/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MITO_VARIANTS — Mitochondrial variant calling with heteroplasmy detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Uses GATK Mutect2 in mitochondrial mode to call variants on chrM,
    including low-frequency heteroplasmic variants (AF < 0.95).

    Four steps:
    1. Extract chrM reads from the BAM (GATK PrintReads; step 20 uses
       samtools view for the same reads)
    2. Run Mutect2 --mitochondria-mode, then the shifted-reference pass:
       the control region (m.16025-16569 and m.1-575) again on chrM rotated
       by 8000 bp (og setup chrm-shift; GATK's BwaSpark realigns the reads,
       as the image has no bwa binary), mapped back and merged with the
       linear calls outside it (og vcf shift-back). Reads across the origin
       align in full there. A failure leaves the linear calls alone.
    3. Filter variants with FilterMutectCalls --mitochondria-mode
    4. Mark possible NuMTs (nuclear copies of chrM) with NuMTFilterTool at
       the median autosomal depth, read from MOSDEPTH's summary and
       global distribution with the awk of scripts/20-mtoolbox.sh. Without
       mosdepth in --tools the depth is 0 and the filter marks nothing, as
       step 20 does without step 16b's output: the task warns, and
       MITO_VARIANTS_ANNOTATE writes the '# numt' line that says so.

    MITO_VARIANTS_ANNOTATE, after it, annotates each call with its ClinVar
    record and gnomAD v3.1's chrM sites (og mito annotate), as step 20's
    last stage does.

    Equivalent to: scripts/20-mtoolbox.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MITO_VARIANTS {
    tag "$meta.id"
    label 'process_medium'
    // Root under Docker: the shifted-reference pass's BwaSpark starts
    // Hadoop's login, which needs a passwd entry for the task's uid, and the
    // GATK image has none for the -u of docker.runOptions (a LoginException,
    // "invalid null input: name"; HADOOP_USER_NAME and SPARK_USER do not
    // help). Step 20 runs that container as root too. A later --user wins.
    // docker.fixOwnership gives the work directory back only after a task
    // that succeeds, so the script gives it back itself on any exit (below),
    // to the directory's owner, also with docker.fixOwnership = false.
    // A process's containerOptions replaces the run's process.containerOptions
    // (the docker profile's '--network none', or any extra option of yours),
    // so '--network none' is repeated here unless the run's docker.runOptions
    // already sets a network (Docker refuses the flag twice). These are
    // Docker's options, for Docker only: Singularity and Apptainer run as you
    // with your passwd entry and need none, and their exec has no --user.
    containerOptions {
        def run_options = ((nextflow.Global.session?.config?.docker as Map)?.runOptions ?: '').toString()
        workflow.containerEngine != 'docker' ? '' :
            '--user 0:0 -e HOME=/root' + (run_options.contains('--net') ? '' : ' --network none')
    }

    publishDir { "${params.outdir}/${meta.id}/mito" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai), path(mosdepth_summary), path(mosdepth_dist)  // [] [] without mosdepth
    path(reference)
    path(reference_fai)
    path(reference_dict)

    output:
    tuple val(meta), path("*_chrM_filtered.vcf.gz"),     emit: mito_vcf
    tuple val(meta), path("*_chrM_filtered.vcf.gz.tbi"), emit: mito_vcf_index, optional: true
    tuple val(meta), path("*_chrM_mutect2.vcf.gz.stats"), emit: stats
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // With mosdepth in --tools the workflow joins its output in; a join that
    // lost it would run the NuMT filter at depth 0 and mark nothing, silently.
    if ((params.tools ?: '').toString().split(',').collect { it.trim() }.contains('mosdepth') && !mosdepth_summary) {
        error "MITO_VARIANTS got no mosdepth output for '${meta.id}' although mosdepth is in --tools: the NuMT filter would run at depth 0"
    }
    def depth_files = mosdepth_summary ? "${mosdepth_summary} ${mosdepth_dist}" : ''
    """
    # Root under Docker (above): give every file back on any exit, also when
    # the task fails, to the owner of the work directory, who made it (you),
    # whatever docker.fixOwnership says.
    OWNER=\$(stat -c %u:%g .)
    trap 'chown -R "\$OWNER" . 2>/dev/null || true' EXIT

    # Step 1: Extract chrM reads
    gatk PrintReads \\
        -I ${bam} \\
        -L chrM \\
        -O ${prefix}_chrM.bam

    # Step 2: Run Mutect2 in mitochondrial mode
    gatk Mutect2 \\
        -R ${reference} \\
        -I ${prefix}_chrM.bam \\
        -L chrM \\
        --mitochondria-mode \\
        --max-mnp-distance 0 \\
        -O ${prefix}_chrM_mutect2.vcf.gz

    # Step 2b: The shifted-reference pass (scripts/20-mtoolbox.sh): the control
    # region again on chrM rotated by 8000 bp, its calls mapped back and merged
    # with the linear ones outside it. A failure leaves the linear calls alone.
    shifted_pass() {
        rm -rf shifted && mkdir shifted || return 1
        # Under Docker the task runs as root (above), and OWNER names you:
        # shifted/ is yours from the start and the large files go to the
        # container's own /tmp, so a task killed outright (no EXIT trap)
        # leaves in shifted/ only plain files, which you can remove.
        # Elsewhere the task runs as you, and they stay in shifted/.
        if [ "\$(id -u)" = 0 ]; then
            chown "\$OWNER" shifted || return 1
            T=\$(mktemp -d) || return 1
        else
            T=shifted/tmp && mkdir "\$T" || return 1
        fi
        og setup chrm-shift --fasta ${reference} \\
            --out shifted/chrM_shifted.fasta --intervals shifted/control.intervals || return 1
        export SPARK_LOCAL_IP=127.0.0.1 SPARK_LOCAL_HOSTNAME=localhost
        gatk CreateSequenceDictionary -R shifted/chrM_shifted.fasta -O shifted/chrM_shifted.dict || return 1
        gatk BwaMemIndexImageCreator -I shifted/chrM_shifted.fasta -O shifted/chrM_shifted.img || return 1
        gatk PrintReads -I ${prefix}_chrM.bam -L chrM -O "\$T/reads.bam" \\
            --read-filter NotDuplicateReadFilter --read-filter NotSecondaryAlignmentReadFilter \\
            --read-filter NotSupplementaryAlignmentReadFilter \\
            --read-filter MateOnSameContigOrNoMappedMateReadFilter || return 1
        # Not --SANITIZE: it re-encodes qualities that are all Q31 or more as
        # if they were phred+64. BwaSpark aligns records two by two as pairs,
        # so awk keeps both mates of a name or neither.
        gatk RevertSam -I "\$T/reads.bam" -O "\$T/unaligned.sam" \\
            --SORT_ORDER queryname --OUTPUT_BY_READGROUP false --RESTORE_ORIGINAL_QUALITIES false \\
            --VALIDATION_STRINGENCY LENIENT --ATTRIBUTE_TO_CLEAR FT --ATTRIBUTE_TO_CLEAR CO || return 1
        awk -F'\\t' '/^@/ { print; next } \$1 == name { print prev; print; name = ""; next } { name = \$1; prev = \$0 }' \\
            "\$T/unaligned.sam" > "\$T/paired.sam" || return 1
        # Unpaired reads (single-end, HiFi or ONT) leave no pair, and the
        # control region would lose its calls: under 90% kept fails the pass.
        n=\$(grep -vc '^@' "\$T/unaligned.sam" || true)
        k=\$(grep -vc '^@' "\$T/paired.sam" || true)
        echo "shifted-reference pass: \${k} of \${n} chrM reads kept as complete pairs" >&2
        if [ "\$k" -eq 0 ] || [ \$((k * 10)) -lt \$((n * 9)) ]; then
            echo "the chrM reads are not paired (under 90% kept): the rotated pass needs read pairs" >&2
            return 1
        fi
        gatk --java-options -Xmx4g BwaSpark -I "\$T/paired.sam" -R shifted/chrM_shifted.fasta \\
            --bwa-mem-index-image shifted/chrM_shifted.img --spark-master "local[${task.cpus}]" \\
            -O "\$T/aligned.bam" || return 1
        gatk SortSam -I "\$T/aligned.bam" -O "\$T/sorted.bam" --SORT_ORDER coordinate \\
            --CREATE_INDEX true || return 1
        gatk Mutect2 -R shifted/chrM_shifted.fasta -I "\$T/sorted.bam" -L shifted/control.intervals \\
            --mitochondria-mode --max-mnp-distance 0 -O shifted/shifted_mutect2.vcf.gz || return 1
        og vcf shift-back --linear ${prefix}_chrM_mutect2.vcf.gz \\
            --shifted shifted/shifted_mutect2.vcf.gz --out shifted/merged_mutect2.vcf.gz || return 1
    }
    T=
    FILTER_IN=${prefix}_chrM_mutect2.vcf.gz
    if shifted_pass; then
        FILTER_IN=shifted/merged_mutect2.vcf.gz
    else
        echo "WARNING: the shifted-reference pass failed: the control region keeps the linear calls alone" >&2
    fi
    [ -z "\$T" ] || rm -rf "\$T"
    rm -f shifted/*.img

    # Step 3: Filter variants (the stats are the linear pass's, over all of chrM)
    gatk FilterMutectCalls \\
        -R ${reference} \\
        -V \${FILTER_IN} \\
        --stats ${prefix}_chrM_mutect2.vcf.gz.stats \\
        --mitochondria-mode \\
        -O ${prefix}_chrM_mutect2_filtered.vcf.gz

    # Step 4: Mark possible NuMTs. The median autosomal depth is the depth at
    # which half of the chr1-22 bases are covered at least that deep, from
    # mosdepth's per-chromosome distribution (scripts/20-mtoolbox.sh, step 5).
    DEPTH_FILES="${depth_files}"
    if [ -n "\${DEPTH_FILES}" ]; then
        AUTOSOMAL_COVERAGE=\$(awk -F'\\t' '
            FNR == NR { if (\$1 ~ /^chr[0-9]+\$/) len[\$1] = \$2; next }
            (\$1 in len) { at_least[\$2] += len[\$1] * \$3 }
            END {
              for (c in len) total += len[c]
              best = 0
              if (total > 0) for (k in at_least) if (at_least[k] / total >= 0.5 && k + 0 > best) best = k + 0
              print best
            }' \${DEPTH_FILES})
        echo "Median autosomal coverage: \${AUTOSOMAL_COVERAGE} (from \${DEPTH_FILES})"
    else
        AUTOSOMAL_COVERAGE=0
        echo "WARNING: mosdepth is not in --tools: NuMTFilterTool runs at the default depth 0 and marks nothing" >&2
    fi
    gatk NuMTFilterTool \\
        -R ${reference} \\
        -V ${prefix}_chrM_mutect2_filtered.vcf.gz \\
        --autosomal-coverage \${AUTOSOMAL_COVERAGE} \\
        -O ${prefix}_chrM_filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // With mosdepth in --tools the workflow joins its output in; a join that
    // lost it would run the NuMT filter at depth 0 and mark nothing, silently.
    if ((params.tools ?: '').toString().split(',').collect { it.trim() }.contains('mosdepth') && !mosdepth_summary) {
        error "MITO_VARIANTS got no mosdepth output for '${meta.id}' although mosdepth is in --tools: the NuMT filter would run at depth 0"
    }
    """
    touch ${prefix}_chrM_filtered.vcf.gz
    touch ${prefix}_chrM_mutect2.vcf.gz.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}


/*
    MITO_VARIANTS_ANNOTATE: each ALT of the filtered chrM calls with its
    ClinVar record (CLNSIG, CLNREVSTAT, CLNDN of the same allele, from
    --mito_clinvar, the whole ClinVar VCF with chr names and its .tbi beside
    it) and gnomAD v3.1's chrM sites (AF_hom, AF_het, max_hl, from
    --gnomad_mito, the file GNOMAD_MITO_URL in versions.env names), into
    <id>_chrM_annotated.tsv. A file not given is named in the table's header
    and its columns are '.', as step 20 does without it. Without mosdepth in
    --tools MITO_VARIANTS had no autosomal depth for the NuMT filter, so the
    table gets the '# numt' line (--numt-default), as step 20 writes it.

    Equivalent to: stage 6 of scripts/20-mtoolbox.sh
*/
process MITO_VARIANTS_ANNOTATE {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/mito" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf)
    path(clinvar)        // [] without --mito_clinvar
    path(clinvar_index)  // staged beside it
    path(gnomad)         // [] without --gnomad_mito

    output:
    tuple val(meta), path("${meta.id}_chrM_annotated.tsv"), emit: annotated
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def cv = clinvar ? "--clinvar ${clinvar}" : ''
    def gn = gnomad ? "--gnomad ${gnomad}" : ''
    // MITO_VARIANTS gets mosdepth's output exactly when mosdepth is in --tools.
    def numt = (params.tools ?: '').toString().split(',').collect { it.trim() }.contains('mosdepth') ? '' : '--numt-default'
    """
    og mito annotate --vcf ${vcf} ${cv} ${gn} ${numt} --out ${meta.id}_chrM_annotated.tsv

    og_tree=\$(og --tree-hash)
    printf '"%s":\\n    ubuntu: %s\\n    owngenome: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" "\$og_tree" > versions.yml
    """

    stub:
    """
    printf 'chrom\\tpos\\tref\\talt\\tfilter\\taf\\tclinvar_id\\tclinvar_significance\\tclinvar_review\\tclinvar_condition\\tgnomad_af_hom\\tgnomad_af_het\\tgnomad_max_hl\\n' > ${meta.id}_chrM_annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ubuntu: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
        owngenome: stub
    END_VERSIONS
    """
}


/*
    MITO_DELETIONS: large mtDNA deletions and the mtDNA copy number of one
    BAM, in the eKLIPse image (it holds samtools and Python too).
    1. samtools idxstats of the whole BAM: chrM reads per base over the
       autosomal reads per base, times 2, is the mtDNA copy number per cell
       (og mito copies, into <id>_mtdna_copies.tsv;
       not assessed under 20 reads per kb of autosome, as on a sliced BAM).
    2. The chrM reads into their own BAM, and eKLIPse on it against the rCRS
       GenBank file its image ships (data/NC_012920.1.gb beside eKLIPse.py;
       GRCh38's chrM is the rCRS). eKLIPse finds deletions from soft-clipped
       reads whose clipped part maps across the deletion
       (og mito deletions turns its CSV into
       <id>_chrM_deletions.tsv). The report names the deletions of
       assets/mito_confirmed.tsv among them.

    Equivalent to: stages 7 and 8 of scripts/20-mtoolbox.sh
*/
process MITO_DELETIONS {
    tag "$meta.id"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/mito" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_chrM_deletions.tsv"), emit: deletions, optional: true
    tuple val(meta), path("${meta.id}_mtdna_copies.tsv"),   emit: copies
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def id = meta.id
    """
    samtools idxstats ${bam} > ${id}.idxstats
    og mito copies --idxstats ${id}.idxstats --out ${id}_mtdna_copies.tsv

    samtools view -b -o ${id}_chrM.bam ${bam} chrM
    samtools index ${id}_chrM.bam
    printf '%s\\t%s\\n' "\$PWD/${id}_chrM.bam" "${id}" > eklipse_input.tsv
    mkdir -p eklipse_out eklipse_tmp
    EK=\$(command -v eKLIPse.py)
    # An eKLIPse error, or a run with no CSV, leaves no deletion table but
    # keeps the copy number, as stage 8 of scripts/20-mtoolbox.sh does; the
    # report then says the deletions were not looked for.
    # A failed run can leave a partial CSV, so it is not read.
    CSV=
    if python3 "\$EK" -in eklipse_input.tsv -ref "\$(dirname "\$EK")/data/NC_012920.1.gb" \\
        -out eklipse_out -tmp eklipse_tmp -thread ${task.cpus}; then
        CSV=\$(find eklipse_out -name eKLIPse_deletions.csv | head -n 1)
    else
        echo "WARNING: eKLIPse failed: no deletion table" >&2
    fi
    if [ -z "\$CSV" ]; then
        echo "WARNING: eKLIPse wrote no eKLIPse_deletions.csv: no deletion table" >&2
    else
        og mito deletions --eklipse "\$CSV" --out ${id}_chrM_deletions.tsv
    fi

    printf '"%s":\\n    eklipse: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    printf 'bkp5\\tbkp3\\tlength\\tlevel_pct\\tlevel_fwd_pct\\tlevel_rev_pct\\tblast_5\\tblast_3\\tdepth_5\\tdepth_3\\trepeat\\n' > ${meta.id}_chrM_deletions.tsv
    printf 'metric\\tvalue\\nstatus\\tnot assessed\\n' > ${meta.id}_mtdna_copies.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eklipse: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
