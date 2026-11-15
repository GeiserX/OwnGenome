/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BLOODTYPE — Blood group antigens (red cell, platelet, neutrophil) predicted from the genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Opt-in ('bloodtype' in --tools). Three processes, as in scripts/41-bloodtype.sh:
      BLOODTYPE_REGIONS og cn bed --locus RHD, once: RHD, its
                       exons and flanks (the RHD row of bin/cn_loci.tsv)
      BLOODTYPE_DEPTH  samtools bedcov -Q 20 -j over the sites of
                       assets/bloodtype_markers.tsv and those regions, for each
                       sample with a BAM: the read depth at every marker site
                       and over RHD and its flanks
      BLOODTYPE        the VCF's records within 100 bp of a marker site,
                       left-aligned against the reference (bcftools norm -c w,
                       multi-allelic records kept whole), so an indel a caller
                       wrote further right in a repeat meets its marker; then
                       og bloodtype report on that file and the depth, reading
                       the marker sites through its index; without a BAM
                       ([] in its place) RhD is not determined and sites with
                       no VCF record are read as reference (confidence medium)

    Equivalent to: scripts/41-bloodtype.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BLOODTYPE_REGIONS {
    label 'process_single'

    output:
    path "rhd_regions.bed", emit: bed
    path "versions.yml",    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    og cn bed --locus RHD > rhd_regions.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch rhd_regions.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process BLOODTYPE_DEPTH {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/bloodtype" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(markers)   // assets/bloodtype_markers.tsv
    path(rhd_bed)   // BLOODTYPE_REGIONS' BED

    output:
    tuple val(meta), path("${meta.id}_bloodtype_depth.tsv"), emit: depth
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The marker table is a BED with extra columns; bedcov reads it once the
    // comment and header lines are gone; the RHD regions follow it. bedcov stops on a contig the BAM has
    // no @SQ line for, so those rows are left out (their sites are then not
    // determined). The same commands as scripts/41-bloodtype.sh.
    """
    samtools view -H ${bam} | awk -F '\\t' '\$1 == "@SQ" {for (i = 2; i <= NF; i++) if (\$i ~ /^SN:/) print substr(\$i, 4)}' > contigs.txt
    { grep -v '^#' ${markers} | tail -n +2; cat ${rhd_bed}; } | awk -F '\\t' 'NR == FNR {c[\$1] = 1; next} \$1 in c' contigs.txt - > targets.bed
    if [ -s targets.bed ]; then samtools bedcov -Q 20 -j targets.bed ${bam} > ${meta.id}_bloodtype_depth.tsv; else : > ${meta.id}_bloodtype_depth.tsv; fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_bloodtype_depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process BLOODTYPE {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/bloodtype" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index), path(depth)   // vcf_index: its .tbi or .csi, or []; depth: BLOODTYPE_DEPTH's, or [] without a BAM
    path(markers)
    path(reference)
    path(reference_fai)

    output:
    tuple val(meta), path("${meta.id}_bloodtype.tsv"),          emit: table
    tuple val(meta), path("${meta.id}_bloodtype_evidence.tsv"), emit: evidence
    path "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def depth_arg = depth ? "--depth ${depth}" : ''
    // og bloodtype report matches a marker by position, REF and ALT, so an
    // indel must sit where the table writes it: left-aligned. The records
    // within 100 bp of a marker site (more than the 38 bp the longest marker
    // indel can shift along its repeat) are left-aligned first. A
    // multi-allelic record stays whole; -c w keeps a record whose REF differs
    // from the reference as it is. Records on a contig the reference lacks go
    // through unchanged. The VCF is read through its index when it has a
    // current one, else once from the start. The same commands as
    // scripts/41-bloodtype.sh; og runs from bin/, on the task's PATH.
    """
    grep -v '^#' ${markers} | tail -n +2 \\
        | awk -F '\\t' -v OFS='\\t' '{s = \$2 - 100; if (s < 0) s = 0; print \$1, s, \$3 + 100; c = \$1; if (sub(/^chr/, "", c)) print c, s, \$3 + 100}' \\
        | sort -k1,1 -k2,2n \\
        | awk -F '\\t' -v OFS='\\t' '\$1 == c && \$2 <= e {if (\$3 > e) e = \$3; next} {if (c != "") print c, s, e; c = \$1; s = \$2; e = \$3} END {if (c != "") print c, s, e}' > bt_windows.bed
    cut -f1 ${reference_fai} > bt_ref_contigs.txt
    awk -F '\\t' 'NR == FNR {c[\$1] = 1; next} \$1 in c' bt_ref_contigs.txt bt_windows.bed > bt_on_ref.bed
    awk -F '\\t' 'NR == FNR {c[\$1] = 1; next} !(\$1 in c)' bt_ref_contigs.txt bt_windows.bed > bt_off_ref.bed
    SEL=-T OVERLAP=--targets-overlap
    for x in tbi csi; do
        if [ -f "${vcf}.\$x" ] && [ ! "${vcf}" -nt "${vcf}.\$x" ]; then SEL=-R OVERLAP=--regions-overlap; fi
    done
    pick() { if [ -s "\$1" ]; then bcftools view "\$SEL" "\$1" "\$OVERLAP" record -Ou ${vcf}; else bcftools view -h -Ou ${vcf}; fi; }
    pick bt_on_ref.bed | bcftools norm -c w -f ${reference} -Ob -o bt_on_ref.bcf
    pick bt_off_ref.bed | bcftools view -Ob -o bt_off_ref.bcf
    bcftools concat -Ou bt_on_ref.bcf bt_off_ref.bcf | bcftools sort -T ./bt_sort -Oz -o bt_norm.vcf.gz
    bcftools index -f -t bt_norm.vcf.gz

    og bloodtype report --sample ${meta.id} --vcf bt_norm.vcf.gz --markers ${markers} ${depth_arg} \\
        --out ${meta.id}_bloodtype.tsv --evidence ${meta.id}_bloodtype_evidence.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf 'system\\tphenotype\\tgenotype\\tconfidence\\tnote\\nABO\\tnot determined\\t.\\tnot determined\\tstub\\n' > ${meta.id}_bloodtype.tsv
    touch ${meta.id}_bloodtype_evidence.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
