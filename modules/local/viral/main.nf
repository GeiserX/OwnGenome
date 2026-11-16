/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VIRAL — Viral DNA in the reads that do not map to the human genome (step 47)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Opt-in ('viral' in --tools, with --viral_ref: setup.sh --viral-data). Two
    processes, as in scripts/47-viral.sh:
      VIRAL_READS   the unmapped reads and the chrEBV reads of the BAM, realigned
                    with minimap2 to the viral genomes; samtools coverage of that
                    BAM and of the sample's BAM (the human depth), samtools
                    depth -aa over every viral genome, the viral alignments
                    (samtools view -F 4) and samtools idxstats of the sample's BAM
      VIRAL_REPORT  og viral report on those files: one row per viral
                    genome, counting only unique, mappable reads, with
                    inherited integrated HHV-6 called by Liu et al. 2020's
                    depth threshold

    Equivalent to: scripts/47-viral.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process VIRAL_READS {
    tag "$meta.id"
    label 'process_medium'

    // The depth, the alignments and the idxstats are og viral report's inputs only.
    publishDir { "${params.outdir}/${meta.id}/viral" }, mode: params.publish_dir_mode,
        pattern: "*_{viral_reads.bam,viral_reads.bam.bai,viral_coverage.tsv,human_coverage.tsv}"

    input:
    tuple val(meta), path(bam), path(bai)
    path(viral_ref)   // reference/viral-<version>/viral_genomes.fa

    output:
    tuple val(meta), path("${meta.id}_viral_coverage.tsv"), path("${meta.id}_viral_depth.tsv"),
          path("${meta.id}_human_coverage.tsv"), path("${meta.id}_viral_hits.tsv"),
          path("${meta.id}_bam_idxstats.tsv"), emit: coverage
    tuple val(meta), path("${meta.id}_viral_reads.bam"), path("${meta.id}_viral_reads.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The same commands as scripts/47-viral.sh: the unmapped reads (a read
    // unmapped beside a mapped mate too), the mapped reads on chrEBV, each
    // read once (samtools fastq leaves out secondary and supplementary
    // records); both coverages count duplicates, so the two depths match.
    // idxstats counts the unplaced unmapped reads; the viral alignments let
    // og viral report count only unique, mappable reads.
    def p = meta.id
    """
    samtools view -@ ${task.cpus} -u -f 4 ${bam} | samtools fastq - > ${p}_unmapped.fq
    if samtools view -H ${bam} | awk -F '\\t' '\$1 == "@SQ" && \$2 == "SN:chrEBV" {f = 1} END {exit !f}'; then
        samtools view -@ ${task.cpus} -u -F 4 ${bam} chrEBV | samtools fastq - >> ${p}_unmapped.fq
    fi
    samtools coverage --ff UNMAP,SECONDARY,QCFAIL ${bam} > ${p}_human_coverage.tsv
    samtools idxstats ${bam} > ${p}_bam_idxstats.tsv
    minimap2 -ax sr --frag=no --secondary=no -t ${task.cpus} -o ${p}_viral.sam ${viral_ref} ${p}_unmapped.fq
    samtools sort -o ${p}_viral_reads.bam ${p}_viral.sam
    samtools index ${p}_viral_reads.bam
    samtools coverage --ff UNMAP,SECONDARY,QCFAIL ${p}_viral_reads.bam > ${p}_viral_coverage.tsv
    samtools depth -aa -G UNMAP,SECONDARY,QCFAIL ${p}_viral_reads.bam > ${p}_viral_depth.tsv
    samtools view -F 4 ${p}_viral_reads.bam > ${p}_viral_hits.tsv
    rm -f ${p}_unmapped.fq ${p}_viral.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2_samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    def p = meta.id
    """
    printf '#rname\\tstartpos\\tendpos\\tnumreads\\tcovbases\\tcoverage\\tmeandepth\\tmeanbaseq\\tmeanmapq\\n' > ${p}_viral_coverage.tsv
    cp ${p}_viral_coverage.tsv ${p}_human_coverage.tsv
    touch ${p}_viral_depth.tsv ${p}_viral_reads.bam ${p}_viral_reads.bam.bai ${p}_viral_hits.tsv
    printf '*\\t0\\t0\\t0\\n' > ${p}_bam_idxstats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2_samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process VIRAL_REPORT {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/viral" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(viral_coverage), path(viral_depth), path(human_coverage), path(viral_hits),
          path(idxstats)
    path(viral_ref)

    output:
    tuple val(meta), path("${meta.id}_viral.tsv"), emit: table
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    og viral report --sample ${meta.id} --viral-fasta ${viral_ref} --viral-coverage ${viral_coverage} \\
        --viral-depth ${viral_depth} --human-coverage ${human_coverage} --viral-hits ${viral_hits} \\
        --idxstats ${idxstats} --out ${meta.id}_viral.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf 'virus\\taccession\\tlength\\treads\\tbreadth\\tmean_depth\\tmode_depth\\thuman_depth\\trelative_depth\\tcopies_per_cell\\tcall\\tnote\\nstub\\t.\\t0\\t0\\t0\\t0\\t0\\t.\\t.\\t.\\tnot found\\tstub\\n' > ${meta.id}_viral.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
