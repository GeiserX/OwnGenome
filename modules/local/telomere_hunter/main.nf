/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TELOMERE_HUNTER — relative telomere content from a WGS BAM (TelomereHunter2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tel_content: intratelomeric reads per million reads with 48-52% GC. A
    higher value means more telomeric sequence. It is a relative research
    measure for comparing samples from one platform, not an age.

    WARNING: Reads the entire BAM (~30-40 GB). Takes 30-60 minutes.

    Equivalent to: scripts/10-telomere-hunter.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process TELOMERE_HUNTER {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/telomere" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(cytoband)  // UCSC GRCh38 chromosome bands: --cytoband, or the copy in assets/reference

    output:
    tuple val(meta), path("${meta.id}"), emit: telomere_results
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The bands are always given: workflows/bam_analysis.nf passes the copy in
    // assets/reference when --cytoband is not set.
    // ulimit: TelomereHunter2 writes one BAM per contig and merges them at
    // once; a BAM with decoy contigs needs more than the 1,024 open files a
    // container gets by default. Raising the soft limit to the hard one needs
    // no privilege, under Docker or Singularity.
    // --plotNone: nothing in the pipeline reads its plots.
    // It has no --version: its own version is the "TelomereHunter2 <v>" line
    // of --help.
    """
    ulimit -n "\$(ulimit -Hn)" 2>/dev/null || true
    telomerehunter2 \\
        -ibt ${bam} \\
        -o ./ \\
        -p ${meta.id} \\
        -c ${task.cpus} \\
        --plotNone \\
        -b ${cytoband}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomerehunter: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
        telomerehunter_reported: \$( { telomerehunter2 --help 2>&1 || true; } | awk '\$1 == "TelomereHunter2" && \$2 ~ /^[0-9]/ { v = \$2; exit } END { print (v != "" ? v : "unknown") }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}
    touch ${meta.id}/${meta.id}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomerehunter: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
