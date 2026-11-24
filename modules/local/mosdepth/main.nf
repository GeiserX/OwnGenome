/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MOSDEPTH — Fast per-base and per-region coverage statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Computes precise per-base depth from actual alignments (not just the index).
    Outputs coverage distributions, threshold summaries, and region BED files:
    500 bp windows for a whole genome, or, with --capture_bed (--model_type
    WES), the capture regions, so an exome's summary gives its on-target
    depth (the script's CAPTURE_BED).

    A second pass over the genes of assets/callability_genes.tsv counts, per
    gene on chr1-22, X and Y, the bases with 10 or more reads of MAPQ 20 or
    more: ${id}_callability.bed.gz, which the report reads to list the genes
    this genome could not read reliably (og summary collect). When the
    tables also hold assets/callability_exons.tsv, the coding exons of those
    genes ride in the same pass into ${id}_callability_exons.bed.gz, so a
    gene with an unread coding exon is not passed on its whole body.

    Equivalent to: scripts/16b-mosdepth.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MOSDEPTH {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/coverage" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(callability_tables)  // assets/callability_genes.tsv, and assets/callability_exons.tsv when exons are measured
    path(capture_bed)         // --capture_bed with --model_type WES, or []

    output:
    tuple val(meta), path("*.mosdepth.summary.txt"),     emit: summary
    tuple val(meta), path("*.mosdepth.global.dist.txt"), emit: global_dist
    tuple val(meta), path("*.regions.bed.gz"),           emit: regions_bed, optional: true
    tuple val(meta), path("*.thresholds.bed.gz"),        emit: thresholds_bed, optional: true
    tuple val(meta), path("*_callability*.bed.gz"),      emit: callability   // the gene table, and the exon table when exons were measured
    path "versions.yml",                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def by     = capture_bed ?: 500
    """
    mosdepth \\
        --by ${by} \\
        --fast-mode \\
        --no-per-base \\
        --threads ${task.cpus} \\
        --thresholds 1,5,10,15,20,30,50 \\
        ${prefix} \\
        ${bam}

    # Gene callability, as scripts/16b-mosdepth.sh: the genes with a position
    # (not the ALT-only ones) and the coding exons of those genes in the exon
    # table, in one pass; reads of MAPQ 20 or more, bases at 10x or more. The
    # rows named GENE:exonN go to the exon table, the others to the gene table.
    # It writes into callability/, so its own summary and region files never
    # match the globs of the outputs above. No --fast-mode: without it the
    # bases where the two mates of a pair overlap count once and the bases of a
    # CIGAR deletion do not count, so a gene is read only when 10 fragments
    # cover it.
    mkdir -p callability
    genes=""; exons=""
    for f in ${callability_tables}; do
        case "\$f" in
            *callability_exons.tsv) exons=\$f ;;
            *) genes=\$f ;;
        esac
    done
    awk -F'\\t' -v OFS='\\t' '!/^#/ && \$1 != "gene" && \$4 ~ /^[0-9]+\$/ {print \$3, \$4, \$5, \$1}' \
        "\$genes" > callability/genes.bed
    : > callability/exons.bed
    if [ -n "\$exons" ]; then
        awk -F'\\t' -v OFS='\\t' 'NR == FNR {g[\$4] = 1; next}
            !/^#/ && \$1 != "exon" && (\$2 in g) {print \$4, \$5, \$6, \$1}' callability/genes.bed "\$exons" > callability/exons.bed
    fi
    cat callability/genes.bed callability/exons.bed > callability/regions.bed
    mosdepth \\
        --by callability/regions.bed \\
        --mapq 20 \\
        --no-per-base \\
        --threads ${task.cpus} \\
        --thresholds 10 \\
        callability/c \\
        ${bam}
    if [ -s callability/exons.bed ]; then
        gzip -dc callability/c.thresholds.bed.gz | awk -F'\\t' \\
            'NR == 1 {print > "callability/g"; print > "callability/e"; next}
             \$4 ~ /:exon[0-9]+\$/ {print > "callability/e"; next} {print > "callability/g"}'
        gzip -nc callability/g > ${prefix}_callability.bed.gz
        gzip -nc callability/e > ${prefix}_callability_exons.bed.gz
    else
        mv callability/c.thresholds.bed.gz ${prefix}_callability.bed.gz
    fi
    rm -rf callability

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    // The stub prints the windows the script passes (tests/e2e/nf-exome-1-stub.sh reads them).
    def prefix = task.ext.prefix ?: "${meta.id}"
    def by     = capture_bed ?: 500
    """
    echo "mosdepth options: --by ${by}"
    touch ${prefix}.mosdepth.summary.txt
    touch ${prefix}.mosdepth.global.dist.txt
    echo | gzip > ${prefix}_callability.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
