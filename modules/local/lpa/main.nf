/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LPA — the LPA KIV-2 repeat count estimated from read depth (Lp(a))
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Opt-in ('lpa' in --tools). Three processes, as in scripts/46-lpa.sh:
      LPA_REGIONS  bin/copy_number.py bed --locus LPA_KIV2, once: the KIV-2
                   and single-copy LPA exons (the LPA_KIV2 row of bin/cn_loci.tsv),
                   and the seven GRCh38 copies of the KIV-2 4925G>A site
                   (bin/traits_report.py --marker-regions LPA_KIV2_4925)
      LPA_DEPTH    samtools bedcov -j over the exons, all MAPQ and MAPQ >= 1,
                   and samtools mpileup at the 4925G>A sites, for each sample
                   with a BAM; empty tables when the BAM has no chr6
      LPA          bin/copy_number.py call --locus LPA_KIV2: the count, its
                   95% interval and the message (laid out by bin/lpa_kiv2.py),
                   then bin/traits_report.py --marker LPA_KIV2_4925: the 4925G>A
                   reading beside it, <id>_lpa_kiv2_4925.tsv

    Equivalent to: scripts/46-lpa.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process LPA_REGIONS {
    label 'process_single'

    output:
    path "lpa_*",        emit: bed       // lpa_regions.bed and lpa_sites.txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    copy_number.py bed --locus LPA_KIV2 > lpa_regions.bed
    traits_report.py --marker-regions LPA_KIV2_4925 > lpa_sites.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch lpa_regions.bed lpa_sites.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process LPA_DEPTH {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/lpa" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(regions)   // LPA_REGIONS' lpa_regions.bed and lpa_sites.txt, read by name

    output:
    tuple val(meta), path("${meta.id}_lpa_kiv2_depth.tsv"), path("${meta.id}_lpa_kiv2_depth_q1.tsv"),
        path("${meta.id}_lpa_kiv2_4925_pileup.tsv"), emit: depth
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // bedcov stops on a contig the BAM has no @SQ line for: without chr6 the
    // tables stay empty. A failed mpileup leaves the count to run: the pileup
    // then holds one line that is not mpileup output, so 4925G>A reads 'not
    // assessed'. The same commands as scripts/46-lpa.sh.
    """
    if samtools view -H ${bam} | awk -F '\\t' '\$1 == "@SQ" {for (i = 2; i <= NF; i++) if (\$i == "SN:chr6") f = 1} END {exit !f}'; then
        samtools bedcov -j lpa_regions.bed ${bam} > ${meta.id}_lpa_kiv2_depth.tsv
        samtools bedcov -j -Q 1 lpa_regions.bed ${bam} > ${meta.id}_lpa_kiv2_depth_q1.tsv
        ( while read -r r; do samtools mpileup -A -B -q 0 -Q 13 -d 100000 -r "\$r" ${bam} || exit 1; done < lpa_sites.txt ) > ${meta.id}_lpa_kiv2_4925_pileup.tsv \\
            || echo '# samtools mpileup failed: 4925G>A not assessed' > ${meta.id}_lpa_kiv2_4925_pileup.tsv
    else
        : > ${meta.id}_lpa_kiv2_depth.tsv
        : > ${meta.id}_lpa_kiv2_depth_q1.tsv
        : > ${meta.id}_lpa_kiv2_4925_pileup.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_lpa_kiv2_depth.tsv ${meta.id}_lpa_kiv2_depth_q1.tsv ${meta.id}_lpa_kiv2_4925_pileup.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process LPA {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/lpa" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(depth), path(depth_q1), path(pileup)

    output:
    tuple val(meta), path("${meta.id}_lpa_kiv2.tsv"),      emit: table
    tuple val(meta), path("${meta.id}_lpa_kiv2_4925.tsv"), emit: null_allele  // KIV-2 4925G>A beside the count
    path "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    copy_number.py call --locus LPA_KIV2 --bedcov ${depth} --bedcov-mapq1 ${depth_q1} --table ${meta.id}_lpa_kiv2.tsv
    traits_report.py --marker LPA_KIV2_4925 --pileup ${pileup} --copy-table ${meta.id}_lpa_kiv2.tsv \\
        --marker-out ${meta.id}_lpa_kiv2_4925.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf 'metric\\tvalue\\nstatus\\tnot assessed\\nmessage\\tstub\\n' > ${meta.id}_lpa_kiv2.tsv
    printf 'metric\\tvalue\\nstatus\\tnot assessed\\nmessage\\tstub\\n' > ${meta.id}_lpa_kiv2_4925.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
