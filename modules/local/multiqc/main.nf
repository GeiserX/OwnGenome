/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MULTIQC — Aggregate QC reports into a single HTML dashboard
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MultiQC auto-discovers supported tool outputs (fastp JSON, mosdepth,
    samtools flagstat/stats, etc.) and produces a combined interactive HTML report.
    MULTIQC          one sample's QC files, published in the sample's own
                     folder (<id>/multiqc, where step 28 writes its report),
                     so a run of another sample into the same folder never
                     replaces it
    MULTIQC_COHORT   every sample's, when the run has more than one sample
                     with QC files, published in multiqc/ beside the sample
                     folders; a run of one sample writes no cohort report

    Equivalent to: scripts/28-multiqc.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// The command of both processes, and their versions.yml.
def multiqc_script(process_name, container) {
    """
    multiqc \\
        . \\
        -f \\
        --no-version-check \\
        --no-ai \\
        -o . \\
        -n "multiqc_report.html" \\
        --title "OwnGenome QC"

    cat <<-END_VERSIONS > versions.yml
    "${process_name}":
        multiqc: ${container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

def multiqc_stub(process_name, container) {
    """
    mkdir -p multiqc_report_data
    touch multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${process_name}":
        multiqc: ${container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process MULTIQC {
    tag "$meta.id"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/multiqc" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(multiqc_files, stageAs: "?/*")

    output:
    tuple val(meta), path("multiqc_report.html"), emit: report
    path "multiqc_report_data",                   emit: data
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    multiqc_script(task.process, task.container)

    stub:
    multiqc_stub(task.process, task.container)
}

process MULTIQC_COHORT {
    label 'process_low'

    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode

    input:
    path(multiqc_files, stageAs: "?/*")

    output:
    path "multiqc_report.html",  emit: report
    path "multiqc_report_data",  emit: data
    path "versions.yml",         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    multiqc_script(task.process, task.container)

    stub:
    multiqc_stub(task.process, task.container)
}
