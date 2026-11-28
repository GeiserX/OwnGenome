/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REPORTING — HTML Report Generation & MultiQC Aggregation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Produces two kinds of report:
    1. Per-sample consolidated HTML report, rendered by og report render
       as the bash report is (QC, sample identity and contamination, variant
       counts, ClinVar, PGx, CPIC, CPSR, clinical filter, slivar, ROH, mito
       haplogroup)
    2. MultiQC dashboards of the QC outputs (fastp, mosdepth, etc.): one per
       sample, in its folder, and one across the samples when there are two or more

    Both modules are gated on params.tools containing their tool name.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HTML_REPORT } from '../modules/local/html_report/main'
include { MULTIQC        } from '../modules/local/multiqc/main'
include { MULTIQC_COHORT } from '../modules/local/multiqc/main'

workflow REPORTING {

    take:
    ch_report_inputs  // channel: [meta, vcf, [output files of the steps that ran for the sample]]
    ch_multiqc_files  // channel: [meta, QC file] (fastp, mosdepth, samtools, etc.)

    main:
    ch_versions     = Channel.empty()
    ch_html_reports = Channel.empty()
    ch_multiqc_html = Channel.empty()

    //
    // MODULE 1: Per-sample consolidated HTML report
    //
    if (params.tools && params.tools.split(',').collect{ it.trim() }.contains('html_report')) {
        HTML_REPORT(ch_report_inputs)
        ch_html_reports = HTML_REPORT.out.html_report
        ch_versions     = ch_versions.mix(HTML_REPORT.out.versions)
    }

    //
    // MODULE 2: MultiQC, one report per sample, and one across the samples
    // when the run has more than one
    //
    if (params.tools && params.tools.split(',').collect{ it.trim() }.contains('multiqc')) {
        // Guard: only run MultiQC when QC inputs exist. VCF-only runs produce
        // no mosdepth summaries and MultiQC would fail with "No analysis
        // results found" and create no output files. The skip is logged, so
        // a run with multiqc selected never ends without a word about it.
        ch_multiqc_files
            .map { meta, f -> [meta.id, f] }
            .collect(flat: false)
            .ifEmpty([])
            .map { rows ->
                if (rows.isEmpty()) {
                    log.warn "multiqc skipped: it reads mosdepth's coverage summaries, and this run produced none " +
                             "(mosdepth needs a BAM in the samplesheet and 'mosdepth' in --tools)."
                }
                rows
            }
            .filter { rows -> !rows.isEmpty() }
            .set { ch_multiqc_gated }

        // Each sample's own report, in its folder: a later run of another
        // sample into the same results folder leaves it alone.
        MULTIQC(
            ch_multiqc_gated
                .flatMap { rows -> rows.groupBy { row -> row[0] }.collect { id, rs -> [[id: id], rs.collect { row -> row[1] }.flatten()] } }
        )
        // The cohort report, only with two samples or more: a run of one
        // sample writes none, so it never replaces a cohort's.
        MULTIQC_COHORT(
            ch_multiqc_gated
                .filter { rows -> rows.collect { row -> row[0] }.unique().size() > 1 }
                .map { rows -> rows.collect { row -> row[1] }.flatten() }
        )
        ch_multiqc_html = MULTIQC.out.report
        ch_versions     = ch_versions.mix(MULTIQC.out.versions, MULTIQC_COHORT.out.versions)
    }

    emit:
    html_reports = ch_html_reports  // channel: [meta, html]
    multiqc_html = ch_multiqc_html  // channel: [meta, html], one per sample
    versions     = ch_versions      // channel: path(versions.yml)
}
