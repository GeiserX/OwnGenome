/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CPSR — Cancer Predisposition Sequencing Reporter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Screens germline VCF for cancer predisposition variants using PCGR 2.x.
    Requires the VEP cache of PCGR_VEP_CACHE_RELEASE in versions.env (115 for
    PCGR 2.3.2), separate from the VEP step's VEP_CACHE_RELEASE cache. CPSR 2.3
    has no --classify_all: it classifies every panel variant by itself.

    CPSR stops on a --sample_id outside 3 to 40 characters (PCGR 2.3.2), so
    it gets bin/cpsr_sample_id's id (S1 -> S1_cpsr, a longer id cut to 40),
    and its files are renamed back to <meta.id>.cpsr.*. CPSR's own HTML
    title then shows that id.

    The cDNA read check on CPSR's P/LP rows (og cdna check --cpsr) runs
    after it, in CPSR_CDNA_CHECK (CLINVAR_CDNA_CHECK included under that name,
    source 'cpsr', workflows/clinical.nf), which gets the BAM and the exons
    and writes cpsr/<id>_cpsr_cdna_check.tsv.

    --cpsr_secondary_findings (default true) passes --secondary_findings:
    CPSR then also reports the P/LP variants in the ACMG SF genes. With false
    it reports the cancer genes only, as scripts/17-cpsr.sh does with
    CPSR_SECONDARY_FINDINGS=false (run-all.sh passes that switch on).

    Equivalent to: scripts/17-cpsr.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// CPSR's options after --panel_id 0 for --cpsr_secondary_findings, which
// must be true or false (a params file can give either as a string). The
// stub prints them too, so a stub run shows what the real task would pass.
def cpsrSecondaryArgs(value) {
    def v = value.toString()
    if (!(v in ['true', 'false'])) {
        error "--cpsr_secondary_findings must be true or false, got '${v}'"
    }
    return v == 'true' ? '--secondary_findings' : ''
}

process CPSR {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/cpsr" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(pcgr_data)
    path(vep_cache_cpsr)

    output:
    tuple val(meta), path("*.cpsr.grch38.html"),                  emit: html_report
    tuple val(meta), path("*.cpsr.grch38.classification.tsv.gz"), emit: tsv_report
    path "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sf_args = cpsrSecondaryArgs(params.cpsr_secondary_findings)
    """
    # Root under Docker (nextflow.config): give every file of the work directory back on any exit,
    # also when the task fails or is killed for memory. docker.fixOwnership gives them back only after
    # a task that succeeds, and only in an image whose chown knows GNU's --from (BusyBox's does not).
    # NXF_OWNER is set under docker.fixOwnership; -h leaves the files the staged links point to alone.
    trap '[ -z "\${NXF_OWNER:-}" ] || chown -Rh "\${NXF_OWNER}" . 2>/dev/null || true' EXIT
    cpsr_id=\$(cpsr_sample_id '${meta.id}')
    cpsr \\
        --input_vcf ${vcf} \\
        --vep_dir ${vep_cache_cpsr} \\
        --refdata_dir ${pcgr_data} \\
        --output_dir ./ \\
        --genome_assembly grch38 \\
        --sample_id "\$cpsr_id" \\
        --panel_id 0 \\
        ${sf_args} \\
        --force_overwrite
    cpsr_sample_id --rename . '${meta.id}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cpsr: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    def sf_args = cpsrSecondaryArgs(params.cpsr_secondary_findings)
    """
    echo "cpsr options: --panel_id 0 ${sf_args} --force_overwrite"
    cpsr_id=\$(cpsr_sample_id '${meta.id}')
    touch "\$cpsr_id".cpsr.grch38.html "\$cpsr_id".cpsr.grch38.classification.tsv.gz
    cpsr_sample_id --rename . '${meta.id}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cpsr: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
