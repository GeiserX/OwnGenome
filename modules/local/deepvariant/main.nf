/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEEPVARIANT — Small variants from the BAM: a VCF and a gVCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    run_deepvariant with the flags of scripts/03-deepvariant.sh: the model
    of --model_type (WGS, or WES for an exome, as the script's MODEL_TYPE),
    the sample name from the samplesheet, one shard per CPU, and a gVCF
    beside the VCF. The gVCF also records where the sample matches the
    reference; PRS and PharmCAT read their hom-ref genotypes from it.

    Sex: for a male sample chrX and chrY are called haploid outside the
    pseudoautosomal regions (assets/par_grch38.bed), so no impossible
    heterozygous call is made there. A female sample is diploid everywhere.
    main.nf requires the sex on every row this process calls. INDEXCOV has
    checked it against the BAM before this runs; an exome's check is
    somalier's (SAMPLE_QC), which runs beside this.

    --intervals (the script's INTERVALS) limits calling to space-separated
    regions, passed to --regions. An exome is called on the regions of
    --capture_bed instead (main.nf refuses both): reads outside the capture
    are off-target strays, and DeepVariant's exome model is trained on the
    targets.

    Equivalent to: scripts/03-deepvariant.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// The model and the regions, written once: the script passes them and the
// stub prints them, so tests/e2e/nf-exome-1-stub.sh checks what runs.
def deepvariant_model() {
    return "--model_type=${params.model_type}"
}

// An exome's capture BED, else --intervals, else the whole genome ('')
def deepvariant_regions(capture_bed) {
    return capture_bed ? "--regions ${capture_bed}" : (params.intervals ? "--regions \"${params.intervals}\"" : '')
}

process DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'

    publishDir { "${params.outdir}/${meta.id}/vcf" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_fai)  // staged beside the FASTA
    path(par_bed)
    path(capture_bed)    // --capture_bed with --model_type WES, or []

    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), path("${meta.id}.vcf.gz.tbi"),     emit: vcf
    tuple val(meta), path("${meta.id}.g.vcf.gz"), path("${meta.id}.g.vcf.gz.tbi"), emit: gvcf
    path "versions.yml",                                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def haploid = meta.sex == 'male' ? "--haploid_contigs=chrX,chrY --par_regions_bed=${par_bed}" : ''
    def regions = deepvariant_regions(capture_bed)
    """
    /opt/deepvariant/bin/run_deepvariant \\
        ${deepvariant_model()} \\
        --ref=${reference} \\
        --reads=${bam} \\
        --output_vcf=${meta.id}.vcf.gz \\
        --output_gvcf=${meta.id}.g.vcf.gz \\
        --intermediate_results_dir=deepvariant_tmp \\
        --sample_name=${meta.id} \\
        --num_shards=${task.cpus} \\
        ${haploid} \\
        ${regions}
    rm -rf deepvariant_tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    // The stub prints the options the script passes (tests/e2e/nf-exome-1-stub.sh reads them).
    """
    echo "run_deepvariant options: ${deepvariant_model()} ${deepvariant_regions(capture_bed)}"
    echo | gzip > ${meta.id}.vcf.gz
    echo | gzip > ${meta.id}.g.vcf.gz
    touch ${meta.id}.vcf.gz.tbi ${meta.id}.g.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
