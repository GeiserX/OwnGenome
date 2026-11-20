/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SURVIVOR_MERGE — SV consensus: the calls two or more callers agree on
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SURVIVOR_PREP   the calls each caller votes with, as plain VCF, one sample
                    column named after the caller (og sv survivor-prep): no
                    record whose FILTER or sample FT is not PASS, none whose
                    genotype is reference only (0/0, 0|0, 0, 0/.), and for
                    CNVpytor none that fails step 05's e-value and q0 rule
    SURVIVOR_MERGE  `SURVIVOR merge LIST 1000 2 1 1 0 50`: both breakpoints
                    within 1,000 bp, same type and strands, >= 50 bp, kept
                    when 2+ callers support it; SUPP and SUPP_VEC say which
    SURVIVOR_SORT   sorted, bgzipped and indexed (bcftools)

    Equivalent to: scripts/22-survivor-merge.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SURVIVOR_PREP {
    tag "$meta.id:$caller"
    label 'process_single'

    input:
    tuple val(meta), val(caller), path(vcf)

    output:
    tuple val(meta), val(caller), path("${caller}.vcf"), emit: vcf
    path "versions.yml",                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // SURVIVOR names its output columns after the input samples: three
    // columns all named after the sample would be one name three times.
    // og (bin/) reads the gzip input itself; scripts/22-survivor-merge.sh
    // runs the same verb.
    """
    og sv survivor-prep --caller ${caller} ${vcf} ${caller}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${caller}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process SURVIVOR_MERGE {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), val(callers), path(vcfs)

    output:
    tuple val(meta), path("${meta.id}_sv_merged.vcf"), emit: vcf
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The list in a fixed caller order, so SUPP_VEC reads the same on every run.
    def order = ['manta', 'delly', 'cnvpytor']
    def listed = callers.sort(false) { c -> order.indexOf(c) }.collect { c -> "${c}.vcf" }.join(' ')
    """
    printf '%s\\n' ${listed} > sv_files.txt
    SURVIVOR merge sv_files.txt 1000 2 1 1 0 50 ${meta.id}_sv_merged.vcf
    # SURVIVOR exits 0 when it cannot open an input, so its output is the check.
    head -n 1 ${meta.id}_sv_merged.vcf | grep -q '^##fileformat=VCF' \\
        || { echo "ERROR: SURVIVOR wrote no VCF" >&2; exit 1; }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_sv_merged.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process SURVIVOR_SORT {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/sv_merged" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(merged)

    output:
    tuple val(meta), path("${meta.id}_sv_consensus.vcf.gz"),     emit: merged_vcf
    tuple val(meta), path("${meta.id}_sv_consensus.vcf.gz.tbi"), emit: merged_vcf_index
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools sort ${merged} -Oz -o ${meta.id}_sv_consensus.vcf.gz
    bcftools index -t ${meta.id}_sv_consensus.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_sv_consensus.vcf.gz
    touch ${meta.id}_sv_consensus.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
