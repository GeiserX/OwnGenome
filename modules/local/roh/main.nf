/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ROH — Runs of Homozygosity (consanguinity screening)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Detects long stretches of homozygous DNA that indicate shared ancestry.
    The summary keeps chr1 to chr22 only and leaves out the part of each segment in a
    centromere or heterochromatin band (masked_bands below), where few
    markers map. With --roh_af_panel (setup.sh installs it under roh_af/)
    bcftools roh reads each site's population frequency from the panel and
    uses the panel's sites only: the global AF, or one continental group's
    when step 26's ancestry table gives that group 0.8 or more. Without it
    every site gets 0.4 (--AF-dflt 0.4).

    Equivalent to: scripts/11-roh-analysis.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ROH {
    tag "$meta.id"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/roh" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index), path(ancestry)   // ancestry: step 26's table, or []
    tuple path(af_panel), path(af_panel_index)                     // --roh_af_panel and its .tbi, or [] and []

    output:
    tuple val(meta), path("*_roh.txt"),         emit: roh_regions
    tuple val(meta), path("*_roh_summary.txt"), emit: roh_summary
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def panel   = af_panel ? af_panel.name : ''
    def anc     = ancestry ? ancestry.name : ''
    """
    # Auto-detect chip data: if FORMAT/PL is absent, use -G30 (genotype-only mode)
    HAS_PL=\$(bcftools view -h ${vcf} | grep -c '##FORMAT=<ID=PL' || true)

    # The allele frequencies, picked and cut as scripts/11-roh-analysis.sh
    # does: the panel's global AF, or one continental group's when step 26
    # gives it 0.8 or more (EUR reads gnomAD's nfe), in a 4-column file of
    # the group's AF of 0.01 to 0.99 that og roh af-file writes and indexes.
    # --AF-file restricts the markers to the panel's sites by itself. Without
    # the panel: 0.4 at every site.
    set -o pipefail
    AF_LINE=""
    if [ -n "${panel}" ]; then
        GROUP=\$( [ ! -s "${anc}" ] || awk -F'\\t' '\$1 ~ /^probability_(AFR|AMR|EAS|EUR|SAS)\$/ && \$2 + 0 >= 0.8 {print substr(\$1, 13), \$2}' "${anc}")
        if [ -z "\${GROUP}" ]; then
            ROH_FLAGS=(--AF-file "${panel}")
            # No table: no ancestry asked for, or PRS failed (it is
            # report-only, so ROH goes on without its table).
            if [ -s "${anc}" ]; then
                AF_DESC="global (panel ${panel}); no step 26 group at 0.8 or more"
            else
                AF_DESC="global (panel ${panel}); no step 26 ancestry table"
            fi
        else
            read -r POP PROB <<< "\${GROUP}"
            case "\${POP}" in AFR) GNOMAD_GROUP=afr ;; AMR) GNOMAD_GROUP=amr ;; EAS) GNOMAD_GROUP=eas ;; EUR) GNOMAD_GROUP=nfe ;; SAS) GNOMAD_GROUP=sas ;; esac
            og roh af-file --panel "${panel}" --group "\${GNOMAD_GROUP}" --out ${meta.id}_roh_af.tsv.gz
            ROH_FLAGS=(--AF-file "${meta.id}_roh_af.tsv.gz")
            AF_DESC="\${GNOMAD_GROUP} (panel ${panel}); step 26 gives \${POP} \${PROB}"
        fi
        AF_LINE="# AF used: \${AF_DESC}"
    else
        ROH_FLAGS=(--AF-dflt 0.4)
    fi
    if [ "\${HAS_PL}" -eq 0 ]; then
        ROH_FLAGS+=(-G30)
    fi

    # Only PASS calls (and records with no filter, as chip VCFs have), as in
    # scripts/11-roh-analysis.sh. pipefail: a failed view must fail the task.
    bcftools view -f PASS,. -Ou ${vcf} \\
        | bcftools roh "\${ROH_FLAGS[@]}" -o ${meta.id}_roh.txt -
    # The table names the frequencies it was found with, when they came from
    # the panel, as the bash step writes it.
    if [ -n "\${AF_LINE}" ]; then
        { printf '%s\\n' "\${AF_LINE}"; cat ${meta.id}_roh.txt; } > ${meta.id}_roh.txt.tmp
        mv -f ${meta.id}_roh.txt.tmp ${meta.id}_roh.txt
    fi

    # Summary: autosomal segments with 5 Mb or more outside GRCh38's
    # centromere and heterochromatin bands (MASKED_BANDS in og's step 11
    # code, the one copy), the threshold of docs/11-roh-analysis.md;
    # length_bp and length_mb are of that part. The bash step runs the same
    # command.
    og roh summary --sample ${meta.id} --roh ${meta.id}_roh.txt --out ${meta.id}_roh_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_roh.txt ${meta.id}_roh_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
