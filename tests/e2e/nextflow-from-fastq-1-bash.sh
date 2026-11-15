#!/usr/bin/env bash
# The bash leg of the parity check (scripts/ci/parity-diff.sh): the single
# step scripts on the fixture reads, as sample HG002P: 01b (fastp), 02
# (minimap2, markdup), 03 (DeepVariant, male), 06, 08 and 36 (T1K's HLA types
# as PharmCAT's outside calls), 20 (chrM calls, PharmCAT's MT-RNR1 call, as
# the pipeline's PHARMCAT reads MITO_VARIANTS'), 07, 11 and 25. The next case runs the
# Nextflow pipeline on the same reads under the same sample name; the E2E
# workflow then compares the two.
#
# PRS scores a synthetic file (the PGS ids of assets/pgs_scores.tsv so nothing
# is downloaded): five chr20 sites where HG002 matches the reference (GIAB
# v4.2.1 truth, no variant within 200 bp; effect allele = reference base) and
# three truth SNVs (effect allele = ALT). The next case scores the same file.
. "$(dirname "$0")/lib.sh"

P=HG002P
D="${GENOME_DIR}/${P}"
rm -rf "$D"
mkdir -p "${D}/fastq"
for r in R1 R2; do
  ln -f "${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_${r}.fastq.gz" "${D}/fastq/${P}_${r}.fastq.gz"
done

run_step 01b-fastp-qc.sh "$P"
check_step_exit 01b-fastp-qc.sh
run_step 02-alignment.sh "$P"
check_step_exit 02-alignment.sh
check "step 02 aligned the trimmed reads" has 'Using trimmed FASTQs from fastp' "$(cat "$STEP_LOG")"
# The fixture's slices only, as case 21 (the rest of the reference has no reads)
INTERVALS=$(awk '{printf "%s%s:%d-%d", (NR > 1 ? " " : ""), $1, $2 + 1, $3}' "${FIXTURE_DIR}/regions.bed")
INTERVALS="$INTERVALS" run_step 03-deepvariant.sh "$P" male
check_step_exit 03-deepvariant.sh
run_step 06-clinvar-screen.sh "$P"
check_step_exit 06-clinvar-screen.sh
# HLA types first: PharmCAT reads them as outside calls (step 36), as the
# pipeline's PHARMCAT reads PGX_CONSENSUS's.
run_step 08-hla-typing.sh "$P"
check_step_exit 08-hla-typing.sh
run_step 36-pgx-consensus.sh "$P"
check_step_exit 36-pgx-consensus.sh
run_step 20-mtoolbox.sh "$P"
check_step_exit 20-mtoolbox.sh
run_step 07-pharmacogenomics.sh "$P"
check_step_exit 07-pharmacogenomics.sh
run_step 11-roh-analysis.sh "$P"
check_step_exit 11-roh-analysis.sh

# --- synthetic PRS scores ---------------------------------------------------------
SITES="${CASE_TMP}/score_sites.tsv"
# Five hom-ref sites, then three of the truth set's SNVs (chrom, pos, effect allele).
"${REPO}/bin/og" dev fixture prs-sites "${FIXTURE_DIR}/HG002_truth_chr20.bed" "${FIXTURE_DIR}/HG002_truth_chr20.vcf.gz" \
  "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta" 5 3 | awk -F'\t' -v OFS='\t' '{print "chr20", $2, $3}' > "$SITES"
cat "$SITES"
check_eq "score sites (5 hom-ref, 3 truth SNVs)" "$(wc -l < "$SITES" | tr -d ' ')" 8
PGS="${E2E_WORK}/parity-pgs"
rm -rf "$PGS" "${GENOME_DIR}/prs_scores"
mkdir -p "$PGS" "${GENOME_DIR}/prs_scores"
# One file per id of assets/pgs_scores.tsv, which step 25 scores.
for id in $(awk -F'\t' '$1 ~ /^PGS[0-9]+$/ {print $1}' "${REPO}/assets/pgs_scores.tsv"); do
  { printf '#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\teffect_weight\n'
    awk -F'\t' -v OFS='\t' '{print "20", $2, $3, 2 ^ (NR - 1)}' "$SITES"; } | gzip -c > "${PGS}/${id}.txt.gz"
done
cp "${PGS}"/*.txt.gz "${GENOME_DIR}/prs_scores/"
run_step 25-prs.sh "$P"
check_step_exit 25-prs.sh
rm -rf "${GENOME_DIR}/prs_scores"

# --- what the parity check reads ------------------------------------------------
for f in "aligned/${P}_sorted.bam" "vcf/${P}.vcf.gz" "vcf/${P}.g.vcf.gz" "clinvar/${P}_clinvar_hits.tsv" \
         "vcf/${P}_roh.txt" "vcf/${P}.report.json" "prs/${P}_prs_summary.tsv"; do
  check "bash leg wrote ${f}" nonempty "${P}/${f}"
done
check_eq "PRS input of the bash leg" "$(awk -F'\t' 'NR == 2 {print $NF}' "${D}/prs/${P}_prs_summary.tsv" 2>/dev/null)" gvcf

finish
