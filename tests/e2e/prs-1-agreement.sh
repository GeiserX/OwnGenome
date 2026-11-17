#!/usr/bin/env bash
# Step 25 runs pgsc_calc; this case ties its sum to a known answer and to the
# plink2 scoring step 25 did before it (cases 21 and 03 left a VCF and a gVCF
# of HG002):
#   - a synthetic score on the chr20 slice: five sites where HG002 matches the
#     reference (GIAB v4.2.1 truth, no variant within 200 bp) with the
#     reference base as effect allele (dosage 2), one such site whose effect
#     allele is another base (dosage 0), and three truth SNVs with their ALT
#     as effect allele, each with its own weight;
#   - pgsc_calc's sum is compared with plink2 --score (cols=+scoresums,
#     no-mean-imputation, the old path) on the genotypes step 25 handed
#     pgsc_calc, and the hom-ref part alone must be 2 x (1+2+4+8+16) = 62;
#   - control: plink2 on the same score with one weight changed must differ
#     from pgsc_calc's sum, so the comparison can fail;
#   - chrX: the score also has three rows on chrX outside the PARs, at
#     positions inside the gVCF's reference blocks (effect allele = reference
#     base). Without a sex step 25 leaves them out (the summary's ChrX says
#     "3 left out", the totals and sums above are the autosomes'); with male
#     pgsc_calc scores them: 12 of 12 rows, the 3 chrX rows matched in its
#     match log, and its sum less the autosomal sum above equals plink2's
#     --score of those rows on the same genotypes with the sex set.
# Writes nothing other cases read; removes what it made. When case prs-3
# measures (a run dispatched on a branch other than main) it also starts that
# case's download of the 1000 Genomes panel in the background, so the 7.4 GB
# arrive while cases prs-1 and prs-2 run.
. "$(dirname "$0")/lib.sh"

if [ "${PGSC_MEASURE:-}" = 1 ] || { [ "${GITHUB_EVENT_NAME:-}" = workflow_dispatch ] && [ "${GITHUB_REF:-}" != refs/heads/main ]; }; then
  rm -f "${E2E_WORK}/panel-download.status"
  # Its own session, so the case's timeout does not stop it; it ends on its own.
  setsid nohup timeout 2700 bash "$(dirname "$0")/prs-3-panel-measure.sh" --download \
    > "${E2E_WORK}/logs/panel-download.log" 2>&1 < /dev/null &
  echo "started the 1000 Genomes panel download for case prs-3"
fi

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
D="${GENOME_DIR}/${SAMPLE}"
check "case 21 wrote a gVCF" vcf_ok "${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz"

# --- sites -----------------------------------------------------------------------
SITES="${CASE_TMP}/sites.tsv"
# chrom, pos, effect allele, other allele, weight, kind: five hom-ref sites,
# a sixth with an effect allele that is neither the reference base nor its
# complement, and three truth SNVs that are not strand-ambiguous (pgsc_calc
# drops A/T and C/G pairs, plink2 would not).
"${REPO}/bin/og" dev fixture prs-sites "${FIXTURE_DIR}/HG002_truth_chr20.bed" "${FIXTURE_DIR}/HG002_truth_chr20.vcf.gz" \
  "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta" 6 3 --no-ambiguous \
  | awk -F'\t' -v OFS='\t' 'BEGIN {split("0.25 -1.5 0.125", w, " "); c["A"] = "T"; c["T"] = "A"; c["C"] = "G"; c["G"] = "C"}
      $1 == "homref" && h < 5 {print "20", $2, $3, "", 2 ^ h, "homref"; h++; next}
      $1 == "homref" {for (i = 1; i <= 4; i++) {o = substr("ACGT", i, 1); if (o != $3 && o != c[$3]) break}
                      print "20", $2, o, $3, "3.5", "homref_effect_not_ref"; next}
      $1 == "snv" {print "20", $2, $3, $4, w[++s], "snv"}' > "$SITES"
cat "$SITES"
check_eq "score sites (5 hom-ref, 1 hom-ref with another effect allele, 3 truth SNVs)" "$(wc -l < "$SITES" | tr -d ' ')" 9

# Three chrX rows outside the PARs: the middle of reference blocks (0/0) of
# 50 bp or more in the gVCF, more than 100 bp from any call with an ALT.
XR=chrX:73700001-74000000
bcf query -r "$XR" -f '%POS\t%INFO/END\t[%GT]\n' "${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz" > "${CASE_TMP}/x_blocks.tsv"
bcf query -r "$XR" -i 'GT="alt"' -f '%POS\n' "${SAMPLE}/vcf/${SAMPLE}.vcf.gz" > "${CASE_TMP}/x_vars.tsv"
X_SITES="${CASE_TMP}/x_sites.tsv"
"${REPO}/bin/og" dev fixture x-ref-sites "${CASE_TMP}/x_blocks.tsv" "${CASE_TMP}/x_vars.tsv" \
  "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta" chrX 3 \
  | awk -F'\t' -v OFS='\t' '{print "X", $1, $2, "", 2 ^ (NR + 4), "chrx"}' > "$X_SITES"
cat "$X_SITES"
check_eq "chrX rows of the score (positions in reference blocks of the gVCF)" "$(wc -l < "$X_SITES" | tr -d ' ')" 3

# The nine ids of assets/pgs_scores.tsv, each the synthetic score (step 25 scores the whole list).
SCORES="${GENOME_DIR}/prs_scores"
rm -rf "$SCORES" "${D}/prs"
mkdir -p "$SCORES"
for id in $(awk -F'\t' '$1 ~ /^PGS[0-9]+$/ {print $1}' "${REPO}/assets/pgs_scores.tsv"); do
  { printf '#pgs_id=%s\n#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n' "$id"
    cut -f1-5 "$SITES" "$X_SITES"; } | gzip -c > "${SCORES}/${id}.txt.gz"
done

# --- pgsc_calc (step 25) ---------------------------------------------------------
ANCESTRY_PANEL=none run_step 25-prs.sh "$SAMPLE"
check_step_exit 25-prs.sh
SUMMARY="${D}/prs/${SAMPLE}_prs_summary.tsv"
col() { awk -F'\t' -v k="$1" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $c["PGS_ID"] == "PGS000018" {print $c[k]}' "$SUMMARY" 2>/dev/null; }
PGSC_SUM=$(col Score_SUM)
check_eq "pgsc_calc matched every score site" "$(col Variants_Matched)" 9
check_eq "no sex: the score's autosomal rows only" "$(col Variants_Total)" 9
check_eq "no sex: the summary says its 3 chrX rows were left out" "$(col ChrX)" "3 left out"
check "no sex: the step log says so" has 'chrX rows of the scores are left out' "$(cat "$STEP_LOG")"
check_eq "input: genotypes from the gVCF" "$(col Input)" gvcf
check_eq "no percentile without a panel" "$(col Percentile)" NA
check "the raw-score line in the step log" has 'Raw score only' "$(cat "$STEP_LOG")"
check "the label comes from assets/pgs_scores.tsv" test "$(col Condition)" = "Coronary artery disease"
TARGET="${D}/prs/pgsc_calc/target.vcf.gz"
check "step 25 kept the genotypes it gave pgsc_calc" test -s "$TARGET"

# --- the old path: plink2 --score on the same genotypes ---------------------------
plink_sum() {  # plink_sum SCORE_TSV: SCORE1_SUM of chr:pos / allele / weight
  docker run --rm --network none --user "$(id -u):$(id -g)" -v "${CASE_TMP}:/w" -v "$(dirname "$TARGET"):/t:ro" -w /w "$PLINK2_IMAGE" \
    sh -c 'plink2 --vcf /t/target.vcf.gz --make-pgen --out g --set-all-var-ids "@:#" --new-id-max-allele-len 100 \
             --chr 1-22 --allow-extra-chr --output-chr chrM >/dev/null &&
           plink2 --pfile g --score "$1" 1 2 3 ignore-dup-ids no-mean-imputation cols=+scoresums --out s --allow-extra-chr >/dev/null' \
    _ "/w/$1" >/dev/null 2>&1 || { echo "plink2 failed"; return; }
  awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} NR == 2 {print $c["SCORE1_SUM"]}' "${CASE_TMP}/s.sscore"
}
awk -F'\t' -v OFS='\t' '{print "chr" $1 ":" $2, $3, $5}' "$SITES" > "${CASE_TMP}/score.tsv"
PLINK_SUM=$(plink_sum score.tsv)
awk -F'\t' -v OFS='\t' '$6 == "homref" {print "chr" $1 ":" $2, $3, $5}' "$SITES" > "${CASE_TMP}/homref.tsv"
HOMREF_SUM=$(plink_sum homref.tsv)
echo "pgsc_calc sum ${PGSC_SUM:-none}; plink2 sum ${PLINK_SUM:-none}; plink2 hom-ref part ${HOMREF_SUM:-none}"
close() { awk -v a="$1" -v b="$2" 'BEGIN { if (a == "" || b == "" || a !~ /^-?[0-9.e+-]+$/ || b !~ /^-?[0-9.e+-]+$/) exit 1; d = a - b; exit !(d < 1e-6 && d > -1e-6) }'; }
check "the hom-ref sites alone sum to 2 x (1+2+4+8+16) = 62 (plink2: ${HOMREF_SUM:-none})" close "${HOMREF_SUM:-}" 62
check "pgsc_calc's sum (${PGSC_SUM:-none}) equals plink2's (${PLINK_SUM:-none})" close "${PGSC_SUM:-}" "${PLINK_SUM:-}"
# Control: one weight changed, so the comparison above is shown able to fail.
awk -F'\t' -v OFS='\t' 'NR == 7 {$3 = $3 + 10} {print}' "${CASE_TMP}/score.tsv" > "${CASE_TMP}/score_changed.tsv"
CHANGED_SUM=$(plink_sum score_changed.tsv)
check "control: a changed weight gives a different sum (${CHANGED_SUM:-none} vs ${PGSC_SUM:-none})" \
  bash -c '! awk -v a="$1" -v b="$2" "BEGIN { d = a - b; exit !(d < 1e-6 && d > -1e-6) }"' _ "${CHANGED_SUM:-0}" "${PGSC_SUM:-0}"
printf '#### PRS: pgsc_calc against the old plink2 path (synthetic nine-site score on chr20)\n\npgsc_calc %s, plink2 %s, hom-ref part %s (want 62), changed-weight control %s\n\n' \
  "${PGSC_SUM:-none}" "${PLINK_SUM:-none}" "${HOMREF_SUM:-none}" "${CHANGED_SUM:-none}" >> "$E2E_NOTES"

# --- chrX with the sex ---------------------------------------------------------------
ANCESTRY_PANEL=none run_step 25-prs.sh "$SAMPLE" male
check_step_exit "25-prs.sh (male)"
PGSC_SUM_X=$(col Score_SUM)
check_eq "male: every row, the 3 on chrX too, is in the total" "$(col Variants_Total)" 12
check_eq "male: pgsc_calc matched all 12" "$(col Variants_Matched)" 12
check_eq "male: the summary says the 3 chrX rows were scored" "$(col ChrX)" "3 scored"
check "male: pgsc_calc's plink2 got the sex" has "update-sex" "$(cat "${D}/prs/pgsc_calc/images.config" 2>/dev/null)"
MATCH_LOG=$(find "${D}/prs/pgsc_calc/results" -name '*_log.csv.gz' 2>/dev/null | head -n 1)
X_MATCHED=$("${REPO}/bin/og" dev fixture csv-count "${MATCH_LOG:-/dev/null}" accession=PGS000018 chr_name=X,23 match_status=matched \
  2>/dev/null || echo "no match log")
check_eq "male: pgsc_calc's match log has the 3 chrX rows matched (${MATCH_LOG:-no log})" "$X_MATCHED" 3
# plink2 on the chrX rows alone, on the genotypes step 25 handed pgsc_calc, with the sex
printf '#IID\tSEX\n%s\t1\n' "$SAMPLE" > "${CASE_TMP}/sex.tsv"
awk -F'\t' -v OFS='\t' '{print "chr" $1 ":" $2, $3, $5}' "$X_SITES" > "${CASE_TMP}/x_score.tsv"
X_PART=$(docker run --rm --network none --user "$(id -u):$(id -g)" -v "${CASE_TMP}:/w" -v "$(dirname "$TARGET"):/t:ro" -w /w "$PLINK2_IMAGE" \
  sh -c 'plink2 --vcf /t/target.vcf.gz --update-sex sex.tsv --make-pgen --out gx --set-all-var-ids "@:#" --new-id-max-allele-len 100 \
           --chr 1-22,X --allow-extra-chr --output-chr chrM >/dev/null &&
         plink2 --pfile gx --score x_score.tsv 1 2 3 ignore-dup-ids no-mean-imputation cols=+scoresums --out sx --allow-extra-chr >/dev/null' \
  >/dev/null 2>&1 && awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} NR == 2 {print $c["SCORE1_SUM"]}' "${CASE_TMP}/sx.sscore")
echo "pgsc_calc sum with chrX ${PGSC_SUM_X:-none}; without ${PGSC_SUM:-none}; plink2 chrX part ${X_PART:-none}"
check "male: pgsc_calc's sum (${PGSC_SUM_X:-none}) is the autosomal sum (${PGSC_SUM:-none}) plus plink2's chrX part (${X_PART:-none})" \
  close "$(awk -v a="${PGSC_SUM_X:-}" -v b="${PGSC_SUM:-}" 'BEGIN { if (a == "" || b == "") exit; print a - b }')" "${X_PART:-}"
check "control: the chrX part is not zero, so the sum above shows the rows were scored" \
  bash -c '! awk -v a="$1" "BEGIN { exit !(a + 0 == 0) }"' _ "${X_PART:-0}"
printf '#### PRS on chrX (3 hom-ref rows, weights 32, 64, 128)\n\nno sex: %s of 9 rows, sum %s; male: %s of 12 rows, sum %s; plink2 chrX part %s\n\n' \
  "9" "${PGSC_SUM:-none}" "$(col Variants_Matched)" "${PGSC_SUM_X:-none}" "${X_PART:-none}" >> "$E2E_NOTES"

rm -rf "$SCORES" "${D:?}/prs"
finish
