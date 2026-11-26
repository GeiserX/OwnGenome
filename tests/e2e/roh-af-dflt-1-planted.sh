#!/usr/bin/env bash
# Step 11 with and without the allele-frequency panel, on a planted chr22
# stretch: 300 sites 10 kb apart, heterozygous at the 30 on each end and
# homozygous ALT at the 240 between them. A separate PASS site is absent
# from the panel. Every panel ALT is common
# (global AF 0.95), so the homozygous stretch is what anyone would carry and
# no run should be found there. With a flat 0.4 at every site (no panel)
# bcftools roh finds one. The panel's European (nfe) column gives the same
# sites AF 0.4, so a step 26 table with EUR at 0.9 brings the run back:
# the group column is the one read. Without the panel, a record's INFO/AC
# and INFO/AN change nothing (--AF-dflt 0.4 is used everywhere), and neither
# the docs nor the step say otherwise.
. "$(dirname "$0")/lib.sh"

G="${CASE_TMP}/genome"
PANEL_NAME="${ROH_AF_PANEL_URL##*/}"
OFF_PANEL_POS=23010000
rm -rf "$G"
mkdir -p "${G}/planted/vcf" "${G}/plantednoac/vcf" "${G}/roh_af"
{
  printf '##fileformat=VCFv4.2\n##contig=<ID=chr22,length=50818468>\n'
  printf '##FILTER=<ID=PASS,Description="All filters passed">\n'
  printf '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">\n'
  printf '##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">\n'
  printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
  printf '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tplanted\n'
  for i in $(seq 1 300); do
    if [ "$i" -le 30 ] || [ "$i" -gt 270 ]; then gt=0/1 pl=30,0,30 ac=1; else gt=1/1 pl=30,30,0 ac=2; fi
    printf 'chr22\t%d\t.\tA\tG\t50\tPASS\tAC=%d;AN=2\tGT:PL\t%s:%s\n' $((20000000 + i * 10000)) "$ac" "$gt" "$pl"
  done
  printf 'chr22\t%d\t.\tA\tG\t50\tPASS\tAC=1;AN=2\tGT:PL\t0/1:30,0,30\n' "$OFF_PANEL_POS"
} > "${G}/planted/vcf/planted.vcf"
grep -v '^##INFO' "${G}/planted/vcf/planted.vcf" | sed 's/AC=[0-9];AN=2/./; s/\tplanted$/\tplantednoac/' \
  > "${G}/plantednoac/vcf/plantednoac.vcf"
{
  printf '# planted panel for the e2e case\n#CHROM\tPOS\tREF,ALT\tAF\tAF_afr\tAF_amr\tAF_eas\tAF_nfe\tAF_sas\n'
  for i in $(seq 1 300); do printf 'chr22\t%d\tA,G\t0.95\t0.95\t0.95\t0.95\t0.40\t.\n' $((20000000 + i * 10000)); done
} > "${G}/roh_af/${PANEL_NAME%.gz}"
docker run --rm --user "$(id -u):$(id -g)" -v "${G}:/g" -w /g "$BCFTOOLS_IMAGE" bash -euo pipefail -c "
  for s in planted plantednoac; do bgzip -f \$s/vcf/\$s.vcf; tabix -f -p vcf \$s/vcf/\$s.vcf.gz; done
  bgzip -f 'roh_af/${PANEL_NAME%.gz}'; tabix -f -s1 -b2 -e2 'roh_af/${PANEL_NAME}'"
check "the planted VCFs and panel are indexed" test -s "${G}/roh_af/${PANEL_NAME}.tbi" -a -s "${G}/plantednoac/vcf/plantednoac.vcf.gz.tbi"

rg() { grep -c '^RG' "${G}/$1/vcf/$1_roh.txt" || true; }
body() { grep -v '^#' "${G}/$1/vcf/$1_roh.txt" | cut -f1,3- ; }
off_panel_sites() {
  awk -v pos="$OFF_PANEL_POS" '$1 == "ST" && $3 == "chr22" && $4 == pos {n++} END {print n+0}' \
    "${G}/planted/vcf/planted_roh.txt"
}

# 1. No panel: a flat 0.4, the same with or without INFO/AC,AN.
mv "${G}/roh_af" "${G}/roh_af.off"
GENOME_DIR="$G" run_step 11-roh-analysis.sh plantednoac
check_step_exit 11-roh-analysis.sh
GENOME_DIR="$G" run_step 11-roh-analysis.sh planted
check_step_exit 11-roh-analysis.sh
NOTE_LOG=$(cat "$STEP_LOG")
check_ge "runs on the planted stretch with AF 0.4 everywhere" "$(rg planted)" 1
check_eq "with flat AF: the off-panel PASS marker has an ST row" "$(off_panel_sites)" 1
check "no panel: the ST and RG lines are the same with and without INFO/AC,AN" \
  cmp -s <(body planted) <(body plantednoac)
check_eq "no panel: the table has no AF line of ours" "$(grep -c '^# AF used:' "${G}/planted/vcf/planted_roh.txt" || true)" 0
check "the step's note says 0.4 is used at every site" has 'gives every site an' "$NOTE_LOG"
check "the step's note no longer says bcftools roh takes INFO/AC" lacks 'takes INFO/AC and INFO/AN where a record' "$NOTE_LOG"
check "docs/11 no longer says bcftools roh reads INFO/AC with --AF-dflt" \
  lacks "reads a site's frequency from INFO/AC" "$(cat "${REPO}/docs/11-roh-analysis.md")"

# 2. The panel, no step 26 table: the global AF (0.95), no run.
mv "${G}/roh_af.off" "${G}/roh_af"
GENOME_DIR="$G" run_step 11-roh-analysis.sh planted
check_step_exit 11-roh-analysis.sh
check_eq "with the panel's global AF: runs on the planted stretch" "$(rg planted)" 0
check_eq "with global AF: the off-panel PASS marker has no ST row" "$(off_panel_sites)" 0
check_eq "with the panel: per-site lines at all 300 panel sites" "$(grep -c '^ST' "${G}/planted/vcf/planted_roh.txt" || true)" 300
check "the table names the global AF" has '^# AF used: global ' "$(head -n 1 "${G}/planted/vcf/planted_roh.txt")"

# 3. Step 26 puts the sample at EUR 0.9: the nfe column (0.4), the run is back.
mkdir -p "${G}/planted/ancestry"
printf 'key\tvalue\nsample\tplanted\npopulation\tEUR\nprobability_AFR\t0.02\nprobability_AMR\t0.05\nprobability_EAS\t0.01\nprobability_EUR\t0.9\nprobability_SAS\t0.02\n' \
  > "${G}/planted/ancestry/planted_ancestry.tsv"
GENOME_DIR="$G" run_step 11-roh-analysis.sh planted
check_step_exit 11-roh-analysis.sh
check_ge "with the nfe AF of a EUR sample: runs on the planted stretch" "$(rg planted)" 1
check_eq "with nfe AF: the off-panel PASS marker has no ST row" "$(off_panel_sites)" 0
check "the table names the nfe AF and step 26's EUR" has '^# AF used: nfe .*EUR 0.9$' "$(head -n 1 "${G}/planted/vcf/planted_roh.txt")"
check_eq "the group's AF file is removed after the run" "$(find "${G}/planted/vcf" -name 'planted_roh_af.tsv.gz*' | wc -l | tr -d ' ')" 0

# 4. EUR at 0.7: no group reaches 0.8, back to the global AF.
sed -i.bak 's/^probability_EUR\t0.9$/probability_EUR\t0.7/' "${G}/planted/ancestry/planted_ancestry.tsv"
GENOME_DIR="$G" run_step 11-roh-analysis.sh planted
check_step_exit 11-roh-analysis.sh
check_eq "no group at 0.8 or more: the global AF, no run" "$(rg planted)" 0
check_eq "with global AF after ancestry fallback: the off-panel PASS marker has no ST row" "$(off_panel_sites)" 0

finish
