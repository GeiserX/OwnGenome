#!/usr/bin/env bash
# Traits, APOE ordering, phase sites and refusal controls on synthetic input.
. "$(dirname "$0")/lib.sh"
OG=${OG:-${REPO}/bin/og}
TABLE="${REPO}/assets/traits.tsv"
RULES="${REPO}/assets/genotype_rules.tsv"
G="${CASE_TMP}/g.tsv"
printf 'chr19\t44908684\tT\tC\tT/C\t40\nchr19\t44908822\tC\tT\tC/T\t40\nchr7\t141973545\tC\tG\tC/G\t30\nchr7\t141972905\tG\tA\tG/A\t31\nchr7\t141972804\tT\tC\tT/C\t32\n' > "$G"
"$OG" traits table --table "$TABLE" --rules "$RULES" --genotypes "$G" --sample S \
  --out "${CASE_TMP}/traits.tsv" --rules-out "${CASE_TMP}/rules.tsv"
check_eq 'traits table exits 0' "$?" 0
check_eq 'APOE double heterozygote is e2/e4' "$(awk -F'\t' '$1 == "APOE" {print $10}' "${CASE_TMP}/rules.tsv")" e2/e4
check_eq 'the taster rule needs phase' "$(awk -F'\t' '$1 == "TAS2R38" {print $11}' "${CASE_TMP}/rules.tsv")" unphased
check_eq 'APOE line occurs exactly once' "$(grep -c '^rs429358+rs7412' "${CASE_TMP}/traits.tsv")" 1
check 'absent LCT has no inferred genotype' awk -F'\t' '$1 == "rs4988235" {ok = ($6 == "." && $7 == "not_covered")} END {exit !ok}' "${CASE_TMP}/traits.tsv"
"$OG" traits phase-sites --rules "$RULES" --genotypes "$G" --sample S --out "${CASE_TMP}/phase.vcf"
check_eq 'phase-sites exits 0' "$?" 0
check_eq 'only the three taster sites need phase' "$(grep -cv '^#' "${CASE_TMP}/phase.vcf")" 3
check_eq 'phase sites are in genome order' "$(grep -v '^#' "${CASE_TMP}/phase.vcf" | cut -f2 | paste -sd, -)" 141972804,141972905,141973545
sed 's/44908684\tT\t/44908684\tG\t/' "$G" > "${CASE_TMP}/bad.tsv"
"$OG" traits table --table "$TABLE" --rules "$RULES" --genotypes "${CASE_TMP}/bad.tsv" --sample S \
  --out "${CASE_TMP}/bad-out.tsv" --rules-out "${CASE_TMP}/bad-rules.tsv" > "${CASE_TMP}/bad.log" 2>&1
check_eq 'control: wrong reference is refused' "$?" 1
check 'control: refusal names the site and GRCh38' has 'rs429358.*GRCh38' "$(cat "${CASE_TMP}/bad.log")"
check 'control: no table is written' test ! -e "${CASE_TMP}/bad-out.tsv"
"$OG" traits check --table "$TABLE" --rules "$RULES"
check_eq 'committed catalogs pass' "$?" 0
finish
