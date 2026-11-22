#!/usr/bin/env bash
# Step 08 with KIR=true (opt-in): a second T1K pass over the KIR genes
# against IPD-KIR release KIR_DB_RELEASE, which `setup.sh --kir-data`
# installs. The fixture has no slice of the KIR region (chr19 leukocyte
# receptor complex), so the step may find too few reads: then the genotype
# file says so. Either way the release is written beside it, and step 24's
# report (and the text report from the same summary) shows the KIR result in
# the HLA section with the caveat that no step interprets it. Red control: a
# summary without the KIR block renders reports that fail that check.
. "$(dirname "$0")/lib.sh"

echo "+ scripts/setup.sh --kir-data"
"${REPO}/scripts/setup.sh" --kir-data "$GENOME_DIR"
check "setup.sh --kir-data installed IPD-KIR ${KIR_DB_RELEASE}" test -s "${GENOME_DIR}/kir/IPD-KIR_${KIR_DB_RELEASE}/kir.dat"

KIR=true run_step 08-hla-typing.sh "$SAMPLE"
check_step_exit 08-hla-typing.sh
K="${GENOME_DIR}/${SAMPLE}/kir_t1k"
cat "${K}/database_release.txt" "${K}/${SAMPLE}_kir_genotype.tsv" 2>/dev/null | head -n 20
check "database_release.txt names IPD-KIR ${KIR_DB_RELEASE}" \
  grep -q "^database: IPD-KIR Release Version ${KIR_DB_RELEASE//./\\.}$" "${K}/database_release.txt"
check "the KIR genotype file is written" test -s "${K}/${SAMPLE}_kir_genotype.tsv"
check "it holds KIR genotypes or says too few reads" \
  grep -Eq $'^KIR[0-9A-Z]+\t|^# KIR not typed: T1K found too few reads' "${K}/${SAMPLE}_kir_genotype.tsv"
check "the HLA pass still wrote its genotypes" test -s "${GENOME_DIR}/${SAMPLE}/hla_t1k/${SAMPLE}_hla_genotype.tsv"

# Both reports: the KIR result, as the genotype file says.
run_step 24-html-report.sh "$SAMPLE"
check_step_exit 24-html-report.sh
HTML="${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.html"
TXT="${CASE_TMP}/${SAMPLE}_report.txt"
"${REPO}/bin/og" report render --summary "${GENOME_DIR}/${SAMPLE}/summary.json" -o "$TXT" >/dev/null
sed -n '/KIR genes/,/^$/p' "$TXT"
# kir_shown TXT HTML: each typed gene, or the reason none was typed, and the caveat, in both reports.
kir_shown() {
  local f g
  for f in "$1" "$2"; do
    if grep -q $'^KIR' "${K}/${SAMPLE}_kir_genotype.tsv"; then
      while IFS=$'\t' read -r g _; do grep -qF "${g}*" "$f" || return 1
      done < <(awk -F'\t' '/^KIR/ && (($3 != "" && $3 != ".") || ($6 != "" && $6 != "."))' "${K}/${SAMPLE}_kir_genotype.tsv")
      g=$(awk -F'\t' '/^KIR/ && ($3 == "" || $3 == ".") && ($6 == "" || $6 == ".") {print $1}' \
        "${K}/${SAMPLE}_kir_genotype.tsv" | paste -sd, - | sed 's/,/, /g')
      [ -z "$g" ] || grep -qF "No allele called (absent, or too few reads): ${g}" "$f" || return 1
    else
      grep -qF "$(sed -n '1s/^# *//p' "${K}/${SAMPLE}_kir_genotype.tsv")" "$f" || return 1
    fi
    grep -qF "No step interprets them" "$f" || return 1
  done
  grep -qF 'KIR genes (T1K, KIR=true):' "$1"
}
not_kir_shown() { ! kir_shown "$@"; }
check "both reports show the KIR result and say no step interprets it" kir_shown "$TXT" "$HTML"
# Red control: the reports of the same summary without its KIR block.
jq 'del(.sections.hla.data.kir)' "${GENOME_DIR}/${SAMPLE}/summary.json" > "${CASE_TMP}/ctl.summary.json"
check "the control summary has no KIR block" \
  jq -e '.sections.hla.data | type == "object" and (has("kir") | not)' "${CASE_TMP}/ctl.summary.json"
"${REPO}/bin/og" report render --summary "${CASE_TMP}/ctl.summary.json" -o "${CASE_TMP}/ctl.txt" -o "${CASE_TMP}/ctl.html" >/dev/null
check "control: without the KIR block the check above rejects both reports" \
  not_kir_shown "${CASE_TMP}/ctl.txt" "${CASE_TMP}/ctl.html"
echo "- KIR on the fixture: $(grep -c $'^KIR' "${K}/${SAMPLE}_kir_genotype.tsv" 2>/dev/null || echo 0) gene rows typed" >> "$E2E_NOTES"

finish
