#!/usr/bin/env bash
# assets/sv_rule_regions.tsv is what scripts/ci/build-sv-rule-regions.sh
# writes from its four pinned downloads, byte for byte: the builder's --check
# mode rebuilds it and compares. A hand edit, a pin that moved without a
# rebuild, or a builder change that moves a row fails here. Control: --check
# on a copy of the committed file with one cell changed must exit 1 and show
# that row. The downloads come from the e2e data cache when the runner filled
# it (entry catalogs, mirror first; scripts/ci/e2e-run.sh), else from their
# sources.
. "$(dirname "$0")/lib.sh"

COMMITTED="${REPO}/assets/sv_rule_regions.tsv"
CACHE="${CASE_TMP}/downloads"
[ -z "${E2E_DATA_CACHE:-}" ] || CACHE="${E2E_DATA_CACHE}/catalogs/sv-rule-regions"

echo "+ scripts/ci/build-sv-rule-regions.sh --check assets/sv_rule_regions.tsv"
SV_RULE_CACHE="$CACHE" bash "${REPO}/scripts/ci/build-sv-rule-regions.sh" --check "$COMMITTED" > "$STEP_LOG" 2>&1
rc=$?
cat "$STEP_LOG"
check_eq "--check on the committed table exits 0" "$rc" 0
check "and says it is what the pinned downloads give" grep -q "is what the pinned downloads give, byte for byte" "$STEP_LOG"
check "with the builder's row counts" grep -q "^Rows: gene_alias 10, ig_tr_locus 6, mhc 1, tandem_repeat 2409$" "$STEP_LOG"
echo "${CASE_NAME}: $(grep -c '^Downloading ' "$STEP_LOG") of 4 downloads not found in ${CACHE}" >> "$E2E_NOTES"

# The committed table's own counts, which docs/05-annotsv.md states.
check_eq "the committed table has 2409 tandem_repeat rows" \
  "$(awk -F'\t' '!/^#/ && $4 == "tandem_repeat"' "$COMMITTED" | wc -l | tr -d ' ')" 2409

# --- control: one cell changed ---------------------------------------------------
EDITED="${CASE_TMP}/edited.tsv"
awk -F'\t' -v OFS='\t' '!/^#/ && $4 == "mhc" {$3 = $3 + 1} 1' "$COMMITTED" > "$EDITED"
check_eq "control: the copy differs from the committed table in one line" \
  "$(diff "$COMMITTED" "$EDITED" | grep -c '^>')" 1
echo "+ scripts/ci/build-sv-rule-regions.sh --check on a copy with the MHC row's end moved by 1"
SV_RULE_CACHE="$CACHE" bash "${REPO}/scripts/ci/build-sv-rule-regions.sh" --check "$EDITED" > "${CASE_TMP}/edited.log" 2>&1
rc=$?
cat "${CASE_TMP}/edited.log"
check_eq "control: --check on the edited copy exits 1" "$rc" 1
check "control: and says to rebuild it" grep -q "is not what the pinned downloads give; rebuild it" "${CASE_TMP}/edited.log"
check "control: and shows the MHC row that differs" grep -qP '^< chr6\t\d+\t\d+\tmhc\t' "${CASE_TMP}/edited.log"
# Cache reuse: only without the e2e data cache, where the first run filled
# CACHE itself. With it the runner filled CACHE before the case, so the first
# run already downloads nothing and this could not go red.
if [ -z "${E2E_DATA_CACHE:-}" ]; then
  check_eq "control: the second run downloaded nothing, it read the first run's cache" \
    "$(grep -c '^Downloading ' "${CASE_TMP}/edited.log")" 0
fi

finish
