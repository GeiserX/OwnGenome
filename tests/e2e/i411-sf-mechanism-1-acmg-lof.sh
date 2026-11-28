#!/usr/bin/env bash
# The clingen_hi and lof columns of assets/acmg_sf.tsv, and TTN's limit cell,
# are what scripts/ci/build-acmg-lof.sh writes from its pinned downloads
# (ClinGen's gene curation list and cardiodb's titin exon table, each checked
# against its sha256) and assets/acmg_lof_overrides.tsv: a hand edit of one
# of those cells fails here. Controls: a copy of the table with DSG2's lof
# turned to no, and one with an end of a TTN region moved by one base, must
# fail the same check and name the gene. The downloads come from the e2e
# data cache when the runner filled it (entry catalogs, mirror first;
# scripts/ci/e2e-run.sh), else from the Internet Archive into the case's own
# folder.
. "$(dirname "$0")/lib.sh"

CACHE="${CASE_TMP}/downloads"
[ -z "${E2E_DATA_CACHE:-}" ] || CACHE="${E2E_DATA_CACHE}/catalogs/acmg-lof"
mkdir -p "$CACHE"
TABLE="${REPO}/assets/acmg_sf.tsv"
echo "+ scripts/ci/build-acmg-lof.sh --check"
bash "${REPO}/scripts/ci/build-acmg-lof.sh" --cache "$CACHE" --check 2>&1 | tee "$STEP_LOG"
check_eq "build-acmg-lof.sh --check exits 0 on the committed table" "${PIPESTATUS[0]}" 0
# shellcheck disable=SC2016  # the positional arguments of bash -c
check "the check printed no difference" bash -c '! grep -q "differs" "$1"' _ "$STEP_LOG"
check_eq "both pinned downloads are in the cache" "$(find "$CACHE" -maxdepth 1 -type f \( -name clingen_gene_curation_list_GRCh38.tsv -o -name cardiodb_titin_transcripts.html \) | wc -l | tr -d ' ')" 2
TTN_REGIONS=$(awk -F'\t' '$1 == "TTN" { n = split($NF, a, ","); print n }' "$TABLE")
check_eq "TTN's limit cell holds the 7 runs of exons with PSI over 90%" "$TTN_REGIONS" 7

# Control 1: DSG2's hand row says yes; the table says no.
EDITED="${CASE_TMP}/edited-lof.tsv"
awk -F'\t' -v OFS='\t' '!/^#/ && NR > 1 && $1 == "DSG2" { for (i = 1; i <= NF; i++) if (h[i] == "lof") $i = "no" }
  !/^#/ && !seen++ { for (i = 1; i <= NF; i++) h[i] = $i } 1' "$TABLE" > "$EDITED"
# shellcheck disable=SC2016  # the positional arguments of bash -c
check "control: the edited copy differs from the table" bash -c '! cmp -s "$1" "$2"' _ "$EDITED" "$TABLE"
bash "${REPO}/scripts/ci/build-acmg-lof.sh" --cache "$CACHE" --table "$EDITED" --check > "${CASE_TMP}/edited-lof.log" 2>&1
check_eq "control: --check fails (exit 1) on a table with DSG2's lof edited" "$?" 1
check "control: the failure names DSG2" grep -Eq "differs for .*\bDSG2\b" "${CASE_TMP}/edited-lof.log"
cat "${CASE_TMP}/edited-lof.log"

# Control 2: the first TTN region one base longer.
EDITED="${CASE_TMP}/edited-ttn.tsv"
awk -F'\t' -v OFS='\t' '$1 == "TTN" && match($NF, /^truncating_in chr2:[0-9]+/) {
    $NF = "truncating_in chr2:" (substr($NF, 20, RLENGTH - 19) - 1) substr($NF, RLENGTH + 1) } 1' "$TABLE" > "$EDITED"
# shellcheck disable=SC2016  # the positional arguments of bash -c
check "control: the TTN copy differs from the table" bash -c '! cmp -s "$1" "$2"' _ "$EDITED" "$TABLE"
bash "${REPO}/scripts/ci/build-acmg-lof.sh" --cache "$CACHE" --table "$EDITED" --check > "${CASE_TMP}/edited-ttn.log" 2>&1
check_eq "control: --check fails (exit 1) on a table with a TTN region moved by one base" "$?" 1
check "control: the failure names TTN" grep -Eq "differs for .*\bTTN\b" "${CASE_TMP}/edited-ttn.log"
cat "${CASE_TMP}/edited-ttn.log"

finish
