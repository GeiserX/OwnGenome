#!/usr/bin/env bash
# assets/callability_exons.tsv is what scripts/ci/build-callability.sh writes
# beside the gene table from its pinned downloads (og dev build callability
# reads GENCODE 50's MANE Select CDS records), byte for byte: a hand edit
# fails here. Control: a copy of the committed table with one exon end moved
# must differ from the rebuild. Every gene of the exon table is a catalog
# gene with a position, and each exon name is GENE:exonN. The downloads come
# from the e2e data cache when the runner filled it (entry catalogs, mirror
# first; scripts/ci/e2e-run.sh), else from their sources: the case runs after
# callability-1-catalog (its needs in scripts/ci/e2e-case-times.tsv), whose
# run fills that entry.
. "$(dirname "$0")/lib.sh"

OUT_DIR="${CASE_TMP}/build"
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"
COMMITTED="${REPO}/assets/callability_exons.tsv"
GENES="${REPO}/assets/callability_genes.tsv"
echo "+ scripts/ci/build-callability.sh"
CACHE="${CASE_TMP}/downloads"
[ -z "${E2E_DATA_CACHE:-}" ] || CACHE="${E2E_DATA_CACHE}/catalogs/callability"
CALLABILITY_CACHE="$CACHE" bash "${REPO}/scripts/ci/build-callability.sh" "${OUT_DIR}/callability_genes.tsv" 2>&1 | tee "$STEP_LOG"
check_eq "build-callability.sh exits 0" "${PIPESTATUS[0]}" 0
OUT="${OUT_DIR}/callability_exons.tsv"
check "the build wrote the exon table beside the gene table" test -s "$OUT"

awk -F'\t' -v OFS='\t' '!/^#/ && $1 != "exon" && !done {$6 = $6 + 1; done = 1} 1' \
  "$COMMITTED" > "${CASE_TMP}/edited.tsv"
if cmp -s "$OUT" "${CASE_TMP}/edited.tsv"; then
  fail "control: the rebuild matches a hand-edited exon table, so this check cannot fail"
else
  pass "control: the rebuild differs from a hand-edited exon table"
fi
if cmp -s "$OUT" "$COMMITTED"; then
  pass "the rebuilt exon table is the committed assets/callability_exons.tsv, byte for byte"
else
  fail "the rebuilt exon table differs from assets/callability_exons.tsv; rebuild it with scripts/ci/build-callability.sh"
  diff "$COMMITTED" "$OUT" | head -20
fi
check "the rebuilt gene table is still the committed one" cmp -s "${OUT_DIR}/callability_genes.tsv" "$GENES"

check_ge "the exon table has 10,000 or more coding exons" "$(awk -F'\t' '!/^#/ && $1 != "exon"' "$COMMITTED" | wc -l | tr -d ' ')" 10000
check_eq "every exon is named GENE:exonN after its gene" \
  "$(awk -F'\t' '!/^#/ && $1 != "exon" && $1 !~ ("^" $2 ":exon[0-9]+$")' "$COMMITTED" | wc -l | tr -d ' ')" 0
check_eq "every gene of the exon table is a catalog gene with a position, on its chromosome" \
  "$(awk -F'\t' 'NR == FNR {if (!/^#/ && $4 ~ /^[0-9]+$/) c[$1] = $3; next}
     !/^#/ && $1 != "exon" && c[$2] != $4' "$GENES" "$COMMITTED" | wc -l | tr -d ' ')" 0

finish
