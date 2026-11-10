#!/usr/bin/env bash
# setup.sh --masked-reference: the pins in versions.env and the masking
# program, without Docker or the network.
#   1. GRC_EXCLUSIONS_SHA256 is the sha256 of the GRC's list as GIAB
#      publishes it (the 12 rows below, read from the file on 2026-11-04);
#      the same check on a copy of versions.env with one changed digit fails.
#      GRC_MASKED_REF_SHA256 is a sha256.
#   2. setup.sh masks with og (og setup mask-fasta, run from the reference's
#      directory with relative paths, og built first), and the Python it
#      replaced is gone from it. With OG set to an og binary (the parity job
#      of guard.yml sets it), og runs on a small FASTA: every base of a row
#      is N, every other byte is the same (line breaks, lower case, the
#      header and a sequence with no row), and it names the row on a
#      sequence the FASTA lacks; the same run with the BED shifted by one
#      base gives another file, so the comparison can fail. The parity cases
#      (tests/parity/cases.tsv) compare og with that Python on these inputs
#      and more. The e2e case fixitfelix-chr21-1 runs the whole option on
#      chr21.
set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
FAILS=0
pass() { echo "[PASS] $*"; }
fail() { echo "[FAIL] $*"; FAILS=$((FAILS + 1)); }
sha() { if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }

# --- 1. the pins ---------------------------------------------------------------------
printf '%s\t%s\t%s\n' \
  chr9 70701806 70719798  chr13 86252979 86269524  chr13 111753855 111789114 \
  chr16 34339329 34456556  chr21 5010000 5166246  chr21 5966593 6161371 \
  chr21 6427258 6580181  chr21 6789085 6934219  chr21 7743700 7865746 \
  chrX 37085011 37099262  chr11_KI270721v1_random 1 100316  chr22_KI270734v1_random 1 165050 \
  > "${WORK}/grc.bed"
GRC_SHA=$(sha "${WORK}/grc.bed")

# pin_ok ENV: GRC_EXCLUSIONS_SHA256 in ENV is the sha256 of the 12 rows.
pin_ok() { [ "$(sed -n 's/^GRC_EXCLUSIONS_SHA256="\([0-9a-f]*\)"$/\1/p' "$1")" = "$GRC_SHA" ]; }
if pin_ok "${REPO}/versions.env"; then pass "GRC_EXCLUSIONS_SHA256 is the sha256 of the GRC list (${GRC_SHA})"
else fail "GRC_EXCLUSIONS_SHA256 in versions.env is not ${GRC_SHA}"; fi
sed 's/^\(GRC_EXCLUSIONS_SHA256="\)./\1x/' "${REPO}/versions.env" > "${WORK}/changed.env"
if ! pin_ok "${WORK}/changed.env" && ! cmp -s "${WORK}/changed.env" "${REPO}/versions.env"; then
  pass "control: a changed GRC_EXCLUSIONS_SHA256 fails the same check"
else fail "control: a changed GRC_EXCLUSIONS_SHA256 passed"; fi
if grep -Eq '^GRC_MASKED_REF_SHA256="[0-9a-f]{64}"$' "${REPO}/versions.env"; then pass "GRC_MASKED_REF_SHA256 is a sha256"
else fail "GRC_MASKED_REF_SHA256 is not a sha256"; fi
if grep -q '^GRC_EXCLUSIONS_URL="https://.*/GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed"$' "${REPO}/versions.env"; then
  pass "GRC_EXCLUSIONS_URL names GIAB's copy of the list"
else fail "GRC_EXCLUSIONS_URL is not GIAB's T2Tv2 list"; fi

# --- 2. the masking program ----------------------------------------------------------
SETUP="${REPO}/scripts/setup.sh"
# mask_checks FILE: the function body of install_masked_reference in FILE,
# checked three ways; prints one word per check that holds: og, built, nopy.
mask_checks() {
  local body
  body=$(awk '/^install_masked_reference\(\) \{$/ {on = 1} on {print} on && /^}$/ {exit}' "$1")
  if grep -q '/pgp-bin/og setup mask-fasta "$(basename "$src")" "$(basename "$bed")" "$(basename "$dest").part"' <<< "$body" \
     && grep -q -- '-w "$(cpath "$dir")"' <<< "$body"; then echo og; fi
  if grep -q 'rust/og.sh" check --required' <<< "$body" && grep -q 'rust/og.sh" install' <<< "$body"; then echo built; fi
  if ! grep -Eq 'MASK_FASTA_PY|python3 -c' <<< "$body" && ! grep -q '^MASK_FASTA_PY=' "$1"; then echo nopy; fi
}
GOT=$(mask_checks "$SETUP")
if grep -qx og <<< "$GOT"; then
  pass "setup.sh masks with og setup mask-fasta, from the reference's directory with relative paths"
else fail "setup.sh does not call og setup mask-fasta from the reference's directory"; fi
if grep -qx built <<< "$GOT"; then
  pass "setup.sh builds og before it masks, when og is missing or stale"
else fail "setup.sh masks without making sure og is built"; fi
if grep -qx nopy <<< "$GOT"; then pass "the Python masking program is gone from setup.sh"
else fail "the Python masking program is still in setup.sh"; fi
# control: the same checks on setup.sh as it was before og (the commit
# tests/parity/oracle.ref names), cut to the lines that matter
cat > "${WORK}/old-setup.sh" <<'OLD'
install_masked_reference() {
  local src=$REF_FASTA dest bed want got dir log
    rm -f "${dest}.part"
    if ! log=$(run_in --rw "$dir" "$PYTHON_IMAGE" python3 -c "$MASK_FASTA_PY" \
          "$(cpath "$src")" "$(cpath "$bed")" "$(cpath "${dest}.part")"); then
      rm -f "${dest}.part"
      return 1
    fi
    mv -f "${dest}.part" "$dest"
}
MASK_FASTA_PY=$(cat <<'PY'
PY
)
OLD
OLD_GOT=$(mask_checks "${WORK}/old-setup.sh")
if ! grep -qx og <<< "$OLD_GOT" && ! grep -qx nopy <<< "$OLD_GOT"; then
  pass "control: on setup.sh before og, 'masks with og' fails and the Python program is found"
else fail "control: the checks passed setup.sh before og (${OLD_GOT//$'\n'/ })"; fi

if [ -n "${OG:-}" ]; then
  printf '>s1 first\nACGTacg\ntACGTAC\nGT\n>s2\nAAAAAAA\nCC\n' > "${WORK}/in.fa"
  printf 'browser position s1:1-3\n#seq\tstart\tend\ns1\t5\t9\ns1\t15\t16\nmissing\t0\t10\n' > "${WORK}/rows.bed"
  # s1 is ACGTacgtACGTACGT: 0-based [5, 9) and [15, 16) are N.
  printf '>s1 first\nACGTaNN\nNNCGTAC\nGN\n>s2\nAAAAAAA\nCC\n' > "${WORK}/want.fa"
  OUT=$("$OG" setup mask-fasta "${WORK}/in.fa" "${WORK}/rows.bed" "${WORK}/out.fa" 2>&1)
  echo "  ${OUT}"
  if cmp -s "${WORK}/out.fa" "${WORK}/want.fa"; then pass "og: the rows are N and every other byte is the same"
  else fail "og: masked FASTA differs from the expected one"; diff "${WORK}/want.fa" "${WORK}/out.fa"; fi
  case "$OUT" in
    "masked 2 of 3 rows (5 bases); rows on sequences the FASTA lacks: missing") pass "og reports 2 of 3 rows and names the missing sequence" ;;
    *) fail "og report: '${OUT}'" ;;
  esac
  printf 's1\t6\t10\ns1\t15\t16\n' > "${WORK}/shifted.bed"
  "$OG" setup mask-fasta "${WORK}/in.fa" "${WORK}/shifted.bed" "${WORK}/shifted.fa" > /dev/null 2>&1
  if [ -s "${WORK}/shifted.fa" ] && ! cmp -s "${WORK}/shifted.fa" "${WORK}/want.fa"; then pass "control: a BED one base off gives another file"
  else fail "control: a BED one base off gave the same file"; fi
else
  echo "  OG is not set: og itself is checked by the parity job (OG=... bash tests/test_fixitfelix.sh)"
fi

if [ "$FAILS" -gt 0 ]; then echo "test_fixitfelix: ${FAILS} check(s) failed"; exit 1; fi
echo "test_fixitfelix: all checks passed"
