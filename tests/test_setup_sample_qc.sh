#!/usr/bin/env bash
# setup.sh, install_sample_qc_data (--sample-qc-data, and every plain run):
# somalier's sites and VerifyBamID2's three panel files for step 33. A file
# already there is kept only when its sha256 is the pinned one. One with
# another sha256 (copied by hand, damaged, or from an older pin) gets a
# warning and is downloaded again; offline, that fails the install instead of
# reporting the wrong file as present.
# The function is read out of setup.sh and _digest out of scripts/lib/common.sh,
# and run with a stub fetch that logs each URL and writes a known body, against
# pins made from those bodies. The first two cases are the red controls of the
# wrong-file cases after them: a pass that never downloads must not mask them.
set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
FAILS=0

{
  awk '/^install_sample_qc_data\(\) \{$/ {on = 1} on {print} on && /^}$/ {exit}' "${REPO}/scripts/setup.sh"
  awk '/^_digest\(\) \{$/ {on = 1} on {print} on && /^}$/ {exit}' "${REPO}/scripts/lib/common.sh"
} > "${WORK}/fn.sh"
for f in install_sample_qc_data _digest; do
  if grep -q "^${f}() {" "${WORK}/fn.sh"; then echo "[PASS] ${f} read from the scripts"
  else echo "[FAIL] ${f} not found in setup.sh or lib/common.sh"; exit 1; fi
done

# body URL: what the stub fetch writes for URL, and what the pins below are made of.
body() { printf 'pinned %s\n' "$(basename "$1")"; }
# shellcheck disable=SC1091
. "${WORK}/fn.sh"
sha() { _digest sha256 <(body "$1"); }
U=https://example.invalid/vb2
SITES_URL=https://example.invalid/somalier/sites.hg38.vcf.gz
PANEL=""
for e in UD mu bed; do PANEL="${PANEL:+${PANEL}$'\n'}${U}/1000g.phase3.100k.b38.vcf.gz.dat.${e} $(sha "${U}/1000g.phase3.100k.b38.vcf.gz.dat.${e}")"; done
SITES_SHA=$(sha "$SITES_URL")

G="${WORK}/genome"
SITES="${G}/reference/somalier/sites.hg38.vcf.gz"
UD="${G}/reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.UD"

# run NAME: install_sample_qc_data on G; its output in WORK/NAME.out, the URLs
# fetched in WORK/NAME.calls, its exit code in RC. OFFLINE=1 makes fetch fail.
run() {
  local CALLS="${WORK}/$1.calls"
  : > "$CALLS"
  (
    # The function sourced below calls it.
    # shellcheck disable=SC2329
    fetch() {
      echo "$1" >> "${CALLS}"
      [ -z "${OFFLINE:-}" ] || return 1
      mkdir -p "$(dirname "$2")"
      body "$1" > "$2"
    }
    # shellcheck disable=SC1091
    . "${WORK}/fn.sh"
    # shellcheck disable=SC2034  # read by the sourced function
    GENOME_DIR=$G SOMALIER_SITES_URL=$SITES_URL SOMALIER_SITES_SHA256=$SITES_SHA VB2_PANEL_FILES=$PANEL
    install_sample_qc_data
  ) > "${WORK}/$1.out" 2>&1
  RC=$?
}

check() {  # check DESC TEST...
  local desc=$1; shift
  if "$@"; then echo "[PASS] ${desc}"; else echo "[FAIL] ${desc}"; FAILS=$((FAILS + 1)); fi
}
calls() { grep -c . "${WORK}/$1.calls" || true; }
has() { grep -Eq "$2" "${WORK}/$1.out"; }
lacks() { ! grep -Eq "$2" "${WORK}/$1.out"; }

# --- nothing there: all four are downloaded -------------------------------------
run fresh
check "fresh: exit 0 (got ${RC})" test "$RC" -eq 0
check "fresh: four downloads (got $(calls fresh))" test "$(calls fresh)" -eq 4
check "fresh: the sites are in place" cmp -s "$SITES" <(body "$SITES_URL")

# --- all four there with their pins: nothing is downloaded ------------------------
run kept
check "kept: exit 0 (got ${RC})" test "$RC" -eq 0
check "kept: no download (got $(calls kept))" test "$(calls kept)" -eq 0
check "kept: the sites are reported present with their sha256" has kept '^\[OK\] somalier sites \(step 33\) already present \(sha256 matches\)\.$'

# --- wrong sites already there: a warning, and only they are downloaded again -----
printf 'a hand-copied file\n' > "$SITES"
run wrong-sites
check "wrong sites: exit 0 (got ${RC})" test "$RC" -eq 0
check "wrong sites: not reported as already present" lacks wrong-sites 'already present'
check "wrong sites: a warning names the file" has wrong-sites "^\[WARN\] .*/somalier/sites\.hg38\.vcf\.gz is not the pinned file \(sha256 differs\): downloading it again\.$"
check "wrong sites: one download, the sites (got $(tr '\n' ' ' < "${WORK}/wrong-sites.calls"))" test "$(cat "${WORK}/wrong-sites.calls")" = "$SITES_URL"
check "wrong sites: the pinned sites replace them" cmp -s "$SITES" <(body "$SITES_URL")

# --- a wrong panel file: a warning, and only it is downloaded again ----------------
printf 'an older panel\n' > "$UD"
run wrong-ud
check "wrong .UD: exit 0 (got ${RC})" test "$RC" -eq 0
check "wrong .UD: a warning names the file" has wrong-ud '^\[WARN\] .*/verifybamid2/1000g\.phase3\.100k\.b38\.vcf\.gz\.dat\.UD is not the pinned file \(sha256 differs\): downloading it again\.$'
check "wrong .UD: one download, the .UD (got $(tr '\n' ' ' < "${WORK}/wrong-ud.calls"))" test "$(cat "${WORK}/wrong-ud.calls")" = "${U}/1000g.phase3.100k.b38.vcf.gz.dat.UD"
check "wrong .UD: the pinned .UD replaces it" cmp -s "$UD" <(body "${U}/1000g.phase3.100k.b38.vcf.gz.dat.UD")

# --- wrong sites and no network: the install fails ---------------------------------
printf 'a damaged file\n' > "$SITES"
OFFLINE=1 run offline
check "offline: exit 1 (got ${RC})" test "$RC" -eq 1
check "offline: not reported as already present" lacks offline 'already present'
check "offline: says the sites could not be installed" has offline "^\[WARN\] Could not install somalier's sites\."

echo
if [ "$FAILS" -eq 0 ]; then echo "RESULT: ALL PASS"; else echo "RESULT: ${FAILS} FAILED"; fi
[ "$FAILS" -eq 0 ]
