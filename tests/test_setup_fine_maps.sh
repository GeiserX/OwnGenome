#!/usr/bin/env bash
# setup.sh's install_fine_maps (step 26b's genetic maps) rebuilds a map an
# older install left cut short. Before og, setup.sh skipped every chromosome
# whose .rfmix.map was not empty, so a map a kill had cut short stayed
# forever. Now maps/DONE, written last, names the map and the chromosomes,
# and without it every map is written again. The test runs the function out
# of setup.sh with og on the host in place of the container (run_in runs the
# command after the image name, with /pgp-bin/og pointing at OG):
#   1. A folder as an install from before og leaves it: the zip, chr22.map
#      whole and chr22.rfmix.map cut short, no maps/DONE. The function
#      rewrites both maps, so they equal a fresh install's, and writes
#      maps/DONE.
#   2. With maps/DONE in place a second call runs og no more (a planted
#      byte stays), so a finished install is not redone.
#   3. A maps/DONE of another chromosome list counts as none: the maps are
#      written again.
#   4. The control: the rule of before (skip a non-empty .rfmix.map), run on
#      the folder of 1, keeps the cut-short map, so check 1 can fail.
# Run: OG=bin/og bash tests/test_setup_fine_maps.sh (skips without og; with
# OG_REQUIRED set, the skip is a failure, as guard.yml's rust job sets it).
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
OG=${OG:-${REPO}/bin/og}
if ! "$OG" version > /dev/null 2>&1; then
  if [ -n "${OG_REQUIRED:-}" ]; then echo "FAIL: og is not built (scripts/setup.sh --og)"; exit 1; fi
  echo "skip: og is not built for this machine (scripts/setup.sh --og)"
  exit 0
fi
PY=$(command -v python3)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
FAILS=0
ok() { echo "ok: $1"; }
bad() { echo "FAIL: $1"; FAILS=$((FAILS + 1)); }

# The function, as setup.sh holds it, and the stand-ins it calls.
awk '/^install_fine_maps\(\) \{$/ {on = 1} on {print} on && /^}$/ {exit}' \
  "${REPO}/scripts/setup.sh" > "${WORK}/fn.sh"
grep -q '^install_fine_maps() {$' "${WORK}/fn.sh" || { echo "FAIL: install_fine_maps is not in scripts/setup.sh"; exit 1; }
# shellcheck source=/dev/null
. "${WORK}/fn.sh"
# shellcheck disable=SC2034  # PGP_ROOT is read by install_fine_maps
PYTHON_IMAGE=python-image PGP_ROOT=$REPO ANCESTRY_FINE_MAP_SHA256=0123abcd
OG_CALLS="${WORK}/og.calls"
need_og() { return 0; }
cpath() { printf '%s\n' "$1"; }
run_in() {
  while [ "$1" != "$PYTHON_IMAGE" ]; do shift; done
  shift
  [ "$1" = /pgp-bin/og ] || { echo "run_in: unexpected command $1" >&2; return 1; }
  shift
  echo "og $*" >> "$OG_CALLS"
  "$OG" "$@"
}

# A zip with two chromosomes, laid out as Beagle's plink.GRCh38.map.zip
# (its members are named plink.chrchrN).
mkzip() {
  "$PY" - "$1" <<'PY'
import sys, zipfile
with zipfile.ZipFile(sys.argv[1], "w", zipfile.ZIP_DEFLATED) as z:
    for c, n in (("chr21", 300), ("chr22", 400)):
        rows = "".join(f"{c}\t.\t{i * 0.0123:.6f}\t{10000000 + i * 997}\n" for i in range(n))
        z.writestr(f"chr_in_chrom_field/plink.chr{c}.GRCh38.map", rows)
PY
}

# The fresh install every other folder is checked against.
FRESH="${WORK}/fresh"
mkdir -p "${FRESH}/maps"
mkzip "${FRESH}/maps/plink.GRCh38.map.zip"
install_fine_maps "$FRESH" "chr21 chr22 " > "${WORK}/fresh.log" 2>&1 || { echo "FAIL: the fresh install failed: $(cat "${WORK}/fresh.log")"; exit 1; }
for f in chr21.map chr21.rfmix.map chr22.map chr22.rfmix.map; do
  [ -s "${FRESH}/maps/${f}" ] || { echo "FAIL: the fresh install wrote no ${f}"; exit 1; }
done

# old_install DIR: the folder an install from before og left after a kill
# in the middle of chr22.rfmix.map.
old_install() {
  mkdir -p "$1/maps"
  cp "${FRESH}/maps/plink.GRCh38.map.zip" "${FRESH}/maps/chr21.map" "${FRESH}/maps/chr21.rfmix.map" "${FRESH}/maps/chr22.map" "$1/maps/"
  head -c 200 "${FRESH}/maps/chr22.rfmix.map" > "$1/maps/chr22.rfmix.map"
}
same_maps() {
  local f
  for f in chr21.map chr21.rfmix.map chr22.map chr22.rfmix.map; do cmp -s "${FRESH}/maps/${f}" "$1/maps/${f}" || return 1; done
}

# 1. The cut-short map is written again.
OLD="${WORK}/old"
old_install "$OLD"
rc=0
install_fine_maps "$OLD" "chr21 chr22 " > "${WORK}/old.log" 2>&1 || rc=$?
if [ "$rc" -ne 0 ]; then
  bad "install_fine_maps failed on the old folder (exit ${rc}): $(tr '\n' ' ' < "${WORK}/old.log")"
elif ! same_maps "$OLD"; then
  bad "the cut-short chr22.rfmix.map was not written again ($(wc -c < "${OLD}/maps/chr22.rfmix.map") bytes, want $(wc -c < "${FRESH}/maps/chr22.rfmix.map"))"
else
  ok "a cut-short chr22.rfmix.map with no maps/DONE is written again, byte for byte as a fresh install"
fi
if [ "$(cat "${OLD}/maps/DONE" 2>/dev/null)" = "map=${ANCESTRY_FINE_MAP_SHA256} chroms=chr21 chr22" ]; then
  ok "maps/DONE names the map and the chromosomes"
else
  bad "maps/DONE is '$(cat "${OLD}/maps/DONE" 2>/dev/null || echo missing)'"
fi

# 2. A finished install is not redone.
: > "$OG_CALLS"
printf 'x' >> "${OLD}/maps/chr21.map"
install_fine_maps "$OLD" "chr21 chr22 " > /dev/null 2>&1 || bad "the second call failed"
if [ ! -s "$OG_CALLS" ] && [ "$(tail -c 1 "${OLD}/maps/chr21.map")" = x ]; then
  ok "with maps/DONE in place a second call runs og no more"
else
  bad "the second call ran og again: $(tr '\n' ' ' < "$OG_CALLS")"
fi

# 3. A maps/DONE of another chromosome list counts as none.
: > "$OG_CALLS"
echo "map=${ANCESTRY_FINE_MAP_SHA256} chroms=chr21" > "${OLD}/maps/DONE"
install_fine_maps "$OLD" "chr21 chr22 " > /dev/null 2>&1 || bad "the call with a stale maps/DONE failed"
if [ -s "$OG_CALLS" ] && same_maps "$OLD"; then
  ok "a maps/DONE of other chromosomes makes og write the maps again"
else
  bad "a stale maps/DONE was taken as done"
fi

# 4. The control: the rule of before keeps the cut-short map.
CTRL="${WORK}/control"
old_install "$CTRL"
for c in chr21 chr22; do
  [ -s "${CTRL}/maps/${c}.rfmix.map" ] && continue
  run_in --rw "$CTRL" "$PYTHON_IMAGE" /pgp-bin/og setup genetic-map --zip "${CTRL}/maps/plink.GRCh38.map.zip" --out "${CTRL}/maps" --rfmix "$c" > /dev/null
done
if same_maps "$CTRL"; then
  bad "control: the rule of before also left whole maps, so check 1 cannot fail"
else
  ok "control: the rule of before (skip a non-empty .rfmix.map) keeps the cut-short map"
fi

if [ "$FAILS" -ne 0 ]; then
  echo "${FAILS} check(s) failed"
  exit 1
fi
echo "all checks passed"
