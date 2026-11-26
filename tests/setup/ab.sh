#!/usr/bin/env bash
# tests/setup/ab.sh: install one setup component twice, with bash (side A)
# and with og (side B), and compare everything the two runs left: the
# install trees file by file (og dev install-manifest), the exit codes, the
# logs with times masked, the missing list, the images and pinned tags, the
# Nextflow plugin folder and bin/.og.
#
# Usage: tests/setup/ab.sh [--bash-vs-bash] COMPONENT [WORK_DIR]
#
# COMPONENT is a setup.sh flag without its dashes (kir-data, viral-data,
# sample-qc-data, archaic, ...), or refresh-clinvar for --refresh clinvar.
#   side A: OG_SETUP=bash scripts/setup.sh --COMPONENT <install folder>
#   side B: OG_SETUP=og og setup install --COMPONENT <install folder>
# Side B calls og directly and never setup.sh, so it cannot fall back to
# bash; when og has no `setup install`, ab.sh stops before side A.
# --bash-vs-bash runs side A's command on both sides: it finds the bytes
# that change from one run to the next (the allowlist's rows), and its
# result says no og side ran.
#
# Both sides install into the same empty folder, WORK_DIR/install, moved to
# WORK_DIR/a or WORK_DIR/b afterwards, so a path written into a file is the
# same on both. WORK_DIR must not hold install, a or b; it defaults to a new
# temporary folder and is kept, with each side's log, exit code and
# manifests, and the compares.
#
# Environment (anything else exported, such as ARCHAIC_REGIONS, reaches
# both sides):
#   AB_SEED       a folder put into each side's empty install folder first,
#                 for a component that needs another one (the reference):
#                 hard links on Linux (a copy on another filesystem), APFS
#                 clones on macOS, never symlinks, which a container cannot
#                 follow out of GENOME_DIR. The seed must still be the same
#                 afterwards.
#   AB_ALLOW      the allowlist (default tests/setup/install-manifest-allow.tsv)
#   OG            the og binary (default <repo>/bin/og)
#   OG_DOWNLOADS  the downloads folder both sides share (default
#                 WORK_DIR/downloads), so side B reuses side A's downloads
#                 where setup reads it.
#   CONTAINER_ENGINE  docker (default) or podman
# Each side gets its own empty NXF_HOME (WORK_DIR/a.nxf, b.nxf).
#
# Every run ends with a control: one bit of the smallest file compared byte
# for byte is flipped in tree B, and the compare must go red and name it.
#
# Exit 0 when every check is green, 1 when one is red, 2 when an argument is
# wrong, side B has no og to run, or the control did not go red.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)

usage() {
  sed -n '9p' "$0" | sed 's/^# //' >&2
  exit 2
}

bvb=false
if [ "${1:-}" = --bash-vs-bash ]; then
  bvb=true
  shift
fi
[ $# -ge 1 ] && [ $# -le 2 ] || usage
comp=$1
case "$comp" in '' | -* | */* | *' '*) usage ;; esac
if [ "$comp" = refresh-clinvar ]; then words=(--refresh clinvar); else words=("--${comp}"); fi
if [ $# -eq 2 ]; then
  work=$2
  mkdir -p "$work"
else
  work=$(mktemp -d "${TMPDIR:-/tmp}/og-ab.XXXXXX")
fi
work=$(cd "$work" && pwd)
for d in install a b; do
  if [ -e "${work}/${d}" ]; then
    echo "ERROR: ${work}/${d} exists; give a WORK_DIR without install, a and b." >&2
    exit 2
  fi
done

OG=${OG:-${ROOT}/bin/og}
ENGINE=${CONTAINER_ENGINE:-docker}
allow=${AB_ALLOW:-${ROOT}/tests/setup/install-manifest-allow.tsv}
export OG_DOWNLOADS=${OG_DOWNLOADS:-${work}/downloads}
mkdir -p "$OG_DOWNLOADS"
if [ -n "${AB_SEED:-}" ]; then
  [ -d "$AB_SEED" ] || { echo "ERROR: AB_SEED=${AB_SEED} is not a folder." >&2; exit 2; }
  AB_SEED=$(cd "$AB_SEED" && pwd)
fi

side_a=(env OG_SETUP=bash "${ROOT}/scripts/setup.sh" "${words[@]}")
if $bvb; then
  side_b=("${side_a[@]}")
else
  side_b=(env OG_SETUP=og "$OG" setup install "${words[@]}")
  # og setup prints its usage for any --help; only a real install verb names itself.
  help=$("$OG" setup install --help 2>&1) || help=""
  if [[ "$help" != *"og setup install"* ]]; then
    echo "ERROR: side B runs ${OG} setup install, and that does not run here. ab.sh never runs bash in its place;" >&2
    echo "       build og from this checkout (scripts/setup.sh --og), or pass --bash-vs-bash to compare bash with bash." >&2
    exit 2
  fi
fi

manifest() { "$OG" dev install-manifest "$@"; }

# images: the engine's tags and image ids, sorted; "unavailable" without one.
images() {
  "$ENGINE" image ls --digests --format '{{.Repository}}:{{.Tag}}\t{{.ID}}' 2>/dev/null | LC_ALL=C sort \
    || echo unavailable
}

# seed DIR: the seed's content into the empty folder DIR, never as symlinks.
seed() {
  [ -n "${AB_SEED:-}" ] || return 0
  if cp --version >/dev/null 2>&1; then
    cp -al "${AB_SEED}/." "$1/" 2>/dev/null || { rm -rf "$1" && mkdir "$1" && cp -pR "${AB_SEED}/." "$1/"; }
  else
    cp -cpR "${AB_SEED}/." "$1/" 2>/dev/null || { rm -rf "$1" && mkdir "$1" && cp -pR "${AB_SEED}/." "$1/"; }
  fi
}

# Untag the localhost/pgp-pinned tags side A added, where the image keeps
# another tag (so only the tag goes): side B must add them again. An image
# side A pulled stays; pulling it again would only cost time.
untag_pinned() {
  local tag id n=0
  while IFS=$'\t' read -r tag id; do
    case "$tag" in localhost/pgp-pinned/* | pgp-pinned/*) ;; *) continue ;; esac
    if [ "$(awk -F'\t' -v i="$id" '$2 == i' "${work}/a.images" | wc -l)" -gt 1 ] \
      && "$ENGINE" rmi "$tag" >/dev/null 2>&1; then
      n=$((n + 1))
    fi
  done < <(LC_ALL=C comm -13 "${work}/before.images" "${work}/a.images")
  echo "   untagged ${n} pinned tag(s) side A added, for side B to add again"
}

images >"${work}/before.images"
og_dir="${ROOT}/bin/.og"
[ -n "${AB_SEED:-}" ] && manifest "$AB_SEED" >"${work}/seed.before.tsv" 2>/dev/null

for side in a b; do
  if [ "$side" = a ]; then cmd=("${side_a[@]}"); else cmd=("${side_b[@]}"); fi
  [ "$side" = b ] && untag_pinned
  mkdir "${work}/install" "${work}/${side}.nxf"
  seed "${work}/install"
  echo "== side ${side}: ${cmd[*]} ${work}/install"
  start=$SECONDS
  rc=0
  NXF_HOME="${work}/${side}.nxf" "${cmd[@]}" "${work}/install" >"${work}/${side}.log" 2>&1 || rc=$?
  echo "$rc" >"${work}/${side}.exit"
  echo "   exit ${rc} after $((SECONDS - start)) s"
  mv "${work}/install" "${work}/${side}"
  manifest "${work}/${side}" >"${work}/${side}.tsv" 2>/dev/null
  manifest "${work}/${side}.nxf" >"${work}/${side}.nxf.tsv" 2>/dev/null
  images >"${work}/${side}.images"
  if [ -d "$og_dir" ]; then manifest "$og_dir" >"${work}/${side}.og.tsv" 2>/dev/null; else echo none >"${work}/${side}.og.tsv"; fi
done

red=0
check() { # NAME COMMAND...: run a check, print its verdict, keep its output.
  local name=$1 rc=0
  shift
  "$@" >"${work}/check.${name}.txt" 2>&1 || rc=$?
  if [ "$rc" -eq 0 ]; then
    echo "[OK]   ${name}"
  else
    red=$((red + 1))
    echo "[RED]  ${name} (${work}/check.${name}.txt):"
    grep -v '^#' "${work}/check.${name}.txt" | head -n 20 | sed 's/^/         /'
  fi
}
# shellcheck disable=SC2329  # called through check
same_file() { cmp "$1" "$2" || { diff "$1" "$2" | head -n 20; return 1; }; }

check exit-codes same_file "${work}/a.exit" "${work}/b.exit"
check tree manifest --compare "${work}/a.tsv" "${work}/b.tsv" --allow "$allow"
check nxf-home manifest --compare "${work}/a.nxf.tsv" "${work}/b.nxf.tsv" --allow "$allow"
check logs manifest --compare-logs "${work}/a.log" "${work}/b.log"
check images same_file "${work}/a.images" "${work}/b.images"
check bin-og same_file "${work}/a.og.tsv" "${work}/b.og.tsv"
if [ -n "${AB_SEED:-}" ]; then
  manifest "$AB_SEED" >"${work}/seed.after.tsv" 2>/dev/null
  check seed-unchanged manifest --compare "${work}/seed.before.tsv" "${work}/seed.after.tsv"
fi
grep '^ALLOWED' "${work}/check.tree.txt" "${work}/check.nxf-home.txt" 2>/dev/null | sed 's/^/         /' || true
if [ "$(cat "${work}/a.exit")" != 0 ]; then
  echo "[WARN] side A exited $(cat "${work}/a.exit"): the trees compared may be partial; read ${work}/a.log."
fi

# The control: one flipped bit in tree B must turn the tree compare red.
path=$(manifest --control "${work}/b.tsv" --allow "$allow") || {
  echo "ERROR: tree B has no file compared byte for byte to plant the control in." >&2
  exit 2
}
# The file may be a hard link into the seed: give B its own copy first.
real="${work}/b/$(printf '%b' "$path")"
cp -p "$real" "${real}.ab-copy" && mv -f "${real}.ab-copy" "$real"
manifest --flip "${work}/b" "$path" >/dev/null
manifest "${work}/b" >"${work}/b.control.tsv" 2>/dev/null
crc=0
manifest --compare "${work}/a.tsv" "${work}/b.control.tsv" --allow "$allow" >"${work}/check.control.txt" || crc=$?
manifest --flip "${work}/b" "$path" >/dev/null
if [ "$crc" -ne 1 ] || ! grep -qF "DIFF	content	${path}	" "${work}/check.control.txt"; then
  echo "ERROR: one flipped bit in ${path} did not turn the compare red (exit ${crc}): the compare cannot fail." >&2
  exit 2
fi
echo "[OK]   control: one flipped bit in ${path} is red, and restored"

if $bvb; then label="bash vs bash (no og side ran)"; else label="bash vs og"; fi
if [ "$red" -eq 0 ]; then
  echo "[OK] ${comp}, ${label}: every check is green. Files: ${work}"
  exit 0
fi
echo "[FAIL] ${comp}, ${label}: ${red} check(s) red. Files: ${work}"
exit 1
