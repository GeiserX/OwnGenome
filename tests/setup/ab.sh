#!/usr/bin/env bash
# tests/setup/ab.sh: install one setup component twice, with bash (side A)
# and with og (side B), and compare everything the two runs left: the
# install trees file by file (og dev install-manifest), the exit codes, the
# logs with times masked, the missing list, the images and pinned tags each
# side added or removed, the Nextflow plugin folder, bin/.og, and who owns
# what each side wrote.
#
# Usage: tests/setup/ab.sh [--bash-vs-bash] COMPONENT [WORK_DIR]
#
# COMPONENT is a setup.sh flag without its dashes (kir-data, viral-data,
# sample-qc-data, archaic, ...), or refresh-clinvar for --refresh clinvar.
#   side A: OG_SETUP=bash scripts/setup.sh --COMPONENT <install folder>
#   side B: OG_SETUP=og og setup install --COMPONENT <install folder>
# Side B calls og directly and never setup.sh, so it cannot fall back to
# bash; when og has no `setup install`, ab.sh stops before side A, so build
# an og that has it first. --bash-vs-bash runs side A's command on both
# sides: it finds the bytes that change from one run to the next (the
# allowlist's rows), and its result says no og side ran.
#
# Both sides install into the same empty folder, WORK_DIR/install, moved to
# WORK_DIR/a or WORK_DIR/b afterwards, so a path written into a file is the
# same on both. WORK_DIR must not hold install, a or b; it defaults to a new
# temporary folder and is kept, with each side's log, exit code, manifests
# (and their errors, NAME.manifest.err), and the compares.
#
# Environment (anything else exported, such as ARCHAIC_REGIONS, reaches
# both sides):
#   AB_SEED       a folder put into each side's empty install folder first,
#                 for a component that needs another one (the reference):
#                 hard links on Linux (a copy on another filesystem), APFS
#                 clones on macOS, never symlinks, which a container cannot
#                 follow out of GENOME_DIR. Everything in it must belong to
#                 the user running ab.sh, and it must still be the same
#                 afterwards.
#   AB_EXPECT_EXIT  the exit code side A must end with (default 0), for a
#                 component that is meant to stop short
#   AB_ALLOW      the allowlist (default tests/setup/install-manifest-allow.tsv)
#   OG            the og binary (default <repo>/bin/og)
#   OG_DOWNLOADS  the downloads folder both sides share (default
#                 WORK_DIR/downloads), so side B reuses side A's downloads
#                 where setup reads it.
#   CONTAINER_ENGINE  docker (default) or podman
# Each side gets its own empty NXF_HOME (WORK_DIR/a.nxf, b.nxf).
#
# Tree A must hold at least one entry the seed did not, or the run proves
# nothing. Every run ends with a control: one bit of the smallest file the
# component wrote and the compare reads byte for byte is flipped in tree B,
# and the compare must go red and name it.
#
# Exit 0 when every check is green, 1 when one is red, 2 when an argument is
# wrong, side B has no og to run, a manifest cannot be taken, tree A holds
# nothing new, or the control did not go red.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)

usage() {
  grep -m 1 '^# Usage:' "$0" | sed 's/^# //' >&2
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
expect=${AB_EXPECT_EXIT:-0}
case "$expect" in '' | *[!0-9]*)
  echo "ERROR: AB_EXPECT_EXIT=${expect} is not an exit code." >&2
  exit 2
  ;;
esac
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
uid=$(id -u)
allow=${AB_ALLOW:-${ROOT}/tests/setup/install-manifest-allow.tsv}
export OG_DOWNLOADS=${OG_DOWNLOADS:-${work}/downloads}
mkdir -p "$OG_DOWNLOADS"
if [ -n "${AB_SEED:-}" ]; then
  [ -d "$AB_SEED" ] || { echo "ERROR: AB_SEED=${AB_SEED} is not a folder." >&2; exit 2; }
  AB_SEED=$(cd "$AB_SEED" && pwd)
  # Seeded files keep their owner: one that is not ours would read as a
  # file the component wrote as another user.
  if [ -n "$(find "$AB_SEED" ! -user "$uid" | head -n 1)" ]; then
    echo "ERROR: AB_SEED holds entries uid ${uid} does not own; run ab.sh as the seed's owner." >&2
    exit 2
  fi
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

# take NAME DIR: the manifest of DIR into WORK/NAME.tsv. Its messages go to
# WORK/NAME.manifest.err, shown when it fails (an unreadable file, say).
take() {
  manifest "$2" >"${work}/$1.tsv" 2>"${work}/$1.manifest.err" || {
    echo "ERROR: og dev install-manifest $2 failed (${work}/$1.manifest.err):" >&2
    sed 's/^/       /' "${work}/$1.manifest.err" >&2
    exit 2
  }
}

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

# named REF: whether side A's log names REF whole, not as part of a longer
# name, tag or path (repo:v1 names neither repo:v10 nor reg.example/repo:v1).
named() {
  awk -v t="$1" '
    {
      s = $0; off = 0
      while ((i = index(s, t)) > 0) {
        a = (off + i > 1) ? substr($0, off + i - 1, 1) : ""
        b = substr($0, off + i + length(t), 1)
        if (a !~ "[[:alnum:]_./-]" && b !~ "[[:alnum:]_./-]") { found = 1; exit }
        off += i; s = substr(s, i + 1)
      }
    }
    END { exit !found }' "${work}/a.log"
}

# undo_images: remove the tags side A added for an image its log names (an
# image with no other tag goes with its last one), so side B starts from the
# images side A started from: it pulls and pins what side A pulled and
# pinned, and logs it the same way. A tag another job on the engine added
# while side A ran is not named in side A's log and is left alone. A tag is
# named by its full NAME:TAG, a pinned tag by its digest, never by its name
# alone: another tag of the same image name is not side A's.
undo_images() {
  local tag ref n=0
  while IFS=$'\t' read -r tag _; do
    case "$tag" in
      */pgp-pinned/*:sha256-*) ref="sha256:${tag##*:sha256-}" ;;
      *) ref=$tag ;;
    esac
    named "$ref" || continue
    "$ENGINE" rmi "$tag" >/dev/null 2>&1 && n=$((n + 1))
  done < <(LC_ALL=C comm -13 "${work}/a.before.images" "${work}/a.images" | grep -v ':<none>' || true)
  echo "   removed ${n} image tag(s) side A added, for side B to add again"
}

# delta SIDE: what the side changed in the image list, "+" added and "-"
# removed; another job on the machine may change it too, but only while a
# side runs.
delta() {
  LC_ALL=C comm -3 "${work}/$1.before.images" "${work}/$1.images" \
    | awk -F'\t' '{ if ($1 != "") print "-\t" $1 "\t" $2; else print "+\t" $2 "\t" $3 }' >"${work}/$1.images.delta"
}

# One run per checkout: each run snapshots bin/.og before side A and puts it
# back before side B, so a second run would snapshot or delete files the
# first one is using.
lock="${ROOT}/bin/.og.ab-lock"
mkdir -p "${ROOT}/bin"
if ! mkdir "$lock" 2>/dev/null; then
  echo "ERROR: another ab.sh run uses this checkout (${lock} exists, pid $(cat "${lock}/pid" 2>/dev/null || echo unknown)). Wait for it, run from another checkout, or remove the folder if no run is left." >&2
  exit 2
fi
echo $$ >"${lock}/pid"
trap 'rm -rf "$lock"' EXIT

# setup.sh builds og when it is missing or stale, and only the first of two
# runs would: build it before side A, so neither side does.
if ! "${ROOT}/rust/og.sh" check --required >/dev/null 2>&1; then
  echo "== og is missing or stale: rust/og.sh install first (log in ${work}/og-build.log)"
  "${ROOT}/rust/og.sh" install >"${work}/og-build.log" 2>&1 || {
    echo "ERROR: could not build og; read ${work}/og-build.log." >&2
    exit 2
  }
fi
og_dir="${ROOT}/bin/.og"
seed_tsv="${work}/seed.none.tsv"
if [ -n "${AB_SEED:-}" ]; then
  take seed.before "$AB_SEED"
  seed_tsv="${work}/seed.before.tsv"
else
  : >"$seed_tsv"
fi
# Both sides start from this bin/.og (og itself included), so side B cannot
# inherit a file side A added there.
if [ -d "$og_dir" ]; then cp -pR "$og_dir" "${work}/og.before"; fi
# og_restore: bin/.og as it was before side A.
og_restore() {
  rm -rf "$og_dir"
  if [ -d "${work}/og.before" ]; then cp -pR "${work}/og.before" "$og_dir"; fi
}

for side in a b; do
  if [ "$side" = a ]; then cmd=("${side_a[@]}"); else cmd=("${side_b[@]}"); fi
  if [ "$side" = b ]; then
    undo_images
    og_restore
  fi
  images >"${work}/${side}.before.images"
  mkdir "${work}/install" "${work}/${side}.nxf"
  seed "${work}/install"
  echo "== side ${side}: ${cmd[*]} ${work}/install"
  start=$SECONDS
  rc=0
  NXF_HOME="${work}/${side}.nxf" "${cmd[@]}" "${work}/install" >"${work}/${side}.log" 2>&1 || rc=$?
  echo "$rc" >"${work}/${side}.exit"
  echo "   exit ${rc} after $((SECONDS - start)) s"
  mv "${work}/install" "${work}/${side}"
  # A file the side left that is not ours (a container run without
  # --user): the next setup run as this user cannot rewrite or delete it.
  find "${work}/${side}" "${work}/${side}.nxf" ! -user "$uid" >"${work}/${side}.foreign" 2>&1 || true
  take "$side" "${work}/${side}"
  take "${side}.nxf" "${work}/${side}.nxf"
  images >"${work}/${side}.images"
  delta "$side"
  if [ -d "$og_dir" ]; then take "${side}.og" "$og_dir"; else echo none >"${work}/${side}.og.tsv"; fi
done

# Tree A must hold something the seed did not: a file or link that is not
# in the seed as it is, or a new folder. Otherwise the run proves nothing.
new=$(awk -F'\t' '
  FILENAME == ARGV[1] { if (FNR > 2) { row[$0] = 1; if ($2 == "d") dir[$1] = 1 } next }
  FNR > 2 && !($0 in row) && !($2 == "d" && ($1 in dir)) { n++ }
  END { print n + 0 }' "$seed_tsv" "${work}/a.tsv")
if [ "$new" -eq 0 ]; then
  echo "ERROR: tree A holds nothing the seed did not (side A exit $(cat "${work}/a.exit")): the compare would prove nothing; read ${work}/a.log." >&2
  exit 2
fi
echo "   tree A: ${new} entries the seed did not hold"

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
    awk '!/^#/ && n++ < 20 {print "         " $0}' "${work}/check.${name}.txt"
  fi
}
# shellcheck disable=SC2329  # called through check
same_file() { cmp "$1" "$2" || { diff "$1" "$2" | head -n 20; return 1; }; }
# shellcheck disable=SC2329  # called through check
side_a_exit() {
  local got
  got=$(cat "${work}/a.exit")
  [ "$got" = "$expect" ] || { echo "side A exited ${got}, not ${expect} (AB_EXPECT_EXIT); read ${work}/a.log"; return 1; }
}
# shellcheck disable=SC2329  # called through check
no_foreign() {
  local side rc=0
  for side in a b; do
    if [ -s "${work}/${side}.foreign" ]; then
      echo "side ${side} left entries uid ${uid} does not own:"
      head -n 20 "${work}/${side}.foreign"
      rc=1
    fi
  done
  return "$rc"
}

check side-a-exit side_a_exit
check exit-codes same_file "${work}/a.exit" "${work}/b.exit"
check tree manifest --compare "${work}/a.tsv" "${work}/b.tsv" --allow "$allow"
check nxf-home manifest --compare "${work}/a.nxf.tsv" "${work}/b.nxf.tsv" --allow "$allow"
check logs manifest --compare-logs "${work}/a.log" "${work}/b.log"
check images same_file "${work}/a.images.delta" "${work}/b.images.delta"
check bin-og same_file "${work}/a.og.tsv" "${work}/b.og.tsv"
check foreign-files no_foreign
if [ -n "${AB_SEED:-}" ]; then
  take seed.after "$AB_SEED"
  check seed-unchanged manifest --compare "${work}/seed.before.tsv" "${work}/seed.after.tsv"
fi
grep '^ALLOWED' "${work}/check.tree.txt" "${work}/check.nxf-home.txt" 2>/dev/null | sed 's/^/         /' || true

# The control: one flipped bit in a file the component wrote in tree B must
# turn the tree compare red.
control_args=(--control "${work}/b.tsv" --allow "$allow")
# With no seed there is nothing to leave out.
if [ -n "${AB_SEED:-}" ]; then control_args+=(--not-in "$seed_tsv"); fi
path=$(manifest "${control_args[@]}") || {
  echo "ERROR: tree B has no file the component wrote that is compared byte for byte, to plant the control in." >&2
  exit 2
}
# The file may still be a hard link into the seed: give B its own copy first.
# printf -v keeps a trailing newline in the name, which $( ) would strip.
printf -v rel '%b' "$path"
real="${work}/b/${rel}"
copy=$(mktemp "${real}.XXXXXX")
if ! cp -p "$real" "$copy" || ! mv -f "$copy" "$real"; then
  rm -f "$copy"
  echo "ERROR: could not give tree B its own copy of ${path}; a flip could reach the seed through a hard link." >&2
  exit 2
fi
# A read-only file gets owner write for each flip only, so tree B keeps its
# mode and the compare sees one difference, the content.
flip() {
  if [ -w "$real" ]; then
    manifest --flip "${work}/b" "$path" >/dev/null
  else
    chmod u+w "$real"
    manifest --flip "${work}/b" "$path" >/dev/null
    chmod u-w "$real"
  fi
}
flip
take b.control "${work}/b"
crc=0
manifest --compare "${work}/a.tsv" "${work}/b.control.tsv" --allow "$allow" >"${work}/check.control.txt" || crc=$?
flip
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
