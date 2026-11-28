#!/usr/bin/env bash
# setup.sh --vep-cache and --clinical-data run by og
# (rust/crates/og-setup/src/install/clinical.rs) against setup.sh's bash
# (OG_SETUP=bash), in each state an install meets: a fresh install, a second
# run, too little disk, a partial download that the disk check counts, a
# stopped extraction, a kept archive that fails its check, a corrupt
# download, no download at all. Each pair gives the same exit code, the same
# messages, the same download calls and the same files with the same modes.
# Step 13 installs its cache at run time through og (common.sh's
# install_vep_cache), and a full setup with too little disk prints the same
# skip lines and the same missing list on both sides.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"
# Each run below says which side it is on (OG_SETUP), whatever the suite was given.
unset OG_SETUP
seed_clinical_downloads
export FETCH_TRIES=1 FETCH_WAIT=0

# An og that notes each `setup install` it is handed, then runs the real one.
mkdir -p "${CASE_WORK}/og-bin"
cat > "${CASE_WORK}/og-bin/og" <<'SHIM'
#!/usr/bin/env bash
if [ "${1:-} ${2:-}" = "setup install" ]; then printf '%s\n' "$*" >> "${CASE_WORK}/og.calls"; fi
exec "${REPO_ROOT}/bin/og" "$@"
SHIM
chmod +x "${CASE_WORK}/og-bin/og"
export OG="${CASE_WORK}/og-bin/og"
: > "${CASE_WORK}/og.calls"

# A df that reports 40 GiB free; a du that reports 22 GiB for a .part under resume.
mkdir -p "${CASE_WORK}/low-bin"
cat > "${CASE_WORK}/low-bin/df" <<'DF'
#!/usr/bin/env bash
echo 'Filesystem 1024-blocks Used Available Capacity Mounted on'
echo 'fakefs 41943040 0 41943040 0% /'
DF
cat > "${CASE_WORK}/low-bin/du" <<'DU'
#!/usr/bin/env bash
case "$2" in
  */resume-*/genome/downloads/*.part) printf '23068672\t%s\n' "$2" ;;
  *) /usr/bin/du "$@" ;;
esac
DU
chmod +x "${CASE_WORK}/low-bin/df" "${CASE_WORK}/low-bin/du"

# SIDE_TREE DIR: every folder, file and link under DIR with its mode, and
# each file's checksum.
side_tree() {
  (cd "$1" && find . -mindepth 1 -print | LC_ALL=C sort | while IFS= read -r f; do
    m=$(stat -c %a "$f" 2>/dev/null || stat -f %Lp "$f")
    case "$f" in */.extract.*) f="${f%%/.extract.*}/.extract.X" ;; esac
    if [ -L "$f" ]; then printf '%s -> %s\n' "$f" "$(readlink "$f")"
    elif [ -d "$f" ]; then printf '%s/ %s\n' "$f" "$m"
    else printf '%s %s %s\n' "$f" "$m" "$(cksum < "$f")"; fi
  done)
}

# compare NAME PREP ARGS...: setup.sh ARGS under OG_SETUP=bash and under og,
# each in its own install folder (named G in what is compared), after PREP G
# (a function that plants a state, or :). LOW=1 puts the low df and du first
# on PATH; FAIL is FAKE_DOWNLOAD_FAIL.
compare() {
  local name=$1 prep=$2 side g start calls want_calls path
  shift 2
  for side in bash og; do
    g="${CASE_WORK}/${name}-${side}/genome"
    mkdir -p "$g"
    "$prep" "$g"
    start=$(wc -l < "$FAKE_DOCKER_LOG")
    calls=$(grep -c '^setup install --clinical-data\|^setup install --vep-cache' "${CASE_WORK}/og.calls" || true)
    path=$PATH
    [ "${LOW:-}" != 1 ] || path="${CASE_WORK}/low-bin:${PATH}"
    PATH=$path OG_SETUP=$side FAKE_DOWNLOAD_FAIL=${FAIL:-} run_rc "${name}-${side}" "${SCRIPTS}/setup.sh" "${@//@G@/$g}"
    echo "$RC" > "${CASE_WORK}/${name}-${side}.rc"
    sed "s|${g}|G|g" "${CASE_WORK}/${name}-${side}.out" > "${CASE_WORK}/${name}-${side}.said"
    tail -n "+$((start + 1))" "$FAKE_DOCKER_LOG" | sed "s|${g}|G|g" > "${CASE_WORK}/${name}-${side}.calls"
    side_tree "$g" > "${CASE_WORK}/${name}-${side}.tree"
    want_calls=$calls
    [ "$side" = og ] && want_calls=$((calls + 1))
    [ "$(grep -c '^setup install --clinical-data\|^setup install --vep-cache' "${CASE_WORK}/og.calls" || true)" -eq "$want_calls" ] \
      || fail "${name}: OG_SETUP=${side} setup.sh $* was$([ "$side" = og ] || echo ' not') to be handed to og"
  done
  for what in rc said calls tree; do
    diff -u "${CASE_WORK}/${name}-bash.${what}" "${CASE_WORK}/${name}-og.${what}" \
      || fail "${name}: og and setup.sh's bash differ in ${what} (the diff above, bash first)"
  done
  echo "[OK] ${name}: og and setup.sh's bash give exit $(cat "${CASE_WORK}/${name}-og.rc"), the same $(wc -l < "${CASE_WORK}/${name}-og.said" | tr -d ' ') lines, $(wc -l < "${CASE_WORK}/${name}-og.calls" | tr -d ' ') calls and $(wc -l < "${CASE_WORK}/${name}-og.tree" | tr -d ' ') entries"
}

"$OG" setup ported --clinical-data > /dev/null || fail "og setup ported: --clinical-data is ported"
"$OG" setup ported --vep-cache > /dev/null || fail "og setup ported: --vep-cache is ported"

# The states, planted the same way on both sides.
installed() { # a whole install, then its second run
  PATH="${CASE_WORK}/og-bin:${PATH}" OG_SETUP=bash "${SCRIPTS}/setup.sh" --clinical-data "$1" > /dev/null 2>&1
}
partial() { # a 22 GiB partial download of VEP 116 (sparse: the du says so)
  mkdir -p "$1/downloads"
  dd if=/dev/null of="$1/downloads/homo_sapiens_vep_116_GRCh38.tar.gz.part" bs=1048576 seek=22528 2>/dev/null
}
stopped() { # an extraction stopped half way, with the checked archives kept
  mkdir -p "$1/vep_cache/.extract.old/homo_sapiens" "$1/pcgr_data/.extract.old/data" "$1/downloads"
  cp "${FAKE_DOWNLOAD_DIR}/homo_sapiens_vep_116_GRCh38.tar.gz" "${FAKE_DOWNLOAD_DIR}/pcgr_ref_data.20260620.grch38.tgz" "$1/downloads/"
}
no_info() { # a cache folder without info.txt
  mkdir -p "$1/vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38"
}
bad_kept() { # a kept PCGR archive with the right layout and the wrong bytes
  mkdir -p "$1/downloads" "${CASE_WORK}/altered/data/grch38"
  printf 'modified reference\n' > "${CASE_WORK}/altered/data/grch38/README"
  tar -czf "$1/downloads/pcgr_ref_data.20260620.grch38.tgz" -C "${CASE_WORK}/altered" data
}
old_place() { # an older setup's archive beside the cache, and a partial one
  mkdir -p "$1/vep_cache" "$1/pcgr_data"
  cp "${FAKE_DOWNLOAD_DIR}/homo_sapiens_vep_116_GRCh38.tar.gz" "$1/vep_cache/"
  head -c 100 "${FAKE_DOWNLOAD_DIR}/pcgr_ref_data.20260620.grch38.tgz" > "$1/pcgr_data/pcgr_ref_data.20260620.grch38.tgz.part"
}

for flag in --clinical-data --vep-cache; do
  n=${flag#--}
  compare "${n}" : "$flag" @G@
  compare "${n}-again" installed "$flag" @G@
  compare "${n}-no-keep" : --no-keep-downloads "$flag" @G@
  LOW=1 compare "${n}-low" : "$flag" @G@
  LOW=1 compare "resume-${n}" partial "$flag" @G@
  FAIL='vep_116|pcgr_ref_data' compare "${n}-stopped" stopped "$flag" @G@
  compare "${n}-no-info" no_info "$flag" @G@
  compare "${n}-old-place" old_place "$flag" @G@
  FAIL='vep_|pcgr_ref_data' compare "${n}-offline" : "$flag" @G@
done
compare clinical-data-bad-kept bad_kept --clinical-data @G@
grep -q '56 GB.*5 GB margin' "${CASE_WORK}/clinical-data-low-og.said" || fail "the low-disk run printed no skip line"
grep -q 'Moved G/vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz, left by an older setup' "${CASE_WORK}/clinical-data-old-place-og.said" \
  || fail "the older setup's archive was not moved to the downloads folder"
grep -q 'exists but has no info.txt' "${CASE_WORK}/vep-cache-no-info-og.said" || fail "a cache folder without info.txt did not stop"

# The corrupt download: the served bundle no longer matches its pin.
cp "${FAKE_DOWNLOAD_DIR}/pcgr_ref_data.20260620.grch38.tgz" "${CASE_WORK}/pcgr.good"
printf 'corrupt\n' >> "${FAKE_DOWNLOAD_DIR}/pcgr_ref_data.20260620.grch38.tgz"
compare clinical-data-corrupt : --clinical-data @G@
grep -q 'sha256 checksum.*expected' "${CASE_WORK}/clinical-data-corrupt-og.said" || fail "the corrupt bundle was not refused"
cp "${CASE_WORK}/pcgr.good" "${FAKE_DOWNLOAD_DIR}/pcgr_ref_data.20260620.grch38.tgz"

# --- a full setup with too little disk: the same skip lines and missing list ---
# The clinical part of a full setup is setup.sh's bash calling common.sh's
# helpers, which hand the work to og unless OG_SETUP=bash.
for side in bash og; do
  g="${CASE_WORK}/full-${side}/genome"
  mkdir -p "$g"
  n=$(grep -c -- '--part' "${CASE_WORK}/og.calls" || true)
  PATH="${CASE_WORK}/low-bin:${PATH}" OG_SETUP=$side run_rc "full-${side}" "${SCRIPTS}/setup.sh" "$g"
  echo "$RC" > "${CASE_WORK}/full-${side}.rc"
  # From the clinical header to the end: the skip lines, the bundle, the summary.
  sed -n '/^Clinical data: about 57 GB/,$p' "${CASE_WORK}/full-${side}.out" | sed "s|${g}|G|g" > "${CASE_WORK}/full-${side}.said"
  [ -s "${CASE_WORK}/full-${side}.said" ] || fail "full setup (${side}) printed no clinical part"
  side_tree "$g/vep_cache" > "${CASE_WORK}/full-${side}.tree"
  side_tree "$g/pcgr_data" >> "${CASE_WORK}/full-${side}.tree"
  m=$(grep -c -- '--part' "${CASE_WORK}/og.calls" || true)
  if [ "$side" = og ]; then
    [ "$m" -gt "$n" ] || fail "a full setup's clinical part did not call og (og setup install --part)"
  else
    [ "$m" -eq "$n" ] || fail "a full setup under OG_SETUP=bash called og for its clinical part"
  fi
done
for what in rc said tree; do
  diff -u "${CASE_WORK}/full-bash.${what}" "${CASE_WORK}/full-og.${what}" \
    || fail "full setup with too little disk: og and bash differ in ${what} (the diff above, bash first)"
done
grep -q '^  - VEP 116 cache (steps 13, 30, 23 and 31)$' "${CASE_WORK}/full-og.said" || fail "the missing list does not name VEP 116"
grep -q '^  - VEP 115 cache (step 17, CPSR)$' "${CASE_WORK}/full-og.said" || fail "the missing list does not name VEP 115"
echo "[OK] full setup with too little disk: og and bash print the same $(wc -l < "${CASE_WORK}/full-og.said" | tr -d ' ') lines from the clinical part on, missing list included"

# --- step 13 installs its cache at run time, through og ----------------------------
for side in bash og; do
  export GENOME_DIR="${CASE_WORK}/step13-${side}/genome"
  seed_reference "$GENOME_DIR"
  seed_sample "$GENOME_DIR" sample1
  use_output_hook
  n=$(grep -c -- '--part vep-cache' "${CASE_WORK}/og.calls" || true)
  OG_SETUP=$side run_expect 0 "step13-${side}" "${SCRIPTS}/13-vep-annotation.sh" sample1
  [ -f "${GENOME_DIR}/vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt" ] \
    || fail "step 13 (${side}) did not install its cache at run time"
  output_has "step13-${side}" "Cache installed at "
  m=$(grep -c -- '--part vep-cache' "${CASE_WORK}/og.calls" || true)
  if [ "$side" = og ]; then [ "$m" -eq $((n + 1)) ] || fail "step 13 did not install its cache through og"
  else [ "$m" -eq "$n" ] || fail "step 13 under OG_SETUP=bash installed its cache through og"; fi
  side_tree "${GENOME_DIR}/vep_cache" > "${CASE_WORK}/step13-${side}.tree"
  sed -n '/not found in/,/^Cache installed at/p' "${CASE_WORK}/step13-${side}.out" | sed "s|${GENOME_DIR}|G|g" > "${CASE_WORK}/step13-${side}.said"
done
unset GENOME_DIR
for what in said tree; do
  diff -u "${CASE_WORK}/step13-bash.${what}" "${CASE_WORK}/step13-og.${what}" \
    || fail "step 13's run-time install: og and bash differ in ${what}"
done
echo "[OK] step 13 installs its cache at run time through og, as the bash did"
