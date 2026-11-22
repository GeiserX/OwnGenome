#!/usr/bin/env bash
# The STRchive catalog of steps 09 and 9b: setup.sh installs it, step 09 runs
# on it by default, step 9b grades with the catalog step 09 ran on, and
# EH_CATALOG=bundled keeps the call of before. The fake download serves the
# two real STRchive files of versions.env, which this case fetches itself.
# setup.sh converts them with og str catalog, which a hook runs for real on
# the host (bin/og; the Guard fake-docker job builds it into bin/.og).
# setup.sh's og check wants both Linux builds, so the case puts a stand-in
# that carries this checkout's rust/ tree hash (the check reads only that)
# for each one bin/.og lacks, and takes them out after. The hook also plays
# Stranger.
#   1. A copy of the ExpansionHunter catalog with one byte flipped: setup.sh
#      --strchive-data exits 1, keeps neither that file nor a catalog, and
#      never runs the converter (the control that the sha256 check fails).
#   2. The real files: setup.sh --strchive-data writes catalog.json with 82
#      loci through og str catalog in OG_IMAGE, run from the catalog's
#      folder with relative paths (so OG_SHADOW can check it), and a second run
#      downloads nothing and converts nothing.
#   2b. A catalog whose stamp no longer matches, with a refetch that fails
#      its sha256 check: setup.sh exits 1 and removes the old catalog.json
#      and str_ranges.tsv, so step 09 does not run on a catalog setup.sh
#      says it could not install.
#   3. Step 09 with no EH_CATALOG runs ExpansionHunter on catalog.json, names
#      the loci it does not grade (MUC1) and copies the catalog beside the
#      VCF; step 9b grades with that copy and counts the motif of the
#      flagged FAME locus with it (og str motif, from expansion_hunter/ with
#      a relative output and OG_SHADOW passed on, never bin/rfc1_motif.py);
#      STRANGER_REPEATS=bundled goes back to
#      Stranger's GRCh38 catalog.
#   4. EH_CATALOG=bundled: the ExpansionHunter call is the one step 09 makes
#      without the STRchive file (the call of before this change), byte for
#      byte, and leaves no catalog copy, so step 9b grades with Stranger's
#      bundled GRCh38 catalog.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
use_output_hook   # writes host-path.sh, which the hook below reads
DIR="${GENOME_DIR}/reference/strchive-${STRCHIVE_VERSION}"
EH_FILE="STRchive-disease-loci-${STRCHIVE_VERSION}.hg38.expansionhunter.json"
LOCI_FILE="STRchive-loci-${STRCHIVE_VERSION}.json"
SERVED="${CASE_WORK}/served" BAD="${CASE_WORK}/served-bad"
mkdir -p "$SERVED" "$BAD"
# The real curl, not the fake one first on PATH.
REAL_CURL=$(PATH=${PATH#"${REPO_ROOT}/scripts/ci/fake-docker:"} command -v curl) || fail "no curl on the system PATH"
for name in "$EH_FILE" "$LOCI_FILE"; do
  "$REAL_CURL" -fsSL --retry 2 --retry-delay 5 --max-time 60 -o "${SERVED}/${name}" \
    "https://github.com/dashnowlab/STRchive/releases/download/${STRCHIVE_VERSION}/${name}" \
    || echo "download of ${name} failed"
done
sha() { if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | awk '{print $1}'; else shasum -a 256 "$1" | awk '{print $1}'; fi; }
[ "$(sha "${SERVED}/${EH_FILE}")" = "$STRCHIVE_EH_SHA256" ] || fail "the served ${EH_FILE} is not the pinned file"
[ "$(sha "${SERVED}/${LOCI_FILE}")" = "$STRCHIVE_LOCI_SHA256" ] || fail "the served ${LOCI_FILE} is not the pinned file"
cp "${SERVED}/${LOCI_FILE}" "${BAD}/"
# The same file with the lowest bit of its middle byte flipped.
flip() {
  local n h b
  n=$(wc -c < "$1" | tr -d ' ')
  h=$((n / 2))
  b=$(od -An -tu1 -j "$h" -N1 "$1" | tr -d ' ')
  { head -c "$h" "$1"; printf "\\$(printf %03o $((b ^ 1)))"; tail -c $((n - h - 1)) "$1"; } > "$2"
}
flip "${SERVED}/${EH_FILE}" "${BAD}/${EH_FILE}"
[ "$(wc -c < "${BAD}/${EH_FILE}")" = "$(wc -c < "${SERVED}/${EH_FILE}")" ] && ! cmp -s "${SERVED}/${EH_FILE}" "${BAD}/${EH_FILE}" \
  || fail "could not flip a byte of ${EH_FILE}"

# The real og runs the conversion; a stand-in for each Linux build bin/.og
# lacks satisfies setup.sh's og check, and goes when the case ends.
"${REPO_ROOT}/bin/og" version > /dev/null || fail "og does not run here: build it into bin/.og (rust/og.sh, the Guard fake-docker job does)"
OGDIR="${REPO_ROOT}/bin/.og"
TREE=$("${REPO_ROOT}/rust/og.sh" tree-hash) || fail "could not compute the rust/ tree hash"
STANDINS=()
trap 'for f in "${STANDINS[@]}"; do rm -f "$f"; done' EXIT
for a in x86_64 aarch64; do
  [ ! -e "${OGDIR}/og-linux-${a}" ] || continue
  printf 'stand-in og-tree-hash:%s;\n' "$TREE" > "${OGDIR}/og-linux-${a}"; chmod 755 "${OGDIR}/og-linux-${a}"
  STANDINS+=("${OGDIR}/og-linux-${a}")
done
"${REPO_ROOT}/rust/og.sh" check --required > /dev/null || fail "bin/.og does not pass setup.sh's og check"
export STR_DIR="$DIR"

cat > "${CASE_WORK}/tools-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
shift   # the image
case " $* " in
  *" /pgp-bin/og str catalog "*)
    # og, in the catalog's folder (setup.sh's -w), with the container paths
    # mapped to the host.
    shift 3   # /pgp-bin/og str catalog
    args=()
    for a in "$@"; do
      case "$a" in
        /*) args+=("$(host_path "$a")") ;;
        *) args+=("$a") ;;
      esac
    done
    cd "$STR_DIR" && exec "${REPO_ROOT}/bin/og" str catalog "${args[@]}" ;;
  *" /pgp-bin/og str motif "*" --flagged "*)
    # --out is relative to the container's -w, which only this call's log
    # line (the last one) shows: the hook gets no docker options.
    wd=$(tail -n 1 "$FAKE_DOCKER_LOG" | sed -n 's/.* -w \([^ ]*\) .*/\1/p')
    prev=""
    for a in "$@"; do [ "$prev" = --out ] && out=$(host_path "${wd}/${a}"); prev=$a; done
    printf 'locus\tverdict\nSAMD12\tpathogenic_motif\n' > "$out" ;;
  *" stranger "*)
    # One FAME record Stranger flags, so step 9b counts its motif.
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n'
    printf 'chr8\t118366913\tSAMD12\tA\t<STR40>\t.\tPASS\tEND=118366918;VARID=SAMD12;REPID=SAMD12;STR_STATUS=pre_mutation\tGT\t0/1\n' ;;
  *) exec "${CASE_WORK}/hook-outputs" image "$@" ;;
esac
HOOK
chmod +x "${CASE_WORK}/tools-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/tools-hook"
fetches() { grep -cE "^(curl|wget) [^ ]*/STRchive-" "$FAKE_DOCKER_LOG" || true; }
converts() { grep -c '/pgp-bin/og str catalog' "$FAKE_DOCKER_LOG" || true; }

# --- 1. a flipped byte: no file, no catalog, no conversion --------------------------------
: > "$FAKE_DOCKER_LOG"
FAKE_DOWNLOAD_DIR="$BAD" run_expect 1 setup-bad "${SCRIPTS}/setup.sh" --strchive-data "$GENOME_DIR"
output_has setup-bad "sha256 checksum of ${EH_FILE//./\\.} is "
output_has setup-bad "\[WARN\] Could not install STRchive's ${EH_FILE//./\\.}"
[ ! -e "${DIR}/${EH_FILE}" ] || fail "setup.sh kept the ExpansionHunter catalog that failed its sha256 check"
[ ! -e "${DIR}/catalog.json" ] || fail "setup.sh wrote a catalog from a file that failed its sha256 check"
[ "$(converts)" -eq 0 ] || fail "setup.sh ran the converter on a file that failed its sha256 check"

# --- 2. the real files: converted once ------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
FAKE_DOWNLOAD_DIR="$SERVED" run_expect 0 setup "${SCRIPTS}/setup.sh" --strchive-data "$GENOME_DIR"
output_has setup "\[OK\] STRchive ${STRCHIVE_VERSION//./\\.} catalog \(steps 09 and 9b\): "
[ "$(fetches)" -eq 2 ] || fail "setup.sh downloaded the STRchive files $(fetches) times, expected 2"
docker_log_has "image=${OG_IMAGE//./\\.} .* -w /genome/reference/strchive-${STRCHIVE_VERSION//./\\.} .*/pgp-bin/og str catalog --eh ${EH_FILE//./\\.} " \
  "setup.sh did not run og str catalog in OG_IMAGE from the catalog's folder"
N=$(jq length "${DIR}/catalog.json" 2>/dev/null || echo 0)
[ "$N" -eq 82 ] || fail "catalog.json has ${N} loci, expected 82"
: > "$FAKE_DOCKER_LOG"
FAKE_DOWNLOAD_DIR="$SERVED" run_expect 0 setup-again "${SCRIPTS}/setup.sh" --strchive-data "$GENOME_DIR"
output_has setup-again "already present"
[ "$(fetches)" -eq 0 ] || fail "the second setup.sh downloaded the STRchive files again"
[ "$(converts)" -eq 0 ] || fail "the second setup.sh converted again"

# --- 2b. a stale catalog whose refetch fails: setup.sh removes it --------------------------
echo stale > "${DIR}/catalog.stamp"
rm -f "${DIR}/${EH_FILE}"
: > "$FAKE_DOCKER_LOG"
FAKE_DOWNLOAD_DIR="$BAD" run_expect 1 setup-stale "${SCRIPTS}/setup.sh" --strchive-data "$GENOME_DIR"
[ ! -e "${DIR}/catalog.json" ] || fail "setup.sh kept the stale catalog.json after its refetch failed, so step 09 would run on it"
[ ! -e "${DIR}/str_ranges.tsv" ] || fail "setup.sh kept the stale str_ranges.tsv after its refetch failed"
FAKE_DOWNLOAD_DIR="$SERVED" run_expect 0 setup-restore "${SCRIPTS}/setup.sh" --strchive-data "$GENOME_DIR"

# --- 3. steps 09 and 9b on the STRchive catalog -------------------------------------------------
EH_DIR="${GENOME_DIR}/sample1/expansion_hunter"
eh_call() { grep -e '--variant-catalog' "$FAKE_DOCKER_LOG" | tail -n 1; }  # the logged docker run of ExpansionHunter
: > "$FAKE_DOCKER_LOG"
run_expect 0 eh-strchive "${SCRIPTS}/09-expansion-hunter.sh" sample1 male
output_has eh-strchive "Variant catalog: STRchive ${STRCHIVE_VERSION//./\\.}, "
output_has eh-strchive "Not graded \(no call from this count\): MUC1: the reference repeat is 3734 bases"
docker_log_has "--variant-catalog /genome/reference/strchive-${STRCHIVE_VERSION//./\\.}/catalog\\.json " \
  "step 09 did not run ExpansionHunter on the STRchive catalog"
cmp -s "${DIR}/catalog.json" "${EH_DIR}/sample1_eh_catalog.json" || fail "step 09 did not copy the catalog it ran on beside the VCF"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n' > "${EH_DIR}/sample1_eh.vcf"
: > "$FAKE_DOCKER_LOG"
run_expect 0 stranger-strchive "${SCRIPTS}/09b-stranger.sh" sample1
output_has stranger-strchive "Repeat catalog: the one step 09 ran on, .*sample1_eh_catalog\.json \(STRchive-${STRCHIVE_VERSION//./\\.}\)"
docker_log_has "stranger --repeats-file /genome/sample1/expansion_hunter/sample1_eh_catalog\\.json " \
  "step 9b did not grade with the catalog step 09 ran on"
docker_log_has "-e OG_SHADOW -w /genome/sample1/expansion_hunter .*/pgp-bin/og str motif .*--catalog /genome/sample1/expansion_hunter/sample1_eh_catalog\\.json --flagged --out sample1_str_motifs\\.tsv\\.part " \
  "step 9b did not count the motif of the flagged FAME locus with og, from expansion_hunter/, with the catalog"
if grep -q 'rfc1_motif\.py' "$FAKE_DOCKER_LOG"; then fail "step 9b still calls bin/rfc1_motif.py"; fi
grep -q $'^SAMD12\tpathogenic_motif' "${EH_DIR}/sample1_str_motifs.tsv" 2>/dev/null || fail "step 9b wrote no sample1_str_motifs.tsv"
rm -f "${EH_DIR}/sample1_eh_stranger.vcf"
: > "$FAKE_DOCKER_LOG"
STRANGER_REPEATS=bundled run_expect 0 stranger-bundled "${SCRIPTS}/09b-stranger.sh" sample1
docker_log_has 'variant_catalog_grch38\.json' "STRANGER_REPEATS=bundled did not pass Stranger's GRCh38 catalog"
[ ! -e "${EH_DIR}/sample1_str_motifs.tsv" ] || fail "with Stranger's bundled catalog step 9b kept the STRchive motif table"

# --- 4. EH_CATALOG=bundled: the call of before ---------------------------------------------------
mv "$DIR" "${DIR}.away"
: > "$FAKE_DOCKER_LOG"
run_expect 0 eh-before "${SCRIPTS}/09-expansion-hunter.sh" sample1 male
output_has eh-before "Variant catalog: the 31 loci inside the image \(reference/strchive-"
BEFORE=$(eh_call)
mv "${DIR}.away" "$DIR"
: > "$FAKE_DOCKER_LOG"
EH_CATALOG=bundled run_expect 0 eh-bundled "${SCRIPTS}/09-expansion-hunter.sh" sample1 male
output_has eh-bundled "Variant catalog: the 31 loci inside the image \(EH_CATALOG=bundled\)"
AFTER=$(eh_call)
echo "without STRchive:    ${BEFORE}"
echo "EH_CATALOG=bundled:  ${AFTER}"
[ -n "$BEFORE" ] && [ "$BEFORE" = "$AFTER" ] || fail "EH_CATALOG=bundled does not make the ExpansionHunter call of before"
case "$AFTER" in
  *"--variant-catalog /usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json "*) ;;
  *) fail "EH_CATALOG=bundled does not pass the catalog inside the image: ${AFTER}" ;;
esac
[ ! -e "${EH_DIR}/sample1_eh_catalog.json" ] || fail "EH_CATALOG=bundled left the STRchive catalog copy beside the VCF"
rm -f "${EH_DIR}/sample1_eh_stranger.vcf"
: > "$FAKE_DOCKER_LOG"
run_expect 0 stranger-after-bundled "${SCRIPTS}/09b-stranger.sh" sample1
docker_log_has 'variant_catalog_grch38\.json' "after EH_CATALOG=bundled step 9b did not grade with Stranger's GRCh38 catalog"
