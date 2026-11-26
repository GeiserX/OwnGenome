#!/usr/bin/env bash
# Raw downloads an older setup left at its own paths (vep_cache/, pcgr_data/)
# move to the downloads folder (adopt in scripts/lib/common.sh), where fetch
# checks them by their checksum:
#   - whole archives are used with no download, a half one resumes, and the
#     installed caches stay where they are; a rerun makes no network call;
#   - a wrong old archive is refused and left in the downloads folder, and
#     nothing is installed from it;
#   - an old copy beside one already in the downloads folder: the copy that
#     passes its check stays, the other goes; a download is never deleted
#     before a check says which is the right one.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
seed_clinical_downloads
export FETCH_TRIES=1 FETCH_WAIT=0
V116=homo_sapiens_vep_116_GRCh38.tar.gz
V115=homo_sapiens_vep_115_GRCh38.tar.gz
PCGR=pcgr_ref_data.20260620.grch38.tgz
D=$FAKE_DOWNLOAD_DIR
downloaded() { grep -E "^(wget|curl) .*${1//./\\.}( |$)" "$FAKE_DOCKER_LOG" | grep -v CHECKSUMS || true; }
installed() {
  for rel in 116 115; do
    [ -s "${1}/vep_cache/homo_sapiens/${rel}_GRCh38/info.txt" ] || fail "VEP ${rel} is not installed in ${1}"
  done
  [ -s "${1}/pcgr_data/20260620/data/grch38/README" ] || fail "the PCGR bundle is not installed in ${1}"
}

# --- the old paths: two whole archives and half of a third ---------------------
G="${CASE_WORK}/old"
mkdir -p "${G}/vep_cache" "${G}/pcgr_data"
cp "${D}/${V116}" "${G}/vep_cache/"
cp "${D}/${PCGR}" "${G}/pcgr_data/"
head -c "$(( $(wc -c < "${D}/${V115}") / 2 ))" "${D}/${V115}" > "${G}/vep_cache/${V115}.part"
: > "$FAKE_DOCKER_LOG"
run_expect 0 old "${SCRIPTS}/setup.sh" --clinical-data "$G"
installed "$G"
[ -z "$(downloaded "$V116")" ] || fail "the old VEP 116 archive was downloaded again: $(downloaded "$V116")"
[ -z "$(downloaded "$PCGR")" ] || fail "the old PCGR bundle was downloaded again: $(downloaded "$PCGR")"
[ -n "$(downloaded "$V115")" ] || fail "the half VEP 115 archive was not resumed"
for f in "$V116" "$V115" "$PCGR"; do
  cmp -s "${G}/downloads/${f}" "${D}/${f}" || fail "${f} is not the whole download in the downloads folder"
done
for f in "vep_cache/${V116}" "vep_cache/${V115}.part" "pcgr_data/${PCGR}"; do
  [ ! -e "${G}/${f}" ] || fail "${f} stayed at the old path"
done
output_has old "Moved ${G//./\\.}/vep_cache/${V116//./\\.}, left by an older setup"
cp "$FAKE_DOCKER_LOG" "${CASE_WORK}/old.log"
run_expect 0 again "${SCRIPTS}/setup.sh" --clinical-data "$G"
cmp -s "$FAKE_DOCKER_LOG" "${CASE_WORK}/old.log" || fail "a rerun made a network call: $(diff "${CASE_WORK}/old.log" "$FAKE_DOCKER_LOG")"

# --- a wrong old archive -----------------------------------------------------
G="${CASE_WORK}/wrong"
mkdir -p "${G}/vep_cache"
printf 'not the VEP cache\n' > "${G}/vep_cache/${V116}"
cp "${G}/vep_cache/${V116}" "${CASE_WORK}/wrong.tgz"
: > "$FAKE_DOCKER_LOG"
run_expect 1 wrong "${SCRIPTS}/setup.sh" --clinical-data "$G"
output_has wrong "is there but is not the file of vep_cache"
cmp -s "${G}/downloads/${V116}" "${CASE_WORK}/wrong.tgz" || fail "the refused old archive was not left as it is in the downloads folder"
[ ! -e "${G}/vep_cache/homo_sapiens/116_GRCh38" ] || fail "VEP 116 was installed from a wrong archive"
[ -z "$(downloaded "$V116")" ] || fail "the refused archive was downloaded over: $(downloaded "$V116")"

# --- an old copy beside the kept one -----------------------------------------
# Both right: the old one goes.
G="${CASE_WORK}/both"
mkdir -p "${G}/vep_cache" "${G}/downloads"
cp "${D}/${V116}" "${G}/vep_cache/"
cp "${D}/${V116}" "${G}/downloads/"
: > "$FAKE_DOCKER_LOG"
run_expect 0 both "${SCRIPTS}/setup.sh" --clinical-data "$G"
installed "$G"
[ ! -e "${G}/vep_cache/${V116}" ] || fail "the second copy of a checked archive stayed"
output_has both "Deleted ${G//./\\.}/vep_cache/${V116//./\\.}, left by an older setup: .* is the same download, checked"
[ -z "$(downloaded "$V116")" ] || fail "VEP 116 was downloaded although two copies were there"
# The kept one wrong, the old one right: the old one replaces it.
G="${CASE_WORK}/kept-wrong"
mkdir -p "${G}/vep_cache" "${G}/downloads"
cp "${D}/${V116}" "${G}/vep_cache/"
printf 'a broken copy\n' > "${G}/downloads/${V116}"
: > "$FAKE_DOCKER_LOG"
run_expect 0 kept_wrong "${SCRIPTS}/setup.sh" --clinical-data "$G"
installed "$G"
cmp -s "${G}/downloads/${V116}" "${D}/${V116}" || fail "the right old copy did not replace the wrong kept one"
for f in "vep_cache/${V116}" "downloads/${V116}.adopt"; do
  [ ! -e "${G}/${f}" ] || fail "${f} was left behind"
done
output_has kept_wrong "failed its check; replaced it with ${G//./\\.}/vep_cache/${V116//./\\.}"
[ -z "$(downloaded "$V116")" ] || fail "VEP 116 was downloaded although a right copy was there"
# Both wrong: neither is deleted, and fetch refuses the kept one.
G="${CASE_WORK}/both-wrong"
mkdir -p "${G}/vep_cache" "${G}/downloads"
printf 'old and broken\n' > "${G}/vep_cache/${V116}"
printf 'kept and broken\n' > "${G}/downloads/${V116}"
run_expect 1 both_wrong "${SCRIPTS}/setup.sh" --clinical-data "$G"
grep -qx 'old and broken' "${G}/vep_cache/${V116}" || fail "the old copy was deleted or changed with no check passed"
grep -qx 'kept and broken' "${G}/downloads/${V116}" || fail "the kept copy was deleted or changed with no check passed"
output_has both_wrong "neither is deleted"
output_has both_wrong "is there but is not the file of vep_cache"
echo "PASS old downloads move to the downloads folder, are checked there, and are never deleted unchecked"
