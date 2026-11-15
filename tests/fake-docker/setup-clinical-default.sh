#!/usr/bin/env bash
# Default clinical downloads, explicit opt-out, disk gates and corrupt input.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
seed_clinical_downloads
export FETCH_TRIES=1 FETCH_WAIT=0
G="${CASE_WORK}/without"
run_expect 0 without "${SCRIPTS}/setup.sh" --no-clinical-data "$G"
output_has without 'Setup complete!'
if grep -E '^(wget|curl).*homo_sapiens_vep_|^(wget|curl).*pcgr_ref_data' "$FAKE_DOCKER_LOG"; then
  fail '--no-clinical-data downloaded clinical data'
fi

G="${CASE_WORK}/clinical"
run_expect 0 clinical "${SCRIPTS}/setup.sh" --clinical-data "$G"
for rel in 116 115; do
  [ -s "${G}/vep_cache/homo_sapiens/${rel}_GRCh38/info.txt" ] || fail "VEP ${rel} missing"
done
[ -s "${G}/pcgr_data/20260620/data/grch38/README" ] || fail 'PCGR bundle missing'
output_lacks clinical 'Phase 1:|Docker Images'
# Installed clinical data turns the three default launch-plan rows on.
seed_reference "$G"
seed_clinvar "$G"
seed_sample "$G" sample1
use_output_hook
GENOME_DIR="$G" SKIP_VALIDATION=true run_expect 0 plan "${SCRIPTS}/run-all.sh" sample1 male
output_has plan '^  13 VEP +runs$'
output_has plan '^  23 Clinical filter +runs$'
output_has plan '^  17 CPSR +runs$'

# A fake df that reports 40 GiB free.
mkdir -p "${CASE_WORK}/low-bin"
cat > "${CASE_WORK}/low-bin/df" <<'DF'
#!/usr/bin/env bash
echo 'Filesystem 1024-blocks Used Available Capacity Mounted on'
echo 'fakefs 41943040 0 41943040 0% /'
DF
chmod +x "${CASE_WORK}/low-bin/df"

# A second install does no network work, even with too little space.
cp "$FAKE_DOCKER_LOG" "${CASE_WORK}/installed.log"
PATH="${CASE_WORK}/low-bin:${PATH}" run_expect 0 again "${SCRIPTS}/setup.sh" --clinical-data "$G"
cmp "$FAKE_DOCKER_LOG" "${CASE_WORK}/installed.log" || fail 'installed clinical data downloaded again'

# A skipped item makes --clinical-data exit 1; the other items still install.
G="${CASE_WORK}/low"
PATH="${CASE_WORK}/low-bin:${PATH}" run_expect 1 low "${SCRIPTS}/setup.sh" --clinical-data "$G"
output_has low '56 GB.*5 GB margin'
output_has low '52 GB.*5 GB margin'
[ ! -e "${G}/vep_cache/homo_sapiens/116_GRCh38" ] || fail 'low-space VEP installed'
[ -s "${G}/pcgr_data/20260620/data/grch38/README" ] || fail 'small bundle skipped with 40 GB free'

# A partial download counts as space already taken: 40 GiB free plus a
# 22 GiB .part covers VEP 116's 56 GB peak and margin, so it resumes.
G="${CASE_WORK}/resume"
mkdir -p "${G}/vep_cache"
dd if=/dev/null of="${G}/vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz.part" bs=1048576 seek=22528 2>/dev/null
PATH="${CASE_WORK}/low-bin:${PATH}" run_expect 1 resume "${SCRIPTS}/setup.sh" --clinical-data "$G"
output_lacks resume 'Skipping VEP 116'
[ -s "${G}/vep_cache/homo_sapiens/116_GRCh38/info.txt" ] || fail 'partial VEP 116 download not resumed with low space'
[ ! -e "${G}/vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz.part" ] || fail 'resumed VEP 116 left its .part'

# An install stopped during extraction keeps its checked archive: the next
# run unpacks it without downloading again and removes the old staging tree.
G="${CASE_WORK}/extract"
mkdir -p "${G}/vep_cache/.extract.old" "${G}/pcgr_data/.extract.old"
cp "${FAKE_DOWNLOAD_DIR}/homo_sapiens_vep_116_GRCh38.tar.gz" "${G}/vep_cache/"
cp "${FAKE_DOWNLOAD_DIR}/pcgr_ref_data.20260620.grch38.tgz" "${G}/pcgr_data/"
FAKE_DOWNLOAD_FAIL='vep_116|pcgr_ref_data' run_expect 0 extract "${SCRIPTS}/setup.sh" --clinical-data "$G"
[ -s "${G}/vep_cache/homo_sapiens/116_GRCh38/info.txt" ] || fail 'checked VEP 116 archive downloaded again'
[ -s "${G}/pcgr_data/20260620/data/grch38/README" ] || fail 'checked PCGR archive downloaded again'
[ ! -e "${G}/vep_cache/.extract.old" ] || fail 'stale VEP staging tree kept'
[ ! -e "${G}/pcgr_data/.extract.old" ] || fail 'stale PCGR staging tree kept'

# The digest belongs to the original fixture. Serve a different valid tarball.
printf 'corrupt\n' >> "${FAKE_DOWNLOAD_DIR}/pcgr_ref_data.20260620.grch38.tgz"
G="${CASE_WORK}/corrupt"
run_expect 1 corrupt "${SCRIPTS}/setup.sh" --clinical-data "$G"
output_has corrupt 'sha256 checksum.*expected'
output_has corrupt '\[WARN\].*PCGR'
[ ! -e "${G}/pcgr_data/20260620" ] || fail 'corrupt bundle published'
[ ! -e "${G}/pcgr_data/pcgr_ref_data.20260620.grch38.tgz.part" ] || fail 'bad bundle kept'
[ -s "${G}/vep_cache/homo_sapiens/115_GRCh38/info.txt" ] || fail 'bundle failure stopped remaining install'

