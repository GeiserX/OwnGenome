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

# A second install does no network work, even with too little space.
cp "$FAKE_DOCKER_LOG" "${CASE_WORK}/installed.log"
run_expect 0 again "${SCRIPTS}/setup.sh" --clinical-data "$G"
cmp "$FAKE_DOCKER_LOG" "${CASE_WORK}/installed.log" || fail 'installed clinical data downloaded again'

mkdir -p "${CASE_WORK}/low-bin"
cat > "${CASE_WORK}/low-bin/df" <<'DF'
#!/usr/bin/env bash
echo 'Filesystem 1024-blocks Used Available Capacity Mounted on'
echo 'fakefs 41943040 0 41943040 0% /'
DF
chmod +x "${CASE_WORK}/low-bin/df"
G="${CASE_WORK}/low"
PATH="${CASE_WORK}/low-bin:${PATH}" run_expect 0 low "${SCRIPTS}/setup.sh" --clinical-data "$G"
output_has low '56 GB.*5 GB margin'
output_has low '52 GB.*5 GB margin'
[ ! -e "${G}/vep_cache/homo_sapiens/116_GRCh38" ] || fail 'low-space VEP installed'
[ -s "${G}/pcgr_data/20260620/data/grch38/README" ] || fail 'small bundle skipped with 40 GB free'

# The digest belongs to the original fixture. Serve a different valid tarball.
printf 'corrupt\n' >> "${FAKE_DOWNLOAD_DIR}/pcgr_ref_data.20260620.grch38.tgz"
G="${CASE_WORK}/corrupt"
run_expect 0 corrupt "${SCRIPTS}/setup.sh" --clinical-data "$G"
output_has corrupt 'sha256 checksum.*expected'
output_has corrupt '\[WARN\].*PCGR'
[ ! -e "${G}/pcgr_data/20260620" ] || fail 'corrupt bundle published'
[ ! -e "${G}/pcgr_data/pcgr_ref_data.20260620.grch38.tgz.part" ] || fail 'bad bundle kept'
[ -s "${G}/vep_cache/homo_sapiens/115_GRCh38/info.txt" ] || fail 'bundle failure stopped remaining install'
