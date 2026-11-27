#!/usr/bin/env bash
# validate-setup.sh and a Nextflow work directory on another disk (-w).
# run-all.sh hands validate-setup.sh the work directory its run gives
# Nextflow, as NXF_WORK: -w or -work-dir, else NXF_WORK, else
# <sample>/nextflow/work, a relative one under <sample>/nextflow. Before,
# validate-setup.sh looked for a resume only at <sample>/nextflow/work and
# measured free space only where GENOME_DIR is, so:
#   1. a resume whose work directory is elsewhere, with GENOME_DIR under the
#      200 GB minimum, was a FAIL and stopped run-all.sh; it is a WARN;
#   2. a work disk under the minimum was never seen; it is a FAIL on a fresh
#      run and a WARN on a resume, and a work directory on GENOME_DIR's disk
#      adds no second line;
#   3. run-all.sh passes -w (absolute or relative) on, so the resume of case 1
#      starts Nextflow instead of stopping at validation.
# A df in front of the fake one answers per path: the work disk's paths get
# FAKE_FREE_WORK_GB on mount /work-disk, every other path FAKE_FREE_GB on /.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
G=$GENOME_DIR
WDISK="${CASE_WORK}/workdisk"
seed_reference "$G"
seed_clinvar "$G"
seed_sample "$G" sample1
mkdir -p "$WDISK"

mkdir -p "${CASE_WORK}/df-bin"
cat > "${CASE_WORK}/df-bin/df" <<SHIM
#!/usr/bin/env bash
p=\${*: -1}
echo "Filesystem 1024-blocks Used Available Capacity Mounted on"
case "\$p" in
  ${WDISK}|${WDISK}/*) echo "workfs 4294967296 0 \$(( \${FAKE_FREE_WORK_GB:?} * 1048576 )) 0% /work-disk" ;;
  *) echo "fakefs 4294967296 0 \$(( \${FAKE_FREE_GB:?} * 1048576 )) 0% /" ;;
esac
SHIM
chmod +x "${CASE_WORK}/df-bin/df"
# with_free GENOME_GB WORK_GB COMMAND...: run COMMAND with that much free on each disk.
with_free() { local g=$1 w=$2; shift 2; PATH="${CASE_WORK}/df-bin:${PATH}" FAKE_FREE_GB=$g FAKE_FREE_WORK_GB=$w "$@"; }

# 1. GENOME_DIR under the minimum, the run's work directory elsewhere and there: a resume.
WD="${WDISK}/owngenome/sample1"
mkdir -p "${WD}/ab"
with_free 150 600 run_expect 0 resume-elsewhere env NXF_WORK="$WD" "${SCRIPTS}/validate-setup.sh" sample1
output_has resume-elsewhere "\[WARN\].*Free disk space: 150 GB, under the 200 GB minimum, but ${WD} exists: a resume"
output_lacks resume-elsewhere '\[FAIL\].*Free disk space'
output_has resume-elsewhere "\[OK\].*Free disk space for Nextflow's work directory ${WD}: 600 GB in /work-disk"
# Control: without NXF_WORK the same run is not a resume, so the old FAIL stands.
with_free 150 600 run_expect 1 resume-unseen "${SCRIPTS}/validate-setup.sh" sample1
output_has resume-unseen '\[FAIL\].*Free disk space: 150 GB, under the 200 GB minimum \(500 GB'

# 2. The work disk under the minimum: a FAIL on a fresh run, a WARN on a resume.
FRESH="${WDISK}/owngenome/fresh"
with_free 600 150 run_expect 1 work-low env NXF_WORK="$FRESH" "${SCRIPTS}/validate-setup.sh" sample1
output_has work-low "\[FAIL\].*Free disk space for Nextflow's work directory ${FRESH}: 150 GB, under the 200 GB minimum"
output_has work-low '\[OK\].*Free disk space: 600 GB in /'
# A full work disk (0 free) is measured too, not skipped.
with_free 600 0 run_expect 1 work-full env NXF_WORK="$FRESH" "${SCRIPTS}/validate-setup.sh" sample1
output_has work-full "\[FAIL\].*Free disk space for Nextflow's work directory ${FRESH}: 0 GB, under the 200 GB minimum"
# An empty work directory, made ahead of a first run, is no resume.
mkdir -p "$FRESH"
with_free 600 150 run_expect 1 work-empty env NXF_WORK="$FRESH" "${SCRIPTS}/validate-setup.sh" sample1
output_has work-empty "\[FAIL\].*Free disk space for Nextflow's work directory ${FRESH}: 150 GB, under the 200 GB minimum"
rmdir "$FRESH"
with_free 600 150 run_expect 0 work-low-resume env NXF_WORK="$WD" "${SCRIPTS}/validate-setup.sh" sample1
output_has work-low-resume "\[WARN\].*Free disk space for Nextflow's work directory ${WD}: 150 GB, under the 200 GB minimum, but ${WD} exists"
# The default work directory sits on GENOME_DIR's disk: one line, as before.
with_free 600 150 run_expect 0 same-disk "${SCRIPTS}/validate-setup.sh" sample1
output_lacks same-disk "work directory"
# A relative NXF_WORK is under <sample>/nextflow, where run-all.sh starts Nextflow.
mkdir -p "${G}/sample1/nextflow/rel/ab"
with_free 150 600 run_expect 0 relative env NXF_WORK=rel "${SCRIPTS}/validate-setup.sh" sample1
output_has relative "but ${G}/sample1/nextflow/rel exists"

# 3. run-all.sh passes -w on: the resume of case 1 reaches Nextflow.
use_output_hook
cat > "${CASE_WORK}/tools-hook" <<'HOOK'
#!/usr/bin/env bash
case "${*:2}" in
  "cat /usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json") echo '[{"LocusId": "HTT"}]' ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/tools-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/tools-hook"
: > "$FAKE_DOCKER_LOG"
with_free 150 600 run_expect 0 run-all-w "${SCRIPTS}/run-all.sh" sample1 male -w "$WD"
output_has run-all-w "\[WARN\].*Free disk space: 150 GB, under the 200 GB minimum, but ${WD} exists"
output_lacks run-all-w 'setup validation failed'
grep -qF -- "-w ${WD}" <(grep '^nextflow :: ' "$FAKE_DOCKER_LOG") || fail "nextflow did not get -w ${WD}"
mkdir -p "${G}/sample1/nextflow/w2/ab"
: > "$FAKE_DOCKER_LOG"
with_free 150 600 run_expect 0 run-all-w-relative "${SCRIPTS}/run-all.sh" sample1 male -work-dir w2
output_has run-all-w-relative "but ${G}/sample1/nextflow/w2 exists"
# Control: the same run without -w is no resume (no <sample>/nextflow/work), so validation stops it.
rm -rf "${G}/sample1/nextflow/rel" "${G}/sample1/nextflow/w2"
with_free 150 600 run_expect 1 run-all-no-w "${SCRIPTS}/run-all.sh" sample1 male
output_has run-all-no-w 'setup validation failed'
echo "PASS validate-setup.sh sees a -w work directory: a resume there is a WARN, its disk is checked, run-all.sh passes it on"
