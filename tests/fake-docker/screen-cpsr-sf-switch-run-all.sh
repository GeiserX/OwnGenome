#!/usr/bin/env bash
# run-all.sh reads CPSR_SECONDARY_FINDINGS as scripts/17-cpsr.sh does, so a
# user can turn CPSR's secondary findings off on the default path too.
# Runs with TOOLS=cpsr and CPSR's data in place:
#   1. unset: no --cpsr_secondary_findings on the nextflow command line (the
#      module's default, true, passes --secondary_findings);
#   2. CPSR_SECONDARY_FINDINGS=false: --cpsr_secondary_findings false reaches
#      nextflow, the console says so, and logs/run_status.tsv keeps the switch
#      in 'meta switches' and its value in 'meta command' (og recheck --rerun
#      replays that command);
#   3. CPSR_SECONDARY_FINDINGS=true: nothing passed, as unset;
#   4. CPSR_SECONDARY_FINDINGS=false with cpsr not in TOOLS: nothing passed;
#   5. CPSR_SECONDARY_FINDINGS=no: exit 2 before nextflow, naming the switch.
# A fake nextflow stands in (scripts/ci/fake-docker/nextflow logs its arguments).
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
ST="${G}/sample1/logs/run_status.tsv"
seed_clinvar "$G"
seed_sample "$G" sample1
use_output_hook   # the reports at the end read the files their containers write
# CPSR's data where run-all.sh looks for it (its plan line for step 17)
mkdir -p "${G}/vep_cache/homo_sapiens/${PCGR_VEP_CACHE_RELEASE}_GRCh38" "${G}/pcgr_data/${PCGR_DATA_BUNDLE}/data"
printf 'placeholder\n' > "${G}/vep_cache/homo_sapiens/${PCGR_VEP_CACHE_RELEASE}_GRCh38/info.txt"

MISSES=0
miss() { echo "ASSERT FAIL: $*" >&2; MISSES=$((MISSES + 1)); }
nf_line() { awk '/^nextflow :: / {l = $0} END {print l}' "$FAKE_DOCKER_LOG"; }
nf_calls() { grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG" || true; }
meta() { awk -F'\t' -v k="$1" '$1 == "meta" && $2 == k {print $3}' "$ST"; }

# 1. unset: the default
run_expect 0 sf-default env TOOLS=cpsr "${SCRIPTS}/run-all.sh" sample1 male
output_has sf-default '^  17 CPSR +runs$'
if nf_line | grep -q -- '--cpsr_secondary_findings'; then miss "sf-default: nextflow was given --cpsr_secondary_findings: $(nf_line)"; fi
output_lacks sf-default 'CPSR_SECONDARY_FINDINGS=false'
meta switches | tr ' ' '\n' | grep -qx CPSR_SECONDARY_FINDINGS || miss "sf-default: 'meta switches' does not name CPSR_SECONDARY_FINDINGS: $(meta switches)"
if meta command | grep -q 'CPSR_SECONDARY_FINDINGS='; then miss "sf-default: 'meta command' sets CPSR_SECONDARY_FINDINGS: $(meta command)"; fi

# 2. false: off on the Nextflow path, and recorded for the replay
run_expect 0 sf-off env TOOLS=cpsr CPSR_SECONDARY_FINDINGS=false "${SCRIPTS}/run-all.sh" sample1 male
output_has sf-off '^  17 CPSR +runs$'
nf_line | grep -q -- ' --cpsr_secondary_findings false' || miss "sf-off: nextflow was not given --cpsr_secondary_findings false: $(nf_line)"
output_has sf-off 'CPSR_SECONDARY_FINDINGS=false: CPSR reports the cancer genes only'
meta switches | tr ' ' '\n' | grep -qx CPSR_SECONDARY_FINDINGS || miss "sf-off: 'meta switches' does not name CPSR_SECONDARY_FINDINGS: $(meta switches)"
meta command | grep -q ' CPSR_SECONDARY_FINDINGS=false [^ ]*/run-all.sh sample1 male$' \
  || miss "sf-off: 'meta command' does not carry CPSR_SECONDARY_FINDINGS=false: $(meta command)"

# 3. true: as unset
run_expect 0 sf-on env TOOLS=cpsr CPSR_SECONDARY_FINDINGS=true "${SCRIPTS}/run-all.sh" sample1 male
if nf_line | grep -q -- '--cpsr_secondary_findings'; then miss "sf-on: nextflow was given --cpsr_secondary_findings: $(nf_line)"; fi

# 4. false, but CPSR does not run: nothing to switch off
run_expect 0 sf-nocpsr env TOOLS=pharmcat CPSR_SECONDARY_FINDINGS=false "${SCRIPTS}/run-all.sh" sample1 male
output_has sf-nocpsr '^  17 CPSR +skipped +\(not in TOOLS\)$'
if nf_line | grep -q -- '--cpsr_secondary_findings'; then miss "sf-nocpsr: nextflow was given --cpsr_secondary_findings without CPSR: $(nf_line)"; fi

# 5. anything else: refused before nextflow
before=$(nf_calls)
run_expect 2 sf-bad env TOOLS=cpsr CPSR_SECONDARY_FINDINGS=no "${SCRIPTS}/run-all.sh" sample1 male
output_has sf-bad "CPSR_SECONDARY_FINDINGS must be true or false, got 'no'"
[ "$(nf_calls)" = "$before" ] || miss "sf-bad: nextflow ran with CPSR_SECONDARY_FINDINGS=no"

[ "$MISSES" -eq 0 ] || fail "${MISSES} result(s) above are wrong"
echo "CPSR_SECONDARY_FINDINGS=false reaches the pipeline as --cpsr_secondary_findings false and is kept for the replay."
