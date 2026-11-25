#!/usr/bin/env bash
# --pgs_academic on the Nextflow path: stub runs of main.nf (profile test, no
# containers) with --tools prs on the stub VCF row. PRS_PREPARE's stub prints
# the og prs format option the parameter sets, so the task's .command.sh
# shows it:
#   1. the default: no --academic (the four academic-only scores are left out);
#   2. --pgs_academic true: --academic;
#   3. --pgs_academic yes: the run stops and names the parameter.
# The stubs only touch files; nothing here needs a database or an image.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

ST="${REPO}/assets/stub"
# stub_run NAME [ARGS...]: one stub run in its own launch directory; sets RC, LOG and PREP (the .command.sh files of PRS_PREPARE)
stub_run() {
  local name=$1 dir
  shift
  dir="${CASE_TMP}/${name}"
  rm -rf "$dir"
  mkdir -p "$dir"
  printf 'sample,vcf,vcf_index,sex\npa1,%s,%s,male\n' "${ST}/sample.vcf.gz" "${ST}/sample.vcf.gz.tbi" > "${dir}/samplesheet.csv"
  echo "+ nextflow run main.nf -profile test -stub --tools prs (${name}) $*"
  ( cd "$dir" && nextflow run "${REPO}/main.nf" -profile test -stub -ansi-log false \
      -work-dir "${dir}/work" --input "${dir}/samplesheet.csv" --outdir "${dir}/out" --tools prs \
      --pgs_scoring "${ST}/db/pgs_scoring" "$@" ) > "${dir}/run.log" 2>&1
  RC=$?
  LOG=$(cat "${dir}/run.log")
  echo "$LOG"
  cp "${dir}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.${name}.nextflow.log" 2>/dev/null
  # PRS_PREPARE's tasks: the ones whose script prints og prs format's options
  mapfile -t PREP < <(find "${dir}/work" -name .command.sh -exec grep -l 'og prs format options:' {} + 2>/dev/null)
}

# --- 1. the default ---------------------------------------------------------------------
stub_run default
check_eq "default: the run exits 0" "$RC" 0
check "default: PRS_PREPARE ran" test "${#PREP[@]}" -gt 0
check "default: og prs format does not get --academic" \
  test -z "$(grep -h -- 'og prs format options:.*--academic' "${PREP[@]}" /dev/null)"

# --- 2. switched on ---------------------------------------------------------------------
stub_run on --pgs_academic true
check_eq "on: the run exits 0" "$RC" 0
check "on: PRS_PREPARE ran" test "${#PREP[@]}" -gt 0
check "on: og prs format gets --academic in every PRS_PREPARE task" \
  test "$(grep -l -- 'og prs format options: --academic' "${PREP[@]}" /dev/null | grep -c .)" -eq "${#PREP[@]}"

# --- 3. a value that is neither true nor false -------------------------------------------
stub_run bad --pgs_academic yes
check "bad: the run stops (exit ${RC})" test "$RC" -ne 0
check "bad: the message names the parameter" has "--pgs_academic takes no value \\(or true / false\\), got 'yes'" "$LOG"

finish
