#!/usr/bin/env bash
# --cpsr_secondary_findings on the Nextflow path: stub runs of main.nf
# (profile test_all, no containers) with --tools cpsr,hla_typing on the stub VCF+BAM row
# (hla_typing because test_all turns on --kir, which needs it; a --kir false
# on the command line reaches main.nf as the string 'false', which is true).
# CPSR's stub prints the options its script passes, so the task's .command.sh
# shows them:
#   1. the default: --secondary_findings is there;
#   2. --cpsr_secondary_findings false: it is not, and --panel_id 0 still is;
#   3. --cpsr_secondary_findings no: the run stops and names the parameter.
# The stubs only touch files; nothing here needs a database or an image.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

ST="${REPO}/assets/stub"
# stub_run NAME [ARGS...]: one stub run in its own launch directory; sets RC, LOG, DIR and CPSR_SH
stub_run() {
  local name=$1
  shift
  DIR="${CASE_TMP}/${name}"
  rm -rf "$DIR"
  mkdir -p "$DIR"
  printf 'sample,bam,bam_index,vcf,vcf_index,sex\nsf1,%s,%s,%s,%s,male\n' \
    "${ST}/sample.bam" "${ST}/sample.bam.bai" "${ST}/sample.vcf.gz" "${ST}/sample.vcf.gz.tbi" > "${DIR}/samplesheet.csv"
  echo "+ nextflow run main.nf -profile test_all -stub --tools cpsr,hla_typing (${name}) $*"
  ( cd "$DIR" && nextflow run "${REPO}/main.nf" -profile test_all -stub -ansi-log false \
      -work-dir "${DIR}/work" --input "${DIR}/samplesheet.csv" --outdir "${DIR}/out" --tools cpsr,hla_typing "$@" ) > "${DIR}/run.log" 2>&1
  RC=$?
  LOG=$(cat "${DIR}/run.log")
  echo "$LOG"
  cp "${DIR}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.${name}.nextflow.log" 2>/dev/null
  # CPSR's task: the one whose script names CPSR's output files
  CPSR_SH=$(find "${DIR}/work" -name .command.sh -exec grep -l 'cpsr.grch38.html' {} + 2>/dev/null | head -n 1)
}

# --- 1. the default ---------------------------------------------------------------------
stub_run default
check_eq "default: the run exits 0" "$RC" 0
check "default: CPSR's task ran" test -n "$CPSR_SH"
check "default: CPSR gets --secondary_findings" grep -q -- '--panel_id 0 --secondary_findings --force_overwrite' "${CPSR_SH:-/dev/null}"

# --- 2. switched off --------------------------------------------------------------------
stub_run off --cpsr_secondary_findings false
check_eq "off: the run exits 0" "$RC" 0
check "off: CPSR's task ran" test -n "$CPSR_SH"
check "off: CPSR still gets --panel_id 0 and --force_overwrite" grep -q -- '--panel_id 0 .*--force_overwrite' "${CPSR_SH:-/dev/null}"
check "off: CPSR does not get --secondary_findings" test -n "$CPSR_SH" -a -z "$(grep -- '--secondary_findings' "${CPSR_SH:-/dev/null}")"

# --- 3. a value that is neither true nor false -------------------------------------------
stub_run bad --cpsr_secondary_findings no
check "bad: the run stops (exit ${RC})" test "$RC" -ne 0
check "bad: the message names the parameter" has "--cpsr_secondary_findings must be true or false, got 'no'" "$LOG"

finish
