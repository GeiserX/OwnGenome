#!/usr/bin/env bash
# Stub runs of main.nf (profile test_all, no containers) on the stub VCF+BAM row.
#   1. Two runs of one sample each, S1 then S2, into one --outdir (what
#      run-all.sh does: --outdir is GENOME_DIR). MultiQC and
#      SOMALIER_RELATE used to publish to multiqc/ and somalier/ beside the
#      sample folders, so S2's run replaced S1's report and tables. Now each
#      sample keeps its own (<id>/multiqc, <id>/qc/somalier) and a run of one
#      sample writes no shared folder.
#   2. One run of both samples into another --outdir: each sample still gets
#      its own MultiQC report and somalier tables, and the cohort's report
#      and somalier tables (every pair) go to multiqc/ and somalier/, which
#      only a run of two samples or more writes.
#   3. PRS failing (planted through beforeScript) with ROH waiting for its
#      ancestry table (--roh_af_panel with prs, --pgs_scoring and
#      --ancestry_ref): PRS is report-only, so the run goes on, ROH runs
#      without the table, the run exits 0, and failed_tasks.tsv names PRS as
#      ignored. It used to stop the run with ROH never started.
# The stubs only touch files; nothing here needs a database or an image.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

ST="${REPO}/assets/stub"
# stub_run NAME OUTDIR IDS [ARGS...]: one stub run of the samples IDS (comma-separated),
# in its own launch directory; sets RC and DIR
stub_run() {
  local name=$1 out=$2 ids=$3 id
  shift 3
  DIR="${CASE_TMP}/${name}"
  rm -rf "$DIR"
  mkdir -p "$DIR"
  echo 'sample,bam,bam_index,vcf,vcf_index,sex' > "${DIR}/samplesheet.csv"
  for id in ${ids//,/ }; do
    printf '%s,%s,%s,%s,%s,male\n' "$id" \
      "${ST}/sample.bam" "${ST}/sample.bam.bai" "${ST}/sample.vcf.gz" "${ST}/sample.vcf.gz.tbi" >> "${DIR}/samplesheet.csv"
  done
  echo "+ nextflow run main.nf -profile test_all -stub (${name}, samples ${ids}, --outdir ${out#"${CASE_TMP}/"}) $*"
  ( cd "$DIR" && nextflow run "${REPO}/main.nf" -profile test_all -stub -ansi-log false \
      -work-dir "${DIR}/work" --input "${DIR}/samplesheet.csv" --outdir "$out" "$@" ) > "${DIR}/run.log" 2>&1
  RC=$?
  tail -n 25 "${DIR}/run.log"
  cp "${DIR}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.${name}.nextflow.log" 2>/dev/null
}
# samples_in FILE: the sample ids of a somalier.samples.tsv, sorted, comma-separated
samples_in() { awk -F'\t' 'NR > 1 {print $2}' "$1" 2>/dev/null | LC_ALL=C sort | paste -sd, -; }
# ran NAME: tasks of process NAME in the run's trace, in any state
ran() {
  local trace
  trace=$(find "$1/pipeline_info" -name 'trace_*.txt' 2>/dev/null | LC_ALL=C sort | awk 'END {print}')
  awk -F'\t' -v p="$2" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    { n = $c["name"]; sub(/ \(.*/, "", n); sub(/.*:/, "", n); if (n == p) k++ }
    END {print k + 0}' "${trace:-/dev/null}" 2>/dev/null
}
# The profile turns KIR on, which needs hla_typing in --tools: a params file
# turns it off (on the command line, --kir false arrives as a string and stays on).
echo '{"kir": false}' > "${CASE_TMP}/no-kir.json"
QC_TOOLS=(--tools "mosdepth,sample_qc,multiqc" -params-file "${CASE_TMP}/no-kir.json")

# --- 1. S1 then S2, one sample each, into one --outdir ------------------------------
SHARED="${CASE_TMP}/shared"
rm -rf "$SHARED"
stub_run s1 "$SHARED" S1 "${QC_TOOLS[@]}"
check_eq "S1's run exits 0" "$RC" 0
touch "${SHARED}/S1/multiqc/multiqc_report.html.s1" 2>/dev/null
stub_run s2 "$SHARED" S2 "${QC_TOOLS[@]}"
check_eq "S2's run exits 0" "$RC" 0
check "S1's MultiQC report is still in S1's folder" test -e "${SHARED}/S1/multiqc/multiqc_report.html"
check "S1's MultiQC folder is S1's run's (S2's run left it alone)" test -e "${SHARED}/S1/multiqc/multiqc_report.html.s1"
check "S2's MultiQC report is in S2's folder" test -e "${SHARED}/S2/multiqc/multiqc_report.html"
check "a run of one sample writes no shared multiqc/" test ! -e "${SHARED}/multiqc"
check_eq "S1's somalier tables hold S1, after S2's run" "$(samples_in "${SHARED}/S1/qc/somalier/somalier.samples.tsv")" S1
check_eq "S2's somalier tables hold S2" "$(samples_in "${SHARED}/S2/qc/somalier/somalier.samples.tsv")" S2
check "a run of one sample writes no shared somalier/" test ! -e "${SHARED}/somalier"
check_eq "one sample: no cohort somalier task" "$(ran "$SHARED" SOMALIER_RELATE_COHORT)" 0
check_eq "one sample: its own somalier task" "$(ran "$SHARED" SOMALIER_RELATE)" 1
check_eq "one sample: no cohort MultiQC task" "$(ran "$SHARED" MULTIQC_COHORT)" 0

# --- 2. S1 and S2 in one run --------------------------------------------------------
COHORT="${CASE_TMP}/cohort"
rm -rf "$COHORT"
stub_run both "$COHORT" S1,S2 "${QC_TOOLS[@]}"
check_eq "the two-sample run exits 0" "$RC" 0
check_eq "two samples: one MultiQC task each" "$(ran "$COHORT" MULTIQC)" 2
check_eq "two samples: one cohort MultiQC task" "$(ran "$COHORT" MULTIQC_COHORT)" 1
check_eq "two samples: one cohort somalier task" "$(ran "$COHORT" SOMALIER_RELATE_COHORT)" 1
check_eq "two samples: one own somalier task each" "$(ran "$COHORT" SOMALIER_RELATE)" 2
for id in S1 S2; do
  check "two samples: ${id}'s own MultiQC report" test -e "${COHORT}/${id}/multiqc/multiqc_report.html"
done
check "two samples: the cohort MultiQC report in multiqc/" test -e "${COHORT}/multiqc/multiqc_report.html"
check_eq "two samples: the somalier tables in somalier/ hold both" "$(samples_in "${COHORT}/somalier/somalier.samples.tsv")" S1,S2
for id in S1 S2; do
  check_eq "two samples: ${id}'s own somalier tables hold ${id} only" "$(samples_in "${COHORT}/${id}/qc/somalier/somalier.samples.tsv")" "$id"
done

# --- 3. PRS fails, ROH waits for its ancestry table ---------------------------------
printf "process { withName: 'PRS' { beforeScript = 'echo planted failure >&2; exit 1' } }\n" > "${CASE_TMP}/plant-prs.config"
PRSOUT="${CASE_TMP}/prs-fails"
rm -rf "$PRSOUT"
stub_run prs "$PRSOUT" S1 --tools prs,roh -params-file "${CASE_TMP}/no-kir.json" -c "${CASE_TMP}/plant-prs.config" --roh_af_panel "${ST}/sample.vcf.gz"
check_eq "PRS failing: the run goes on and exits 0" "$RC" 0
check_eq "PRS ran (and failed)" "$(ran "$PRSOUT" PRS)" 1
check_eq "ROH ran without the ancestry table" "$(ran "$PRSOUT" ROH)" 1
check "ROH's summary is published" test -e "${PRSOUT}/S1/roh/S1_roh_summary.txt"
# shellcheck disable=SC2016  # the awk program is awk's, not the shell's
check "failed_tasks.tsv names PRS as ignored" \
  awk -F'\t' '$1 ~ /(^|:)PRS$/ && $2 == "ignored" {f = 1} END {exit !f}' "${DIR}/failed_tasks.tsv"

rm -rf "$SHARED" "$COHORT" "$PRSOUT" "${CASE_TMP}/s1" "${CASE_TMP}/s2" "${CASE_TMP}/both" "${CASE_TMP}/prs"
finish
