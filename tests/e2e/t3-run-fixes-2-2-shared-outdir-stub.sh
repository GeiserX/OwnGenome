#!/usr/bin/env bash
# Stub runs of main.nf (profile test_all, no containers) on the stub VCF+BAM row.
#   1. Two runs of one sample each, SAMPLE_A then SAMPLE_B, into one
#      --outdir (what run-all.sh does: --outdir is GENOME_DIR). MultiQC and
#      SOMALIER_RELATE used to publish to multiqc/ and somalier/ beside the
#      sample folders, so SAMPLE_B's run replaced SAMPLE_A's report and
#      tables. Now each sample keeps its own (<id>/multiqc, and
#      <id>/qc/somalier/<id>.samples.tsv and .pairs.tsv, the names step 33
#      gives them) and a run of one sample writes no shared folder.
#   2. One run of both samples into another --outdir: each sample still gets
#      its own MultiQC report and somalier tables, and the cohort's report
#      and somalier tables (every pair) go to multiqc/ and somalier/, which
#      only a run of two samples or more writes.
#   3. PRS failing (planted through beforeScript) with ROH waiting for its
#      ancestry table (--roh_af_panel with prs, --pgs_scoring and
#      --ancestry_ref): PRS is report-only, so the run goes on, ROH runs
#      without the table, the run exits 0, and failed_tasks.tsv names PRS as
#      ignored. It used to stop the run with ROH never started.
#   4. A run of one sample into a results folder that a run before 1.0
#      wrote into: the old multiqc/ and somalier/ are left as they are, and
#      the log says they are old.
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
# samples_in FILE: the sample ids of a somalier relate samples table, sorted, comma-separated
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

# --- 1. SAMPLE_A then SAMPLE_B, one sample each, into one --outdir ------------------------------
SHARED="${CASE_TMP}/shared"
rm -rf "$SHARED"
stub_run s1 "$SHARED" SAMPLE_A "${QC_TOOLS[@]}"
check_eq "SAMPLE_A's run exits 0" "$RC" 0
touch "${SHARED}/SAMPLE_A/multiqc/multiqc_report.html.s1" 2>/dev/null
stub_run s2 "$SHARED" SAMPLE_B "${QC_TOOLS[@]}"
check_eq "SAMPLE_B's run exits 0" "$RC" 0
check "SAMPLE_A's MultiQC report is still in SAMPLE_A's folder" test -e "${SHARED}/SAMPLE_A/multiqc/multiqc_report.html"
check "SAMPLE_A's MultiQC folder is SAMPLE_A's run's (SAMPLE_B's run left it alone)" test -e "${SHARED}/SAMPLE_A/multiqc/multiqc_report.html.s1"
check "SAMPLE_B's MultiQC report is in SAMPLE_B's folder" test -e "${SHARED}/SAMPLE_B/multiqc/multiqc_report.html"
check "a run of one sample writes no shared multiqc/" test ! -e "${SHARED}/multiqc"
check_eq "SAMPLE_A's somalier tables hold SAMPLE_A, after SAMPLE_B's run" "$(samples_in "${SHARED}/SAMPLE_A/qc/somalier/SAMPLE_A.samples.tsv")" SAMPLE_A
check_eq "SAMPLE_B's somalier tables hold SAMPLE_B" "$(samples_in "${SHARED}/SAMPLE_B/qc/somalier/SAMPLE_B.samples.tsv")" SAMPLE_B
check_eq "SAMPLE_A's somalier tables have the names step 33 gives them, and no other" \
  "$(find "${SHARED}/SAMPLE_A/qc/somalier" -maxdepth 1 -type f \( -name '*.tsv' -o -name '*.html' \) -exec basename {} \; 2>/dev/null | LC_ALL=C sort | paste -sd, -)" \
  SAMPLE_A.pairs.tsv,SAMPLE_A.samples.tsv
check "a run of one sample writes no shared somalier/" test ! -e "${SHARED}/somalier"
check_eq "one sample: no cohort somalier task" "$(ran "$SHARED" SOMALIER_RELATE_COHORT)" 0
check_eq "one sample: its own somalier task" "$(ran "$SHARED" SOMALIER_RELATE)" 1
check_eq "one sample: no cohort MultiQC task" "$(ran "$SHARED" MULTIQC_COHORT)" 0

# --- 2. SAMPLE_A and SAMPLE_B in one run --------------------------------------------------------
COHORT="${CASE_TMP}/cohort"
rm -rf "$COHORT"
stub_run both "$COHORT" SAMPLE_A,SAMPLE_B "${QC_TOOLS[@]}"
check_eq "the two-sample run exits 0" "$RC" 0
check_eq "two samples: one MultiQC task each" "$(ran "$COHORT" MULTIQC)" 2
check_eq "two samples: one cohort MultiQC task" "$(ran "$COHORT" MULTIQC_COHORT)" 1
check_eq "two samples: one cohort somalier task" "$(ran "$COHORT" SOMALIER_RELATE_COHORT)" 1
check_eq "two samples: one own somalier task each" "$(ran "$COHORT" SOMALIER_RELATE)" 2
for id in SAMPLE_A SAMPLE_B; do
  check "two samples: ${id}'s own MultiQC report" test -e "${COHORT}/${id}/multiqc/multiqc_report.html"
done
check "two samples: the cohort MultiQC report in multiqc/" test -e "${COHORT}/multiqc/multiqc_report.html"
check_eq "two samples: the somalier tables in somalier/ hold both" "$(samples_in "${COHORT}/somalier/somalier.samples.tsv")" SAMPLE_A,SAMPLE_B
for id in SAMPLE_A SAMPLE_B; do
  check_eq "two samples: ${id}'s own somalier tables hold ${id} only" "$(samples_in "${COHORT}/${id}/qc/somalier/${id}.samples.tsv")" "$id"
done

# --- 3. PRS fails, ROH waits for its ancestry table ---------------------------------
printf "process { withName: 'PRS' { beforeScript = 'echo planted failure >&2; exit 1' } }\n" > "${CASE_TMP}/plant-prs.config"
PRSOUT="${CASE_TMP}/prs-fails"
rm -rf "$PRSOUT"
stub_run prs "$PRSOUT" SAMPLE_A --tools prs,roh -params-file "${CASE_TMP}/no-kir.json" -c "${CASE_TMP}/plant-prs.config" --roh_af_panel "${ST}/sample.vcf.gz"
check_eq "PRS failing: the run goes on and exits 0" "$RC" 0
check_eq "PRS ran (and failed)" "$(ran "$PRSOUT" PRS)" 1
check_eq "ROH ran without the ancestry table" "$(ran "$PRSOUT" ROH)" 1
check "ROH's summary is published" test -e "${PRSOUT}/SAMPLE_A/roh/SAMPLE_A_roh_summary.txt"
# shellcheck disable=SC2016  # the awk program is awk's, not the shell's
check "failed_tasks.tsv names PRS as ignored" \
  awk -F'\t' '$1 ~ /(^|:)PRS$/ && $2 == "ignored" {f = 1} END {exit !f}' "${DIR}/failed_tasks.tsv"

# --- 4. A results folder a run before 1.0 wrote into ---------------------------------
# Such a run published multiqc/ and somalier/ (with a versions.yml) beside the
# sample folders. A run of one sample leaves them alone now, and the log says
# they are old; the runs above, with no such folder, say nothing.
OLD="${CASE_TMP}/old-outdir"
rm -rf "$OLD"
mkdir -p "${OLD}/multiqc" "${OLD}/somalier"
echo "old report" > "${OLD}/multiqc/multiqc_report.html"
printf '"SOMALIER_RELATE":\n    somalier: 0.2.19\n' > "${OLD}/somalier/versions.yml"
echo "old table" > "${OLD}/somalier/somalier.samples.tsv"
before=$(cd "$OLD" && cat multiqc/multiqc_report.html somalier/versions.yml somalier/somalier.samples.tsv | cksum)
stub_run old "$OLD" SAMPLE_A "${QC_TOOLS[@]}"
check_eq "a run of one sample into an old results folder exits 0" "$RC" 0
check "the log says multiqc/ and somalier/ are left from a run before 1.0" \
  grep -qF "${OLD}/somalier are left from a run before 1.0" "${DIR}/run.log"
check_eq "the old multiqc/ and somalier/ are left as they were" \
  "$(cd "$OLD" && cat multiqc/multiqc_report.html somalier/versions.yml somalier/somalier.samples.tsv | cksum)" "$before"
check "the sample's own MultiQC report is in its folder" test -e "${OLD}/SAMPLE_A/multiqc/multiqc_report.html"
# shellcheck disable=SC2016  # bash -c expands its own "$@"
check "a results folder with no old somalier/ gets no such line" \
  bash -c '! grep -qF "are left from a run before 1.0" "$@"' _ "${CASE_TMP}/s1/run.log" "${CASE_TMP}/s2/run.log" "${CASE_TMP}/both/run.log"

rm -rf "$SHARED" "$COHORT" "$PRSOUT" "$OLD" "${CASE_TMP}/s1" "${CASE_TMP}/s2" "${CASE_TMP}/both" "${CASE_TMP}/prs" "${CASE_TMP}/old"
finish
