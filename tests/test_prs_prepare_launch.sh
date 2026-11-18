#!/usr/bin/env bash
# PRS_PREPARE formats the score files from the samplesheet alone, so it must
# start at launch, while the samples are still aligned and called, and not
# wait for DeepVariant. This stub run (no container engine) has two rows that
# DeepVariant calls (FASTQ, BAM) and two VCF-only rows without a gVCF, one
# with the sex and one without, so the run needs both score sets. From the
# trace it checks that:
#   - PRS_PREPARE ran once per set, with_x and autosomes;
#   - each PRS_PREPARE task was submitted before any DEEPVARIANT task
#     completed (before the change, the first set waited for the first VCF);
#   - every sample still got its PRS_SUMMARY, so each joined its set.
#
# Needs nextflow on PATH. The Guard unit-test job has none, so there it says
# SKIP; nextflow.yml runs it with PRS_LAUNCH_REQUIRE=1, where a missing
# nextflow is a failure.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
if ! command -v nextflow > /dev/null 2>&1; then
  if [ "${PRS_LAUNCH_REQUIRE:-}" = 1 ]; then
    echo "FAIL: nextflow is not on PATH and PRS_LAUNCH_REQUIRE=1"; exit 1
  fi
  echo "SKIP: nextflow is not on PATH (nextflow.yml runs this test)"; exit 0
fi

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
S="${REPO}/assets/stub"

cat > "${WORK}/samplesheet.csv" <<EOF
sample,fastq_1,fastq_2,bam,bam_index,vcf,vcf_index,sex
dv_fastq,${S}/sample_R1.fastq.gz,${S}/sample_R2.fastq.gz,,,,,male
dv_bam,,,${S}/sample.bam,${S}/sample.bam.bai,,,female
vcf_sex,,,,,${S}/sample.vcf.gz,${S}/sample.vcf.gz.tbi,male
vcf_nosex,,,,,${S}/sample.vcf.gz,${S}/sample.vcf.gz.tbi,
EOF

cat > "${WORK}/trace.config" <<EOF
trace {
    enabled   = true
    overwrite = true
    raw       = true
    file      = '${WORK}/trace.txt'
    fields    = 'task_id,name,status,submit,complete'
}
// A stub task ends in about a second, and on a small runner the memory the
// first tasks ask for can hold PRS_PREPARE back for about as long. Each
// DEEPVARIANT task waits 20 s first, so the check measures the order the
// channels allow, not which stub happens to be quicker. Before the change
// PRS_PREPARE waits for a VCF, so it still comes after.
process {
    withName: 'DEEPVARIANT' {
        beforeScript = 'sleep 20'
    }
}
EOF

# The test profile with the stub scoring files; only prs is selected, so only
# the upstream processes and the four PRS processes run.
rc=0
(cd "$REPO" && nextflow -q run main.nf -profile test -stub -ansi-log false \
  -c "${WORK}/trace.config" -work-dir "${WORK}/work" --outdir "${WORK}/out" \
  --input "${WORK}/samplesheet.csv" --tools prs --pgs_scoring "${S}/db/pgs_scoring") > "${WORK}/nextflow.log" 2>&1 || rc=$?
if [ "$rc" -ne 0 ] || [ ! -s "${WORK}/trace.txt" ]; then
  cat "${WORK}/nextflow.log"
  echo "FAIL: the stub run exited ${rc}"; exit 1
fi

# The trace can list a finished task twice (docs/lessons-learned.md, "The
# trace file can list a finished task twice"): one entry per task_id, with the
# earliest submit and the latest complete time its rows give.
awk -F'\t' '
  NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
  {
    id = $col["task_id"]
    if (!(id in name)) { name[id] = $col["name"]; order[++n] = id; sub_[id] = "-"; done_[id] = "-" }
    t = $col["submit"]; if (t != "" && t != "-" && (sub_[id] == "-" || t + 0 < sub_[id] + 0)) sub_[id] = t
    t = $col["complete"]; if (t != "" && t != "-" && (done_[id] == "-" || t + 0 > done_[id] + 0)) done_[id] = t
  }
  # proc NAME: the process of a task name, "PRS:PRS_PREPARE (with_x)" -> PRS_PREPARE
  function proc(s) { sub(/ \(.*/, "", s); sub(/.*:/, "", s); return s }
  function ms(id, t, field) { if (t == "-") { print "FAIL: " name[id] " has no " field " time in the trace"; exit 1 } return t + 0 }
  END {
    for (k = 1; k <= n; k++) {
      id = order[k]; p = proc(name[id])
      if (p == "PRS_PREPARE") prep[++np] = id
      else if (p == "DEEPVARIANT") dv[++nd] = id
      else if (p == "PRS_SUMMARY") ns++
    }
    printf "tasks: %d, PRS_PREPARE %d, DEEPVARIANT %d, PRS_SUMMARY %d\n", n, np, nd, ns
    sets = ""; na = nx = 0
    for (k = 1; k <= np; k++) {
      s = name[prep[k]]; sub(/.*\(/, "", s); sub(/\)+$/, "", s)
      sets = sets (k > 1 ? ", " : "") s; na += (s == "autosomes"); nx += (s == "with_x")
    }
    if (np != 2 || na != 1 || nx != 1) bad[++nb] = "PRS_PREPARE ran for [" sets "], expected autosomes and with_x once each"
    if (nd != 2) bad[++nb] = (nd + 0) " DEEPVARIANT tasks, expected 2 (the FASTQ and the BAM row)"
    if (ns != 4) bad[++nb] = (ns + 0) " PRS_SUMMARY tasks, expected one per sample (4)"
    if (nd && np) {
      first = ""
      for (k = 1; k <= nd; k++) { t = ms(dv[k], done_[dv[k]], "complete"); if (first == "" || t < first) first = t }
      for (k = 1; k <= np; k++) {
        lag = (ms(prep[k], sub_[prep[k]], "submit") - first) / 1000
        printf "%s: submitted %.1f s %s the first DEEPVARIANT task completed\n", name[prep[k]], (lag < 0 ? -lag : lag), (lag >= 0 ? "after" : "before")
        if (lag >= 0) bad[++nb] = name[prep[k]] " waited for DeepVariant"
      }
    }
    for (k = 1; k <= nb; k++) print "FAIL:", bad[k]
    if (nb) exit 1
    print "OK: PRS_PREPARE started at launch, once per set, and every sample got its PRS summary"
  }' "${WORK}/trace.txt"
