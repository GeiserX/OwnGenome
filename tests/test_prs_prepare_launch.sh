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

python3 - "${WORK}/trace.txt" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))

def tasks(proc):
    return [r for r in rows if r["name"].split(" (")[0].split(":")[-1] == proc]

def ms(r, field):
    v = r[field]
    if not v or v == "-":
        sys.exit(f"FAIL: {r['name']} has no {field} time in the trace")
    return int(v)

prepare, dv, summary = tasks("PRS_PREPARE"), tasks("DEEPVARIANT"), tasks("PRS_SUMMARY")
print(f"tasks: {len(rows)}, PRS_PREPARE {len(prepare)}, DEEPVARIANT {len(dv)}, PRS_SUMMARY {len(summary)}")
bad = []
sets = sorted(r["name"].split("(")[-1].rstrip(")") for r in prepare)
if sets != ["autosomes", "with_x"]:
    bad.append(f"PRS_PREPARE ran for {sets}, expected autosomes and with_x once each")
if len(dv) != 2:
    bad.append(f"{len(dv)} DEEPVARIANT tasks, expected 2 (the FASTQ and the BAM row)")
if len(summary) != 4:
    bad.append(f"{len(summary)} PRS_SUMMARY tasks, expected one per sample (4)")
if dv and prepare:
    first_dv = min(ms(r, "complete") for r in dv)
    for r in prepare:
        lag = (ms(r, "submit") - first_dv) / 1000
        word = "after" if lag >= 0 else "before"
        print(f"{r['name']}: submitted {abs(lag):.1f} s {word} the first DEEPVARIANT task completed")
        if lag >= 0:
            bad.append(f"{r['name']} waited for DeepVariant")
for b in bad:
    print("FAIL:", b)
if bad:
    sys.exit(1)
print("OK: PRS_PREPARE started at launch, once per set, and every sample got its PRS summary")
PY
