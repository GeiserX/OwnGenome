#!/usr/bin/env bash
# run-all.sh with HLA_QC=true passes --hla_qc true to the pipeline when
# hla_typing is selected, and records step 08 (HLA typing) from OptiType's two
# processes as well: HLA_QC_READS and HLA_QC are not named after the step's
# tool (hla_typing), so each must be counted under it, or an OptiType task that
# failed leaves step 08 'ok'. Runs with TOOLS=hla_typing and the pipeline exit 3:
#   1. HLA_QC=true, every task completed: --hla_qc true reaches nextflow and
#      step 08 is 'ok' (the control);
#   2. HLA_QC=true, HLA_QC failed: step 08 is 'failed';
#   3. HLA_QC=true, HLA_QC_READS failed: step 08 is 'failed';
#   4. HLA_QC unset, then HLA_QC=false: no --hla_qc on the command line.
# A fake nextflow stands in; a wrapper in front of it writes the two files
# main.nf's completion handler would (run-all-kir-step.sh is the KIR twin).
# Then step 08 itself with HLA_QC=true, through the fake docker:
#   5. OptiType types the read pairs T1K extracted (hla_t1k/*_hla_candidate_1
#      and _2.fq: the reads over the HLA genes and their mates, the unmapped
#      pairs and the ALT and decoy reads that match HLA sequence), not a
#      samtools cut of the BAM around HLA-A, -B and -C, and the step prints
#      their number;
#   6. a single-end BAM: T1K writes one candidate file, so OptiType does not
#      run, and paired files an earlier run left are not typed instead.
#   7. T1K's pairs hold class II reads only (no allele of HLA-A, -B or -C
#      called): OptiType, with nothing to align, does not run; a no-call.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
ST="${G}/sample1/logs/run_status.tsv"
seed_clinvar "$G"
seed_sample "$G" sample1
use_output_hook   # the reports at the end read the files their containers write
# Step 08's data where run-all.sh looks for it
for name in hla_dat gencode_genes; do
  # shellcheck source=../../scripts/lib/common.sh
  path=$(. "${REPO_ROOT}/scripts/lib/common.sh"; data_file "$name" || true)
  mkdir -p "$(dirname "$path")"; printf 'placeholder\n' > "$path"
done

mkdir -p "${CASE_WORK}/nfbin"
cat > "${CASE_WORK}/nfbin/nextflow" <<'NF'
#!/usr/bin/env bash
[ -z "${PROCESS_STATUS:-}" ] || printf 'process\tsucceeded\tcached\tfailed\tignored\tunfinished\n%b' "$PROCESS_STATUS" > process_status.tsv
[ -z "${FAILED_TASKS:-}" ] || printf 'process\tstatus\n%b' "$FAILED_TASKS" > failed_tasks.tsv
exec "${REPO_ROOT}/scripts/ci/fake-docker/nextflow" "$@"
NF
chmod +x "${CASE_WORK}/nfbin/nextflow"
NFPATH="${CASE_WORK}/nfbin:${PATH}"

# row PROCESS SUCCEEDED CACHED FAILED IGNORED UNFINISHED, as one escaped line for printf %b
row() { printf 'BAM_ANALYSIS:%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "$@"; }
MISSES=0
miss() { echo "ASSERT FAIL: $*" >&2; MISSES=$((MISSES + 1)); }
# nf_line NAME: the nextflow command line of the run NAME (the last one logged)
nf_line() { awk '/^nextflow :: / {l = $0} END {print l}' "$FAKE_DOCKER_LOG"; }
# qc NAME HLA_QC WANT FAILED_TASKS ROWS: a run that stops with ROWS as its process table; step 08 must be WANT
qc() {
  local name=$1 hq=$2 want=$3 failed=$4 rows=$5 got
  rm -rf "${G}/sample1/logs"
  run_expect 3 "$name" env PATH="$NFPATH" FAKE_NEXTFLOW_RC=3 TOOLS=hla_typing ${hq:+HLA_QC=$hq} \
    PROCESS_STATUS="$rows" FAILED_TASKS="$failed" "${SCRIPTS}/run-all.sh" sample1 male
  output_has "$name" '^  08 HLA typing \(T1K\) +runs$'
  got=$(awk -F'\t' '$2 == "08" {print $3}' "$ST")
  [ "$got" = "$want" ] || miss "${name}: step 08 is '${got}', want '${want}'"
}
T1K="$(row T1K_BUILD 1 0 0 0 0)$(row HLA_TYPING 1 0 0 0 0)"

qc qc-done true ok '' "${T1K}$(row HLA_QC_READS 1 0 0 0 0)$(row HLA_QC 1 0 0 0 0)"
nf_line | grep -q -- ' --hla_qc true' || miss "qc-done: nextflow was not given --hla_qc true"
output_has qc-done 'HLA_QC=true: OptiType checks'
qc qc-failed true failed 'BAM_ANALYSIS:HLA_QC\tfailed\n' "${T1K}$(row HLA_QC_READS 1 0 0 0 0)$(row HLA_QC 0 0 1 0 0)"
grep -q 'hla_typing (HLA_QC) failed' "${CASE_WORK}/qc-failed.out" || miss "qc-failed: the console does not name hla_typing for HLA_QC"
qc qc-reads-failed true failed 'BAM_ANALYSIS:HLA_QC_READS\tfailed\n' "${T1K}$(row HLA_QC_READS 0 0 1 0 0)$(row HLA_QC 0 0 0 0 0)"
for hq in '' false; do
  qc "qc-off${hq:+-$hq}" "$hq" ok '' "$T1K"
  if nf_line | grep -q -- '--hla_qc'; then miss "qc-off${hq:+-$hq}: nextflow was given --hla_qc without HLA_QC=true"; fi
  output_lacks "qc-off${hq:+-$hq}" 'OptiType checks'
done

# --- step 08 with HLA_QC=true: what OptiType reads ---------------------------------
seed_reference "$G"
cat > "${CASE_WORK}/hla-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
args=" ${*:2} "
opt() { sed -n "s/.* $1 \([^ ]*\) .*/\1/p" <<<"$args"; }
fq() { local i; for i in $(seq "$2"); do printf '@r%s\nACGTACGTAC\n+\nIIIIIIIIII\n' "$i"; done > "$1"; }
case "$args" in
  *" t1k-build.pl "*)
    o=$(host_path "$(opt -o)")
    mkdir -p "$o"
    printf '>HLA-A*01:01:01:01\nACGT\n' > "${o}/hla_dna_seq.fa"
    for g in HLA-A HLA-B HLA-C HLA-DRB1 HLA-DQB1 HLA-DPB1; do
      printf '>%s*01:01 chr6 29942532 29945870 +\nACGT\n' "$g"
    done > "${o}/hla_dna_coord.fa"
    exit 0 ;;
  *" run-t1k "*)
    # What run-t1k leaves: its genotypes and the reads its extractor took.
    p="$(host_path "$(opt --od)")$(opt -o)"
    if [ "${FAKE_T1K_CLASS2:-}" = 1 ]; then
      printf 'HLA-A\t0\t.\t0\t-1\t.\t0\t-1\t.\nHLA-DRB1\t2\tHLA-DRB1*15:01:01\t80\t60\tHLA-DRB1*03:01:01\t70\t60\t.\n'
    else
      printf 'HLA-A\t2\tHLA-A*01:01:01\t100\t60\tHLA-A*26:01:01\t90\t60\n'
    fi > "${p}_genotype.tsv"
    if [ "${FAKE_T1K_SINGLE:-}" = 1 ]; then fq "${p}_candidate.fq" 5
    else fq "${p}_candidate_1.fq" 7; fq "${p}_candidate_2.fq" 7; fi
    exit 0 ;;
  *" optitype run "*)
    o="$(host_path "$(opt -o)")/run1"
    mkdir -p "$o"
    printf '\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n0\tA*01:01\tA*26:01\tB*08:01\tB*07:02\tC*07:01\tC*07:02\t7\t1.0\n' \
      > "${o}/$(opt -p)_result.tsv"
    exit 0 ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/hla-hook"
T1K_DIR="${G}/sample1/hla_t1k"
OT_RES="${G}/sample1/hla_optitype/sample1_optitype_result.tsv"

: > "$FAKE_DOCKER_LOG"
FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hla-hook" HLA_QC=true run_expect 0 step08-qc "${SCRIPTS}/08-hla-typing.sh" sample1
docker_log_has '^run image=[^ ]*t1k.* run-t1k -b /genome/sample1/aligned/sample1_sorted\.bam ' \
  "step08-qc: T1K did not extract the HLA reads from the BAM itself"
docker_log_has '^run image=[^ ]*optitype.* optitype run -i /genome/sample1/hla_t1k/sample1_hla_candidate_1\.fq -i /genome/sample1/hla_t1k/sample1_hla_candidate_2\.fq ' \
  "step08-qc: OptiType did not type the read pairs T1K extracted"
if grep -q '^run image=[^ ]*samtools' "$FAKE_DOCKER_LOG"; then
  miss "step08-qc: step 08 cut OptiType's reads out of the BAM with samtools: $(grep '^run image=[^ ]*samtools' "$FAKE_DOCKER_LOG" | head -n 1)"
fi
output_has step08-qc "^HLA read pairs T1K extracted, OptiType's input: 7$"
grep -q $'^0\tA\\*01:01' "$OT_RES" 2>/dev/null || miss "step08-qc: no OptiType result in ${OT_RES}"

# A single-end BAM, with the paired files of the run above still there.
: > "$FAKE_DOCKER_LOG"
FAKE_T1K_SINGLE=1 FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hla-hook" HLA_QC=true run_expect 0 step08-qc-single "${SCRIPTS}/08-hla-typing.sh" sample1
[ -s "${T1K_DIR}/sample1_hla_candidate.fq" ] || miss "step08-qc-single: the fake run-t1k wrote no single-end candidate file"
[ ! -e "${T1K_DIR}/sample1_hla_candidate_1.fq" ] || miss "step08-qc-single: the last run's paired candidate reads were left in hla_t1k/"
if grep -q '^run image=[^ ]*optitype' "$FAKE_DOCKER_LOG"; then
  miss "step08-qc-single: OptiType ran without read pairs from this run"
fi
output_has step08-qc-single "^HLA read pairs T1K extracted, OptiType's input: 0$"
grep -q '^# OptiType not run: no read pairs' "$OT_RES" 2>/dev/null || miss "step08-qc-single: the result file does not say OptiType was not run"

# Read pairs, but T1K called no allele of HLA-A, -B or -C (class II reads only).
: > "$FAKE_DOCKER_LOG"
rm -f "$OT_RES"
FAKE_T1K_CLASS2=1 FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hla-hook" HLA_QC=true run_expect 0 step08-qc-class2 "${SCRIPTS}/08-hla-typing.sh" sample1
output_has step08-qc-class2 "^HLA read pairs T1K extracted, OptiType's input: 7$"
if grep -q '^run image=[^ ]*optitype' "$FAKE_DOCKER_LOG"; then
  miss "step08-qc-class2: OptiType ran on reads with no HLA-A, -B or -C allele called"
fi
grep -q '^# OptiType not run: no read pairs' "$OT_RES" 2>/dev/null || miss "step08-qc-class2: the result file does not say OptiType was not run"

[ "$MISSES" -eq 0 ] || fail "${MISSES} result(s) above are wrong"
echo "HLA_QC=true reaches the pipeline as --hla_qc true, step 08 counts OptiType's two processes, and OptiType types the reads T1K extracted."
