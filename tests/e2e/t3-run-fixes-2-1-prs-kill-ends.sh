#!/usr/bin/env bash
# A pgsc_calc task killed for memory ends the PRS task, which Nextflow then
# retries; the run never hangs on it.
#
# PRS runs pgsc_calc, a Nextflow pipeline of its own. Its config used to give
# FORMAT_SCOREFILES errorStrategy 'finish', and pgsc_calc's own default is
# 'finish' for any other failure. With 'finish' the nested Nextflow waits for
# its pending tasks and then for its dataflow operators. In pgsc_calc's
# ancestry branch a task that ends minutes after the failure leaves an
# operator that never ends, so the nested run never exits, the PRS task
# stays RUNNING and the outer run waits on it for ever. PRS now gives
# pgsc_calc 'terminate': a failure there kills its pending tasks at once,
# pgsc_calc exits, and PRS exits 137, so nextflow.config's retry rule runs it
# again with twice the memory. PRS is report-only, so a PRS that fails again
# is ignored and the run goes on.
#
# The real PRS module script (no -stub) runs the pgsc_calc copy setup.sh
# installs, in Docker, on a synthetic sample, score file and ancestry panel.
# A launch-folder nextflow.config that the PRS task's beforeScript copies in
# (Nextflow reads it beside pgsc_calc's own, under PRS's -c images.config)
# gives FORMAT_SCOREFILES an 8 MB container, so the kernel kills it for
# memory (exit 137), and holds EXTRACT_DATABASE for 120 s before it starts,
# so a task is still pending when the kill comes. On every attempt.
# case-timeout: 900
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
PGSC="${GENOME_DIR}/tools/pgsc_calc-${PGSC_CALC_VERSION}"
check "setup.sh installed pgsc_calc ${PGSC_CALC_VERSION}" test -s "${PGSC}/main.nf"
[ -s "${PGSC}/main.nf" ] || finish

T="${CASE_TMP}/nf"
NP="${T}/project"
rm -rf "$T"
mkdir -p "${NP}" "${T}/pgs" "${T}/panel"
# The run gives up after BOUND seconds: a hung PRS shows as a timeout.
BOUND=600
# EXTRACT_DATABASE is held this long: a pgsc_calc that waits on it runs
# longer than this, one that ends at the kill does not.
HOLD=120

# --- inputs -------------------------------------------------------------------
{
  printf '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n'
  printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tkilled\n'
  for p in 100 200 300 400 500; do printf 'chr1\t%s\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n' "$p"; done
} | gzip -c > "${T}/killed.vcf.gz"
printf '#pgs_id=PGS000001\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n1\t100\tG\tA\t1\n' \
  | gzip -c > "${T}/pgs/PGS000001.txt.gz"
# The panel EXTRACT_DATABASE unpacks: pgsc_calc checks meta.txt's version only.
( cd "${T}/panel" && touch GRCh38_x_ALL.pgen GRCh38_x_ALL.psam GRCh38_x_ALL.pvar.zst GRCh38_x.king.cutoff.out.id \
  && sed -n "s/^ *ref_format_version *= *['\"]\([^'\"]*\)['\"].*/\1/p" "${PGSC}/nextflow.config" > meta.txt \
  && tar -cf "${T}/panel.tar" GRCh38_x_ALL.pgen GRCh38_x_ALL.psam GRCh38_x_ALL.pvar.zst GRCh38_x.king.cutoff.out.id meta.txt )
check "the panel names pgsc_calc's reference format" test -s "${T}/panel/meta.txt"

# What the PRS task puts beside pgsc_calc: the kill and the pending task.
cat > "${T}/nested.config" <<NFCONF
process {
    withName: 'FORMAT_SCOREFILES' { containerOptions = '--memory 8m --memory-swap 8m' }
    withName: 'EXTRACT_DATABASE' { beforeScript = 'sleep ${HOLD}' }
}
NFCONF

cat > "${NP}/main.nf" <<NF
include { PRS } from '${REPO}/modules/local/prs/main'

workflow {
    PRS(channel.of([[id: 'killed'], file('${T}/killed.vcf.gz'), 'vcf', file('${T}/pgs')]), file('${T}/panel.tar'))
}
NF
# Every attempt fits the local executor at once (docs/lessons-learned.md, "The
# local executor can lose track of free memory"). PRS asks for 1 GB, twice that on
# its retry: pgsc_calc runs its tasks on this machine with a local executor of
# its own, which refuses a task that asks for more than the machine has, and PRS
# gives every pgsc_calc task its own memory as --max_memory.
cat > "${T}/run.config" <<NFCONF
executor.cpus = 64
executor.memory = '512.GB'
params.pgsc_calc = '${PGSC}'
process {
    withName: 'PRS' {
        memory = { 1.GB * task.attempt }
        beforeScript = "cp '${T}/nested.config' nextflow.config"
    }
}
trace.raw = true
trace.overwrite = true
trace.fields = 'task_id,name,attempt,exit,memory,realtime,status'
timeline.enabled = false
report.enabled = false
dag.enabled = false
NFCONF

# attempts TASK: "attempt:exit:GB" for each attempt of TASK, in attempt order,
# one row per task_id. TASK is matched as a plain string, not a regex: in a -v
# value gawk turns "\(" into "(", and mawk keeps it.
attempts() {
  awk -F'\t' -v task="$1" 'NR > 1 && $2 == task {
      if (!($1 in ex) || ex[$1] == "-") { ex[$1] = $4; row[$1] = $3 ":" $4 ":" int($5 / 1073741824) } }
    END { for (k in row) print row[k] }' "${T}/run.trace" 2>/dev/null | sort -t: -k1,1n | paste -sd, -
}
# longest TASK: the longest realtime of TASK's attempts, in whole seconds.
longest() {
  awk -F'\t' -v task="$1" 'NR > 1 && $2 == task && $6 ~ /^[0-9]+$/ { s = int($6 / 1000); if (s > m) m = s }
    END { print m + 0 }' "${T}/run.trace" 2>/dev/null
}

echo "+ nextflow run (PRS, FORMAT_SCOREFILES killed for memory, EXTRACT_DATABASE held ${HOLD} s), at most ${BOUND} s"
SECONDS=0
rc=0
(cd "$T" && timeout "$BOUND" nextflow run "${NP}/main.nf" -ansi-log false -c "${REPO}/nextflow.config" -c "${T}/run.config" \
  -with-trace "${T}/run.trace" --outdir "${T}/out" -work-dir "${T}/work" --max_cpus 4 --max_memory 8.GB) > "${T}/run.log" 2>&1 || rc=$?
took=$SECONDS
echo "nextflow exited ${rc} after ${took} s"
tail -n 30 "${T}/run.log"
cat "${T}/run.trace" 2>/dev/null

check "the run ends before ${BOUND} s (a hung pgsc_calc runs into the timeout)" test "$rc" -ne 124
check_eq "the run goes on past the failed PRS: nextflow exits 0" "$rc" 0
check_eq "PRS exits 137 on both attempts, with twice the memory on the retry" "$(attempts 'PRS (killed)')" "1:137:1,2:137:2"
check "each PRS attempt ends within ${HOLD} s: pgsc_calc did not wait for its pending task" \
  test "$(longest 'PRS (killed)')" -lt "$HOLD"
# The nested run of the first attempt: the planted kill, then terminate.
NLOG=$(grep -rl --include=pgsc_calc.nextflow.log 'FORMAT_SCOREFILES' "${T}/work" 2>/dev/null | head -n 1)
check "the kill is the kernel's: FORMAT_SCOREFILES exits 137 in pgsc_calc" \
  grep -qF 'FORMAT_SCOREFILES (1)` terminated with an error exit status (137)' "${NLOG:-/dev/null}"
check "pgsc_calc ends at the kill (Session aborted), not after its pending tasks" \
  grep -q 'Session aborted' "${NLOG:-/dev/null}"
check "pgsc_calc never waits for its pending tasks (no 'Execution cancelled')" \
  bash -c '! grep -q "Execution cancelled" "$1"' _ "${NLOG:-/dev/null}"
check "PRS names the memory kill" \
  grep -rqF 'a pgsc_calc task was killed for memory' --include=.command.err "${T}/work"

rm -rf "$T"
finish
