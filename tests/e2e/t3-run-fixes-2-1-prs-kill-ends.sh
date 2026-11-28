#!/usr/bin/env bash
# A pgsc_calc task killed for memory ends the PRS task, which Nextflow then
# retries; the run never hangs on it.
#
# PRS runs pgsc_calc, a Nextflow pipeline of its own. Its config used to give
# FORMAT_SCOREFILES errorStrategy 'finish', and pgsc_calc's own default is
# 'finish' for any other failure; MATCH_VARIANTS and PLINK2_SCORE set
# 'finish' themselves. With 'finish' the nested Nextflow waits for its
# pending tasks and then for its dataflow operators. In pgsc_calc's ancestry
# branch a task that ends minutes after the failure leaves an operator that
# never ends, so the nested run never exits, the PRS task stays RUNNING and
# the outer run waits on it for ever. PRS now gives every pgsc_calc task
# 'terminate' (a withName line for the processes that set errorStrategy
# themselves, which outrank a generic line): a failure there kills its
# pending tasks at once, pgsc_calc exits, and PRS exits 137, so
# nextflow.config's retry rule runs it again with twice the memory. PRS is
# report-only, so a PRS that fails again is ignored and the run goes on.
#
# The real PRS module script (no -stub) runs the pgsc_calc copy setup.sh
# installs, in Docker, on a synthetic sample, score file and ancestry panel.
# A launch-folder nextflow.config that the PRS task's beforeScript copies in
# (Nextflow reads it beside pgsc_calc's own, under PRS's -c images.config)
# gives one task an 8 MB container, so the kernel kills it for memory (exit
# 137), and holds EXTRACT_DATABASE for 120 s before it starts, so a task is
# still pending when the kill comes. On every attempt. Two runs: the kill in
# FORMAT_SCOREFILES, then in MATCH_VARIANTS, which sets 'finish' itself (the
# score file is one pgscatalog-format reads, so pgsc_calc gets that far).
# case-timeout: 1200
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
PGSC="${GENOME_DIR}/tools/pgsc_calc-${PGSC_CALC_VERSION}"
check "setup.sh installed pgsc_calc ${PGSC_CALC_VERSION}" test -s "${PGSC}/main.nf"
[ -s "${PGSC}/main.nf" ] || finish

# Each run gives up after BOUND seconds: a hung PRS shows as a timeout.
BOUND=420
# EXTRACT_DATABASE is held this long: a pgsc_calc that waits on it runs
# longer than this, one that ends at the kill does not.
HOLD=120

# run_kill KILLED: run PRS with the task KILLED killed for memory, and check
# that every part of it ends at the kill.
run_kill() {
  local killed=$1
  local T="${CASE_TMP}/${killed}"
  local NP="${T}/project"
  rm -rf "$T"
  mkdir -p "${NP}" "${T}/pgs" "${T}/panel"

  # --- inputs -----------------------------------------------------------------
  {
    printf '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n'
    printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tkilled\n'
    for p in 100 200 300 400 500; do printf 'chr1\t%s\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n' "$p"; done
  } | gzip -c > "${T}/killed.vcf.gz"
  # The header PRS_PREPARE writes (og prs format), which pgscatalog-format
  # reads, with the score's five positions.
  {
    printf '#pgs_id=PGS000001\n#pgs_name=PGS000001\n#trait_reported=planted trait\n#genome_build=GRCh38\n'
    printf 'chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n'
    for p in 100 200 300 400 500; do printf '1\t%s\tG\tA\t1\n' "$p"; done
  } | gzip -c > "${T}/pgs/PGS000001.txt.gz"
  # The panel EXTRACT_DATABASE unpacks: pgsc_calc checks meta.txt's version only.
  ( cd "${T}/panel" && touch GRCh38_x_ALL.pgen GRCh38_x_ALL.psam GRCh38_x_ALL.pvar.zst GRCh38_x.king.cutoff.out.id \
    && sed -n "s/^ *ref_format_version *= *['\"]\([^'\"]*\)['\"].*/\1/p" "${PGSC}/nextflow.config" > meta.txt \
    && tar -cf "${T}/panel.tar" GRCh38_x_ALL.pgen GRCh38_x_ALL.psam GRCh38_x_ALL.pvar.zst GRCh38_x.king.cutoff.out.id meta.txt )
  check "${killed}: the panel names pgsc_calc's reference format" test -s "${T}/panel/meta.txt"

  # What the PRS task puts beside pgsc_calc: the kill and the pending task.
  # The hold ends at the TERM Nextflow sends when it stops the task: a plain
  # sleep defers the trap, and Nextflow's own TERM handler does not exit, so
  # the held task would start its container minutes after the case ended.
  cat > "${T}/nested.config" <<NFCONF
process {
    withName: '${killed}' { containerOptions = '--memory 8m --memory-swap 8m' }
    withName: 'EXTRACT_DATABASE' { beforeScript = 'sleep ${HOLD} & trap "kill \$! 2>/dev/null; exit 143" TERM; wait \$!' }
}
NFCONF

  cat > "${NP}/main.nf" <<NF
include { PRS } from '${REPO}/modules/local/prs/main'

workflow {
    PRS(channel.of([[id: 'killed'], file('${T}/killed.vcf.gz'), 'vcf', file('${T}/pgs')]), file('${T}/panel.tar'))
}
NF
  # Every attempt fits the local executor at once (docs/lessons-learned.md,
  # "The local executor can lose track of free memory"). PRS asks for 1 GB,
  # twice that on its retry: pgsc_calc runs its tasks on this machine with a
  # local executor of its own, which refuses a task that asks for more than
  # the machine has, and PRS gives every pgsc_calc task its own memory as
  # --max_memory.
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

  echo "+ nextflow run (PRS, ${killed} killed for memory, EXTRACT_DATABASE held ${HOLD} s), at most ${BOUND} s"
  SECONDS=0
  local rc=0
  (cd "$T" && timeout "$BOUND" nextflow run "${NP}/main.nf" -ansi-log false -c "${REPO}/nextflow.config" -c "${T}/run.config" \
    -with-trace "${T}/run.trace" --outdir "${T}/out" -work-dir "${T}/work" --max_cpus 4 --max_memory 8.GB) > "${T}/run.log" 2>&1 || rc=$?
  local took=$SECONDS
  echo "nextflow exited ${rc} after ${took} s"
  tail -n 30 "${T}/run.log"
  cat "${T}/run.trace" 2>/dev/null

  check "${killed}: the run ends before ${BOUND} s (a hung pgsc_calc runs into the timeout)" test "$rc" -ne 124
  check_eq "${killed}: the run goes on past the failed PRS: nextflow exits 0" "$rc" 0
  check_eq "${killed}: PRS exits 137 on both attempts, with twice the memory on the retry" \
    "$(attempts "$T" 'PRS (killed)')" "1:137:1,2:137:2"
  # Only an attempt that ended has a realtime; a hung one has '-'. So both
  # attempts must have ended, each within HOLD.
  check_eq "${killed}: each of the 2 PRS attempts ended within ${HOLD} s: pgsc_calc did not wait for its pending task" \
    "$(ended_within "$T" 'PRS (killed)' "$HOLD")" "2 of 2"
  # The nested run of the first attempt: the planted kill, then terminate.
  local nlog
  nlog=$(grep -rlF --include=pgsc_calc.nextflow.log "${killed} " "${T}/work" 2>/dev/null | head -n 1)
  check "${killed}: the kill is the kernel's: ${killed} exits 137 in pgsc_calc" \
    grep -qE "${killed} \\([^)]*\\)\` terminated with an error exit status \\(137\\)" "${nlog:-/dev/null}"
  # 'Session aborted' comes also when the outer run's timeout stops a hung
  # pgsc_calc, so it must come within HOLD of the kill, by the log's clock.
  check_eq "${killed}: pgsc_calc ends at the kill: 'Session aborted' within ${HOLD} s of the 137" \
    "$(aborted_after_kill "${nlog:-/dev/null}" "$killed" "$HOLD")" "aborted at the kill"
  # shellcheck disable=SC2016  # bash -c expands its own "$1"
  check "${killed}: pgsc_calc never waits for its pending tasks (no 'Execution cancelled')" \
    bash -c 'test -s "$1" && ! grep -q "Execution cancelled" "$1"' _ "${nlog:-/dev/null}"
  check "${killed}: PRS names the memory kill" \
    grep -rqF 'a pgsc_calc task was killed for memory' --include=.command.err "${T}/work"

  # Nothing of the run is left running: no process works inside the case folder.
  sleep 3
  local left=0 p
  for p in /proc/[0-9]*; do
    case "$(readlink "${p}/cwd" 2>/dev/null)" in "${T}"/*) left=$((left + 1)) ;; esac
  done
  check_eq "${killed}: no task of the run is still running after it ends" "$left" 0
  rm -rf "$T"
}

# attempts DIR TASK: "attempt:exit:GB" for each attempt of TASK in DIR's
# trace, in attempt order, one row per task_id. TASK is matched as a plain
# string, not a regex: in a -v value gawk turns "\(" into "(", and mawk
# keeps it.
attempts() {
  awk -F'\t' -v task="$2" 'NR > 1 && $2 == task {
      if (!($1 in ex) || ex[$1] == "-") { ex[$1] = $4; row[$1] = $3 ":" $4 ":" int($5 / 1073741824) } }
    END { for (k in row) print row[k] }' "${1}/run.trace" 2>/dev/null | sort -t: -k1,1n | paste -sd, -
}
# ended_within DIR TASK S: "N of M": of TASK's M attempts, the N that ended
# (a realtime in the trace; a running or killed-by-timeout one has none)
# within S seconds. One row per task_id, its last.
ended_within() {
  awk -F'\t' -v task="$2" -v s="$3" 'NR > 1 && $2 == task { rt[$1] = $6 }
    END { for (k in rt) { m++; if (rt[k] ~ /^[0-9]+$/ && rt[k] / 1000 < s) n++ }
          print n + 0 " of " m + 0 }' "${1}/run.trace" 2>/dev/null
}
# aborted_after_kill LOG TASK S: "aborted at the kill" when LOG's first
# 'Session aborted' line comes after TASK's exit 137 line and within S
# seconds of it, by the clock of Nextflow's log lines (a line without one,
# such as the 137 line of an error report, has the time of the line before
# it); else what it saw.
aborted_after_kill() {
  awk -v task="$2" -v s="$3" '
    BEGIN { kill = -1; abort = -1; now = -1 }
    /^[A-Z][a-z][a-z]-[0-9][0-9] [0-9][0-9]:[0-9][0-9]:[0-9][0-9]/ {
      now = substr($0, 8, 2) * 3600 + substr($0, 11, 2) * 60 + substr($0, 14, 2) }
    kill < 0 && index($0, task " (") && index($0, "terminated with an error exit status (137)") { kill = now }
    kill >= 0 && /Session aborted/ { abort = now - kill; if (abort < 0) abort += 86400; exit }
    END { if (kill < 0) print "no 137 line"
          else if (abort < 0) print "no Session aborted after the 137"
          else if (abort > s) print "Session aborted " abort " s after the 137"
          else print "aborted at the kill" }' "$1" 2>/dev/null
}

run_kill FORMAT_SCOREFILES
run_kill MATCH_VARIANTS
finish
