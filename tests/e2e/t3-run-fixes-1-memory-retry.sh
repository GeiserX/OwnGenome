#!/usr/bin/env bash
# LOCITYPER and CPSR under the pipeline's retry rule, on the modules' real
# scripts (no -stub, no Docker), with stand-ins for locityper, og and cpsr on
# PATH. Nextflow retries a task that exits 137 or 143 once, with twice the
# memory, capped by --max_memory.
#
# LOCITYPER used to catch Locityper's exit code and write its table with it,
# so a Locityper killed for memory left a green task and no retry. Now a kill
# (137 or 143) on the first attempt ends the task with that code, and
# nextflow.config's rule for the report-only tools retries attempt 1 only; the
# attempt after it writes the table. Samples: "once" is killed on its first
# attempt only, "always" on every attempt, "term" is stopped by SIGTERM (143)
# on its first attempt only, "ok" never.
#
# CPSR used to ask for process_medium's 8 GB, and on a whole genome it is
# killed under that. It asks for 16 GB now, still capped by --max_memory:
# sample "once" is killed on its first attempt (16 GB, then 32 GB), and a
# second run with --max_memory 8.GB must give it 8 GB.
#
# The stand-ins kill themselves with SIGKILL, as the kernel does on out of
# memory, except "term". Rows of the trace are keyed by task_id: Nextflow can write a
# finished task twice (docs/lessons-learned.md, "The trace file can list a
# finished task twice").
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

T="${CASE_TMP}/nf"
NP="${T}/project"
rm -rf "$T"
mkdir -p "${NP}/bin" "${T}/marks" "${T}/db/db" "${T}/pcgr" "${T}/vep115"
cp "${REPO}/bin/cpsr_sample_id" "${NP}/bin/"

# --- stand-ins ----------------------------------------------------------------
# locityper: preproc passes; genotype is killed as the sample (the BAM's name) says.
cat > "${NP}/bin/locityper" <<STANDIN
#!/bin/bash
cmd=\$1; shift
out="" bam=""
while [ \$# -gt 0 ]; do case "\$1" in -o) out=\$2; shift ;; -a) bam=\$2; shift ;; esac; shift; done
s=\$(basename "\$bam" .bam)
if [ "\$cmd" = genotype ]; then
    case "\$s" in
        always) kill -9 \$\$ ;;
        once) [ -e "${T}/marks/locityper-once" ] || { touch "${T}/marks/locityper-once"; kill -9 \$\$; } ;;
        term) [ -e "${T}/marks/locityper-term" ] || { touch "${T}/marks/locityper-term"; kill -15 \$\$; } ;;
    esac
fi
mkdir -p "\$out"
STANDIN
# og locityper report: one line with the exit code it was given.
cat > "${NP}/bin/og" <<'STANDIN'
#!/bin/bash
[ "$1 $2" = "locityper report" ] || { echo "unexpected og $*" >&2; exit 2; }
shift 2; rc="" out=""
while [ $# -gt 0 ]; do case "$1" in --exit-code) rc=$2; shift ;; --out) out=$2; shift ;; esac; shift; done
printf 'exit_code\t%s\n' "$rc" > "$out"
STANDIN
# cpsr: killed on sample once's first attempt, else writes the two files CPSR writes.
cat > "${NP}/bin/cpsr" <<STANDIN
#!/bin/bash
id="" dir=""
while [ \$# -gt 0 ]; do case "\$1" in --sample_id) id=\$2; shift ;; --output_dir) dir=\$2; shift ;; esac; shift; done
case "\$id" in
    once) [ -e "${T}/marks/cpsr-once" ] || { touch "${T}/marks/cpsr-once"; kill -9 \$\$; } ;;
esac
touch "\$dir/\$id.cpsr.grch38.html" "\$dir/\$id.cpsr.grch38.classification.tsv.gz"
STANDIN
chmod +x "${NP}/bin/locityper" "${NP}/bin/og" "${NP}/bin/cpsr"

# --- inputs -------------------------------------------------------------------
printf '>chr1\nACGT\n' > "${T}/ref.fa"
printf 'chr1\t4\t6\t4\t5\n' > "${T}/ref.fa.fai"
printf 'chr6\t1\t2\tHLA-DQB1\n' > "${T}/db/loci.bed"
: > "${T}/db/counts.jf"
echo "fai=$(sha256sum "${T}/ref.fa.fai" | cut -d' ' -f1) fasta=$(sha256sum "${T}/ref.fa" | cut -d' ' -f1) " > "${T}/db/DONE"
for s in once always ok term; do : > "${T}/${s}.bam"; : > "${T}/${s}.bam.bai"; done
for s in once ok; do : > "${T}/${s}.vcf.gz"; : > "${T}/${s}.vcf.gz.tbi"; done

cat > "${NP}/main.nf" <<NF
include { LOCITYPER } from '${REPO}/modules/local/locityper/main'
include { CPSR } from '${REPO}/modules/local/cpsr/main'

params.samples = 'once,always,ok,term'
params.cpsr_samples = 'once'

workflow {
    def loc = params.samples.toString().tokenize(',')
    def cps = params.cpsr_samples.toString().tokenize(',')
    LOCITYPER(channel.fromList(loc).map { s -> [[id: s], file("${T}/\${s}.bam"), file("${T}/\${s}.bam.bai")] },
              file('${T}/ref.fa'), file('${T}/ref.fa.fai'), file('${T}/db'))
    CPSR(channel.fromList(cps).map { s -> [[id: s], file("${T}/\${s}.vcf.gz"), file("${T}/\${s}.vcf.gz.tbi")] },
         file('${T}/pcgr'), file('${T}/vep115'))
}
NF
# Every attempt of the run fits the local executor at once, so no task waits
# for room (docs/lessons-learned.md, "The local executor can lose track of
# free memory").
cat > "${T}/run.config" <<'NFCONF'
executor.cpus = 64
executor.memory = '512.GB'
trace.raw = true
trace.overwrite = true
trace.fields = 'task_id,name,attempt,exit,memory'
timeline.enabled = false
report.enabled = false
dag.enabled = false
NFCONF

# nf_go NAME OPTIONS...: one run with the checkout's nextflow.config (labels,
# withName requests, errorStrategy, resourceLimits); its trace is NAME.trace.
nf_go() {
  local name=$1 rc=0
  shift
  echo "+ nextflow run (${name}) $*"
  (cd "$T" && nextflow run "${NP}/main.nf" -ansi-log false -c "${REPO}/nextflow.config" -c "${T}/run.config" \
    -with-trace "${T}/${name}.trace" --outdir "${T}/${name}-out" -work-dir "${T}/${name}-work" "$@") > "${T}/${name}.log" 2>&1 || rc=$?
  echo "nextflow exited ${rc}"
  [ "$rc" -eq 0 ] || tail -n 40 "${T}/${name}.log"
  cat "${T}/${name}.trace" 2>/dev/null
}
# attempts NAME TASK: "attempt:exit:GB" for each attempt of the task named TASK,
# in attempt order, one row per task_id. TASK is matched as a plain string, not
# a regex: in a -v value gawk turns "\(" into "(", and mawk keeps it.
attempts() {
  awk -F'\t' -v task="$2" 'NR > 1 && $2 == task {
      if (!($1 in ex) || ex[$1] == "-") { ex[$1] = $4; row[$1] = $3 ":" $4 ":" int($5 / 1073741824) } }
    END { for (k in row) print row[k] }' "${T}/$1.trace" 2>/dev/null | sort -t: -k1,1n | paste -sd, -
}

nf_go big --max_cpus 8 --max_memory 64.GB
check_eq "LOCITYPER killed once: attempt 1 exits 137 at 32 GB, attempt 2 passes at 64 GB" \
  "$(attempts big 'LOCITYPER (once)')" "1:137:32,2:0:64"
check_eq "LOCITYPER killed on every attempt: the retry at 64 GB writes the table" \
  "$(attempts big 'LOCITYPER (always)')" "1:137:32,2:0:64"
check_eq "LOCITYPER not killed: one attempt at 32 GB" "$(attempts big 'LOCITYPER (ok)')" "1:0:32"
check_eq "LOCITYPER stopped by SIGTERM once: attempt 1 exits 143 at 32 GB, attempt 2 passes at 64 GB" \
  "$(attempts big 'LOCITYPER (term)')" "1:143:32,2:0:64"
check_eq "the term table is the retry's: Locityper ended well" \
  "$(cut -f2 "${T}/big-out/term/locityper/term_locityper.tsv" 2>/dev/null)" 0
check_eq "the once table is the retry's: Locityper ended well" \
  "$(cut -f2 "${T}/big-out/once/locityper/once_locityper.tsv" 2>/dev/null)" 0
check_eq "the always table holds the last kill's exit code" \
  "$(cut -f2 "${T}/big-out/always/locityper/always_locityper.tsv" 2>/dev/null)" 137
check "the killed attempt says why it ends" \
  grep -rqF "Locityper was killed (exit 137) on attempt 1: retrying with more memory" --include=.command.err "${T}/big-work"
check_eq "CPSR asks for 16 GB, and 32 GB on the retry after a kill" "$(attempts big 'CPSR (once)')" "1:137:16,2:0:32"
check "CPSR's report is published" test -e "${T}/big-out/once/cpsr/once.cpsr.grch38.html"

nf_go capped --max_cpus 8 --max_memory 8.GB --samples ok --cpsr_samples ok
check_eq "--max_memory 8.GB still caps CPSR's 16 GB" "$(attempts capped 'CPSR (ok)')" "1:0:8"

rm -rf "$T"
finish
