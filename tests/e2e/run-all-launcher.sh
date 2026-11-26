#!/usr/bin/env bash
# scripts/run-all.sh, the launcher, from the fixture's FASTQ pair as sample
# HG002P, in a GENOME_DIR of its own that holds only what case
# nextflow-from-fastq-2-nextflow gives the pipeline directly: the reference
# without its minimap2 index, ClinVar, the HLA data and the PRS score files.
# TOOLS is that case's list without vcfanno (run-all.sh runs it only with VEP)
# and html_report (run-all.sh renders the full report with step 24 instead).
#   1. It publishes the same files as the direct `nextflow run`, and the seven
#      items of scripts/ci/parity-diff.sh agree with the bash steps' outputs.
#   2. A second run starts from the BAM, VCF and gVCF the first one published,
#      not from the FASTQ its samplesheet named: no FASTP, minimap2 or
#      DeepVariant task runs. A third run finds every task in Nextflow's cache
#      (-resume): no task runs again, and it ends in minutes.
#   3. The reports, logs/run_status.tsv and OwnHealth's facts file are written.
#   4. It runs as the guide documents it, without THREADS or --max_memory:
#      the launcher passes the runner's CPU count and RAM, and Nextflow starts
#      the 8-CPU and 32 GB tasks capped to them instead of refusing them.
#   5. Hard-link publishing. The tasks' files are root's (the docker profile),
#      and this runner's user may link them only while fs.protected_hardlinks
#      is 0. The first run sets it to 0: run-all.sh then publishes by hard link,
#      the BAM is one file with two names (work/ and aligned/), and no published
#      file is a symbolic link the direct (copy) run lacks. The second run has
#      the runner's own setting back (1): run-all.sh passes no link mode.
#   6. A container engine with less memory and fewer CPUs than the host, as
#      Docker Desktop and colima have (a docker shim that reports 8 GB and
#      2 CPUs, and a stand-in nextflow that records its arguments): run-all.sh
#      passes the engine's as --max_memory and --max_cpus, and the pool
#      Nextflow fills, executor.config, is the engine's. Each limit also
#      alone; THREADS, an explicit --max_memory and a -c of the user's win,
#      and a cap the user sets above the engine's raises the pool to it.
# Each comparison gets a negative control in the same run: a planted missing
# file must be reported, the first run's trace must show the alignment and
# calling tasks, and the first run's trace (tasks COMPLETED) must fail the
# all-cached check.
# SKIP_VALIDATION=true: the fixture's reference is a slice that validate-setup.sh
# rejects (case 10); the fake-docker cases check that run-all.sh calls it.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
P=HG002P
G2="${E2E_WORK}/genome-runall"
DIRECT="${GENOME_DIR}/nf-fastq/${P}"
check "case nextflow-from-fastq-2-nextflow published its outputs" test -d "$DIRECT"
rm -rf "$G2"
mkdir -p "${G2}/reference" "${G2}/clinvar" "${G2}/prs_scores" "${G2}/${P}/fastq"
for e in fasta fasta.fai dict; do ln -f "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.${e}" "${G2}/reference/"; done
ln -f "${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz" "${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz.tbi" "${G2}/clinvar/"
for r in R1 R2; do ln -f "${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_${r}.fastq.gz" "${G2}/${P}/fastq/${P}_${r}.fastq.gz"; done
ln -f "${E2E_WORK}/parity-pgs/"*.txt.gz "${G2}/prs_scores/"
for d in hla_dat gencode_genes; do
  from=$(GENOME_DIR="$GENOME_DIR" bash -c '. "$1/scripts/lib/common.sh" && data_file "$2"' _ "$REPO" "$d")
  to=$(GENOME_DIR="$G2" bash -c '. "$1/scripts/lib/common.sh" && data_file "$2"; true' _ "$REPO" "$d")
  mkdir -p "$(dirname "$to")" && ln -f "$from" "$to"
done
INTERVALS=$(awk '{printf "%s%s:%d-%d", (NR > 1 ? " " : ""), $1, $2 + 1, $3}' "${FIXTURE_DIR}/regions.bed")
TOOLS=pharmcat,cpic,roh,prs,mito_haplogroup,telomere_hunter,mosdepth,mito_variants,multiqc,clinvar,hla_typing

launch() {  # launch NAME: run-all.sh, output to the case log and logs/<case>.<NAME>.log; exit code in RC, seconds in SECS
  local t0=$SECONDS
  echo "+ scripts/run-all.sh ${P} male --sex_check warn, THREADS unset (${1})"
  env -u THREADS GENOME_DIR="$G2" SKIP_VALIDATION=true TOOLS="$TOOLS" INTERVALS="$INTERVALS" \
    "${REPO}/scripts/run-all.sh" "$P" male --sex_check warn 2>&1 | tee "${E2E_WORK}/logs/${CASE_NAME}.${1}.log"
  RC=${PIPESTATUS[0]} SECS=$((SECONDS - t0))
  echo "+ exit ${RC} after ${SECS}s"
}
latest_trace() { find "${G2}/pipeline_info" -name 'trace_*.txt' 2>/dev/null | LC_ALL=C sort | awk 'END {print}'; }
# status_count TRACE STATUS: tasks with that status
status_count() {
  awk -F'\t' -v s="$2" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $c["status"] == s {k++} END {print k + 0}' "$1" 2>/dev/null || echo 0
}
all_cached() { [ "$(status_count "$1" CACHED)" -ge 10 ] && [ "$(status_count "$1" COMPLETED)" -eq 0 ] && [ "$(status_count "$1" FAILED)" -eq 0 ]; }
not_all_cached() { ! all_cached "$1"; }
# upstream_tasks TRACE: the FASTP, minimap2 and DeepVariant tasks the run started (not cached)
upstream_tasks() {
  awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    $c["status"] != "CACHED" && $c["name"] ~ /(^|:)(FASTP|MINIMAP2_INDEX|ALIGN_MINIMAP2|DEEPVARIANT)( |$)/ {k++} END {print k + 0}' "$1" 2>/dev/null || echo 0
}
# published DIR: the files under DIR, relative, without the launcher's own (inputs, launch directory,
# logs, manifest, bash reports, the run-to-run delta and its history copies, OwnHealth's facts
# file) and the two outputs
# the direct case makes and run-all.sh does not. The PRS summary's dated copy is named by the
# time of its run, so its stamp is masked.
published() {
  (cd "$1" && find . -type f | sed -E 's|^\./||; s/_prs_summary\.[0-9]{8}-[0-9]{4}(-[0-9]+)?\.tsv$/_prs_summary.TIME.tsv/' | LC_ALL=C sort) \
    | grep -vE "^(fastq|nextflow|logs|vep)/|^run_manifest\.tsv$|^summary\.json$|^${P}_report\.(html|txt)$|^${P}_delta\.(html|txt|json)$|^history/|^facts/${P}_genome_facts\.json$"
}

# --- 6. a container engine smaller than the host --------------------------------------------
# Its own GENOME_DIR with an empty VCF: the stand-in nextflow runs nothing and exits 1, so
# run-all.sh then writes the reports of nothing. The negative control is an engine as big as
# the host, as on native Linux: the host's caps, and no executor.config passed or left behind.
ENG="${CASE_TMP}/engine" G3="${E2E_WORK}/genome-engine" E=ENGINE1
rm -rf "$ENG" "$G3"; mkdir -p "${ENG}/bin" "${G3}/reference" "${G3}/${E}/vcf"
: > "${G3}/${E}/vcf/${E}.vcf.gz"; : > "${G3}/${E}/vcf/${E}.vcf.gz.tbi"
cat > "${ENG}/bin/docker" <<'SH'
#!/usr/bin/env bash
# docker info prints ENGINE_INFO; anything else goes to the docker on PATH without this folder. The
# folder leaves PATH too: the e2e job's own shim (tests/e2e/bin/docker) skips only its own folder, so
# it would find this one again, and the two would exec each other forever.
if [ "${1:-}" = info ]; then echo "$ENGINE_INFO"; exit 0; fi
self=$(cd "$(dirname "$0")" && pwd) rest=""
IFS=: read -r -a dirs <<< "$PATH"
for d in "${dirs[@]}"; do [ "$d" = "$self" ] || rest+="${rest:+:}${d}"; done
PATH=$rest exec docker "$@"
SH
cat > "${ENG}/bin/nextflow" <<'SH'
#!/usr/bin/env bash
# Records its arguments on one line and a copy of the executor.config it was given; runs nothing.
echo " $* " > "${ENGINE_REC}.args"
prev=""; for a in "$@"; do if [ "$prev" = -c ] && [ "${a##*/}" = executor.config ]; then cp "$a" "${ENGINE_REC}.config"; fi; prev=$a; done
exit 1
SH
chmod +x "${ENG}/bin/docker" "${ENG}/bin/nextflow"
# The host's CPUs and RAM in whole GB, as run-all.sh reads them; the small engine has 8 GB (less on a
# host with 8 GB or less) and 2 CPUs (1 on a host with 2 or fewer).
HC=$(getconf _NPROCESSORS_ONLN)
HM=$(awk '/^MemTotal:/ {print int($2 / 1048576)}' /proc/meminfo 2>/dev/null || echo $(( $(sysctl -n hw.memsize) / 1073741824 )))
EG=$(( HM > 8 ? 8 : HM / 2 )) ECPU=$(( HC > 2 ? 2 : 1 ))
engine_run() {  # engine_run NAME 'MEMTOTAL NCPU' [OPTION...]: run-all.sh with docker info answering that, THREADS=$ET when
  local n=$1 info=$2; shift 2  # ET is set; its log is logs/<case>.engine-NAME.log
  rm -f "${ENG}/${n}".*
  echo "+ ${ET:+THREADS=${ET} }scripts/run-all.sh ${E} female $*, docker info reporting '${info}' (${n})"
  env -u THREADS ${ET:+"THREADS=${ET}"} PATH="${ENG}/bin:${PATH}" ENGINE_INFO="$info" ENGINE_REC="${ENG}/${n}" GENOME_DIR="$G3" \
    SKIP_VALIDATION=true TOOLS=roh EH_CATALOG="${G3}/none.json" "${REPO}/scripts/run-all.sh" "$E" female "$@" \
    > "${E2E_WORK}/logs/${CASE_NAME}.engine-${n}.log" 2>&1
  echo "+ exit $? (the stand-in nextflow exits 1)"
}
EXC="${G3}/${E}/nextflow/executor.config"
engine_run small "$((EG * 1073741824)) ${ECPU}"
SMALL=$(cat "${ENG}/small.args" 2>/dev/null)
echo "+ nextflow${SMALL}"
check "the stand-in nextflow ran (small engine)" test -n "$SMALL"
check "nextflow got --max_memory with the engine's ${EG} GB, not the host's ${HM} GB" has " --max_memory ${EG}\.GB " "$SMALL"
check "nextflow got --max_cpus with the engine's ${ECPU} CPUs, not the host's ${HC}" has " --max_cpus ${ECPU} " "$SMALL"
check "nextflow got -c ${EXC}" has " -c ${EXC} " "$SMALL"
check "executor.config sets the pool's memory to ${EG} GB" has "^ *memory = '${EG} GB'\$" "$(cat "${ENG}/small.config" 2>/dev/null)"
check "executor.config sets the pool's CPUs to ${ECPU}" has "^ *cpus += ${ECPU}\$" "$(cat "${ENG}/small.config" 2>/dev/null)"
check "run-all.sh says it books the engine's size" has '^\[Resources\] the container engine has less than this host' \
  "$(cat "${E2E_WORK}/logs/${CASE_NAME}.engine-small.log" 2>/dev/null)"
engine_run big "$(( (HM + 4) * 1073741824 )) $((HC + 4))"
BIG=$(cat "${ENG}/big.args" 2>/dev/null)
echo "+ nextflow${BIG}"
check "negative control: an engine as big as the host gives the host's --max_memory ${HM}.GB" has " --max_memory ${HM}\.GB " "$BIG"
check "negative control: and the host's --max_cpus ${HC}" has " --max_cpus ${HC} " "$BIG"
check "negative control: and no executor.config" lacks "executor\.config" "$BIG"
check "negative control: the small engine's executor.config is not left behind" test ! -e "$EXC"
# Each limit alone, and what the user passes wins: THREADS and an explicit --max_memory stay, and a
# -c of the user's comes after the engine's executor.config, so each setting in its executor block wins.
UC="${G3}/user.config"
ET=3 engine_run memory "$((EG * 1073741824)) $((HC + 4))" -c "$UC"
MEMO=$(cat "${ENG}/memory.args" 2>/dev/null) MEMC=$(cat "${ENG}/memory.config" 2>/dev/null)
echo "+ nextflow${MEMO}"
check "less memory alone: --max_memory ${EG}.GB" has " --max_memory ${EG}\.GB " "$MEMO"
check "less memory alone: THREADS=3 stays --max_cpus 3" has " --max_cpus 3 " "$MEMO"
check "less memory alone: the pool keeps the host's ${HC} CPUs" has "^ *cpus += ${HC}\$" "$MEMC"
check "less memory alone: the pool has ${EG} GB" has "^ *memory = '${EG} GB'\$" "$MEMC"
check "the user's -c comes after the engine's executor.config" has " -c ${EXC} (.* )?-c ${UC} " "$MEMO"
engine_run cpus "$(( (HM + 4) * 1073741824 )) ${ECPU}" --max_memory 6.GB
CPUO=$(cat "${ENG}/cpus.args" 2>/dev/null) CPUC=$(cat "${ENG}/cpus.config" 2>/dev/null)
echo "+ nextflow${CPUO}"
check "fewer CPUs alone: --max_cpus ${ECPU}" has " --max_cpus ${ECPU} " "$CPUO"
check_eq "fewer CPUs alone: the user's --max_memory 6.GB is the only one" "$(grep -o ' --max_memory [^ ]*' <<<"$CPUO" | tr -d '\n')" " --max_memory 6.GB"
check "fewer CPUs alone: the pool has ${ECPU} CPUs" has "^ *cpus += ${ECPU}\$" "$CPUC"
check "fewer CPUs alone: the pool keeps the host's ${HM} GB" has "^ *memory = '${HM} GB'\$" "$CPUC"
# A cap the user sets above the engine's figure: the pool must hold it, or Nextflow refuses each task
# that asks for more than the pool (Process requirement exceeds available CPUs/memory) and the run stops.
AT=$((ECPU + 2)) AM=$((EG + 4))
ET=$AT engine_run above "$((EG * 1073741824)) ${ECPU}" --max_memory "${AM}.GB"
ABO=$(cat "${ENG}/above.args" 2>/dev/null) ABC=$(cat "${ENG}/above.config" 2>/dev/null)
echo "+ nextflow${ABO}"
check "THREADS=${AT} above the engine's ${ECPU} CPUs: --max_cpus ${AT}" has " --max_cpus ${AT} " "$ABO"
check "THREADS=${AT} above the engine's ${ECPU} CPUs: the pool has ${AT} CPUs, not ${ECPU}" has "^ *cpus += ${AT}\$" "$ABC"
check "--max_memory ${AM}.GB above the engine's ${EG} GB: the pool has ${AM} GB" has "^ *memory = '${AM} GB'\$" "$ABC"
check "run-all.sh says the pool follows the user's cap" has '^\[Resources\] the pool follows a task cap you set' \
  "$(cat "${E2E_WORK}/logs/${CASE_NAME}.engine-above.log" 2>/dev/null)"
engine_run above2 "$((EG * 1073741824)) ${ECPU}" "--max_cpus=$((AT + 1))" --max_memory 9000.MB
AB2=$(cat "${ENG}/above2.config" 2>/dev/null)
check "--max_cpus=$((AT + 1)) above the engine's: the pool has $((AT + 1)) CPUs" has "^ *cpus += $((AT + 1))\$" "$AB2"
check "--max_memory 9000.MB, not whole GB: the pool leaves memory to Nextflow" lacks "memory" "$AB2"
rm -rf "$G3" "$ENG"

# --- first run, with fs.protected_hardlinks 0 ---------------------------------------------
PH=$(cat /proc/sys/fs/protected_hardlinks 2>/dev/null)
restore_ph() { [ -z "$PH" ] || sudo -n sysctl -q -w "fs.protected_hardlinks=${PH}" > /dev/null 2>&1; }
trap restore_ph EXIT
check "fs.protected_hardlinks is set to 0 for the first run (it was ${PH:-unknown})" sudo -n sysctl -q -w fs.protected_hardlinks=0
launch first
restore_ph
check_eq "run-all.sh exits 0" "$RC" 0
FIRST_TRACE=$(latest_trace)
check "the first run wrote a trace" test -s "${FIRST_TRACE:-/dev/null}"
check_ge "tasks completed in the first run" "$(status_count "${FIRST_TRACE:-/dev/null}" COMPLETED)" 10
check "the samplesheet is a FASTQ row" has "^${P},${G2}/${P}/fastq/${P}_R1.fastq.gz,${G2}/${P}/fastq/${P}_R2.fastq.gz,,,,,male\$" \
  "$(sed -n 2p "${G2}/${P}/nextflow/samplesheet.csv")"

# --- 4. the caps are the runner's ---------------------------------------------------------
NFCMD=$(grep -m1 -F '$> nextflow run' "${G2}/${P}/nextflow/.nextflow.log" 2>/dev/null || true)
echo "+ ${NFCMD#*\$> }"
check "nextflow got --max_cpus with the runner's CPU count" has " --max_cpus $(getconf _NPROCESSORS_ONLN) " "${NFCMD} "
check "nextflow got --max_memory with the runner's RAM" has " --max_memory $(awk '/^MemTotal:/ {print int($2 / 1048576)}' /proc/meminfo)\.GB " "${NFCMD} "

# --- 5. published by hard link ------------------------------------------------------------
check "nextflow got --publish_dir_mode link" has " --publish_dir_mode link " "${NFCMD} "
PBAM="${G2}/${P}/aligned/${P}_sorted.bam"
check_eq "the published BAM has two links (work/ and aligned/)" "$(stat -c %h "$PBAM" 2>/dev/null)" 2
cp "$PBAM" "${CASE_TMP}/copied.bam" 2>/dev/null
check_eq "negative control: a copy of it has one link" "$(stat -c %h "${CASE_TMP}/copied.bam" 2>/dev/null)" 1
rm -f "${CASE_TMP}/copied.bam"
# symlinks DIR: the symbolic links under DIR (not its nextflow/ launch directory), relative
symlinks() { (cd "$1" 2>/dev/null && find . -path ./nextflow -prune -o -type l -print | LC_ALL=C sort); }
SYMLINKS=$(LC_ALL=C comm -13 <(symlinks "$DIRECT") <(symlinks "${G2}/${P}"))
[ -z "$SYMLINKS" ] || printf 'symbolic links the direct (copy) run does not have:\n%s\n' "$SYMLINKS"
check_eq "published symbolic links the direct run does not have" "$(grep -c . <<<"$SYMLINKS" || true)" 0

# --- 1. the same published files as the direct run, and parity with the bash steps -------
published "$DIRECT" > "${CASE_TMP}/direct.txt"
published "${G2}/${P}" > "${CASE_TMP}/launcher.txt"
check_ge "files the direct run published" "$(wc -l < "${CASE_TMP}/direct.txt" | tr -d ' ')" 20
DIFF=$(diff "${CASE_TMP}/direct.txt" "${CASE_TMP}/launcher.txt")
[ -z "$DIFF" ] || printf 'published files that differ (< direct, > launcher):\n%s\n' "$DIFF"
check_eq "published files that differ from the direct run" "$(grep -c '^[<>]' <<<"$DIFF" || true)" 0
sed 1d "${CASE_TMP}/launcher.txt" > "${CASE_TMP}/planted.txt"
check "negative control: a planted missing file is reported" \
  test "$(diff "${CASE_TMP}/direct.txt" "${CASE_TMP}/planted.txt" | grep -c '^<' || true)" -eq 1
echo "+ scripts/ci/parity-diff.sh ${GENOME_DIR}/${P} ${G2}/${P} ${P}"
GITHUB_STEP_SUMMARY=/dev/null "${REPO}/scripts/ci/parity-diff.sh" "${GENOME_DIR}/${P}" "${G2}/${P}" "$P"
check_eq "parity of the launcher's outputs with the bash steps (exit)" "$?" 0

# --- 3. reports and step status -----------------------------------------------------------
check "the HTML report (step 24) is written" test -s "${G2}/${P}/${P}_report.html"
check "the text report is written" test -s "${G2}/${P}/${P}_report.txt"
check "summary.json validates" "${REPO}/bin/og" dev schema "${REPO}/tests/schema/summary.schema.json" "${G2}/${P}/summary.json"
check "OwnHealth's facts file validates against its schema" "${REPO}/bin/og" dev schema \
  "${REPO}/assets/ownhealth/genome-facts.v1.json" "${G2}/${P}/facts/${P}_genome_facts.json"
STATUS=$(cat "${G2}/${P}/logs/run_status.tsv" 2>/dev/null)
check "run_status.tsv records step 07 ok" has $'^step\t07\tok$' "$STATUS"
check "run_status.tsv records step 13 skipped" has $'^step\t13\tskipped ' "$STATUS"

# --- 2. the second run starts from the BAM and VCF, the third finds every task cached ------
launch second
NFCMD2=$(grep -m1 -F '$> nextflow run' "${G2}/${P}/nextflow/.nextflow.log" 2>/dev/null || true)
if [ "$PH" = 1 ]; then
  check "the second run, with fs.protected_hardlinks ${PH}, passes no link mode" lacks " --publish_dir_mode link " "${NFCMD2} "
fi
check_eq "the second run exits 0" "$RC" 0
SECOND_TRACE=$(latest_trace)
check "the second run wrote its own trace" test "${SECOND_TRACE:-}" != "${FIRST_TRACE:-}"
RA="${G2}/${P}/aligned/${P}_sorted.bam" RV="${G2}/${P}/vcf/${P}.vcf.gz" RG="${G2}/${P}/vcf/${P}.g.vcf.gz"
ROW2="${P},,,${RA},${RA}.bai,${RV},${RV}.tbi,male,${RG},${RG}.tbi"
check_eq "the second run's samplesheet names the BAM, VCF and gVCF the first run published" \
  "$(sed -n 2p "${G2}/${P}/nextflow/samplesheet.csv")" "$ROW2"
# A trace that is missing, empty or holds only its header would give 0 above for nothing
check_ge "tasks in the second run's trace" "$(awk 'NR > 1' "${SECOND_TRACE:-/dev/null}" 2>/dev/null | grep -c . || true)" 10
check_eq "FASTP, minimap2 and DeepVariant tasks the second run started" "$(upstream_tasks "${SECOND_TRACE:-/dev/null}")" 0
check_ge "negative control: the first run's trace shows them" "$(upstream_tasks "${FIRST_TRACE:-/dev/null}")" 3
launch third
check_eq "the third run exits 0" "$RC" 0
THIRD_TRACE=$(latest_trace)
check "the third run wrote its own trace" test "${THIRD_TRACE:-}" != "${SECOND_TRACE:-}"
echo "third run: $(status_count "${THIRD_TRACE:-/dev/null}" CACHED) cached, $(status_count "${THIRD_TRACE:-/dev/null}" COMPLETED) completed, ${SECS}s"
check "every task of the third run came from the cache" all_cached "${THIRD_TRACE:-/dev/null}"
check "negative control: the first run's trace fails the all-cached check" not_all_cached "${FIRST_TRACE:-/dev/null}"
check "the third run takes under 10 minutes (${SECS}s)" test "$SECS" -lt 600
check_eq "the samplesheet is the same in the third run" "$(sed -n 2p "${G2}/${P}/nextflow/samplesheet.csv")" "$ROW2"

finish
