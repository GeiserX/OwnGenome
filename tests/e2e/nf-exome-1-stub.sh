#!/usr/bin/env bash
# --model_type WES on the Nextflow path: stub runs of main.nf (profile
# test_all, no containers) on the stub BAM-only row, so DEEPVARIANT calls it.
# The stubs of DEEPVARIANT, MOSDEPTH and MANTA print the options their
# scripts pass, so each task's .command.sh shows them:
#   1. WES: DeepVariant gets the exome model and the capture BED as its
#      regions, mosdepth the capture BED as its windows, Manta --exome;
#      INDEXCOV does not run and somalier (sample_qc) does; CNVpytor, ROH,
#      TelomereHunter, archaic, MoChA and viral do not run and the log says
#      "not run: exome"; mito_variants runs, MITO_DELETIONS (the mtDNA copy
#      number) does not;
#      survivor_merge keeps its two callers, Manta and Delly;
#   2. WES with Manta and CNVpytor as the only callers: survivor_merge goes
#      with CNVpytor, and the run still starts;
#   3. the default (WGS): the WGS model, no regions, 500 bp windows, no
#      --exome, INDEXCOV and the whole-genome steps run;
#   4. the stops: WES without --capture_bed, --capture_bed without WES, an
#      unknown model, WES without sample_qc, WES with --intervals.
# The stubs only touch files; nothing here needs a database or an image.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

ST="${REPO}/assets/stub"
BED="${CASE_TMP}/capture_regions.bed"
printf 'chr1\t1000\t2000\tEXON1\nchr2\t5000\t5300\tEXON2\n' > "$BED"
WES_TOOLS=mosdepth,manta,delly,survivor_merge,sample_qc,hla_typing,cnvpytor,roh,telomere_hunter,archaic,mocha,viral,mito_variants
# stub_run NAME [ARGS...]: one stub run in its own launch directory; sets RC, LOG and DIR
stub_run() {
  local name=$1
  shift
  DIR="${CASE_TMP}/${name}"
  rm -rf "$DIR"
  mkdir -p "$DIR"
  printf 'sample,bam,bam_index,sex\nex1,%s,%s,male\n' "${ST}/sample.bam" "${ST}/sample.bam.bai" > "${DIR}/samplesheet.csv"
  echo "+ nextflow run main.nf -profile test_all -stub (${name}) $*"
  ( cd "$DIR" && nextflow run "${REPO}/main.nf" -profile test_all -stub -ansi-log false \
      -work-dir "${DIR}/work" --input "${DIR}/samplesheet.csv" --outdir "${DIR}/out" "$@" ) > "${DIR}/run.log" 2>&1
  RC=$?
  LOG=$(cat "${DIR}/run.log")
  echo "$LOG"
  cp "${DIR}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.${name}.nextflow.log" 2>/dev/null
}
# task_sh TEXT: the .command.sh of the run's task whose script holds TEXT
task_sh() { find "${DIR}/work" -name .command.sh -exec grep -l -- "$1" {} + 2>/dev/null | head -n 1; }
# exome_line: the run log's "not run: exome" line
exome_line() { grep -- '--model_type WES: not run: exome:' "${DIR}/run.log" | head -n 1; }

# --- 1. an exome -------------------------------------------------------------------------
stub_run wes --model_type WES --capture_bed "$BED" --tools "$WES_TOOLS"
check_eq "WES: the run exits 0" "$RC" 0
DV_SH=$(task_sh 'run_deepvariant options:')
MD_SH=$(task_sh 'mosdepth options:')
MT_SH=$(task_sh 'configManta.py options:')
check "WES: DeepVariant's task ran" test -n "$DV_SH"
check "WES: DeepVariant runs the exome model on the capture BED" \
  grep -q -- 'run_deepvariant options: --model_type=WES --regions capture_regions.bed' "${DV_SH:-/dev/null}"
check "WES: mosdepth measures the capture BED" grep -q -- 'mosdepth options: --by capture_regions.bed' "${MD_SH:-/dev/null}"
check "WES: Manta runs with --exome" grep -q -- 'configManta.py options: .*--exome' "${MT_SH:-/dev/null}"
check "WES: indexcov does not run" test -z "$(find "${DIR}/work" -name .command.sh -exec grep -l 'goleft' {} + 2>/dev/null)"
check "WES: no sex check file from indexcov" test ! -e "${DIR}/out/ex1/ex1_sex_check.tsv"
check "WES: somalier's sex check runs (SAMPLE_QC)" test -e "${DIR}/out/ex1/qc/ex1_sample_qc.tsv"
L=$(exome_line)
check "WES: the log says which steps are not run: exome" test -n "$L"
for t in cnvpytor roh telomere_hunter archaic mocha viral; do
  check "WES: the not-run line names ${t}" has "[ ,]${t}([,.]|\$)" "$L"
done
check "WES: survivor_merge keeps Manta and Delly, so it is not named" lacks 'survivor_merge' "$L"
for p in CNVPYTOR ROH TELOMERE_HUNTER ARCHAIC_GENOTYPES MOCHA_PREP VIRAL_READS MITO_DELETIONS INDEXCOV; do
  check "WES: no ${p} task" test -z "$(grep -E "(^|[ :>])${p} " "${DIR}/run.log")"
done
check "WES: MITO_VARIANTS still runs" test -n "$(grep -E '[ :>]MITO_VARIANTS ' "${DIR}/run.log")"
check "WES: SURVIVOR_MERGE still runs" test -n "$(grep -E '[ :>]SURVIVOR_MERGE ' "${DIR}/run.log")"
check "WES: the mito line says what is not run" has 'mito_variants: not run: exome: the mtDNA copy number' "$LOG"
check "WES: no mtDNA copy number is published" test -z "$(find "${DIR}/out/ex1" -name '*_mtdna_copies.tsv' 2>/dev/null)"

# --- 2. an exome where CNVpytor was one of survivor_merge's two callers ------------------
stub_run wes-two --model_type WES --capture_bed "$BED" --tools mosdepth,manta,cnvpytor,survivor_merge,sample_qc,hla_typing
check_eq "WES, Manta and CNVpytor: the run exits 0" "$RC" 0
L=$(exome_line)
check "WES, Manta and CNVpytor: the not-run line names survivor_merge" has 'survivor_merge' "$L"
check "WES, Manta and CNVpytor: no SURVIVOR_MERGE task" test -z "$(grep -E '[ :>]SURVIVOR_MERGE ' "${DIR}/run.log")"
check "WES, Manta and CNVpytor: Manta still runs" test -n "$(task_sh 'configManta.py options:')"

# --- 3. a whole genome, the default --------------------------------------------------------
stub_run wgs --tools "$WES_TOOLS"
check_eq "WGS: the run exits 0" "$RC" 0
DV_SH=$(task_sh 'run_deepvariant options:')
MD_SH=$(task_sh 'mosdepth options:')
MT_SH=$(task_sh 'configManta.py options:')
check "WGS: DeepVariant runs the WGS model" grep -q -- 'run_deepvariant options: --model_type=WGS' "${DV_SH:-/dev/null}"
check "WGS: DeepVariant calls the whole genome (no --regions)" test -n "$DV_SH" -a -z "$(grep -- '--regions' "${DV_SH:-/dev/null}")"
check "WGS: mosdepth uses 500 bp windows" grep -q -- 'mosdepth options: --by 500' "${MD_SH:-/dev/null}"
check "WGS: Manta runs without --exome" test -n "$MT_SH" -a -z "$(grep -- '--exome' "${MT_SH:-/dev/null}")"
check "WGS: indexcov runs" test -e "${DIR}/out/ex1/ex1_sex_check.tsv"
check "WGS: no not-run line" test -z "$(exome_line)"
for p in CNVPYTOR ROH TELOMERE_HUNTER ARCHAIC_GENOTYPES MOCHA_PREP VIRAL_READS MITO_DELETIONS; do
  check "WGS: ${p} runs" test -n "$(grep -E "[ :>]${p} " "${DIR}/run.log")"
done

# --- 4. the stops ------------------------------------------------------------------------
stub_run no-bed --model_type WES --tools mosdepth,sample_qc,hla_typing
check "WES without --capture_bed: the run stops (exit ${RC})" test "$RC" -ne 0
check "WES without --capture_bed: the message names it" has '--model_type WES needs --capture_bed' "$LOG"

stub_run bed-only --capture_bed "$BED" --tools mosdepth,hla_typing
check "--capture_bed without WES: the run stops (exit ${RC})" test "$RC" -ne 0
check "--capture_bed without WES: the message says to add the model" has '--capture_bed holds the capture regions of an exome: add --model_type WES' "$LOG"

stub_run bad-model --model_type PACBIO --tools mosdepth,hla_typing
check "--model_type PACBIO: the run stops (exit ${RC})" test "$RC" -ne 0
check "--model_type PACBIO: the message names the two values" has "--model_type must be WGS .* or WES .*, got 'PACBIO'" "$LOG"

stub_run no-qc --model_type WES --capture_bed "$BED" --tools mosdepth,hla_typing
check "WES without sample_qc: the run stops (exit ${RC})" test "$RC" -ne 0
check "WES without sample_qc: the message says to add it" has '--model_type WES takes the sex check from somalier: add sample_qc' "$LOG"

stub_run intervals --model_type WES --capture_bed "$BED" --intervals chr1:1-5000 --tools mosdepth,sample_qc,hla_typing
check "WES with --intervals: the run stops (exit ${RC})" test "$RC" -ne 0
check "WES with --intervals: the message says why" has '--intervals and --capture_bed both set' "$LOG"
for n in no-bed bed-only bad-model no-qc intervals; do
  check_eq "${n}: no task ran" "$(find "${CASE_TMP}/${n}/work" -name .command.run 2>/dev/null | wc -l | tr -d ' ')" 0
done

finish
