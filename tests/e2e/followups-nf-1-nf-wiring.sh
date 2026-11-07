#!/usr/bin/env bash
# The Nextflow report and couple screen get the step 35 tables, the gVCF and
# the release the bash path already uses; CPSR's read check counts under step
# 17; a Locityper database without DONE stops the run before any task.
#   1. HTML_REPORT, run for real in its image through a two-line pipeline, on
#      HG002's VCF, a gVCF under a samplesheet-style name (meta.gvcf), and
#      made-up CYP21A2 and gene-dosage tables: the gVCF and its index are
#      linked as vcf/<id>.g.vcf.gz(.tbi), summary.json carries the CYP21A2
#      deletion and the HP count, and its manifest has the run
#      pipeline_version row with the manifest version of the pipeline's
#      config. Red controls: a copy of the module without the two place
#      alternatives stops on the CYP21A2 table; a copy without the gVCF and
#      release lines links no gVCF and writes no release.
#   2. main.nf passes them on: a stub run of the stub pedigree (two partners,
#      both called by DEEPVARIANT) with clinvar, parascopy, cpsr and
#      html_report stages the CYP21A2 and dosage tables in every HTML_REPORT
#      task that got the HBA table, each partner's CYP21A2 table in
#      COUPLE_CARRIER, and DEEPVARIANT's gVCF in the HTML_REPORT task of each
#      called sample. Its process table has a CPSR_CDNA_CHECK row, which
#      bin/run_state.py maps to step 17. Red control: a copy of the pipeline
#      without the new lines stages none of them.
#   3. --tools locityper with a --locityper_data folder that has no DONE stops
#      a -preview run (the workflow is built, no task runs) with the setup.sh
#      command; with DONE it
#      passes. Red control: a copy of main.nf without the check takes the
#      folder in silence.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
check "an earlier case left ${VCF#"$GENOME_DIR"/}" nonempty "${VCF#"$GENOME_DIR"/}"
DECLARED=$(sed -nE "/^manifest[[:space:]]*\{/,/^\}/ s/^[[:space:]]*version[[:space:]]*=[[:space:]]*'([^']+)'.*/\1/p" \
  "${REPO}/nextflow.config")
check "nextflow.config declares a manifest version (${DECLARED})" test -n "$DECLARED"

# --- 1. HTML_REPORT ----------------------------------------------------------------
S=NFF
T="${CASE_TMP}/tables"
rm -rf "$T"; mkdir -p "$T"
ln -f "$VCF" "${T}/given.gvcf.gz" 2>/dev/null || cp "$VCF" "${T}/given.gvcf.gz"
ln -f "${VCF}.tbi" "${T}/given.gvcf.gz.tbi" 2>/dev/null || cp "${VCF}.tbi" "${T}/given.gvcf.gz.tbi"
printf 'metric\tvalue\nstatus\tcalled\ncyp21a2_copies\t1\nmodule_total\t3\nexpected_total\t4\nestimate\t1.02\nfinding\tdeletion\nreason\t.\nreading\tplanted\n' \
  > "${T}/${S}_cyp21a2_psv.tsv"
printf 'locus\tstatus\tcopies\testimate\treason\treading\nHP\tcalled\t3\t3.04\t.\tplanted\n' > "${T}/${S}_cn_dosage.tsv"

NP="${CASE_TMP}/project"
rm -rf "$NP"; mkdir -p "$NP"
cp -R "${REPO}/bin" "${NP}/bin"
nfconf() {  # nfconf OUT
  printf "params.outdir = '%s'\nparams.publish_dir_mode = 'copy'\nmanifest.version = '%s'\ndocker.enabled = true\ndocker.runOptions = '--user %s:%s --network none'\n" \
    "$1" "$DECLARED" "$(id -u)" "$(id -g)" > "${NP}/nextflow.config"
}
nfrun() {  # nfrun NAME: the pipeline in NP/main.nf, log CASE_TMP/NAME.log
  echo "+ nextflow run ($1)"
  ( cd "$CASE_TMP" && nextflow run "${NP}/main.nf" -ansi-log false -c "${REPO}/conf/containers.config" \
      -work-dir "${CASE_TMP}/work-$1" ) > "${CASE_TMP}/$1.log" 2>&1
  NF_RC=$?
  [ "$NF_RC" = 0 ] || grep -E 'ERROR|Error executing process|Command error' -A 8 "${CASE_TMP}/$1.log" | head -40
}
report_run() {  # report_run NAME MODULE
  rm -rf "${CASE_TMP:?}/$1" "${CASE_TMP:?}/work-$1"
  nfconf "${CASE_TMP}/$1"
  cat > "${NP}/main.nf" <<NF
include { HTML_REPORT } from '$2'

workflow {
    HTML_REPORT(channel.of(
        [[id: '${S}', sex: 'male', gvcf: ['given.gvcf.gz', 'given.gvcf.gz.tbi']], file('${VCF}'),
         [file('${T}/${S}_cyp21a2_psv.tsv'), file('${T}/${S}_cn_dosage.tsv'),
          file('${T}/given.gvcf.gz'), file('${T}/given.gvcf.gz.tbi')]]))
}
NF
  nfrun "$1"
}
# shellcheck disable=SC2329  # called through check
# linked NAME FILE: the layout link FILE of the task in work-NAME resolves to a file named "given.gvcf.gz..."
linked() {
  local l
  l=$(find "${CASE_TMP}/work-$1" -path "*/layout/${S}/vcf/$2" 2>/dev/null | head -1)
  [ -n "$l" ] && [ -s "$l" ] && case "$(basename "$(readlink -f "$l")")" in "given.gvcf.gz"*) true ;; *) false ;; esac
}
summary() {  # summary NAME EXPR: a Python expression over d, the summary.json of run NAME
  python3 -c "import json, sys
try:
    d = json.load(open(sys.argv[1]))
except (OSError, ValueError):
    print('no summary'); sys.exit(0)
print($2)" "${CASE_TMP}/$1/${S}/summary.json" 2>&1
}
report_run report "${REPO}/modules/local/html_report/main"
check_eq "HTML_REPORT with the gVCF and both tables exits 0" "$NF_RC" 0
check "the gVCF is linked as vcf/${S}.g.vcf.gz" linked report "${S}.g.vcf.gz"
check "and its index as vcf/${S}.g.vcf.gz.tbi" linked report "${S}.g.vcf.gz.tbi"
PARA="(d['sections'].get('parascopy') or {}).get('data') or {}"
check_eq "summary.json: the CYP21A2 finding and copies" \
  "$(summary report "[($PARA).get('cyp21a2', {}).get(k) for k in ('status', 'finding', 'cyp21a2_copies')]")" \
  "['called', 'deletion', '1']"
check_eq "summary.json: the HP dosage count" \
  "$(summary report "[(r['locus'], r['copies']) for r in ($PARA).get('dosage') or []]")" "[('HP', '3')]"
check_eq "summary.json: the manifest names the release" \
  "$(summary report "[r[2] for r in d['manifest'] if r[:2] == ['run', 'pipeline_version']]")" "['${DECLARED}']"
check "the report itself shows the release" grep -q "pipeline_version</td><td>${DECLARED}<" "${CASE_TMP}/report/${S}/${S}_report.html"
check "the run manifest is not published (run-all.sh writes the full one there)" test ! -e "${CASE_TMP}/report/${S}/run_manifest.tsv"

mkdir -p "${NP}/noplace" "${NP}/nogvcf"
# shellcheck disable=SC2016  # the dollar signs are sed text
sed 's/|\${id}_cn_dosage\.tsv|\${id}_cyp21a2_psv\.tsv)/)/' "${REPO}/modules/local/html_report/main.nf" > "${NP}/noplace/main.nf"
check_eq "control copy: the module has no place for the two tables" \
  "$(grep -cE '_cn_dosage\.tsv\||_cyp21a2_psv\.tsv\)' "${NP}/noplace/main.nf")" 0
report_run report-noplace "${NP}/noplace/main"
check "red: the module without the place alternatives stops (exit ${NF_RC})" test "$NF_RC" -ne 0
check "and names the table it has no place for" grep -qE "no place for: ${S}_(cyp21a2_psv|cn_dosage)\.tsv" -r "${CASE_TMP}/work-report-noplace" --include=.command.err
sed -e 's/^\( *\)ln -s .*\.g\.vcf\.gz.*$/\1true/' -e '/run_manifest\.tsv"$/d' "${REPO}/modules/local/html_report/main.nf" > "${NP}/nogvcf/main.nf"
check_eq "control copy: the gVCF links and the release line are gone" \
  "$(grep -cE 'g\.vcf\.gz"|run_manifest\.tsv"$' "${NP}/nogvcf/main.nf")" 0
report_run report-nogvcf "${NP}/nogvcf/main"
check_eq "control: that module still runs" "$NF_RC" 0
check "red: it links no gVCF" test -z "$(find "${CASE_TMP}/work-report-nogvcf" -path "*/layout/${S}/vcf/${S}.g.vcf.gz" 2>/dev/null)"
check_eq "red: and its manifest has no release" \
  "$(summary report-nogvcf "[r[2] for r in d['manifest'] if r[:2] == ['run', 'pipeline_version']]")" "[]"

# --- 2. main.nf: stub run of the stub pedigree ---------------------------------------
staged_in() {  # staged_in WORK PROCESS PATTERN: tasks of PROCESS in WORK that stage a file matching PATTERN (ERE)
  find "$1" -name .command.run -exec grep -l "^### name: .*$2" {} + 2>/dev/null \
    | while read -r f; do grep -lE -- "$3" "$f"; done | wc -l | tr -d ' '
}
for n in stub stub-nowire; do
  rm -rf "${CASE_TMP:?}/${n}"; mkdir -p "${CASE_TMP}/${n}"
  git -C "$REPO" ls-files -z | (cd "$REPO" && tar --null -T - -cf -) | tar -xf - -C "${CASE_TMP}/${n}"
done
W1='BAM_ANALYSIS.out.cn_dosage.map'; W2='BAM_ANALYSIS.out.cyp21a2.map'
W3='SV.out.cnvpytor_vcf, BAM_ANALYSIS.out.cyp21a2$'; W4='.join(ch_gvcf.map { meta, gvcf, idx'
check_eq "main.nf joins the dosage and CYP21A2 tables into the report inputs once each" \
  "$(grep -cF -e "$W1" -e "$W2" "${CASE_TMP}/stub/main.nf")" 2
check_eq "main.nf hands the CYP21A2 table to the couple screen once" "$(grep -c "$W3" "${CASE_TMP}/stub/main.nf")" 1
check_eq "main.nf joins the gVCF into the report inputs once" "$(grep -cF "$W4" "${CASE_TMP}/stub/main.nf")" 1
sed -i -e "/${W1}/d" -e "/${W2}/d" -e 's/SV.out.cnvpytor_vcf, BAM_ANALYSIS.out.cyp21a2$/SV.out.cnvpytor_vcf/' \
  -e 's/\.join(ch_gvcf\.map { meta, gvcf, idx/.join(ch_gvcf.filter { false }.map { meta, gvcf, idx/' \
  "${CASE_TMP}/stub-nowire/main.nf"
check_eq "control copy: the new lines are gone" \
  "$(grep -cF -e "$W1" -e "$W2" -e "$W4" "${CASE_TMP}/stub-nowire/main.nf")$(grep -c "$W3" "${CASE_TMP}/stub-nowire/main.nf")" 00
echo '{"kir": false}' > "${CASE_TMP}/no-kir.json"
for n in stub stub-nowire; do
  echo "+ nextflow run main.nf -profile test_all -stub --pedigree assets/stub/pedigree.tsv (${n})"
  ( cd "${CASE_TMP}/${n}" && nextflow run main.nf -profile test_all -stub -ansi-log false --outdir "${CASE_TMP}/${n}/out" \
      -work-dir "${CASE_TMP}/${n}/work" --pedigree assets/stub/pedigree.tsv --tools clinvar,parascopy,cpsr,html_report \
      -params-file "${CASE_TMP}/no-kir.json" ) > "${CASE_TMP}/${n}.log" 2>&1
  rc=$?
  tail -n 6 "${CASE_TMP}/${n}.log"
  check_eq "stub run (${n}) exits 0" "$rc" 0
done
SW="${CASE_TMP}/stub/work"; S0="${CASE_TMP}/stub-nowire/work"
HBA=$(staged_in "$SW" HTML_REPORT '_hba_copy_number[.]tsv')
check_ge "HTML_REPORT tasks with the HBA table" "$HBA" 1
check_eq "every one of them stages the CYP21A2 table" "$(staged_in "$SW" HTML_REPORT '_cyp21a2_psv[.]tsv')" "$HBA"
check_eq "and the dosage table" "$(staged_in "$SW" HTML_REPORT '_cn_dosage[.]tsv')" "$HBA"
check_eq "both partners' CYP21A2 tables reach COUPLE_CARRIER" \
  "$(find "$SW" -name .command.run -exec grep -l '^### name: .*COUPLE_CARRIER' {} + 2>/dev/null | while read -r f; do grep -ho 'stub_[a-z]*_cyp21a2_psv.tsv' "$f"; done | sort -u | wc -l | tr -d ' ')" 2
for id in stub_fastq stub_bam; do
  check_eq "the HTML_REPORT task of ${id} stages DEEPVARIANT's gVCF" \
    "$(staged_in "$SW" "HTML_REPORT .${id}." "${id}[.]g[.]vcf[.]gz")" 1
done
check_eq "red: without the lines no HTML_REPORT task stages the CYP21A2 or dosage table" \
  "$(staged_in "$S0" HTML_REPORT '_cyp21a2_psv[.]tsv|_cn_dosage[.]tsv')" 0
check_ge "control: the stub run without the lines still ran COUPLE_CARRIER" "$(staged_in "$S0" COUPLE_CARRIER '_clinvar_hits[.]vcf')" 1
check_eq "red: and COUPLE_CARRIER stages no CYP21A2 table" "$(staged_in "$S0" COUPLE_CARRIER '_cyp21a2_psv[.]tsv')" 0
check_eq "red: and no HTML_REPORT task stages a gVCF" "$(staged_in "$S0" HTML_REPORT '[.]g[.]vcf[.]gz')" 0
# CPSR's read check under step 17: the run's process table, read the way run-all.sh and the live view read it
STEP17=$(python3 - "${CASE_TMP}/stub/process_status.tsv" "${REPO}" <<'PY'
import sys
sys.path.insert(0, sys.argv[2] + "/bin")
import run_state
ps = run_state.read_process_status(sys.argv[1]) or {}
step_of = {p: s["step"] for s in run_state.read_steps(sys.argv[2] + "/assets/steps.tsv") for p in s["processes"]}
r = ps.get("CPSR_CDNA_CHECK") or {}
print(step_of.get("CPSR_CDNA_CHECK"), sum(v for v in r.values() if isinstance(v, int)))
PY
)
echo "CPSR_CDNA_CHECK: step and tasks: ${STEP17}"
check "the process table has CPSR_CDNA_CHECK tasks, and run_state.py maps them to step 17" \
  has '^17 [1-9]' "$STEP17"

# --- 3. Locityper's DONE, before any task ------------------------------------------
LD="${CASE_TMP}/locityper-db"
rm -rf "$LD"; mkdir -p "$LD"
cp -R "${CASE_TMP}/stub" "${CASE_TMP}/nocheck" 2>/dev/null || true
rm -rf "${CASE_TMP:?}/nocheck/work" "${CASE_TMP:?}/nocheck/out"
# shellcheck disable=SC2016  # the dollar signs are sed text
sed -i '/!file("${params.locityper_data}\/DONE").exists()) {/,/^    }$/d' "${CASE_TMP}/nocheck/main.nf"
sed -i '/tools_list.contains(.locityper.) && params.locityper_data && !workflow.stubRun$/d' "${CASE_TMP}/nocheck/main.nf"
check_eq "control copy: the DONE check is gone" "$(grep -c 'locityper_data} has no DONE' "${CASE_TMP}/nocheck/main.nf")" 0
preview() {  # preview NAME DIR: a -preview run of DIR/main.nf with locityper and the database LD
  ( cd "$2" && nextflow run main.nf -profile test -preview -ansi-log false --outdir "${CASE_TMP}/$1-out" \
      -work-dir "${CASE_TMP}/$1-work" --tools locityper --locityper_data "$LD" ) > "${CASE_TMP}/$1.log" 2>&1
  PV_RC=$?
  tail -n 4 "${CASE_TMP}/$1.log"
}
preview nodone "${CASE_TMP}/stub"
check "a database without DONE stops the run (exit ${PV_RC})" test "$PV_RC" -ne 0
check "with the setup.sh command" grep -q "locityper_data ${LD} has no DONE: run scripts/setup.sh --locityper-data" "${CASE_TMP}/nodone.log"
preview nocheck "${CASE_TMP}/nocheck"
check_eq "red: without the check the same folder passes the preview" "$PV_RC" 0
touch "${LD}/DONE"
preview with-done "${CASE_TMP}/stub"
check_eq "control: with DONE the preview passes" "$PV_RC" 0

rm -rf "${CASE_TMP:?}/stub" "${CASE_TMP:?}/stub-nowire" "${CASE_TMP:?}/nocheck" "${CASE_TMP:?}"/work-* "${CASE_TMP:?}"/*-work
finish
