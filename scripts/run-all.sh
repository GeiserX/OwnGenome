#!/usr/bin/env bash
# run-all.sh: the whole pipeline for one sample: validate-setup.sh, then main.nf with -resume (a rerun redoes only what changed).
# Usage: GENOME_DIR=/data ./scripts/run-all.sh <sample> <male|female> [nextflow options]
#   Options after the sex go to `nextflow run` as given (--sex_check warn), except -bg. --max_cpus (THREADS) and --max_memory
#   cap each task (default: the host's); they do not limit how much runs at once: Nextflow fills the machine's CPUs and RAM.
#   When the container engine has less (Docker Desktop and colima run the containers in a VM), the defaults are the
#   engine's, and <sample>/nextflow/executor.config makes Nextflow fill only that.
#   A full run on 8 CPUs can take more than a day: start it inside tmux or screen, or with nohup. A closed terminal stops
#   it, and -resume then restarts the task that was running (DeepVariant is one task).
# Needs Docker, bash 4.4+, Java 17+ and Nextflow (NEXTFLOW_VERSION in versions.env). Results: GENOME_DIR/<sample>/.
# Input, first match: aligned/<sample>_sorted.bam with .bai (plus vcf/<sample>.vcf.gz with .tbi: not called
#   again, and the gVCF beside it, vcf/<sample>.g.vcf.gz with .tbi, goes to PharmCAT and PRS); fastq/<sample>_R1.fastq.gz
#   and _R2; the VCF alone. Kept in <sample>/nextflow/samplesheet.csv while its files exist, its BAM is this call's
#   and it is no lower in that order or its run did not end: after a run that ended, a FASTQ row gives way to the BAM
#   and VCF that run published.
# Switches: SKIP_VALIDATION=true; THREADS=N (--max_cpus); SKIP_TRIM=true (--skip_trim);
#   INTERVALS="chr20 chr22" (--intervals); ALIGN_DIR=dir (the BAM from <sample>/dir/); TOOLS=a,b
#   (only these --tools names; parascopy, paraphase, y_haplogroup, y_str, consumer_export, bloodtype, ancient, lpa, vntyper, viral, chip, mei, mocha and locityper run only when named; paraphase needs the BAM
#   to be PacBio HiFi, PL:PACBIO, or its step fails and the rest goes on; cyrius runs when setup.sh --cyrius installed it, whose
#   licence is non-commercial: TOOLS without cyrius leaves it out); KIR=true (KIR genes with HLA);
#   HLA_QC=true (OptiType beside T1K: HLA-A and HLA-B reach PharmCAT only when both agree);
#   CPSR_SECONDARY_FINDINGS=false (--cpsr_secondary_findings false: step 17's CPSR reports the cancer genes only, not
#   the ACMG SF genes; the report's own secondary-findings card from step 23 stays);
#   HPO=HP:0001250,HP:0001263 (the symptoms as HPO ids: step 48 ranks the variants that could explain them with Exomiser,
#   after setup.sh --exomiser-data);
#   EVERYTHING=true (every optional step of assets/everything.tsv: every --tools name, KIR, HLA_QC, GRIDSS, TIDDIT,
#   IMPUTATION, SOMATIC, BENCHMARK and every EXTRA_CALLERS; TOOLS or a switch set to another value stops it with exit 2.
#   A step that would be skipped because its data is not installed, or because ANCESTRY_PANEL, ANCESTRY_FINE_DIR or
#   ARCHAIC_DIR is none, stops it with exit 2 before Nextflow, naming the setup command; a skip that comes from the
#   input stays a skip: no BAM, no gVCF, a female sample for the Y steps, a BAM that is not PacBio HiFi for Paraphase);
#   REF_FASTA, EH_CATALOG (a file, or bundled), MANTA_CALL_REGIONS, PARASCOPY_POPULATION, PARAPHASE_REGIONS, ANCESTRY_PANEL (--ancestry_ref; none: raw
#   scores), ANCESTRY_FINE_DIR (--ancestry_fine_panel, step 26b,
#   always 100 bootstrap replicates; none: skip), ARCHAIC_DIR and ARCHAIC_POP (--archaic_data and --archaic_pop,
#   step 43; ARCHAIC_DIR=none: skip), PGSC_CALC_DIR (--pgsc_calc), LOCITYPER_DATA (--locityper_data, step 53's
#   database) and LOCITYPER_BG_REGION (--locityper_bg_region) as for the single steps; the panels, the databases and
#   pgsc_calc are passed when installed.
# Run as scripts after the pipeline: TIDDIT=true (04a), GRIDSS=true (04b), IMPUTATION=true (14),
#   SOMATIC=true (29), EXTRA_CALLERS=gatk,freebayes,strelka2,octopus (03a-03d), BENCHMARK=true (needs
#   EXTRA_CALLERS or a second caller VCF); then the HTML report (24) and the text report, also after a failed pipeline.
# FAMILIES=families.tsv (who is related to whom, assets/families.example.tsv): when the sample is in a family, the sex
#   must agree with the file (exit 2 before the pipeline). After the pipeline, before the reports, og family
#   plan writes GENOME_DIR/<family>/family_plan.tsv from what the members have, and step 38 (38-family.sh) and step 39
#   (39-couple.sh) run when it says so, or are recorded as skipped with the reason (waiting: a member has no gVCF, or for step 39 no ClinVar hits, yet).
#   The other members are never run here: run run-all.sh for each, and the last one runs the family steps.
#   A step without its data or BAM is skipped. A report-only tool whose task fails is skipped, and its step is marked
#   failed (logs/run_status.tsv). Results are published by hard link when work/ and GENOME_DIR share a filesystem,
#   you passed no -w, -c or --publish_dir_mode, and you run as root or Linux's fs.protected_hardlinks is 0 (then a user
#   may link the tasks' root-owned files).
# logs/run_status.tsv keeps this command with its switches: bin/og recheck --rerun replays it when a newer
#   ClinVar or catalog changes a result (docs/recheck.md).
set -euo pipefail
case "${1:-}" in -h|--help) sed -n '2,/^set -euo/p' "$0" | sed -e '$d' -e 's/^# \{0,1\}//'; exit 0 ;; esac
SAMPLE=${1:-} SEX=${2:-}
[ -n "$SAMPLE" ] && [[ "$SEX" =~ ^(male|female)$ ]] || { { echo "Usage: $0 <sample> <male|female> [nextflow options]"
  [ -z "$SEX" ] || echo "ERROR: sex must be 'male' or 'female', got '${SEX}'."; } >&2; exit 2; }
((BASH_VERSINFO[0] * 100 + BASH_VERSINFO[1] >= 404)) || { echo "ERROR: run-all.sh needs bash 4.4 or later (this is ${BASH_VERSION}); on macOS: brew install bash" >&2; exit 2; }
shift 2; SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd) USER_THREADS=${THREADS:-}
# This command as given, for bin/og recheck --rerun (logs/run_status.tsv, 'meta command'): the switches
# below that are set, then the arguments, each shell-quoted, run from this directory so relative paths still hold.
# 'meta switches' names them all, so a replay drops one set now that the run did not have. The replay reads this
# shape back word by word and starts run-all.sh with those arguments, never through a shell; it refuses any other.
SWITCHES=(GENOME_DIR SKIP_VALIDATION THREADS SKIP_TRIM INTERVALS ALIGN_DIR TOOLS KIR HLA_QC HPO CPSR_SECONDARY_FINDINGS REF_FASTA EH_CATALOG
  MANTA_CALL_REGIONS PARASCOPY_POPULATION PARAPHASE_REGIONS ANCESTRY_PANEL ANCESTRY_FINE_DIR ARCHAIC_DIR ARCHAIC_POP
  PGSC_CALC_DIR ANCIENT_PANEL_DIR VIRAL_REF MEGANE_KMERS MOCHA_DATA GRIDSS TIDDIT IMPUTATION SOMATIC EXTRA_CALLERS BENCHMARK
  FAMILIES OWNGENOME_LIVE_URL LOCITYPER_DATA LOCITYPER_BG_REGION EVERYTHING NXF_VER NXF_WORK NXF_HOME NXF_OFFLINE)
RERUN_CMD="cd $(printf '%q' "$PWD") && "; for v in "${SWITCHES[@]}"; do [ -z "${!v+x}" ] || RERUN_CMD+="${v}=$(printf '%q' "${!v}") "; done
RERUN_CMD+=$(printf '%q ' "${SCRIPT_DIR}/run-all.sh" "$SAMPLE" "$SEX" "$@"); RERUN_CMD=${RERUN_CMD% }
# With -bg nextflow returns at once, and the steps would be recorded and the reports written before the run ends.
for a in "$@"; do [ "$a" != -bg ] || { echo "ERROR: run-all.sh does not take -bg: it waits for the pipeline to record each step and write the reports. Run it inside tmux or screen, or with nohup." >&2; exit 2; }; done
export GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "${SCRIPT_DIR}/lib/common.sh"
validate_sample "$SAMPLE"
# EVERYTHING=true: every tool row of assets/everything.tsv goes into TOOLS (html_report aside: this script writes
# the report itself, step 24) and every switch row is set to its value. EVHOW keeps each row's install command.
EV="" EVMISS=(); declare -A EVHOW=()
if [[ "${EVERYTHING:-}" =~ ^(true|1)$ ]]; then
  EV=1 EVT="" EVF="${PGP_ROOT}/assets/everything.tsv"
  [ -s "$EVF" ] || { echo "ERROR: EVERYTHING=true reads ${EVF}, which is missing" >&2; exit 2; }
  [ -z "${TOOLS:-}" ] || { echo "ERROR: EVERYTHING=true names every step itself: unset TOOLS (TOOLS=${TOOLS})" >&2; exit 2; }
  while IFS=$'\t' read -r kind name value _; do
    case "$kind" in
      tool) [ "$name" = html_report ] || EVT+="${EVT:+,}${name}"; EVHOW[$name]=$value ;;
      input) EVHOW[$name]=$value ;;
      switch) if [ -n "${!name:-}" ] && [ "${!name}" != "$value" ] && ! { [ "$value" = true ] && [ "${!name}" = 1 ]; }; then
          echo "ERROR: EVERYTHING=true sets ${name}=${value}: unset ${name} (it is ${!name})" >&2; exit 2; fi
        printf -v "$name" '%s' "$value" ;;
    esac
  done < <(grep -v '^#' "$EVF" | tail -n +2)
  TOOLS=$EVT
  echo "[Everything] EVERYTHING=true: every optional step (${EVF}); one that lacks its data stops the run before Nextflow"
fi
G=$GENOME_DIR S="${GENOME_DIR}/${SAMPLE}" LOG_DIR="${GENOME_DIR}/${SAMPLE}/logs" T=${TOOLS:-} C=${EXTRA_CALLERS:-}
on() { [[ "${!1:-}" =~ ^(true|1)$ ]]; }
SF=${CPSR_SECONDARY_FINDINGS:-true}  # as scripts/17-cpsr.sh reads it: unset or empty is true, else true or false only
[[ "$SF" =~ ^(true|false)$ ]] || { echo "ERROR: CPSR_SECONDARY_FINDINGS must be true or false, got '${SF}'" >&2; exit 2; }
# The ClinVar files are shared by every sample of GENOME_DIR. bin/og recheck --refresh-clinvar holds this
# lock while it replaces them; the run holds it shared on fd 7 until it ends, so a run waits for a refresh going on
# and a refresh refuses while a run goes on. Without flock (macOS: brew install flock) nothing waits here.
CVLOCK="${G}/clinvar/.refresh.lock"
if [ -d "${G}/clinvar" ] && command -v flock >/dev/null 2>&1; then
  [ -e "$CVLOCK" ] || : 2>/dev/null >> "$CVLOCK" || true
  if [ -r "$CVLOCK" ]; then
    exec 7<"$CVLOCK"
    flock -s -n 7 || { echo "Waiting for the ClinVar refresh that holds ${CVLOCK} to end ..."; flock -s 7; }
  fi
fi
java_major=$(java -version 2>&1 | awk -F'"' '/version "/ {split($2, v, "."); print (v[1] == 1 ? v[2] : v[1]); exit}') || true
if ! command -v nextflow >/dev/null 2>&1 || [[ ! "${java_major:-}" =~ ^[0-9]+$ ]] || [ "$java_major" -lt 17 ]; then
  printf '%s\n' "ERROR: run-all.sh starts the Nextflow pipeline and needs Java 17 or later (found: ${java_major:-none}) and Nextflow ($(command -v nextflow || echo 'not on PATH'))." \
    "  Install the release CI validates:  curl -s https://get.nextflow.io | NXF_VER=${NEXTFLOW_VERSION} bash && sudo mv nextflow /usr/local/bin/" \
    "  Or run the steps one by one: ./scripts/<step>.sh ${SAMPLE} (docs/getting-started.md)." >&2; exit 2; fi
on SKIP_VALIDATION || "${SCRIPT_DIR}/validate-setup.sh" "$SAMPLE" \
  || { echo "ERROR: setup validation failed; fix the items above. To bypass: SKIP_VALIDATION=true $0 ${SAMPLE} ${SEX}" >&2; exit 1; }
# og (bin/.og) built from another rust/ than this checkout's stops the run, SKIP_VALIDATION or not (rust/og.sh); before the family plan, which runs og.
[ ! -f "${PGP_ROOT}/rust/og.sh" ] || OG_MSG=$(bash "${PGP_ROOT}/rust/og.sh" check 2>&1) || { echo "$OG_MSG" >&2; exit 2; }
# The sample's family (FAMILIES): the file is checked and the sex must be the one it gives, before anything runs.
FAM=""
famplan() {  # famplan [run_in options] -- [plan options]: og family plan on the FAMILIES file
  local o=(); while [ "$1" != -- ]; do o+=("$1"); shift; done; shift
  run_in ${o[@]+"${o[@]}"} --cpus 1 --memory 1g -v "${FAMILIES}:/pgp/families.tsv:ro" -v "${PGP_ROOT}/bin:/pgp-bin:ro" \
    "$PYTHON_IMAGE" /pgp-bin/og family plan --families /pgp/families.tsv "$@"
}
if [ -n "${FAMILIES:-}" ]; then
  [ -f "$FAMILIES" ] || { echo "ERROR: FAMILIES=${FAMILIES}: no such file (assets/families.example.tsv shows the format)" >&2; exit 2; }
  FAMILIES="$(cd "$(dirname "$FAMILIES")" && pwd)/$(basename "$FAMILIES")"
  FAM=$(famplan -- --sample "$SAMPLE" --sex "$SEX" --which) || { echo "ERROR: FAMILIES=${FAMILIES}: fix the file or the sex (above)." >&2; exit 2; }
  FAM=$(tail -n 1 <<< "$FAM")
  if [ -n "$FAM" ]; then echo "[Family] ${SAMPLE} is in family ${FAM} (${FAMILIES}): steps 38 and 39 are planned after the pipeline"
  else echo "[Family] ${SAMPLE} is in no family of ${FAMILIES}"; fi
fi
mkdir -p "${S}/nextflow" "$LOG_DIR"
GENOME_DIR="$G" bash "${PGP_ROOT}/bin/write_manifest.sh" "$SAMPLE" run-all.sh "$SEX" || echo "WARNING: could not write ${S}/run_manifest.tsv"
STATUS="${LOG_DIR}/run_status.tsv"  # bin/collect_summary.py marks results older than this run stale
# EVERYTHING=true can still refuse before Nextflow: keep the last run's status to put back then.
[ -z "$EV" ] || [ ! -f "$STATUS" ] || cp "$STATUS" "${STATUS}.prev"
printf '# run-all.sh: when this run started and how each step ended\nmeta\tstarted_epoch\t%s\nmeta\tstarted_utc\t%s\nmeta\tdeclared_sex\t%s\nmeta\tcommand\t%s\nmeta\tswitches\t%s\n' \
  "$(date +%s)" "$(date -u '+%Y-%m-%d %H:%M:%S UTC')" "$SEX" "$RERUN_CMD" "${SWITCHES[*]}" > "$STATUS"
# As run-all.sh ends, however it ends: a step still 'planned' did not finish, and the time is recorded
# (bin/run_state.py reads both). A planned row gives way to the step's result, so each step keeps one row.
finish_status() {
  local rc=$?
  live_post "$SAMPLE" "{\"event\":\"pgp_run_finished\",\"exit\":${rc},\"utcTime\":\"$(date -u +%Y-%m-%dT%H:%M:%SZ)\"}"
  awk -F'\t' -v OFS='\t' 'NR == FNR { if ($1 == "step" && $3 != "planned") ended[$2] = 1; next }
    $1 == "step" && $3 == "planned" { if ($2 in ended) next; $3 = "not finished" } { print }' "$STATUS" "$STATUS" \
    > "${STATUS}.tmp" && mv "${STATUS}.tmp" "$STATUS" || true
  printf 'meta\tfinished_utc\t%s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')" >> "$STATUS" || true
}
trap finish_status EXIT
NF=() SEL=() RUNS=() KNOWN="" OPTIN=()  # every step of a default run: it runs, or it is skipped with the reason
declare -A TOOL_OF=()                   # step number -> its --tools name, for the steps that run
need() { local f; for f in "$@"; do [ -e "$f" ] || { echo "data not installed: ${f#"$G"/}"; return; }; done; }
optin() { [[ ",${T// /}," == *",$1,"* ]] && return 1; echo "opt-in: add $1 to TOOLS"; }  # prints why it is skipped
plan() {  # plan "STEP Label" TOOL [REASON]: 0 when TOOL runs
  local r=${3:-} t=${TOOLS:-}; KNOWN+=" $2"
  [[ -n "$NOBAM" || " 16 16b 10 20 21 04 19 15 22 08 09 09b 18 05 28 32 33 35 35b 37 45 49 " != *" ${1%% *} "* ]] || r="no BAM"
  [ -z "$t" ] || [[ ",${t// /}," == *",$2,"* ]] || r=${r:-not in TOOLS}
  if [ -z "$r" ]; then SEL+=("$2") RUNS+=("${1%% *}") TOOL_OF[${1%% *}]=$2; printf '  %-28s runs\n' "$1"; return 0; fi
  printf '  %-28s skipped    (%s)\n' "$1" "$r"; printf 'step\t%s\tskipped (%s)\n' "${1%% *}" "$r" >> "$STATUS"
  evmiss "$1" "$2" "$r"; return 1
}
# evmiss "STEP Label" TOOL REASON: with EVERYTHING=true, a skip for data that is not installed, or for X=none, stops
# the run before Nextflow (the list prints after the plan). A skip that comes from the input is not one of them.
evmiss() {
  local h=${EVHOW[$2]:-setup.sh}; [ "$h" != - ] || h=setup.sh
  [ -n "$EV" ] || return 0
  case "$3" in
    *=none) EVMISS+=("$1: $3: unset ${3%%=*}") ;;
    "data not installed"*|"not installed"*|"opt-in"*) EVMISS+=("$1: $3; install: $h") ;;
  esac
}
arg() { [ ! -e "$2" ] || NF+=("$1" "$2"); }  # arg --param FILE: pass FILE when it exists
CV="${G}/clinvar/clinvar_pathogenic_chr.vcf.gz" A="${G}/annotations" CAT=${EH_CATALOG:-}
# The catalog of steps 09 and 9b: STRchive's when setup.sh built it, as step 09
# picks it; EH_CATALOG=bundled or no STRchive file: the one inside the image,
# which step 09 reads there and the pipeline takes as a file. Its copy is
# named after the image tag, so an image bump extracts a fresh one.
STRCHIVE_CAT="${G}/reference/strchive-${STRCHIVE_VERSION}/catalog.json"
if [ -z "$CAT" ] && [ -s "$STRCHIVE_CAT" ]; then CAT=$STRCHIVE_CAT; fi
EH_TAG=${EXPANSIONHUNTER_IMAGE##*:} EH_TAG=${EH_TAG//\//_}
if [ -z "$CAT" ] || [ "$CAT" = bundled ]; then CAT="${G}/reference/expansionhunter_variant_catalog-${EH_TAG}.json"; fi
# Each run extracts to its own temporary file, so two runs after an image bump never share one.
[ "${EH_CATALOG:-bundled}" != bundled ] || [ "$CAT" = "$STRCHIVE_CAT" ] || [ -s "$CAT" ] || { CT=$(mktemp "${CAT}.XXXXXX") \
  && { { run_in "$EXPANSIONHUNTER_IMAGE" cat /usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json > "$CT" \
  && [ -s "$CT" ] && chmod 644 "$CT" && mv "$CT" "$CAT"; } || rm -f "$CT"; }; } || :
VEPW=$(need "${G}/vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt") SCORES=()
for s in whole_genome_SNVs.tsv.gz:cadd_snv gnomad.genomes.r4.0.indel.tsv.gz:cadd_indel spliceai_scores.raw.snv.hg38.vcf.gz:spliceai_snv \
         spliceai_scores.raw.indel.hg38.vcf.gz:spliceai_indel revel_grch38.tsv.gz:revel AlphaMissense_hg38.tsv.gz:alphamissense; do
  f="${A}/${s%%:*}"; [ -e "${f}.tbi" ] || f=${f/.raw./.masked.}
  [ ! -e "${f}.tbi" ] || SCORES+=("--${s#*:}" "$f" "--${s#*:}_index" "${f}.tbi")
done
# The samplesheet: the first match below, or the last run's while every file it names exists, its BAM is this call's
# (ALIGN_DIR) and either it starts no earlier in that list or the run that wrote it did not end (so -resume continues it).
# A run ended when main.nf wrote process_status.tsv with no failed, ignored or unfinished task: after it, a FASTQ row
# gives way to the BAM and VCF that run produced.
SHEET="${S}/nextflow/samplesheet.csv" row="" first="" ended="" PT="${S}/nextflow/process_status.tsv"
[ ! -s "$PT" ] || ! awk -F'\t' 'NR > 1 && ($4 + $5 + $6) > 0 {bad = 1} END {exit bad}' "$PT" || ended=1
B="${S}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam" V="${S}/vcf/${SAMPLE}.vcf.gz" F="${S}/fastq/${SAMPLE}_R"
rank() { local f1 b v; IFS=, read -r f1 _ b _ v _ <<< "$1"  # 1 BAM+VCF, 2 BAM, 3 FASTQ, 4 VCF: the order of the first match
  if [ -n "$b" ] && [ -n "$v" ]; then echo 1; elif [ -n "$b" ]; then echo 2; elif [ -n "$f1" ]; then echo 3; else echo 4; fi; }
[ ! -f "$SHEET" ] || row=$(awk -F, 'NR == 2 {print $2","$3","$4","$5","$6","$7}' "$SHEET")
IFS=, read -r -a cols <<< "$row"; for f in "${cols[@]}"; do [ -z "$f" ] || [ -e "$f" ] || row=""; done
[ "${cols[2]:-}" = "$B" ] || [ -z "${cols[2]:-}${ALIGN_DIR:-}" ] || row=""
if [ -f "$B" ] && [ -f "${B}.bai" ] && [ -f "$V" ] && [ -f "${V}.tbi" ]; then first=",,${B},${B}.bai,${V},${V}.tbi"
elif [ -f "$B" ] && [ -f "${B}.bai" ]; then first=",,${B},${B}.bai,,"
elif [ -f "${F}1.fastq.gz" ] && [ -f "${F}2.fastq.gz" ]; then first="${F}1.fastq.gz,${F}2.fastq.gz,,,,"
elif [ -f "$V" ] && [ -f "${V}.tbi" ]; then first=",,,,${V},${V}.tbi"; fi
[ -n "$row" ] && { [ -z "$first" ] || [ -z "$ended" ] || [ "$(rank "$row")" -le "$(rank "$first")" ]; } || row=$first
[ -n "$row" ] || { echo "ERROR: no input for ${SAMPLE}: no ${B}, no ${F}1/2.fastq.gz and no ${V}, each with its index." >&2; exit 1; }
IFS=, read -r f1 _ fb _ fv _ <<< "$row"; NOBAM=${f1}${fb} GV="${S}/vcf/${SAMPLE}.g.vcf.gz" GCOL="" GROW=""
# A BAM+VCF row is not called again: the gVCF beside the VCF, when there is one with its index, gives PharmCAT and PRS the reference calls.
if [ -n "$fb" ] && [ -n "$fv" ]; then
  if [ -f "$GV" ] && [ -f "${GV}.tbi" ]; then GCOL=",gvcf,gvcf_index" GROW=",${GV},${GV}.tbi"; echo "NOTE: starting from the existing VCF; PharmCAT and PRS read the gVCF beside it (${GV})."
  else echo "NOTE: starting from the existing VCF, with no gVCF with its index beside it (${GV}): PharmCAT and PRS read its variant sites only. To call again with a gVCF, remove the VCF and its index."; fi
fi
printf 'sample,fastq_1,fastq_2,bam,bam_index,vcf,vcf_index,sex%s\n%s,%s,%s%s\n' "$GCOL" "$SAMPLE" "$row" "$SEX" "$GROW" > "$SHEET"
echo "[Input] ${SHEET}: ${row}"
if [ -n "$NOBAM" ]; then RUNS=(16); else printf 'step\t16\tskipped (no BAM)\n' >> "$STATUS"; fi
echo "[Steps] ${SAMPLE} (${SEX})"
for p in "07 PharmCAT:pharmcat" "27 CPIC lookup:cpic" "11 ROH:roh" "12 Mito haplogroup:mito_haplogroup" "16b mosdepth:mosdepth" \
         "10 TelomereHunter:telomere_hunter" "20 Mito variants (Mutect2):mito_variants" "04 Manta:manta" \
         "19 Delly:delly" "15 duphold:duphold" "22 SV consensus merge:survivor_merge" "28 MultiQC:multiqc"; do plan "${p%:*}" "${p##*:}" || true; done
plan "06 ClinVar screen" clinvar "$(need "$CV" "${CV}.tbi")" && NF+=(--clinvar "$CV" --clinvar_index "${CV}.tbi")
H=$(data_file hla_dat || true) GT=$(data_file gencode_genes || true)
plan "08 HLA typing (T1K)" hla_typing "$(need "$H" "$GT")" && NF+=(--hla_dat "$H" --hla_genes "$GT") \
  && if on KIR; then K="${G}/kir/IPD-KIR_${KIR_DB_RELEASE}/kir.dat"; if [ -s "$K" ]; then NF+=(--kir true --kir_dat "$K")
    else echo "  KIR typing skipped: data not installed: ${K#"$G"/} (setup.sh --kir-data)"
      [ -z "$EV" ] || EVMISS+=("08 KIR typing (KIR=true): data not installed: ${K#"$G"/}; install: setup.sh --kir-data"); fi; fi
[[ " ${SEL[*]} " != *" hla_typing "* ]] || ! on HLA_QC || { NF+=(--hla_qc true); echo "    (HLA_QC=true: OptiType checks T1K's HLA-A and HLA-B)"; }
plan "09 ExpansionHunter" expansion_hunter "$(need "$CAT")" && NF+=(--expansion_catalog "$CAT")
plan "09b Stranger" stranger "$(need "$CAT")" || true
plan "13 VEP" vep "$VEPW" && NF+=(--vep_cache "${G}/vep_cache") && arg --gnomad_constraint "${A}/gnomad_v4.1_constraint.tsv"
NOSCORE=$([ "${#SCORES[@]}" -gt 0 ] || echo 'data not installed: no score file in annotations/')
plan "30 vcfanno" vcfanno "${VEPW:+needs VEP, }${VEPW:-$NOSCORE}" && NF+=("${SCORES[@]}")
plan "23 Clinical filter" clinical_filter "${VEPW:+needs VEP, ${VEPW}}" || true
plan "31 slivar" slivar "${VEPW:+needs VEP, ${VEPW}}" || true
plan "17 CPSR" cpsr "$(need "${G}/vep_cache/homo_sapiens/${PCGR_VEP_CACHE_RELEASE}_GRCh38/info.txt" "${G}/pcgr_data/${PCGR_DATA_BUNDLE}/data")" \
  && NF+=(--pcgr_data "${G}/pcgr_data/${PCGR_DATA_BUNDLE}" --vep_cache_cpsr "${G}/vep_cache")
[[ " ${SEL[*]} " != *" cpsr "* ]] || [ "$SF" = true ] \
  || { NF+=(--cpsr_secondary_findings false); echo "    (CPSR_SECONDARY_FINDINGS=false: CPSR reports the cancer genes only)"; }
plan "18 CNVpytor" cnvpytor "$(need "${G}/reference/cnvpytor/gc_hg38.pytor")" && NF+=(--cnvpytor_resources "${G}/reference/cnvpytor")
plan "05 AnnotSV" annotsv "$(need "${G}/annotsv_annotations/Annotations_Human/Genes/GRCh38")" && NF+=(--annotsv_annotations "${G}/annotsv_annotations")
plan "32 pypgx" pypgx "$(need "${G}/reference/pypgx-bundle")" && NF+=(--pypgx_bundle "${G}/reference/pypgx-bundle")
# Step 33 (sample identity, sex from the reads, nuclear contamination) with the two files setup.sh installs.
SQ="${G}/reference/somalier/sites.hg38.vcf.gz" VB="${G}/reference/verifybamid2" VP=1000g.phase3.100k.b38.vcf.gz.dat
plan "33 Sample identity and QC" sample_qc "$(need "$SQ" "${VB}/${VP}.UD" "${VB}/${VP}.mu" "${VB}/${VP}.bed")" \
  && NF+=(--somalier_sites "$SQ" --verifybamid2_panel "$VB")
CY="${G}/tools/cyrius-${CYRIUS_VERSION}" PS="${G}/reference/parascopy-${PARASCOPY_DATA_VERSION}"
CYSTAMP="python=${PYTHON_IMAGE} lock=$(_digest sha256 "${SCRIPT_DIR}/cyrius-constraints.txt")"
# Step 21 runs when setup.sh --cyrius installed this version (whoever ran it accepted Cyrius's non-commercial licence).
CYOK=$([ "$(cat "${CY}/INSTALLED" 2>/dev/null)" != "$CYSTAMP" ] || echo yes)
plan "21 Cyrius CYP2D6" cyrius "$([ -n "$CYOK" ] || if [ -e "${CY}/INSTALLED" ] || [[ ",${T// /}," == *",cyrius,"* ]]; then
  echo "data not installed: ${CY#"$G"/} for this version (setup.sh --cyrius)"; else echo "not installed: setup.sh --cyrius (non-commercial licence)"; fi)" \
  && NF+=(--cyrius_install "$CY") && { [ -n "$T" ] || echo "    (planned because setup.sh --cyrius installed it; TOOLS without cyrius leaves it out)"; }
plan "35 Parascopy SMN1/SMN2" parascopy "$(optin parascopy || need "${PS}/homology_table/GRCh38.bed.gz")" && NF+=(--parascopy_data "$PS") \
  && NF+=(--parascopy_population "${PARASCOPY_POPULATION:-EUR}")
# Step 35b reads PacBio HiFi reads only. og paraphase check-bam, which PARAPHASE_CHECK runs too, says whether every
# read group of the BAM is PL:PACBIO (exit 3: not). Another BAM is a skip, not a failed step; so is a FASTQ row, which
# the pipeline aligns as short reads. A check that cannot run leaves the step to PARAPHASE_CHECK.
hifi() {
  [ -n "$NOBAM" ] || return 0
  [ -n "$fb" ] || { echo 'not PacBio HiFi: a FASTQ row is aligned as short reads'; return 0; }
  local rc=0
  run_in --cpus 1 --memory 1g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYTHON_IMAGE" /pgp-bin/og paraphase check-bam "/genome/${fb#"$G"/}" \
    > /dev/null 2>&1 || rc=$?
  [ "$rc" -ne 3 ] || echo 'not PacBio HiFi: a read group of the BAM is not PL:PACBIO'
}
plan "35b Paraphase (PacBio HiFi)" paraphase "$(optin paraphase || hifi)" && { [ -z "${PARAPHASE_REGIONS:-}" ] || NF+=(--paraphase_regions "$PARAPHASE_REGIONS"); }
plan "25 PRS" prs "$(need "$(compgen -G "${G}/prs_scores/*.txt.gz" | head -n 1 || echo "${G}/prs_scores/<PGS id>.txt.gz")")" && NF+=(--pgs_scoring "${G}/prs_scores") \
  && { PC=${PGSC_CALC_DIR:-${G}/tools/pgsc_calc-${PGSC_CALC_VERSION}}; [ ! -f "${PC}/main.nf" ] || NF+=(--pgsc_calc "$PC"); }
# Step 26 is pgsc_calc's projection inside the PRS run, with its panel and the site list setup.sh writes beside it.
PANEL=${ANCESTRY_PANEL:-${G}/reference/pgsc_calc/${PGSC_PANEL}.tar.zst}
plan "26 Ancestry (pgsc_calc)" ancestry "$(if [[ " ${SEL[*]} " != *" prs "* ]]; then echo 'needs PRS'; elif [ "$PANEL" = none ]; then echo 'ANCESTRY_PANEL=none'
  else need "$PANEL" "${PANEL%.tar.zst}_GRCh38_sites.tsv"; fi)" && NF+=(--ancestry_ref "$PANEL")
# Step 26b against the panel setup.sh --ancestry-fine-panel installs (model/DONE is written last).
AF=${ANCESTRY_FINE_DIR:-${G}/reference/ancestry_fine/${ANCESTRY_FINE_PANEL}}
plan "26b Fine-grained ancestry" ancestry_fine "$(if [ "$AF" = none ]; then echo 'ANCESTRY_FINE_DIR=none'; else need "${AF}/model/DONE"; fi)" \
  && NF+=(--ancestry_fine_panel "$AF")
YD="${G}/reference/yleaf-${YLEAF_DATA_VERSION}/data"
plan "37 Y haplogroup (Yleaf)" y_haplogroup "$(optin y_haplogroup || if [ "$SEX" = male ]; then need "${YD}/hg38/new_positions.txt"
  else echo 'female sample: no Y chromosome'; fi)" && NF+=(--yleaf_data "$YD")
# Step 40 reads the gVCF: a row that starts from a VCF without one beside it has none.
CS="${G}/reference/consumer/23andme_v5_sites_dbsnp${CONSUMER_DBSNP_BUILD}.tsv.gz"
plan "40 Raw-data export" consumer_export "$(optin consumer_export || if [ -n "$fv" ] && [ -z "$GROW" ]; then echo 'no gVCF beside the VCF'; else need "$CS"; fi)" \
  && NF+=(--consumer_sites "$CS")
# Step 41 (blood groups) is opt-in and listed only when TOOLS names it; without a BAM its RhD is not determined.
if [[ ",${T// /}," == *",bloodtype,"* ]]; then plan "41 Blood groups" bloodtype || true; else KNOWN+=" bloodtype"; fi
# Step 46 (Lp(a), the LPA KIV-2 count) is opt-in and listed only when TOOLS names it; it reads the BAM.
if [[ ",${T// /}," == *",lpa,"* ]]; then plan "46 Lp(a) KIV-2 count" lpa "$([ -n "$NOBAM" ] || echo 'no BAM')" || true; else KNOWN+=" lpa"; fi
# Step 52 (MUC1 VNTR frameshifts, VNtyper) is opt-in and listed only when TOOLS names it; it reads the BAM.
if [[ ",${T// /}," == *",vntyper,"* ]]; then plan "52 MUC1 VNTR (VNtyper)" vntyper "$([ -n "$NOBAM" ] || echo 'no BAM')" || true; else KNOWN+=" vntyper"; fi
# Step 49 (CHIP) is opt-in and listed only when TOOLS names it: Mutect2 on the BAM, with step 29's gnomAD, panel of normals and common sites when installed.
CG="${G}/reference/ensembl/Homo_sapiens.GRCh38.${ENSEMBL_GFF_RELEASE}.chr.gff3.gz"
if [[ ",${T// /}," == *",chip,"* ]]; then plan "49 CHIP (clonal blood cells)" chip "$(need "$CG")" && NF+=(--chip_gff "$CG") \
  && for r in chip_germline_resource:af-only-gnomad chip_pon:1000g_pon chip_common_sites:small_exac_common_3; do f="${G}/somatic/${r#*:}.hg38.vcf.gz"; [ ! -e "${f}.tbi" ] || arg "--${r%%:*}" "$f"; done
else KNOWN+=" chip"; fi
# Step 42 reads the gVCF only: DeepVariant writes one for a FASTQ or BAM row, a BAM+VCF row has one when it lies beside the VCF.
plan "42 Traits panel" traits "$([ -z "$fv" ] || [ -n "$GROW" ] || echo "no gVCF: ${GV#"$G"/}")" || true
# Step 43 against the table setup.sh --archaic builds (DONE is written last); the genes of its segments come from GENCODE.
# It reads the gVCF too (a site with no call is the reference only inside a reference block): with the table installed,
# a row that starts from a VCF without one beside it is skipped, as for step 42.
AR=${ARCHAIC_DIR:-${G}/reference/archaic/${ARCHAIC_DATA_VERSION}}
plan "43 Archaic ancestry (IBDmix)" archaic "$(if [ "$AR" = none ]; then echo 'ARCHAIC_DIR=none'
  elif [ -e "${AR}/DONE" ] && [ -n "$fv" ] && [ -z "$GROW" ]; then echo "no gVCF: ${GV#"$G"/}"; else need "${AR}/DONE"; fi)" \
  && NF+=(--archaic_data "$AR" --archaic_pop "${ARCHAIC_POP:-ALL}") && { [[ " ${NF[*]} " == *" --hla_genes "* ]] || [ ! -s "$GT" ] || NF+=(--hla_genes "$GT"); }
# Step 44 (ancient populations) is opt-in and listed only when TOOLS names it; it reads the gVCF and the panel
# setup.sh --ancient-panel builds (DONE is written last). ANCIENT_PANEL_DIR picks another panel.
AN=${ANCIENT_PANEL_DIR:-${G}/reference/ancient/aadr_${AADR_VERSION}}
if [[ ",${T// /}," == *",ancient,"* ]]; then
  plan "44 Ancient populations" ancient "$(if [ -n "$fv" ] && [ -z "$GROW" ]; then echo "no gVCF: ${GV#"$G"/}"; else need "${AN}/DONE"; fi)" \
    && NF+=(--ancient_panel "$AN")
else KNOWN+=" ancient"; fi
# Step 45 (Y-STR repeat counts) is opt-in and listed only when TOOLS names it; it reads chrY, so a female sample is skipped.
if [[ ",${T// /}," == *",y_str,"* ]]; then plan "45 Y-STR markers (HipSTR)" y_str "$([ "$SEX" = male ] || echo 'female sample: no Y chromosome')" || true
else KNOWN+=" y_str"; fi
# Step 47 (viral DNA in the unmapped reads) is opt-in and listed only when TOOLS names it; it reads the BAM.
VR=${VIRAL_REF:-${G}/reference/viral-${VIRAL_REF_VERSION}/viral_genomes.fa}
if [[ ",${T// /}," == *",viral,"* ]]; then plan "47 Viral DNA" viral "$(if [ -z "$NOBAM" ]; then echo 'no BAM'; else need "$VR"; fi)" \
  && NF+=(--viral_ref "$VR") || true; else KNOWN+=" viral"; fi
# Step 48 (Exomiser) runs when HPO names the symptoms, against the data setup.sh --exomiser-data installs.
ED="${G}/reference/exomiser-${EXOMISER_DATA_VERSION}"
if [ -n "${HPO:-}" ] || [[ ",${T// /}," == *",exomiser,"* ]]; then
  # A set HPO turns step 48 on even when TOOLS leaves exomiser out, as --hpo does in main.nf.
  if TOOLS=$([ -n "${HPO:-}" ] && echo exomiser || echo "${TOOLS:-}") plan "48 Exomiser (HPO terms)" exomiser "$([ -n "${HPO:-}" ] || echo 'HPO not set: HPO=HP:0001250,...'; [ -z "${HPO:-}" ] \
       || need "${ED}/${EXOMISER_DATA_VERSION}_hg38" "${ED}/${EXOMISER_DATA_VERSION}_phenotype")"; then
    NF+=(--hpo "$HPO" --exomiser_data "$ED")
  fi
else KNOWN+=" exomiser"; fi
# Step 50 (Alu, L1 and SVA insertions, MEGAnE) is opt-in and listed only when TOOLS names it; it reads the BAM.
MN=$(basename "$REF_FASTA"); MN=${MN%.gz}; MN=${MN%.*}; MK=${MEGANE_KMERS:-${G}/reference/megane/${MN}.mk}
case "$MK" in /*) ;; *) MK="${G}/${MK}" ;; esac
if [[ ",${T// /}," == *",mei,"* ]]; then plan "50 Mobile element insertions" mei "$(if [ -z "$NOBAM" ]; then echo 'no BAM'; else need "$MK" "${MK%.mk}.mi"; fi)" \
  && NF+=(--megane_kmers "$MK") || true; else KNOWN+=" mei"; fi
# Step 51 (mosaic chromosomal alterations, MoChA) is opt-in and listed only when TOOLS names it; it reads the VCF.
MD=${MOCHA_DATA:-${G}/reference/mocha-${MOCHA_PANEL_VERSION}}
if [[ ",${T// /}," == *",mocha,"* ]]; then plan "51 Mosaic alterations (MoChA)" mocha "$(need "${MD}/DONE")" \
  && NF+=(--mocha_data "$MD") || true; else KNOWN+=" mocha"; fi
# Step 53 (a second opinion at hard polymorphic genes, Locityper) is opt-in and listed only when TOOLS names it; it reads
# the BAM and the database setup.sh --locityper-data builds (DONE is written last). LOCITYPER_BG_REGION picks preproc's region.
LD=${LOCITYPER_DATA:-${G}/reference/locityper-${LOCITYPER_PANEL_VERSION}}
case "$LD" in /*) ;; *) LD="${G}/${LD}" ;; esac
if [[ ",${T// /}," == *",locityper,"* ]]; then plan "53 Hard loci (Locityper)" locityper "$(if [ -z "$NOBAM" ]; then echo 'no BAM'; else need "${LD}/DONE"; fi)" \
  && NF+=(--locityper_data "$LD") && { [ -z "${LOCITYPER_BG_REGION:-}" ] || NF+=(--locityper_bg_region "$LOCITYPER_BG_REGION"); } || true; else KNOWN+=" locityper"; fi
for t in ${T//,/ }; do [[ "${KNOWN} " == *" ${t} "* ]] || { echo "ERROR: unknown step '${t}' in TOOLS. Known:${KNOWN}" >&2; exit 2; }; done
# Step 4b (GRIDSS=true) runs after the pipeline, and GRIDSS needs the classic BWA index beside the reference, which
# setup.sh does not build: check it now, so the run stops here and not after hours of pipeline.
if [ -n "$EV" ] && [ -n "$NOBAM" ]; then for e in amb ann bwt pac sa; do [ -f "${REF_FASTA}.${e}" ] || {
  EVMISS+=("04b GRIDSS (GRIDSS=true): classic BWA index not installed: ${REF_FASTA#"$G"/}.${e}; install: bwa index (docs/04b-gridss.md)"); break; }; done; fi
if [ "${#EVMISS[@]}" -gt 0 ]; then
  { echo "ERROR: EVERYTHING=true runs every optional step, and ${#EVMISS[@]} cannot run (Nextflow did not start):"; printf '  %s\n' "${EVMISS[@]}"
    echo "Install what is missing, or run without EVERYTHING=true to skip them."; } >&2
  # Nothing ran: the last run's status stays as it was (or there is none), and no run-finished event goes out.
  trap - EXIT; if [ -f "${STATUS}.prev" ]; then mv "${STATUS}.prev" "$STATUS"; else rm -f "$STATUS"; fi; exit 2
fi
rm -f "${STATUS}.prev"
arg --cytoband "$(data_file cytoband || true)"; ! GEX=$(data_file gencode_exons) || NF+=(--gencode_exons "$GEX"); arg --delly_exclude "$(data_file delly_exclude || true)"; arg --manta_call_regions "${MANTA_CALL_REGIONS:-}"
# Step 20's annotation of the chrM calls (MITO_VARIANTS_ANNOTATE): ClinVar with its index, and gnomAD's chrM sites.
MC="${G}/clinvar/clinvar_chr.vcf.gz"; if [ -s "$MC" ] && [ -s "${MC}.tbi" ]; then NF+=(--mito_clinvar "$MC"); fi; arg --gnomad_mito "${G}/gnomad_mito/${GNOMAD_MITO_URL##*/}"
[ -z "${INTERVALS:-}" ] || NF+=(--intervals "$INTERVALS")
# The CPUs and RAM the tasks get: the host's, or the container engine's when docker info reports less. Docker Desktop
# and colima run the containers in a VM, often smaller than the host; Nextflow books the host's RAM by default, and
# a task the VM cannot hold is killed (exit 137). Then the default caps are the engine's, and executor.config sets the
# pool Nextflow fills to it, so a task waits for room instead. Each executor setting in a later -c of the user's wins.
NC=$(getconf _NPROCESSORS_ONLN) EM="" EC="" EXC="${S}/nextflow/executor.config" EX=""
M=$(awk '/^MemTotal:/ {print int($2 / 1048576)}' /proc/meminfo 2>/dev/null || echo $(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1073741824 )))
read -r EM EC _ < <(docker info --format '{{.MemTotal}} {{.NCPU}}' 2>/dev/null) || true
if [[ "$EM" =~ ^[0-9]+$ ]] && [ $((EM / 1073741824)) -ge 1 ] && { [ "${M:-0}" -lt 1 ] || [ $((EM / 1073741824)) -lt "$M" ]; }; then M=$((EM / 1073741824)) EX=1; fi
if [[ "$EC" =~ ^[0-9]+$ ]] && [ "$EC" -ge 1 ] && [ "$EC" -lt "$NC" ]; then NC=$EC EX=1; fi
[[ " $* " == *" --max_cpus"* ]] || NF+=(--max_cpus "${USER_THREADS:-$NC}")
[[ " $* " == *" --max_memory"* ]] || [ "${M:-0}" -lt 1 ] || NF+=(--max_memory "${M}.GB")
if [ -n "$EX" ]; then
  # The pool never goes below a cap the user set (THREADS, --max_cpus, --max_memory): Nextflow refuses a task that asks
  # for more than the pool. A cap it cannot read (not whole CPUs, not whole GB) leaves that line out, as on the host.
  XC=$NC XM=$M UC=${USER_THREADS:-} UM="" prev="" GB_RE='^([0-9]+)[[:space:]]*\.?[[:space:]]*[Gg][Bb]?$'
  for a in "$@"; do
    case "$prev" in --max_cpus) UC=$a ;; --max_memory) UM=$a ;; esac
    case "$a" in --max_cpus=*) UC=${a#*=} ;; --max_memory=*) UM=${a#*=} ;; esac; prev=$a
  done
  if [ -n "$UC" ]; then if [[ "$UC" =~ ^[0-9]+$ ]]; then [ "$UC" -le "$XC" ] || XC=$UC; else XC=""; fi; fi
  if [ -n "$UM" ]; then if [[ "$UM" =~ $GB_RE ]]; then [ "${BASH_REMATCH[1]}" -le "${XM:-0}" ] || XM=${BASH_REMATCH[1]}; else XM=""; fi; fi
  { printf 'executor {\n'; [ -z "$XC" ] || printf '    cpus   = %s\n' "$XC"; [ "${XM:-0}" -lt 1 ] || printf "    memory = '%s GB'\n" "$XM"; printf '}\n'; } > "$EXC"
  XCS="${XC} CPUs" XMS="${XM} GB"; [ -n "$XC" ] || XCS="the host's CPUs"; [ "${XM:-0}" -ge 1 ] || XMS="the host's RAM"
  NF+=(-c "$EXC"); echo "[Resources] the container engine has less than this host: Nextflow books at most ${XCS} and ${XMS} (${EXC})"
  [ "$XC" = "$NC" ] && [ "$XM" = "$M" ] || echo "[Resources] the pool follows a task cap you set above the engine's ${NC} CPUs and ${M} GB, or one not in whole CPUs or GB: a task can then be killed inside the engine (exit 137)"
else rm -f "$EXC"; fi
# Publish by hard link when it is sure to work: the BAM and every other output are then stored once, not twice.
# Nextflow does not fall back to a copy when a link fails; the run stops. The tasks run as root (the docker profile)
# and their files can stay root's, and Linux lets another user link a file it does not own only when
# fs.protected_hardlinks is 0. So: root or protected_hardlinks 0, and a probe link from work/ to GENOME_DIR that works.
# The probe only sees the default work/: with -w, -work-dir, NXF_WORK, or a -c/-config that may set workDir, keep the copy.
if [[ " $* " != *" --publish_dir_mode"* && " $* " != *" -w "* && " $* " != *" -work-dir"* && " $* " != *" -c "* \
   && " $* " != *" -config "* ]] && [ -z "${NXF_WORK:-}" ] \
   && { [ "$(id -u)" -eq 0 ] || [ "$(cat /proc/sys/fs/protected_hardlinks 2>/dev/null)" = 0 ]; }; then
  LP="${S}/nextflow/work/.pgp-link-probe.$$"
  mkdir -p "${S}/nextflow/work" && : > "$LP" && ln "$LP" "${S}/.pgp-link-probe.$$" 2>/dev/null && ln "$LP" "${G}/.pgp-link-probe.$$" 2>/dev/null \
    && NF+=(--publish_dir_mode link)
  rm -f "$LP" "${S}/.pgp-link-probe.$$" "${G}/.pgp-link-probe.$$"
fi
if on SKIP_TRIM; then NF+=(--skip_trim true); fi
for v in TIDDIT:04a-tiddit GRIDSS:04b-gridss IMPUTATION:14-imputation-prep SOMATIC:29-mutect2-somatic; do if on "${v%%:*}"; then OPTIN+=("${v#*:}.sh"); fi; done
for c in ${C//,/ }; do OPTIN+=("$(cd "$SCRIPT_DIR" && compgen -G "03[a-d]-${c}*.sh")") || { echo "ERROR: unknown caller '${c}' in EXTRA_CALLERS (gatk, freebayes, strelka2, octopus)" >&2; exit 2; }; done
if ! on BENCHMARK; then :; elif [ -n "$C" ] || compgen -G "${S}/vcf_*/${SAMPLE}.vcf.gz" >/dev/null || [ -f "${S}/vcf_strelka2/results/variants/variants.vcf.gz" ]; then OPTIN+=(benchmark-variants.sh)
else echo "  benchmark-variants skipped (only one caller VCF: set EXTRA_CALLERS, or run a 03a-03d script first)"; fi
[ -z "${MAX_JOBS:-}" ] || echo "NOTE: MAX_JOBS is no longer read. --max_cpus (THREADS) and --max_memory cap each task; Nextflow fills the machine's CPUs and RAM."
[ -z "${ANCESTRY:-}" ] || echo "NOTE: ANCESTRY is no longer read: step 26 runs in the pipeline whenever the ancestry panel is installed (setup.sh --ancestry-panel)."
export NXF_VER=${NXF_VER:-$NEXTFLOW_VERSION}
# The live site (docs/live-site.md). Nextflow's weblog is the nf-weblog plugin, which Nextflow 26.10 loads only when
# a config names it, so the block lives in a config of this run's own, passed only when OWNGENOME_LIVE_URL is set: a run
# without it neither downloads the plugin nor warns about an unknown weblog scope. The pin is NF_WEBLOG_PLUGIN in versions.env.
if [ -n "${OWNGENOME_LIVE_URL:-}" ]; then
  # The ingest listener speaks plain HTTP and serves /ingest/<run> at its root, so https or a path would lose every post.
  [[ "$OWNGENOME_LIVE_URL" =~ ^http://[A-Za-z0-9.:_-]+/?$ ]] \
    || { echo "ERROR: OWNGENOME_LIVE_URL must be a plain http://host:port URL, got '${OWNGENOME_LIVE_URL}'" >&2; OWNGENOME_LIVE_URL=""; exit 2; }
  export OWNGENOME_LIVE_RUN=$SAMPLE
  printf "plugins {\n    id '%s'\n}\nweblog {\n    enabled = true\n    url     = '%s'\n}\n" "$NF_WEBLOG_PLUGIN" \
    "${OWNGENOME_LIVE_URL%/}/ingest/${SAMPLE}" > "${S}/nextflow/live.config"
  NF+=(-c "${S}/nextflow/live.config")
  echo "Live site: messages go to ${OWNGENOME_LIVE_URL%/}/ingest/${SAMPLE} (${NF_WEBLOG_PLUGIN}, downloaded by Nextflow on first use)"
  # The address the browser opens (the site's own, not the ingest port), when the launcher knows it.
  [ -z "${OWNGENOME_LIVE_PAGE_URL:-}" ] || echo "Live site: this run's page is ${OWNGENOME_LIVE_PAGE_URL%/}/runs/${SAMPLE}/"
else
  rm -f "${S}/nextflow/live.config"
fi
live_step() { live_post "$SAMPLE" "{\"event\":\"pgp_step\",\"step\":\"$1\",\"state\":\"$2\",\"utcTime\":\"$(date -u +%Y-%m-%dT%H:%M:%SZ)\"}"; }
echo "[Nextflow ${NXF_VER}] launch directory ${S}/nextflow: .nextflow.log, and every task's files under work/"
FT="${S}/nextflow/failed_tasks.tsv" PT="${S}/nextflow/process_status.tsv"; rm -f "$FT" "$PT"  # main.nf writes both as the run ends
# Every step that will run, before it starts: 'planned' until its result row (bin/run_state.py shows the run as it goes).
for n in "${RUNS[@]}" "${OPTIN[@]}"; do n=${n%%-*}; [[ "$n" != [0-9]* ]] || printf 'step\t%s\tplanned\n' "$n" >> "$STATUS"; done
# An earlier run's OwnHealth file goes before any result changes, so a run cut short never leaves it beside new results.
rm -f "${S}/facts/${SAMPLE}_genome_facts.json"
rc=0; (cd "${S}/nextflow" && nextflow run "${PGP_ROOT}/main.nf" -profile docker -resume --input "$SHEET" --reference "$REF_FASTA" \
  --outdir "$G" --tools "$(IFS=,; echo "${SEL[*]}")" "${NF[@]}" "$@") || rc=$?
# A step's result is what Nextflow says of its tool's processes in this run. 'failed': a task failed (a report-only
# tool's failure is ignored by the pipeline, which goes on; any other stops it). 'ok': the pipeline ended well, or it
# stopped and process_status.tsv shows every process of the tool with a task that completed or came from the cache
# (-resume: the same inputs and code, so the result is this run's even when its file is older) and none left
# unfinished. 'not finished': the pipeline stopped before that. Stranger reads ExpansionHunter's calls, so it does not
# run when ExpansionHunter failed.
# tool_of PROCESS: the selected tool the process belongs to, by name: the tool's own (PRS: prs) or the tool's plus _<part>
# (PRS_SUMMARY, CPIC_LOOKUP: cpic). Y_STR is step 45's, Yleaf's three are the other Y_*, the SV merge's SURVIVOR_*, the mito haplogroup's chrM
# extraction and haplocheck are step 12's, T1K's index build, the KIR pass (KIR=true) and OptiType (HLA_QC=true) are step 08's, somalier's and
# VerifyBamID2's are step 33's (sample_qc); indexcov is step 16.
tool_of() {
  local p=${1##*:} t best=""; p=${p,,}
  case "$p" in y_str) echo y_str; return ;; y_*) echo y_haplogroup; return ;; survivor_*) echo survivor_merge; return ;;
    mito_extract_chrm|haplocheck) echo mito_haplogroup; return ;; t1k_build|kir_build|kir_typing|hla_qc|hla_qc_reads) echo hla_typing; return ;;
    somalier|somalier_relate|verifybamid2) echo sample_qc; return ;; esac
  for t in "${SEL[@]}" indexcov; do if [[ "$p" == "$t" || "$p" == "${t}_"* ]] && [ "${#t}" -gt "${#best}" ]; then best=$t; fi; done
  echo "${best:-$p}"
}
declare -A BAD=() FIN=() OPEN=(); failed=0
while IFS=$'\t' read -r p st; do
  [ "$p" != process ] && [ -n "$p" ] || continue
  t=$(tool_of "$p") p=${p##*:}
  [[ " ${SEL[*]} " == *" ${t} "* ]] || { echo "  ${p}: ${st}" >&2; continue; }
  BAD[$t]=1 failed=1 rest=$(IFS=,; echo "${SEL[*]}"); rest=",${rest},"; rest=${rest/,${t},/,}; rest=${rest#,}
  echo "  ${t} (${p}) failed$([ "$st" != ignored ] || echo ', skipped'): see ${S}/nextflow/.nextflow.log. To run without it: TOOLS=${rest%,}" >&2
done < <(cat "$FT" 2>/dev/null)
# A process with no task had not started, except haplocheck without step 20's Mutect2 calls: it never gets one then.
while IFS=$'\t' read -r p st; do
  t=$(tool_of "$p"); case "$st" in bad) BAD[$t]=1 ;; fin) FIN[$t]=1 ;; open) OPEN[$t]=1 ;;
    *) [[ "${p##*:}" == HAPLOCHECK && ( " ${SEL[*]} " != *" mito_variants "* || -n "${BAD[mito_variants]:-}" ) ]] || OPEN[$t]=1 ;; esac
done < <(awk -F'\t' 'NR > 1 && NF == 6 {print $1 "\t" ($4 + $5 > 0 ? "bad" : $6 > 0 ? "open" : $2 + $3 > 0 ? "fin" : "none")}' "$PT" 2>/dev/null)
for n in "${RUNS[@]}"; do
  t=${TOOL_OF[$n]:-indexcov}; [ "$t" != ancestry ] || t="prs"  # step 26 is computed inside the PRS run
  r="not finished"; if [ "$rc" -eq 0 ] || { [ -n "${FIN[$t]:-}" ] && [ -z "${OPEN[$t]:-}" ]; }; then r=ok; fi
  [ -z "${BAD[$t]:-}" ] || r=failed
  [ "$t" != stranger ] || [ -z "${BAD[expansion_hunter]:-}" ] || r="skipped (expansion_hunter failed)"
  printf 'step\t%s\t%s\n' "$n" "$r" >> "$STATUS"
done
# The reports are written after a failed pipeline too: they show what finished and mark older results stale.
[ "$rc" -eq 0 ] || { echo "ERROR: the pipeline failed (exit ${rc}): see above and ${S}/nextflow/.nextflow.log. Fix the cause and run the same command; -resume reruns only what did not finish. Writing the reports of what finished." >&2; OPTIN=(); }
# Steps 38 and 39 for the sample's family, after the pipeline and the opt-in scripts and before the reports (so they
# show this run's couple table), as og family plan decides from what is under GENOME_DIR now: each member's
# gVCF (with its index) and step 06's ClinVar hits, and step 26b's panel for the shared segments.
family_steps() {
  [ -n "$FAM" ] && [ "$rc" -eq 0 ] || return 0
  HG=() HC=()
  for f in "$G"/*/vcf/*.g.vcf.gz; do m=${f#"$G"/}; m=${m%%/*}; [ "$f" != "${G}/${m}/vcf/${m}.g.vcf.gz" ] || [ ! -f "${f}.tbi" ] || HG+=("$m"); done
  for f in "$G"/*/clinvar/*_clinvar_hits.vcf; do m=${f#"$G"/}; m=${m%%/*}; [ "$f" != "${G}/${m}/clinvar/${m}_clinvar_hits.vcf" ] || [ ! -f "$f" ] || HC+=("$m"); done
  IBD=no; [ "$AF" = none ] || { [ -s "${AF}/panel/DONE" ] && [ -d "${AF}/maps" ] && IBD=yes; } || true
  FP="${G}/${FAM}/family_plan.tsv"; mkdir -p "${G}/${FAM}"
  echo "[Family] ${FAM}: ${FP}"
  famplan --rw "${G}/${FAM}" -- --family "$FAM" --have-gvcf ${HG[@]+"${HG[@]}"} --have-clinvar ${HC[@]+"${HC[@]}"} \
    --ibd-panel "$IBD" --out "/genome/${FAM}/family_plan.tsv" || {
    # never act on an older family_plan.tsv: the family steps wait for a plan that worked
    echo "  FAILED: og family plan; steps 38 and 39 not run" >&2; failed=1
    printf 'step\t%s\tskipped (og family plan failed)\n' 38 39 >> "$STATUS"; return 0; }
  # decision<TAB>reason of step 38 (the joint call) and step 39 (a couple that runs, else the first couple)
  d38=$(awk -F'\t' '$2 == "joint_call" {print $4 "\t" $5; exit}' "$FP" 2>/dev/null)
  d39=$(awk -F'\t' '$2 == "couple_carrier" {r = r ? r : $4 "\t" $5; if ($4 == "runs") {r = $4 "\t" $5; exit}} END {print r}' "$FP" 2>/dev/null)
  famstep() {  # famstep STEP SCRIPT: run the family step's script on the family, record ok or failed
    local r=ok; printf 'step\t%s\tplanned\n' "$1" >> "$STATUS"; echo "[${2%.sh}] ${FAM}, log: ${LOG_DIR}/${2%.sh}.log"
    bash "${SCRIPT_DIR}/$2" "$FAM" "$FAMILIES" > "${LOG_DIR}/${2%.sh}.log" 2>&1 || { r=failed failed=1; echo "  FAILED"; }
    printf 'step\t%s\t%s\n' "$1" "$r" >> "$STATUS"
  }
  case "${d38%%$'\t'*}" in
    runs) famstep 38 38-family.sh ;;
    waiting) printf 'step\t38\tskipped (waiting: %s)\n' "${d38#*$'\t'}" >> "$STATUS"; echo "  38 Family: waiting: ${d38#*$'\t'}" ;;
    *) printf 'step\t38\tskipped (%s)\n' "${d38#*$'\t'}" >> "$STATUS"; echo "  38 Family: skipped: ${d38#*$'\t'}" ;;
  esac
  case "${d39%%$'\t'*}" in
    runs) famstep 39 39-couple.sh ;;
    waiting) printf 'step\t39\tskipped (waiting: %s)\n' "${d39#*$'\t'}" >> "$STATUS"; echo "  39 Couple carrier screen: waiting: ${d39#*$'\t'}" ;;
    *) printf 'step\t39\tskipped (%s)\n' "${d39#*$'\t'}" >> "$STATUS"; echo "  39 Couple carrier screen: skipped: ${d39#*$'\t'}" ;;
  esac
}
# The OwnHealth file reads summary.json, so it is written only after step 24 wrote a fresh one.
report=failed facts=skipped
for sc in "${OPTIN[@]}" 24-html-report.sh generate-report.sh ownhealth-facts.sh; do
  [ "$sc" != ownhealth-facts.sh ] || [ "$report" = ok ] || { echo "[ownhealth-facts] skipped: 24-html-report.sh failed, so summary.json may be an earlier run's; an earlier run's OwnHealth file was removed"; continue; }
  [ "$sc" != 24-html-report.sh ] || family_steps
  [ "$sc" != 24-html-report.sh ] || GENOME_DIR="$G" HLA_DAT="${H:-}" bash "${PGP_ROOT}/bin/write_manifest.sh" "$SAMPLE" run-all.sh "$SEX" || true
  if [[ "$sc" == [0-9]*-*.sh && "$sc" != 24-* ]]; then live_step "${sc%%-*}" started; fi
  r=ok; echo "[${sc%.sh}] log: ${LOG_DIR}/${sc%.sh}.log"; bash "${SCRIPT_DIR}/${sc}" "$SAMPLE" > "${LOG_DIR}/${sc%.sh}.log" 2>&1 || { r=failed failed=1; echo "  FAILED"; }
  case $sc in 24-html-report.sh) report=$r ;; ownhealth-facts.sh) facts=$r ;; esac
  if [[ "$sc" == [0-9]*-*.sh && "$sc" != 24-* ]]; then printf 'step\t%s\t%s\n' "${sc%%-*}" "$r" >> "$STATUS"; live_step "${sc%%-*}" "$r"; fi
done
[ "$rc" -eq 0 ] || exit "$rc"
[ "$facts" = ok ] && oh=", the OwnHealth file ${S}/facts/${SAMPLE}_genome_facts.json" || oh="; the OwnHealth file was not written (${facts})"
echo "Done: results in ${S}/, the HTML report ${S}/${SAMPLE}_report.html, the text report ${S}/${SAMPLE}_report.txt${oh}"
exit "$failed"
