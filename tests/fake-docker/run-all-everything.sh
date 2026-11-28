#!/usr/bin/env bash
# EVERYTHING=true runs every optional step of assets/everything.tsv.
#   all data      every data folder stubbed, a male BAM sample: every tool row
#                 but html_report (run-all.sh writes the report itself) is in
#                 --tools and listed as runs, with --kir true and --hla_qc
#                 true, and the scripts of TIDDIT, GRIDSS, IMPUTATION, SOMATIC,
#                 every EXTRA_CALLERS and BENCHMARK run after the pipeline
#   no Yleaf      exit 2 before Nextflow, naming setup.sh --yleaf-data, and
#                 the last run's logs/run_status.tsv is left as it was
#   no BWA index  exit 2 before Nextflow, naming the classic BWA index that
#                 script 04b (GRIDSS) needs and bwa index
#   female        still runs; the two Y steps are skipped (female sample)
#   not HiFi      a BAM whose read groups are not PL:PACBIO skips Paraphase
#                 and the run goes on
#   conflicts     TOOLS, or a switch set to another value, stops it with exit 2
# Without EVERYTHING the opt-in steps stay out, as before.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
seed_reference "$G"
seed_clinvar "$G"
use_output_hook
seed_sample "$G" sample1
rm -f "${G}/sample1/vcf/sample1.vcf.gz" "${G}/sample1/vcf/sample1.vcf.gz.tbi"   # a BAM row: DeepVariant writes the gVCF
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"
HLA=$(grep -E '^HLA_DB_RELEASE=' "${SCRIPTS}/lib/common.sh" | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')
GC=$(sed -n 's/^GENCODE_RELEASE=\([0-9]*\)$/\1/p' "${SCRIPTS}/lib/common.sh")
[ -n "$HLA" ] && [ -n "$GC" ] || fail "cannot read HLA_DB_RELEASE and GENCODE_RELEASE from scripts/lib/common.sh"

# Every data folder run-all.sh looks for, as setup.sh leaves it.
f() { local p; for p in "$@"; do mkdir -p "$(dirname "${G}/${p}")"; printf 'placeholder\n' > "${G}/${p}"; done; }
f "hla/IPD-IMGT-HLA_${HLA}/hla.dat" "reference/gencode.v${GC}.basic.genes.gtf" "kir/IPD-KIR_${KIR_DB_RELEASE}/kir.dat" \
  "reference/strchive-${STRCHIVE_VERSION}/catalog.json" \
  "vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt" "annotations/revel_grch38.tsv.gz" "annotations/revel_grch38.tsv.gz.tbi" \
  "vep_cache/homo_sapiens/${PCGR_VEP_CACHE_RELEASE}_GRCh38/info.txt" "pcgr_data/${PCGR_DATA_BUNDLE}/data/README" \
  "reference/cnvpytor/gc_hg38.pytor" "annotsv_annotations/Annotations_Human/Genes/GRCh38/README" "reference/pypgx-bundle/README" \
  "reference/somalier/sites.hg38.vcf.gz" \
  "reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.UD" "reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.mu" \
  "reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.bed" \
  "reference/parascopy-${PARASCOPY_DATA_VERSION}/homology_table/GRCh38.bed.gz" "prs_scores/PGS000018.txt.gz" \
  "prs_scores/PGS000014.txt.gz" "prs_scores/PGS000016.txt.gz" "prs_scores/PGS000017.txt.gz" "prs_scores/PGS000027.txt.gz" \
  "reference/pgsc_calc/${PGSC_PANEL}.tar.zst" "reference/pgsc_calc/${PGSC_PANEL}_GRCh38_sites.tsv" \
  "reference/ancestry_fine/${ANCESTRY_FINE_PANEL}/model/DONE" "reference/yleaf-${YLEAF_DATA_VERSION}/data/hg38/new_positions.txt" \
  "reference/consumer/23andme_v5_sites_dbsnp${CONSUMER_DBSNP_BUILD}.tsv.gz" \
  "reference/ensembl/Homo_sapiens.GRCh38.${ENSEMBL_GFF_RELEASE}.chr.gff3.gz" "reference/archaic/${ARCHAIC_DATA_VERSION}/DONE" \
  "reference/ancient/aadr_${AADR_VERSION}/DONE" "reference/viral-${VIRAL_REF_VERSION}/viral_genomes.fa" \
  "reference/megane/GRCh38_no_alt_analysis_set.mk" "reference/megane/GRCh38_no_alt_analysis_set.mi" \
  "reference/mocha-${MOCHA_PANEL_VERSION}/DONE" "reference/locityper-${LOCITYPER_PANEL_VERSION}/DONE"
# GRIDSS (script 04b) reads the classic BWA index beside the reference.
BWAX="amb ann bwt pac sa"
for e in $BWAX; do f "reference/GRCh38_no_alt_analysis_set.fasta.${e}"; done
# Cyrius counts as installed when its stamp names this lock file and image (setup.sh --cyrius).
if command -v sha256sum > /dev/null 2>&1; then LOCK=$(sha256sum "${SCRIPTS}/cyrius-constraints.txt"); else LOCK=$(shasum -a 256 "${SCRIPTS}/cyrius-constraints.txt"); fi
mkdir -p "${G}/tools/cyrius-${CYRIUS_VERSION}"
printf 'python=%s lock=%s\n' "$PYTHON_IMAGE" "${LOCK%% *}" > "${G}/tools/cyrius-${CYRIUS_VERSION}/INSTALLED"

last() { grep '^nextflow :: ' "$FAKE_DOCKER_LOG" | tail -1; }
nfcalls() { grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG" || true; }
has_args() { grep -qF -- "$(printf '%q ' "$@")" <<<"$(last)" || fail "nextflow lacks '$*': $(last)"; }
tools_arg() { grep -oE -- '--tools [^ ]+' <<<"$(last)" | cut -d' ' -f2 | tr -d '\134' | tr ',' '\n' | LC_ALL=C sort; }
STATUS="${G}/sample1/logs/run_status.tsv"
WANT=$(awk -F'\t' '$1 == "tool" && $2 != "html_report" {print $2}' "${REPO_ROOT}/assets/everything.tsv" | LC_ALL=C sort)
[ "$(grep -c . <<<"$WANT")" -ge 40 ] || fail "assets/everything.tsv has fewer than 40 tool rows: ${WANT}"

# Every data folder there, a male BAM sample: every tool row runs, and every script-only step is queued.
run_rc all env EVERYTHING=true "${SCRIPTS}/run-all.sh" sample1 male
[ "$(nfcalls)" -eq 1 ] || fail "nextflow was not called once"
[ "$(tools_arg)" = "$WANT" ] || fail "--tools is not every tool row of assets/everything.tsv but html_report:
$(diff <(echo "$WANT") <(tools_arg))"
output_lacks all '^  [0-9]+b? .* skipped '
has_args --kir true
has_args --hla_qc true
has_args --pgs_academic true
output_has all '^\[Everything\] EVERYTHING=true'
for s in 04a-tiddit 04b-gridss 14-imputation-prep 29-mutect2-somatic 03a-gatk-haplotypecaller 03b-freebayes \
         03c-strelka2-germline 03d-octopus benchmark-variants; do
  [ -f "${G}/sample1/logs/${s}.log" ] || fail "scripts/${s}.sh did not run after the pipeline"
done
grep -qE $'^step\t04a\t(ok|failed)$' "$STATUS" || fail "step 04a (TIDDIT=true) has no status line"
grep -qE $'^meta\tcommand\t.*EVERYTHING=true ' "$STATUS" || fail "the recorded command lacks EVERYTHING=true: $(grep command "$STATUS")"
if grep -qE $'^meta\tcommand\t.*TOOLS=' "$STATUS"; then fail "the recorded command names TOOLS: a replay would refuse it"; fi

# Without Yleaf's tables: exit 2 before Nextflow, naming the install command.
mv "${G}/reference/yleaf-${YLEAF_DATA_VERSION}" "${CASE_WORK}/yleaf-away"
n=$(nfcalls)
cp "$STATUS" "${CASE_WORK}/status-before"
run_expect 2 no-yleaf env EVERYTHING=true "${SCRIPTS}/run-all.sh" sample1 male
[ "$(nfcalls)" -eq "$n" ] || fail "nextflow started without Yleaf's tables"
output_has no-yleaf 'ERROR: EVERYTHING=true runs every optional step, and 1 cannot run \(Nextflow did not start\)'
output_has no-yleaf '^  37 Y haplogroup \(Yleaf\): data not installed: reference/yleaf-.*; install: setup.sh --yleaf-data$'
cmp -s "${CASE_WORK}/status-before" "$STATUS" || fail "a refused EVERYTHING=true run changed the last run's status:
$(diff "${CASE_WORK}/status-before" "$STATUS")"
[ ! -e "${STATUS}.prev" ] || fail "a refused EVERYTHING=true run left ${STATUS}.prev behind"

# A female sample still runs without them: the Y steps are skipped for the sample, not for want of data.
run_rc female env EVERYTHING=true IMPUTATION=true "${SCRIPTS}/run-all.sh" sample1 female
[ "$(nfcalls)" -eq $((n + 1)) ] || fail "nextflow did not run for a female sample"
output_has female '^  37 Y haplogroup \(Yleaf\) +skipped +\(female sample: no Y chromosome\)$'
output_has female '^  45 Y-STR markers \(HipSTR\) +skipped +\(female sample: no Y chromosome\)$'
if tools_arg | grep -qE '^y_(haplogroup|str)$'; then fail "a Y step reached --tools for a female sample: $(last)"; fi
mv "${CASE_WORK}/yleaf-away" "${G}/reference/yleaf-${YLEAF_DATA_VERSION}"

# Without the classic BWA index: GRIDSS (script 04b, after the pipeline) would fail hours in, so the run stops first.
for e in $BWAX; do mv "${G}/reference/GRCh38_no_alt_analysis_set.fasta.${e}" "${CASE_WORK}/"; done
run_expect 2 no-bwa env EVERYTHING=true "${SCRIPTS}/run-all.sh" sample1 male
[ "$(nfcalls)" -eq $((n + 1)) ] || fail "nextflow started without the classic BWA index"
output_has no-bwa '^  04b GRIDSS \(GRIDSS=true\): classic BWA index not installed: reference/GRCh38_no_alt_analysis_set\.fasta\.amb; install: \./scripts/setup\.sh --bwa-index \(docs/04b-gridss\.md\)$'
for e in $BWAX; do mv "${CASE_WORK}/GRCh38_no_alt_analysis_set.fasta.${e}" "${G}/reference/"; done

# Without the academic-only PGS scores (PRS_ACADEMIC=true, setup.sh --prs-academic): the run stops first.
mv "${G}/prs_scores/PGS000016.txt.gz" "${CASE_WORK}/"
run_expect 2 no-academic env EVERYTHING=true "${SCRIPTS}/run-all.sh" sample1 male
[ "$(nfcalls)" -eq $((n + 1)) ] || fail "nextflow started without the academic-only PGS scores"
output_has no-academic '^  25 PRS academic-only scores \(PRS_ACADEMIC=true\): data not installed: PGS000016; install: setup.sh --prs-academic$'
mv "${CASE_WORK}/PGS000016.txt.gz" "${G}/prs_scores/"

# A BAM that is not PacBio HiFi: og paraphase check-bam exits 3, and step 35b is a skip, not a stop.
cat > "${CASE_WORK}/not-hifi" <<'HOOK'
#!/usr/bin/env bash
case "$*" in *"og paraphase check-bam "*) exit 3 ;; esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/not-hifi"
run_rc not-hifi env EVERYTHING=true FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/not-hifi" "${SCRIPTS}/run-all.sh" sample1 male
[ "$(nfcalls)" -eq $((n + 2)) ] || fail "nextflow did not run with a BAM that is not PacBio HiFi"
output_has not-hifi '^  35b Paraphase \(PacBio HiFi\) +skipped +\(not PacBio HiFi: a read group of the BAM is not PL:PACBIO\)$'
if tools_arg | grep -qx paraphase; then fail "paraphase reached --tools for a BAM that is not PacBio HiFi"; fi
docker_log_has 'og paraphase check-bam /genome/sample1/aligned/sample1_sorted\.bam' "run-all.sh did not check the BAM with og paraphase check-bam"

# EVERYTHING names the steps and sets the switches itself: another value is refused before anything runs.
run_expect 2 with-tools env EVERYTHING=true TOOLS=pharmcat "${SCRIPTS}/run-all.sh" sample1 male
output_has with-tools 'EVERYTHING=true names every step itself: unset TOOLS'
run_expect 2 kir-off env EVERYTHING=true KIR=false "${SCRIPTS}/run-all.sh" sample1 male
output_has kir-off 'EVERYTHING=true sets KIR=true: unset KIR'
run_expect 2 none env EVERYTHING=true ARCHAIC_DIR=none "${SCRIPTS}/run-all.sh" sample1 male
output_has none '^  43 Archaic ancestry \(IBDmix\): ARCHAIC_DIR=none: unset ARCHAIC_DIR$'
[ "$(nfcalls)" -eq $((n + 2)) ] || fail "nextflow started on a refused EVERYTHING=true run"

# Without EVERYTHING a default run leaves the opt-in steps out, as before.
run_expect 0 default "${SCRIPTS}/run-all.sh" sample1 male
output_has default '^  37 Y haplogroup \(Yleaf\) +skipped +\(opt-in: add y_haplogroup to TOOLS\)$'
if grep -qE -- ' --kir true' <<<"$(last)"; then fail "a default run passed --kir true"; fi
if grep -qE -- ' --pgs_academic true' <<<"$(last)"; then fail "a default run passed --pgs_academic true"; fi
echo "EVERYTHING=true runs every optional step, and stops before Nextflow when one lacks its data."
