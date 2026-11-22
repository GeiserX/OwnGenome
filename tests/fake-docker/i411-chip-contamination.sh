#!/usr/bin/env bash
# Step 49 (CHIP) gives FilterMutectCalls the contamination estimate when step
# 29's common-sites VCF is installed with its .tbi: GetPileupSummaries at every
# site of that VCF (not only inside the CHIP exons), CalculateContamination, then
# FilterMutectCalls with --contamination-table and --tumor-segmentation.
# Without the VCF, or without its .tbi, there is no estimate and the log says
# so; with it but without gnomAD's germline resource the log warns that the
# estimate changes few calls. The Nextflow CHIP_MUTECT2 runs the same commands: its script block, run
# here with a fake gatk, passes the table when it gets --chip_common_sites and
# not without it; main.nf hands it the file, and run-all.sh passes the VCF
# under somatic/ to the pipeline only with its .tbi.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"

export GENOME_DIR="${CASE_WORK}/genome"
G=$GENOME_DIR
seed_reference "$G"
seed_sample "$G" sample1
printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n' > "${G}/reference/GRCh38_no_alt_analysis_set.dict"
GFF="${G}/reference/ensembl/Homo_sapiens.GRCh38.${ENSEMBL_GFF_RELEASE}.chr.gff3.gz"
mkdir -p "$(dirname "$GFF")"
printf 'placeholder\n' > "$GFF"
COMMON="${G}/somatic/small_exac_common_3.hg38.vcf.gz"

# og chip prepare writes a BED with one exon and og chip report an empty table,
# so the step goes on past its checks; every other call writes its -O output.
use_output_hook
cat > "${CASE_WORK}/hook-chip" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
args=" $* "
# arg OPTION: the value after OPTION in this call
arg() { sed -E "s/.* $1 ([^ ]+) .*/\\1/" <<<"$args"; }
case "$args" in
  *"/og chip prepare "*) printf 'chr2\t25234268\t25234430\tDNMT3A\n' > "$(host_path "$(arg --out-bed)")" ;;
  *"/og chip report "*) : > "$(host_path "$(arg --out)")" ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/hook-chip"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hook-chip"
C=/genome/sample1/chip/sample1_chip

# --- the bash step ---------------------------------------------------------------------
# The common-sites VCF with its .tbi: the estimate, and FilterMutectCalls reads it.
mkdir -p "${G}/somatic"
printf 'placeholder\n' > "$COMMON"
printf 'placeholder\n' > "${COMMON}.tbi"
: > "$FAKE_DOCKER_LOG"
run_expect 0 chip-common "${SCRIPTS}/49-chip.sh" sample1
docker_log_has "^run image=[^ ]*gatk.* FilterMutectCalls .*--contamination-table ${C}_contamination\\.table --tumor-segmentation ${C}_segments\\.table" \
  "FilterMutectCalls did not get the contamination table"
docker_log_has "^run image=[^ ]*gatk.* GetPileupSummaries -I /genome/sample1/aligned/sample1_sorted\\.bam -V /genome/somatic/small_exac_common_3\\.hg38\\.vcf\\.gz -L /genome/somatic/small_exac_common_3\\.hg38\\.vcf\\.gz -O ${C}_pileups\\.table" \
  "step 49 did not summarise the pileups at every common site"
docker_log_has "^run image=[^ ]*gatk.* CalculateContamination -I ${C}_pileups\\.table --tumor-segmentation ${C}_segments\\.table -O ${C}_contamination\\.table" \
  "step 49 did not estimate contamination from those pileups"
# The estimate comes before the filter that reads it.
awk '/ CalculateContamination / { c = NR } / FilterMutectCalls / { f = NR } END { exit !(c && f > c) }' "$FAKE_DOCKER_LOG" \
  || fail "FilterMutectCalls ran before CalculateContamination"
output_lacks chip-common 'no contamination estimate'
output_has chip-common 'contamination estimate: common sites '
output_has chip-common 'WARNING: the contamination estimate without gnomAD'
# With gnomAD's germline resource too: Mutect2 reads it and the warning goes.
printf 'placeholder\n' > "${G}/somatic/af-only-gnomad.hg38.vcf.gz"
printf 'placeholder\n' > "${G}/somatic/af-only-gnomad.hg38.vcf.gz.tbi"
: > "$FAKE_DOCKER_LOG"
run_expect 0 chip-common-gnomad "${SCRIPTS}/49-chip.sh" sample1
docker_log_has '^run image=[^ ]*gatk.* Mutect2 .*--germline-resource /genome/somatic/af-only-gnomad\.hg38\.vcf\.gz' "Mutect2 did not get gnomAD"
docker_log_has "^run image=[^ ]*gatk.* FilterMutectCalls .*--contamination-table ${C}_contamination\\.table" "FilterMutectCalls with gnomAD did not get the contamination table"
output_lacks chip-common-gnomad 'without gnomAD'
rm -f "${G}/somatic/af-only-gnomad.hg38.vcf.gz" "${G}/somatic/af-only-gnomad.hg38.vcf.gz.tbi"

# The VCF without its .tbi is not used.
rm -f "${COMMON}.tbi"
: > "$FAKE_DOCKER_LOG"
run_expect 0 chip-common-no-tbi "${SCRIPTS}/49-chip.sh" sample1
if grep -q -- '--contamination-table' "$FAKE_DOCKER_LOG"; then fail "step 49 used the common-sites VCF without its .tbi"; fi
output_has chip-common-no-tbi 'no common-sites VCF .*: no contamination estimate'

# No common-sites VCF: no estimate, and the log says so.
rm -f "$COMMON"
: > "$FAKE_DOCKER_LOG"
run_expect 0 chip-no-common "${SCRIPTS}/49-chip.sh" sample1
docker_log_has '^run image=[^ ]*gatk.* FilterMutectCalls .*--ob-priors ' "step 49 did not run FilterMutectCalls"
if grep -qE 'GetPileupSummaries|CalculateContamination|--contamination-table' "$FAKE_DOCKER_LOG"; then
  fail "step 49 estimated contamination without the common-sites VCF"
fi
output_has chip-no-common 'no common-sites VCF .*: no contamination estimate'
printf 'placeholder\n' > "$COMMON"
printf 'placeholder\n' > "${COMMON}.tbi"

# --- CHIP_MUTECT2 (Nextflow) -------------------------------------------------------------
# The shell of the module's script block, with its Groovy values filled in.
MOD="${REPO_ROOT}/modules/local/chip/main.nf"
# shellcheck disable=SC2016  # the dollar signs are Groovy placeholders, matched literally
render_mutect2() {
  awk '/^process CHIP_MUTECT2/ { p = 1 } p && /^    script:/ { on = 1; next } on && /^    """/ { n++; if (n == 2) exit; next } on && n == 1' "$MOD" \
    | sed -e 's/\${reference}/ref.fa/g' -e 's/\${bam}/in.bam/g' -e 's/\${bed}/chip.bed/g' -e 's/\${task\.cpus}/2/g' \
          -e 's/\${resources}//g' -e 's/\${meta\.id}/s1/g' -e "s|\\\${common_sites}|$1|g" -e '/versions\.yml/d' \
          -e 's/\\\$/$/g' -e 's/\\\\/\\/g'
}
NFW="${CASE_WORK}/chip-nf"
mkdir -p "${NFW}/bin"
cat > "${NFW}/bin/gatk" <<'FAKE'
#!/usr/bin/env bash
printf 'gatk %s\n' "$*" >> "${NFW_LOG:?}"
FAKE
chmod +x "${NFW}/bin/gatk"
# run_module COMMON_SITES: run the rendered block in $NFW, with the gatk calls in gatk.log
run_module() {
  render_mutect2 "$1" > "${NFW}/script.sh"
  : > "${NFW}/gatk.log"
  (cd "$NFW" && PATH="${NFW}/bin:${PATH}" NFW_LOG="${NFW}/gatk.log" bash -euo pipefail script.sh) > "${NFW}/out.txt" 2>&1 \
    || { cat "${NFW}/script.sh" "${NFW}/out.txt"; fail "CHIP_MUTECT2's script block failed with common sites '$1'"; }
  cat "${NFW}/gatk.log" "${NFW}/out.txt"
}
run_module common.vcf.gz
grep -q '^gatk GetPileupSummaries -I in\.bam -V common\.vcf\.gz -L common\.vcf\.gz -O pileups\.table$' "${NFW}/gatk.log" \
  || fail "CHIP_MUTECT2 did not summarise the pileups at every common site, as step 49 does"
grep -q '^gatk CalculateContamination -I pileups\.table --tumor-segmentation segments\.table -O contamination\.table$' "${NFW}/gatk.log" \
  || fail "CHIP_MUTECT2 did not estimate contamination"
grep -q '^gatk FilterMutectCalls .*--contamination-table contamination\.table --tumor-segmentation segments\.table .*-O s1_chip_filtered\.vcf\.gz$' "${NFW}/gatk.log" \
  || fail "CHIP_MUTECT2's FilterMutectCalls did not get the contamination table"
run_module ""
grep -q '^gatk FilterMutectCalls .*-O s1_chip_filtered\.vcf\.gz$' "${NFW}/gatk.log" || fail "CHIP_MUTECT2 without common sites did not filter"
if grep -qE 'GetPileupSummaries|CalculateContamination|--contamination-table' "${NFW}/gatk.log"; then
  fail "CHIP_MUTECT2 estimated contamination without --chip_common_sites"
fi
grep -q 'no --chip_common_sites' "${NFW}/out.txt" || fail "CHIP_MUTECT2 without common sites does not say it has no estimate"
# main.nf warns on the console (a task's stderr stays in its work directory)
# without the common sites, and with them but without gnomAD.
grep -qF 'if (!common) log.warn "chip: no --chip_common_sites, so FilterMutectCalls runs without a contamination estimate' "${REPO_ROOT}/main.nf" \
  || fail "main.nf does not warn when chip runs without --chip_common_sites"
grep -qF 'else if (!germ) log.warn "chip: --chip_common_sites without --chip_germline_resource' "${REPO_ROOT}/main.nf" \
  || fail "main.nf does not warn when the contamination estimate runs without gnomAD"
# main.nf hands the module the file and its index, in the order of its inputs.
grep -qF 'def chip_common       = common ? file(common, checkIfExists: true) : []' "${REPO_ROOT}/main.nf" \
  || fail "main.nf does not read --chip_common_sites"
grep -qF 'chip_pon, chip_pon_index, chip_common, chip_common_index)' "${REPO_ROOT}/main.nf" \
  || fail "main.nf does not pass the common sites to CHIP_MUTECT2"
# The rendered script substitutes the VCF directly, so check the positional
# binding of the VCF and its index in the module's input declarations too.
awk '
  /^process CHIP_MUTECT2[[:space:]]*\{/ { in_process = 1 }
  in_process && /^}/ { exit }
  in_process && /^[[:space:]]*input:/ { in_input = 1; next }
  in_input && /^[[:space:]]*output:/ { exit }
  in_input {
    sub(/[[:space:]]*\/\/.*/, "")
    gsub(/^[[:space:]]+|[[:space:]]+$/, "")
    if (two_back == "path(pon_index)" && previous == "path(common)" && $0 == "path(common_index)") {
      found = 1
      exit
    }
    two_back = previous
    previous = $0
  }
  END { exit !found }
' "$MOD" || fail "CHIP_MUTECT2 does not take the common sites and their index after the panel of normals"

# --- run-all.sh --------------------------------------------------------------------------
# TOOLS=chip passes the installed common-sites VCF to the pipeline, only with its .tbi.
export SKIP_VALIDATION=true
seed_clinvar "$G"
: > "$FAKE_DOCKER_LOG"
TOOLS=chip run_expect 0 run-all-common "${SCRIPTS}/run-all.sh" sample1 male
docker_log_has "^nextflow :: .* --chip_common_sites ${COMMON//./\\.} " "run-all.sh did not pass the installed common-sites VCF"
rm -f "${COMMON}.tbi"
: > "$FAKE_DOCKER_LOG"
TOOLS=chip run_expect 0 run-all-no-tbi "${SCRIPTS}/run-all.sh" sample1 male
docker_log_has '^nextflow :: .* --tools ' "run-all.sh did not run the pipeline"
if grep -q -- '--chip_common_sites' "$FAKE_DOCKER_LOG"; then fail "run-all.sh passed the common-sites VCF without its .tbi"; fi
echo "Step 49, CHIP_MUTECT2 and run-all.sh give FilterMutectCalls the contamination estimate with the common sites, and only then."
