#!/usr/bin/env bash
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
use_output_hook

sample=s1
root="${GENOME_DIR}/${sample}"
bam="${root}/aligned/${sample}_sorted.bam"
vcf="${root}/vcf/${sample}.vcf.gz"
gvcf="${root}/vcf/${sample}.g.vcf.gz"
sheet="${root}/nextflow/samplesheet.csv"
header='sample,fastq_1,fastq_2,bam,bam_index,vcf,vcf_index,sex'
launcher=("${SCRIPTS}/run-all.sh" "$sample" female)
expect_sheet() {
  local want
  want=$(printf '%s\n%s\n' "$1" "$2")
  [ "$(cat "$sheet")" = "$want" ] || fail "unexpected samplesheet: $(cat "$sheet")"
}

seed_sample "$GENOME_DIR" "$sample"
touch "$gvcf" "${gvcf}.tbi"
run_expect 0 auto_gvcf "${launcher[@]}"
expect_sheet "${header},gvcf,gvcf_index" "${sample},,,${bam},${bam}.bai,${vcf},${vcf}.tbi,female,${gvcf},${gvcf}.tbi"
note=$(grep '^NOTE: .*gVCF' "${CASE_WORK}/auto_gvcf.out")
[[ "$note" == *"$gvcf"* ]] || fail "the NOTE line must name the gVCF"

seed_sample "$GENOME_DIR" s2
rm -f "${GENOME_DIR}/s2/aligned/s2_sorted.bam" "${GENOME_DIR}/s2/aligned/s2_sorted.bam.bai"
touch "${GENOME_DIR}/s2/vcf/s2.g.vcf.gz" "${GENOME_DIR}/s2/vcf/s2.g.vcf.gz.tbi"
run_expect 0 vcf_only "${SCRIPTS}/run-all.sh" s2 female
[ "$(head -n 1 "${GENOME_DIR}/s2/nextflow/samplesheet.csv")" = "$header" ] || fail "VCF-only input has gVCF columns"

mkdir -p "${root}/fastq"
r1="${root}/fastq/${sample}_R1.fastq.gz" r2="${root}/fastq/${sample}_R2.fastq.gz"
printf 'new reads\n' > "$r1"
printf 'new reads\n' > "$r2"
run_expect 0 fastq env INPUT=fastq "${launcher[@]}"
expect_sheet "$header" "${sample},${r1},${r2},,,,,female"
output_has fastq '^Starting from FASTQ: aligning and calling$'
grep -qF 'INPUT=fastq ' "${root}/logs/run_status.tsv" || fail "recheck would lose INPUT=fastq"

run_expect 0 bam env INPUT=bam "${launcher[@]}"
expect_sheet "$header" "${sample},,,${bam},${bam}.bai,,,female"
output_has bam '^Starting from the existing BAM: calling variants$'

mkdir -p "${CASE_WORK}/archaic"
touch "${CASE_WORK}/archaic/DONE"
run_expect 0 vcf env INPUT=vcf TOOLS=consumer_export,traits,archaic,ancient \
  ARCHAIC_DIR="${CASE_WORK}/archaic" "${launcher[@]}"
expect_sheet "$header" "${sample},,,,,${vcf},${vcf}.tbi,female"
output_has vcf '^Starting from the existing VCF: reusing variant calls$'

run_expect 0 auto_resume "${launcher[@]}"
expect_sheet "$header" "${sample},,,,,${vcf},${vcf}.tbi,female"
output_has auto_resume '^Starting from the existing VCF: reusing variant calls$'

for step in 40 42 43 44; do
  output_has vcf "^[[:space:]]+${step} .*skipped [[:space:]]+\\(INPUT=vcf: the gVCF beside the VCF is not read\\)$"
done

# Mark the last run complete so auto takes the highest-ranked input.
printf 'process\tsucceeded\tcached\tfailed\tignored\tunfinished\nROH\t1\t0\t0\t0\t0\n' > "${root}/nextflow/process_status.tsv"
run_expect 0 auto "${launcher[@]}"
expect_sheet "${header},gvcf,gvcf_index" "${sample},,,${bam},${bam}.bai,${vcf},${vcf}.tbi,female,${gvcf},${gvcf}.tbi"
output_has auto '^Starting from the existing BAM and VCF: alignment and calling are reused$'
output_lacks auto '^\[Input\].*samplesheet\.csv:'

run_expect 0 empty_input env INPUT= "${launcher[@]}"
expect_sheet "${header},gvcf,gvcf_index" "${sample},,,${bam},${bam}.bai,${vcf},${vcf}.tbi,female,${gvcf},${gvcf}.tbi"
output_has empty_input '^Starting from the existing BAM and VCF: alignment and calling are reused$'
for f in "$bam" "${bam}.bai" "$vcf" "${vcf}.tbi" "$gvcf" "${gvcf}.tbi"; do
  [ -f "$f" ] || fail "the launcher removed $f"
done

calls=$(grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG")
run_expect 2 missing_bam env INPUT=bam "${SCRIPTS}/run-all.sh" s2 female
grep -qF "${GENOME_DIR}/s2/aligned/s2_sorted.bam" "${CASE_WORK}/missing_bam.out" || fail "missing BAM was not named"
run_expect 2 invalid env INPUT=other "${launcher[@]}"
output_has invalid 'INPUT must be auto, fastq, bam, or vcf'
[ "$(grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG")" -eq "$calls" ] || fail "Nextflow ran for an invalid or missing input"

if grep -niE 'does not read a gVCF beside it|gVCF.*never (reads?|read)' "${REPO_ROOT}/docs/getting-started.md"; then
  fail "Getting started says an existing VCF never gets its gVCF"
fi

# Model a run recorded before INPUT was added to the switch list.
printf '##fileformat=VCFv4.1\n##fileDate=20260901\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' \
  | gzip -c > "${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz"
run_expect 0 replay_source "${launcher[@]}"
status="${root}/logs/run_status.tsv"
awk -F '\t' -v OFS='\t' '$1 == "meta" && $2 == "switches" {sub(/(^| )INPUT( |$)/, " ", $3); sub(/^ /, "", $3)} {print}' \
  "$status" > "${status}.old"
mv "${status}.old" "$status"
calls=$(grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG")
run_expect 0 rerun env INPUT=vcf "${REPO_ROOT}/bin/og" recheck --sample "$sample" \
  --offline "${REPO_ROOT}/tests/fixtures/recheck/offline-nothing-newer.json" --rerun
[ "$(grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG")" -eq "$((calls + 1))" ] || fail "Replay did not start Nextflow"
expect_sheet "${header},gvcf,gvcf_index" "${sample},,,${bam},${bam}.bai,${vcf},${vcf}.tbi,female,${gvcf},${gvcf}.tbi"

cat > "${CASE_WORK}/validation-hook" <<'HOOK'
#!/usr/bin/env bash
case "$*" in
  *"samtools view -H "*)
    printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:%s\n@RG\tID:fixture\tSM:s1\n' "${BAM_CHR1_LENGTH:-248956422}" ;;
  *"bcftools view -h "*)
    printf '##fileformat=VCFv4.2\n##contig=<ID=%s,length=%s>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n' \
      "${VCF_CHR1:-chr1}" "${VCF_CHR1_LENGTH:-248956422}" ;;
  *) exec "${CASE_WORK}/hook-outputs" "$@" ;;
esac
HOOK
chmod +x "${CASE_WORK}/validation-hook"
validate=(env -u SKIP_VALIDATION "FAKE_DOCKER_RUN_HOOK=${CASE_WORK}/validation-hook")
calls=$(grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG")
run_expect 0 fastq_wrong_build_bam "${validate[@]}" INPUT=fastq BAM_CHR1_LENGTH=249250621 "${launcher[@]}"
expect_sheet "$header" "${sample},${r1},${r2},,,,,female"
output_has fastq_wrong_build_bam 'INPUT=fastq: the existing BAM is not read, so it is not checked'
output_has fastq_wrong_build_bam 'INPUT=fastq: the existing VCF is not read, so it is not checked'
[ "$(grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG")" -eq "$((calls + 1))" ] || fail "FASTQ run did not start Nextflow"

calls=$(grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG")
run_expect 0 bam_wrong_build_vcf "${validate[@]}" INPUT=bam VCF_CHR1=1 VCF_CHR1_LENGTH=249250621 "${launcher[@]}"
expect_sheet "$header" "${sample},,,${bam},${bam}.bai,,,female"
output_has bam_wrong_build_vcf 'INPUT=bam: the existing VCF is not read, so it is not checked'
[ "$(grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG")" -eq "$((calls + 1))" ] || fail "BAM run did not start Nextflow"

calls=$(grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG")
run_expect 1 auto_wrong_build_bam "${validate[@]}" BAM_CHR1_LENGTH=249250621 "${launcher[@]}"
output_has auto_wrong_build_bam '\[FAIL\].*BAM genome build: GRCh37/hg19'
[ "$(grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG")" -eq "$calls" ] || fail "Auto started Nextflow with a wrong-build BAM"

# An explicit INPUT still checks the file it reads, and only skips the others.
calls=$(grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG")
run_expect 1 bam_wrong_build_bam "${validate[@]}" INPUT=bam BAM_CHR1_LENGTH=249250621 "${launcher[@]}"
output_has bam_wrong_build_bam '\[FAIL\].*BAM genome build: GRCh37/hg19'
run_expect 1 vcf_nochr_vcf "${validate[@]}" INPUT=vcf VCF_CHR1=1 VCF_CHR1_LENGTH=249250621 "${launcher[@]}"
output_has vcf_nochr_vcf "chromosome names without the 'chr' prefix"
[ "$(grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG")" -eq "$calls" ] || fail "Nextflow ran with a wrong-build input the run reads"
run_expect 0 vcf_wrong_build_bam "${validate[@]}" INPUT=vcf BAM_CHR1_LENGTH=249250621 "${launcher[@]}"
output_has vcf_wrong_build_bam 'INPUT=vcf: the BAM is not read, so it is not checked'
output_has vcf_wrong_build_bam 'INPUT=vcf: the BAM is not read — BAM-dependent steps \(4, 10, 15, 16, 18, 19, 20\) will be skipped\.'
output_lacks vcf_wrong_build_bam 'No BAM found'
[ "$(grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG")" -eq "$((calls + 1))" ] || fail "VCF run did not start Nextflow"
