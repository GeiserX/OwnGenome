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

run_expect 0 vcf env INPUT=vcf "${launcher[@]}"
expect_sheet "$header" "${sample},,,,,${vcf},${vcf}.tbi,female"
output_has vcf '^Starting from the existing VCF: reusing variant calls$'

# Mark the last run complete so auto takes the highest-ranked input.
printf 'process\tsucceeded\tcached\tfailed\tignored\tunfinished\nROH\t1\t0\t0\t0\t0\n' > "${root}/nextflow/process_status.tsv"
run_expect 0 auto "${launcher[@]}"
expect_sheet "${header},gvcf,gvcf_index" "${sample},,,${bam},${bam}.bai,${vcf},${vcf}.tbi,female,${gvcf},${gvcf}.tbi"
output_has auto '^Starting from the existing BAM and VCF: alignment and calling are reused$'
output_lacks auto '^\[Input\].*samplesheet\.csv:'
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
