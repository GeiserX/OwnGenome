#!/usr/bin/env bash
# Step 43 (archaic ancestry, IBDmix) reads the gVCF: a site with no call is the
# reference only inside a reference block. With the archaic table installed,
# run-all.sh lists step 43 as skipped with 'no gVCF: <sample>/vcf/<sample>.g.vcf.gz'
# for a row that starts from a VCF with no gVCF and index beside it, as it does
# for step 42, and records the same reason in logs/run_status.tsv:
#   1. BAM+VCF with the gVCF and its index beside the VCF: step 43 runs and is
#      'ok' (the control: the table counts as installed);
#   2. BAM+VCF, no gVCF: skipped, and nextflow gets no --archaic_data;
#   3. a VCF-only row: skipped for the same reason.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
seed_reference "$G"
seed_clinvar "$G"
use_output_hook
AV=$(sed -n 's/^ARCHAIC_DATA_VERSION="\(.*\)".*/\1/p' "${REPO_ROOT}/versions.env")
[ -n "$AV" ] || fail "cannot read ARCHAIC_DATA_VERSION from versions.env"
mkdir -p "${G}/reference/archaic/${AV}"
printf 'placeholder\n' > "${G}/reference/archaic/${AV}/DONE"   # setup.sh --archaic writes DONE last

MISSES=0
miss() { echo "ASSERT FAIL: $*" >&2; MISSES=$((MISSES + 1)); }
# step43 NAME SAMPLE: step 43's result in that run's run_status.tsv (its last row)
step43() { awk -F'\t' '$1 == "step" && $2 == "43" {r = $3} END {print r}' "${G}/$2/logs/run_status.tsv"; }

# 1. The control: the gVCF and its index beside the VCF.
seed_sample "$G" s1
printf 'placeholder\n' > "${G}/s1/vcf/s1.g.vcf.gz"
printf 'placeholder\n' > "${G}/s1/vcf/s1.g.vcf.gz.tbi"
run_expect 0 with-gvcf "${SCRIPTS}/run-all.sh" s1 male
output_has with-gvcf '^  43 Archaic ancestry \(IBDmix\) +runs$'
[ "$(step43 with-gvcf s1)" = ok ] || miss "with-gvcf: step 43 is '$(step43 with-gvcf s1)', want 'ok'"
docker_log_has "^nextflow :: .* --archaic_data ${G}/reference/archaic/${AV} " "with-gvcf: nextflow was not given --archaic_data"

# 2. BAM+VCF with no gVCF beside the VCF.
seed_sample "$G" s2
N=$(wc -l < "$FAKE_DOCKER_LOG")   # this run's lines of the docker log start after line N
run_expect 0 no-gvcf "${SCRIPTS}/run-all.sh" s2 male
output_has no-gvcf '^  43 Archaic ancestry \(IBDmix\) +skipped +\(no gVCF: s2/vcf/s2\.g\.vcf\.gz\)$'
output_has no-gvcf '^  42 Traits panel +skipped +\(no gVCF: s2/vcf/s2\.g\.vcf\.gz\)$'
output_lacks no-gvcf '^  43 .*data not installed'
[ "$(step43 no-gvcf s2)" = 'skipped (no gVCF: s2/vcf/s2.g.vcf.gz)' ] || miss "no-gvcf: step 43 is '$(step43 no-gvcf s2)'"
# One awk, no pipe: under pipefail, grep -q closing a pipe early would hide the match.
if awk -v n="$N" 'NR > n && / --archaic_data / { hit = 1 } END { exit !hit }' "$FAKE_DOCKER_LOG"; then
  miss "no-gvcf: nextflow was given --archaic_data"
fi

# 3. A VCF-only row.
mkdir -p "${G}/s3/vcf"
printf 'placeholder\n' > "${G}/s3/vcf/s3.vcf.gz"
printf 'placeholder\n' > "${G}/s3/vcf/s3.vcf.gz.tbi"
run_expect 0 vcf-only "${SCRIPTS}/run-all.sh" s3 male
output_has vcf-only '^  43 Archaic ancestry \(IBDmix\) +skipped +\(no gVCF: s3/vcf/s3\.g\.vcf\.gz\)$'
[ "$(step43 vcf-only s3)" = 'skipped (no gVCF: s3/vcf/s3.g.vcf.gz)' ] || miss "vcf-only: step 43 is '$(step43 vcf-only s3)'"

[ "$MISSES" -eq 0 ] || fail "${MISSES} result(s) above are wrong"
echo "With the archaic table installed, step 43 runs on a row with its gVCF and is skipped as 'no gVCF' on a BAM+VCF or VCF-only row without one."
