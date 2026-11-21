#!/usr/bin/env bash
# The single-process callers scatter: on a reference with chr1, chr2 and an
# unplaced contig, steps 03a (GATK), 03b (FreeBayes) and 29 with
# INTERVALS=genome (Mutect2) each run three units, chr1, chr2 and the rest,
# at most SCATTER_JOBS at a time, and join them: bcftools concat for GATK, one
# raw VCF for FreeBayes, MergeVcfs, MergeMutectStats and every unit's
# orientation counts for Mutect2. SCATTER=false runs one process. A contig in
# INTERVALS that the reference lacks stops the step before any container.
# A rerun of 03a or 03b keeps a finished VCF called with the same INTERVALS,
# SCATTER, BAM, reference and image, and starts no caller; another
# INTERVALS or SCATTER, a realigned BAM, a BAM replaced under the same size
# and time, or a missing record calls again; so does a replaced BAM for step
# 29. A region off its contig, a start of 0, an end before its start and two
# regions that overlap stop the step before any container; regions that only
# touch are two units. With THREADS=16 on the fake docker's 16 GB, two 8 GB
# units run at once, not eight; SCATTER_JOBS still sets the number.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" THREADS=4
seed_reference "$GENOME_DIR"
REF="${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set"
printf 'chr1\t248956422\t6\t60\t61\nchr2\t242193529\t253105752\t60\t61\nchrUn_KI270302v1\t2274\t499335640\t60\t61\n' > "${REF}.fasta.fai"
printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n' > "${REF}.dict"
seed_sample "$GENOME_DIR" sample1
use_output_hook
cat > "${CASE_WORK}/scatter-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
# bgzf HOST_PATH: a finished bgzipped VCF (ends with the BGZF end-of-file block)
bgzf() {
  mkdir -p "$(dirname "$1")"
  { printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n' | gzip -c
    printf '\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'; } > "$1"
}
case " ${*:2} " in
  *" freebayes "*)
    t=$(sed -n 's/.* --targets \([^ ]*\) .*/\1/p' <<<" ${*:2} ")
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n'
    printf 'chr1\t%s\t.\tA\tG\t50\t.\t.\tGT\t0/1\n' "${#t}" ;;
  *" bcftools concat "*|*" bcftools sort "*)
    bgzf "$(host_path "$(sed -n 's/.* -o \([^ ]*\) .*/\1/p' <<<" ${*:2} ")")"
    exit 0 ;;
  *" bcftools index "*)
    printf 'TBI\001' > "$(host_path "${*: -1}").tbi"
    exit 0 ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/scatter-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/scatter-hook"
count() { grep -cE -- "$1" "$FAKE_DOCKER_LOG" || true; }

# --- 03a GATK HaplotypeCaller ------------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 gatk "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1
output_has gatk '3 unit\(s\), 2 at a time'
[ "$(count ' HaplotypeCaller .*-L /genome/sample1/vcf_gatk/scatter/00[123]\.bed ')" -eq 3 ] || fail "03a did not run HaplotypeCaller once per unit"
docker_log_has 'bcftools concat -a -D .*scatter/001\.vcf\.gz .*scatter/002\.vcf\.gz .*scatter/003\.vcf\.gz' "03a did not join the three units in order"
[ ! -e "${GENOME_DIR}/sample1/vcf_gatk/scatter" ] || fail "03a kept its scatter folder"
G="${GENOME_DIR}/sample1/vcf_gatk"
grep -q '^INTERVALS= scatter=true bam=.* reference=.* image=' "${G}/sample1.run" 2>/dev/null || fail "03a wrote no run record beside its VCF"
# A rerun of the same call keeps the VCF and starts no container for GATK.
: > "$FAKE_DOCKER_LOG"
run_expect 0 gatk-again "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1
output_has gatk-again 'Output already exists'
[ "$(count ' HaplotypeCaller ')" -eq 0 ] || fail "a rerun of 03a called again a VCF it had finished"
# A BAM replaced by another of the same size and modification time calls
# again: the record holds a checksum of the BAM, not only its size and time.
BAM="${GENOME_DIR}/sample1/aligned/sample1_sorted.bam"
touch -r "$BAM" "${CASE_WORK}/bam-time"
printf 'placeholdeR\n' > "$BAM"
touch -r "${CASE_WORK}/bam-time" "$BAM"
: > "$FAKE_DOCKER_LOG"
run_expect 0 gatk-same-size-bam "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1
output_lacks gatk-same-size-bam 'Output already exists'
[ "$(count ' HaplotypeCaller ')" -eq 3 ] || fail "03a reused a VCF called from a BAM replaced under the same size and time"
# Another INTERVALS calls again.
: > "$FAKE_DOCKER_LOG"
INTERVALS=chr1 run_expect 0 gatk-intervals "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1
output_lacks gatk-intervals 'Output already exists'
[ "$(count ' HaplotypeCaller ')" -eq 1 ] || fail "03a reused a whole-genome VCF for INTERVALS=chr1"
grep -q '^INTERVALS=chr1 ' "${G}/sample1.run" || fail "03a did not record INTERVALS=chr1"
# A realigned BAM (another modification time) calls again.
touch -t 209901010000 "${GENOME_DIR}/sample1/aligned/sample1_sorted.bam"
: > "$FAKE_DOCKER_LOG"
INTERVALS=chr1 run_expect 0 gatk-new-bam "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1
[ "$(count ' HaplotypeCaller ')" -eq 1 ] || fail "03a reused a VCF called from an older BAM"
# SCATTER=false after a scattered run calls again, in one process.
: > "$FAKE_DOCKER_LOG"
SCATTER=false run_expect 0 gatk-one "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1
[ "$(count ' HaplotypeCaller ')" -eq 1 ] || fail "SCATTER=false did not run one HaplotypeCaller"
if grep -q 'bcftools concat' "$FAKE_DOCKER_LOG"; then fail "SCATTER=false joined units"; fi
: > "$FAKE_DOCKER_LOG"
run_expect 1 gatk-bad-contig env INTERVALS="chr1 chr9" "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1
output_has gatk-bad-contig 'contig chr9 \(INTERVALS\) is not in the reference'
[ "$(count ' HaplotypeCaller ')" -eq 0 ] || fail "03a started GATK with an unknown contig in INTERVALS"
# A region past the end of its contig (chrUn_KI270302v1 has 2,274 bases).
: > "$FAKE_DOCKER_LOG"
run_expect 1 gatk-past-end env INTERVALS="chr1 chrUn_KI270302v1:1-3000" "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1
output_has gatk-past-end 'region chrUn_KI270302v1:1-3000 \(INTERVALS\) is not on chrUn_KI270302v1 \(1-2274\)'
[ "$(count ' HaplotypeCaller ')" -eq 0 ] || fail "03a started GATK on a region past the end of its contig"
# Two regions that overlap.
: > "$FAKE_DOCKER_LOG"
run_expect 1 freebayes-overlap env INTERVALS="chr2:1-500 chr1 chr2:400-900" "${SCRIPTS}/03b-freebayes.sh" sample1
output_has freebayes-overlap 'regions chr2:1-500 and chr2:400-900 \(INTERVALS\) overlap'
[ "$(count ' freebayes ')" -eq 0 ] || fail "03b started FreeBayes on two regions that overlap"
# scatter_beds on its own, on a 1 kb chr20: a start of 0, an end past the
# contig, an end before the start, a whole contig with one of its regions,
# and the same region twice are refused; regions that only touch are two units.
mkdir -p "${GENOME_DIR}/small"
printf 'chr20\t1000\t7\t60\t61\n' > "${GENOME_DIR}/small/ref.fasta.fai"
beds() {
  REF_FASTA="${GENOME_DIR}/small/ref.fasta" \
    bash -c '. "$1/scripts/lib/common.sh"; scatter_beds "$2" "$3"' _ "$REPO_ROOT" "${CASE_WORK}/small-units" "$1"
}
for bad in "chr20:0-2000" "chr20:0-500" "chr20:1-1001" "chr20:600-500" "chr20 chr20:10-20" "chr20:5-9 chr20:5-9"; do
  if beds "$bad" > "${CASE_WORK}/beds.out" 2>&1; then fail "scatter_beds took INTERVALS='${bad}' on a 1 kb chr20"; fi
  grep -qE 'ERROR: (region .* is not on chr20 \(1-1000\)|regions .* overlap)' "${CASE_WORK}/beds.out" \
    || fail "scatter_beds refused INTERVALS='${bad}' without saying why: $(cat "${CASE_WORK}/beds.out")"
done
TOUCH=$(beds "chr20:501-1000 chr20:1-500" | while read -r f; do cut -f1-3 "$f"; done | tr '\t\n' ':,')
[ "$TOUCH" = "chr20:0:500,chr20:500:1000," ] || fail "scatter_beds did not take two regions that only touch as two units (${TOUCH})"
# Memory: THREADS=16 would run 8 units of 8 GB at once on the fake docker's
# 16 GB; it runs 2, and says why. SCATTER_JOBS still wins.
rm -rf "$G"
: > "$FAKE_DOCKER_LOG"
THREADS=16 run_expect 0 gatk-memory "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1
output_has gatk-memory '3 unit\(s\), 2 at a time'
output_has gatk-memory 'NOTE: 2 unit\(s\) at a time, not 8: docker has 16 GB and each unit may use 8 GB'
rm -rf "$G"
THREADS=16 SCATTER_JOBS=3 run_expect 0 gatk-jobs "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1
output_has gatk-jobs '3 unit\(s\), 3 at a time'

# --- 03b FreeBayes -----------------------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 freebayes "${SCRIPTS}/03b-freebayes.sh" sample1
[ "$(count ' freebayes .*--targets /genome/sample1/vcf_freebayes/scatter/00[123]\.bed ')" -eq 3 ] || fail "03b did not run FreeBayes once per unit"
docker_log_has 'bcftools sort .*vcf_freebayes/sample1_raw\.vcf' "03b did not sort the joined raw VCF"
grep -q '^INTERVALS= scatter=true bam=.* image=' "${GENOME_DIR}/sample1/vcf_freebayes/sample1.run" 2>/dev/null || fail "03b wrote no run record beside its VCF"
: > "$FAKE_DOCKER_LOG"
run_expect 0 freebayes-again "${SCRIPTS}/03b-freebayes.sh" sample1
output_has freebayes-again 'Output already exists'
[ "$(count ' freebayes ')" -eq 0 ] || fail "a rerun of 03b called again a VCF it had finished"
# A VCF without its record (an older version of this step) is called again.
rm -f "${GENOME_DIR}/sample1/vcf_freebayes/sample1.run"
: > "$FAKE_DOCKER_LOG"
run_expect 0 freebayes-no-record "${SCRIPTS}/03b-freebayes.sh" sample1
[ "$(count ' freebayes ')" -eq 3 ] || fail "03b reused a VCF without knowing how it was called"
INTERVALS="chr2 chr1:1-100" run_expect 0 freebayes-regions "${SCRIPTS}/03b-freebayes.sh" sample1
output_has freebayes-regions '2 unit\(s\), 2 at a time'
# The units, and so the joined VCF, follow the reference, not INTERVALS.
ORDER=$(GENOME_DIR="$GENOME_DIR" bash -c '. "$1/scripts/lib/common.sh"; for f in $(scatter_beds "$2" "chr2 chr1:1-100"); do cut -f1,2 "$f"; done' \
  _ "$REPO_ROOT" "${CASE_WORK}/units" | tr '\t\n' ': ')
[ "$ORDER" = "chr1:0 chr2:0 " ] || fail "scatter_beds did not put INTERVALS in reference order (${ORDER})"

# --- 29 Mutect2, INTERVALS=genome ----------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
INTERVALS=genome run_expect 0 mutect2 "${SCRIPTS}/29-mutect2-somatic.sh" sample1
output_has mutect2 'Scattered: 3 units, 2 at a time'
[ "$(count ' Mutect2 .*--f1r2-tar-gz /genome/sample1/somatic/scatter/00[123]_f1r2\.tar\.gz .*-L /genome/sample1/somatic/scatter/00[123]\.bed')" -eq 3 ] \
  || fail "29 did not run Mutect2 once per unit with its own orientation counts"
docker_log_has ' MergeVcfs -I .*scatter/001\.vcf\.gz -I .*scatter/002\.vcf\.gz -I .*scatter/003\.vcf\.gz -O /genome/sample1/somatic/sample1_somatic_unfiltered\.vcf\.gz' \
  "29 did not merge the unit VCFs"
docker_log_has ' MergeMutectStats --stats .*001\.vcf\.gz\.stats --stats .*002\.vcf\.gz\.stats --stats .*003\.vcf\.gz\.stats -O ' "29 did not merge the unit statistics"
docker_log_has ' LearnReadOrientationModel -I .*001_f1r2\.tar\.gz -I .*002_f1r2\.tar\.gz -I .*003_f1r2\.tar\.gz ' "the orientation model did not read every unit"
# A rerun keeps the result; a BAM replaced under the same size and time calls again.
SO="${GENOME_DIR}/sample1/somatic"
{ printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n' | gzip -c
  printf '\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'; } > "${SO}/sample1_somatic_filtered.vcf.gz"
: > "$FAKE_DOCKER_LOG"
INTERVALS=genome run_expect 0 mutect2-again "${SCRIPTS}/29-mutect2-somatic.sh" sample1
output_has mutect2-again 'Output already exists'
touch -r "$BAM" "${CASE_WORK}/bam-time"
printf 'placeholdER\n' > "$BAM"
touch -r "${CASE_WORK}/bam-time" "$BAM"
: > "$FAKE_DOCKER_LOG"
INTERVALS=genome THREADS=16 run_expect 0 mutect2-same-size-bam "${SCRIPTS}/29-mutect2-somatic.sh" sample1
output_lacks mutect2-same-size-bam 'Output already exists'
[ "$(count ' Mutect2 ')" -eq 3 ] || fail "29 reused a result called from a BAM replaced under the same size and time"
output_has mutect2-same-size-bam 'Scattered: 3 units, 2 at a time'
echo "Scattered callers run one container per unit and join them."
