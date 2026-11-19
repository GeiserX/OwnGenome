#!/usr/bin/env bash
# Steps 03 and 03e take the sample's sex:
#   03 male    DeepVariant gets --haploid_contigs=chrX,chrY and the PAR BED;
#              female and no sex do not; any other value stops with exit 2.
#   03         writes a VCF and a gVCF with their indexes under temporary
#              names, honours THREADS, DV_MEM, VCF_OUT_DIR and keeps its
#              intermediate files in the sample directory, removed after; a
#              run that leaves no gVCF fails and keeps the earlier VCF.
#   03e male   Clair3 gets --gender=male and the PAR BED; CLAIR3_MODEL wins.
#   03e ont    the model comes from the BAM's @RG DS basecall_model= or from
#              BASECALL_MODEL, through og longread clair3-model run for real
#              (bin/og) against a fake /opt/models; no basecaller, one Clair3
#              has no model for, or a model the image lacks stops the step
#              with the image's models listed, before Clair3 runs.
#   02b ont    the basecaller model (BASECALL_MODEL, or the unaligned BAM's
#              own @RG DS) goes into the aligned BAM's @RG as DS.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
mkdir -p "${GENOME_DIR}/sample1/aligned_longread"
cp "${GENOME_DIR}/sample1/aligned/"* "${GENOME_DIR}/sample1/aligned_longread/"

use_output_hook
cat > "${CASE_WORK}/tools-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
shift   # the image
args=" $* "
# eq NAME: the value of --NAME=value in the command.
eq() { sed -n "s/.* --$1=\([^ ]*\) .*/\1/p" <<<"$args"; }
bgzf() {
  mkdir -p "$(dirname "$1")"
  { printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\nchr1\t1000\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n# %s\n' "${FAKE_TAG:-1}" | gzip -c
    printf '\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'; } > "$1"
  printf 'TBI\001' > "${1}.tbi"
}
case "$args" in
  *"/opt/deepvariant/bin/run_deepvariant "*)
    [ "${FAKE_DV:-ok}" = fail ] && exit 1
    v=$(host_path "$(eq output_vcf)")
    bgzf "$v"
    [ -d "$(host_path "$(eq intermediate_results_dir)")" ] || { echo "fake DeepVariant: no intermediate directory" >&2; exit 1; }
    : > "$(host_path "$(eq intermediate_results_dir)")/make_examples.tfrecord.gz"
    [ "${FAKE_DV:-ok}" = nogvcf ] || bgzf "$(host_path "$(eq output_gvcf)")" ;;
  *"/opt/bin/run_clair3.sh "*)
    bgzf "$(host_path "$(eq output)")/merge_output.vcf.gz" ;;
  " /pgp-bin/og longread clair3-model "*)
    # The real og, with the image's /opt/models played by CASE_WORK/models.
    run=()
    for a in "$@"; do
      case "$a" in
        /pgp-bin/*) run+=("${REPO_ROOT}/bin/${a#/pgp-bin/}") ;;
        /genome*) run+=("$(host_path "$a")") ;;
        /opt/models*) run+=("${CASE_WORK}/models${a#/opt/models}") ;;
        *) run+=("$a") ;;
      esac
    done
    "${run[@]}" | sed "s#^${CASE_WORK}/models#/opt/models#" ;;
  *" samtools index "*)
    for a in "$@"; do last=$a; done
    : > "$(host_path "$last").bai" ;;
  *" bcftools stats "*)
    printf 'SN\t0\tnumber of records:\t1\n' ;;
  *) exec "${CASE_WORK}/hook-outputs" image "$@" ;;
esac
HOOK
chmod +x "${CASE_WORK}/tools-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/tools-hook"

V="${GENOME_DIR}/sample1/vcf"
leftovers() { find "${GENOME_DIR}/sample1" -mindepth 1 \( -name '*.part.*' -o -name deepvariant_tmp \) 2>/dev/null; }
dv_line() { awk '/run_deepvariant/' "$FAKE_DOCKER_LOG" | tail -n 1; }

# --- 03 male, THREADS, DV_MEM ----------------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 dv-male env THREADS=2 DV_MEM=20g "${SCRIPTS}/03-deepvariant.sh" sample1 male
L=$(dv_line)
# The docker log quotes each word with printf %q, which writes the comma as \,
grep -qE -- '--haploid_contigs=chrX\\?,chrY ' <<<"$L" || fail "DeepVariant call lacks --haploid_contigs=chrX,chrY: ${L}"
for want in '--cpus 2 ' '--memory 20g ' '--num_shards=2 ' \
            '--par_regions_bed=/pgp/par_grch38.bed ' '/assets/par_grch38.bed:/pgp/par_grch38.bed:ro ' \
            '--output_vcf=/genome/sample1/vcf/sample1.part.vcf.gz ' \
            '--output_gvcf=/genome/sample1/vcf/sample1.part.g.vcf.gz ' \
            '--intermediate_results_dir=/genome/sample1/vcf/deepvariant_tmp '; do
  grep -qF -- "$want" <<<"$L" || fail "DeepVariant call lacks '${want}': ${L}"
done
for f in sample1.vcf.gz sample1.vcf.gz.tbi sample1.g.vcf.gz sample1.g.vcf.gz.tbi; do
  [ -f "${V}/${f}" ] || fail "no ${V}/${f} after step 03"
done
bgzf_complete() { [ "$(tail -c 28 "$1" | od -An -v -tx1 | tr -d ' \n')" = 1f8b08040000000000ff0600424302001b0003000000000000000000 ]; }
bgzf_complete "${V}/sample1.g.vcf.gz" || fail "the gVCF is not the file DeepVariant wrote"
[ -z "$(leftovers)" ] || fail "files left after step 03: $(leftovers)"
output_has dv-male 'chrX and chrY haploid outside the PARs'

# --- 03 female, no sex, a wrong sex ---------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 dv-female "${SCRIPTS}/03-deepvariant.sh" sample1 female
L=$(dv_line)
if grep -q -- '--haploid_contigs' <<<"$L"; then fail "female: DeepVariant got --haploid_contigs"; fi
grep -qF -- '--cpus 8 ' <<<"$L" || fail "the default THREADS is not 8: ${L}"
grep -qF -- '--memory 32g ' <<<"$L" || fail "the default DV_MEM is not 32g: ${L}"
: > "$FAKE_DOCKER_LOG"
run_expect 0 dv-nosex "${SCRIPTS}/03-deepvariant.sh" sample1
L=$(dv_line)
if grep -q -- '--haploid_contigs' <<<"$L"; then fail "no sex: DeepVariant got --haploid_contigs"; fi
output_has dv-nosex 'chrX and chrY are called diploid'
run_expect 2 dv-badsex "${SCRIPTS}/03-deepvariant.sh" sample1 m

# --- VCF_OUT_DIR keeps the primary VCF ------------------------------------------
before=$(od -An -tx1 < "${V}/sample1.vcf.gz" | tr -d ' \n')
run_expect 0 dv-outdir env ALIGN_DIR=aligned_longread VCF_OUT_DIR=vcf_longread FAKE_TAG=lr "${SCRIPTS}/03-deepvariant.sh" sample1
[ -f "${GENOME_DIR}/sample1/vcf_longread/sample1.g.vcf.gz" ] || fail "VCF_OUT_DIR: no gVCF in vcf_longread/"
[ "$(od -An -tx1 < "${V}/sample1.vcf.gz" | tr -d ' \n')" = "$before" ] || fail "VCF_OUT_DIR: the primary VCF changed"
run_expect 2 dv-outdir-bad env VCF_OUT_DIR=../x "${SCRIPTS}/03-deepvariant.sh" sample1

# --- DeepVariant leaves no gVCF: the step fails, the old VCF stays ---------------
before=$(od -An -tx1 < "${V}/sample1.vcf.gz" | tr -d ' \n')
run_rc dv-nogvcf env FAKE_DV=nogvcf FAKE_TAG=new "${SCRIPTS}/03-deepvariant.sh" sample1
expect_rc dv-nogvcf 1
[ "$(od -An -tx1 < "${V}/sample1.vcf.gz" | tr -d ' \n')" = "$before" ] || fail "a run without a gVCF replaced the VCF"
[ -z "$(leftovers)" ] || fail "files left after the failed run: $(leftovers)"

# --- 03e Clair3 ---------------------------------------------------------------------
"${REPO_ROOT}/bin/og" version > /dev/null || fail "og does not run here: build it into bin/.og (rust/og.sh, the Guard fake-docker job does)"
for m in r1041_e82_400bps_sup_v500 r1041_e82_400bps_sup_v520 r1041_e82_400bps_hac_v410 r941_prom_sup_g5014 hifi_revio; do
  mkdir -p "${CASE_WORK}/models/${m}"
done
# mkbam FILE HEADER: a BAM that holds only HEADER (gzip reads as BGZF here).
mkbam() {
  local n=${#2}
  { printf 'BAM\001'
    # shellcheck disable=SC2059  # the format is the four length bytes
    printf "$(printf '\\x%02x\\x%02x\\x%02x\\x%02x' $((n & 255)) $((n >> 8 & 255)) $((n >> 16 & 255)) $((n >> 24 & 255)))"
    printf '%s' "$2"
    printf '\000\000\000\000'; } | gzip -c > "$1"
}
LRBAM="${GENOME_DIR}/sample1/aligned_longread/sample1_sorted.bam"
mkbam "$LRBAM" $'@HD\tVN:1.6\n@RG\tID:sample1\tSM:sample1\tPL:ONT\n'
clair3_ran() { grep -q 'run_clair3' "$FAKE_DOCKER_LOG"; }
: > "$FAKE_DOCKER_LOG"
run_expect 0 clair3-male env PLATFORM=ont CLAIR3_MODEL=/opt/models/r1041_e82_400bps_sup_v520 "${SCRIPTS}/03e-clair3.sh" sample1 male
docker_log_has 'run_clair3\.sh .*--model_path=/opt/models/r1041_e82_400bps_sup_v520 .*--gender=male --par_regions_bed=/pgp/par_grch38\.bed' \
  "Clair3 did not get the model override, --gender=male and the PAR BED"
: > "$FAKE_DOCKER_LOG"
run_expect 0 clair3-nosex env PLATFORM=hifi "${SCRIPTS}/03e-clair3.sh" sample1
if awk '/run_clair3/ && /--gender/ { bad = 1 } END { exit !bad }' "$FAKE_DOCKER_LOG"; then
  fail "Clair3 got --gender with no sex given"
fi
docker_log_has 'run_clair3\.sh .*--model_path=/opt/models/hifi_revio ' "Clair3 HiFi model path changed"
run_expect 2 clair3-badsex env PLATFORM=ont "${SCRIPTS}/03e-clair3.sh" sample1 M

# --- 03e ONT: the model follows the basecaller -------------------------------------
mkbam "$LRBAM" $'@HD\tVN:1.6\n@RG\tID:sample1\tSM:sample1\tPL:ONT\tDS:basecall_model=dna_r10.4.1_e8.2_400bps_hac@v4.1.0\n'
: > "$FAKE_DOCKER_LOG"
run_expect 0 clair3-from-rg env PLATFORM=ont "${SCRIPTS}/03e-clair3.sh" sample1
docker_log_has 'run_clair3\.sh .*--model_path=/opt/models/r1041_e82_400bps_hac_v410 ' \
  "Clair3 did not get the model of the BAM's basecaller (hac v4.1.0)"
: > "$FAKE_DOCKER_LOG"
run_expect 0 clair3-basecall-env env PLATFORM=ont BASECALL_MODEL=dna_r9.4.1_e8_sup@v3.3 "${SCRIPTS}/03e-clair3.sh" sample1
docker_log_has 'run_clair3\.sh .*--model_path=/opt/models/r941_prom_sup_g5014 ' \
  "BASECALL_MODEL did not choose the R9.4.1 SUP model over the BAM's read group"

# No basecaller anywhere, a fast basecaller, a model the image lacks, a
# CLAIR3_MODEL the image lacks: each stops before Clair3 runs.
mkbam "$LRBAM" $'@HD\tVN:1.6\n@RG\tID:sample1\tSM:sample1\tPL:ONT\n'
: > "$FAKE_DOCKER_LOG"
run_expect 1 clair3-no-basecaller env PLATFORM=ont "${SCRIPTS}/03e-clair3.sh" sample1
output_has clair3-no-basecaller 'names no basecaller model.*BASECALL_MODEL.*CLAIR3_MODEL'
output_has clair3-no-basecaller 'models in .*: hifi_revio r1041_e82_400bps_hac_v410 r1041_e82_400bps_sup_v500 r1041_e82_400bps_sup_v520 r941_prom_sup_g5014'
! clair3_ran || fail "Clair3 ran with no basecaller model known"
: > "$FAKE_DOCKER_LOG"
run_expect 1 clair3-fast env PLATFORM=ont BASECALL_MODEL=dna_r10.4.1_e8.2_400bps_fast@v5.0.0 "${SCRIPTS}/03e-clair3.sh" sample1
output_has clair3-fast 'Clair3 has no model for the basecaller model dna_r10\.4\.1_e8\.2_400bps_fast@v5\.0\.0'
! clair3_ran || fail "Clair3 ran on a fast basecaller"
: > "$FAKE_DOCKER_LOG"
run_expect 1 clair3-not-in-image env PLATFORM=ont BASECALL_MODEL=dna_r10.4.1_e8.2_400bps_sup@v4.3.0 "${SCRIPTS}/03e-clair3.sh" sample1
output_has clair3-not-in-image 'r1041_e82_400bps_sup_v430 .*is not in the image'
! clair3_ran || fail "Clair3 ran with a model the image lacks"
: > "$FAKE_DOCKER_LOG"
run_expect 1 clair3-bad-override env PLATFORM=ont CLAIR3_MODEL=/opt/models/no_such_model "${SCRIPTS}/03e-clair3.sh" sample1
output_has clair3-bad-override 'no_such_model \(set by CLAIR3_MODEL\) is not in the image'
! clair3_ran || fail "Clair3 ran with a CLAIR3_MODEL the image lacks"

# --- 02b ONT carries the basecaller model into the aligned BAM ----------------------
mkdir -p "${GENOME_DIR}/sample1/fastq"
printf 'placeholder\n' | gzip -c > "${GENOME_DIR}/sample1/fastq/sample1.fastq.gz"
: > "$FAKE_DOCKER_LOG"
run_expect 0 longread-basecall-env env PLATFORM=ont BASECALL_MODEL=dna_r10.4.1_e8.2_400bps_sup@v5.0.0 "${SCRIPTS}/02b-alignment-longread.sh" sample1
docker_log_has 'minimap2 .*LB:sample1[^ ]*tDS:basecall_model=dna_r10\.4\.1_e8\.2_400bps_sup@v5\.0\.0 ' \
  "step 02b did not put BASECALL_MODEL into the read group"
mkbam "${GENOME_DIR}/sample1/fastq/sample1.bam" $'@HD\tVN:1.6\n@RG\tID:run1\tPL:ONT\tDS:runid=r1 basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.2.0 modbase_models=x\n'
: > "$FAKE_DOCKER_LOG"
run_expect 0 longread-ubam env PLATFORM=ont INPUT="${GENOME_DIR}/sample1/fastq/sample1.bam" "${SCRIPTS}/02b-alignment-longread.sh" sample1
docker_log_has 'minimap2 .*LB:sample1[^ ]*tDS:basecall_model=dna_r10\.4\.1_e8\.2_400bps_sup@v5\.2\.0 ' \
  "step 02b did not carry the unaligned BAM's basecaller model into the read group"
: > "$FAKE_DOCKER_LOG"
run_expect 0 longread-hifi env PLATFORM=hifi INPUT="${GENOME_DIR}/sample1/fastq/sample1.bam" BASECALL_MODEL=dna_r10.4.1_e8.2_400bps_sup@v5.0.0 "${SCRIPTS}/02b-alignment-longread.sh" sample1
if grep -q 'DS:basecall_model' "$FAKE_DOCKER_LOG"; then fail "step 02b wrote a basecaller model for HiFi reads"; fi
run_expect 2 longread-bad-basecaller env PLATFORM=ont BASECALL_MODEL=$'sup\tx' "${SCRIPTS}/02b-alignment-longread.sh" sample1
