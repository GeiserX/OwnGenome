#!/usr/bin/env bash
# Step 16b's gene callability counts the two mates of a pair once where they
# overlap, on both paths, with the real mosdepth image.
#
# A made-up BAM on a 3,000-base chr20 holds three 10-base genes:
#   TEST_OVERLAP  5 pairs whose mates both cover the gene: 5 fragments, so no
#                 base reaches 10x. With --fast-mode mosdepth counts each base
#                 twice, 10x, and the gene read as read.
#   TEST_PAIRS    10 such pairs: 10 fragments, every base at 10x (control: the
#                 pass still counts each pair once, not zero times).
#   TEST_DEL      10 single reads with a 10-base deletion over the gene: no read
#                 has a base there. With --fast-mode the deletion counted 10x.
# scripts/16b-mosdepth.sh (CALLABILITY_GENES) and the MOSDEPTH module, run alone
# through a one-line pipeline, must both write 0, 10 and 0 bases at 10x. A run
# of the same image with --fast-mode on the same BAM must write 10, 10 and 10,
# so the BAM does hold what the gene pass must not count twice.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

S="${SAMPLE}-ovl"
SD="${GENOME_DIR}/${S}"
rm -rf "${SD:?}"
mkdir -p "${SD}/aligned"

CAT="${CASE_TMP}/callability_genes.tsv"
{
  printf 'gene\thgnc_id\tchrom\tstart\tend\treasons\towner\n'
  printf 'TEST_OVERLAP\tHGNC:0\tchr20\t1000\t1010\tsegdup\t.\n'
  printf 'TEST_PAIRS\tHGNC:1\tchr20\t1500\t1510\tsegdup\t.\n'
  printf 'TEST_DEL\tHGNC:2\tchr20\t2000\t2010\tsegdup\t.\n'
} > "$CAT"

# The reads: 60-base mates 20 bases apart (insert 80), both over the gene, and
# 80-base reads with 40M10D40M whose deletion is the gene's 10 bases.
SAM="${SD}/aligned/in.sam"
s60=$(printf 'A%.0s' $(seq 60)); q60=$(printf 'I%.0s' $(seq 60))
s80=$(printf 'A%.0s' $(seq 80)); q80=$(printf 'I%.0s' $(seq 80))
{
  printf '@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr20\tLN:3000\n'
  for i in 1 2 3 4 5; do
    printf 'ov%s\t99\tchr20\t961\t60\t60M\t=\t981\t80\t%s\t%s\n' "$i" "$s60" "$q60"
    printf 'ov%s\t147\tchr20\t981\t60\t60M\t=\t961\t-80\t%s\t%s\n' "$i" "$s60" "$q60"
  done
  for i in $(seq 10); do
    printf 'pr%s\t99\tchr20\t1461\t60\t60M\t=\t1481\t80\t%s\t%s\n' "$i" "$s60" "$q60"
    printf 'pr%s\t147\tchr20\t1481\t60\t60M\t=\t1461\t-80\t%s\t%s\n' "$i" "$s60" "$q60"
  done
  for i in $(seq 10); do
    printf 'dl%s\t0\tchr20\t1961\t60\t40M10D40M\t*\t0\t0\t%s\t%s\n' "$i" "$s80" "$q80"
  done
} > "$SAM"
docker run --rm --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -w /genome "$SAMTOOLS_IMAGE" sh -c \
  "samtools sort -o '${S}/aligned/${S}_sorted.bam' '${S}/aligned/in.sam' && samtools index '${S}/aligned/${S}_sorted.bam'"
BAM="${SD}/aligned/${S}_sorted.bam"
check "the made-up BAM and its index exist" test -s "${BAM}.bai"
check_eq "the made-up BAM holds 40 reads" "$(sam view -c "${S}/aligned/${S}_sorted.bam" 2>/dev/null)" 40

# at10 TABLE GENE: the gene's bases at 10x or more (column 10X of the table)
at10() { gzip -dc "$1" 2>/dev/null | awk -F'\t' -v g="$2" '$4 == g {print $5}'; }

# --- the input: --fast-mode counts the overlap and the deletion -----------------------
awk -F'\t' -v OFS='\t' 'NR > 1 {print $3, $4, $5, $1}' "$CAT" > "${SD}/fast.bed"
docker run --rm --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -w "/genome/${S}" "$MOSDEPTH_IMAGE" \
  mosdepth --by fast.bed --mapq 20 --fast-mode --no-per-base --thresholds 10 fast "aligned/${S}_sorted.bam"
F="${SD}/fast.thresholds.bed.gz"
check_eq "input: --fast-mode counts TEST_OVERLAP's 5 pairs as 10x" "$(at10 "$F" TEST_OVERLAP)" 10
check_eq "input: --fast-mode counts TEST_DEL's deletion as 10x" "$(at10 "$F" TEST_DEL)" 10
check_eq "input: --fast-mode counts TEST_PAIRS at 10x" "$(at10 "$F" TEST_PAIRS)" 10

# --- step 16b ---------------------------------------------------------------------------
CALLABILITY_GENES="$CAT" THREADS=2 run_step 16b-mosdepth.sh "$S"
check_step_exit 16b-mosdepth.sh
T="${SD}/mosdepth/${S}_callability.bed.gz"
echo "step 16b's table:"; gzip -dc "$T" 2>/dev/null
check_eq "bash: TEST_OVERLAP's 5 overlapping pairs reach 10x on no base" "$(at10 "$T" TEST_OVERLAP)" 0
check_eq "bash: TEST_PAIRS's 10 pairs reach 10x on all 10 bases" "$(at10 "$T" TEST_PAIRS)" 10
check_eq "bash: TEST_DEL's deleted bases reach 10x on no base" "$(at10 "$T" TEST_DEL)" 0

# --- the MOSDEPTH module ----------------------------------------------------------------
NP="${CASE_TMP}/project"
rm -rf "$NP" "${CASE_TMP}/out" "${CASE_TMP}/work"
mkdir -p "$NP"
cat > "${NP}/main.nf" <<NF
include { MOSDEPTH } from '${REPO}/modules/local/mosdepth/main'

workflow {
    MOSDEPTH(channel.of([[id: '${S}'], file('${BAM}'), file('${BAM}.bai')]), file('${CAT}'))
}
NF
cat > "${NP}/nextflow.config" <<NFCONF
params.outdir = '${CASE_TMP}/out'
params.publish_dir_mode = 'copy'
docker.enabled = true
docker.runOptions = '--user $(id -u):$(id -g) --network none'
NFCONF
echo "+ nextflow run (MOSDEPTH)"
(cd "$CASE_TMP" && nextflow run "${NP}/main.nf" -ansi-log false -c "${REPO}/conf/containers.config" \
   -work-dir "${CASE_TMP}/work") > "${CASE_TMP}/nf.log" 2>&1
RC=$?
check_eq "nextflow run (MOSDEPTH) exits 0" "$RC" 0
[ "$RC" = 0 ] || grep -E 'ERROR|Error executing process|Command error' -A 8 "${CASE_TMP}/nf.log" | head -40
N="${CASE_TMP}/out/${S}/coverage/${S}_callability.bed.gz"
check_eq "Nextflow: TEST_OVERLAP's 5 overlapping pairs reach 10x on no base" "$(at10 "$N" TEST_OVERLAP)" 0
check_eq "Nextflow: TEST_PAIRS's 10 pairs reach 10x on all 10 bases" "$(at10 "$N" TEST_PAIRS)" 10
check_eq "Nextflow: TEST_DEL's deleted bases reach 10x on no base" "$(at10 "$N" TEST_DEL)" 0
check_eq "MOSDEPTH's table is step 16b's" "$(gzip -dc "$N" 2>/dev/null | md5sum)" "$(gzip -dc "$T" 2>/dev/null | md5sum)"

rm -rf "${SD:?}"
finish
