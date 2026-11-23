#!/usr/bin/env bash
# Step 16b measures each coding exon of the catalog genes in the same mosdepth
# pass as the genes, on both paths, with the real mosdepth image, and a gene
# read on its body with an unread coding exon reads "exon not read".
#
# A made-up BAM on a 3,000-base chr20: 50-base single reads start every 10 bases
# from 900 to 2650, 12 at each start, except the starts that would cover
# 1500-1520, so those 20 bases have no read and every other base of 950-2650 has
# 12 or more. The catalog has two genes and the exon table three exons:
#   TEST_GAP  1000-2000, 980 of 1000 bases at 10x (98%): read on its body.
#             TEST_GAP:exon1 1100-1140 is read, TEST_GAP:exon2 1500-1520 is not.
#   TEST_OK   2200-2600, every base at 10x; TEST_OK:exon1 2300-2400 is read.
# An exon row of a gene the catalog lacks (OTHER:exon1) is not measured.
# scripts/16b-mosdepth.sh (CALLABILITY_GENES, CALLABILITY_EXONS) and the
# MOSDEPTH module, run alone through a one-line pipeline with both tables, must
# write the same gene table (two genes, no exon row) and the same exon table
# (three exons: 40, 0 and 100 bases at 10x). A separate mosdepth run over the
# exons alone gives the same counts, so the shared pass changes no number.
# og summary collect on step 16b's folder lists TEST_GAP as "exon not read"
# naming TEST_GAP:exon2 and TEST_OK as read. Controls: without the exon table
# (an exon table of other genes) step 16b writes no exon table and TEST_GAP
# reads read; MOSDEPTH given the gene table alone writes the gene table alone.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

S="${SAMPLE}-exons"
SD="${GENOME_DIR}/${S}"
rm -rf "${SD:?}"
mkdir -p "${SD}/aligned"

CAT_DIR="${CASE_TMP}/assets"
rm -rf "$CAT_DIR"
mkdir -p "$CAT_DIR"
CAT="${CAT_DIR}/callability_genes.tsv"
EXONS="${CAT_DIR}/callability_exons.tsv"
{
  printf '# a test catalog\n'
  printf 'gene\thgnc_id\tchrom\tstart\tend\treasons\towner\n'
  printf 'TEST_GAP\tHGNC:0\tchr20\t1000\t2000\tsegdup\t.\n'
  printf 'TEST_OK\tHGNC:1\tchr20\t2200\t2600\tsegdup\t.\n'
} > "$CAT"
{
  printf '# a test exon table\n'
  printf 'exon\tgene\ttranscript\tchrom\tstart\tend\n'
  printf 'TEST_GAP:exon1\tTEST_GAP\tT1\tchr20\t1100\t1140\n'
  printf 'TEST_GAP:exon2\tTEST_GAP\tT1\tchr20\t1500\t1520\n'
  printf 'TEST_OK:exon1\tTEST_OK\tT2\tchr20\t2300\t2400\n'
  printf 'OTHER:exon1\tOTHER\tT3\tchr20\t100\t200\n'
} > "$EXONS"
OTHER_EXONS="${CASE_TMP}/other_exons.tsv"
printf 'exon\tgene\ttranscript\tchrom\tstart\tend\nOTHER:exon1\tOTHER\tT3\tchr20\t100\t200\n' > "$OTHER_EXONS"

# The reads: SAM positions are 1-based, so a read at POS p covers bases p-1 to
# p+48 (0-based). Starts whose read would reach a base of 1500-1519 are left out.
SAM="${SD}/aligned/in.sam"
s50=$(printf 'A%.0s' $(seq 50)); q50=$(printf 'I%.0s' $(seq 50))
{
  printf '@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr20\tLN:3000\n'
  for start in $(seq 900 10 2650); do
    if [ "$start" -gt 1450 ] && [ "$start" -lt 1520 ]; then continue; fi
    for i in $(seq 12); do
      printf 'r%s_%s\t0\tchr20\t%s\t60\t50M\t*\t0\t0\t%s\t%s\n' "$start" "$i" "$((start + 1))" "$s50" "$q50"
    done
  done
} > "$SAM"
docker run --rm --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -w /genome "$SAMTOOLS_IMAGE" sh -c \
  "samtools sort -o '${S}/aligned/${S}_sorted.bam' '${S}/aligned/in.sam' && samtools index '${S}/aligned/${S}_sorted.bam'"
BAM="${SD}/aligned/${S}_sorted.bam"
check "the made-up BAM and its index exist" test -s "${BAM}.bai"
check_eq "the made-up BAM holds 2,040 reads" "$(sam view -c "${S}/aligned/${S}_sorted.bam" 2>/dev/null)" 2040

# at10 TABLE REGION: the region's bases at 10x or more (column 10X of the table)
at10() { gzip -dc "$1" 2>/dev/null | awk -F'\t' -v g="$2" '$4 == g {print $5}'; }
# rows TABLE: the table's data rows
rows() { gzip -dc "$1" 2>/dev/null | grep -vc '^#'; }

# --- the input: mosdepth over the exons alone ---------------------------------------------
awk -F'\t' -v OFS='\t' '!/^#/ && $1 != "exon" && $2 != "OTHER" {print $4, $5, $6, $1}' "$EXONS" > "${SD}/exons.bed"
docker run --rm --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -w "/genome/${S}" "$MOSDEPTH_IMAGE" \
  mosdepth --by exons.bed --mapq 20 --no-per-base --thresholds 10 alone "aligned/${S}_sorted.bam"
A="${SD}/alone.thresholds.bed.gz"
check_eq "input: TEST_GAP:exon1 has 40 bases at 10x" "$(at10 "$A" TEST_GAP:exon1)" 40
check_eq "input: TEST_GAP:exon2 has no base at 10x" "$(at10 "$A" TEST_GAP:exon2)" 0
check_eq "input: TEST_OK:exon1 has 100 bases at 10x" "$(at10 "$A" TEST_OK:exon1)" 100

# --- step 16b ---------------------------------------------------------------------------
CALLABILITY_GENES="$CAT" CALLABILITY_EXONS="$EXONS" THREADS=2 run_step 16b-mosdepth.sh "$S"
check_step_exit 16b-mosdepth.sh
T="${SD}/mosdepth/${S}_callability.bed.gz"
E="${SD}/mosdepth/${S}_callability_exons.bed.gz"
echo "step 16b's gene table:"; gzip -dc "$T" 2>/dev/null
echo "step 16b's exon table:"; gzip -dc "$E" 2>/dev/null
check_eq "bash: the gene table has the two genes and no exon row" "$(rows "$T")" 2
check_eq "bash: TEST_GAP has 980 bases at 10x" "$(at10 "$T" TEST_GAP)" 980
check_eq "bash: TEST_OK has 400 bases at 10x" "$(at10 "$T" TEST_OK)" 400
check_eq "bash: the exon table has the three exons of catalog genes" "$(rows "$E")" 3
for x in TEST_GAP:exon1 TEST_GAP:exon2 TEST_OK:exon1; do
  check_eq "bash: ${x} has the bases at 10x of a pass over the exons alone" "$(at10 "$E" "$x")" "$(at10 "$A" "$x")"
done

# reading SAMPLE DIR: "GENE STATUS WORST_EXON" of each gene og summary collect
# lists as not read, or "none", from the sample folder DIR with the test catalog
# as the assets.
reading() {
  local t
  t=$(mktemp -d "${CASE_TMP}/reading.XXXXXX")
  mkdir -p "${t}/assets"
  cp "$CAT" "${t}/assets/callability_genes.tsv"
  "${REPO}/bin/og" summary collect --sample "$1" --sample-dir "$2" --assets "${t}/assets" \
    --out "${t}/summary.json" --only a > /dev/null
  jq -r '.sections.callability.data.not_read
    | if length == 0 then "none" else .[] | "\(.gene) \(.status) \(.worst_exon.exon // "-")" end' "${t}/summary.json"
}
R=$(reading "$S" "$SD" 2>&1); echo "og summary collect, not read: ${R}"
check_eq "bash: TEST_GAP reads exon not read, naming TEST_GAP:exon2; TEST_OK is read" "$R" "TEST_GAP exon not read TEST_GAP:exon2"

# --- control: an exon table that names none of the catalog genes ------------------------
# Its own sample folder, with a copy of the BAM, so step 16b's tables above stay.
C="${S}-ctl"
CD="${GENOME_DIR}/${C}"
rm -rf "${CD:?}"
mkdir -p "${CD}/aligned"
cp "$BAM" "${CD}/aligned/${C}_sorted.bam"
cp "${BAM}.bai" "${CD}/aligned/${C}_sorted.bam.bai"
CALLABILITY_GENES="$CAT" CALLABILITY_EXONS="$OTHER_EXONS" THREADS=2 run_step 16b-mosdepth.sh "$C"
check_step_exit "16b-mosdepth.sh (control)"
check_eq "control: step 16b wrote the gene table" "$(at10 "${CD}/mosdepth/${C}_callability.bed.gz" TEST_GAP)" 980
check "control: with no exon of a catalog gene, step 16b writes no exon table" \
  test ! -e "${CD}/mosdepth/${C}_callability_exons.bed.gz"
R0=$(reading "$C" "$CD" 2>&1); echo "control, not read: ${R0}"
check_eq "control: without the exon table TEST_GAP reads read" "$R0" "none"

# --- the MOSDEPTH module ----------------------------------------------------------------
NP="${CASE_TMP}/project"
nf_mosdepth() { # NAME TABLES: MOSDEPTH run alone with the tables (a Groovy list)
  rm -rf "$NP" "${CASE_TMP}/out-$1" "${CASE_TMP}/work-$1"
  mkdir -p "$NP"
  cat > "${NP}/main.nf" <<NF
include { MOSDEPTH } from '${REPO}/modules/local/mosdepth/main'

workflow {
    MOSDEPTH(channel.of([[id: '${S}'], file('${BAM}'), file('${BAM}.bai')]), $2)
}
NF
  cat > "${NP}/nextflow.config" <<NFCONF
params.outdir = '${CASE_TMP}/out-$1'
params.publish_dir_mode = 'copy'
docker.enabled = true
docker.runOptions = '--user $(id -u):$(id -g) --network none'
NFCONF
  echo "+ nextflow run (MOSDEPTH, $1)"
  (cd "$CASE_TMP" && nextflow run "${NP}/main.nf" -ansi-log false -c "${REPO}/conf/containers.config" \
     -work-dir "${CASE_TMP}/work-$1") > "${CASE_TMP}/nf-$1.log" 2>&1
  local rc=$?
  check_eq "nextflow run (MOSDEPTH, $1) exits 0" "$rc" 0
  [ "$rc" = 0 ] || grep -E 'ERROR|Error executing process|Command error' -A 8 "${CASE_TMP}/nf-$1.log" | head -40
}
nf_mosdepth both "[file('${CAT}'), file('${EXONS}')]"
NT="${CASE_TMP}/out-both/${S}/coverage/${S}_callability.bed.gz"
NE="${CASE_TMP}/out-both/${S}/coverage/${S}_callability_exons.bed.gz"
check_eq "Nextflow: the gene table has the two genes and no exon row" "$(rows "$NT")" 2
check_eq "Nextflow: the exon table has the three exons of catalog genes" "$(rows "$NE")" 3
for x in TEST_GAP:exon1 TEST_GAP:exon2 TEST_OK:exon1; do
  check_eq "Nextflow: ${x} has the bases at 10x of a pass over the exons alone" "$(at10 "$NE" "$x")" "$(at10 "$A" "$x")"
done
check_eq "Nextflow: TEST_GAP has 980 bases at 10x" "$(at10 "$NT" TEST_GAP)" 980

check_eq "MOSDEPTH's gene table is step 16b's" "$(gzip -dc "$NT" 2>/dev/null | md5sum)" "$(gzip -dc "$T" 2>/dev/null | md5sum)"
check_eq "MOSDEPTH's exon table is step 16b's" "$(gzip -dc "$NE" 2>/dev/null | md5sum)" "$(gzip -dc "$E" 2>/dev/null | md5sum)"

nf_mosdepth genes "file('${CAT}')"
check "control: MOSDEPTH given the gene table alone writes it" test -s "${CASE_TMP}/out-genes/${S}/coverage/${S}_callability.bed.gz"
check "control: and writes no exon table" test ! -e "${CASE_TMP}/out-genes/${S}/coverage/${S}_callability_exons.bed.gz"

rm -rf "${SD:?}" "${CD:?}"
finish
