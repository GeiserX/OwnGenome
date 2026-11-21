#!/usr/bin/env bash
# Run all six synthetic score tracks through both annotation paths and feed
# the Nextflow output into CLINICAL_FILTER. Only the high scores select hits.
# shellcheck source=lib.sh
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
WORK=$(mktemp -d "${CASE_TMP}/run.XXXXXX") || { fail "cannot create the case directory"; finish; }
GENOME_DIR="${WORK}/genome"
S=VCFANNO_TEST
P="${WORK}/project"
VEP="${S}/vep/${S}_vep.vcf"
ANNOT="${GENOME_DIR}/annotations"

(
  set -e
  mkdir -p "${GENOME_DIR}/${S}/vep" "$ANNOT" "$P/bin/.og"
  cp "${REPO}/bin/og" "$P/bin/"
  cp "${REPO}"/bin/.og/og-linux-* "$P/bin/.og/"

  # Use synthetic CSQ entries with LOW impact so that an impact tier cannot
  # pick any record by itself.
  {
    cat <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=CSQ,Number=.,Type=String,Description="VEP annotations. Format: Allele|Consequence|IMPACT|SYMBOL|Feature|MAX_AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
VCF
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$S"
    printf 'chr1\t%s\t%s\tA\t%s\t.\tPASS\tCSQ=%s|%s|LOW|TESTGENE|TX1|0.001\tGT\t0/1\n' \
      101 cadd_snv C C synonymous_variant \
      202 cadd_indel AT T intron_variant \
      303 revel G G missense_variant \
      404 alphamissense T T missense_variant \
      505 cadd_below AAG AG intron_variant \
      550 missense_below G G missense_variant \
      606 spliceai_snv C C intron_variant \
      707 spliceai_indel AT T intron_variant
  } > "${GENOME_DIR}/${VEP}"

  printf '#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n1\t101\tA\tC\t3.1\t25\n' \
    > "$ANNOT/whole_genome_SNVs.tsv"
  printf '#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n1\t202\tA\tAT\t3.2\t26\n1\t505\tA\tAAG\t1.1\t10\n' \
    > "$ANNOT/gnomad.genomes.r4.0.indel.tsv"
  printf '%s\n' \
    $'#chr\tpos\tref\talt\tREVEL' \
    $'chr1\t303\tA\tG\t0.9' \
    $'chr1\t404\tA\tT\t0.1' \
    $'chr1\t550\tA\tG\t0.2' > "$ANNOT/revel_grch38.tsv"
  printf '%s\n' \
    $'#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class' \
    $'chr1\t303\tA\tG\thg38\tPTEST\tTX1\tA1G\t0.1\tlikely_benign' \
    $'chr1\t404\tA\tT\thg38\tPTEST\tTX1\tA2T\t0.9\tlikely_pathogenic' \
    $'chr1\t550\tA\tG\thg38\tPTEST\tTX1\tA3G\t0.2\tlikely_benign' > "$ANNOT/AlphaMissense_hg38.tsv"

  splice_track() {
    printf '##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n'
    printf '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAI prediction">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    printf '1\t%s\t.\tA\t%s\t.\tPASS\tSpliceAI=%s|TESTGENE|0.3|0|0|0|1|0|0|0\n' "$1" "$2" "$2"
  }
  splice_track 606 C > "$ANNOT/spliceai_scores.masked.snv.hg38.vcf"
  splice_track 707 AT > "$ANNOT/spliceai_scores.masked.indel.hg38.vcf"

  docker run --rm --user "$(id -u):$(id -g)" --network none \
    -v "${GENOME_DIR}:/genome" -w /genome "$GATK_IMAGE" bash -ec '
      for file in annotations/*.tsv; do
        bgzip -c "$file" > "$file.gz"
        tabix -s 1 -b 2 -e 2 "$file.gz"
      done
      for file in annotations/*.vcf "$1"; do
        bgzip -c "$file" > "$file.gz"
        tabix -p vcf "$file.gz"
      done
    ' bash "$VEP"
)
setup_rc=$?

[ "$setup_rc" -eq 0 ] || { fail "could not prepare the synthetic VCF and score tables"; finish; }

cat > "$P/main.nf" <<NF
include { VCFANNO; VCFANNO_INDEX } from '${REPO}/modules/local/vcfanno/main'
include { CLINICAL_FILTER } from '${REPO}/modules/local/clinical_filter/main'

workflow {
    def input = file('${GENOME_DIR}/${VEP}.gz')
    def tracks = [
        'whole_genome_SNVs.tsv.gz',
        'gnomad.genomes.r4.0.indel.tsv.gz',
        'spliceai_scores.masked.snv.hg38.vcf.gz',
        'spliceai_scores.masked.indel.hg38.vcf.gz',
        'revel_grch38.tsv.gz',
        'AlphaMissense_hg38.tsv.gz'
    ].collect { name -> file("${ANNOT}/\$name") }
    VCFANNO(channel.of([[id: '${S}'], input, file("\${input}.tbi")]),
        tracks[0], file("\${tracks[0]}.tbi"), tracks[1], file("\${tracks[1]}.tbi"),
        tracks[2], file("\${tracks[2]}.tbi"), tracks[3], file("\${tracks[3]}.tbi"),
        tracks[4], file("\${tracks[4]}.tbi"), tracks[5], file("\${tracks[5]}.tbi"))
    VCFANNO_INDEX(VCFANNO.out.vcf)
    CLINICAL_FILTER(VCFANNO_INDEX.out.vcf.join(VCFANNO_INDEX.out.vcf_index))
}
NF
cat > "$P/nextflow.config" <<NFCONF
params.outdir = '${GENOME_DIR}/nextflow'
params.publish_dir_mode = 'copy'
docker.enabled = true
docker.runOptions = '--user $(id -u):$(id -g) --network none'
NFCONF

(
  cd "$WORK" && nextflow run "$P/main.nf" -ansi-log false \
    -c "${REPO}/conf/containers.config" -work-dir "${WORK}/work"
) 2>&1 | tee "$STEP_LOG"
NF_RC=${PIPESTATUS[0]}
check_eq "VCFANNO, VCFANNO_INDEX and CLINICAL_FILTER exit 0" "$NF_RC" 0
[ "$NF_RC" -eq 0 ] || finish

clinical="nextflow/${S}/clinical/${S}"
check_eq "CLINICAL_FILTER selects the six high-score records and excludes the two low-score records" \
  "$(bcf query -f '%ID\n' "${clinical}_clinical.vcf.gz" | paste -sd, -)" \
  cadd_snv,cadd_indel,revel,alphamissense,spliceai_snv,spliceai_indel
check_eq "CLINICAL_FILTER reports one row per high-score hit" \
  "$(awk 'END {print NR - 1}' "${GENOME_DIR}/${clinical}_clinical_summary.tsv")" 6

run_step 30-vcfanno.sh "$S"
check_step_exit 30-vcfanno.sh

info_types() {
  bcf view -h "$1" | sed -nE 's/^##INFO=<ID=(CADD_PHRED(_indel)?|SpliceAI(_indel)?|REVEL|AM_pathogenicity|AM_class),Number=[^,]+,Type=([^,]+),.*/\1=\4/p' | LC_ALL=C sort
}
expected_types=$(printf '%s\n' AM_class=String AM_pathogenicity=Float CADD_PHRED=Float \
  CADD_PHRED_indel=Float REVEL=Float SpliceAI=String SpliceAI_indel=String)
check_eq "Nextflow score INFO types" \
  "$(info_types "nextflow/${S}/vep/${S}_annotated.vcf.gz")" "$expected_types"
check_eq "step 30 score INFO types" \
  "$(info_types "${S}/vep/${S}_annotated.vcf.gz")" "$expected_types"

finish
