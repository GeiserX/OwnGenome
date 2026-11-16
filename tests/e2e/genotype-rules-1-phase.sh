#!/usr/bin/env bash
# Step 42's genotype rules (assets/genotype_rules.tsv) on synthetic data, with
# the real bcftools, WhatsHap and Python images, and the Nextflow mirror
# (the made-up chr7, gVCFs and reads are og dev fixture tas2r38):
#   - scripts/ci/check-rule-alleles.sh --online: every rule site matches
#     dbSNP's GRCh38 record, and the dbSNP table pinned in the script is what
#     dbSNP's VCF of the pinned build holds today; controls: a rules copy with
#     rs1726866's ref and alt swapped fails, naming it, and a copy of the
#     script whose pinned table differs from dbSNP fails the online compare;
#   - TAS2R38 on a made-up chr7 of GRCh38's length (N everywhere except a 2.6 kb window around the
#     gene, its three sites carrying the GRCh38 bases), one gVCF heterozygous at
#     the three sites, four samples:
#       GRSPAN    read pairs that span the three sites, from PAV and AVI copies
#                 -> WhatsHap phases them in one block: PAV/AVI, read_backed
#       GRSWAP    the same gVCF, read pairs from G,A,C and C,G,T copies
#                 -> other_haplotype, so the result rests on the reads
#       GRNOSPAN  the same gVCF, single reads that each cover one site
#                 -> unphased: no read or read pair joins the sites
#       GRREF     a gVCF that is one reference block (GQ 42), no BAM
#                 -> every site homozygous reference with the block's GQ;
#                    the reference haplotype AAI is not one the rule reads
#     and the traits table keeps one row per table row plus the APOE line;
#   - nextflow run main.nf --tools traits on GRSPAN (VCF, gVCF and BAM) writes
#     the same rules table, traits table and genotypes as the bash step, byte
#     for byte; control: the comparison tells GRSPAN's rules table from
#     GRNOSPAN's.
# Writes four sample folders and a reference folder under GENOME_DIR, all
# removed at the end; no other case reads them.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
CHECK="${REPO}/scripts/ci/check-rule-alleles.sh"
RULES="${REPO}/assets/genotype_rules.tsv"
TABLE="${REPO}/assets/traits.tsv"

# --- the dbSNP orientation check, online ---------------------------------------------
bash "$CHECK" --online > "${CASE_TMP}/dbsnp.log" 2>&1
RC=$?
cat "${CASE_TMP}/dbsnp.log"
check_eq "check-rule-alleles.sh --online exits 0" "$RC" 0
check "the pinned dbSNP table is what dbSNP holds today" has '^ok: the pinned table is what dbSNP' "$(cat "${CASE_TMP}/dbsnp.log")"
awk -F'\t' -v OFS='\t' '$1 == "site" && $3 == "rs1726866" { t = $6; $6 = $7; $7 = t } { print }' "$RULES" > "${CASE_TMP}/swapped.tsv"
check_eq "control: the swapped copy differs from the rules in one line" \
  "$(diff "$RULES" "${CASE_TMP}/swapped.tsv" | grep -c '^>' || true)" 1
bash "$CHECK" --online "${CASE_TMP}/swapped.tsv" > "${CASE_TMP}/swapped.log" 2>&1
check "control: the swapped copy fails (exit $?), naming rs1726866" has 'FAIL: rs1726866: ref and alt are swapped' \
  "$(cat "${CASE_TMP}/swapped.log")"
mkdir -p "${CASE_TMP}/root/scripts/ci"
cp "${REPO}/versions.env" "${CASE_TMP}/root/"
sed 's/^#: rs7412\tchr19\t44908822\tC\tT$/#: rs7412\tchr19\t44908822\tC\tG/' "$CHECK" > "${CASE_TMP}/root/scripts/ci/check-rule-alleles.sh"
check "control: the script copy's pinned table differs in one line" \
  test "$(diff "$CHECK" "${CASE_TMP}/root/scripts/ci/check-rule-alleles.sh" | grep -c '^>')" -eq 1
bash "${CASE_TMP}/root/scripts/ci/check-rule-alleles.sh" --online "$RULES" > "${CASE_TMP}/table.log" 2>&1
check "control: a pinned table that differs from dbSNP fails the online compare" \
  has '^FAIL: the pinned table differs from dbSNP' "$(cat "${CASE_TMP}/table.log")"

# --- synthetic chr7, reads and gVCFs --------------------------------------------------
RD="${G}/reference-genotype-rules"
rm -rf "$RD"
mkdir -p "$RD"
for s in GRSPAN GRSWAP GRNOSPAN GRREF; do rm -rf "${G:?}/${s}"; mkdir -p "${G}/${s}/vcf" "${G}/${s}/aligned"; done
"${REPO}/bin/og" dev fixture tas2r38 "$RD" "$G"
check "the synthetic chr7 was written" test -s "${RD}/chr7.fa"
sam faidx reference-genotype-rules/chr7.fa
check_eq "chr7.fa.fai: one contig of GRCh38 chr7's length" "$(cut -f1,2 "${RD}/chr7.fa.fai" 2>/dev/null | tr '\t' ' ')" "chr7 159345973"
for s in GRSPAN GRSWAP GRNOSPAN GRREF; do
  in_genome "$BCFTOOLS_IMAGE" bash -c "bcftools view -Oz -o '${s}/vcf/${s}.g.vcf.gz' '${s}/vcf/${s}.g.vcf' && bcftools index -t '${s}/vcf/${s}.g.vcf.gz'"
  check "${s}: the gVCF is indexed" test -s "${G}/${s}/vcf/${s}.g.vcf.gz.tbi"
done
for s in GRSPAN GRSWAP GRNOSPAN; do
  in_genome "$BCFTOOLS_IMAGE" bash -c "bcftools view -Oz -o '${s}/vcf/${s}.vcf.gz' '${s}/vcf/${s}.vcf' && bcftools index -t '${s}/vcf/${s}.vcf.gz'"
  in_genome "$SAMTOOLS_IMAGE" bash -c "samtools sort -o '${s}/aligned/${s}_sorted.bam' '${s}/aligned/${s}.sam' && samtools index '${s}/aligned/${s}_sorted.bam'"
  check "${s}: the BAM is indexed" test -s "${G}/${s}/aligned/${s}_sorted.bam.bai"
done
check_eq "GRSPAN: 80 reads in the BAM" "$(sam view -c GRSPAN/aligned/GRSPAN_sorted.bam 2>/dev/null)" 80
check_eq "GRSPAN: the gVCF genotypes at the three sites" \
  "$(bcf query -r chr7:141972804,chr7:141972905,chr7:141973545 -f '[%GT]\n' GRSPAN/vcf/GRSPAN.g.vcf.gz 2>/dev/null | tr '\n' ' ')" "0/1 0/1 0/1 "

# --- the bash step ---------------------------------------------------------------------
export REF_FASTA="${RD}/chr7.fa"
rule() { awk -F'\t' -v r="$2" -v k="$3" 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
  $1 == r { print $c[k]; exit }' "${G}/$1/traits/$1_genotype_rules.tsv" 2>/dev/null; }
ROWS=$(grep -v '^#' "$TABLE" | awk 'NR > 1' | wc -l | tr -d ' ')
for s in GRSPAN GRSWAP GRNOSPAN GRREF; do
  run_step 42-traits.sh "$s"
  check_step_exit 42-traits.sh
  cp "$STEP_LOG" "${CASE_TMP}/${s}.log"
  check_eq "${s}: one traits row per table row plus the APOE line" \
    "$(awk 'NR > 1' "${G}/${s}/traits/${s}_traits.tsv" 2>/dev/null | wc -l | tr -d ' ')" "$((ROWS + 1))"
  check_eq "${s}: one rules line per rule" \
    "$(awk 'NR > 1' "${G}/${s}/traits/${s}_genotype_rules.tsv" 2>/dev/null | wc -l | tr -d ' ')" \
    "$(awk -F'\t' '$1 == "rule"' "$RULES" | wc -l | tr -d ' ')"
  check_eq "${s}: APOE, whose sites are not on the synthetic reference, is not covered" "$(rule "$s" APOE status)" not_covered
done
unset REF_FASTA
for s in GRSPAN GRSWAP GRNOSPAN; do
  check_eq "${s}: the three TAS2R38 genotypes, reference allele first" "$(rule "$s" TAS2R38 genotypes)" "C/G,G/A,T/C"
  check_eq "${s}: the TAS2R38 traits rows read het" \
    "$(awk -F'\t' '$2 == "TAS2R38" { printf "%s:%s ", $6, $7 }' "${G}/${s}/traits/${s}_traits.tsv" 2>/dev/null)" \
    "C/G:called G/A:called T/C:called "
done
check_eq "GRSPAN: WhatsHap phased the pairs: PAV/AVI" "$(rule GRSPAN TAS2R38 result) $(rule GRSPAN TAS2R38 status)" "PAV/AVI called"
check_eq "GRSPAN: phase read_backed" "$(rule GRSPAN TAS2R38 phase)" read_backed
check_eq "GRSPAN: the GQ of the three records" "$(rule GRSPAN TAS2R38 gq)" "35,35,35"
check_eq "GRSWAP: the same genotypes phased into G,A,C and C,G,T: other_haplotype" \
  "$(rule GRSWAP TAS2R38 status) $(rule GRSWAP TAS2R38 phase)" "other_haplotype read_backed"
check "GRSWAP: the reading names both copies" has 'G,A,C and C,G,T' "$(rule GRSWAP TAS2R38 reading)"
check_eq "GRNOSPAN: no read joins the sites: unphased" "$(rule GRNOSPAN TAS2R38 status) $(rule GRNOSPAN TAS2R38 result)" "unphased ."
check "GRNOSPAN: the reading says no read or read pair covered them" has 'no read or read pair covered them together' \
  "$(rule GRNOSPAN TAS2R38 reading)"
check_eq "GRREF: reference blocks only: homozygous reference at the three sites" "$(rule GRREF TAS2R38 genotypes)" "C/C,G/G,T/T"
check_eq "GRREF: each site with the block's GQ" "$(rule GRREF TAS2R38 gq)" "42,42,42"
check_eq "GRREF: the reference haplotype AAI is not one the rule reads" \
  "$(rule GRREF TAS2R38 status) $(rule GRREF TAS2R38 phase)" "other_haplotype not_needed"
check_eq "GRREF: the TAS2R38 traits rows read homozygous reference" \
  "$(awk -F'\t' '$2 == "TAS2R38" { printf "%s:%s ", $6, $7 }' "${G}/GRREF/traits/GRREF_traits.tsv" 2>/dev/null)" \
  "C/C:called G/G:called T/T:called "
check "GRREF: the step says there is no BAM to phase from" has 'No BAM' "$(cat "${CASE_TMP}/GRREF.log")"
echo "TAS2R38 lines:"
for s in GRSPAN GRSWAP GRNOSPAN GRREF; do
  echo "  ${s}: $(rule "$s" TAS2R38 result) $(rule "$s" TAS2R38 status) $(rule "$s" TAS2R38 phase)"
done

# --- the Nextflow mirror ---------------------------------------------------------------
if command -v nextflow >/dev/null; then
  NF_OUT="${G}/nf-genotype-rules"
  SHEET="${CASE_TMP}/samplesheet.csv"
  S="${G}/GRSPAN"
  printf 'sample,vcf,vcf_index,gvcf,gvcf_index,bam,bam_index\nGRSPAN,%s,%s,%s,%s,%s,%s\n' \
    "${S}/vcf/GRSPAN.vcf.gz" "${S}/vcf/GRSPAN.vcf.gz.tbi" "${S}/vcf/GRSPAN.g.vcf.gz" "${S}/vcf/GRSPAN.g.vcf.gz.tbi" \
    "${S}/aligned/GRSPAN_sorted.bam" "${S}/aligned/GRSPAN_sorted.bam.bai" > "$SHEET"
  (
    cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
      -work-dir "${E2E_WORK}/nf-work-genotype-rules" \
      --input "$SHEET" --reference "${RD}/chr7.fa" --tools traits --sex_check warn \
      --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB
  ) 2>&1 | tee "$STEP_LOG"
  NF_RC=${PIPESTATUS[0]}
  check_eq "nextflow run --tools traits on GRSPAN exits 0" "$NF_RC" 0
  check "the run had a TRAITS_PHASE task" has 'TRAITS_PHASE' "$(cat "$STEP_LOG")"
  for f in GRSPAN_genotype_rules.tsv GRSPAN_traits.tsv GRSPAN_traits_genotypes.tsv; do
    check "Nextflow published traits/${f}" test -s "${NF_OUT}/GRSPAN/traits/${f}"
    check "traits/${f}: Nextflow and the bash step agree byte for byte" \
      cmp -s "${NF_OUT}/GRSPAN/traits/${f}" "${S}/traits/${f}"
  done
  check "control: the comparison tells GRSPAN's rules table from GRNOSPAN's" \
    test "$(cmp -s "${NF_OUT}/GRSPAN/traits/GRSPAN_genotype_rules.tsv" "${G}/GRNOSPAN/traits/GRNOSPAN_genotype_rules.tsv"; echo $?)" -ne 0
  in_genome "$BCFTOOLS_IMAGE" rm -rf nf-genotype-rules
else
  fail "nextflow is not on PATH"
fi

in_genome "$BCFTOOLS_IMAGE" rm -rf GRSPAN GRSWAP GRNOSPAN GRREF reference-genotype-rules
finish
