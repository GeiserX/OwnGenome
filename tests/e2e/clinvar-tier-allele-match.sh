#!/usr/bin/env bash
# Step 23's ClinVar tier takes a VEP record only when its allele is a step 6
# hit, not every record at a hit's position. A synthetic sample (no fixture
# file is read) carries, on chr22:
#   30000037 AGC>A and 30000042 T>TGGG  the two parts of a combined hit
#   30000100 C>T                        a plain hit
#   30000200 T>A and T>C                one position, two alleles (a split
#                                       T>A,C call); only T>A is a hit
# The step 6 hits are the combined allele (30000038 GCCAT>CATGGG, a position
# neither part has, with PGP_JOINED_FROM naming the parts), 30000100 C>T and
# 30000200 T>A. The tier must hold the two parts, C>T and T>A: four records,
# and not T>C. Every record is common and MODIFIER, so no other tier keeps
# one, and the summary must not list T>C either.
. "$(dirname "$0")/lib.sh"

T="${SAMPLE}am"
rm -rf "${GENOME_DIR:?}/${T}"
mkdir -p "${GENOME_DIR}/${T}/vep" "${GENOME_DIR}/${T}/clinvar"

CSQ_FMT='Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|CANONICAL|MANE_SELECT|MAX_AF'
csq() { echo "CSQ=${1}|intergenic_variant|MODIFIER||||||||0.3"; }
row() { printf 'chr22\t%s\t.\t%s\t%s\t50\tPASS\t%s\tGT:GQ:DP:AD\t0/1:50:30:15,15\n' "$1" "$2" "$3" "$(csq "$4")"; }
{
  printf '##fileformat=VCFv4.2\n'
  printf '##FILTER=<ID=PASS,Description="All filters passed">\n'
  printf '##contig=<ID=chr22,length=50818468>\n'
  printf '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: %s">\n' "$CSQ_FMT"
  printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
  printf '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
  printf '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
  printf '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$T"
  row 30000037 AGC A -
  row 30000042 T TGGG GGG
  row 30000100 C T T
  row 30000200 T A A
  row 30000200 T C C
} > "${GENOME_DIR}/${T}/vep/${T}_vep.vcf"

{
  printf '##fileformat=VCFv4.2\n'
  printf '##contig=<ID=chr22,length=50818468>\n'
  printf '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="ClinVar significance">\n'
  printf '##INFO=<ID=PGP_JOINED_FROM,Number=1,Type=String,Description="Step 06: the combined allele of these two neighbouring sample records, CHROM:POS:REF:ALT joined by |">\n'
  printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$T"
  printf 'chr22\t30000038\t9900003\tGCCAT\tCATGGG\t.\tPASS\tCLNSIG=Pathogenic;PGP_JOINED_FROM=chr22:30000037:AGC:A|chr22:30000042:T:TGGG\tGT\t0/1\n'
  printf 'chr22\t30000100\t9900009\tC\tT\t.\tPASS\tCLNSIG=Pathogenic\tGT\t0/1\n'
  printf 'chr22\t30000200\t9900011\tT\tA\t.\tPASS\tCLNSIG=Likely_pathogenic\tGT\t1/.\n'
} > "${GENOME_DIR}/${T}/clinvar/${T}_clinvar_hits.vcf"

run_step 23-clinical-filter.sh "$T"
check_step_exit 23-clinical-filter.sh
check "the log names the step 6 hits as the tier's source" has 'ClinVar tier: step 6 hits' "$(cat "$STEP_LOG")"
TIER="${T}/clinical/${T}_clinvar_pathogenic.vcf.gz"
check "the ClinVar tier VCF is readable" vcf_ok "$TIER"
ALLELES=$(bcf query -f '%CHROM:%POS:%REF:%ALT\n' "$TIER" 2>/dev/null)
echo "tier alleles: $(tr '\n' ' ' <<< "$ALLELES")"
has_allele() { grep -cx -- "$1" <<< "$ALLELES"; }
check_eq "the classified allele chr22:30000200 T>A is in the tier" "$(has_allele chr22:30000200:T:A)" 1
check_eq "the other allele at that position, T>C, is not" "$(has_allele chr22:30000200:T:C)" 0
check_eq "the plain hit chr22:30000100 C>T is in the tier" "$(has_allele chr22:30000100:C:T)" 1
check_eq "the combined hit's first part chr22:30000037 AGC>A is in the tier" "$(has_allele chr22:30000037:AGC:A)" 1
check_eq "the combined hit's second part chr22:30000042 T>TGGG is in the tier" "$(has_allele chr22:30000042:T:TGGG)" 1
check_eq "the tier holds exactly those four records" "$(vcf_count "$TIER")" 4
check "the log counts four ClinVar records" has 'ClinVar pathogenic/LP: +4 ' "$(cat "$STEP_LOG")"
SUMMARY="${GENOME_DIR}/${T}/clinical/${T}_clinical_summary.tsv"
check "the summary exists" test -s "$SUMMARY"
check_eq "the summary has no row for chr22:30000200 T>C" \
  "$(awk -F'\t' '$1 == "chr22" && $2 == 30000200 && $3 == "T" && $4 == "C"' "$SUMMARY" | wc -l | tr -d ' ')" 0
check_eq "the summary lists the four tier records" "$(awk 'NR > 1' "$SUMMARY" | wc -l | tr -d ' ')" 4

finish
