#!/usr/bin/env bash
# Step 22 (SURVIVOR merge, then og sv consensus-check): a consensus record is
# one event two or more callers agree on. To the synthetic set of
# tests/fixtures/sv/ (its README has the known answer: one record, at
# 10,101,001) this case adds six events. SURVIVOR 1.0.7 merges each of them
# into one record of two or three callers; the case shows that first, by
# running SURVIVOR alone on the calls the step let vote. Then:
#
# - two deletions of 100 bp that start 900 bp apart and share no base: gone;
# - a Manta deletion and a Delly duplication at the same place (SURVIVOR
#   wrote one DUP of two callers): gone;
# - two insertions of 1,000 and 300 bp, 200 bp apart: gone (sizes too far
#   apart). Both use Delly's layout (<INS>, END = POS + 1, SVLEN, CT=NtoN).
#   Delly's CT=NtoN gives ++ strands. Manta uses +- strands.
# - two insertions of 1,000 and 600 bp, 200 bp apart: kept, support 2;
# - a Manta inversion and a Delly breakend pair (BND, 3to3) of it: kept
#   (SURVIVOR's rule counts BND and INV as one type);
# - three 600 bp deletions 300 and 400 bp apart, one per caller: the first
#   and the last share no base, so the record keeps Manta's and Delly's
#   calls (support 2, SUPP_VEC 110) at Manta's coordinates, and drops
#   CNVpytor's.
#
# Red before the fix: main's step 22 wrote SURVIVOR's records as they were,
# so the consensus held 7 records, not 4, and the chain read 3 callers.
. "$(dirname "$0")/lib.sh"

T="${SAMPLE}svcc"
FIX="${REPO}/tests/fixtures/sv"
rm -rf "${GENOME_DIR:?}/${T}"
mkdir -p "${GENOME_DIR}/${T}/manta/results/variants" "${GENOME_DIR}/${T}/delly" "${GENOME_DIR}/${T}/cnvpytor" "${CASE_TMP}/in"
cp "${FIX}/manta.vcf.in" "${FIX}/delly.vcf.in" "${FIX}/cnvpytor.vcf.in" "${CASE_TMP}/in/"
rec() {  # rec FILE POS ID ALT INFO GT: one chr20 record added to the copy of FILE
  printf 'chr20\t%s\t%s\tN\t%s\t60\tPASS\t%s\tGT\t%s\n' "$2" "$3" "$4" "$5" "$6" >> "${CASE_TMP}/in/$1"
}
rec manta.vcf.in 10800000 manta_e1 '<DEL>' 'SVTYPE=DEL;END=10800100;SVLEN=-100' 0/1
rec delly.vcf.in 10800900 delly_e1 '<DEL>' 'SVTYPE=DEL;END=10801000;SVLEN=-100;CT=3to5' 0/1
rec manta.vcf.in 10900000 manta_e2 '<DEL>' 'SVTYPE=DEL;END=10903000;SVLEN=-3000' 0/1
rec delly.vcf.in 10900000 delly_e2 '<DUP>' 'SVTYPE=DUP;END=10903000;SVLEN=3000;CT=5to3' 0/1
rec manta.vcf.in 11000000 manta_e3 '<INS>' 'SVTYPE=INS;END=11000001;SVLEN=1000;CT=NtoN' 0/1
rec delly.vcf.in 11000200 delly_e3 '<INS>' 'SVTYPE=INS;END=11000201;SVLEN=300;CT=NtoN' 0/1
rec manta.vcf.in 11100000 manta_e4 '<INS>' 'SVTYPE=INS;END=11100001;SVLEN=1000;CT=NtoN' 0/1
rec delly.vcf.in 11100200 delly_e4 '<INS>' 'SVTYPE=INS;END=11100201;SVLEN=600;CT=NtoN' 0/1
rec manta.vcf.in 11200000 manta_e5 '<INV>' 'SVTYPE=INV;END=11205000;SVLEN=5000' 0/1
rec delly.vcf.in 11200100 delly_e5 'N]chr20:11205100]' 'SVTYPE=BND;CHR2=chr20;END=11205100;CT=3to3' 0/1
rec manta.vcf.in 11300000 manta_e6 '<DEL>' 'SVTYPE=DEL;END=11300600;SVLEN=-600' 0/1
rec delly.vcf.in 11300300 delly_e6 '<DEL>' 'SVTYPE=DEL;END=11300900;SVLEN=-600;CT=3to5' 0/1
rec cnvpytor.vcf.in 11300700 cnv_e6 '<DEL>' 'SVTYPE=DEL;END=11301300;SVLEN=-600;pytorP1=1e-08;pytorQ0=0;pytorPN=0' './.'
bgz() {  # bgz NAME DEST: the VCF NAME, bgzipped and indexed at GENOME_DIR/DEST
  docker run --rm -i --user "$(id -u):$(id -g)" -v "${CASE_TMP}/in:/in:ro" -v "${GENOME_DIR}:/genome" "$BCFTOOLS_IMAGE" \
    sh -c "bcftools view -Oz -o /genome/$2 /in/$1 && bcftools index -t /genome/$2"
}
bgz manta.vcf.in "${T}/manta/results/variants/diploidSV.vcf.gz"
bgz delly.vcf.in "${T}/delly/${T}_sv.vcf.gz"
bgz cnvpytor.vcf.in "${T}/cnvpytor/${T}_cnvs.vcf.gz"

run_step 22-survivor-merge.sh "$T"
check_step_exit 22-survivor-merge.sh
VCF="${T}/sv_merged/${T}_sv_consensus.vcf.gz"
check "consensus VCF is readable" vcf_ok "$VCF"

# SURVIVOR alone, on the calls the step let vote: the merges the check undoes.
IN="${T}/sv_merged/inputs"
printf '/genome/%s/%s.vcf\n' "$IN" manta "$IN" delly "$IN" cnvpytor > "${GENOME_DIR}/${IN}/raw_files.txt"
docker run --rm --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" "$SURVIVOR_IMAGE" \
  SURVIVOR merge "/genome/${IN}/raw_files.txt" 1000 2 1 1 0 50 "/genome/${IN}/raw.vcf" >/dev/null 2>&1 || true
raw() {  # POS SUPP SUPP_VEC SVTYPE of SURVIVOR's own records
  grep -v '^#' "${GENOME_DIR}/${IN}/raw.vcf" 2>/dev/null | awk -F'\t' '{
    n = split($8, i, ";"); s = v = t = "."
    for (k = 1; k <= n; k++) {
      if (i[k] ~ /^SUPP=/) s = substr(i[k], 6)
      if (i[k] ~ /^SUPP_VEC=/) v = substr(i[k], 10)
      if (i[k] ~ /^SVTYPE=/) t = substr(i[k], 8)
    }
    print $2 "\t" s "\t" v "\t" t }' | sort -n
}
echo "SURVIVOR alone (POS SUPP SUPP_VEC SVTYPE):"
raw | tee "${CASE_TMP}/raw.tsv"
inraw() { awk -F'\t' -v a="$1" -v b="$2" '$1 >= a && $1 <= b' "${CASE_TMP}/raw.tsv"; }
check_eq "SURVIVOR alone: records" "$(grep -c . "${CASE_TMP}/raw.tsv" || true)" 7
check_eq "SURVIVOR alone: the two disjoint deletions are one record" "$(cut -f2,3 <<<"$(inraw 10799000 10802000)")" $'2\t110'
check_eq "SURVIVOR alone: the deletion and the duplication are one DUP" "$(cut -f2,4 <<<"$(inraw 10899000 10904000)")" $'2\tDUP'
check_eq "SURVIVOR alone: the 1,000 and 300 bp insertions are one record" "$(cut -f2 <<<"$(inraw 10999000 11001000)")" 2
check_eq "SURVIVOR alone: the chain is one record of 3 callers" "$(cut -f2,3 <<<"$(inraw 11299000 11302000)")" $'3\t111'

echo "Consensus records (POS END SUPP SUPP_VEC SVTYPE CIPOS):"
bcf query -f '%POS\t%INFO/END\t%INFO/SUPP\t%INFO/SUPP_VEC\t%INFO/SVTYPE\t%INFO/CIPOS\n' "$VCF" 2>/dev/null | tee "${CASE_TMP}/consensus.tsv"
at() { awk -F'\t' -v a="$1" -v b="$2" '$1 >= a && $1 <= b' "${CASE_TMP}/consensus.tsv"; }
check_eq "the fixture's deletion is still one record of Manta and Delly" "$(cut -f3,4 <<<"$(at 10100999 10101001)")" $'2\t110'
check_eq "no record for the two deletions 900 bp apart" "$(at 10799000 10802000 | grep -c . || true)" 0
check_eq "no record for the deletion and the duplication" "$(at 10899000 10904000 | grep -c . || true)" 0
check_eq "no record for the 1,000 and 300 bp insertions" "$(at 10999000 11001000 | grep -c . || true)" 0
check_eq "the 1,000 and 600 bp insertions stay one record" "$(cut -f3,4,5 <<<"$(at 11099000 11101000)")" $'2\t110\tINS'
check_eq "the inversion and the breakend pair stay one record" "$(cut -f3,4,5 <<<"$(at 11199000 11201000)")" $'2\t110\tINV'
check_eq "the chain keeps Manta's and Delly's calls at Manta's coordinates" \
  "$(at 11299000 11302000)" $'11300000\t11300600\t2\t110\tDEL\t0,300'
check_eq "consensus records in all" "$(vcf_count "$VCF")" 4
check "the header states the check" has '^##og_sv_consensus_check=' "$(bcf view -h "$VCF" 2>/dev/null)"
check "the log counts what the check did" has \
  'og sv consensus-check: 7 records, 4 written; 3 one event as merged, 4 rechecked: split into 1 records of 2\+ callers, 3 dropped with fewer than 2 callers left; 0 malformed records dropped, 0 unresolved member votes dropped' \
  "$(cat "$STEP_LOG")"

# The round-one review calls go through the same step and bcftools sort.
T="${SAMPLE}svccr"
rm -rf "${GENOME_DIR:?}/${T}"
mkdir -p "${GENOME_DIR}/${T}/manta/results/variants" "${GENOME_DIR}/${T}/delly" "${GENOME_DIR}/${T}/cnvpytor"
cp "${FIX}"/consensus-review-{manta,delly,cnvpytor}.vcf.in "${CASE_TMP}/in/"
bgz consensus-review-manta.vcf.in "${T}/manta/results/variants/diploidSV.vcf.gz"
bgz consensus-review-delly.vcf.in "${T}/delly/${T}_sv.vcf.gz"
bgz consensus-review-cnvpytor.vcf.in "${T}/cnvpytor/${T}_cnvs.vcf.gz"
run_step 22-survivor-merge.sh "$T"
check_step_exit 22-survivor-merge.sh
VCF="${T}/sv_merged/${T}_sv_consensus.vcf.gz"
bcf query -f '%POS\t%ID\t%ALT\t%QUAL\t%INFO/SUPP_VEC[\t%ID]\n' "$VCF" >"${CASE_TMP}/review.tsv"
review() { awk -F'\t' -v m="$1" 'int($1 / 1000000) == m' "${CASE_TMP}/review.tsv"; }
check_eq "short insertion from the same caller joins the matching insertion" \
  "$(cut -f1,2,4,5 <<<"$(review 1)")" $'1000200\tsmall\t55\t110'
check_eq "the older call keeps its own orientation" "$(review 2)" ""
check_eq "breakends to different contigs cannot agree" "$(review 3)" ""
check_eq "two same-caller deletions leave the compatible pair" \
  "$(cut -f1,2,4,5 <<<"$(review 4)")" $'4000250\tmiddle\t58\t110'
check_eq "the discarded source ID and QUAL leave with its call" \
  "$(cut -f1,2,4,5,6,7 <<<"$(review 5)")" $'5000050\tmemberD\t59\t110\t.\tmemberD'
check_eq "explicit INV and TRA cannot agree" "$(review 6)" ""
check_eq "the split translocation keeps the retained mate, QUAL, and callers" \
  "$(cut -f1,2,3,4,5 <<<"$(review 7)")" $'7000000\tkept-bnd\tN]chr21:7005000]\t55\t110'
check_eq "consensus records across the review probes" "$(vcf_count "$VCF")" 5
finish
