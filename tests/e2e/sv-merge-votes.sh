#!/usr/bin/env bash
# Step 22 (SURVIVOR merge): a caller votes only with the calls it makes. On the
# synthetic set of tests/fixtures/sv/ (its README has the known answer), three
# Manta-only deletions each have a second record beside them that is no call:
# a PASS Delly record genotyped 0/0, a Delly call whose sample filter (FT) is
# LowQual, and a CNVpytor segment with e-value 0.5. None of them may reach the
# consensus, while a CNVpytor call that passes the quality rule and has
# genotype ./. still pairs with Manta's call (support 2, SUPP_VEC 101). Then the
# same calls come from CNVpytor's table instead of its VCF: the e-value and q0
# of the table decide the same way, also when an earlier version of the step
# left the table's VCF without them.
#
# Red before the fix: main's step 22 kept every PASS record, so the three
# Manta-only deletions read SUPP=2 and the consensus held 5 records, not 2.
. "$(dirname "$0")/lib.sh"

T="${SAMPLE}svvote"
FIX="${REPO}/tests/fixtures/sv"
rm -rf "${GENOME_DIR:?}/${T}"
mkdir -p "${GENOME_DIR}/${T}/manta/results/variants" "${GENOME_DIR}/${T}/delly" "${GENOME_DIR}/${T}/cnvpytor" "${CASE_TMP}/in"
cp "${FIX}/manta.vcf.in" "${FIX}/delly.vcf.in" "${FIX}/cnvpytor.vcf.in" "${CASE_TMP}/in/"
# The positive control, a 10 kb deletion both Manta and CNVpytor call, CNVpytor
# with genotype ./. as the VCF made from its table writes it.
printf 'chr20\t10700001\tmanta_6\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=10710001;SVLEN=-10000\tGT\t0/1\n' >> "${CASE_TMP}/in/manta.vcf.in"
printf 'chr20\t10700201\tcnv_3\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=10710201;SVLEN=-10000;pytorP1=1e-08;pytorQ0=0.01;pytorPN=0\tGT\t./.\n' \
  >> "${CASE_TMP}/in/cnvpytor.vcf.in"
bgz() {  # bgz NAME DEST: the VCF NAME, bgzipped and indexed at GENOME_DIR/DEST
  docker run --rm -i --user "$(id -u):$(id -g)" -v "${CASE_TMP}/in:/in:ro" -v "${GENOME_DIR}:/genome" "$BCFTOOLS_IMAGE" \
    sh -c "bcftools view -Oz -o /genome/$2 /in/$1 && bcftools index -t /genome/$2"
}
bgz manta.vcf.in "${T}/manta/results/variants/diploidSV.vcf.gz"
bgz delly.vcf.in "${T}/delly/${T}_sv.vcf.gz"
bgz cnvpytor.vcf.in "${T}/cnvpytor/${T}_cnvs.vcf.gz"

VCF="${T}/sv_merged/${T}_sv_consensus.vcf.gz"
IN="${GENOME_DIR}/${T}/sv_merged/inputs"
consensus() {  # POS SUPP SUPP_VEC of every consensus record
  bcf query -f '%POS\t%INFO/SUPP\t%INFO/SUPP_VEC\n' "$VCF" 2>/dev/null | tee "${CASE_TMP}/consensus.tsv"
}
at() {  # at FROM TO: the consensus records that start in FROM..TO
  awk -F'\t' -v a="$1" -v b="$2" '$1 >= a && $1 <= b' "${CASE_TMP}/consensus.tsv"
}
ids() { grep -v '^#' "$1" 2>/dev/null | cut -f3 | tr '\n' ' '; }
votes() {  # the checks of one run; $1 names the CNVpytor input
  echo "Consensus records (POS SUPP SUPP_VEC):"
  consensus
  check_eq "$1: one record for the deletion Manta and Delly both call" "$(at 10100999 10101001)" $'10101001\t2\t110'
  check_eq "$1: no record for the Manta call beside Delly's 0/0 record" "$(at 10399000 10405000 | grep -c . || true)" 0
  check_eq "$1: no record for the Manta call beside Delly's FT=LowQual call" "$(at 10499000 10504000 | grep -c . || true)" 0
  check_eq "$1: no record for the Manta call beside a CNVpytor e-value of 0.5" "$(at 10599000 10611000 | grep -c . || true)" 0
  check_eq "$1: the ./. CNVpytor call that passes its quality rule still votes" "$(cut -f2,3 <<<"$(at 10699000 10711000)")" $'2\t101'
  check_eq "$1: consensus records in all" "$(vcf_count "$VCF")" 2
  check_eq "$1: the Delly records that vote" "$(ids "${IN}/delly.vcf")" "delly_1 delly_2 "
  check "$1: the log counts Delly's left-out records" \
    has 'og sv survivor-prep delly: 4 records, 2 vote; no vote: 0 FILTER not PASS, 1 sample FT not PASS, 1 genotyped reference' "$(cat "$STEP_LOG")"
  check "$1: the log counts the CNVpytor call below the rule" has 'og sv survivor-prep cnvpytor: .*, 1 below CNVpytor.s quality rule' "$(cat "$STEP_LOG")"
}

run_step 22-survivor-merge.sh "$T"
check_step_exit 22-survivor-merge.sh
votes "CNVpytor VCF"
check_eq "CNVpytor VCF: the CNVpytor records that vote" "$(ids "${IN}/cnvpytor.vcf")" "cnv_1 cnv_3 "

# The same three CNVpytor calls from its table (step 18's _cnvs.txt: type,
# region, size, level, e-val1..4, q0, pN, dG) when step 18 left no VCF.
rm -f "${GENOME_DIR}/${T}/cnvpytor/${T}_cnvs.vcf.gz" "${GENOME_DIR}/${T}/cnvpytor/${T}_cnvs.vcf.gz.tbi"
printf 'deletion\tchr20:%s\t10001\t0.5\t%s\t0\t0\t0\t%s\t0\t1000\n' \
  10300000-10310000 1e-08 0 10600001-10610001 0.5 0 10700201-10710201 1e-08 0.01 \
  > "${GENOME_DIR}/${T}/cnvpytor/${T}_cnvs.txt"
run_step 22-survivor-merge.sh "$T"
check_step_exit 22-survivor-merge.sh
check "the log says the table was turned into a VCF" has 'CNVpytor: converting' "$(cat "$STEP_LOG")"
check "the table's e-value reaches INFO" has 'pytorP1=0.5;pytorQ0=0;pytorPN=0' \
  "$(bcf view -H "${T}/cnvpytor/${T}_cnvs.vcf.gz" 2>/dev/null)"
votes "CNVpytor table"
check_eq "CNVpytor table: the CNVpytor calls that vote" "$(grep -v '^#' "${IN}/cnvpytor.vcf" 2>/dev/null | cut -f2 | tr '\n' ' ')" "10300000 10700201 "

# An earlier version of step 22 left the table's VCF without the e-value tags.
# Kept as it is, none of its calls would vote; the step makes it again.
# shellcheck disable=SC2016  # $1 expands inside the container
docker run --rm --network none --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" "$BCFTOOLS_IMAGE" bash -euo pipefail -c '
  f="/genome/$1/cnvpytor/$1_cnvs.vcf.gz"
  bcftools view "$f" | sed -e "/^##INFO=<ID=pytor/d" -e "s/;pytorP1=[^\t]*//" | bcftools view -Oz -o "$f.old" -
  mv "$f.old" "$f" && bcftools index -f -t "$f"' _ "$T"
check "the old table VCF has no e-value" lacks 'pytorP1' "$(bcf view "${T}/cnvpytor/${T}_cnvs.vcf.gz" 2>/dev/null)"
run_step 22-survivor-merge.sh "$T"
check_step_exit 22-survivor-merge.sh
check "the log says the old table VCF is made again" has 'has no e-value \(pytorP1\); making it again' "$(cat "$STEP_LOG")"
votes "old table VCF"

finish
