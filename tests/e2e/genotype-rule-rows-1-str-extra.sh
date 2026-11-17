#!/usr/bin/env bash
# The repeats STRchive leaves out (assets/str_extra_loci.json) with the real tools:
#   1. STR_EXTRA_LOCI=yes setup.sh --strchive-data appends them to STRchive's
#      catalog (82 + 2 entries, the asset's sha256 in catalog.stamp); the same
#      folder without the opt-in converts again and is back to 82 entries.
#   2. ExpansionHunter 5.0.0 (step 09) on made-up reads, against a made-up
#      reference whose chr7 and chr2 are N except a window around each tract,
#      the GRCh38 tract at its GRCh38 place (CFTR TG11 T7, UGT1A1 A(TA)6TAA):
#        EXA  CFTR TG12-T5 on one copy and TG11-T7 on the other, UGT1A1 7 and
#             8 TA units (TA6/TA7, *1/*28)
#             -> CFTR_TG 11/12, CFTR_T 5 and 7, UGT1A1 7/8 (the counts are compared
#                as a set: ExpansionHunter printed CFTR_T as 7/5, not sorted)
#        EXB  TG11-T7 on both copies, UGT1A1 7 units on both
#             -> CFTR_TG 11/11, CFTR_T 7/7 (no 5T: the red control), UGT1A1 7/7
#      and step 09 prints both loci as not graded, with the asset's reason.
# What this cannot show: reads of a real person at these tracts, or a reading
# of the counts (the 5T line and *28 are not written by any step yet).
# Writes two sample folders and a reference folder under GENOME_DIR, all removed
# at the end; no other case reads them.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
EXTRA="${REPO}/assets/str_extra_loci.json"
SETUP="${CASE_TMP}/setup"
rm -rf "$SETUP"
mkdir -p "$SETUP"
CAT="${SETUP}/reference/strchive-${STRCHIVE_VERSION}/catalog.json"
count() { jq length "$1" 2>/dev/null || echo 0; }

# --- 1. setup, with and without the opt-in ------------------------------------------------
STR_EXTRA_LOCI=yes "${REPO}/scripts/setup.sh" --strchive-data "$SETUP" > "${CASE_TMP}/setup.log" 2>&1
check_eq "STR_EXTRA_LOCI=yes setup.sh --strchive-data exits 0" "$?" 0
cat "${CASE_TMP}/setup.log"
check_eq "the catalog has STRchive's 82 entries and the 2 extra loci" "$(count "$CAT")" 84
check "catalog.stamp names the extra catalog's sha256" has "extra=$(sha256sum "$EXTRA" | cut -d' ' -f1)" \
  "$(cat "${CAT%/*}/catalog.stamp" 2>/dev/null)"
jq 'map(select(.Source != "STRchive"))' "$CAT" > "${G}/reference-str-extra.json"
check_eq "the extra loci of the merged catalog" "$(jq -r 'map(.LocusId) | join(",")' "${G}/reference-str-extra.json")" "CFTR,UGT1A1"
"${REPO}/scripts/setup.sh" --strchive-data "$SETUP" > "${CASE_TMP}/setup-off.log" 2>&1
check_eq "control: setup.sh again without the opt-in exits 0" "$?" 0
check_eq "control: without the opt-in the catalog is back to 82 entries" "$(count "$CAT")" 82
check "control: catalog.stamp says extra=off" has 'extra=off$' "$(cat "${CAT%/*}/catalog.stamp" 2>/dev/null)"

# --- 2. made-up reference and reads ---------------------------------------------------------
RD="${G}/reference-str-extra"
rm -rf "$RD"
mkdir -p "$RD"
for s in EXA EXB; do rm -rf "${G:?}/${s}"; mkdir -p "${G}/${s}/aligned"; done
"${REPO}/bin/og" dev fixture str-extra "${RD}/ref.fa" EXA "${G}/EXA/aligned/EXA.sam" EXB "${G}/EXB/aligned/EXB.sam"
check "the made-up reference was written" test -s "${RD}/ref.fa"
sam faidx reference-str-extra/ref.fa
check_eq "ref.fa.fai: chr7 and chr2" "$(cut -f1 "${RD}/ref.fa.fai" 2>/dev/null | tr '\n' ' ')" "chr7 chr2 "
check_eq "the CFTR tract on the made-up chr7 is GRCh38's TG11 T7" \
  "$(sam faidx reference-str-extra/ref.fa chr7:117548601-117548641 2>/dev/null | tail -n +2 | tr -d '\n')" \
  "TTTTGATGTGTGTGTGTGTGTGTGTGTGTTTTTTTAACAGG"
for s in EXA EXB; do
  in_genome "$SAMTOOLS_IMAGE" bash -c "samtools sort -o '${s}/aligned/${s}_sorted.bam' '${s}/aligned/${s}.sam' && samtools index '${s}/aligned/${s}_sorted.bam'"
  check "${s}: the BAM is indexed" test -s "${G}/${s}/aligned/${s}_sorted.bam.bai"
done

# --- 3. step 09 on the extra loci -----------------------------------------------------------
repcn() {  # repcn <sample> <VARID>: REPCN of that record
  # shellcheck disable=SC2016  # awk program
  awk -F'\t' -v id="$2" '!/^#/ {
      n = split($8, kv, ";"); v = ""
      for (i = 1; i <= n; i++) if (kv[i] ~ /^VARID=/) v = substr(kv[i], 7)
      if (v != id) next
      nk = split($9, k, ":"); split($10, x, ":"); r = "."
      for (i = 1; i <= nk; i++) if (k[i] == "REPCN") r = x[i]
      print r
    }' "${G}/$1/expansion_hunter/$1_eh.vcf" 2>/dev/null
}
export REF_FASTA="${RD}/ref.fa"
for s in EXA EXB; do
  EH_CATALOG="${G}/reference-str-extra.json" run_step 09-expansion-hunter.sh "$s" female
  check_step_exit 09-expansion-hunter.sh
  cp "$STEP_LOG" "${CASE_TMP}/${s}.log"
  check "${s}: step 09 prints CFTR as not graded, with the asset's reason" \
    has 'Not graded.*CFTR: the CFTR intron 9 TG and poly-T tract' "$(cat "$STEP_LOG")"
  check "${s}: step 09 prints UGT1A1 as not graded, with the asset's reason" \
    has 'Not graded.*UGT1A1: the UGT1A1 promoter TA repeat' "$(cat "$STEP_LOG")"
  check_eq "${s}: one ExpansionHunter record per repeat (CFTR_TG, CFTR_T, UGT1A1)" \
    "$(grep -vc '^#' "${G}/${s}/expansion_hunter/${s}_eh.vcf" 2>/dev/null || echo 0)" 3
done
unset REF_FASTA
echo "counts: EXA $(repcn EXA CFTR_TG) $(repcn EXA CFTR_T) $(repcn EXA UGT1A1); EXB $(repcn EXB CFTR_TG) $(repcn EXB CFTR_T) $(repcn EXB UGT1A1)"
check_eq "EXA: CFTR TG counts" "$(repcn EXA CFTR_TG)" "11/12"
check_eq "EXA: CFTR poly-T counts, one 5T copy and one 7T copy" "$(repcn EXA CFTR_T | tr '/' '\n' | sort -n | paste -sd/ -)" "5/7"
check_eq "EXA: UGT1A1 TA units 7/8 (TA6/TA7, *1/*28)" "$(repcn EXA UGT1A1)" "7/8"
check_eq "EXB: CFTR TG counts" "$(repcn EXB CFTR_TG)" "11/11"
check_eq "EXB (red control): CFTR poly-T 7/7, no 5T copy" "$(repcn EXB CFTR_T)" "7/7"
check_eq "EXB: UGT1A1 TA units 7/7 (*1/*1)" "$(repcn EXB UGT1A1)" "7/7"

in_genome "$SAMTOOLS_IMAGE" rm -rf EXA EXB reference-str-extra reference-str-extra.json
finish
