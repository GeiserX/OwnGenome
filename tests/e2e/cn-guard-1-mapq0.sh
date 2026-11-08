#!/usr/bin/env bash
# The MAPQ 0 guard of the CYP21A2 row of bin/cn_loci.tsv on real reads: HG002
# from the 35x BWA-MEM GRCh38 no-ALT BAM of DeepVariant's case study, over the
# row's segments (CYP21A2, CYP21A1P, the two flanks). Two BAMs:
#   - unmodified: the MAPQ 0 share over the two targets is over 0.05 (HG002's
#     reads give about 0.1), and the row is never not assessed for it;
#   - every MAPQ 0 read removed (samtools view -q 1), as an aligner that leaves
#     out the reads fitting both copies writes it: not assessed,
#     low_mapq0_share, no copies, no module total, no finding, and a message
#     that says the BAM has no multi-mapping reads.
# Red control: the BAM without MAPQ 0 reads gives a lower estimate (the
# table keeps it when not assessed), so the guard stops a real drift.
# The block of step 35 runs on each (cut from scripts/35-paralogs.sh, with
# scripts/lib/common.sh's run_in), and the PARASCOPY module's block, cut from
# its main.nf, writes the same table on the BAM without MAPQ 0 reads, byte for
# byte.
. "$(dirname "$0")/lib.sh"

SRC="https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"
P=HG002cnguard
W="${GENOME_DIR}/${P}"
rm -rf "$W"; mkdir -p "${W}/aligned" "${W}/paralogs"
python3 "${REPO}/bin/copy_number.py" bed --locus CYP21A2 | sort -k1,1 -k2,2n > "${W}/fetch.bed"
check_ge "segments of the CYP21A2 row to fetch" "$(wc -l < "${W}/fetch.bed" | tr -d ' ')" 10

retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && { echo "  try ${i} of 4 failed; again in $((i * 15)) s"; sleep $((i * 15)); }; done; return 1; }
fetch_reads() {
  rm -f "${W}/unmodified.bam" "${W}/$(basename "$SRC").bai"
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${GENOME_DIR}:/genome" -w "/genome/${P}" "$SAMTOOLS_IMAGE" \
    samtools view -b -M -L "/genome/${P}/fetch.bed" -o "/genome/${P}/unmodified.bam" "$SRC"
}
check "HG002's reads over the CYP21A2 row download" retry fetch_reads
rm -f "${W}/$(basename "$SRC").bai"
sam index "${P}/unmodified.bam"
sam view -b -q 1 -o "${P}/no_mapq0.bam" "${P}/unmodified.bam" && sam index "${P}/no_mapq0.bam"
all=$(sam view -c "${P}/unmodified.bam" 2>/dev/null || echo 0)
kept=$(sam view -c "${P}/no_mapq0.bam" 2>/dev/null || echo 0)
check_ge "HG002 reads over the CYP21A2 row" "$all" 5000
check_eq "no read of MAPQ 0 is left in the BAM without them" \
  "$(sam view -c -e 'mapq == 0' "${P}/no_mapq0.bam" 2>/dev/null | tr -d ' ')" 0
check "the BAM without MAPQ 0 reads lost some (${kept} of ${all} kept)" awk -v a="$all" -v k="$kept" 'BEGIN {exit !(k > 0 && k < a)}'

col() { awk -F'\t' -v k="$2" '$1 == k {print $2}' "$1" 2>/dev/null; }

# --- step 35's CYP21A2 block, with common.sh's run_in ---------------------------------
sed -n '/^CYP21_OUT=/,/^rm -rf "\$CYP21_WORK"$/p' "${REPO}/scripts/35-paralogs.sh" > "${CASE_TMP}/step_body.sh"
check_ge "step 35's CYP21A2 block was cut from the script" "$(wc -l < "${CASE_TMP}/step_body.sh" | tr -d ' ')" 20
check "step 35's CYP21A2 block runs the MAPQ 1 pass" grep -q -- '--bedcov-mapq1' "${CASE_TMP}/step_body.sh"
step_block() {  # step_block NAME: the block on ${P}/<NAME>.bam, its table copied to CASE_TMP/<NAME>.tsv
  { echo 'set -euo pipefail'
    printf 'SAMPLE=%q GENOME_DIR=%q\n' "$P" "$GENOME_DIR"
    printf '. %q\n' "${REPO}/scripts/lib/common.sh"
    printf 'OUTDIR=%q BAM=%q\n' "${W}/paralogs" "${W}/$1.bam"
    cat "${CASE_TMP}/step_body.sh"; } > "${CASE_TMP}/block.sh"
  echo "+ step 35's CYP21A2 block ($1)"
  bash "${CASE_TMP}/block.sh" 2>&1 | tail -n 5
  cp "${W}/paralogs/${P}_cyp21a2_psv.tsv" "${CASE_TMP}/$1.tsv" 2>/dev/null || : > "${CASE_TMP}/$1.tsv"
  cat "${CASE_TMP}/$1.tsv"
}
for v in unmodified no_mapq0; do step_block "$v"; done
U="${CASE_TMP}/unmodified.tsv" N="${CASE_TMP}/no_mapq0.tsv"

check "unmodified: the MAPQ 0 share ($(col "$U" mapq0_share)) is over 0.05" \
  awk -v s="$(col "$U" mapq0_share)" 'BEGIN {exit !(s != "" && s != "." && s + 0 > 0.05)}'
check "unmodified: never not assessed for the MAPQ 0 share" lacks 'low_mapq0_share' "$(col "$U" reason)"
check_eq "without MAPQ 0 reads: not assessed (low_mapq0_share), share 0.000, no copies, module total or finding" \
  "$(col "$N" status)/$(col "$N" reason)/$(col "$N" mapq0_share)/$(col "$N" cyp21a2_copies)/$(col "$N" module_total)/$(col "$N" finding)" \
  "not assessed/low_mapq0_share/0.000/././."
check "without MAPQ 0 reads: the message says the BAM has no multi-mapping reads" \
  has '^not assessed: the BAM has no multi-mapping reads here' "$(col "$N" message)"

# Red control: the estimate the guard withheld. Removing the MAPQ 0 reads lowers
# the module total estimate (4.22 to 3.81 on these reads, measured 2026-11-08),
# the drift that, on a genome whose depth runs lower, reads as one module fewer.
check "red control: without MAPQ 0 reads the estimate ($(col "$N" estimate)) is at least 0.2 under the unmodified one ($(col "$U" estimate))" \
  awk -v n="$(col "$N" estimate)" -v u="$(col "$U" estimate)" 'BEGIN {exit !(n != "" && n != "." && u != "" && u != "." && n + 0.2 <= u + 0)}'
echo "cn-guard: HG002 MAPQ 0 share $(col "$U" mapq0_share) ($(col "$U" status), module total $(col "$U" module_total)," \
  "estimate $(col "$U" estimate)); without MAPQ 0 reads ${kept} of ${all} reads, $(col "$N" status) ($(col "$N" reason));" \
  "estimate $(col "$N" estimate)" >> "$E2E_NOTES"

# --- the PARASCOPY module's CYP21A2 block, in its image -------------------------------
MOD="${REPO}/modules/local/parascopy/main.nf"
sed -n '/# CYP21A2 from depth and PSVs (as scripts\/35-paralogs.sh)/,/rm -f cyp21a2_regions.bed/p' "$MOD" \
  | sed -e 's/\\\\/\\/g' -e 's/\\\$/$/g' > "${CASE_TMP}/module_body.sh"
check_ge "the module's CYP21A2 block was cut from main.nf" "$(wc -l < "${CASE_TMP}/module_body.sh" | tr -d ' ')" 15
M="${W}/module"
rm -rf "$M"; mkdir -p "$M"
sed -e "s|\${bam}|/genome/${P}/no_mapq0.bam|g" -e "s|\${meta.id}|${P}|g" "${CASE_TMP}/module_body.sh" > "${M}/block.sh"
docker run --rm --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -v "${REPO}/bin:/pgp-bin:ro" -w "/genome/${P}/module" \
  "$PARASCOPY_IMAGE" sh -c 'PATH="/pgp-bin:$PATH"; sh -e block.sh' 2>&1 | tail -n 3
check "the module's block writes step 35's table on the BAM without MAPQ 0 reads, byte for byte" \
  cmp -s "$N" "${M}/${P}_cyp21a2_psv.tsv"

rm -rf "$W"
finish
