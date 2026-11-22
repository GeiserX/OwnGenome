#!/usr/bin/env bash
# Step 06's cDNA read check (og cdna check) on the fixture. The fixture has
# no GENCODE file, so the case writes an exon file of its own: a made-up exon
# of 60 bp centred on every hit case 30 found, so each hit's reads cross exon
# edges on both sides, and a second exon 3 kb after the planted hit.
#   1. The real sample: the hits table gets the cdna_check column and no hit
#      is flagged; the planted hit, a SNV on genomic reads, is clear.
#   2. A copy of the sample whose reads at the planted hit are replaced by
#      genomic reads with the ref allele and reads with the alt allele clipped
#      at the exon end, the rest of each read on the second exon (SA tag), as
#      spliced cDNA aligns: the planted hit is flagged, and both reports print
#      the flag beside it.
. "$(dirname "$0")/lib.sh"

GENCODE_RELEASE=$(sed -n 's/^GENCODE_RELEASE=//p' "${REPO}/scripts/lib/common.sh")
EXONS="${GENOME_DIR}/reference/gencode.v${GENCODE_RELEASE}.basic.exons.bed"
trap 'rm -f "$EXONS"' EXIT
CHR=$(planted chrom)
POS=$(planted pos)
REF=$(planted ref)
ALT=$(planted alt)
P0=$((POS - 1))
E=$((P0 + 30))     # the end of the planted hit's exon (0-based, exclusive)
B0=$((P0 + 3000))  # the start of the second exon
HITS_TSV="${GENOME_DIR}/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.tsv"
check "case 30 left a hits table" test -s "$HITS_TSV"
{
  awk -F'\t' -v OFS='\t' 'NR > 1 {print $1, $2 - 31, $2 + 29}' "$HITS_TSV"
  printf '%s\t%s\t%s\n' "$CHR" "$B0" "$((B0 + 150))"
} | sort -k1,1 -k2,2n > "$EXONS"
echo "exons: $(wc -l < "$EXONS") lines; planted ${CHR}:${POS} ${REF}>${ALT}, exon end ${E}, next exon ${B0}"

# --- 1. the real sample ---------------------------------------------------------
run_step 06-clinvar-screen.sh "$SAMPLE"
check_step_exit 06-clinvar-screen.sh
T1="${GENOME_DIR}/${SAMPLE}/clinvar/${SAMPLE}_clinvar_cdna_check.tsv"
check_eq "the cDNA check table has the cdna_check column" \
  "$(head -n 1 "$T1" 2>/dev/null | tr '\t' '\n' | grep -cx cdna_check)" 1
check_eq "one row per hit" "$(awk 'NR > 1' "$T1" 2>/dev/null | wc -l | tr -d ' ')" \
  "$(awk 'NR > 1' "$HITS_TSV" | wc -l | tr -d ' ')"
check_eq "no real hit is flagged" "$(awk -F'\t' 'NR > 1 && $9 == "flagged"' "$T1" 2>/dev/null | wc -l | tr -d ' ')" 0
check_eq "the planted hit on genomic reads is checked and clear" \
  "$(awk -F'\t' -v c="$CHR" -v p="$POS" '$1 == c && $2 == p {print $9}' "$T1" 2>/dev/null)" clear
check_ge "the planted hit's alt reads are counted" \
  "$(awk -F'\t' -v c="$CHR" -v p="$POS" '$1 == c && $2 == p {print $5}' "$T1" 2>/dev/null)" 5

# --- 2. a copy with spliced alt reads at the planted hit -------------------------
T="${SAMPLE}cdna"
rm -rf "${GENOME_DIR:?}/${T}"
mkdir -p "${GENOME_DIR}/${T}/vcf" "${GENOME_DIR}/${T}/aligned"
ln "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${GENOME_DIR}/${T}/vcf/${T}.vcf.gz"
ln "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" "${GENOME_DIR}/${T}/vcf/${T}.vcf.gz.tbi"
FROM=$((P0 - 2000))
TO=$((P0 + 5000))
in_user() { docker run --rm -u "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -w /genome "$SAMTOOLS_IMAGE" samtools "$@"; }
in_user view -h "${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${CHR}:$((FROM + 1))-${TO}" > "${CASE_TMP}/region.sam"
in_user faidx reference/GRCh38_no_alt_analysis_set.fasta "${CHR}:$((FROM + 1))-${TO}" > "${CASE_TMP}/region.fa"
# The region's reads, less those over the planted position, plus 10
# genomic reads with the ref allele running on into the intron and 12
# spliced alt reads: clipped at the exon end, the rest on the next exon.
awk -F'\t' -v chrom="$CHR" -v frm="$FROM" -v p0="$P0" -v ref="$REF" -v alt="$ALT" -v e="$E" -v b0="$B0" \
    -v out="${GENOME_DIR}/${T}/aligned/planted.sam" '
  function base(a, b) { return substr(g, a - frm + 1, b - a) }
  function q(n,   s) { s = ""; while (n-- > 0) s = s "I"; return s }
  FILENAME == ARGV[1] { if (!/^>/) { sub(/^[[:space:]]+/, ""); sub(/[[:space:]]+$/, ""); g = g toupper($0) } next }
  /^@/ { lines[++n] = $0; next }
  {
    span = 0; c = $6
    while (match(c, /[0-9]+[MIDNSHP=X]/)) {
      op = substr(c, RSTART + RLENGTH - 1, 1)
      if (op ~ /[MDN=X]/) span += substr(c, RSTART, RLENGTH - 1)
      c = substr(c, RSTART + RLENGTH)
    }
    start = $4 - 1
    if (start <= p0 && p0 < start + span) { dropped++; next }
    lines[++n] = $0
  }
  END {
    if (base(p0, p0 + 1) != toupper(ref)) {
      print "the reference has " base(p0, p0 + 1) " at the planted position, not " ref > "/dev/stderr"
      exit 1
    }
    for (i = 0; i < 10; i++) {
      s = p0 - 50 + 3 * i
      lines[++n] = "ref" i "\t0\t" chrom "\t" s + 1 "\t60\t100M\t*\t0\t0\t" base(s, s + 100) "\t" q(100)
    }
    for (i = 0; i < 12; i++) {
      aligned = 70 + i; s = e - aligned; clip = 100 - aligned
      seq = base(s, p0) toupper(alt) base(p0 + 1, e) base(b0, b0 + clip)
      lines[++n] = "cdna" i "\t0\t" chrom "\t" s + 1 "\t60\t" aligned "M" clip "S\t*\t0\t0\t" seq "\t" q(100) \
        "\tSA:Z:" chrom "," b0 + 1 ",+," aligned "S" clip "M,60,0;"
    }
    for (i = 1; i <= n; i++) print lines[i] > out
    print "replaced " dropped + 0 " reads at the planted position with 10 ref and 12 spliced alt reads"
  }' "${CASE_TMP}/region.fa" "${CASE_TMP}/region.sam"
in_user sort -o "${T}/aligned/${T}_sorted.bam" "${T}/aligned/planted.sam"
in_user index "${T}/aligned/${T}_sorted.bam"
rm -f "${GENOME_DIR}/${T}/aligned/planted.sam"

run_step 06-clinvar-screen.sh "$T"
check_step_exit 06-clinvar-screen.sh
T2="${GENOME_DIR}/${T}/clinvar/${T}_clinvar_cdna_check.tsv"
ROW=$(awk -F'\t' -v c="$CHR" -v p="$POS" '$1 == c && $2 == p' "$T2" 2>/dev/null)
echo "planted row: ${ROW}"
check_eq "the spliced alt reads at the planted hit are flagged" "$(cut -f9 <<< "$ROW")" flagged
check_eq "12 of 12 alt reads and 0 of 10 ref reads clipped at the exon edge" "$(cut -f5-8 <<< "$ROW" | tr '\t' ' ')" "12 12 10 0"
check "the step's log names the flagged hit" has "FLAGGED ${CHR}:${POS} " "$(cat "$STEP_LOG")"
check_eq "the other hits stay unflagged" \
  "$(awk -F'\t' -v c="$CHR" -v p="$POS" 'NR > 1 && $9 == "flagged" && !($1 == c && $2 == p)' "$T2" 2>/dev/null | wc -l | tr -d ' ')" 0

# --- 3. both reports print the flag beside the hit --------------------------------
FLAG="possible cDNA/transgene contamination: check the reads"
"${REPO}/bin/og" report render --sample "$T" --sample-dir "${GENOME_DIR}/${T}" --json "${CASE_TMP}/t.json" \
  -o "${CASE_TMP}/t.txt" -o "${CASE_TMP}/t.html"
"${REPO}/bin/og" report render --sample "$SAMPLE" --sample-dir "${GENOME_DIR}/${SAMPLE}" \
  --json "${CASE_TMP}/s.json" -o "${CASE_TMP}/s.txt" -o "${CASE_TMP}/s.html"
check "the text report prints the flag" has "reads: ${FLAG}" "$(cat "${CASE_TMP}/t.txt" 2>/dev/null)"
check "the HTML hits table has the cDNA read check column" has '<th>cDNA read check</th>' "$(cat "${CASE_TMP}/t.html" 2>/dev/null)"
check "the HTML hits table prints the flag" has "<td>${FLAG}" "$(cat "${CASE_TMP}/t.html" 2>/dev/null)"
check "the real sample's HTML report has the column" has '<th>cDNA read check</th>' "$(cat "${CASE_TMP}/s.html" 2>/dev/null)"
check "and no flag" lacks "$FLAG" "$(cat "${CASE_TMP}/s.html" "${CASE_TMP}/s.txt" 2>/dev/null)"

finish
