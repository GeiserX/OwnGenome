#!/usr/bin/env bash
# Step 35's SMN1 g.27134T>G reading (bin/traits_report.py --marker
# SMN1_27134) on real reads:
#   - dbSNP orientation, online: the record pinned for the marker is what
#     dbSNP's GRCh38 VCF of the build versions.env pins holds for rs143838139;
#     control: the same compare with REF and ALT swapped fails;
#   - HG002's Illumina reads over the marker's SMN1 and SMN2 copies (the 35x
#     BWA-MEM GRCh38 no-ALT BAM of DeepVariant's case study), three BAMs:
#       unmodified  HG002's own reading, whatever it is, goes to the notes
#       planted     every 4th read with T at a copy (a read of one M block)
#                   gets G there: 'present'
#       reference   every read with G at a copy gets T, and reads that start up
#                   to 500 bases before a copy with any other CIGAR are dropped:
#                   'not seen'
#     Running Parascopy on these windows would measure nothing, so the SMN
#     table beside them is a stand-in written here: two SMN1 copies at exon 7
#     (PASS, quality 30) and agCN 4. Step 35's marker block runs on each (cut
#     from scripts/35-paralogs.sh, with scripts/lib/common.sh's run_in), and
#     the PARASCOPY module's block, cut from its main.nf, writes the same table
#     in its image, byte for byte, on the planted BAM.
#   - every message says the marker does not prove 2+0.
. "$(dirname "$0")/lib.sh"

SRC="https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"
P=HG002smn
W="${GENOME_DIR}/${P}"
rm -rf "$W"; mkdir -p "${W}/aligned" "${W}/paralogs"
retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && { echo "  try ${i} of 4 failed; again in $((i * 15)) s"; sleep $((i * 15)); }; done; return 1; }
col() { awk -F'\t' -v k="$2" '$1 == k {print $2}' "$1" 2>/dev/null; }
SITES=$(python3 -c "import sys; sys.path.insert(0, '${REPO}/bin'); import traits_report as t; \
print(','.join(map(str, t.MARKERS['SMN1_27134']['copies'])))")
check_eq "two copies of the marker site (SMN1, SMN2)" "$SITES" "70952074,70076654"

# --- dbSNP orientation, online --------------------------------------------------------
URL="https://ftp.ncbi.nih.gov/snp/archive/${CONSUMER_DBSNP_BUILD:?}/VCF/GCF_000001405.40.gz"
PINNED=$(python3 -c "import sys; sys.path.insert(0, '${REPO}/bin'); import traits_report as t; \
d = t.MARKERS['SMN1_27134']['dbsnp']; print(d[0], d[2], d[3], d[4], sep='\t')")
same_record() {
  awk -F'\t' 'NR == 1 { id = $1; pos = $2; ref = $3; alt = $4; next }
              $1 == id && $2 == pos && $3 == ref && index("," $4 ",", "," alt ",") { ok = 1 }
              END { exit ok ? 0 : 1 }' <(echo "$1") <(echo "$2")
}
LIVE=$(retry docker run --rm "${BCFTOOLS_HTTPS_IMAGE:?}" bcftools query -r NC_000005.10:70952074-70952074 \
  -i 'ID="rs143838139"' -f '%ID\t%POS\t%REF\t%ALT\n' "$URL" 2>"${CASE_TMP}/dbsnp.err")
echo "pinned: ${PINNED}; dbSNP ${CONSUMER_DBSNP_BUILD}: ${LIVE:-<nothing>}"
check "the pinned SMN1 g.27134T>G record is dbSNP's" same_record "$PINNED" "$LIVE"
SWAPPED=$(awk -F'\t' -v OFS='\t' '{ t = $3; $3 = $4; $4 = t; print }' <<< "$PINNED")
check "control: the record with REF and ALT swapped does not match dbSNP" \
  test "$(same_record "$SWAPPED" "$LIVE"; echo $?)" -ne 0

# --- HG002's reads over both copies ---------------------------------------------------
fetch_reads() {
  rm -f "${W}/win.bam" "${W}/$(basename "$SRC").bai"
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${GENOME_DIR}:/genome" -w "/genome/${P}" "$SAMTOOLS_IMAGE" \
    samtools view -b -o "/genome/${P}/win.bam" "$SRC" chr5:70075600-70077700 chr5:70951000-70953100
}
check "HG002's reads over the SMN1 and SMN2 copies download" retry fetch_reads
rm -f "${W}/$(basename "$SRC").bai"
sam view -h "${P}/win.bam" > "${W}/orig.sam"
rm -f "${W}/win.bam"

plant() {
  awk -F'\t' -v OFS='\t' -v sites="$SITES" -v mode="$1" '
    BEGIN { n = split(sites, s, ",") }
    /^@/ { print; next }
    $3 == "chr5" {
      for (i = 1; i <= n; i++) {
        if (s[i] < $4 || s[i] > $4 + 500) continue
        if ($6 !~ /^[0-9]+M$/) { if (mode == "reference") { drop = 1; break }; continue }
        if (s[i] >= $4 + length($10)) continue
        o = s[i] - $4 + 1; b = substr($10, o, 1)
        if (mode == "planted" && b == "T" && (k++ % 4) == 0) { $10 = substr($10, 1, o - 1) "G" substr($10, o + 1); f++ }
        if (mode == "reference" && b == "G") { $10 = substr($10, 1, o - 1) "T" substr($10, o + 1); f++ }
      }
      if (drop) { drop = 0; d++; next }
    }
    { print }
    END { print (f + 0) " bases changed, " (d + 0) " reads dropped" > "/dev/stderr" }' "${W}/orig.sam"
}
plant planted > "${W}/planted.sam" 2> "${CASE_TMP}/planted.log"
plant reference > "${W}/reference.sam" 2> "${CASE_TMP}/reference.log"
cp "${W}/orig.sam" "${W}/unmodified.sam"
echo "planted: $(cat "${CASE_TMP}/planted.log"); reference: $(cat "${CASE_TMP}/reference.log")"
check_ge "planted: bases changed to G" "$(awk '{print $1}' "${CASE_TMP}/planted.log")" 5
for v in unmodified planted reference; do
  sam sort -o "${P}/${v}.bam" "${P}/${v}.sam" && sam index "${P}/${v}.bam"
done
rm -f "${W}"/*.sam
check_ge "HG002 reads over both copies" "$(sam view -c "${P}/unmodified.bam" 2>/dev/null || echo 0)" 200

# The stand-in SMN table: two SMN1 copies at exon 7, agCN 4.
printf 'chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions\n' \
  > "${W}/paralogs/${P}_smn_copy_number.tsv"
printf 'chr5\t70925000\t70953000\tSMN1\tPASS\t4\t50\tPASS\t2,2\t30,30\tchr5:70049000-70078000\n' \
  >> "${W}/paralogs/${P}_smn_copy_number.tsv"

# --- step 35's marker block, with common.sh's run_in ----------------------------------
sed -n '/^SMN1M_OUT=/,/^rm -rf "\$SMN1M_WORK"$/p' "${REPO}/scripts/35-paralogs.sh" > "${CASE_TMP}/step_body.sh"
check_ge "step 35's marker block was cut from the script" "$(wc -l < "${CASE_TMP}/step_body.sh" | tr -d ' ')" 15
step_block() {  # step_block NAME: the block on ${P}/<NAME>.bam, its table copied to CASE_TMP/<NAME>.tsv
  { echo 'set -euo pipefail'
    printf 'SAMPLE=%q GENOME_DIR=%q\n' "$P" "$GENOME_DIR"
    printf '. %q\n' "${REPO}/scripts/lib/common.sh"
    printf 'OUTDIR=%q BAM=%q\n' "${W}/paralogs" "${W}/$1.bam"
    cat "${CASE_TMP}/step_body.sh"; } > "${CASE_TMP}/block.sh"
  echo "+ step 35's marker block ($1)"
  bash "${CASE_TMP}/block.sh" 2>&1 | tail -n 5
  cp "${W}/paralogs/${P}_smn1_marker.tsv" "${CASE_TMP}/$1.tsv" 2>/dev/null || : > "${CASE_TMP}/$1.tsv"
  cat "${CASE_TMP}/$1.tsv"
}
for v in unmodified planted reference; do step_block "$v"; done
U="${CASE_TMP}/unmodified.tsv" L="${CASE_TMP}/planted.tsv" R="${CASE_TMP}/reference.tsv"

check "unmodified: a reading" has '^(present|not seen|uncertain)$' "$(col "$U" status)"
check_ge "unmodified: 10 reads or more at the two copies" "$(col "$U" depth)" 10
check_eq "the SMN1 psCN comes from the SMN table" "$(col "$U" smn1_pscn)" "2 (quality 30, filter PASS)"
check_eq "planted: present" "$(col "$L" status)" present
check "planted: the fraction ($(col "$L" fraction)) is about a quarter of the T reads (0.15 to 0.45)" \
  awk -v f="$(col "$L" fraction)" 'BEGIN {exit !(f >= 0.15 && f <= 0.45)}'
check "planted: copies = fraction x agCN" \
  awk -v f="$(col "$L" fraction)" -v c="$(col "$L" copies)" 'BEGIN {d = f * 4 - c; exit !(c > 0 && d < 0.06 && d > -0.06)}'
check "planted: with two SMN1 copies it raises the chance of 2+0" has "higher chance" "$(col "$L" message)"
check_eq "reference: not seen" "$(col "$R" status)" "not seen"
check_eq "reference: no read with G" "$(col "$R" alt_reads)" 0
for v in U L R; do
  check "${v}: the message says the marker does not prove 2+0" has "does not prove 2\+0" "$(col "${!v}" message)"
done
check "control: the planted and the reference tables differ" bash -c '! cmp -s "$1" "$2"' _ "$L" "$R"
echo "genotype-rules-2-smn1-marker: HG002 SMN1 g.27134T>G $(col "$U" status), $(col "$U" alt_reads) of" \
  "$(col "$U" depth) reads (fraction $(col "$U" fraction)); planted $(col "$L" fraction)" >> "$E2E_NOTES"

# --- the PARASCOPY module's marker block, in its image --------------------------------
MOD="${REPO}/modules/local/parascopy/main.nf"
sed -n '/# SMN1 g.27134T>G beside the SMN1 psCN (as scripts\/35-paralogs.sh)/,/rm -f smn1_marker_sites.txt smn1_marker_pileup.tsv/p' "$MOD" \
  | sed -e 's/\\\\/\\/g' -e 's/\\\$/$/g' > "${CASE_TMP}/module_body.sh"
check_ge "the module's marker block was cut from main.nf" "$(wc -l < "${CASE_TMP}/module_body.sh" | tr -d ' ')" 12
M="${W}/module"
rm -rf "$M"; mkdir -p "$M"
cp "${W}/paralogs/${P}_smn_copy_number.tsv" "$M/"
sed -e "s|\${bam}|/genome/${P}/planted.bam|g" -e "s|\${meta.id}|${P}|g" "${CASE_TMP}/module_body.sh" > "${M}/block.sh"
docker run --rm --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -v "${REPO}/bin:/pgp-bin:ro" -w "/genome/${P}/module" \
  "$PARASCOPY_IMAGE" sh -c 'PATH="/pgp-bin:$PATH"; sh -e block.sh' 2>&1 | tail -n 3
check "the module's block writes step 35's table, byte for byte" cmp -s "$L" "${M}/${P}_smn1_marker.tsv"

rm -rf "$W"
finish
