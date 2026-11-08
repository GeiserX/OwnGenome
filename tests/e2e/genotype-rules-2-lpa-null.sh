#!/usr/bin/env bash
# Step 46's KIV-2 4925G>A reading (bin/lpa_kiv2.py null-allele, MARKERS in
# bin/traits_report.py) on real reads:
#   - dbSNP orientation, online: the record pinned for the marker is what
#     dbSNP's GRCh38 VCF of the build versions.env pins holds for rs1484919544
#     (same position and REF, the pinned ALT among dbSNP's); control: the same
#     compare with REF and ALT swapped fails;
#   - HG002's Illumina reads over LPA (chr6:160.53-160.67 Mb, the 35x BWA-MEM
#     GRCh38 no-ALT BAM of DeepVariant's case study, as case lpa-1 uses),
#     carried onto the fixture's header. Three samples, each run through
#     scripts/46-lpa.sh:
#       unmodified  HG002's own reading, whatever it is, goes to the notes; its
#                   reads at the seven copies are counted
#       planted     every 8th read with C at one of the seven copies (a read
#                   of one M block) gets T there: 'present', with the copies
#                   that carry it about fraction x (N + 2)
#       reference   every read with T at a copy gets C, and reads that start
#                   up to 500 bases before a copy with any other CIGAR are
#                   dropped: 'not seen'
#     The reading goes to lpa/<sample>_lpa_kiv2_4925.tsv; step 46's own table
#     keeps its rows.
#   - the LPA module's pileup and 4925G>A lines, cut from its main.nf, write the
#     same pileup and table in their images, byte for byte, on the planted BAM.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR" T="$CASE_TMP"
SRC="https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"
BASE="${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
dk() { docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$1" "${@:2}"; }
retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && { echo "  try ${i} of 4 failed; again in $((i * 15)) s"; sleep $((i * 15)); }; done; return 1; }
metric() { awk -F'\t' -v k="$2" '$1 == k {print $2}' "${G}/$1/lpa/$1_lpa_kiv2.tsv" 2>/dev/null; }
null() { awk -F'\t' -v k="$2" '$1 == k {print $2}' "${G}/$1/lpa/$1_lpa_kiv2_4925.tsv" 2>/dev/null; }
sample() {
  rm -rf "${G:?}/$1"; mkdir -p "${G}/$1/aligned"
  cp "$2" "${G}/$1/aligned/$1_sorted.bam"; cp "$2.bai" "${G}/$1/aligned/$1_sorted.bam.bai"
}
SITES=$(python3 -c "import sys; sys.path.insert(0, '${REPO}/bin'); import traits_report as t; \
print(','.join(map(str, t.MARKERS['LPA_KIV2_4925']['copies'])))")
check_eq "seven copies of the 4925G>A site" "$(tr ',' '\n' <<< "$SITES" | wc -l | tr -d ' ')" 7

# --- dbSNP orientation, online --------------------------------------------------------
URL="https://ftp.ncbi.nih.gov/snp/archive/${CONSUMER_DBSNP_BUILD:?}/VCF/GCF_000001405.40.gz"
PINNED=$(python3 -c "import sys; sys.path.insert(0, '${REPO}/bin'); import traits_report as t; \
d = t.MARKERS['LPA_KIV2_4925']['dbsnp']; print(d[0], d[2], d[3], d[4], sep='\t')")
same_record() {
  awk -F'\t' 'NR == 1 { id = $1; pos = $2; ref = $3; alt = $4; next }
              $1 == id && $2 == pos && $3 == ref && index("," $4 ",", "," alt ",") { ok = 1 }
              END { exit ok ? 0 : 1 }' <(echo "$1") <(echo "$2")
}
LIVE=$(retry docker run --rm "${BCFTOOLS_HTTPS_IMAGE:?}" bcftools query -r NC_000006.12:160640667-160640667 \
  -i 'ID="rs1484919544"' -f '%ID\t%POS\t%REF\t%ALT\n' "$URL" 2>"${T}/dbsnp.err")
echo "pinned: ${PINNED}; dbSNP ${CONSUMER_DBSNP_BUILD}: ${LIVE:-<nothing>}"
check "the pinned 4925G>A record is dbSNP's" same_record "$PINNED" "$LIVE"
SWAPPED=$(awk -F'\t' -v OFS='\t' '{ t = $3; $3 = $4; $4 = t; print }' <<< "$PINNED")
check "control: the record with REF and ALT swapped does not match dbSNP" \
  test "$(same_record "$SWAPPED" "$LIVE"; echo $?)" -ne 0

# --- HG002's reads over LPA, on the fixture's header ---------------------------------
fetch_reads() {
  rm -f "${T}/win.bam" "${T}/$(basename "$SRC").bai"
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$SAMTOOLS_IMAGE" \
    samtools view -b -o "${T}/win.bam" "$SRC" chr6:160530000-160670000
}
check "case 20's BAM is there" test -s "${BASE}.bai"
check "HG002's reads over LPA download" retry fetch_reads
{ dk "$SAMTOOLS_IMAGE" samtools view -H "$BASE" | grep -v -e '^@RG' -e '^@PG'
  dk "$SAMTOOLS_IMAGE" samtools view -H "${T}/win.bam" | grep '^@RG'
  dk "$SAMTOOLS_IMAGE" samtools view "${T}/win.bam" | awk -F'\t' '$7 == "=" || $7 == "*" || $7 == "chr6"'
} > "${T}/orig.sam"
rm -f "${T}/win.bam"

# plant MODE: the SAM with the reads at the seven copies changed (see the header).
plant() {
  awk -F'\t' -v OFS='\t' -v sites="$SITES" -v mode="$1" '
    BEGIN { n = split(sites, s, ",") }
    /^@/ { print; next }
    $3 == "chr6" {
      for (i = 1; i <= n; i++) {
        if (s[i] < $4 || s[i] > $4 + 500) continue
        if ($6 !~ /^[0-9]+M$/) { if (mode == "reference") { drop = 1; break }; continue }
        if (s[i] >= $4 + length($10)) continue
        o = s[i] - $4 + 1; b = substr($10, o, 1)
        if (mode == "planted" && b == "C" && (k++ % 8) == 0) { $10 = substr($10, 1, o - 1) "T" substr($10, o + 1); f++ }
        if (mode == "reference" && b == "T") { $10 = substr($10, 1, o - 1) "C" substr($10, o + 1); f++ }
      }
      if (drop) { drop = 0; d++; next }
    }
    { print }
    END { print (f + 0) " bases changed, " (d + 0) " reads dropped" > "/dev/stderr" }' "${T}/orig.sam"
}
plant planted > "${T}/planted.sam" 2> "${T}/planted.log"
plant reference > "${T}/reference.sam" 2> "${T}/reference.log"
echo "planted: $(cat "${T}/planted.log"); reference: $(cat "${T}/reference.log")"
check_ge "planted: bases changed to T" "$(awk '{print $1}' "${T}/planted.log")" 10

H="${SAMPLE}nul"
for v in unmodified planted reference; do
  src="${T}/${v}.sam"; [ "$v" = unmodified ] && src="${T}/orig.sam"
  dk "$SAMTOOLS_IMAGE" samtools sort -o "${T}/${v}.bam" "$src"
  dk "$SAMTOOLS_IMAGE" samtools index "${T}/${v}.bam"
  sample "${H}${v:0:1}" "${T}/${v}.bam"
  run_step 46-lpa.sh "${H}${v:0:1}"; check_step_exit 46-lpa.sh
  cat "${G}/${H}${v:0:1}/lpa/${H}${v:0:1}_lpa_kiv2.tsv" "${G}/${H}${v:0:1}/lpa/${H}${v:0:1}_lpa_kiv2_4925.tsv" 2>/dev/null
done
U="${H}u" P="${H}p" R="${H}r"
check_ge "HG002 reads over LPA" "$(dk "$SAMTOOLS_IMAGE" samtools view -c "${T}/unmodified.bam" 2>/dev/null || echo 0)" 20000

check_eq "unmodified: a count" "$(metric "$U" status)" estimated
check "unmodified: a 4925G>A reading of the seven copies" has '^(present|not seen|uncertain)$' "$(null "$U" status)"
check_ge "unmodified: 100 reads or more at the seven copies" "$(null "$U" depth)" 100
check "unmodified: the pileup is published, one line per copy with reads" \
  test "$(wc -l < "${G}/${U}/lpa/${U}_lpa_kiv2_4925_pileup.tsv" | tr -d ' ')" -eq 7
check "unmodified: the count's message is the count's alone" bash -c '! grep -q 4925G "$1"' _ "${G}/${U}/lpa/${U}_lpa_kiv2.tsv"
check_eq "unmodified: the 4925G>A table carries the count" "$(null "$U" kiv2_copies)" "$(metric "$U" kiv2_copies)"

check_eq "planted: present" "$(null "$P" status)" present
check "planted: the fraction ($(null "$P" fraction)) is about an 8th of the C reads (0.08 to 0.25)" \
  awk -v f="$(null "$P" fraction)" 'BEGIN {exit !(f >= 0.08 && f <= 0.25)}'
check "planted: copies that carry it = fraction x (N + 2)" \
  awk -v f="$(null "$P" fraction)" -v n="$(null "$P" kiv2_copies)" -v c="$(null "$P" copies)" \
    'BEGIN {d = f * (n + 2) - c; exit !(c > 0 && d < 0.06 && d > -0.06)}'
check "planted: the sentence names what a short-read fraction cannot say" \
  has "cannot say which chromosome or which repeat" "$(null "$P" message)"
check_eq "planted: the count is the unmodified one (only bases changed)" "$(metric "$P" kiv2_copies)" "$(metric "$U" kiv2_copies)"

check_eq "reference: not seen" "$(null "$R" status)" "not seen"
check_eq "reference: no read with T" "$(null "$R" alt_reads)" 0
check "control: the planted and the reference tables differ" bash -c '! cmp -s "$1" "$2"' _ \
  "${G}/${P}/lpa/${P}_lpa_kiv2_4925.tsv" "${G}/${R}/lpa/${R}_lpa_kiv2_4925.tsv"
echo "genotype-rules-2-lpa-null: HG002 KIV-2 4925G>A $(null "$U" status)," \
  "$(null "$U" alt_reads) of $(null "$U" depth) reads (fraction $(null "$U" fraction))," \
  "KIV-2 count $(metric "$U" kiv2_copies); planted $(null "$P" fraction), copies $(null "$P" copies)" >> "$E2E_NOTES"

# --- the LPA module's pileup and 4925G>A lines, in their images ----------------------
MOD="${REPO}/modules/local/lpa/main.nf"
M="${T}/module"
rm -rf "$M"; mkdir -p "$M"
unesc() { sed -e 's/\\\\/\\/g' -e 's/\\\$/$/g' -e "s|\${bam}|${G}/${P}/aligned/${P}_sorted.bam|g" -e "s|\${meta.id}|${P}|g" \
  -e "s|\${pileup}|${P}_lpa_kiv2_4925_pileup.tsv|g"; }
grep -E '^ *traits_report\.py --marker-regions LPA_KIV2_4925 > lpa_sites\.txt$' "$MOD" | unesc > "${M}/sites.sh"
sed -n '/^ *( while read -r r; do samtools mpileup /,/|| echo /p' "$MOD" | unesc > "${M}/pileup.sh"
sed -n '/^ *traits_report\.py --marker LPA_KIV2_4925 --pileup/,/--marker-out/p' "$MOD" | unesc > "${M}/null.sh"
check_eq "the module's sites, pileup and 4925G>A lines were cut from main.nf" \
  "$(cat "${M}/sites.sh" "${M}/pileup.sh" "${M}/null.sh" | wc -l | tr -d ' ')" 5
cat "${M}/sites.sh" "${M}/pileup.sh" "${M}/null.sh"
cp "${G}/${P}/lpa/${P}_lpa_kiv2.tsv" "$M/"
mod() { docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${E2E_WORK}:${E2E_WORK}" -v "${REPO}/bin:/pgp-bin:ro" \
  -w "$M" "$1" sh -c 'PATH="/pgp-bin:$PATH"; set -e; . "./$1"' _ "$2"; }
mod "$PYTHON_IMAGE" sites.sh && mod "$SAMTOOLS_IMAGE" pileup.sh && mod "$PYTHON_IMAGE" null.sh 2>&1 | tail -n 3
check "the module writes step 46's pileup, byte for byte" \
  cmp -s "${G}/${P}/lpa/${P}_lpa_kiv2_4925_pileup.tsv" "${M}/${P}_lpa_kiv2_4925_pileup.tsv"
check "the module writes step 46's 4925G>A table, byte for byte" \
  cmp -s "${G}/${P}/lpa/${P}_lpa_kiv2_4925.tsv" "${M}/${P}_lpa_kiv2_4925.tsv"
check "the module leaves step 46's table byte for byte" cmp -s "${G}/${P}/lpa/${P}_lpa_kiv2.tsv" "${M}/${P}_lpa_kiv2.tsv"

rm -rf "${G:?}/${U}" "${G:?}/${P}" "${G:?}/${R}" "${M:?}" "${T:?}"/*.sam "${T:?}"/*.bam "${T:?}"/*.bai
finish
