#!/usr/bin/env bash
# 41-bloodtype.sh — Blood group antigens (red cell, platelet, neutrophil) predicted from the genome
# Usage: ./scripts/41-bloodtype.sh <sample_name>
#
# Opt-in. Reads the sample's genotype at the defining variants of ABO (A1, A2,
# B, O), RhCE, Kell, Duffy (with the GATA box), Kidd, MNS, Lutheran, Diego,
# Dombrock, Colton, the platelet antigens HPA-1 to -6 and -15 and the
# neutrophil antigens HNA-1, -3, -4 and -5 (assets/bloodtype_markers.tsv),
# and RhD from read depth (the RHD row of bin/cn_loci.tsv): the common
# RhD-negative allele deletes the whole RHD gene, which a VCF does not show.
# The VCF's records around the marker sites are left-aligned against the
# reference first (bcftools norm, multi-allelic records kept whole), so an
# indel a caller wrote further right in a repeat still meets its marker. og
# bloodtype report then applies the rules, reading the marker sites through
# that file's .tbi; it is the same code as the BLOODTYPE Nextflow processes.
# A prediction from the DNA, not a serological typing.
#
# Input:  vcf/<sample>.vcf.gz; the reference (REF_FASTA) with its .fai;
#         aligned/<sample>_sorted.bam with its .bai
#         (ALIGN_DIR as for the other BAM steps), optional: without it RhD is
#         not determined, sites the VCF has no record for are read as
#         reference with confidence medium, and a system with no record at
#         any of its sites is not assessed.
# Output: bloodtype/<sample>_bloodtype.tsv           one row per system
#         bloodtype/<sample>_bloodtype_evidence.tsv  what the VCF and the reads said per site
#         bloodtype/<sample>_bloodtype_depth.tsv     samtools bedcov over the sites and RHD (with a BAM)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image OG_IMAGE SAMTOOLS_IMAGE BCFTOOLS_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
VCF="${S}/vcf/${SAMPLE}.vcf.gz"
BAM="${S}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam"
OUTDIR="${S}/bloodtype"
TABLE="${OUTDIR}/${SAMPLE}_bloodtype.tsv"
EVIDENCE="${OUTDIR}/${SAMPLE}_bloodtype_evidence.tsv"
DEPTH="${OUTDIR}/${SAMPLE}_bloodtype_depth.tsv"
NORM="${OUTDIR}/${SAMPLE}_bloodtype_norm.vcf.gz"

echo "=== Blood groups (predicted from the genome): ${SAMPLE} ==="
for f in "$VCF" "$REF_FASTA" "${REF_FASTA}.fai"; do
  [ -f "$f" ] || { echo "ERROR: File not found: ${f}" >&2; exit 1; }
done
echo "VCF: ${VCF}"
mkdir -p "$OUTDIR"
rm -f "$TABLE" "$EVIDENCE" "$DEPTH" "$NORM" "${NORM}.tbi" "${OUTDIR}/targets.bed" "${OUTDIR}/contigs.txt" "${OUTDIR}/rhd.bed"
O=$(cpath "$OUTDIR")
ARGS=()
if [ -f "$BAM" ] && [ -f "${BAM}.bai" ]; then
  echo "BAM: ${BAM}"
  echo "[1/3] Read depth at the marker sites and over RHD (samtools bedcov, MAPQ >= 20)..."
  # The marker table is a BED with extra columns; bedcov reads it once the
  # comment and header lines are gone. The RHD regions come from
  # og cn bed (BLOODTYPE_REGIONS in Nextflow) as a 4-column BED.
  # bedcov stops on a contig the BAM has no @SQ line for, so those rows are
  # left out: their sites get no depth row and are not determined. The
  # BLOODTYPE_DEPTH process runs the same commands.
  # shellcheck disable=SC2016  # $1 belongs to the inner shell
  run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 1g "$OG_IMAGE" \
    sh -c '/pgp-bin/og cn bed --locus RHD > "$1"' _ "${O}/rhd.bed"
  # shellcheck disable=SC2016  # $1 to $6 belong to the inner bash, the two programs to awk
  run_in -v "${PGP_ROOT}/assets:/pgp-assets:ro" --cpus 1 --memory 2g "$SAMTOOLS_IMAGE" bash -euo pipefail -c '
      samtools view -H "$3" | awk -F "\t" "$5" > "$2/contigs.txt"
      { grep -v "^#" "$1" | tail -n +2; cat "$2/rhd.bed"; } | awk -F "\t" "$6" "$2/contigs.txt" - > "$2/targets.bed"
      if [ -s "$2/targets.bed" ]; then samtools bedcov -Q 20 -j "$2/targets.bed" "$3" > "$4"; else : > "$4"; fi' \
    _ /pgp-assets/bloodtype_markers.tsv "$O" "$(cpath "$BAM")" "$(cpath "$DEPTH")" \
    '$1 == "@SQ" {for (i = 2; i <= NF; i++) if ($i ~ /^SN:/) print substr($i, 4)}' \
    'NR == FNR {c[$1] = 1; next} $1 in c'
  rm -f "${OUTDIR}/targets.bed" "${OUTDIR}/contigs.txt" "${OUTDIR}/rhd.bed"
  ARGS+=(--depth "$(cpath "$DEPTH")")
else
  echo "[1/3] No BAM (${BAM}): RhD needs read depth and is reported as not determined."
fi
echo "[2/3] The VCF's records around the marker sites, left-aligned (bcftools norm)..."
# og bloodtype report matches a marker by position, REF and ALT, so an indel
# must sit where the table writes it: left-aligned, as bcftools norm writes it.
# A caller can write it further right in a repeat (DeepVariant does without
# --normalize_reads, and a VCF from any caller is accepted). So the records
# within 100 bp of a marker site (more than the 38 bp the longest marker indel
# can shift along its repeat) are left-aligned against the reference first.
# A multi-allelic record stays whole: split, a 1/2 call becomes 1/0 and 0/1
# records and a null row would lose one of its two copies. -c w keeps a record
# whose REF differs from the reference as it is (og then reads it as a no
# call). Records on a contig the reference lacks (a VCF named 1, 2, ... against
# a chr-named reference) cannot be aligned and go through unchanged. The VCF is
# read through its index when the one htslib loads (.csi before .tbi) is
# current, else once from the start. Read that way, a VCF with no ##contig line
# for a window's contig first gets one, because bcftools writes no record on an
# undeclared contig: the new header, then the records as text, so a VCF in plain
# gzip (not bgzip) works too. Through an index, htslib takes the names from the
# index. bcftools writes no BCF record that uses a FILTER, INFO or FORMAT tag
# the header does not declare (a hand-made VCF may lack even the GT line), so
# the picked records go out as text and get a line for each such tag first.
# The BLOODTYPE process runs the same commands before og.
# shellcheck disable=SC2016  # the program's $ belong to the inner bash and awk
NORM_PROG='
  MARKERS=$1 VCF=$2 REF=$3 W=$4 OUT=$5 DECLARE=$9
  grep -v "^#" "$MARKERS" | tail -n +2 \
    | awk -F "\t" -v OFS="\t" "$6" | sort -k1,1 -k2,2n | awk -F "\t" -v OFS="\t" "$7" > "$W/bt_windows.bed"
  cut -f1 "$REF.fai" > "$W/bt_ref_contigs.txt"
  awk -F "\t" "NR == FNR {c[\$1] = 1; next} \$1 in c" "$W/bt_ref_contigs.txt" "$W/bt_windows.bed" > "$W/bt_on_ref.bed"
  awk -F "\t" "NR == FNR {c[\$1] = 1; next} !(\$1 in c)" "$W/bt_ref_contigs.txt" "$W/bt_windows.bed" > "$W/bt_off_ref.bed"
  SEL=-T OVERLAP=--targets-overlap HDR=$VCF
  for x in csi tbi; do
    if [ -f "$VCF.$x" ]; then
      if [ ! "$VCF" -nt "$VCF.$x" ]; then SEL=-R OVERLAP=--regions-overlap; fi
      break
    fi
  done
  if [ "$SEL" = -R ]; then echo "  The VCF is read through its index."
  else echo "  The VCF is read once from the start (no index, or one older than the VCF)."; fi
  if [ "$SEL" = -T ]; then
    bcftools head "$VCF" > "$W/bt_header.vcf"
    cut -f1 "$W/bt_windows.bed" | awk "$8" "$W/bt_header.vcf" - > "$W/bt_contigs.txt"
    if [ -s "$W/bt_contigs.txt" ]; then
      { grep -v "^#CHROM" "$W/bt_header.vcf"; cat "$W/bt_contigs.txt"; grep "^#CHROM" "$W/bt_header.vcf"; } > "$W/bt_header2.vcf"
      HDR=$W/bt_header2.vcf
    fi
  fi
  pick() {
    if [ ! -s "$1" ]; then bcftools view -h -Ov "$HDR"
    elif [ "$HDR" != "$VCF" ]; then { cat "$HDR"; gzip -cdf "$VCF" | awk "!/^#/"; } | bcftools view -T "$1" --targets-overlap record -Ov -
    else bcftools view "$SEL" "$1" "$OVERLAP" record -Ov "$VCF"; fi | awk -F "\t" "$DECLARE"
  }
  pick "$W/bt_on_ref.bed" | bcftools norm -c w -f "$REF" -Ob -o "$W/bt_on_ref.bcf"
  pick "$W/bt_off_ref.bed" | bcftools view -Ob -o "$W/bt_off_ref.bcf"
  bcftools concat -Ou "$W/bt_on_ref.bcf" "$W/bt_off_ref.bcf" | bcftools sort -T "$W/bt_sort" -Oz -o "$OUT"
  bcftools index -f -t "$OUT"
  rm -f "$W/bt_windows.bed" "$W/bt_ref_contigs.txt" "$W/bt_on_ref.bed" "$W/bt_off_ref.bed" "$W/bt_on_ref.bcf" "$W/bt_off_ref.bcf" \
    "$W/bt_header.vcf" "$W/bt_contigs.txt" "$W/bt_header2.vcf"'
# Each marker row 100 bp each side, under its chr name and without the chr
# (og reads a VCF named either way); then the overlapping windows merged.
# shellcheck disable=SC2016  # awk programs
WINDOWS_AWK='{s = $2 - 100; if (s < 0) s = 0; print $1, s, $3 + 100; c = $1; if (sub(/^chr/, "", c)) print c, s, $3 + 100}'
# shellcheck disable=SC2016
MERGE_AWK='$1 == c && $2 <= e {if ($3 > e) e = $3; next} {if (c != "") print c, s, e; c = $1; s = $2; e = $3} END {if (c != "") print c, s, e}'
# The window contigs (on stdin) the header (the first file) declares no
# ##contig line for, as ##contig lines.
# shellcheck disable=SC2016
CONTIG_AWK='NR == FNR {if (sub(/^##contig=<ID=/, "")) {sub(/[,>].*/, ""); d[$0] = 1}; next} !($1 in d) && !s[$1]++ {print "##contig=<ID=" $1 ">"}'
# The picked VCF text, with a header line for each FILTER, INFO and FORMAT tag
# its records use and its header lacks, before #CHROM. Text the header already
# covers comes out unchanged.
# shellcheck disable=SC2016
DECLARE_AWK='function need(k, t) {if (!(k in d)) {d[k] = 1; a[++na] = "##" k t ",Description=" q "Not declared in the input VCF" q ">"}}
  BEGIN {q = sprintf("%c", 34)}
  /^##/ {h[++nh] = $0; if (match($0, /^##(FILTER|INFO|FORMAT)=<ID=[^,>]*/)) d[substr($0, 3, RLENGTH - 2)] = 1; next}
  /^#/ {c = $0; next}
  {r[++nr] = $0
   n = split($7, f, ";"); for (i = 1; i <= n; i++) if (f[i] != "." && f[i] != "PASS") need("FILTER=<ID=" f[i], "")
   n = split($8, f, ";"); for (i = 1; i <= n; i++) if (f[i] != ".") {k = f[i]; v = sub(/=.*/, "", k); need("INFO=<ID=" k, v ? ",Number=.,Type=String" : ",Number=0,Type=Flag")}
   n = split($9, f, ":"); for (i = 1; i <= n; i++) need("FORMAT=<ID=" f[i], f[i] == "GT" ? ",Number=1,Type=String" : ",Number=.,Type=String")}
  END {for (i = 1; i <= nh; i++) print h[i]; for (i = 1; i <= na; i++) print a[i]; if (c != "") print c; for (i = 1; i <= nr; i++) print r[i]}'
run_in -v "${PGP_ROOT}/assets:/pgp-assets:ro" --cpus 1 --memory 2g "$BCFTOOLS_IMAGE" bash -euo pipefail -c "$NORM_PROG" \
  _ /pgp-assets/bloodtype_markers.tsv "$(cpath "$VCF")" "$REF_FASTA_C" "$O" "$(cpath "$NORM")" "$WINDOWS_AWK" "$MERGE_AWK" \
  "$CONTIG_AWK" "$DECLARE_AWK"
[ -s "$NORM" ] || { echo "ERROR: no left-aligned VCF written (${NORM})" >&2; exit 1; }
echo "[3/3] Blood group rules (og bloodtype report)..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" -v "${PGP_ROOT}/assets:/pgp-assets:ro" --cpus 1 --memory 2g "$OG_IMAGE" \
  /pgp-bin/og bloodtype report --sample "$SAMPLE" --vcf "$(cpath "$NORM")" \
    --markers /pgp-assets/bloodtype_markers.tsv ${ARGS[@]+"${ARGS[@]}"} \
    --out "$(cpath "$TABLE")" --evidence "$(cpath "$EVIDENCE")"
rm -f "$NORM" "${NORM}.tbi"

[ -s "$TABLE" ] || { echo "ERROR: no table written (${TABLE})" >&2; exit 1; }
echo ""
echo "Table: ${TABLE}"
echo "Evidence per site: ${EVIDENCE}"
