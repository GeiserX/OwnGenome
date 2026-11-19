#!/usr/bin/env bash
# build-refget-table.sh: write assets/ownhealth/grch38_refget.tsv, the GA4GH refget identifier
# (sha512t24u of the upper-case sequence, "SQ.") and MD5 of each GRCh38 chromosome. og facts
# names the sequence of every VRS identifier it writes by these (rust/crates/og-facts/src/vrs.rs).
# The no-ALT analysis set the pipeline runs masks parts of chr5, chr19, chr22 and chrY, so the
# digest of the reference on disk is not GRCh38's for those four; the RefSeq sequences are.
#
# Usage: scripts/ci/build-refget-table.sh [OUT]   (default assets/ownhealth/grch38_refget.tsv)
# CHROMS="21 MT" limits it to those chromosomes, to check a few rows against the table; it then
# needs an OUT other than the checked-in table, which it would otherwise cut to those rows.
# Downloads NCBI's GRCh38 (GCF_000001405.26) chromosome FASTA, about 950 MB, one at a time.
# Needs curl, gzip, openssl, base64 and md5sum or md5.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd -P)"
TABLE=${REPO}/assets/ownhealth/grch38_refget.tsv
OUT=${1:-$TABLE}
# The physical path (pwd -P), so a symlinked folder cannot hide the table.
if [ -n "${CHROMS:-}" ] && { [ $# -eq 0 ] || [ "$(cd "$(dirname "$OUT")" 2>/dev/null && pwd -P)/$(basename "$OUT")" = "$TABLE" ]; }; then
  echo "ERROR: CHROMS writes only some chromosomes: give an OUT other than ${TABLE}" >&2
  exit 2
fi
[ ! -d "$OUT" ] || { echo "ERROR: ${OUT} is a folder; give the table's file name" >&2; exit 2; }
BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38
STRUCT=${BASE}/GCF_000001405.26_GRCh38_assembly_structure
WORK=$(mktemp -d)
PART=""
trap 'rm -rf "$WORK"; [ -z "$PART" ] || rm -f "$PART"' EXIT
# This run's own staging file beside OUT, so two runs never write one file.
PART=$(mktemp "${OUT}.XXXXXX")

md5_of() { if command -v md5sum >/dev/null; then md5sum | cut -d' ' -f1; else md5 -q; fi; }

curl -sfL --retry 3 "${BASE}/GCF_000001405.26_GRCh38_assembly_report.txt" -o "${WORK}/report.txt"
{
  echo "# assets/ownhealth/grch38_refget.tsv: the GA4GH refget identifier and MD5 of each GRCh38 chromosome."
  echo "# og facts writes VRS identifiers on these sequences. Written by scripts/ci/build-refget-table.sh from"
  echo "# NCBI's RefSeq GRCh38 chromosome FASTA (GCF_000001405.26), upper-cased; refget is sha512t24u (the"
  echo "# first 24 bytes of SHA-512, base64url), md5 is the M5 tag a .dict gives. Lengths from the assembly report."
  printf 'chrom\taccession\tlength\trefget\tmd5\n'
  for c in ${CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT}; do
    if [ "$c" = MT ]; then url=${STRUCT}/non-nuclear/assembled_chromosomes/FASTA/chrMT.fna.gz; name=chrM
    else url=${STRUCT}/Primary_Assembly/assembled_chromosomes/FASTA/chr${c}.fna.gz; name=chr${c}; fi
    curl -sfL --retry 3 "$url" -o "${WORK}/chr.fna.gz"
    gzip -dc "${WORK}/chr.fna.gz" | grep -v '^>' | tr -d '\n' | LC_ALL=C tr '[:lower:]' '[:upper:]' > "${WORK}/seq"
    acc=$(gzip -dc "${WORK}/chr.fna.gz" | sed -n '1s/^>\([^ ]*\).*/\1/p')
    len=$(awk -F'\t' -v a="$acc" '$7 == a {print $9}' "${WORK}/report.txt")
    [ "$len" = "$(wc -c < "${WORK}/seq" | tr -d ' ')" ] || { echo "ERROR: ${name}: length differs from the assembly report" >&2; exit 1; }
    openssl dgst -sha512 -binary < "${WORK}/seq" > "${WORK}/sha512"
    refget=$(head -c 24 "${WORK}/sha512" | base64 | tr '+/' '-_')
    md5=$(md5_of < "${WORK}/seq")
    [ "${#md5}" -eq 32 ] || { echo "ERROR: ${name}: no MD5" >&2; exit 1; }
    printf '%s\t%s\t%s\tSQ.%s\t%s\n' "$name" "$acc" "$len" "$refget" "$md5"
  done
} > "$PART"
chmod 644 "$PART"
mv "$PART" "$OUT"
echo "wrote ${OUT}"
