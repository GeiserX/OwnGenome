#!/usr/bin/env bash
# build-sv-rule-regions.sh: rebuild assets/sv_rule_regions.tsv, the regions
# where og sv postfilter overrides AnnotSV's class 4 or 5.
#
# Usage: scripts/ci/build-sv-rule-regions.sh [OUT]   (default: assets/sv_rule_regions.tsv)
#        scripts/ci/build-sv-rule-regions.sh --check [FILE]   exit 1 when FILE
#            (default: assets/sv_rule_regions.tsv) is not what the downloads
#            give, byte for byte, and show the first lines that differ
#        scripts/ci/build-sv-rule-regions.sh --list   the downloads, one
#            "NAME URL SHA256" per line
#
# Downloads four pinned files (about 145 MB, each checked against its sha256)
# into SV_RULE_CACHE, where a file with the right sha256 is read as it is, or
# else into a temporary folder that is deleted, and writes OUT:
#
#   TRExplorer v2.1's GRCh38 tandem repeat catalog (broadinstitute/trexplorer-catalog)
#   GENCODE 50's basic annotation (the release scripts/lib/common.sh pins)
#   HGNC's complete set of 2026-10-06 (monthly archive)
#   the GRCh38 assembly regions report of NCBI (GCA_000001405.15)
#
# and keeps, all on GRCh38 chr1-22, X and Y:
#
#   tandem_repeat  every catalog repeat of 50 bp or more (the shortest SV)
#                  that overlaps a coding exon of a MANE Select transcript;
#                  name GENE:MOTIF (a motif over 12 bases as its length)
#   ig_tr_locus    the span of the IG and TR genes of GENCODE (gene_type IG_*
#                  or TR_*) of each of IGH, IGK, IGL, TRA (with TRD inside it),
#                  TRB and TRG: per locus the largest group of genes on one
#                  chromosome with gaps under 1 Mb, so the orphon copies on
#                  other chromosomes stay out
#   mhc            the MHC region of the assembly report
#   gene_alias     the GENCODE span of each HGNC pseudogene whose alias or
#                  previous symbol is the approved symbol of a gene with a
#                  GenCC disease curation; name ALIAS|PSEUDOGENE|its HGNC id|
#                  the HGNC id of the gene that holds ALIAS as its symbol
#
# A gene's span is that of the first GENCODE record of its HGNC id, and the
# _PAR_Y copies of the pseudoautosomal genes are skipped, as
# build-callability.sh does. The rules are in og dev build sv-rule-regions
# (rust/crates/og-dev/src/build/sv_rules.rs). Needs curl, og
# (scripts/setup.sh --og) and sha256sum or shasum. It reads about 3 GB of
# text in a few seconds and keeps little of it.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
CHECK=""
case "${1:-}" in
  --list) ;;
  --check) CHECK=${2:-${REPO}/assets/sv_rule_regions.tsv} ;;
  *) OUT=${1:-${REPO}/assets/sv_rule_regions.tsv} ;;
esac

TRX_URL=https://github.com/broadinstitute/trexplorer-catalog/releases/download/v2.1/TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.bed.gz
# Upstream rebuilt this v2.1 asset in place on 2026-10-30 (it had dropped
# loci); the sha256 is that of the rebuilt file, the only one still served.
TRX_SHA256=df9c23c4c66ae798759a00f1b042a255819dd95e6b5f2ac68b37e5072117642b
GENCODE_URL=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz
GENCODE_SHA256=5b9ad9aa71f6ab5d55123c51f74723644a895bc0cbaa9cee8fbc11bb9dec0bf7
HGNC_URL=https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/hgnc_complete_set_2026-10-06.txt
HGNC_SHA256=8bf6f686e640dbc864a346484626953621b0500d900d26cd085a35c2782d725d
REGIONS_URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_regions.txt
REGIONS_SHA256=6b91d2c7961a322936db8fd09b7fb4a9143a590882f28c0a628a4ef26897f6c7

if [ "${1:-}" = --list ]; then
  echo "trx.bed.gz ${TRX_URL} ${TRX_SHA256}"
  echo "gencode.gtf.gz ${GENCODE_URL} ${GENCODE_SHA256}"
  echo "hgnc.txt ${HGNC_URL} ${HGNC_SHA256}"
  echo "regions.txt ${REGIONS_URL} ${REGIONS_SHA256}"
  exit 0
fi

# shellcheck source=scripts/ci/lib/fetch-pinned.sh
. "${REPO}/scripts/ci/lib/fetch-pinned.sh"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
DIR=${SV_RULE_CACHE:-$WORK}
[ -z "$CHECK" ] || OUT="${WORK}/rebuilt.tsv"

fetch_pinned "$DIR" trx.bed.gz "$TRX_URL" "$TRX_SHA256" || exit 1
fetch_pinned "$DIR" gencode.gtf.gz "$GENCODE_URL" "$GENCODE_SHA256" || exit 1
fetch_pinned "$DIR" hgnc.txt "$HGNC_URL" "$HGNC_SHA256" || exit 1
fetch_pinned "$DIR" regions.txt "$REGIONS_URL" "$REGIONS_SHA256" || exit 1

"${REPO}/bin/og" dev build sv-rule-regions --work "$DIR" --out "$OUT"
if [ -n "$CHECK" ]; then
  if cmp -s "$OUT" "$CHECK"; then
    echo "${CHECK} is what the pinned downloads give, byte for byte" >&2
    exit 0
  fi
  echo "ERROR: ${CHECK} is not what the pinned downloads give; rebuild it with scripts/ci/build-sv-rule-regions.sh. The first differences (< the file, > the rebuild):" >&2
  diff "$CHECK" "$OUT" | head -n 20 >&2 || true
  exit 1
fi
echo "Wrote ${OUT}" >&2
