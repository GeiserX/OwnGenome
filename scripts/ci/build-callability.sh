#!/usr/bin/env bash
# build-callability.sh: rebuild assets/callability_genes.tsv, the genes where
# short reads often cannot be placed, so a "no finding" there needs a reading
# of its own (step 16b measures each one, the report lists the ones it could
# not read).
#
# Usage: scripts/ci/build-callability.sh [OUT]   (default: assets/callability_genes.tsv)
#        scripts/ci/build-callability.sh --list    the downloads, one "NAME URL SHA256" per line
#
# Downloads seven pinned files (about 85 MB, each checked against its
# sha256) into CALLABILITY_CACHE, where a file with the right sha256 is read
# as it is, or else into a temporary folder that is deleted, and writes OUT:
#
#   GIAB CMRG v1.00's GRCh38 gene coordinates (Wagner 2022, PMID 35132260):
#     the 273 challenging medically relevant genes
#   GIAB genome stratifications v3.6, GRCh38: segmental duplications, the
#     false duplications of GRCh38 (both copies; Behera 2023, PMID 36810122)
#     and the collapsed duplications
#   GENCODE 50's basic annotation with the patches and ALT scaffolds
#   the GRCh38 assembly report of NCBI (GCA_000001405.15), which names the
#     ALT scaffolds
#
# A gene is a gene of GENCODE 50 with an HGNC id, on chr1-22, X or Y unless it
# is ALT-only, and protein-coding unless it is a CMRG gene (a few of those are
# T-cell receptor segments). It is listed with every reason that applies:
#
#   CMRG               a CMRG gene (matched by symbol, or by the protein-coding
#                      gene that covers most of its CMRG span when GENCODE
#                      has renamed it)
#   segdup             half or more of the gene body in a segmental duplication
#   false_duplication  the gene body overlaps a false duplication of GRCh38
#                      (either copy): reads of the real gene get MAPQ 0
#   collapsed          the gene body overlaps a collapsed duplication: two
#                      copies in people, one in GRCh38
#   ALT-only           the gene is on an ALT scaffold and nowhere on chr1-22,
#                      X, Y or M, so the default no-ALT reference cannot hold it
#
# and the step that calls it with a method of its own (owner): parascopy for
# the genes of assets/paralog_loci.tsv, pypgx for the bam genes of
# assets/pgx_genes.tsv that step 32 can call on the no-ALT reference,
# cn_loci:<locus> for a gene a gene row of bin/cn_loci.tsv names (its locus,
# the locus up to its first '_', or a target name), in that order; '.' for
# none. A change to one of those three files changes the owners, so rebuild
# after it.
#
# assets/paralog_loci.tsv names owners only; it is not a source of rows. The
# rows say where short reads fail by the measures above, gene body by gene
# body, and a Parascopy model may cover only part of a gene. Seventeen catalog
# genes are not rows: seven are pseudogenes (CELP, GBAP1, NCF1B, NCF1C,
# OTOAP1, PMS2CL, STRCP1), left out as not protein-coding, and no reason
# above applies to the other ten (ABCC6, CFH, CFHR4, CHRNA7, DEFA3, KRT86,
# NEB, PIK3CA, PLG, VWF). Step 35's own report section reads all of them.
# tests/test_paralog_loci.py holds that list and goes red when it changes.
#
# The rules are in og dev build callability
# (rust/crates/og-dev/src/build/callability.rs). Needs curl, og
# (scripts/setup.sh --og) and sha256sum or shasum. The e2e case
# tests/e2e/callability-1-catalog.sh runs it and fails on any difference with
# the committed file.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
OUT=${1:-${REPO}/assets/callability_genes.tsv}

GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release
STRAT=${GIAB}/genome-stratifications/v3.6/GRCh38@all
FILES=(
  "cmrg.bed|${GIAB}/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SupplementaryFiles/GRCh38_CMRG_benchmark_gene_coordinates.bed|67ce028e3acc0f4946104b8a6e0b76ce1ba9da8eb8457e3d7185bd1e81cbf203"
  "segdups.bed.gz|${STRAT}/SegmentalDuplications/GRCh38_segdups.bed.gz|09ac029a7d0d681a9f0ed6581a6bc10fc708a834bc3f91c46a9bef504f6055c5"
  "falsedup_correct.bed.gz|${STRAT}/OtherDifficult/GRCh38_false_duplications_correct_copy.bed.gz|44d07f373f6841a99fb85bd712e535e3e1364a033fd37cfd42e9e17988277717"
  "falsedup_incorrect.bed.gz|${STRAT}/OtherDifficult/GRCh38_false_duplications_incorrect_copy.bed.gz|d2de1abe40a76cb1d29dfa5ae4c0a9ae05f2301641721f4a78f6b72c3a0b5bcd"
  "collapsed.bed.gz|${STRAT}/OtherDifficult/GRCh38_collapsed_duplication_FP_regions.bed.gz|32b6f04e25afeb12078112a771f76e3f880f3a1d1bb2f454f60c5462644bf8af"
  "gencode.gtf.gz|https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.chr_patch_hapl_scaff.basic.annotation.gtf.gz|e9199c3ef075e54f339411a090385dd65ffe6ff59725403230cf1688ce9e8524"
  "report.txt|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_report.txt|aa733ce92719f6339b3a4540137d4322894434363314c4b7704ba6691e9eab66"
)

if [ "${1:-}" = --list ]; then
  for f in "${FILES[@]}"; do
    IFS='|' read -r name url sum <<<"$f"
    echo "${name} ${url} ${sum}"
  done
  exit 0
fi

# shellcheck source=scripts/ci/lib/fetch-pinned.sh
. "${REPO}/scripts/ci/lib/fetch-pinned.sh"
DIR=${CALLABILITY_CACHE:-}
if [ -z "$DIR" ]; then
  DIR=$(mktemp -d)
  trap 'rm -rf "$DIR"' EXIT
fi
for f in "${FILES[@]}"; do
  IFS='|' read -r name url sum <<<"$f"
  fetch_pinned "$DIR" "$name" "$url" "$sum" || exit 1
done

"${REPO}/bin/og" dev build callability --work "$DIR" --repo "$REPO" --out "$OUT"
echo "Wrote ${OUT}" >&2
