#!/usr/bin/env bash
# build-cn-intervals.sh: rewrite the interval rows of bin/cn_loci.tsv, the
# copy-number checks step 35 runs with one samtools bedcov pass
# (og cn intervals).
#
# Usage: scripts/ci/build-cn-intervals.sh [--catalog FILE] [--cache DIR] [--check]
#        scripts/ci/build-cn-intervals.sh --list   the downloads, one "NAME URL SHA256" per line
#
# Rewrites FILE (default bin/cn_loci.tsv) from its marker line on, or with
# --check exits 1 when it is not what the downloads give. The four downloads
# (about 150 KB, each checked against its sha256) are kept in --cache when
# given, else in a temporary folder that is deleted:
#
#   ClinGen Dosage Sensitivity's recurrent CNV BED, version 2.1, GRCh38
#     (CLINGEN_RECURRENT_CNV_URL and _SHA256 in versions.env, CC0)
#   GIAB genome stratifications v3.6, GRCh38: segmental duplications,
#     assembly gaps with 15 kb on each side, satellites
#
# with assets/parascopy_models_v1.7.txt. Every region row of the BED is a
# row, then the rows the BED lacks: 1p36 and Xq28 (MECP2) from ClinGen's
# curations, five carrier deletions (GJB6-D13S1830, CFTRdele2,3, three
# beta-globin deletions) and the Y chromosome's AZF regions. Targets are cut
# to pieces of at most 500 kb and lose the BED's breakpoint and LCR rows,
# segmental duplications (not on the AZF rows), gaps, satellites and
# Parascopy models; controls are two 100 kb flanks clear of all of that. The
# rules and the hand rows are in og dev build cn-intervals
# (rust/crates/og-dev/src/build/cn_intervals.rs); docs/cn-loci.md says what
# each row reads.
#
# Needs curl, og (scripts/setup.sh --og) and sha256sum or shasum. The e2e
# case tests/e2e/cn-intervals-1-catalog.sh runs --check.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
# shellcheck source=scripts/ci/lib/fetch-pinned.sh
. "${REPO}/scripts/ci/lib/fetch-pinned.sh"

CATALOG=${REPO}/bin/cn_loci.tsv
CACHE=""
CHECK=""
LIST=""
while [ $# -gt 0 ]; do
  case $1 in
    --catalog) CATALOG=${2:?--catalog needs a file}; shift 2 ;;
    --cache) CACHE=${2:?--cache needs a folder}; shift 2 ;;
    --check) CHECK=1; shift ;;
    --list) LIST=1; shift ;;
    *) echo "Usage: $0 [--catalog FILE] [--cache DIR] [--check] | --list" >&2; exit 2 ;;
  esac
done
if [ -n "$LIST" ] && [ -n "$CHECK" ]; then
  echo "ERROR: --list and --check do not go together" >&2
  exit 2
fi

# versions_env KEY: the value of KEY in versions.env.
versions_env() {
  local v
  v=$(sed -n "s/^$1=\"\([^\"]*\)\".*/\1/p" "${REPO}/versions.env" | head -n 1)
  [ -n "$v" ] || { echo "ERROR: no $1 in versions.env" >&2; return 1; }
  echo "$v"
}

GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all
clingen_url=$(versions_env CLINGEN_RECURRENT_CNV_URL) || exit 2
clingen_sha=$(versions_env CLINGEN_RECURRENT_CNV_SHA256) || exit 2
FILES=(
  "clingen_recurrent_cnv.bed|${clingen_url}|${clingen_sha}"
  "segdups.bed.gz|${GIAB}/SegmentalDuplications/GRCh38_segdups.bed.gz|09ac029a7d0d681a9f0ed6581a6bc10fc708a834bc3f91c46a9bef504f6055c5"
  "gaps.bed.gz|${GIAB}/OtherDifficult/GRCh38_gaps_slop15kb.bed.gz|b5fdbfc7e5b9cf4fa356f23841fd7b8a74fa239fd8cbc5b5c3463c9e7bdc0902"
  "satellites.bed.gz|${GIAB}/LowComplexity/GRCh38_satellites_slop5.bed.gz|10207aa2dac0e1189ea553a374abbdcb4de25a81378d88dc4ba03aea8718b5a2"
)

if [ -n "$LIST" ]; then
  for f in "${FILES[@]}"; do
    IFS='|' read -r name url sum <<< "$f"
    echo "${name} ${url} ${sum}"
  done
  exit 0
fi

if [ -z "$CACHE" ]; then
  CACHE=$(mktemp -d)
  trap 'rm -rf "$CACHE"' EXIT
fi
for f in "${FILES[@]}"; do
  IFS='|' read -r name url sum <<< "$f"
  fetch_pinned "$CACHE" "$name" "$url" "$sum" || exit 2
done

# Not exec: the EXIT trap deletes a temporary download folder.
rc=0
"${REPO}/bin/og" dev build cn-intervals \
  --clingen "${CACHE}/clingen_recurrent_cnv.bed" --segdups "${CACHE}/segdups.bed.gz" \
  --gaps "${CACHE}/gaps.bed.gz" --satellites "${CACHE}/satellites.bed.gz" \
  --parascopy "${REPO}/assets/parascopy_models_v1.7.txt" --catalog "$CATALOG" ${CHECK:+--check} || rc=$?
exit "$rc"
