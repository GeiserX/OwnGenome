#!/usr/bin/env bash
# build-acmg-lof.sh: rebuild the clingen_hi and lof columns of
# assets/acmg_sf.tsv, which say whether a predicted loss-of-function call
# counts as a secondary finding in each ACMG gene (og summary collect).
#
# Usage: scripts/ci/build-acmg-lof.sh [--table FILE] [--cache DIR] [--check]
#        scripts/ci/build-acmg-lof.sh --list   the download, as "NAME URL SHA256"
#
# Rewrites the two cells of each gene row of FILE (default
# assets/acmg_sf.tsv) and nothing else, or with --check exits 1 when they
# are not what the download gives. The download is ClinGen's gene curation
# list for GRCh38 of 2026-09-07 (CC0, about 250 KB). ClinGen rewrites that
# file every day, so the pin is the Internet Archive's copy of that day,
# checked against its sha256. It is kept in --cache when given, else in a
# temporary folder that is deleted. The rules are in og dev build acmg-lof
# (rust/crates/og-dev/src/build/acmg_lof.rs); the table's header says what
# each column holds.
#
# Needs curl, og (scripts/setup.sh --og) and sha256sum or shasum.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
# shellcheck source=scripts/ci/lib/fetch-pinned.sh
. "${REPO}/scripts/ci/lib/fetch-pinned.sh"

TABLE=${REPO}/assets/acmg_sf.tsv
CACHE=""
CHECK=""
LIST=""
while [ $# -gt 0 ]; do
  case $1 in
    --table) TABLE=${2:?--table needs a file}; shift 2 ;;
    --cache) CACHE=${2:?--cache needs a folder}; shift 2 ;;
    --check) CHECK=1; shift ;;
    --list) LIST=1; shift ;;
    *) echo "Usage: $0 [--table FILE] [--cache DIR] [--check] | --list" >&2; exit 2 ;;
  esac
done
if [ -n "$LIST" ] && [ -n "$CHECK" ]; then
  echo "ERROR: --list and --check do not go together" >&2
  exit 2
fi

NAME=clingen_gene_curation_list_GRCh38.tsv
URL=https://web.archive.org/web/20260908054916id_/https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
SHA=a493366d01c1a0550a43dcd20e09da4807186749b71a39d6f97ad81ee6235de8

if [ -n "$LIST" ]; then
  echo "${NAME} ${URL} ${SHA}"
  exit 0
fi

if [ -z "$CACHE" ]; then
  CACHE=$(mktemp -d)
  trap 'rm -rf "$CACHE"' EXIT
fi
fetch_pinned "$CACHE" "$NAME" "$URL" "$SHA" || exit 2

# Not exec: the EXIT trap deletes a temporary download folder.
rc=0
"${REPO}/bin/og" dev build acmg-lof --clingen "${CACHE}/${NAME}" --table "$TABLE" ${CHECK:+--check} || rc=$?
exit "$rc"
