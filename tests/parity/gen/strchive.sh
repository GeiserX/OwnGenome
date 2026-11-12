#!/usr/bin/env bash
# strchive.sh: the two STRchive files versions.env pins (STRCHIVE_VERSION),
# downloaded into the current directory as eh.json and loci.json (fixed
# names, so a new release needs no case edit) and checked against their
# sha256, and assets/str_extra_loci.json beside them. scripts/ci/oracle-diff.sh runs
# it for the str-catalog cases on the real release.
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
# shellcheck source=../../../versions.env
. "${ROOT}/versions.env"
sha() { if command -v sha256sum > /dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }
for row in "eh.json STRchive-disease-loci-${STRCHIVE_VERSION}.hg38.expansionhunter.json ${STRCHIVE_EH_SHA256}" \
           "loci.json STRchive-loci-${STRCHIVE_VERSION}.json ${STRCHIVE_LOCI_SHA256}"; do
  read -r out name want <<< "$row"
  curl -sSLf --retry 3 -o "$out" "https://github.com/dashnowlab/STRchive/releases/download/${STRCHIVE_VERSION}/${name}"
  [ "$(sha "$out")" = "$want" ] || { echo "strchive.sh: ${name} does not have the sha256 versions.env pins" >&2; exit 1; }
done
# Named apart from the planted extra.json of tests/parity/inputs/str-catalog.
cp "${ROOT}/assets/str_extra_loci.json" str_extra_loci.json
