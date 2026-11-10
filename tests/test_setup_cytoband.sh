#!/usr/bin/env bash
# setup installs the chromosome bands from the copy the repository ships
# (assets/reference/cytoBand.hg38.txt), so step 10, --cytoband and the genome
# map read one file. The test runs offline: curl and wget on PATH fail and
# log each call.
#   1. install_data_file cytoband writes bytes equal to the asset and calls
#      neither curl nor wget. Before this change it downloaded UCSC's
#      cytoBand.txt.gz, so here it failed.
#   2. validate-setup.sh passes the installed copy.
#   3. A planted copy that differs from the asset by one byte makes
#      validate-setup.sh print a FAIL for the bands and exit 1.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
ASSET="${REPO}/assets/reference/cytoBand.hg38.txt"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
FAILS=0
ok() { echo "ok: $1"; }
bad() { echo "FAIL: $1"; FAILS=$((FAILS + 1)); }

mkdir -p "${WORK}/bin" "${WORK}/genome"
for tool in curl wget; do
  printf '#!/bin/sh\necho "%s $*" >> "%s/downloads.log"\nexit 7\n' "$tool" "$WORK" > "${WORK}/bin/${tool}"
  chmod +x "${WORK}/bin/${tool}"
done
export PATH="${WORK}/bin:${PATH}" GENOME_DIR="${WORK}/genome" FETCH_WAIT=0 FETCH_TRIES=1

[ -s "$ASSET" ] || { echo "FAIL: ${ASSET} is missing"; exit 1; }

# 1. Install offline.
rc=0
out=$(bash -c '. "$1/scripts/lib/common.sh" && install_data_file cytoband && data_file cytoband' _ "$REPO" 2>"${WORK}/install.err") || rc=$?
if [ "$rc" -ne 0 ]; then
  bad "install_data_file cytoband failed offline (exit ${rc}): $(tr '\n' ' ' < "${WORK}/install.err")"
elif [ "$out" != "${GENOME_DIR}/reference/cytoBand.hg38.txt" ]; then
  bad "data_file cytoband printed '${out}'"
elif ! cmp -s "$ASSET" "$out"; then
  bad "the installed bands differ from ${ASSET}"
else
  ok "setup installed the bands byte for byte from the asset"
fi
if [ -s "${WORK}/downloads.log" ]; then
  bad "setup tried to download the bands: $(tr '\n' ' ' < "${WORK}/downloads.log")"
else
  ok "no download for the bands"
fi

# validate LABEL: run validate-setup.sh and keep the line about the bands.
validate() {
  local rc=0
  bash "${REPO}/scripts/validate-setup.sh" > "${WORK}/validate.$1" 2>&1 || rc=$?
  echo "$rc"
}

# 2. The installed copy passes.
if [ -s "${GENOME_DIR}/reference/cytoBand.hg38.txt" ]; then
  validate good > /dev/null
  if grep -q '\[OK\].*GRCh38 chromosome bands: present' "${WORK}/validate.good"; then
    ok "validate-setup passes the installed bands"
  else
    bad "validate-setup did not pass the installed bands: $(grep -i 'chromosome bands' "${WORK}/validate.good" || echo 'no line')"
  fi
fi

# 3. A planted copy one byte off fails.
mkdir -p "${GENOME_DIR}/reference"
sed '1s/p36\.33/p36.32/' "$ASSET" > "${GENOME_DIR}/reference/cytoBand.hg38.txt"
if cmp -s "$ASSET" "${GENOME_DIR}/reference/cytoBand.hg38.txt"; then
  bad "the planted copy is the same as the asset; the sed no longer matches its first row"
else
  rc=$(validate planted)
  if grep -q '\[FAIL\].*GRCh38 chromosome bands.*not the copy in assets/reference' "${WORK}/validate.planted" && [ "$rc" -eq 1 ]; then
    ok "validate-setup fails a copy that differs from the asset (exit ${rc})"
  else
    bad "validate-setup did not refuse the planted copy (exit ${rc}): $(grep -i 'chromosome bands' "${WORK}/validate.planted" || echo 'no line')"
  fi
fi

if [ "$FAILS" -ne 0 ]; then
  echo "${FAILS} check(s) failed"
  exit 1
fi
echo "all checks passed"
