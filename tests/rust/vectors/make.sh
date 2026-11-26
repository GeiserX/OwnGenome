#!/usr/bin/env bash
# RECORD (this commit only). Writes tests/rust/vectors/pycompat/cases.tsv.gz
# and tree.tsv.gz from tests/rust/pycompat_oracle.py in PYTHON_IMAGE: the
# hand-picked hard cases of every operation, then 1000 random cases per
# operation with seed 20261126, and the tree the tree operations walked,
# listed by GNU find. This is step 1 and the pycompat part of step 3 of
# make.sh at commit 00596ca8; chars.tsv.gz still comes from that commit.
# Needs docker; no network.
#
# Usage: tests/rust/vectors/make.sh WORK   (an empty scratch folder)
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
OUT="${ROOT}/tests/rust/vectors"
WORK=$(cd "${1:?Usage: $0 WORK}" && pwd)
# shellcheck source=../../../versions.env
. "${ROOT}/versions.env"
mkdir -p "${OUT}/pycompat"

docker run --rm --network none -e PYCOMPAT_TREE=/tree "$PYTHON_IMAGE" sh -c \
  'mkdir /tree && python3 -I -c "$1" 1000 20261126 > /tmp/pairs && find /tree -mindepth 1 -printf "%y\t%P\t%l\t%T@\n" > /tmp/tree && find /tree -maxdepth 0 -printf "%y\t\t\t%T@\n" >> /tmp/tree && cat /tmp/pairs && echo @@TREE@@ && cat /tmp/tree' \
  _ "$(cat "${ROOT}/tests/rust/pycompat_oracle.py")" > "${WORK}/oracle.out"
LC_ALL=C sed '/^@@TREE@@$/,$d' "${WORK}/oracle.out" > "${WORK}/pairs.txt"
LC_ALL=C sed '1,/^@@TREE@@$/d' "${WORK}/oracle.out" > "${WORK}/tree.txt"

LC_ALL=C awk 'NR % 2 == 1 { chars = ($0 ~ /^chars\t/) } { print (NR % 2 == 0 && chars) ? "@CHARS@" : $0 }' \
  "${WORK}/pairs.txt" | gzip -9 -n > "${OUT}/pycompat/cases.tsv.gz"
gzip -9 -n < "${WORK}/tree.txt" > "${OUT}/pycompat/tree.tsv.gz"
ls -l "${OUT}/pycompat"
