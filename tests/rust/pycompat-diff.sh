#!/usr/bin/env bash
# pycompat-diff.sh: the pycompat differential test. CPython (python3 -I in
# PYTHON_IMAGE) draws N random cases per operation and writes its result for
# each; `og pycompat-probe` must write the same results byte for byte.
#
# Usage: tests/rust/pycompat-diff.sh OG_BINARY [N] [SEED] [OP ...]
#   N defaults to 100000, SEED to a random one (printed, so a failure can be
#   replayed). PYTHON=python3 runs the oracle on the host instead of in the
#   image, for a quick local loop; CI always uses the image.
# Exit 0 when every result matches, 1 with the first differences otherwise.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
OG=${1:?Usage: $0 OG_BINARY [N] [SEED] [OP ...]}
N=${2:-100000} SEED=${3:-$RANDOM$RANDOM}
shift $(($# < 3 ? $# : 3))
# shellcheck source=../../versions.env
. "${ROOT}/versions.env"
TMP=$(mktemp -d "${TMPDIR:-/tmp}/pycompat-diff.XXXXXX")
trap 'rm -rf "$TMP"' EXIT
echo "pycompat-diff: ${N} cases per operation, seed ${SEED}, oracle ${PYTHON:-${PYTHON_IMAGE}}"
t0=$(date +%s)
# The oracle goes in as an argument and its output comes back on stdout, so
# nothing is mounted (Docker on macOS shares only some folders).
ORACLE=$(cat "${ROOT}/tests/rust/pycompat_oracle.py")
if [ -n "${PYTHON:-}" ]; then
  "$PYTHON" -I -c "$ORACLE" "$N" "$SEED" "$@" > "${TMP}/both.txt"
else
  docker run --rm --network none "$PYTHON_IMAGE" python3 -I -c "$ORACLE" "$N" "$SEED" "$@" > "${TMP}/both.txt"
fi
awk 'NR % 2 == 1' "${TMP}/both.txt" > "${TMP}/cases.tsv"
awk 'NR % 2 == 0' "${TMP}/both.txt" > "${TMP}/expected.txt"
t1=$(date +%s)
"$OG" pycompat-probe < "${TMP}/cases.tsv" > "${TMP}/actual.txt"
t2=$(date +%s)
cases=$(wc -l < "${TMP}/cases.tsv" | tr -d ' ')
[ "$cases" -gt 0 ] || { echo "FAIL: the oracle wrote no cases"; exit 1; }
echo "pycompat-diff: ${cases} cases; CPython $((t1 - t0)) s, og $((t2 - t1)) s"
if cmp -s "${TMP}/expected.txt" "${TMP}/actual.txt"; then
  echo "PASS: og pycompat-probe matches CPython on all ${cases} cases"
  exit 0
fi
paste "${TMP}/cases.tsv" "${TMP}/expected.txt" "${TMP}/actual.txt" | awk -F'\t' '
  { want = $(NF - 1); got = $NF; if (want == got) next; bad[$1]++; n++
    if (n <= 20) { c = $1; for (i = 2; i < NF - 1; i++) c = c "\t" $i; printf "case:     %s\npython:   %s\nog:       %s\n\n", c, want, got } }
  END { for (op in bad) printf "FAIL: %s differs on %d cases\n", op, bad[op] }'
echo "FAIL: og pycompat-probe differs from CPython (seed ${SEED}; rerun with that seed to replay)"
exit 1
