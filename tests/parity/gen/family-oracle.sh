#!/usr/bin/env bash
# The oracle of og family: bin/family_report.py imports bin/family_segments.py,
# and oracle-diff.sh materialises one file per tool. This stages the second
# file, at the commit tests/parity/oracle.ref names, and oracle.sh, which
# puts both beside each other and runs the one og's verb names.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
sha=$(awk -F'\t' '$1 == "family" { print $2 }' "${ROOT}/tests/parity/oracle.ref")
[ -n "$sha" ] || { echo "family-oracle.sh: no family row in oracle.ref" >&2; exit 1; }
git -C "$ROOT" cat-file -e "${sha}^{commit}" 2> /dev/null || git -C "$ROOT" fetch -q --depth=1 origin "$sha"
git -C "$ROOT" show "${sha}:bin/family_segments.py" > family_segments.py
cat > oracle.sh <<'SH'
# sh ./oracle.sh PROG ARGS: PROG is family_report.py; family_segments.py is staged here.
prog=$1
shift
d=$(mktemp -d)
cp "$prog" "${d}/family_report.py"
cp ./family_segments.py "${d}/family_segments.py"
script=family_report.py
if [ "${1:-}" = segments ]; then
  shift
  script=family_segments.py
fi
rc=0
PYTHONDONTWRITEBYTECODE=1 python3 "${d}/${script}" "$@" || rc=$?
rm -rf "$d"
exit "$rc"
SH
