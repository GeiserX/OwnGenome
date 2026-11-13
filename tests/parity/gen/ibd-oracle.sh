#!/usr/bin/env bash
# The oracle of og family ibd: bin/ibd_report.py imports bin/family_report.py,
# and oracle-diff.sh materialises one file per tool. This stages the second
# file, at the commit tests/parity/oracle.ref names, and ibd_oracle.sh, which
# puts both beside each other and runs the first.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
sha=$(awk -F'\t' '$1 == "family-ibd" { print $2 }' "${ROOT}/tests/parity/oracle.ref")
[ -n "$sha" ] || { echo "ibd-oracle.sh: no family-ibd row in oracle.ref" >&2; exit 1; }
git -C "$ROOT" cat-file -e "${sha}^{commit}" 2> /dev/null || git -C "$ROOT" fetch -q --depth=1 origin "$sha"
git -C "$ROOT" show "${sha}:bin/family_report.py" > family_report.py
cat > ibd_oracle.sh <<'SH'
# sh ./ibd_oracle.sh PROG ARGS: PROG is ibd_report.py; family_report.py is staged here.
prog=$1
shift
d=$(mktemp -d)
cp "$prog" "${d}/ibd_report.py"
cp ./family_report.py "${d}/family_report.py"
rc=0
PYTHONDONTWRITEBYTECODE=1 python3 "${d}/ibd_report.py" "$@" || rc=$?
rm -rf "$d"
exit "$rc"
SH
