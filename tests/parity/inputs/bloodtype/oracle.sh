#!/bin/sh
# The oracle of og bloodtype (tests/parity/oracle.ref): bin/bloodtype_report.py
# with the copy_number.py it imported (staged by gen/bloodtype-oracle.sh).
#   sh ./oracle.sh PROG report [ARGS]
d=$(mktemp -d) || exit 70
trap 'rm -rf "$d"' EXIT
cp "$1" "${d}/bloodtype_report.py"
cp copy_number.py "$d"/
shift
[ "$1" = report ] || { echo "oracle.sh: only report has an oracle" >&2; exit 2; }
shift
python3 -B "${d}/bloodtype_report.py" "$@"
