#!/usr/bin/env bash
# Extract the archived heredoc, not a second maintained implementation.
set -euo pipefail
mode=$1 prog=$2; shift 2
case "$mode" in
  diplotypes) start='^import json, sys$' ;;
  called-lines) start='^import json, sys$' ;;
  called) start='^import json, sys$' ;;
  csv-cell) start='^import csv, sys$' ;;
  mini-ref) start='^import random, sys$' ;;
  long-reads) start='^import gzip, sys$' ;;
  eh-catalog) start='^import json, re, sys$' ;;
  requests) start='^import csv, sys$' ;;
  attempts) start='^import csv, sys$' ;;
  *) exit 2 ;;
esac
# The workflow's first inline reader is requests, the second is attempts.
occurrence=1
[ "$mode" != attempts ] || occurrence=2
awk -v re="$start" -v want="$occurrence" '
  { sub(/^          /, "") }
  $0 ~ re { n++; if (n == want) on=1 }
  on && /^PY$/ { exit }
  on { print }
' "$prog" > oracle-body
case "$mode" in
  called|called-lines) set -- "$1" ;;
  requests) set -- "$1" "$5" "$7" ;;
  attempts) set -- "$1" "$7" ;;
esac
rc=0
python3 -I oracle-body "$@" || rc=$?
rm -f oracle-body
exit "$rc"
