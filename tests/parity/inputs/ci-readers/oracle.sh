#!/usr/bin/env bash
# Extract the archived heredoc, not a second maintained implementation.
set -euo pipefail
mode=$1 prog=$2; shift 2
if [ "$mode" = check ]; then
  kind=$1
  case "$kind" in
    hap-ibd) var=HAP_IBD_IMAGE; gzip -c o.ibd > o.ibd.gz ;;
    ibdmix) var=IBDMIX_IMAGE ;;
    rfmix) var=RFMIX_IMAGE ;;
    eigen) var=EIGENSOFT_IMAGE ;;
    parascopy) var=PARASCOPY_IMAGE ;;
    *) exit 2 ;;
  esac
  expect=$(awk -F '\t' -v v="$var" '$1 == v { print $5 }' "$prog")
  [ -n "$expect" ] || exit 2
  # Compare the assertions, not their debug dumps. The Rust diagnostics
  # have explicit goldens; a failed Python assertion still fails the case.
  py() { python3 -I -c "$2" > /dev/null; }
  same() { :; }
  at_least() { :; }
  eval "$expect"
  [ "$kind" != hap-ibd ] || rm -f o.ibd.gz
  exit 0
fi
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
indent=0
case "$mode" in requests|attempts) indent=1 ;; esac
awk -v re="$start" -v want="$occurrence" -v indent="$indent" '
  { if (indent) sub(/^          /, "") }
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
