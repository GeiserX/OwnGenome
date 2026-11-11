#!/bin/sh
# The oracle of og table constraint-join (tests/parity/oracle.ref):
# bin/constraint_join.awk, given og's options as the awk variables they
# replace.
#   sh ./oracle.sh PROG --gene-col NAME [--constrained] CONSTRAINT TABLE
prog=$1
shift
set -- "$@" --end
while [ "$1" != --end ]; do
  case "$1" in
    --gene-col) set -- "$@" "gene_col=$2"; shift ;;
    --constrained) set -- "$@" constrained=1 ;;
    *) set -- "$@" "$1" ;;
  esac
  shift
done
shift
LC_ALL=C exec awk -f "$prog" "$@"
