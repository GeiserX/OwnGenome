#!/bin/sh
# The oracle of og roh summary (tests/parity/oracle.ref): the summary block
# of scripts/11-roh-analysis.sh as it was before og took it over, its
# MASKED_BANDS line and its awk, given og's options.
#   sh ./oracle.sh PROG --sample SAMPLE --roh ROH [--out OUT]
prog=$1
out=/dev/stdout
shift
while [ $# -gt 1 ]; do
  case $1 in
    --sample) SAMPLE=$2 ;;
    --roh) OUTPUT=$2 ;;
    --out) out=$2 ;;
  esac
  shift 2
done
MASKED_BANDS=$(sed -n 's/^MASKED_BANDS="\(.*\)"$/\1/p' "$prog")
block=$(awk '/^SUMMARY=/ { on = 1; next } on && /^} > / { print "}"; exit } on { print }' "$prog")
[ -n "$MASKED_BANDS" ] && [ -n "$block" ] || { echo "oracle.sh: no summary block in $prog" >&2; exit 2; }
export SAMPLE OUTPUT MASKED_BANDS
eval "$block" > "$out"
