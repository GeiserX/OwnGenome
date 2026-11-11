#!/bin/sh
# The oracle of og clinvar deep-intronic (tests/parity/oracle.ref):
# bin/clinvar_hits.awk in its sites mode, given og's options as the awk
# variables it took.
#   sh ./oracle.sh PROG [--no-gvcf] [--min-gq N] SITES VCF_CALLS GVCF_CALLS
prog=$1
shift
gvcf=yes min_gq=30
for a in "$@"; do
  shift
  case "$a" in
    --no-gvcf) gvcf=no ;;
    --min-gq) min_gq=next ;;
    *) if [ "$min_gq" = next ]; then min_gq=$a; else set -- "$@" "$a"; fi ;;
  esac
done
exec awk -v mode=sites -v "gvcf=${gvcf}" -v "min_gq=${min_gq}" -f "$prog" "$@"
