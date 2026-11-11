#!/bin/sh
# The oracle of og vcf csq-rank (tests/parity/oracle.ref): bin/csq_rank.awk,
# given og's options as the awk variables they replace.
#   sh ./oracle.sh PROG [--keep K] [--min-gq N] ... [--in FILE] [--out FILE]
prog=$1
shift
in=- out=-
set -- "$@" --end
while [ "$1" != --end ]; do
  case "$1" in
    --keep) set -- "$@" -v "keep=$2" ;;
    --rare-af) set -- "$@" -v "rare_af=$2" ;;
    --min-gq) set -- "$@" -v "min_gq=$2" ;;
    --min-dp) set -- "$@" -v "min_dp=$2" ;;
    --het-lo) set -- "$@" -v "het_lo=$2" ;;
    --het-hi) set -- "$@" -v "het_hi=$2" ;;
    --hom-lo) set -- "$@" -v "hom_lo=$2" ;;
    --in) in=$2 ;;
    --out) out=$2 ;;
    *) echo "oracle.sh: unknown option $1" >&2; exit 2 ;;
  esac
  shift 2
done
shift
[ "$in" = - ] || exec < "$in"
[ "$out" = - ] || exec > "$out"
LC_ALL=C exec awk "$@" -f "$prog"
