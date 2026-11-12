#!/bin/sh
# The oracle of og cn (tests/parity/oracle.ref): bin/copy_number.py with the
# modules it imported (staged by gen/cn-oracle.sh), each og cn verb turned
# into the old command line.
#   sh ./oracle.sh PROG VERB [ARGS]
#   call ... --marker LPA_KIV2_4925 --marker-pileup P --marker-out O
#       copy_number.py call ..., then lpa_kiv2.py null-allele --table T
#   marker-regions NAME      traits_report.py --marker-regions NAME
#   marker-read NAME --pileup P [--copy-table T] --out O
#       traits_report.py --marker NAME ... --marker-out O
d=$(mktemp -d) || exit 70
trap 'rm -rf "$d"' EXIT
cp "$1" "${d}/copy_number.py"
cp hba_copy_number.py lpa_kiv2.py traits_report.py cn_loci.tsv "$d"/
shift
verb=$1
shift
case "$verb" in
  marker-regions) python3 -B "${d}/traits_report.py" --marker-regions "$@"; exit $? ;;
  marker-read)
    name=$1
    shift
    set -- "$@" --end
    while [ "$1" != --end ]; do
      case "$1" in
        --out) set -- "$@" --marker-out "$2" ;;
        --catalog) ;;
        *) set -- "$@" "$1" "$2" ;;
      esac
      shift 2
    done
    shift
    python3 -B "${d}/traits_report.py" --marker "$name" "$@"
    exit $?
    ;;
  call)
    pileup="" out="" table=""
    set -- "$@" --end
    while [ "$1" != --end ]; do
      case "$1" in
        --marker) ;;
        --marker-pileup) pileup=$2 ;;
        --marker-out) out=$2 ;;
        --table) table=$2; set -- "$@" "$1" "$2" ;;
        *) set -- "$@" "$1" "$2" ;;
      esac
      shift 2
    done
    shift
    python3 -B "${d}/copy_number.py" call "$@" || exit $?
    [ -z "$out" ] || python3 -B "${d}/lpa_kiv2.py" null-allele --pileup "$pileup" --table "$table" --out "$out"
    exit $?
    ;;
  *) python3 -B "${d}/copy_number.py" "$verb" "$@"; exit $? ;;
esac
