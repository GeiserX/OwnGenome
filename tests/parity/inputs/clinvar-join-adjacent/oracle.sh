#!/bin/sh
# The oracle of og clinvar join-adjacent (tests/parity/oracle.ref):
# bin/clinvar_join_adjacent.awk, given og's options as the awk variables it
# took. --mode build is the pipe step 06 ran: prepare, then
# bcftools +fill-from-fasta, then finish. tag read its file twice.
#   sh ./oracle.sh PROG --mode M [--window N] [--fasta REF] [VCF]
prog=$1
shift
mode= window= fasta= file=
while [ $# -gt 0 ]; do
  case "$1" in
    --mode) mode=$2; shift ;;
    --window) window=$2; shift ;;
    --fasta) fasta=$2; shift ;;
    *) file=$1 ;;
  esac
  shift
done
set -- -f "$prog"
[ -z "$window" ] || set -- -v "window=${window}" "$@"
case "$mode" in
  tag) exec awk -v mode=tag "$@" "$file" "$file" ;;
  build)
    exec bash -o pipefail -c 'awk -v mode=prepare -v "window=$1" -f "$0" ${3:+"$3"} \
      | bcftools +fill-from-fasta --no-version -Ov - -- -c INFO/PGP_JOIN_SEQ -f "$2" \
      | awk -v mode=finish -f "$0"' "$prog" "$window" "$fasta" "$file" ;;
  *) exec awk -v "mode=${mode}" "$@" ${file:+"$file"} ;;
esac
