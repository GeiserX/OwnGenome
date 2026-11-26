#!/usr/bin/env bash
# What samtools view returns on each region of tests/data/regions.tsv, in the
# line format of `og bam-fetch-probe` and of tests/data/fetch.pysam.tsv: the
# region, then the read's name, flag, 0-based start, end (None when it has
# no CIGAR), MAPQ, CIGAR (None when it has none) and its aligned blocks, one
# per M, = or X operation, a D or N moving the reference on between them.
#
#   rust/crates/og-bam/tests/samtools_fetch.sh SAMTOOLS_IMAGE > fetch.tsv
#
# Run from the repository root. A region is 0-based and half-open, as
# pysam's fetch takes it, so samtools gets start+1 to stop.
# OG_FETCH_PLANT=deletion plants one mistake (a D no longer moves the
# reference on), for the workflow's control.
set -euo pipefail
d=rust/crates/og-bam/tests/data
img=${1:?Usage: $0 SAMTOOLS_IMAGE}
args=()
while IFS=$'\t' read -r c s e; do
  case "$c" in "" | "#"*) continue ;; esac
  args+=("$c:$s-$e" "$c:$((s + 1))-$e")
done < "${d}/regions.tsv"
[ "${#args[@]}" -gt 0 ] || { echo "ERROR: ${d}/regions.tsv has no region" >&2; exit 1; }
moves="DN"
[ "${OG_FETCH_PLANT:-}" != deletion ] || moves="N"
# shellcheck disable=SC2016 # the inner script expands in the container
docker run --rm --network none -v "${PWD}:/w:ro" -w /w "$img" sh -c '
  b=$1; shift
  while [ $# -gt 0 ]; do echo "@region $1"; samtools view "$b" "$2"; shift 2; done
' _ "${d}/cyp2d6.bam" "${args[@]}" < /dev/null | awk -F '\t' -v OFS='\t' -v moves="$moves" '
  /^@region / { region = substr($0, 9); next }
  {
    start = $4 - 1; cigar = $6
    if (cigar == "*") { print region, $1, $2, start, "None", $5, "None", ""; next }
    p = start; blocks = ""; rest = cigar
    while (match(rest, /^[0-9]+[MIDNSHP=X]/)) {
      n = substr(rest, 1, RLENGTH - 1) + 0; op = substr(rest, RLENGTH, 1)
      rest = substr(rest, RLENGTH + 1)
      if (op == "M" || op == "=" || op == "X") {
        blocks = blocks (blocks == "" ? "" : ",") p "-" (p + n); p += n
      } else if (index(moves, op)) p += n
    }
    if (rest != "") { print "ERROR: a CIGAR this script cannot read: " cigar > "/dev/stderr"; exit 1 }
    print region, $1, $2, start, p, $5, cigar, blocks
  }'
