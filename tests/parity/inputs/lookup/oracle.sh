#!/usr/bin/env bash
# oracle.sh PROG ARGS: the oracle of `og lookup query`. It does what
# scripts/owngenome-lookup.sh did with bin/genome_lookup.py (PROG) before
# og: for each query, `resolve`, then `report` for each region with the
# whole files (bcftools only cut slices of them; the report kept the same
# records either way), a blank line after each answer and after each query
# that does not resolve, and the exit code of the last one that failed.
set -uo pipefail
prog=$1
shift
gvcf="" ann="" summary="" sites="" genes=""
tables=()
while [ $# -gt 0 ]; do
  case "$1" in
    --gvcf) gvcf=$2; shift 2 ;;
    --annotated) ann=$2; shift 2 ;;
    --summary) summary=$2; shift 2 ;;
    --sites) sites=$2; shift 2 ;;
    --genes) genes=$2; shift 2 ;;
    --table) tables+=(--table "$2"); shift 2 ;;
    --) shift; break ;;
    *) break ;;
  esac
done
note() {
  if [ ! -s "$1.tbi" ] && [ ! -s "$1.csi" ]; then
    echo "NOTE: $1 has no .tbi index; reading the whole file (bcftools index -t makes this fast)." >&2
  fi
}
rc=0
for q in "$@"; do
  res=(${tables[@]+"${tables[@]}"})
  [ -z "$sites" ] || res+=(--sites "$sites")
  [ -z "$genes" ] || res+=(--genes "$genes")
  [ -z "$ann" ] || res+=(--annotated "$ann")
  qrc=0
  regions=$(python3 "$prog" resolve --query "$q" "${res[@]}") || qrc=$?
  if [ "$qrc" -ne 0 ]; then
    rc=$qrc
    echo ""
    continue
  fi
  while IFS=$'\t' read -r -u 3 kind chrom start end ref source; do
    note "$gvcf"
    rep=(--gvcf "$gvcf")
    if [ -n "$ann" ]; then
      note "$ann"
      rep+=(--annotated "$ann")
    fi
    [ -z "$summary" ] || rep+=(--summary "$summary")
    python3 "$prog" report --query "$q" --kind "$kind" --chrom "$chrom" --start "$start" --end "$end" \
      --ref "$ref" --source "$source" "${rep[@]}" || rc=$?
    echo ""
  done 3<<< "$regions"
done
exit "$rc"
