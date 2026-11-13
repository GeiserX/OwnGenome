#!/usr/bin/env bash
# Public catalogs and synthetic genotype grids for all three traits verbs.
set -euo pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd)
cp "$root/assets/traits.tsv" traits.tsv
cp "$root/assets/genotype_rules.tsv" rules.tsv
cat > oracle.sh <<'EOF'
prog=$1; shift
mode=$1; shift
case "$mode" in
 table) exec python3 "$prog" "$@" ;;
 check) exec python3 "$prog" "${@/--table/--check-table}" ;;
 phase-sites)
   args=()
   while [ "$#" -gt 0 ]; do
     if [ "$1" = --out ]; then args+=(--phase-sites "$2"); shift 2
     else args+=("$1"); shift; fi
   done
   exec python3 "$prog" "${args[@]}" ;;
esac
EOF
: > empty.tsv
awk -F'\t' '$1 !~ /^#/ && $1 != "chrom" {print $1 "\t" $2 "\t" $5 "\t" $6 "\t" $5 "/" $5 "\t40"}' traits.tsv > reference.tsv
for a in 0 1 2; do
 for b in 0 1 2; do
   case $a in 0) ga=T/T;; 1) ga=T/C;; 2) ga=C/C;; esac
   case $b in 0) gb=C/C;; 1) gb=C/T;; 2) gb=T/T;; esac
   printf 'chr19\t44908684\tT\tC\t%s\t40\nchr19\t44908822\tC\tT\t%s\t39\n' "$ga" "$gb" > "apoe-$a$b.tsv"
 done
done
printf 'chr19\t44908684\tT\tC\tT/G\t30\nchr19\t44908822\tC\tT\tC/C\t40\n' > other.tsv
printf 'chr19\t44908684\tT\tC\tT/C\t30\n' > missing.tsv
printf 'chr19\t44908684\tTT\tT\tTT/T\t10\nchr19\t44908684\tT\tC\tC/T\t40\n' > duplicate.tsv
printf 'chr19\t44908684\tG\tC\tG/C\t30\n' > wrong-ref.tsv
printf 'chr19\t44908684\tT\tC\t./C\t30\n' > half.tsv
printf 'chr7\t141973545\tC\tG\tC/G\t30\nchr7\t141972905\tG\tA\tG/A\t31\nchr7\t141972804\tT\tC\tT/C\t32\n' > triple.tsv
head -1 triple.tsv > single.tsv
: > empty.vcf
printf 'chr7\t141973545\t.\tC\tG\t.\tPASS\t.\tGT:PS\t1|0:7\nchr7\t141972905\t.\tG\tA\t.\tPASS\t.\tGT:PS\t0|1:7\nchr7\t141972804\t.\tT\tC\t.\tPASS\t.\tGT:PS\t1|0:7\n' > phased.vcf
sed '1s/:7$/:8/' phased.vcf > split.vcf
awk -F'\t' -v OFS='\t' '$1 == "combination" && $2 == "APOE" && $3 == "e3/e3" {next} {print}' rules.tsv > bad-rules.tsv
awk -F'\t' -v OFS='\t' '$3 == "rs671" {$6=$5} {print}' traits.tsv > bad-table.tsv
