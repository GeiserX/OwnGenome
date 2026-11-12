#!/bin/sh
# Written by og ancestry-fine sample-bed. Usage: sh run_admixture.sh THREADS
set -eu
t=${1:-1}
while read -r level k; do
  admixture -P -j"$t" "S_${level}.bed" "$k" > "S_${level}.log" 2>&1 \
    || { tail -n 20 "S_${level}.log"; exit 1; }
  : > "S_${level}.boot.Q"
  r=1
  while [ "$r" -le 7 ]; do
    awk -v r="$r" -v nb=6 -v k="$k" '
      FNR == 1 { f++ }
      f == 1 { g[n++] = $0; next }
      f == 2 { m[n2++] = $0; next }
      f == 3 { p[n3++] = $0; next }
      END {
        srand(r); bs = int((n + nb - 1) / nb); nblk = int((n + bs - 1) / bs); printf "S S 0 0 0 -9" > "rep.ped"
        for (b = 0; b < nb; b++) {
          s = int(rand() * nblk) * bs
          for (i = s; i < s + bs && i < n; i++) { printf " %s", g[i] > "rep.ped"; print m[i] > "rep.map"; print p[i] > ("rep." k ".P.in") }
        }
        printf "\n" > "rep.ped"
      }' "S.geno" "S.map" "S_${level}.${k}.P.in"
    admixture -P "rep.ped" "$k" > rep.log 2>&1 || { tail -n 20 rep.log; exit 1; }
    cat "rep.${k}.Q" >> "S_${level}.boot.Q"
    rm -f rep.ped rep.map "rep.${k}.P.in" "rep.${k}.P" "rep.${k}.Q"
    r=$((r + 1))
  done
done < levels.txt
