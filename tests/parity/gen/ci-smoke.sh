#!/usr/bin/env bash
set -euo pipefail
awk 'BEGIN {print ">chr20"; for(i=0;i<500000;i++) {printf "%s",(i>=21000 && i<21024 ? (i%2==0 ? "C":"A") : "G"); if(i%60==59)printf "\n"} printf "\n"}' > slice.fa
for h in 1 2; do awk 'BEGIN {print ">chr20";for(i=0;i<11000;i++){printf "A";if(i%60==59)printf "\n"}printf "\n"}' > "hap${h}.fa"; done
