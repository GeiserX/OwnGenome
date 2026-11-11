#!/usr/bin/env bash
# sex-loss-genome.sh: writes genome.bed.gz in the current directory, the
# same every time: 500 bp windows over GRCh38's chr1-22, X and Y (about 6.2
# million lines, the size of a real mosdepth regions file), for the timing
# of the qc-sex-loss parity case genome. scripts/ci/oracle-diff.sh runs it.
set -euo pipefail
export LC_ALL=C
# A whole genome of 500 bp windows: GRCh38's chr1-22, X and Y at 30x, with
# one Y in 80% of the cells (a loss of 0.20); depths from a fixed LCG whose
# products stay below 2^53, so every awk computes the same file.
awk 'BEGIN {
  n = split("248956422 242193529 198295559 190214555 181538259 170805979 159345973 145138636 138394717 133797422 135086622 133275309 114364328 107043718 101991189 90338345 83257441 80373285 58617616 64444167 46709983 50818468", len, " ")
  seed = 12345
  for (c = 1; c <= 24; c++) {
    if (c <= 22) { name = "chr" c; L = len[c]; mean = 30 }
    else if (c == 23) { name = "chrX"; L = 156040895; mean = 15 }
    else { name = "chrY"; L = 57227415; mean = 12 }
    for (s = 0; s < L; s += 500) {
      seed = (seed * 22695477 + 1) % 16777216
      e = s + 500 > L ? L : s + 500
      printf "%s\t%d\t%d\t%.2f\n", name, s, e, mean * (0.85 + 0.3 * seed / 16777216)
    }
  }
}' | gzip -1 > genome.bed.gz
