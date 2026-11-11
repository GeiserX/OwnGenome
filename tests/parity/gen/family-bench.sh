#!/usr/bin/env bash
# A synthetic joint VCF for og family report at scale, five samples: GF and
# GM, their son DAD, MOM, and their daughter KID, so two children are
# mapped, KID's father's copy told apart by the grandparents. FAMILY_SITES
# sites (default 200,000) over chr1-22 and chrX, every 100 bases; KID's
# grandparental copy switches at random (about one site in 25), which gives
# the smoothing many short runs to absorb. Genotypes come from awk's rand
# (seed 1), so they differ between awk builds; both sides of a case read the
# same file. Writes bench.vcf.gz.
set -euo pipefail
awk -v n="${FAMILY_SITES:-200000}" 'BEGIN {
  srand(1); OFS = "\t"
  per = int(n / 23)
  print "##fileformat=VCFv4.2"
  for (c = 1; c <= 23; c++) printf "##contig=<ID=%s,length=%d>\n", name(c), per * 100 + 1000
  print "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GF", "GM", "DAD", "MOM", "KID"
  for (c = 1; c <= 23; c++) {
    copy = 0
    for (i = 1; i <= per; i++) {
      for (k = 0; k < 8; k++) h[k] = rand() < 0.4 ? 1 : 0
      if (rand() < 0.04) copy = 1 - copy
      d0 = h[0]; d1 = h[3]
      kd = copy ? d1 : d0
      print name(c), i * 100, ".", "A", "G", "50", "PASS", ".", "GT:DP:GQ:AD", g(h[0], h[1]), g(h[2], h[3]), g(d0, d1), g(h[4], h[5]), g(kd, h[5])
    }
  }
}
function name(c) { return c == 23 ? "chrX" : "chr" c }
function g(a, b) { return a "/" b ":30:50:15,15" }' | gzip -n -1 > bench.vcf.gz
