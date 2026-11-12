#!/usr/bin/env bash
# A three-generation family for og family's chromosome map (the fixture of
# tests/test_family_segments.py): grandparents gpa and gma, their son dad,
# mum, and their son kid. chr1 has a crossover in the father's meiosis at
# 20 Mb and three planted genotype errors, chr2 a 10 Mb gap with no site,
# chrX the son's single X from his mother. Writes 3gen.vcf (its genotypes
# drawn with awk's rand, so they differ between awk builds; both sides of a
# case read the same file) and 3gen.vcf.gz.
set -euo pipefail
awk 'BEGIN {
  srand(38); OFS = "\t"
  print "##fileformat=VCFv4.2"
  print "##contig=<ID=chr1,length=248956422>"
  print "##contig=<ID=chr2,length=242193529>"
  print "##contig=<ID=chrX,length=156040895>"
  print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
  print "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "gpa", "gma", "dad", "mum", "kid"
  span("chr1", 1000000, 41000000); span("chr2", 1000000, 11000000); span("chr2", 21000000, 31000000)
  span("chrX", 5000000, 15000000)
}
function gt(a, b) { return a <= b ? a "/" b : b "/" a }
function r() { return rand() < 0.5 ? 0 : 1 }
function span(chrom, lo, hi,   pos, gpa1, gpa2, gma1, gma2, mum1, mum2, planted, from_dad) {
  for (pos = lo; pos < hi; pos += 2000) {
    gpa1 = r(); gpa2 = r(); gma1 = r(); gma2 = r(); mum1 = r(); mum2 = r()
    planted = chrom == "chr1" && (pos == 5000000 || pos == 9000000 || pos == 13000000)
    if (planted) { gpa1 = 0; gpa2 = 0; gma1 = 1; gma2 = 1; mum1 = 0; mum2 = 0 }
    if (chrom == "chrX") {
      print chrom, pos, ".", "A", "G", "50", "PASS", ".", "GT", gpa1, gt(gma1, gma2), gma1, gt(mum1, mum2), mum1
      continue
    }
    from_dad = (chrom == "chr1" && pos < 20000000) ? gpa1 : gma1
    if (planted) from_dad = gma1
    print chrom, pos, ".", "A", "G", "50", "PASS", ".", "GT", gt(gpa1, gpa2), gt(gma1, gma2), gt(gpa1, gma1), gt(mum1, mum2), gt(from_dad, mum1)
  }
}' > 3gen.vcf
gzip -n -c 3gen.vcf > 3gen.vcf.gz
