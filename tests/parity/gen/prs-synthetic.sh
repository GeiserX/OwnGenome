#!/usr/bin/env bash
# prs-synthetic.sh: about 1.5 million scoring-file rows for og prs format,
# written into the current folder (tests/parity/cases.tsv, prs-format).
#
#   PGS9000NN.txt.gz         eight GRCh38-harmonised scoring files of 1,000 to
#                           400,000 rows. Their positions come from one pool
#                           of 300,000 sites, so most alleles are shared, as
#                           in the large PGS Catalog scores; rows on chrX,
#                           chrY and MT, rows the catalog could not place,
#                           strand-ambiguous pairs, indels, lower case, a
#                           chr prefix, other alleles as lists or missing
#                           (then hm_inferOtherAllele), positions with a
#                           leading zero and positions past 32 bits (each
#                           apart from every other, so no two rows tie in
#                           the sort; Python orders a tie by its set).
#   synth.fa, synth.fa.fai  chr1 to chr22 and chrX, 1 Mb each, the base at
#                           each position a fixed function of it (some soft
#                           masked). PGS900001's single-base rows all carry
#                           the reference base, so it passes the strand
#                           check; the others fail it.
#   synth-labels.tsv        a label, a sex and evaluations for some ids.
#
# A fixed Park-Miller generator (products stay under 2^53), so mawk and gawk write the same files.
set -euo pipefail
awk 'BEGIN {
  B = "ACGT"
  for (c = 1; c <= 23; c++) {
    name = (c == 23) ? "chrX" : "chr" c
    print ">" name > "synth.fa"
    line = ""
    for (p = 1; p <= 1000000; p++) {
      b = substr(B, (p * 7 + c * 13 + int(p / 97)) % 4 + 1, 1)
      if (p % 1000 < 50) b = tolower(b)
      line = line b
      if (length(line) == 60) { print line > "synth.fa"; line = "" }
    }
    if (line != "") print line > "synth.fa"
  }
}'
# 1,000,000 bases in lines of 60: 16,667 lines of 61 bytes, the last one 41 bytes.
awk 'BEGIN { off = 0
  for (c = 1; c <= 23; c++) {
    name = (c == 23) ? "chrX" : "chr" c
    off += length(name) + 2
    printf "%s\t1000000\t%d\t60\t61\n", name, off
    off += 16666 * 61 + 41
  }
}' > synth.fa.fai

awk -v OFS='\t' '
function rnd(n) { seed = (seed * 16807) % 2147483647; return seed % n }
function ref(c, p) { return substr("ACGT", (p * 7 + c * 13 + int(p / 97)) % 4 + 1, 1) }
BEGIN {
  seed = 42
  split("1000 400000 300000 200000 100000 50000 10000 400000", size, " ")
  split("A C G T AT GCC TA a c g t N", al, " ")
  comp["A"] = "T"; comp["T"] = "A"; comp["C"] = "G"; comp["G"] = "C"
  for (i = 1; i <= 300000; i++) {
    k = rnd(100)
    pc[i] = (k < 92) ? 1 + rnd(22) : (k < 97) ? 23 : (k < 99) ? 24 : 25
    pp[i] = 1 + rnd(1000000)
  }
  for (s = 1; s <= 8; s++) {
    f = sprintf("PGS9000%02d.txt.gz", s)
    cmd = "gzip -c > " f
    printf "#pgs_id=PGS9000%02d\n#trait_reported=Synthetic trait %d\n#HmPOS_build=GRCh38\n", s, s | cmd
    print "rsID", "chr_name", "chr_position", "effect_allele", "other_allele", "effect_weight", "hm_chr", "hm_pos", "hm_inferOtherAllele" | cmd
    for (r = 1; r <= size[s]; r++) {
      i = 1 + rnd(300000)
      c = pc[i]; p = pp[i]
      chrom = (c <= 22) ? c : (c == 23) ? "X" : (c == 24) ? "Y" : "MT"
      if (rnd(10) == 0) chrom = "chr" chrom
      pos = p
      k = rnd(1000)
      if (k == 0) pos = "0" (p + 1000000)
      else if (k == 1) pos = 4294967296 + p
      else if (k == 2) { chrom = ""; pos = "" }
      ea = al[1 + rnd(12)]; oa = al[1 + rnd(12)]
      if (rnd(20) == 0) oa = comp[toupper(ea)] == "" ? oa : comp[toupper(ea)]
      if (s == 1 && c <= 23 && length(ea) == 1 && length(oa) == 1) {
        ea = ref(c, p)
        if (comp[ea] == oa) oa = (oa == "A" || oa == "T") ? "C" : "A"
      }
      inf = ""
      k = rnd(10)
      if (k == 0) { inf = oa; oa = "" }
      else if (k == 1) { inf = al[1 + rnd(4)]; oa = oa "/" al[1 + rnd(4)] }
      else if (k == 2) oa = "."
      w = (rnd(2) ? "-" : "") "0.00" rnd(100000)
      print "rs" r, "1", r, ea, oa, w, chrom, pos, inf | cmd
    }
    close(cmd)
  }
}'
printf 'pgs_id\ttrait_reported\tapplies_to\tppm_id\tancestry\teffect\nPGS900001\tStrand passes\t\tPPM000001\tEUR\tOR 1.2\nPGS900002\tFor women\tfemale\tPPM000002;PPM000003\tEUR;SAS\tHR 1.1;HR 1.3\nPGS900003\t\tmale\t\t\t\n' > synth-labels.tsv
