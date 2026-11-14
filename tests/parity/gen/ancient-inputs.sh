#!/usr/bin/env bash
# The inputs of the ancient cases of step 44 (ancient-select, ancient-lift,
# ancient-person-geno, ancient-align, ancient-report), all synthetic,
# written with awk:
#   aadr.anno, aadr.ind, few.ind
#             an AADR-style annotation (long header cells, a quoted one,
#             ASSESSMENT WARNINGS before ASSESSMENT) with present-day people
#             of the reference groups, one individual under two genetic
#             ids, ancient people inside and outside West Eurasia, below the
#             site floor, undated, failing the assessment; the .ind file
#             (unknown ids, short lines) and one that matches too few
#   hg19.chain, hg19.chain.gz, minus-source.chain, orphan.chain, hg38.fa(.fai), stale.fa(.fai)
#             a chain with a gap, a minus-strand block, two chains over one
#             stretch and a block that moves to another chromosome; a
#             reference with lower-case lines; one whose .fai is wrong
#   aadr.snp  sites in blocks, in gaps, on block edges, on the minus
#             strand, off the autosomes, with equal or non-base alleles and
#             with a fractional position
#   sites.tsv, genotypes.tsv, panel.snp, panel.geno, short.geno,
#   long.geno, wide.geno, nowhere.tsv
#             person-geno: two SNPs at one position, other alleles, missing
#             and haploid calls; a .geno that is short, long or uneven
#   panel.evec, run.evec, flipped.evec, bad.evec, few.evec, ancients.tsv, stats.tsv
#             two smartpca runs on one panel (the second with PC2's sign
#             flipped), one whose axes differ, one with too few present-day
#             people, and the ancient people's table (names to escape, a
#             long name, undated people, dates on both sides of year 1)
set -euo pipefail

# --- select -------------------------------------------------------------------------
awk 'BEGIN {
  OFS = "\t"
  print "Genetic ID", "Master ID", "Individual ID (PubMed)", \
    "Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date]", \
    "Date standard deviation in BP [OxCal sigma]", "Full Date One of two formats", \
    "Group ID", "\"Locality\"", "Political Entity", "Latitude", "Longitude", "Data type", \
    "SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)", \
    "ASSESSMENT WARNINGS (Xcontam interval is listed if lower bound is >0.005)", " ASSESSMENT " > "aadr.anno"
  split("CEU French Sardinian Basque GBR IBS TSI Druze Armenian Polish Mozabite Saami", M, " ")
  for (i = 1; i <= 90; i++) {
    g = M[i % 12 + 1]
    a = (i % 9 == 0) ? "QUESTIONABLE" : (i % 10 == 0 ? " pass (low coverage)" : "PASS")
    dt = (i % 15 == 0) ? "1240K" : "Shotgun.diploid"
    print "M" i, "M" i, "Ind_M" i, (i % 17 == 0 ? "-0.0" : "0"), "0", "present", g, "Paris\tFrance", "48.8", "2.3", dt, 900000 + i, "", a >> "aadr.anno"
  }
  print "M91", "M91", "Ind_M5", "0", "0", "present", "CEU", "Utah", "USA", "40", "-111", "Shotgun.diploid", 999999999, "", "PASS" >> "aadr.anno"
  split("Germany_CordedWare England_BellBeaker Turkey_N <Mixed&Group> Italy_Imperial_Rome Hungary_LBA_with_a_very_long_name_that_runs_past_45 Türkiye_Neolithic", G, " ")
  for (i = 1; i <= 60; i++) {
    g = G[i % 7 + 1]
    bp = (i % 13 == 0) ? ".." : 500 + 137 * i
    lat = (i % 11 == 0) ? 10 : 30 + i % 40
    lon = (i % 19 == 0) ? "" : -20 + i
    snps = (i % 8 == 0) ? 99999 : 100000 + 1000 * i
    a = (i % 14 == 0) ? "FAIL" : (i % 6 == 0 ? "PROVISIONAL_PASS" : "MERGE_PASS")
    ind = (i % 10 == 0) ? "" : "Ind_A" int(i / 2)
    print "A" i, "A" i, ind, bp, "50", "  2500-2000   BCE ", g, "Site " (i % 3), (i % 4 == 0 ? "" : "Country" (i % 2)), lat, lon, "1240K", snps, "x", a >> "aadr.anno"
  }
  print "short", "row" >> "aadr.anno"
  for (i = 1; i <= 91; i++) print "M" i, (i % 2 ? "M" : "F"), "Pop" > "aadr.ind"
  for (i = 1; i <= 60; i++) print "  A" i, "U", "Pop" >> "aadr.ind"
  print "X1 M" >> "aadr.ind"
  print "Unknown", "M", "Pop" >> "aadr.ind"
  for (i = 1; i <= 20; i++) print "M" i, "M", "Pop" > "few.ind"
}'

# --- lift ---------------------------------------------------------------------------
cat > hg19.chain <<'EOF'
chain 1000 chr1 100000 + 1000 1600 chr1 3000 + 2000 2600 1
200	100	100
300

chain 900 chr2 50000 + 500 900 chr2 3000 - 100 500 2
400

chain 800 chr3 50000 + 100 400 chr3 3000 + 100 400 3
300

chain 700 chr3 50000 + 300 500 chr3 3000 + 1300 1500 4
200

chain 600 chr4 50000 + 0 500 chr5 3000 + 0 500 5
500

chain 500 chr5 50000 + 100 1100 chr5 3000 + 1100 2100 6
1000
EOF
gzip -n -c hg19.chain > hg19.chain.gz
printf 'chain 1 chr1 100 - 0 10 chr1 100 + 0 10 9\n10\n' > minus-source.chain
printf 'chain 1 chr1 100 + 0 10 chr1 100 + 0 10 9\n10\n5\n' > orphan.chain
awk 'BEGIN {
  s = "ACGTTGCAGC"; off = 0
  for (c = 1; c <= 5; c++) {
    h = ">chr" c; printf "%s\n", h > "hg38.fa"; off += length(h) + 1
    printf "chr%d\t3000\t%d\t60\t61\n", c, off > "hg38.fa.fai"
    printf "chr%d\t3000\t%d\t60\t61\n", c, off + 1 > "stale.fa.fai"
    line = ""
    for (i = 0; i < 3000; i++) {
      line = line substr(s, (i * 7 + c * 3) % 10 + 1, 1)
      if (length(line) == 60) {
        if (int(i / 60) % 3 == 1) line = tolower(line)
        printf "%s\n", line > "hg38.fa"; off += 61; line = ""
      }
    }
  }
}'
cp hg38.fa stale.fa
awk 'BEGIN {
  OFS = " "
  split("A C A G A T C G C T G T", P, " ")
  n = 0
  pos[1] = "999 1000 1001 1100 1200 1201 1250 1300 1301 1450 1600 1601"
  pos[2] = "500 501 600 750 900 901"
  pos[3] = "100 101 250 300 301 350 400 401 450 500 501"
  pos[4] = "1 250 500"
  pos[5] = "101 500 1100 1101"
  for (c = 1; c <= 5; c++) {
    k = split(pos[c], Q, " ")
    for (j = 1; j <= k; j++) {
      for (p = 1; p <= 6; p++) {
        n++
        a1 = P[2 * p - 1]; a2 = P[2 * p]
        if (n % 5 == 0) { t = a1; a1 = a2; a2 = t }
        if (n % 23 == 0) { a1 = tolower(a1) }
        q = (n % 29 == 0) ? Q[j] ".0" : Q[j]
        printf "  rs%d\t%d   0.0123 %s %s %s\n", n, c, q, a1, a2 > "aadr.snp"
      }
    }
  }
  print "rsA", 1, "0", 1150, "A", "A" > "aadr.snp"
  print "rsN", 1, "0", 1150, "N", "G" > "aadr.snp"
  print "rsX", 23, "0", 1150, "A", "G" > "aadr.snp"
  print "rsY", "X", "0", 1150, "A", "G" > "aadr.snp"
  print "rs0", 0, "0", 1150, "A", "G" > "aadr.snp"
  print "short", 1, "0", 1150, "A" > "aadr.snp"
  print "rsDup", 1, "0", 1150, "A", "G" > "aadr.snp"
  print "rsDup0", 1, "0", 1150, "A", "G" > "aadr.snp"
}'

# --- person-geno --------------------------------------------------------------------
awk 'BEGIN {
  OFS = "\t"
  split("A G C T A C G T", P, " ")
  print "snp", "chrom", "pos", "allele1", "allele2" > "sites.tsv"
  for (i = 1; i <= 200; i++) {
    pos = (i == 51) ? 500 : 10 * i
    print "rs" i, "chr1", pos, P[i % 4 * 2 + 1], P[i % 4 * 2 + 2] > "sites.tsv"
  }
  print "rs999", "chr1" > "sites.tsv"
  split("A/G A|A G/G ./. C/T A a/g G|A A/G/G", T, " ")
  for (i = 1; i <= 220; i += 3) {
    a1 = P[i % 4 * 2 + 1]; a2 = P[i % 4 * 2 + 2]
    g = T[i % 10 + 1]
    gsub(/A/, "1", g); gsub(/G/, "2", g); gsub(/a/, "3", g); gsub(/g/, "4", g)
    gsub(/1/, a1, g); gsub(/2/, a2, g); gsub(/3/, tolower(a1), g); gsub(/4/, tolower(a2), g)
    print "chr1", 10 * i, a1, a2, g > "genotypes.tsv"
  }
  print "chr1", 500, "A", "G", "A/G" > "genotypes.tsv"
  print "chr1", 510, "C", "T", "T/T" > "genotypes.tsv"
  print "chr1", 510, "C", "T", "C/T" > "genotypes.tsv"
  print "1", 20, "A", "G", "A/G" > "genotypes.tsv"
  print "chr1", 30 > "genotypes.tsv"
  print "chr2", 10, "A", "G", "A/G" > "nowhere.tsv"
  for (i = 1; i <= 230; i++) {
    if (i == 100) print "" > "panel.snp"
    print "rs" i, 1, "0.0", 10 * i, "A", "G" > "panel.snp"
    g = sprintf("%d%d%d%d%d", i % 3, (i + 1) % 3, 9, i % 2, 2)
    print g > "panel.geno"
    if (i < 230) print g > "short.geno"
    print g > "long.geno"
    print (i == 77 ? g "0" : g) > "wide.geno"
  }
  print "" > "long.geno"
}'

# --- align and report ---------------------------------------------------------------
awk 'BEGIN {
  split("French Sardinian CEU Basque GBR", M, " ")
  printf "     #eigvals:     5.512     3.201     2.100     1.400     0.900    -0.100\n" > "panel.evec"
  printf "#eigvals: 5.512 3.201 2.100 1.400 0.900 -0.100\n" > "run.evec"
  printf "#eigvals: 5.512 3.201 2.100 1.400 0.900 -0.100\n" > "bad.evec"
  printf "#eigvals: 5.512 3.201\n" > "few.evec"
  for (i = 1; i <= 40; i++) {
    p1 = sin(i * 1.7) * 0.05; p2 = cos(i * 0.9) * 0.04; p3 = sin(i * 0.31) * 0.03; p4 = cos(i * 2.3) * 0.02; p5 = sin(i) * 0.01
    lab = M[i % 5 + 1]
    printf "%12s %10.6f %10.6f %10.6f %10.6f %10.6f %15s\n", "M" i, p1, p2, p3, p4, p5, lab > "panel.evec"
    printf "%s %.6f %.6f %.6f %.6f %.6f %s\n", "M" i, p1, -p2, p3, -p4, p5, lab > "run.evec"
    printf "%s %.6f %.6f %.6f %.6f %.6f %s\n", "M" i, (i % 2 ? p1 : -p1), p2, p3, p4, p5, lab > "bad.evec"
    if (i <= 9) printf "%s %.6f %.6f %s\n", "M" i, p1, p2, lab > "few.evec"
  }
  for (i = 1; i <= 60; i++) {
    printf "%s %.6f %.6f %.6f %.6f %.6f %s\n", "A" i, sin(i * 0.77) * 0.06, cos(i * 1.31) * 0.05, sin(i * 0.5) * 0.02, cos(i) * 0.02, 0.001, "Ancient" > "panel.evec"
  }
  printf "%s %.6f %.6f %.6f %.6f %.6f %s\n", "M3", 0.01, 0.02, 0.03, 0.04, 0.05, "Other" > "panel.evec"
  printf "%s %.6f %.6f %.6f %.6f %.6f %s\n", "S1", 0.0123, -0.0045, 0.002, 0.001, 0.0, "Person" > "run.evec"
  printf "%s %.6f %.6f %.6f %.6f %.6f %s\n", "S1", 0.0123, -0.0045, 0.002, 0.001, 0.0, "Person" > "bad.evec"
  printf "%s %.6f %.6f %s\n", "S1", 0.0123, -0.0045, "Person" > "few.evec"
  printf "%s %.6f %.6f %.6f %.6f %.6f %s\n", "S2", 0.0123, -0.0045, 0.002, 0.001, 0.0, "Ancient" > "run.evec"
  OFS = "\t"
  print "id", "group", "date_bp", "date_sd", "full_date", "locality", "country", "lat", "lon", "snps" > "ancients.tsv"
  split("Germany_CordedWare England_BellBeaker Turkey_N <Mixed&Group> Italy_Imperial_Rome Hungary_LBA_with_a_very_long_name_that_runs_past_45 Türkiye_Neolithic", G, " ")
  split("1650 1950 1949 2450.5 4250 7300 9960 6000 3000 .", D, " ")
  for (i = 1; i <= 60; i++) {
    loc = (i % 3 == 0) ? "." : "Site " (i % 4)
    print "A" i, G[i % 7 + 1], D[i % 10 + 1], "50", ".", loc, "Country" (i % 2), "45", "10", 200000 > "ancients.tsv"
  }
  print "A99", "Nowhere", "5000", "50", ".", ".", ".", "45", "10", 200000 > "ancients.tsv"
  print "A7", "Germany_CordedWare", "5100", "50", ".", "Halle", "Germany" > "ancients.tsv"
  print "" > "ancients.tsv"
  printf "panel_sites\t230\nsites_used\t61\nother\n" > "stats.tsv"
}'
cp run.evec flipped.evec
