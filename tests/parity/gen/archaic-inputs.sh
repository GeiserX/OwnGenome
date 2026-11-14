#!/usr/bin/env bash
# The inputs of the archaic cases of step 43 (archaic-table, archaic-strict,
# archaic-genotypes, archaic-report), all synthetic, written with awk:
#   modern.tsv, altai.tsv, vindija.tsv, denisova.tsv
#             bcftools query lines of one GRCh37 chromosome (22): 1000 Genomes
#             SNVs with their six frequencies (indels, a split multi-allelic
#             site, frequencies that need :.6g's exponent form, NaN, out of
#             range and unreadable values) and the three genomes' calls
#             (missing, partial, multi-allelic, gVCF <*> and REF-mismatch
#             calls, and sites only an archaic genome carries), also gzipped
#   altai.bed, vindija.bed, excl.bed, empty.bed
#             the masks: chromosome 22 named with and without chr, unsorted,
#             overlapping, reversed and empty intervals, track lines
#   unsorted.tsv  a modern file that goes back
#   rows.tsv, rows-unsorted.tsv, strict.bed
#             lifted GRCh38 rows (repeated positions) and the strict mask
#   chr21.tsv.gz, chr22.tsv.gz, chrX.tsv.gz
#             the GRCh38 table step 43 reads (chrX is not a table name),
#             with masked codes and frequencies that round half to even
#   s.vcf, s.vcf.gz, two.vcf, back.vcf, none.vcf, badcode/
#             the sample's VCF: split records, clashes, indels, a
#             multi-allelic record, RefCall, missing GT, a chromosome the
#             table lacks; a two-sample VCF, one whose chromosome comes
#             back, one with no table chromosome
#   chr21.Altai.txt ... chr22.Denisova.txt, spans.tsv, genes.gtf
#             IBDmix's segments (short, weak, overlapping and touching
#             ones, other IDs) and a GENCODE-style gene GTF
set -euo pipefail

# --- table (GRCh37 inputs) ----------------------------------------------------------
awk 'BEGIN {
  OFS = "\t"
  split("A C G T", B, " ")
  split("0.5 1e-05 0.000123 . nan 1.5 -0.1 0.123456789 1 0 2.5e-7 x", F, " ")
  for (i = 1; i <= 60; i++) {
    pos = 1000 + 10 * i
    r = B[i % 4 + 1]; a = B[(i + 1) % 4 + 1]
    if (i % 11 == 0) { print pos, r "T", a, F[1], F[2], F[3], F[4], F[5], F[6] > "modern.tsv"; continue }
    if (i % 13 == 0) {
      print pos, r, a, "0.1", "0.2", "0.3", "0.4", "0.5", "0.6" > "modern.tsv"
      print pos, r, B[(i + 2) % 4 + 1], "0.1", "0.2", "0.3", "0.4", "0.5", "0.6" > "modern.tsv"
      continue
    }
    if (i % 7 == 3) continue
    print pos, r, a, F[i % 12 + 1], F[(i + 1) % 12 + 1], F[(i + 2) % 12 + 1], F[(i + 3) % 12 + 1], F[(i + 4) % 12 + 1], F[(i + 5) % 12 + 1] > "modern.tsv"
  }
  split("0/0 0/1 1/1 ./. 1/2 0|1 1 0/1/1 . 1|1", G, " ")
  for (i = 1; i <= 64; i++) {
    pos = 1000 + 10 * i + (i % 9 == 0 ? 5 : 0)
    r = B[i % 4 + 1]; a = B[(i + 1) % 4 + 1]
    if (i % 3 != 0) print pos, r, (i % 10 == 4 ? "<*>" : a), G[i % 10 + 1] > "altai.tsv"
    if (i % 4 != 1) print pos, (i % 17 == 0 ? "C" : r), (i % 6 == 0 ? "." : a), G[(i + 3) % 10 + 1] > "vindija.tsv"
    if (i % 5 != 2) print pos, r, (i % 8 == 0 ? B[(i + 2) % 4 + 1] : a), G[(i + 7) % 10 + 1] > "denisova.tsv"
  }
}'
for f in modern altai vindija denisova; do gzip -n -c "${f}.tsv" > "${f}.tsv.gz"; done
printf '1200\tA\tG\t0.1\t0.1\t0.1\t0.1\t0.1\t0.1\n1100\tA\tG\t0.1\t0.1\t0.1\t0.1\t0.1\t0.1\n' > unsorted.tsv
cat > altai.bed <<'EOF'
track name=altai
# a comment
22	1300	1500
22	1000	1100
22	1080	1250
22	1600	1600
22	1700	1650
chr21	0	100000
EOF
printf 'browser position chr22\nchr22\t1000\t1400\nchr22\t1450\t1700\n' > vindija.bed
printf 'chr22\t1195\t1215\n22\t1500\t1520\n' > excl.bed
printf 'chr21\t0\t1000\n' > empty.bed

# --- strict (lifted rows) -----------------------------------------------------------
awk 'BEGIN {
  OFS = "\t"
  split("0 1 2 9 10 11 12", C, " ")
  for (i = 1; i <= 40; i++) {
    pos = 5000 + 7 * i
    row = "chr22" OFS pos OFS "A" OFS "G" OFS C[i % 7 + 1] OFS C[(i + 2) % 7 + 1] OFS C[(i + 4) % 7 + 1] OFS "0.1" OFS "0.2" OFS "0.3" OFS "0.4" OFS "0.5" OFS "0.6"
    print row > "rows.tsv"
    if (i % 9 == 0) print "chr22" OFS pos OFS "C" OFS "T" OFS "1" OFS "1" OFS "1" OFS "0" OFS "0" OFS "0" OFS "0" OFS "0" OFS "0" > "rows.tsv"
  }
}'
{ head -5 rows.tsv; printf 'chr22\t5001\tA\tG\t1\t1\t1\t0\t0\t0\t0\t0\t0\n'; tail -3 rows.tsv; } > rows-unsorted.tsv
printf 'chr22\t5010\t5100\n22\t5150\t5200\nchr22\t5240\t5260\n' > strict.bed

# --- genotypes (the GRCh38 table and the sample) ------------------------------------
HEAD='chrom\tpos\tref\talt\tAltai\tVindija\tDenisova\tAF_ALL\tAF_AFR\tAF_AMR\tAF_EAS\tAF_EUR\tAF_SAS'
table() { # CHROM START N
  awk -v c="$1" -v s="$2" -v n="$3" 'BEGIN {
    OFS = "\t"
    split("A C G T", B, " ")
    split("0 1 2 9 10 11 12", C, " ")
    split("0.0025 0.0075 0.5 1e-05 . nan 0.99 1.2 -0.3 0.0125 0.333333", F, " ")
    for (i = 1; i <= n; i++) {
      pos = s + 10 * i
      print c, pos, B[i % 4 + 1], B[(i + 1) % 4 + 1], C[i % 7 + 1], C[(i + 3) % 7 + 1], C[(i * 5) % 7 + 1], \
        F[i % 11 + 1], F[(i + 1) % 11 + 1], F[(i + 2) % 11 + 1], F[(i + 3) % 11 + 1], F[(i + 4) % 11 + 1], F[(i + 5) % 11 + 1]
    }
  }'
}
{ printf '%b\n' "$HEAD"; table chr21 1000 30; } | gzip -n > chr21.tsv.gz
{ printf '%b\n' "$HEAD"; table chr22 2000 50; } | gzip -n > chr22.tsv.gz
{ printf '%b\n' "$HEAD"; table chrX 3000 5; } | gzip -n > chrX.tsv.gz
mkdir -p badcode emptydir
{ printf '%b\n' "$HEAD"; table chr22 2000 4 | awk -F'\t' -v OFS='\t' 'NR == 3 { $5 = 3 } { print }'; } | gzip -n > badcode/chr22.tsv.gz
vcf_head() {
  printf '##fileformat=VCFv4.2\n##contig=<ID=chr22>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$1"
}
{
  vcf_head S
  awk 'BEGIN {
    OFS = "\t"
    split("A C G T", B, " ")
    split("0/1 1/1 0/0 ./. 0|1 1|1 0/2 1", G, " ")
    print "chr1", 50, ".", "A", "G", ".", "PASS", ".", "GT", "0/1"
    for (i = 1; i <= 30; i += 2) {
      pos = 1000 + 10 * i; r = B[i % 4 + 1]; a = B[(i + 1) % 4 + 1]
      print "chr21", pos, ".", r, a, ".", (i % 5 == 0 ? "RefCall" : "PASS"), ".", "GT:GQ", G[i % 8 + 1] ":30"
    }
    for (i = 1; i <= 50; i++) {
      pos = 2000 + 10 * i; r = B[i % 4 + 1]; a = B[(i + 1) % 4 + 1]; o = B[(i + 2) % 4 + 1]
      k = i % 12
      if (k == 0) print "chr22", pos, ".", r, a, ".", ".", ".", "GT", G[i % 8 + 1]
      else if (k == 1) { print "chr22", pos, ".", r, a, ".", "PASS", ".", "GT", "0/1"; print "chr22", pos, ".", r, o, ".", "PASS", ".", "GT", "0/0" }
      else if (k == 2) { print "chr22", pos, ".", r, a, ".", "PASS", ".", "GT", "0/0"; print "chr22", pos, ".", r, o, ".", "PASS", ".", "GT", "0/0" }
      else if (k == 3) { print "chr22", pos, ".", r, a, ".", "PASS", ".", "GT", "0/1"; print "chr22", pos, ".", r, o, ".", "PASS", ".", "GT", "0/1"; print "chr22", pos, ".", r, a, ".", "PASS", ".", "GT", "1/1" }
      else if (k == 4) { print "chr22", pos, ".", r, a "," o, ".", "PASS", ".", "GT", "1/2" }
      else if (k == 5) { print "chr22", pos, ".", r "T", r, ".", "PASS", ".", "GT", "0/1"; print "chr22", pos, ".", r, a, ".", "PASS", ".", "GT", "1/1" }
      else if (k == 6) { print "chr22", pos, ".", r, o, ".", "PASS", ".", "GT", "1/1" }
      else if (k == 7) { print "chr22", pos, ".", r, a, ".", "PASS", ".", "GQ", "20" }
      else if (k == 8) { print "chr22", pos, ".", r, a, ".", "PASS", ".", "GT:DP", "./.:3" }
      else if (k == 9) { print "chr22", pos, ".", r, a, ".", "PASS", ".", "GT", "0/1"; print "chr22", pos, ".", r, a, ".", "PASS", ".", "GT", "./." }
      else if (k == 10) { print "chr22", pos, ".", r, ".", ".", "PASS", ".", "GT", "0/0" }
      else print "chr22", pos, ".", r, a, ".", "PASS", ".", "GT:AD", "1/1:0,9"
    }
    print "chr22", 9999, ".", "A", "G", ".", "PASS", ".", "GT", "0/1"
  }'
} > s.vcf
gzip -n -c s.vcf > s.vcf.gz
{ printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n'; printf 'chr22\t2010\t.\tC\tG\t.\tPASS\t.\tGT\t0/1\t0/0\n'; } > two.vcf
{ vcf_head S; printf 'chr22\t2010\t.\tC\tG\t.\tPASS\t.\tGT\t0/1\nchr21\t1010\t.\tC\tG\t.\tPASS\t.\tGT\t0/1\nchr22\t2020\t.\tG\tT\t.\tPASS\t.\tGT\t0/1\n'; } > back.vcf
{ vcf_head S; printf 'chr1\t2010\t.\tC\tG\t.\tPASS\t.\tGT\t0/1\nchrX\t3010\t.\tC\tG\t.\tPASS\t.\tGT\t0/1\n'; } > none.vcf

# --- report (IBDmix output) ---------------------------------------------------------
printf 'chrom\tfirst\tlast\tsites\tsites_with_your_alt\tsites_other_allele\nchr21\t1000000\t46000000\t3000\t400\t2\nchr22\t16000000\t50800000\t2500\t300\t1\n' > spans.tsv
seg() { printf 'ID\tchrom\tstart\tend\tslod\n'; cat; }
seg > chr21.Altai.txt <<'EOF'
SAMPLE	21	1000000	1200000	12.5
SAMPLE	21	1150000	1400000	7.25
SAMPLE	21	2000000	2040000	30
SAMPLE	21	3000000	3100000	3.99
P001	21	5000000	9000000	99
SAMPLE	21	9000000	9500000	4
EOF
seg > chr21.Vindija.txt <<'EOF'
SAMPLE	21	1400000	1600000	8.5
SAMPLE	21	9000000	9300000	45.5
SAMPLE	chr21	20000000	20060000	4.5
EOF
seg > chr21.Denisova.txt <<'EOF'
SAMPLE	21	1100000	1300000	6
SAMPLE	21	30000000	30500000	14.5
SAMPLE	21	30400000	30700000	5.5
EOF
seg > chr22.Altai.txt <<'EOF'
SAMPLE	22	17000000	17800000	22.5
SAMPLE	22	40000000	40050000	9
EOF
seg > chr22.Vindija.txt <<'EOF'
SAMPLE	22	17500000	18100000	10.5
EOF
seg > chr22.Denisova.txt <<'EOF'
SAMPLE	22	25000000	25900000	8.5
SAMPLE	22	17000000	17050000	15
EOF
cat > genes.gtf <<'EOF'
##description: synthetic
chr21	HAVANA	gene	1000	1000000	.	+	.	gene_id "G1"; gene_type "protein_coding"; gene_name "TOUCHES_START";
chr21	HAVANA	gene	1199999	1250000	.	+	.	gene_id "G2"; gene_type "protein_coding"; gene_name "INSIDE";
chr21	HAVANA	gene	1100000	1300000	.	+	.	gene_id "G3"; gene_type "lncRNA"; gene_name "NOT_CODING";
chr21	HAVANA	transcript	1100000	1300000	.	+	.	gene_id "G2"; gene_type "protein_coding"; gene_name "INSIDE";
chr21	HAVANA	gene	1600000	1700000	.	+	.	gene_id "G4"; gene_type "protein_coding"; gene_name "AT_END";
chr21	HAVANA	gene	30500000	30600000	.	-	.	gene_id "G5"; gene_type "protein_coding"; gene_name "DENI";
chr21	HAVANA	gene	1300000	1350000	.	-	.	gene_id "G6"; gene_type "protein_coding"; gene_name "INSIDE";
chr22	HAVANA	gene	17700000	17900000	.	+	.	gene_id "G7"; gene_type "protein_coding"; gene_name "ON22";
chr22	HAVANA	gene	17000000	17100000	.	+	.	gene_id "G8"; gene_type "protein_coding"; gene_id_only "x";
EOF
