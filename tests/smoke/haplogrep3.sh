# shellcheck shell=sh
# HAPLOGREP3_IMAGE row. Start in an empty writable task directory. The pinned
# image finds its bundled config and trees beside the jar, without a network
# or a change into the binary directory (which can be read-only).
set -eu
WORK=$(mktemp -d /out/haplogrep3.XXXXXX)
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
echo "running in $(pwd)"
haplogrep3 classify --tree phylotree-fu-rcrs@1.2 --hetLevel 0.9 \
  --input /in/mito.vcf.gz --output /out/haplogroup.txt --extend-report
cat /out/haplogroup.txt

# Mutect2-style GT:AF calls: below, on and above the chosen threshold.
# A homozygous ALT call stays in the profile even with a low AF.
cat > fractions.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=chrM,length=16569>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chrM	750	.	A	G	50	PASS	.	GT:AF	1/1:0.05
chrM	3243	.	A	G	50	PASS	.	GT:AF	0/1:0.05
chrM	5000	.	C	T	50	PASS	.	GT:AF	0/1:0.95
chrM	9000	.	T	C	50	PASS	.	GT:AF	0/1:0.60
chrM	10000	.	A	G	50	PASS	.	GT:AF	0/1:0.90
VCF
haplogrep3 classify --tree phylotree-fu-rcrs@1.2 --hetLevel 0.9 \
  --input fractions.vcf --output fractions.txt --extend-report
awk -F '\t' '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      gsub(/"/, "", $i)
      if ($i == "Input_Sample") col = i
    }
    next
  }
  NR == 2 && col {
    profile = " " $col " "
    gsub(/"/, "", profile)
    ok = profile ~ / 750G / && profile ~ / 5000T / && profile ~ / 10000G / \
      && profile !~ / 3243G / && profile !~ / 9000C /
  }
  END {
    if (!ok) {
      print "FAIL: haplogrep3 profile does not obey hetLevel 0.9" > "/dev/stderr"
      exit 1
    }
    print "PASS: AF 0.90 and 0.95 retained; AF 0.05 and 0.60 excluded; homozygous ALT retained"
  }
' fractions.txt
