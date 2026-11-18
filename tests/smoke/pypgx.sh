# shellcheck shell=sh
# PYPGX_IMAGE row: the CYP2D6 path of step 32 (depth of coverage, VDR control
# statistics, run-ngs-pipeline) on the GIAB slice, with the pypgx-bundle tag
# PYPGX_BUNDLE_VERSION names. pypgx makes the input VCF itself.
set -e
export PYPGX_BUNDLE=/in/pypgx-bundle
ln -sfn /in/pypgx-bundle "${HOME}/pypgx-bundle"
BAM=/in/HG002_slice.bam
pypgx create-input-vcf input.vcf.gz /in/ref.fa "$BAM" --assembly GRCh38 --genes CYP2D6
pypgx prepare-depth-of-coverage depth.zip "$BAM" --assembly GRCh38 --genes CYP2D6
pypgx compute-control-statistics VDR control.zip "$BAM" --assembly GRCh38
pypgx run-ngs-pipeline CYP2D6 CYP2D6 --variants input.vcf.gz --depth-of-coverage depth.zip \
  --control-statistics control.zip --assembly GRCh38
# The Genotype of the first row of results.zip's data.tsv.
DATA=$(unzip -l CYP2D6/results.zip | awk '$NF ~ /data\.tsv$/ {print $NF; exit}')
[ -n "$DATA" ] || { echo "CYP2D6/results.zip holds no data.tsv"; exit 1; }
unzip -p CYP2D6/results.zip "$DATA" | awk -F'\t' '
  NR == 1 { for (i = 1; i <= NF; i++) if ($i == "Genotype") c = i; next }
  c { print $c; found = 1; exit }
  END { if (!found) { print "no Genotype row in data.tsv" > "/dev/stderr"; exit 1 } }' > CYP2D6.genotype.txt
cat CYP2D6.genotype.txt
