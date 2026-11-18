# shellcheck shell=sh
# HIPSTR_IMAGE row: what step 45 runs, on the fixture's GIAB slice BAM. Eleven
# Y-STR regions of HipSTR's annotated marker set (the source of
# assets/ystr_panel.tsv) lie in the fixture's chrY slices; three of them
# (DYS570, DYS388, DYS442) are FamilyTreeDNA markers. HipSTR runs with the
# step's flags, and awk reads the bgzipped VCF that og y str-report turns
# into step 45's table. The image's python3 --version is kept for the row's
# check.
set -e
printf 'chrY\t%s\t%s\t%s\t%s\t%s\n' \
  2954700 2954743 4 11 DYS453   2974166 2974201 4 9 DYS579    6993190 6993257 4 17 DYS570 \
  12635604 12635639 3 12 DYS388 12649170 12649229 4 15 DYS442 14611544 14611591 4 12 DYS558 \
  14635532 14635551 5 4 DYS571.1 14635563 14635602 4 10 DYS571.2 19647924 19647959 3 12 DYS489 \
  20607709 20607741 3 11 DYS470.1 20607773 20607784 3 4 DYS470.2 > regions.bed
HipSTR --bams /in/HG002_slice.bam --fasta /in/ref.fa --regions regions.bed --str-vcf calls.vcf.gz \
  --bam-samps HG002 --bam-libs HG002 --haploid-chrs chrY --def-stutter-model --min-reads 5 \
  --max-str-len 150 --output-filters --log hipstr.log
python3 --version > python_version.txt 2>&1
# called.tsv: id, GT, GB, DP and Q of every call with a genotype and those
# four fields and DSTUTTER.
gzip -dc calls.vcf.gz | awk -F'\t' -v OFS='\t' '
  /^#/ { next }
  {
    split("", f); k = split($9, key, ":"); m = split($10, val, ":")
    for (i = 1; i <= k && i <= m; i++) f[key[i]] = val[i]
    if (("GT" in f) && f["GT"] != "." && ("GB" in f) && ("DP" in f) && ("Q" in f) && ("DSTUTTER" in f))
      print $3, f["GT"], f["GB"], f["DP"], f["Q"]
  }' > called.tsv
cat called.tsv
