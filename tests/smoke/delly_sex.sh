# shellcheck shell=sh
# DELLY_IMAGE row for the --sex option step 19 and the Nextflow DELLY pass.
# xmini.bam is mini.bam moved onto a chrX of GRCh38's length, past PAR1, so
# the planted deletion has no reads on a male chrX. Delly applies --sex only
# there, and ignores it without a word on a chrX of any other length: the
# genotype must be haploid (1) for a male and 1/1 for a female, and the log
# must name the sex. A value Delly does not know must fail the call.
set -e
for sex in male female; do
  delly sr --sex "$sex" -g /in/xmini.fa -o "${sex}.bcf" /in/xmini.bam > "${sex}.log" 2>&1 || { cat "${sex}.log"; exit 1; }
  grep -i 'sex=' "${sex}.log" || true
done
if delly sr --sex bogus -g /in/xmini.fa -o bogus.bcf /in/xmini.bam > bogus.log 2>&1; then
  echo accepted > bogus.txt
else
  echo rejected > bogus.txt
fi
