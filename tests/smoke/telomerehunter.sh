# shellcheck shell=sh
# TELOMEREHUNTER_IMAGE row: TelomereHunter2 on telo.bam (the slice plus 1,100
# planted telomeric reads, input telobam) with the GRCh38 bands, as step 10
# runs it. TelomereHunter2 writes one BAM per contig and merges them at once:
# the slice's 2,580 contigs need more than the default 1,024 open files, so
# the soft limit goes up to the hard one, as the Nextflow module does.
set -e
# shellcheck disable=SC3045  # the image's sh is bash, which has ulimit -n and -H
ulimit -n "$(ulimit -Hn)"
telomerehunter2 -ibt /in/telo.bam -o th -p HG002 -c 4 --plotNone -b /in/cytoBand.hg38.txt
cat th/HG002/HG002_summary.tsv
