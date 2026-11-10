"""What pysam's fetch returns, in the line format of `og bam-fetch-probe`.

    python3 -I pysam_fetch.py BAM REGIONS > fetch.pysam.tsv

REGIONS has one `contig<TAB>start<TAB>stop` per line, 0-based and half-open
as pysam takes them; `#` lines are skipped. One line per read, in fetch
order: the region, then the read's name, flag, reference_start,
reference_end, mapping_quality, cigarstring and get_blocks().
"""

import sys

import pysam


def main(bam_path, regions_path):
    out = sys.stdout
    with pysam.AlignmentFile(bam_path) as bam, open(regions_path) as regions:
        for line in regions:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            contig, start, stop = line.split("\t")
            region = f"{contig}:{start}-{stop}"
            for r in bam.fetch(contig, int(start), int(stop)):
                blocks = ",".join(f"{s}-{e}" for s, e in r.get_blocks())
                out.write(
                    f"{region}\t{r.query_name}\t{r.flag}\t{r.reference_start}\t{r.reference_end}"
                    f"\t{r.mapping_quality}\t{r.cigarstring}\t{blocks}\n"
                )


if __name__ == "__main__":
    main(*sys.argv[1:])
