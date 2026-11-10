# GENCODE 50 slice

A small public slice of GENCODE release 50's basic annotation, laid out as
setup installs it under `GENOME_DIR`, so this folder works as `--root` for
`og map genes`:

- `reference/gencode.v50.basic.genes.gtf`: 346 gene lines.
- `reference/gencode.v50.basic.exons.bed`: 3,489 exon lines (chrom, 0-based
  start, end), sorted and unique.

The Rust tests of `rust/crates/og-map` (`tests/genes.rs`) read it, and the
demo can draw genes from it.

## Source

https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz
(md5 `11e77cf1d4b78b43f57d608e3f7f2202`, as listed in that folder's `MD5SUMS`;
downloaded 2026-11-10). GENCODE's data is free to use with no restriction.

## How it was cut

1. The two full files, written the way `scripts/lib/common.sh`
   (`_install_gencode`, `_gtf_exons`) writes them:

   ```bash
   gzip -cd gencode.v50.basic.annotation.gtf.gz | awk -F'\t' '$3 == "gene"' > gencode.v50.basic.genes.gtf
   gzip -cd gencode.v50.basic.annotation.gtf.gz \
     | awk -F'\t' -v OFS='\t' '$3 == "exon" {print $1, $4 - 1, $5}' \
     | LC_ALL=C sort -k1,1 -k2,2n -k3,3n -u > gencode.v50.basic.exons.bed
   ```

2. The gene lines that overlap a window of 100 kb on each side of these
   genes: GBA1, UGT1A1, DNMT3A, CYP2B6, APOE, CYP2C19, CYP2C9, CYP3A5, CFTR,
   TPMT, HLA-A, SLCO1B1, NUDT15, BRCA2, BRCA1, CYP2D6, HBB, DMD, SRY and
   MT-RNR1. The chrM genes stay in on purpose: the map skips them.

3. The exon lines that overlap the span of a kept gene, in their original
   order.
