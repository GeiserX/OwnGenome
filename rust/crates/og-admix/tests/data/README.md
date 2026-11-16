# og-admix test data

- `panel.vcf.gz`: a bgzipped phased panel of five people on chr22 with the
  records the panel reader must skip (a multi-allelic site, an indel, a
  site with no ALT, a `*` allele) and the ones it must keep (two SNVs at
  one position, extra FORMAT fields, a missing or unphased genotype in a
  person nobody asked for, a failed FILTER).
- `*.bcftools.tsv`: what bcftools 1.25
  (`quay.io/biocontainers/bcftools:1.25--h118bc1c_1`) prints for the
  biallelic SNVs of the named people, in the order named:

  ```
  bcftools query -s P4,P1,P3 -i 'TYPE="snp" && N_ALT=1' \
    -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' panel.vcf.gz
  bcftools query -s mum,dad -i 'TYPE="snp" && N_ALT=1' \
    -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' tests/parity/inputs/family/phased.vcf
  ```

  `tests/panel_bcftools.rs` reads the same files with the panel reader and
  compares its haplotypes with these lines.
