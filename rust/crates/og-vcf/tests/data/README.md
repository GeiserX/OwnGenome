# Test slices

Four small slices of public HG002 files (Genome in a Bottle's Ashkenazi son), so the round-trip tests read real VCF layouts: many INFO items per record, a VEP `CSQ` annotation, long alleles and gVCF reference blocks. Each slice is the file's whole header and its first records, compressed with `bgzip` 1.24 (`staphb/bcftools:1.24`). No private genome is here. The repository's `.gitignore` ignores `*.vcf.gz`, so a new slice is added with `git add -f`.

| File | Records | Taken from |
|---|---|---|
| `hg002_benchmark.vcf.gz` | first 300 | GIAB v4.2.1 benchmark, `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` (https://giab.s3.amazonaws.com/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/), as the e2e fixture release `fixture-v7` carries it |
| `hg002_smvar_trf.vcf.gz` | first 100 with `TRFdiff` in INFO | GIAB v5.0q, `HG002_GRCh38_v5.0q_smvar.vcf.gz` (https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/), as `fixture-v7` carries it |
| `hg002_vep.vcf.gz` | first 12 | `HG002_vep.vcf` of `fixture-v7`: GIAB truth variants annotated by VEP 116 (`scripts/ci/build-fixture.sh`) |
| `hg002_chr20.g.vcf.gz` | first 400 | DeepVariant 1.10.0 case-study gVCF, https://storage.googleapis.com/deepvariant/case-study-outputs/1.10.0/deeptrio/wgs-chr20/HG002.child.g.vcf.gz |

How a slice was made (the benchmark one; the others differ in the file and the record filter):

```bash
{ gzip -dc HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | grep '^#'
  gzip -dc HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | grep -v '^#' | head -300; } > hg002_benchmark.vcf
bgzip hg002_benchmark.vcf
```
