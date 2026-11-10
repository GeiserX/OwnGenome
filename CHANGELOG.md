# Changelog

## Unreleased

- Step 33's estimate of mosaic loss of chrY or of one chrX now runs in the Rust `og` binary: `og qc sex-loss` replaces `bin/mosaic_sex_loss.py`, which is gone. The `mosaic_*` rows it adds to `qc/<sample>_sample_qc.tsv` are the same, byte for byte, down to the last digit of CPython's exact variance, and it reads a whole genome's mosdepth windows about three times faster. Step 33 and the Nextflow `SAMPLE_QC` process now need `og`: run `scripts/setup.sh --og` once. Three things change on purpose:
  - A depth of `nan` or `inf` in the mosdepth windows is now named in `mosaic_reason` as unreadable. Before, the median of a block with `nan` depended on where it sat.
  - A VCF site whose `AD` adds up past 4,294,967,295 reads is now skipped. Before, it stopped the tool with an `OverflowError`, and the table got no `mosaic_*` rows.
  - A regions file or VCF that is not UTF-8 is named in `mosaic_reason` with its line number, not Python's byte offset into its read buffer.

- Genome map, first pieces, in the Rust `og` binary. `og map build` writes one person's map document: the chromosome bands, the runs of homozygosity, and every finding with a place on the genome (its own position, or its gene's span). A result the page keeps closed (a finding to act on, and APOE on every page) is tagged with its id on every row that could reveal it, so the page's own closed-result rule hides those rows too. `og map genes` builds the gene and exon files once per GENCODE release, offline. No pipeline output, option or page changes yet; the live page starts drawing the map in a later change. The format is in `rust/crates/og-map/README.md`.
