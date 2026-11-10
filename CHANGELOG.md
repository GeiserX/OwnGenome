# Changelog

## Unreleased

- Genome map, first pieces, in the Rust `og` binary. `og map build` writes one person's map document: the chromosome bands, the runs of homozygosity, and every finding with a place on the genome (its own position, or its gene's span). A result the page keeps closed (a finding to act on, and APOE on every page) is tagged with its id on every row that could reveal it, so the page's own closed-result rule hides those rows too. `og map genes` builds the gene and exon files once per GENCODE release, offline. No pipeline output, option or page changes yet; the live page starts drawing the map in a later change. The format is in `rust/crates/og-map/README.md`.
