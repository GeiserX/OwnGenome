# Changelog

Each release's notes. From 1.0.0 on, a major or minor release has hand-written notes under `docs/releases/`, and the release workflow uses them as the release body. Older releases list their commits on GitHub.

## Unreleased

- Genome map, first pieces, in the Rust `og` binary. `og map build` writes one person's map document: the chromosome bands, the runs of homozygosity, and every finding with a place on the genome (its own position, or its gene's span). A result the page keeps closed (a finding to act on, and APOE on every page) is tagged with its id on every row that could reveal it, so the page's own closed-result rule hides those rows too. `og map genes` builds the gene and exon files once per GENCODE release, offline. No pipeline output, option or page changes yet; the live page starts drawing the map in a later change. The format is in `rust/crates/og-map/README.md`.

## Releases

- 1.0.0 (not released yet): [docs/releases/v1.0.0.md](docs/releases/v1.0.0.md)
- 0.2.0 to 0.8.2: https://github.com/GeiserX/OwnGenome/releases
