# Changelog

Each release's notes. From 1.0.0 on, a major or minor release has hand-written notes under `docs/releases/`, and the release workflow uses them as the release body. Older releases list their commits on GitHub.

## Unreleased

- The CYP2D6 depth check and the RFC1 motif count, ported to the Rust `og` binary: `og pgx cyp2d6-bed`, `og pgx cyp2d6-measure` and `og pgx cyp2d6-check` (from `bin/cyp2d6_depth_check.py`) and `og str motif` (from `bin/rfc1_motif.py`) read BAMs with htslib, without pysam. On 67 cases they write what the Python writes, byte for byte, with one fix: a mosdepth regions file that lists one region twice with two depths is now refused (exit 2); the Python averaged the two. Steps 21, 32 and 9b still run the Python for now; no option, output file or format changes.
- Genome map, first pieces, in the Rust `og` binary. `og map build` writes one person's map document: the chromosome bands, the runs of homozygosity, and every finding with a place on the genome (its own position, or its gene's span). A result the page keeps closed (a finding to act on, and APOE on every page) is tagged with its id on every row that could reveal it, so the page's own closed-result rule hides those rows too. `og map genes` builds the gene and exon files once per GENCODE release, offline. No pipeline output, option or page changes yet; the live page starts drawing the map in a later change. The format is in `rust/crates/og-map/README.md`.
- Rust building blocks for the ports that look things up by position. `og-io` reads the lines of a region from a bgzipped file through its `.tbi` or `.csi` index, as `tabix` prints them, reading only the blocks the index points at. `og-intervals` reads bases from an uncompressed FASTA through its `.fai`, as `samtools faidx` prints them, and refuses a compressed FASTA or a `.fai` that does not match the file. It also adds sorted interval sets (union, intersection, lookups by binary search, BED files such as the PAR regions) and reciprocal overlap. Nothing calls them yet, so no pipeline output, option or file changes.

## Releases

- 1.0.0 (not released yet): [docs/releases/v1.0.0.md](docs/releases/v1.0.0.md)
- 0.2.0 to 0.8.2: https://github.com/GeiserX/OwnGenome/releases
