# rust/

OwnGenome's own tools, ported from `bin/*.py` one at a time into one static binary, `og`. Each port must write what the Python wrote, byte for byte. The full guide is [docs/rust.md](../docs/rust.md).

| Path | What it is |
|---|---|
| `crates/pycompat` | Python-exact `repr`, format types, `round`, `float()`/`int()`, `math.fsum`, `sum`, `json.dumps`/`loads`, `csv`, `str` methods, truthiness and `str()`/`repr()` of JSON values, float mtimes and `glob.glob` order |
| `crates/og-io` | plain, gzip and BGZF input; Python's universal newlines; atomic writes |
| `crates/og-vcf` | VCF records read in place: INFO kept in order, GT and VEP's CSQ fields, the header passed through |
| `crates/og-bam` | BAM reading as pysam does it (fetch, `reference_end`, `cigarstring`, `get_blocks`), on htslib linked statically, and the tools that read BAMs: the reads at a variant (`reads`: one CIGAR walk per read, SA pieces, filters per tool), `og cdna check|scan`, `og comphet phase`, `og bam contigs`, the CYP2D6 depth check and the RFC1 motif count |
| `crates/og-pedigree` | the pedigree rules of step 38: families files, members, kinship, declared relationships, PARs |
| `crates/og-qc` | `og qc sex-loss`: step 33's mosaic loss of Y or X, with CPython's medians, `fmean` and exact `variance` |
| `crates/og-chip` | step 49, clonal haematopoiesis: `og chip prepare`, `og chip report`, `og chip check-table` |
| `crates/og-cn` | copy number from depth, steps 32, 35, 41 and 46: `og cn bed`, `og cn call`, `og cn intervals`, `og cn check-bam`, `og cn marker-regions`, `og cn marker-read` |
| `crates/og-reports` | the opt-in step reports: `og y str-bed`, `og y str-report`, `og viral report`, `og exomiser prepare`, `og exomiser report`, `og mei report`, `og mei annotsv-input`, `og mocha report`, `og vntyper report`, `og locityper report` |
| `crates/og-filters` | the VCF filters of steps 04, 05 and 07: `og sv svlen-fix`, `og sv inv-filter`, `og pgx refblock-fill`, `og pgx mito-merge` |
| `crates/og-setup` | what `scripts/setup.sh` runs once: `og str catalog` (the STRchive catalog of steps 09 and 9b), `og setup unzip` (Exomiser's zips), `og setup genetic-map` (the genetic map of steps 26b and 51) |
| `crates/og-ancestry` | step 26b, fine-grained ancestry: `og ancestry-fine panel-map`, `check-map`, `order`, `sample-bed`, `chroms`, `report` |
| `og` | the binary: `og <tool> [args]`, `og version`, `og --tree-hash` |
| `og.sh` | build or download `bin/.og/og-linux-{x86_64,aarch64}`, and check they match this `rust/` |
| `rust-toolchain.toml` | the Rust release, the one `OG_BUILDER` (`versions.env`) ships |

```bash
cargo test                                   # unit tests (needs the pinned toolchain)
cargo clippy --all-targets -- -D warnings
../scripts/setup.sh --og                     # the binaries the pipeline runs, built in Docker
cargo build --release && ../tests/rust/pycompat-diff.sh target/release/og   # pycompat against CPython
```

A binary from a plain `cargo build` is stamped `unknown`, and `rust/og.sh check` refuses it: only `scripts/setup.sh --og` stamps the tree hash the pipeline checks.
