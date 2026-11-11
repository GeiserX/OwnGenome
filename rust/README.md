# rust/

OwnGenome's own tools, ported from `bin/*.py` one at a time into one static binary, `og`. Each port must write what the Python wrote, byte for byte. The full guide is [docs/rust.md](../docs/rust.md).

| Path | What it is |
|---|---|
| `crates/pycompat` | Python-exact `repr`, format types, `round`, `float()`/`int()`, `math.fsum`, `sum`, `json.dumps`/`loads`, `csv`, `str` methods, truthiness and `str()`/`repr()` of JSON values, float mtimes and `glob.glob` order |
| `crates/og-io` | plain, gzip and BGZF input; Python's universal newlines; atomic writes |
| `crates/og-vcf` | VCF records read in place: INFO kept in order, GT and VEP's CSQ fields, the header passed through |
| `crates/og-bam` | BAM reading as pysam does it (fetch, `reference_end`, `cigarstring`, `get_blocks`), on htslib linked statically |
| `crates/og-pedigree` | the pedigree rules of step 38: families files, members, kinship, declared relationships, PARs |
| `crates/og-chip` | step 49, clonal haematopoiesis: `og chip prepare`, `og chip report`, `og chip check-table` |
| `crates/og-filters` | the VCF filters of steps 04, 05 and 07: `og sv svlen-fix`, `og sv inv-filter`, `og pgx refblock-fill`, `og pgx mito-merge` |
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
