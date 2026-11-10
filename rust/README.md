# rust/

OwnGenome's own tools, ported from `bin/*.py` one at a time into one static binary, `og`. Each port must write what the Python wrote, byte for byte. The full guide is [docs/rust.md](../docs/rust.md).

| Path | What it is |
|---|---|
| `crates/pycompat` | Python-exact `repr`, format types, `round`, `float()`/`int()`, `math.fsum`, `sum`, `json.dumps`/`loads` |
| `crates/og-io` | plain, gzip and BGZF input; Python's universal newlines; atomic writes |
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
