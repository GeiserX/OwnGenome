# The Rust Core

This page is for contributors. OwnGenome's own tools are moving from Python to Rust, one tool per pull request, into a single binary, `og`. This page says how the Rust code is laid out, how the binary reaches the pipeline, how it stays byte-for-byte compatible with the Python it replaces, and what CI checks.

Steps 04, 05 and 07 call `og` (the filters of `og sv` and `og pgx`), steps 23 and 31 and their Nextflow modules call it (`og vcf csq-rank` and `og table constraint-join`), step 27 and `CPIC_LOOKUP` call it (`og pgx cpic-report`), and `setup.sh` calls it to mask the reference. The pieces below are the base every port builds on.

## Layout

```
rust/
  Cargo.toml            the workspace (edition 2024); Cargo.lock pins every crate
  rust-toolchain.toml   the Rust release, the one OG_BUILDER ships
  og.sh                 build, download and check the binary
  crates/pycompat/      Python-exact formatting and parsing
  crates/og-io/         plain, gzip and BGZF input, universal newlines, atomic writes
  crates/og-vcf/        VCF records read in place, INFO kept in order, GT and CSQ fields; csq-rank
  crates/og-bam/        BAM reading as pysam does it, on a statically linked htslib
  crates/og-pgx/        og pgx: PharmCAT's report read into typed gene calls (step 27)
  crates/og-filters/    the VCF filters of steps 04, 05 and 07
  crates/og-steps/      the step helpers: constraint-join
  crates/og-chip/       step 49, clonal haematopoiesis: og chip prepare, report, check-table
  crates/og-reports/    the opt-in step reports of steps 45, 47, 48 and 50 to 53
  og/                   the binary: og <tool> [args]
bin/og                  the shim tasks call; runs bin/.og/og-linux-<arch> (og-darwin-<arch> on a Mac)
tests/rust/             the differential test against CPython, the image test
```

Each ported tool becomes a subcommand named after its file in `bin/` without `.py`, in kebab-case. Today `og` has `og version`, `og --tree-hash`, `og setup mask-fasta`, `og run state|tasks|eta|headroom|delta`, `og pgx cyp2d6-bed|cyp2d6-measure|cyp2d6-check`, `og str motif` and two hidden ones the tests drive: `og pycompat-probe` and `og bam-fetch-probe`. The filters of steps 04, 05 and 07 are grouped by step: `og sv svlen-fix` and `og sv inv-filter` (were `sv_svlen_fix.py` and `manta_inv_filter.py`), `og pgx refblock-fill` and `og pgx mito-merge` (were `pgx_refblock_fill.py` and `pgx_mito_merge.py`). Steps 23 and 31 use `og vcf csq-rank` and `og table constraint-join` (ports of the awk programs `bin/csq_rank.awk` and `bin/constraint_join.awk`, named by what they work on). Step 27 uses `og pgx cpic-report` and `og pgx step36-outside` (were the two commands of `pgx_parse.py`). Step 49 runs `og chip prepare`, `og chip report` and `og chip check-table` (the commands of `chip_report.py`, which keeps only the detection limit and the column names the reports import). The opt-in steps read their tool's output with `og y str-bed|str-report` (45), `og viral report` (47), `og exomiser prepare|report` (48), `og mei report|annotsv-input` (50), `og mocha report` (51), `og vntyper report` (52) and `og locityper report` (53); they replaced `ystr_report.py`, `viral_report.py`, `exomiser_report.py`, `mei_report.py`, `mocha_report.py`, `vntyper_report.py` and `locityper_report.py`.

The rules: no `unsafe` (the workspace forbids it), `cargo clippy -D warnings`, `cargo fmt`, and few, well-known crates (`indexmap`, `flate2`, `zlib-rs`, `memchr`, `smallvec` and `rust-htslib` so far, declared in `rust/Cargo.toml` and locked in `Cargo.lock`).

## Building

```bash
./scripts/setup.sh --og
```

`scripts/setup.sh <genome_dir>` does the same when `og` is missing or stale, after it pulls the images. `--og` builds only `og`. This puts four binaries in `bin/.og/` (ignored by git):

- `og-linux-x86_64` and `og-linux-aarch64`, two static musl builds for the pipeline's containers. Both are kept because an Apple Silicon host runs amd64 images under emulation, where `uname -m` says `x86_64`.
- `og-darwin-arm64` and `og-darwin-x86_64`, for a Mac that hosts the live site (`scripts/owngenome-live.sh`), whose genome map runs `og map` on the host. They are built with `--no-default-features`, which leaves out the cargo feature `bam`: that host never reads a BAM, and BAM reading needs htslib, which is C and would have to cross-compile for Apple. The BAM tools (`og pgx cyp2d6-*`, `og str motif`, `og bam-fetch-probe`) say so and exit 2 there. Every Linux build has the feature, since it is on by default.

- **From a checkout of `main` or a branch,** it builds them from `rust/` inside `OG_BUILDER`, a pinned [cargo-zigbuild](https://github.com/rust-cross/cargo-zigbuild) image (`versions.env`). Docker is the only thing it needs; no Rust on the host. The cargo cache lives in `rust/target/og-build/`.
- **From a checkout of a release tag** (`vX.Y.Z`), it downloads the release's four binaries, checks them against the release's `SHA256SUMS` and, when `gh` is installed and logged in, against the build attestation (`gh attestation verify`). A release without binaries falls back to a build.

For day-to-day work on the Rust code a host toolchain is quicker: `cd rust && cargo test`. A plain `cargo build` stamps the binary `unknown`, which the check below refuses, so use `scripts/setup.sh --og` for a binary the pipeline will run.

The image's Rust release (1.93.0) and `rust/rust-toolchain.toml` must agree, or the build downloads a toolchain first. `OG_BUILDER` has no `_IMAGE` suffix on purpose: no pipeline step runs it, so it has no smoke row and Renovate does not bump it. Bump the image and the toolchain file together, by hand.

## How a task runs it

A Nextflow task has `bin/` on its PATH, so a module calls `og <tool>` by its bare name. A bash step calls `/pgp-bin/og <tool>` inside the image it already uses. Either way it starts `bin/og`, a `#!/bin/sh` shim that runs `bin/.og/og-linux-$(uname -m)`, or exits 127 with `ERROR: og is not built for <arch>; run scripts/setup.sh --og`. On macOS the same shim runs `bin/.og/og-darwin-$(uname -m)` (`ERROR: og is not built for macOS <arch>` when it is missing); any other system gets exit 127.

The binary is static, so it runs in any image that has `/bin/sh`: the Python image, the bcftools image, the Manta and Parascopy images. `tests/rust/og-in-images.sh` checks those four on every platform each one publishes.

## The stale-binary check

`bin/.og` is not in git, so after a pull or an edit it can hold a binary built from other code. Results would then come from code the checkout does not show. To stop that:

- The build stamps the binary with the git tree hash of `rust/`: `git rev-parse HEAD:rust` on a clean checkout, or for an edited `rust/` the tree those edits would commit. Any edit, committed or not, changes the hash.
- `rust/og.sh check` reads the stamp from each binary with `grep`, without running it, and compares it with the checkout's. The macOS pair is optional, because no pipeline step runs it, but a stale one is refused like a Linux one.
- `run-all.sh` (even with `SKIP_VALIDATION=true`), `main.nf` and `validate-setup.sh` run that check. A binary built from another `rust/` stops them with `run scripts/setup.sh --og`.
- A missing binary is still a note, which names steps 04, 05 and 07. Those steps, and steps 23 and 31 with their modules, stop without it: `bin/og` says `ERROR: og is not built for <arch>; run scripts/setup.sh --og` (exit 127). `REQUIRED` in `rust/og.sh` stays `false` until the fake-docker and Nextflow stub tests, which start runs without `og`, build it first; then a missing binary stops a run before any step.
- `run_manifest.tsv` gets a row `run og_tree_hash`, so a result names the Rust code that made it (`not built` before the first build).

## Matching Python byte for byte

A port must write exactly what the Python tool wrote. Most differences hide in number formatting, so the shared primitives live in `pycompat`, each one written to return what CPython 3.15 returns:

| Module | What it covers |
|---|---|
| `fmt` | `repr(float)` (also `str` and `f"{x}"`); the `f`, `e`, `g` and `%` format types at any precision; `{:,}` and `{:,.Nf}`; `{:0Nd}`; `%-Ns` and `%Ns`; `round(x)` and `round(x, n)`, both rounding the exact binary value half to even (`round(0.35, 1) == 0.3`, `round(2.675, 2) == 2.67`) |
| `num` | `float(str)` and `int(str)` with Python's whitespace, signs and underscores; `math.fsum`; the builtin `sum`, which since Python 3.12 is a compensated sum and so differs from a plain loop |
| `json` | `json.dumps` with `indent` (none or a number), `separators`, `ensure_ascii` (surrogate pairs for astral characters), `sort_keys` and `allow_nan`; `json.loads` with `NaN` and `Infinity`, where a duplicate key keeps its first place and takes the last value |
| `csv` | `csv.writer`, `csv.reader` and `csv.DictReader` for the default dialect, with `QUOTE_MINIMAL` or `QUOTE_NONE` and any delimiter; a port of CPython's `_csv.c` state machine |
| `text` | `str.split`, `strip`/`lstrip`/`rstrip`, `isspace`/`isdigit`/`isdecimal`/`isalpha`, `upper`/`lower` with full case mapping, `s[start:stop]` in code points, and `html.escape` |
| `value` | over a value `json.loads` read: `bool(v)`, `str(v)`, `repr(v)`, `v or w` and `d.get(k) or default` |
| `mtime` | `os.path.getmtime` as a float (`sec + nsec * 1e-9`, so times a few hundred nanoseconds apart can be equal), `os.stat().st_mtime_ns`, `glob.glob` in its own order (`**` with `recursive`), `glob.escape`, `fnmatch`, `max(paths, key=os.path.getmtime)` with its tie rule, and `sorted()` of paths; paths are bytes, and a byte that is not UTF-8 compares as its `surrogateescape` character. Unix only |

`og-io` reads a path or `-` as plain text, gzip or BGZF (sniffed by the magic bytes; NUL padding after a gzip member is skipped, as Python's `gzip` does), splits lines as Python's text mode does (`\n`, `\r\n` and a lone `\r`), decodes UTF-8 strictly or with `replace`, and writes a file under a temporary name before renaming it into place.

`og-bam` reads an indexed BAM through `rust-htslib`, which builds htslib from source and links it into `og` on both musl targets, so the binary stays static. Its curl, S3, GCS, bzip2 and lzma features are off: `og` reads local BAM files only, and only CRAM needs bzip2 and lzma. `fetch` is htslib's own region query, the call pysam makes, so a region returns the same reads; `reference_end`, `cigarstring` and `get_blocks` are rebuilt from the CIGAR the way pysam computes them. htslib adds about 700 KB to each binary. The fetch smoke checks it: `crates/og-bam/tests/data` holds a 2 kb slice of the test fixture's BAM at CYP2D6, a list of regions with the edge cases (a read's first and last base, a base inside a deletion, a placed unmapped read, a contig without reads) and what pysam fetches there (`tests/pysam_fetch.py` writes it). `cargo test` and `og bam-fetch-probe` must write the same lines.

Known limits, each an error in Rust where Python would go on: ints beyond 128 bits, a lone surrogate in a JSON string, non-ASCII digits in `float()` and `int()`, a string `indent`, JSON nested deeper than 1000 levels (`RecursionError`; CPython 3.14 reads about 100,000 on a default stack), and `repr()` of a string that holds a character above U+00FF, which also covers `str()` of a list or dict holding one (`value` returns `Unsupported`; `str()` of the string itself works). The first five never appear in what the pipeline reads. The U+00FF limit can: a port that turns a list or dict into text must handle that error, and the PR that ports it says how.

### The differential test

`tests/rust/pycompat-diff.sh <og binary> [N] [SEED]` checks every `pycompat` operation against CPython. `tests/rust/pycompat_oracle.py` runs with `python3 -I` in `PYTHON_IMAGE`, draws N random cases per operation and writes its result for each. `og pycompat-probe` reads the same cases and must write the same lines, byte for byte. The values include the ones that catch a formatting bug: binary ties (0.125, 2.675, 0.35), 1e-05 and 1e16 where `repr` switches layout, -0.0, NaN, subnormals, random bit patterns, numbers exactly halfway between two shortest `repr` candidates, astral characters, CRLF and lone CR (both as lines without their endings and as `readline()` results with them), and broken JSON. A failure prints the first cases that differ and the seed, so `pycompat-diff.sh <og> <N> <seed>` replays it. `PYTHON=python3` runs the oracle on the host for a quick loop.

To add a primitive: write it in `pycompat` with unit tests, add an operation to `og/src/probe.rs` and the same operation to the oracle, then run the differential test.

## CI

- **Guard, job `rust`** (every pull request, required through `ci-ok`): `cargo fmt --check`; `cargo clippy -D warnings`; a planted failing test (`--features planted-failure`) that must fail; `cargo test`; the differential test at 100,000 cases per operation; the four builds through `scripts/setup.sh --og` (both musl ones and the macOS pair); and the stale check on the real binaries, which must refuse them after an edit to `rust/`, in `rust/og.sh check` and in `nextflow run main.nf -stub`, before any task. Then the tests of the `og sv` and `og pgx` filters and the step 05 tests run with `bin/og` and `OG_REQUIRED=1`; in `unit-tests`, which has no `og`, they skip with a note.
- **Guard, job `unit-tests`**: `tests/test_og_check.sh` runs the check on planted trees (no binary, a stale one, an edit that is uncommitted, untracked or committed, a file git ignores, a binary stamped `unknown`, a fresh and a stale macOS pair, a tree hash that cannot be computed, a release download with a bad checksum or attestation, without the macOS pair, with `gh` not logged in, with a network error, without binaries, or whose copy into `bin/.og` fails), and checks that `run-all.sh`, `validate-setup.sh` and `write_manifest.sh` act on it.
- **Rust** (`rust.yml`, on a pull request that touches the Rust code, the shim or `versions.env`): times a cold and a cached build, runs `og` by bare name in the four images on amd64 and arm64 (QEMU), runs that check once more with a wrong-arch binary planted, which must fail, runs the fetch smoke with both binaries (arm64 under QEMU) and with pysam in `PYPGX_IMAGE`, which must all match `fetch.pysam.tsv`, and keeps the binaries as the artifact `og-binaries`. Every night it runs the differential test at 10,000,000 cases per operation with a new seed.
- **Release** (`release.yml`): a tag builds the four binaries, writes `SHA256SUMS`, attests them and attaches all five files to the release. A pull request that changes `rust/` or the workflow builds them as a dry run.

[Testing](testing.md) has the rest of CI.
