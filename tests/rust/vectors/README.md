# CPython answers the Rust tests replay

Some Rust tests check og against the Python it replaced or against CPython itself. They used to run Python, or old Python tools out of git history, on every test run. They now read the answers that Python gave, recorded once and committed here. No test installs or runs Python.

| File | Read by | What it holds |
|---|---|---|
| `pycompat/cases.tsv.gz` | `rust/og/src/probe.rs`, module `vectors` | the hard cases (every hand-picked float of the oracle and its negation for each float operation, at precisions 0, 1, 2, 3 and 17 and ndigits -1 to 3; the fixed number, integer and JSON texts and each whitespace character), then 1000 random cases per `og pycompat-probe` operation (50 for the five that walk a folder), each case line followed by CPython's answer line. A `chars` answer is `@CHARS@`. |
| `pycompat/chars.tsv.gz` | the same | the `chars` answer of each code point: `R lo hi F` for a run of code points that are their own upper and lower case, `X cp F upper lower` for the rest |
| `pycompat/tree.tsv.gz` | the same | the folder the glob, first, newest, mtime and mtimens cases walked (`/tree`), listed by GNU find: type, path, link target, mtime |
| `scanner/oracle.jsonl.gz` | `rust/crates/og-dev/src/personal_data.rs` | what the Python personal-data scanner of commit 4104a93f printed on each test case, keyed by a fingerprint of the arguments, the files and the git index |
| `scanner/spans.jsonl.gz`, `scanner/model.jsonl.gz` | `rust/crates/og-dev/src/lookaround.rs`, `personal_data.rs` | CPython 3.15's `re.finditer` spans on each pattern and text (`model` holds a fingerprint of the spans of each pattern over all texts) |
| `scanner/regex_case.jsonl.gz` | `personal_data.rs` | `re.IGNORECASE` relations of CPython 3.9, 3.12, 3.15 and the image's own 3.13 (3.13.5): the pairs that differ, the singleton relation graph, and each interpreter's spans on each fold class |
| `scanner/unicode.jsonl.gz` | `personal_data.rs` | the scalars whose `str` classes or case mappings differ between those four interpreters, and the scalars Unicode 17 leaves unassigned |
| `ci_readers.jsonl.gz` | `rust/crates/og-dev/tests/ci_readers.rs` | whether each replaced CPython smoke assertion of `tests/smoke/commands.tsv` at commit 669a8196 passed on each control input |
| `runlog.jsonl.gz` | `rust/crates/og-runlog/tests/parity.rs` | exit code, stdout, stderr and written files of `bin/run_state.py`, `run_headroom.py` and `run_delta.py` of commit 630292bd on each case |
| `live_ops.jsonl.gz` | `rust/crates/og-live/tests/ops.rs` | the card and metrics text of `bin/owngenome_live_ops.py` of commit 1ec63069 on each case |

A `.jsonl.gz` file holds one `[key, answer]` JSON array per line, sorted, gzipped with `gzip -9 -n`. Paths in an answer are written as `{TMP}`, `{TMP:<name>}` or `{REPO}`, and the test puts its own folders back. The scanner and `ci_readers` keys are made from the input itself or a fingerprint of it: a test whose input changed finds no answer for it and fails with a message that says so. The pycompat cases sit beside their answers. The `runlog` and `live_ops` keys are the case labels: when a case's input changes, record its answer again, because the test still finds the old one.

## How they were made

Commit 00596ca8 ("test(rust): record the CPython answers the oracle tests compare against") still runs the Python and writes each answer it compares against when `OG_RECORD_VECTORS` is set. Its `tests/rust/vectors/make.sh` ran:

1. the pycompat oracle (`tests/rust/pycompat_oracle.py`) with `python3 -I` in `PYTHON_IMAGE`, 1000 cases per operation and seed 20261126, then GNU find over the folder it built;
2. the tests above in the official Rust image (Debian trixie) of the release in `rust/rust-toolchain.toml`, with CPython 3.9, 3.12 and 3.15 from uv 0.12.18 and the image's own `/usr/bin/python3` (3.13.5), the old tools fetched from their commits;
3. `sort -u` and `gzip -9 -n` into this folder.

Commit 55691a5e ("test(rust): record the pycompat hard cases before the random draws") recorded `pycompat/cases.tsv.gz` again, with the hard cases before the same random cases. It is that oracle with `hard_cases()` added, and a `tests/rust/vectors/make.sh` that runs only step 1 and writes `cases.tsv.gz` and `tree.tsv.gz`. The random cases and the tree came out byte for byte as before. `chars.tsv.gz` and every other file are still the ones 00596ca8 wrote.

To record a new pycompat answer, check out 55691a5e, change `tests/rust/pycompat_oracle.py` there, run `tests/rust/vectors/make.sh WORK` with Docker (an empty scratch folder; it takes seconds), and bring `cases.tsv.gz` to the current branch. For any other answer, do the same from 00596ca8 and its full `make.sh`. Both commits stay in the history of main, because pull requests land as merge commits. `rust/og/src/probe.rs` checks that every hand-picked float is a case of each float operation.
