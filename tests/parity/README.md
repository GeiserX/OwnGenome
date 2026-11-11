# tests/parity: og against the Python it replaces

Each OwnGenome tool moves from Python to the Rust binary `og` in one pull request. That pull request deletes the Python, so the Python must stay runnable to prove the port writes the same thing. This folder holds what `scripts/ci/oracle-diff.sh` needs for that proof. The Python is called the *oracle*.

| File | What it holds |
|---|---|
| `oracle.ref` | One row per ported tool: the commit that still had the Python, where it was (`heredoc:FILE:NAME` or `file:FILE`), its interpreter, and the `og` words that replace it. |
| `cases.tsv` | One row per case: tool, case name, inputs, arguments, the outputs it must write, the compare kind for its output files, and the `versions.env` image the oracle runs in. |
| `inputs/<tool>/` | The case inputs, staged into each run by their base name. A staged input named `stdin` is the standard input of both sides; without one it is empty. `.gitattributes` keeps their line endings as they are. |
| `gen/` | Generators for inputs too large to commit, named in `inputs` as `gen:<script>`. |
| `golden/<tool>/<case>/` | What `og` must write where it fixes the oracle on purpose: one file per item (`exit`, `stdout`, `stderr` or an output name), or `<item>.absent` when `og` must write nothing. |
| `known.tsv` | Differences that cannot be stable: zip and deflate bytes, clocks, mtimes, Python exception text. |
| `shadow-self-test.sh` | The self-test of `bin/og`'s shadow mode. |

## Running it

```bash
scripts/ci/oracle-diff.sh                  # every case; builds rust/target/release/og when needed
scripts/ci/oracle-diff.sh setup-mask-fasta/overlap
scripts/ci/oracle-diff.sh --self-test      # planted differences must each be reported
OG=path/to/og scripts/ci/oracle-diff.sh    # another og binary
```

Each case stages the same inputs into two directories. `og` runs in one directory. The oracle runs in the other, inside its image with no network. Both get the same arguments. The harness then compares the exit code, stdout, stderr, the outputs the case names and every new file either side wrote. A staged input that a tool changes in place is compared only when the case names it as an output, so name such a file in the `outputs` column. Files use the case's compare kind. The exit code, stdout and stderr compare as bytes. `og dev parity-compare` does each comparison. Its kinds are `bytes`, `gz-content`, `content` (`gz-content` when the oracle's file is gzip, else `bytes`, for a case whose outputs mix the two), `zip-content`, `png-pixels`, `tsv[:keys=A,B][:unordered]`, `kv`, `vcf`, `json[:paths=..][:mask=..][:tol=..]`, `mask:REGEX` and `stderr-first-line`. Its own `--self-test` plants a difference of each kind.

Each difference has one of three outcomes:

- A *quirk* is a byte that cannot be stable. It gets a `known.tsv` row with the reason.
- A *deliberate fix* gets a golden file and a line in the pull request.
- Anything else is a bug in the port, and the run is red.

A `known.tsv` row whose difference is gone fails. A golden that the oracle also matches fails too. So neither list can go stale.

## Shadow mode in e2e

`bin/og` has a shadow mode for CI. `OG_SHADOW` names tools of `oracle.ref`, comma or space separated. `scripts/ci/oracle-diff.sh materialise` writes their oracles into `bin/.og/oracle/`. An `og` call that starts with such a tool's words then runs `og` as the real output. It also runs the oracle in a scratch directory of symlinks to what the working directory held before, with the same arguments. The exit code, stdout and every new top-level file each side wrote must be the same, byte for byte, or the call fails. The call also fails when either side changes a file that was there before or writes into a directory that was there before: the oracle would write those through its symlinks into the same tree, so the shadow cannot compare them. stderr is not compared. A call that agrees prints `og shadow (<tool>): og and the Python it replaces agree` on stderr, so a log shows the shadow ran. The e2e workflow sets `OG_SHADOW` for ported tools, so every case that calls one checks it on real data. Shadow mode compares bytes only and knows no goldens, so a tool leaves `OG_SHADOW` once a later change makes `og` differ from its oracle on purpose. The tool must be called with paths relative to its working directory.

## Adding a port

1. Write the `og` subcommand.
2. Add a row to `oracle.ref` with the parent commit of your pull request. Add the cases: the old unit tests' inputs, and the traps the old code had. When the old program took its options in another form (awk's `-v name=value`), stage a small adapter with the case inputs and name it as the interpreter, as `sh ./oracle.sh` does for `vcf-csq-rank` and `table-constraint-join`: it turns `og`'s options into the old ones. Such an oracle runs only in the harness, not in shadow mode, which has no adapter in the task's directory.
3. Run `scripts/ci/oracle-diff.sh`. Turn each difference into a fix, a golden or a `known.tsv` row.
4. Switch the callers to `og`, delete the Python, and add the tool to `OG_SHADOW` in `.github/workflows/e2e.yml`.
5. Link a red run in the pull request: the cases with one planted change in the port must fail.
