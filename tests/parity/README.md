# Rust parity cases

The parity harness runs each case in a fresh working directory and compares
its exit code, stdout, stderr, named outputs, and new files with committed
goldens. The goldens came from the archived implementations and the existing
goldens that recorded deliberate changes in the Rust ports.

| Path | Contents |
|---|---|
| `cases.tsv` | Tool, case, inputs, arguments, named outputs, comparison kind, and the words of the og command. |
| `inputs/` | Files staged by base name. If a staged file is named `stdin`, the tool reads it as standard input. |
| `gen/` | Generators that stage inputs. The large synthetic VEP, ClinVar, consumer, and IBD generators use fixed seeds. |
| `golden/fixtures/` | Small fixture archives captured from the retired generators and historical catalogs. |
| `golden/TOOL/CASE/` | Expected exit, stdout, stderr, and output files. |
| `known.tsv` | Record of archived oracle differences and deliberate golden updates. |

Inputs are comma-separated. Use `gen:SCRIPT` for a generator, `fixture:NAME` for
a frozen input archive, or `-` for an empty input list. Name an in-place edit in
the outputs column so the harness compares it. The comparison kind for exit,
stdout, and stderr is bytes. An output uses the case's comparison kind or its
matching file glob. The [Rust comparator](../../rust/crates/og-dev/src/parity.rs) lists the supported comparison kinds.

An `ITEM.absent` file records that the tool must leave ITEM absent. An
`ITEM.sha256` file holds the SHA-256 of the canonical stream for that item's
comparison kind. Large output files and selected public or synthetic fixture
outputs use checksums. The checksum for a ZIP covers each member's name
and uncompressed contents. Each case can store up to 1,000,000 bytes across
its golden files. The Yleaf data-folder cases replace the temporary working
directory with `<case-dir>` in their diagnostic. The two known Paraphase
exception differences keep the diagnostic from og.

## Run and inspect cases

Run these commands from the repository root.

	scripts/ci/oracle-diff.sh
	scripts/ci/oracle-diff.sh setup-mask-fasta/overlap
	scripts/ci/oracle-diff.sh --self-test
	scripts/ci/oracle-diff.sh --check-goldens

Set `OG` to the absolute path of a prebuilt binary to use it. Set `PARITY_WORK` to
an absolute scratch directory to retain the case's files and actual results.
When a change intentionally alters an expectation, review the case output in
that directory and update the affected golden in the same pull request. Record
the changed case, output, and reason in `known.tsv` when the current expectation
intentionally differs from the archived oracle. To recompute a checksum, pass
the comparison kind and the actual filename to
`og dev parity-compare --canonical` and pipe the result to `sha256sum`.

The `cn-engine` goldens come from the per-gene tools at
605645249e6762c52a84ac58a28968f269b16536. The E2E case checks HBA regions and
tables, the RhD rows, and the LPA table and interval against those results.
Its reads are public HG002 reference data at full depth and one eighth depth.
