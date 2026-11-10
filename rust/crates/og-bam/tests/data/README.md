# og-bam fetch smoke data

| File | What it is |
|---|---|
| `cyp2d6.bam`, `.bai` | the 357 reads of the test fixture's `HG002_slice.bam` (public GIAB HG002, built by `scripts/ci/build-fixture.sh`) that pysam's `fetch("chr22", 42125999, 42128000)` returns, at CYP2D6, with the header unchanged |
| `regions.tsv` | the regions the smoke fetches: the whole slice and the edge cases, each with a comment |
| `fetch.pysam.tsv` | what `../pysam_fetch.py` writes for them with the pysam in `PYPGX_IMAGE` (`versions.env`) |

To make them again, from the fixture directory, in `PYPGX_IMAGE`:

```python
import pysam
src = pysam.AlignmentFile("HG002_slice.bam")
with pysam.AlignmentFile("cyp2d6.bam", "wb", template=src) as out:
    for r in src.fetch("chr22", 42125999, 42128000):
        out.write(r)
pysam.index("cyp2d6.bam")
```

then `python3 -I pysam_fetch.py cyp2d6.bam regions.tsv > fetch.pysam.tsv`. On the fixture of `tests/fixtures/VERSION` this gives the same bytes as the files here.

# Parity data of `og pgx cyp2d6-*` and `og str motif`

`parity/` is what `../parity.rs` checks og against: the Python tools it replaced (`bin/cyp2d6_depth_check.py`, `bin/rfc1_motif.py`) run on the same inputs.

| Path | What it is |
|---|---|
| `parity/in/` | the inputs `../make_parity.py` makes: BAMs of made-up reads at RFC1 and at a made-up second motif locus (the reads of `tests/test_rfc1_motif.py`), the fixture slice above renamed to `22` and to another contig, made-up reads at CYP2D6 for each MAPQ floor and flag, VCFs, catalogs and mosdepth region files with each defect the check refuses |
| `parity/cases.tsv` | one row per case: name, tool, arguments (`{d}` this folder, `{r}` the repository, `{o}` the output folder) |
| `parity/want/<case>/` | what the Python wrote: `exit`, `stdout` and `files/` |
| `parity/og/<case>/` | where og differs on purpose, what og must write instead |

To make them again, from the repository root of a commit that still has the two Python tools, in `PYPGX_IMAGE`:

```bash
docker run --rm --network none -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD:/repo" -w /repo \
  "$PYPGX_IMAGE" python3 -I rust/crates/og-bam/tests/make_parity.py
```

The `og/` folders are written by hand, and `make_parity.py` leaves them as they are.
