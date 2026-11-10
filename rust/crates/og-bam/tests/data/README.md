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
