# og-bam fetch smoke data

| File | What it is |
|---|---|
| `cyp2d6.bam`, `.bai` | the 357 reads of the test fixture's `HG002_slice.bam` (public GIAB HG002, built by `scripts/ci/build-fixture.sh`) that pysam's `fetch("chr22", 42125999, 42128000)` returns, at CYP2D6, with the header unchanged |
| `regions.tsv` | the regions the smoke fetches: the whole slice and the edge cases, each with a comment |
| `fetch.pysam.tsv` | what pysam's fetch returned for them, with the pysam in `PYPGX_IMAGE` (`versions.env`): per read the region, name, flag, `reference_start`, `reference_end`, `mapping_quality`, `cigarstring` and `get_blocks()` |

`cargo test` (`../fetch_pysam.rs`) holds og to `fetch.pysam.tsv`. The Rust workflow also checks that samtools in `SAMTOOLS_IMAGE` returns the same reads there: `../samtools_fetch.sh` writes the same lines from `samtools view`, run from the repository root.

The BAM and `fetch.pysam.tsv` were made with pysam; this README at commit 182811a7 has the recipe, and `../pysam_fetch.py` there wrote the lines.

# Parity data of `og pgx cyp2d6-*`, `og str motif`, `og cdna`, `og comphet phase` and `og bam contigs`

`parity/` and `reads/` are what `../parity.rs` checks og against: what the Python tools it replaced (`bin/cyp2d6_depth_check.py`, `bin/rfc1_motif.py`, `bin/clip_check.py`, `bin/cdna_scan.py`, `bin/comphet_phase.py`) wrote on the same inputs. Those tools are gone, so these answers are frozen.

| Path | What it is |
|---|---|
| `parity/in/` | BAMs of made-up reads at RFC1 and at a made-up second motif locus, the fixture slice above renamed to `22` and to another contig, made-up reads at CYP2D6 for each MAPQ floor and flag, VCFs, catalogs and mosdepth region files with each defect the check refuses |
| `reads/in/` | BAMs, tables and VCFs of made-up reads: clipped, spliced and scanned cDNA reads, and the comphet fixture plain and `--swap` |
| `*/cases.tsv` | one row per case: name, tool, arguments (`{d}` this folder, `{r}` the repository, `{o}` the output folder) |
| `*/want/<case>/` | what the Python wrote: `exit`, `stdout` and `files/` |
| `*/og/<case>/` | where og differs on purpose, what og must write instead (written by hand) |

The inputs and answers were made with pysam by `../make_parity.py` and `../make_reads_parity.py`, run on a commit that still had the Python tools. Both scripts are at commit 182811a7 and earlier.
