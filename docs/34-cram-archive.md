# Step 34: CRAM Archive

## What This Does
Writes the sample's alignments as CRAM beside the BAM and checks their read content. It deletes the BAM only when you ask, after the checks pass. `--restore` writes the BAM back from the CRAM.

## Why
A 30x genome keeps 80 to 120 GB of BAM after the analysis is done. A CRAM stores each read as its difference from the reference and is about half that size or less. The CI fixture has 68 MB of CRAM for 133 MB of BAM. Read content is checked before the original can be deleted. Most people keep their alignments to rerun a step when a tool or database improves. A CRAM keeps that option for half the disk.

The price is the reference: a CRAM can only be read with a FASTA that holds the same sequences as the one it was written with. Keep `reference/GRCh38_no_alt_analysis_set.fasta` (or your `REF_FASTA`) as long as you keep the CRAM. samtools checks each contig's MD5 when it reads a CRAM and stops with an error on a different file, so a wrong reference cannot give silently wrong reads.

## Tool
- **samtools** (`view -C`, `index`, `quickcheck`, `flagstat`, `checksum`)

## Docker Image
- `SAMTOOLS_IMAGE`

Pinned in `versions.env`. [Image versions](versions.md) lists the current tag.

## Command
```bash
./scripts/34-cram-archive.sh your_name               # write and check the CRAM; the BAM is kept
./scripts/34-cram-archive.sh your_name --delete-bam  # the same, then delete the BAM and its index
./scripts/34-cram-archive.sh your_name --restore     # write the BAM back from the CRAM
```

The script runs, in the samtools container:
```bash
samtools view -C --reference "${REF_FASTA}" --output-fmt-option store_md=1 --output-fmt-option store_nm=1 -o aligned/${SAMPLE}_sorted.part.cram aligned/${SAMPLE}_sorted.bam
samtools quickcheck -v aligned/${SAMPLE}_sorted.part.cram
samtools index aligned/${SAMPLE}_sorted.part.cram
samtools flagstat aligned/${SAMPLE}_sorted.bam
samtools flagstat --input-fmt-option reference="${REF_FASTA}" aligned/${SAMPLE}_sorted.part.cram
samtools checksum -PCMc -b 0xfff -f0 -F0 -z all,cigarx -t '*' aligned/${SAMPLE}_sorted.bam
samtools checksum -PCMc -b 0xfff -f0 -F0 -z all,cigarx -t '*' --input-fmt-option reference="${REF_FASTA}" --input-fmt-option decode_md=0 aligned/${SAMPLE}_sorted.part.cram
```

The CRAM must pass `samtools quickcheck`, the two `flagstat` reports must match, and the content checksums must match. The checksum reads every record, including secondary and supplementary alignments. It covers read names, sequences, qualities, flags, positions, CIGARs, mate fields and all auxiliary tags. A changed base quality or tag can leave `flagstat` counts unchanged and still fail this check. A decoding error also fails the step.

The checksums use samtools' order-independent aggregation. The comparison ignores the filename header and sorts the report rows, so read and tag order do not affect it. samtools sanitizes format-level differences, such as an unmapped read's MAPQ and equivalent CIGARs. Explicit `MD` and `NM` tags are stored in the CRAM. For a CRAM this step wrote, they are not regenerated during the check or restore, so missing tags stay missing. An older or external CRAM, written without stored `MD`/`NM`, keeps only the tags that differ from the reference. Restoring it regenerates `MD` and `NM` from the reference, and the check decodes the CRAM the same way. These CRC-based content checks detect conversion faults. They are not a cryptographic identity proof. See the [samtools 1.24 checksum manual](https://www.htslib.org/doc/1.24/samtools-checksum.html).

Only after all checks pass does the CRAM lose its `.part` name. A mismatch or a failed or empty checksum command removes the partial CRAM, keeps the BAM and its index, and exits non-zero. This also applies with `--delete-bam`. An older published CRAM stays in place if a new attempt fails. Every archive attempt writes a fresh CRAM from the current BAM before it can delete anything. The checksum adds a full read of each file to the archive and restore checks.

`--restore` does the same in reverse: the BAM is written under a `.part` name, checked against the CRAM, and moved into place. It refuses to write over a BAM that exists. A failed check removes the partial BAM and keeps the CRAM.

One archive or restore of a sample runs at a time. Archive, `--delete-bam` and `--restore` all take the same lock, `aligned/${SAMPLE}_sorted.lock`, with `flock` before they read or write any of the sample's files. A second run on the same sample, from another terminal or a Nextflow `CRAM_ARCHIVE` task, stops at once with `ERROR: another archive or restore of <sample> is running`. It changes nothing. Run it again when the first has finished. Without the lock, two runs would delete or overwrite the `.part` files the other is writing and checking, and `--delete-bam` could remove a BAM the other run is still reading. The lock goes when its run ends, including on an error, Ctrl-C or `kill -9`. The empty `.lock` file stays and means nothing on its own. `flock` comes with util-linux on Linux. On macOS install it with `brew install flock`, or the step warns and runs without the lock.

In Nextflow, `cram_archive` in `--tools` runs the same conversion and check (`CRAM_ARCHIVE`) for every BAM of the run, except the rows given as CRAM, which have one already, and publishes `<outdir>/<sample>/aligned/<sample>_sorted.cram`. It never deletes anything: delete the BAM yourself once the CRAM is there, or run this script on it. It takes the same lock, next to the BAM the samplesheet row names (`<BAM name without .bam>.lock`), and fails the task while this script holds it. This script and a `CRAM_ARCHIVE` task exclude each other on Linux, where the container shares the host kernel. On macOS with Docker Desktop, only two runs of this script exclude each other. When it cannot write the lock file there, it warns and reads the BAM without the lock.

## Output Files
| File | Description |
|---|---|
| `aligned/${SAMPLE}_sorted.cram` | The alignments |
| `aligned/${SAMPLE}_sorted.cram.crai` | Its index |
| `aligned/${SAMPLE}_sorted.cram.flagstat`, `aligned/${SAMPLE}_sorted.bam.flagstat` | The two `flagstat` outputs the check compared |
| `aligned/${SAMPLE}_sorted.cram.checksum`, `aligned/${SAMPLE}_sorted.bam.checksum` | Normalized content checksum reports, without the filename header. Nextflow publishes only the CRAM report. |
| `aligned/${SAMPLE}_sorted.lock` | The empty file the lock is taken on |

## Rerunning steps from a CRAM

| How you run | From a CRAM |
|---|---|
| Nextflow | Every step. Put the CRAM in the samplesheet's `cram` and `crai` columns (with the `vcf` columns, or a `sex` column to call it again) and pass the reference it was written with as `--reference`. `CRAM_TO_BAM` writes it out once as a BAM in the work directory, checked against the CRAM, and every BAM step reads that. The BAM needs the disk space of the original until the run's work directory is removed |
| The bash steps that read alignments | None directly: steps 03, 04, 04b, 08, 09, 10, 15, 16, 16b, 18, 19, 20, 21, 28 (its samtools flagstat), 29, 32 and 33, and the alternative callers 03a to 03d and 04a, read `aligned/${SAMPLE}_sorted.bam` (03e and 04c read `aligned_longread/`, which this step does not archive). Run `./scripts/34-cram-archive.sh your_name --restore` first, then the step |
| VCF-only steps | Every one: steps 06, 07, 11, 12, 13, 14, 17, 22, 23, 25, 26, 27, 30 and 31 read the VCF, which this step never touches |

## Runtime
Writing a CRAM reads the whole BAM and compresses it again. The checks read both files for flagstat and again for the checksum. On a synthetic 100,000-read BAM in the pinned image, three trials took 4.44 to 4.99 seconds with 44.1 to 45.3 MiB peak RSS. The count-only checks took 3.39 to 3.99 seconds with 36.3 to 36.5 MiB. These measurements use one thread under amd64 emulation and exclude container startup. A whole genome has not been timed here.

## Notes
- Keep the reference FASTA (and its `.fai`) with the CRAM. A FASTA decodes the CRAM only when every contig the CRAM uses has the same sequence. Keep the file you wrote it with.
- The CRAM is written with samtools' default CRAM version and compression, which every current htslib, Picard and GATK read.
- A CRAM made elsewhere (nf-core/sarek writes them) can go straight into the Nextflow samplesheet. For the bash steps, `--restore` needs it at `aligned/${SAMPLE}_sorted.cram` with its `.crai`.
