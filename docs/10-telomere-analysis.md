# Step 10: Telomere Length Estimation

## What This Does
Estimates relative telomere content from WGS BAM files by quantifying telomeric repeat reads (TTAGGG/CCCTAA).

## Why
Telomere length correlates with cellular aging at a population level. Comparing telomere content between individuals of similar age, sequenced on the same platform, provides a rough relative comparison. However, telomere length alone is not established as a clinically important standalone risk marker for individuals — it provides only a rough estimate of aging rate and is influenced by many non-age factors (genetics, cell type, technical variables).

## Tool
- **TelomereHunter2** 1.0.12 (German Cancer Research Center; GPL-3.0), the maintained Python 3 successor of TelomereHunter

## Docker Image
- `TELOMEREHUNTER_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

> The image is the Bioconda build of TelomereHunter2. It runs as you, not as root, and with `--plotNone`: nothing in the pipeline reads its plots.

### Which TelomereHunter, and what changed

Step 10 runs TelomereHunter2 1.0.12. Until now it ran TelomereHunter 1.1.0, Python 2.7 code with its plots off. On the image test's input, the HG002 fixture BAM (fixture-v8) plus 1,100 planted unmapped telomeric reads (600 `TTAGGG`, 400 `CCCTAA`, 100 with one `TCAGGG`), both were run with the GRCh38 bands. The second run drops two planted reads:

| Input | Image | intratelomeric reads | tel_content | TCAGGG per intratelomeric read |
|---|---|---|---|---|
| 1,100 planted | TelomereHunter 1.1.0 | 1100 | 5347.255643 | 0.0909 |
| 1,100 planted | TelomereHunter2 1.0.12 (this step) | 1100 | 5347.255643 | 0.0909 |
| 1,098 planted | TelomereHunter 1.1.0 | 1098 | 5337.585253 | 0.0911 |
| 1,098 planted | TelomereHunter2 1.0.12 (this step) | 1098 | 5337.585253 | 0.0911 |

So a value from an earlier run of this step stays comparable, and both versions move by the same amount when reads change. The summary file keeps its name, but TelomereHunter2 writes other columns in another order: `tel_content` is the third column, `intratel_reads` is now `intratelomeric_reads`, `tel_reads` is `tel_read_count`, `read_length` is `read_lengths`, and it adds the junction-spanning, subtelomeric and intrachromosomal read counts. The reports and `summary.json` read both layouts by column name, so they do not change. A script of your own that reads the summary by column position must read it by name.

## Command
```bash
export GENOME_DIR=/path/to/your/data
./scripts/10-telomere-hunter.sh your_sample
```

`ALIGN_DIR` (default `aligned`) picks the BAM, for example `ALIGN_DIR=aligned_bwamem2`. `THREADS` (default 4) caps the container's CPUs and is TelomereHunter2's worker count (`-c`).

TelomereHunter2 sorts telomeric reads into intratelomeric, subtelomeric, junction-spanning and intrachromosomal classes by chromosome band. The BAM is GRCh38, so the script always passes UCSC's GRCh38 bands, `reference/cytoBand.hg38.txt` (chr1-22, X and Y; see [reference setup](00-reference-setup.md#small-pinned-data-files)). The repository ships that file in `assets/reference`, and `setup.sh` installs it. When it is not installed, the script installs it first and says so; it never runs without the bands. (TelomereHunter 1.1.0 used its own hg19 band ends without `-b`, and TelomereHunter2 leaves the band classes out.) The Nextflow pipeline passes `--cytoband` when it is set, and otherwise the copy in `assets/reference`.

TelomereHunter2 writes one BAM per contig and merges them at once, so a BAM with decoy contigs (the HG002 test BAM has 2,580) needs more than the 1,024 open files a container gets by default. The script raises the container's limit with `--ulimit nofile=65536:65536`; the Nextflow module raises its soft limit to the hard one. Without that the run stops with "Too many open files". What the script runs:

```bash
source versions.env   # from the repository root
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

docker run --rm --network none --user "$(id -u):$(id -g)" \
  --cpus 4 --memory 4g --ulimit nofile=65536:65536 \
  -v ${GENOME_DIR}:/genome \
  "${TELOMEREHUNTER_IMAGE}" \
  telomerehunter2 \
    -ibt /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
    -o /genome/${SAMPLE}/telomere/${SAMPLE} \
    -p ${SAMPLE} \
    -c 4 \
    --plotNone \
    -b /genome/reference/cytoBand.hg38.txt

# Output: telomere content report in ${GENOME_DIR}/${SAMPLE}/telomere/${SAMPLE}/
```

TelomereHunter2 uses `-ibt` (input BAM tumor) for a single-sample analysis. No `--tumor_only` flag is needed: with only `-ibt` (no `-ibc` for a matched control BAM), it runs in single-sample mode.

## Key Metric
- **`tel_content`** — intratelomeric reads per million reads with 48-52% GC, the GC content of telomeric repeats (TelomereHunter2 divides the intratelomeric read count by `total_reads_with_tel_gc`, as TelomereHunter 1.1.0 did). It is a relative telomere content, not a telomere length; the report labels it "Telomere content (relative)"
- Higher values indicate longer/more abundant telomeres
- Compare only between samples sequenced on the same platform. The value is not an age

### Read-count check

The value rests on the intratelomeric reads (`intratelomeric_reads` in TelomereHunter2's summary, `intratel_reads` in a summary TelomereHunter 1.1.0 wrote). The reports print that count beside the value and warn "too few telomeric reads for a stable estimate" when it is below 2,000 per 30x of depth. The depth is the one TelomereHunter counted: `total_reads` times the read length (`read_lengths`, the largest when there are several) over 3.1 Gb. A typical 30x short-read genome has tens of thousands of intratelomeric reads, so the floor sits far below that: a sample under it has a library or platform problem, or reads too short for TelomereHunter's repeat threshold, and its value is not comparable. The floor is a first choice, not one fitted to many samples. A summary without `total_reads` or the read length gets the count and no check, and the report says so.

## Important Notes
- The container runs as you (`--user "$(id -u):$(id -g)"`), so the outputs are yours. A telomere directory an older version wrote as root has to be taken back first: `sudo chown -R "$(id -u):$(id -g)" "${GENOME_DIR}/${SAMPLE}/telomere"`.
- Short-read WGS (150bp reads) systematically underestimates true telomere length because reads cannot span long repetitive regions
- Results are useful as a **relative comparison** between samples, NOT as an absolute telomere length measurement
- Long-read sequencing (PacBio/ONT) provides more accurate telomere length if absolute values are needed
