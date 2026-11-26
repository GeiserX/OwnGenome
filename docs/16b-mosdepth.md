# Step 16b: mosdepth Coverage Statistics

Fast per-base and per-region depth statistics from BAM. Produces coverage distributions, threshold reports, and a genome-wide summary.

---

## What It Does

Reads the BAM file and computes:
1. **Mean coverage** per chromosome and genome-wide
2. **Coverage distribution** — cumulative fraction of bases at each depth level
3. **Region coverage** — mean depth in 500bp windows across the genome
4. **Threshold report** — fraction of bases at 1x, 5x, 10x, 15x, 20x, 30x, 50x depth
5. **Gene callability** — for each gene of `assets/callability_genes.tsv` on chr1-22, X and Y, and for each coding exon of those genes in `assets/callability_exons.tsv`, the bases with 10 or more reads of mapping quality 20 or more (a second mosdepth pass, which counts the two mates of a pair once where they overlap)

## Why

- **Validates sequencing depth** — confirms your WGS actually has 30X (or whatever depth was ordered)
- **Identifies low-coverage regions** — gaps in coverage cause missed variant calls
- **Detects sex chromosome anomalies** — chrX depth relative to autosomes reveals XX vs XY
- **Says which hard genes this genome could not read** — in a gene inside a segmental duplication, a false duplication of GRCh38 or on an ALT contig, short reads often get mapping quality 0, so no caller sees them and a "no finding" there means nothing. The report lists those genes (see [Gene callability](#gene-callability))
- **More precise than indexcov** — indexcov (step 16) reads only the BAM index for a quick estimate; mosdepth reads actual alignments

## Tool

**mosdepth** v0.3.14 — fast BAM/CRAM depth calculation.

- Paper: Pedersen & Quinlan, Bioinformatics 2018 (doi:10.1093/bioinformatics/btx699)
- Source: [github.com/brentp/mosdepth](https://github.com/brentp/mosdepth)

## Docker Image

- `MOSDEPTH_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Command

```bash
export GENOME_DIR=/path/to/data
./scripts/16b-mosdepth.sh <sample_name>
```

## Output

| File | Description |
|---|---|
| `mosdepth/<sample>.mosdepth.summary.txt` | Per-chromosome mean coverage + genome total |
| `mosdepth/<sample>.mosdepth.global.dist.txt` | Cumulative coverage distribution (for plotting) |
| `mosdepth/<sample>.regions.bed.gz` | Mean depth per 500bp window |
| `mosdepth/<sample>.thresholds.bed.gz` | Fraction of bases at 1x/5x/10x/15x/20x/30x/50x |
| `mosdepth/<sample>_callability.bed.gz` | Per gene of `assets/callability_genes.tsv` with a position: chrom, start, end, gene, and the bases with 10 or more reads of MAPQ 20 or more (mosdepth's thresholds table, column `10X`) |
| `mosdepth/<sample>_callability_exons.bed.gz` | The same per coding exon of those genes in `assets/callability_exons.tsv`, the region named `GENE:exonN`; only when the exon table names a gene of the gene table |

The Nextflow MOSDEPTH task writes the same files into `<outdir>/<sample>/coverage/`; it measures exons when it is given `assets/callability_exons.tsv` with the gene table, as the pipeline does.

## Interpreting Results

### Summary file

```
chrom   length      bases          mean    min  max
chr1    248956422   7394254892     29.70   0    312
chr2    242193529   7185433827     29.67   0    290
...
total   3088286401  91784821430    29.72   0    312
```

- **mean ~30**: You have 30X whole-genome coverage (good for germline calling)
- **mean < 15**: Variant calling accuracy degrades significantly
- **chrX mean ~15 (male) or ~30 (female)**: Expected sex-linked depth

### Threshold report

Shows what fraction of each region is covered at key depths:
- **>95% at 10x**: Adequate for most germline variant callers
- **>90% at 20x**: Good quality WGS
- **>80% at 30x**: High-quality WGS

### Gene callability

`assets/callability_genes.tsv` lists 1,572 genes where short reads are often not placed with confidence. `scripts/ci/build-callability.sh` writes it from pinned public files, and a CI case rebuilds it and fails on any difference:

| Reason | Genes | Source |
|---|---|---|
| `CMRG` | 273 | GIAB's challenging medically relevant genes, CMRG v1.00 (Wagner et al., Nat Biotechnol 2022, PMID 35132260) |
| `segdup` | 1,276 | protein-coding genes with half or more of the gene body in a segmental duplication (GIAB genome stratifications v3.6, GRCh38) |
| `false_duplication` | 22 | genes over a false duplication of GRCh38, either copy: KCNE1, CBS and CRYAA among them (Behera et al., Genome Biol 2023, PMID 36810122; the same GIAB stratifications) |
| `collapsed` | 50 | genes over a collapsed duplication: two copies in people, one in GRCh38 (GIAB stratifications) |
| `ALT-only` | 19 | protein-coding genes only on an ALT scaffold of GRCh38, so the default no-ALT reference cannot hold them: GSTT1, HLA-DRB3, HLA-DRB4 and others (GENCODE 50, NCBI's GRCh38 assembly report) |

A gene can have several reasons. Gene positions are GENCODE 50 gene spans. The `owner` column names the step that calls the gene with a method of its own: `parascopy` (step 35, the genes of `assets/paralog_loci.tsv`), `pypgx` (step 32, the bam genes of `assets/pgx_genes.tsv` it can call on the no-ALT reference, so not GSTT1) or `cn_loci:<locus>` (the depth engine, `bin/cn_loci.tsv`: the genes a row names, HBA1, HBA2, LPA and RHD, not a gene that only holds one of its targets); `.` for none.

A gene counts as **read** when 90% or more of its bases have 10 or more reads of MAPQ 20 or more. The gene pass runs mosdepth without `--fast-mode`, so a base where the two mates of one pair overlap counts one read, not two, and a base inside a deletion of a read's CIGAR counts no read. With `--fast-mode` mosdepth counts every base from a read's start to its end once per read: short-insert pairs, whose mates overlap, and small deletions lift a base to 10x when fewer than 10 fragments read it, and the gene reads as read. The genome-wide pass (the 500 bp windows and their thresholds) keeps `--fast-mode`: it describes the run's depth, and nothing is decided per gene from it. Both reports have a section "Genes this genome could not read reliably" in the Quality chapter. It lists every gene with no owner that was not read, with the share of its bases that reached that depth, and every ALT-only gene with no owner as "not on this reference". Genes with an owner that were not read are listed apart: their own step's section says what it called. In `summary.json` each of them (`not_read_owned`) also names that step (`owner_step`: `35` for Parascopy and the HBA, CYP21A2, HP and CCL3L1 rows of the depth engine, `32` for pypgx, `46` for LPA, `41` for RHD) and says whether it ran for this sample (`owner_status`: `ran`, `not run` or `from an earlier run`). A step ran when it wrote the gene's file in the sample folder: for a Parascopy gene its own table, `paralogs/<sample>_<table>_copy_number.tsv` with the `table` of its row in `assets/paralog_loci.tsv` (the `smn` table for SMN1 and SMN2), `pypgx/<sample>_pypgx_summary.tsv`, `paralogs/<sample>_hba_copy_number.tsv`, `_cyp21a2_psv.tsv` or `_cn_dosage.tsv`, `lpa/<sample>_lpa_kiv2.tsv`, or `bloodtype/<sample>_bloodtype.tsv`. A file older than the start of the last run, from a step that was not `ok` in that run (`logs/run_status.tsv`), is `from an earlier run`: the rule that makes a section of `summary.json` stale. A gene whose own step did not run has no answer from it, only its coverage here. The share is the share of bases at the depth, to four decimals. On a partial BAM (the CI fixture, an exome) most genes are not read, because they have no reads at all.

A gene also needs each of its **coding exons** read. 90% of a whole gene body can reach 10x while one exon has no read at all: the introns are most of the body, so a whole exon can fail and the gene still pass. `assets/callability_exons.tsv` lists the coding exons of each catalog gene with a position: the CDS records of the gene's MANE Select transcript in GENCODE 50, each with 2 bases more on both sides for the splice site, named `GENE:exonN` by the transcript's exon number. `scripts/ci/build-callability.sh` writes it beside the gene table from the same pinned GTF, and a CI case rebuilds it and fails on any difference. It has 11,888 exons of 1,501 genes; 52 catalog genes with a position have no MANE Select CDS on chr1-22, X or Y and are judged on their body alone: a gene whose MANE Select transcript codes no protein (H19), one whose MANE Select transcript sits on a patch or ALT scaffold of GRCh38 (ABO, NBPF1, GSTT2), and one GENCODE gives no MANE Select tag (MUC1, the T-cell receptor segment TRBV9). The exons ride in the gene pass: mosdepth reads the BAM once and measures each region on its own, and step 16b writes the rows named `GENE:exonN` to `<sample>_callability_exons.bed.gz` and the gene rows to `<sample>_callability.bed.gz`. The two tables stay apart because the report keys gene rows by name. A gene whose body reached the share but one of whose coding exons did not is **not read**: `summary.json` lists it with `status` `exon not read`, the share of its body, and `worst_exon` (`exon`, `share`), its exon with the lowest share. The text and HTML reports list it with the genes not read, with the share of its body; the exon it failed on is named in `summary.json`. A gene not read on its body also names its worst exon when one is under the share. `exons_checked` counts the exon rows of catalog genes the table held. A folder without the exon table (one from before it existed) gives the summary it gave before: no gene is judged by its exons and neither key appears.

`CALLABILITY_GENES` names another table with the same columns (the CI case uses a two-gene one), and `CALLABILITY_EXONS` another exon table; only the exons of the genes of the gene table are measured, so a gene table none of whose genes the exon table names writes no exon table. A table from an older catalog, one that lacks a gene of the current catalog, makes the section unreadable: run step 16b again. When only `<sample>_callability.bed.gz` is missing, or only `<sample>_callability_exons.bed.gz` while the exon table names a gene of the gene table, step 16b runs the gene pass alone.

The CMRG gene list and the stratification BEDs are only inputs of the build script; a run needs only the committed table.

## Runtime

| Dataset | Threads | Time | Memory |
|---|---|---|---|
| The genome-wide pass, 30-40x WGS, measured in three runs ([runtime](hardware-requirements.md#runtime-per-step)) | 4 | 2-3 min | 1.9 GB |
| The gene callability pass, public HG002 30x, measured in two runs (below) | 4 | 3 min 26 s - 3 min 50 s | 2.0 GB |
| The same pass with the coding exons, public HG002 30x, measured in two runs against two runs without them (below) | 4 | 2-7 s more (about 3%) | not measured |
| chr22 BAM (estimate) | 2 | < 30 sec | < 1 GB |

The gene pass was measured on Google's public HG002 NovaSeq PCR-free 30x BAM on GRCh38 (`HG002.novaseq.pcr-free.30x.dedup.grch38.bam`, 43.7 GB, from `storage.googleapis.com/brain-genomics-public`), with mosdepth 0.3.14 and 4 threads, on a 10-core Mac mini that was running other jobs. Without `--fast-mode` it took 3 min 26 s and 3 min 50 s; with it, as before, 1 min 46 s and 2 min 29 s. The CPU time is the steadier figure: about 7 min without `--fast-mode` against 5 min with it, 37% more. The genome-wide pass took 3 min 15 s on the same machine. Peak memory was 2.0 GB in every run.

On that genome 934 of the 1,553 catalog genes with a position are read, against 938 with `--fast-mode`. The share of bases at 10x fell by 0.2 points on average and by 7 points at most. Four genes moved from read to not read: DHX40, SPATA31C2 and GOLGA6L10 were at 90.1-90.2% and are now at 89.6-89.9%, and UGT2B17 (step 32 calls it) went from 91.7% to 89.9%. No gene moved the other way.

The coding exons were measured on the same public BAM, with mosdepth 0.3.14 and 4 threads, on a 20-core Intel test machine that was running other jobs. The gene pass alone took 2 min 6 s and 2 min 9 s (5 min 37 s and 5 min 51 s of CPU); with the 11,888 exons in the same pass it took 2 min 13 s and 2 min 11 s (6 min 0 s and 5 min 52 s of CPU). The exons add 2 to 7 seconds, about 3%: mosdepth reads the BAM once either way, and the extra regions cost only their own sums. The gene rows of the shared pass are byte for byte those of the gene pass alone.

On that genome 71 of the 934 catalog genes read on their body have a coding exon under 90% at 10x, so they now read `exon not read`. The gene, its worst exon and that exon's share: ACOT1 exon3 88.5%, ACR exon5 83.1%, AGAP5 exon7 28.6%, ANAPC1 exon32 0.0%, ANKRD30B exon23 0.0%, ANKRD36C exon58 0.0%, ARSD exon9 27.0%, ASAH2B exon5 0.0%, BCAP31 exon2 0.0%, BMP8A exon4 0.0%, BMP8B exon3 0.0%, CASTOR2 exon1 0.0%, CCDC144A exon10 45.2%, CDK11B exon20 0.0%, CEL exon11 88.6%, CYP4Z1 exon3 0.0%, FAM90A20 exon1 0.0%, FAM9C exon5 44.5%, FBXW10 exon14 71.9%, FCGR1A exon2 24.0%, GAGE10 exon5 0.0%, GOLGA6L2 exon2 33.1%, GOLGA8H exon4 0.0%, GOLGA8M exon15 30.2%, GPR89A exon6 0.0%, GPR89B exon6 0.0%, GYG2 exon6 0.8%, GYPA exon2 31.1%, HERC2 exon45 83.7%, HLA-A exon6 0.0%, HP exon5 81.1%, HRNR exon3 88.2%, HYDIN exon83 0.0%, IQCK exon7 0.0%, KIR3DL3 exon5 71.5%, KRTAP5-10 exon1 82.3%, MAGEB6 exon2 85.6%, MAGEB6B exon1 87.8%, MUC12 exon2 85.1%, NBPF11 exon15 10.7%, NBPF12 exon36 0.0%, NBPF15 exon19 26.8%, NBPF3 exon14 0.0%, NBPF9 exon29 0.0%, NOMO1 exon6 0.0%, OTOA exon23 0.0%, PMS2 exon15 82.8%, POM121C exon13 77.7%, PRAMEF10 exon2 78.3%, PRKX exon6 0.0%, PSG9 exon5 83.4%, RFPL4AL1 exon3 76.5%, RHD exon8 4.8%, SAGE1 exon19 76.5%, SLC35E2B exon3 65.3%, SLC6A8 exon1 59.4%, SPANXB1 exon1 79.8%, SPANXN1 exon1 72.2%, SPANXN2 exon2 81.5%, SPATA31A6 exon4 87.7%, SPIN2A exon2 37.8%, SRGAP2 exon9 0.0%, SSX1 exon4 33.0%, SSX3 exon3 27.7%, SSX7 exon7 64.7%, TRIM49 exon5 6.4%, UGT2B15 exon3 89.7%, USP18 exon11 0.0%, WASHC2A exon7 13.6%, WASHC2C exon7 0.0%, XG exon10 87.6%. Nine of them have an owner (CEL, FCGR1A, GYPA, HYDIN, OTOA and PMS2 for Parascopy, HP and RHD for the depth engine, UGT2B15 for pypgx); the other 62 join the genes with no owner that were not read. 27 of the 71 have an exon with no base at 10x: a whole coding exon no short read of mapping quality 20 placed.

## MultiQC Integration

mosdepth output is automatically detected by MultiQC. The `mosdepth.global.dist.txt` and `mosdepth.summary.txt` files are consumed to generate coverage distribution plots and summary statistics in the aggregated report.

## Notes

- The genome-wide pass runs with `--fast-mode`: mosdepth reads only each read's start and end, so it neither corrects for overlapping mates nor reads the CIGAR. That is faster and close enough for the depth summaries. `--no-per-base`, not `--fast-mode`, is what leaves out the large per-base BED. The gene callability pass runs without `--fast-mode` (see [Gene callability](#gene-callability))
- Supports `ALIGN_DIR` variable to use alternative alignments: `ALIGN_DIR=aligned_bwamem2 ./scripts/16b-mosdepth.sh sample`
- Supports `CAPTURE_BED` for WES on-target coverage: `CAPTURE_BED=${GENOME_DIR}/captures/my_panel.bed ./scripts/16b-mosdepth.sh sample`. **The BED file must be inside `GENOME_DIR`** — the Docker container only mounts `GENOME_DIR` as `/genome/`, so paths outside it are invisible to the container
- mosdepth uses 1 main thread + N decompression threads, so `--threads 4` actually uses 5 threads total
