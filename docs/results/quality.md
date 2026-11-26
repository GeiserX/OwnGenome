# Data Quality and Reports

Part of [Interpreting Your Results](../interpreting-results.md). Read [Before You Look](../interpreting-results.md#before-you-look-what-the-pipeline-can-tell-you) first.

## Reads, Alignment, Variant Calling and Coverage (Steps 1b, 2, 3, 16, 16b)

These steps check that the data is good enough for everything else. Look at them first if a later result looks strange.

- **fastp (step 1b):** `${SAMPLE}/fastq_trimmed/${SAMPLE}_fastp.html`. A few percent of reads trimmed or dropped is normal; a large share points at a library or upload problem.
- **Alignment and variant calling (steps 2 and 3):** the [example output](../interpreting-results.md#variant-calling-step-3) shows what a 30X genome looks like: about 4.5 to 5.5 million variants and a Ti/Tv ratio of 2.0 to 2.1.
- **Coverage (step 16b, mosdepth):** `${SAMPLE}/mosdepth/${SAMPLE}.mosdepth.summary.txt`. The `total` row should show a mean near the depth you paid for (about 30). Below 15, small-variant calls lose accuracy. See [step 16b](../16b-mosdepth.md#interpreting-results).
- **Sex check (step 16, indexcov):** `${SAMPLE}/indexcov/` estimates the copy number of chrX and chrY from the BAM index. About 1 and 1 is XY, about 2 and 0 is XX. A result that does not match the sex you gave the pipeline means a sample swap, a mislabelled file or a real sex-chromosome difference; [step 16](../16-indexcov.md#interpretation) lists the patterns. Re-check the input before reading anything else.

## What Part of the Genome a Run Covered

Some results are a share, a sum or a placement over the whole genome: callability, the polygenic scores, fine-grained ancestry, telomere content, runs of homozygosity and archaic ancestry. A run on part of the genome (`--intervals`, or a VCF or BAM cut down to a few chromosomes) can also produce these results. When the summary records a scope mark for a result, both reports show it in the result section and list it in the run header. The live site does not show scope marks yet.

So `summary.json` records the run's scope under `scope`. It comes from the chromosomes that the records of the sample VCF are on (`chr20` and `20` count alike):

- `autosomes_covered`: how many of chr1 to chr22 have calls, and `autosomes`: which ones.
- `chrx`: whether chrX has calls. chrX alone does not make a run partial, because these results read the autosomes.
- `intervals`: the intervals the run was limited to, from the `run intervals` row of `run_manifest.tsv`, else `null`. `run-all.sh` writes that row when it was given `INTERVALS` (its `--intervals`), and so does the Nextflow report for `--intervals`; a whole-genome run has no such row. So a run limited to parts of all 22 autosomes is marked too.
- `partial` and `mark`: the run is partial when it covers fewer than 22 autosomes (`partial: covers 2 of 22 autosomes`) or was limited to intervals (`partial: limited to the run's intervals (...)`).

When the run is partial, each of those six sections that has a result of this run gets the same text in its `partial` field. A stale section holds a result from an earlier run, and that run's scope was not recorded: it gets `partial: not known, from an earlier run`, whatever this run covered. Without a sample VCF that can be read to its end, the scope is not known: `summary.json` has no `scope`, and no section is marked. The scope says which chromosomes have calls. It does not say how much of each chromosome was read well; callability does that.

## Sample Identity and Contamination (Step 33)

The Quality Control block in the text and HTML reports and the live QC card show when step 33 calls sample identity or the FREEMIX estimate uncertain. The reason appears when step 33 supplies one. A similarity signal under uncertain identity appears as a possible match.

## Reports (Steps 24 and 28)

- **HTML report (step 24):** `${SAMPLE}/${SAMPLE}_report.html` collects the headline results of the other steps on one page, including the ClinVar hits table. It contains health findings: share it only as you would share a medical record.
- **Selected steps with no result:** When the summary identifies one, the text and HTML report headers show the step and the reason recorded in the summary.
- **MultiQC (step 28):** `${SAMPLE}/multiqc/multiqc_report.html` puts the QC of fastp, samtools, mosdepth and the other tools in one page. It is about data quality, not about your health.
