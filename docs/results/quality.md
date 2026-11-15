# Data Quality and Reports

Part of [Interpreting Your Results](../interpreting-results.md). Read [Before You Look](../interpreting-results.md#before-you-look-what-the-pipeline-can-tell-you) first.

## Reads, Alignment, Variant Calling and Coverage (Steps 1b, 2, 3, 16, 16b)

These steps check that the data is good enough for everything else. Look at them first if a later result looks strange.

- **fastp (step 1b):** `${SAMPLE}/fastq_trimmed/${SAMPLE}_fastp.html`. A few percent of reads trimmed or dropped is normal; a large share points at a library or upload problem.
- **Alignment and variant calling (steps 2 and 3):** the [example output](../interpreting-results.md#variant-calling-step-3) shows what a 30X genome looks like: about 4.5 to 5.5 million variants and a Ti/Tv ratio of 2.0 to 2.1.
- **Coverage (step 16b, mosdepth):** `${SAMPLE}/mosdepth/${SAMPLE}.mosdepth.summary.txt`. The `total` row should show a mean near the depth you paid for (about 30). Below 15, small-variant calls lose accuracy. See [step 16b](../16b-mosdepth.md#interpreting-results).
- **Sex check (step 16, indexcov):** `${SAMPLE}/indexcov/` estimates the copy number of chrX and chrY from the BAM index. About 1 and 1 is XY, about 2 and 0 is XX. A result that does not match the sex you gave the pipeline means a sample swap, a mislabelled file or a real sex-chromosome difference; [step 16](../16-indexcov.md#interpretation) lists the patterns. Re-check the input before reading anything else.

## Reports (Steps 24 and 28)

- **HTML report (step 24):** `${SAMPLE}/${SAMPLE}_report.html` collects the headline results of the other steps on one page, including the ClinVar hits table. It contains health findings: share it only as you would share a medical record.
- **MultiQC (step 28):** `${SAMPLE}/multiqc/multiqc_report.html` puts the QC of fastp, samtools, mosdepth and the other tools in one page. It is about data quality, not about your health.
