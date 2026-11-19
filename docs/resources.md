# Recommended Resources

Free resources for understanding your genomic data, organized from most accessible to most technical.

---

## Start Here (No Background Needed)

### Understanding Genetics Basics

- **[Your Genome](https://www.yourgenome.org/)** — Wellcome Connecting Science. Excellent visual explainers on DNA, genes, inheritance, and genetic testing. Start here if you are new to genetics.
- **[MedlinePlus Genetics](https://medlineplus.gov/genetics/)** — NIH's plain-language genetics reference. Look up any gene or condition mentioned in your results.
- **[OMIM (Online Mendelian Inheritance in Man)](https://www.omim.org/)** — The definitive catalog of human genes and genetic disorders. Search for any gene from your ClinVar report to understand the associated conditions.

### Understanding Your Specific Results

- **[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)** — Look up any variant ID (rsXXXXX) from your ClinVar screen to see its clinical significance, supporting evidence, and review status. Pay attention to the star rating and the number of submitters.
- **[gnomAD Browser](https://gnomad.broadinstitute.org/)** — Look up any variant to see how common it is in the general population. If a "pathogenic" variant has a gnomAD AF > 1%, the classification may be wrong.
- **[PharmGKB](https://www.pharmgkb.org/)** — The pharmacogenomics knowledge base. Look up any gene from your PharmCAT report to understand which drugs are affected and what the clinical recommendations are.
- **[CPIC Guidelines](https://cpicpgx.org/guidelines/)** — Clinical Pharmacogenetics Implementation Consortium. The clinical guidelines that PharmCAT's recommendations are based on. Organized by gene-drug pair.

---

## Intermediate (Some Biology Background Helpful)

### Variant Interpretation

- **[ClinGen](https://www.clinicalgenome.org/)** — The Clinical Genome Resource. Maintains gene-disease validity assessments and variant interpretation standards. Useful for understanding how variants are classified.
- **[ACMG SF v3.3 Gene List](https://doi.org/10.1016/j.gim.2025.101454)** — The 84 genes that ACMG recommends reporting incidental findings for (Lee et al., Genet Med 2025;27(8):101454, doi:10.1016/j.gim.2025.101454). This is the list the pipeline's secondary-findings section applies (`assets/acmg_sf.tsv`). Note: CPSR uses its own cancer-focused panels (Genomics England PanelApp), not the ACMG SF list — CPSR covers 500+ cancer genes but omits non-cancer ACMG genes (cardiac, metabolic). A negative CPSR result does not constitute a complete ACMG incidental-findings screen.
- **[MitoMap](https://www.mitomap.org/)** — The human mitochondrial genome database, with the curated disease status of mtDNA variants. Step 20 already puts ClinVar's record and gnomAD v3.1's frequency beside each call; MITOMAP is a manual lookup in a browser, as it answers scripted requests with a 403.
- **[GeneReviews](https://www.ncbi.nlm.nih.gov/books/NBK1116/)** — Expert-authored, peer-reviewed descriptions of genetic conditions. If you find a pathogenic variant, the GeneReviews entry for the associated condition is the best single resource for understanding it.

### Structural Variants and CNVs

- **[DECIPHER](https://www.deciphergenomics.org/)** — Database of genomic variation and associated phenotype. Useful for checking if a CNV or SV overlaps with known pathogenic regions.
- **[DGV (Database of Genomic Variants)](http://dgv.tcag.ca/)** — Catalog of structural variants found in healthy populations. If your SV is in DGV, it is likely benign.
- **[ClinGen Dosage Sensitivity](https://dosage.clinicalgenome.org/)** — Which genes are sensitive to having too few copies (haploinsufficiency) or too many copies (triplosensitivity). Essential for interpreting CNV results.

---

## Advanced (For Researchers and Power Users)

### Tools and Methods

- **[GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)** — The gold standard for variant calling pipelines. Useful for understanding the theory behind this pipeline's approach.
- **[DeepVariant paper (2018)](https://www.nature.com/articles/nbt.4235)** — The original publication describing the deep learning variant caller used in step 3.
- **[Manta paper (2016)](https://academic.oup.com/bioinformatics/article/32/8/1220/1743909)** — How Manta detects structural variants from paired-end reads.
- **[VEP documentation](https://www.ensembl.org/info/docs/tools/vep/index.html)** — Complete reference for understanding VEP's annotation fields and consequence types.

### Population Genetics

- **[1000 Genomes Project](https://www.internationalgenome.org/)** — Catalog of human genetic variation from 2,504 individuals across 26 populations.
- **[gnomAD paper (2020)](https://www.nature.com/articles/s41586-020-2308-7)** — Understanding allele frequencies and how gnomAD's constraint metrics (pLI, LOEUF) predict gene essentiality.

### Cohort data terms

`setup.sh` downloads or builds each population panel on your own machine, from its source. The repository ships none of them. These are the terms of each source:

- **1000 Genomes** (the phase 3 sites of step 43, the high-coverage panel of step 51 (MoChA), and pgsc_calc's panel of steps 25 and 26): available without embargo. Cite the project (The 1000 Genomes Project Consortium, Nature 2015, doi:10.1038/nature15393; https://www.internationalgenome.org/IGSR_disclaimer). For the high-coverage data, cite Byrska-Bishop et al., Cell 2022 (doi:10.1016/j.cell.2022.08.004) and add: "These data were generated at the New York Genome Center with funds provided by NHGRI Grant 3UM1HG008901-03S1." Its cell lines come from the NIGMS Human Genetic Cell Repository at the Coriell Institute for Medical Research.
- **gnomAD HGDP + 1000 Genomes** (step 26b): released under gnomAD's terms, CC0 1.0 (https://gnomad.broadinstitute.org/terms).
- **HPRC release 2** (step 53, Locityper): public domain (https://humanpangenome.org/data-use/).
- **HAPNEST** (the PGS Catalog's small synthetic panel the tests use; `ANCESTRY_PANEL_NAME` can pick it for steps 25 and 26): synthetic genotypes made from 1000 Genomes and HGDP, published by EMBL-EBI (BioStudies S-BSST936, https://www.ebi.ac.uk/biostudies/studies/S-BSST936) under EMBL-EBI's terms of use (https://www.ebi.ac.uk/about/terms-of-use/).

The test-data prerelease `e2e-data-1` (https://github.com/GeiserX/OwnGenome/releases/tag/e2e-data-1) mirrors some of these files and the table builders' downloads unchanged, for CI only. Its release notes name each file's source and terms.

### Online Courses (Free)

- **[Coursera: Genomic Data Science](https://www.coursera.org/specializations/genomic-data-science)** — Johns Hopkins. 8-course specialization covering the full analysis pipeline from FASTQ to interpretation.
- **[MIT OpenCourseWare: Computational Biology](https://ocw.mit.edu/courses/6-047-computational-biology-fall-2015/)** — MIT 6.047. More theoretical but excellent for understanding algorithms behind variant calling and annotation.

---

## Useful Command-Line References

### bcftools Cheat Sheet

```bash
# Count variants by type
bcftools stats file.vcf.gz | grep "^SN"

# Extract specific fields
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' file.vcf.gz

# Filter to PASS variants only
bcftools view -f PASS file.vcf.gz

# Filter to a specific region
bcftools view -r chr22:16000000-17000000 file.vcf.gz

# Filter to heterozygous variants only
bcftools view -g het file.vcf.gz

# Filter to homozygous ALT variants only
bcftools view -g hom file.vcf.gz

# Count variants per chromosome
bcftools view -f PASS file.vcf.gz | grep -v '^#' | cut -f1 | sort | uniq -c | sort -rn
```

### samtools Cheat Sheet

```bash
# Quick alignment summary
samtools flagstat file.bam

# Per-chromosome read counts
samtools idxstats file.bam

# Extract reads from a specific region
samtools view -b file.bam chr22:16000000-17000000 > region.bam

# View BAM header (reference info, read groups)
samtools view -H file.bam

# Calculate average depth
samtools depth -a file.bam | awk '{sum+=$3; n++} END {print sum/n}'
```

---

## Community

- **[r/bioinformatics](https://www.reddit.com/r/bioinformatics/)** — General bioinformatics discussion
- **[r/genomics](https://www.reddit.com/r/genomics/)** — Genomics-focused community
- **[BioStars](https://www.biostars.org/)** — Q&A forum for bioinformatics (like StackOverflow for genomics)
- **[SEQanswers](http://seqanswers.com/)** — Sequencing technology discussion forum
