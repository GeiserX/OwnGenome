# Step 3: Variant Calling (BAM to VCF and gVCF)

## What This Does
Identifies all positions where the sample's DNA differs from the reference genome: SNPs (single nucleotide changes) and small indels (insertions/deletions <50bp).

It writes two files:

- `vcf/${SAMPLE}.vcf.gz`: the variant sites only.
- `vcf/${SAMPLE}.g.vcf.gz` (a gVCF): the same calls plus "reference blocks", stretches where the sample matches the reference with a genotype quality. A position inside a block is a confirmed 0/0. A position outside any covered block was not covered. The variant-only VCF cannot tell those two apart. PharmCAT (step 7), the PRS (step 25) and imputation prep with panel sites (step 14) read hom-ref genotypes from the gVCF.

## Why
The VCF file is the foundation for ALL downstream analyses: ClinVar screening, pharmacogenomics, PRS, ROH, etc.

## Tool
- **DeepVariant** v1.10.0 — Google's deep learning variant caller (state-of-the-art accuracy)

## Docker Image
- `DEEPVARIANT_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Prerequisites
- Sorted, indexed BAM file
- GRCh38 reference genome (FASTA + FAI)

## Command
```bash
./scripts/03-deepvariant.sh your_sample male    # or female; see "Sex chromosomes" below
```

What the script runs:

```bash
source versions.env   # from the repository root
REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta   # see 00-reference-setup.md#the-reference-path-on-every-page
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data
THREADS=8

docker run --rm \
  --cpus ${THREADS} --memory 32g \
  -v ${GENOME_DIR}:/genome \
  -v "$PWD/assets/par_grch38.bed:/pgp/par_grch38.bed:ro" \
  "${DEEPVARIANT_IMAGE}" \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref="/genome/${REF_FASTA}" \
    --reads=/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
    --output_vcf=/genome/${SAMPLE}/vcf/${SAMPLE}.part.vcf.gz \
    --output_gvcf=/genome/${SAMPLE}/vcf/${SAMPLE}.part.g.vcf.gz \
    --intermediate_results_dir=/genome/${SAMPLE}/vcf/deepvariant_tmp \
    --sample_name="${SAMPLE}" \
    --num_shards=${THREADS} \
    --haploid_contigs=chrX,chrY \
    --par_regions_bed=/pgp/par_grch38.bed   # these two lines for a male sample only

# The script then renames the .part files (and their .tbi indexes) to
# ${SAMPLE}.vcf.gz and ${SAMPLE}.g.vcf.gz and removes deepvariant_tmp/.

# For WES data, use MODEL_TYPE=WES:
# MODEL_TYPE=WES ./scripts/03-deepvariant.sh your_sample male

# To call only some regions, set INTERVALS (space-separated, passed to --regions):
# INTERVALS="chr20:10000001-10500000" ./scripts/03-deepvariant.sh your_sample male

# A second BAM (BWA-MEM2, long reads) into its own directory, so vcf/ is kept:
# ALIGN_DIR=aligned_bwamem2 VCF_OUT_DIR=vcf_bwamem2 ./scripts/03-deepvariant.sh your_sample male

# Output: ~93MB VCF with ~5.5M total variants (~4.6M PASS), and the gVCF
```

Settings, all optional:

| Variable | Default | What it does |
|---|---|---|
| `THREADS` | 8 | CPUs of the container and DeepVariant's `--num_shards` |
| `DV_MEM` | `32g` | container memory |
| `ALIGN_DIR` | `aligned` | where the BAM is, inside the sample directory |
| `VCF_OUT_DIR` | `vcf` | where the VCF and gVCF go, inside the sample directory |
| `MODEL_TYPE` | `WGS` | `WGS`, `WES`, `PACBIO` or `ONT_R104` |
| `INTERVALS` | whole genome | regions to call |

DeepVariant writes its outputs under `.part` names, and the script renames them only when the VCF, the gVCF and both indexes are complete, so a killed run leaves no VCF that looks finished. Its intermediate files go to `deepvariant_tmp/` in the output directory, not into the container's own disk, and are removed at the end of a run.

## Sex chromosomes

A male sample has one X and one Y. Outside the pseudoautosomal regions (PARs, the ends of X and Y that pair with each other), a heterozygous call on chrX or chrY cannot be real. With `male`, DeepVariant calls chrX and chrY haploid outside the PARs listed in `assets/par_grch38.bed` (the GRCh38 PAR1 and PAR2 of chrX and chrY), so those calls are homozygous. The PAR coordinates come from the [Genome Reference Consortium's GRCh38 page](https://www.ncbi.nlm.nih.gov/grc/human): PAR1 is chrX and chrY 10,001-2,781,479; PAR2 is chrX 155,701,383-156,030,895 and chrY 56,887,903-57,217,415 (1-based). The BED holds the same intervals with 0-based starts. It has no header line, because some of its readers take every line as an interval. `female` and no sex call every contig diploid; no sex prints a note saying so. `run-all.sh` passes the sex it was given.

The alternative callers (03a GATK, 03b FreeBayes, 03c Strelka2, 03d Octopus) take the same `male` or `female` second argument, and with `male` call chrX and chrY with ploidy 1 outside the same PARs:

| Step | How it calls a male's chrX and chrY outside the PARs |
|---|---|
| 03a GATK | the chrX and chrY units are cut at the PARs; the parts outside run as units of their own with `--sample-ploidy 1` |
| 03b FreeBayes | the same cut; the parts outside run with `--ploidy 1` |
| 03c Strelka2 | `--ploidy` with `<sample>.ploidy.vcf.gz` in `vcf_strelka2/`: one `CN=1` record per stretch outside the PARs, its column named after the BAM's read group sample (`@RG SM`) |
| 03d Octopus | Octopus sets ploidy per contig, so the parts outside run in a second Octopus run with `--contig-ploidies chrX=1 chrY=1`, and the two VCFs are joined. The first run passes `--contig-ploidies MT=1 chrM=1`, Octopus' default without its `chrY=1`, so chrY's PARs are diploid like chrX's |

`female` calls every contig diploid (Octopus keeps its default: chrY and chrM haploid). With no second argument, the steps take the sex `run-all.sh` was given for the sample, which it writes as `declared_sex` in `logs/run_status.tsv` before any step runs, and print a note saying so. `run-all.sh` starts the `EXTRA_CALLERS` steps with the sample name only, so that is how they get the sex. With no sex in either place they call as before and print a note. A male run records `sex=male` in the step's run record, so a result called with the other ploidy is called again; female and no sex call the same way and share a record. None of these steps has a Nextflow module: they run only as these scripts. Clair3 (03e, long reads) takes the same `male` or `female` argument (`--gender`).

## Alternative callers: scatter

GATK HaplotypeCaller (03a) multithreads only its PairHMM and FreeBayes (03b) runs on one thread, so each would leave most of the CPUs idle for a whole genome. Both split the work into units: chr1-22, chrX, chrY and chrM one each and every other contig together, or each region of `INTERVALS` (space-separated, for example `INTERVALS="chr20 chr22"`). `SCATTER_JOBS` units (default THREADS/2, fewer when their 8 GB each do not fit in the memory the container engine reports, with a note that says so) run at once, GATK with 2 CPUs and 8 GB each, FreeBayes with 1 CPU and 8 GB, and the parts are joined in reference order. `SCATTER=false` runs one process over everything (32 GB); for a male sample, two in turn, the haploid parts of chrX and chrY and the rest. On chr20 and chr22 of the e2e fixture the scattered and the single run of each caller give the same records (`tests/e2e/sv-mito-telomere-steps-6-scatter.sh`).

Each region of `INTERVALS` must lie on its contig: a start of 1 or more, an end no larger than the contig's length in the reference's `.fai`, and the start not after the end. Two regions must not overlap, because the bases they share would be called twice and joined twice. Regions that only touch, such as `chr20:1-1000 chr20:1001-2000`, are fine. Any of these mistakes stops the step before a caller starts. Step 29 (Mutect2) scatters `INTERVALS=genome` the same way.

## Alternative callers: reruns

`run-all.sh` runs the callers of `EXTRA_CALLERS` (03a-03d) and step 29 on every invocation, so each keeps a finished result only when it was called the way this run would call it. Each step writes a run record beside its output: `<sample>.run` in `vcf_gatk/`, `vcf_freebayes/`, `vcf_strelka2/` and `vcf_octopus/`, and `<sample>_somatic_filtered.run` in `somatic/`. The record holds `INTERVALS` (not for 03c, which calls the whole genome; `SCATTER` too for 03a and 03b; `sex=male` for a male sample in 03a-03d; the optional resources for 29), the BAM's path, the reference and the image, then the BAM's size, modification time and the sha256 of its first and last MiB. A realigned or replaced BAM, even under the same name, size and time, is called again; so are other `INTERVALS`, and a result called with the other ploidy (male, or female and no sex). A result written before this change is called once more on the next run: 03a, 03b and 29 records changed form, and 03c and 03d had none. Octopus also needs the VCF's `.tbi`: it indexes the VCF itself when Octopus did not, and a VCF without its index is called again. For Strelka2, a run of another BAM, finished or not, is removed and configured again. Delete the output folder to call again by hand.


## Resource Requirements
- CPU: `THREADS` (default 8) sets `--cpus` and `--num_shards`. `make_examples` runs one process per shard; `call_variants` on CPU scales sub-linearly ([upstream](https://github.com/google/deepvariant/blob/r1.10/docs/deepvariant-details.md#call_variants)). `run-all.sh` runs DeepVariant with at most 8 shards, fewer when `--max_cpus` (`THREADS`) is lower; [CPU requirements](hardware-requirements.md#cpu-requirements) shows how to give it more
- RAM: 32GB recommended (`DV_MEM`)
- Disk: the intermediate files in `deepvariant_tmp/` need free space in the sample directory while the step runs
- GPU: neither entry point has a GPU option; only `call_variants` could use one (see [troubleshooting](troubleshooting.md#step-3-deepvariant-and-a-gpu))
- Time: see [Hardware and storage requirements](hardware-requirements.md#runtime-per-step)

## Output Interpretation
- **PASS** variants: high-confidence calls (~4.6M per 30X WGS sample)
- **RefCall**: site looks like reference (not a variant)
- **LowQual**: low confidence
- GQ (Genotype Quality): higher = more confident
- DP (Read Depth): typical 25-35x for a 30X WGS sample

## Notes
- DeepVariant is the gold standard for SNPs/indels but does NOT detect:
  - Structural variants >50bp (use Manta, step 4)
  - Repeat expansions (use ExpansionHunter, step 9)
  - Copy number variants (use CNVpytor or GATK gCNV)
- bcftools can also call variants but is significantly less accurate than DeepVariant
- **The small model.** DeepVariant 1.8.0 added a small model, on by default in `run_deepvariant`, that calls the sites it finds easy and passes the rest to the full model. DeepVariant's release notes give a speed-up of about 1.7x (40% less runtime) for WGS with similar or better accuracy, and the flag `--disable_small_model=true` to turn it off. This pipeline runs DeepVariant 1.10.0 with the default, so the small model is on. On one real ~30x genome, a site whose alt reads came from contaminating cDNA was `RefCall` under DeepVariant 1.6 and `PASS` at a modest GQ under 1.10, with the same allele depths. We have not measured whether the small model caused that change, and the fixture has no such site. Turning the small model off would cost about the speed-up above, so the default stays. When such a site is a ClinVar hit, step 06's [cDNA read check](06-clinvar-screen.md#the-cdna-read-check) can flag it from the reads, whatever made it PASS
