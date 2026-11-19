# Long-Read Sequencing Guide

## What Is Long-Read Sequencing?

Standard WGS (from Illumina, Novogene, Dante Labs, etc.) produces **short reads** — millions of DNA fragments typically 150 base pairs long. Two companies offer a fundamentally different approach: reading much longer stretches of DNA in a single pass.

### Oxford Nanopore Technology (ONT)

ONT sequencers (MinION, PromethION, GridION) push a strand of DNA through a tiny protein pore. As each nucleotide passes through, it changes the electrical current. The device reads these current changes in real time to determine the sequence.

- **Read lengths:** 1,000 to 100,000+ bp (some reads exceed 1 million bp)
- **Accuracy:** ~99% with R10.4.1 chemistry and super-accuracy basecalling (Q20+)
- **Common data format:** FASTQ (basecalled from FAST5/POD5 raw signal)
- **Vendors/services:** ONT direct (MinION, PromethION), some clinical labs

### PacBio HiFi (High-Fidelity)

PacBio sequencers (Sequel II, Revio) use a different approach: a polymerase enzyme copies a circular DNA template while a camera records fluorescent signals from each incorporated nucleotide. By reading the same circle multiple times, errors cancel out.

- **Read lengths:** 10,000-25,000 bp
- **Accuracy:** >99.9% (Q30+, comparable to Illumina)
- **Common data format:** Unaligned BAM or FASTQ (HiFi/CCS reads)
- **Vendors/services:** PacBio direct, Revio service providers, some clinical labs

### Why Use Long Reads?

Long reads solve problems that short reads cannot:

| Capability | Short Reads (150bp) | Long Reads (10-100Kbp) |
|---|---|---|
| SNPs and small indels | Excellent | Excellent (HiFi) / Good (ONT) |
| Structural variants (>50bp) | Limited | Excellent |
| Repeat expansions | Indirect (ExpansionHunter) | Direct measurement |
| Phasing (which allele is on which chromosome) | Statistical | Read-based, much more accurate |
| Segmental duplications / tandem repeats | Mostly invisible | Resolvable |
| Methylation | Requires bisulfite conversion | Native (ONT) |

---

## Pipeline Compatibility

**What the long-read branch covers:** alignment (script 2b), small variants with Clair3 (3e) and structural variants with Sniffles2 (4c). After that, the VCF-based steps below work on the result. Nothing else in this pipeline is long-read specific: no repeat genotyping, no paralog or CYP2D6 resolution, no methylation and no phasing. For those, use the vendor workflows, which run the long-read tools named on this page: PacBio's [HiFi-human-WGS-WDL](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL) and Oxford Nanopore's [wf-human-variation](https://github.com/epi2me-labs/wf-human-variation).

Not every step in this pipeline works with long-read data. Here is the full breakdown.

### Works As-Is (No Changes Needed)

These steps take a VCF or BAM and work identically regardless of read technology:

| Step | Tool | Notes |
|---|---|---|
| 5 | AnnotSV | Takes any SV VCF — point it at Sniffles2 output |
| 6 | ClinVar Screen | VCF-only, technology-independent |
| 7 | PharmCAT | VCF-only — requires VCF + `.tbi` symlinked to `vcf/` (see below) |
| 11 | ROH Analysis | VCF-only via plink2 |
| 12 | Mito Haplogroup | VCF-only via haplogrep3 |
| 13 | VEP Annotation | VCF-only |
| 14 | Imputation Prep | VCF-only |
| 17 | CPSR | VCF-only |
| 22 | SURVIVOR Merge | Takes any SV VCF set |
| 23 | Clinical Filter | VCF-only |
| 24 | HTML Report | Aggregates existing outputs |
| 25 | PRS | VCF-only via pgsc_calc |
| 26 | Ancestry | VCF-only via pgsc_calc (with the ancestry panel) |
| 27 | CPIC Lookup | Reads PharmCAT JSON output |

### Needs Long-Read Specific Scripts

| Step | Short-Read Tool | Long-Read Alternative | Script |
|---|---|---|---|
| 2 | minimap2 (`-x sr`) | minimap2 (`-x map-ont` / `-x map-hifi`) | `scripts/02b-alignment-longread.sh` |
| 3 | DeepVariant (WGS model) | **Clair3** or DeepVariant (ONT/PACBIO model) | `scripts/03e-clair3.sh` |
| 4 | Manta | **Sniffles2** | `scripts/04c-sniffles2.sh` |

### Does NOT Work with Long Reads

These tools are specifically designed for short-read data and will produce incorrect results or crash with long-read input:

| Step | Tool | Why It Fails | Alternative |
|---|---|---|---|
| 1b | fastp | Expects paired-end Illumina reads; no long-read mode | NanoPlot, LongQC, or chopper for long-read QC |
| 2a | BWA-MEM2 | Short-read aligner only; cannot handle reads >500bp | Use minimap2 with long-read preset |
| 3a | GATK HaplotypeCaller | Active-region model assumes short reads; poor with long reads | Clair3 or DeepVariant long-read mode |
| 3b | FreeBayes | Bayesian model tuned for short reads | Clair3 or DeepVariant long-read mode |
| 3c | Strelka2 | Scoring model trained on BWA-MEM short-read data | Clair3 or DeepVariant long-read mode |
| 3d | Octopus | Haplotype model designed for short reads | Clair3 or DeepVariant long-read mode |
| 4 | Manta | Illumina-specific insert size model | Sniffles2 |
| 4a | TIDDIT | Short-read coverage/discordance model | Sniffles2 |
| 4b | GRIDSS | Assembly-based, short-read specific | Sniffles2 |
| 9 | ExpansionHunter | Illumina short-read graph model; expects paired-end data | TRGT (PacBio) or Straglr (ONT), outside this pipeline (vendor workflows) |
| 15 | duphold | Adds depth fold-change tags (DHFFC, DHBFC) to Manta's calls and filters on them; reads Manta's VCF and the short-read BAM only | Not needed: Sniffles2 writes its own QUAL and GQ |
| 18 | CNVpytor | Read-depth model calibrated for short reads | Sniffles2 detects CNVs natively |
| 19 | Delly | Paired-end and split-read model | Sniffles2 |
| 21 | Cyrius (CYP2D6) | Short-read depth-based star allele caller | Paraphase (long-read CYP2D6 resolver), outside this pipeline (vendor workflows) |

---

## Step-by-Step: Long-Read Analysis

### Prerequisites

Same reference data as the short-read pipeline:

```bash
export GENOME_DIR=/path/to/your/data
# Reference genome and index should already exist from initial setup
# See docs/00-reference-setup.md
```

### Step 1: Prepare Your Data

Place your long-read data in the sample directory:

```bash
# Create sample directory
mkdir -p ${GENOME_DIR}/your_sample/fastq/

# Copy or symlink your long-read FASTQ
# ONT: typically one FASTQ per run (basecalled from FAST5/POD5)
cp /path/to/ont_reads.fastq.gz ${GENOME_DIR}/your_sample/fastq/your_sample.fastq.gz

# PacBio: typically one FASTQ or unaligned BAM per cell
cp /path/to/hifi_reads.fastq.gz ${GENOME_DIR}/your_sample/fastq/your_sample.fastq.gz
# OR for unaligned BAM:
cp /path/to/hifi_reads.bam ${GENOME_DIR}/your_sample/fastq/your_sample.bam
```

Note: long-read data is **single-end** (no R1/R2 pairs). The alignment script expects a single file, not paired files.

### Step 2: Align Reads

```bash
# For Oxford Nanopore
PLATFORM=ont ./scripts/02b-alignment-longread.sh your_sample

# For PacBio HiFi
PLATFORM=hifi ./scripts/02b-alignment-longread.sh your_sample
```

Output: `${GENOME_DIR}/your_sample/aligned_longread/your_sample_sorted.bam`

The script writes the BAM under a temporary name, indexes it, checks it with `samtools quickcheck` and only then renames it, so a killed run never puts a partial BAM under `your_sample_sorted.bam` for Clair3 or Sniffles2 to read. A BAM and index from an earlier complete run stay in place, and those steps still read them. `samtools sort` uses `THREADS` threads (default 8) with 1 GB each, and spills to `aligned_longread/your_sample.sort_tmp/`, which is removed at the end.

The script auto-detects input files in this order:
1. `fastq/your_sample.fastq.gz`
2. `fastq/your_sample_lr.fastq.gz`
3. `fastq/your_sample.fq.gz`
4. `fastq/your_sample.bam` (unaligned BAM)

Override with: `INPUT=${GENOME_DIR}/your_sample/fastq/custom_reads.fastq.gz PLATFORM=ont ./scripts/02b-alignment-longread.sh your_sample`

For ONT reads the script also records the basecaller model in the aligned BAM's read group (`DS:basecall_model=...`), where step 3e reads it to pick the Clair3 model. It takes `BASECALL_MODEL` when set, else the `basecall_model=` that Dorado writes into an unaligned BAM's read group. A FASTQ does not carry it, so set it for FASTQ input:

```bash
BASECALL_MODEL=dna_r10.4.1_e8.2_400bps_sup@v5.0.0 PLATFORM=ont ./scripts/02b-alignment-longread.sh your_sample
```

> **Note:** The INPUT path must be **physically** inside `GENOME_DIR` because Docker only mounts that directory. Symlinks that point outside `GENOME_DIR` will not work — Docker cannot follow them. Copy or move the file instead.

### Step 3: Variant Calling (SNPs + Indels)

You have two options:

#### Option A: Clair3 (recommended for long reads)

```bash
# For ONT
PLATFORM=ont ./scripts/03e-clair3.sh your_sample

# For PacBio HiFi
PLATFORM=hifi ./scripts/03e-clair3.sh your_sample
```

Output: `${GENOME_DIR}/your_sample/vcf_clair3/your_sample.vcf.gz`

Clair3 uses deep-learning models trained on long-read error profiles, one per platform and, for ONT, one per basecaller model. HiFi uses `hifi_revio`. For ONT the step picks the model that matches the reads' basecaller (`og longread clair3-model`):

1. `CLAIR3_MODEL`, a model directory, when you set it.
2. Else `BASECALL_MODEL`, a Dorado model name such as `dna_r10.4.1_e8.2_400bps_sup@v5.0.0`.
3. Else the `basecall_model=` in the BAM's read group (step 2b writes it; a BAM aligned by Dorado has it too).

| Basecaller model | Clair3 model |
|---|---|
| `dna_r10.4.1_e8.2_400bps_sup@vX.Y.Z` | `r1041_e82_400bps_sup_vXYZ` (for example v5.0.0 gives `r1041_e82_400bps_sup_v500`) |
| `dna_r10.4.1_e8.2_400bps_hac@vX.Y.Z` | `r1041_e82_400bps_hac_vXYZ` |
| R9.4.1 SUP (Dorado `dna_r9.4.1_e8_sup@v3.3`, Guppy `dna_r9.4.1_450bps_sup`) | `r941_prom_sup_g5014` |
| R9.4.1 HAC | `r941_prom_hac_g360+g422` |
| fast models, other pores or speeds, RNA | none |

The pinned image holds the R10.4.1 SUP models v4.1.0, v5.0.0 and v5.2.0 and the HAC models v4.1.0, v5.0.0, v5.2.0 and v6.0.0. Before Clair3 runs, the step checks that the chosen directory exists in the image, `CLAIR3_MODEL` included. With no basecaller model known, one with no Clair3 model (fast basecalls), or a model the image lacks (for example SUP v4.3.0), the step stops with exit 1 and lists the image's models. It never falls back to another model: a model trained on other basecalls still writes a VCF, but every call is less accurate and nothing in the output says so. Set `CLAIR3_MODEL` to use a model by hand.

#### Option B: DeepVariant (long-read mode)

DeepVariant 1.10.0+ also supports long reads with dedicated models:

```bash
source versions.env   # from the repository root
REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta   # see 00-reference-setup.md#the-reference-path-on-every-page
SAMPLE=your_sample

docker run --rm \
  --cpus 8 --memory 32g \
  -v "${GENOME_DIR}:/genome" \
  "${DEEPVARIANT_IMAGE}" \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=ONT_R104 \
    --ref="/genome/${REF_FASTA}" \
    --reads="/genome/${SAMPLE}/aligned_longread/${SAMPLE}_sorted.bam" \
    --output_vcf="/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
    --num_shards=8

# For PacBio, use --model_type=PACBIO instead of ONT_R104
```

**Clair3 vs DeepVariant for long reads:**

| Aspect | Clair3 | DeepVariant |
|---|---|---|
| ONT accuracy | Excellent (purpose-built) | Good (newer models improving) |
| PacBio HiFi accuracy | Excellent | Excellent |
| Speed | Generally faster | Slower (more compute-heavy) |
| Indel calling | Very good | Very good |
| Community support | Active (HKU) | Active (Google) |
| Model updates | Frequent (matches new chemistries) | Periodic |

For PacBio HiFi data, both tools perform comparably. For ONT data, Clair3 tends to have an edge because its models are updated more frequently to match new ONT basecalling improvements.

### Step 4: Structural Variant Calling

```bash
# Default: uses aligned_longread/ directory
./scripts/04c-sniffles2.sh your_sample

# Or explicitly:
ALIGN_DIR=aligned_longread ./scripts/04c-sniffles2.sh your_sample
```

Output: `${GENOME_DIR}/your_sample/sv_sniffles/your_sample_sv.vcf.gz`

Before Sniffles2 runs, the step reads the first 10,000 primary reads of the BAM with the pinned samtools image. When their mean length is below 1,000 bp, it stops with exit 1 and names the long-read BAM to use. Sniffles2 finds no SV in short reads and still exits 0, so without this check a short-read BAM gave an empty VCF that looked like a finished run. Short-read SVs come from steps 04 (Manta), 04a, 04b and 19.

`setup.sh` installs a tandem-repeat map for Sniffles2: every GRCh38 tandem repeat in GIAB's genome stratifications v3.6 (`reference/GRCh38_AllTandemRepeats.bed`, pinned by sha256 in `TANDEM_REPEATS_URL` and `TANDEM_REPEATS_SHA256`; NIST data, free to reuse). When the map is there, the step passes it as `--tandem-repeats`, which gives SVs inside repeats better breakpoints and merging. When it is not, the step prints a note and runs without it. We do not use the GRCh38 file that ships with Sniffles: it comes from PacBio's pbsv under a licence that allows PacBio data only.

What we measured: the step on GIAB's public HG002 PacBio HiFi BAM (Sequel II, 15-20 kb, 10x), chr20:10-20 Mb, scored with truvari against GIAB's v5.0q SV benchmark inside its benchmark regions (chr20:10.1-19.9 Mb, PASS calls, 1 kb matching distance):

| | Calls on chr20 | Recall | Precision |
|---|---|---|---|
| Without the tandem-repeat map | 62 | 0.77 (49 of 64) | 0.96 |
| With the map | 63 | 0.78 (50 of 64) | 0.96 |

That is one 10 Mb window at low depth, not a genome-wide benchmark. For Sniffles2's own accuracy figures, see Smolka et al., Nature Biotechnology 2024 (https://doi.org/10.1038/s41587-023-02024-y). No CI case runs this step on real long reads; [Testing](testing.md) lists what checks it.

### Step 5: Continue with Standard Pipeline Steps

After alignment and variant calling, most downstream steps work normally. Point them at the long-read outputs:

```bash
SAMPLE=your_sample

# ClinVar screen (uses VCF from Clair3)
VCF_DIR=vcf_clair3 ./scripts/06-clinvar-screen.sh "$SAMPLE"

# PharmCAT (uses VCF from Clair3)
# PharmCAT needs the VCF + .tbi index at the default vcf/ location.
# Symlink or copy both files:
mkdir -p "${GENOME_DIR}/${SAMPLE}/vcf"
ln -sf "../vcf_clair3/${SAMPLE}.vcf.gz" "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
ln -sf "../vcf_clair3/${SAMPLE}.vcf.gz.tbi" "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi"
./scripts/07-pharmacogenomics.sh "$SAMPLE"

# VEP annotation (uses VCF from Clair3)
VCF_DIR=vcf_clair3 ./scripts/13-vep-annotation.sh "$SAMPLE"

# Annotate Sniffles2 SVs
SV_VCF="${GENOME_DIR}/${SAMPLE}/sv_sniffles/${SAMPLE}_sv.vcf.gz" ./scripts/05-annotsv.sh "$SAMPLE"
```

---

## Docker Images Used

| Tool | Image | Size |
|---|---|---|
| minimap2 | `MINIMAP2_IMAGE` | ~30 MB |
| samtools | `SAMTOOLS_IMAGE` | ~200 MB |
| Clair3 | `CLAIR3_IMAGE` | ~3 GB (includes all models) |
| Sniffles2 | `SNIFFLES_IMAGE` | ~200 MB |
| DeepVariant | `DEEPVARIANT_IMAGE` | ~5 GB |

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

Pre-pull images before your first run:

```bash
source versions.env   # from the repository root
docker pull "${MINIMAP2_IMAGE}"
docker pull "${SAMTOOLS_IMAGE}"
docker pull "${CLAIR3_IMAGE}"
docker pull "${SNIFFLES_IMAGE}"
```

---

## ONT-Specific Notes

### Basecalling Matters

ONT raw data (FAST5/POD5) must be basecalled before entering this pipeline. Use the latest Dorado basecaller with the **super accuracy (SUP)** model for best variant calling:

```bash
# This happens BEFORE the pipeline — on a GPU machine
dorado basecaller sup pod5_dir/ > reads.bam
dorado demux reads.bam  # if barcoded
```

Step 3e picks the Clair3 model from the basecaller model (see Step 3 above), so keep the `basecall_model=` that Dorado writes into the BAM, or set `BASECALL_MODEL`. Clair3 has no model for fast basecalls.

### Read Length Distribution

ONT generates a wide range of read lengths. Very short reads (<1kb) add noise. If your N50 is below 5kb, consider filtering:

No step of this pipeline filters reads. If you want to drop reads under 1 kb, run a read filter such as chopper on the FASTQ yourself before script 2b.

### Methylation

ONT natively detects methylation (5mC, 6mA) during basecalling with Dorado. Methylation tags are stored in the BAM as MM/ML tags. This pipeline does not process methylation data. When the input is an unaligned BAM, script 2b carries the MM/ML tags into the aligned BAM, where the vendor workflows above can read them.

---

## PacBio-Specific Notes

### HiFi vs CLR

This pipeline supports **HiFi (CCS)** reads only, not older CLR (continuous long reads). HiFi reads are generated by PacBio Sequel II, Sequel IIe, and Revio systems using circular consensus sequencing.

If you received CLR data, you'll need to generate CCS reads first using `ccs` or `pbccs` (PacBio tool, outside this pipeline).

### Kinetic Tags

PacBio HiFi BAMs may contain kinetic modification tags (fi, fp, ri, rp) for methylation detection. Like ONT methylation, these are preserved in alignment but not processed by this pipeline.

---

## Hybrid Analysis (Short + Long Reads)

If you have both short-read and long-read data for the same sample, you can run both pipelines and combine results:

1. **Short-read pipeline:** standard `02-alignment.sh` + `03-deepvariant.sh` + `04-manta.sh`
2. **Long-read pipeline:** `02b-alignment-longread.sh` + `03e-clair3.sh` + `04c-sniffles2.sh`
3. **Merge variant calls:**
   - SNPs/indels: use `benchmark-variants.sh` to compare, then take the intersection for high-confidence calls
   - SVs: use `22-survivor-merge.sh` with both Manta and Sniffles2 VCFs for consensus SVs

This is the gold standard approach when both data types are available.

---

## Troubleshooting

### "minimap2 runs out of memory"

Long-read alignment uses more memory than short-read because the index is loaded into RAM. Increase Docker memory limit:

```bash
# In the script, change --memory 16g to --memory 24g or higher
```

### "names no basecaller model", "has no model" or "is not in the image"

Step 3e stopped before Clair3 ran, because it could not match a model to the reads' basecaller. The message lists the models in `CLAIR3_IMAGE` (`/opt/models/`).

1. No basecaller model: set `BASECALL_MODEL` to the Dorado model the reads were basecalled with. Dorado's BAM names it in its `@RG` line (`samtools view -H reads.bam | grep '^@RG'`), and its FASTQ headers carry `basecall_model_version_id=`.
2. No Clair3 model or not in the image: set `CLAIR3_MODEL` to the listed model closest to your basecaller, for example `CLAIR3_MODEL=/opt/models/r1041_e82_400bps_sup_v520`, or use a newer `CLAIR3_IMAGE` that ships the matching model.

### "Sniffles2 needs long reads" or "Sniffles2 produces 0 SVs"

- Step 04c stops when the first 10,000 primary reads of the BAM average below 1,000 bp. Point `ALIGN_DIR` at the long-read BAM (step 02b writes `aligned_longread/<sample>_sorted.bam`).
- A mean above that but still low (below about 5,000 bp for HiFi) usually means a filtered or fragmented run: check `samtools stats your.bam | grep "average length"`
- Check that the BAM is sorted and indexed

### "DeepVariant long-read model fails"

- Use `--model_type=ONT_R104` for ONT (not `ONT` alone)
- Use `--model_type=PACBIO` for PacBio HiFi
- DeepVariant 1.10.0 requires the exact model type string — check the DeepVariant release notes for your version

### "VCF downstream steps fail on Clair3 output"

Steps 6 (ClinVar), 13 (VEP), and 17 (CPSR) accept a `VCF_DIR` env var:
```bash
VCF_DIR=vcf_clair3 ./scripts/06-clinvar-screen.sh your_sample
VCF_DIR=vcf_clair3 ./scripts/13-vep-annotation.sh your_sample
```
For steps that don't support `VCF_DIR` (e.g., PharmCAT), symlink the Clair3 VCF:
```bash
ln -s ../vcf_clair3/${SAMPLE}.vcf.gz ${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz
ln -s ../vcf_clair3/${SAMPLE}.vcf.gz.tbi ${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi
```

---

## Expected Output Comparison

Rough numbers for a 30X human genome. Each runtime names its source or says it was not measured:

| Metric | Short-Read (Illumina) | ONT (R10.4.1 SUP) | PacBio HiFi (Revio) |
|---|---|---|---|
| SNPs | ~4.5M | ~4.3M | ~4.5M |
| Indels | ~1.0M | ~0.8M | ~1.0M |
| PASS rate | ~93% | ~85-90% | ~93% |
| SVs (total) | ~7-9K (Manta) | ~20-25K (Sniffles2) | ~18-22K (Sniffles2) |
| SV insertions | ~500 (poor) | ~8-10K (excellent) | ~7-9K (excellent) |
| Alignment time | minimap2: 3 h 51 min - 5 h 21 min with duplicate marking, measured on 30-40x ([runtime table](hardware-requirements.md#runtime-per-step)) | not measured here | not measured here |
| Variant calling time | DeepVariant: 8 h 50 min - 13 h 30 min with 8 shards, three measured 30-40x runs ([runtime table](hardware-requirements.md#runtime-per-step)) | Clair3: not measured here | Clair3: not measured here |

The biggest advantage of long reads is SV detection — especially insertions, which short reads systematically miss.
