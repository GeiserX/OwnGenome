# Whole Exome Sequencing (WES) Guide

If your data comes from whole exome sequencing instead of whole genome sequencing, most of this pipeline still works -- but some steps need different settings, and others should be skipped entirely. This guide explains what WES is, which pipeline steps apply, and how to run them: the Nextflow pipeline with `--model_type WES` sets them for you ([Run the Nextflow Pipeline](#run-the-nextflow-pipeline)), and the numbered scripts take them one by one.

---

## What WES Is (and How It Differs from WGS)

Whole exome sequencing targets the **exons** of your genome -- the ~22,000 protein-coding gene regions that make up roughly 2% of total DNA. A capture kit (set of probes) binds to these regions, pulling them out of a DNA sample for sequencing. Everything else is discarded before sequencing begins.

| | WGS (30X) | WES (100-150X) |
|---|---|---|
| **Genome covered** | ~100% | ~2% (exons only) |
| **Raw data size** | 60-90 GB FASTQ | 5-10 GB FASTQ |
| **BAM size** | 80-120 GB | 5-15 GB |
| **Variants called** | ~4-5 million | ~30,000-50,000 |
| **Mean depth (on-target)** | 30X | 100-150X |
| **Cost** | $200-1,000 | $100-500 |
| **Detects coding SNPs/indels** | Yes | Yes (higher depth) |
| **Detects intronic/intergenic variants** | Yes | No |
| **Detects structural variants** | Yes (large SVs) | Limited |
| **Detects repeat expansions** | Yes | Only within captured regions |
| **Telomere length** | Yes | No |
| **CNV calling (genome-wide)** | Yes | No |

**Bottom line:** WES is optimized for finding coding variants -- the mutations most likely to directly affect protein function. It has higher per-base depth in exonic regions than typical WGS, which can mean better sensitivity for rare coding variants. But it is completely blind to non-coding regions, structural variants, telomeres, and most repeat expansions.

---

## What You Need: The Capture BED File

WES data is only meaningful in the regions targeted by the capture kit. Outside those regions, coverage drops to near zero and any reads that map there are off-target noise.

You must provide a **BED file** defining the capture regions. This file tells the pipeline where your data actually has coverage. Without it, coverage statistics will be misleading (showing ~2X average across the whole genome instead of ~100-150X on-target) and some tools will waste time analyzing off-target regions.

### How to Set It

Pass the BED file path to the Nextflow pipeline as `--capture_bed`, or to the scripts via the `CAPTURE_BED` environment variable:

```bash
export CAPTURE_BED="${GENOME_DIR}/reference/capture_regions.bed"
```

### Where to Get Your BED File

Contact your sequencing provider or download from the capture kit vendor's website:

| Capture Kit | Vendor | Exonic Targets | BED Download |
|---|---|---|---|
| SureSelect Human All Exon V8 | Agilent | ~35 Mb | [Agilent SureDesign](https://earray.chem.agilent.com/suredesign/) (free account) |
| SureSelect Clinical Research Exome V3 | Agilent | ~54 Mb | Same portal |
| Twist Human Core Exome | Twist Bioscience | ~33 Mb | [Twist support](https://www.twistbioscience.com/resources/data-files) |
| xGen Exome Hyb Panel v2 | IDT | ~34 Mb | [IDT download page](https://www.idtdna.com/pages/products/next-generation-sequencing/hybridization-capture/lockdown-panels/xgen-exome-panel) |
| KAPA HyperExome | Roche | ~43 Mb | [Roche Sequencing Solutions](https://sequencing.roche.com/) |
| Nextera Rapid Capture Exome | Illumina | ~37 Mb | [Illumina support](https://support.illumina.com/) |

> **Build matters.** Make sure the BED file matches your reference genome build (GRCh38/hg38). If you only have an hg19/GRCh37 BED, convert it with UCSC liftOver or Picard LiftoverVcf before using it.

> **If you cannot find your BED file:** Ask your sequencing provider directly. The lab that performed your sequencing knows exactly which capture kit was used and can provide the target regions file. Without it, coverage QC and some variant-calling adjustments cannot be applied correctly.

---

## Per-Step WES Compatibility

### Works Unchanged

These steps operate on VCF data and do not depend on genome-wide BAM coverage. They work identically for WES and WGS.

| Step | Name | Why It Works |
|---|---|---|
| **6** | ClinVar Screen | Checks your variants (whatever you have) against ClinVar |
| **8** | HLA Typing (T1K) | HLA genes are included in all major exome capture kits |
| **12** | Mito Haplogroup | Works if your capture kit includes mitochondrial targets (most clinical kits do) |
| **13** | VEP Annotation | Annotates whatever variants you supply |
| **17** | CPSR | Cancer predisposition screening on your VCF. All ACMG SF v3.2 genes are coding |
| **25** | PRS | Matches scoring file variants against your VCF. Expect ~5-15% variant matching (far fewer than from WGS) because most PRS scoring files include non-coding variants |

### Runs, but Some Answers Are Incomplete

These steps run as they are, but an exome does not cover every position they read. Read [Drug-Gene Answers on an Exome](#drug-gene-answers-on-an-exome) before you use their results.

| Step | Name | What Is Missing |
|---|---|---|
| **7** | PharmCAT | Some defining variants of CYP2C19, CYP3A5, CYP3A4, UGT1A1, DPYD, CYP2B6, VKORC1, IFNL3, CFTR and F2 sit outside the coding exons, where a typical capture does not read |
| **27** | CPIC Recommendations | Reads PharmCAT's calls, so it carries the same gaps |

### Needs Adjustment

These steps require different parameters when processing WES data.

| Step | Name | What Changes | Why |
|---|---|---|---|
| **3** | DeepVariant | `MODEL_TYPE=WES ./scripts/03-deepvariant.sh sample` | DeepVariant has separate trained models for WGS and WES. The WES model is calibrated for the higher depth, sharper coverage boundaries, and different error profiles of exome data. Using the WGS model on WES data inflates false positive rates at capture region boundaries. |
| **4** | Manta (SVs) | Add `--exome` flag to `configManta.py` (the pipeline does, with `--model_type WES`) | Without `--exome`, Manta applies its whole-genome depth filters, which drop calls in the deep coverage of the capture regions. The `--exome` flag turns those depth filters off. It does not restrict calls to the capture regions; `--manta_call_regions` does that, with a bgzipped, indexed BED, and `--model_type WES` does not set it. |
| **9** | ExpansionHunter | Limited to captured loci only | ExpansionHunter estimates repeat lengths from reads spanning the repeat region. For WES, only loci that fall within capture regions have sufficient coverage. Most clinically important repeat expansions (Huntington's HTT, Fragile X FMR1, ALS C9orf72) are in UTRs or introns and are **not captured** by standard exome kits. Expect very few usable results. |
| **16b** | mosdepth | `--by ${CAPTURE_BED}` instead of genome-wide | Calculates on-target coverage (mean depth within capture regions) rather than genome-wide coverage. Without the BED file, mosdepth reports ~2X mean coverage across the whole genome, which is technically correct but completely misleading for WES. |

### Skip These Steps

These steps require genome-wide data and produce incorrect or meaningless results with WES.

| Step | Name | Why It Fails with WES |
|---|---|---|
| **4b** | GRIDSS | Assembly-based SV caller that needs genome-wide read pairs. Off-target WES regions have no coverage, causing GRIDSS to generate thousands of false positive breakpoints at capture boundaries. |
| **10** | TelomereHunter | Estimates telomere length from telomeric repeat reads (TTAGGG). Exome capture kits do not target telomeres, so there are essentially zero telomeric reads. Output will show near-zero telomere content regardless of actual telomere length. |
| **11** | ROH Analysis | `bcftools roh` reads runs of homozygous genotypes. An exome leaves gaps of megabases between the captured exons, so each run and the share of the genome in runs rest on a few sites. Large runs (>5 Mb) may still show, but their edges and the totals are not reliable. |
| **16** | indexcov (Coverage QC and sex check) | Designed for whole-genome BAMs. Reports normalized coverage across 16kb windows, which is meaningless when 98% of the genome has zero coverage by design. Use mosdepth with `--by ${CAPTURE_BED}` instead (step 16b), and take the sex check from somalier (step 33). |
| **18** | CNVpytor | Detects CNVs from read-depth signal across the genome. WES has coverage in ~2% of the genome, making genome-wide read-depth analysis impossible. For exome CNV calling, dedicated tools like ExomeDepth or XHMM are needed (not currently in this pipeline). |
| **19** | Delly | Relies on discordant read pairs and split reads across the genome. WES read pairs that span capture boundaries look like structural variant evidence to Delly, causing massive false positive rates. |
| **20** | Mito analysis (GATK Mutect2) | Mutect2 mitochondrial mode expects reads covering the full mtDNA circle. Most exome kits capture few or no mitochondrial regions, so heteroplasmy detection fails. If your kit does capture mtDNA, step 12 (haplogrep3) still works. The mtDNA copy number weighs chrM's reads against the autosomes', which an exome reads only inside its capture. |
| **21** | Cyrius (CYP2D6) | Cyrius is a whole-genome caller: it calls CYP2D6 copy number from the depth at CYP2D6 and CYP2D7 against the sample's depth elsewhere, and an exome's capture depth gives neither. A wrong call can reach the drug answers through the CYP2D6 consensus. |
| **26** | Ancestry projection | pgsc_calc places the sample on a reference panel's principal components from common SNVs spread over the whole genome. An exome covers few of them, and its coding variants are biased toward conserved regions, so the placement is unreliable. |
| **35** | Parascopy (SMN1/SMN2 and other paralogs) | Parascopy weighs the depth of each duplicated gene against background windows spread across GRCh38, which an exome does not read. The copy numbers, among them SMN1's for spinal muscular atrophy carrier status, would rest on that missing background. |
| **43** | Archaic ancestry (IBDmix) | IBDmix finds Neanderthal and Denisovan segments from runs of sites along the genome; an exome holds too few of them per segment. |
| **47** | Viral DNA | Reads the unmapped reads of a whole genome. An exome's capture keeps few reads that are not human, so a virus present in the blood can show as absent. |
| **50** | Mobile element insertions (MEGAnE) | MEGAnE is made for paired-end whole genomes and its own documentation says it does not support exomes: it calls insertions from reads anywhere in the genome and genotypes them against the genome's depth. |
| **51** | Mosaic alterations (MoChA) | Reads allele balance along whole chromosome arms, from heterozygous sites spread over the genome. An exome has too few of them to see a mosaic event. |
| **53** | Hard loci (Locityper) | Locityper is a whole-genome genotyper: `locityper preproc` models the read depth of a background region, and an exome reads the target genes only where the capture has probes. |

---

## Expected Variant Counts

If you have run WGS before, the variant counts from WES will look very different:

| Metric | WGS (30X) | WES (100-150X) |
|---|---|---|
| Total variants | ~4-5 million | ~30,000-50,000 |
| PASS variants | ~4.5 million | ~25,000-40,000 |
| Coding variants | ~20,000-25,000 | ~20,000-25,000 |
| ClinVar hits (pathogenic + LP) | 0-10 typical | 0-10 typical |
| PharmCAT gene calls | 20+ genes | 20+ genes, some without a defining position ([see below](#drug-gene-answers-on-an-exome)) |
| VEP HIGH impact | 100-200 | 80-150 |

The coding variant counts are nearly identical because WES specifically targets those regions. Where WES differs is the total count: you are missing ~99% of variants, almost all of which are in non-coding regions. For clinical interpretation of coding mutations, WES and WGS produce comparable results.

---

## Quick Start for WES Data

### 1. Gather Your Files

You need:
- **BAM or FASTQ** from your WES provider
- **Capture BED file** from your provider or the kit vendor (see table above)
- The standard pipeline reference data ([step 00](00-reference-setup.md))

Place your BED file in the reference directory:
```bash
cp /path/to/your/capture_regions.bed "${GENOME_DIR}/reference/capture_regions.bed"
```

### 2. Set Environment Variables

```bash
export GENOME_DIR=/path/to/your/data
export SAMPLE=your_name
export CAPTURE_BED="${GENOME_DIR}/reference/capture_regions.bed"
```

### Run the Nextflow Pipeline

The pipeline runs an exome with `--model_type WES` and the capture BED. It needs `sample_qc` in `--tools`, with the data `scripts/setup.sh --sample-qc-data` installs, because the sex check of an exome comes from somalier:

```bash
nextflow run main.nf -profile docker \
  --input samplesheet.csv --reference /path/to/GRCh38_no_alt_analysis_set.fasta --outdir ./results \
  --model_type WES --capture_bed "$CAPTURE_BED" \
  --tools pharmcat,cpic,vcfanno,prs,mito_haplogroup,mosdepth,mito_variants,sample_qc,traits,html_report,multiqc \
  --somalier_sites "${GENOME_DIR}/reference/somalier/sites.hg38.vcf.gz" \
  --verifybamid2_panel "${GENOME_DIR}/reference/verifybamid2"
```

With `--model_type WES` the pipeline:

- calls with DeepVariant's exome model on the regions of `--capture_bed`. Step 03 with `MODEL_TYPE=WES` calls the whole BAM; the calls inside the capture are the ones to read;
- measures coverage over the capture regions, as step 16b with `CAPTURE_BED`;
- runs Manta with `--exome`;
- skips indexcov and takes the sex check from somalier (`sample_qc`, step 33). A sex that disagrees with the samplesheet stops the run, as on a whole genome, but DeepVariant does not wait for it;
- does not run `cnvpytor`, `roh`, `telomere_hunter`, `archaic`, `mocha`, `viral`, `cyrius`, `parascopy`, `locityper` or `mei` when they are in `--tools`, nor the mtDNA copy number and large deletions of `mito_variants`. The log says "not run: exome" and names them, and the report lists them as not run. `survivor_merge` goes too when CNVpytor was one of its two callers. The NuMT filter of `mito_variants` gets the exome's median autosomal depth, 0 or close to it, so it marks little or nothing.
- reports the on-target mean depth: the `total_region` row of `<outdir>/<sample>/coverage/<sample>.mosdepth.summary.txt`, the mean over the capture regions, not the genome-wide `total` row (a few X for an exome). The report's QC lines, the live page and the genome facts' depth floor all read it. The report takes the `total_region` row whenever it is under half as long as the `total` row (500 bp windows cover the whole genome, a capture kit 1-3% of it), so a report built again with `24-html-report.sh`, a run through `run-all.sh` and step 16b with `CAPTURE_BED` read it too.

The pipeline refuses `--model_type WES` without `--capture_bed`, `--capture_bed` without `--model_type WES`, `--intervals` with `--capture_bed` (both set DeepVariant's regions: cut the capture BED to your regions instead), and a capture BED none of whose contigs is in the reference, such as a kit BED named `1` instead of `chr1`: use the kit's GRCh38 BED, or add `chr` to the names. Delly and the ancestry projection still run when you select them; leave them out, as the table above says. `run-all.sh` passes the options after the sex to Nextflow, so `run-all.sh <sample> <sex> --model_type WES --capture_bed <bed>` runs an exome under the same rules. It selects `sample_qc` only when step 33's data is installed; without it the run stops.

### Or Run the Scripts One by One

`MODEL_TYPE=WES` for DeepVariant and `CAPTURE_BED` for mosdepth are the only settings the scripts read for exome data.

### 3. Align (if starting from FASTQ)

```bash
./scripts/02-alignment.sh $SAMPLE
```

Alignment works identically for WES and WGS data. minimap2 does not need different settings.

### 4. Call Variants

```bash
MODEL_TYPE=WES ./scripts/03-deepvariant.sh $SAMPLE
```

### 5. Run VCF-Based Analysis

These steps run with no exome setting:

```bash
./scripts/06-clinvar-screen.sh $SAMPLE     # ClinVar pathogenic variants
./scripts/07-pharmacogenomics.sh $SAMPLE   # PharmCAT drug-gene interactions
./scripts/13-vep-annotation.sh $SAMPLE     # Functional annotation
./scripts/17-cpsr.sh $SAMPLE               # Cancer predisposition screening
./scripts/27-cpic-lookup.sh $SAMPLE        # Drug recommendations
```

Steps 7 and 27 run as they are, but some of their answers are incomplete on an exome: read [Drug-Gene Answers on an Exome](#drug-gene-answers-on-an-exome) first.

### 6. Run BAM-Based Steps (with Adjustments)

```bash
# Manta: scripts/04-manta.sh has no exome option; run configManta.py by hand with --exome
# mosdepth: pass capture BED for on-target coverage stats
CAPTURE_BED=${GENOME_DIR}/reference/your_capture.bed ./scripts/16b-mosdepth.sh $SAMPLE
```

### 7. Skip These Steps

Do **not** run: ROH (11), TelomereHunter (10), indexcov (16), CNVpytor (18), Delly (19), GRIDSS (4b), Mito analysis (20), Cyrius (21), Ancestry projection (26), Parascopy (35), Archaic ancestry (43), Viral DNA (47), Mobile element insertions (50), MoChA (51), Locityper (53).

---

## Drug-Gene Answers on an Exome

PharmCAT (step 7) calls a gene from a fixed list of positions. Some of the variants that define a star allele sit in a promoter, an intron or a UTR, and a typical exome capture (the coding exons plus about 50 bases on each side) does not read them. PharmCAT drops a position it has no call for. It cannot call an allele whose defining variant sits there, and it calls the gene from the positions it did read. When those match the reference, the gene can come out as \*1/\*1 with normal function, even if you carry the allele. A few stray reads outside the capture can also give such a position a weak reference call, which PharmCAT reads like any other.

These are the PharmCAT 3.4.0 positions that a capture of the coding exons with 50 bases of padding does not read (exons of the RefSeq Select transcript, GRCh38):

| Gene | Position | Where | Alleles it defines | What a missed allele looks like |
|---|---|---|---|---|
| CYP2C19 | rs12248560 (c.-806C>T) | promoter | \*17 (and \*4, still called from its start-codon variant) | a rapid or ultrarapid metabolizer reads as normal |
| CYP3A5 | rs776746 (\*3, 6986A>G) | intron 3 | \*3 | the most common allele in Europeans, a non-expresser, reads as \*1 (expresser): this changes the tacrolimus dose |
| CYP3A4 | rs35599367 | intron 6 | \*22 and \*37 (a \*37 then matches \*3) | decreased function reads as normal |
| UGT1A1 | rs887829 (c.-364C>T) | promoter | \*80 and its combinations with \*28 and \*37 | decreased function reads as normal |
| UGT1A1 | rs3064744, the TA repeat | promoter, about 55 bases before the start codon | \*28, \*36, \*37 | decreased function reads as normal, unless the kit pads more than 55 bases and the repeat has depth at the capture edge |
| DPYD | rs75017182 (c.1129-5923C>G) | deep intron | HapB3 | HapB3 is called from c.1236G>A alone, the coding variant that travels with it; c.1129-5923C>G on its own is not seen |
| CYP2B6 | rs34223104 (c.-82T>C) | promoter, 82 bases before the start codon | \*22 and four rarer alleles | increased function reads as normal, unless the kit pads more than 82 bases |
| VKORC1 | rs9923231 (c.-1639G>A) | promoter | the warfarin-sensitive allele (T) | the gene has no call, and step 27's single-position fallback reports the position as not callable |
| IFNL3 | rs12979860 | in IFNL4, about 3.2 kb from IFNL3's coding exons | the variant allele (T) | the gene has no call |
| CFTR | rs75039782 (3849+10kbC>T) | deep intron | an ivacaftor-responsive variant | a carrier of this variant reads as having none |
| F2 | rs1799963 (20210G>A) | 3' UTR, 97 bases after the stop codon | prothrombin 20210A | the gene has no call, unless the kit targets the UTRs |

Every position of ABCG2, CACNA1S, CYP2C9, CYP4F2, F5, G6PD, NAT2, NUDT15, RYR1, SLCO1B1 and TPMT lies within 50 bases of a coding exon, so a capture reads them where it has depth. CYP2D6, HLA-A, HLA-B and MT-RNR1 do not come from these positions: see [step 7's limitations](07-pharmacogenomics.md#limitations).

Your kit's BED file and the reads at its edges decide the exact list. PharmCAT writes the positions it had no call for to `${SAMPLE}.missing_pgx_var.vcf`. Check that file before you act on a call in one of the genes above, and treat a normal-function call there as unconfirmed. [Step 7's limitations](07-pharmacogenomics.md#limitations) say what the reports show for these genes.

---

## WES Coverage QC: What to Look For

With WGS, you check for ~30X mean genome-wide coverage. With WES, the key metrics are different:

| Metric | Good WES | Concerning | How to Check |
|---|---|---|---|
| **Mean on-target depth** | >80X | <50X | mosdepth with `--by` BED |
| **Fraction of targets at 20X** | >95% | <85% | mosdepth thresholds file |
| **Fraction of targets at 10X** | >98% | <90% | mosdepth thresholds file |
| **On-target rate** | >70% | <50% | On-target reads / total reads |
| **Duplication rate** | <20% | >30% | fastp or Picard MarkDuplicates |

Low on-target rates mean the capture was inefficient and many reads are wasted on non-target regions. This does not affect variant calling in captured regions but means you paid for sequencing that produced no usable data.

---

## WES Vendors and Capture Kits

Most clinical and research WES providers use one of a few standard capture kits. Knowing which kit was used helps you find the correct BED file:

| Provider | Typical Kit | Notes |
|---|---|---|
| Invitae, GeneDx (clinical) | Custom panels | Usually target specific gene lists, not the full exome |
| Ambry Genetics | Agilent SureSelect CRE | Clinical-grade exome |
| Blueprint Genetics | Twist or IDT | Research and clinical |
| Novogene, BGI (research) | Agilent SureSelect V6/V7/V8 | Common for research WES |
| Macrogen | IDT xGen or Agilent | Varies by order |
| Illumina (TruSeq/Nextera) | Nextera Rapid Capture or Illumina Exome Panel | Older kits; BED files on Illumina support site |

> **Clinical vs research WES:** Clinical exome labs often use curated gene panels that cover fewer regions than a full research exome kit. If your report says "clinical exome" or names specific gene panels, the target region may be smaller than a standard exome capture. The pipeline still works -- just expect fewer total variants.

---

## Limitations

1. **No structural variant calling.** WES cannot reliably detect deletions, duplications, inversions, or translocations. The exome SV callers that exist (ExomeDepth, XHMM, CoNIFER) detect exon-level CNVs from read depth, but they are not currently in this pipeline.

2. **No telomere or repeat expansion analysis.** These require reads from non-coding regions that WES does not capture.

3. **Reduced PRS accuracy.** Most polygenic risk scoring files include many non-coding variants. WES typically matches only 5-15% of a scoring file's variants, far fewer than WGS. Imputation from exome data does not help much because the target regions are too sparse for accurate imputation of flanking non-coding variants.

4. **No runs of homozygosity.** Exome targets are spread across all chromosomes, but the gaps of megabases between them leave each run, and the share of the genome in runs, resting on a few sites. The pipeline does not run ROH on an exome.

5. **Some drug-gene answers are incomplete.** Defining variants of CYP2C19, CYP3A5, CYP3A4, UGT1A1, DPYD and other genes sit outside the coding exons, so PharmCAT can call a normal-function allele where you carry one it could not see. See [Drug-Gene Answers on an Exome](#drug-gene-answers-on-an-exome).

6. **Ancestry estimation does not work.** The projection onto a reference panel (step 26) needs common variants spread across the genome, which exome data cannot provide.

7. **Deep intronic and regulatory variants are invisible.** Some clinically significant variants are in deep intronic regions (e.g., splicing regulators) or promoters/enhancers. WES will not detect these. If you have a clinical suspicion that WES missed, WGS is the appropriate follow-up.

---

## Should You Get WES or WGS?

If you are choosing between the two:

| Factor | WES Advantage | WGS Advantage |
|---|---|---|
| **Cost** | Cheaper ($100-500) | More comprehensive ($200-1,000) |
| **Storage** | Smaller files (5-15 GB BAM) | Much larger (80-120 GB BAM) |
| **Coding variant detection** | Higher per-base depth in exons | Same variants, slightly lower depth |
| **Non-coding variants** | None | Full coverage |
| **Structural variants** | Very limited | Comprehensive |
| **Pharmacogenomics** | Misses defining variants in promoters and introns (CYP2C19 \*17, CYP3A5 \*3, UGT1A1 \*28 and others) | Reads them |
| **Repeat expansions** | Very limited | Full |
| **ClinVar screening** | ~90% of ClinVar pathogenic variants are coding | Catches the other ~10% too |
| **Future reanalysis value** | Limited to exons | Full genome for future discoveries |

**If you already have WES data:** This pipeline extracts substantial value from it. Run the compatible steps and you will get ClinVar screening, cancer predisposition and functional annotation that is comparable to WGS for coding regions, and drug-gene answers with the gaps listed in [Drug-Gene Answers on an Exome](#drug-gene-answers-on-an-exome).

**If you are ordering new sequencing:** WGS is generally recommended because the price difference has shrunk dramatically while WGS captures everything WES does plus non-coding regions, structural variants, and repeat expansions. The extra data storage is the main downside.
