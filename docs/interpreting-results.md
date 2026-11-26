# Interpreting Your Results

You've run the pipeline. Now you have directories full of VCFs, TSVs, and HTML reports. This guide explains what to look at first and what it all means — no bioinformatics degree required.

## Before You Look: What the Pipeline Can Tell You

Decide what you want to know before you open the reports. A default run looks further than most people expect:

- **Secondary findings are on by default.** Step 17 runs CPSR with `--secondary_findings`. Besides the cancer genes, CPSR then reports pathogenic and likely pathogenic variants in the genes of the ACMG SF list of secondary findings. That list includes genes for inherited heart conditions that can cause sudden death (cardiomyopathies, arrhythmias), familial hypercholesterolaemia and some metabolic diseases. A finding there can be serious and actionable even though you never asked about it.
- **APOE and Alzheimer's disease.** Step 25 computes a score for late-onset Alzheimer's disease (PGS000334). Its two largest weights are the APOE variants rs429358 and rs7412, which define the ε2, ε3 and ε4 alleles, and the summary report lists that score with the others. The ε4 allele raises the risk; it does not say who will get the disease. Many people choose not to learn their APOE status. Decide before you open the step 25 output.
- **Your relatives.** You share half your DNA with each parent, child and sibling. A pathogenic variant in a dominant gene means each of them has a 50% chance of carrying it too, so a result about you is also information about them. Comparing two genomes (see [multi-sample](multi-sample.md)) can also reveal unexpected family relationships.
- **Insurance.** The rules depend on where you live. In the United States, GINA stops health insurers and employers from using genetic information, but it does not cover life, disability or long-term-care insurance. Elsewhere the rules differ. Some insurers ask whether you have had a genetic test, and a result entered in your medical record can count. Check what applies to you before you act on a finding.
- **Some genes cannot be read well with short reads.** In a gene inside a segmental duplication, over a false duplication of GRCh38 or only on an ALT contig, short reads often fit two places equally and get mapping quality 0, so no caller uses them. A "no finding" in such a gene says nothing. The report's section "Genes this genome could not read reliably" (Quality chapter, from step 16b) lists the genes of `assets/callability_genes.tsv` where fewer than 90% of the bases had 10 or more reads of mapping quality 20 or more, with the share each reached, and the genes the default reference cannot hold at all ("not on this reference": GSTT1, HLA-DRB3, HLA-DRB4, several KIR genes and others). Genes with a step of their own (SMN1 and SMN2 in step 35, the pypgx genes in step 32, HBA1, HBA2, LPA and RHD in the depth engine) are listed apart; read that step's section for them. GSTT1 is not one of them: step 32 cannot call it on the default reference. [Step 16b](16b-mosdepth.md#gene-callability) says where the list comes from.
- **Switching secondary findings off.** Set `CPSR_SECONDARY_FINDINGS=false` before `run-all.sh` or `scripts/17-cpsr.sh` (Nextflow: `--cpsr_secondary_findings false`). CPSR then reports the cancer panel only. The report's own secondary-findings card (step 23) is computed on every run and does not follow this switch: [Turning CPSR's secondary findings off](results/secondary-findings.md#cpsr-report-step-17).

---

## Before You Panic: What Every Genome Looks Like

If this is your first time looking at your own genomic data, the numbers can be alarming. Here is what a **completely normal, healthy person's genome** looks like:

| Finding | Normal Range | Why It Seems Scary |
|---|---|---|
| Total variants | 4.5-5.5 million | Sounds like millions of "mutations" — but >99.9% are normal human variation |
| ClinVar pathogenic hits (step 6) | 0-10 | Step 6 screens against Pathogenic/Likely_pathogenic only. Most hits are recessive carriers — you need TWO copies to be affected |
| HIGH impact variants (VEP) | 100-150 | Most are heterozygous in non-essential genes. For recessive genes, one copy is typically tolerated (but see haploinsufficiency). |
| Structural variants | 5,000-10,000 | Most are in non-coding regions. Your parents had them too. |
| Heteroplasmic mitochondrial variants | 20-40 | Low-level heteroplasmy (<5%) is universal and age-related |
| VUS (Variants of Uncertain Significance) | 20-200+ | "Uncertain" means **not enough data yet**, not "probably bad" |

---

## Start Here: The Three Most Important Outputs

1. [ClinVar Screen (Step 6)](results/clinvar.md#clinvar-screen-step-6){ #1-clinvar-screen-step-6 }
2. [PharmCAT Report (Step 7)](results/pharmacogenomics.md#pharmcat-report-step-7){ #2-pharmcat-report-step-7 }
3. [CPSR Report (Step 17)](results/secondary-findings.md#cpsr-report-step-17){ #3-cpsr-report-step-17 }

---

## The Report Chapters

Each chapter of the report has its own page in this guide. The links below keep the old addresses of the sections on this page working.

- [Secondary Findings and Carrier Status](results/secondary-findings.md)
    - [Carrier Status Is Not Disease](results/secondary-findings.md#carrier-status-is-not-disease){ #carrier-status-is-not-disease }
    - [Clinical Filter (Step 23)](results/secondary-findings.md#clinical-filter-step-23){ #clinical-filter-step-23 }
    - [Variant Prioritization (Step 31)](results/secondary-findings.md#variant-prioritization-step-31){ #variant-prioritization-step-31 }
- [ClinVar and Variant Annotation](results/clinvar.md)
    - [The VUS Trap](results/clinvar.md#the-vus-trap){ #the-vus-trap }
    - [ClinVar Star Ratings](results/clinvar.md#clinvar-star-ratings){ #clinvar-star-ratings }
    - [VEP Annotation (Step 13)](results/clinvar.md#vep-annotation-step-13){ #vep-annotation-step-13 }: [What It Adds](results/clinvar.md#what-it-adds){ #what-it-adds }, [gnomAD Frequency: Your Best Sanity Check](results/clinvar.md#gnomad-frequency-your-best-sanity-check){ #gnomad-frequency-your-best-sanity-check }, [Filtering Strategy](results/clinvar.md#filtering-strategy){ #filtering-strategy }, [What "HIGH Impact" Means](results/clinvar.md#what-high-impact-means){ #what-high-impact-means }
    - [Pathogenicity Scores (Step 30 — vcfanno)](results/clinvar.md#pathogenicity-scores-step-30-vcfanno){ #pathogenicity-scores-step-30-vcfanno }: [CADD](results/clinvar.md#cadd-combined-annotation-dependent-depletion){ #cadd-combined-annotation-dependent-depletion }, [REVEL](results/clinvar.md#revel-rare-exome-variant-ensemble-learner){ #revel-rare-exome-variant-ensemble-learner }, [AlphaMissense](results/clinvar.md#alphamissense){ #alphamissense }, [SpliceAI](results/clinvar.md#spliceai){ #spliceai }, [gnomAD Gene Constraint (Step 23 summary)](results/clinvar.md#gnomad-gene-constraint-step-23-summary){ #gnomad-gene-constraint-step-23-summary }
    - [Quick Variant Filtering Recipes](results/clinvar.md#quick-variant-filtering-recipes){ #quick-variant-filtering-recipes }
- [Pharmacogenomics and HLA](results/pharmacogenomics.md)
    - [HLA Typing (Step 8)](results/pharmacogenomics.md#hla-typing-step-8){ #hla-typing-step-8 }
    - [More Pharmacogenomics: Cyrius, CPIC, pypgx and the Consensus (Steps 21, 27, 32, 36)](results/pharmacogenomics.md#more-pharmacogenomics-cyrius-cpic-pypgx-and-the-consensus-steps-21-27-32-36){ #more-pharmacogenomics-cyrius-cpic-pypgx-and-the-consensus-steps-21-27-32-36 }
- [Structural Variants, Copy Number and Repeats](results/structural-and-repeats.md)
    - [Structural Variants (Steps 4, 5, 15, 18, 19, 22)](results/structural-and-repeats.md#structural-variants-steps-4-5-15-18-19-22){ #structural-variants-steps-4-5-15-18-19-22 }: [What Are Structural Variants?](results/structural-and-repeats.md#what-are-structural-variants){ #what-are-structural-variants }, [A Note About BND (Breakend) Calls](results/structural-and-repeats.md#a-note-about-bnd-breakend-calls){ #a-note-about-bnd-breakend-calls }, [How Many Is Normal?](results/structural-and-repeats.md#how-many-is-normal){ #how-many-is-normal }, [Which Callers to Trust?](results/structural-and-repeats.md#which-callers-to-trust){ #which-callers-to-trust }, [SV Consensus (Step 22)](results/structural-and-repeats.md#sv-consensus-step-22){ #sv-consensus-step-22 }, [AnnotSV Output](results/structural-and-repeats.md#annotsv-output){ #annotsv-output }
    - [STR Expansions (Step 9)](results/structural-and-repeats.md#str-expansions-step-9){ #str-expansions-step-9 }: [What Are Repeat Expansions?](results/structural-and-repeats.md#what-are-repeat-expansions){ #what-are-repeat-expansions }, [How to Read ExpansionHunter Output](results/structural-and-repeats.md#how-to-read-expansionhunter-output){ #how-to-read-expansionhunter-output }, [Stranger Annotation (Step 9b)](results/structural-and-repeats.md#stranger-annotation-step-9b){ #stranger-annotation-step-9b }
    - [CNVpytor Results (Step 18)](results/structural-and-repeats.md#cnvpytor-results-step-18){ #cnvpytor-results-step-18 }
    - [Delly Results (Step 19)](results/structural-and-repeats.md#delly-results-step-19){ #delly-results-step-19 }
    - [SMN1 and SMN2 Copy Number (Step 35, opt-in)](results/structural-and-repeats.md#smn1-and-smn2-copy-number-step-35-opt-in){ #smn1-and-smn2-copy-number-step-35-opt-in }
    - [Lp(a) and the LPA KIV-2 Count (Step 46, opt-in)](results/structural-and-repeats.md#lpa-and-the-lpa-kiv-2-count-step-46-opt-in){ #lpa-and-the-lpa-kiv-2-count-step-46-opt-in }
    - [Mobile Element Insertions (Step 50, opt-in)](results/structural-and-repeats.md#mobile-element-insertions-step-50-opt-in){ #mobile-element-insertions-step-50-opt-in }
- [Ancestry, Lineage and Polygenic Scores](results/ancestry-and-lineage.md)
    - [ROH Analysis (Step 11)](results/ancestry-and-lineage.md#roh-analysis-step-11){ #roh-analysis-step-11 }: [What It Means](results/ancestry-and-lineage.md#what-it-means){ #what-it-means_1 }, [How to Read](results/ancestry-and-lineage.md#how-to-read){ #how-to-read }, [Centromeric Artifacts](results/ancestry-and-lineage.md#centromeric-artifacts){ #centromeric-artifacts }
    - [Mitochondrial Haplogroup (Step 12)](results/ancestry-and-lineage.md#mitochondrial-haplogroup-step-12){ #mitochondrial-haplogroup-step-12 }: [What It Means](results/ancestry-and-lineage.md#what-it-means_1){ #what-it-means_2 }, [Common European Haplogroups](results/ancestry-and-lineage.md#common-european-haplogroups){ #common-european-haplogroups }, [Medical Relevance](results/ancestry-and-lineage.md#medical-relevance){ #medical-relevance }, [Exploring Your Branch](results/ancestry-and-lineage.md#exploring-your-branch){ #exploring-your-branch }
    - [Mitochondrial Variants (Step 20)](results/ancestry-and-lineage.md#mitochondrial-variants-step-20){ #mitochondrial-variants-step-20 }
    - [Ancient Populations (Step 44, opt-in)](results/ancestry-and-lineage.md#ancient-populations-step-44-opt-in){ #ancient-populations-step-44-opt-in }
    - [Polygenic Risk Scores (Step 25)](results/ancestry-and-lineage.md#polygenic-risk-scores-step-25){ #polygenic-risk-scores-step-25 }
- [Telomeres and Somatic Variants](results/traits-and-somatic.md)
    - [Telomere Length (Step 10)](results/traits-and-somatic.md#telomere-length-step-10){ #telomere-length-step-10 }: [What It Means](results/traits-and-somatic.md#what-it-means){ #what-it-means }, [How to Interpret](results/traits-and-somatic.md#how-to-interpret){ #how-to-interpret }, [Limitations](results/traits-and-somatic.md#limitations){ #limitations }
    - [Somatic Variants (Step 29) [EXPERIMENTAL]](results/traits-and-somatic.md#somatic-variants-step-29-experimental){ #somatic-variants-step-29-experimental }
- [Data Quality and Reports](results/quality.md)
    - [Reads, Alignment, Variant Calling and Coverage (Steps 1b, 2, 3, 16, 16b)](results/quality.md#reads-alignment-variant-calling-and-coverage-steps-1b-2-3-16-16b){ #reads-alignment-variant-calling-and-coverage-steps-1b-2-3-16-16b }
    - [Reports (Steps 24 and 28)](results/quality.md#reports-steps-24-and-28){ #reports-steps-24-and-28 }

---

## What to Do Next

1. **Share your PharmCAT report** with your prescribing physician or pharmacist
2. **Review ClinVar pathogenic hits** — check if any are in dominant genes or if you're homozygous for recessive genes
3. **Read the CPSR HTML report** — it's designed for clinical interpretation and will highlight anything that needs attention
4. **If you find something concerning:** Don't panic. Discuss with a genetic counselor. Many "pathogenic" variants have incomplete penetrance (not everyone with the variant gets the disease)
5. **For carrier status findings:** Relevant mainly for family planning. If both partners carry the same recessive condition, each child has a 25% chance of being affected

---

## Re-running with Updated Databases

Genomic databases are updated continuously. Variants classified as VUS today may be reclassified next year. Periodic re-analysis is one of the most valuable things you can do.

### What to Update and When

| Database | Update Frequency | Pipeline Steps Affected | How to Update |
|---|---|---|---|
| ClinVar | Monthly (the release `setup.sh` downloads; NCBI also builds weekly files) | Step 6 (ClinVar screen) | Re-download from NCBI FTP (see [00-reference-setup.md](00-reference-setup.md)) |
| VEP cache | Every 6 months | Step 13 (VEP annotation) | Download new release from Ensembl FTP |
| PCGR/CPSR data | Every 6-12 months | Step 17 (CPSR) | Download new bundle from PCGR GitHub releases |
| PharmCAT | Every few months | Step 7 (pharmacogenomics) | Bump `PHARMCAT_IMAGE` in `versions.env`, then pull it |

### Recommended Re-analysis Schedule

- **Every 6 months:** Re-run steps 6 (ClinVar) and 17 (CPSR) with updated databases. These are quick steps (about 20-40 minutes together, measured; CPSR is most of it) and the most likely to have new classifications.
- **Every 12 months:** Re-run step 13 (VEP) with updated cache for new gnomAD frequencies and consequence predictions.
- **After major database releases:** ClinVar periodically reclassifies large batches of variants. Follow [@ClinVarUpdates](https://twitter.com/ClinVarUpdates) or check the NCBI blog for announcements.

### What You Do NOT Need to Re-run

- Steps 2-3 (alignment + variant calling): Your variants don't change. Only re-run if a major DeepVariant version is released with improved accuracy.
- Steps 4, 18, 19 (SV callers): Structural variant calling is compute-intensive and results don't change with database updates.
- Step 10 (telomere): Telomere content doesn't change with database updates.

---

## Example Outputs: What Correct Results Look Like

Invented examples in the real output formats, so you know what to expect. Every number and every call below is made up; none comes from a real genome.

### Variant Calling (Step 3)

```
bcftools stats output:
SN  0  number of samples:     1
SN  0  number of records:     5500000
SN  0  number of SNPs:        4200000
SN  0  number of indels:      1300000
SN  0  number of multiallelic sites:  45000

# PASS variants only: 4,650,000-4,700,000
# Ti/Tv ratio: 2.05-2.10 (if < 1.8, something is wrong)
```

### ClinVar Screen (Step 6)

```
Pathogenic/Likely Pathogenic hits: 3

  chr<N>:<pos> <REF>><ALT> (rs<id>)   — <GENE> carrier (<condition>, recessive)
  chr<N>:<pos> <REF>><ALT> (rs<id>)   — <GENE> carrier (<condition>, recessive)
  chr<N>:<pos> <REF>><ALT> (rs<id>)   — <GENE> carrier (<condition>, recessive)
```

A handful of heterozygous (0/1) hits in recessive genes is the usual result and means carrier status only.

### PharmCAT (Step 7)

The HTML report will show a table like:

```
Gene        Diplotype           Phenotype              Affected Drugs
CYP2C19    *x/*y               <phenotype>             PPIs, SSRIs, clopidogrel
CYP2C9     *x/*y               <phenotype>             Warfarin, NSAIDs
NAT2       *x/*y               <phenotype>             Isoniazid, hydralazine
DPYD       *x/*y               <phenotype>             5-FU
SLCO1B1    *x/*y               <phenotype>             Statins
```

The phenotype column reads Poor, Intermediate, Normal, Rapid or Ultrarapid Metabolizer (Normal or Decreased Function for transporters such as SLCO1B1).

Typically 18-21 of PharmCAT's 23 genes will have confident calls. CYP2D6 may be "Inconclusive" from short-read WGS (known limitation).

### ExpansionHunter (Step 9)

```json
{
  "LocusResults": {
    "HTT": { "Genotype": "<repeats>/<repeats>" },
    "FMR1": { "Genotype": "<repeats>" },
    "C9orf72": { "Genotype": "<repeats>/<repeats>" },
    "ATXN1": { "Genotype": "<repeats>/<repeats>" },
    "DMPK": { "Genotype": "<repeats>/<repeats>" }
  }
}
```

Each number is the repeat count on one allele. The report shows these counts and, with step 9b (Stranger), lists the loci it marks `pre_mutation` or `full_mutation`. It gives no all-clear verdict.

### CPSR (Step 17)

The HTML report counts the variants in each class:

```
Pathogenic:                  0 variants
Likely pathogenic:           0 variants
VUS:                        22 variants
Likely benign:             <n> variants
Benign:                    <n> variants
Secondary findings:          0 variants
```

No Pathogenic or Likely pathogenic variant, in the cancer panel or in the secondary findings, means nothing actionable was found in those genes. The VUS count varies widely (20-200+) and is not cause for concern.

### ROH (Step 11)

```
Autosomal ROH >5MB (potential consanguinity signal):

NOTE: Centromeric ROH (chr1:125-143MB, chr9:42-60MB, chr18:15-20MB) are technical artifacts, not real.
```

No line under the heading means no segment of 5 Mb or more: no evidence of related parents. If there are some, add them up and compare the total with the table in [step 11](11-roh-analysis.md#total-roh-and-parental-relationship).

### Telomere Length (Step 10) { #telomere-length-step-10_1 }

```
tel_content: <value>
```

No universal "normal" range — compare between samples of the same age, sequenced on the same platform.

---

## Investigating Specific Variants

Found something interesting? These free tools help you dig deeper:

### Visual Inspection

- **[IGV Web](https://igv.org/app/)** — Load your BAM file (or a region of it) to visually inspect read-level evidence for a variant. Essential for confirming structural variants and checking for sequencing artifacts.
- **[gene.iobio](https://gene.iobio.io/)** — Clinically-driven variant interrogation tool. Load your VCF and BAM, search for specific genes, and see coverage, variant calls, and population frequency in one view.
- **[UCSC Genome Browser](https://genome.ucsc.edu/)** — Search for any genomic coordinate to see the surrounding genes, conservation, regulatory elements, and known variants.

### Database Lookups

- **[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)** — Search by rsID, gene name, or genomic position. Check the review status (stars) and submission history.
- **[gnomAD](https://gnomad.broadinstitute.org/)** — Search any variant to see its population frequency across 800,000+ individuals (v4). If common in gnomAD, almost certainly benign.
- **[OMIM](https://www.omim.org/)** — The definitive catalog of genetic disorders. Search by gene name to understand what conditions it causes and the inheritance pattern.
- **[GeneReviews](https://www.ncbi.nlm.nih.gov/books/NBK1116/)** — Expert-written disease descriptions for genetic conditions. The single best resource for understanding a specific genetic disease.

---

## Annotation Tool Disagreement

An important caveat: **different annotation tools may classify the same variant differently**.

VEP (used in step 13), SnpEff, and ANNOVAR are the three most common variant annotation tools. Studies have shown that they disagree on consequence predictions for ~5-10% of variants, particularly at splice sites and multi-transcript genes.

**What this means for you:**
- If VEP says a variant is "HIGH impact" but ClinVar says it is benign, **trust ClinVar** (human-reviewed evidence > computational prediction)
- If VEP and ClinVar agree on pathogenicity, this strengthens the interpretation
- If you find a potentially significant variant using VEP that is NOT in ClinVar, search gnomAD for its population frequency before drawing conclusions
- For the most important findings, consider running a second annotation tool as validation

The pipeline uses VEP because it is the most widely used and well-maintained tool, with direct gnomAD frequency integration. But no single tool is perfect.

---

## What This Pipeline Does Not Assess

A clean result from these steps says nothing about the following. Each needs a different test or a different kind of data. [What the pipeline does not report](not-reported.md) has the full list, with the results left out on purpose and a source for each item.

- **Copy-number changes in genes with a near-identical copy**, except the loci the opt-in step 35 counts: SMN1/SMN2, GBA1, PMS2, STRC, NCF1, C4A/C4B, the CFH cluster and the rest of its [paralog catalog](35-paralogs.md#the-paralog-catalog), the number of alpha-globin genes, the dosage loci HP, LCE3B/LCE3C and CCL3L1 ([copy number loci](cn-loci.md)), and CYP21A2 deletions and possible gene conversions, read from depth and six sites where CYP21A2 differs from its pseudogene ([the CYP21A2 row](cn-loci.md#the-cyp21a2-row)). For every catalog locus but SMN1, step 35 counts the gene and its paralog together: it can show a copy lost or gained, not which copy. Most alpha-thalassaemia changes (HBA1/HBA2 deletions) sit in regions where short reads cannot tell the gene from its paralog. The CYP21A2 call is a screen: it does not say which chromosome carries a change, and a two-copy result does not rule out two copies on one chromosome and none on the other. The VCF-based screens do not see them, and without step 35 neither is SMA carrier status. Step 35 counts HBA1 and HBA2 together from read depth: it can show a gene lost, not which one or on which chromosome. The count has been checked on one real genome only (HG002, four genes), so a count under four is marked unconfirmed and needs a clinical test. When the estimate is one gene or more, a BAM whose aligner gives no mapping-quality-0 reads over the two genes (novoalign by default) is reported as not assessed, because the reads that fit both genes are missing from it. Step 35 estimates SMN1 copy number but cannot see a "2+0" carrier; SMA carrier status still needs a clinical SMN1 test.
- **Mobile-element insertions**, unless the opt-in step 50 ran. Manta reports insertions but does not classify them as Alu, LINE-1 or SVA insertions; step 50 (MEGAnE) does, and reads three known disease insertions.
- **Methylation and phasing from long reads.** The long-read branch stops at alignment, small variants and structural variants; see the [long-read guide](long-read-guide.md).
- **Mosaic copy-neutral loss of heterozygosity**, present in only some cells (common in blood with age). Step 11 sees runs of homozygosity that are in every cell; it cannot see a change carried by a fraction of them.
- **Repeat expansions outside STRchive's 82 disease loci** (ExpansionHunter's 31 without `setup.sh`'s STRchive catalog), the seven of them whose repeat short reads cannot size or STRchive gives no threshold for, and accurate sizing of very large expansions.

## Important Caveats

- **This is not a clinical diagnosis.** These tools use the same algorithms as clinical labs, but the pipeline has not been clinically validated.
- **False positives exist.** Short-read WGS has limitations in repetitive regions, homologous genes (CYP2D6, HLA), and structural variants.
- **False negatives exist.** Some pathogenic variants are in regions that short reads can't cover (deep intronic, repeat expansions beyond read length, large structural variants).
- **VUS (Variants of Uncertain Significance)** are not actionable. They may be reclassified in the future as more data accumulates.
- **ClinVar classifications can change.** A variant classified as pathogenic today may be reclassified as benign (or vice versa) as new evidence emerges. Re-run step 6 periodically with updated ClinVar databases.
