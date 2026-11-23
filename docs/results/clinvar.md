# ClinVar and Variant Annotation

Part of [Interpreting Your Results](../interpreting-results.md). Read [Before You Look](../interpreting-results.md#before-you-look-what-the-pipeline-can-tell-you) first.

## The VUS Trap

The single biggest source of unnecessary anxiety in personal genomics is **VUS — Variants of Uncertain Significance**. These are variants where:

- There is not enough scientific evidence to classify them as either pathogenic or benign
- The vast majority will eventually be reclassified as **benign** as more data accumulates
- They are **not actionable** — no clinical decision should be made based on a VUS
- CPSR may report dozens or hundreds of VUS. This is normal and expected.

**Rule of thumb:** If a variant is classified as VUS, treat it the same as if it were not tested. Do not change screening or management based on a VUS. Check back in 1-2 years when ClinVar may have reclassified it.

**One nuance:** If you have a strong family history of a condition AND a VUS appears in the relevant high-penetrance gene (e.g., BRCA1/2, TP53, MLH1/MSH2), it may be worth mentioning to a genetic counselor — not to act on the VUS, but because the family history itself may warrant enhanced screening regardless of the variant's classification.

## ClinVar Star Ratings

Not all ClinVar classifications are equally reliable. Each entry has a **review status** indicated by stars:

| Stars | Review Status | Reliability |
|---|---|---|
| 0 | No assertion criteria | Low — submitter did not explain their reasoning |
| 1 | Single submitter, criteria provided | Moderate — one lab's interpretation |
| 2 | Two or more submitters, no conflict | Good — multiple labs agree |
| 3 | Expert panel reviewed | High — reviewed by specialists |
| 4 | Practice guideline | Highest — established clinical standard |

**Focus on 2+ star entries.** Single-submitter (1-star) pathogenic calls are sometimes reclassified. If you find a scary-looking pathogenic variant with 0-1 stars, check the ClinVar entry directly at [ncbi.nlm.nih.gov/clinvar](https://www.ncbi.nlm.nih.gov/clinvar/) — look at the "Review status" and "Last evaluated" date.

---

## ClinVar Screen (Step 6)

**What it tells you:** Known pathogenic variants in your genome, as classified by ClinVar (NCBI's public database of clinically significant variants).

**Where to look:** `${SAMPLE}/clinvar/`

**How to read it:**
- Each line in the output VCF is a variant in your genome that matches a known ClinVar entry
- Step 6 intersects against the **pathogenic-only** ClinVar subset (Pathogenic + Likely_pathogenic). Every hit in this output is at a position ClinVar classifies as disease-associated — benign/VUS entries are excluded at the database level
- The `CLNSIG` field confirms the classification

**What to expect:**
- 0-10 pathogenic/likely pathogenic hits is typical for a 30X WGS
- Most pathogenic hits are **carrier status** (heterozygous) for recessive conditions — you carry one copy but aren't affected
- A heterozygous pathogenic variant in a recessive gene (like GJB2 for hearing loss) means you're a **carrier**, not affected
- A homozygous pathogenic variant, or a heterozygous variant in a dominant gene, needs attention

**When to worry:**
- Pathogenic variant in a **dominant** gene (one copy is enough to cause disease)
- **Two** pathogenic variants in the same recessive gene (one from each parent)
- Any variant in cancer predisposition genes (BRCA1, BRCA2, MLH1, MSH2, etc.)

**When "no hit" means less:** no hit in a gene says nothing when the reads could not see the gene. The Quality chapter lists the hard-to-read genes this genome did not read ([step 16b](../16b-mosdepth.md#gene-callability)). A gene marked `exon not read` was read over most of its length, but one coding exon (named, as `GENE:exonN`) had fewer than 90% of its bases at 10 reads of mapping quality 20. A pathogenic variant in that exon would not be found, so a ClinVar screen with no hit there is not a clean result for that gene.


---

## VEP Annotation (Step 13)

### What It Adds

VEP annotates every variant with:
- **Consequence type:** missense, nonsense, synonymous, splice site, etc.
- **SIFT score:** Predicts if amino acid change is tolerated (>0.05) or damaging (<0.05)
- **PolyPhen score:** Predicts if change is benign (<0.15), possibly damaging (0.15-0.85), or probably damaging (>0.85)
- **gnomAD frequencies:** How common this variant is in gnomAD's exomes (`gnomADe_AF`) and genomes (`gnomADg_AF`)

### gnomAD Frequency: Your Best Sanity Check

The single most useful annotation VEP adds is the **gnomAD allele frequency** — how common a variant is in the general population. Step 13 runs VEP with `--everything`, which turns on both `--af_gnomade` and `--af_gnomadg`, so each variant gets two frequencies: `gnomADe_AF` from gnomAD's exomes and `gnomADg_AF` from its genomes. A non-coding variant outside exome capture regions has no exome frequency but can still have a genome frequency, so check both before calling a variant absent from gnomAD. The recipes below do that.

**Key principle:** A variant that is common in healthy people is almost certainly benign, regardless of what any prediction tool says.

| gnomAD AF | Interpretation | Action |
|---|---|---|
| > 5% (0.05) | Common polymorphism | Benign. Ignore. |
| 1-5% | Low-frequency variant | Almost certainly benign |
| 0.1-1% | Uncommon | Probably benign, but check ClinVar |
| 0.01-0.1% | Rare | Worth investigating if in a disease gene |
| < 0.01% | Very rare | Potentially significant. Check ClinVar + literature |
| Absent | Novel or ultra-rare | Could be significant OR a sequencing artifact. Verify with a second method |

**If a variant is "pathogenic" in ClinVar but has gnomAD AF > 1%:** The ClinVar entry may be outdated or wrong. Truly pathogenic variants for severe diseases are almost always rare (< 0.1%) because natural selection removes them from the population.

### Filtering Strategy

Step 23 (the [clinical filter](secondary-findings.md#clinical-filter-step-23)) already applies the filters most people want. For your own queries, use the [Quick Variant Filtering Recipes](#quick-variant-filtering-recipes) below. Do not `grep` the VEP VCF: VEP writes its values by position inside the pipe-separated `CSQ` field, never as `name=value`, so a grep for `gnomAD_AF` or `SYMBOL=` matches nothing, and a grep for `HIGH` or `1/1` can match text in other fields. `bcftools +split-vep` reads the fields by name.

### What "HIGH Impact" Means

VEP classifies variant impact as:

| Impact | Types | Interpretation |
|---|---|---|
| HIGH | Stop gained, frameshift, splice donor/acceptor | Likely breaks the protein |
| MODERATE | Missense, in-frame insertion/deletion | Changes the protein, may or may not matter |
| LOW | Synonymous, splice region | Probably no functional effect |
| MODIFIER | Intronic, intergenic, UTR | Usually non-functional |

**Everyone has ~100-150 HIGH impact variants.** Most are in one copy (heterozygous) of non-essential genes. Don't panic at the number.

---

## Pathogenicity Scores (Step 30 — vcfanno)

If you ran step 30, your VCF now includes quantitative pathogenicity scores beyond VEP's qualitative SIFT/PolyPhen predictions. These scores are widely used as computational evidence in variant interpretation (see ClinGen's PP3/BP4 calibration framework).

### CADD (Combined Annotation Dependent Depletion)

Scores **all** variant types (coding, non-coding, splice, regulatory). Uses a PHRED-like scale where higher = more deleterious.

| CADD PHRED | Interpretation | Context |
|---|---|---|
| < 10 | Likely benign | Bottom 90% of genome variation |
| 10-20 | Uncertain | Top 10%, but most are still benign |
| 20-25 | Potentially deleterious | Top 1% — investigate if in a disease gene |
| 25-30 | Likely deleterious | Top 0.3% — strong candidate for pathogenicity |
| > 30 | Highly deleterious | Top 0.1% — likely damaging if in a constrained gene |

**When to use CADD:** Best for non-coding and splice-region variants where SIFT/PolyPhen don't apply. For missense variants, REVEL and AlphaMissense are more specific.

### REVEL (Rare Exome Variant Ensemble Learner)

Scores **missense variants only**. Combines 13 individual tools into a single 0-1 score. Recommended by ClinGen for ACMG PP3/BP4 evidence.

The levels below are the REVEL row of Table 2 in [Pejaver et al. 2022](https://doi.org/10.1016/j.ajhg.2022.10.013), the ClinGen calibration of PP3/BP4:

| REVEL Score | ClinGen Evidence Level | Interpretation |
|---|---|---|
| <= 0.003 | BP4_Very Strong | Very strong evidence of benign |
| > 0.003 to <= 0.016 | BP4_Strong | Strong evidence of benign |
| > 0.016 to <= 0.183 | BP4_Moderate | Moderate evidence of benign |
| > 0.183 to <= 0.290 | BP4_Supporting | Supporting evidence of benign |
| > 0.290 to < 0.644 | No evidence | Uncertain significance |
| >= 0.644 to < 0.773 | PP3_Supporting | Supporting evidence of pathogenicity |
| >= 0.773 to < 0.932 | PP3_Moderate | Moderate evidence of pathogenicity |
| >= 0.932 | PP3_Strong | Strong evidence of pathogenicity |

There is no PP3 Very Strong level for REVEL: no score reached it in the calibration.

**When to use REVEL:** First-line score for evaluating missense variants. If REVEL >= 0.644, investigate the variant seriously.

### AlphaMissense

DeepMind's protein-structure-informed **missense** classifier. Uses AlphaFold2 protein structure predictions to assess amino acid substitution impact.

| am_pathogenicity | am_class | Interpretation |
|---|---|---|
| < 0.34 | likely_benign | Predicted benign by protein structure analysis |
| 0.34-0.564 | ambiguous | Uncertain — use other evidence |
| > 0.564 | likely_pathogenic | Predicted damaging based on protein structure |

These cut-offs are AlphaMissense's own class boundaries, set by its authors. They are not ACMG/ClinGen evidence levels.

**When to use AlphaMissense:** Complements REVEL as a second opinion. If they disagree, investigate further. Agreement does not add up to stronger PP3 evidence: ClinGen's calibration recommends one tool, used genome-wide, for PP3/BP4, so picking whichever tool scores a variant highest would bias the result.

### SpliceAI

Deep learning model predicting **splice-altering** variants. Scores four types of splice disruption: acceptor gain (AG), acceptor loss (AL), donor gain (DG), donor loss (DL).

| Max Delta Score | Interpretation |
|---|---|
| < 0.2 | Unlikely to affect splicing |
| 0.2-0.5 | May affect splicing — investigate |
| 0.5-0.8 | Likely affects splicing |
| > 0.8 | Strong evidence of splice disruption |

**When to use SpliceAI:** VEP already flags canonical splice site variants (GT/AG dinucleotides). SpliceAI catches **cryptic** splice variants — intronic or exonic variants that create new splice sites or disrupt existing ones through more subtle mechanisms.

### gnomAD Gene Constraint (Step 23 summary)

These are per-gene metrics (not per-variant) added to the clinical filter summary TSV. They measure how intolerant a gene is to different types of mutations.

| Metric | Threshold | Meaning |
|---|---|---|
| LOEUF < 0.35 | Constrained for loss-of-function | Gene is intolerant to LoF mutations — a HIGH impact variant here is more likely to cause disease |
| pLI >= 0.9 | Loss-of-function intolerant | Same as LOEUF but older metric. LOEUF is preferred. |
| mis_Z > 3.09 | Constrained for missense | Gene is intolerant to missense mutations — a REVEL-high missense here is more concerning |

**Combining scores:** A rare variant (gnomAD AF < 0.01%) with CADD > 25, in a constrained gene (LOEUF < 0.35), with ClinVar pathogenic classification, is a high-confidence pathogenic finding. Any one of these alone is insufficient.

---

### Quick Variant Filtering Recipes

Copy-paste these commands to extract the most clinically relevant variants. They read step 13's output, `${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz`, with `bcftools +split-vep` from the pipeline's pinned bcftools image, which picks each value out of the `CSQ` field by name (the same way step 23 does). `-s worst` keeps the most severe consequence of each variant; `-d` keeps every transcript.

```bash
# Run from the repository root, with GENOME_DIR and SAMPLE set
source versions.env
bcf() { docker run --rm -i -v "${GENOME_DIR}:/genome" -w /genome "${BCFTOOLS_IMAGE}" bcftools "$@"; }
VEP_VCF="${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz"   # relative to GENOME_DIR

# Rare: below 0.1% in gnomAD exomes and genomes, or absent from them
RARE='(gnomADe_AF="." || gnomADe_AF<0.001) && (gnomADg_AF="." || gnomADg_AF<0.001)'
FIELDS='%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%gnomADe_AF\t%gnomADg_AF\t[%GT]\n'

# 1. Homozygous high-impact variants (stop gained, frameshift, splice donor/acceptor)
bcf view -i 'GT="AA"' "$VEP_VCF" |
  bcf +split-vep - -s worst -i 'IMPACT="HIGH"' -f "$FIELDS"

# 2. Rare high-impact variants
bcf +split-vep "$VEP_VCF" -s worst -i "IMPACT=\"HIGH\" && ${RARE}" -f "$FIELDS"

# 3. Compound heterozygous candidates: genes with two or more rare heterozygous
#    HIGH or MODERATE variants. The phase is unknown: both can sit on the same
#    copy of the gene, so this is a list to curate, not a finding.
bcf view -i 'GT="het"' "$VEP_VCF" |
  bcf +split-vep - -s worst -i "(IMPACT=\"HIGH\" || IMPACT=\"MODERATE\") && ${RARE}" -f '%SYMBOL\n' |
  sort | uniq -c | awk '$1 >= 2' | sort -rn

# 4. HIGH or MODERATE variants in genes you name, matched on the exact symbol.
#    The list here is a few cancer genes for illustration; CPSR (step 17) already
#    reports the full ACMG secondary-findings list. For pharmacogenes, use
#    GENES="CYP2D6 CYP2C19 CYP2C9 DPYD UGT1A1 SLCO1B1 TPMT NUDT15".
GENES="BRCA1 BRCA2 MLH1 MSH2 MSH6 PMS2 APC MUTYH TP53"
bcf +split-vep "$VEP_VCF" -d -i 'IMPACT="HIGH" || IMPACT="MODERATE"' -f "$FIELDS" |
  awk -F'\t' -v genes="$GENES" 'BEGIN { n = split(genes, g, " "); for (i = 1; i <= n; i++) want[g[i]] = 1 } $5 in want' |
  sort -u

# 5. Rare missense variants that SIFT and PolyPhen both call damaging
bcf +split-vep "$VEP_VCF" -s worst \
  -i "Consequence~\"missense_variant\" && SIFT~\"^deleterious[(]\" && PolyPhen~\"^probably_damaging\" && ${RARE}" \
  -f "$FIELDS"

# 6. Missense variants absent from both gnomAD sets (novel or ultra-rare)
bcf +split-vep "$VEP_VCF" -s worst \
  -i 'Consequence~"missense_variant" && gnomADe_AF="." && gnomADg_AF="."' -f "$FIELDS"
```

If the VCF came from a VEP run without gnomAD frequencies, `+split-vep` stops with `the tag "gnomADe_AF" is not defined`: drop the `RARE` term and the two gnomAD columns from `FIELDS`. PharmCAT (step 7) misses some alleles, and CYP2D6 needs the BAM-based callers (steps 21 and 32), so recipe 4 with the pharmacogene list is a cross-check, not a replacement.

**Important:** These are starting points, not definitive screens. Any interesting finding should be cross-referenced with ClinVar and ideally confirmed by a second method (Sanger sequencing or a clinical lab).
