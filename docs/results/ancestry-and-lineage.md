# Ancestry, Lineage and Polygenic Scores

Part of [Interpreting Your Results](../interpreting-results.md). Read [Before You Look](../interpreting-results.md#before-you-look-what-the-pipeline-can-tell-you) first.

## ROH Analysis (Step 11)

### What It Means

Runs of Homozygosity (ROH) are long stretches where both copies of your DNA are identical. Everyone has some ROH, but extensive ROH can indicate:
- Parental relatedness (consanguinity)
- Uniparental disomy
- Population bottleneck effects

### How to Read

Add up the autosomal segments of 5 Mb or more and compare the total with the table in [step 11](../11-roh-analysis.md#total-roh-and-parental-relationship), which gives the expected total for each parental relationship (about 45 Mb for second cousins, 180 Mb for first cousins). In short:

- **No segments of 5 Mb or more, or only a few:** no sign that your parents are related
- **Many small segments (1-5 Mb):** population-level background, typical of population isolates (Ashkenazi, Finnish, etc.)
- **One chromosome carrying most of the long runs:** possible uniparental disomy rather than related parents. The reports flag a chromosome when its segments of 5 Mb or more add up to at least 20 Mb and to more than half of all such autosomal segments ([step 11](../11-roh-analysis.md#one-chromosome-carrying-most-of-the-runs-possible-uniparental-disomy))

### Centromeric Artifacts

Some apparent ROH near centromeres are artifacts of low-coverage sequencing in repetitive regions. These can be ignored.

---

## Mitochondrial Haplogroup (Step 12)

### What It Means

Your mitochondrial haplogroup traces your maternal ancestry lineage. It's determined by the specific set of variants in your mitochondrial DNA (inherited only from your mother).

### Common European Haplogroups

| Haplogroup | Origin | Notes |
|---|---|---|
| H | Western Europe | Most common in Europe (~40%) |
| U | Northern/Eastern Europe | Second most common (~15%) |
| T | Near East / Mediterranean | ~10% of Europeans |
| K | Near East | ~6% of Europeans, Ashkenazi ~30% |
| J | Near East | ~8% of Europeans |
| V | Iberian Peninsula / Scandinavia | ~4% |
| I | Near East / Europe | ~3% |

### Medical Relevance

Mitochondrial haplogroups have weak associations with some diseases (Parkinson's, diabetes, longevity), but these are population-level statistics, not individual predictions. The main clinical value is in step 20 (GATK Mutect2 mitochondrial mode), which detects disease-causing mitochondrial variants and heteroplasmy.

### Exploring Your Branch

Both reports link your haplogroup to FamilyTreeDNA Discover (free, no upload) and to YFull's tree. Discover shows the branch's age, where its testers live and the ancient DNA samples on it. For the Y haplogroup of step 37 the reports add the same two links. The free and paid routes, what each one takes and what an upload shares: [step 12](../12-mito-haplogroup.md#explore-your-branch-further) and [step 37](../37-y-haplogroup.md#explore-your-branch-further). The pipeline never uploads anything.

---

## Mitochondrial Variants (Step 20)

GATK Mutect2 in mitochondrial mode detects variants with heteroplasmy fractions — the proportion of your mitochondria carrying each variant.

**Where to look:** `${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz`

**Key field:** `AF` (allele fraction) indicates heteroplasmy level:

| AF Level | Meaning |
|---|---|
| >0.95 | Homoplasmic — fixed in all mitochondria (haplogroup-defining) |
| 0.10-0.95 | Heteroplasmic — clinically significant range |
| 0.03-0.10 | Low-level heteroplasmy — often age-related somatic |
| <0.03 | Near detection limit |

**What to expect:**
- 50-70 PASS variants total
- 25-35 homoplasmic (haplogroup variants)
- 25-35 low-level heteroplasmic (mostly <5%)
- Length variants in the D-loop poly-C tracts (m.303-315 and m.16184-16193) in most samples: the reports list them apart and do not count them as heteroplasmy

**When to investigate further:**
- Heteroplasmic variant at a known disease position: the report lists each heteroplasmic call, and each call ClinVar classes as anything but benign, with its ClinVar class and gnomAD v3.1 frequency (`mito/<sample>_chrM_annotated.tsv`). [MITOMAP](https://www.mitomap.org/) has the curated status; look it up by hand, as it refuses scripted requests
- m.3243A>G (MELAS) or m.8344A>G (MERRF) at detectable heteroplasmy levels
- Any position in MT-ATP6, MT-ND genes with AF >0.10

**Important caveats about heteroplasmy thresholds:**
- There is **no single absolute heteroplasmy threshold** that determines clinical significance. Thresholds vary by variant and by tissue (ClinGen/MSeqDR mtDNA interpretation specifications)
- **Blood and saliva underrepresent heteroplasmy** for many mitochondrial diseases. WGS from saliva, a cheek swab or blood may show lower heteroplasmy levels than affected tissues (muscle, nerve). m.3243A>G in particular shows different clinical phenotypes at very different heteroplasmy levels across tissues
- The AF values from this pipeline reflect the tissue your sample came from. A low or absent heteroplasmy level there does **not** rule out clinically significant heteroplasmy in other tissues
- For any detected pathogenic mtDNA variant, discuss with a specialist who can order tissue-specific testing if warranted

**Cross-reference:** Compare with step 12 (haplogrep3) — your homoplasmic variants should match your assigned haplogroup.

---

## Ancient Populations (Step 44, opt-in)

**Where to look:** `${SAMPLE}/ancient/${SAMPLE}_ancient.tsv`, the 20 ancient groups nearest to you on a PCA of present-day West Eurasians, with their dates, places and number of people (`n`), and `${SAMPLE}/ancient/${SAMPLE}_ancient_pca.svg`, the plot with you marked.

Closeness on the PCA is genetic similarity, not descent: being near an ancient group does not mean you descend from those people. A group with `n` of 1 or 2 has a noisy position. The axes describe West Eurasia only, so a genome with most of its ancestry from elsewhere lands at the edge and its nearest groups say little. `sites_used` in `${SAMPLE}_ancient_info.tsv` tells how many panel sites your gVCF covered; a whole genome covers nearly all of them (1,133,951 of the 1,150,106 that lift to GRCh38, for a gVCF built from the public GIAB HG002 genotypes). [Step 44](../44-ancient.md) explains the columns.

## Polygenic Risk Scores (Step 25)

**Where to look:** `${SAMPLE}/prs/${SAMPLE}_prs_summary.tsv`, one row per condition, and the PRS section of the report.

With the ancestry panel installed (`scripts/setup.sh --ancestry-panel`), each score has a **Percentile** and an **Ancestry_Group**: the 80th percentile means your score is higher than 80% of the reference samples whose genetic ancestry is most similar to yours (step 26 names that group). A percentile is still not a probability: most people at a high percentile never get the condition, and lifestyle and family history often weigh more. Percentiles of different conditions are not comparable either.

Without the panel the report says "Raw score only": the sum cannot be compared with anyone, so read nothing into its size. A score with `NA` as its sum matched under 75% of its variants and was not scored; `Matched_Pct` says how many. When the input lists variants only (a VCF without its gVCF), the sums leave out the sites where you match the reference ([step 25](../25-prs.md#hom-ref-sites-come-from-the-gvcf) explains why). Each score then keeps its raw sum but gets no percentile, `Percentile_Note` says `percentile withheld: the input lists variants only`, and the report says "Percentile withheld" in place of "Raw score only". A chip VCF lists its reference calls and keeps its percentiles. When `ChrX` says `left out`, the sample's sex was not given and the score's chrX rows are not in its sum ([chrX and the sex](../25-prs.md#chrx-and-the-sex)). The breast cancer score is for women and the prostate cancer score for men: when you declared the other sex, `Sex_Check` says `not applicable to the declared sex` and the score has no percentile. The blood pressure and lipid scores rank a genetic tendency only; a measured blood pressure or lipid panel says far more. The Alzheimer's score includes APOE: read [Before you look](../interpreting-results.md#before-you-look-what-the-pipeline-can-tell-you) first.
