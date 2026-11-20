# Telomeres and Somatic Variants

Part of [Interpreting Your Results](../interpreting-results.md). Read [Before You Look](../interpreting-results.md#before-you-look-what-the-pipeline-can-tell-you) first.

## Telomere Length (Step 10)

### What It Means

Telomeres are protective caps at the ends of chromosomes that shorten with cell division. TelomereHunter measures `tel_content`: intratelomeric reads per million reads with 48-52% GC. It is a relative telomere content, not a telomere length, and the report labels it "Telomere content (relative)".

### How to Interpret

- **Higher `tel_content` = more telomeric reads** (generally correlates with longer telomeres)
- **Lower `tel_content` = fewer telomeric reads** (generally correlates with shorter telomeres)
- There is no universal "normal" range — compare between samples of similar age, sequenced on the same platform
- Typical `tel_content` for 30X WGS: 300-800 (varies by sequencing platform and coverage)

### Limitations

- **Not a clinical telomere length measurement.** TelomereHunter was developed for cancer genome analysis, not as a validated healthy-population aging assay
- **Not an age.** Telomere length follows age only across a population. For one person it is a rough estimate, and it is not an established standalone risk marker (Vaiserman and Krasnienkov 2021, Frontiers in Genetics 11:630186, https://doi.org/10.3389/fgene.2020.630186)
- Short-read WGS systematically underestimates telomere length compared to dedicated assays (TRF, FlowFISH)
- Useful only for **relative comparisons** between samples run on the same platform — not absolute measurements and not individual health predictions

---

## Somatic Variants (Step 29) [EXPERIMENTAL]

Mutect2 in tumor-only mode looks for somatic mutations -- variants acquired during your lifetime rather than inherited. Consumer WGS is usually made from saliva or a cheek swab, a mix of cheek cells and white blood cells; the main category of interest is **clonal hematopoiesis (CHIP)**, which lives in the blood cells.

**Where to look:** `${SAMPLE}/somatic/${SAMPLE}_somatic_filtered.vcf.gz`

**Key field:** `AF` (allele fraction) indicates the clone size:

| AF Range | Likely Source |
|---|---|
| 0.45-0.55 | Heterozygous germline (false positive) |
| ~1.0 | Homozygous germline (false positive) |
| 0.10-0.40 | Could be a large somatic clone, mosaic, or noisy germline |
| below 0.10 | At 30X this is one to three reads: mostly noise, see below |

**What to expect:**
- Thousands of PASS calls in a healthy individual -- the vast majority are germline false positives
- At 30X a variant in 2% of the reads has less than one supporting read on average, and one in 10% about three. Only variants at roughly 10% allele fraction or more can be told apart from noise. That is an approximate rule about reads, not a clone size: a heterozygous variant at 10% allele fraction sits in about 20% of the sampled cells, and copy number and the mix of cell types in the sample shift it. Most CHIP (defined from 2% allele fraction) is invisible here, and a clean result does not rule it out
- Without a matched normal sample, germline variants that are rare in gnomAD will often pass all filters

**CHIP genes to check:** DNMT3A, TET2, ASXL1, TP53, JAK2, SF3B1, SRSF2, PPM1D, CBL. CHIP prevalence increases with age and is associated with elevated cardiovascular risk and risk of hematologic malignancies.

**Important caveats:**
- This step has a much higher false positive rate than any other step in the pipeline
- Do not interpret PASS variants as confirmed somatic without cross-referencing with the germline VCF (step 3) and gnomAD frequencies (step 13)
- If a variant is also called at ~50% AF by DeepVariant in step 3, it is almost certainly germline
