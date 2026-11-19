# Structural Variants, Copy Number and Repeats

Part of [Interpreting Your Results](../interpreting-results.md). Read [Before You Look](../interpreting-results.md#before-you-look-what-the-pipeline-can-tell-you) first.

## Structural Variants (Steps 4, 5, 15, 18, 19, 22)

### What Are Structural Variants?

Unlike SNPs (single letter changes), structural variants are large rearrangements:
- **Deletions (DEL):** A chunk of DNA is missing
- **Duplications (DUP):** A chunk is copied extra times
- **Inversions (INV):** A chunk is flipped backwards
- **Translocations (BND):** A chunk moved to a different chromosome
- **Insertions (INS):** New DNA inserted

### A Note About BND (Breakend) Calls

If you run Manta or Delly, you will see many **BND** calls — often hundreds or thousands. BND indicates a "breakend" where one end of a read pair maps to a different chromosome or a distant location. This sounds alarming ("translocation!") but:

- **Most BND calls are artifacts** of repetitive regions, segmental duplications, or mobile elements
- A typical genome has 1,000-3,000 BND calls from Manta and 5,000+ from Delly
- **Fewer than 5 are likely real** inter-chromosomal translocations in a healthy genome
- BND calls require **multi-caller support** (called by both Manta and Delly at overlapping breakpoints) to be considered credible
- Unless a BND disrupts a known disease gene AND is confirmed by a second caller, it can be safely ignored

### How Many Is Normal?

A typical human genome has:
- ~5,000-10,000 structural variants total
- Most are in non-coding regions and harmless
- ~5-20 may affect genes
- 0-2 may be clinically significant

### Which Callers to Trust?

If you ran multiple SV callers:
- **Called by 2+ callers (Manta + Delly, or Manta + CNVpytor):** Lower false-positive rate
- **Called by 1 caller only:** Lower confidence, may be false positive
- **duphold DHFFC < 0.7 for deletions:** High confidence (depth drops as expected)
- **duphold DHBFC > 1.3 for duplications:** High confidence (depth rises as expected)

### SV Consensus (Step 22)

`${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf.gz` keeps the SVs that two or more callers found (Manta, Delly, CNVpytor and, when they ran, Sniffles2 and TIDDIT): SURVIVOR pairs two calls when both breakpoints lie within 1 kb of each other and their type and strands agree. `SUPP` says how many callers found each one and `SUPP_VEC` which. GRIDSS is left out (its breakend records never match). It is the short list to read first, but it drops real SVs that only one caller found; see [step 22](../22-survivor-merge.md#limitations).

### AnnotSV Output

The AnnotSV TSV (step 5) adds clinical annotations to each SV:
- `ACMG_class`: 1 (benign) to 5 (pathogenic)
- `Overlapped_CDS_percent`: How much of a gene is affected
- Focus on SVs with `ACMG_class` 4 or 5 that overlap known disease genes


---

## STR Expansions (Step 9)

### What Are Repeat Expansions?

Some regions of DNA have short sequences repeated many times (e.g., CAG CAG CAG...). When the number of repeats exceeds a threshold, it can cause disease.

### How to Read ExpansionHunter Output

The output VCF lists each tested locus with the number of repeats found. Key loci:

| Locus | Gene | Normal | Intermediate / Premutation | Pathogenic | Disease |
|---|---|---|---|---|---|
| HTT | HTT | <=26 | 27-35 (mutable normal); 36-39 (reduced penetrance) | >=40 (full penetrance) | Huntington's disease |
| FMR1 | FMR1 | <45 | 45-54 (intermediate); 55-200 (premutation) | >200 | Fragile X syndrome |
| ATXN1 | ATXN1 | <33 | — | >39 | Spinocerebellar ataxia 1 |
| C9orf72 | C9orf72 | <24 | 24-30 (gray zone, lab cutoffs vary) | Typically hundreds-thousands | ALS / Frontotemporal dementia |
| DMPK | DMPK | <35 | 35-49 (premutation) | >=50 | Myotonic dystrophy type 1 |

**HTT 36-39 repeats (reduced penetrance):** Individuals in this range may or may not develop Huntington's disease. The risk increases with repeat length but is not certain. GeneReviews classifies >=40 as full penetrance and 36-39 as reduced penetrance. Alleles of 27-35 ("mutable normal") do not cause disease but may expand in offspring.

**C9orf72 gray zone:** The exact pathogenic threshold for C9orf72 is not established. Laboratory cutoffs are discordant (JNNP 2021 review). Clearly pathogenic expansions are typically hundreds to thousands of repeats. Short-read WGS has limited ability to size very large expansions accurately.

**FMR1 intermediate zone (45-54 repeats):** Not affected, but repeats may expand in offspring. Carriers should receive genetic counseling. Premutation (55-200) carries risk of FXTAS (males >50) and FXPOI.

The report prints no all-clear verdict for these loci. It lists every locus in the ExpansionHunter output, the five above first, with its repeat count, Stranger's status and the catalog's normal maximum (`STR_NORMAL_MAX`), and, when step 9b (Stranger) ran, every locus Stranger marks `pre_mutation` or `full_mutation`; without step 9b, compare the counts with the table yourself. It also says:

- **Not sized.** A locus whose call rests on flanking reads only prints "size not determined (flanking reads only)"; one with in-repeat reads only (longer than a read) prints its count as at least the lower end of ExpansionHunter's interval. An allele counts as in-repeat when its `SO` says `INREPEAT`, or when ExpansionHunter found no spanning and no flanking read for it but at least one in-repeat read (`ADSP` and `ADFL` 0, `ADIR` above 0), whatever its `SO` says; one with fewer than 2 spanning reads, or an interval that reaches the catalog's pathogenic minimum, is marked low confidence. Such a locus is listed apart and never shown as normal, because an unsized repeat can hide an expansion.
- **Intermediate range.** A sized count inside the HTT, FMR1 or DMPK intermediate range of the table above prints the range, a note that it does not cause disease but may expand in offspring, and the GeneReviews chapter. The report reads these ranges from `bin/str_intermediate_ranges.tsv`, and `rust/crates/og-setup/tests/genome_str_ranges.rs` checks that the table above prints the same ones.
- **Reduced penetrance.** Stranger grades HTT with one pathogenic minimum, 36 on the STRchive catalog, so a count of 36 to 39 reads as plain `full_mutation`. It is the reduced-penetrance range of the table above. `summary.json` marks a sized HTT allele of 36 to 39 with `penetrance` (`36-39`), `penetrance_label` and `penetrance_citation` (the GeneReviews chapter) on the locus, and on its entry in the flagged list. The range comes from `assets/str_penetrance_ranges.tsv`. An allele that was not sized gets no such mark: its point estimate is not a count. Stranger's status stays as Stranger wrote it. The report does not print this mark yet; it is in `summary.json` only, and only from step 24 (`run-all.sh`): the Nextflow report does not stage the table yet, so a Nextflow run has no such mark.
- **RFC1.** A `full_mutation` flag at RFC1 carries the motif verdict of step 9b ([RFC1 note](../09b-stranger.md#notes)).

A count under every threshold is not a guarantee: repeat-expansion calls from short reads can be wrong at some loci, and a long expansion cannot be sized. Intermediate-range results should be discussed with a genetic counselor.

### Stranger Annotation (Step 9b)

Step 9b adds a `STR_STATUS` field to each locus in the VCF so you do not need to look up thresholds manually:

| STR_STATUS | Meaning |
|---|---|
| `normal` | Repeat count is within the established normal range |
| `pre_mutation` | Elevated repeat count; not currently disease-causing but carries risk of expansion in offspring or late-onset carrier effects (e.g. FXTAS, FXPOI for FMR1) |
| `full_mutation` | Repeat count exceeds the pathogenic threshold for this locus |

The `Disease`, `InheritanceMode`, `SourceDisplay`, `STR_NORMAL_MAX`, and `STR_PATHOLOGIC_MIN` INFO fields give the clinical context for each locus (a field the catalog has no value for prints as `.`). Use a VCF viewer or `bcftools query` to extract them:

```bash
bcftools query -f '%INFO/STR_STATUS\t%INFO/Disease\t%INFO/InheritanceMode\t%INFO/SourceDisplay\t%INFO/STR_NORMAL_MAX\t%INFO/STR_PATHOLOGIC_MIN\n' \
    expansion_hunter/<sample>_eh_stranger.vcf
```

`full_mutation` at any locus warrants follow-up with a clinical geneticist. Short-read WGS has limited sizing accuracy for very large expansions (>150 repeats), so confirmatory testing may be recommended.

---

## CNVpytor Results (Step 18)

CNVpytor detects **copy number variants** using read depth analysis — complementary to Manta's paired-end/split-read approach.

**Where to look:** `${SAMPLE}/cnvpytor/${SAMPLE}_cnvs.txt`

**Format:** Each line has: type, region, size, normalized_RD, e-value1, e-value2, e-value3, e-value4, q0

**What to expect:**
- 3,000-4,000 total CNVs (mostly deletions)
- 1,500-2,000 significant (e-value < 0.01)
- Calls at chromosome starts (chr1:1-10000) are telomeric artifacts — ignore them

**Filtering:**
```bash
# Significant CNVs only (e-value < 0.01)
awk '$5 < 0.01' ${SAMPLE}_cnvs.txt

# Large deletions (>100kb, potentially clinically relevant)
awk '$1 == "deletion" && $3 > 100000 && $5 < 0.01' ${SAMPLE}_cnvs.txt

# Large duplications
awk '$1 == "duplication" && $3 > 100000 && $5 < 0.01' ${SAMPLE}_cnvs.txt
```

**Multi-caller overlap:** CNVs found by both Manta AND CNVpytor have lower false-positive rates. Cross-reference by checking if the same genomic region appears in both output files.

---

## Delly Results (Step 19)

Delly is a third structural variant caller, detecting deletions, duplications, inversions, and translocations.

**Where to look:** `${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz`

**Quick summary:**
```bash
# Count SVs by type
bcftools query -f '%INFO/SVTYPE\n' ${SAMPLE}_sv.vcf.gz | sort | uniq -c

# Filter PASS variants only
bcftools view -f PASS ${SAMPLE}_sv.vcf.gz | grep -cv '^#'
```

**What to expect:**
- 5,000-15,000 total SV calls
- Most are small deletions (<1kb)
- PASS filter reduces count significantly

**Multi-caller overlap:** SVs detected by Manta + Delly + CNVpytor (or any 2 of 3) have substantially lower false-positive rates. Single-caller calls, especially large ones, should be viewed with caution.

---

## SMN1 and SMN2 Copy Number (Step 35, opt-in)

**Where to look:** `${SAMPLE}/paralogs/${SAMPLE}_smn_copy_number.tsv`, one row per stretch of the SMN1/SMN2 locus with the copy number of the two genes together (agCN) and of each (psCN), each with a quality.

One SMN1 copy suggests SMA carrier status; two do not rule it out (two copies on one chromosome, none on the other). Trust a value only with quality 20 or more and filter `PASS`, and confirm anything that matters with a clinical SMN1 test. [Step 35](../35-paralogs.md#interpreting-results) explains the columns.

**The silent-carrier marker:** `${SAMPLE}/paralogs/${SAMPLE}_smn1_marker.tsv` reads SMN1 g.27134T>G (rs143838139), and both reports and the live person page print it right under the SMN1 copies: `present`, `uncertain`, `not seen`, `too little depth` or `not assessed`. With two SMN1 copies, `present` raises the chance that both sit on one chromosome (2+0, a silent carrier) and `not seen` lowers it; it proves neither, and the figures depend on ancestry. A carrier test in a clinical laboratory confirms it. [The SMN1 silent-carrier marker](../35-paralogs.md#the-smn1-silent-carrier-marker) gives the thresholds.

## Lp(a) and the LPA KIV-2 Count (Step 46, opt-in)

**Where to look:** `${SAMPLE}/lpa/${SAMPLE}_lpa_kiv2.tsv`, the KIV-2 repeat count over both chromosomes with its 95% interval and class, and `${SAMPLE}/lpa/${SAMPLE}_lpa_kiv2_4925.tsv`, the KIV-2 4925G>A splice variant. Both reports put the count beside the two LPA risk variants of step 42 (rs10455872, rs3798220), and print the 4925G>A line right under the count: `present` or `uncertain` with the reads that carry it and the copies, `not seen`, `too little depth` or `not assessed`.

Fewer repeats or a risk allele point to a higher Lp(a); a 4925G>A copy lowers the Lp(a) made from the chromosome that carries it. Neither tells which chromosome carries what, and `not seen` is not proof of absence at low depth. The count is a research estimate: the Lp(a) blood test is the measure, once in every adult. [Step 46](../46-lpa.md) explains the columns.

## Mobile Element Insertions (Step 50, opt-in)

**Where to look:** `${SAMPLE}/mei/${SAMPLE}_mei_sites.tsv` for the three known disease insertions (SVA insertions in FKTN, TAF1 and BBS1), and `${SAMPLE}/mei/${SAMPLE}_mei.tsv` for every Alu, L1 and SVA insertion MEGAnE called, with one or two copies.

Expect many calls: a genome carries about 1,200 such insertions, nearly all common and harmless. Read the catalog first. **Seen** means an insertion of the same element was called inside the site's window; one copy of the FKTN or BBS1 insertion is carrier status for a recessive condition, and the TAF1 insertion on a male's single X is what causes X-linked dystonia-parkinsonism. **Not seen** is not proof of absence: short reads find an SVA's edges, not its body. **Not callable** means too few reads there (under 10x) to tell. Confirm anything that matters in a clinical laboratory. [Step 50](../50-mei.md) explains the columns.
