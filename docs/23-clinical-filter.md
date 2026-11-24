# Step 23: Clinical Variant Filter

## What This Does

Extracts the small subset of clinically interesting variants from your VEP-annotated VCF. Instead of manually searching through 4-5 million variants, this step produces a focused list of a few hundred variants that are rare and functionally impactful, plus the known ClinVar pathogenic ones.

## Why

The biggest challenge after running VEP annotation is: "I have millions of variants, what do I look at?" This step solves that by applying conservative filters to surface the variants most likely to be medically relevant.

## Tool

bcftools + `bcftools +split-vep` plugin (parses VEP CSQ fields structurally — no grep)

## Docker Image

- `BCFTOOLS_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Input

- The step 30 output `${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_annotated.vcf.gz` when it is newer than the VEP output, else the VEP output of step 13 (`${SAMPLE}_vep.vcf.gz`, or `${SAMPLE}_vep.vcf`, which is compressed first). A derived file older than its source is ignored with a notice, so a re-run of step 13 is never hidden behind an old copy.
- Optional: the step 6 hits `${GENOME_DIR}/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.vcf` and the gnomAD v4.1 constraint table `${GENOME_DIR}/annotations/gnomad_v4.1_constraint.tsv`.

## Command

```bash
./scripts/23-clinical-filter.sh your_name
```

## What Gets Filtered

Every tier starts from the PASS records. "Rare" means VEP's `MAX_AF` (the highest allele frequency in any 1000 Genomes, gnomAD exome or gnomAD genome population) is below 1%. Without `MAX_AF` the step uses `gnomADe_AF` and `gnomADg_AF`: rare when the values it has are below 1%. A variant common in gnomAD genomes but absent from the exomes is therefore not rare. When the VEP output has none of these fields (VEP run without `--everything`, `--max_af` or `--af_gnomadg`), no tier is filtered by frequency, every MODERATE variant is kept, and the step prints a notice saying so.

A variant with no frequency value at all is **frequency unknown**, not rare. VEP gives a variant a frequency only when it finds the same allele in gnomAD, written the same way. A multi-base indel written right-shifted, or inside a multiallelic record, finds nothing and gets no value even when it is common. Step 13 normalises its input for that reason, but a VCF annotated elsewhere may not be. HIGH and MODERATE variants of unknown frequency go to their own tier (`_freq_unknown.vcf.gz`) and their own count, never into the rare HIGH or MODERATE counts. The CADD, SpliceAI and REVEL/AlphaMissense tiers keep both rare and unknown-frequency variants; the summary's `FREQ_TIER` column says which.

The gene, impact and consequence of a variant come from one transcript: its MANE Select transcript (VEP's `MANE_SELECT`, written with `--everything`), else the canonical one (`CANONICAL=YES`), else the most severe over every transcript. A variant with a frameshift on a minor transcript and an intronic consequence on the MANE Select one is ranked intronic, so it is not in the rare HIGH tier. The summary keeps the worst consequence beside the ranked one. `og vcf csq-rank` (OwnGenome's Rust binary, run in the bcftools image) applies this rule, and the Nextflow module and the reports use the same one. Step 31 keeps the worst consequence for its tiers and shows the MANE Select one beside it.

### Genotype check

Each variant's call is checked against the sample's `GQ`, `DP` and `AD`. A call is **weak** when GQ is below 30 or DP below 10 (step 6's floor), when a heterozygote's alt reads are outside 20% to 80% of the depth, or when a homozygous alt call's alt reads are below 80%. A 1/1 call at GQ 11 with 12 reference and 5 alt reads is weak: the reads support a heterozygote or no call, not a homozygote. Weak calls stay in their tier and are marked (`CALL_CHECK`, `CALL_WHY`); the step prints how many there are, and the reports list them apart from the counted findings. The floors are settable: `CLINICAL_MIN_GQ`, `CLINICAL_MIN_DP`, `CLINICAL_HET_AB_MIN`, `CLINICAL_HET_AB_MAX`, `CLINICAL_HOM_AB_MIN` for the script, and `ext.min_gq`, `ext.min_dp`, `ext.het_ab_min`, `ext.het_ab_max`, `ext.hom_ab_min` on the `CLINICAL_FILTER` process. The clinical VCF records the floors it used (`##PGP_CALL_FLOOR`), and the reports read them from there.

### Rare HIGH impact
Stop-gained, frameshift, splice donor/acceptor, start-lost, on the ranked transcript.

### Rare MODERATE impact
Missense variants and in-frame insertions/deletions, on the ranked transcript.

### HIGH or MODERATE impact, frequency unknown
The same impacts with no frequency value. Their own tier, not rare (see above).

### ClinVar pathogenic/likely pathogenic, at any frequency
- Preferred source: the step 6 hits file, built from the ClinVar file in `clinvar/` that `setup.sh` refreshes. The tier holds the VEP records whose allele is a hit: the same CHROM, POS, REF and ALT.
- Another allele at a hit's position is not in the tier. A call of T>A,C splits into two records at one position; when ClinVar classified only T>A, the T>C record stays out.
- Step 6 and step 13 both split and left-align the sample against the same reference, so one allele is written the same way in the hits and in the VEP output, one ALT per record. A VEP record with several ALTs (an input step 13 did not split) is kept when any of its ALTs is a hit.
- A combined hit of step 6 (two neighbouring records that make one allele, `PGP_JOINED_FROM` in the hits) puts both of its parts in the tier. VEP annotated the parts, not the joined allele. The tier's count then counts each part as one record.
- Next, when step 6 has not run: `ClinVar_CLNSIG`, the same ClinVar file step 13 annotated with `--custom` (pathogenic or likely pathogenic, not conflicting). It follows a ClinVar refresh once step 13 runs again.
- Last: VEP's `CLIN_SIG`, from the VEP cache release, so a ClinVar refresh never reaches it. The step says which source it used.
- A common pathogenic allele (for example HFE p.C282Y) stays: this tier has no frequency filter.

### Rare high CADD (step 30)
CADD PHRED >= 20 for variants that are not HIGH or MODERATE on the ranked transcript. Rare or of unknown frequency.

### Rare cryptic splice (step 30)
A SpliceAI delta score >= 0.2 for any gene of the value. SpliceAI writes one entry per gene, joined by commas; every entry is tested.

### Rare deleterious missense (step 30)
REVEL >= 0.644 (ClinGen's PP3 Supporting threshold) or AlphaMissense >= 0.564 (AlphaMissense's own likely_pathogenic class boundary, not an ACMG evidence level).

## Output

| File | Contents |
|---|---|
| `${SAMPLE}_clinical.vcf.gz` | All tiers merged |
| `${SAMPLE}_clinical_summary.tsv` | One row per variant: `CHROM`, `POS`, `REF`, `ALT`, `GT`, `IMPACT`, `GENE`, `Consequence`, `MAX_AF`, `CADD_PHRED`, `REVEL`, `AM_CLASS`, `TRANSCRIPT`, `WORST_IMPACT`, `WORST_Consequence`, `FREQ_TIER`, `GQ`, `DP`, `AD`, `CALL_CHECK`, `CALL_WHY`, and with the constraint table `LOEUF`, `pLI`, `mis_z` |
| `${SAMPLE}_high_impact.vcf.gz` | Rare HIGH impact |
| `${SAMPLE}_rare_moderate.vcf.gz` | Rare MODERATE impact |
| `${SAMPLE}_freq_unknown.vcf.gz` | HIGH or MODERATE impact of unknown frequency (when the VEP output has a frequency field) |
| `${SAMPLE}_clinvar_pathogenic.vcf.gz` | ClinVar P/LP (when step 6 hits, `ClinVar_CLNSIG` or `CLIN_SIG` exist) |
| `${SAMPLE}_cadd_high.vcf.gz`, `${SAMPLE}_spliceai_high.vcf.gz`, `${SAMPLE}_missense_deleterious.vcf.gz` | The step 30 tiers, when their scores exist |

`IMPACT`, `GENE` and `Consequence` are those of the ranked transcript; `TRANSCRIPT` says which set it came from (`MANE`, `CANONICAL` or `ANY`), and `WORST_IMPACT` and `WORST_Consequence` give the most severe consequence over every transcript. `GENE` is `.` for an intergenic variant. `FREQ_TIER` is `rare`, `unknown`, `common` (a ClinVar record) or `none` (no frequency field). `CALL_CHECK` is `ok`, `weak` or `unknown` (no GQ, DP or AD to judge), and `CALL_WHY` names what made a call weak, for example `GQ11,AB0.29`. A score, frequency or FORMAT field the input does not carry is written as `.`. The records of every tier carry the same marks as `PGP_*` INFO tags in the VCFs.

The constraint columns come from `og table constraint-join`, the loader step 31 and the Nextflow slivar module run too: only canonical transcripts count, the Ensembl row wins over the RefSeq one, and `mis_z` is gnomAD v4.1's `mis.z_score`. When the table is present and rows carry gene symbols but not one matches it, the step fails instead of writing `.` everywhere.

The Nextflow `CLINICAL_FILTER` module applies the same tiers, ranking, frequency tier and genotype check (it runs `og vcf csq-rank` too) and writes the same summary columns. It takes the ClinVar tier from `ClinVar_CLNSIG` (the `--clinvar` file, which the VEP module adds with `--custom`), else from VEP's `CLIN_SIG`, and adds no constraint columns.

### Secondary-findings genes (ACMG SF v3.3)

Both reports list, beside the clinical filter's counts, the secondary findings in the 84 genes of the ACMG SF v3.3 list (Lee et al., Genet Med 2025;27(8):101454; the genes in which the ACMG recommends reporting pathogenic variants found by chance, such as BRCA1, LDLR, MYBPC3 or TTN): the step 6 ClinVar hits in those genes, the clinical filter's records with a HIGH-impact consequence on the ranked transcript in them, and step 5's AnnotSV class 4 and 5 structural variants in the genes whose rule allows them. Each gene's rule from Table 1 of the list is a row of `assets/acmg_sf.tsv`: one variant in a dominant gene, two copies hit in a recessive one (MUTYH, BTD, ATP7B and others: a homozygote, two variants that read-backed phase puts on different copies, or two variants of unknown phase, which `summary.json` marks `"pair": "possibly biallelic"`; two variants phased on one copy are listed apart, a combined allele of step 6 counts once with its parts, and a record with no alt copy called counts nothing), HFE only as homozygous p.C282Y, GLA and OTC also in a heterozygous female, ABCD1 hemizygous, homozygous or with two heterozygous variants. A variant the rule does not report is listed apart with the reason; a pseudodeficiency allele is not counted, and a BTD finding with the p.Asp444His hypomorph says partial, not profound. Of the reported HIGH-impact records the reports count as "rare HIGH impact" only the rare ones with a sound call; a weak call or an unknown frequency is listed with its mark and counted apart, and a common one is left out. The list version is `ACMG_SF_VERSION` in `og summary collect` (`rust/crates/og-summary`). It is a list to review with a clinician: a rare HIGH-impact variant is not a ClinVar classification. [Interpreting results](interpreting-results.md#carrier-status-is-not-disease) gives every rule. Step 17 (CPSR) reports the ACMG secondary findings with its own classification.

summary.json holds the card twice. `sections.clinical.data.acmg_sf` is the card as both reports read it. It is there whenever this step's summary table or clinical VCF can be read, also one left by an earlier run (then, when this step did not run ok, the clinical section is stale). `sections.acmg_sf` is the same card as a section of step 6: it is written whenever step 6's hits are there, so a run without the VEP cache still gets the ClinVar and SV findings. Without this step's clinical VCF it has `high_impact_checked: false` and a `high_impact_note` that names why (this step skipped, with the run's reason, failed, or its VCF unreadable) and what that can miss: a predicted loss-of-function variant that is not in ClinVar, alone or as the second hit of a recessive gene. With the clinical VCF the two cards are equal, except when that VCF is from an earlier run (this step not ok in this run): then `sections.acmg_sf` adds a `high_impact_note` with the earlier VCF's date. A clinical VCF that cannot be read, or whose HIGH-impact records cannot be parsed, makes the clinical section's card unavailable, while `sections.acmg_sf` keeps its ClinVar rows with `high_impact_checked: false` and says the VCF could not be read. The two cards then differ.

The `table_text` column of `assets/acmg_sf.tsv` quotes each gene's rule from Table 1 of the paper (doi:10.1016/j.gim.2025.101454), such as "All P and LP" or "P and LP (2 variants)". summary.json gives that text as the rule of each finding, and the reports use it in the reason a variant is not reported. These are short cited quotations of a factual reporting rule: six phrases, the longest 44 characters. "(See text)" after a rule points to the paper's section on the dilated cardiomyopathy (DCM) genes. The `note` column says what the table's footnotes add, in our own words.

A ClinVar hit reaches its gene's rule through GENEINFO, which names the gene, not the disease. A RET variant that ClinVar asserts only for Hirschsprung disease would fill the MEN2 slot like a MEN2 variant. So each ClinVar finding of the card keeps the condition ClinVar names for the allele, which step 6 copies from `CLNDN` and `CLNDISDB`: `clinvar_condition` (the names, without "not provided" and "not specified") and `clinvar_condition_ids` (the database ids, such as `OMIM:142623`, also those of a condition `CLNDISDB` gives and `CLNDN` does not name). Beside them, `acmg_phenotype` is the disease Table 1 lists the gene for (the `phenotype` column). `condition_not_acmg` is `true` when every condition ClinVar names is another disease. A condition is the ACMG disease when one of its ids is a key of the gene's `phenotype_keys` column (Table 1's MIM disorder numbers first), or when a key phrase starts a word of its name. "not provided", "not specified", panel names such as "Hereditary cancer-predisposing syndrome", "Inborn genetic diseases" or "Cardiac arrhythmia", and a name that only names the hit's own gene ("RYR1-related disorder" on a RYR1 hit; also "-related disorders", "-related condition(s)" and "-related cancer predisposition") count neither way. Another gene's such name ("FBXO11-related disorders" on an MSH6 hit) is another disease. A condition with ids and no name never sets the mark either. A marked hit stays a finding: the mark tells the reader to check which disease the variant was reported for, it never drops a row. A hits VCF without `CLNDN` and `CLNDISDB` (written before step 6 kept them) gives none of these four fields. Each carrier list row keeps `clinvar_condition` and `clinvar_condition_ids` too, beside the list's own `condition` for the gene.

## Runtime

2.5 to 6.5 minutes, measured in three 30-40x WGS runs ([runtime](hardware-requirements.md#runtime-per-step)) (I/O-bound, reading the large VEP VCF)

## How to Use the Output

### Quick look at the summary

```bash
# View the most important variants
column -t ${GENOME_DIR}/${SAMPLE}/clinical/${SAMPLE}_clinical_summary.tsv | head -20
```

### Cross-reference with ClinVar

```bash
source versions.env   # from the repository root
# Find which clinical variants are also in ClinVar
docker run --rm -v "${GENOME_DIR}:/genome" "${BCFTOOLS_IMAGE}" \
  bcftools isec -n=2 -w1 \
    /genome/${SAMPLE}/clinical/${SAMPLE}_clinical.vcf.gz \
    /genome/clinvar/clinvar.vcf.gz \
    -Oz -o /genome/${SAMPLE}/clinical/${SAMPLE}_clinical_clinvar.vcf.gz
```

### Load in a genome browser

The `_clinical.vcf.gz` file is small enough to load in [IGV Web](https://igv.org/app/) or [gene.iobio](https://gene.iobio.io/) for visual inspection.

## Limitations

- This is a **computational filter**, not a clinical interpretation
- Some pathogenic variants are LOW impact (e.g., synonymous variants affecting splicing, regulatory variants) and will be missed by this filter
- Frequency filtering depends on VEP having written `MAX_AF` or the gnomAD fields, and on VEP finding the allele: an unmatched one is "frequency unknown"
- Always cross-reference findings with ClinVar and consult a professional for clinical decisions

## Notes

- No additional Docker images required — uses the same bcftools image as other steps
- Uses `bcftools +split-vep` to parse VEP's pipe-delimited CSQ annotation structurally (not grep)
- The VEP VCF is compressed and indexed automatically if needed
- PASS filter is applied to exclude low-quality variant calls
- Available CSQ subfields and INFO scores are detected from the header; a tier whose data is missing is skipped and the step says so
