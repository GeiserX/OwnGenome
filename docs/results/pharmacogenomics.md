# Pharmacogenomics and HLA

Part of [Interpreting Your Results](../interpreting-results.md). Read [Before You Look](../interpreting-results.md#before-you-look-what-the-pipeline-can-tell-you) first.

## PharmCAT Report (Step 7)

**What it tells you:** How your genes affect drug metabolism. These results are clinically relevant and should be shared with your prescribing physician.

**Where to look:** `${SAMPLE}/pharmcat/` (Nextflow) or `${SAMPLE}/vcf/` (bash scripts) — PharmCAT writes its reports there. Open the HTML report in a browser.

**Key genes to check:**

| Gene | Affects | Common Impact |
|---|---|---|
| CYP2C19 | PPIs, clopidogrel, SSRIs, voriconazole | Rapid metabolizers burn through drugs too fast |
| CYP2D6 | Codeine, tramadol, tamoxifen, many psych meds | Poor metabolizers get little pain relief from codeine and tramadol, less benefit from tamoxifen, and a toxic buildup of many psych meds; ultrarapid metabolizers risk opioid toxicity from codeine and tramadol |
| CYP2C9 | Warfarin, NSAIDs, phenytoin | Dose adjustment needed |
| DPYD | 5-fluorouracil (cancer drug) | Poor metabolizers can die from standard doses |
| SLCO1B1 | Statins (simvastatin, atorvastatin) | Increased myopathy risk |
| NAT2 | Isoniazid (TB), hydralazine | Slow acetylators have more side effects |
| UGT1A1 | Irinotecan, atazanavir | Two \*28 alleles are associated with Gilbert syndrome (elevated bilirubin) |

**What to do:** Share the PharmCAT report with your prescribing physician or pharmacist. PharmCAT is a research tool — its authors explicitly note that missing positions, unphased input, and undetected structural variation (especially CYP2D6) can affect genotype and phenotype calls. The report is a valuable starting point for pharmacogenomic-guided prescribing, but clinical confirmation may be warranted before making medication changes, especially for high-risk drugs (DPYD, CYP2D6-dependent opioids).

---

## HLA Typing (Step 8)

**Where to look:** `${SAMPLE}/hla_t1k/${SAMPLE}_hla_genotype.tsv`, two alleles per HLA gene.

The main use is drug safety: a few HLA alleles predict severe reactions to specific drugs, for example HLA-B\*57:01 with abacavir and HLA-B\*58:01 with allopurinol. The report marks each allele of `bin/hla_alleles.tsv` (the drug alleles, and disease-associated ones such as HLA-B\*27 and HLA-C\*06:02) present or absent, and derives the coeliac DQ types (DQ2.5, DQ2.2, DQ8) and the number of DQB1\*02 copies from the DQA1 and DQB1 rows; [step 8](../08-hla-typing.md#what-the-report-derives) gives the rules and the sources. The genotypes are unphased, so DQ2.5 reads "likely cis" or "possible trans". These lines are information about alleles, not a diagnosis: carrying one does not mean the condition or reaction will happen, and not carrying one does not rule it out. A quality-0 gene with many tied candidates reads "no confident call", and the HLA-HFE row is a sequence allele, not the HFE C282Y/H63D genotype, which comes from the variant calls: the ClinVar hits list each change, and the HFE rule of the traits step reads the two together ([step 42](../42-traits.md)). Typing from short-read WGS is approximate. HLA-A and HLA-B also reach PharmCAT (step 36), so the PharmCAT report and the CPIC recommendations (step 27) give the drug guidance for them; the CPIC file marks them `[outside call]`. When T1K's tied candidates, or a null allele, leave open whether one of the drug alleles is there, step 36 does not pass that gene: it reads `indeterminate` with a reason such as `ambiguous: B*57:01 possible`, and the CPIC file says why. If one of those alleles appears, or is missing and you are about to take the drug, ask for a clinical HLA test; do not rely on this output for transplant matching. With `HLA_QC=true` (opt-in) the HLA section also says, for HLA-A, HLA-B and HLA-C, whether OptiType agrees with T1K at two fields, with both calls; a gene they disagree on, or did not compare, is not passed to PharmCAT when it is HLA-A or HLA-B. The one exception is a T1K null that OptiType calls without its suffix: OptiType never calls a null, so step 36 passes that pair, and the line's note says so. With `KIR=true` (opt-in) the KIR genes are typed too, in `${SAMPLE}/kir_t1k/`, and the HLA section lists the genes T1K typed; no step interprets them.

## More Pharmacogenomics: Cyrius, CPIC, pypgx and the Consensus (Steps 21, 27, 32, 36)

- **CPIC lookup (step 27):** `${SAMPLE}/cpic/${SAMPLE}_cpic_recommendations.txt` turns PharmCAT's calls into the drugs with CPIC guidance. Only genes whose phenotype names a changed function or a positive test get drug entries. Genes PharmCAT called with no such phenotype (`n/a`, the reference VKORC1 genotype `-1639 GG`, `Uncertain Susceptibility`) are listed in their own section without drugs; PharmCAT's HTML report has any guidance for them. Genes PharmCAT could not call are listed separately at the end; their absence from the drug list does not mean normal function. Its section "Calls From Other Tools" says which of HLA-A, HLA-B and CYP2D6 reached PharmCAT from the BAM-based callers, and why the others did not.
- **PGx consensus (step 36):** `${SAMPLE}/pgx_consensus/${SAMPLE}_pgx_consensus.tsv` shows, for HLA-A, HLA-B and CYP2D6, what each caller said and what was passed to PharmCAT. CYP2D6 is passed only when pypgx and Cyrius give the same diplotype and the depth at CYP2D6 passed its check; otherwise it reads `indeterminate`. HLA-A and HLA-B read `indeterminate` too when T1K's call leaves a drug allele open (`ambiguous: B*57:01 possible`).
- **pypgx (step 32):** `${SAMPLE}/pypgx/${SAMPLE}_pypgx_summary.tsv` calls the 29 genes of `assets/pgx_genes.tsv`, 13 of them (the `bam` rows) with copy number from the BAM (GSTT1, which GRCh38 has only on the ALT contig `chr22_KI270879v1_alt`, is reported as not callable on the default no-ALT reference; a BAM aligned to a reference with ALT or HLA contigs stops the step with a `REFUSED` line, because its depth splits between a gene and its ALT copy), so it sees gene deletions and duplications PharmCAT cannot. `${SAMPLE}_pharmcat_comparison.tsv`, written into the same folder by step 27 (CPIC lookup), shows where the two tools agree.
- **Cyrius (step 21, opt-in, non-commercial licence):** `${SAMPLE}/cyrius/${SAMPLE}_cyp2d6.tsv` gives a second CYP2D6 call from the BAM. Without it, CYP2D6 stays `indeterminate`.

The text and HTML reports and the live CPIC card count and name genes marked `incomplete`. PharmCAT did not read one or more of the positions that define an allele for those genes, so step 27 leaves their drug guidance out. The HTML report shows drug contraindication warnings even when every gene has a normal result.

The CYP2D6 consensus shows weak pypgx evidence even when the callers agree. The CYP2D6 depth-check line shows the mean depth over CYP2D6 and its flanks, for all reads and for MAPQ >= 1, beside the mapping-quality shares, and says once when no share was measured. The HLA section flags T1K results withheld from PharmCAT because of tied candidates or an expression suffix.

CYP2D6 is the hard gene: a nearby pseudogene and frequent copy-number changes confuse short reads. The pipeline gives CYP2D6 drug guidance only when two callers agree; even then, take any result that would change a prescription to a pharmacist or a certified pharmacogenomics test first.
