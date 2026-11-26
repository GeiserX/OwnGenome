# Couple carrier screen: A and B

A (female) and B (male), partners in the pedigree. The screen covers 3292 genes of autosomal recessive and X-linked conditions (assets/carrier_genes.tsv: the ACMG and ACOG carrier screening genes and the recessive and X-linked genes of Orphanet).

This is a screen, not a diagnosis, and not medical advice. Before acting on a result, or on its absence, see a genetic counsellor: a clinical carrier test confirms a finding, and the counsellor can say what a negative screen leaves out for your ancestry.

## Shared risk: both partners carry a change in the same gene

No gene of the list where both partners carry a change.

## X-linked findings

- **FMR1** (Fragile X syndrome; ACMG,ACOG): B: daughters carriers
  - B: FMR1 CGG repeat 250, one allele full mutation (more than 200) (ExpansionHunter estimate)
  - B has a change in FMR1 on his X chromosome (X-linked: Fragile X syndrome). Every daughter inherits it and is a carrier; a son gets his Y chromosome, not this X.
- **GLA** (Fabry disease; ACMG): A: sons 1 in 2
  - A: chrX:101400642 A>G het, Pathogenic, 2 stars, ClinVar 90000
  - A carries a change in GLA on one X chromosome (X-linked: Fabry disease). Each son has a 1 in 2 chance of inheriting it and can have the condition; each daughter has a 1 in 2 chance of being a carrier.
- **OTC** (Ornithine transcarbamylase deficiency; ACMG): A: every son, every daughter a carrier
  - A: chrX:38376115 A>G hom, Pathogenic, 2 stars, ClinVar 90000
  - A carries the change on both X chromosomes: every son inherits it and can have Ornithine transcarbamylase deficiency; every daughter is a carrier.

## Genes where only one partner carries a change

None.

## What was assessed

| Partner | ClinVar hits (in listed genes / all) | Deletions | SMN1 | FMR1 and FXN repeats |
|---|---|---|---|---|
| A | 2 / 2 | not assessed (no Delly calls, step 19) | not assessed (no Parascopy, step 35) | not assessed (no ExpansionHunter, step 09) |
| B | 0 / 0 | not assessed (no Delly calls, step 19) | not assessed (no Parascopy, step 35) | 1 expanded |

Deletions shorter than 50 bp are left to the small-variant caller. Common deletions: not filtered, no benign SV set was given (AnnotSV's annotations, step 05): a deletion many people carry can show as a finding.

## What this screen does not see

- Deletions of A: no Delly calls (step 19).
- SMN1 copy number of A: no Parascopy (step 35).
- CYP21A2 deletion and conversions of A: step 35's CYP21A2 table was not passed to the couple screen (--cyp21a2).
- FMR1 and FXN repeats of A: no ExpansionHunter (step 09).
- Deletions of B: no Delly calls (step 19).
- SMN1 copy number of B: no Parascopy (step 35).
- CYP21A2 deletion and conversions of B: step 35's CYP21A2 table was not passed to the couple screen (--cyp21a2).
- Changes ClinVar does not list as pathogenic or likely pathogenic: a new loss-of-function change, or one of uncertain significance, is not counted.
- Genes outside the carrier list (assets/carrier_genes.tsv).
- SMN1 silent carriers: two copies on one chromosome and none on the other read as two copies.
- Changes in genes with a near-identical copy (HBA1 and HBA2, CYP21A2, GBA1, SMN1): short reads place them poorly, and a gene conversion is not seen. A deletion there counts only when Delly calls it; for CYP21A2, step 35's table (--cyp21a2) adds the 30 kb chimeric deletion and conversions at six of its pseudogene variant sites, nothing else.
- Repeat expansions other than FMR1 and FXN, and the exact size of a long FMR1 or FXN allele: short reads do not span it.
- Deletions Delly does not call, deletions longer than 1000 bp that read depth (Delly's RDCN) does not confirm, duplications, inversions and insertions.
- Delly deletions shorter than 50 bp: they count only when the small-variant caller also calls them and ClinVar lists them, so one it misses is missed here.
- With AnnotSV's benign set, a deletion that matches a common one there (1% or more, 50% reciprocal overlap) is not counted: a rare pathogenic deletion at the same place as a common one is missed.
- Low-level mosaicism, and whether two changes of one person in one gene sit on the same copy.

A screen that finds no shared gene lowers the chance of a child with one of these conditions; it does not remove it.
