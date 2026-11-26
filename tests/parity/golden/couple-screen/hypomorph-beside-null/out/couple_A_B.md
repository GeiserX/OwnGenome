# Couple carrier screen: A and B

A (female) and B (male), partners in the pedigree. The screen covers 3292 genes of autosomal recessive and X-linked conditions (assets/carrier_genes.tsv: the ACMG and ACOG carrier screening genes and the recessive and X-linked genes of Orphanet).

This is a screen, not a diagnosis, and not medical advice. Before acting on a result, or on its absence, see a genetic counsellor: a clinical carrier test confirms a finding, and the counsellor can say what a negative screen leaves out for your ancestry.

## Shared risk: both partners carry a change in the same gene

- **RBM8A** (Thrombocytopenia-absent radius syndrome; Orphanet): 1 in 4 per pregnancy
  - A: chr1:145927447 C>T het, c.-21G>A (rs139428292), a low-expression allele: TAR syndrome with a null allele such as the 1q21.1 deletion on the other copy, a hypomorph allele
  - B: chr1:145924501 A>G het, Pathogenic, 2 stars, ClinVar 90000
  - Both partners carry a change in RBM8A, linked to Thrombocytopenia-absent radius syndrome (autosomal recessive). A child who inherits the change of each parent can have the condition: 1 in 4 per pregnancy. The change of A is a low-function (hypomorph) allele, not a full loss of function: with the other parent's loss-of-function change a child can have the condition, for some genes in a milder form; two such alleles alone are set aside.

## X-linked findings

No change in an X-linked gene of the list.

## Genes where only one partner carries a change

- **HEXA** (Tay-Sachs disease; ACMG,ACOG): none shared
  - A: chr15:72358680 A>G het, Pathogenic, 2 stars, ClinVar 900011
  - Only A carries a change in HEXA (Tay-Sachs disease, autosomal recessive). A child can be a carrier like A, not affected, unless B carries a change this screen does not see.

## What was assessed

| Partner | ClinVar hits (in listed genes / all) | Deletions | SMN1 | FMR1 and FXN repeats |
|---|---|---|---|---|
| A | 1 / 3 | not assessed (no Delly calls, step 19) | not assessed (no Parascopy, step 35) | not assessed (no ExpansionHunter, step 09) |
| B | 1 / 1 | not assessed (no Delly calls, step 19) | not assessed (no Parascopy, step 35) | not assessed (no ExpansionHunter, step 09) |

Deletions shorter than 50 bp are left to the small-variant caller. Common deletions: not filtered, no benign SV set was given (AnnotSV's annotations, step 05): a deletion many people carry can show as a finding.

## What this screen does not see

- Set aside for A, not a carrier finding on its own: chr15:72350584 G>A het in HEXA: c.739C>T p.Arg247Trp, a pseudodeficiency allele.
- Set aside for A, not a carrier finding on its own: chr1:145927447 C>T het in RBM8A: c.-21G>A (rs139428292), a low-expression allele: TAR syndrome with a null allele such as the 1q21.1 deletion on the other copy, a hypomorph allele.
- Deletions of A: no Delly calls (step 19).
- SMN1 copy number of A: no Parascopy (step 35).
- CYP21A2 deletion and conversions of A: step 35's CYP21A2 table was not passed to the couple screen (--cyp21a2).
- FMR1 and FXN repeats of A: no ExpansionHunter (step 09).
- Deletions of B: no Delly calls (step 19).
- SMN1 copy number of B: no Parascopy (step 35).
- CYP21A2 deletion and conversions of B: step 35's CYP21A2 table was not passed to the couple screen (--cyp21a2).
- FMR1 and FXN repeats of B: no ExpansionHunter (step 09).
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
