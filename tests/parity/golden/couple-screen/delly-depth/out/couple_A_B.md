# Couple carrier screen: A and B

A (female) and B (male), partners in the pedigree. The screen covers 3292 genes of autosomal recessive and X-linked conditions (assets/carrier_genes.tsv: the ACMG and ACOG carrier screening genes and the recessive and X-linked genes of Orphanet).

This is a screen, not a diagnosis, and not medical advice. Before acting on a result, or on its absence, see a genetic counsellor: a clinical carrier test confirms a finding, and the counsellor can say what a negative screen leaves out for your ancestry.

## Shared risk: both partners carry a change in the same gene

No gene of the list where both partners carry a change.

## X-linked findings

- **OTC** (Ornithine transcarbamylase deficiency; ACMG): B: daughters carriers
  - B: deletion chrX:38373616-38378615 (5.0 kb) hemizygous, 1 exon of the gene (Delly)
  - B has a change in OTC on his X chromosome (X-linked: Ornithine transcarbamylase deficiency). Every daughter inherits it and is a carrier; a son gets his Y chromosome, not this X.

## Genes where only one partner carries a change

- **CFTR** (Cystic fibrosis; ACMG,ACOG): none shared
  - A: deletion chr7:117501246-117501845 (600 bp) het, 1 exon of the gene (Delly)
  - Only A carries a change in CFTR (Cystic fibrosis, autosomal recessive). A child can be a carrier like A, not affected, unless B carries a change this screen does not see.
- **PAH** (Phenylketonuria; ACMG): none shared
  - A: deletion chr12:77897650-127897649 (50000.0 kb) het, 1 exon of the gene (Delly)
  - Only A carries a change in PAH (Phenylketonuria, autosomal recessive). A child can be a carrier like A, not affected, unless B carries a change this screen does not see.

## What was assessed

| Partner | ClinVar hits (in listed genes / all) | Deletions | SMN1 | FMR1 and FXN repeats |
|---|---|---|---|---|
| A | 0 / 0 | 2 over exons of listed genes; 2 longer than 1000 bp set aside, read depth does not confirm the loss | not assessed (no Parascopy, step 35) | not assessed (no ExpansionHunter, step 09) |
| B | 0 / 0 | 1 over exons of listed genes; 1 longer than 1000 bp set aside, read depth does not confirm the loss | not assessed (no Parascopy, step 35) | not assessed (no ExpansionHunter, step 09) |

Deletions shorter than 50 bp are left to the small-variant caller. Common deletions: not filtered, no benign SV set was given (AnnotSV's annotations, step 05): a deletion many people carry can show as a finding.

## What this screen does not see

- Common deletions of A: not filtered, no benign SV set was given (AnnotSV's annotations, step 05): a deletion many people carry can show as a finding.
- SMN1 copy number of A: no Parascopy (step 35).
- CYP21A2 deletion and conversions of A: step 35's CYP21A2 table was not passed to the couple screen (--cyp21a2).
- FMR1 and FXN repeats of A: no ExpansionHunter (step 09).
- Common deletions of B: not filtered, no benign SV set was given (AnnotSV's annotations, step 05): a deletion many people carry can show as a finding.
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
