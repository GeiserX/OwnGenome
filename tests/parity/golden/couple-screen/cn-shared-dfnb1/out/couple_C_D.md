# Couple carrier screen: C and D

C (female) and D (sex unknown), partners in the pedigree. The screen covers 3292 genes of autosomal recessive and X-linked conditions (assets/carrier_genes.tsv: the ACMG and ACOG carrier screening genes and the recessive and X-linked genes of Orphanet).

This is a screen, not a diagnosis, and not medical advice. Before acting on a result, or on its absence, see a genetic counsellor: a clinical carrier test confirms a finding, and the counsellor can say what a negative screen leaves out for your ancestry.

## Shared risk: both partners carry a change in the same gene

No gene of the list where both partners carry a change in the gene table; read depth (step 35) adds a possible shared risk, below.

Possible shared risk from read depth (step 35's carrier deletion alleles, not in the gene table):

- **GJB2/GJB6 (DFNB1)**: C: del(GJB6-D13S1830) present (read depth, step 35: 1 of 2 copies); D: chr13:20190221 A>G het, Pathogenic, 2 stars, ClinVar 90000. Both partners may carry a change in GJB2/GJB6 (DFNB1): the gene table does not count read-depth alleles, so a clinical test has to confirm them before the chance for a child can be given.

## X-linked findings

No change in an X-linked gene of the list.

## Genes where only one partner carries a change

- **GJB2** (Nonsyndromic hearing loss, recessive (DFNB1A); ACMG): none shared
  - D: chr13:20190221 A>G het, Pathogenic, 2 stars, ClinVar 90000
  - Only D carries a change in GJB2 (Nonsyndromic hearing loss, recessive (DFNB1A), autosomal recessive). A child can be a carrier like D, not affected, unless C carries a change this screen does not see.
- **UGT1A1** (Crigler-Najjar syndrome type 1; and 1 other; Orphanet): none shared
  - C: chr2:233766787 A>G het, Pathogenic, 2 stars, ClinVar 90000
  - Only C carries a change in UGT1A1 (Crigler-Najjar syndrome type 1; and 1 other, autosomal recessive). A child can be a carrier like C, not affected, unless D carries a change this screen does not see.

## Carrier deletion alleles and AZF regions (read depth, step 35)

Each allele is present, absent or not assessable from the read depth over it. A present deletion needs a clinical test. The AZF regions are read on a male partner only. The three beta-globin rows overlap: the Sicilian deletion covers both HBB rows, and the 619 bp and exon 3 rows share 487 bp, so one deletion can read as present on two or three of them; only a clinical test tells which allele it is.

- C, del(GJB6-D13S1830): **present** (read depth, step 35: 1 of 2 copies)
- C, CFTRdele2,3: **not assessable** (read depth, step 35: no row in step 35's interval table)
- C, HBB exon 3 deletion: **not assessable** (read depth, step 35: no row in step 35's interval table)
- C, HBB 619 bp deletion: **not assessable** (read depth, step 35: no row in step 35's interval table)
- C, Sicilian (delta-beta)0 deletion: **not assessable** (read depth, step 35: no row in step 35's interval table)
- C, AZF regions: not applicable (no Y chromosome)

## What was assessed

| Partner | ClinVar hits (in listed genes / all) | Deletions | SMN1 | FMR1 and FXN repeats |
|---|---|---|---|---|
| C | 1 / 1 | not assessed (no Delly calls, step 19) | not assessed (no Parascopy, step 35) | not assessed (no ExpansionHunter, step 09) |
| D | 1 / 1 | not assessed (no Delly calls, step 19) | not assessed (no Parascopy, step 35) | not assessed (no ExpansionHunter, step 09) |

Deletions shorter than 50 bp are left to the small-variant caller. Common deletions: not filtered, no benign SV set was given (AnnotSV's annotations, step 05): a deletion many people carry can show as a finding.

## What this screen does not see

- Deletions of C: no Delly calls (step 19).
- SMN1 copy number of C: no Parascopy (step 35).
- CYP21A2 deletion and conversions of C: step 35's CYP21A2 table was not passed to the couple screen (--cyp21a2).
- FMR1 and FXN repeats of C: no ExpansionHunter (step 09).
- Carrier deletion alleles of D: no interval table (step 35).
- Deletions of D: no Delly calls (step 19).
- SMN1 copy number of D: no Parascopy (step 35).
- CYP21A2 deletion and conversions of D: step 35's CYP21A2 table was not passed to the couple screen (--cyp21a2).
- FMR1 and FXN repeats of D: no ExpansionHunter (step 09).
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
