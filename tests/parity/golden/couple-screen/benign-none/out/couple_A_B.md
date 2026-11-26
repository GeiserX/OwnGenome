# Couple carrier screen: A and B

A (female) and B (male), partners in the pedigree. The screen covers 3292 genes of autosomal recessive and X-linked conditions (assets/carrier_genes.tsv: the ACMG and ACOG carrier screening genes and the recessive and X-linked genes of Orphanet).

This is a screen, not a diagnosis, and not medical advice. Before acting on a result, or on its absence, see a genetic counsellor: a clinical carrier test confirms a finding, and the counsellor can say what a negative screen leaves out for your ancestry.

## Shared risk: both partners carry a change in the same gene

No gene of the list where both partners carry a change.

## X-linked findings

No change in an X-linked gene of the list.

## Genes where only one partner carries a change

- **CLN3** (Juvenile CLN3 disease; and 1 other; Orphanet): none shared
  - A: deletion chr16:28482344-28487343 (5.0 kb) het, 1 exon of the gene (Delly)
  - Only A carries a change in CLN3 (Juvenile CLN3 disease; and 1 other, autosomal recessive). A child can be a carrier like A, not affected, unless B carries a change this screen does not see.
- **GALC** (Adult Krabbe disease; and 2 others; Orphanet): none shared
  - A: deletion chr14:87915718-87915767 (50 bp) het, 1 exon of the gene (Delly)
  - Only A carries a change in GALC (Adult Krabbe disease; and 2 others, autosomal recessive). A child can be a carrier like A, not affected, unless B carries a change this screen does not see.
- **GJB2** (Nonsyndromic hearing loss, recessive (DFNB1A); ACMG): none shared
  - A: deletion chr13:20187722-20192721 (5.0 kb) het, 1 exon of the gene (Delly)
  - Only A carries a change in GJB2 (Nonsyndromic hearing loss, recessive (DFNB1A), autosomal recessive). A child can be a carrier like A, not affected, unless B carries a change this screen does not see.
- **HEXA** (Tay-Sachs disease; ACMG,ACOG): none shared
  - A: deletion chr15:72356170-72361169 (5.0 kb) het, 1 exon of the gene (Delly)
  - Only A carries a change in HEXA (Tay-Sachs disease, autosomal recessive). A child can be a carrier like A, not affected, unless B carries a change this screen does not see.
- **MCOLN1** (Mucolipidosis type IV; ACMG,ACOG): none shared
  - A: deletion chr19:7525811-7530810 (5.0 kb) het, 1 exon of the gene (Delly)
  - Only A carries a change in MCOLN1 (Mucolipidosis type IV, autosomal recessive). A child can be a carrier like A, not affected, unless B carries a change this screen does not see.
- **NPC1** (Niemann-Pick disease type C, adult neurologic onset; and 4 others; Orphanet): none shared
  - A: deletion chr18:23544077-23549076 (5.0 kb) het, 1 exon of the gene (Delly)
  - Only A carries a change in NPC1 (Niemann-Pick disease type C, adult neurologic onset; and 4 others, autosomal recessive). A child can be a carrier like A, not affected, unless B carries a change this screen does not see.
- **SMPD1** (Niemann-Pick disease types A and B; ACMG,ACOG): none shared
  - A: deletion chr11:6390190-6395189 (5.0 kb) het, 1 exon of the gene (Delly)
  - Only A carries a change in SMPD1 (Niemann-Pick disease types A and B, autosomal recessive). A child can be a carrier like A, not affected, unless B carries a change this screen does not see.

## What was assessed

| Partner | ClinVar hits (in listed genes / all) | Deletions | SMN1 | FMR1 and FXN repeats |
|---|---|---|---|---|
| A | 0 / 0 | 7 over exons of listed genes; 2 shorter than 50 bp set aside | not assessed (no Parascopy, step 35) | not assessed (no ExpansionHunter, step 09) |
| B | 0 / 0 | not assessed (no Delly calls, step 19) | not assessed (no Parascopy, step 35) | not assessed (no ExpansionHunter, step 09) |

Deletions shorter than 50 bp are left to the small-variant caller. Common deletions: not filtered, no benign SV set was given (AnnotSV's annotations, step 05): a deletion many people carry can show as a finding.

## What this screen does not see

- Common deletions of A: not filtered, no benign SV set was given (AnnotSV's annotations, step 05): a deletion many people carry can show as a finding.
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
