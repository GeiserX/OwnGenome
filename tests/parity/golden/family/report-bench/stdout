# Family analysis

Sequenced: DAD, GF, GM, KID, MOM. In the pedigree, not sequenced: none.

## Relationships

Relatedness is twice the kinship coefficient: about the share of the genome the two have inherited from common ancestors (1 the same person, 0.5 parent and child or full siblings, 0.25 second degree, 0 unrelated). somalier measures it at its sites; the pedigree gives the expected value.

| Pair | Declared | Expected | Measured | Opposite homozygotes | Sites | Verdict |
|---|---|---|---|---|---|---|

## DAD

Father: GF (sequenced). Mother: GM (sequenced).

Sites of the joint VCF where DAD has a genotype and a sequenced parent does too: 199936. `DAD_inheritance.tsv.gz` lists them.

- the allele from the father is known at 176839 (88.4%)
- the allele from the mother is known at 176890 (88.5%)
- both parents heterozygous like the child, so no origin from the genotypes: 23046 (11.5%)
- inconsistent with the declared parents (a genotype error, a deletion, or a new mutation): 0 (0.0%)
- left out on chrX or chrY outside the pseudoautosomal regions (sex unknown, or a heterozygous call on a single X or Y): 49

De novo candidates (an allele neither parent has): 0, 0 of them with GQ >= 20 and depth >= 10 in the child and in each parent who passes the site on (on a son's X outside the pseudoautosomal regions the mother, on his Y the father), fewer than 2 reads and under 5% of the reads on the new allele in each such parent and, for a heterozygous call, an allele balance of 0.25 to 0.75 (`de_novo.tsv`). Most candidates of a genome are genotype errors; a candidate that matters needs confirming by another method.

## KID

Father: DAD (sequenced). Mother: MOM (sequenced).

Sites of the joint VCF where KID has a genotype and a sequenced parent does too: 199985. `KID_inheritance.tsv.gz` lists them.

- the allele from the father is known at 176935 (88.5%)
- the allele from the mother is known at 176935 (88.5%)
- both parents heterozygous like the child, so no origin from the genotypes: 23050 (11.5%)
- inconsistent with the declared parents (a genotype error, a deletion, or a new mutation): 0 (0.0%)

De novo candidates (an allele neither parent has): 0, 0 of them with GQ >= 20 and depth >= 10 in the child and in each parent who passes the site on (on a son's X outside the pseudoautosomal regions the mother, on his Y the father), fewer than 2 reads and under 5% of the reads on the new allele in each such parent and, for a heterozygous call, an allele balance of 0.25 to 0.75 (`de_novo.tsv`). Most candidates of a genome are genotype errors; a candidate that matters needs confirming by another method.

## Compound heterozygotes

No candidate variants were given (the clinical filter, step 23, did not run for a child), so no pair is reported.

## Chromosome map

For each child, the stretches of each chromosome that came from each parent, from the sites where that parent is heterozygous and the genotypes give the allele passed on. When a parent's own parent is sequenced too, each stretch also says which of the parent's two copies it is (the grandfather's or the grandmother's), and a change between them is a crossover. `<child>_segments.tsv` lists them, `<child>_dnapainter.csv` is the same list in DNA Painter's import format (GRCh38 positions), `<child>_segments.svg` draws it.

- DAD: from the father (GF): 22 stretches over 19.1 Mb; no grandparent on this side is sequenced, so the two copies of the father are not told apart
- DAD: from the mother (GM): 23 stretches over 19.1 Mb; no grandparent on this side is sequenced, so the two copies of the mother are not told apart
- KID: from the father (DAD): 22 stretches over 19.1 Mb, told apart as the copy from GF or GM; 0 crossovers, 25582 discordant sites absorbed
- KID: from the mother (MOM): 22 stretches over 19.1 Mb; no grandparent on this side is sequenced, so the two copies of the mother are not told apart

A stretch starts and ends at its first and last informative site: a crossover lies somewhere between two stretches, not at a base. Short reads give no long-range phase, so without a sequenced grandparent the map says only which parent a stretch came from.
