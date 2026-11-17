# Family analysis

Sequenced: dad, gma, gpa, kid, mum. In the pedigree, not sequenced: none.

## Relationships

Relatedness is twice the kinship coefficient: about the share of the genome the two have inherited from common ancestors (1 the same person, 0.5 parent and child or full siblings, 0.25 second degree, 0 unrelated). somalier measures it at its sites; the pedigree gives the expected value.

| Pair | Declared | Expected | Measured | Opposite homozygotes | Sites | Verdict |
|---|---|---|---|---|---|---|

## dad

Father: gpa (sequenced). Mother: gma (sequenced).

Sites of the joint VCF where dad has a genotype and a sequenced parent does too: 35000. `dad_inheritance.tsv.gz` lists them.

- the allele from the father is known at 26209 (74.9%)
- the allele from the mother is known at 31209 (89.2%)
- both parents heterozygous like the child, so no origin from the genotypes: 3791 (10.8%)
- inconsistent with the declared parents (a genotype error, a deletion, or a new mutation): 0 (0.0%)

De novo candidates (an allele neither parent has): 0, 0 of them with GQ >= 20 and depth >= 10 in the child and in each parent who passes the site on (on a son's X outside the pseudoautosomal regions the mother, on his Y the father), fewer than 2 reads and under 5% of the reads on the new allele in each such parent and, for a heterozygous call, an allele balance of 0.25 to 0.75 (`de_novo.tsv`). Most candidates of a genome are genotype errors; a candidate that matters needs confirming by another method.

## kid

Father: dad (sequenced). Mother: mum (sequenced).

Sites of the joint VCF where kid has a genotype and a sequenced parent does too: 35000. `kid_inheritance.tsv.gz` lists them.

- the allele from the father is known at 26289 (75.1%)
- the allele from the mother is known at 31289 (89.4%)
- both parents heterozygous like the child, so no origin from the genotypes: 3711 (10.6%)
- inconsistent with the declared parents (a genotype error, a deletion, or a new mutation): 0 (0.0%)

De novo candidates (an allele neither parent has): 0, 0 of them with GQ >= 20 and depth >= 10 in the child and in each parent who passes the site on (on a son's X outside the pseudoautosomal regions the mother, on his Y the father), fewer than 2 reads and under 5% of the reads on the new allele in each such parent and, for a heterozygous call, an allele balance of 0.25 to 0.75 (`de_novo.tsv`). Most candidates of a genome are genotype errors; a candidate that matters needs confirming by another method.

## Compound heterozygotes

No candidate variants were given (the clinical filter, step 23, did not run for a child), so no pair is reported.

## Chromosome map

For each child, the stretches of each chromosome that came from each parent, from the sites where that parent is heterozygous and the genotypes give the allele passed on. When a parent's own parent is sequenced too, each stretch also says which of the parent's two copies it is (the grandfather's or the grandmother's), and a change between them is a crossover. `<child>_segments.tsv` lists them, `<child>_dnapainter.csv` is the same list in DNA Painter's import format (GRCh38 positions), `<child>_segments.svg` draws it.

- dad: from the father (gpa): 3 stretches over 60.0 Mb; no grandparent on this side is sequenced, so the two copies of the father are not told apart
- dad: from the mother (gma): 4 stretches over 70.0 Mb; no grandparent on this side is sequenced, so the two copies of the mother are not told apart
- kid: from the father (dad): 4 stretches over 59.9 Mb, told apart as the copy from gma or gpa; 1 crossover, 3 discordant sites absorbed
- kid: from the mother (mum): 4 stretches over 70.0 Mb; no grandparent on this side is sequenced, so the two copies of the mother are not told apart

A stretch starts and ends at its first and last informative site: a crossover lies somewhere between two stretches, not at a base. Short reads give no long-range phase, so without a sequenced grandparent the map says only which parent a stretch came from.
