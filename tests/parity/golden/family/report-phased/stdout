# Family analysis

Sequenced: dad, kid, mum. In the pedigree, not sequenced: grandpa.

## Relationships

Relatedness is twice the kinship coefficient: about the share of the genome the two have inherited from common ancestors (1 the same person, 0.5 parent and child or full siblings, 0.25 second degree, 0 unrelated). somalier measures it at its sites; the pedigree gives the expected value.

| Pair | Declared | Expected | Measured | Opposite homozygotes | Sites | Verdict |
|---|---|---|---|---|---|---|
| dad and kid | parent-child | 0.500 | 0.490 | 0.0% | 5000 | ok: first degree, as declared |
| dad and mum | partners | 0.000 | 0.010 | 12.0% | 5000 | ok: unrelated, as declared |
| kid and mum | parent-child | 0.500 | 0.510 | 0.0% | 5000 | ok: first degree, as declared |

## kid

Father: dad (sequenced). Mother: mum (sequenced).

Sites of the joint VCF where kid has a genotype and a sequenced parent does too: 3. `kid_inheritance.tsv.gz` lists them.

- the allele from the father is known at 3 (100.0%)
- the allele from the mother is known at 3 (100.0%)
- both parents heterozygous like the child, so no origin from the genotypes: 0 (0.0%)
- inconsistent with the declared parents (a genotype error, a deletion, or a new mutation): 0 (0.0%)

De novo candidates (an allele neither parent has): 0, 0 of them with GQ >= 20 and depth >= 10 in the child and in each parent who passes the site on (on a son's X outside the pseudoautosomal regions the mother, on his Y the father), fewer than 2 reads and under 5% of the reads on the new allele in each such parent and, for a heterozygous call, an allele balance of 0.25 to 0.75 (`de_novo.tsv`). Most candidates of a genome are genotype errors; a candidate that matters needs confirming by another method.

WhatsHap phased 3 heterozygous sites of kid by inheritance. At 3 of them the genotypes also give the origin, and the phase agrees with it at 2 (66.7%).

## Compound heterozygotes

No candidate variants were given (the clinical filter, step 23, did not run for a child), so no pair is reported.

## Chromosome map

For each child, the stretches of each chromosome that came from each parent, from the sites where that parent is heterozygous and the genotypes give the allele passed on. When a parent's own parent is sequenced too, each stretch also says which of the parent's two copies it is (the grandfather's or the grandmother's), and a change between them is a crossover. `<child>_segments.tsv` lists them, `<child>_dnapainter.csv` is the same list in DNA Painter's import format (GRCh38 positions), `<child>_segments.svg` draws it.


A stretch starts and ends at its first and last informative site: a crossover lies somewhere between two stretches, not at a base. Short reads give no long-range phase, so without a sequenced grandparent the map says only which parent a stretch came from.
