# Family analysis

Sequenced: kid, mum. In the pedigree, not sequenced: dad, grandpa.

## Relationships

Relatedness is twice the kinship coefficient: about the share of the genome the two have inherited from common ancestors (1 the same person, 0.5 parent and child or full siblings, 0.25 second degree, 0 unrelated). somalier measures it at its sites; the pedigree gives the expected value.

| Pair | Declared | Expected | Measured | Opposite homozygotes | Sites | Verdict |
|---|---|---|---|---|---|---|
| kid and mum | parent-child | 0.500 | 0.500 | 0.0% | 5000 | ok: first degree, as declared |

## kid

Father: dad (not sequenced). Mother: mum (sequenced).

Sites of the joint VCF where kid has a genotype and a sequenced parent does too: 3. `kid_inheritance.tsv.gz` lists them.

- the allele from the father is known at 2 (66.7%)
- the allele from the mother is known at 2 (66.7%)
- both parents heterozygous like the child, so no origin from the genotypes: 1 (33.3%)
- inconsistent with the declared parents (a genotype error, a deletion, or a new mutation): 0 (0.0%)

The father is not sequenced. At 2 of these sites (66.7%) the genotypes give the allele the father passed on to kid: that is the half of the father's genome kid carries. The other half was not passed on and is not seen here.

No de novo candidates are called: that needs both parents sequenced.

## Compound heterozygotes

No candidate variants were given (the clinical filter, step 23, did not run for a child), so no pair is reported.

## Chromosome map

For each child, the stretches of each chromosome that came from each parent, from the sites where that parent is heterozygous and the genotypes give the allele passed on. When a parent's own parent is sequenced too, each stretch also says which of the parent's two copies it is (the grandfather's or the grandmother's), and a change between them is a crossover. `<child>_segments.tsv` lists them, `<child>_dnapainter.csv` is the same list in DNA Painter's import format (GRCh38 positions), `<child>_segments.svg` draws it.


A stretch starts and ends at its first and last informative site: a crossover lies somewhere between two stretches, not at a base. Short reads give no long-range phase, so without a sequenced grandparent the map says only which parent a stretch came from.
