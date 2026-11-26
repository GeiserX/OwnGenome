# Shared DNA segments

hap-ibd looked for the stretches of DNA each pair of members shares by descent, on the joint call phased by Beagle against the reference panel. Analysed: 2 chromosomes (chr1, chr2), 145.0 cM of genetic map. Only stretches of 2 cM or more count: shorter ones cannot be told apart from chance sharing.

Fewer than the 22 autosomes were analysed, so the degree below is a rough estimate.

| Pair | Declared | Shared | Of the analysed map | Longest | Stretches | Both copies | Estimate |
|---|---|---|---|---|---|---|---|
| A and B | parent-child | 145.0 cM | 100.0% | 98.0 cM | 2 | 0.0 cM | first degree (parent-child) |
| A and C | full siblings | 98.0 cM | 67.6% | 98.0 cM | 1 | 24.5 cM | first degree (full siblings) |
| A and D | not in the pedigree | 0.0 cM | 0.0% | 0.0 cM | 0 | 0.0 cM | no shared segment |
| B and C | parent-child | 0.0 cM | 0.0% | 0.0 cM | 0 | 0.0 cM | no shared segment |
| B and D | not in the pedigree | 3.0 cM | 2.1% | 3.0 cM | 1 | 0.0 cM | distant or unrelated (short segments only) |
| C and D | not in the pedigree | 0.0 cM | 0.0% | 0.0 cM | 0 | 0.0 cM | no shared segment |

Every stretch is in ibd_segments.tsv (chromosome, start, end, cM, and the cM where both copies are shared). Two stretches a phasing error split are joined when the markers between them hold at most one opposite homozygote (one 0/0, the other 1/1), which two people who share a copy never are bar a genotype error. A parent and child share one copy along the whole genome; full siblings share both copies over about a quarter of it.
