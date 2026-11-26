# Shared DNA segments

hap-ibd looked for the stretches of DNA each pair of members shares by descent, on the joint call phased by Beagle against the reference panel. Analysed: 1 chromosome (chr3), 98.0 cM of genetic map. Only stretches of 2 cM or more count: shorter ones cannot be told apart from chance sharing.

Fewer than the 22 autosomes were analysed, so the degree below is a rough estimate.

| Pair | Declared | Shared | Of the analysed map | Longest | Stretches | Both copies | Estimate |
|---|---|---|---|---|---|---|---|
| P and Q | not in the pedigree | 88.0 cM | 89.8% | 59.0 cM | 2 | 0.0 cM | first degree (parent-child) |
| P and R | not in the pedigree | 45.0 cM | 45.9% | 45.0 cM | 1 | 0.0 cM | second degree |
| P and S | not in the pedigree | 98.0 cM | 100.0% | 98.0 cM | 1 | 0.0 cM | first degree (parent-child) |
| Q and R | not in the pedigree | 0.0 cM | 0.0% | 0.0 cM | 0 | 0.0 cM | no shared segment |
| Q and S | not in the pedigree | 0.0 cM | 0.0% | 0.0 cM | 0 | 0.0 cM | no shared segment |
| R and S | not in the pedigree | 0.0 cM | 0.0% | 0.0 cM | 0 | 0.0 cM | no shared segment |

Every stretch is in ibd_segments.tsv (chromosome, start, end, cM, and the cM where both copies are shared). Two stretches a phasing error split are joined when the markers between them hold at most one opposite homozygote (one 0/0, the other 1/1), which two people who share a copy never are bar a genotype error. A parent and child share one copy along the whole genome; full siblings share both copies over about a quarter of it.
