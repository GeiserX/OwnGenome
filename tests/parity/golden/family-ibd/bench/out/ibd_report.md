# Shared DNA segments

hap-ibd looked for the stretches of DNA each pair of members shares by descent, on the joint call phased by Beagle against the reference panel. Analysed: 22 chromosomes (chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22), 1100.0 cM of genetic map. Only stretches of 2 cM or more count: shorter ones cannot be told apart from chance sharing.

| Pair | Declared | Shared | Of the analysed map | Longest | Stretches | Both copies | Estimate |
|---|---|---|---|---|---|---|---|
| DAD and KID | parent-child | 1031.0 cM | 93.7% | 50.0 cM | 50 | 0.0 cM | first degree (parent-child) |
| DAD and MUM | unrelated | 0.0 cM | 0.0% | 0.0 cM | 0 | 0.0 cM | no shared segment |
| DAD and OTHER | not in the pedigree | 0.0 cM | 0.0% | 0.0 cM | 0 | 0.0 cM | no shared segment |
| DAD and SIB | parent-child | 1004.0 cM | 91.3% | 45.0 cM | 51 | 0.0 cM | first degree (parent-child) |
| KID and MUM | parent-child | 1028.0 cM | 93.5% | 49.0 cM | 53 | 0.0 cM | first degree (parent-child) |
| KID and OTHER | not in the pedigree | 0.0 cM | 0.0% | 0.0 cM | 0 | 0.0 cM | no shared segment |
| KID and SIB | full siblings | 1086.0 cM | 98.7% | 50.0 cM | 30 | 953.0 cM | same person or identical twins |
| MUM and OTHER | not in the pedigree | 0.0 cM | 0.0% | 0.0 cM | 0 | 0.0 cM | no shared segment |
| MUM and SIB | parent-child | 0.0 cM | 0.0% | 0.0 cM | 0 | 0.0 cM | no shared segment |
| OTHER and SIB | not in the pedigree | 0.0 cM | 0.0% | 0.0 cM | 0 | 0.0 cM | no shared segment |

Every stretch is in ibd_segments.tsv (chromosome, start, end, cM, and the cM where both copies are shared). Two stretches a phasing error split are joined when the markers between them hold at most one opposite homozygote (one 0/0, the other 1/1), which two people who share a copy never are bar a genotype error. A parent and child share one copy along the whole genome; full siblings share both copies over about a quarter of it.
