# Synthetic SV calls for step 22

Three caller VCFs on chr20 with a known answer, read by
`tests/e2e/sv-mito-telomere-steps-1-sv-merge.sh`, `tests/e2e/sv-merge-votes.sh`
and `tests/e2e/sv-merge-consensus-check.sh` (which adds six events of its own
and lists them); the `SURVIVOR_IMAGE` row of
`tests/smoke/commands.tsv` uses the first two events. The records are made up; they describe no
person.

| Event | manta.vcf.in | delly.vcf.in | cnvpytor.vcf.in | Consensus |
|---|---|---|---|---|
| one 2 kb deletion, called 2 bp apart across a 1 kb boundary | 10,100,999-10,102,999 | 10,101,001-10,103,001 | | one record, SUPP=2 |
| two deletions starting in the same 1 kb window, 600 bp and 5 kb long | 10,200,100-10,200,700 | 10,200,300-10,205,300 | | none: the ends are 4.6 kb apart |
| a deletion one caller saw | | | 10,300,000-10,310,000 | none |
| a deletion beside a Delly record genotyped 0/0 | 10,400,001-10,404,001 | 10,400,011-10,404,011, `0/0` | | none: Delly found the reference |
| a deletion beside a Delly call whose sample filter (FT) is LowQual | 10,500,001-10,503,001 | 10,500,011-10,503,011, FT `LowQual` | | none: Delly filtered it |
| a deletion beside a CNVpytor segment with e-value 0.5 | 10,600,001-10,610,001 | | 10,600,001-10,610,001, `./.` | none: CNVpytor's quality rule |

The position binning step 22 used before (chromosome, `int(POS/1000)` and
SVTYPE) gets the first two wrong: it splits the first event over two windows
and counts the second pair as agreement.

In each of the last three events only the Manta call votes. Step 22 kept every
PASS record before `og sv survivor-prep` (`tests/e2e/sv-merge-votes.sh`), so
each read as two callers and the consensus held 4 records instead of 1.

The consensus-review caller inputs exercise repeated calls within a caller,
a difference in mate contig or orientation, explicit INV and TRA, and split
records with missing or discarded IDs. The pinned VCF captures SURVIVOR
1.0.7 output from those three synthetic inputs in Manta, Delly, CNVpytor
order. The Rust tests read that capture. The end-to-end case runs the pinned
merge and checker through step 22.
