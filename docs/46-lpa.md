# Step 46: Lp(a), the LPA KIV-2 Repeat Count

## What This Does

Estimates how many kringle IV type 2 (KIV-2) repeats the LPA gene carries on both chromosomes together, from read depth, with a 95% interval for read sampling plus the error measured against assembled genomes. Beside the count it reads the KIV-2 4925G>A splice variant, which lowers the Lp(a) made from the chromosome that carries it. The report's Lp(a) section puts that count beside the two common LPA risk variants that [step 42](42-traits.md) reads, rs10455872 and rs3798220, and flags a genetically high Lp(a) when either points that way. Opt-in.

It is a research estimate, not a diagnostic test. The test that counts is a blood Lp(a) measurement.

## Why

Lp(a) is a lipoprotein particle whose level is mostly inherited, and a high level is a causal risk factor for atherosclerotic cardiovascular disease and is linked to aortic valve stenosis. The 2022 European Atherosclerosis Society consensus statement advises measuring it at least once in every adult (Kronenberg et al., Eur Heart J 2022, PMID 36036785).

Two parts of LPA predict the level:

- **The KIV-2 repeat count.** The KIV-2 repeat is a 5.5 kb unit, repeated a variable number of times. Fewer repeats give a smaller apo(a) and, on average, a higher Lp(a). In a Danish population study the total over both chromosomes ranged from 6 to 99, and mean Lp(a) was 56 mg/dL in the quarter with the fewest repeats against 15 mg/dL in the quarter with the most; the count explained 21% to 27% of the variation in Lp(a) level (Kamstrup et al., JAMA 2009, PMID 19509380).
- **Two common variants.** rs10455872 and rs3798220 were each associated with a higher Lp(a) level, a smaller apo(a) and coronary disease, with odds ratios of 1.70 and 1.92; one risk allele at either went with an odds ratio of 1.51, two or more with 2.57 (Clarke et al., NEJM 2009, PMID 20032323). The study was of people of European descent.

- **The KIV-2 4925G>A splice variant.** It sits at the donor splice site of the second KIV-2 exon and is usually in one or two of a person's KIV-2 repeats. In a European population 22% carried it, most of them with isoforms of 19 to 25 kringle IV repeats; carriers had a lower Lp(a) than non-carriers across the whole isoform range, and the variant explained 20.6% of the variation in Lp(a) among people with small isoforms (Coassin et al., Eur Heart J 2017, PMID 28444229).

A VCF cannot show the repeat count or the 4925G>A variant: the six KIV-2 units in GRCh38 are near-identical, so short reads from any copy pile up on them. That pile-up is what this step measures.

## Tool

- samtools `bedcov -j` for the read depth over the exons, `mpileup` for the reads at the 4925G>A sites
- `og cn` (OwnGenome's Rust binary, [docs/rust.md](rust.md)), the copy-number engine, on the LPA_KIV2 row of `bin/cn_loci.tsv` for the count, its interval, the reading and the table, and in the same call the 4925G>A reading. See [Copy number from depth](cn-loci.md)

## Docker Image

- `SAMTOOLS_IMAGE` for the depth and the pileup, `OG_IMAGE` for the count and the 4925G>A reading

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags. The step adds no image and no reference data.

## Input

- `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` with its `.bai` (`ALIGN_DIR` as for the other BAM steps), aligned to GRCh38 with `chr` contig names.
- A failure of the 4925G>A reader never costs the count: step 46 and the `LPA` process then write a 4925G>A table that says `not assessed, the check failed`, print a warning, and go on.
- The aligner must report a read that fits several places equally well, with mapping quality 0. minimap2 (step 2) and BWA-MEM do. novoalign by default leaves such reads unaligned: then the step finds almost no MAPQ 0 reads over KIV-2 and reports the count as not assessed rather than a count that comes out low.

## Command

```bash
./scripts/46-lpa.sh your_name
```

With `run-all.sh`: `TOOLS=...,lpa`. It is opt-in and needs no data, so a run lists step 46 only when `TOOLS` names it, and it is skipped for a sample with no BAM. In the Nextflow pipeline: `lpa` in `--tools` (processes `LPA_REGIONS`, `LPA_DEPTH` and `LPA`, for each sample with a BAM, with the same commands). The two risk variants come from step 42 (`traits`, on by default).

## How the count is made

The method is the one published short-read estimators use: depth over the KIV-2 repeat against depth over single-copy sequence (Zekavat et al., Nat Commun 2018, PMID 29973585; Behera et al., BMC Med Genomics 2024, PMID 39449055). Behera et al. give the KIV-2 region in GRCh38 as chr6:160,611,568-160,646,868, six units, and its reads reach into the last exon of KIV-1 and the first exon of KIV-3, which resemble KIV-2.

The regions are exons of the MANE Select transcript ENST00000316300.10 (Ensembl, checked 2026-10-25), numbered from its first exon. They are the LPA_KIV2 row of `bin/cn_loci.tsv`, and `og cn bed --locus LPA_KIV2` prints them:

| Group | Exons | What they are |
|---|---|---|
| KIV-2 block | 3 to 16 | the six KIV-2 units (exons 4 to 15), plus exon 3 (KIV-1) and exon 16 (KIV-3): seven 160 bp and seven 182 bp exons, 2,394 bp |
| single-copy | 19 to 33 | kringle exons from the second exon of KIV-4 to KV, each under 90% identical to any KIV-2 exon |

Exons 17 and 18 (about 95% identical to KIV-2) are in neither group, nor are exons 1 and 2 or the protease exons, which resemble plasminogen.

samtools `bedcov -j` sums the depth over each exon, every mapping quality counted and duplicates, secondary and failed reads left out; a second pass counts reads of MAPQ 1 or more only. With D the mean depth over the single-copy exons and M the mean depth over the KIV-2 block, the block holds N + 2 copies of a 342 bp exon pair (N KIV-2 units plus the KIV-1/KIV-3 pair on each chromosome) spread over seven reference pairs:

```
N = 14 x M / D - 2      (KIV-2 repeats on both chromosomes)
```

GRCh38's own six units on each chromosome give M = D and N = 12.

The 95% interval has two parts, added in quadrature and applied to N + 2 on the log scale:

- read sampling, from how much the depth varies between exons: sqrt(cv_KIV2² / 14 + cv_single² / 15), where cv is each group's coefficient of variation of the per-exon depth;
- a systematic error of 0.035 (3.5% of N + 2), the `cv:0.035` in the `interval` column of the LPA_KIV2 row, measured against the KIV-2 repeats of 96 assembled genomes ([below](#how-the-interval-was-measured)). Read sampling alone treats the exons as independent and leaves out any bias of an aligner or a library, which moves every exon of the block together.

The table's `message` says so: `about N KIV-2 repeats in total over both copies (95% interval LOW to HIGH, read sampling plus the error measured against assembled genomes)`, then the class reading. The `ci95_low` and `ci95_high` columns keep their names.

The class compares that interval with the range of most people: about 70% of the 1000 Genomes samples had 30 to 45 KIV-2 copies with DRAGEN's caller (Behera et al. 2024; mean 37.5 in Europeans, 35.7 in Africans, 44.6 in East Asians).

| Class | When |
|---|---|
| `low` | the whole interval is under 30: a low count goes with small apo(a) isoforms and a higher Lp(a) |
| `high` | the whole interval is over 45: on average a lower Lp(a) |
| `typical` | otherwise |

The count is not made when:

| Status | When |
|---|---|
| `too little depth` | the single-copy exons have under 10x |
| `not assessed` | reads of MAPQ 0 make up under 20% of the depth over the KIV-2 block (the aligner left out the reads that fit several copies), the KIV-2 depth is too low for any repeat count, or the BAM has no chr6 |

## How the interval was measured

We checked the count against genomes whose KIV-2 repeats can be counted directly, in assemblies of each chromosome, and widened the interval by the error we found.

**The assembled count.** For each haplotype assembly, we found the contig that covers LPA from its alignment to GRCh38, took a window of 450 kb on each side, and aligned two GRCh38 anchors to it (chr6:160,555,001-160,565,000 and chr6:160,660,001-160,670,000). Between the anchors we counted the copies of the two KIV-2 exons (182 bp and 160 bp) at 97% identity or more. GRCh38 itself gives 7 of each (its six KIV-2 units plus the KIV-1 and KIV-3 exon that match them), so the KIV-2 count of a haplotype is that number minus 1. A haplotype was kept only when the two exons gave the same number and the span between the anchors agreed with it within half a unit (5,545 bp per unit).

**The estimate.** For the same genomes, we read the 1000 Genomes 30x reads (PCR-free, BWA-MEM, GRCh38) over chr6:160,530,000-160,670,000, ran the two `samtools bedcov -j` passes of this step over the 29 exons, and ran `og cn call`. HG002 used its 35x reads of the Genome in a Bottle DeepVariant release.

**The genomes.** HG002 from its Q100 assembly (24 + 14 = 38), and the 95 samples of the HPRC release 2 that are also in 1000 Genomes with 30x reads and gave a count for both haplotypes: 96 genomes in all. Four had no haplotype alignment over LPA (HG00272, HG00280, HG02922, NA20850), one had exon counts that disagreed (NA19682), and one had a truncated CRAM index (HG01786); they are left out. The assembled totals range from 16 to 57 (median 37), the flank depth from 28x to 57x.

**What we found.**

- The estimate is above the assembled count in all 96 genomes: by a factor of 1.114 on N + 2 on average (4.6 repeats), with a spread of 0.038 on the log scale.
- With read sampling alone, the interval held the assembled count in 88 of 96 genomes (91.7%), under the 95% it claims.
- With a systematic term of 0.035, the smallest in steps of 0.005 that reaches 95%, it holds it in 93 of 96 (96.9%). The interval gets 7% wider on average (0.381 to 0.406 of N + 2), and no class changes. The three it misses are HG00097, HG00146 and NA18982, each with the estimate high.
- We widen the interval and leave the estimate as it is. The bias is the same sign in every genome, so a later change can correct the scale; that needs reads from this pipeline's own aligner (minimap2), not BWA-MEM.

`og dev evidence lpa-truth` prints the assembled counts below, and `og dev evidence lpa-bench CALLS.tsv OUT.tsv MIN` compares a list of calls with them and fails under MIN. The `lpa-benchmark` job of the public-evidence workflow repeats the whole comparison on fresh slices of the 1000 Genomes reads, gates at 95%, and checks that the gate fails when each genome is given another genome's estimate.

**Sources.** HPRC release 2 assemblies and their GRCh38 alignments: https://github.com/human-pangenomics/hprc_intermediate_assembly (`data_tables/assemblies_release2_v1.0.index.csv` and `annotation/alignments_to_ref/alignments_to_ref_grch38_winnowmap_hprc_r2_v1.1.index.csv`). HG002 Q100: `s3://human-pangenomics/working/HPRC_PLUS/HG002/assemblies/Q100/pansn/` (hg002v1.1). 1000 Genomes 30x reads: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index, which lists each CRAM's md5. Reference: GRCh38_full_analysis_set_plus_decoy_hla, chr6 MD5 5691468a67c7e7a7b5f2a3a683792c29.

| Genome | Assembled (hap 1 + hap 2) | Estimate | 95% interval | Holds |
|---|---|---|---|---|
| HG00097 | 28 (10 + 18) | 34.4 | 28.9 to 40.9 | no |
| HG00099 | 20 (10 + 10) | 22.2 | 18.5 to 26.5 | yes |
| HG00126 | 25 (14 + 11) | 28.2 | 24.0 to 33.0 | yes |
| HG00128 | 46 (23 + 23) | 52.1 | 43.8 to 61.9 | yes |
| HG00133 | 57 (27 + 30) | 62.1 | 47.2 to 81.5 | yes |
| HG00140 | 36 (20 + 16) | 39.8 | 32.0 to 49.4 | yes |
| HG00146 | 41 (21 + 20) | 51.0 | 43.1 to 60.4 | no |
| HG002 | 38 (24 + 14) | 43.1 | 37.4 to 49.7 | yes |
| HG00232 | 40 (18 + 22) | 45.0 | 37.9 to 53.4 | yes |
| HG00235 | 30 (18 + 12) | 30.1 | 23.9 to 37.7 | yes |
| HG00253 | 41 (19 + 22) | 47.2 | 40.2 to 55.3 | yes |
| HG00290 | 29 (11 + 18) | 32.5 | 27.2 to 38.8 | yes |
| HG00320 | 39 (18 + 21) | 43.1 | 34.1 to 54.4 | yes |
| HG00321 | 28 (15 + 13) | 32.5 | 25.6 to 41.2 | yes |
| HG00323 | 27 (8 + 19) | 29.1 | 22.9 to 36.8 | yes |
| HG00329 | 37 (15 + 22) | 41.6 | 32.0 to 54.0 | yes |
| HG00344 | 36 (22 + 14) | 41.0 | 32.4 to 51.7 | yes |
| HG00350 | 42 (23 + 19) | 48.2 | 41.5 to 56.0 | yes |
| HG01167 | 53 (29 + 24) | 60.6 | 48.6 to 75.5 | yes |
| HG01530 | 31 (19 + 12) | 35.6 | 28.9 to 43.7 | yes |
| HG01784 | 38 (20 + 18) | 41.5 | 33.3 to 51.7 | yes |
| HG02040 | 24 (15 + 9) | 25.3 | 20.5 to 31.1 | yes |
| HG02155 | 49 (22 + 27) | 55.4 | 45.2 to 67.9 | yes |
| HG02165 | 31 (11 + 20) | 34.9 | 27.9 to 43.4 | yes |
| HG02178 | 42 (19 + 23) | 47.0 | 39.6 to 55.7 | yes |
| HG02391 | 51 (28 + 23) | 58.4 | 48.8 to 70.0 | yes |
| HG02392 | 46 (21 + 25) | 51.1 | 43.1 to 60.4 | yes |
| HG02583 | 30 (17 + 13) | 32.5 | 26.3 to 40.0 | yes |
| HG02976 | 32 (14 + 18) | 38.6 | 29.8 to 49.9 | yes |
| HG03130 | 37 (16 + 21) | 41.2 | 32.0 to 52.9 | yes |
| HG03139 | 33 (18 + 15) | 37.5 | 29.2 to 47.9 | yes |
| HG03195 | 37 (17 + 20) | 41.4 | 32.1 to 53.3 | yes |
| HG03209 | 36 (19 + 17) | 42.8 | 33.7 to 54.3 | yes |
| HG03225 | 36 (15 + 21) | 41.4 | 32.2 to 53.1 | yes |
| HG03270 | 34 (19 + 15) | 36.4 | 28.4 to 46.6 | yes |
| HG03369 | 42 (20 + 22) | 45.0 | 36.0 to 56.1 | yes |
| HG03470 | 38 (20 + 18) | 43.4 | 35.6 to 52.9 | yes |
| HG03521 | 38 (20 + 18) | 42.6 | 34.0 to 53.2 | yes |
| HG03583 | 28 (9 + 19) | 31.9 | 25.4 to 40.0 | yes |
| HG03742 | 39 (20 + 19) | 44.3 | 36.5 to 53.6 | yes |
| HG03784 | 34 (22 + 12) | 39.8 | 30.7 to 51.4 | yes |
| HG03874 | 43 (21 + 22) | 48.3 | 39.7 to 58.6 | yes |
| NA18505 | 45 (27 + 18) | 49.3 | 38.4 to 63.1 | yes |
| NA18508 | 52 (20 + 32) | 56.1 | 47.2 to 66.6 | yes |
| NA18522 | 36 (26 + 10) | 38.6 | 32.3 to 46.1 | yes |
| NA18565 | 40 (24 + 16) | 44.3 | 37.5 to 52.4 | yes |
| NA18570 | 35 (22 + 13) | 36.8 | 28.5 to 47.2 | yes |
| NA18608 | 31 (19 + 12) | 33.6 | 28.7 to 39.4 | yes |
| NA18620 | 30 (15 + 15) | 34.3 | 28.6 to 41.0 | yes |
| NA18747 | 44 (21 + 23) | 48.8 | 41.0 to 58.1 | yes |
| NA18879 | 37 (22 + 15) | 42.3 | 34.3 to 52.2 | yes |
| NA18940 | 39 (20 + 19) | 45.8 | 39.0 to 53.7 | yes |
| NA18943 | 47 (20 + 27) | 51.2 | 43.3 to 60.6 | yes |
| NA18944 | 48 (22 + 26) | 50.9 | 39.4 to 65.5 | yes |
| NA18945 | 42 (18 + 24) | 48.5 | 41.7 to 56.4 | yes |
| NA18948 | 28 (9 + 19) | 31.8 | 24.7 to 40.7 | yes |
| NA18952 | 55 (25 + 30) | 63.6 | 53.5 to 75.4 | yes |
| NA18959 | 49 (25 + 24) | 56.7 | 43.5 to 73.7 | yes |
| NA18960 | 41 (19 + 22) | 47.0 | 36.4 to 60.5 | yes |
| NA18967 | 51 (26 + 25) | 55.4 | 46.6 to 65.9 | yes |
| NA18970 | 37 (14 + 23) | 41.2 | 34.3 to 49.3 | yes |
| NA18971 | 48 (25 + 23) | 55.5 | 47.2 to 65.3 | yes |
| NA18974 | 31 (17 + 14) | 32.2 | 25.8 to 40.0 | yes |
| NA18976 | 35 (14 + 21) | 38.9 | 30.1 to 50.3 | yes |
| NA18982 | 35 (22 + 13) | 42.6 | 36.2 to 50.1 | no |
| NA18983 | 47 (23 + 24) | 50.8 | 39.0 to 65.9 | yes |
| NA19036 | 36 (14 + 22) | 39.5 | 32.4 to 48.1 | yes |
| NA19043 | 16 (1 + 15) | 18.6 | 16.0 to 21.7 | yes |
| NA19087 | 45 (20 + 25) | 45.4 | 38.0 to 54.1 | yes |
| NA19159 | 35 (22 + 13) | 37.5 | 29.7 to 47.3 | yes |
| NA19185 | 36 (22 + 14) | 40.2 | 31.9 to 50.4 | yes |
| NA19338 | 39 (20 + 19) | 47.1 | 36.6 to 60.5 | yes |
| NA19391 | 24 (9 + 15) | 26.8 | 21.2 to 33.7 | yes |
| NA19443 | 27 (20 + 7) | 31.3 | 24.5 to 39.9 | yes |
| NA19468 | 32 (22 + 10) | 37.4 | 29.0 to 48.2 | yes |
| NA19700 | 44 (17 + 27) | 48.4 | 38.4 to 61.0 | yes |
| NA19776 | 32 (7 + 25) | 36.3 | 29.4 to 44.7 | yes |
| NA19835 | 43 (19 + 24) | 47.7 | 36.6 to 62.0 | yes |
| NA19909 | 33 (20 + 13) | 37.2 | 28.7 to 48.1 | yes |
| NA20282 | 41 (21 + 20) | 43.4 | 33.3 to 56.4 | yes |
| NA20346 | 44 (22 + 22) | 47.9 | 38.2 to 60.0 | yes |
| NA20503 | 41 (22 + 19) | 47.7 | 38.8 to 58.6 | yes |
| NA20752 | 33 (14 + 19) | 37.6 | 29.8 to 47.3 | yes |
| NA20762 | 53 (32 + 21) | 60.4 | 47.1 to 77.2 | yes |
| NA20799 | 31 (20 + 11) | 36.3 | 29.3 to 44.9 | yes |
| NA20805 | 34 (14 + 20) | 38.4 | 30.3 to 48.5 | yes |
| NA20806 | 50 (23 + 27) | 59.0 | 45.5 to 76.3 | yes |
| NA20809 | 35 (16 + 19) | 38.5 | 30.8 to 48.0 | yes |
| NA20827 | 41 (22 + 19) | 43.1 | 33.8 to 54.9 | yes |
| NA20870 | 39 (28 + 11) | 40.2 | 34.0 to 47.5 | yes |
| NA20905 | 39 (16 + 23) | 44.7 | 35.9 to 55.6 | yes |
| NA21093 | 31 (17 + 14) | 33.1 | 26.4 to 41.2 | yes |
| NA21102 | 43 (24 + 19) | 49.6 | 38.8 to 63.2 | yes |
| NA21106 | 27 (11 + 16) | 30.2 | 25.6 to 35.6 | yes |
| NA21110 | 43 (19 + 24) | 45.7 | 39.0 to 53.5 | yes |
| NA21144 | 37 (17 + 20) | 40.0 | 33.1 to 48.2 | yes |

## How the 4925G>A variant is read

The variant is the last base of the second KIV-2 exon, G>A on LPA's minus strand, so C>T on the plus strand. GRCh38 holds that base at the end of each of its six KIV-2 units (exons 5, 7, ..., 15; `NM_005577.4:c.733G>A` at the first, the one dbSNP names as rs1484919544) and at the near-identical last exon of KIV-1 (exon 3). The 21 bases around it are the same at all seven, so a read from any copy can land on any of them. `samtools mpileup` without a reference reads every read at the seven positions, every mapping quality, duplicates left out, bases of quality 13 or more, and `og cn call --locus LPA_KIV2 ... --marker LPA_KIV2_4925` adds them up into its own table, in the same call that makes the count. The share of reads with T is a fraction of the copies, not a genotype. With a count, the copies that carry it are that fraction times N + 2 (the KIV-2 repeats plus one KIV-1 per chromosome).

| `status` | When |
|---|---|
| `present` | 5 reads or more with T, and at least 0.8% of the reads (half a copy among 60) |
| `not seen` | 2 reads with T or fewer (one or two can be sequencing errors) |
| `uncertain` | between the two |
| `too little depth` | under 100 reads over the seven positions together |
| `not assessed` | the pileup could not be read (also when `samtools mpileup` failed; the count still runs), or the count found that the aligner left out the reads of MAPQ 0, which are the reads a carrying repeat gives. The marker reads the count's own MAPQ 0 share, not the rounded value in the count table, so the two never disagree at the 20% line |

What a short-read fraction cannot say: which chromosome or which repeat carries the variant, so whether it sits on the small or the large isoform. A variant in one repeat among many gives few reads, so `not seen` is not proof of absence at low depth.

## Output

| File | Contents |
|---|---|
| `lpa/${SAMPLE}_lpa_kiv2.tsv` | `metric`, `value`: `status`, `kiv2_copies`, `ci95_low`, `ci95_high`, `class`, `kiv2_depth`, `flank_depth`, `mapq0_share`, `message`, `advice` |
| `lpa/${SAMPLE}_lpa_kiv2_4925.tsv` | `metric`, `value`: `marker`, `rsid`, `status`, `alt_reads`, `depth`, `fraction`, `copies` (the copies that carry it, `.` without a count or a reading), `kiv2_copies`, `message` |
| `lpa/${SAMPLE}_lpa_kiv2_depth.tsv` | `samtools bedcov -j` over the regions, every mapping quality |
| `lpa/${SAMPLE}_lpa_kiv2_depth_q1.tsv` | the same, reads of MAPQ 1 or more |
| `lpa/${SAMPLE}_lpa_kiv2_4925_pileup.tsv` | `samtools mpileup` at the seven 4925G>A positions |

Both reports have an Lp(a) section: the count with its interval and class, the two risk variants from step 42 with their risk allele counts (rs10455872 G, rs3798220 C), a flag when there is at least one risk allele or the count is low, and the advice to measure Lp(a) in blood. Without step 42's table the two variants are shown as not read. Right under the count, both reports and the live person page print the 4925G>A line: `present` or `uncertain` with the reads that carry it and the copies, `not seen` with the reads, `too little depth`, or `not assessed`; then the reader's message and the caveat that a null allele lowers the Lp(a) made from the copies it sits on and that the Lp(a) blood test is the measure. Without the 4925G>A table (a run from before this reader) or with a table that has no known status, the line says `not assessed` with the reason. The live card also lists the state as a fact. The Nextflow `HTML_REPORT` gets the table from the `LPA` process the same way.

## Runtime

Seconds: `bedcov` reads only the 29 exons and `mpileup` the seven 4925G>A positions, through the BAM index.

## Notes

- Research estimate from read depth, not a diagnostic test, and not a measurement of Lp(a). Lp(a) also depends on variants this step does not read.
- The count is the total over both chromosomes. Which chromosome carries the shorter allele matters for the level, and depth cannot tell.
- Exome data can work in principle, since the regions are exons, but the step has only been run on whole-genome reads.
- The 4925G>A reading and the count are separate: a sample with too little depth for the count can still have a 4925G>A reading, then without the copies that carry it.
- On the e2e fixture, a slice of HG002 with no LPA, the step reports the count as not assessed or as too little depth. The e2e case `lpa-1` runs it on HG002's real reads over LPA, checks that the interval holds the 38 KIV-2 repeats of HG002's Q100 assembly, and checks the arithmetic on a copy with half the KIV-2 reads. The case `genotype-rules-2-lpa-null` plants T in a share of HG002's reads at the seven 4925G>A positions (`present`) and sets them all to C (`not seen`).
