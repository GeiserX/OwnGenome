# Step 46: Lp(a), the LPA KIV-2 Repeat Count

## What This Does

Estimates how many kringle IV type 2 (KIV-2) repeats the LPA gene carries on both chromosomes together, from read depth, with a 95% interval. Beside the count it reads the KIV-2 4925G>A splice variant, which lowers the Lp(a) made from the chromosome that carries it. The report's Lp(a) section puts that count beside the two common LPA risk variants that [step 42](42-traits.md) reads, rs10455872 and rs3798220, and flags a genetically high Lp(a) when either points that way. Opt-in.

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
- `bin/copy_number.py` (standard library Python), the copy-number engine, on the LPA_KIV2 row of `bin/cn_loci.tsv` for the count, its interval and the reading; `bin/lpa_kiv2.py` lays out the table. See [Copy number from depth](cn-loci.md)

## Docker Image

- `SAMTOOLS_IMAGE` for the depth and the pileup, `PYTHON_IMAGE` for the count and the 4925G>A reading

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags. The step adds no image and no reference data.

## Input

- `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` with its `.bai` (`ALIGN_DIR` as for the other BAM steps), aligned to GRCh38 with `chr` contig names.
- The aligner must report a read that fits several places equally well, with mapping quality 0. minimap2 (step 2) and BWA-MEM do. novoalign by default leaves such reads unaligned: then the step finds almost no MAPQ 0 reads over KIV-2 and reports the count as not assessed rather than a count that comes out low.

## Command

```bash
./scripts/46-lpa.sh your_name
```

With `run-all.sh`: `TOOLS=...,lpa`. It is opt-in and needs no data, so a run lists step 46 only when `TOOLS` names it, and it is skipped for a sample with no BAM. In the Nextflow pipeline: `lpa` in `--tools` (processes `LPA_REGIONS`, `LPA_DEPTH` and `LPA`, for each sample with a BAM, with the same commands). The two risk variants come from step 42 (`traits`, on by default).

## How the count is made

The method is the one published short-read estimators use: depth over the KIV-2 repeat against depth over single-copy sequence (Zekavat et al., Nat Commun 2018, PMID 29973585; Behera et al., BMC Med Genomics 2024, PMID 39449055). Behera et al. give the KIV-2 region in GRCh38 as chr6:160,611,568-160,646,868, six units, and its reads reach into the last exon of KIV-1 and the first exon of KIV-3, which resemble KIV-2.

The regions are exons of the MANE Select transcript ENST00000316300.10 (Ensembl, checked 2026-10-25), numbered from its first exon. They are the LPA_KIV2 row of `bin/cn_loci.tsv`, and `bin/copy_number.py bed --locus LPA_KIV2` prints them:

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

The 95% interval comes from how much the depth varies between exons: the relative standard error is sqrt(cv_KIV2² / 14 + cv_single² / 15), where cv is each group's coefficient of variation of the per-exon depth, applied to N + 2 on the log scale. It covers sampling noise only, not a bias of an aligner or a library.

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

## How the 4925G>A variant is read

The variant is the last base of the second KIV-2 exon, G>A on LPA's minus strand, so C>T on the plus strand. GRCh38 holds that base at the end of each of its six KIV-2 units (exons 5, 7, ..., 15; `NM_005577.4:c.733G>A` at the first, the one dbSNP names as rs1484919544) and at the near-identical last exon of KIV-1 (exon 3). The 21 bases around it are the same at all seven, so a read from any copy can land on any of them. `samtools mpileup` without a reference reads every read at the seven positions, every mapping quality, duplicates left out, bases of quality 13 or more, and `bin/traits_report.py --marker LPA_KIV2_4925` adds them up into its own table (`bin/lpa_kiv2.py` lays it out). The share of reads with T is a fraction of the copies, not a genotype. With a count, the copies that carry it are that fraction times N + 2 (the KIV-2 repeats plus one KIV-1 per chromosome).

| `status` | When |
|---|---|
| `present` | 5 reads or more with T, and at least 0.8% of the reads (half a copy among 60) |
| `not seen` | 2 reads with T or fewer (one or two can be sequencing errors) |
| `uncertain` | between the two |
| `too little depth` | under 100 reads over the seven positions together |
| `not assessed` | the pileup could not be read (also when `samtools mpileup` failed; the count still runs), or the count found that the aligner left out the reads of MAPQ 0, which are the reads a carrying repeat gives |

What a short-read fraction cannot say: which chromosome or which repeat carries the variant, so whether it sits on the small or the large isoform. A variant in one repeat among many gives few reads, so `not seen` is not proof of absence at low depth.

## Output

| File | Contents |
|---|---|
| `lpa/${SAMPLE}_lpa_kiv2.tsv` | `metric`, `value`: `status`, `kiv2_copies`, `ci95_low`, `ci95_high`, `class`, `kiv2_depth`, `flank_depth`, `mapq0_share`, `message`, `advice` |
| `lpa/${SAMPLE}_lpa_kiv2_4925.tsv` | `metric`, `value`: `marker`, `rsid`, `status`, `alt_reads`, `depth`, `fraction`, `copies` (the copies that carry it, `.` without a count or a reading), `kiv2_copies`, `message` |
| `lpa/${SAMPLE}_lpa_kiv2_depth.tsv` | `samtools bedcov -j` over the regions, every mapping quality |
| `lpa/${SAMPLE}_lpa_kiv2_depth_q1.tsv` | the same, reads of MAPQ 1 or more |
| `lpa/${SAMPLE}_lpa_kiv2_4925_pileup.tsv` | `samtools mpileup` at the seven 4925G>A positions |

Both reports have an Lp(a) section: the count with its interval and class, the two risk variants from step 42 with their risk allele counts (rs10455872 G, rs3798220 C), a flag when there is at least one risk allele or the count is low, and the advice to measure Lp(a) in blood. Without step 42's table the two variants are shown as not read. The reports do not read the 4925G>A table yet; step 46 prints its message beside the count.

## Runtime

Seconds: `bedcov` reads only the 29 exons and `mpileup` the seven 4925G>A positions, through the BAM index.

## Notes

- Research estimate from read depth, not a diagnostic test, and not a measurement of Lp(a). Lp(a) also depends on variants this step does not read.
- The count is the total over both chromosomes. Which chromosome carries the shorter allele matters for the level, and depth cannot tell.
- Exome data can work in principle, since the regions are exons, but the step has only been run on whole-genome reads.
- The 4925G>A reading and the count are separate: a sample with too little depth for the count can still have a 4925G>A reading, then without the copies that carry it.
- On the e2e fixture, a slice of HG002 with no LPA, the step reports the count as not assessed or as too little depth. The e2e case `lpa-1` runs it on HG002's real reads over LPA and checks the arithmetic on a copy with half the KIV-2 reads. The case `genotype-rules-2-lpa-null` plants T in a share of HG002's reads at the seven 4925G>A positions (`present`) and sets them all to C (`not seen`).
