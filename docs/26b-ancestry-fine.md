# Step 26b: Fine-Grained Ancestry (admixture fractions and a chromosome painting)

## What This Does

Estimates which reference groups your genome is most similar to, as fractions with a bootstrap range, at two levels: continental groups and regional ones (5 and 17 on the default panel, which leaves Oceania out; see [The groups](#the-groups)). It also paints each chromosome: each segment of each of your two copies gets the continental group it is most similar to. Opt-in.

Without the panel the step prints one line and exits 0:

```
Step 26b skipped: no fine-grained ancestry panel at <genome_dir>/reference/ancestry_fine/hgdp1kgp_v2; install it with scripts/setup.sh --ancestry-fine-panel <genome_dir>
```

## Why

[Step 26](26-ancestry.md) names the most similar of five continental groups of the 1000 Genomes Project. It gives no fractions and cannot tell a mix of two groups from a person between them. The methods behind the fractions and the painting of commercial ancestry tests are open; what differs is the reference panel. This step uses an open panel of 80 populations and open tools, so you can see the result and also where it stops being supported.

## Tools

- **Beagle** 5.5 (Browning et al., AJHG 2021) phases your genotypes against the panel.
- **og ancestry-fine**, OwnGenome's Rust tool, builds the model, projects the fractions, paints the phased chromosomes, and writes the report. `scripts/setup.sh --og` installs it. The bootstrap resamples 5 cM blocks and reference people.

Step 26b uses og for painting. RFMix is licensed for academic research and appears in the benchmark for comparison. [RFMix licence notice](https://github.com/slowkoni/rfmix)

## The Reference Panel

gnomAD's HGDP + 1000 Genomes callset (Koenig et al., Genome Research 2024), phased with SHAPEIT5: 4,091 people of 80 populations sequenced to about 30x on GRCh38. It is part of gnomAD's genome data, which gnomAD releases under the [Creative Commons Zero dedication](https://gnomad.broadinstitute.org/policies): open, no access agreement. `setup.sh --ancestry-fine-panel` reads it from gnomAD's public bucket (`resources/hgdp_1kg/phased_haplotypes_v2`), a region at a time, and keeps:

- the biallelic SNVs with a minor allele frequency of 5% or more, one per 2 kb (`ANCESTRY_FINE_THIN`), for everyone in the panel: phasing reference and painting reference;
- the people's groups (`samples.tsv`), from gnomAD's metadata of the release;
- after LD pruning (r2 above 0.1 within 500 kb removed, among the unrelated reference people), the SNVs ADMIXTURE uses, and each group's allele frequencies at both levels;
- Beagle's GRCh38 genetic map (HapMap, lifted from GRCh37). Beagle's server is its only host; when it does not answer, the install takes the Internet Archive's copy of the same file, which has the same pinned sha256.

The bucket reports a crc32c checksum per file. `versions.env` pins the sha256 of those checksums and sizes (`ANCESTRY_FINE_PANEL_SHA256`), and of the metadata and the map, so a file gnomAD replaces stops the install rather than changing the panel silently.

### The groups

A group name is ours: a set of reference populations we put together and named, not a nationality. Admixed populations (between continental groups) are in no group and are never a reference: ACB, ASW, CLM, MXL, PEL, PUR, Hazara and Uygur. Related people are no reference either: a person is related when gnomAD's metadata flags it (`sample_filters.all_samples_related`).

**A group needs 10 unrelated reference people.** A regional group with fewer leaves the model: its people are no reference at either level and no painting reference. A continental group leaves with it when none of its regional groups is left. The install log and every report name the groups left out. Why 10:

- RFMix fits its CRF weight on people it simulates from `floor(n / 10) + 1` reference people of each group, and leaves those people out of the reference while it does. With one reference person a group then has no reference haplotype at all: RFMix logs `mm::allocate() called at setup_ref_haplotypes():686 (random-forest.cpp) requesting zero or negative size (0)` and can stop with a segmentation fault, on some chromosomes and not others. Two people is RFMix's hard floor.
- Above that floor, a group is meant to stand for a population. Below 10 people (20 haplotypes), the group's supervised ADMIXTURE frequencies and RFMix's reference are a handful of people's genotypes.

On hgdp1kgp_v2 this leaves out Oceania only: gnomAD flags 29 of its 30 people as related (second-degree or ambiguous, in three small populations), so it has one unrelated person. The next smallest group, Southern Africa, has 12. gnomAD's metadata says which kind of relative each person has but not with whom, so an unrelated subset of the 30 cannot be chosen from it. An Oceanian ancestry is then spread over the groups that are modelled, and the reports say Oceania is not one of them.

| Continental group | Regional groups (reference populations) |
|---|---|
| Africa | West Africa (YRI, ESN, GWD, MSL, Yoruba, Mandenka); East Africa (LWK, BantuKenya); Central African rainforest (Biaka, Mbuti); Southern Africa (San, BantuSouthAfrica) |
| West Eurasia (Europe, Middle East, North Africa) | Northwestern Europe (CEU, GBR, Orcadian, French); Northeastern Europe (FIN, Russian); Iberia (IBS, Basque); Italy (TSI, Tuscan, BergamoItalian); Sardinia (Sardinian); Caucasus (Adygei); Levant (Druze, Palestinian); Bedouin (Bedouin); North Africa (Mozabite) |
| Central and South Asia | Central and South Asia (GIH, PJL, BEB, STU, ITU, Balochi, Brahui, Makrani, Pathan, Sindhi, Burusho, Kalash) |
| East Asia | East Asia (CHB, CHS, CDX, JPT, KHV, Han, NorthernHan, Japanese, Dai, Daur, Hezhen, Lahu, Miao, Mongolian, Naxi, Oroqen, She, Tu, Tujia, Xibo, Yi, Cambodian); Siberia (Yakut) |
| Americas (Native American) | Americas (Karitiana, Maya, Pima, Surui, Colombian) |
| Oceania | Oceania (Bougainville, PapuanHighlands, PapuanSepik); left out of the model on hgdp1kgp_v2, see above |

Europe, the Middle East and North Africa form one continental group: at that level their populations are not well apart, and a painting with them split flips between them along every European or Middle Eastern genome. The regional level splits them.

`ANCESTRY_FINE_SAMPLE_MAP=<file>` installs a `samples.tsv` of your own (the same six columns) and fits the frequencies again; the panel's genotypes are kept. Your table is used as given, so give every regional group with reference people 10 or more, and every continental group with reference people at least 2 painting references, each of them a reference person, and every reference person both groups: the install stops on a table that does not, since a group of one makes RFMix fail as above.

The model is fitted again whenever `samples.tsv` changes, so after an update that changes this rule, run `setup.sh --ancestry-fine-panel` once more: it rewrites `samples.tsv` and refits the model, and reuses the panel's genotypes.

### Not in this panel

The Allen Ancient DNA Resource (AADR, release v66.p1, CC0 1.0 on the [Harvard Dataverse](https://doi.org/10.7910/DVN/FFIDCW)) has present-day people of many more regions, Greek and regional Spanish among them. It is not used here: its present-day data are array genotypes on GRCh37 in the packed EIGENSTRAT format, so they need a liftover and a converter no pinned image here has, and they are unphased, so they could not feed the painting. qpAdm models of ancient sources are not part of this step either.

## Docker Images

`ADMIXTURE_IMAGE`, `BEAGLE_IMAGE`, `RFMIX_IMAGE`, with `BCFTOOLS_IMAGE`, `OG_IMAGE` and, for the install, `BCFTOOLS_HTTPS_IMAGE` (the panel is read over https; `BCFTOOLS_IMAGE` has no CA certificates). All are pinned in `versions.env`; [Image versions](versions.md) lists the tags. `ADMIXTURE_IMAGE`, `BEAGLE_IMAGE`, `RFMIX_IMAGE` and `BCFTOOLS_HTTPS_IMAGE` are `# optional`: `setup.sh --ancestry-fine-panel` pulls them, a plain `setup.sh` does not.

## Input

- VCF from DeepVariant (step 3) and its gVCF beside it: the panel's SNVs are genotyped from the gVCF, so the sites where you match the reference count.
- Without a gVCF, a VCF that holds reference calls (0/0), such as the output of `chip-to-vcf.sh`. A variant-only VCF stops the step: from the sites where you differ from the reference alone, the fractions would lean toward the groups where those alleles are common. The rule: at least one in five of the genotyped panel SNVs must be 0/0.
- One sample per file. The step reads the first sample column only, so a gVCF or VCF that holds more than one stops it, in the script and in the pipeline (keep one with `bcftools view -s <sample>`).
- The panel, installed once: `./scripts/setup.sh --ancestry-fine-panel <genome_dir>`.

## Command

```bash
./scripts/setup.sh --ancestry-fine-panel /path/to/genome_dir   # once
./scripts/26b-ancestry-fine.sh your_name
```

`run-all.sh` runs it in the pipeline whenever the panel is installed (`ancestry_fine` in `--tools`, `--ancestry_fine_panel` set); `ANCESTRY_FINE_DIR=none` skips it. In the Nextflow pipeline the step is `ancestry_fine` in `--tools` with `--ancestry_fine_panel <genome_dir>/reference/ancestry_fine/hgdp1kgp_v2`. Settings: `THREADS` (default 4), `ANCESTRY_FINE_BOOTSTRAP` (replicates, default 100; 0 for none), `ANCESTRY_FINE_DIR` (another panel folder inside the genome directory). `ANCESTRY_FINE_BOOTSTRAP` applies to `26b-ancestry-fine.sh` only: the pipeline, and so `run-all.sh`, always runs 100 replicates.

## What the Script Does Internally

1. Genotypes the panel's SNVs from the gVCF (`bcftools convert --gvcf2vcf`, as step 25 does) and keeps a record only with the panel's two alleles. Both gVCF and VCF inputs keep only `FILTER=PASS` or `.`. LowQual, RefCall and other filtered records stay missing. Unfiltered reference blocks still count as reference. This can lower the genotype and reference counts. The VCF rule requiring one reference call in five uses those filtered counts.
2. Projects your genotypes on each level's reference frequencies with `og ancestry-fine project` and runs the bootstrap replicates.
3. Chooses chromosomes with at least 100 genotyped SNVs spanning at least 5 cM.
4. Phases those chromosomes with Beagle against the panel (`impute=false`, `seed=1`).
5. Paints them with `og ancestry-fine paint` against 200 unrelated reference people per continental group at most, taken in turn from each population (`ANCESTRY_FINE_PAINT_PER_GROUP` at install, at least 10: a lower cap stops the install, since it would undo the minimum of 10 reference people per group). The model has two levels. First, for each continental group alone, a copying model over that group's reference haplotypes (Li and Stephens, as FLARE uses it) gives each site's evidence for the group, read both ways along the chromosome. Then an ancestry model over the groups calls each site's group from that evidence, weighted by 0.2, since neighbouring sites copied from one reference haplotype are not independent evidence. Its proportions stay at the continental fractions of part 2, and only the rate of ancestry changes is fitted, over all the chromosomes. It writes `work/paint/<chrom>.msp.tsv` and `work/paint/params.tsv`. A chromosome it cannot read (its phased file, panel or map missing or broken, a panel or phased file of another chromosome, a map with no marker, or no SNV shared with the panel) goes to `work/paint/failed.txt` with the reason in the log, and the others are still painted; the reports list it as not painted.
6. Writes the outputs below.

The Nextflow module (`modules/local/ancestry_fine`) runs five processes. The end-to-end case compares its fractions and painting with the Bash step.

## Output

| File | Contents |
|---|---|
| `ancestry_fine/${SAMPLE}_ancestry_fine.tsv` | one row per group: `level` (`continent`, `region`, or `painting` for the share of the painted genome), `code`, `group`, `fraction`, `se` (bootstrap standard error), `low` and `high` (2.5 and 97.5 percentiles of the bootstrap). Also a `left_out` row for each group the model leaves out, and a `not_painted` row (reason `painting failed`) for each chromosome the painting failed on, both with `.` in the number columns |
| `ancestry_fine/${SAMPLE}_painting.bed` | the painting: chromosome, start, end, haplotype (1 or 2), continental group code |
| `ancestry_fine/${SAMPLE}_painting.png` | chromosomes 1 to 22 to scale, two bars each, coloured by continental group; grey is not painted |
| `ancestry_fine/${SAMPLE}_ancestry_fine.txt` | the same in plain text, with what it means |
| `ancestry_fine/work/` | the genotypes, the fractions and their fit, the phased haplotypes, and the painting (`paint/`) |

The HTML and text reports have a section with the fractions, the painted shares and the picture.

## Runtime

Historical end-to-end timings in CI on the reduced panel: the fixture's slices plus chr22:20-35 Mb, 4 CPUs.

| Part | Time | Notes |
|---|---|---|
| `setup.sh --ancestry-fine-panel` on that panel | 1 min 43 s | 154 MB on disk, 7,789 SNVs of 4,089 people, 1,003 of them for ADMIXTURE |
| The step on the fixture's slices from a gVCF (no chromosome long enough to paint) | 4 s | 2,136 panel SNVs genotyped |
| The step with the painting of chr22:20-35 Mb | 1 min 19 s | 100 bootstrap replicates, 236 painting references (235 since Oceania left the model), measured at the earlier cap of 50 per group |
| The same through the Nextflow pipeline, with the HTML report | 1 min 38 s | earlier engine |

### Painting accuracy and runtime

The benchmark used five simulated genomes per scenario, covering chromosomes 9, 15, and 21. Both engines painted with the same references from a panel capped at 200 per group. The simulator used a fixed set of 19 held-out people who each project at least 99% to their own continental group. Three are Karitiana or Surui references for the Americas, and each of the other groups contributes four. The held-out people were removed from the phasing panel, painting map, and frequency model.

| Scenario | og accuracy | RFMix accuracy | og dosage R² | RFMix dosage R² | og seconds | RFMix seconds | og peak MiB | RFMix peak MiB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| S0 | 0.9981 | 0.9927 | 0.9988 | 0.9954 | 14 | 2663 | 140 | 1280 |
| S1 | 0.9614 | 0.9657 | 0.9704 | 0.9738 | 17 | 2850 | 207 | 1287 |
| S2 | 0.8764 | 0.8983 | 0.8829 | 0.9004 | 22 | 2779 | 216 | 1309 |
| S3 | 0.7227 | 0.6267 | 0.7469 | 0.6618 | 15 | 2123 | 182 | 1272 |
| S4 | 0.7817 | 0.6622 | 0.8049 | 0.6899 | 9 | 2880 | 206 | 1280 |

S0 is unadmixed. S1 and S3 are two-way mixes at different ages. For recent three-way admixture (S2), og was about 2 percentage points behind RFMix in accuracy. For older admixture (S4), og was about 12 points ahead. The comparison ran on a 4-CPU Linux host and measures the painting step per simulated genome. The 1.1 follow-up will work on the gap in recent three-way admixture. See `docs/research/ancestry-engine-benchmark.md` in the repository for the simulation and scoring details.

## Interpreting Results

- **The fractions** say which reference groups your genome is most similar to, and in what proportion. They are not percentages of your ancestors from a country.
- **The bootstrap range** is how much a fraction moves when the SNVs are resampled. A regional fraction with a wide range is not supported by the data; read it as "somewhere in this range", never as one exact percentage.
- **The two levels** are two answers to one question. When the regional fractions inside West Eurasia add up to the continental one, the regional split is how the step divides it, not a second measurement.
- **The painting** shows where along each chromosome each continental group is most similar. A short segment of another group can be a real recent ancestor or noise; long segments are the stronger evidence.

## Limitations

- A label such as Iberia is a cluster we name after the reference populations in it, not a nationality or an ethnicity.
- Continental-level similarity is well supported. Country-level similarity is supported only where the panel has that country and its neighbours are genetically apart. Below that the step says little.
- Splits within Iberia are weak for every method; the panel has IBS (Spanish) and Basque only.
- Greek and Southern Italian ancestry overlap, and the panel has no Greek population at all: a Greek genome is described with Italy, the Levant and the Caucasus.
- Groups the panel samples poorly (Native American, North African) are estimated less well. Oceania is not modelled at all on the default panel (see [The groups](#the-groups)).
- A population the panel lacks, Ashkenazi Jewish for one, gets no group of its own: its genomes are described with the nearest groups.
- The painting uses continental groups only; local ancestry within Europe is not reliable with these tools.

## Notes

- The e2e case (`tests/e2e/ancestry-fine-grained-1.sh`) checks HG002 from its gVCF, HG002's GIAB genotypes on chr22:20-35 Mb, and a planted person whose first half comes from a Yoruba and second half from a Han Chinese person (both left out of the panel): the planted person must come out near half and half with a painting that switches at the junction. A red control in the same case moves every reference person's label to the next group and checks that the HG002 checks then fail. The case also checks that Oceania is left out and named, and that og writes its fitted painting parameters for chr22. A two-sample gVCF must stop the step and the Nextflow module with the reason and no fractions.

## Links

- [gnomAD downloads](https://gnomad.broadinstitute.org/downloads), where the HGDP + 1000 Genomes resource is listed, and its [data policy](https://gnomad.broadinstitute.org/policies)
- [ADMIXTURE](https://dalexander.github.io/admixture/)
- [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html)
- [RFMix](https://github.com/slowkoni/rfmix)
