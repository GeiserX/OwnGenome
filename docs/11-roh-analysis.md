# Step 11: Runs of Homozygosity (ROH) Analysis

## What This Does
Detects long stretches of homozygous genotypes (autozygous segments) in the genome. These arise when both copies of a chromosomal region are inherited from a common ancestor.

## Why
ROH analysis screens for consanguinity and uniparental disomy (UPD). Long ROH segments increase the risk of autosomal recessive disease by unmasking deleterious variants. ROH patterns also provide population-level ancestry information.

## Tool
- **bcftools roh** (samtools/bcftools)

## Docker Image
- `BCFTOOLS_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Command
```bash
export GENOME_DIR=/path/to/your/data
./scripts/11-roh-analysis.sh your_sample
```

What the script runs. Only PASS records (and records with no filter, as chip VCFs have) go into `bcftools roh`: DeepVariant's `RefCall` and other filtered records are not genotypes to count on. With the allele-frequency panel installed (see [Allele frequencies](#allele-frequencies)), bcftools roh reads each site's population frequency from it; this is the command with the panel's global AF. The Nextflow ROH module reads the same records with the same flags.

```bash
source versions.env   # from the repository root
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data
PANEL=roh_af/${ROH_AF_PANEL_URL##*/}

docker run --rm \
  -v ${GENOME_DIR}/${SAMPLE}/vcf:/data \
  -v ${GENOME_DIR}/roh_af:/genome/roh_af:ro \
  "${BCFTOOLS_IMAGE}" \
  bash -euo pipefail -c "bcftools view -f PASS,. -Ou /data/${SAMPLE}.vcf.gz \
    | bcftools roh --AF-file /genome/${PANEL} -o /data/${SAMPLE}_roh.txt -"
```

Without the panel the script runs `bcftools roh --AF-dflt 0.4` instead. For chip data (no `FORMAT/PL`), the script adds `-G30`.

## Allele frequencies

`bcftools roh` weighs each homozygous call by how common its allele is: a homozygous call at an allele most people carry says little about a run, one at a rare allele says more. So it needs each site's population allele frequency (AF).

- **With the panel** (`scripts/setup.sh` installs it as `roh_af/roh_af_gnomad_v4.1.tsv.gz`, with its `.tbi`): the biallelic PASS SNVs of chr1 to chr22 whose global AF in gnomAD v4.1 genomes (76,215 genomes, released under CC0) is 0.01 to 0.99, with the global AF and the AF of five continental groups. `--AF-file` makes bcftools roh use the panel's sites only, and only where the alleles match; a site of the VCF that is not in the panel is skipped. (bcftools refuses `-T` beside `--AF-file`: the panel already sets the sites.) The step uses the global AF, or the AF of one group when step 26's ancestry table (`ancestry/<sample>_ancestry.tsv`) gives that group a probability of 0.8 or more: AFR reads gnomAD's `afr`, AMR `amr`, EAS `eas`, EUR `nfe` (non-Finnish European) and SAS `sas`. For a group, the step writes a temporary 4-column file of the sites whose group AF is 0.01 to 0.99, the global panel's range. Run step 26 before step 11 to have it; the Nextflow pipeline waits for step 26's table when it will come (`--roh_af_panel` with `prs`, `--pgs_scoring` and `--ancestry_ref`).
- **Without the panel**: `--AF-dflt 0.4`, so every site gets an alternate allele frequency of 0.4. This holds also where a record carries `INFO/AC` and `INFO/AN`: with `--AF-dflt`, bcftools roh does not read them. A homozygous call at a common allele then counts as rare evidence, so the short segments come out too many and too long. The step's note says so.

The first line of `${SAMPLE}_roh.txt` names the frequencies used when they came from the panel, for example `# AF used: nfe (panel roh_af_gnomad_v4.1.tsv.gz); step 26 gives EUR 0.93`. Without the panel there is no such line, and bcftools roh's own header line gives its command with `--AF-dflt 0.4`.

How the panel is built and published: [Reference setup](00-reference-setup.md#roh-allele-frequency-panel).

### Measured on two public genomes

We ran `bcftools roh` 1.24 with the step's flags (PASS records, `-G30`: neither input has `FORMAT/PL`) on two public genomes. Both inputs hold variant records only, as step 11's VCF does:

- **HG002**: the GIAB v4.2.1 benchmark calls, GRCh38, chr1 to chr22.
- **HG00096**: the 1000 Genomes phase 3 genotypes on GRCh38 positions (low coverage, imputed), the sample's non-reference sites on chr1 to chr22.

"Sites" counts the records bcftools roh used. "Runs of 5 Mb or more" is what `${SAMPLE}_roh_summary.txt` lists, with the centromere and heterochromatin bands left out. Both genomes are European, so the group column is gnomAD's `nfe`.

| Genome | Allele frequencies | Sites | Segments | Total | Runs of 5 Mb or more |
|---|---|---:|---:|---:|---|
| HG002 | 0.4 at every site | 4,000,561 | 6,346 | 483.6 Mb | none |
| HG002 | panel, global AF | 1,746,930 | 2,608 | 362.8 Mb | none |
| HG002 | panel, `nfe` AF | 1,684,685 | 2,387 | 328.4 Mb | none |
| HG00096 | 0.4 at every site | 3,997,277 | 6,347 | 449.2 Mb | none |
| HG00096 | panel, global AF | 1,762,917 | 2,570 | 361.5 Mb | chr22:26,579,132-32,429,957, 5.9 Mb |
| HG00096 | panel, `nfe` AF | 1,709,455 | 2,283 | 317.9 Mb | chr22:26,579,132-32,429,957, 5.9 Mb |

What it shows:

- The flat 0.4 inflates the short segments. It gives 2.4 to 2.8 times as many segments as the panel, and 24% to 47% more total length.
- HG002 has no run of 5 Mb or more with any of the three.
- HG00096 has one run of 5 Mb or more with the panel and none with 0.4. With 0.4, bcftools roh split that stretch of chr22 into 13 segments of 1.1 Mb or less. The stretch holds 259 heterozygous calls. 223 of them are indels or SNVs that are not in the panel (rare in gnomAD, or not PASS there), so the panel run does not read them. The other 36 are panel sites, beside 2,387 homozygous panel sites. These counts come from the same panel markers that found the run, so they explain the difference but do not prove the run independently: read it as a candidate that the panel's markers find and the flat 0.4 splits.
- The group's AF gives fewer and shorter segments than the global AF. The runs of 5 Mb or more are the same.

## Output
- `${SAMPLE}/vcf/${SAMPLE}_roh.txt` — bcftools roh output: per-site states (`ST` lines) and segments (`RG` lines)
- `${SAMPLE}/vcf/${SAMPLE}_roh_summary.txt` — the autosomal segments with **5 Mb or more outside the centromere and heterochromatin bands** (chrom, start, end, and `length_bp` and `length_mb` of the part outside the bands). The script and the Nextflow module write the same file with the same threshold and bands.

### Centromere and heterochromatin bands

Few markers map in GRCh38's centromeres (cytoband stain `acen`), its large heterochromatin blocks (`gvar`, such as 1q12, 9q12 and 16q11.2) and the acrocentric stalks (`stalk`), so `bcftools roh` reads long homozygous runs there that are artefacts. The script, the Nextflow module and the reports take the part of each segment that lies in one of these bands out of its length; a segment wholly inside one counts as nothing. The bands are the `acen`, `gvar` and `stalk` rows of UCSC's `cytoBand.hg38.txt` (the file `setup.sh` installs for step 10). Neighbouring rows are merged around each centromere, with two spans on chrY. `MASKED_BANDS` in `rust/crates/og-summary/src/collect/group_c/roh.rs` is their only copy. The step and the Nextflow module call `og roh summary`, which reads that same constant. On chr1 the span is 121.7-143.2 Mb, on chr9 42.2-61.5 Mb, on chr18 15.4-21.5 Mb.

The reports count autosomes only: Segments, Total and Largest keep chr1 to chr22, so they leave out chrX and chrY (a man's single X reads as one long run) and any other contig (chrM, chrEBV, chrUn_*), and the part of each segment in a band. They say how many chrX and chrY segments and how many megabases in bands were left out.

## Interpretation

### Total ROH and parental relationship

When parents are related, a child is autozygous (both copies from the same ancestor) over a fraction F of the genome, the inbreeding coefficient. The expected total of long ROH is F times the length of the autosomes, about 2,900 Mb. The table counts **autosomal segments of 5 Mb or more** (the ones the script prints), without their part in the centromere and heterochromatin bands above. Those bands cover about 217 Mb of chr1 to chr22 (7.5% of 2,875 Mb), and the table below uses the full autosomes, so expect a reported total up to about 7.5% under its row. Shorter segments come mostly from distant shared ancestry and population history, not from the parents' relationship.

| Parents' relationship | F | Expected total of ROH segments of 5 Mb or more |
|---|---|---|
| Not related | ~0 | None or a few segments |
| Second cousins | 1/64 | ~45 Mb |
| First cousins once removed | 1/32 | ~90 Mb |
| First cousins | 1/16 | ~180 Mb |
| Half siblings, uncle and niece, double first cousins | 1/8 | ~360 Mb |
| Parent and child, full siblings | 1/4 | ~720 Mb |

These are averages. Inheritance is random, so one person's total can be well above or below the value for their parents' relationship, and the ranges of neighbouring rows overlap. A total near one row is consistent with it; it does not prove it.

`og summary collect` writes **F_ROH**, the measure of McQuillan et al. 2008 ([PMID 18760389](https://pubmed.ncbi.nlm.nih.gov/18760389/)). Its numerator is the length of autosomal runs of 5 Mb or more outside the bands. Its denominator is the GRCh38 length of the chromosomes this run assessed, with their bands removed. With all 22 autosomes, that is 2,875,001,522 bp less 217,000,000 bp, or 2,658,001,522 bp. With chr20 alone it is 64,444,167 bp less 4,700,000 bp, or 59,744,167 bp. A 21 Mb run then gives F_ROH 0.3515, rather than 0.0079.

The collector streams the per-site `ST` rows to identify assessed autosomes. A chromosome with sites counts even when it has no ROH. No marker-count floor is imposed, so sparse exome and chip input remain supported. This is a chromosome-length denominator, **not a count of callable bases**. One site does not establish coverage across a chromosome. A capture panel or a small interval on every autosome can still look complete at this resolution.

The new `roh.data.scope` object records `autosomes` in numeric order, `unmasked_bp`, `source`, `partial` and a note. Older segment-only files have `source: segments_only`; their observed chromosomes give only a lower bound on scope. An empty file has no denominator and `f_roh: null`. The segment list and length totals remain available.

A relationship label or a one-chromosome UPD flag needs sites on all 22 autosomes and the full 2,658,001,522 bp denominator. Below that floor, or without per-site evidence, `relationship` and `upd_flag` are null. A subset F_ROH is a descriptive fraction of the listed chromosomes and must not be mapped to the parental relationship table. The report and live-view wording for this scope is handled separately.

For complete chromosome scope, the relationship is the row whose F is nearest on a log scale. Boundaries are the geometric means of neighbouring rows; for example 1/16 and 1/8 meet at about 0.088. An F_ROH under 1/128 reads "no sign of related parents". A one-chromosome flag suppresses that relationship label.
To compute the total the table means from the summary file, which takes the band parts out before it keeps the segments of 5 Mb or more, as the reports do:

```bash
awk -F'\t' 'NR > 3 {sum += $4; n++} END {printf "%d segments, %.0f Mb\n", n, sum / 1e6}' \
  "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}_roh_summary.txt"
```

### One chromosome carrying most of the runs: possible uniparental disomy

In uniparental disomy (UPD) both copies of a chromosome come from one parent. When they are two copies of the same parental chromosome (isodisomy), that chromosome reads as one long run of homozygosity, or a few, while the rest of the genome has few long runs. Nakka et al. 2019 (PMID 31607426) found 675 cases of UPD this way, from runs of homozygosity and segments shared with relatives, among 4,400,363 23andMe research participants and 431,094 UK Biobank participants.

The reports flag a chromosome when its segments of 5 Mb or more (band parts left out) add up to at least 20 Mb **and** to more than half of the length of all such autosomal segments. Both numbers are our thresholds, not a published cut-off: 20 Mb on one chromosome is the long-standing sign of isodisomy (see the uniparental disomy note under [Important Notes](#important-notes)), and the share rule keeps out a child of related parents, whose runs spread over many chromosomes. The flag says "may": a run that long can also be chance in a child of related parents.

A flagged chromosome is looked up in [assets/imprinted_regions.tsv](https://github.com/GeiserX/OwnGenome/blob/main/assets/imprinted_regions.tsv). That file is the only imprinting table and is embedded into `og` at build time. Its header credits [Eggermann et al. 2023](https://doi.org/10.1038/s41572-023-00443-4), with per-row PMIDs for the reference facts.

The flag lists the chromosome, run length, share, region names and citations. It names no syndrome. The existing `maternal_upd` and `paternal_upd` JSON keys stay as `"-"` for reader compatibility. The asset retains their reference facts. Those facts cannot be inferred from a region match in one genome.

All six chromosomes use the same chromosome-wide lookup, including chr11 and chr15. Mixed iso/heterodisomy can leave an imprinted locus outside the observed homozygous runs ([Robinson 2000](https://pubmed.ncbi.nlm.nih.gov/10797485/)). A region hint therefore does not claim that a run overlaps the imprinted locus. Region names are broad cytoband labels, not diagnostic locus intervals or a specialist-reviewed assay. No interval was narrowed on the basis of a ROH alone.

The flag's `note` describes long chromosome-focused ROH and suggests investigating possible isodisomy. Its limits are:

- A large deletion of one copy or acquired loss of heterozygosity can give the same pattern. Check the copy-number results and the MoChA results. The collector does not establish copy neutrality or perform that cross-check.
- Shared ancestry, including related parents, can also produce long runs. A flag is a screening hint, not a UPD diagnosis.
- A single genome cannot determine parental origin. The parents' genotypes can establish it. Methylation testing assesses an imprinting defect and alone does not establish UPD.
- Heterodisomy can leave no ROH. Mosaic or short segmental UPD can fall below the thresholds.
- Sites on all autosomes establish chromosome presence only. They do not prove assay completeness within any chromosome.

### Single segments

| Individual ROH Segment | Interpretation |
|---|---|
| <1 Mb | Common, population-level background |
| 1-5 Mb | Distant shared ancestry; many of them is typical of population isolates |
| 5-10 Mb | Counted in the total above; a few can occur even when the parents are not related |
| >10 Mb | Recent shared ancestry; possible uniparental disomy if confined to one chromosome (the reports flag one chromosome from 20 Mb, see above) |

## Important Notes
- The script auto-detects chip data (no FORMAT/PL tag) and adds `-G30` for genotype-only mode
- With the panel, bcftools roh uses gnomAD frequencies at the panel's sites only (see [Allele frequencies](#allele-frequencies)); without it, `--AF-dflt 0.4` gives every site an AF of 0.4, whatever the record's `INFO/AC` and `INFO/AN` say. That inflates the short segments; the long segments are the signal.
- The centromere and heterochromatin bands (above) are taken out of every count; the step and the reports do it, so there is nothing to subtract by hand
- Long runs concentrated on one chromosome may indicate uniparental disomy; the reports flag a chromosome whose segments of 5 Mb or more add up to at least 20 Mb and to more than half of the length of all such autosomal segments (above). The flag requires sites on all 22 autosomes. Check copy number and acquired LOH, and confirm with the parents' genotypes. A methylation test assesses an imprinting defect
- ROH within known disease gene regions warrants checking for homozygous pathogenic variants
