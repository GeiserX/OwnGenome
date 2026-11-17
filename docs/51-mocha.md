# Step 51: Mosaic Chromosomal Alterations (MoChA)

## What This Does

Looks for large mosaic chromosomal alterations (mCAs) on chr1 to chr22: a gain, a loss or a copy-neutral loss of heterozygosity (CN-LOH) of millions of bases that a share of the blood cells carries, not every cell. It reports each one with its type, where it is, the chromosome arm it covers and the share of cells that carry it (the cell fraction). Opt-in.

Mosaic loss of Y and of X are [step 33](33-sample-qc.md)'s. This step reads the autosomes only.

## Why

mCAs in blood become more common with age. [Loh et al. 2018](https://pubmed.ncbi.nlm.nih.gov/29995854/) (PMID 29995854) called 8,342 mCAs in 151,202 UK Biobank participants, and people who carried one had a higher risk of a later blood cancer. Most carriers never get one. The same paper introduced the method this step runs: in the cells that carry the event one parental copy of the stretch is lost or gained, so at every heterozygous site of the stretch the allele balance leans the same way along one phased haplotype. Following the phase makes a small lean over many sites visible.

## Tool

- [MoChA](https://github.com/freeseek/mocha) (`bcftools +mocha` and `bcftools +mochatools`, MIT licence), the authors' implementation of that method, version `MOCHA_VERSION` (2025-08-19, commit `MOCHA_COMMIT`)
- Beagle to phase the VCF against the 1000 Genomes high-coverage panel, as the MoChA README does
- `og mocha report` (OwnGenome's Rust binary, [og](rust.md)) for the floor and the calls reported

## Docker Image

- `MOCHA_IMAGE` (bcftools 1.22 with the MoChA plugins, the authors' image, pinned by digest), `BEAGLE_IMAGE` and `PYTHON_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags. MoChA has no tagged release, so `MOCHA_COMMIT` and `MOCHA_VERSION` move with the image, by hand. The e2e case checks that the image prints `MOCHA_VERSION` and that `mocha.c` at `MOCHA_COMMIT` defines it.

## Reference data

`./scripts/setup.sh --mocha-data <genome_dir>` installs `reference/mocha-<MOCHA_PANEL_VERSION>/`:

- `panel/chrN.vcf.gz`: the 1000 Genomes high-coverage phased panel of 3,202 people (release 20220422, EBI), one chromosome at a time. Each file is checked against the md5 of the release's manifest, and the manifest against `MOCHA_PANEL_MANIFEST_SHA256`. Only the biallelic SNVs seen twice or more are kept, with genotypes only; the downloaded file is then deleted.
- `maps/chrN.map`: the PLINK GRCh38 genetic map of step 26b (`ANCESTRY_FINE_MAP_URL`, its mirror and its sha256).
- `cnps.bed`: the panel's common germline deletions and duplications, MoChA's `--cnp` list. Each chromosome file gives its structural variants seen twice or more and over 10 kb (`SVTYPE` CNV, DEL or DUP) before it is deleted, the recipe of the MoChA README (GRCh38, "List of common germline duplications and deletions"). A panel installed before this list existed is rebuilt by a rerun of `--mocha-data`.

The 22 downloads come to about 30 GB (29.8 GB by the server's sizes on 2026-10-27). They are fetched one at a time and each is deleted once cut, so what stays on disk is much smaller. `MOCHA_REGIONS=<bed>` keeps only the regions of a BED file (the e2e case cuts the panel to one 500 kb slice of chr20 this way). `MOCHA_DATA` overrides the folder.

## Input

- `${SAMPLE}/vcf/${SAMPLE}.vcf.gz` with its `.tbi`, with GT, AD, DP and GQ (DeepVariant writes all four)
- the panel above
- `MOCHA_PHASED_VCF=<file inside GENOME_DIR>`: a VCF already phased (GT with `|`), read instead; no Beagle run and no panel needed. The bash step only (`scripts/51-mocha.sh`); the Nextflow processes always phase with Beagle. Without the panel folder's `cnps.bed`, MoChA then runs with no `--cnp` list and says so

## Command

```bash
./scripts/setup.sh --mocha-data /path/to/genome_dir
GENOME_DIR=/path/to/genome_dir ./scripts/51-mocha.sh your_name
```

With `run-all.sh`: `TOOLS=...,mocha`; a run lists step 51 only when `TOOLS` names it. In the Nextflow pipeline: `mocha` in `--tools` with `--mocha_data /path/to/genome_dir/reference/mocha-20220422` (processes `MOCHA_PREP`, `MOCHA_PHASE`, `MOCHA_CALL` and `MOCHA_REPORT`).

## How it works

1. The VCF's biallelic PASS SNVs on chr1 to chr22. A genotype under 10 reads or with GQ under 20 is set missing, as the MoChA README does for WGS. `mochatools` adds the GC content around each site from the reference.
2. Beagle phases the genotypes against the panel, one chromosome at a time (`impute=false`, `seed=1`). A chromosome the panel has no file for is not phased, and the step prints a warning that names it.
3. The phased genotypes go back onto the sites of step 1, on every chromosome, and MoChA runs in its WGS mode (allele depths from AD) with `--genome GRCh38` and `cnps.bed` as `--cnp`. The MHC and KIR regions are left out of its model.
4. `og mocha report` reads MoChA's tables and decides which calls are reported.

A call is reported when all of these hold:

- it is on chr1 to chr22;
- MoChA's own filter of the UK Biobank callset (MoChA README, "Filter callset") keeps it: not a germline copy-number polymorphism, not a likely constitutional event (a BAF deviation of 0.1 or more with sparse heterozygous sites), not a weak call (no BAF deviation estimate and a BAF-and-phase LOD of 10 or less), not a likely germline duplication;
- MoChA estimated its cell fraction, and it is at or above the floor.

**The floor** is the larger of two numbers:

- 0.05. This one is our choice: below a cell fraction of 0.05 the pipeline reports nothing, whatever the depth.
- the read-sampling floor at the sample's depth, `6 / sqrt(D x H)`, with `D` MoChA's median depth and `H` the heterozygous sites per 10 Mb. This is the cell fraction at which a 10 Mb loss or gain of one copy moves the allele balance by three standard errors of read sampling alone. A loss or gain of one copy in a fraction `c` of the cells moves the balance by about `c / 4`; CN-LOH moves it by `c / 2`, so its signal shows at half that floor. A CN-LOH is still reported only at or above the same floor as a loss or gain.

At 30x with about 3,000 heterozygous sites per 10 Mb the read-sampling floor is about 0.02, so the floor is 0.05. At 5x with 300 sites per 10 Mb it is about 0.15. The floor is a lower bound: phasing errors and read bias are not in it, so an event above it is not always found.

The sample is **not assessed** when MoChA's sample filter of the same README fails (a call rate under 0.97, or a phased BAF autocorrelation over 0.03), when MoChA gives no value for one of the two, when MoChA gives no median depth (so no floor can be set), when no heterozygous site on chr1 to chr22 reached MoChA, or when an autosome has under 0.80 of its heterozygous sites phased. Then no call is reported, and the reports say why. The message names each autosome under the minimum, with its phased and heterozygous sites.

**The call rate on WGS.** MoChA reads a WGS VCF by its allele depths (AD). In that mode it drops a site with a missing genotype before it counts the sites, so the call rate in its stats is 1 on every WGS VCF, whatever step 1's depth and GQ filter set missing. The call-rate half of the sample filter therefore never fires on this step's input. It applies only to a stats table that carries a real call rate. On WGS the phased BAF autocorrelation, the phased share below and the floor are what gate the sample.

**The phased share.** MoChA follows a haplotype only across phased heterozygous sites. An unphased site adds nothing to its phased model and adds 0 to its phased BAF autocorrelation, so a chromosome that is mostly unphased passes that filter and shows no event, even one above the floor. Beagle phases only the sites the panel holds, so the share of an autosome's heterozygous sites that end up phased is the share the panel holds. A panel cut with `MOCHA_REGIONS`, or a `MOCHA_PHASED_VCF` that is only partly phased, leaves whole stretches unphased. An autosome the panel has no file for (one a `MOCHA_REGIONS` file does not name) keeps its sites unphased and its row in the sites table, so this check names it with 0 phased sites and the sample is not assessed.

We measured that share on two public genomes, with step 1's filters (biallelic PASS SNVs, at least 10 reads and GQ 20) and the panel's sites (biallelic SNVs seen twice or more, the same position and alleles):

- HG002, of Ashkenazi Jewish descent: the DeepVariant 1.0 call set of a 30x Illumina whole genome (Google's DeepVariant WGS case study), every region.
- HG005, of Han Chinese descent: the Genome in a Bottle v4.2.1 benchmark call set, which holds only its benchmark regions.

| Autosome | HG002 heterozygous sites | HG002 share in the panel | HG005 heterozygous sites | HG005 share in the panel |
|---|---|---|---|---|
| chr1 | 167816 | 0.9705 | 149656 | 0.9614 |
| chr2 | 168167 | 0.9698 | 162468 | 0.9684 |
| chr3 | 157028 | 0.9732 | 138142 | 0.9752 |
| chr4 | 158076 | 0.9729 | 145939 | 0.9717 |
| chr5 | 148594 | 0.9733 | 124577 | 0.9709 |
| chr6 | 145168 | 0.9598 | 131235 | 0.9680 |
| chr7 | 131005 | 0.9692 | 115320 | 0.9676 |
| chr8 | 124187 | 0.9712 | 106123 | 0.9685 |
| chr9 | 95946 | 0.9614 | 81566 | 0.9570 |
| chr10 | 113721 | 0.9706 | 95840 | 0.9641 |
| chr11 | 103275 | 0.9723 | 89849 | 0.9730 |
| chr12 | 107852 | 0.9720 | 92700 | 0.9739 |
| chr13 | 81021 | 0.9747 | 68967 | 0.9679 |
| chr14 | 73041 | 0.9679 | 63691 | 0.9733 |
| chr15 | 64200 | 0.9555 | 55996 | 0.9566 |
| chr16 | 66355 | 0.9683 | 61311 | 0.9561 |
| chr17 | 58942 | 0.9554 | 48337 | 0.9647 |
| chr18 | 59704 | 0.9747 | 56497 | 0.9742 |
| chr19 | 53980 | 0.9517 | 44114 | 0.9695 |
| chr20 | 51614 | 0.9471 | 40867 | 0.9751 |
| chr21 | 29933 | 0.9335 | 24548 | 0.9733 |
| chr22 | 31687 | 0.9410 | 24923 | 0.9651 |

Over the whole genome the share is 0.967 (HG002) and 0.968 (HG005). The lowest share on one autosome is 0.933 (HG002, chr21). The minimum is that lowest share less 0.10, rounded down to a multiple of 0.05: 0.80. The 0.10 margin is for a genome whose ancestry the panel holds fewer people of, which has more variants the panel lacks, and for a caller with more false heterozygous calls. A whole genome phased against the full panel stays well above it; a chromosome that the phasing left mostly unphased falls far below it.

Below the minimum, MoChA's unphased model can still see a large event, but the step reports no call, the same as when nothing is phased. The check counts genotypes written with `|`. A VCF phased only in short read-backed blocks (WhatsHap, for example) passes it, but MoChA needs phase along the whole chromosome: give `MOCHA_PHASED_VCF` phase from a population panel or a trio.

## Output

| File | Contents |
|---|---|
| `mosaic/${SAMPLE}_mocha.tsv` | one row per MoChA call: `chrom`, `start`, `end`, `length_mb`, `arms` (for example "20p, to the p telomere"), `type` (loss, gain, CN-LOH, undetermined), `cell_fraction`, `bdev`, `rel_cov`, `n_hets`, `lod_baf_phase`, `lod_lrr_baf`, `reported` (yes or no) and `reason` |
| `mosaic/${SAMPLE}_mocha_summary.tsv` | `metric`, `value`: status (assessed or not assessed) and why, depth, heterozygous sites (all and phased), sites per 10 Mb, the floor and its two parts, the calls MoChA made and the calls reported, and how the VCF was phased |
| `mosaic/${SAMPLE}_mocha_calls_raw.tsv`, `_mocha_stats.tsv` | MoChA's own `--calls` and `--stats` tables |
| `mosaic/${SAMPLE}_mocha_sites.tsv` | per chromosome: first and last position, sites, heterozygous sites, phased heterozygous sites |
| `mosaic/${SAMPLE}_mocha.log` | MoChA's log |

Both reports show a Mosaic Chromosomal Alterations (MoChA) section with the reported calls and the floor. Without step 51 the report lists autosomal mosaic alterations as not assessed.

## Runtime

Not measured on a whole genome. Beagle's phasing of the 22 chromosomes takes most of the time; it runs with `THREADS` threads and a 14 GB heap (16 GB for the container). MoChA itself reads one phased BCF.

## Notes

- Research result, not a diagnosis. An mCA in blood is common with age and most carriers never get a blood cancer. A clinical test for clonal blood cells is done in a clinical laboratory.
- Blood only. The method assumes the DNA comes from blood cells; saliva holds a mix of cells, which changes what a cell fraction means.
- A constitutional event (one every cell carries, such as a germline deletion) has a cell fraction near 1 and is left out by MoChA's filter; this step is not a test for those.
- A germline copy-number polymorphism is left out by type: MoChA types a call over a region of `cnps.bed` as `CNP_Loss` or `CNP_Gain`. A germline CNV missing from that list (rare, or under 10 kb) gets no such type. Short events (under a few Mb) are hard to tell from these and are reported only when MoChA's filter keeps them.
- Calls are labelled by chromosome arm, not by gene: an mCA spans millions of bases.
- The e2e case (`tests/e2e/mocha-1.sh`) runs the step on HG002's slices with the panel cut to the chr20 slice (MoChA makes no call on chr20; the fixture's other autosomes have no panel file, so they stay unphased and the sample is not assessed, with each of them named; a call MoChA makes on one of those slices, such as a germline CNV like the chr4 UGT2B17 deletion, is untyped because the cut panel's `cnps.bed` holds chr20 only, and is reported "no" because the sample is not assessed), and on a synthetic phased VCF with a 20 Mb CN-LOH on 20p in 15% of the cells (one CN-LOH call, cell fraction 0.10 to 0.20) beside the same VCF with no shift (no call). Nothing in the case is a finding about the fixture's sample.
