# Step 49: Clonal Haematopoiesis (CHIP)

## What This Does

Looks for mutations that blood stem cells acquired during life and that now sit in a clone of the blood: clonal haematopoiesis of indeterminate potential (CHIP). It calls low-VAF changes in 13 CHIP genes (DNMT3A, TET2, ASXL1, JAK2, TP53, PPM1D, SF3B1, SRSF2, U2AF1, CBL, GNB1, IDH1, IDH2), keeps the ones in each gene's published classes, and says how small a clone it could have found at the sample's depth. Opt-in.

This is not a diagnosis. CHIP becomes common with age, and most people who have it never get a blood cancer.

## Why

A whole genome from blood holds the blood's acquired mutations as well as the inherited ones. Jaiswal et al. found such mutations in 9.5% of people aged 70 to 79 in exome data of 17,182 people, most of them in DNMT3A, TET2 and ASXL1, and they were linked to a higher risk of blood cancer, of coronary heart disease and of death (PMID 25426837; Genovese et al., PMID 25426838, found the same). Bick et al. found CHIP in 4,229 of 97,691 whole genomes of the TOPMed programme (PMID 33057201). The term CHIP and its threshold, a VAF of at least 2% in a gene recurrently mutated in blood cancers, are from Steensma et al. (PMID 25931582).

The germline caller of step 3 is built for 50% and 100% alleles and does not report a change carried by a few reads. Step 29 runs Mutect2 on these genes but leaves its calls unread. This step reads them.

## Tool

- GATK Mutect2, tumour-only, then LearnReadOrientationModel and FilterMutectCalls; GetPileupSummaries and CalculateContamination for the contamination estimate
- bcftools csq, which names each change on the gene's MANE Select transcript, and bcftools mpileup for the depth
- `og chip prepare` and `og chip report` (the Rust binary `og`, `rust/crates/og-chip`), applying the gene table `assets/chip_genes.tsv`; `og chip check-table` checks that table

## Docker Image

- `GATK_IMAGE`, `BCFTOOLS_IMAGE` and `OG_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags.

## Input

- `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` with its `.bai` (`ALIGN_DIR` as for the other BAM steps)
- the reference with its `.fai` and `.dict`
- Ensembl's GRCh38 gene annotation, `reference/ensembl/Homo_sapiens.GRCh38.<release>.chr.gff3.gz` (`ENSEMBL_GFF_RELEASE` in `versions.env`, about 110 MB). Install it once with `./scripts/setup.sh --chip-data "${GENOME_DIR}"`; it is checked against `ENSEMBL_GFF_SHA256`.
- Optional, and the same files as [step 29](29-mutect2-somatic.md): gnomAD's AF-only germline resource `somatic/af-only-gnomad.hg38.vcf.gz` and a panel of normals `somatic/1000g_pon.hg38.vcf.gz`, each with its `.tbi`. Without gnomAD, inherited variants are left out by their VAF only.
- Optional, and recommended together with gnomAD: step 29's common-sites VCF `somatic/small_exac_common_3.hg38.vcf.gz` with its `.tbi` ([step 29](29-mutect2-somatic.md) says how to download it). With it the step estimates how much of the sample is another person's DNA and FilterMutectCalls reads that estimate. Without it the log says there is no contamination estimate, and a contaminant's alleles can pass as a small clone. With it but without gnomAD the log warns that the estimate changes few calls (step 2 below says why).

## Command

```bash
./scripts/setup.sh --chip-data "${GENOME_DIR}"   # once
./scripts/49-chip.sh your_name
```

With `run-all.sh`: `TOOLS=...,chip`. It is opt-in, so a run lists step 49 only when `TOOLS` names it, and it is skipped without a BAM or without the annotation. In the Nextflow pipeline: `chip` in `--tools` with `--chip_gff` (and optionally `--chip_germline_resource`, `--chip_pon` and `--chip_common_sites`; `run-all.sh` passes each one installed under `somatic/`), processes `CHIP_PREPARE`, `CHIP_MUTECT2`, `CHIP_ANNOTATE` and `CHIP_REPORT`. `CHIP_GENES=FILE` runs the script with another table in the same layout.

## How it decides

1. For each gene of `assets/chip_genes.tsv` the step takes the transcript Ensembl tags `MANE_Select`, and calls inside its coding exons, 10 bases wider on each side so the splice donor and acceptor sites are in. That is 162 exons and about 35 kb with Ensembl 116.
2. Mutect2 calls in tumour-only mode, with the read orientation model, and with gnomAD and the panel of normals when they are installed. With the common-sites VCF, GetPileupSummaries counts the reads at every one of its 59,295 sites across the genome, CalculateContamination estimates from them how much of the sample is another person's DNA, and FilterMutectCalls gets that estimate (`--contamination-table`). Another person's DNA adds its own alleles at a low VAF, the VAF of a small clone; with the estimate, FilterMutectCalls marks a call that contamination explains `contamination`, and the call reads `filtered`. The pileups do not stop at the CHIP exons, as step 29's do at its intervals: the exons hold 8 of the 59,295 sites, too few for an estimate. FilterMutectCalls weighs each allele by its population frequency, which Mutect2 takes from gnomAD, so the estimate does its work only with gnomAD installed too. On public data (GIAB's HG002 reads with HG005's reads mixed in, 19,365 of 209,006 reads, over DNMT3A, SF3B1 and IDH1 and 450 windows around common sites on chr2), CalculateContamination read 0.147 (error 0.023) for the mix and 0.006 for HG002 alone. With gnomAD, the mix had 35 PASS calls without the estimate and 4 with it, as many as HG002 alone; without gnomAD the estimate removed none of them. Step 33's VerifyBamID2 gives a second, separate estimate (FREEMIX).
3. Each call that changes the protein is put in the gene's class from the table:

| Gene | Counted | Source (PMID) |
|---|---|---|
| DNMT3A | truncating anywhere; missense at R882 | 25426837, 21067377 |
| TET2 | truncating anywhere | 25426837, 19474426 |
| ASXL1 | truncating in the last two exons, 12 and 13 (11 and 12 in the older transcript numbering of the CHIP studies) | 25426837, 19388938 |
| JAK2 | V617F | 25426837, 15858187 |
| TP53 | truncating anywhere; missense in the DNA-binding core domain, residues 102 to 292 | 25426837, 8023157 |
| PPM1D | truncating in the last exon (6) | 23242139, 29954749 |
| SF3B1 | missense at K700 | 21995386, 23335386 |
| SRSF2 | missense at P95 | 21909114, 23335386 |
| U2AF1 | missense at S34 and Q157 | 22158538, 23335386 |
| CBL | missense in exons 8 and 9 (codons 366 to 477, the linker and RING finger) | 19620960, 19276253 |
| GNB1 | missense at K57 | 25485910 |
| IDH1 | missense at R132 | 19657110, 26016821 |
| IDH2 | missense at R140 and R172 | 20171147, 26016821 |

"Truncating" is a frameshift, a stop gained, or a change of a splice donor or acceptor site. Every codon and exon is numbered on the MANE Select transcript; the residue at each hotspot was checked against Ensembl's protein sequence of that transcript. These are the best-known classes only, narrower than the variant lists of the large CHIP studies: a change outside them is listed as `other_change` and not counted.

4. Each change then gets a status:

| Status | Rule |
|---|---|
| `chip` | passed FilterMutectCalls, in the gene's class, VAF from 2% to under 35%, at least 3 reads with the change |
| `likely_germline` | VAF 35% or more. An inherited variant sits near 50% or 100%, and without a matched normal sample it cannot be told apart from a large clone |
| `filtered` | FilterMutectCalls did not pass it; its filters are in the `note` column |
| `low_support` | fewer than 3 reads with the change |
| `low_vaf` | VAF under 2%, the VAF that defines CHIP (Steensma et al., PMID 25931582). Not counted |
| `other_change` | a protein change outside the gene's class |

The rules are tested in that order, so a call is `low_vaf` only when it passed the filters and has 3 reads. Synonymous, intronic and UTR changes are left out.

A `chip` row is a possible finding, not a diagnosis, and its `note` starts with `possible clonal haematopoiesis (not a diagnosis)`, then the class rule it matched. CHIP itself is a clinical term: it also needs a normal blood count and no blood cancer, which a genome cannot show.

5. A `chip` call whose reads also fit an inherited heterozygous variant keeps its status and gets one more sentence in its `note`: `10 of 30 reads also fit an inherited heterozygous variant`. The test is a binomial count at 50%: the reads fit when at least 1 in 100 inherited heterozygous variants at that depth would show as few reads with the change, or fewer. At 20x that is 5 reads of 20 and up (25%), at 30x 9 of 30 and up (30%), at 40x 13 of 40 and up (33%). From 60x the fit starts at the 35% cut, so no counted call gets the sentence there. These calls stay counted, because a real clone of 25% to 35% looks the same, and the step's log line counts them apart: `2 CHIP calls in 13 genes (possible clonal haematopoiesis, not a diagnosis; 1 of them also fits an inherited variant; ...)`.

## The read-support limit

The step counts the depth of each coding base with bcftools mpileup, at Mutect2's default read filters (mapping quality 20 and up, base quality 10 and up, duplicates left out). A base with no pileup line has depth 0.

Each gene then gets a read-support limit: the smallest VAF at which a clone shows at least 3 reads with the change 8 times in 10, averaged over the gene's coding bases, each at its own depth (a binomial count per base). A base the run did not read counts as a base where no clone can show. With every base at one depth, the limit is about 20% at 20x, 14% at 30x, 10% at 40x, 7% at 60x and 4% at 100x. Uneven depth raises it: half the bases at 20x and half at 40x give about 16%, where the mean of 30x would say 14%.

It is a read count, not the caller's sensitivity. Mutect2's model and FilterMutectCalls can still reject a change that 3 reads carry, so they find fewer clones than the limit says.

The status of each gene:

| Status | Rule |
|---|---|
| `covered` | some VAF reaches 8 in 10; `detection_limit_vaf` holds it |
| `partly_covered` | some coding bases have 3 reads or more, but too few: no VAF reaches 8 in 10 (a fifth of the bases or more under 3 reads is enough). No limit |
| `too_shallow` | no coding base has 3 reads. No limit |
| `not_on_reference` | none of its coding bases is on the reference |

`detection_limit_vaf` is rounded up at the third decimal, so a clone at the written VAF reaches 8 in 10 (0.138 at 30x, where the limit is 0.1373). The stated limit of the run is the median of the covered genes' written limits. The step's log line, the summary (`median_limit`) and both reports give that one figure, with the median of the covered genes' mean depths.

Two more columns help to find a gene with a region the run did not read:

- `callable_fraction`: the share of the gene's coding bases with 10 reads or more.
- `assessable_fraction`: the share of its coding bases where a clone at the stated limit, as written, shows 3 reads 8 times in 10. The limit is an average over bases, so a gene read near the median depth has only part of its bases at full power, often under half. Compare genes by it: a gene far below the others has bases the run did not read. Empty when no gene is covered.

A coding base that two genes of the BED share counts for the limit and both shares of each gene. Its depth counts in the mean depth of only the last of them in the BED, as before. The CHIP genes' exons do not overlap, so this matters only for a BED of your own with overlapping or padded rows.

So at the depth of a usual 30x genome, a clone has to hold about 14% of the alleles to show its reads reliably. CHIP starts at 2%, so a clone between 2% and the limit is likely missed: a result with no call does not rule out CHIP. A gene whose depth is far below the genome's (reads that map with low quality there) gets a higher limit, `partly_covered` when part of it has no reads, or `too_shallow`.

## Output

| File | Contents |
|---|---|
| `chip/${SAMPLE}_chip.tsv` | one row per protein change in a CHIP gene: `gene`, `chrom`, `pos`, `ref`, `alt`, `consequence`, `protein_change` (for example `R882H`, `Q1200fs`), `exon`, `vaf`, `alt_reads`, `depth`, `filter`, `status`, `class` (truncating or missense, for a change in the gene's class), `note` |
| `chip/${SAMPLE}_chip_depth.tsv` | one row per gene: `transcript`, `chrom`, `coding_bases`, `mean_depth`, `detection_limit_vaf` (the read-support limit), `status` (covered, partly_covered, too_shallow, not_on_reference), `callable_fraction`, `assessable_fraction` |
| `chip/${SAMPLE}_chip_filtered.vcf.gz` | Mutect2's calls after FilterMutectCalls |
| `chip/${SAMPLE}_chip.gff3`, `chip/${SAMPLE}_chip.bed` | the transcripts and exons it called on (the bash step) |
| `chip/${SAMPLE}_chip_contamination.table` | CalculateContamination's estimate, with the common-sites VCF (the bash step) |

Both reports show a "Clonal Haematopoiesis (CHIP)" section: the calls, the stated read-support limit, the genes that are not covered, and how many other changes were listed.

## Runtime

Not measured on a whole genome. Mutect2 reads about 35 kb of the BAM, which takes minutes, not hours. With the common-sites VCF, GetPileupSummaries also reads the BAM at the VCF's 59,295 sites.

## Notes

- A `chip` row is a possible finding worth showing a haematologist, together with a blood count. It does not mean a blood cancer, and the risks the studies report are averages over many people.
- Tumour-only calling cannot prove a change is acquired. The VAF rule removes most inherited variants, but a change at 35% or more is never counted, and a real large clone there is missed too. A counted call with the germline-fit sentence may be inherited: a test of a second tissue (a cheek swab or hair) or of a parent tells the two apart.
- A clone under 2% is not CHIP by its definition, so a call under 2% is `low_vaf` and not counted. Such a call needs 3 reads with the change, so more than 150 reads over the site: it shows only where the depth is far above a usual 30x genome.
- **U2AF1 on the default reference.** GRCh38 added to the short arm of chr21 a near copy of the U2AF1 locus that real genomes do not have, and the no-ALT analysis set keeps it. Reads from U2AF1 then fit both copies equally well and get mapping quality 0. Miller et al. 2022 (PMID 35041928) found that GRCh38 pipelines missed U2AF1 mutations for this reason, and that they were found again once the copy was masked. This step counts depth and Mutect2 reads only from mapping quality 20 up, so on the default reference expect U2AF1 to show `too_shallow` or `partly_covered` in `chip_depth.tsv`, or a read-support limit far above the other genes'. Either way, no call at S34 or Q157 does not mean there is none. The fix is the masked reference, at the same coordinates: `setup.sh --masked-reference` builds it, with the false copies set to N as FixItFelix (Behera et al. 2023, PMID 36810122) and GIAB do, and a sample realigned to it gets its U2AF1 reads back ([the masked reference](realignment.md#the-masked-reference)). On chr21 of GIAB's HG002, step 16b reads U2AF1 on the masked reference and not on the default one (CI case `fixitfelix-chr21-1`); this step has not been measured on a whole genome either way.
- Steps 29 and 49 call the same genes in two ways: step 29 calls gene bodies and leaves the calls unread, step 49 calls the coding exons and reads them.
