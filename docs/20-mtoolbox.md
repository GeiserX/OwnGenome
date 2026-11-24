# Step 20: Mitochondrial Variant Calling and Heteroplasmy Detection

## What This Does
Calls mitochondrial DNA variants with heteroplasmy fractions — detecting variants present in only a fraction of your mitochondrial copies. Uses GATK Mutect2 in mitochondrial mode. It also looks for single large mtDNA deletions with eKLIPse, counts the mtDNA copies per cell, and names the confirmed disease variants of a short curated list among the calls.

## Why
The main VCF (step 03) calls chrM as if it were a diploid nuclear contig. This step goes deeper, and step 12 then reads its calls for the haplogroup and the contamination check (run this step first):
- **Heteroplasmy detection**: Identifies variants present in only a fraction of mtDNA copies (clinically important for mitochondrial diseases)
- **Dedicated mitochondrial calling**: Mutect2's mitochondrial mode handles the unique properties of mtDNA (high copy number, circular genome, no recombination)
- **Somatic-grade sensitivity**: Detects variants at allele fractions as low as 1-3%

## Tool
- **GATK Mutect2** (Broad Institute) in `--mitochondria-mode`
- **eKLIPse** ([dooguypapua/eKLIPse](https://github.com/dooguypapua/eKLIPse), AGPL-3.0, Goudenège et al. 2019, PMID 30393377) for large deletions, run as its own tool in its own container

> **Note:** This step was originally planned for MToolBox, but no working Docker image exists for MToolBox (see [lessons-learned.md](lessons-learned.md#mtoolbox-no-working-docker-image-exists)). GATK Mutect2 is the standard clinical alternative. The script is still called `scripts/20-mtoolbox.sh`; the name is historical.

## Docker Image
- `GATK_IMAGE`
- `EKLIPSE_IMAGE`: the BioContainers build of Bioconda's `eklipse` 1.8, pinned by digest. eKLIPse has no GitHub release and is no longer maintained; the Bioconda recipe builds upstream commit `3606cb2edac983d2623ddc667b49206c3d01373c` (2022-02-03) with three patches that port it to Python 3. The e2e case `tests/e2e/mito-extended-1-deletions.sh` reads the recipe inside the pinned build and checks that commit.
- `SAMTOOLS_IMAGE` for depth. `og` runs in `OG_IMAGE` (stock Ubuntu) for the copy number and the tables.

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Command
```bash
# Extract chrM reads
samtools view -b sorted.bam chrM > chrM.bam
samtools index chrM.bam

# Call variants in mitochondrial mode
gatk Mutect2 \
  -R reference.fasta \
  -I chrM.bam \
  -L chrM \
  --mitochondria-mode \
  --max-mnp-distance 0 \
  -O chrM_mutect2.vcf.gz

# The control region again on chrM rotated by 8000 bp (below)
og setup chrm-shift --fasta reference.fasta --out shifted/chrM_shifted.fasta \
  --intervals shifted/control.intervals
#   ... realign the chrM reads to it with GATK's BwaSpark, then Mutect2
#   -L shifted/control.intervals into shifted/shifted_mutect2.vcf.gz
og vcf shift-back --linear chrM_mutect2.vcf.gz \
  --shifted shifted/shifted_mutect2.vcf.gz --out shifted/merged_mutect2.vcf.gz

# Filter
gatk FilterMutectCalls \
  -R reference.fasta \
  -V shifted/merged_mutect2.vcf.gz \
  --stats chrM_mutect2.vcf.gz.stats \
  --mitochondria-mode \
  -O chrM_mutect2_filtered.vcf.gz

# Mark possible NuMTs, given the median autosomal depth
gatk NuMTFilterTool \
  -R reference.fasta \
  -V chrM_mutect2_filtered.vcf.gz \
  --autosomal-coverage 30 \
  -O chrM_filtered.vcf.gz
```

The script runs these steps: `./scripts/20-mtoolbox.sh your_sample`. Then it annotates every call with ClinVar and gnomAD (below).

### The shifted-reference pass for the control region

chrM is a circle, but GRCh38 writes it as a line from m.1 to m.16569. A read across the origin, where m.16569 joins m.1, aligns to neither end in full, so the control region around the origin is called from fewer and clipped reads. As Broad's mitochondria pipeline does, the step calls the control region (m.16025-16569 and m.1-575) a second time on chrM rotated by 8000 bp, where the origin lies in the middle:

1. `og setup chrm-shift` writes chrM of the reference rotated by 8000 bp (its base 1 is m.8001), with its `.fai`, and the control region in its positions (`chrM:8025-9144`). It runs in each sample's step, in under a second, so no setup step or extra reference file is needed. It refuses a chrM that is not the 16569 bp of the rCRS.
2. GATK realigns the chrM reads to it (`PrintReads` without duplicates, secondary or supplementary records, `RevertSam` grouped by name, then `BwaSpark` on local Spark: the GATK image has no `bwa` binary), and Mutect2 calls the control region on it. BwaSpark aligns the records two by two as pairs, so a read whose mate is not there is left out first. `RevertSam --SANITIZE` would do that too, but it re-encodes base qualities that are all Q31 or more as if they were the old phred+64 scale, so it is not used.
3. `og vcf shift-back` maps those calls back with `(pos - 1 + 8000) mod 16569 + 1`, so no chain file is needed, and merges them with the linear calls outside the control region. A rotated call whose reference bases cross the origin cannot be written on the linear chrM and is dropped; the step log counts them. The merged VCF carries a `##og_shifted_pass=` header line, which the filters below keep.

FilterMutectCalls and NuMTFilterTool then run on the merged calls, with the linear pass's stats. `mito/shifted/` keeps the rotated calls and the merged ones; the BAMs are removed. Should any part of the pass fail, the step prints a `WARNING` and filters the linear calls alone, as before; the Nextflow `MITO_VARIANTS` does the same.

### ClinVar and gnomAD beside each call

The last stage (`og mito annotate`, in `OG_IMAGE`) writes `${SAMPLE}_chrM_annotated.tsv`: one row per allele of the filtered VCF, with

- ClinVar's record of the same allele from `clinvar/clinvar_chr.vcf.gz`, the ClinVar file `setup.sh` installs: `clinvar_id`, `clinvar_significance` (CLNSIG), `clinvar_review` (CLNREVSTAT) and `clinvar_condition` (CLNDN);
- gnomAD v3.1's chrM sites of the same allele (`gnomad_mito/`, pinned in `versions.env`): `gnomad_af_hom` (the share of gnomAD's samples with the allele at 95% or more), `gnomad_af_het` (at 10% to 95%) and `gnomad_max_hl` (the highest level gnomAD saw).

Alleles are matched by position, REF and ALT after the bases they share are trimmed, so Mutect2's multi-allelic indels meet ClinVar's and gnomAD's single-allele records. An indel in a repeat can be written at any place in it, so an indel that finds no record that way is moved to the left end of its repeat on the rCRS (GRCh38's chrM) and matched there, on both sides. One C more in the m.311-315 tract then meets its record whether either side writes it m.310T>TC or m.315C>CC. An exact spelling wins over a moved one, and an allele whose REF, as written, is not the rCRS at its position is not moved. The table writes each call's position, REF and ALT trimmed as above (Mutect2's m.302 AC>ACC reads m.302 A>AC), never the moved spelling. A file that is not installed leaves its columns `.`, and the table's first lines say which. The reports list each PASS call that is heteroplasmic, or that ClinVar classes as anything but benign (a record with no classification, `not_provided`, does not count), with these columns, and count the ClinVar pathogenic or likely pathogenic ones.

Step 13's VEP annotates the main VCF, whose chrM records come from the nuclear caller: a heteroplasmic call only Mutect2 makes has no VEP record, so this step annotates its own calls. In the Nextflow pipeline `MITO_VARIANTS_ANNOTATE` does the same after `MITO_VARIANTS`, with the files of `--mito_clinvar` and `--gnomad_mito`; without them its table names each file as not given and leaves its columns `.`.

[MITOMAP](https://www.mitomap.org/) has the curated disease status of mtDNA variants. It stays a manual lookup: it answers scripted requests with a 403 and its redistribution terms are not settled, so no step reads it, and none of its tables is shipped.

### Confirmed disease variants

`assets/mito_confirmed.tsv` is a short list of mtDNA variants whose link to a disease is settled, one row each with its ClinVar variation id and the PubMed ids of the primary papers. The reports name every call of the list with its level (AF) and FILTER, whatever the filter says, so a confirmed variant is never hidden. The first rows:

| Variant | Gene | Disease | ClinVar | First described |
|---|---|---|---|---|
| m.3243A>G | MT-TL1 | MELAS; maternally inherited diabetes and deafness | 9589 | Goto 1990, PMID 2102678 |
| m.8344A>G | MT-TK | MERRF | 9579 | Shoffner 1990, PMID 2112427 |
| m.11778G>A | MT-ND4 | LHON | 9708 | Wallace 1988, PMID 3201231 |
| m.3460G>A | MT-ND1 | LHON | 9722 | Howell 1991, PMID 1928099 |
| m.14484T>C | MT-ND6 | LHON | 9688 | Johns 1992, PMID 1417830 |
| m.8470_13446del (4977 bp) | MT-ATP8 to MT-ND5 | Kearns-Sayre syndrome, CPEO, Pearson syndrome | 430680 | Schon 1989, PMID 2711184 |

`scripts/ci/build-mito-confirmed.sh --clinvar clinvar/clinvar.vcf.gz` writes the file: the curated part (the disease, the papers, a note on heteroplasmy) is a table in the script, and ClinVar gives the variation id, the gene and the HGVS name of the same allele. A variant ClinVar no longer holds, or no longer classes pathogenic, stops the build. Only those stable ClinVar fields go in, so any ClinVar release that holds the records writes the same file, and `--check` fails when the committed file is not what ClinVar gives. ClinVar writes the common deletion at the 3' end of its 13 bp repeat (m.8483_13459del); the file gives the left-aligned bases of ClinVar's VCF record (m.8470 to m.13446). To add a variant, add its row to the script's table and rerun it.

### Large deletions (eKLIPse)

A single large deletion removes several kilobases of the mtDNA in part of the molecules. It causes Kearns-Sayre syndrome, chronic progressive external ophthalmoplegia (CPEO) and Pearson syndrome; the commonest one, 4977 bp long, lies between two 13 bp direct repeats (Schon 1989, PMID 2711184). Short reads that cross the deletion's junction align on one side and are soft-clipped on the other. eKLIPse BLASTs the clipped parts against the mtDNA, keeps the breakpoints found from both sides, and gives each deletion's level: the share of the molecules that carry it, in %.

Stage 8 runs eKLIPse on the chrM reads of stage 1, against the rCRS GenBank file its image ships (GRCh38's chrM is the rCRS), with eKLIPse's default settings, and writes `${SAMPLE}_chrM_deletions.tsv` (`og mito deletions`): per deletion `bkp5` and `bkp3` (the last base kept before it and the first one kept after it, as eKLIPse writes them), the deleted length, the level and its forward and reverse strand estimates, the BLAST hits and depths at each breakpoint, and the repeat around it. eKLIPse's own files and its circular plot stay in `mito/eklipse/out/`.

The reports name a deletion of `assets/mito_confirmed.tsv` when both breakpoints lie within 50 bp of the listed ones, since a deletion between two repeats can be written anywhere along them. They list every named deletion and any other at a level of 1% or more, and count the rest.

**Blood often lacks the deletion.** In Kearns-Sayre syndrome and CPEO a single deletion is often found in muscle and not in blood; in Pearson syndrome it is in blood (GeneReviews, PMID 20301382; Moraes 1989, PMID 2541333). A blood genome with no deletion does not rule one out; a muscle biopsy is the usual sample to test.

### mtDNA copy number

Stage 7 reads the read counts of the BAM index (`samtools idxstats`): chrM reads per base over the reads per base of chr1-22, times 2 for the two nuclear copies, is the number of mtDNA copies per cell (`og mito copies`, `${SAMPLE}_mtdna_copies.tsv`; it reads chrM under any of the names `chrM`, `MT`, `chrMT` and `M`). It is the depth ratio that whole-genome studies use (Ding 2015, PMID 26172475; Longchamps 2020, PMID 32004343). The counts include duplicates and reads of any mapping quality, on chrM and the autosomes alike. Under 20 reads per kb of autosome (about 3x at 150 bp reads) the BAM is a slice, a panel or a shallow genome, and the number is not assessed. An exome can pass that floor and is not told apart: the number is meant for whole-genome BAMs. The copy number depends on the tissue, and in blood on the mix of cells, and on the library, so it compares samples of one kind only.

The copy-number catalog (`bin/cn_loci.tsv`, [docs/cn-loci.md](cn-loci.md)) has no row for it yet: its engine reads targets and controls on one contig, and this ratio needs controls on the autosomes.

### NuMTs

NuMTs are copies of mitochondrial DNA in the nuclear genome. Their reads can map to chrM and look like low-level heteroplasmy. GATK's NuMTFilterTool marks an allele `possible_numt` when its depth is no more than such copies could give at the sample's autosomal depth. The script reads the median autosomal depth from step 16b's mosdepth output (`mosdepth/<sample>.mosdepth.summary.txt` and `.global.dist.txt`), so run step 16b first; `AUTOSOMAL_COVERAGE=30` sets it by hand. Without either, the filter runs at depth 0 and marks nothing. The step then prints a `WARNING`, and the annotation table gets a `# numt` header line that reads `NuMT filter used a default depth` (`og mito annotate --numt-default`). The Nextflow `MITO_VARIANTS` runs the same filter with the same depth, read from `MOSDEPTH`'s output when `mosdepth` is in `--tools` (a default tool). Without it the depth is 0: the run prints a warning on the console when it starts, and `MITO_VARIANTS_ANNOTATE` writes the same `# numt` line. The reports do not show that line yet.

## Output
- `${SAMPLE}_chrM_mutect2.vcf.gz` — Raw mitochondrial variant calls on the linear chrM, with its `.stats`
- `shifted/shifted_mutect2.vcf.gz` — raw calls of the control region on the rotated chrM, in its positions; `shifted/merged_mutect2.vcf.gz` — those calls mapped back, merged with the linear ones outside the control region (the bash step only; Nextflow keeps them in its work directory)
- `${SAMPLE}_chrM_mutect2_filtered.vcf.gz` — after FilterMutectCalls
- `${SAMPLE}_chrM_filtered.vcf.gz` — after FilterMutectCalls and NuMTFilterTool, with PASS or the reasons a call failed (`possible_numt` among them); the file the reports and step 12 (haplogroup, haplocheck) read
- Each variant includes an `AF` (allele fraction) field indicating heteroplasmy level
- `${SAMPLE}_chrM_annotated.tsv` — each allele of the filtered VCF with its FILTER, AF, ClinVar record and gnomAD v3.1 frequencies (above)
- `${SAMPLE}_mtdna_copies.tsv` — the mtDNA copy number per cell, or why it was not assessed (from `${SAMPLE}_idxstats.txt`)
- `${SAMPLE}_chrM_deletions.tsv` — the large deletions eKLIPse found; `eklipse/out/` holds eKLIPse's own output. Should eKLIPse exit 0 without writing its `eKLIPse_deletions.csv`, the step warns and writes no table, and the report says the deletions were not looked for

## Interpreting Heteroplasmy
| AF Level | Meaning |
|---|---|
| >0.95 | Homoplasmic — effectively fixed variant |
| 0.10-0.95 | Heteroplasmic — mixed population, clinically significant threshold varies |
| 0.03-0.10 | Low-level heteroplasmy — may be age-related somatic |
| <0.03 | Near detection limit |

Step 07 reads these calls for the aminoglycoside risk gene MT-RNR1: a PASS allele PharmCAT names (m.1555A>G, m.1494C>T and others) at AF 0.05 or more is passed to PharmCAT, heteroplasmic or not, because CPIC sets no heteroplasmy level at which aminoglycosides are safe ([Pharmacogenomics](07-pharmacogenomics.md#input-the-preprocessor-or-pharmcat-refuses)).

Count only PASS calls: a call marked `possible_numt` or with another filter is not evidence of heteroplasmy. Without the shifted-reference pass, calls near the two ends of the control region (m.1-500 and m.16000-16569) are less reliable: the reference is circular, so reads that span its end align poorly at both edges.

PASS alone does not make a call reliable, so `summary.json` flags two kinds of heteroplasmic call and keeps both in the count:

- `origin_edge`: the bases the call changes touch m.1-500 or m.16000-16569, the edges above. Only when the calls have no `##og_shifted_pass=` line, so only when the shifted-reference pass did not run: after it, no call sits next to a reference end, and step 20's section of `summary.json` names the pass in `shifted_pass`.
- `low_af`: the AF is under 0.10, the low-level range above. NuMT reads that the filter missed look like calls in this range.

Each entry of `heteroplasmic_sites` carries both flags, and `heteroplasmic_origin_edge` and `heteroplasmic_low_af` count the flagged ones beside the `heteroplasmic` total. The rows of the annotation list carry the same two flags, so a ClinVar class on a flagged call reads as a class on a less reliable call. On those rows `origin_edge` is set for any listed call, homoplasmic or not, and `low_af` only for a heteroplasmic one. A flag never removes a call: the ClinVar pathogenic count and the confirmed variants still name it. A PASS m.3243A>G at AF 0.07 carries `low_af` and stays in both.

Length variants in the D-loop poly-C tracts (m.303-315 and m.16184-16193), one C more or less, are found in most people and often read as a mix of lengths. The reports list them as "length heteroplasmy in the D-loop poly-C tracts" and leave them out of the heteroplasmic count. A call counts as in a tract when the bases it changes (an indel: the bases after its anchor base) overlap one.

## Runtime
1 to 2 minutes per sample through the pipeline, measured in three 30-40x WGS runs ([runtime](hardware-requirements.md#runtime-per-step)), before eKLIPse and the shifted-reference pass were added; eKLIPse's run on a whole genome's chrM reads has not been timed yet. The shifted-reference pass took about 3 minutes on 4 CPUs for the 2.5 million chrM reads of the public HG002 35x genome (about 18,000x on chrM), where the linear Mutect2 took about 6.

In the Nextflow pipeline `MITO_DELETIONS` runs stages 7 and 8 on the same BAM, beside `MITO_VARIANTS`, whenever `mito_variants` is in `--tools`. As in the script, an eKLIPse error or a run with no CSV leaves no deletion table and keeps the copy number; the report then says the deletions were not looked for.

## Notes
- `--mitochondria-mode` disables several filters inappropriate for mtDNA: no germline filtering, adjusted LOD thresholds, handles high copy number.
- `--max-mnp-distance 0` prevents merging nearby variants into multi-nucleotide polymorphisms.
- The GATK Docker image is large (~2.2 GB) but well-maintained and versioned.
- Disease annotation: ClinVar and gnomAD come with the step (`${SAMPLE}_chrM_annotated.tsv`); [MITOMAP](https://www.mitomap.org/) is a manual lookup in a browser, as it refuses scripted requests.
- Some mitochondrial diseases require heteroplasmy above a tissue-specific threshold (e.g., m.3243A>G MELAS requires >60% in blood).

## What is not done, and why

- **mtDNA-Server 2 (or mutserve) in place of this step.** Both were last released in December 2024 (mtDNA-Server 2.1.16, mutserve 2.0.3). The cheaper fix for the main source of false heteroplasmy, NuMT reads, is in place (NuMTFilterTool above), and haplocheck, the contamination check of the mtDNA-Server tools, runs in step 12 on these calls.

