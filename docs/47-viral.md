# Step 47: Viral DNA in the Unmapped Reads

## What This Does

Takes the reads that did not map to the human genome, realigns them to 24 viral genomes, and reports per virus how many unique, mappable reads it got, how much of its genome they cover and how deep, against the human depth of the same sample. Opt-in.

Its main finding is **inherited chromosomally integrated HHV-6** (iciHHV-6). In under 1% of people in the USA and the UK, and about 1% in other studies, a whole HHV-6A or HHV-6B genome sits in a telomere of one chromosome and is inherited like any other part of the genome ([Pellett et al. 2012](https://doi.org/10.1002/rmv.715); [Liu et al. 2020](https://doi.org/10.1371/journal.pgen.1008915)). Every cell then carries one viral copy, so the reads cover the whole viral genome at about half the depth of a human chromosome, which has two copies. Knowing it matters because a blood PCR test for HHV-6 then reads high, which can be taken for an active infection ([Pellett et al. 2012](https://doi.org/10.1002/rmv.715)).

The other viruses are reported as found, as below the threshold, or as not found. **This is not a test for an infection.** Most adults carry herpesviruses that stay in their cells for life, and laboratory reagents carry viral DNA of their own: in 8,240 blood genomes, [Moustafa et al. 2017](https://doi.org/10.1371/journal.ppat.1006292) found sequences of 19 human viruses in 42% of the people, and judged most of the other 75 viruses they saw to be contamination of reagents and the environment.

## Why

The reads a human alignment leaves unmapped are already in the BAM and nothing else in the pipeline reads them. Moustafa et al. 2017 and Liu et al. 2020 both took exactly these reads from whole-genome sequencing of blood to find viral DNA and integrated HHV-6.

## Tool

- samtools to take the reads out of the BAM and to measure the depths
- minimap2 (`-ax sr --frag=no`) to realign them to the viral genomes
- `og viral report` (OwnGenome's Rust binary, [og](rust.md)) for the calls

## Docker Image

- `SAMTOOLS_IMAGE`, `MINIMAP2_IMAGE` and `OG_IMAGE` in the bash step; `MINIMAP2_SAMTOOLS_IMAGE` (the same minimap2 and samtools releases in one image) and `OG_IMAGE` in the Nextflow processes

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags. The step adds no image.

## Reference data

`./scripts/setup.sh --viral-data <genome_dir>` downloads the 24 NCBI RefSeq records of `VIRAL_REF_ACCESSIONS` (versions.env) as one FASTA from NCBI's E-utilities and checks it against `VIRAL_REF_SHA256`; it is about 1.6 MB, stored as `reference/viral-<VIRAL_REF_VERSION>/viral_genomes.fa`. An accession with its version number names one fixed sequence. The download was checked on 2026-10-25: two fetches gave the same file, with 24 sequences.

| Group | Genomes (RefSeq) |
|---|---|
| HHV-6 | HHV-6A (NC_001664.4), HHV-6B (NC_000898.1) |
| Other herpesviruses | EBV (NC_007605.1), CMV (NC_006273.2), HSV-1 (NC_001806.2), HSV-2 (NC_001798.2), VZV (NC_001348.1), HHV-7 (NC_001716.2), HHV-8 (NC_009333.1) |
| Papillomaviruses | HPV-16 (NC_001526.4), HPV-18 (NC_001357.1), HPV-6b (NC_001355.1), HPV-11 (NC_001525.1) |
| Hepatitis and retroviruses | HBV (NC_003977.2), HCV (NC_004102.1), HIV-1 (NC_001802.1), HTLV-1 (NC_001436.1) |
| Polyomaviruses | BK (NC_001538.1), JC (NC_001699.1), Merkel cell (NC_010277.2) |
| Others | adenovirus 5 (AC_000008.1), Torque teno virus 1 (NC_002076.2), parvovirus B19 (NC_000883.2) |
| Lab control | phiX174 (NC_001422.1), which Illumina adds to sequencing runs |

HCV is an RNA virus with no DNA stage; its row stays because Moustafa et al. 2017 saw HCV DNA once, which they read as an integration.

## Input

- `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` with its `.bai` (`ALIGN_DIR` as for the other BAM steps). A VCF alone has no reads, so a sample without a BAM is left out.
- the viral FASTA above (`VIRAL_REF` overrides its path, inside `GENOME_DIR`)

## Command

```bash
./scripts/47-viral.sh your_name
```

With `run-all.sh`: `TOOLS=...,viral`; a run lists step 47 only when `TOOLS` names it. In the Nextflow pipeline: `viral` in `--tools` with `--viral_ref reference/viral-1/viral_genomes.fa` (processes `VIRAL_READS` and `VIRAL_REPORT`).

## How it works

1. The reads with the unmapped flag (`samtools view -f 4`), which includes a read unmapped beside a mapped mate, and the reads mapped to chrEBV, the EBV genome the GRCh38 analysis set carries. `samtools fastq` leaves out secondary and supplementary records, so each read is realigned once.
2. minimap2 aligns them to the viral genomes as single reads (`--frag=no`, so two mates that sit next to each other in the FASTQ are not paired), with no secondary alignments.
3. `samtools coverage` gives each viral genome's length and its aligned records, `samtools depth -aa` the depth at every base, zero included, and `samtools view -F 4` the alignments themselves.
4. `samtools coverage` over the sample's BAM gives the human depth: chr1 to chr22, over the bases with at least one read. Both coverages count duplicates and leave out unmapped, secondary and QC-failed reads, so the two depths are measured the same way.
5. `samtools idxstats` of the sample's BAM gives its unplaced unmapped reads, the pairs where neither read mapped. A whole-genome BAM always has many. None means the BAM was cut to regions or filtered (`samtools view -F 4`), so the step had nothing to look in.

`og viral report` counts a read as evidence for a virus only when all of these hold:

- it is the read's primary alignment, not a supplementary, secondary or QC-failed one;
- its mapping quality is 20 or more, so it does not sit equally well on two genomes, such as the segments HHV-6A and HHV-6B share, or on a repeat;
- it is not low complexity: its DUST score is 7 or less. The score is PRINSEQ's ([Schmieder and Edwards 2011](https://doi.org/10.1093/bioinformatics/btr026)): the DUST score of [Morgulis et al. 2006](https://doi.org/10.1089/cmb.2006.13.1028) over 64-base windows every 32 bases, scaled from 0 to 100. A random sequence scores about 1.5, a dinucleotide repeat about 49 and a run of one base 100; 7 is the cut-off PRINSEQ's authors suggest. Poly-A tails, telomere repeats (TTAGGG) and other simple repeats do not count;
- it is not a duplicate: no counted read on the same genome and strand starts at the same unclipped 5' position. Copies of one DNA molecule made by PCR count once.

A supplementary record of counted quality adds to the breadth, not to the reads, so a read across the joining point of a circular genome counts once and still covers both ends. The integration call does not use this selection: it reads the depth of every primary alignment, as the human depth does.

Per genome, `og viral report` writes:

| Column | Meaning |
|---|---|
| `reads` | counted reads, as above; every aligned record is in `_viral_coverage.tsv` |
| `breadth` | the fraction of the genome covered by a counted read |
| `mean_depth`, `mode_depth` | the mean depth, and the most common depth over the genome's bases, from every primary alignment |
| `relative_depth` | `mode_depth` / human depth |
| `copies_per_cell` | for an integration, 2 x `mode_depth` / human depth (a cell holds two copies of each autosome); for viral DNA found and phiX, 2 x `mean_depth` / human depth; `.` for the other calls, where the depth is a pile on a few bases and no count per cell |
| `call` | see below |

The calls:

- **not assessed**: every genome, when the BAM has no unplaced unmapped read. Nothing could be found, so nothing is ruled out.
- **inherited integrated HHV-6 likely**: HHV-6A or HHV-6B with `relative_depth` of 0.30 or more. This is the rule of [Liu et al. 2020](https://doi.org/10.1371/journal.pgen.1008915): the mode of the per-base depth over the viral genome against the genome's sequencing depth, here the mean depth over chr1 to chr22, with integration at 0.30 or more. One integrated copy gives about 0.5; viraemia gives 10 to 1000 times less, by the studies they cite. The mode is used, not the mean, so a pile of reads on a short repeat (HHV-6's ends carry the same TTAGGG repeat as human telomeres) cannot lift a genome that has no reads elsewhere. When both species pass, both are called: one integrated genome does not bring the other species near this depth, as the reads they share are a small part of each genome. A laboratory test tells them apart.
- **lab control (phiX)**: any counted phiX read. It was added in the lab and says nothing about the sample.
- **viral DNA found**: at least 10 counted reads, over at least 10% of the genome.
- **a few reads**: at least 10 counted reads, over less than 10% of the genome.
- **below threshold**: 1 to 9 counted reads. The note gives the count. Too few to tell the virus from contamination at this depth.
- **not found**: no counted read. When the genome has aligned records that did not count, the note says how many.

The 10 reads and 10% are this pipeline's choices, not published values. The 10% is a share of each genome's own length, so the reads it takes grow with the genome. With 150-base reads that do not overlap, 10 reads cover 10% of a genome up to 15 kb: HPV, HBV, HIV-1, the polyomaviruses, parvovirus B19 and TTV need only the 10 reads. Adenovirus 5 (36 kb) needs about 24, and the herpesviruses (125 to 236 kb) about 84 (VZV) to 158 (CMV). On a herpesvirus, "a few reads" can therefore stand for over a hundred reads; its `reads` and `breadth` columns say how many.

The step prints one line per sample: each genome with a call other than not found, with its count when it is below the threshold, or "no reads above the threshold at this depth".

## Output

| File | Contents |
|---|---|
| `viral/${SAMPLE}_viral.tsv` | one row per viral genome: `virus`, `accession`, `length`, `reads`, `breadth`, `mean_depth`, `mode_depth`, `human_depth`, `relative_depth`, `copies_per_cell`, `call`, `note` |
| `viral/${SAMPLE}_viral_reads.bam` (`.bai`) | the realigned reads, on the viral genomes; open it in IGV with the viral FASTA |
| `viral/${SAMPLE}_viral_coverage.tsv` | `samtools coverage` of that BAM |
| `viral/${SAMPLE}_human_coverage.tsv` | `samtools coverage` of the sample's BAM |

Both reports show the HHV-6 call and every genome with reads in their "Viral DNA in the unmapped reads" section.

## Runtime

Not measured on a whole genome. samtools reads the BAM twice from end to end, once for the unmapped reads and once for the human depth; minimap2 then aligns only the unmapped reads, to 24 small genomes.

## Notes

- Not a diagnosis. A virus found means its DNA was in the sequenced sample, which can come from latent virus in blood cells, a past or current infection, or a reagent.
- What integrated HHV-6 means for health is not settled. [Pellett et al. 2012](https://doi.org/10.1002/rmv.715) list what is known: it can make blood tests read as an active infection, transplant recipients with it may be at higher risk, and its role in disease needs more research. Each child of a carrier has about a one-in-two chance of inheriting it, as for any variant on one chromosome.
- A finding is confirmed by a test in a clinical laboratory; Liu et al. 2020 used droplet digital PCR.
- [Liu et al. 2020](https://doi.org/10.1371/journal.pgen.1008915) also describe people who carry only one direct repeat (DR) of HHV-6B, a "solo-DR": reads on the genome's ends only. Such a pattern gives a low mode depth here and is not called integrated; look at the BAM when HHV-6 has reads but no call.
- A negative is narrower than "no virus". "Not found" and "below threshold" mean no counted read, or too few, at this sequencing depth, among the 24 genomes, at the strains of their RefSeq records. A strain far from its record can be missed.
- Only the reads the human alignment left unmapped are realigned, and the reads on chrEBV. Viral sequence that maps to the human genome stays there: the human end of an integration junction, and HHV-6's TTAGGG ends, which map to the telomeres.
- A BAM cut to a few regions, or filtered to its mapped reads, has no unplaced unmapped read, and the step says "not assessed" for every genome.
- On the e2e fixture the case plants synthetic reads drawn from the HHV-6B genome at half the human depth as unmapped reads, and the call is integration; at a tenth of the human depth over the same genome it is not. It plants CMV at a low level and phiX too, and three decoys that must not count: three HBV molecules read four times each, which give 3 counted reads and "below threshold"; simple repeats (poly-A, TTAGGG and AT); and the HHV-6A reads an HHV-6B genome brings with it, which read "not found". A copy of the BAM without its unmapped reads reads "not assessed". Nothing in the case is a finding about the fixture's own sample.
