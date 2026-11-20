# Step 41: Blood Groups (predicted)

## What This Does

Predicts the red cell, platelet and neutrophil antigens the genome encodes: ABO (with A1 and A2), RhD, RhCE (C/c, E/e), Kell (K/k, Kp, Js, and the Kell-null K0), Duffy (with the GATA-box variant), Kidd (with the Kidd-null Jk(a-b-)), MNS (M/N, S/s), Lutheran, Diego, Dombrock, Colton, Vel, Langereis (Lan), Junior (Jra), H (with Bombay and para-Bombay), secretor status, the platelet antigens HPA-1 to HPA-35, and the neutrophil antigens HNA-1 (a, b, c), HNA-3, HNA-4 and HNA-5. Opt-in, and the last optional step.

The result is a prediction from the DNA, not a serological typing. It is information about the sample, and the table states only what the variants predict. For a platelet or neutrophil antigen the genome predicts the person lacks, the row's note also says what lacking it means (see [What a platelet or neutrophil antigen means](#what-a-platelet-or-neutrophil-antigen-means)). That is information about the antigen, not advice.

A person who lacks an antigen almost everyone carries (Vel, Lan, Jra, H, k and the others the table marks as high-prevalence) is rare, and their row comes first in the table.

## Why

Most blood group antigens differ by one known variant, so the genome predicts them. Nothing else in the pipeline reports them, and working out even ABO by hand from a VCF means knowing which reference allele GRCh38 carries (it carries the O allele's c.261delG, so an A or a B shows up as an insertion).

## Tool

- `og bloodtype report` (OwnGenome's Rust binary, `rust/crates/og-bloodtype`), applying the marker table `assets/bloodtype_markers.tsv`. It reads the marker sites through the VCF's `.tbi` or `.csi`; without a usable index (none, or one older than the VCF) it reads the VCF once, line by line, with the same result. `og bloodtype check-markers` checks that the table and the rules agree
- bcftools `norm` to left-align the VCF's records near the marker sites against the reference first (see [Indels written another way](#indels-written-another-way))
- samtools `bedcov` for the read depth over the marker sites and RHD; the RHD regions are the RHD row of `bin/cn_loci.tsv`, and og cn's engine turns their depth into the copy number (see [Copy number from depth](cn-loci.md))

The variants and the antigens they define are the ISBT allele definitions, as RBCeq2 2.4.4 packages them (database v2.5.1, MIT licence, Australian Red Cross Lifeblood). The table's header names the release and the checksum of the file they were read from. Each position and REF base was checked against Ensembl GRCh38 and its dbSNP id. The HPA positions come from dbSNP through Ensembl, because that RBCeq2 database lists them at their GRCh37 positions.

The null alleles (an allele that makes no antigen) are the single-variant null alleles RBCeq2 lists for SMIM1 (Vel), ABCB6 (Lan), ABCG2 (Jra), FUT1 (H), FUT2 (secretor), KEL (K0) and SLC14A1 (Jk(a-b-)): 205 rows. Null alleles at one position with the same REF share a row, its ALT then a comma list. Left out are the alleles of two or more variants, the two deletions over 50 bp, one KEL allele that shares its position with another under a different REF, and three alleles too common to be a null on their own (see Notes). HPA-7 to HPA-35 are RBCeq2's HPA alleles. That database gives them at GRCh37 positions, so each was mapped to GRCh38 with Ensembl's assembly map. Every new row was then checked on GRCh38 on 2026-11-03: Ensembl VEP on the row's REF and ALT gives the allele's c. change on its RefSeq transcript, and the REF is Ensembl's sequence. Where RBCeq2's text has a slip, the table keeps VEP's change: HPA-13 is ITGA2 c.2483C>T, HPA-16 is ITGB3 c.497C>T and HPA-35 is ITGB3 c.1514G>A. The IPD-HPA database (https://www.ebi.ac.uk/ipd/hpa/) names the same antigens.

The neutrophil antigens are not in RBCeq2. Their variants are those of the ISBT HNA allele table (Flesch et al. 2016, Transfusion 56:1477, https://doi.org/10.1111/trf.13575). Each was mapped to GRCh38 with Ensembl VEP on the gene's transcript, which also checks the REF base (2026-10-25). dbSNP places some FCGR3B variant ids on FCGR3A, its near copy (rs447536, for one, maps to FCGR3A on GRCh38), so the FCGR3B positions in the table are FCGR3B's own, in exon 3.

Why not run RBCeq2 or bloodTyper itself: RBCeq2 is on PyPI only, with no Bioconda package and no container image, and its README says it does not type RH from ordinary short-read calls (it needs long reads or a caller that writes copy numbers). bloodTyper is on neither PyPI nor Bioconda; the one Docker Hub image is a personal upload. A small table applied in the pinned Python image needs no new image and reads RhD from depth.

## Docker Image

- `PYTHON_IMAGE` for the rules, `SAMTOOLS_IMAGE` for the depth, `BCFTOOLS_IMAGE` for the left-alignment. In the Nextflow pipeline the `BLOODTYPE` process runs both the left-alignment and the rules in `BCFTOOLS_IMAGE` (`og` is a static binary)

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags.

## Input

- `${SAMPLE}/vcf/${SAMPLE}.vcf.gz`: one sample, from any caller. Its records near the marker sites are left-aligned first, so an indel need not be written the way the table writes it. Split and whole multi-allelic records, all-sites VCFs (a reference record and then the variant at one position) and gVCF reference blocks are read too: every record that covers a site counts, whatever their order (see below)
- the reference (`REF_FASTA`, `reference/GRCh38_no_alt_analysis_set.fasta` by default) with its `.fai`, for the left-alignment; `--reference` in the Nextflow pipeline
- `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` with its `.bai` (`ALIGN_DIR` as for the other BAM steps). Optional: without it RhD is not determined, and a site the VCF has no record for is read as the reference allele at confidence medium. A system with no VCF record at any of its sites is then not assessed (see [Without the BAM](#without-the-bam)).

## Command

```bash
./scripts/41-bloodtype.sh your_name
```

With `run-all.sh`: `TOOLS=...,bloodtype`. It is opt-in and needs no data, so a run lists step 41 only when `TOOLS` names it. In the Nextflow pipeline: `bloodtype` in `--tools` (processes `BLOODTYPE_REGIONS`, `BLOODTYPE_DEPTH` and `BLOODTYPE`).

## How each system is read

| System | Variants | Rule |
|---|---|---|
| ABO | c.261delG, c.802G>A, c.803G>C, c.796C>A, c.1061delC | GRCh38 carries c.261delG (O.01): each insertion there is a functional allele. c.802G>A makes one of them O.02, c.803G>C makes one B, the rest are A; c.1061delC makes an A an A2. A1 with A2 reads A1 |
| RhD | read depth over RHD; c.487-19_504dup with c.609G>A | Copies = 2 x RHD depth / depth of the two flanks outside the RHD deletion (0 under 0.5, 1 under 1.5, else 2). The RHD pseudogene (RHD*08N.01) is an inactive copy. Exons 3 to 7 without RHD reads in a gene that is there is the pattern of an RHD-RHCE hybrid gene: RhD is then not determined |
| RhCE | c.307C>T (C/c), c.676G>C (E/e) | one allele per copy of the variant |
| Kell | c.578C>T (K/k), c.841C>T (Kp), c.1790T>C (Js) | the same |
| Duffy | c.125G>A (Fya/Fyb), c.-67T>C (GATA box), c.265C>T (FY*X) | The GATA box stops expression on red cells. With one Fyb and one GATA box, both phased calls in one phase set (the same `PS`, or no `PS` on either), it goes on the chromosome the VCF gives: on the Fya chromosome that is Fy(a-b+). Else it is placed on a Fyb chromosome first, where it usually lies (FY*02N.01). FY*X gives a weak Fyb (`b+w`) |
| Kidd, MNS, Lutheran, Diego, Dombrock, Colton | one variant each pair (c.838G>A; c.59T>C and c.143C>T; c.230G>A; c.2561C>T; c.793G>A; c.134C>T) | one allele per copy |
| HPA-1, -2, -3, -4, -5, -6, -15 | ITGB3 c.176T>C, GP1BA c.482C>T, ITGA2B c.2621T>G, ITGB3 c.506G>A, ITGA2 c.1600G>A, ITGB3 c.1544G>A, CD109 c.2108C>A | one allele per copy; GRCh38 carries HPA-15b. HPA-6b is a low-frequency antigen, written HPA-6bw in the ISBT table |
| HNA-1 | FCGR3B c.108G>C, c.114C>T, c.244G>A, c.316G>A (a or b); c.194A>G (checked); c.233C>A (c) | FCGR3B\*01 (HNA-1a) and \*02 (HNA-1b) differ at five bases. GRCh38 carries the \*02 base at four of them and the \*01 base at c.194, so the copies of HNA-1a are counted on the four, which must agree (else not determined), and c.194 is only checked: a c.194 that does not match is named in the note. FCGR3B\*03 (HNA-1c) is \*02 with c.233C>A and carries HNA-1b too. Without a call at c.233 an HNA-1b allele reads \*02 or \*03. Confidence is always low: FCGR3B and FCGR3A are near copies, and FCGR3B has from none to four or more copies, which the step does not count |
| Vel, Langereis, Junior | SMIM1 c.64_80del; 36 ABCB6 null rows; 23 ABCG2 null rows | Read from the null rows alone: two null copies lack the antigen (Vel-, Lan-, Jr(a-)). One copy is a carrier, still positive. Two different null alleles, each once, are taken to sit on different chromosomes (no phase), so the antigen is lacking at confidence low |
| H, Secretor | 34 FUT1 null rows; 14 FUT2 null rows (c.461G>A among them, often written c.428G>A) | The same rule. Two FUT1 null copies give H- on red cells: `H- (Bombay, Oh)` with a non-secretor, `H- (para-Bombay)` with a secretor. The ABO row then says the red cells carry little or no A or B and usually group as O, whatever the ABO alleles are; an ABO row that reads `not assessed` or `not determined` has no alleles to qualify and gets no such note |
| Kell null (K0), Kidd null (Jk(a-b-)) | 57 KEL null rows; 40 SLC14A1 null rows, c.342-1G>A and c.871T>C among them | Each null copy silences one chromosome. With no phase it is placed on a chromosome with the antigens of the row's `null_on` column (k, Kpb and Jsb for KEL; Jka or Jkb, as the ISBT table finds it): one KEL null with K/k reads K+k-. Two null copies give K-k- Kp(a-b-) Js(a-b-) (lacking Ku) or Jk(a-b-) (lacking Jk3), also when an antigen site has no call. Two different nulls with different `null_on` alleles are each named with their own in the note |
| HPA-7 to -14, -16 to -35 | one variant each in ITGB3, ITGA2B, ITGA2, GP1BB or GP9 | one allele per copy; the a antigen is on the reference allele, the low-prevalence b antigen on the other |
| HNA-3, -4, -5 | SLC44A2 c.461G>A, ITGAM c.230G>A, ITGAL c.2372G>C | one allele per copy; GRCh38 carries HNA-3b |

A site the VCF has no record for is the reference allele when the BAM has at least 8 reads there at MAPQ 20, and not determined below that. Without the BAM it is read as reference, and the confidence says so (see [Without the BAM](#without-the-bam)). The step leaves out of `samtools bedcov` the rows on a contig the BAM has no `@SQ` line for (a BAM cut to a few chromosomes), because `bedcov` would stop on them. Those sites then have no depth row: the ones with no VCF record are not determined, and RhD is not determined without chr1. A record at a site whose REF is not the table's REF is a no call, because it says nothing about the table's ALT. There are two exceptions. A 0/0 call over a longer REF that starts with the table's REF is the reference allele. A 0/0 record over a shorter REF (an all-sites VCF's `C . 0/0`, or a gVCF block that starts at a deletion marker such as ABO c.1061delC) says nothing past its own base, so the site reads as one with no record: the reads decide. A gVCF reference block (ALT `.`, `<*>` or `<NON_REF>`) whose `INFO/END` reaches the table's last REF base is a reference call there, wherever the block starts, with the block's `GQ` and `DP` (else `MIN_DP`) (see [Weak calls](#weak-calls)).

### Without the BAM

A VCF lists where the genome differs from the reference. A site with no record is either the reference allele or a site the VCF does not cover, and only the reads tell the two apart. So without the BAM:

- A system with a record at one of its sites at least is read, and its sites with no record count as the reference allele at confidence medium. The record shows that the VCF covers the gene.
- A system with no record at any of its sites is `not assessed` (phenotype and confidence), and the note names the sites' chromosome. A VCF with no chr9 record would otherwise read ABO as O, the group GRCh38 carries.

A person who matches GRCh38 at every site of a system (group O at ABO, HPA-1a/1a, most null rows) also has no record there in a variants-only VCF, so that system is not assessed either. A gVCF covers every site with a record or a reference block, and needs no BAM for this. RhD needs the BAM in every case.

### Weak calls

A call keeps its genotype, but it is weak when its genotype quality (FORMAT `GQ`) is under 20 or its read depth (FORMAT `DP`, else a gVCF block's `MIN_DP`) is under 8, the same depth a site with no record needs. A system that uses a weak call has confidence low, and its note names the site and the values (`KEL c.578C>T: weak call (GQ 0 under 20, DP 1 under 8)`). The evidence file adds the same values to the site's note. A record without `GQ` or `DP` is not checked. The allele balance of a call (the share of reads that carry the ALT) is not checked.

### Indels written another way

The table writes each indel left-aligned, the way `bcftools norm` writes it: ABO c.1061delC is `CG C` at chr9:133255669, at the start of a run of G. A caller can write the same deletion further right in the run, as `GG G` at chr9:133255671 (DeepVariant does so in rare cases without `--normalize_reads`, which the pipeline does not pass, and a VCF from another caller may do so too). Matched by position as written, that record would miss the marker, and with reads there the site would read as the reference allele at confidence high. So the step first takes the VCF's records within 100 bp of each marker site and left-aligns them against the reference (`bcftools norm -c w -f`). 100 bp is more than the 38 bp the longest marker indel can shift along its repeat on GRCh38. The rules then read that file.

- A multi-allelic record stays one record (no `-m -any`). A `C A,G 1/2` at ABO c.803 keeps its one B copy, and a null row's `1/2` of two of its ALTs keeps both copies.
- A multi-allelic indel is left-aligned as a whole, so its REF can stay longer than the marker's (`TCC TC,T` against `TC T`). Its REF then differs from the table's, and the site is a no call, as the paragraph above says.
- A record whose REF differs from the reference is kept as it is (`-c w`), with a warning in the log.
- A record on a contig the reference lacks (a VCF named `1`, `2`, ... against the chr-named reference) cannot be aligned. It goes through unchanged.
- bcftools reads the VCF through its index when the index htslib loads (the `.csi` if there is one, else the `.tbi`) is not older than the VCF. Otherwise it reads the VCF once from the start. Read that way, a VCF whose header has no `##contig` line for a marker's contig first gets one, because bcftools writes no record on an undeclared contig. That path also reads a VCF in plain gzip, not bgzip.
- The left-aligned file is a work file: the bash step deletes it after the rules run, and the pipeline does not publish it.

### Several records at one site

When several records cover one site, each one is read and their order does not matter:

- A call of the table's ALT wins over the other records at the site's position. An all-sites VCF writes `T . 0/0` and then `T TC 0/1` at ABO c.261, so the insertion counts. Another allele that a second record calls there (a split `T G 1/0`) is named in the note.
- When one null marker lists several ALTs, calls of distinct listed ALTs add their copies. Splitting the record keeps the same evidence and phenotype. Repeated calls of the same ALT count once, even if its allele index changes. Different dosages for the same ALT, or a total above two copies, make the site a no call with an incompatible-dosage note.
- A no call of a record that names the table's ALT (`./.`, or a failed FILTER) comes next: the ALT may be there, so a reference record beside it does not make the site reference.
- Then a reference call, then any other no call.
- A record that starts upstream and spans the site with an allele that changes the site's base (a deletion over it, or a padded record such as `CG CA`) makes the site a no call, whatever the site's own records say. A reference call or a no call upstream changes nothing, and neither does an allele that keeps the site's base.

## Output

| File | Contents |
|---|---|
| `bloodtype/${SAMPLE}_bloodtype.tsv` | one row per system: `system`, `phenotype` (for example `O`, `D+`, `Fy(a+b-)`, `K-k+ Kp(a-b+) Js(a-b+)`), `genotype` (the alleles, for example `O.01/O.01`, `RHD/RHD deleted (RHD copies by read depth 0.91)`), `confidence`, `note` |
| `bloodtype/${SAMPLE}_bloodtype_evidence.tsv` | one row per site (status: called, covered, assumed, low depth or no call; alt copies; read depth) and per RHD region (mean depth and copies) |
| `bloodtype/${SAMPLE}_bloodtype_depth.tsv` | samtools `bedcov -Q 20 -j` over the marker sites and the RHD regions of `bin/cn_loci.tsv` (with a BAM). A marker row keeps the marker table's columns, name in column 6, then the depth sum; an RHD row is a 4-column BED (chrom, start, end, name) plus the depth sum |

Confidence:

- **high**: every site was called in the VCF or covered by reads
- **medium**: a site was read as reference with no BAM depth to confirm it, or one of the system's null rows was not determined (for example too few reads there)
- **low**: a site in a region where short reads are unreliable, named in the note (RHCE exon 2 for C/c, which is RHD-like on a C allele; GYPA, GYPB and GYPE, which are near copies), markers that disagree, another allele at a site of a two-antigen pair, two different null alleles taken to sit on different chromosomes, a weak call (see [Weak calls](#weak-calls)), or an ABO A allele (no variant shows A: it is a functional allele with neither c.802G>A nor c.803G>C, so a rare O allele that keeps c.261 or a weak A subgroup reads A too)
- **not determined**: a needed site had no call or too few reads; RhD without the BAM, with flanks under 10x, or with the hybrid pattern
- **not assessed**: without the BAM, no VCF record at any of the system's sites

Row order: first the systems where the genome predicts the person lacks a high-prevalence antigen, then the others in the table's order. The high-prevalence antigens are the ones in the marker table's `high_prevalence` column: ISBT's high-prevalence antigens Vel, Lan, Jra, H, k, Kpb, Jsb, Ku, Jk3, Lub, Dib and Coa, and for platelets HPA-1a and the a antigen of HPA-4 and of HPA-6 to HPA-35 except HPA-15, whose b antigen is a low-prevalence one. Such a row's note ends with `lacks <antigen>, carried by most people (high prevalence)`. A genome with only the reference alleles lacks none of them, so its rows keep the table's order.

Both reports show the table as their "Blood Groups (predicted)" section: system, phenotype, genotype and confidence. The note, with the meaning, is in the TSV.

## What a platelet or neutrophil antigen means

When the genome predicts a person lacks one of these antigens (both copies carry the other allele, for example HPA-1b/1b), the note of that row starts with `lacks HPA-1a:` and the meaning from the table's `meaning` column:

- **Platelet antigens (HPA).** A person can make antibodies to a platelet antigen they lack, after a pregnancy or a transfusion. Those antibodies cause fetal and neonatal alloimmune thrombocytopenia (FNAIT: the mother's antibodies destroy the baby's platelets), post-transfusion purpura and refractoriness to platelet transfusions. HPA-1a is the antigen most often behind FNAIT and post-transfusion purpura in white populations, so an HPA-1a-negative (HPA-1b/1b) mother can make anti-HPA-1a against an HPA-1a-positive baby. Source: Curtis and McFarland, Human platelet antigens 2013, Vox Sanguinis 106:93, https://doi.org/10.1111/vox.12085
- **Neutrophil antigens (HNA).** The same holds for the antigens on neutrophils: antibodies after a pregnancy or a transfusion cause neonatal alloimmune neutropenia and transfusion-related acute lung injury (TRALI). Sources: Flesch et al. 2016 (above) and Esmaeili et al. 2022, J Immunol Res 2022:4387555, https://doi.org/10.1155/2022/4387555

Lacking an antigen is common (most people lack HPA-1b, for example). It means a person could make that antibody in those situations, not that they have it.

## Runtime

Not measured on a whole genome. bcftools reads the VCF's records near the marker sites through its index (the VCF once from the start without one) and left-aligns them; og reads the marker positions of that small file; samtools reads the BAM at 284 small regions.

## Notes

- RhD here is the number of RHD copies, not their kind: weak D and partial D (RHD missense alleles) are not assessed, and neither are RHD null alleles other than RHD*08N.01, the RHD-CE-D hybrids that keep exons 3 to 7, and DEL alleles. The common RhD-negative allele, the deletion of the whole gene (RHD*01N.01), is what the depth sees. A D+ row's note says so.
- Phase is read for Duffy only (see the table). ABO assumes A and B sit on different chromosomes, and two different null alleles of one system are taken to sit on different chromosomes.
- Rare and null alleles outside the table (null alleles of two or more variants, ABO's weak subgroups, cis-AB) are not assessed. A high-prevalence antigen reads positive when none of the table's null alleles is found, so a null allele the table does not hold goes unseen: the note of each system read from null rows (Vel, Langereis, Junior, H, secretor, Kell, Kidd) says so while the antigen reads present. The ABO row's note names the ABO alleles outside its five variants. Neither are HNA-2 (CD177), HNA-1d and HNA-1null (no FCGR3B gene) assessed.
- Every null row's allele was looked up in gnomAD v4 (genomes and exomes, 2026-11-03): 214 alleles, 53 of them not in gnomAD at all. Three were dropped because they are far more common than their null phenotype allows: FUT2 c.772G>A (rs602662, about 0.5 of chromosomes) and c.335C>T (rs200157007, 0.14 in South Asians) sit on the same haplotype as the common non-secretor allele c.461G>A, so with no phase they turned Se/se carriers into non-secretors; SLC14A1 c.810G>A (rs17675299, 0.08 in Finns) is far too common for a Kidd null. The rows kept that pass 1% in some population are known nulls of that population: FUT2 c.461G>A (the common non-secretor allele), SMIM1 c.64_80del (Vel), ABCG2 c.376C>T (Jr(a-), East Asians), SLC14A1 c.871T>C (Finns), and FUT2 c.812delC (rs1799761, 0.015 in African genomes).
- Another ALT at a null row's position (for example G>T at the position of KEL c.1546C>T, rs8176034) is not that null allele: the note names it, and the confidence stays as the other sites set it.
- Secretor status: weak secretor alleles (such as FUT2 c.418A>T, often written c.385A>T) and a deletion of the whole FUT2 gene are not assessed; the row's note says so.
- Lewis (FUT3) is not read. GRCh38 carries FUT3 c.202T>C and c.314C>T, the variants of a common Lewis-null allele, so a working FUT3 shows up as alternative alleles whose phase a VCF does not give.
- On the e2e fixture the case reads HG002's reads at the marker sites from GIAB's 60x BAM and HG002's genotypes from GIAB's v4.2.1 benchmark, because the fixture's slices hold none of these genes. Those published calls give O (no c.261 insertion), and so do the reads on their own (under 10% of them carry the insertion there); the reads give two RHD copies (D+). A copy of the same reads without any RHD read gives D-, and with half of them D+ with one copy deleted.
