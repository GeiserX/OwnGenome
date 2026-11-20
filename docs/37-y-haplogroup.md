# Step 37: Y-Chromosome Haplogroup

## What This Does

Assigns the Y-chromosome haplogroup of a male sample: the branch of the paternal line's family tree its Y variants place it on. Opt-in.

## Why

The mitochondrial haplogroup (step 12) traces the maternal line only. The Y chromosome passes from father to son almost unchanged, so its haplogroup traces the paternal line. Like step 12, it is ancestry, not health: no Y haplogroup is a diagnosis.

## Tool

- **Yleaf** 4.1.7 (Ralf et al., Mol Biol Evol 2018; GPL-3.0), Erasmus MC. It reads the BAM's pileup at its Y markers and predicts the haplogroup they support on its default tree, YFull v14.01 (about 183,000 marker positions on GRCh38).

- **The current YFull YTree** (release 14.06.0 of 22 September 2026, CC BY 4.0, [github.com/YFullTeam/YTree](https://github.com/YFullTeam/YTree)), walked by `og y tree-place` (`rust/crates/og-reports/src/ytree.rs`). It is newer than Yleaf's v14.01, and the second position table below places SNPs Yleaf's own table does not have, so the branch it reports can sit below the one Yleaf reports. It also writes the path from the root with the calls of every level. The SNPs on YFull's tree have names only; their GRCh38 positions and bases come from two tables of Yleaf 4.1.6 (GPL-3.0): `ftdna_positions_hg38.txt` (made for FamilyTreeDNA's tree, the widest of its tables) first, then `new_positions.txt` (the one Yleaf 4 uses for its YFull v14 tree). SNP names are shared between the trees, so a name is looked up as it stands. Together they place 85% of the tree's 459,468 SNPs (76% from the first table alone), counted on 2026-10-25.

The step reports both results, each with the tree it comes from. Both reports (step 24, text and HTML) show them as two sections side by side: "Y-Chromosome Haplogroup (Yleaf)" with Yleaf's call, and "Y Placement on the Current Public Tree" with this one, its level, the derived, ancestral and no-call counts along the path, any ancestral SNP inside the terminal branch or outvoted above it, and one row per level of the path.

## Docker Image

- `YLEAF_IMAGE`, and `SAMTOOLS_IMAGE` for the pileup

The step runs in three parts: Yleaf's marker positions, written by `og y yleaf-positions` in its image, `samtools idxstats` and `samtools mpileup -l <positions> -ABQ20q1` (Yleaf 4's own flags at its default quality 20; samtools reads `20q1` as the `-Q` value, so that pileup has no mapping-quality filter; it stays as Yleaf writes it) in `SAMTOOLS_IMAGE`, in the same run as the chrY BAM (step 37 makes the current tree's pileup in a second run, with its flags written out: `-A -Q20 -q20`, so reads of mapping quality below 20 do not count; Nextflow's `Y_PILEUP` task makes all three in one task), then Yleaf's own command line on that pileup. The step puts the pileup where Yleaf keeps its own, `<out>/<bam name>/temp_haplogroup_pileup.pu`, and runs `Yleaf -bam <bam> -rg hg38 --ref-fasta <reference> -tree yfull -q 20 -t 1 -rp -force --no-update-check`; `-rp` makes Yleaf read that file instead of running samtools. The image test (the `YLEAF_IMAGE` row of `tests/smoke/commands.tsv`) checks that the image still has these options, that pileup file name and those samtools flags, and that its marker table and tree are the ones in the archive `setup.sh --yleaf-data` installs. Yleaf stops on an empty pileup (no read at any marker) and writes no row; the step then writes the row Yleaf gives for too few markers (`Hg` NA) with `og y yleaf-no-reads`.

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Input

- `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` with its `.bai` (`ALIGN_DIR` as for the other BAM steps)
- `${SAMPLE}/vcf/${SAMPLE}.vcf.gz` from step 03, for the private-SNP list only. Without it the rest of the step runs and that one file is skipped with a line saying so.
- `${SAMPLE}/indexcov/indexcov-indexcov.ped` from step 16: the sex indexcov infers from the reads. The step runs only when it is male, and asks for no sex of its own. A female or undetermined sample is skipped with one line and exit 0.

In the Nextflow pipeline (`y_haplogroup` in `--tools`) the process runs for the rows whose samplesheet sex is male. `INDEXCOV` has checked that sex against the reads before any BAM step starts, and a mismatch stops the run (unless `--sex_check warn`).

## Command

```bash
./scripts/setup.sh --yleaf-data /path/to/genome_dir   # once: Yleaf's marker tables and trees, and the current tree (about 210 MB to download)
./scripts/37-y-haplogroup.sh your_name
```

`setup.sh --yleaf-data` takes Yleaf's data folder (the GRCh38 marker positions in `hg38/` and the haplogroup trees in `hg_prediction_tables/`, the same files the image ships) from GitHub's archive of the same release, tag `v<version>`, checked against `YLEAF_DATA_SHA256` in `versions.env`, into `reference/yleaf-<version>/data`; the pipeline takes it as `--yleaf_data`. The other builds in the archive (hg19, T2T) are left out. The same command then puts the current tree in `data/ytree`: YFull's `current_tree.json` at its release commit, and the two Yleaf 4.1.6 position tables at the commit of that tag, each checked against its sha256 in `versions.env` (`YTREE_*`). An install made before the current tree existed is completed by running the command again. If the current tree cannot be downloaded, the command warns and still exits 0, because Yleaf's own data is in place. Without `data/ytree` the step prints `Current tree: not installed` and writes Yleaf's result only.

With `run-all.sh`: `TOOLS=...,y_haplogroup` (it is opt-in, so a default run lists it as skipped; without the data it is skipped with the reason).

Yleaf downloads the whole hg38 FASTA on its first run unless it is given one. The step passes the pipeline's reference with `--ref-fasta` (Yleaf checks its chrY length against hg38 in the `.fai`; a BAM never needs the sequence itself) and `--no-update-check` turns off its check for a newer release, so nothing is downloaded and every container runs without a network. `--ref-fasta` takes only a file ending in `.fa`, `.fasta` or `.fna`, as `setup.sh` installs it, so step 37 stops with a message on any other `REF_FASTA` (a compressed FASTA, say).

## Output

| File | Contents |
|---|---|
| `y_haplogroup/${SAMPLE}_y_haplogroup.txt` | Yleaf's prediction: `Hg` (the haplogroup), `Hg_marker`, `Total_reads`, `Valid_markers` (Y markers with haplogroup information and enough reads), `QC-score` |
| `y_haplogroup/yleaf/` | Yleaf's working files: the markers it read with their alleles and its log |
| `y_haplogroup/positions.txt` | the marker positions the pileup was made at |
| `y_haplogroup/${SAMPLE}_y_tree.txt` | the placement on the current tree: `Hg` (the terminal branch, YFull's name), `Level` (its depth below the root), `Path_derived`, `Path_ancestral`, `Path_no_call` (the SNP calls of every branch on the path), `Ancestral_in_terminal`, `SNPs_called`, `Min_depth`, `Tree` and `Positions` (the tree release and the position tables it used), and `Outvoted_ancestral` (below; `.` when none). Then how much of the evidence was readable and what the ranking set aside: `Informative_callable` (the tree SNPs on the tables whose site has at least `Min_depth` reads, called or not), `Path_breadth` (that share of the path's SNPs on the tables, `.` without a path), `Conflict` (the positive branch on another lineage with the most support of its own from the common branch, above 0; `.` when none), `Conflict_support` and `Terminal_support` (derived minus ancestral calls of that branch and of the terminal one, each counted from their common branch down) |
| `y_haplogroup/${SAMPLE}_y_tree_path.tsv` | one row per branch from the root to the terminal one: `SNPs` (its defining SNPs), `On_tables` (those the tables place), `Derived`, `Ancestral`, `No_call`, and the names of its ancestral SNPs |
| `y_haplogroup/${SAMPLE}_y_private_snps.tsv` | the chrY SNVs no public SNP table names: `Pos`, `Ref`, `Alt`, `Depth`, `Alt_fraction`, `GQ`, `Genotype`, with the gates and the dropped counts in the header lines ("Private SNPs" below) |
| `y_haplogroup/${SAMPLE}_chrY.bam` and `.bai` | the BAM's chrY reads only, the file the services below read; the Nextflow pipeline publishes the same file there (`Y_PILEUP`) |

The chrY BAM is written every time this step runs on a male sample. Expect it to be about 1 to 2 percent of the analysis BAM, so a few hundred megabytes to a gigabyte on a 30x genome. That is an estimate from chrY's share of the reads; we have not measured it on a full genome. If you do not need it, delete it and its `.bai`.

How the placement works. A SNP is called when at least 10 reads cover it and 90% of them carry one of its two bases, derived or ancestral (Yleaf's own defaults, `-r 10` and `-b 90`); anything else is a no-call. A branch is positive when it has a derived SNP and more derived than ancestral ones, and negative when it has an ancestral SNP and no more derived than ancestral ones. A negative branch blocks the branches below it, with one exception: a branch with a single ancestral SNP is outvoted when the branches below it on the path carry derived SNPs at two or more distinct positions (two names for one site count once). Its SNP is then listed in `Outvoted_ancestral` and the step line names it. That is one site, not one read: a call already needs 10 reads and 90% of them on one base, so what it guards against is a wrong call at one site (a mismapped stretch, an error in the position table, a recurrent mutation) stopping a placement many SNPs below support. A branch with two or more ancestral SNPs still blocks, however much sits below it, and so does a second single-ancestral-SNP branch on the same path: one odd site is outvoted, two are not. A branch with no call at all neither helps nor blocks, so a branch whose SNPs no table places does not stop the walk. SNPs in the pseudoautosomal regions are left out: the analysis set masks them on chrY.

Which positive branch is the terminal one. A positive branch with no blocking branch above it is a candidate. Candidates rank by their support: derived minus ancestral calls on the path from the root. For two candidates on different lineages, that compares the support of each from their common branch down, because the shared part counts the same for both. Ties go to the deeper branch, then to more derived calls on the path, then to the smaller name. On one lineage the deeper candidate never has less support, so the walk still goes as deep as the calls allow. Across lineages, an outvoted branch no longer wins on depth alone: a clean branch with 3 derived calls below the common branch beats a deeper one with 2 derived calls and its outvoted ancestral one (support 1). The candidate on another lineage with the most support from its common branch with the terminal one, above 0, is named in `Conflict`. That support is counted from the common branch down, so a lineage that splits off low gets no credit for the calls it shares, and the step line names it too. Two lineages with support are a sign of a wrong call, a mismapped stretch or an error in the tree; compare the two supports before you trust the placement.

`Ancestral_in_terminal` lists the SNPs of the terminal branch that read ancestral. The sample then shares only part of that branch's block of SNPs, which is a sign of a sub-branch the tree does not have yet; it can also be a read error, so look at the depth before reading much into it. In `_y_tree.txt`, `Hg` is `NA` when no branch has a derived call.

In Yleaf's `_y_haplogroup.txt`, `Hg` is `NA` when too few markers had reads for a call; the step then prints `Y haplogroup: insufficient markers` with the marker count, and both reports say "insufficient markers".

With a call, both reports add the links of the next section, built from the branch name (for `R-M269*(xR-L23)`, Yleaf's mark for "R-M269 but none of these sub-branches", the links use `R-M269`).

## Explore your branch further

A haplogroup is a branch of the paternal family tree. Other people's tests tell you more about that branch: how old it is, where its men live today, which ancient remains carry it. The routes below are free first. Each was checked on its own pages on 2026-10-25; prices and policies change, so read the service's page before you send anything.

**This pipeline never uploads anything.** An upload shares the file with that service, under its own terms. A BAM of chrY holds the reads of your paternal line only, but it is still your DNA, and through the paternal line it is your father's, brothers' and sons' too. It also keeps the BAM's header (the sample name, the read groups and the command lines of the tools that made it) and, for each read, where its mate aligned. Look at the header with `samtools view -H` before you send it.

| Route | Cost | What you send | What you get | Checked 2026-10-25 |
|---|---|---|---|---|
| FamilyTreeDNA Discover, `https://discover.familytreedna.com/y-dna/<branch>` | free, no account | no file: the branch name goes in the link | the branch's age (when it formed, and its most recent common ancestor), the countries its testers name, ancient DNA samples on it, a migration map, its sub-branches | the page loads its data from `https://discover.familytreedna.com/resources/y-dna/<branch>.json`, which answered for `R-M269` with no login; its search (`/HaplogroupSearch/<name>`) found `R-M269` exactly |
| YFull YTree, `https://www.yfull.com/tree/<branch>/` | free to browse | no file: the branch name goes in the link | the branch on YFull's tree, with its age estimates and the samples placed under it | `/tree/R-M269/` answers 200, a branch YFull does not have answers 404 |
| Decoding Us Navigator, `https://decoding-us.com/download` | free | nothing: the app runs on your computer | haplogroup placement and coverage metrics from your BAM or CRAM; you choose later whether to publish anonymised summaries | the site says it accepts no file uploads and that raw data never leaves your computer; the Navigator is an alpha (v0.1.0-beta.1, 2026-09-03). Its trees are built on T2T-CHM13; we have not run it on a GRCh38 BAM |
| YFull Y analysis, `https://www.yfull.com/order/` | €44 plus VAT, paid after the analysis | a link to a shared BAM or CRAM (YFull does not accept FASTQ; a VCF only for Big Y) | your terminal branch, your novel SNPs, STRs, age estimates, and comparison with YFull's samples; the mtDNA analysis is then free in the same account | YFull's pages do not say whether a BAM with chrY only is enough; we have not tried it |

The two volunteer projects older guides name:

- **James Kane's YDNA Warehouse**: `ydna-warehouse.org` now redirects to Decoding Us (the row above).
- **Alex Williamson's The Big Tree**: `thebigtree.org` had no DNS record on 2026-10-25, so we could not check whether it takes uploads. Treat it as unavailable.

Names differ between trees. Yleaf names branches the YFull way; FamilyTreeDNA keeps its own tree. If Discover shows no page for your branch, type the branch name, or the SNP after the dash, into Discover's search: it offers the nearest branch it has (a search for YFull's `A0-T` returns FamilyTreeDNA's `A-L1098`).

**Private SNPs**, the Y mutations you carry that are on no public tree, are what would place you below your branch, next to the men closest to you. The step lists them in `y_haplogroup/${SAMPLE}_y_private_snps.tsv` (`og y private-snps`), from your own chrY calls minus every change the two public SNP tables name. YFull's analysis reports the same kind of thing as novel SNPs; the file is a starting point for a comparison with other men, not a result in itself.

How that list is built, and what each gate removes:

- one man: the VCF's only sample. A VCF with more than one sample (a joint call) is read only when the step's sample name matches one of its columns (`--sample`), and refused otherwise, so another man's calls never land in the list.
- chrY only, outside the pseudoautosomal regions.
- one REF base and one ALT base, with no second ALT allele at the site: a private SNP is a single-base change, and a multiallelic chrY record is a mapping artefact on a haploid chromosome.
- the genotype carries the ALT allele and no REF allele (`1`, `1/1`, `1|1`). Outside PAR the Y is haploid, so a heterozygous call there is an artefact.
- `FILTER` PASS or `.`, depth 10 or more, genotype quality 20 or more, and 90% or more of the reads carrying ALT (`--min-depth`, `--min-gq`, `--min-fraction`; the depth and fraction defaults are Yleaf's own `-r 10` and `-b 90`).
- inside the accessible part of chrY: the base is marked `P` (passed) in the chrY strict accessibility mask of the 1000 Genomes Project on GRCh38 (release 20160622, about 8.0 Mb of chrY). `setup.sh --yleaf-data` installs it as `ytree/callable_mask_chrY.fasta.gz`, checked against its sha256 in `versions.env`; `--callable-mask` points at another one. Each run of `setup.sh --yleaf-data` checks the installed mask against that sha256. A mask that matches stays, with no download, also when a new tree pin installs the tables again or the tree download fails; one that fails the check, or is missing, is downloaded again on its own (`tests/test_y_private_mask_install.py`). The mask leaves out every base where the short reads of 2,691 samples had too many reads of mapping quality 0, a low average mapping quality, or depth outside a narrow band. That removes the ampliconic, palindromic, X-transposed and heterochromatic sequence, where a call is mostly a difference between near-identical copies rather than a mutation. The sample's own reads cannot find those places: in near-identical copies the aligner gives reads mapping quality 0, the caller drops them, and the calls that remain look normal.
- depth no more than twice the median depth of the sample's own callable chrY calls (`--max-depth-factor`, 0 turns it off). This catches only a collapsed repeat this sample carries more copies of than the mask's samples did, where reads pile up on one place.
- not on a public SNP list: a call is removed when `positions_ftdna_hg38.txt` or `positions_yfull_hg38.txt`, the two tables `setup.sh --yleaf-data` installs, has a row at its position with the same two bases, in either order (the reference can carry the derived base). A row at that position whose bases are not usable (not two different single bases, each A, C, G or T) removes it too. Another base at a known position is a new mutation at that site, so it stays. That subtraction is what makes the list private; without it the file would be mostly the markers the haplogroup call is already built on.

The header lines of the file record the gates, the mask, the sample's median chrY depth and how many records each gate dropped, so a surprising count can be traced. **Each row is unconfirmed until another man is found who carries the same mutation.** Some rows are sequencing or mapping errors. Some are real mutations only this line carries, which no tree can place. The count is not a measure of anything on its own: it rises with coverage and falls with every tree release that names one of them.

The scan needs step 03's VCF (`${SAMPLE}/vcf/${SAMPLE}.vcf.gz`) and the current tree's tables and mask. Without them, the step says so in one line and writes the rest of its output as before.

## How well it places

We measured both placements against Y haplogroups other people published for the same men. The men are 58 men of the 1000 Genomes Project. Poznik et al. 2016 (Nature Genetics 48:593) published their haplogroups in `1000Y.sample_level_haplogroup_calls.ver5.txt` on the 1000 Genomes FTP; we use its column 3, the most derived known SNP the authors saw (`R1b-DF27`). The reads are each man's chrY from his 30x CRAM of the New York Genome Center (NYGC), the same file the public-evidence run reads.

How the men were chosen, fixed before any result was seen. We take the men of that file who have a 30x CRAM, 1,233 of its 1,244. A man counts only when his published SNP sits on one branch of the current tree. Where YFull gives one SNP name to two branches (a recurrent mutation: `M356` is on `C-K98` and on `T-Y31540`), the clade letter of the published name settles it (`C5-M356` is `C-K98`). For each of the file's high-level haplogroups we take its most common published branches: one branch, two for a haplogroup of 10 men or more, three for 50 or more. Branches with the same count are taken in the sort order of their published name. From each branch we take its first man by ID. That gives 59 men in 35 haplogroups, from A0 to T. One branch was left out because its SNP is on two R branches (`R2-L294`). One man was left out because EBI serves a truncated CRAM index for him (HG03846, 4,096 bytes), and his branch had no other man.

Each man ran through step 37's own commands: the two pileups with the step's flags, then Yleaf 4.1.7 with `-rp` on its pileup and `og y tree-place` on the other. `og dev evidence y-bench` then places the published SNP and each tool's call on YFull YTree 14.06.0. It names the call by what it is against the published branch:

| Result | The tool's branch is |
|---|---|
| exact | the published branch |
| below | below it: the same line, deeper |
| ancestor only | above it, on its line |
| sibling | on another line that splits off at the published branch's parent |
| wrong | on a line that splits off higher |
| no call, unknown name | no branch, or a name the tree cannot place |

The published calls come from the ISOGG tree of 2014, which is much shallower than YFull's, so at 30x a correct call is nearly always "below". The benchmark checks the line down to the published branch. It cannot check the levels below it.

| Tool | Men | Exact | Below | Ancestor only | Sibling | Wrong | No call | Unknown name | Published branch or below |
|---|---|---|---|---|---|---|---|---|---|
| Yleaf 4.1.7 | 58 | 1 | 57 | 0 | 0 | 0 | 0 | 0 | 58 of 58 (100.0%) |
| og y tree-place | 58 | 1 | 57 | 0 | 0 | 0 | 0 | 0 | 58 of 58 (100.0%) |

Both tools put every man on his published line, and both go 7 levels below the published branch on average (og 7.4, Yleaf 7.1). The two named the same branch for 44 men. For 13 more both are on one line, and og's branch is deeper for 12 of them (one to three levels). For one man, HG00530, they split into two sister branches below the published one: og's `O-K481` against Yleaf's `O-Y174537`, and og names `O-Y174537` in `Conflict` with support 1 against 4 for its own branch. The published call does not say which of the two is right.

The outvote rule fired for one man (HG01583: one ancestral SNP, `BY37185`, outvoted). That placement is on the published line, 9 levels below it. No man was placed off his published line by the rule or for lack of it, so the data give no reason to change it, and it stays as it is. `Conflict` named a branch for 4 men, each with support 1 against 4 to 31 for the terminal branch. For three R1a men it was the same branch, `R-Y40`.

The public-evidence workflow's `y-benchmark` job runs the step itself again (`scripts/37-y-haplogroup.sh`, with `setup.sh --yleaf-data`) on 6 of these men by default (input `y_samples`), and runs `og dev evidence y-bench` on its calls. Its gate fails when a tool puts fewer than 95% of the men on the published branch or below it. Its red control gives each tool's calls to the man half the list away, and that table must fail the gate. Before it reads any CRAM, it refuses a `y_samples` list in which so many men share a line with that partner that the handed calls could still pass the gate.

<details>
<summary>The 58 men</summary>

| Man | Population | Published (Poznik et al. 2016), its branch on YFull 14.06.0 | Yleaf 4.1.7 | og y tree-place |
|---|---|---|---|---|
| HG01890 | ACB | `A0-L1038`, `A-L1035` | `A-L1035`, exact | `A-L1035`, exact |
| HG02613 | GWD | `A1a-M31`, `A1a` | `A-Y32399`, 3 levels below | `A-Y32399`, 3 levels below |
| HG02588 | GWD | `B-M181`, `B` | `B-Y311987*(xB-ZS382)`, 3 levels below | `B-ZS382`, 4 levels below |
| HG00553 | PUR | `B2-M109`, `B-Y10267` | `B-C137898`, 7 levels below | `B-C137898`, 7 levels below |
| NA18971 | JPT | `C1-P121`, `C-M8` | `C-Y161239`, 3 levels below | `C-Y161239`, 3 levels below |
| HG03917 | BEB | `C3-Z1300`, `C-Z1300` | `C-F3777*(xC-F9966)`, 3 levels below | `C-F3777`, 3 levels below |
| HG03850 | STU | `C5-P92`, `C-P92` | `C-Y33575*(xC-FT358545,C-FTB33293,C-Y82010,C-FTA79799)`, 5 levels below | `C-Y33575`, 5 levels below |
| HG03228 | PJL | `C5-M356`, `C-K98` | `C-FT327336`, 2 levels below | `C-FT327336`, 2 levels below |
| NA18952 | JPT | `D2-M179`, `D-M64.1` | `D-CTS2826`, 5 levels below | `D-CTS2826`, 5 levels below |
| NA18966 | JPT | `D2-M125`, `D-M125` | `D-Y439573`, 2 levels below | `D-Y439573`, 2 levels below |
| HG02461 | GWD | `E1a-L94`, `E-L94` | `E-Z5986*(xE-Z5987,E-Y478064,E-Z5989)`, 2 levels below | `E-FT333060`, 3 levels below |
| HG01094 | PUR | `E1a-CTS2808`, `E-M132` | `E-FT296631`, 6 levels below | `E-FT296631`, 6 levels below |
| HG01556 | CLM | `E1b-P252`, `E-U174` | `E-Z5968*(xE-FTA5825,E-Z5969,E-Z37526)`, 5 levels below | `E-Z5968`, 5 levels below |
| HG01204 | PUR | `E1b-P277`, `E-M4254` | `E-FT295369`, 14 levels below | `E-FT295369`, 14 levels below |
| HG02947 | ESN | `E1b-CTS8030`, `E-Z1704` | `E-Z37814`, 1 level below | `E-Z37814`, 1 level below |
| NA19035 | LWK | `E2-M200`, `E-M200` | `E-Y25288`, 3 levels below | `E-Y25288`, 3 levels below |
| HG02040 | KHV | `F-M213`, `F` | `F4`, 2 levels below | `F4`, 2 levels below |
| NA20858 | GIH | `G1-L833`, `G-M342` | `G-Y262772*(xG-Y99680)`, 9 levels below | `G-Y262772`, 9 levels below |
| NA19670 | MXL | `G2-CTS4803`, `G-CTS4803` | `G-S2784`, 5 levels below | `G-S2784`, 5 levels below |
| HG02224 | IBS | `G2-L30`, `G-L30` | `G-PF3345*(xG-CTS342,G-Y46211,G-FTA56155,G-L497,G-Y507487,G-BY149069,G-Y65645,G-Z39367,G-BY1433)`, 8 levels below | `G-PF3345`, 8 levels below |
| HG03809 | BEB | `H-M370`, `H-M69` | `H-Z34843*(xH-Y168231)`, 9 levels below | `H-Y166956`, 10 levels below |
| HG02684 | PJL | `FH-M2713`, `H` | `H-PH1063`, 5 levels below | `H-M4771`, 6 levels below |
| HG01583 | PJL | `H1-M52`, `H-M52` | `H-Z5885*(xH-Z4695,H-Z12600,H-FT282446,H-Y291749)`, 6 levels below | `H-FT167357`, 9 levels below |
| HG04198 | ITU | `H2-Apt`, `H-Z13973` | `H-Z13969*(xH-Y373420)`, 3 levels below | `H-C116648`, 4 levels below |
| HG00186 | FIN | `I1-L258`, `I-L258` | `I-Y131678`, 7 levels below | `I-Y131678`, 7 levels below |
| HG00117 | GBR | `I1-L118`, `I1` | `I-Y172395`, 7 levels below | `I-Y172395`, 7 levels below |
| HG01527 | IBS | `I2-CTS9183`, `I-CTS616` | `I-Y167082`, 9 levels below | `I-Y167082`, 9 levels below |
| HG01101 | PUR | `I2-L160`, `I-L160` | `I-Y2181`, 9 levels below | `I-Y2181`, 9 levels below |
| HG01069 | PUR | `J1-P58`, `J-P58` | `J-Z2289`, 18 levels below | `J-Z2289`, 18 levels below |
| HG01253 | CLM | `J1-L321`, `J1` | `J-M10096`, 7 levels below | `J-M10096`, 7 levels below |
| HG02236 | IBS | `J2-F4326`, `J-L26` | `J-Y33322*(xJ-Y87966)`, 18 levels below | `J-Y33322`, 18 levels below |
| HG02690 | PJL | `J2-M241`, `J-M241` | `J-FT19136`, 11 levels below | `J-Z8372`, 10 levels below |
| HG01412 | PUR | `J2-L397`, `J-L70` | `J-Z45235*(xJ-Y275293)`, 5 levels below | `J-Z45235`, 5 levels below |
| HG03742 | ITU | `KxLT-M526`, `K2` | `K-Y28301`, 6 levels below | `K-Y28301`, 6 levels below |
| HG02654 | PJL | `L1-M27`, `L-M27` | `L-Y28521*(xL-MF68713)`, 10 levels below | `L-FTA9535`, 12 levels below |
| HG03672 | STU | `L1-L1307`, `L-L1307` | `L-Y35083*(xL-M5513,L-Y106928)`, 6 levels below | `L-Y35083`, 6 levels below |
| HG02138 | KHV | `N-M2283`, `N-Z4762` | `N-Y24355*(xN-Y61250,N-Y65390,N-Y64234)`, 16 levels below | `N-Y24355`, 16 levels below |
| HG00185 | FIN | `N1-Z1940`, `N-Z1940` | `N-Y141790`, 4 levels below | `N-Y141790`, 4 levels below |
| HG00182 | FIN | `N1-CTS10035`, `N-Z1925` | `N-BY142665`, 9 levels below | `N-BY142665`, 9 levels below |
| HG00580 | CHS | `O1-F589`, `O-M119` | `O-FT323420*(xO-MF37402)`, 11 levels below | `O-FT323420`, 11 levels below |
| HG02396 | CDX | `O1-M119`, `O-M119` | `O-Z23482*(xO-Y252345)`, 18 levels below | `O-Y252345`, 19 levels below |
| HG00701 | CHS | `O2-M88`, `O-M111` | `O-CTS3745`, 9 levels below | `O-CTS3745`, 9 levels below |
| HG00442 | CHS | `O2-M95`, `O-M95` | `O-CTS10007`, 2 levels below | `O-CTS10007`, 2 levels below |
| HG00530 | CHS | `O2-F1942`, `O-F1658` | `O-Y174537*(xO-MF81782)`, 5 levels below | `O-K481`, 5 levels below |
| HG00622 | CHS | `O3-M188`, `O-M188` | `O-FT315752`, 5 levels below | `O-FT315752`, 5 levels below |
| HG00524 | CHS | `O3-M1561`, `O-P164` | `O-CTS1174`, 19 levels below | `O-CTS1174`, 19 levels below |
| HG00403 | CHS | `O3-F11`, `O-F11` | `O-ACT3882`, 3 levels below | `O-ACT3882`, 3 levels below |
| HG01124 | CLM | `Q1a-M3`, `Q-M3` | `Q-M848*(xQ-Y126607,Q-Z5906,Q-M925,Q-Z35841,Q-Y592325,Q-CTS2731,Q-MPB139,Q-Y166361,Q-Z19357,Q-Y27993,Q-Z19301,Q-MPB096,Q-Y202441,Q-Z35597,Q-Z19432)`, 1 level below | `Q-M848`, 1 level below |
| HG01944 | PEL | `Q1a-M120`, `Q-M120` | `Q-FT8312*(xQ-Y280718)`, 6 levels below | `Q-FT9713`, 7 levels below |
| HG03652 | PJL | `Q1b-F711`, `Q-L275` | `Q-Y478758`, 5 levels below | `Q-Y478758`, 5 levels below |
| HG02600 | PJL | `R1a-Y7`, `R-Y7` | `R-C10046`, 9 levels below | `R-C10046`, 9 levels below |
| HG02597 | PJL | `R1a-L657`, `R-L657` | `R-C131139`, 13 levels below | `R-Y672703`, 14 levels below |
| HG03636 | PJL | `R1a-Z2123`, `R-Z2123` | `R-YP6547*(xR-FT74091)`, 6 levels below | `R-YP6547`, 6 levels below |
| HG00107 | GBR | `R1b-DF27`, `R-DF27` | `R-A17453*(xR-Y235960,R-Y32852,R-A16155)`, 5 levels below | `R-A17453`, 5 levels below |
| HG00129 | GBR | `R1b-PF6570`, `R-U152` | `R-Y222042`, 3 levels below | `R-Y222042`, 3 levels below |
| HG00126 | GBR | `R1b-M167`, `R-Z201` | `R-BY75410`, 10 levels below | `R-BY75410`, 10 levels below |
| HG01586 | PJL | `R2-L266`, `R2` | `R-Y154917*(xR-Y555693,R-Y155742)`, 18 levels below | `R-FT302409`, 19 levels below |
| HG01133 | CLM | `T-Z709`, `T-Z709` | `T-CTS9882*(xT-Y169564,T-Y15648,T-Y4981,T-Y273415)`, 12 levels below | `T-CTS6901`, 13 levels below |

</details>

## Runtime

Not measured on a full genome; the end test runs the step on the fixture's chrY slices only. An estimate is a few minutes: Yleaf reads the pileup at about 183,000 Y marker positions, and the current tree's pileup is at about 369,000 more. The placement reads the tree and the tables in about a second.

## Notes

- Yleaf 4 places the sample on the YFull tree v14.01 and names the haplogroup by clade and defining marker, for example `R-M269`. Its `QC-score` (0 to 1) is how consistently the markers on the path to that branch agree; the default acceptance threshold is 0.95.
- Short reads cover only the unique parts of the Y; that is enough for a haplogroup at 30x.
- The two results can differ in depth and in name: YFull renames and splits branches between releases, and the two use different position tables. When they disagree on the lineage itself, not just the depth, look at `Conflict` and the per-level table: the lineage with less support (derived minus ancestral calls from the common branch), an outvoted one included, is the weaker one.
- On the e2e fixture (six stretches of chrY, 1.3 Mb in all) the step runs on a copy of HG002 with a male indexcov row, because indexcov reads the fixture's slices as female.
