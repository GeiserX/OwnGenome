# Step 0: Reference Data Setup

One-time downloads required before running the pipeline. Each heading below gives the size of its download; [Hardware and storage requirements](hardware-requirements.md#shared-reference-data-one-time) adds them up (about 73 GB for a default run, 248 GB with the optional annotation databases).

> **Estimated time:** 1-3 hours depending on internet speed. The two VEP caches (26 GB for step 13, 24 GB for step 17) are the largest downloads of a default run.

A plain `./scripts/setup.sh "$GENOME_DIR"` installs both VEP caches and the PCGR/CPSR bundle after AnnotSV.
These clinical downloads add about 57 GB to transfer and 65 GB on disk. Allow 1 to 3 hours.
Setup checks each install separately: VEP 116 needs 56 GB peak, the bundle 15 GB, and VEP 115 needs 52 GB.
Each check adds a 5 GB margin and credits the disk space occupied by a partial download.
Stale extraction files are removed before the check. Retained PCGR archives are checked against the pinned SHA-256 before reuse. An item that does not fit is skipped with a warning. Setup continues.

```bash
./scripts/setup.sh --no-clinical-data "$GENOME_DIR"  # skip all three clinical downloads
./scripts/setup.sh --clinical-data "$GENOME_DIR"     # install or retry only these three
```

An interrupted transfer resumes. Existing caches and bundles are kept.
`--clinical-data` exits 1 when any of the three is still missing after it runs.
Without these data, `run-all.sh` skips the affected steps and records the reason.
A bare `nextflow run` still needs its tools and data parameters. Setup does not change that selection.

See [Licences for downloads](licences.md) for the terms that apply to individual bundles and tools.

## The downloads folder

Setup keeps every raw download, the archive or file as the server sent it, in one folder: `<genome_dir>/downloads`, or the folder `OG_DOWNLOADS` names (a relative one is taken from the working directory). A rerun, a `--refresh`, or a second install with the same `OG_DOWNLOADS` then downloads nothing it already has. That matters when a server is slow or limits how much you can download. `assets/data-sources.tsv` lists every source: its path, URL, checksum, size and licence.

- **Before you start**, `bin/og setup plan --for setup` (or `--all` for `--everything`, or `--for <flag>`) lists each source the chosen options need, whether it is already there, its size, the total to download, the raw files kept, and the free disk. `setup.sh` runs the same check first and stops before any download when the raw downloads it still needs do not fit in the downloads folder's disk. A download whose install is already done is not counted, so a rerun after `bin/og setup prune` needs no extra space. `SETUP_DISK_CHECK=false` skips that check. A source kept as a folder of many files (the archaic genomes, Locityper's region) counts as there only once its install is done; until then the files it holds are taken off its size, and a lock, a `.part` or another source's file in it counts for nothing. An install that setup records with a stamp (the fine ancestry and MoChA panels, the archaic tables, Locityper, gnomAD-SV and MANE for AnnotSV) counts as done only while the stamp names the pins setup compares (the release, its checksum, the image), so a reinstall that a new pin forces is counted.
- **A file already there** that matches its pinned checksum is used with no download. When the checksum is a list on the server (VEP's `CHECKSUMS`, an md5 list), that small list is read first. The checksum that matched is saved beside the file as `<name>.checksum`, and when the list cannot be read (the server is down), a kept file that matches that saved checksum is used. A file with the wrong checksum is refused with a message that names the source and its row of the table, and it is never overwritten: delete it or replace it yourself.
- **An older install**: a raw file or `.part` that an older setup left at its old place (in `vep_cache/`, `pcgr_data/`, `reference/archaic/<version>/download/` and the like) is moved into the downloads folder before the fetch, which checks it there. When a copy is in both places, the one that passes its check stays and the other goes; when no check says which one is right, both stay and setup says so.
- **og does the downloads**: `setup.sh` builds `bin/og` before its first download and stops when it cannot (`setup.sh --og` shows why). A step script that calls `fetch` needs `bin/og` too.
- **A failed unpack** of a download that matched its pinned checksum keeps the download: the cause is the disk or the tool, not the file. Only a download with no checksum is deleted, so the next run fetches it again.
- **Bring your own file**: copy it into the downloads folder under its name in the table (`bin/og setup source <id> path` prints the full path). Setup checks it like a download.
- **A cut transfer** resumes from `<name>.part`. A failed try is tried again `FETCH_TRIES` times (default 3), `FETCH_WAIT` seconds apart (default 5, then 10, ...). A transfer under 1 KB/s for `FETCH_STALL` seconds (default 120) counts as failed. curl, or wget when there is no curl, comes from `PATH`. Two setups that share the folder take turns on the same file: the second one waits, then uses the file the first one checked.
- **The cost**: the kept raw files take about 69 GB for a plain setup and about 344 GB for `--everything` (the 150 GB of archaic genomes are most of it). `bin/og setup plan` prints the figure for your options.
- **To keep less**: `setup.sh --no-keep-downloads` (or `SETUP_KEEP_DOWNLOADS=false`) deletes each raw download once it is installed. `bin/og setup prune` deletes the kept downloads later, and nothing else. It refuses a file whose folder resolves outside the downloads folder (a subfolder moved away and linked back).

## Install data for an everything run

```bash
./scripts/setup.sh --everything "$GENOME_DIR"
```

This runs base setup, then installs the data listed in [the everything manifest](https://github.com/GeiserX/OwnGenome/blob/main/assets/everything.tsv) in order. It includes the extra STR loci, even if `STR_EXTRA_LOCI` is unset. Setup prints the restricted-use terms before downloading. Review the terms and the [hardware requirements](hardware-requirements.md#shared-reference-data-one-time) before starting. MEGAnE alone needs about 50 GB of memory for the full reference.

A failed installer goes on the missing list with a retry command, and setup continues through the manifest. An incomplete setup exits with status 3. After retrying, complete downloads and indexes stay in place. User inputs such as HPO terms and family relationships are still supplied when you run the pipeline.

Cyrius uses the PolyForm Strict 1.0.0 licence for non-commercial use. The following PGS Catalog scores are freely available to the academic community for research use. Review each score's terms before using it.

- [PGS000014](https://www.pgscatalog.org/score/PGS000014/)
- [PGS000016](https://www.pgscatalog.org/score/PGS000016/)
- [PGS000017](https://www.pgscatalog.org/score/PGS000017/)
- [PGS000027](https://www.pgscatalog.org/score/PGS000027/)

## GRCh38 Reference Genome

The foundation for everything. All tools need this. The pipeline uses NCBI's **GRCh38 no-ALT analysis set**: chr1-22, X, Y and M, the unplaced and unlocalized scaffolds and the EBV genome, 195 sequences, with UCSC names (`chr1`, `chrM`). `setup.sh` downloads it, checks it and stores it as `reference/GRCh38_no_alt_analysis_set.fasta`:

```bash
./scripts/setup.sh ${GENOME_DIR}
```

By hand, the same steps:

```bash
export GENOME_DIR=/path/to/your/data
mkdir -p ${GENOME_DIR}/reference
cd ${GENOME_DIR}/reference
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids

# Download (~0.9 GB compressed, ~3.2 GB unpacked) and its index
wget -c ${NCBI}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
wget -c ${NCBI}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
wget ${NCBI}/md5checksums.txt

# Verify both against the md5 NCBI lists for them; each line must say OK
grep -E ' \./GCA_000001405\.15_GRCh38_no_alt_analysis_set\.fna\.(gz|fai)$' md5checksums.txt | md5sum -c -

# Unpack under the name the scripts read
gzip -dc GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > GRCh38_no_alt_analysis_set.fasta
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai GRCh38_no_alt_analysis_set.fasta.fai
wc -l GRCh38_no_alt_analysis_set.fasta.fai   # 195
```

The md5 checks the download, and `validate-setup.sh` checks the sequence dictionary `setup.sh` makes from the unpacked FASTA (`GRCh38_no_alt_analysis_set.dict`, below) on every run. It does not read the FASTA itself. It reads the MD5 of each sequence's bases (`M5`) from the `.dict`, and compares the 195 names, lengths and MD5s with the values pinned in the script for NCBI's file and for [the masked copy](#the-masked-reference-optional). A FASTA whose `.dict` shows the same names and lengths but one other base fails with `REF_FASTA has the names and lengths of NCBI's GRCh38 no-ALT analysis set but other bases`. A masked FASTA without `_maskedGRC` in its name fails too, as does a `.dict` whose names or lengths are not the `.fai`'s, and a `.dict` older than the FASTA gets a warning: a FASTA changed after its `.dict` was made gets that warning, not a failure, because `setup.sh` keeps an existing `.dict`. A reference with other names or lengths, such as the decoy set below, has nothing pinned: an `[INFO]` line says so.

### Why the no-ALT analysis set

GRCh38 also has ALT contigs: second copies of regions that vary a lot between people, such as the MHC (the HLA genes), KIR and the CYP2D6 locus. A reference with them (the Broad `hg38` FASTA, 3,366 sequences with ALT, HLA and decoy contigs) only helps an aligner that is run ALT-aware, and none of the aligners here is. With them, a read that matches the primary copy and an ALT copy equally well gets mapping quality 0, and callers ignore it: depth thins at exactly those loci, and depth-based callers can report a deletion that is not there. The [measured loss](realignment.md#how-much-depth-alt-contigs-cost) is in the realignment page.

The no-ALT analysis set has each of those regions once. It is the reference DeepVariant's case studies use, the one the GIAB v4.2.1 benchmark is defined on, and the one PharmCAT normalises against. chr1-22, X, Y and M are the same sequence as in the Broad file, so ClinVar, the VEP cache and the score files stay valid.

Two other references were considered and are not the default:

- **The no-ALT set plus the hs38d1 decoys** (`GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz` in the same NCBI directory, 2,580 sequences). The decoys catch reads from sequence missing in the primary assembly, which removes some false calls for callers less robust than DeepVariant. None of the tools here expects them, and its 2,385 extra contigs are more that CNVpytor and Delly have to be kept away from. It works through `REF_FASTA` (below).
- **The same set with GRCh38's false duplications masked.** `setup.sh --masked-reference` builds it from the default reference; it is optional, and the next section says what it changes.

### The masked reference (optional)

GRCh38 holds some stretches twice where people have one copy: false duplications, errors of the assembly (Aganezov et al., Science 2022, PMID 35357935). The largest are five stretches of chr21's short arm, about 0.77 Mb, which repeat parts of the long arm: KCNE1, CBS, U2AF1 and CRYAA among them. A read from the real gene fits both copies equally well, gets mapping quality 0 and is skipped by the callers, so on the default reference those genes are not read at all. Behera et al. (FixItFelix, Genome Biol 2023, PMID 36810122) showed that masking the false copy brings their reads and calls back, without changing any coordinate.

```bash
./scripts/setup.sh --masked-reference ${GENOME_DIR}
```

It writes `reference/GRCh38_no_alt_analysis_set_maskedGRC.fasta` beside the default reference, with its `.fai` and `.dict`: the same 195 sequences, the same lengths, and the 12 regions of the GRC's list of false duplications (`GRC_EXCLUSIONS_URL` in `versions.env`, GIAB's copy, checked against its sha256) set to N, about 1.2 Mb. It needs the default reference installed first and takes a few minutes and about 3.2 GB more. The FASTA it writes must have the sha256 pinned in `versions.env` (`GRC_MASKED_REF_SHA256`); any other file is deleted. Every one of its 195 sequences has the same bases (upper and lower case aside) as in GIAB's `GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2` file (compared once, 2026-11-04), and chr21 has the same masked stretches as FixItFelix's own reference, which also adds decoy contigs. This one adds none, so the sequence list stays the one `validate-setup.sh` expects.

Because the coordinates do not change, ClinVar, the VEP cache, the score files and the other downloads stay valid. The reads do not: a BAM aligned to the default reference keeps its MAPQ 0 reads, so a sample has to be aligned again to gain anything, as [Realigning: the masked reference](realignment.md#the-masked-reference) shows.

The CI case `tests/e2e/fixitfelix-chr21-1.sh` measures it on chr21 with GIAB's HG002 and the CMRG v1.00 benchmark (Wagner et al., Nat Biotechnol 2022, PMID 35132260), see [the realignment page](realignment.md#the-masked-reference) for the numbers.

GSTT1 lies on an ALT contig in GRCh38, so pypgx (step 32) cannot call it on the default reference; the step says so and calls the other genes. No other step calls it, so the report lists it with the genes no step reads. It is one of the ALT-only genes of `assets/callability_genes.tsv`, the genes short reads often cannot place, which step 16b measures and the report lists ([gene callability](16b-mosdepth.md#gene-callability)); that table is committed, so it needs no download. Steps 32 (pypgx) and 35 (Parascopy and the HBA1/HBA2 count) refuse a BAM aligned to a reference with ALT or HLA contigs, with a `REFUSED` line, even when `ALLOW_ALT_REFERENCE=true` lets the other steps run on it: their depth would split between a gene and its ALT copy. Run one by one, each of the two scripts stops with exit 1 and the other steps are not touched. When the check cannot run (the BAM header cannot be read, or `og cn check-bam` fails), each stops with exit 2 instead of going on with a BAM nobody checked. On the Nextflow path (`run-all.sh`) the refusal fails the run: the tasks already running finish, no new task starts, and `run-all.sh` writes the reports of what finished. The RHD count of step 41 and the LPA KIV-2 count of step 46 do not check the BAM for ALT contigs yet.

Every sample aligned to another reference has to be aligned again. [Realigning after a reference change](realignment.md) says what to redo, what to keep and how to check a BAM.

### The reference path on every page

The step scripts read the reference from `${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta`. The commands on the other pages of these docs write that path as `${REF_FASTA}`, relative to `GENOME_DIR`: on the host it is `${GENOME_DIR}/${REF_FASTA}`, and inside a container that mounts `GENOME_DIR` at `/genome` it is `/genome/${REF_FASTA}`. Set it once in the shell where you paste those commands:

```bash
export REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta
```

The scripts read `REF_FASTA` too, as a path relative to `GENOME_DIR` or as an absolute path inside it, so the same variable points both at another reference. `setup.sh` then also needs `REF_FASTA_URL` (a `.gz` URL is unpacked) and `REF_FASTA_MD5` (an md5, the URL of a checksum file that lists the download, or empty), and `REF_FAI_MD5` the same way for the `.fai` published beside it (empty builds the index with samtools). For the decoy variant:

```bash
export REF_FASTA=reference/GRCh38_no_alt_plus_hs38d1_analysis_set.fasta
export REF_FASTA_URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
./scripts/setup.sh ${GENOME_DIR}   # checks both files against NCBI's md5checksums.txt
```

`validate-setup.sh` fails when the reference has ALT or HLA contigs. To keep such a reference on purpose, set `ALLOW_ALT_REFERENCE=true` as well; the check then warns instead.

### Why GRCh38?

This pipeline uses **GRCh38** (also called hg38) exclusively. It's the current standard genome build with:
- Corrected mitochondrial sequence (rCRS). Step 20 also calls the chrM control region on chrM rotated by 8000 bp; it writes that rotated copy from this reference in each run (`og setup chrm-shift`), so setup installs no extra file for it ([step 20](20-mtoolbox.md#the-shifted-reference-pass-for-the-control-region))
- Better representation of centromeres and telomeres
- `chr` prefix naming (chr1, chr2, ..., chrX, chrY, chrM)

If your data is on **GRCh37/hg19**, extract FASTQ from BAM and re-align. See [vendor-guide.md](vendor-guide.md#genome-build-grch37-hg19-vs-grch38-hg38).

## ClinVar Database

Updated monthly by NCBI. Contains known pathogenic/benign variant classifications. `setup.sh` downloads it, checks it against NCBI's published md5, builds the two files the steps read from it and records its release date:

| File | What it is |
|---|---|
| `clinvar/clinvar.vcf.gz` (+ `.tbi`) | NCBI's file as published, chromosomes named `1, 2 ... MT` |
| `clinvar/clinvar_chr.vcf.gz` (+ `.tbi`) | the same records with `chr1, chr2 ... chrM`, the names the BAMs and VCFs use |
| `clinvar/clinvar_pathogenic_chr.vcf.gz` (+ `.tbi`) | the Pathogenic and Likely_pathogenic records, read by step 6 |
| `clinvar/RELEASE` | the release date, from the file's `##fileDate` line |

To move to the current release, run:

```bash
./scripts/setup.sh --refresh clinvar ${GENOME_DIR}
```

It downloads the new file under `clinvar/.refresh/`, checks the md5, builds both derived files there, and only then replaces all six files and `RELEASE`. Step 6's normalised copy (`clinvar_pathogenic_chr.norm.vcf.gz`) is removed, so step 6 rebuilds it from the new release. A failed download or build leaves the installed release as it was. If moving the new files in fails part way, the refresh exits with an error and leaves no `RELEASE`, so `validate-setup.sh` reports the date as unknown until a refresh completes. Do not re-download with a plain `wget`: it writes `clinvar.vcf.gz.1` beside the old file, and the derived files stay built from the old one.

`validate-setup.sh` prints the release date and warns when it is more than 35 days old.

> **Tip:** Refresh ClinVar monthly for the latest classifications. ClinVar adds ~1000 new pathogenic variants per month. Rerun step 6 afterwards.

## Small pinned data files

`setup.sh` also installs these small files, each from a fixed commit or release and checked before it is stored (`./scripts/setup.sh --sample-qc-data ${GENOME_DIR}` installs only the last two). The chromosome bands are the exception: the repository ships them, and setup copies them without a download. A failed download does not stop setup; the next run tries again, and `validate-setup.sh` lists what is missing.

| File | Source | Used by | Without it |
|---|---|---|---|
| `reference/delly_human.hg38.excl.tsv` | Delly's GRCh38 exclude map at a pinned commit (sha256 checked) | step 19 (`delly sr -x`) | Delly runs without it: slower, and with calls in centromeres, telomeres and the unplaced scaffolds |
| `reference/cytoBand.hg38.txt` | UCSC's GRCh38 chromosome bands, chr1-22, X and Y, copied from `assets/reference/cytoBand.hg38.txt` (its README says how it was made). `validate-setup.sh` fails an installed copy that differs from it | step 10 (`telomerehunter2 -b`), Nextflow `--cytoband`; the genome map reads the asset itself | Step 10 installs it from `assets/reference` before it runs; Nextflow passes the asset itself |
| `hla/IPD-IMGT-HLA_<release>/hla.dat` | IPD-IMGT/HLA release `HLA_DB_RELEASE` (3.66.0) from the IMGTHLA repository (md5 checked) | step 8 | step 8 is skipped |
| `reference/gencode.v50.basic.genes.gtf` | the gene lines of GENCODE 50's basic annotation (md5 checked) | step 8 (gene positions for T1K) | step 8 is skipped |
| `reference/gencode.v50.basic.exons.bed` | the exon intervals of the same GENCODE 50 basic annotation, written from the same download as the gene lines (md5 checked) | steps 6 and 17 and Nextflow `--gencode_exons` with `clinvar` (cDNA read check); step 39 and `--gencode_exons` with `--pedigree` (deletions over an exon) | steps 6 and 17 and `CLINVAR_CDNA_CHECK` write every hit as `not_checked`; step 39 does not use deletions |
| `gnomad_mito/gnomad.genomes.v3.1.sites.chrM.vcf.bgz` | gnomAD v3.1's mitochondrial sites, about 5 MB (`GNOMAD_MITO_URL`, sha256 `GNOMAD_MITO_SHA256` in `versions.env`, checked by `setup.sh` and `validate-setup.sh`) | step 20 (gnomAD AF_hom, AF_het and max_hl beside each chrM call) | step 20 writes its annotation table with the gnomAD columns empty, and `validate-setup.sh` warns (it fails on a file whose sha256 is not the pinned one) |
| `roh_af/roh_af_gnomad_v4.1.tsv.gz` and its `.tbi` | step 11's allele-frequency panel: the 7,439,999 biallelic PASS SNVs of chr1 to chr22 with a global AF of 0.01 to 0.99 in gnomAD v4.1 genomes (CC0), with the global AF and five continental group AFs, about 208 MB. A release asset of this repository built by `og dev build roh-af` ([below](#roh-allele-frequency-panel)); `ROH_AF_PANEL_URL`, sha256 `ROH_AF_PANEL_SHA256` in `versions.env`, checked by `setup.sh`, which then indexes it with `tabix` | step 11 and Nextflow `--roh_af_panel` (`bcftools roh --AF-file`) | step 11 gives every site an allele frequency of 0.4 (`--AF-dflt 0.4`) and says so; `validate-setup.sh` does not check it |
| `reference/GRCh38_AllTandemRepeats.bed.gz` and its unpacked copy `reference/GRCh38_AllTandemRepeats.bed` | every GRCh38 tandem repeat of GIAB's genome stratifications v3.6 (TRF and RepeatMasker simple repeats, 5 bp added on each side), about 9 MB, NIST data free to reuse (`TANDEM_REPEATS_URL`, sha256 `TANDEM_REPEATS_SHA256` in `versions.env`, checked by `setup.sh`). Sniffles2 reads a plain BED only, so setup keeps the unpacked copy beside the download | step 04c (`sniffles --tandem-repeats`, long reads only) | step 04c prints a note and runs Sniffles2 without it: SVs inside repeats get less precise breakpoints and merging |
| `reference/somalier/sites.hg38.vcf.gz` | somalier's GRCh38 sites, 17,766 SNPs (sha256 checked, also when the file is already there: one with another sha256 gets a warning and is downloaded again) | step 33 (`somalier extract`), Nextflow `--somalier_sites` | step 33 stops and says to install it |
| `reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.{UD,mu,bed}` | VerifyBamID2's 1000 Genomes panel of 100,000 markers, from the VerifyBamID v2.0.3 release (sha256 checked, also when a file is already there, as for the sites) | step 33 (`verifybamid2`), Nextflow `--verifybamid2_panel` | step 33 stops and says to install it |

`setup.sh` also writes the reference's sequence dictionary (`GRCh38_no_alt_analysis_set.dict`), which GATK, Picard and `chip-to-vcf.sh` need. It holds the MD5 of each sequence, which `validate-setup.sh` compares with the pinned values ([above](#grch38-reference-genome)).

## ROH allele-frequency panel

Step 11 (`bcftools roh`) needs each site's population allele frequency: see [Allele frequencies](11-roh-analysis.md#allele-frequencies). gnomAD publishes them in its genome sites VCFs, 563 GB for chr1 to chr22, far more than step 11 needs. `og dev build roh-af` keeps the biallelic PASS SNVs of chr1 to chr22 whose global AF is 0.01 to 0.99, with the global AF and the AFs of the African (`afr`), admixed American (`amr`), East Asian (`eas`), non-Finnish European (`nfe`) and South Asian (`sas`) groups, copied as gnomAD writes them (`.` where a group has no allele number). A site is biallelic when gnomAD saw one alternate allele at its position (`n_alt_alleles=1`). The result is one tab-separated file, `CHROM POS REF,ALT AF AF_afr AF_amr AF_eas AF_nfe AF_sas`, in BGZF; its first four columns are the `--AF-file` format of `bcftools roh`, which reads column 4 only.

`setup.sh` downloads the published file (`ROH_AF_PANEL_URL`), checks its sha256 and writes its index with `tabix -s1 -b2 -e2` from `GATK_IMAGE`. To build it again from gnomAD, byte for byte, with `og` built from this checkout: each VCF streams through `og` and is not kept, and each stream's MD5 is compared with the one gnomAD's bucket lists for the file (a stream cut short makes `og` stop with "the gzip input ends inside a member").

```bash
#!/usr/bin/env bash
# A failed download, a failed og run or an MD5 that differs stops the build
# before the panel is joined.
set -euo pipefail
mkdir -p rows
for c in $(seq 1 22); do
  f=gnomad.genomes.v4.1.sites.chr${c}.vcf.bgz
  want=$(curl -sS --fail "https://storage.googleapis.com/storage/v1/b/gcp-public-data--gnomad/o/release%2F4.1%2Fvcf%2Fgenomes%2F${f}?fields=md5Hash" \
    | sed -n 's/.*"md5Hash": *"\([^"]*\)".*/\1/p')
  curl -sS --fail "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/${f}" | tee >(openssl dgst -md5 -binary | base64 > "rows/chr${c}.md5") \
    | bin/og dev build roh-af rows --vcf - --out "rows/chr${c}.tsv"
  sleep 2
  if [ -z "$want" ] || [ "$(cat "rows/chr${c}.md5")" != "$want" ]; then
    echo "chr${c}: MD5 differs from gnomAD's, build it again" >&2
    exit 1
  fi
done
bin/og dev build roh-af panel --rows $(printf 'rows/chr%s.tsv ' $(seq 1 22)) --out roh_af_gnomad_v4.1.tsv.gz
sha256sum roh_af_gnomad_v4.1.tsv.gz   # ROH_AF_PANEL_SHA256
```

Run it as a script, not pasted into a shell: `exit 1` stops it. `og dev build roh-af rows` also fails on an input with no `#CHROM` line, such as an empty download.

`og dev build roh-af panel ... --check` compares a panel with what the rows give and exits 1 when they differ. The published panel was built this way, three chromosomes at a time, in 1 h 45 min. That is 563 GB read at about 90 MB/s. It holds 7,439,999 sites and is 207,888,195 bytes. It is published as a release asset of the tag `data-roh-af-gnomad-v4.1`, the file `ROH_AF_PANEL_URL` names.

## AnnotSV Annotations (~5 GB)

AnnotSV's annotation bundle includes third-party data with their own terms. Review the sources in the [AnnotSV manual](https://www.lbgi.fr/AnnotSV/Documentation/README.AnnotSV_latest.pdf) before downloading.

Required for step 5. `setup.sh` downloads them (about 5 GB, about 20 GB unpacked) into `${GENOME_DIR}/annotsv_annotations/`; its source server is slow, so this can take an hour or two, and an interrupted download resumes when `setup.sh` runs again. `./scripts/setup.sh --annotsv-data ${GENOME_DIR}` installs only these (the bundle when it is missing) and the two files below.

| File | Source | Used by | Without it |
|---|---|---|---|
| `annotsv_annotations/Annotations_Human/SVincludedInFt/BenignSV/GRCh38/gnomad-sv.installed` | gnomAD-SV `GNOMAD_SV_VERSION` (4.1) sites, sha256 checked; its common SVs are merged into AnnotSV's benign set, the download stays in the downloads folder, and this file records the release | step 5 (population frequencies for the benign criteria) | step 5 and the Nextflow `ANNOTSV` stop and point at `--annotsv-data` |
| `annotsv_annotations/mane_select_refseq.txt` | the MANE Select RefSeq ids of NCBI's MANE `MANE_RELEASE` (1.5) summary, sha256 checked | step 5 (`-txFile`: the MANE Select transcript of each gene) | step 5 and the Nextflow `ANNOTSV` stop and point at `--annotsv-data` |

[Step 5](05-annotsv.md#gnomad-sv-frequencies-and-mane-select-transcripts) says what each one changes.

## VEP Cache (~26 GB)

Ensembl Variant Effect Predictor annotation database. Required for step 13 and the pipeline's VEP. Setup installs it by default. Without it, a run skips VEP, vcfanno, the clinical filter and slivar. To install only this cache:

```bash
./scripts/setup.sh --vep-cache ${GENOME_DIR}
```

It downloads the release-116 cache, checks it against Ensembl's `CHECKSUMS` file, unpacks it into a temporary folder and moves it to `${GENOME_DIR}/vep_cache/homo_sapiens/116_GRCh38/` only when it is complete. The tarball (~26 GB) and the unpacked cache (~30 GB) are both on disk; the tarball stays in the downloads folder (`--no-keep-downloads` deletes it once the cache is in place). An interrupted download resumes when you run the command again. It runs on the host (curl or wget, then tar), not in a container. Most of its time is the download: 26 GB at 100 Mbit/s is about 35 minutes, plus the unpacking.

By hand, without the checksum check and the move from a temporary folder:

```bash
mkdir -p ${GENOME_DIR}/vep_cache/tmp
cd ${GENOME_DIR}/vep_cache/tmp

# Download (manual wget is more reliable than VEP INSTALL.pl)
# The -c flag enables resume if the download is interrupted
wget -c https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz

# Extract to parent directory (~30 GB extracted)
cd ${GENOME_DIR}/vep_cache
tar xzf tmp/homo_sapiens_vep_116_GRCh38.tar.gz
# Creates: ${GENOME_DIR}/vep_cache/homo_sapiens/116_GRCh38/

# Optional: delete the tarball to save 26 GB
# rm tmp/homo_sapiens_vep_116_GRCh38.tar.gz
```

> **Warning:** The VEP `INSTALL.pl` script downloads to a temporary directory that may lack write permissions inside Docker. Always download manually with `wget -c`. See [lessons-learned.md](lessons-learned.md) for details.

## PCGR/CPSR Ref Data Bundle (~7 GB)

Required for step 17 (CPSR cancer predisposition screening). Includes ClinVar, gnomAD, CancerMine, and other databases. PCGR 2.x uses a separate, smaller ref data bundle — VEP cache is mounted independently.

> **Important:** PCGR 2.3.2 bundles VEP 115, which requires a **release-115** cache — different from the release-116 cache used by step 13 above. See the next section for the VEP 115 download. The bundle date (`PCGR_DATA_BUNDLE`) and the cache release (`PCGR_VEP_CACHE_RELEASE`) are in `versions.env` and move with `PCGR_IMAGE`.

PCGR/CPSR reference data includes third-party databases with their own terms. The [PCGR bundle report](https://sigven.github.io/pcgrdb/grch38/) lists the sources and licences of the newest bundle release. Setup installs release 20260620 (`PCGR_DATA_BUNDLE`), so that list can differ from the installed bundle. Setup prints the [PCGR resource list](https://sigven.github.io/pcgr/articles/annotation_resources.html) before it downloads the bundle.

Setup installs it by default. To add the clinical data to an existing setup:

```bash
./scripts/setup.sh --clinical-data "$GENOME_DIR"
```

The installer checks the archive against `PCGR_DATA_SHA256` in `versions.env`, then extracts it into a temporary folder.
The server publishes no checksum, so the size (`PCGR_DATA_BYTES`) and the digest are measured from a full download. `PCGR_DATA_MEASURED_FOR` names the image tag and bundle date they belong to, and `scripts/ci/check-coupled-versions.sh` fails when `PCGR_IMAGE` or `PCGR_DATA_BUNDLE` moves without it.
Only a complete tree is moved to `pcgr_data/20260620/data/`. The archive stays in the downloads folder (`--no-keep-downloads` deletes it after installation).

## VEP 115 Cache for CPSR (~24 GB)

PCGR 2.3.2 (step 17) bundles VEP 115 internally, which needs the **release-115** cache. This is separate from the release-116 cache used by step 13. Both coexist in the same `vep_cache/` directory under different subdirectories (`116_GRCh38/` and `115_GRCh38/`). A `113_GRCh38/` cache left from PCGR 2.2.5 is no longer read and can be deleted.

Setup installs this cache by default, after the PCGR bundle. `--clinical-data` retries it using Ensembl's `CHECKSUMS` file and a temporary extraction folder.

> If you only run step 13 (VEP annotation) and skip step 17 (CPSR), you only need the release-116 cache. If you only run step 17, you only need release-115.

## T1K HLA Reference (Optional)

Only needed for step 8 (HLA typing). Step 8 builds its T1K index itself, the first time it runs (12-26 seconds in the measured runs), from the `hla.dat` and GENCODE gene lines `setup.sh` installs (see "Small pinned data files" above). The index goes to `t1k_idx/t1k-<T1K version>_imgt-<release>_gencode-<release>/`, so a new T1K image or a new database release builds a new index instead of reusing the old one. Step 8 writes the release it typed against to `hla_t1k/database_release.txt`.

To type against another IPD-IMGT/HLA release, set `HLA_DB_RELEASE` (for example `HLA_DB_RELEASE=3.64.0`) for both `setup.sh` and step 8.

The coordinate file takes each HLA gene's GRCh38 position from the GENCODE annotation (`t1k-build.pl -g`, as T1K's README describes). Built from the FASTA or its `.fai` instead, every gene gets `-1 -1` coordinates and T1K extracts no reads; step 8 stops when a typed gene has no coordinates.

> **Note:** HLA typing from WGS is challenging. For clinical HLA typing, dedicated lab assays are more reliable. See [lessons-learned.md](lessons-learned.md#t1k-coordinate-file-with-wrong-values).

## CNVpytor GC/Mask Resources (Optional, for Step 18)

The CNVpytor 1.3.2 container needs its GC and mask files from the upstream release. Setup verifies each SHA-256 and writes the files under `reference/cnvpytor/`.

```bash
./scripts/setup.sh --cnvpytor-data "$GENOME_DIR"
```

## PyPGx Data Bundle (Optional, for Step 32)

Setup clones the matching bundle release from upstream into `reference/pypgx-bundle/` and checks its commit against the pin. Upstream publishes no licence terms for the data bundle. See [Step 32](32-pypgx.md#pypgx-bundle-required-one-time-download) for details.

```bash
./scripts/setup.sh --pypgx-data "$GENOME_DIR"
```

## Somatic Calling Resources (Optional, for Steps 29 and 49)

Install the GATK Best Practices hg38 resources for Mutect2. Setup checks the SHA-256 of each VCF and its Tabix index before publishing it to `somatic/`. The gnomAD AF-only VCF, 1000 Genomes panel of normals, and ExAC common sites take about 3 GiB together. See [Step 29](29-mutect2-somatic.md#recommended-resources-reduce-false-positives) for what each does.

```bash
./scripts/setup.sh --somatic-data "$GENOME_DIR"
```

## MEGAnE k-mer Set (Optional, for Step 50)

Only needed for step 50, the mobile element insertions ([50-mei.md](50-mei.md)). MEGAnE drops reads that map to many places by a list of every 32-mer seen more than once in the reference, which it builds from the reference itself:

```bash
./scripts/setup.sh --megane-kmers ${GENOME_DIR}
# Creates: ${GENOME_DIR}/reference/megane/GRCh38_no_alt_analysis_set.mk and .mi
```

It pulls `MEGANE_IMAGE` (about 6.4 GB) and runs MEGAnE's `build_kmerset` on `REF_FASTA`. The builder holds two 8-byte k-mers per reference base in memory, so GRCh38 needs about 50 GB of RAM (MEGAnE's authors give 48 GB and about 10 minutes); `setup.sh` checks the memory Docker has first and stops with a message when it is short. On a smaller machine, build it where the memory is and copy `reference/megane/` over. A set newer than the reference is kept; a new reference needs a new set.

## Locityper Database (Optional, for Step 53)

Only needed for step 53, the second opinion at HLA-DRB1, HLA-DQB1, KIR, CYP2D6 and FCGR ([53-locityper.md](53-locityper.md)):

```bash
./scripts/setup.sh --locityper-data ${GENOME_DIR}
# Creates: ${GENOME_DIR}/reference/locityper-hprc-v2.1-mc-grch38/ (counts.jf, db/, loci.bed, DONE)
```

It pulls `LOCITYPER_IMAGE`, counts the 25-mers of `REF_FASTA` with jellyfish, reads the HPRC release 2 pangenome VCF (`LOCITYPER_PANEL_URL`) over https only around the loci of `assets/locityper_loci.bed`, through its `.tbi` checked against `LOCITYPER_PANEL_TBI_SHA256`, and builds the database with `locityper target` and `locityper augment`. A new reference, image, panel or locus list builds again.

## Annotation Databases (Optional, for Steps 30-31)

Deep pathogenicity scoring and variant prioritization. These databases power vcfanno (step 30) and slivar (step 31). All are optional — step 30 gracefully skips any missing track.

Install all four score sets into `annotations/` with about 175 GiB of free disk for the downloads, plus temporary space for the REVEL sort.

```bash
./scripts/setup.sh --annotation-scores "$GENOME_DIR"
```

Setup checks for enough space to finish the missing CADD and SpliceAI score downloads, plus a 5 GiB margin. It resumes interrupted transfers and skips masked SpliceAI SNV or indel scores when the corresponding raw file and its Tabix index are already present.

It checks CADD scores against the publisher's MD5 and validates the masked SpliceAI scores against the pinned byte sizes, SHA-256 of the first and last MiB, and gzip integrity. Both publishers' Tabix indexes, the REVEL ZIP, and the AlphaMissense download use full SHA-256 pins in `versions.env`.


### CADD v1.7 Pre-scored (SNVs + Indels) — ~83 GB

The most widely used deleteriousness metric. Scores all possible SNVs genome-wide plus gnomAD indels.

> **License:** Free for non-commercial/academic use. Commercial use requires a license from the University of Washington. The pipeline does not redistribute CADD data — users download directly from the source.

> **Note:** CADD TSVs use chromosome names without `chr` prefix (1, 2, 3...). vcfanno handles this via column mapping in the TOML config — no manual renaming needed.

### SpliceAI Pre-scored — ~91 GB

Deep learning splice-site variant prediction. Catches pathogenic intronic variants that VEP's rule-based splice prediction misses.

The pipeline uses the **masked** score files. SpliceAI's authors recommend the masked scores for variant interpretation and the raw scores for alternative-splicing research. In the masked files, a gain at an annotated splice site and a loss at a site that is not an annotated splice site are set to 0, because neither changes how a variant is interpreted.

> **License:** The precomputed SpliceAI scores are free for academic and not-for-profit use only; any other use requires a commercial license from Illumina. They are not Apache 2.0. The pipeline does not redistribute them.

Step 30 and `validate-setup.sh` find the masked files under these names. They also accept the raw files (`spliceai_scores.raw.*`) if you already have them; when both sets are present, step 30 uses the raw ones.

> **Alternative source:** Illumina BaseSpace at `https://basespace.illumina.com/s/otSPW8hnhaZR` (requires free account).

### REVEL v1.3 — ~0.6 GB

Ensemble pathogenicity scoring for missense variants, combining 13 individual tools. Recommended by ClinGen for missense variant classification.

The REVEL authors make the scores freely available for non-commercial use. The [Zenodo record](https://doi.org/10.5281/zenodo.7072866) lists ODbL 1.0. For other uses, contact the maintainers through the [REVEL download page](https://sites.google.com/site/revelgenomics/downloads) before downloading.

Setup projects the release onto GRCh38 with `og setup revel ZIP TSV`, then compresses and indexes `revel_grch38.tsv.gz`.

> **Thresholds:** ClinGen's calibration (Pejaver et al. 2022) gives REVEL >= 0.644 as PP3_Supporting, >= 0.773 as PP3_Moderate and >= 0.932 as PP3_Strong for missense variants. There is no PP3 Very Strong level. See [interpreting-results.md](interpreting-results.md#revel-rare-exome-variant-ensemble-learner) for the BP4 (benign) levels.

### AlphaMissense — ~613 MB

DeepMind's protein-structure-informed missense classifier. Complements REVEL with structural context. Covers all ~71 million possible human missense variants.

> **License:** CC BY-NC-SA 4.0 (non-commercial, share-alike). The pipeline does not redistribute these scores.

Setup indexes the published BGZF file as `AlphaMissense_hg38.tsv.gz.tbi`.

> **Thresholds:** am_pathogenicity < 0.34 = likely benign, > 0.564 = likely pathogenic. These are AlphaMissense's own class boundaries, not ACMG evidence levels.

### gnomAD v4.1 Constraint Metrics — ~91 MB

Per-gene loss-of-function intolerance scores. Essential for interpreting novel variants in constrained genes.

```bash
cd ${GENOME_DIR}/annotations

# Gene-level constraint table
wget -c -O gnomad_v4.1_constraint.tsv \
  https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv
```

> **Key metrics:** LOEUF (oe_lof_upper) < 0.35 = constrained for LoF; pLI >= 0.9 = intolerant to LoF; missense Z > 3.09 = constrained for missense.

## Ensembl Gene Annotation (Optional, for Step 49, ~110 MB)

Only needed for step 49, the clonal haematopoiesis (CHIP) screen. It reads each CHIP gene's MANE Select transcript from Ensembl's GRCh38 gene annotation:

```bash
./scripts/setup.sh --chip-data ${GENOME_DIR}
# Creates: ${GENOME_DIR}/reference/ensembl/Homo_sapiens.GRCh38.<release>.chr.gff3.gz
```

The release is `ENSEMBL_GFF_RELEASE` in `versions.env`, and the download is checked against `ENSEMBL_GFF_SHA256`; a file already there with another sha256 is downloaded again. The Nextflow pipeline takes the same file as `--chip_gff`. See [Step 49](49-chip.md).

## Exomiser Data (Optional, for Step 48, ~37 GB)

Exomiser's hg38 and phenotype bundles include third-party data with their own terms. Review Exomiser's [data source notes](https://exomiser.readthedocs.io/en/15.1.0/advanced_analysis.html#pathogenicitysources) before installing.

Only needed for step 48, the symptom-driven ranking that runs when you give HPO terms. `./scripts/setup.sh --exomiser-data ${GENOME_DIR}` downloads Exomiser's data release `EXOMISER_DATA_VERSION` (hg38 and phenotype archives, 37 GB, each checked against its sha256 in `versions.env`), unpacks it into `${GENOME_DIR}/reference/exomiser-<release>/` (58 GB; about 72 GB free while it runs), deletes the archives and pulls `EXOMISER_IMAGE`. An interrupted download resumes when the command runs again. See [Step 48](48-exomiser.md).

## AADR Ancient Panel (Optional, for Step 44, ~1.05 GB kept, 9.6 GiB at the peak)

UCSC allows free non-commercial use of its hg19ToHg38 chain by independent researchers and nonprofit organisations. Other uses require a commercial licence. See the [UCSC chain download terms](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/).

Only needed for step 44, the placement among ancient genomes. `./scripts/setup.sh --ancient-panel ${GENOME_DIR}` downloads the AADR release `AADR_VERSION` (1240K, the 7.1 GB `.geno` and three small files, each checked against its md5 in `versions.env`) and UCSC's hg19ToHg38 chain, then builds the panel in `${GENOME_DIR}/reference/ancient/aadr_<release>/` with convertf and two smartpca runs. On a 4-CPU GitHub runner it took 61 minutes, peaked at 9.6 GiB of extra disk and kept 1.05 GiB; its smartpca container gets 12 GB, and the setup's containers held 7.6 GiB of anonymous memory at most (any setup container, convertf included). See [Step 44](44-ancient.md#runtime).

## Docker Images — Pre-Pull All

`setup.sh` pulls every image in `versions.env` except the lines marked `# optional`, and `scripts/setup.sh --pull-only` pulls the same list without the rest of setup. A failed pull is tried three times (`FETCH_TRIES` and `FETCH_WAIT` change that), and Docker's own message is printed: `toomanyrequests` is Docker Hub's rate limit, `manifest unknown` a tag that does not exist. `--parascopy-data`, `--yleaf-data` and `--megane-kmers` pull their step's optional image along with its data. Each image pinned by digest (`name:tag@sha256:...`) also gets a local tag, `localhost/pgp-pinned/<name>:sha256-<digest>`: Docker gives an image pulled by digest no tag, and a dangling-image cleanup such as `docker image prune` would delete it ([Troubleshooting](troubleshooting.md#images-that-were-pulled-are-missing-after-a-docker-cleanup)). Tasks still run the digest reference. To pull every image in advance, the optional ones included, from the repository root in bash (then run `scripts/setup.sh --pull-only` once to tag the digest-pinned ones):

```bash
source versions.env
for var in $(grep -oE '^[A-Z0-9_]+_IMAGE=' versions.env | tr -d '='); do
  docker pull "${!var}"
done
```

## Alternative Callers (Optional, for Benchmarking)

Only needed if you plan to run alternative variant callers. See [benchmarking.md](benchmarking.md).

The alternative aligners and callers and hap.py (`benchmark-variants.sh`) are the lines marked `# optional` in `versions.env`. `setup.sh` skips them, and each one is pulled the first time its script runs. [Image versions](versions.md) lists them under "Alternative aligners and callers (optional)".

### GATK Sequence Dictionary

GATK HaplotypeCaller requires a `.dict` file alongside the reference FASTA. `setup.sh` creates it; to create it by hand:

```bash
source versions.env   # from the repository root
REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta   # see "The reference path on every page" above
docker run --rm --user root \
  -v ${GENOME_DIR}:/genome \
  "${GATK_IMAGE}" \
  gatk CreateSequenceDictionary \
    -R /genome/${REF_FASTA}
```

### BWA-MEM2 Index

BWA-MEM2 requires its own index files (different from minimap2's `.mmi`). Building them needs a lot of memory: the bwa-mem2 README states 28 GB per Gbp of reference, which is **about 90 GB of RAM** for the 3.2 Gbp GRCh38 FASTA. The finished index is about 10 GB on disk, and aligning with it needs about 10 GB of RAM, so only the one-time build is the problem.

The practical route is to build the index once on a machine (or a rented cloud instance) with at least 96 GB of RAM, then copy the five index files next to the FASTA on your own machine. With less memory the build is killed (exit code 137).

```bash
REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta   # see "The reference path on every page" above
source versions.env   # from the repository root: BWAMEM2_IMAGE, the image step 02a uses
docker run --rm --user root \
  --cpus 8 --memory 96g \
  -v ${GENOME_DIR}:/genome \
  "${BWAMEM2_IMAGE}" \
  bwa-mem2 index "/genome/${REF_FASTA}"
# Creates: .0123, .amb, .ann, .bwt.2bit.64, .pac alongside the FASTA
```

GRIDSS (step 4b) needs the classic BWA index instead (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`); the two are not interchangeable. Run `./scripts/setup.sh --bwa-index "$GENOME_DIR"` to build it. See [GRIDSS prerequisites](04b-gridss.md#prerequisites).

### GIAB Truth Set (for hap.py Benchmarking)

Download a GIAB truth set for benchmarking variant callers. **HG002** (Ashkenazi Jewish male) is preferred because its truth set covers more difficult genomic regions. HG001 (NA12878) is an alternative used in the [quick test](quick-test.md).

**Important:** Truth set benchmarking only works when the query VCF comes from the **same biological sample** as the truth set. You must sequence HG002 (or HG001) DNA, not your own sample. See [benchmarking.md](benchmarking.md) for details.

```bash
mkdir -p ${GENOME_DIR}/giab
cd ${GENOME_DIR}/giab

# HG002 truth set (recommended, GRCh38 v4.2.1)
# Pinned to the NISTv4.2.1 directory: GIAB's latest/ directory moves to each new
# release (it now holds v5.0q), so latest/ URLs for v4.2.1 files return 404.
wget -c https://giab.s3.amazonaws.com/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget -c https://giab.s3.amazonaws.com/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
wget -c https://giab.s3.amazonaws.com/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed

# Alternative: HG001/NA12878 (used in quick-test.md)
# wget -c https://giab.s3.amazonaws.com/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
# wget -c https://giab.s3.amazonaws.com/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
# wget -c https://giab.s3.amazonaws.com/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed
```

**Total Docker image size:** ~10-15 GB (compressed, after layer deduplication). Alternative tools add ~3-5 GB.

## Disk Space Summary

[Hardware and storage requirements](hardware-requirements.md#shared-reference-data-one-time) has the one table of every download above, with the totals.

> **Tip:** If disk space is tight, use `setup.sh --no-clinical-data` to skip both VEP caches and the PCGR bundle initially. The core pipeline (steps 2-3-6-7) only needs the reference FASTA and ClinVar (~3.5 GB unpacked).

## Verifying Your Setup

After all downloads, verify everything is in place:

```bash
REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta   # see "The reference path on every page" above
echo "Checking reference setup..."
[ -f "${GENOME_DIR}/${REF_FASTA}" ] && echo "  GRCh38 FASTA: OK" || echo "  GRCh38 FASTA: MISSING"
[ -f "${GENOME_DIR}/${REF_FASTA}.fai" ] && echo "  FASTA index: OK" || echo "  FASTA index: MISSING"
[ -f "${GENOME_DIR}/clinvar/clinvar_chr.vcf.gz" ] && echo "  ClinVar (chr): OK" || echo "  ClinVar: MISSING"
[ -d "${GENOME_DIR}/vep_cache/homo_sapiens/116_GRCh38" ] && echo "  VEP cache: OK" || echo "  VEP cache: MISSING"
[ -d "${GENOME_DIR}/pcgr_data/20260620/data" ] && echo "  PCGR data: OK" || echo "  PCGR data: MISSING"
[ -d "${GENOME_DIR}/vep_cache/homo_sapiens/115_GRCh38" ] && echo "  VEP 115 cache (CPSR): OK" || echo "  VEP 115 cache (CPSR): MISSING"

echo "Annotation databases (optional, for steps 30-31):"
for DB_PAIR in \
  "whole_genome_SNVs.tsv.gz:CADD SNVs" \
  "gnomad.genomes.r4.0.indel.tsv.gz:CADD indels" \
  "spliceai_scores.masked.snv.hg38.vcf.gz:SpliceAI SNVs" \
  "spliceai_scores.masked.indel.hg38.vcf.gz:SpliceAI indels" \
  "revel_grch38.tsv.gz:REVEL" \
  "AlphaMissense_hg38.tsv.gz:AlphaMissense" \
  "gnomad_v4.1_constraint.tsv:gnomAD constraint"; do
  DB_FILE="${DB_PAIR%%:*}"
  DB_NAME="${DB_PAIR#*:}"
  if [ -f "${GENOME_DIR}/annotations/${DB_FILE}" ]; then
    # Check for .tbi index (not needed for plain TSV files)
    if [[ "$DB_FILE" == *.vcf.gz ]] || [[ "$DB_FILE" == *.tsv.gz && "$DB_FILE" != *constraint* ]]; then
      [ -f "${GENOME_DIR}/annotations/${DB_FILE}.tbi" ] && echo "  ${DB_NAME}: OK (+index)" || echo "  ${DB_NAME}: OK (WARNING: .tbi index missing)"
    else
      echo "  ${DB_NAME}: OK"
    fi
  else
    echo "  ${DB_NAME}: not downloaded"
  fi
done
echo "Done."
```
