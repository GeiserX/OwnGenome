#!/usr/bin/env bash
# setup.sh — One-stop setup: download references, pull Docker images, validate
# Usage: ./scripts/setup.sh <genome_dir>
#        ./scripts/setup.sh --everything <genome_dir>
#        ./scripts/setup.sh --no-clinical-data <genome_dir>
#        ./scripts/setup.sh --clinical-data <genome_dir>
#        ./scripts/setup.sh --pull-only      pull the Docker images and exit
#        ./scripts/setup.sh --og             put OwnGenome's og binaries in bin/.og
#                                            (Linux, and macOS for a live-site host):
#                                            built from rust/ in OG_BUILDER, or on
#                                            a release tag downloaded and checked
#                                            (rust/og.sh install; docs/rust.md)
#        ./scripts/setup.sh --refresh clinvar <genome_dir>
#                                            replace ClinVar with NCBI's current
#                                            release and rebuild the files made from it
#        ./scripts/setup.sh --cyrius <genome_dir>
#        ./scripts/setup.sh --parascopy-data <genome_dir>
#        ./scripts/setup.sh --kir-data <genome_dir>
#        ./scripts/setup.sh --ancestry-panel <genome_dir>
#        ./scripts/setup.sh --ancestry-fine-panel <genome_dir>
#        ./scripts/setup.sh --yleaf-data <genome_dir>
#        ./scripts/setup.sh --consumer-sites <genome_dir>
#        ./scripts/setup.sh --ancient-panel <genome_dir>
#        ./scripts/setup.sh --viral-data <genome_dir>
#        ./scripts/setup.sh --masked-reference <genome_dir>
#        ./scripts/setup.sh --exomiser-data <genome_dir>
#        ./scripts/setup.sh --strchive-data <genome_dir>
#        ./scripts/setup.sh --chip-data <genome_dir>
#        ./scripts/setup.sh --megane-kmers <genome_dir>
#        ./scripts/setup.sh --mocha-data <genome_dir>
#        ./scripts/setup.sh --locityper-data <genome_dir>
#                                            what an opt-in step needs, and nothing
#                                            else: Cyrius (step 21), Parascopy's
#                                            homology table and models plus its
#                                            image (step 35), IPD-KIR (step 08
#                                            with KIR=true), pgsc_calc's ancestry
#                                            panel (step 26, and percentiles in
#                                            step 25), the fine-grained ancestry
#                                            panel plus its images (step 26b),
#                                            Yleaf's marker tables plus its image
#                                            (step 37), the site table of the
#                                            raw-data export (step 40), the
#                                            AADR panel of step 44, the
#                                            viral genomes of step 47, the
#                                            reference with GRCh38's false
#                                            duplications masked, Ensembl's
#                                            gene annotation (step 49), Exomiser's
#                                            hg38 and phenotype data plus its image
#                                            (step 48), MEGAnE's k-mer set of the
#                                            reference plus its image (step 50), the
#                                            STRchive catalog of steps 09 and 9b
#                                            (also part of a full setup), MoChA's
#                                            phasing panel plus its images (step 51),
#                                            Locityper's database of the loci of
#                                            assets/locityper_loci.bed plus its
#                                            image (step 53)
#        ./scripts/setup.sh --archaic <genome_dir>
#                                            the archaic genomes and 1000 Genomes
#                                            frequencies of step 43, as one table
#        ./scripts/setup.sh --annotsv-data <genome_dir>
#                                            AnnotSV's data for step 5 (the bundle
#                                            when missing, gnomAD-SV's common SVs
#                                            and the MANE Select list) and exit
#        ./scripts/setup.sh --vep-cache <genome_dir>
#                                            the VEP cache of VEP_CACHE_RELEASE
#                                            (VEP, vcfanno, clinical filter and
#                                            slivar need it), before the first run
#
# This script downloads everything needed to run the pipeline:
#   1. GRCh38 reference genome + index: NCBI's GRCh38 no-ALT analysis set
#      (~0.9 GB download, ~3.2 GB unpacked)
#   2. ClinVar database (~200 MB), with its release date in clinvar/RELEASE
#   3. All Docker images (~10-15 GB)
#   4. The reference's sequence dictionary and small pinned data files: Delly's
#      exclude map, GRCh38 chromosome bands, an IPD-IMGT/HLA release, the
#      GENCODE gene coordinates T1K needs (~350 MB) and GENCODE's exons for
#      the cDNA read check of steps 06 and 17, both from one ~80 MB download,
#      and somalier's sites and VerifyBamID2's marker panel for step 33 (~10 MB),
#      and gnomAD v3.1's mitochondrial sites for step 20 (~5 MB), the
#      STRchive repeat catalog of steps 09 and 9b (~0.5 MB), and GIAB's
#      tandem-repeat map for Sniffles2 in step 04c (~9 MB, ~35 MB unpacked)
#   5. AnnotSV annotation data for step 5 (~5.3 GB download, ~20 GB unpacked),
#      plus gnomAD-SV's common SVs (~1.7 GB download, kept as a few MB) and the
#      MANE Select transcript list
#   6. Step 25: the PGS Catalog scores of assets/pgs_scores.tsv (~400 MB), a
#      checkout of pgsc_calc and the Nextflow plugin it needs, so step 25 and
#      the PRS process run without the network. The four scores free for
#      academic research use only (licence academic-research) are left out:
#      --prs-academic installs them (PRS_ACADEMIC=true scores them)
#
#   7. VEP 116, PCGR/CPSR data and VEP 115 (~57 GB download, ~65 GB unpacked).
#      --no-clinical-data skips them; --clinical-data installs only these three.
#
# Exit status of a full setup: 0 when everything is installed ("Setup
# complete!"); 3 when the core is ready but an item of 4 to 7 failed or was
# skipped ("Setup finished: core ready, N items missing", each listed with the
# command that retries it); 1 when the reference, ClinVar or an image failed.
#
# Every download goes through fetch (scripts/lib/common.sh): it is written to
# <file>.part, checked, and only then renamed, so a file that exists is whole.
# A failed image pull is tried FETCH_TRIES times too, FETCH_WAIT apart, and
# Docker's own message (a rate limit, a missing tag) is printed as it comes.
#
# The reference can be replaced by another GRCh38 build: set REF_FASTA (where
# it is stored), REF_FASTA_URL (a .gz URL is unpacked) and REF_FASTA_MD5 (an
# md5, the URL of a checksum file that lists the download, or empty for no
# check), and REF_FAI_MD5 the same way for the .fai published beside it (empty
# builds the index with samtools). docs/00-reference-setup.md lists the builds
# that work and docs/realignment.md what a change of reference means for
# existing samples.
set -euo pipefail

# Keep an explicit checksum override for fixture tests and verified mirrors.
PCGR_SHA_OVERRIDE=${PCGR_DATA_SHA256:-}
GNOMAD_SV_SHA_OVERRIDE=${GNOMAD_SV_SHA256:-}
MANE_SHA_OVERRIDE=${MANE_SUMMARY_SHA256:-}

# These notices also print before the corresponding downloads.
ANNOTSV_TERMS="AnnotSV's annotation bundle includes third-party data with their own terms."
EXOMISER_TERMS="Exomiser's hg38 and phenotype bundles include third-party data with their own terms."
CYRIUS_TERMS="Cyrius: PolyForm Strict 1.0.0, non-commercial use only."
PRS_ACADEMIC_TERMS=(
  'PGS000014: freely available to the academic community for research use.'
  'PGS000016: freely available to the academic community for research use.'
  'PGS000017: freely available to the academic community for research use.'
  'PGS000027: freely available to the academic community for research use.'
)
UCSC_CHAIN_TERMS="UCSC allows free non-commercial use of its hg19ToHg38 chain by independent researchers and nonprofit organisations. Other uses require a commercial licence."

# The usage text below names data releases, so read them before it prints.
# shellcheck source=../versions.env
. "$(dirname "$0")/../versions.env"

EVERYTHING_SETUP=false
PULL_ONLY=false
REFRESH=""
SAMPLE_QC_ONLY=false
OPT_IN=""
CLINICAL_DATA=true
# --no-keep-downloads, anywhere among the arguments, is SETUP_KEEP_DOWNLOADS=false:
# each raw download is deleted once installed instead of kept in the
# downloads folder (OG_DOWNLOADS, default <genome_dir>/downloads).
SETUP_KEEP_DOWNLOADS=${SETUP_KEEP_DOWNLOADS:-true}
ARGS=()
for a in "$@"; do
  if [ "$a" = --no-keep-downloads ]; then SETUP_KEEP_DOWNLOADS=false; else ARGS+=("$a"); fi
done
set -- ${ARGS[@]+"${ARGS[@]}"}
case "$SETUP_KEEP_DOWNLOADS" in
  true|false) export SETUP_KEEP_DOWNLOADS ;;
  *) echo "ERROR: SETUP_KEEP_DOWNLOADS must be true or false, got '${SETUP_KEEP_DOWNLOADS}'" >&2; exit 2 ;;
esac
case "${1:-}" in
  --everything)
    EVERYTHING_SETUP=true
    shift ;;
  --og)
    exec "$(dirname "$0")/../rust/og.sh" install ;;
  --cyrius|--prs-academic|--parascopy-data|--kir-data|--ancestry-panel|--ancestry-fine-panel|--yleaf-data|--consumer-sites|--ancient-panel|--viral-data|--masked-reference|--exomiser-data|--chip-data|--strchive-data|--megane-kmers|--mocha-data|--locityper-data|--vep-cache|--annotsv-data|--annotation-scores|--somatic-data|--bwa-index|--cnvpytor-data|--pypgx-data|--archaic)
    OPT_IN=${1#--}
    shift ;;
  --clinical-data)
    OPT_IN=clinical-data
    shift ;;
  --no-clinical-data)
    CLINICAL_DATA=false
    shift ;;
  --pull-only)
    PULL_ONLY=true
    shift ;;
  --sample-qc-data)
    SAMPLE_QC_ONLY=true
    shift ;;
  --refresh)
    REFRESH=${2:-}
    if [ "$REFRESH" != clinvar ]; then
      echo "Usage: $0 --refresh clinvar <genome_dir>" >&2
      echo "  (clinvar is the one database this flag refreshes)" >&2
      exit 1
    fi
    shift 2 ;;
  --*) echo "ERROR: unknown setup flag: $1" >&2; exit 2 ;;
esac

GENOME_DIR=${1:-${GENOME_DIR:-""}}
if [ -z "$GENOME_DIR" ] && ! $PULL_ONLY; then
  echo "Usage: $0 <genome_dir>"
  echo "       $0 --everything <genome_dir>"
  echo "       $0 --no-clinical-data <genome_dir>"
  echo "       $0 --clinical-data <genome_dir>"
  echo "       $0 --pull-only"
  echo "       $0 --og"
  echo "       $0 --refresh clinvar <genome_dir>"
  echo "       $0 --sample-qc-data <genome_dir>"
  echo "       $0 --cyrius | --prs-academic | --parascopy-data | --kir-data | --ancestry-panel | --ancestry-fine-panel | --yleaf-data | --consumer-sites | --ancient-panel | --viral-data | --masked-reference | --exomiser-data | --chip-data | --strchive-data | --megane-kmers | --mocha-data | --locityper-data | --archaic | --vep-cache | --annotsv-data <genome_dir>"
  echo ""
  echo "  <genome_dir>  Where to store reference data and sample outputs."
  echo "                Needs at least 500 GB free space per sample."
  echo "  --everything  Run full setup and every installer in assets/everything.tsv."
  echo "  --annotation-scores  Install CADD, masked SpliceAI, REVEL and AlphaMissense (~175 GiB)."
  echo "  --somatic-data       Install GATK somatic resources for steps 29 and 49."
  echo "  --bwa-index          Build the classic BWA index for GRIDSS."
  echo "  --cnvpytor-data      Install CNVpytor reference files (step 18)."
  echo "  --pypgx-data         Install the upstream PyPGx bundle (step 32)."
  echo "  --no-clinical-data  Skip both VEP caches and the PCGR/CPSR bundle."
  echo "  --no-keep-downloads  Delete each raw download once it is installed (SETUP_KEEP_DOWNLOADS=false)."
  echo "                Without it every raw download stays in the downloads folder, OG_DOWNLOADS"
  echo "                (default <genome_dir>/downloads), so a rerun downloads nothing again; this keeps"
  echo "                about 300 GB more for --everything. bin/og setup plan lists what a setup downloads,"
  echo "                and bin/og setup prune deletes the kept downloads. Goes with any other flag."
  echo "  --clinical-data     Install only those three downloads for an existing setup."
  echo "  --pull-only   Pull the Docker images listed in versions.env and exit."
  echo "  --og          Put OwnGenome's og binaries in bin/.og, for Linux and for a macOS live-site host"
  echo "                (built from rust/, or a release's) and exit."
  echo "  --refresh clinvar"
  echo "                Download NCBI's current ClinVar, rebuild the files made from it"
  echo "                and record its release date in <genome_dir>/clinvar/RELEASE."
  echo "  --sample-qc-data"
  echo "                Install only somalier's sites and VerifyBamID2's panel (step 33) and exit."
  echo "  --cyrius      Install Cyrius for the opt-in step 21 (PyPI, hash-locked) and exit."
  echo "                Cyrius is under the PolyForm Strict licence: non-commercial use only."
  echo "  --prs-academic"
  echo "                Install the four PGS Catalog scores of step 25 that are free for academic research use only"
  echo "                (PGS000014, PGS000016, PGS000017, PGS000027; PRS_ACADEMIC=true scores them) and exit."
  echo "  --parascopy-data"
  echo "                Install Parascopy's GRCh38 homology table and models (step 35, ~50 MB), pull its image and exit."
  echo "  --kir-data    Install the IPD-KIR ${KIR_DB_RELEASE:-} database (KIR=true in step 08, ~40 MB) and exit."
  echo "  --ancestry-panel"
  echo "                Install pgsc_calc's ancestry reference panel ${PGSC_PANEL:-} (~7 GB download): step 26, and"
  echo "                percentiles instead of raw scores in step 25; then exit."
  echo "  --ancestry-fine-panel"
  echo "                Install the reference of step 26b, gnomAD's HGDP + 1000 Genomes panel (~19 GB is read, a few GB"
  echo "                stored, hours on a home connection), and its images; then exit."
  echo "  --yleaf-data  Install Yleaf ${YLEAF_DATA_VERSION:-}'s Y marker tables and trees and the current YFull tree (step 37, about 210 MB to download, 220 MB on disk), pull its image and exit."
  echo "  --consumer-sites"
  echo "                Build the site table of the raw-data export (step 40): the 23andMe v5 site list joined with"
  echo "                dbSNP ${CONSUMER_DBSNP_BUILD:-}'s GRCh38 VCF (~30 GB download, kept in the downloads folder), with UCSC's"
  echo "                hg19ToHg38 chain for each site's strand; then exit."
  echo "  --archaic     Build the archaic table of step 43: the Altai and Vindija Neanderthals and the Altai Denisovan"
  echo "                (MPI-EVA, ~165 GB, kept in the downloads folder) with 1000 Genomes frequencies, lifted to GRCh38"
  echo "                (hours on a home connection; --no-keep-downloads keeps 13 GB at most at a time); then exit."
  echo "  --ancient-panel"
  echo "                Build the panel of step 44 from the Allen Ancient DNA Resource ${AADR_VERSION:-} (1240K, ~7.3 GB download,"
  echo "                kept in the downloads folder; needs about 12 GB of free RAM) and pull EIGENSOFT; then exit."
  echo "  --viral-data  Install the human viral genomes and phiX of step 47 (NCBI RefSeq, VIRAL_REF_ACCESSIONS, ~2 MB) and exit."
  echo "  --masked-reference"
  echo "                Build the default reference with the false duplications of GRCh38 masked (KCNE1, CBS, U2AF1, CRYAA"
  echo "                and others become readable; ~3.2 GB beside it, a few minutes) and exit. Samples must be realigned to it."
  echo "  --exomiser-data"
  echo "                Install Exomiser ${EXOMISER_DATA_VERSION:-}'s hg38 and phenotype data (step 48: about 37 GB to download, 58 GB"
  echo "                unpacked, about 72 GB free while it unpacks), pull its image and exit."
  echo "  --chip-data   Install Ensembl ${ENSEMBL_GFF_RELEASE:-}'s GRCh38 gene annotation (step 49, the CHIP genes' MANE Select"
  echo "                transcripts; ~110 MB) and exit."
  echo "  --strchive-data"
  echo "                Install only the STRchive ${STRCHIVE_VERSION:-} repeat catalog of steps 09 and 9b (~0.5 MB) and exit."
  echo "  --megane-kmers"
  echo "                Build MEGAnE's k-mer set of the reference for step 50 (Alu, L1 and SVA insertions; ~1.8 GB,"
  echo "                a few minutes and about 5 GB of RAM, with og), pull MEGAnE's image (~6.4 GB) and exit."
  echo "  --mocha-data  Install MoChA's phasing panel for step 51 (1000 Genomes high coverage, ${MOCHA_PANEL_VERSION:-}: ~30 GB read,"
  echo "                one chromosome at a time) and the genetic map, pull its images and exit."
  echo "  --locityper-data"
  echo "                Build Locityper's database for step 53: the k-mer counts of the reference (several GB of RAM"
  echo "                once, a few GB on disk), the HPRC ${LOCITYPER_PANEL_VERSION:-} haplotypes of the loci in"
  echo "                assets/locityper_loci.bed (read over https, tens of MB), their pairwise alignments; pull its image and exit."
  echo "  --annotsv-data"
  echo "                Install AnnotSV's data for step 5 and exit: the bundle when it is missing (~5.3 GB), then"
  echo "                gnomAD-SV ${GNOMAD_SV_VERSION:-}'s common SVs (~1.7 GB download) and MANE ${MANE_RELEASE:-}'s Select transcripts."
  echo "  --vep-cache   Install the VEP ${VEP_CACHE_RELEASE:-} cache (step 13 and the Nextflow VEP: ~26 GB download, ~30 GB"
  echo "                unpacked, both on disk during the install) and exit. Without it a run skips VEP, vcfanno,"
  echo "                the clinical filter and slivar."
  echo ""
  echo "Download and use terms (see docs/licences.md for the sources):"
  echo "  ${ANNOTSV_TERMS}"
  echo "  PCGR/CPSR reference data includes third-party databases with their own terms."
  echo "  ${EXOMISER_TERMS}"
  echo "  The REVEL authors make the scores freely available for non-commercial use."
  echo "  ${UCSC_CHAIN_TERMS}"
  echo "  Illumina's DRAGEN ORA Decompression Software requires acceptance of its end-user licence."
  echo "  The PacBio Software License Agreement permits Paraphase use only with data generated on a PacBio instrument or provided by PacBio."
  echo ""
  echo "Example:"
  echo "  ./scripts/setup.sh /data/genomics"
  echo "  ./scripts/setup.sh ~/genome_data"
  exit 1
fi

export GENOME_DIR

# Image tags, data versions and the helpers (run_in, fetch, pipeline_images)
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
if [ -n "$PCGR_SHA_OVERRIDE" ]; then PCGR_DATA_SHA256=$PCGR_SHA_OVERRIDE; fi
if [ -n "$GNOMAD_SV_SHA_OVERRIDE" ]; then GNOMAD_SV_SHA256=$GNOMAD_SV_SHA_OVERRIDE; fi
if [ -n "$MANE_SHA_OVERRIDE" ]; then MANE_SUMMARY_SHA256=$MANE_SHA_OVERRIDE; fi

# Every download goes through og (og setup fetch), so og is built before the
# first one. Then, when the raw downloads are kept, the sources this run
# downloads must fit in the downloads folder: og setup plan --check refuses
# up front instead of filling the disk halfway (--no-keep-downloads, or
# OG_DOWNLOADS on a larger disk, are the ways out it names). SETUP_DISK_CHECK=false
# skips that check: the test harnesses set it, as their sources are small
# fixtures, not the sizes the table records.
if ! $PULL_ONLY; then
  # OG names another og (the tests set it); otherwise bin/og is built when it
  # is stale or does not run here. The downloads need only this machine's og,
  # so a missing Linux pair is need_og's matter, for the installers that run
  # og in a container.
  if [ -n "${OG:-}" ]; then
    [ -x "$OG" ] || { echo "ERROR: OG=${OG} is not an executable og" >&2; exit 1; }
  elif { ! "${PGP_ROOT}/rust/og.sh" check > /dev/null 2>&1 || ! "${PGP_ROOT}/bin/og" version > /dev/null 2>&1; } \
      && ! "${PGP_ROOT}/rust/og.sh" install; then
    echo "ERROR: could not build og, which does every download of setup: run ${0} --og and look at its error." >&2
    exit 1
  fi
  if [ "$SETUP_KEEP_DOWNLOADS" = true ] && [ "${SETUP_DISK_CHECK:-true}" != false ] && [ -z "$REFRESH" ] && [ "${SETUP_EVERYTHING_SCOPE:-}" != row ]; then
    if $EVERYTHING_SETUP; then PLAN_FOR=(--all)
    elif [ -n "$OPT_IN" ]; then PLAN_FOR=(--for "$OPT_IN")
    elif $SAMPLE_QC_ONLY; then PLAN_FOR=(--for sample-qc-data)
    elif $CLINICAL_DATA; then PLAN_FOR=(--for setup)
    else PLAN_FOR=(--for setup --skip clinical-data); fi
    og_setup plan "${PLAN_FOR[@]}" --check > /dev/null || exit 1
  fi
fi

# What a full setup could not install, one "what|retry command" entry each.
# The summary lists them, and a run with any left exits SETUP_PARTIAL_EXIT.
MISSING=()
SETUP_PARTIAL_EXIT=3
# The retry commands quote the script path, so a checkout path with a space
# still copies and pastes as one word.
SETUP_CMD=$(printf '%q' "$0")
note_missing() { MISSING+=("$1|$2"); }

print_missing() {
  local item
  echo ""
  echo "Missing (the steps named run without it or say it is missing; the warnings above say why):"
  for item in "${MISSING[@]}"; do
    echo "  - ${item%%|*}"
    echo "      retry: ${item#*|}"
  done
}

# NCBI's GRCh38 analysis set without ALT contigs: chr1-22, X, Y, M, the
# unplaced and unlocalized scaffolds and chrEBV, 195 sequences. NCBI lists the
# md5 of every file of the directory in md5checksums.txt, which is read at
# download time, and publishes the .fai beside the FASTA.
REF_FASTA_URL=${REF_FASTA_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz}
REF_FAI_URL="${REF_FASTA_URL%.gz}.fai"
if [ -z "${REF_FASTA_MD5+x}" ]; then REF_FASTA_MD5="$(dirname "$REF_FASTA_URL")/md5checksums.txt"; fi
if [ -z "${REF_FAI_MD5+x}" ]; then REF_FAI_MD5="$(dirname "$REF_FASTA_URL")/md5checksums.txt"; fi

# --- Step 33: somalier's sites and VerifyBamID2's marker panel ----------------------
# Pinned files, each checked against the sha256 recorded here. Step 33
# (scripts/33-sample-qc.sh) and validate-setup.sh read them under these names.
SOMALIER_SITES_URL=https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz
SOMALIER_SITES_SHA256=d1a853b8bb2e5f1a520bc67c8303be699543d225d785bce51e8335a2420489b5
# The panel of the VerifyBamID release VERIFYBAMID2_IMAGE packages (v2.0.3):
# one line per file, its URL and its sha256.
VB2_PANEL_FILES="https://raw.githubusercontent.com/Griffan/VerifyBamID/v2.0.3/resource/1000g.phase3.100k.b38.vcf.gz.dat.UD 259e320123756bb702542e3b3c4d766b481426d96cd63abd4b9f9fced035f375
https://raw.githubusercontent.com/Griffan/VerifyBamID/v2.0.3/resource/1000g.phase3.100k.b38.vcf.gz.dat.mu f7e8b4fad17cc433d887f18bf183fa4b89597d5a375f22133f44f39c5e35ffc8
https://raw.githubusercontent.com/Griffan/VerifyBamID/v2.0.3/resource/1000g.phase3.100k.b38.vcf.gz.dat.bed 0025d782137e5906bbc7afc553c49d2c82c5b8f0bfa0649bc415887442f42869"

# install_sample_qc_data: install both under GENOME_DIR/reference unless they
# are there with their pinned sha256. A file already there is checked again on
# every run: one with another sha256 (copied by hand, damaged, or from an older
# pin) is downloaded again. Returns 1 when a download fails its check.
install_sample_qc_data() {
  local dest="${GENOME_DIR}/reference" url sha rc=0 f
  f="${dest}/somalier/sites.hg38.vcf.gz"
  if [ -s "$f" ] && [ "$(_digest sha256 "$f")" = "$SOMALIER_SITES_SHA256" ]; then
    echo "[OK] somalier sites (step 33) already present (sha256 matches)."
  else
    if [ -s "$f" ]; then
      echo "[WARN] ${f} is not the pinned file (sha256 differs): downloading it again."
    fi
    if fetch "$SOMALIER_SITES_URL" "$f" sha256 "$SOMALIER_SITES_SHA256"; then
      echo "[OK] somalier sites (step 33): ${f}"
    else
      echo "[WARN] Could not install somalier's sites. Step 33 needs them; re-run setup.sh."
      rc=1
    fi
  fi
  while read -r url sha; do
    f="${dest}/verifybamid2/$(basename "$url")"
    if [ -s "$f" ] && [ "$(_digest sha256 "$f")" = "$sha" ]; then
      continue
    fi
    if [ -s "$f" ]; then
      echo "[WARN] ${f} is not the pinned file (sha256 differs): downloading it again."
    fi
    if ! fetch "$url" "$f" sha256 "$sha"; then
      echo "[WARN] Could not install $(basename "$url") of VerifyBamID2's panel. Step 33 needs it; re-run setup.sh."
      rc=1
    fi
  done <<<"$VB2_PANEL_FILES"
  [ "$rc" -eq 0 ] && echo "[OK] VerifyBamID2 panel (step 33): ${dest}/verifybamid2/"
  return "$rc"
}

if $SAMPLE_QC_ONLY; then
  install_sample_qc_data
  exit $?
fi

# --- Opt-in steps: Cyrius, Parascopy, KIR ------------------------------------------
# None of these is installed by a plain setup.sh run: each belongs to a step a
# default run leaves out, and each is asked for by its own flag.

# install_cyrius: Cyrius (step 21) into GENOME_DIR/tools/cyrius-<version>,
# installed by pip in PYTHON_IMAGE, the image step 21 runs it in, from
# scripts/cyrius-constraints.txt: every package pinned with the sha256 of its
# wheels (--require-hashes), nothing resolved (--no-deps) and nothing built
# (--only-binary). The network is used here, never when step 21 runs. Cyrius
# is under the PolyForm Strict licence 1.0.0 (non-commercial use only).
install_cyrius() {
  local dest="${GENOME_DIR}/tools/cyrius-${CYRIUS_VERSION}" lock="${PGP_ROOT}/scripts/cyrius-constraints.txt" stamp
  stamp="python=${PYTHON_IMAGE} lock=$(_digest sha256 "$lock")"
  if [ "$(cat "${dest}/INSTALLED" 2>/dev/null)" = "$stamp" ]; then
    echo "[OK] Cyrius ${CYRIUS_VERSION} (step 21) already installed: ${dest}"
    return 0
  fi
  echo "Installing Cyrius ${CYRIUS_VERSION} from ${lock} (hash-locked) into ${dest}"
  rm -rf "${dest}.part"
  mkdir -p "${dest}.part"
  # net: pip installs the locked wheels from PyPI, each checked by its hash
  # (a package install, not a data source).
  if ! run_in --net --rw "${GENOME_DIR}/tools" -v "${lock}:/lock.txt:ro" "${PYTHON_IMAGE}" \
      pip install --no-cache-dir --disable-pip-version-check -q --require-hashes --no-deps \
        --only-binary :all: --target "$(cpath "${dest}.part")" -r /lock.txt; then
    rm -rf "${dest}.part"
    echo "[WARN] Could not install Cyrius. Step 21 needs it; run: $0 --cyrius ${GENOME_DIR}"
    return 1
  fi
  echo "$stamp" > "${dest}.part/INSTALLED"
  rm -rf "$dest"
  mv "${dest}.part" "$dest"
  echo "[OK] Cyrius ${CYRIUS_VERSION} (step 21): ${dest}"
}

# install_parascopy_data: Parascopy's precomputed GRCh38 homology table and
# the model parameters of 1000 Genomes populations (step 35), Zenodo record
# PARASCOPY_DATA_RECORD (row parascopy_data of assets/data-sources.tsv reads
# it from here), each archive checked against its md5 there.
# shellcheck disable=SC2034  # og setup source reads it (and og dev freshness)
PARASCOPY_DATA_RECORD=15019940
PARASCOPY_FILES="GRCh38 a95bf674f43317d3a4c1b8ddbb140945
models_GRCh38_1KGP 244110d8fa883cf11334527ec2383498"
install_parascopy_data() {
  local dest="${GENOME_DIR}/reference/parascopy-${PARASCOPY_DATA_VERSION}" part md5 url raw
  if [ -s "${dest}/homology_table/GRCh38.bed.gz" ] && [ -s "${dest}/models_GRCh38_1KGP/EUR/SMN1.gz" ]; then
    echo "[OK] Parascopy ${PARASCOPY_DATA_VERSION} homology table and models (step 35) already present."
    return 0
  fi
  while read -r part md5; do
    url=$(source_of parascopy_data --part "$part" url) && raw=$(source_of parascopy_data --part "$part" path) || return 1
    adopt "${dest}.part/${raw##*/}" "$raw" "$url" md5 "$md5"
  done <<<"$PARASCOPY_FILES"
  rm -rf "${dest}.part"
  mkdir -p "${dest}.part"
  while read -r part md5; do
    url=$(source_of parascopy_data --part "$part" url) && raw=$(source_of parascopy_data --part "$part" path) || return 1
    if ! fetch "$url" "$raw" md5 "$md5" || ! tar -xzf "$raw" -C "${dest}.part"; then
      rm -rf "${dest}.part"
      echo "[WARN] Could not install Parascopy's data. Step 35 needs it; run: $0 --parascopy-data ${GENOME_DIR}"
      return 1
    fi
    done_with "$raw"
  done <<<"$PARASCOPY_FILES"
  rm -rf "$dest"  # not kept: the older install, replaced
  mv "${dest}.part" "$dest"
  echo "[OK] Parascopy ${PARASCOPY_DATA_VERSION} homology table and models (step 35): ${dest}"
}

# install_kir_data: kir.dat of IPD-KIR release KIR_DB_RELEASE (step 08 with
# KIR=true), from the release branch of the IPD-KIR repository (2.15.0 is
# branch 2150), checked against the md5 that release lists for it.
install_kir_data() {
  local branch=${KIR_DB_RELEASE//./} dest="${GENOME_DIR}/kir/IPD-KIR_${KIR_DB_RELEASE}/kir.dat" url md5
  if [ -s "$dest" ]; then
    echo "[OK] IPD-KIR ${KIR_DB_RELEASE} (step 08, KIR=true) already present."
    return 0
  fi
  url="https://raw.githubusercontent.com/ANHIG/IPDKIR/${branch}"
  # md5checksum.txt lines read "MD5 (kir.dat) = <md5>".
  # net: checksum, the md5 of kir.dat that the release lists
  md5=$(_get "${url}/md5checksum.txt" - 2>/dev/null | awk '$2 == "(kir.dat)" {print $NF}') || md5=""
  if [ -z "$md5" ] || ! fetch "${url}/kir.dat" "$dest" md5 "$md5"; then
    echo "[WARN] Could not install IPD-KIR ${KIR_DB_RELEASE} (is it a release with a branch ${branch}?)."
    return 1
  fi
  if ! grep -q "IPD-KIR Release Version ${KIR_DB_RELEASE}" "$dest"; then
    rm -f "$dest"  # not kept: the wrong release, which step 08 must not read
    echo "[WARN] ${url}/kir.dat does not say IPD-KIR Release Version ${KIR_DB_RELEASE}; removed it."
    return 1
  fi
  echo "[OK] IPD-KIR ${KIR_DB_RELEASE} (step 08, KIR=true): ${dest}"
}

# install_yleaf_data: the hg38 and hg_prediction_tables folders of yleaf/data
# in Yleaf YLEAF_DATA_VERSION (step 37): its GRCh38 marker positions and
# haplogroup trees, the data YLEAF_IMAGE ships too. GitHub's archive of the tag
# v<version>, checked against YLEAF_DATA_SHA256; the other builds (hg19, T2T)
# are left out. Then the current Y tree beside it, in data/ytree
# (og y tree-place): YFull's YTree YTREE_VERSION and the GRCh38 SNP tables
# of Yleaf YTREE_POSITIONS_VERSION, each at a pinned commit and sha256, and
# the chrY strict accessibility mask of the 1000 Genomes Project
# (YCALLABLE_MASK_RELEASE) that the private SNP list keeps to. A failed
# download of the current tree is a warning, not a failure. The mask is
# checked against YCALLABLE_MASK_SHA256 on every run: one that matches is
# kept with no download, also when a new tree pin installs the tables again,
# and one that fails is downloaded again on its own.
install_yleaf_data() {
  local dest="${GENOME_DIR}/reference/yleaf-${YLEAF_DATA_VERSION}" tgz url
  if [ -s "${dest}/data/hg38/new_positions.txt" ] && [ -s "${dest}/data/hg_prediction_tables/tree.json" ]; then
    echo "[OK] Yleaf ${YLEAF_DATA_VERSION} marker tables (step 37) already present."
  else
    url=$(source_of yleaf_data url) && tgz=$(source_of yleaf_data path) || return 1
    adopt "${dest}.tar.gz" "$tgz" "$url" sha256 "$YLEAF_DATA_SHA256"
    if ! fetch "$url" "$tgz" sha256 "$YLEAF_DATA_SHA256"; then
      echo "[WARN] Could not download Yleaf ${YLEAF_DATA_VERSION}."
      return 1
    fi
    rm -rf "${dest}.part"
    mkdir -p "${dest}.part"
    tar -xzf "$tgz" -C "${dest}.part" --strip-components 2 \
      "Yleaf-${YLEAF_DATA_VERSION}/yleaf/data/hg38" "Yleaf-${YLEAF_DATA_VERSION}/yleaf/data/hg_prediction_tables"
    if [ ! -s "${dest}.part/data/hg38/new_positions.txt" ] || [ ! -s "${dest}.part/data/hg_prediction_tables/tree.json" ]; then
      echo "[WARN] The Yleaf ${YLEAF_DATA_VERSION} archive has no yleaf/data/hg38/new_positions.txt or hg_prediction_tables/tree.json."
      rm -rf "${dest}.part"
      return 1
    fi
    done_with "$tgz"
    rm -rf "$dest"  # not kept: the older install, replaced
    mv "${dest}.part" "$dest"
    echo "[OK] Yleaf ${YLEAF_DATA_VERSION} marker tables and trees (step 37): ${dest}/data"
  fi

  # version.txt names every pin (releases, commits, sha256), so a re-pin of
  # any of them, at the same tree release too, installs the tables again.
  local yt="${dest}/data/ytree" label="YFull YTree ${YTREE_VERSION} (CC BY 4.0)" want mask_ok=false
  local mask_url tree_url ftdna_url yfull_url
  mask_url=$(source_of y_callable_mask url) && tree_url=$(source_of ytree url) \
    && ftdna_url=$(source_of ytree_ftdna_positions url) && yfull_url=$(source_of ytree_yfull_positions url) || return 1
  local mask="${yt}/callable_mask_chrY.fasta.gz"
  if [ -s "$mask" ] && [ "$(_digest sha256 "$mask")" = "$YCALLABLE_MASK_SHA256" ]; then
    mask_ok=true
  fi
  want="$(printf 'tree\t%s\npositions\tYleaf %s GRCh38 SNP tables (GPL-3.0)\npins\t%s %s %s %s %s %s %s' \
    "$label" "$YTREE_POSITIONS_VERSION" "$YTREE_COMMIT" "$YTREE_SHA256" \
    "$YTREE_POSITIONS_COMMIT" "$YTREE_FTDNA_SHA256" "$YTREE_YFULL_SHA256" \
    "$YCALLABLE_MASK_RELEASE" "$YCALLABLE_MASK_SHA256")"
  if [ -s "${yt}/version.txt" ] && [ "$(cat "${yt}/version.txt")" = "$want" ]; then
    if $mask_ok; then
      echo "[OK] Current Y tree (${label}, step 37) already present."
      return 0
    fi
    # The tables are current but the mask is missing or fails its sha256:
    # fetch it alone. Without it step 37 skips the private SNP list.
    echo "[WARN] ${mask} is missing or fails its sha256; downloading it again."
    rm -f "$mask"  # not kept: an installed mask that fails its sha256
    if fetch "$mask_url" "$mask" sha256 "$YCALLABLE_MASK_SHA256"; then
      echo "[OK] Current Y tree (${label}, step 37) present, chrY callable mask downloaded again."
    else
      echo "[WARN] Could not download the chrY callable mask; step 37 will skip the private SNP list. Run this again to retry."
    fi
    return 0
  fi
  rm -rf "${yt}.part"
  mkdir -p "${yt}.part"
  # A mask that matches its pin is kept: a new tree pin does not download it again.
  # The two position tables are not: they are fetched again with the tree, even
  # when YTREE_POSITIONS_COMMIT did not move.
  if $mask_ok; then
    cp "$mask" "${yt}.part/callable_mask_chrY.fasta.gz"
  fi
  if ! fetch "$tree_url" "${yt}.part/current_tree.json" sha256 "$YTREE_SHA256" \
     || ! fetch "$ftdna_url" "${yt}.part/positions_ftdna_hg38.txt" sha256 "$YTREE_FTDNA_SHA256" \
     || ! fetch "$yfull_url" "${yt}.part/positions_yfull_hg38.txt" sha256 "$YTREE_YFULL_SHA256" \
     || { ! $mask_ok && ! fetch "$mask_url" "${yt}.part/callable_mask_chrY.fasta.gz" sha256 "$YCALLABLE_MASK_SHA256"; }; then
    # Yleaf's own data is in place, so step 37 still runs: it says the
    # current tree is not installed and writes Yleaf's result only. An older
    # tree from an earlier pin goes too, so no step places on tables that
    # versions.env no longer names. A mask that matches its pin stays: with
    # no version.txt beside it no step reads it, and the next run keeps it.
    echo "[WARN] Could not download the current Y tree (${label}); step 37 will place the sample on Yleaf's tree only. Run this again to retry."
    rm -rf "${yt}.part"
    if $mask_ok; then
      find "$yt" -mindepth 1 -maxdepth 1 ! -name callable_mask_chrY.fasta.gz -exec rm -rf {} +  # not kept: the installed tree of an older pin (its position tables are fetched again with the tree); the mask that matches its pin stays
    else
      rm -rf "$yt"  # not kept: the installed tree of an older pin; its position tables are fetched again with the tree
    fi
    return 0
  fi
  printf '%s\n' "$want" > "${yt}.part/version.txt"
  rm -rf "$yt"  # not kept: the installed tree of an older pin, replaced; its position tables were fetched again with the tree
  mv "${yt}.part" "$yt"
  echo "[OK] Current Y tree (${label}) and Yleaf ${YTREE_POSITIONS_VERSION}'s SNP positions (step 37): ${yt}"
}

# install_viral_data: the viral genomes of step 47, VIRAL_REF_ACCESSIONS in
# that order as one FASTA from NCBI's E-utilities, checked against
# VIRAL_REF_SHA256. minimap2 indexes a file this small in a second at each run.
install_viral_data() {
  local dest="${GENOME_DIR}/reference/viral-${VIRAL_REF_VERSION}/viral_genomes.fa" url
  if [ -s "$dest" ] && [ "$(_digest sha256 "$dest")" = "$VIRAL_REF_SHA256" ]; then
    echo "[OK] Viral genomes ${VIRAL_REF_VERSION} (step 47) already present."
    return 0
  fi
  rm -f "$dest"  # not kept: an installed file that fails its sha256
  url=$(source_of viral_genomes url) || return 1
  if ! fetch "$url" "$dest" sha256 "$VIRAL_REF_SHA256"; then
    echo "[WARN] Could not install the viral genomes of step 47 (NCBI E-utilities)."
    return 1
  fi
  echo "[OK] Viral genomes ${VIRAL_REF_VERSION} (step 47, $(grep -c '^>' "$dest") sequences): ${dest}"
}

# need_og WHAT: build og (bin/.og) when it is missing or stale; when that
# fails, warn that og WHAT and return 1.
need_og() {
  "${PGP_ROOT}/rust/og.sh" check --required > /dev/null 2>&1 && return 0
  "${PGP_ROOT}/rust/og.sh" install && return 0
  echo "[WARN] Could not build og, which ${1}: run ${0} --og and look at its error."
  return 1
}

# install_masked_reference (--masked-reference): REF_FASTA (the default
# no-ALT analysis set) with the false duplications of GRCh38 set to N, in
# <REF_FASTA without extension>_maskedGRC.fasta beside it, with its .fai and
# .dict. A false duplication is a stretch GRCh38 holds twice where people
# have one copy (chr21 holds second copies of KCNE1, CBS, U2AF1 and CRYAA):
# reads of the real gene fit both copies and get MAPQ 0. Masked, the false
# copy takes no reads. Sequence names, lengths and coordinates do not change.
# The list is GRC_EXCLUSIONS_URL (sha256 checked); the FASTA it writes must
# have GRC_MASKED_REF_SHA256, or MASKED_REF_SHA256 when that is set (a source
# other than NCBI's file, as the CI case's chr21), or it is deleted. The
# masking is og's `og setup mask-fasta` (rust/crates/og-intervals), run in
# OG_IMAGE with bin/ mounted, so og is built first when it is missing or
# stale. The command uses relative paths from the reference directory.
install_masked_reference() {
  local src=$REF_FASTA dest bed want got dir log
  case "$src" in
    *.gz) echo "[WARN] --masked-reference needs an unpacked FASTA; REF_FASTA is ${src}."; return 1 ;;
    *_maskedGRC.*) echo "[WARN] REF_FASTA (${src}) is already a masked reference. Unset it, so the mask is built from the default reference."; return 1 ;;
  esac
  dest="${src%.*}_maskedGRC.fasta"
  dir=$(dirname "$dest")
  want=${MASKED_REF_SHA256:-$GRC_MASKED_REF_SHA256}
  bed="${dir}/${GRC_EXCLUSIONS_URL##*/}"
  if [ -s "$dest" ] && [ -s "${dest}.fai" ] && [ -s "${dest%.*}.dict" ] && [ "$(_digest sha256 "$dest")" = "$want" ]; then
    echo "[OK] Masked reference already present (sha256 matches): ${dest}"
  else
    for f in "$src" "${src}.fai"; do
      [ -s "$f" ] || { echo "[WARN] ${f} is missing: install the reference first (${0} ${GENOME_DIR})."; return 1; }
    done
    if [ -s "$bed" ] && [ "$(_digest sha256 "$bed")" = "$GRC_EXCLUSIONS_SHA256" ]; then
      :
    elif ! fetch "$GRC_EXCLUSIONS_URL" "$bed" sha256 "$GRC_EXCLUSIONS_SHA256"; then
      echo "[WARN] Could not download the GRC's list of false duplications (${GRC_EXCLUSIONS_URL})."
      return 1
    fi
    # Only the partial output goes now: a masked reference already there stays
    # until its replacement has passed the sha256 check.
    rm -f "${dest}.part"
    if ! "${PGP_ROOT}/rust/og.sh" check --required > /dev/null 2>&1 && ! "${PGP_ROOT}/rust/og.sh" install; then
      echo "[WARN] Could not build og, which masks the reference: run ${0} --og and look at its error."
      return 1
    fi
    echo "Masking the false duplications of $(basename "$src") (under a minute for GRCh38)..."
    if ! log=$(run_in --rw "$dir" -v "${PGP_ROOT}/bin:/pgp-bin:ro" -w "$(cpath "$dir")" "$OG_IMAGE" \
          /pgp-bin/og setup mask-fasta "$(basename "$src")" "$(basename "$bed")" "$(basename "$dest").part"); then
      rm -f "${dest}.part"
      echo "[WARN] Could not write the masked reference."
      return 1
    fi
    echo "  ${log}"
    got=$(_digest sha256 "${dest}.part")
    if [ "$got" != "$want" ]; then
      rm -f "${dest}.part"
      echo "[WARN] The masked reference has sha256 ${got}, not ${want}: $(basename "$src") is not the file it is pinned for (NCBI's no-ALT analysis set). Removed it."
      return 1
    fi
    mv -f "${dest}.part" "$dest"
    rm -f "${dest}.fai" "${dest%.*}.dict"
    # The index and the dictionary go next to the FASTA, so its directory is writable here.
    if ! run_in --rw "$dir" "$SAMTOOLS_IMAGE" samtools faidx "$(cpath "$dest")" \
       || ! run_in --rw "$dir" "$SAMTOOLS_IMAGE" samtools dict -o "$(cpath "${dest%.*}.dict")" "$(cpath "$dest")"; then
      rm -f "${dest}.fai" "${dest%.*}.dict"
      echo "[WARN] Could not index the masked reference. Kept the FASTA (its sha256 matches) and removed its partial index; run setup again."
      return 1
    fi
    echo "[OK] Masked reference (sha256 matches): ${dest}"
  fi
  echo "  Use it for a sample by realigning it (docs/realignment.md#the-masked-reference):"
  printf '    export REF_FASTA=%q   # bash steps; the Nextflow run takes --reference %q\n' "${dest#"${GENOME_DIR}"/}" "$dest"
}

# install_exomiser_data: Exomiser's hg38 and phenotype data of release
# EXOMISER_DATA_VERSION (step 48) into reference/exomiser-<release>/, the
# folder Exomiser reads as its data directory. Each archive is checked against
# its sha256 in versions.env, unpacked by og (`og setup unzip`, so no unzip
# on the host is needed; each member streams to disk), and kept in the downloads folder. A
# folder is moved into place only once whole.
install_exomiser_data() {
  local dest="${GENOME_DIR}/reference/exomiser-${EXOMISER_DATA_VERSION}" part sha zip url
  if [ -d "${dest}/${EXOMISER_DATA_VERSION}_hg38" ] && [ -d "${dest}/${EXOMISER_DATA_VERSION}_phenotype" ]; then
    echo "[OK] Exomiser ${EXOMISER_DATA_VERSION} data (step 48) already present: ${dest}"
    return 0
  fi
  echo "  ${EXOMISER_TERMS} See docs/48-exomiser.md."
  mkdir -p "$dest"
  while read -r part sha; do
    [ -d "${dest}/${EXOMISER_DATA_VERSION}_${part}" ] && continue
    url=$(source_of exomiser_data --part "$part" url) && zip=$(source_of exomiser_data --part "$part" path) || return 1
    adopt "${dest}/${zip##*/}" "$zip" "$url" sha256 "$sha"
    # A zip kept by an unpack that failed is used again once its sha256 still matches.
    if [ -f "$zip" ] && [ "$(_digest sha256 "$zip")" = "$sha" ]; then
      echo "Using the checked ${zip##*/} kept from the last run."
    elif echo "Downloading Exomiser's ${EXOMISER_DATA_VERSION}_${part}.zip (resumes when run again)..." \
        && ! fetch "$url" "$zip" sha256 "$sha"; then
      echo "[WARN] Could not download Exomiser's ${part} data. Step 48 needs it; run: $0 --exomiser-data ${GENOME_DIR}"
      return 1
    fi
    rm -rf "${dest}/unpack.part"
    need_og "unpacks Exomiser's data" || return 1
    if ! run_in --rw "$dest" --cpus 1 --memory 2g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$OG_IMAGE" \
        /pgp-bin/og setup unzip "$(cpath "$zip")" "$(cpath "${dest}/unpack.part")" \
        || [ ! -d "${dest}/unpack.part/${EXOMISER_DATA_VERSION}_${part}" ]; then
      rm -rf "${dest}/unpack.part"
      echo "[WARN] Could not unpack ${zip}; it is kept: run this command again."
      return 1
    fi
    mv "${dest}/unpack.part/${EXOMISER_DATA_VERSION}_${part}" "${dest}/${EXOMISER_DATA_VERSION}_${part}"
    rm -rf "${dest}/unpack.part"
    done_with "$zip"
  done <<EOF2
hg38 ${EXOMISER_HG38_SHA256}
phenotype ${EXOMISER_PHENOTYPE_SHA256}
EOF2
  echo "[OK] Exomiser ${EXOMISER_DATA_VERSION} hg38 and phenotype data (step 48): ${dest}"
}

# install_chip_data: Ensembl ENSEMBL_GFF_RELEASE's GRCh38 GFF3 (step 49 reads
# the CHIP genes' MANE Select transcripts from it), checked against
# ENSEMBL_GFF_SHA256. A file that is there with another sha256 is downloaded again.
install_chip_data() {
  local dest url
  dest=$(source_of ensembl_gff path) && url=$(source_of ensembl_gff url) || return 1
  if [ -s "$dest" ] && [ "$(_digest sha256 "$dest")" = "$ENSEMBL_GFF_SHA256" ]; then
    echo "[OK] Ensembl ${ENSEMBL_GFF_RELEASE} gene annotation (step 49) already present."
    return 0
  fi
  rm -f "$dest"  # not kept: an installed file that fails its sha256
  if ! fetch "$url" "$dest" sha256 "$ENSEMBL_GFF_SHA256"; then
    echo "[WARN] Could not download Ensembl ${ENSEMBL_GFF_RELEASE}'s gene annotation. Step 49 needs it; run setup.sh --chip-data again."
    return 1
  fi
  echo "[OK] Ensembl ${ENSEMBL_GFF_RELEASE} gene annotation (step 49): ${dest}"
}

# megane_kmers_path: where step 50 reads MEGAnE's k-mer set of REF_FASTA, the
# path scripts/50-mei.sh and run-all.sh build the same way (MEGANE_KMERS
# overrides it): reference/megane/<reference file name without extension>.mk,
# with MEGAnE's .mi count file beside it.
megane_kmers_path() {
  local n p
  n=$(basename "$REF_FASTA"); n=${n%.gz}; n=${n%.*}
  p=${MEGANE_KMERS:-${GENOME_DIR}/reference/megane/${n}.mk}
  # A relative MEGANE_KMERS is inside GENOME_DIR, as scripts/50-mei.sh reads it.
  case "$p" in /*) ;; *) p="${GENOME_DIR}/${p}" ;; esac
  printf '%s\n' "$p"
}

# install_megane_kmers: MEGAnE's k-mer set of the reference (step 50): every
# 32-mer seen more than once in it, on either strand, which MEGAnE uses to
# drop reads that map to many places. It is built from the reference itself,
# so a new reference needs a new set; one newer than the reference is kept.
# og builds it (`og setup megane-kmers`, in OG_IMAGE with bin/ mounted), byte
# for byte the files MEGAnE's own build_kmerset writes, in og's 5 GB budget:
# the genome at 4 bits per base, then one range of k-mers at a time.
# MEGAnE's own builder peaked at about 67 GB for GRCh38, RAM plus swap.
install_megane_kmers() {
  local mk dir need=6 have log
  mk=$(megane_kmers_path); dir=$(dirname "$mk")
  if [ ! -s "$REF_FASTA" ] || [ ! -s "${REF_FASTA}.fai" ]; then
    echo "[WARN] Step 50's k-mer set is built from the reference, and ${REF_FASTA} (with its .fai) is not there yet: run ./scripts/setup.sh ${GENOME_DIR} first."
    return 1
  fi
  if [ -s "$mk" ] && [ -s "${mk%.mk}.mi" ] && [ "$mk" -nt "$REF_FASTA" ]; then
    echo "[OK] MEGAnE k-mer set (step 50) already present: ${mk}"
    return 0
  fi
  # Docker names the engine's memory .MemTotal; podman has it under .Host.
  have=$("$CONTAINER_ENGINE" info --format '{{.MemTotal}}' 2>/dev/null) || have=""
  case "$have" in ''|*[!0-9]*) have=$("$CONTAINER_ENGINE" info --format '{{.Host.MemTotal}}' 2>/dev/null) || have="" ;; esac
  have=$(printf '%s\n' "$have" | awk '$1 ~ /^[0-9]+$/ { printf "%d", $1 / 1e9 }')
  if [ -n "$have" ] && [ "$have" -gt 0 ] && [ "$have" -lt "$need" ]; then
    echo "[WARN] MEGAnE's k-mer set needs about ${need} GB of RAM and ${CONTAINER_ENGINE} has ${have} GB."
    echo "       Build it on a machine with more memory (./scripts/setup.sh --megane-kmers <genome_dir> there) and copy ${dir#"$GENOME_DIR"/}/ here."
    return 1
  fi
  need_og "builds MEGAnE's k-mer set" || return 1
  mkdir -p "$dir"
  echo "Building MEGAnE's k-mer set of $(basename "$REF_FASTA") (step 50; a few minutes, about 5 GB of RAM)..."
  # og writes each file under a hidden temporary name (.<name>.tmp-<pid>-<n>)
  # and renames it at the end; a build killed on the way (out of memory,
  # Ctrl-C) leaves that file behind, up to 1.2 GB for GRCh38.
  rm -f "$mk" "${mk%.mk}.mi" "${dir}/.$(basename "$mk").tmp-"* "${dir}/.$(basename "${mk%.mk}").mi.tmp-"*
  if ! log=$(run_in --rw "$dir" -v "${PGP_ROOT}/bin:/pgp-bin:ro" --memory "${need}g" "$OG_IMAGE" \
        /pgp-bin/og setup megane-kmers --fa "$REF_FASTA_C" --out "$(cpath "$dir")" \
        --prefix "$(basename "${mk%.mk}")" --mem 5) || [ ! -s "$mk" ] || [ ! -s "${mk%.mk}.mi" ]; then
    rm -f "$mk" "${mk%.mk}.mi"
    echo "[WARN] Could not build MEGAnE's k-mer set (step 50); og's error is above."
    return 1
  fi
  echo "  ${log}"
  echo "[OK] MEGAnE k-mer set (step 50, $(cat "${mk%.mk}.mi") repeated 32-mers): ${mk}"
}

# install_strchive_data: STRchive's ExpansionHunter catalog and loci file
# (STRCHIVE_VERSION), each checked against its sha256 in versions.env, joined
# by og (`og str catalog`) into reference/strchive-<version>/catalog.json,
# the catalog steps 09 and 9b and the pipeline use by default, and
# str_ranges.tsv beside it, the size ranges of each locus. A file that
# fails its check is removed and nothing is converted; the steps then use
# the catalogs inside their images. catalog.stamp names the inputs and the
# converter's output format, so a changed one converts again: bump format
# below whenever og str catalog writes something different.
# og uses relative paths from the catalog's directory in OG_IMAGE.
# STR_EXTRA_LOCI=yes (opt-in) also appends assets/str_extra_loci.json, the
# repeats STRchive leaves out (CFTR TG/poly-T, UGT1A1 TA): ExpansionHunter
# sizes them and Stranger leaves them ungraded (docs/09-str-expansions.md).
# The stamp names that file's sha256, or extra=off.
install_strchive_data() {
  local dir="${GENOME_DIR}/reference/strchive-${STRCHIVE_VERSION}" url name sha stamp format=og-str-catalog-1
  local extra_json="${PGP_ROOT}/assets/str_extra_loci.json" extra_sha=off extra_args=()
  case "${STR_EXTRA_LOCI:-no}" in
    yes) extra_sha=$(_digest sha256 "$extra_json") extra_args=(--extra /pgp-extra/str_extra_loci.json) ;;
    no) ;;
    *) echo "[WARN] STR_EXTRA_LOCI is '${STR_EXTRA_LOCI}', not yes or no: the extra repeats are left out." ;;
  esac
  # The URLs are rows strchive_eh and strchive_loci of assets/data-sources.tsv,
  # which scripts/ci/check-links.sh asks too.
  local eh_url loci_url
  eh_url=$(source_of strchive_eh url) && loci_url=$(source_of strchive_loci url) || return 1
  local eh=${eh_url##*/} loci=${loci_url##*/}
  stamp="release=${STRCHIVE_VERSION} eh=${STRCHIVE_EH_SHA256} loci=${STRCHIVE_LOCI_SHA256} converter=${format} extra=${extra_sha}"
  if [ -s "${dir}/catalog.json" ] && [ -s "${dir}/str_ranges.tsv" ] && [ "$(cat "${dir}/catalog.stamp" 2>/dev/null)" = "$stamp" ]; then
    echo "[OK] STRchive ${STRCHIVE_VERSION} catalog (steps 09 and 9b) already present."
    return 0
  fi
  # og first: when it cannot be built, the catalog already there stays.
  need_og "converts the STRchive catalog" || return 1
  # A stale catalog goes next, so a failed refetch leaves none for step 09.
  rm -f "${dir}/catalog.json" "${dir}/str_ranges.tsv" "${dir}/catalog.stamp"
  for url in "${eh_url}|${STRCHIVE_EH_SHA256}" "${loci_url}|${STRCHIVE_LOCI_SHA256}"; do
    sha=${url#*|} url=${url%%|*} name=${url##*/}
    [ -s "${dir}/${name}" ] && [ "$(_digest sha256 "${dir}/${name}")" = "$sha" ] && continue
    rm -f "${dir}/${name}"  # not kept: an installed file that fails its sha256
    if ! fetch "$url" "${dir}/${name}" sha256 "$sha"; then
      echo "[WARN] Could not install STRchive's ${name}. Steps 09 and 9b use the catalogs inside their images; re-run setup.sh to try again."
      return 1
    fi
  done
  if ! run_in --rw "$dir" -v "${PGP_ROOT}/bin:/pgp-bin:ro" -v "${PGP_ROOT}/assets:/pgp-extra:ro" \
       -w "$(cpath "$dir")" "$OG_IMAGE" /pgp-bin/og str catalog \
       --eh "$eh" --loci "$loci" --release "$STRCHIVE_VERSION" \
       --out catalog.json --ranges str_ranges.tsv ${extra_args[@]+"${extra_args[@]}"} \
       || [ ! -s "${dir}/catalog.json" ] || [ ! -s "${dir}/str_ranges.tsv" ]; then
    rm -f "${dir}/catalog.json" "${dir}/str_ranges.tsv"
    echo "[WARN] Could not convert the STRchive ${STRCHIVE_VERSION} catalog (og str catalog). Steps 09 and 9b use the catalogs inside their images."
    return 1
  fi
  printf '%s\n' "$stamp" > "${dir}/catalog.stamp"
  echo "[OK] STRchive ${STRCHIVE_VERSION} catalog (steps 09 and 9b): ${dir}/catalog.json"
}

# install_tandem_repeats: the tandem-repeat map step 04c passes to Sniffles2
# (--tandem-repeats): TANDEM_REPEATS_URL, checked against
# TANDEM_REPEATS_SHA256, and its unpacked copy beside it, because Sniffles2
# reads a plain BED only. The copy is written under a temporary name first,
# so a file under the final name is whole. Without the map step 04c runs as
# before and says so.
install_tandem_repeats() {
  local gz="${GENOME_DIR}/reference/${TANDEM_REPEATS_URL##*/}"
  local bed="${gz%.gz}"
  if [ -s "$bed" ] && [ -s "$gz" ] && [ "$(_digest sha256 "$gz")" = "$TANDEM_REPEATS_SHA256" ]; then
    echo "[OK] GIAB tandem-repeat map (step 04c) already present."
    return 0
  fi
  rm -f "$bed"
  if ! fetch "$TANDEM_REPEATS_URL" "$gz" sha256 "$TANDEM_REPEATS_SHA256"; then
    echo "[WARN] Could not install GIAB's tandem-repeat map (${TANDEM_REPEATS_URL}). Step 04c runs Sniffles2 without it; re-run setup.sh to try again."
    return 1
  fi
  if ! gzip -dc "$gz" > "${bed}.part" || [ ! -s "${bed}.part" ]; then
    rm -f "${bed}.part"
    echo "[WARN] Could not unpack ${gz}. Step 04c runs Sniffles2 without the tandem-repeat map; re-run setup.sh to try again."
    return 1
  fi
  mv -f "${bed}.part" "$bed"
  echo "[OK] GIAB tandem-repeat map (step 04c): ${bed}"
}

# install_consumer_sites: the site table of the raw-data export (step 40). The
# 23andMe v5 site list (GRCh37 positions and reference bases) of
# CONSUMER_SITES_COMMIT, checked against CONSUMER_SITES_SHA256, is joined by
# rsID with dbSNP CONSUMER_DBSNP_BUILD's GRCh38 VCF (about 30 GB, checked
# against CONSUMER_DBSNP_MD5) in og's `og consumer build-sites`
# (rust/crates/og-steps), which takes the GRCh38 position from dbSNP instead
# of a liftover, each site's strand from UCSC's hg19ToHg38 chain (the one
# step 44 uses, checked against HG19_TO_HG38_CHAIN_SHA256), and lists every
# site it drops with the reason. It reads the
# BGZF file on 4 cores and looks at the ID column of a line before anything
# else, so no grep pass runs first. It runs in OG_IMAGE with bin/
# mounted, so og is built first when it is missing or stale.
# The dbSNP file stays in the downloads folder (done_with). The table built is then
# compared with CONSUMER_SITES_TABLE_SHA256: a difference is a warning.
install_consumer_sites() {
  local dir="${GENOME_DIR}/reference/consumer" tpl db out dropped url chain
  out="${dir}/23andme_v5_sites_dbsnp${CONSUMER_DBSNP_BUILD}.tsv.gz"
  dropped="${dir}/23andme_v5_dropped_dbsnp${CONSUMER_DBSNP_BUILD}.tsv"
  chain="${GENOME_DIR}/reference/liftover/hg19ToHg38.over.chain.gz"
  if [ -s "$out" ]; then
    consumer_sites_pin_check "$out" "already present"
    return 0
  fi
  tpl=$(source_of consumer_sites_template path) && url=$(source_of consumer_sites_template url) || return 1
  if ! fetch "$url" "$tpl" sha256 "$CONSUMER_SITES_SHA256"; then
    echo "[WARN] Could not download the 23andMe v5 site list."
    return 1
  fi
  [ -s "$chain" ] || echo "  ${UCSC_CHAIN_TERMS} See docs/40-export-consumer.md."
  if ! fetch "$HG19_TO_HG38_CHAIN_URL" "$chain" sha256 "$HG19_TO_HG38_CHAIN_SHA256"; then
    echo "[WARN] Could not download UCSC's hg19ToHg38 chain; run: $0 --consumer-sites ${GENOME_DIR}"
    return 1
  fi
  # In the downloads folder: fetch reuses a copy whose md5 matches the pin
  # (a rerun after a failed build, another install) and refuses one that does
  # not, so the table is never built from another file.
  db=$(source_of dbsnp path) && url=$(source_of dbsnp url) || return 1
  adopt "${dir}/${db##*/}" "$db" "$url" md5 "$CONSUMER_DBSNP_MD5"
  echo "dbSNP ${CONSUMER_DBSNP_BUILD} for GRCh38 (about 30 GB, kept in ${OG_DOWNLOADS}; a copy already there with the pinned md5 is used)..."
  if ! fetch "$url" "$db" md5 "$CONSUMER_DBSNP_MD5"; then
    echo "[WARN] Could not download dbSNP ${CONSUMER_DBSNP_BUILD}. Run this command again: the download resumes."
    return 1
  fi
  if ! "${PGP_ROOT}/rust/og.sh" check --required > /dev/null 2>&1 && ! "${PGP_ROOT}/rust/og.sh" install; then
    echo "[WARN] Could not build og, which builds the site table: run ${0} --og and look at its error."
    return 1
  fi
  echo "Building the site table (reads the whole dbSNP file once)..."
  # The table is written to reference/consumer/, outside any sample directory.
  if ! run_in --rw "$dir" -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 4 --memory 4g "$OG_IMAGE" \
       /pgp-bin/og consumer build-sites --template "$(cpath "$tpl")" --dbsnp "$(cpath "$db")" \
       --chain "$(cpath "$chain")" --out "$(cpath "$out")" --dropped "$(cpath "$dropped")" \
       --source "23andme2vcf ${CONSUMER_SITES_COMMIT} 23andme_v5_hg19_ref.txt.gz; dbSNP ${CONSUMER_DBSNP_BUILD} GCF_000001405.40; UCSC hg19ToHg38.over.chain.gz"; then
    rm -f "$out"
    echo "[WARN] Could not build the site table of step 40; the dbSNP download is kept: run this command again."
    return 1
  fi
  done_with "$db"
  consumer_sites_pin_check "$out" "built (dropped sites and their reasons: ${dropped})"
}

# consumer_sites_pin_check TABLE WHAT: compare the site table with
# CONSUMER_SITES_TABLE_SHA256 and say so. A difference is a warning: both
# inputs are checked against their own pins, so it means build-sites changed
# since the pin was taken, and step 40 can still use the table.
consumer_sites_pin_check() {
  local got
  got=$(consumer_sites_digest "$1")
  if [ "$got" = "$CONSUMER_SITES_TABLE_SHA256" ]; then
    echo "[OK] Site table of step 40 $2: $1, the same as the pinned build"
  else
    echo "[WARN] Site table of step 40 $2: $1 differs from the pinned build (content sha256 ${got}, pinned ${CONSUMER_SITES_TABLE_SHA256})."
    echo "       It was built by another version of og consumer build-sites. Step 40 uses it; to rebuild it, remove it and run this command again."
  fi
}

# --- Step 44: the ancient populations panel (opt-in) ---------------------------------
# The Allen Ancient DNA Resource's 1240K files (AADR_VERSION, by Dataverse file
# id, each checked against the md5 Dataverse publishes) become the panel of
# step 44 in reference/ancient/aadr_<release>:
#   panel.ind, poplist.txt, ancients.tsv  who builds the axes (the present-day
#        West Eurasians of og ancient's MODERN_GROUPS) and who is
#        projected (dated West Eurasians with 100,000 sites or more);
#   sites_grch38.tsv  the autosomal sites on GRCh38 (UCSC's hg19ToHg38 chain,
#        each kept only when one of its alleles is the GRCh38 base);
#   moderns.geno/.snp/.ind  the present-day people in EIGENSTRAT text, which
#        step 44 adds the person to;
#   panel.evec, panel.eval  smartpca's axes with the ancient people projected
#        (lsqproject), read by step 44's report;
#   DONE, written last, once a second smartpca run on moderns.* alone gives
#        the same axes, which is how step 44 runs.
# smartpca holds the whole genotype file in memory (7.1 GB for the 23,089
# people of v66.p1), so the install gives it a 12 GB container and needs that
# much free RAM. The AADR downloads (7.3 GB) stay in the downloads folder
# (done_with). ANCIENT_SOURCE_DIR=<folder> takes the four files (.anno, .ind, .snp,
# .geno) from that folder instead, unchecked (the tests' synthetic panel),
# staged in raw_<version> beside the panel and never in the downloads folder.
HG19_TO_HG38_CHAIN_URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
install_ancient_panel() {
  local base="${GENOME_DIR}/reference/ancient" dir part raw old chain e url md5 f src=${ANCIENT_SOURCE_DIR:-}
  dir="${base}/aadr_${AADR_VERSION}" part="${base}/.aadr_${AADR_VERSION}.part" old="${base}/raw_${AADR_VERSION}"
  chain="${GENOME_DIR}/reference/liftover/hg19ToHg38.over.chain.gz"
  if [ -s "${dir}/DONE" ]; then
    echo "[OK] ancient populations panel AADR ${AADR_VERSION} (step 44) already present: ${dir}"
    return 0
  fi
  for f in "$REF_FASTA" "${REF_FASTA}.fai"; do
    [ -s "$f" ] || { echo "[WARN] ${f} is missing: install the reference first (${0} ${GENOME_DIR})."; return 1; }
  done
  if [ -z "${ANCIENT_CHAIN:-}" ] && [ ! -s "$chain" ]; then
    echo "  ${UCSC_CHAIN_TERMS} See docs/44-ancient.md."
  fi
  echo "Installing the ancient populations panel, AADR ${AADR_VERSION} 1240K (step 44) in ${dir}"
  ensure_image EIGENSOFT_IMAGE 44 || return 1
  need_og "chooses the ancient panel's people and places its sites" || return 1
  # ANCIENT_SOURCE_DIR's unchecked files are staged in raw_<version> here,
  # never in the downloads folder, where they would pass for AADR's.
  if [ -n "$src" ]; then raw=$old; else raw=$(dl "aadr-${AADR_VERSION}"); fi
  mkdir -p "$raw"
  for e in anno ind snp geno; do
    f="${raw}/aadr.${e}"
    if [ -n "$src" ]; then
      cp "$(compgen -G "${src}/*.${e}" | head -n 1)" "$f" 2>/dev/null \
        || { echo "[WARN] ANCIENT_SOURCE_DIR=${src} has no .${e} file."; return 1; }
      continue
    fi
    # One URL per file, spelled out so scripts/ci/check-links.sh asks each pinned id.
    case "$e" in
      anno) url="https://dataverse.harvard.edu/api/access/datafile/${AADR_ANNO_ID}" md5=$AADR_ANNO_MD5 ;;
      ind) url="https://dataverse.harvard.edu/api/access/datafile/${AADR_IND_ID}" md5=$AADR_IND_MD5 ;;
      snp) url="https://dataverse.harvard.edu/api/access/datafile/${AADR_SNP_ID}" md5=$AADR_SNP_MD5 ;;
      geno) url="https://dataverse.harvard.edu/api/access/datafile/${AADR_GENO_ID}" md5=$AADR_GENO_MD5 ;;
    esac
    # An older setup downloaded into raw_<version>: a file there that
    # matches its md5 moves to the downloads folder (a staged copy does not).
    if [ -s "${old}/aadr.${e}" ] && [ "$(_digest md5 "${old}/aadr.${e}")" = "$md5" ]; then
      adopt "${old}/aadr.${e}" "$f" "$url" md5 "$md5"
    fi
    [ "$e" != geno ] || echo "  Downloading the AADR genotypes (7.1 GB; an interrupted download resumes)..."
    if ! fetch "$url" "$f" md5 "$md5"; then
      echo "[WARN] Could not download the AADR ${AADR_VERSION} .${e} file; run: $0 --ancient-panel ${GENOME_DIR}"
      return 1
    fi
  done
  if [ -n "${ANCIENT_CHAIN:-}" ]; then
    chain=$ANCIENT_CHAIN
  elif ! fetch "$HG19_TO_HG38_CHAIN_URL" "$chain" sha256 "$HG19_TO_HG38_CHAIN_SHA256"; then
    echo "[WARN] Could not download UCSC's hg19ToHg38 chain; run: $0 --ancient-panel ${GENOME_DIR}"
    return 1
  fi
  rm -rf "$part"
  mkdir -p "$part"
  local P R C
  P=$(cpath "$part") R=$(cpath "$raw")
  C=$(cpath "$chain") || { echo "[WARN] the chain file ${chain} must be inside GENOME_DIR."; return 1; }
  echo "  [1/5] Choosing the present-day reference people and the ancient people to project..."
  echo "  [2/5] Placing the 1240K sites on GRCh38..."
  # shellcheck disable=SC2016  # $1 to $4 belong to the inner bash
  if ! run_in --rw "$part" -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 4g "$OG_IMAGE" \
         bash -euo pipefail -c '
           /pgp-bin/og ancient select --anno "$1/aadr.anno" --ind "$1/aadr.ind" \
             --out-ind "$2/panel.ind" --poplist "$2/poplist.txt" --ancients "$2/ancients.tsv"
           /pgp-bin/og ancient lift --snp "$1/aadr.snp" --chain "$3" --fasta "$4" --out "$2/sites_grch38.tsv"' \
         _ "$R" "$P" "$C" "$REF_FASTA_C"; then
    echo "[WARN] Could not choose the panel's people or place its sites (the lines above say why)."
    return 1
  fi
  # The parameter files name container paths; EIGENSOFT reads one per run.
  printf '%s\n' "genotypename: ${R}/aadr.geno" "snpname: ${R}/aadr.snp" "indivname: ${P}/panel.ind" \
    "poplistname: ${P}/poplist.txt" "outputformat: EIGENSTRAT" "genooutfilename: ${P}/moderns.geno" \
    "snpoutfilename: ${P}/moderns.snp" "indoutfilename: ${P}/moderns.ind" "hashcheck: NO" > "${part}/convertf.par"
  ancient_smartpca_par "${R}/aadr.geno" "${R}/aadr.snp" "${P}/panel.ind" "${P}/poplist.txt" "${P}/panel" > "${part}/panel.par"
  ancient_smartpca_par "${P}/moderns.geno" "${P}/moderns.snp" "${P}/moderns.ind" "${P}/poplist.txt" "${P}/check" > "${part}/check.par"
  echo "  [3/5] The present-day people in EIGENSTRAT format (convertf)..."
  echo "  [4/5] The axes, with the ancient people projected (smartpca, lsqproject)..."
  echo "  [5/5] The same axes again from the present-day people alone, as step 44 runs..."
  # shellcheck disable=SC2016  # $1 belongs to the inner shell
  if ! run_in --rw "$part" --cpus "$THREADS" --memory 12g "$EIGENSOFT_IMAGE" sh -euc '
         convertf -p "$1/convertf.par" > "$1/convertf.log"
         smartpca -p "$1/panel.par" > "$1/panel.log"
         smartpca -p "$1/check.par" > "$1/check.log"' _ "$P"; then
    tail -n 20 "${part}"/*.log 2>/dev/null || true
    echo "[WARN] EIGENSOFT failed on the AADR files (logs in ${part}); run: $0 --ancient-panel ${GENOME_DIR}"
    return 1
  fi
  if ! run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 2g "$OG_IMAGE" \
         /pgp-bin/og ancient align --panel-evec "${P}/panel.evec" --run-evec "${P}/check.evec"; then
    echo "[WARN] The two smartpca runs on the same present-day people gave different axes (logs in ${part})."
    return 1
  fi
  rm -f "${part}/check.evec" "${part}/check.eval"
  printf 'release\t%s\nsites_grch38\t%s\nbuilt\t%s\n' "$AADR_VERSION" "$(($(wc -l < "${part}/sites_grch38.tsv") - 1))" \
    "$(date -u '+%Y-%m-%d')" > "${part}/DONE"
  rm -rf "$dir"  # not kept: the older panel, replaced
  mv "$part" "$dir"
  # The staged copy of ANCIENT_SOURCE_DIR is no download: it goes at once,
  # as does what an older setup left in raw_<version> (its AADR files were
  # moved to the downloads folder above).
  rm -rf "$old"  # not kept: the staged copy of ANCIENT_SOURCE_DIR, or an older setup's leftovers
  [ -n "$src" ] || done_with "$raw"
  echo "[OK] ancient populations panel AADR ${AADR_VERSION} (step 44): ${dir} ($(wc -l < "${dir}/moderns.ind" | tr -d ' ') present-day people, $(($(wc -l < "${dir}/ancients.tsv") - 1)) ancient people)"
}

# ancient_smartpca_par GENO SNP IND POPLIST OUT: the smartpca parameters of
# setup.sh and step 44 (scripts/44-ancient.sh writes the same lines): axes from
# the present-day groups of POPLIST, everyone else projected by least squares
# (lsqproject, for samples with missing sites), no outlier removal, so both
# runs build their axes from the same people.
ancient_smartpca_par() {
  printf '%s\n' "genotypename: $1" "snpname: $2" "indivname: $3" "poplistname: $4" \
    "evecoutname: $5.evec" "evaloutname: $5.eval" "lsqproject: YES" "numoutevec: 10" \
    "numoutlieriter: 0" "hashcheck: NO" "numthreads: ${THREADS}"
}

# --- Step 25 and 26: pgsc_calc, its scores and its ancestry panel ---------------------
PGS_BASE_URL=${PGS_BASE_URL:-https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores}

# install_pgsc_calc DIR: GitHub's archive of pgsc_calc PGSC_CALC_VERSION,
# checked against PGSC_CALC_SHA256, unpacked into DIR. Step 25 holds the same.
install_pgsc_calc() {
  local dir=$1 tgz url
  url=$(source_of pgsc_calc url) && tgz=$(source_of pgsc_calc path) || return 1
  adopt "${dir}.tar.gz" "$tgz" "$url" sha256 "$PGSC_CALC_SHA256"
  fetch "$url" "$tgz" sha256 "$PGSC_CALC_SHA256" || return 1
  rm -rf "${dir}.part"
  mkdir -p "${dir}.part"
  tar -xzf "$tgz" -C "${dir}.part" --strip-components 1 || { rm -rf "${dir}.part"; return 1; }
  done_with "$tgz"
  # An incomplete folder left by an earlier run would receive the new one inside it.
  rm -rf "$dir"  # not kept: an incomplete folder of an earlier run
  mv "${dir}.part" "$dir"
}

# fetch_prs_scores LICENCE: the GRCh38-harmonised file of every score of
# assets/pgs_scores.tsv whose licence column is LICENCE (an empty cell counts
# as author), each checked against the md5 the PGS Catalog publishes beside
# it. Exit 1 when one download failed.
fetch_prs_scores() {
  local id url dest rc=0 list="${PGP_ROOT}/assets/pgs_scores.tsv"
  if [ ! -r "$list" ]; then
    echo "[WARN] Could not read the PGS score list ${list}"
    return 1
  fi
  mkdir -p "${GENOME_DIR}/prs_scores" || return 1
  while read -r id; do
    dest="${GENOME_DIR}/prs_scores/${id}.txt.gz"
    [ -s "$dest" ] && continue
    url="${PGS_BASE_URL}/${id}/ScoringFiles/Harmonized/${id}_hmPOS_GRCh38.txt.gz"
    if ! fetch "$url" "$dest" md5 "${url}.md5"; then
      echo "[WARN] Could not download ${id} (step 25 tries again when it runs): ${url}"
      rc=1
    fi
  done < <(awk -F'\t' -v want="$1" '!h && $1 == "pgs_id" { for (i = 1; i <= NF; i++) c[$i] = i; h = 1; next }
    h && $1 ~ /^PGS[0-9]+$/ { l = c["licence"] ? $c["licence"] : ""; if (l == "") l = "author"; if (l == want) print $1 }' "$list")
  return "$rc"
}

# install_prs: the scores of assets/pgs_scores.tsv under the PGS Catalog's
# default terms (licence author; the four free for academic research use
# only come with --prs-academic), pgsc_calc PGSC_CALC_VERSION in tools/, and
# the nf-schema plugin its nextflow.config names. With them
# step 25 and the PRS process need no network. A failure is a warning: step
# 25 fetches what is missing the first time it runs.
install_prs() {
  local rc=0 calc="${GENOME_DIR}/tools/pgsc_calc-${PGSC_CALC_VERSION}"
  fetch_prs_scores author || rc=1
  [ "$rc" -eq 0 ] && echo "[OK] PGS Catalog scores (step 25): ${GENOME_DIR}/prs_scores/"
  if [ -f "${calc}/main.nf" ]; then
    echo "[OK] pgsc_calc ${PGSC_CALC_VERSION} (step 25) already present."
  elif install_pgsc_calc "$calc"; then
    echo "[OK] pgsc_calc ${PGSC_CALC_VERSION} (step 25): ${calc}"
  else
    echo "[WARN] Could not fetch pgsc_calc ${PGSC_CALC_VERSION}; step 25 tries again when it runs."
    rc=1
  fi
  if ! command -v nextflow >/dev/null 2>&1; then
    echo "[WARN] Nextflow is not installed: run-all.sh and step 25 need it (docs/nextflow.md)."
    rc=1
  elif [ -d "${NXF_HOME:-${HOME}/.nextflow}/plugins/nf-schema-${PGSC_CALC_NF_SCHEMA}" ]; then
    echo "[OK] Nextflow plugin nf-schema ${PGSC_CALC_NF_SCHEMA} (pgsc_calc) already present."
  elif nextflow plugin install "nf-schema@${PGSC_CALC_NF_SCHEMA}" >/dev/null; then
    echo "[OK] Nextflow plugin nf-schema ${PGSC_CALC_NF_SCHEMA} (pgsc_calc)."
  else
    echo "[WARN] Could not install the nf-schema ${PGSC_CALC_NF_SCHEMA} plugin; step 25 tries again when it runs."
    rc=1
  fi
  if [ "${SETUP_EVERYTHING_SCOPE:-}" = base ]; then return "$rc"; fi
  return 0
}

# install_ancestry_panel: pgsc_calc's reference panel (PGSC_PANEL, the 1000
# Genomes database; ANCESTRY_PANEL_NAME picks another file of the PGS
# Catalog's resources folder, the tests use its small synthetic one), checked
# against the md5 that folder lists, and beside it the panel's common GRCh38
# biallelic SNVs (chrN, position, REF, ALT), from the panel's own genotypes:
# step 25 genotypes them from the gVCF, so the projection counts the sites
# where the sample matches the reference.
install_ancestry_panel() {
  local name=${ANCESTRY_PANEL_NAME:-$PGSC_PANEL} dir="${GENOME_DIR}/reference/pgsc_calc" panel sites url md5s
  panel="${dir}/${name}.tar.zst"
  sites="${dir}/${name}_GRCh38_sites.tsv"
  if [ -s "$panel" ] && [ -s "$sites" ]; then
    echo "[OK] ancestry panel ${name} (steps 25 and 26) already present."
    return 0
  fi
  [ -s "$panel" ] || echo "Downloading pgsc_calc's ancestry panel ${name} (the 1000 Genomes panel is about 7 GB; an interrupted download resumes)..."
  # Row pgsc_ancestry_panel; ANCESTRY_PANEL_NAME is another file of its folder.
  url=$(source_of pgsc_ancestry_panel url) && md5s=$(source_of pgsc_ancestry_panel value) || return 1
  url="${url%/*}/${name}.tar.zst"
  # fetch stores a download only once its md5 matches, so a panel that is there is whole.
  if [ ! -s "$panel" ] && ! fetch "$url" "$panel" md5 "$md5s"; then
    echo "[WARN] Could not download the ancestry panel; run: $0 --ancestry-panel ${GENOME_DIR}"
    return 1
  fi
  # Only the SNVs pgsc_calc can use for the projection: its FILTER_VARIANTS
  # keeps biallelic ACGT SNVs on the autosomes with a panel frequency of 5%
  # or more (its maf_ref). The whole GRCh38 table holds about 62 million
  # SNVs; genotyping all of them from a gVCF would take step 25 hours.
  echo "  Listing the panel's common GRCh38 SNVs (plink2 on the panel's genotypes)..."
  local tmp="${dir}/.${name}.sites.tmp" prefix
  rm -rf "$tmp"
  mkdir -p "$tmp"
  # shellcheck disable=SC2016  # $1 and $2 belong to the inner sh
  if ! run_in --rw "$dir" "$PGSC_ZSTD_IMAGE" sh -c 'cd "$2" && tar -xf "$1" --wildcards "GRCh38_*_ALL.pgen" "GRCh38_*_ALL.psam" "GRCh38_*_ALL.pvar.zst"' \
      _ "$(cpath "$panel")" "$(cpath "$tmp")"; then
    rm -rf "$tmp"
    echo "[WARN] Could not unpack the GRCh38 genotypes of ${panel}; run: $0 --ancestry-panel ${GENOME_DIR}"
    return 1
  fi
  prefix=$(find "$tmp" -maxdepth 1 -name 'GRCh38_*_ALL.pgen' | head -1)
  prefix=${prefix%.pgen}
  # One chromosome at a time: plink2 then loads a 22nd of the panel's variant
  # table (the whole table, 61.6 million variants, did not fit a 16 GB runner).
  # --allow-extra-chr: the panel names contigs plink2 does not know, as
  # pgsc_calc's own plink2 calls on it allow.
  # plink2 still indexes the whole table: 4 GB of workspace ran out of memory
  # even for chromosome 1. 8 GB, the most pgsc_calc's own call on the panel
  # gets here, or three quarters of a smaller machine's RAM.
  local c ok=true mem_mb
  if [ -r /proc/meminfo ]; then
    mem_mb=$(awk '/^MemTotal:/ {print int($2 * 3 / 4 / 1024)}' /proc/meminfo)
  else
    mem_mb=$(( $(sysctl -n hw.memsize 2>/dev/null || echo 17179869184) * 3 / 4 / 1048576 ))
  fi
  [ "$mem_mb" -le 8000 ] || mem_mb=8000
  [ -n "$prefix" ] || ok=false
  for c in $(seq 1 22); do
    $ok || break
    # A panel without a common SNV on a chromosome ("No variants remaining") is not an error.
    if ! run_in --rw "$dir" "$PLINK2_IMAGE" plink2 --pfile "$(cpath "$prefix")" vzs --chr "$c" \
          --allow-extra-chr --snps-only just-acgt --max-alleles 2 --maf 0.05 \
          --make-just-pvar --threads "$THREADS" --memory "$mem_mb" --out "$(cpath "${tmp}/common_${c}")" > "${tmp}/plink2.log" 2>&1 \
       && ! grep -q 'No variants remaining' "${tmp}/plink2.log"; then
      head -n 40 "${tmp}/plink2.log"
      ok=false
    fi
  done
  if ! $ok || ! awk -F'\t' -v OFS='\t' '
          /^##/ { next }
          /^#/ { for (i = 1; i <= NF; i++) c[$i] = i; next }
          {
            chr = $c["#CHROM"]; sub(/^chr/, "", chr)
            r = $c["REF"]; a = $c["ALT"]
            if (r ~ /^[ACGT]$/ && a ~ /^[ACGT]$/) print "chr" chr, $c["POS"], r, a
          }' "${tmp}"/common_*.pvar | LC_ALL=C sort -u -k1,1 -k2,2n -k3,3 -k4,4 > "${sites}.part" || [ ! -s "${sites}.part" ]; then
    rm -rf "$tmp" "${sites}.part"
    echo "[WARN] Could not list the common GRCh38 SNVs of ${panel}; run: $0 --ancestry-panel ${GENOME_DIR}"
    return 1
  fi
  rm -rf "$tmp"
  mv "${sites}.part" "$sites"
  echo "[OK] ancestry panel ${name} (steps 25 and 26): ${panel} ($(wc -l < "$sites" | tr -d ' ') common GRCh38 SNVs in ${sites##*/})"
}

# --- Step 26b: the fine-grained ancestry panel (opt-in) -----------------------------
# gnomAD's HGDP + 1000 Genomes callset, phased with SHAPEIT5 (Koenig et al.
# 2024): 4,091 people of 80 populations, open access. setup.sh streams the
# common biallelic SNVs of chr1 to chr22 from gnomAD's public bucket (about
# 19 GB is read; only the kept SNVs are stored) and keeps one per 2 kb.
# The 22 files are row ancestry_fine_bcf of assets/data-sources.tsv
# (scripts/ci/check-links.sh asks each one).
ANCESTRY_FINE_META_URL=https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc/gnomad_meta_updated.tsv
# Beagle's GRCh38 PLINK genetic map (HapMap, lifted from GRCh37), for phasing and painting.
# Beagle's server is its only host and was unreachable for hours on 2026-10-15;
# the Internet Archive's copy is byte for byte the same file (same sha256), so
# the install falls back to it, checked against the same pin.
ANCESTRY_FINE_MAP_URL=https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
ANCESTRY_FINE_MAP_MIRROR_URL=https://web.archive.org/web/20260306143720id_/https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip

# ancestry_fine_pins: one line per chromosome file of the panel: chrN, the
# crc32c Google reports for it and its size, from a HEAD request each. Their
# sha256 is ANCESTRY_FINE_PANEL_SHA256, so a file the bucket replaces changes
# it. bcftools reads the panel a region at a time, so the bytes cannot be
# checked as they stream; the bucket's own checksum of each file can.
ancestry_fine_pins() {
  local c u
  command -v curl >/dev/null 2>&1 || { echo "ERROR: --ancestry-fine-panel needs curl to check the panel's files" >&2; return 1; }
  for c in $(seq 1 22); do
    u=$(source_of ancestry_fine_bcf --part "$c" url) || return 1
    curl -fsI --connect-timeout 20 --max-time 60 --retry 3 --retry-delay 5 "$u" | tr -d '\r' \
      | awk -v c="chr${c}" 'tolower($1) == "x-goog-hash:" && $2 ~ /^crc32c=/ {h = substr($2, 8)}
                            tolower($1) == "content-length:" {n = $2}
                            END {if (h == "" || n == "") exit 1; print c "\t" h "\t" n}' \
      || { echo "ERROR: no crc32c or size for the panel's chr${c} file (${u})" >&2; return 1; }
  done
}

# install_fine_maps DIR CHROMS: DIR/maps/ of step 26b, the genetic map of each
# chromosome of CHROMS (space separated) that Beagle and og's painting read,
# cut from DIR/maps/plink.GRCh38.map.zip by og (`og setup genetic-map`), which
# writes each file whole or not at all.
# maps/DONE, written last, names the map and the chromosomes: without it (an
# install from before og, whose maps a kill could have left cut short), or
# with a map file missing, every map is written again. The chrN.rfmix.map
# files an older install wrote for RFMix stay on disk, read by nothing.
# tests/test_setup_fine_maps.sh runs it.
install_fine_maps() {
  local dir=$1 chroms=$2
  local gstamp="map=${ANCESTRY_FINE_MAP_SHA256} chroms=${chroms% }" c whole=yes
  for c in $chroms; do
    [ -s "${dir}/maps/${c}.map" ] || whole=no
  done
  [ "$whole" = yes ] && [ "$(cat "${dir}/maps/DONE" 2>/dev/null)" = "$gstamp" ] && return 0
  rm -f "${dir}/maps/DONE"
  need_og "unpacks the genetic map" || return 1
  # SC2086: one argument per chromosome.
  # shellcheck disable=SC2086
  if ! run_in --rw "$dir" -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$OG_IMAGE" /pgp-bin/og setup genetic-map \
      --zip "$(cpath "${dir}/maps/plink.GRCh38.map.zip")" --out "$(cpath "${dir}/maps")" $chroms; then
    echo "[WARN] Could not unpack the genetic map"
    return 1
  fi
  echo "$gstamp" > "${dir}/maps/DONE"
}

# install_ancestry_fine_panel: the reference of step 26b, in
# reference/ancestry_fine/ANCESTRY_FINE_PANEL (ANCESTRY_FINE_DIR: another
# folder inside GENOME_DIR). Four parts, each redone only when its inputs change:
#   samples.tsv  every person of the panel with population, continental and
#                regional group, reference and painting flags (og ancestry-fine
#                panel-map, rust/crates/og-ancestry); ANCESTRY_FINE_SAMPLE_MAP gives
#                a table of your own instead;
#   panel/       chrN.vcf.gz: their phased genotypes at the biallelic SNVs with a
#                minor allele frequency of 5% or more, one per 2 kb (bcftools
#                +prune); sites.tsv (chrN, position, REF, ALT) and positions.tsv;
#   maps/        chrN.map (Beagle and og's painting);
#   model/       ref.bim, ref.fam and ref.geno: the reference individuals' genotypes
#                at the SNVs left after LD pruning (r2 above 0.1 within 500 kb
#                removed); for each level (continent, region) each group's allele
#                frequencies, <level>.P, counted by og ancestry-fine freqs, with
#                their column order in levels.tsv; and paint_map.tsv, the
#                reference individuals og paints with.
# ANCESTRY_FINE_REGIONS (a BED file) keeps only those regions, and
# ANCESTRY_FINE_EXCLUDE (sample ids) leaves people out of the panel: the e2e
# case uses both, for the fixture's slices and the two people of its
# synthetic mix. ANCESTRY_FINE_THIN (default 2kb) is the window that keeps
# one SNV (the e2e case keeps more of its few megabases with 200bp), and
# ANCESTRY_FINE_PAINT_PER_GROUP (default 200, at least 10) caps the painting
# references of each continental group.
install_ancestry_fine_panel() {
  local dir=${ANCESTRY_FINE_DIR:-${GENOME_DIR}/reference/ancestry_fine/${ANCESTRY_FINE_PANEL}}
  local regions=${ANCESTRY_FINE_REGIONS:-} exclude=${ANCESTRY_FINE_EXCLUDE:-} smap=${ANCESTRY_FINE_SAMPLE_MAP:-}
  local thin=${ANCESTRY_FINE_THIN:-2kb}
  local v c e k m fit rstamp="" pstamp mstamp tmp chroms="" src bcf net chr22 folder
  cpath "$dir" >/dev/null || return 1
  chr22=$(source_of ancestry_fine_bcf --part 22 url) || return 1
  folder=${chr22%/*}
  echo "Installing the fine-grained ancestry panel ${ANCESTRY_FINE_PANEL} (step 26b) in ${dir}"
  echo "  Step 26b finds the fractions and paints the chromosomes with og."
  for v in OG_IMAGE BCFTOOLS_HTTPS_IMAGE BEAGLE_IMAGE; do ensure_image "$v" 26b || return 1; done
  # The sample table and the column order come from og (og ancestry-fine).
  if ! "${PGP_ROOT}/rust/og.sh" check --required > /dev/null 2>&1 && ! "${PGP_ROOT}/rust/og.sh" install; then
    echo "[WARN] Could not build og, which writes the panel's sample table: run ${0} --og and look at its error."
    return 1
  fi
  mkdir -p "${dir}/panel" "${dir}/maps" "${dir}/model"
  echo "  Checking the panel's files against ANCESTRY_FINE_PANEL_SHA256..."
  if ! ancestry_fine_pins > "${dir}/pins.tsv"; then
    echo "[WARN] Could not read the panel's files at ${folder}; run: $0 --ancestry-fine-panel ${GENOME_DIR}"
    return 1
  fi
  if [ "$(_digest sha256 "${dir}/pins.tsv")" != "$ANCESTRY_FINE_PANEL_SHA256" ]; then
    echo "[WARN] The panel's files at ${folder} are not the ones ANCESTRY_FINE_PANEL_SHA256 pins (see ${dir}/pins.tsv)."
    echo "       The bucket replaced a file: read what changed before moving the pin in versions.env."
    return 1
  fi
  # The metadata and the map, unless they are there with their pinned sha256;
  # the map's mirror (the fourth field) is tried when its own host fails.
  for v in "${ANCESTRY_FINE_META_URL}|${dir}/meta.tsv|${ANCESTRY_FINE_META_SHA256}|" \
           "${ANCESTRY_FINE_MAP_URL}|${dir}/maps/plink.GRCh38.map.zip|${ANCESTRY_FINE_MAP_SHA256}|${ANCESTRY_FINE_MAP_MIRROR_URL}"; do
    IFS='|' read -r c e k m <<< "$v"
    [ -s "$e" ] && [ "$(_digest sha256 "$e")" = "$k" ] && continue
    fetch "$c" "$e" sha256 "$k" && continue
    if [ -n "$m" ]; then
      rm -f "${e}.part"
      echo "  Trying the mirror of ${c##*/}: ${m}"
      fetch "$m" "$e" sha256 "$k" && continue
    fi
    echo "[WARN] Could not download ${c}; run: $0 --ancestry-fine-panel ${GENOME_DIR}"
    return 1
  done
  if [ -n "$regions" ]; then
    [ -s "$regions" ] || { echo "[WARN] ANCESTRY_FINE_REGIONS=${regions} is not a file"; return 1; }
    cp "$regions" "${dir}/regions.bed"
    rstamp=$(_digest sha256 "${dir}/regions.bed")
    chroms=$(awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {sub(/^chr/, "", $1); print $1}' "${dir}/regions.bed" | sort -un | sed 's/^/chr/' | tr '\n' ' ')
  else
    rm -f "${dir}/regions.bed"
    for c in $(seq 1 22); do chroms+="chr${c} "; done
  fi

  # samples.tsv: from the metadata, or the table given.
  # net: slice, the header only: the sample list of a panel file in gnomAD's
  # bucket, saved as bcf_samples.txt. SC2016: $1 belongs to the inner sh.
  # shellcheck disable=SC2016
  if ! run_in --net --rw "$dir" "$BCFTOOLS_HTTPS_IMAGE" sh -c 'bcftools query -l "$1" 2>/dev/null' \
        _ "$chr22" > "${dir}/bcf_samples.txt.part" \
     || [ ! -s "${dir}/bcf_samples.txt.part" ]; then
    rm -f "${dir}/bcf_samples.txt.part"
    echo "[WARN] Could not read the panel's sample list; run: $0 --ancestry-fine-panel ${GENOME_DIR}"
    return 1
  fi
  mv "${dir}/bcf_samples.txt.part" "${dir}/bcf_samples.txt"
  if [ -n "$smap" ]; then
    head -n 1 "$smap" | grep -q $'^sample\tpopulation\tcontinent\tregion\treference\tpaint$' \
      || { echo "[WARN] ANCESTRY_FINE_SAMPLE_MAP=${smap} does not start with the header of samples.tsv"; return 1; }
    cp "$smap" "${dir}/samples.tsv.part"
    if ! run_in --rw "$dir" -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$OG_IMAGE" \
          /pgp-bin/og ancestry-fine check-map --map "$(cpath "${dir}/samples.tsv.part")"; then
      rm -f "${dir}/samples.tsv.part"
      echo "[WARN] ANCESTRY_FINE_SAMPLE_MAP=${smap} has a group too small to model"
      return 1
    fi
  elif ! run_in --rw "$dir" -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$OG_IMAGE" \
        /pgp-bin/og ancestry-fine panel-map --meta "$(cpath "${dir}/meta.tsv")" \
          --bcf-samples "$(cpath "${dir}/bcf_samples.txt")" --exclude "$exclude" \
          --paint-per-group "${ANCESTRY_FINE_PAINT_PER_GROUP:-200}" --out "$(cpath "${dir}/samples.tsv.part")"; then
    rm -f "${dir}/samples.tsv.part"
    echo "[WARN] Could not write the panel's sample table"
    return 1
  fi
  mv "${dir}/samples.tsv.part" "${dir}/samples.tsv"
  awk -F'\t' 'NR > 1 {print $1}' "${dir}/samples.tsv" > "${dir}/panel/keep.txt"

  # panel/: redone when the pins, the regions or the people kept change; a
  # chromosome done by an interrupted run with the same inputs is kept.
  pstamp="pins=${ANCESTRY_FINE_PANEL_SHA256} regions=${rstamp:-all} thin=${thin} keep=$(_digest sha256 "${dir}/panel/keep.txt") bcftools=${BCFTOOLS_HTTPS_IMAGE}"
  if [ "$(cat "${dir}/panel/DONE" 2>/dev/null)" != "$pstamp" ]; then
    rm -f "${dir}/panel/DONE" "${dir}/model/DONE"
    for c in $chroms; do
      [ ! -s "${dir}/panel/${c}.vcf.gz.tbi" ] || [ "$(cat "${dir}/panel/${c}.stamp" 2>/dev/null)" != "$pstamp" ] || continue
      # A whole panel: each chromosome's file is downloaded to the downloads
      # folder, checked by the crc32c pins.tsv holds, and read there (kept,
      # done_with). With ANCESTRY_FINE_REGIONS only those regions are read
      # over https (htslib stores the remote index in the working directory,
      # so the container works in panel/ and deletes it after).
      src=$(source_of ancestry_fine_bcf --part "${c#chr}" url) || return 1
      net=(--net) bcf=""
      if [ -z "$rstamp" ]; then
        bcf=$(source_of ancestry_fine_bcf --part "${c#chr}" path) || return 1
        if ! fetch "$src" "$bcf" crc32c "$(awk -F'\t' -v c="$c" '$1 == c {print $2}' "${dir}/pins.tsv")"; then
          echo "[WARN] Could not download ${c} of the panel; run: $0 --ancestry-fine-panel ${GENOME_DIR}"
          return 1
        fi
        src=$(cpath "$bcf") net=()
      fi
      echo "  ${c}: reading the panel's common SNVs..."
      # bcftools warns about the files' PP header line on every read; those
      # lines go.
      # net: slice, the regions of ANCESTRY_FINE_REGIONS only (a whole panel is read from the downloads folder).
      # shellcheck disable=SC2016  # $1 to $5 belong to the inner bash
      if ! run_in ${net[@]+"${net[@]}"} --rw "$dir" "$BCFTOOLS_HTTPS_IMAGE" bash -euo pipefail -c '
            cd "$(dirname "$3")"; trap "rm -f ./*.bcf.csi" EXIT
            r=(); [ -z "$4" ] || r=(-R "$4" --regions-overlap pos)
            { bcftools view "${r[@]}" --force-samples -S "$2" -m2 -M2 -v snps -Ou "$1" \
                | bcftools view -q 0.05:minor -Ou \
                | bcftools annotate -x ID,INFO,^FORMAT/GT -Ou \
                | bcftools +prune -n 1 -w "$5" -N maxAF -Oz -o "$3.part"; } 2> >(grep -v bcf_hdr_check_sanity >&2)
            bcftools index -t "$3.part"
            mv "$3.part.tbi" "$3.tbi"; mv "$3.part" "$3"' \
            _ "$src" \
            "$(cpath "${dir}/panel/keep.txt")" "$(cpath "${dir}/panel/${c}.vcf.gz")" \
            "$([ -z "$rstamp" ] || cpath "${dir}/regions.bed")" "$thin"; then
        echo "[WARN] Could not read ${c} of the panel; run: $0 --ancestry-fine-panel ${GENOME_DIR}"
        return 1
      fi
      echo "$pstamp" > "${dir}/panel/${c}.stamp"
      [ -z "$bcf" ] || done_with "$bcf"
    done
    # sites.tsv: chrN, position, REF and ALT of every SNV; positions.tsv: the
    # first two columns, the form bcftools -R and -T read.
    # shellcheck disable=SC2016  # $1 and $2 belong to the inner bash
    if ! run_in --rw "$dir" "$BCFTOOLS_IMAGE" bash -euo pipefail -c \
          'cd "$1"; for c in $2; do bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" "$c.vcf.gz"; done' \
          _ "$(cpath "${dir}/panel")" "$chroms" > "${dir}/panel/sites.tsv.part" || [ ! -s "${dir}/panel/sites.tsv.part" ]; then
      echo "[WARN] Could not list the panel's sites"
      return 1
    fi
    cut -f1,2 "${dir}/panel/sites.tsv.part" > "${dir}/panel/positions.tsv"
    mv "${dir}/panel/sites.tsv.part" "${dir}/panel/sites.tsv"
    echo "$pstamp" > "${dir}/panel/DONE"
  fi
  echo "  panel: $(wc -l < "${dir}/panel/sites.tsv" | tr -d ' ') SNVs of $(wc -l < "${dir}/panel/keep.txt" | tr -d ' ') people"

  install_fine_maps "$dir" "$chroms" || return 1

  # model/: redone when samples.tsv or the panel changes. Bump freqs=og-N
  # when the output of og ancestry-fine freqs or og_admix::geno (ref.geno) changes.
  mstamp="panel=$(_digest sha256 "${dir}/panel/sites.tsv") samples=$(_digest sha256 "${dir}/samples.tsv") freqs=og-1"
  if [ "$(cat "${dir}/model/DONE" 2>/dev/null)" != "$mstamp" ]; then
    tmp="${dir}/model/tmp"
    rm -rf "$tmp" "${dir}/model/DONE"
    mkdir -p "$tmp"
    awk -F'\t' 'NR > 1 && $5 == 1 {print $1}' "${dir}/samples.tsv" > "${tmp}/ref.txt"
    echo "  model: LD pruning on $(wc -l < "${tmp}/ref.txt" | tr -d ' ') reference individuals, then each group's allele frequencies..."
    # shellcheck disable=SC2016  # $1 to $4 belong to the inner bash
    if ! run_in --rw "$dir" "$BCFTOOLS_IMAGE" bash -euo pipefail -c '
          for c in $4; do bcftools view -S "$2" -Ou "$1/$c.vcf.gz" | bcftools +prune -m 0.1 -w 500kb -Ob -o "$3/$c.bcf"; done
          for c in $4; do echo "$3/$c.bcf"; done > "$3/bcfs.txt"
          bcftools concat -f "$3/bcfs.txt" -Oz -o "$3/adm.vcf.gz"' \
          _ "$(cpath "${dir}/panel")" "$(cpath "${tmp}/ref.txt")" "$(cpath "$tmp")" "$chroms" \
       || ! fit=$(run_in --rw "$dir" -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$OG_IMAGE" /pgp-bin/og ancestry-fine freqs \
              --vcf "$(cpath "${tmp}/adm.vcf.gz")" --samples "$(cpath "${dir}/samples.tsv")" --out "$(cpath "$tmp")"); then
      echo "[WARN] Could not count the reference groups' allele frequencies"
      return 1
    fi
    echo "  model: ${fit}"
    awk -F'\t' -v OFS='\t' 'NR > 1 && $6 == 1 {print $1, $3}' "${dir}/samples.tsv" > "${tmp}/paint_map.tsv"
    for e in ref.bim ref.fam ref.geno continent.P region.P levels.tsv paint_map.tsv; do mv -f "${tmp}/${e}" "${dir}/model/${e}"; done
    # ref.bed: the model's genotypes as an older install wrote them, read by nothing now.
    rm -rf "$tmp" "${dir}/model/ref.bed"
    echo "$mstamp" > "${dir}/model/DONE"
  fi
  echo "[OK] fine-grained ancestry panel ${ANCESTRY_FINE_PANEL} (step 26b): ${dir} ($(wc -l < "${dir}/model/ref.bim" | tr -d ' ') SNVs in the model, $(wc -l < "${dir}/model/paint_map.tsv" | tr -d ' ') painting references)"
}

# install_mocha_data: the phasing panel of step 51 in
# reference/mocha-MOCHA_PANEL_VERSION (MOCHA_DATA: another folder inside
# GENOME_DIR). For each of chr1 to chr22 the 1000 Genomes high-coverage phased
# panel's VCF is downloaded, checked against the md5 of the release's manifest
# (itself checked against MOCHA_PANEL_MANIFEST_SHA256), cut to its biallelic
# SNVs seen at least twice with their genotypes only (panel/chrN.vcf.gz, what
# Beagle reads), and kept in the downloads folder. From it, its common germline
# deletions and duplications go to panel/chrN.cnp.bed (the MoChA README's
# GRCh38 recipe: AC>1, over 10 kb, SVTYPE CNV, DEL or DUP), joined into
# cnps.bed, the --cnp list of step 51. maps/chrN.map is the GRCh38 genetic map
# of step 26b. About 30 GB is read; a chromosome done is kept by a rerun. MOCHA_REGIONS (a BED file) keeps only those regions and the
# chromosomes they name: the e2e case uses it for the fixture's chr20 slice.
install_mocha_data() {
  local dir=${MOCHA_DATA:-${GENOME_DIR}/reference/mocha-${MOCHA_PANEL_VERSION}}
  local regions=${MOCHA_REGIONS:-} chroms="" c e f md5 stamp rstamp=all
  cpath "$dir" >/dev/null || return 1
  echo "Installing MoChA's phasing panel ${MOCHA_PANEL_VERSION} (step 51) in ${dir}"
  for c in MOCHA_IMAGE BEAGLE_IMAGE; do ensure_image "$c" 51 || return 1; done
  local raw x row url
  raw=$(dl "mocha-${MOCHA_PANEL_VERSION}")
  mkdir -p "${dir}/panel" "${dir}/maps" "$raw"
  # An older setup downloaded into raw/ here: those files move to the
  # downloads folder, where fetch checks them.
  for x in "${dir}/raw/"*; do
    [ ! -e "$x" ] || adopt "$x" "${raw}/${x##*/}"
  done
  rmdir "${dir}/raw" 2>/dev/null || true
  if ! { [ -s "${dir}/manifest.txt" ] && [ "$(_digest sha256 "${dir}/manifest.txt")" = "$MOCHA_PANEL_MANIFEST_SHA256" ]; } \
     && ! fetch "${MOCHA_PANEL_URL}/20220804_manifest.txt" "${dir}/manifest.txt" sha256 "$MOCHA_PANEL_MANIFEST_SHA256"; then
    echo "[WARN] Could not read the panel's manifest (${MOCHA_PANEL_URL}); run: $0 --mocha-data ${GENOME_DIR}"
    return 1
  fi
  if [ -n "$regions" ]; then
    [ -s "$regions" ] || { echo "[WARN] MOCHA_REGIONS=${regions} is not a file"; return 1; }
    cp "$regions" "${dir}/regions.bed"
    rstamp=$(_digest sha256 "${dir}/regions.bed")
    chroms=$(awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {sub(/^chr/, "", $1); print $1}' "${dir}/regions.bed" | sort -un | sed 's/^/chr/' | tr '\n' ' ')
  else
    rm -f "${dir}/regions.bed"
    for c in $(seq 1 22); do chroms+="chr${c} "; done
  fi
  stamp="manifest=${MOCHA_PANEL_MANIFEST_SHA256} regions=${rstamp} bcftools=${MOCHA_IMAGE} cnp=1"
  [ "$(cat "${dir}/DONE" 2>/dev/null)" = "$stamp" ] && { echo "[OK] MoChA's phasing panel (step 51) already installed: ${dir}"; return 0; }
  rm -f "${dir}/DONE"
  for c in $chroms; do
    [ -s "${dir}/panel/${c}.vcf.gz.tbi" ] && [ -e "${dir}/panel/${c}.cnp.bed" ] \
      && [ "$(cat "${dir}/panel/${c}.stamp" 2>/dev/null)" = "$stamp" ] && continue
    f=$(source_of mocha_panel --part "${c#chr}" path) || return 1
    f=${f##*/}
    for e in "" .tbi; do
      row=mocha_panel; [ -z "$e" ] || row=mocha_panel_index
      url=$(source_of "$row" --part "${c#chr}" url) || return 1
      md5=$(awk -v n="${f}${e}" '$1 == n {print $2}' "${dir}/manifest.txt")
      [ -n "$md5" ] || { echo "[WARN] ${f}${e} is not in the panel's manifest"; return 1; }
      [ -s "${raw}/${f}${e}" ] && [ "$(_digest md5 "${raw}/${f}${e}")" = "$md5" ] && continue
      fetch "$url" "${raw}/${f}${e}" md5 "$md5" \
        || { echo "[WARN] Could not download ${f}${e}; run: $0 --mocha-data ${GENOME_DIR}"; return 1; }
    done
    echo "  ${c}: the biallelic SNVs seen twice or more, genotypes only, and the common CNVs..."
    # shellcheck disable=SC2016  # $1 to $4 belong to the inner bash
    if ! run_in --rw "$dir" "$MOCHA_IMAGE" bash -euo pipefail -c '
          r=(); [ -z "$3" ] || r=(-R "$3")
          bcftools query "${r[@]}" -i "AC>1 && END-POS+1>10000 && (SVTYPE==\"CNV\" || SVTYPE==\"DEL\" || SVTYPE==\"DUP\")" \
            -f "%CHROM\t%POS0\t%END\t%SVTYPE\n" "$1" > "$4.part"
          bcftools view --no-version "${r[@]}" -m2 -M2 -v snps -c 2 -Ou "$1" \
            | bcftools annotate --no-version -x ID,QUAL,FILTER,INFO,^FMT/GT -Oz -o "$2.part"
          bcftools index -f -t "$2.part"
          mv "$4.part" "$4"; mv "$2.part.tbi" "$2.tbi"; mv "$2.part" "$2"' \
          _ "$(cpath "${raw}/${f}")" "$(cpath "${dir}/panel/${c}.vcf.gz")" "$([ -z "$regions" ] || cpath "${dir}/regions.bed")" \
          "$(cpath "${dir}/panel/${c}.cnp.bed")"; then
      echo "[WARN] Could not cut ${c} of the panel"
      return 1
    fi
    echo "$stamp" > "${dir}/panel/${c}.stamp"
    done_with "${raw}/${f}" "${raw}/${f}.tbi"
  done
  # The genetic map: the PLINK GRCh38 map of step 26b, from its host or the mirror.
  local map="${dir}/maps/plink.GRCh38.map.zip"
  if ! { [ -s "$map" ] && [ "$(_digest sha256 "$map")" = "$ANCESTRY_FINE_MAP_SHA256" ]; }; then
    fetch "$ANCESTRY_FINE_MAP_URL" "$map" sha256 "$ANCESTRY_FINE_MAP_SHA256" \
      || { rm -f "${map}.part"; echo "  Trying the mirror: ${ANCESTRY_FINE_MAP_MIRROR_URL}"
           fetch "$ANCESTRY_FINE_MAP_MIRROR_URL" "$map" sha256 "$ANCESTRY_FINE_MAP_SHA256"; } \
      || { echo "[WARN] Could not download the genetic map; run: $0 --mocha-data ${GENOME_DIR}"; return 1; }
  fi
  # og writes each map whole or not at all. SC2086: one argument per chromosome.
  need_og "unpacks the genetic map" || return 1
  # shellcheck disable=SC2086
  if ! run_in --rw "$dir" -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$OG_IMAGE" /pgp-bin/og setup genetic-map \
      --zip "$(cpath "$map")" --out "$(cpath "${dir}/maps")" $chroms; then
    echo "[WARN] Could not unpack the genetic map"
    return 1
  fi
  for c in $chroms; do cat "${dir}/panel/${c}.cnp.bed"; done > "${dir}/cnps.bed.part" && mv "${dir}/cnps.bed.part" "${dir}/cnps.bed"
  echo "$stamp" > "${dir}/DONE"
  echo "[OK] MoChA's phasing panel (step 51, ${chroms% }): ${dir} ($(du -sh "$dir" | cut -f1), $(wc -l < "${dir}/cnps.bed" | tr -d ' ') common CNVs)"
}

# install_locityper_data: the database of step 53 in
# reference/locityper-LOCITYPER_PANEL_VERSION (LOCITYPER_DATA: another folder
# inside GENOME_DIR), as the Locityper documentation builds one from a
# pangenome VCF (locityper.vercel.app, "Preparing targets"):
#   counts.jf    jellyfish's canonical 25-mer counts of REF_FASTA, the
#                documentation's command (counts of 2 or more); rebuilt when
#                the reference's .fai changes
#   panel.vcf.gz the HPRC VCF of LOCITYPER_PANEL_URL, read over https 250 kb
#                either side of each locus through its .tbi, whose sha256
#                LOCITYPER_PANEL_TBI_SHA256 pins (Locityper extends a locus by
#                up to 200 kb to clear the variants on its edges)
#   db/          `locityper target` on those loci (GRCh38 is the panel's
#                reference haplotype), then `locityper augment`: the pairwise
#                alignments of the haplotypes, from which Locityper's weighted
#                distance is an edit distance
#   loci.bed     the loci it was built for, without comments
# LOCITYPER_LOCI (a BED file of chrom, start, end, name) replaces
# assets/locityper_loci.bed, for a smaller set of loci.
# DONE names the inputs, so a changed one builds again.
install_locityper_data() {
  local dir=${LOCITYPER_DATA:-${GENOME_DIR}/reference/locityper-${LOCITYPER_PANEL_VERSION}}
  local loci=${LOCITYPER_LOCI:-${PGP_ROOT}/assets/locityper_loci.bed} stamp jstamp n regions bad name hsize key saved
  case "$dir" in /*) ;; *) dir="${GENOME_DIR}/${dir}" ;; esac
  cpath "$dir" >/dev/null || return 1
  if [ ! -s "$REF_FASTA" ] || [ ! -s "${REF_FASTA}.fai" ]; then
    echo "[WARN] Step 53's database is built from the reference, and ${REF_FASTA} (with its .fai) is not there yet: run ./scripts/setup.sh ${GENOME_DIR} first."
    return 1
  fi
  [ -s "$loci" ] || { echo "[WARN] ${loci} is not a file"; return 1; }
  echo "Building Locityper's database (step 53, panel ${LOCITYPER_PANEL_VERSION}) in ${dir}"
  ensure_image LOCITYPER_IMAGE 53 || return 1
  mkdir -p "$dir"
  # Each row: chrom, start, end, name, with start < end inside a contig of the reference's .fai.
  bad=$(awk -F'\t' 'NR == FNR {len[$1] = $2; next}
    !/^#/ && NF && (NF < 4 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $4 == "" || $2 + 0 >= $3 + 0 || !($1 in len) || $3 + 0 > len[$1] + 0) {
      print FILENAME ":" FNR; exit}' "${REF_FASTA}.fai" "$loci")
  [ -z "$bad" ] || { echo "[WARN] ${bad}: want chrom, start, end, name (tab separated), start below end, inside a contig of ${REF_FASTA}"; return 1; }
  awk -F'\t' -v OFS='\t' '!/^#/ && NF {print $1, $2, $3, $4}' "$loci" > "${dir}/loci.bed.new"
  n=$(wc -l < "${dir}/loci.bed.new" | tr -d ' ')
  [ "$n" -gt 0 ] || { rm -f "${dir}/loci.bed.new"; echo "[WARN] ${loci} lists no locus"; return 1; }
  # The reference's .fai and its sequence (a new FASTA with the same contigs is a new reference too).
  jstamp="fai=$(_digest sha256 "${REF_FASTA}.fai") fasta=$(_digest sha256 "$REF_FASTA") image=${LOCITYPER_IMAGE}"
  stamp="${jstamp} panel=${LOCITYPER_PANEL_TBI_SHA256} loci=$(_digest sha256 "${dir}/loci.bed.new")"
  if [ "$(cat "${dir}/DONE" 2>/dev/null)" = "$stamp" ] && [ -d "${dir}/db" ] && [ -s "${dir}/counts.jf" ] && [ -s "${dir}/loci.bed" ]; then
    rm -f "${dir}/loci.bed.new"
    echo "[OK] Locityper's database (step 53, ${n} loci) already built: ${dir}"
    return 0
  fi
  rm -f "${dir}/DONE"

  if [ ! -s "${dir}/counts.jf" ] || [ "$(cat "${dir}/counts.stamp" 2>/dev/null)" != "$jstamp" ]; then
    echo "  k-mer counts of $(basename "$REF_FASTA") (jellyfish, several GB of RAM for GRCh38)..."
    rm -f "${dir}/counts.jf" "${dir}/counts.stamp"
    # The hash holds one entry per distinct k-mer, at most one per base, so it
    # is sized to the reference's length (about 3.1G for GRCh38), never below
    # 10M. A fixed 3G made a 235 Mb reference take 36 minutes on 16 GB.
    hsize=$(awk '{n += $2} END {n = n < 10000000 ? 10000000 : n; printf "%.0f", n}' "${REF_FASTA}.fai")
    if ! run_in --rw "$dir" "$LOCITYPER_IMAGE" jellyfish count --canonical --lower-count 2 --out-counter-len 2 \
        --mer-len 25 --threads "$THREADS" --size "$hsize" --output "$(cpath "${dir}/counts.jf.part")" "$REF_FASTA_C"; then
      rm -f "${dir}/counts.jf.part"
      echo "[WARN] jellyfish could not count the reference's k-mers"
      return 1
    fi
    mv "${dir}/counts.jf.part" "${dir}/counts.jf"
    echo "$jstamp" > "${dir}/counts.stamp"
  fi

  fetch "${LOCITYPER_PANEL_URL}.tbi" "${dir}/panel.remote.tbi" sha256 "$LOCITYPER_PANEL_TBI_SHA256" \
    || { echo "[WARN] Could not download the panel's index; run: $0 --locityper-data ${GENOME_DIR}"; return 1; }
  # The loci 250 kb either side, merged where they meet, as chrom:start-end.
  regions=$(awk -F'\t' -v OFS='\t' '{s = $2 - 250000; if (s < 0) s = 0; print $1, s, $3 + 250000}' "${dir}/loci.bed.new" \
    | sort -k1,1 -k2,2n \
    | awk -F'\t' '$1 == c && $2 <= e {if ($3 > e) e = $3; next} {if (c != "") print c ":" s + 1 "-" e; c = $1; s = $2; e = $3}
                  END {if (c != "") print c ":" s + 1 "-" e}' | tr '\n' ' ')
  # What the regions read: kept in the downloads folder by the panel URL, its
  # pin and the regions (with its sha256 beside it), so a rerun over the same
  # loci reads nothing again.
  key=$(printf '%s %s %s\n' "$LOCITYPER_PANEL_URL" "$LOCITYPER_PANEL_TBI_SHA256" "${regions% }" | _digest sha256 -)
  saved=$(dl "locityper-${LOCITYPER_PANEL_VERSION}/panel-${key:0:16}.vcf.gz")
  if [ -s "$saved" ] && [ -s "${saved}.tbi" ] && [ "$(_digest sha256 "$saved")" = "$(cat "${saved}.sha256" 2>/dev/null)" ]; then
    echo "  The panel's haplotypes over these loci are already in ${saved} (sha256 matches); not read again."
    cp "$saved" "${dir}/panel.vcf.gz.part" && cp "${saved}.tbi" "${dir}/panel.vcf.gz.tbi" \
      && mv "${dir}/panel.vcf.gz.part" "${dir}/panel.vcf.gz" || return 1
  else
  echo "  The panel's haplotypes over ${regions% }..."
  # net: slice, the regions of the loci only, saved below. SC2016: $1 to $3
  # belong to the inner bash; SC2086: one argument per region.
  # shellcheck disable=SC2016,SC2086
  if ! run_in --net --rw "$dir" "$LOCITYPER_IMAGE" bash -euo pipefail -c '
        url=$1 idx=$2 out=$3; shift 3
        tabix -h "${url}##idx##${idx}" "$@" | bgzip > "${out}.part"
        tabix -f -p vcf "${out}.part"
        mv "${out}.part.tbi" "${out}.tbi"; mv "${out}.part" "$out"' \
      _ "$LOCITYPER_PANEL_URL" "$(cpath "${dir}/panel.remote.tbi")" "$(cpath "${dir}/panel.vcf.gz")" $regions; then
    rm -f "${dir}/panel.vcf.gz.part" "${dir}/panel.vcf.gz.part.tbi"
    echo "[WARN] Could not read the panel over https; run: $0 --locityper-data ${GENOME_DIR}"
    return 1
  fi
  if [ "${SETUP_KEEP_DOWNLOADS:-true}" != false ]; then
    mkdir -p "$(dirname "$saved")" \
      && cp "${dir}/panel.vcf.gz" "${saved}.part" && cp "${dir}/panel.vcf.gz.tbi" "${saved}.tbi" \
      && _digest sha256 "${saved}.part" > "${saved}.sha256" && mv "${saved}.part" "$saved" \
      || echo "[WARN] Could not keep the panel's haplotypes in ${saved}; the next run reads them again."
  fi
  fi

  echo "  locityper target and augment on ${n} loci (the pairwise alignments take the longest)..."
  rm -rf "${dir}/db" "${dir}/db.part"
  # shellcheck disable=SC2016  # $1 to $5 belong to the inner bash
  if ! run_in --rw "$dir" "$LOCITYPER_IMAGE" bash -euo pipefail -c '
        locityper target -d "$1/db.part" -r "$2" -j "$1/counts.jf" -v "$1/panel.vcf.gz" -g GRCh38 -L "$3" -@ "$4" > "$5" 2>&1
        locityper augment -d "$1/db.part" -@ "$4" >> "$5" 2>&1' \
      _ "$(cpath "$dir")" "$REF_FASTA_C" "$(cpath "${dir}/loci.bed.new")" "$THREADS" "$(cpath "${dir}/build.log")"; then
    echo "[WARN] Locityper could not build the database; the last lines of ${dir}/build.log:"
    tail -n 20 "${dir}/build.log" || true
    return 1
  fi
  bad=""
  while IFS=$'\t' read -r _ _ _ name; do
    [ -e "${dir}/db.part/loci/${name}/success" ] || bad+="${name} "
  done < "${dir}/loci.bed.new"
  if [ -n "$bad" ]; then
    echo "[WARN] Locityper built no database for: ${bad% } (see ${dir}/build.log)"
    return 1
  fi
  mv "${dir}/db.part" "${dir}/db"
  mv "${dir}/loci.bed.new" "${dir}/loci.bed"
  echo "$stamp" > "${dir}/DONE"
  echo "[OK] Locityper's database (step 53, ${n} loci, panel ${LOCITYPER_PANEL_VERSION}): ${dir} ($(du -sh "$dir" | cut -f1))"
}

# pull_image IMAGE: docker pull, tried FETCH_TRIES times (default 3),
# FETCH_WAIT seconds apart (default 5, then 10, ...), the knobs fetch uses.
# Docker's own message goes to the terminal: a rate limit (toomanyrequests)
# or a registry error passes with a retry, a wrong tag (manifest unknown)
# does not.
pull_image() {
  local img=$1 tries=${FETCH_TRIES:-3} wait=${FETCH_WAIT:-} i
  case "$tries" in ''|*[!0-9]*) tries=0 ;; *) tries=$((10#$tries)) ;; esac
  [ "$tries" -gt 0 ] || { echo "ERROR: FETCH_TRIES must be a whole number above 0, got '${FETCH_TRIES}'" >&2; return 1; }
  case "$wait" in *[!0-9]*) echo "ERROR: FETCH_WAIT must be a whole number of seconds, got '${wait}'" >&2; return 1 ;; esac
  for ((i = 1; i <= tries; i++)); do
    "$CONTAINER_ENGINE" pull "$img" && return 0
    echo "  Pull attempt ${i}/${tries} failed: ${img} (Docker's message is above)" >&2
    [ "$i" -lt "$tries" ] && sleep "${wait:-$((i * 5))}"
  done
  return 1
}

# keep_image IMAGE: tag a digest-pinned image (tag_pinned_image in
# scripts/lib/common.sh), so a dangling-image cleanup does not delete it.
# When the tag fails, what Docker lists decides the warning, as in
# validate-setup.sh (image_is_dangling): on the classic image store an image
# with no tag is dangling and the cleanup deletes it; on the containerd image
# store, or with another local tag, it is not; or Docker did not say.
keep_image() {
  local rc=0 msg
  tag_pinned_image "$1" && return 0
  image_is_dangling "$1" || rc=$?
  case $rc in
    0) msg="Docker lists it as dangling, so a dangling-image cleanup would delete it" ;;
    1) msg="it has no pgp-pinned tag, but Docker does not list it as dangling, so a dangling-image cleanup keeps it" ;;
    *) msg="Docker did not list its dangling images, so a dangling-image cleanup may delete it" ;;
  esac
  echo "  WARNING: Could not tag ${1} as $(pinned_tag "$1"); ${msg}." >&2
  return 1
}

# --- Step 5: AnnotSV's annotation data, gnomAD-SV and MANE ----------------------------
# The AnnotSV image holds code only; without its data AnnotSV exits with an
# error. Two files are added to the bundle, each pinned in versions.env:
#   - gnomAD-SV GNOMAD_SV_VERSION's sites (GNOMAD_SV_SHA256): its common
#     deletions, duplications, insertions and inversions, written where AnnotSV
#     reads a gnomAD-SV v4 source of its benign set (ACMG 2F/4O, B_*_AFmax);
#   - MANE MANE_RELEASE's summary (MANE_SUMMARY_SHA256): the MANE Select RefSeq
#     transcripts, the list step 05 passes as -txFile, so the Tx, Frameshift,
#     Location and exon columns follow MANE Select.
ANNOTSV_DIR="${GENOME_DIR:+${GENOME_DIR}/annotsv_annotations}"
ANNOTSV_NAME="Annotations_Human_${ANNOTSV_ANNOTATIONS_VERSION}.tar.gz"
ANNOTSV_URL="https://www.lbgi.fr/~geoffroy/Annotations/${ANNOTSV_NAME}"
ANNOTSV_TARBALL="${OG_DOWNLOADS:+${OG_DOWNLOADS}/${ANNOTSV_NAME}}"
GNOMAD_SV_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/${GNOMAD_SV_VERSION}/genome_sv/gnomad.v${GNOMAD_SV_VERSION}.sv.sites.vcf.gz"
MANE_SUMMARY_URL="https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_${MANE_RELEASE}/MANE.GRCh38.v${MANE_RELEASE}.summary.txt.gz"

# install_annotsv_bundle: AnnotSV's annotation bundle, downloaded under a .part
# name, extracted into a temporary directory and moved into place only when
# complete.
install_annotsv_bundle() {
  if [ -d "${ANNOTSV_DIR}/Annotations_Human/Genes/GRCh38" ]; then
    echo "[OK] AnnotSV annotation data already present."
    return 0
  fi
  echo "  ${ANNOTSV_TERMS} See docs/05-annotsv.md."
  echo "Downloading AnnotSV ${ANNOTSV_ANNOTATIONS_VERSION} annotation data (~5.3 GB)..."
  echo "  The AnnotSV server is slow (about 0.8 MB/s measured from a GitHub runner), so this can take"
  echo "  1-2 hours. An interrupted download resumes when setup.sh runs again."
  adopt "${GENOME_DIR}/${ANNOTSV_NAME}" "$ANNOTSV_TARBALL" "$ANNOTSV_URL"
  if ! fetch "$ANNOTSV_URL" "$ANNOTSV_TARBALL"; then
    echo "[WARN] AnnotSV annotation download failed. Step 5 (AnnotSV) will be skipped until it is in place."
    echo "  Re-run setup.sh, or download it by hand:"
    echo "    curl -fL -C - -o ${ANNOTSV_TARBALL} ${ANNOTSV_URL}"
    echo "    mkdir -p ${ANNOTSV_DIR} && tar -xzf ${ANNOTSV_TARBALL} -C ${ANNOTSV_DIR}"
    return 1
  fi
  echo "  Extracting..."
  rm -rf "${ANNOTSV_DIR}.part"
  mkdir -p "${ANNOTSV_DIR}.part"
  if ! tar -xzf "$ANNOTSV_TARBALL" -C "${ANNOTSV_DIR}.part"; then
    rm -rf "${ANNOTSV_DIR}.part"
    echo "[WARN] Could not extract ${ANNOTSV_TARBALL}. Step 5 (AnnotSV) will be skipped until it is in place."
    return 1
  fi
  if [ ! -d "${ANNOTSV_DIR}.part/Annotations_Human/Genes/GRCh38" ]; then
    rm -rf "${ANNOTSV_DIR}.part"
    rm -f "$ANNOTSV_TARBALL"  # not kept: it holds no AnnotSV bundle, so the next run downloads it again
    echo "[WARN] ${ANNOTSV_TARBALL} holds no Annotations_Human/Genes/GRCh38; removed it. Step 5 (AnnotSV) will be skipped until it is in place."
    return 1
  fi
  rm -rf "$ANNOTSV_DIR"  # not kept: the older bundle, replaced
  mv "${ANNOTSV_DIR}.part" "$ANNOTSV_DIR"
  done_with "$ANNOTSV_TARBALL"
  echo "[OK] AnnotSV annotation data: ${ANNOTSV_DIR}/ ($(du -sh "$ANNOTSV_DIR" | cut -f1))"
}

# install_annotsv_gnomad_sv: gnomAD-SV's sites, reduced to the records AnnotSV
# keeps for its benign set (FILTER PASS; DEL, DUP, INS or INV; AN >= 1000,
# AF >= 0.001 and N_HOMALT >= 5; AnnotSV then applies its -benignAF, 1%, when
# it annotates). AnnotSV 3.5 reads gnomAD-SV v4 from files named
# gnomad.v4.0.sv.chr*.vcf.gz in BenignSV/GRCh38, and from the INFO keys END,
# SVTYPE, controls_and_biobanks_AN, _AF and _N_HOMALT in that order. The v4.1
# sites file has no controls-and-biobanks subset, so its AN, AF and N_HOMALT of
# all genomes are written under those keys, in a file of that name. One AnnotSV
# run (annotsv_merge_benign) merges the file into its benign_*_SV_GRCh38.sorted.bed
# and deletes it; only then does gnomad-sv.installed record what was added.
install_annotsv_gnomad_sv() {
  local benign="${ANNOTSV_DIR}/Annotations_Human/SVincludedInFt/BenignSV/GRCh38"
  local stamp="${benign}/gnomad-sv.installed"
  local want="gnomAD-SV ${GNOMAD_SV_VERSION} sha256 ${GNOMAD_SV_SHA256}"
  local src
  src=$(dl "gnomad.v${GNOMAD_SV_VERSION}.sv.sites.vcf.gz")
  local out="${benign}/gnomad.v4.0.sv.chrAll.from-v${GNOMAD_SV_VERSION}-sites.vcf.gz"
  local n
  if [ "$(cat "$stamp" 2>/dev/null)" = "$want" ]; then
    echo "[OK] gnomAD-SV ${GNOMAD_SV_VERSION} frequencies for AnnotSV already present."
    return 0
  fi
  if [ -s "$stamp" ]; then
    echo "  ${stamp} records '$(cat "$stamp")'; adding ${want}."
    echo "  The rows of the earlier release stay in AnnotSV's benign files; reinstall the bundle to drop them."
  fi
  mkdir -p "$benign"
  # fetch checks a file already there and downloads only a missing one; a
  # wrong file in the downloads folder is refused.
  adopt "${ANNOTSV_DIR}/${src##*/}" "$src" "$GNOMAD_SV_URL" sha256 "$GNOMAD_SV_SHA256"
  [ -s "$src" ] || echo "Downloading gnomAD-SV ${GNOMAD_SV_VERSION} sites (~1.7 GB)..."
  if ! fetch "$GNOMAD_SV_URL" "$src" sha256 "$GNOMAD_SV_SHA256"; then
    echo "[WARN] Could not download gnomAD-SV ${GNOMAD_SV_VERSION}. Step 5 stops until it is in place; run setup.sh again."
    return 1
  fi
  echo "  Keeping the records AnnotSV uses for its benign set..."
  if ! gzip -dc "$src" | awk -F'\t' -v OFS='\t' -v ver="$GNOMAD_SV_VERSION" '
      BEGIN {
        print "##fileformat=VCFv4.2"
        print "##source=gnomAD-SV v" ver " sites, all genomes: AN, AF and N_HOMALT under the controls_and_biobanks_* keys AnnotSV 3.5 reads"
        print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
      }
      /^#/ || $7 != "PASS" { next }
      {
        end = ""; t = ""; an = ""; af = ""; hom = ""
        n = split($8, kv, ";")
        for (i = 1; i <= n; i++) {
          eq = index(kv[i], "=")
          if (!eq) continue
          k = substr(kv[i], 1, eq - 1); v = substr(kv[i], eq + 1)
          if (k == "END") end = v
          else if (k == "SVTYPE") t = v
          else if (k == "AN") an = v
          else if (k == "AF") af = v
          else if (k == "N_HOMALT") hom = v
        }
        if (t != "DEL" && t != "DUP" && t != "INS" && t != "INV") next
        if (end == "" || an == "" || af == "" || hom == "") next
        if (an + 0 < 1000 || af + 0 < 0.001 || hom + 0 < 5) next
        print $1, $2, $3, "N", "<" t ">", ".", "PASS", "END=" end ";SVTYPE=" t ";controls_and_biobanks_AN=" an ";controls_and_biobanks_AF=" af ";controls_and_biobanks_N_HOMALT=" hom
      }' | gzip -c > "${out}.part"; then
    rm -f "${out}.part"
    echo "[WARN] Could not read ${src}; removed nothing. Run setup.sh again."
    return 1
  fi
  n=$(gzip -dc "${out}.part" | grep -vc '^#' || true)
  if [ "$n" -eq 0 ]; then
    rm -f "${out}.part"
    # The file matched its pinned sha256, so the conversion failed (a full
    # disk, the tool), not the download: it stays for the next run.
    echo "[WARN] Converting ${src} gave no record AnnotSV can use, though it matches its pinned sha256; it is kept. Free some disk if it is full, and run setup.sh again."
    return 1
  fi
  mv "${out}.part" "$out"
  if ! annotsv_merge_benign "$out"; then
    rm -f "$out"
    echo "[WARN] AnnotSV did not merge gnomAD-SV ${GNOMAD_SV_VERSION} into its benign set; kept ${src}. Run setup.sh again."
    return 1
  fi
  printf '%s\n' "$want" > "$stamp"
  done_with "$src"
  echo "[OK] gnomAD-SV ${GNOMAD_SV_VERSION}: ${n} common SVs merged into AnnotSV's benign set (${benign})"
}

# annotsv_merge_benign ADAPTER: AnnotSV merges a gnomad.v4.0.sv.chr*.vcf.gz it
# finds in BenignSV/GRCh38 into its shared benign_*_SV_GRCh38.sorted.bed at
# the start of a run, then deletes it. One AnnotSV run here, as the caller, on
# a one-deletion VCF does that merge at setup: sample runs then find nothing to
# merge, so parallel Nextflow tasks never rewrite those files at once, and
# tasks running as root never take them over. Succeeds when AnnotSV exits 0,
# ADAPTER is gone and the benign loss file exists.
annotsv_merge_benign() {
  local adapter=$1 probe="${ANNOTSV_DIR}/.setup-merge" rc=0
  local loss="${ANNOTSV_DIR}/Annotations_Human/SVincludedInFt/BenignSV/GRCh38/benign_Loss_SV_GRCh38.sorted.bed"
  rm -rf "$probe"
  mkdir -p "${probe}/out"
  printf '%s\n' '##fileformat=VCFv4.2' \
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End">' \
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">' \
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">' \
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' > "${probe}/probe.vcf"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tprobe\n' >> "${probe}/probe.vcf"
  printf 'chr20\t1000001\tsetup_probe\tN\t<DEL>\t100\tPASS\tEND=1001000;SVTYPE=DEL;SVLEN=-1000\tGT\t0/1\n' >> "${probe}/probe.vcf"
  # The image first, with pull_image's retries: left to `docker run`, one
  # broken pull (a short read from the registry) failed the merge in two
  # seconds with exit 125.
  if ! "$CONTAINER_ENGINE" image inspect "$ANNOTSV_IMAGE" &>/dev/null && ! pull_image "$ANNOTSV_IMAGE"; then
    echo "[WARN] Could not pull ${ANNOTSV_IMAGE} for AnnotSV's merge run; Docker's message is above."
    rm -rf "$probe"
    return 1
  fi
  echo "  AnnotSV merges the records into its benign set (one run, a few minutes)..."
  # The annotations directory is writable: AnnotSV rewrites its benign files there.
  run_in --rw "$ANNOTSV_DIR" --cpus 2 --memory 8g "$ANNOTSV_IMAGE" \
    AnnotSV -SVinputFile "$(cpath "${probe}/probe.vcf")" \
      -outputFile "$(cpath "${probe}/out/probe.tsv")" -outputDir "$(cpath "${probe}/out")" \
      -genomeBuild GRCh38 -annotationsDir "$(cpath "$ANNOTSV_DIR")" > "${probe}.log" 2>&1 || rc=$?
  if [ "$rc" -ne 0 ] || [ -e "$adapter" ] || [ ! -s "$loss" ]; then
    echo "[WARN] AnnotSV's merge run: exit ${rc}, $([ -e "$adapter" ] && echo "${adapter##*/} still there" || echo "${adapter##*/} consumed"). Its log: ${probe}.log"
    tail -n 20 "${probe}.log" >&2
    rm -rf "$probe"
    return 1
  fi
  rm -rf "$probe" "${probe}.log"
}

# install_annotsv_mane: the MANE Select RefSeq transcripts of MANE MANE_RELEASE,
# one per line, read by step 05 and the Nextflow ANNOTSV as -txFile.
# mane_select_refseq.release records the release and checksum.
install_annotsv_mane() {
  local out="${ANNOTSV_DIR}/mane_select_refseq.txt"
  local stamp="${ANNOTSV_DIR}/mane_select_refseq.release"
  local want="MANE ${MANE_RELEASE} sha256 ${MANE_SUMMARY_SHA256}"
  local src
  src=$(dl "MANE.GRCh38.v${MANE_RELEASE}.summary.txt.gz")
  local n
  if [ -s "$out" ] && [ "$(cat "$stamp" 2>/dev/null)" = "$want" ]; then
    echo "[OK] MANE ${MANE_RELEASE} Select transcripts for AnnotSV already present."
    return 0
  fi
  adopt "${ANNOTSV_DIR}/${src##*/}" "$src" "$MANE_SUMMARY_URL" sha256 "$MANE_SUMMARY_SHA256"
  if ! fetch "$MANE_SUMMARY_URL" "$src" sha256 "$MANE_SUMMARY_SHA256"; then
    echo "[WARN] Could not download the MANE ${MANE_RELEASE} summary. Step 5 stops until it is in place; run setup.sh again."
    return 1
  fi
  # Columns by name: RefSeq_nuc holds NM_ ids with a version, MANE_status
  # "MANE Select" or "MANE Plus Clinical".
  if ! gzip -dc "$src" | awk -F'\t' '
      NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
      $col["MANE_status"] == "MANE Select" && $col["RefSeq_nuc"] ~ /^NM_/ { print $col["RefSeq_nuc"] }' > "${out}.part"; then
    rm -f "${out}.part"
    echo "[WARN] Could not read ${src}. Run setup.sh again."
    return 1
  fi
  n=$(wc -l < "${out}.part" | tr -d ' ')
  if [ "$n" -eq 0 ]; then
    rm -f "${out}.part"
    # The file matched its pinned sha256, so reading it failed (a full disk,
    # the tool), not the download: it stays for the next run.
    echo "[WARN] Reading ${src} gave no MANE Select RefSeq transcript, though it matches its pinned sha256; it is kept. Free some disk if it is full, and run setup.sh again."
    return 1
  fi
  mv "${out}.part" "$out"
  printf '%s\n' "$want" > "$stamp"
  done_with "$src"
  echo "[OK] MANE ${MANE_RELEASE}: ${n} MANE Select transcripts for AnnotSV (${out})"
}

# install_annotsv_data: the bundle, then the two files step 05 needs with it.
install_annotsv_data() {
  install_annotsv_bundle || return 1
  local rc=0
  install_annotsv_gnomad_sv || rc=1
  install_annotsv_mane || rc=1
  return "$rc"
}

# ensure_image VAR STEP: the image of an opt-in step (`# optional` in
# versions.env, so the main setup leaves it out), pulled with its data, so
# the run does not pull it in the middle.
ensure_image() {
  local img=${!1}
  if "$CONTAINER_ENGINE" image inspect "$img" &>/dev/null; then
    echo "[OK] ${img} (step $2) already pulled."
  elif pull_image "$img"; then
    echo "[OK] ${img} (step $2) pulled."
  else
    echo "[WARN] Could not pull ${img}. Step $2 needs it; run this command again."
    return 1
  fi
  keep_image "$img"
}

# --- Step 43: the archaic data (opt-in) ----------------------------------------------
# The snpAD genotypes of three high-coverage archaic genomes (Prüfer et al.
# 2017, the Vindija 33.19 Neanderthal, with the Altai Neanderthal and the
# Altai Denisovan genotyped the same way), one GRCh37 file per chromosome and
# genome on MPI-EVA's server; the 1000 Genomes phase 3 sites file (GRCh37),
# whose allele frequencies stand for the modern population; UCSC's
# hg19ToHg38 chain. Each archaic file is checked against the md5 its folder's
# md5sums.txt lists, and each md5sums.txt against the sha256 in versions.env.
# Two masks keep the table to sequence the method can read, as Chen et al.
# 2020 do before IBDmix: each archaic genome's own FilterBed from MPI-EVA
# (coverage, mapability 35, MQ25, no tandem repeats, no indels; GRCh37), and
# the 1000 Genomes GRCh38 strict mask for the modern side. MPI-EVA publishes
# no checksum beside the FilterBed files, so versions.env pins the sha256 of
# each of the 66 files.
ARCHAIC_URL=https://cdna.eva.mpg.de/neandertal/Vindija/VCF
ARCHAIC_FILTERBED_URL=https://cdna.eva.mpg.de/neandertal/Vindija/FilterBed
ARCHAIC_1KG_SITES_URL=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz
ARCHAIC_STRICT_MASK_URL=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/StrictMask/20160622.allChr.mask.bed
ARCHAIC_SEGDUPS_URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz
ARCHAIC_CHAIN_URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# archaic_filterbed_sha NAME CHROM: the sha256 versions.env pins for one
# MPI-EVA FilterBed file, from that genome's "CHROM:sha256 ..." list.
archaic_filterbed_sha() {
  local list v
  case "$1" in
    Altai) list=$ARCHAIC_FILTERBED_SHA256_ALTAI ;;
    Vindija) list=$ARCHAIC_FILTERBED_SHA256_VINDIJA ;;
    Denisova) list=$ARCHAIC_FILTERBED_SHA256_DENISOVA ;;
    *) echo "ERROR: no FilterBed pins for '${1}'" >&2; return 1 ;;
  esac
  for v in $list; do
    [ "${v%%:*}" = "$2" ] || continue
    printf '%s' "${v#*:}"
    return 0
  done
  echo "ERROR: versions.env pins no FilterBed sha256 for ${1} chromosome ${2}" >&2
  return 1
}

# The version of the table og archaic table and strict write: a change to
# their bytes changes it, and every install rebuilds its table.
ARCHAIC_TABLE_FORMAT=archaic-table-v1
# What code= held when bin/archaic_report.py built the table: the sha256 of
# its part above "# --- genotypes (step 43)", at its last version.
ARCHAIC_TABLE_PYTHON_SHA256=6bdc4c99a8d80d96179551e37022a3cd5324f650577b9fd64e27d84383e0f1d7

# install_archaic: the table of step 43 in reference/archaic/ARCHAIC_DATA_VERSION
# (ARCHAIC_DIR: another folder inside GENOME_DIR), one chrN.tsv.gz per
# autosome, GRCh38: position, REF, ALT, the ALT count of the Altai, Vindija
# and Denisova genomes (9: not called there) and the 1000 Genomes ALT
# frequency of all and of each continental group. The rows are the 1000
# Genomes biallelic SNVs where an archaic genome was called, and the sites
# where one of the three carries a non-reference allele: the sites IBDmix's
# own generate_gt would keep from the same files, each archaic genome read
# only inside its own FilterBed and every row inside the 1000 Genomes GRCh38
# strict mask. Per chromosome: the three
# archaic files are downloaded (about 165 GB in all) to the downloads folder
# and kept there (done_with; with SETUP_KEEP_DOWNLOADS=false each is deleted
# once the chromosome is built, 13 GB at most at a time); an interrupted
# install resumes at the chromosome it stopped in. The 1000 Genomes file
# (1.5 GB) is kept the same way. ARCHAIC_REGIONS (a GRCh37 BED, chromosomes named 1 to
# 22) keeps only those regions: the e2e case builds chr22 alone with it.
install_archaic() {
  local dir=${ARCHAIC_DIR:-${GENOME_DIR}/reference/archaic/${ARCHAIC_DATA_VERSION}}
  local regions=${ARCHAIC_REGIONS:-} rstamp=all chroms="" c g v k stamp wk n
  local genomes="Altai:Altai Vindija33.19:Vindija Denisova:Denisova"
  cpath "$dir" >/dev/null || return 1
  { [ -s "${dir}/hg19ToHg38.over.chain.gz" ] && [ "$(_digest sha256 "${dir}/hg19ToHg38.over.chain.gz")" = "$ARCHAIC_CHAIN_SHA256" ]; } \
    || echo "  ${UCSC_CHAIN_TERMS} See docs/43-archaic.md."
  echo "Installing the archaic data ${ARCHAIC_DATA_VERSION} (step 43) in ${dir}"
  for v in BCFTOOLS_IMAGE PICARD_IMAGE OG_IMAGE IBDMIX_IMAGE; do ensure_image "$v" 43 || return 1; done
  need_og "builds the archaic table" || return 1
  [ -s "$REF_DICT" ] || { echo "[WARN] ${REF_DICT} is missing: run $0 ${GENOME_DIR} first (it writes the sequence dictionary)"; return 1; }
  local adl x
  adl=$(dl "archaic-${ARCHAIC_DATA_VERSION}")
  mkdir -p "${dir}/md5" "$adl" "${dir}/work" "${dir}/mask" "${dir}/exclude"
  # An older setup downloaded into download/ here: those files move to the
  # downloads folder, where fetch checks them.
  for x in "${dir}/download/"*; do
    [ ! -e "$x" ] || adopt "$x" "${adl}/${x##*/}"
  done
  rmdir "${dir}/download" 2>/dev/null || true
  for v in "Altai|${ARCHAIC_ALTAI_MD5SUMS_SHA256}" "Vindija33.19|${ARCHAIC_VINDIJA_MD5SUMS_SHA256}" "Denisova|${ARCHAIC_DENISOVA_MD5SUMS_SHA256}"; do
    IFS='|' read -r g k <<< "$v"
    [ -s "${dir}/md5/${g}.txt" ] && [ "$(_digest sha256 "${dir}/md5/${g}.txt")" = "$k" ] && continue
    if ! fetch "${ARCHAIC_URL}/${g}/md5sums.txt" "${dir}/md5/${g}.txt" sha256 "$k"; then
      echo "[WARN] The md5 list of ${g} at ${ARCHAIC_URL}/${g}/ is not the one versions.env pins: read what changed before moving the pin."
      return 1
    fi
  done
  if ! { [ -s "${dir}/hg19ToHg38.over.chain.gz" ] && [ "$(_digest sha256 "${dir}/hg19ToHg38.over.chain.gz")" = "$ARCHAIC_CHAIN_SHA256" ]; } \
     && ! fetch "$ARCHAIC_CHAIN_URL" "${dir}/hg19ToHg38.over.chain.gz" sha256 "$ARCHAIC_CHAIN_SHA256"; then
    echo "[WARN] Could not download UCSC's hg19ToHg38 chain; run: $0 --archaic ${GENOME_DIR}"
    return 1
  fi
  if ! { [ -s "${dir}/strict_GRCh38.bed" ] && [ "$(_digest sha256 "${dir}/strict_GRCh38.bed")" = "$ARCHAIC_STRICT_MASK_SHA256" ]; } \
     && ! fetch "$ARCHAIC_STRICT_MASK_URL" "${dir}/strict_GRCh38.bed" sha256 "$ARCHAIC_STRICT_MASK_SHA256"; then
    echo "[WARN] Could not download the 1000 Genomes GRCh38 strict mask; run: $0 --archaic ${GENOME_DIR}"
    return 1
  fi
  if ! { [ -s "${dir}/segdups_hg19.txt.gz" ] && [ "$(_digest sha256 "${dir}/segdups_hg19.txt.gz")" = "$ARCHAIC_SEGDUPS_SHA256" ]; } \
     && ! fetch "$ARCHAIC_SEGDUPS_URL" "${dir}/segdups_hg19.txt.gz" sha256 "$ARCHAIC_SEGDUPS_SHA256"; then
    echo "[WARN] Could not download UCSC's segmental duplications; run: $0 --archaic ${GENOME_DIR}"
    return 1
  fi
  if [ -n "$regions" ]; then
    [ -s "$regions" ] || { echo "[WARN] ARCHAIC_REGIONS=${regions} is not a file"; return 1; }
    cp "$regions" "${dir}/regions.bed"
    rstamp=$(_digest sha256 "${dir}/regions.bed")
    chroms=$(awk '$1 ~ /^([1-9]|1[0-9]|2[0-2])$/ {print $1}' "${dir}/regions.bed" | sort -un | tr '\n' ' ')
  else
    rm -f "${dir}/regions.bed"
    chroms=$(seq 1 22 | tr '\n' ' ')
  fi
  [ -n "${chroms// /}" ] || { echo "[WARN] ARCHAIC_REGIONS=${regions} names no autosome (1 to 22)"; return 1; }
  stamp="data=${ARCHAIC_DATA_VERSION} regions=${rstamp} bcftools=${BCFTOOLS_IMAGE} picard=${PICARD_IMAGE} masks=$(_digest sha256 <(printf '%s\n%s\n%s\n%s\n%s\n' "$ARCHAIC_STRICT_MASK_SHA256" "$ARCHAIC_SEGDUPS_SHA256" "$ARCHAIC_FILTERBED_SHA256_ALTAI" "$ARCHAIC_FILTERBED_SHA256_VINDIJA" "$ARCHAIC_FILTERBED_SHA256_DENISOVA")) code=${ARCHAIC_TABLE_FORMAT}"
  # Until og archaic replaced bin/archaic_report.py, code= held the sha256 of
  # that file's table half. The last Python wrote the same table bytes as
  # og archaic table and strict (tests/parity), so a table it built counts as
  # built here: its stamp is accepted once and rewritten. Remove this one
  # release after 1.0.0.
  local legacy="${stamp% code=*} code=${ARCHAIC_TABLE_PYTHON_SHA256}"
  local have=1
  for c in $chroms; do [ -s "${dir}/chr${c}.tsv.gz" ] || have=0; done
  for v in "${dir}/DONE" "${dir}"/chr*.stamp; do
    [ "$(cat "$v" 2>/dev/null)" = "$legacy" ] && echo "$stamp" > "$v"
  done
  if [ "$have" = 1 ] && [ "$(cat "${dir}/DONE" 2>/dev/null)" = "$stamp" ]; then
    echo "[OK] archaic data ${ARCHAIC_DATA_VERSION} (step 43): ${dir}"
    return 0
  fi
  rm -f "${dir}/DONE"
  local sites="${adl}/1kg_sites.vcf.gz"
  for c in $chroms; do
    [ -s "${dir}/chr${c}.tsv.gz" ] && [ "$(cat "${dir}/chr${c}.stamp" 2>/dev/null)" = "$stamp" ] && continue
    if [ ! -s "${sites}.tbi" ]; then
      # A whole copy already there (an earlier run that stopped at the index) is kept.
      { [ -s "$sites" ] && [ "$(_digest md5 "$sites")" = "$ARCHAIC_1KG_SITES_MD5" ]; } \
        || fetch "$ARCHAIC_1KG_SITES_URL" "$sites" md5 "$ARCHAIC_1KG_SITES_MD5" \
        || { echo "[WARN] Could not download the 1000 Genomes sites file; run: $0 --archaic ${GENOME_DIR}"; return 1; }
      run_in --rw "$adl" "$BCFTOOLS_IMAGE" bcftools index -t -f "$(cpath "$sites")" \
        || { echo "[WARN] Could not index ${sites}"; return 1; }
    fi
    wk="${dir}/work/${c}"
    rm -rf "$wk"; mkdir -p "$wk"
    if [ -n "$regions" ]; then awk -v c="$c" '$1 == c' "${dir}/regions.bed" > "${wk}/regions.bed"; fi
    echo "  chr${c}: the 1000 Genomes SNVs and their frequencies..."
    # shellcheck disable=SC2016  # $1 to $4 belong to the inner bash
    if ! run_in --rw "$dir" "$BCFTOOLS_IMAGE" bash -euo pipefail -c '
          r=(-r "$2"); [ ! -s "$4" ] || r=(-R "$4")
          bcftools view "${r[@]}" -m2 -M2 -v snps -Ou "$1" \
            | bcftools query -f "%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AFR_AF\t%INFO/AMR_AF\t%INFO/EAS_AF\t%INFO/EUR_AF\t%INFO/SAS_AF\n" > "$3/modern.tsv"
          bcftools view -r "$2" -v indels -Ou "$1" \
            | bcftools query -f "%POS\t%REF\t%ALT\n" > "$3/indels.tsv"' \
          _ "$(cpath "$sites")" "$c" "$(cpath "$wk")" "$(cpath "${wk}/regions.bed")"; then
      echo "[WARN] Could not read chr${c} of the 1000 Genomes sites"; return 1
    fi
    # The regions no genome is read at: UCSC's segmental duplications, and
    # each 1000 Genomes indel with 5 bp on each side, the window reaching past
    # the longer of its REF and ALT as the published script writes it (Chen et
    # al.'s own excluded beds, minus the CpG one, which needs a GRCh37
    # reference and three primate SNP sets this pipeline does not install;
    # docs/43-archaic.md says so). The indels are read from the whole
    # chromosome even when ARCHAIC_REGIONS cuts the table down, so an indel
    # just outside a region still excludes the sites next to it.
    { awk -v c="$c" 'BEGIN {OFS = "\t"} $2 == "chr" c {print c, $3, $4}' <(gzip -dc "${dir}/segdups_hg19.txt.gz")
      awk -F'\t' -v c="$c" 'BEGIN {OFS = "\t"}
        { n = length($2); m = split($3, a, ","); for (i = 1; i <= m; i++) if (length(a[i]) > n) n = length(a[i])
          s = $1 - 6; if (s < 0) s = 0
          print c, s, $1 + 5 + n }' "${wk}/indels.tsv"
    } | sort -k2,2n -k3,3n | gzip -c > "${dir}/exclude/chr${c}.bed.gz"
    n=$(gzip -dc "${dir}/exclude/chr${c}.bed.gz" | awk 'END {print NR + 0}')
    if [ "$n" -lt 100 ]; then
      echo "[WARN] the exclusion list of chr${c} has ${n} regions: the segmental duplications or the 1000 Genomes indels did not arrive"; return 1
    fi
    echo "  chr${c}: ${n} excluded regions (segmental duplications and indel neighbourhoods)"
    for v in $genomes; do
      g=${v%%:*}
      fetch "${ARCHAIC_URL}/${g}/chr${c}_mq25_mapab100.vcf.gz" "${adl}/${g}.chr${c}.vcf.gz" md5 "file://${dir}/md5/${g}.txt" \
        || { echo "[WARN] Could not download chr${c} of ${g}; run: $0 --archaic ${GENOME_DIR} (it resumes)"; return 1; }
      k=$(archaic_filterbed_sha "${v##*:}" "$c") || return 1
      { [ -s "${dir}/mask/${g}.chr${c}.bed.gz" ] && [ "$(_digest sha256 "${dir}/mask/${g}.chr${c}.bed.gz")" = "$k" ]; } \
        || fetch "${ARCHAIC_FILTERBED_URL}/${g}/chr${c}_mask.bed.gz" "${dir}/mask/${g}.chr${c}.bed.gz" sha256 "$k" \
        || { echo "[WARN] Could not download the FilterBed of chr${c} of ${g}; run: $0 --archaic ${GENOME_DIR} (it resumes)"; return 1; }
    done
    echo "  chr${c}: the archaic genotypes at those SNVs and at every archaic non-reference site..."
    # The sites to keep: the 1000 Genomes SNVs, and every position where one of
    # the three genomes has a non-reference call. Then each genome's called
    # genotypes at them. The downloads have no index: -T (targets) streams
    # through a file, where -R would need one. They have no ##contig lines
    # either, so they are read as text only: htslib adds the contig while it
    # parses, but a BCF stream (-Ou) carries the header written before that
    # and the next reader stops on "Invalid CONTIG id".
    # shellcheck disable=SC2016  # $1 to $4 belong to the inner bash
    if ! run_in --rw "$dir" "$BCFTOOLS_IMAGE" bash -euo pipefail -c '
          r=(); [ ! -s "$4" ] || r=(-T "$4")
          for g in Altai Vindija33.19 Denisova; do
            bcftools query "${r[@]}" -i "GT=\"alt\"" -f "%POS\n" "$1/$g.chr$2.vcf.gz"
          done > "$3/nonref.pos"
          { cut -f1 "$3/modern.tsv"; cat "$3/nonref.pos"; } | sort -nu | awk -v c="$2" "{print c \"\t\" \$1}" > "$3/keep.tsv"
          for g in Altai Vindija33.19 Denisova; do
            bcftools query -T "$3/keep.tsv" -i "GT!=\"mis\"" -f "%POS\t%REF\t%ALT\t[%GT]\n" "$1/$g.chr$2.vcf.gz" > "$3/$g.tsv"
          done' _ "$(cpath "$adl")" "$c" "$(cpath "$wk")" "$(cpath "${wk}/regions.bed")"; then
      echo "[WARN] Could not read the archaic genotypes of chr${c}"; return 1
    fi
    if ! run_in --rw "$dir" -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$OG_IMAGE" \
          /pgp-bin/og archaic table --chrom "$c" --modern "$(cpath "${wk}/modern.tsv")" \
            --archaic "Altai=$(cpath "${wk}/Altai.tsv")" --archaic "Vindija=$(cpath "${wk}/Vindija33.19.tsv")" \
            --archaic "Denisova=$(cpath "${wk}/Denisova.tsv")" \
            --mask "Altai=$(cpath "${dir}/mask/Altai.chr${c}.bed.gz")" \
            --mask "Vindija=$(cpath "${dir}/mask/Vindija33.19.chr${c}.bed.gz")" \
            --mask "Denisova=$(cpath "${dir}/mask/Denisova.chr${c}.bed.gz")" \
            --exclude "$(cpath "${dir}/exclude/chr${c}.bed.gz")" \
            --out "$(cpath "${wk}/hg19.vcf")"; then
      echo "[WARN] Could not merge chr${c}"; return 1
    fi
    echo "  chr${c}: lifting to GRCh38 (Picard LiftoverVcf, UCSC's chain)..."
    if ! run_in --rw "$dir" --cpus 2 --memory 10g "$PICARD_IMAGE" java -Xmx8g -jar /usr/picard/picard.jar LiftoverVcf \
          I="$(cpath "${wk}/hg19.vcf")" O="$(cpath "${wk}/hg38.vcf")" CHAIN="$(cpath "${dir}/hg19ToHg38.over.chain.gz")" \
          R="$REF_FASTA_C" REJECT="$(cpath "${wk}/rejected.vcf")" WARN_ON_MISSING_CONTIG=true \
          TMP_DIR="$(cpath "$wk")" > "${wk}/liftover.log" 2>&1; then
      tail -n 20 "${wk}/liftover.log"
      echo "[WARN] Picard could not lift chr${c}"; return 1
    fi
    # The rows that land on the same chromosome, once per position. A row
    # outside the 1000 Genomes GRCh38 strict mask is kept with its archaic
    # calls marked masked (as the FilterBed ones are), not dropped: the
    # published workflow hands IBDmix that mask, and IBDmix still scores a
    # mismatch at a masked site.
    # shellcheck disable=SC2016  # $1 and $2 belong to the inner bash
    if ! run_in --rw "$dir" "$BCFTOOLS_IMAGE" bash -euo pipefail -c \
           'bcftools query -i "CHROM==\"chr$2\"" -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/Altai\t%INFO/Vindija\t%INFO/Denisova\t%INFO/AF_ALL\t%INFO/AF_AFR\t%INFO/AF_AMR\t%INFO/AF_EAS\t%INFO/AF_EUR\t%INFO/AF_SAS\n" "$1" > "$3/lifted.tsv"' \
           _ "$(cpath "${wk}/hg38.vcf")" "$c" "$(cpath "$wk")"; then
      echo "[WARN] Could not read the lifted chr${c}"; return 1
    fi
    awk -v c="chr${c}" '$1 == c' "${dir}/strict_GRCh38.bed" > "${wk}/strict.bed"
    if [ ! -s "${wk}/strict.bed" ]; then
      echo "[WARN] the GRCh38 strict mask has no interval on chr${c}: ${dir}/strict_GRCh38.bed"; return 1
    fi
    if ! run_in --rw "$dir" -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$OG_IMAGE" \
          /pgp-bin/og archaic strict --chrom "$c" --strict "$(cpath "${wk}/strict.bed")" \
            --rows "$(cpath "${wk}/lifted.tsv")" --out "$(cpath "${wk}/marked.tsv")"; then
      echo "[WARN] Could not mark the strict mask on chr${c}"; return 1
    fi
    if ! { printf 'chrom\tpos\tref\talt\tAltai\tVindija\tDenisova\tAF_ALL\tAF_AFR\tAF_AMR\tAF_EAS\tAF_EUR\tAF_SAS\n'
           cat "${wk}/marked.tsv"
         } | gzip -c > "${dir}/chr${c}.tsv.gz.part"; then
      echo "[WARN] Could not write chr${c} of the table"; return 1
    fi
    n=$(gzip -dc "${dir}/chr${c}.tsv.gz.part" | awk 'END {print NR - 1}')
    if [ "$n" -lt 1 ]; then echo "[WARN] chr${c} of the table has no row"; rm -f "${dir}/chr${c}.tsv.gz.part"; return 1; fi
    mv "${dir}/chr${c}.tsv.gz.part" "${dir}/chr${c}.tsv.gz"
    echo "$stamp" > "${dir}/chr${c}.stamp"
    echo "  chr${c}: ${n} sites ($(awk '!/^#/ {n++} END {print n+0}' "${wk}/rejected.vcf" 2>/dev/null || echo 0) not lifted)"
    for v in $genomes; do done_with "${adl}/${v%%:*}.chr${c}.vcf.gz"; done
    rm -rf "$wk"
  done
  # A table of other chromosomes left by an earlier install with other regions goes.
  for v in "${dir}"/chr*.tsv.gz; do
    [ -e "$v" ] || continue
    c=$(basename "$v" .tsv.gz); c=${c#chr}
    [[ " ${chroms} " == *" ${c} "* ]] || rm -f "$v" "${dir}/chr${c}.stamp" "${dir}/mask/"*".chr${c}.bed.gz" "${dir}/exclude/chr${c}.bed.gz"  # not kept: what was installed for a chromosome no longer chosen
  done
  rm -rf "${dir}/work"
  done_with "$sites" "${sites}.tbi"
  echo "$stamp" > "${dir}/DONE"
  echo "[OK] archaic data ${ARCHAIC_DATA_VERSION} (step 43): ${dir} ($(echo "$chroms" | wc -w | tr -d ' ') chromosomes)"
}

# install_vep_cache_opt_in: the VEP cache step 13 and the Nextflow VEP read,
# installed before the first run; without it run-all.sh skips VEP, vcfanno,
# the clinical filter and slivar. The tarball (~26 GB) and the unpacked tree
# (~30 GB) are both on disk until the install ends.
install_vep_cache_opt_in() {
  local dir="${GENOME_DIR}/vep_cache"
  if [ -f "${dir}/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt" ]; then
    echo "[OK] VEP ${VEP_CACHE_RELEASE} cache already present: ${dir}/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38"
    return 0
  fi
  clinical_disk_ok "$dir" 56 "VEP ${VEP_CACHE_RELEASE} cache" "$(source_of vep_cache path)" || return 1
  echo "Installing the VEP ${VEP_CACHE_RELEASE} cache: ~26 GB download, ~30 GB unpacked, both on disk until it ends."
  if ! install_vep_cache "$dir" "$VEP_CACHE_RELEASE"; then
    echo "[WARN] Could not install the VEP ${VEP_CACHE_RELEASE} cache. Run this command again: the download resumes."
    return 1
  fi
  echo "[OK] VEP ${VEP_CACHE_RELEASE} cache: ${dir}/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38"
}

# Each missing item gets its own disk check. A skipped or failed item must
# not prevent the smaller bundle or the other cache from installing; it makes
# the return status 1.
install_clinical_data() {
  local dir="${GENOME_DIR}/vep_cache" rel peak missing=0
  echo "Clinical data: about 57 GB to download, 65 GB unpacked; allow 1 to 3 hours."
  for rel in "$VEP_CACHE_RELEASE" bundle "$PCGR_VEP_CACHE_RELEASE"; do
    if [ "${SETUP_EVERYTHING_SCOPE:-}" = row ] && [ "$rel" = "$VEP_CACHE_RELEASE" ]; then continue; fi
    if [ "$rel" = bundle ]; then
      [ ! -d "${GENOME_DIR}/pcgr_data/${PCGR_DATA_BUNDLE}/data" ] || { echo "[OK] PCGR/CPSR bundle already present."; continue; }
      if clinical_disk_ok "${GENOME_DIR}/pcgr_data" 15 'PCGR/CPSR bundle' "$(source_of pcgr_bundle path)"; then
        install_pcgr_bundle "${GENOME_DIR}/pcgr_data" || { missing=1; echo "[WARN] Could not install PCGR/CPSR. Step 17 needs it; run: ./scripts/setup.sh --clinical-data \"${GENOME_DIR}\""; }
      else
        missing=1
      fi
      [ -d "${GENOME_DIR}/pcgr_data/${PCGR_DATA_BUNDLE}/data" ] \
        || note_missing "PCGR/CPSR bundle ${PCGR_DATA_BUNDLE} (step 17)" "${SETUP_CMD} --clinical-data \"${GENOME_DIR}\""
    else
      [ ! -f "${dir}/homo_sapiens/${rel}_GRCh38/info.txt" ] || { echo "[OK] VEP ${rel} cache already present."; continue; }
      peak=56
      [ "$rel" != "$PCGR_VEP_CACHE_RELEASE" ] || peak=52
      if clinical_disk_ok "$dir" "$peak" "VEP ${rel} cache" "$(source_of "$(vep_cache_row "$rel")" path)"; then
        if ! install_vep_cache "$dir" "$rel"; then
          missing=1
          echo "[WARN] Could not install VEP ${rel}. Run: ./scripts/setup.sh --clinical-data \"${GENOME_DIR}\" (interrupted downloads resume)."
          if [ "$rel" = "$VEP_CACHE_RELEASE" ]; then
            echo "  To retry only this cache: ./scripts/setup.sh --vep-cache \"${GENOME_DIR}\""
          fi
        fi
      else
        missing=1
      fi
      if [ ! -f "${dir}/homo_sapiens/${rel}_GRCh38/info.txt" ]; then
        if [ "$rel" = "$VEP_CACHE_RELEASE" ]; then
          note_missing "VEP ${rel} cache (steps 13, 30, 23 and 31)" "${SETUP_CMD} --vep-cache \"${GENOME_DIR}\""
        else
          note_missing "VEP ${rel} cache (step 17, CPSR)" "${SETUP_CMD} --clinical-data \"${GENOME_DIR}\""
        fi
      fi
    fi
  done
  return "$missing"
}

# Expand the pin names listed in versions.env to file, digest type and value.
setup_pin_rows() {
  local key keys=()
  read -r -a keys <<< "$1"
  for key in "${keys[@]}"; do printf '%s\n' "${!key}"; done
}

# Install a published data/index pair. Keep existing pairs (including manual
# installs). Before adding an index to a lone data file, verify that the data
# is the matching release. A failed index fetch leaves the data for a retry.
install_indexed_downloads() (
  local dir=$1 base=$2 files=$3 rows file kind pin index index_kind index_pin dest rc=0
  rows=$(setup_pin_rows "$files") || exit 1
  mkdir -p "$dir" || exit 1
  lock_acquire "${dir}/.setup-indexed.lock" || exit 1
  trap 'lock_release "${dir}/.setup-indexed.lock"' EXIT
  while read -r file kind pin && read -r index index_kind index_pin; do
    dest="${dir}/${file}"
    if [ -s "$dest" ] && [ -s "${dir}/${index}" ]; then
      continue
    fi
    rm -f "${dir}/${index}" || exit 1  # not kept: an index without its data, fetched again below
    # fetch keeps a data file that matches its pin and downloads it otherwise.
    if ! fetch "${base}/${file}" "$dest" "$kind" "$pin"; then
      rc=1
      continue
    fi
    fetch "${base}/${index}" "${dir}/${index}" "$index_kind" "$index_pin" || rc=1
  done <<< "$rows"
  exit "$rc"
)

annotation_score_terms() {
  printf '  %s\n' \
    'CADD: free for non-commercial and academic use. Commercial use requires a licence from the University of Washington.' \
    'SpliceAI precomputed scores: academic and not-for-profit use only. Other use requires a commercial licence from Illumina.' \
    'REVEL: the authors allow non-commercial use. Zenodo lists ODbL 1.0. Contact the maintainers for other uses.' \
    'AlphaMissense: CC BY-NC-SA 4.0 (non-commercial, share-alike).'
}

install_revel() (
  local dir=$1 dest="${1}/revel_grch38.tsv.gz" zip stage tsv bgz
  [ -s "$dest" ] && [ -s "${dest}.tbi" ] && exit 0
  need_og "converts the REVEL scores" || exit 1
  ensure_image GATK_IMAGE 30 || exit 1
  # The zip and its pin are row revel of assets/data-sources.tsv.
  zip=$(source_of revel path) || exit 1
  og_setup fetch --source revel > /dev/null || exit 1
  stage=$(mktemp -d "${dir}/.revel.XXXXXX") || exit 1
  trap 'rm -rf "$stage"' EXIT
  tsv=$(cpath "${stage}/revel_grch38.tsv") || exit 1
  bgz="${tsv}.gz"
  run_in --rw "$stage" -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$OG_IMAGE" \
    /pgp-bin/og setup revel "$(cpath "$zip")" "$tsv" || exit 1
  run_in --rw "$stage" "$GATK_IMAGE" bgzip -f "$tsv" || exit 1
  run_in --rw "$stage" "$GATK_IMAGE" tabix -s 1 -b 2 -e 2 "$bgz" || exit 1
  [ -s "${stage}/${dest##*/}" ] && [ -s "${stage}/${dest##*/}.tbi" ] || exit 1
  rm -f "${dest}.tbi" || exit 1
  mv -f "${stage}/${dest##*/}" "$dest" && mv -f "${stage}/${dest##*/}.tbi" "${dest}.tbi" || exit 1
  done_with "$zip"
)

install_alphamissense() (
  local dir=$1 dest="${1}/AlphaMissense_hg38.tsv.gz" stage url sha
  [ -s "$dest" ] && [ -s "${dest}.tbi" ] && exit 0
  ensure_image GATK_IMAGE 30 || exit 1
  # The URL and its pin are row alphamissense of assets/data-sources.tsv.
  url=$(source_of alphamissense url) && sha=$(source_of alphamissense value) || exit 1
  if [ ! -s "$dest" ] || [ "$(_digest sha256 "$dest")" != "$sha" ]; then
    rm -f "${dest}.tbi" || exit 1
    fetch "$url" "$dest" sha256 "$sha" || exit 1
  fi
  stage=$(mktemp -d "${dir}/.alphamissense.XXXXXX") || exit 1
  trap 'rm -rf "$stage"' EXIT
  # The published file is BGZF. A hard link lets tabix write its index in
  # staging without another copy of the scores.
  ln "$dest" "${stage}/${dest##*/}" || cp -f "$dest" "${stage}/${dest##*/}" || exit 1
  run_in --rw "$stage" "$GATK_IMAGE" tabix -s 1 -b 2 -e 2 "$(cpath "${stage}/${dest##*/}")" || exit 1
  [ -s "${stage}/${dest##*/}.tbi" ] && mv -f "${stage}/${dest##*/}.tbi" "${dest}.tbi"
)

annotation_score_disk_ok() {
  local dir=$1 files=$2 rows file pin index dest bytes size part remaining=0 free retry
  rows=$(setup_pin_rows "$files") || return 1
  while read -r file _ pin && read -r index _; do
    dest="${dir}/${file}"
    if [ -s "$dest" ] && [ -s "${dir}/${index}" ]; then continue; fi
    case "$file" in
      whole_genome_SNVs.tsv.gz) bytes=$CADD_SNV_BYTES ;;
      gnomad.genomes.r4.0.indel.tsv.gz) bytes=$CADD_INDEL_BYTES ;;
      *) bytes=${pin%%:*} ;;  # SpliceAI pins start with the byte count.
    esac
    # Credit a size-matching file on an index retry. The installer verifies
    # its content before reusing it.
    if [ -s "$dest" ]; then
      size=$(stat -L -c %s "$dest" 2>/dev/null || stat -L -f %z "$dest") || return 1
      if [ "$size" = "$bytes" ]; then continue; fi
    fi
    if [ -s "${dest}.part" ]; then
      part=$(stat -L -c %s "${dest}.part" 2>/dev/null || stat -L -f %z "${dest}.part") || return 1
      bytes=$((bytes - part))
    fi
    if [ "$bytes" -gt 0 ]; then remaining=$((remaining + bytes)); fi
  done <<< "$rows"
  [ "$remaining" -gt 0 ] || return 0

  retry="${SETUP_CMD} --annotation-scores $(printf '%q' "$GENOME_DIR")"
  free=$(df -Pk "$dir" 2>/dev/null | awk 'NR == 2 {print $4}') || free=""
  case "$free" in
    ''|*[!0-9]*)
      echo "[WARN] Cannot read free space for annotation scores. Retry: ${retry}" >&2
      return 1 ;;
  esac
  if [ "$free" -lt "$(( (remaining + 1023) / 1024 + 5 * 1048576 ))" ]; then
    echo "[WARN] Skipping annotation scores: needs ${remaining} more bytes for CADD and SpliceAI plus a 5 GiB margin." >&2
    echo "  Install later: ${retry}" >&2
    return 1
  fi
}

install_annotation_scores() (
  local dir="${GENOME_DIR}/annotations" rc=0 i key file raw masked_files=""
  local pin_keys=()
  mkdir -p "$dir" || exit 1
  lock_acquire "${dir}/.setup-scores.lock" || exit 1
  trap 'lock_release "${dir}/.setup-scores.lock"' EXIT
  read -r -a pin_keys <<< "$SPLICEAI_MASKED_FILES"
  for ((i = 0; i < ${#pin_keys[@]}; i += 2)); do
    key=${pin_keys[i]}
    read -r file _ <<< "${!key}"
    raw="${dir}/${file/.masked./.raw.}"
    if [ -s "$raw" ] && [ -s "${raw}.tbi" ]; then
      echo "[OK] SpliceAI raw scores already present: ${raw} and ${raw}.tbi"
    else
      masked_files+=" ${key} ${pin_keys[i+1]}"
    fi
  done
  annotation_score_disk_ok "$dir" "$CADD_FILES $masked_files" || exit 1
  install_indexed_downloads "$dir" "$CADD_BASE_URL" "$CADD_FILES" || rc=1
  if [ -n "$masked_files" ]; then
    install_indexed_downloads "$dir" "$SPLICEAI_MASKED_BASE_URL" "$masked_files" || rc=1
  fi
  install_revel "$dir" || { echo '[WARN] Could not install REVEL.' >&2; rc=1; }
  install_alphamissense "$dir" || { echo '[WARN] Could not install AlphaMissense.' >&2; rc=1; }
  exit "$rc"
)

install_somatic_data() {
  install_indexed_downloads "${GENOME_DIR}/somatic" "$SOMATIC_DATA_BASE_URL" "$SOMATIC_DATA_FILES"
}

has_bwa_index() {
  local ext
  for ext in amb ann bwt pac sa; do
    [ -s "${1}.${ext}" ] || return 1
  done
}

install_bwa_index() (
  local dir stage prefix ext
  dir=$(dirname "$REF_FASTA")
  [ -s "$REF_FASTA" ] || { echo "ERROR: install the reference first: ${REF_FASTA}" >&2; exit 1; }
  lock_acquire "${REF_FASTA}.bwa-index.lock" || exit 1
  stage=""
  trap 'if [ -n "$stage" ]; then rm -rf "$stage"; fi; lock_release "${REF_FASTA}.bwa-index.lock"' EXIT
  if has_bwa_index "$REF_FASTA"; then
    echo "[OK] Classic BWA index already present: ${REF_FASTA}"
    exit 0
  fi
  ensure_image BWA_IMAGE 04b || exit 1
  stage=$(mktemp -d "${dir}/.bwa-index.XXXXXX") || exit 1
  prefix="${stage}/index"
  echo 'Building the classic BWA index for GRIDSS (about 5 GB on disk and 5.4 GB RAM)...'
  run_in --rw "$stage" "$BWA_IMAGE" bwa index -p "$(cpath "$prefix")" "$REF_FASTA_C" || exit 1
  has_bwa_index "$prefix" || { echo 'ERROR: bwa index did not write all five files.' >&2; exit 1; }
  # Publish .sa last so a partial move cannot look complete to GRIDSS.
  rm -f "${REF_FASTA}.sa" || exit 1
  for ext in amb ann bwt pac sa; do mv -f "${prefix}.${ext}" "${REF_FASTA}.${ext}" || exit 1; done
  echo "[OK] Classic BWA index: ${REF_FASTA}"
)

install_cnvpytor_data() (
  local dir="${GENOME_DIR}/reference/cnvpytor" rows file kind pin rc=0
  rows=$(setup_pin_rows "$CNVPYTOR_DATA_FILES") || exit 1
  mkdir -p "$dir" || exit 1
  lock_acquire "${dir}/.setup.lock" || exit 1
  trap 'lock_release "${dir}/.setup.lock"' EXIT
  while read -r file kind pin; do
    [ -s "${dir}/${file}" ] && continue
    fetch "${CNVPYTOR_DATA_BASE_URL}/${file}" "${dir}/${file}" "$kind" "$pin" || rc=1
  done <<< "$rows"
  exit "$rc"
)

install_pypgx_data() (
  local dir="${GENOME_DIR}/reference" dest stage=""
  dest="${dir}/pypgx-bundle"
  mkdir -p "$dir" || exit 1
  lock_acquire "${dest}.setup.lock" || exit 1
  trap 'if [ -n "$stage" ]; then rm -rf "$stage"; fi; lock_release "${dest}.setup.lock"' EXIT
  if [ -d "$dest" ]; then
    [ "$(git -c safe.directory="$dest" -C "$dest" rev-parse HEAD 2>/dev/null)" = "$PYPGX_BUNDLE_COMMIT" ] && exit 0
    echo "ERROR: move the existing PyPGx bundle aside before installing ${PYPGX_BUNDLE_VERSION}: ${dest}" >&2
    exit 1
  fi
  stage=$(mktemp -d "${dir}/.pypgx.XXXXXX") || exit 1
  git clone --depth 1 --branch "$PYPGX_BUNDLE_VERSION" "$PYPGX_BUNDLE_URL" "${stage}/bundle" || exit 1
  if [ "$(git -C "${stage}/bundle" rev-parse HEAD)" != "$PYPGX_BUNDLE_COMMIT" ]; then
    echo "ERROR: PyPGx bundle tag ${PYPGX_BUNDLE_VERSION} differs from its pinned commit." >&2
    exit 1
  fi
  mv "${stage}/bundle" "$dest"
)

install_everything() {
  local manifest="${PGP_ROOT}/assets/everything.tsv" kind flag reason base_retry
  local dir_arg flags=0
  [ -s "$manifest" ] || { echo "ERROR: missing installer manifest: ${manifest}" >&2; return 1; }
  dir_arg=$(printf '%q' "$GENOME_DIR")
  base_retry="${SETUP_CMD} --everything ${dir_arg}"
  echo '=== Download and use terms (docs/00-reference-setup.md) ==='
  annotation_score_terms
  printf '  %s\n' "$CYRIUS_TERMS"
  printf '  %s\n' "${PRS_ACADEMIC_TERMS[@]}" \
    'PyPGx bundle: upstream publishes no licence terms for the data bundle.'

  echo '=== Everything base setup ==='
  if ! SETUP_EVERYTHING_SCOPE=base "$0" "$GENOME_DIR"; then
    note_missing 'base setup (see the base setup messages above)' "$base_retry"
  fi
  while IFS=$'\t' read -r kind flag _ _ reason; do
    [ "$kind" = setup ] || continue
    flags=$((flags + 1))
    echo "=== Everything installer: ${flag} ==="
    if [[ "$flag" != --* || "$flag" = --everything ]]; then
      note_missing "invalid setup row: ${flag}" "$base_retry"
      continue
    fi
    # Each invocation has its own errexit state. A failed installer cannot
    # end the loop or consume further manifest rows.
    if ! SETUP_EVERYTHING_SCOPE=row STR_EXTRA_LOCI=yes "$0" "$flag" "$GENOME_DIR" </dev/null; then
      note_missing "${flag} (${reason})" "STR_EXTRA_LOCI=yes ${SETUP_CMD} ${flag} ${dir_arg}"
    fi
  done < "$manifest"
  if [ "$flags" -eq 0 ]; then
    echo "ERROR: no installer rows in ${manifest}" >&2
    return 1
  fi
  if [ "${#MISSING[@]}" -eq 0 ]; then
    echo 'Everything setup complete!'
    return 0
  fi
  print_missing
  echo "Everything setup exits with status ${SETUP_PARTIAL_EXIT}."
  return "$SETUP_PARTIAL_EXIT"
}

# install_prs_academic: the four scores of assets/pgs_scores.tsv free for
# academic research use only (licence academic-research), which step 25 and
# the PRS process score with PRS_ACADEMIC=true (Nextflow --pgs_academic true).
install_prs_academic() {
  if fetch_prs_scores academic-research; then
    echo "[OK] Academic-only PGS Catalog scores (step 25 with PRS_ACADEMIC=true): ${GENOME_DIR}/prs_scores/"
  else
    return 1
  fi
}

case "$OPT_IN" in
  prs-academic)
    if [ "${SETUP_EVERYTHING_SCOPE:-}" != row ]; then printf '  %s\n' "${PRS_ACADEMIC_TERMS[@]}"; fi
    install_prs_academic; exit $? ;;
  annotation-scores)
    if [ "${SETUP_EVERYTHING_SCOPE:-}" != row ]; then annotation_score_terms; fi
    install_annotation_scores; exit $? ;;
  somatic-data) install_somatic_data; exit $? ;;
  bwa-index) install_bwa_index; exit $? ;;
  cnvpytor-data) install_cnvpytor_data; exit $? ;;
  pypgx-data) install_pypgx_data; exit $? ;;
  ancestry-panel) install_ancestry_panel; exit $? ;;
  ancestry-fine-panel) install_ancestry_fine_panel; exit $? ;;
  cyrius)
    if [ "${SETUP_EVERYTHING_SCOPE:-}" != row ]; then printf '  %s\n' "$CYRIUS_TERMS"; fi
    ensure_image PYTHON_IMAGE 21 && install_cyrius; exit $? ;;
  parascopy-data) install_parascopy_data && ensure_image PARASCOPY_IMAGE 35; exit $? ;;
  kir-data) install_kir_data; exit $? ;;
  yleaf-data) install_yleaf_data && ensure_image YLEAF_IMAGE 37; exit $? ;;
  consumer-sites) install_consumer_sites; exit $? ;;
  ancient-panel) install_ancient_panel; exit $? ;;
  viral-data) install_viral_data; exit $? ;;
  masked-reference) install_masked_reference; exit $? ;;
  exomiser-data) install_exomiser_data && ensure_image EXOMISER_IMAGE 48; exit $? ;;
  chip-data) install_chip_data; exit $? ;;
  strchive-data) install_strchive_data; exit $? ;;
  megane-kmers) ensure_image MEGANE_IMAGE 50 && install_megane_kmers; exit $? ;;
  mocha-data) install_mocha_data; exit $? ;;
  locityper-data) install_locityper_data; exit $? ;;
  clinical-data) install_clinical_data; exit $? ;;
  vep-cache) install_vep_cache_opt_in; exit $? ;;
  annotsv-data) install_annotsv_data; exit $? ;;
  archaic) install_archaic; exit $? ;;
esac

# pull_images: pull every image setup pre-pulls (versions.env, minus the
# lines marked `# optional`), and tag the digest-pinned ones, those already
# present included. Sets PULLED, SKIPPED, FAILED and TAGGED (pinned tags added).
pull_images() {
  local img
  PULLED=0
  SKIPPED=0
  FAILED=0
  TAGGED=0
  while IFS= read -r img; do
    if "$CONTAINER_ENGINE" image inspect "$img" &>/dev/null; then
      SKIPPED=$((SKIPPED + 1))
    else
      echo "  Pulling: ${img}..."
      if pull_image "$img"; then
        PULLED=$((PULLED + 1))
      else
        echo "  WARNING: Failed to pull ${img}; Docker's message is above."
        FAILED=$((FAILED + 1))
        continue
      fi
    fi
    pinned_tag_present "$img" && continue
    if keep_image "$img"; then TAGGED=$((TAGGED + 1)); else FAILED=$((FAILED + 1)); fi
  done < <(pipeline_images)
  echo "[OK] Docker images: ${PULLED} pulled, ${SKIPPED} already present, ${FAILED} failed, ${TAGGED} pgp-pinned tag(s) added."
  if [ "$FAILED" -gt 0 ]; then
    echo ""
    echo "ERROR: ${FAILED} Docker image(s) failed to pull or tag. Fix the issues above and re-run setup."
    return 1
  fi
}

# --- ClinVar ---------------------------------------------------------------------
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CLINVARDIR="${GENOME_DIR:+${GENOME_DIR}/clinvar}"
# The raw download and the three files built from it, each with its index.
CLINVAR_SET="clinvar.vcf.gz clinvar.vcf.gz.tbi clinvar_chr.vcf.gz clinvar_chr.vcf.gz.tbi
  clinvar_pathogenic_chr.vcf.gz clinvar_pathogenic_chr.vcf.gz.tbi
  clinvar_conflicting_plp_chr.vcf.gz clinvar_conflicting_plp_chr.vcf.gz.tbi"

# clinvar_build_derived RAW_DIR OUT_DIR: from RAW_DIR/clinvar.vcf.gz build, in
# OUT_DIR (both under GENOME_DIR), clinvar_chr.vcf.gz (NCBI's 1, 2 ... MT
# renamed to chr1, chr2 ... chrM) and clinvar_pathogenic_chr.vcf.gz (the
# Pathogenic and Likely_pathogenic records), each with its .tbi.
clinvar_build_derived() {
  local raw out
  raw=$(cpath "$1") && out=$(cpath "$2") || return 2
  awk 'BEGIN { for (i = 1; i <= 22; i++) print i, "chr" i; print "X chrX"; print "Y chrY"; print "MT chrM" }' \
    > "${2}/chr_rename.txt"
  echo "  Creating chr-prefixed ClinVar..."
  run_in --rw "$2" "$BCFTOOLS_IMAGE" \
    bcftools annotate --rename-chrs "${out}/chr_rename.txt" "${raw}/clinvar.vcf.gz" \
      -Oz -o "${out}/clinvar_chr.vcf.gz" || return 1
  run_in --rw "$2" "$BCFTOOLS_IMAGE" bcftools index -f -t "${out}/clinvar_chr.vcf.gz" || return 1
  echo "  Creating pathogenic/likely pathogenic subset..."
  run_in --rw "$2" "$BCFTOOLS_IMAGE" \
    bcftools view -i 'CLNSIG~"Pathogenic" || CLNSIG~"Likely_pathogenic"' "${out}/clinvar_chr.vcf.gz" \
      -Oz -o "${out}/clinvar_pathogenic_chr.vcf.gz" || return 1
  run_in --rw "$2" "$BCFTOOLS_IMAGE" bcftools index -f -t "${out}/clinvar_pathogenic_chr.vcf.gz" || return 1
  rm -f "${2}/chr_rename.txt"
}

# clinvar_build_conflicting SRC_DIR OUT_DIR: from SRC_DIR/clinvar_chr.vcf.gz build,
# in OUT_DIR (both under GENOME_DIR), clinvar_conflicting_plp_chr.vcf.gz with its .tbi: the records
# ClinVar labels conflicting whose CLNSIGCONF counts at least one Pathogenic or
# Likely_pathogenic submission. og clinvar conflicting-subset selects them; the
# pathogenic subset holds none of them. Step 06 lists the ones a sample carries
# apart from its hits, for a clinician to review. Needs og (need_og).
clinvar_build_conflicting() {
  local src d
  src=$(cpath "$1") && d=$(cpath "$2") || return 2
  echo "  Creating the conflicting subset (records some labs call pathogenic)..."
  run_in --rw "$2" -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$BCFTOOLS_IMAGE" bash -o pipefail -c "set -e
    /pgp-bin/og clinvar conflicting-subset '${src}/clinvar_chr.vcf.gz' \\
      | bcftools view --no-version -Oz -o '${d}/clinvar_conflicting_plp_chr.vcf.gz' -
    bcftools index -f -t '${d}/clinvar_conflicting_plp_chr.vcf.gz'" || return 1
}

# clinvar_install DIR: move every file of CLINVAR_SET that DIR holds into
# clinvar/, then remove DIR. The data files go first and their indexes after.
# Returns 1 at the first move that fails (set -e does not apply inside a
# function called as `f || ...`), keeping DIR.
clinvar_install() {
  local f
  for f in $CLINVAR_SET; do
    case "$f" in *.tbi) continue ;; esac
    if [ -f "${1}/${f}" ]; then mv -f "${1}/${f}" "${CLINVARDIR}/${f}" || return 1; fi
  done
  for f in $CLINVAR_SET; do
    case "$f" in *.tbi) ;; *) continue ;; esac
    if [ -f "${1}/${f}" ]; then mv -f "${1}/${f}" "${CLINVARDIR}/${f}" || return 1; fi
  done
  rm -rf "$1"
}

# clinvar_record_release: write the ##fileDate of clinvar/clinvar.vcf.gz
# (YYYY-MM-DD) to clinvar/RELEASE, which validate-setup.sh and step 06 print.
clinvar_record_release() {
  local d
  d=$(gzip -cd "${CLINVARDIR}/clinvar.vcf.gz" 2>/dev/null | head -n 100 \
      | awk -F= '/^##fileDate=/ { print $2; exit }') || true
  case "$d" in
    [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]) d="${d:0:4}-${d:4:2}-${d:6:2}" ;;
  esac
  case "$d" in
    [0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9])
      printf '%s\n' "$d" > "${CLINVARDIR}/RELEASE.tmp"
      mv -f "${CLINVARDIR}/RELEASE.tmp" "${CLINVARDIR}/RELEASE"
      echo "[OK] ClinVar release ${d} (recorded in ${CLINVARDIR}/RELEASE)." ;;
    *)
      echo "[WARN] ${CLINVARDIR}/clinvar.vcf.gz has no ##fileDate line; its release date is unknown." ;;
  esac
}

# refresh_clinvar: download NCBI's current ClinVar into clinvar/.refresh/,
# check it against NCBI's md5, build the three derived files from it there, and
# only then replace all eight files. The normalised copy step 06 keeps is removed, so
# step 06 builds it again from the new release. RELEASE and that copy go before
# the files move: a refresh that stops between two moves leaves no release date
# on a set that may mix two releases (validate-setup.sh then says the date is
# unknown) and fails, so the next refresh starts again. When og cannot build the
# conflicting subset, the other six files still refresh and the old subset is
# removed, so it never sits beside a newer release: step 06 then says those
# records were not checked, and the next setup.sh run builds the subset.
refresh_clinvar() {
  local new="${CLINVARDIR}/.refresh"
  # A download cut short last time (*.part) is resumed; anything else is redone.
  mkdir -p "$new"
  find "$new" -mindepth 1 ! -name '*.part' -exec rm -rf {} +  # not kept: a stopped refresh's staging files, redone; its .part resumes
  echo "Downloading the current ClinVar release..."
  fetch "$CLINVAR_URL" "${new}/clinvar.vcf.gz" md5 "${CLINVAR_URL}.md5" || return 1
  fetch "${CLINVAR_URL}.tbi" "${new}/clinvar.vcf.gz.tbi" || return 1
  clinvar_build_derived "$new" "$new" || {
    echo "ERROR: could not build the ClinVar files from the new release; the old ones are kept." >&2
    return 1
  }
  if ! need_og "builds the conflicting ClinVar subset of step 06" || ! clinvar_build_conflicting "$new" "$new"; then
    rm -f "${new}/clinvar_conflicting_plp_chr.vcf.gz.tbi" "${new}/clinvar_conflicting_plp_chr.vcf.gz" \
      "${CLINVARDIR}/clinvar_conflicting_plp_chr.vcf.gz.tbi" "${CLINVARDIR}/clinvar_conflicting_plp_chr.vcf.gz"
    echo "[WARN] Could not build the conflicting ClinVar subset of the new release; step 06 says those records were not checked until setup.sh builds it." >&2
  fi
  rm -f "${CLINVARDIR}/RELEASE" \
    "${CLINVARDIR}/clinvar_pathogenic_chr.norm.vcf.gz" "${CLINVARDIR}/clinvar_pathogenic_chr.norm.vcf.gz.tbi"
  clinvar_install "$new" || {
    echo "ERROR: could not move the new ClinVar files into ${CLINVARDIR}; they may mix two releases. Run setup.sh --refresh clinvar again." >&2
    return 1
  }
  clinvar_record_release
}

if $EVERYTHING_SETUP; then
  install_everything
  exit $?
fi

echo "============================================"
echo "  OwnGenome — Setup"
echo "  Data directory: ${GENOME_DIR:-(none: --pull-only)}"
echo "============================================"
echo ""

# Check Docker
if ! command -v "$CONTAINER_ENGINE" &>/dev/null; then
  echo "ERROR: Docker is not installed."
  echo "  Install Docker: https://docs.docker.com/get-docker/"
  exit 1
fi
if ! "$CONTAINER_ENGINE" info &>/dev/null; then
  echo "ERROR: Docker daemon is not running."
  echo "  Start Docker Desktop or run: sudo systemctl start docker"
  exit 1
fi
echo "[OK] Docker is running."

if $PULL_ONLY; then
  echo ""
  echo "=== Docker Images (~10-15 GB total) ==="
  pull_images || exit 1
  exit 0
fi

if [ "$REFRESH" = clinvar ]; then
  mkdir -p "$CLINVARDIR"
  echo ""
  echo "=== Refreshing ClinVar ==="
  if [ -f "${CLINVARDIR}/RELEASE" ]; then echo "  Current release: $(cat "${CLINVARDIR}/RELEASE")"; fi
  refresh_clinvar || exit 1
  if [ -f "${CLINVARDIR}/clinvar_conflicting_plp_chr.vcf.gz.tbi" ]; then
    echo "[OK] ClinVar refreshed: raw file, chr-prefixed file, pathogenic and conflicting subsets replaced."
  else
    echo "[OK] ClinVar refreshed: raw file, chr-prefixed file and pathogenic subset replaced; the conflicting subset was not built (see the warning above)."
  fi
  exit 0
fi

# Genomes are private data: only the owner may enter the data directory.
mkdir -p "$GENOME_DIR"
chmod 700 "$GENOME_DIR"

# Check disk space (df -Pk works on Linux and macOS; GENOME_DIR may not
# exist yet, so measure its nearest existing parent)
DF_DIR="$GENOME_DIR"
while [ ! -d "$DF_DIR" ]; do
  DF_DIR=$(dirname "$DF_DIR")
done
AVAIL_KB=$(df -Pk "$DF_DIR" | awk 'NR==2 {print $4}')
AVAIL_GB=$(( ${AVAIL_KB:-0} / 1048576 ))
if [ "$AVAIL_GB" -lt 50 ]; then
  echo "WARNING: Only ${AVAIL_GB} GB free in ${GENOME_DIR}. Need at least 50 GB for references."
fi
echo "[OK] ${AVAIL_GB} GB free in ${GENOME_DIR}."

###############################################################################
# Phase 1: Reference Genome
###############################################################################
echo ""
echo "=== Phase 1: Reference Genome (~0.9 GB download, ~3.2 GB unpacked) ==="

FASTA="$REF_FASTA"
FAI="${REF_FASTA}.fai"
REFDIR=$(dirname "$FASTA")
mkdir -p "$REFDIR"

if [ -f "$FASTA" ] && [ -f "$FAI" ]; then
  echo "[OK] Reference genome already downloaded."
else
  echo "Downloading GRCh38 reference genome..."
  echo "  Source: ${REF_FASTA_URL}"

  if [ ! -f "$FASTA" ]; then
    case "$REF_FASTA_URL" in
      *.gz) DL=$(dl "${REF_FASTA_URL##*/}")
        # shellcheck disable=SC2046  # the checksum arguments are dropped when REF_FASTA_MD5 is empty
        adopt "${FASTA}.download.gz" "$DL" "$REF_FASTA_URL" $([ -z "$REF_FASTA_MD5" ] || printf 'md5 %s' "$REF_FASTA_MD5") ;;
      *) DL="$FASTA" ;;
    esac
    # shellcheck disable=SC2046  # the checksum arguments are dropped when REF_FASTA_MD5 is empty
    fetch "$REF_FASTA_URL" "$DL" $([ -z "$REF_FASTA_MD5" ] || printf 'md5 %s' "$REF_FASTA_MD5") || {
      echo "  The download is resumed when setup.sh runs again."
      exit 1
    }
    if [ "$DL" != "$FASTA" ]; then
      echo "  Unpacking $(basename "$REF_FASTA_URL")..."
      if gzip -dc "$DL" > "${FASTA}.tmp"; then
        mv -f "${FASTA}.tmp" "$FASTA"
        done_with "$DL"
      else
        rm -f "${FASTA}.tmp"
        if [ -n "$REF_FASTA_MD5" ]; then
          # It matched its md5, so the unpack failed (a full disk, gzip), not
          # the download: it stays for the next run.
          echo "ERROR: could not unpack $(basename "$REF_FASTA_URL"), though it matches its md5; it is kept in ${DL}. Free some disk if it is full, and run setup.sh again." >&2
        else
          rm -f "$DL"  # not kept: no checksum vouches for it and it does not unpack, so the next run downloads it again
          echo "ERROR: could not unpack $(basename "$REF_FASTA_URL"); removed it, run setup.sh again." >&2
        fi
        exit 1
      fi
    fi
  fi

  if [ ! -f "$FAI" ]; then
    if [ -z "$REF_FAI_MD5" ] || ! fetch "$REF_FAI_URL" "$FAI" md5 "$REF_FAI_MD5"; then
      echo "  Generating index with samtools..."
      # The index goes next to the FASTA, so its directory is writable here.
      run_in --rw "$REFDIR" "$SAMTOOLS_IMAGE" \
        samtools faidx "$REF_FASTA_C"
    fi
  fi
  echo "[OK] Reference genome downloaded."
fi
if [ -f "$FAI" ]; then echo "[OK] Reference: ${FASTA} ($(grep -c . "$FAI") sequences)"; fi

###############################################################################
# Phase 2: ClinVar Database
###############################################################################
echo ""
echo "=== Phase 2: ClinVar Database (~200 MB) ==="

mkdir -p "$CLINVARDIR"

CLINVAR="${CLINVARDIR}/clinvar.vcf.gz"
CLINVAR_TBI="${CLINVARDIR}/clinvar.vcf.gz.tbi"

if [ -f "$CLINVAR" ] && [ -f "$CLINVAR_TBI" ]; then
  echo "[OK] ClinVar database already downloaded. To replace it with the current release:"
  echo "  ./scripts/setup.sh --refresh clinvar ${GENOME_DIR}"
else
  echo "Downloading ClinVar database..."
  if [ ! -f "$CLINVAR" ]; then
    # NCBI publishes the md5 next to the file.
    fetch "$CLINVAR_URL" "$CLINVAR" md5 "${CLINVAR_URL}.md5" || exit 1
  fi
  if [ ! -f "$CLINVAR_TBI" ]; then
    fetch "${CLINVAR_URL}.tbi" "$CLINVAR_TBI" || exit 1
  fi
fi

# The two derived files are built whenever one is missing, also on a rerun
# with the raw download already present. They are built in clinvar/.build/
# and moved in when both are complete, so a half-built file is never in place.
if [ ! -f "${CLINVARDIR}/clinvar_chr.vcf.gz.tbi" ] || [ ! -f "${CLINVARDIR}/clinvar_pathogenic_chr.vcf.gz.tbi" ]; then
  rm -rf "${CLINVARDIR}/.build"
  mkdir -p "${CLINVARDIR}/.build"
  clinvar_build_derived "$CLINVARDIR" "${CLINVARDIR}/.build" || exit 1
  clinvar_install "${CLINVARDIR}/.build"
fi
if [ ! -f "${CLINVARDIR}/RELEASE" ]; then clinvar_record_release; fi
echo "[OK] ClinVar database and its chr-prefixed and pathogenic subsets are ready."

###############################################################################
# Phase 3: Docker Images
###############################################################################
echo ""
echo "=== Phase 3: Docker Images (~10-15 GB total) ==="

# The list is every *_IMAGE line of versions.env not marked `# optional`.
pull_images || exit 1

# Steps 04, 05, 07 and 25 run OwnGenome's og binary (bin/.og): built from rust/
# in OG_BUILDER, or downloaded on a release tag, when it is missing or stale.
if "${PGP_ROOT}/rust/og.sh" check --required > /dev/null 2>&1; then
  echo "[OK] og (bin/.og) already built from this checkout's rust/."
elif "${PGP_ROOT}/rust/og.sh" install; then
  echo "[OK] og (bin/.og) built."
else
  echo "[WARN] Could not build og. Steps 04, 05, 07, 25 and 33 stop without it; run: $0 --og"
  note_missing "og, OwnGenome's own binary (steps 04, 05, 07, 25 and 33)" "${SETUP_CMD} --og"
fi

# ClinVar's conflicting subset (step 06 lists the records of it a sample carries,
# apart from the hits) is built by og, so it comes after og. It is built when
# it is missing or older than clinvar_chr.vcf.gz, in clinvar/.build-conflicting/,
# and moved in once complete. Without it step 06 writes "not checked".
if [ -f "${CLINVARDIR}/clinvar_chr.vcf.gz" ] \
   && { [ ! -f "${CLINVARDIR}/clinvar_conflicting_plp_chr.vcf.gz.tbi" ] \
        || [ "${CLINVARDIR}/clinvar_chr.vcf.gz" -nt "${CLINVARDIR}/clinvar_conflicting_plp_chr.vcf.gz" ]; }; then
  rm -rf "${CLINVARDIR}/.build-conflicting"
  mkdir -p "${CLINVARDIR}/.build-conflicting"
  if need_og "builds the conflicting ClinVar subset of step 06" \
     && clinvar_build_conflicting "$CLINVARDIR" "${CLINVARDIR}/.build-conflicting" \
     && clinvar_install "${CLINVARDIR}/.build-conflicting"; then
    echo "[OK] ClinVar conflicting subset (records some labs call pathogenic, step 06) built."
  else
    rm -rf "${CLINVARDIR}/.build-conflicting"
    echo "[WARN] Could not build the conflicting ClinVar subset; step 06 says the records were not checked."
    note_missing "the conflicting ClinVar subset (step 06)" "${SETUP_CMD} --refresh clinvar \"${GENOME_DIR}\""
  fi
fi

###############################################################################
# Phase 4: Sequence dictionary and small pinned data files
###############################################################################
echo ""
echo "=== Phase 4: Sequence dictionary and pinned data files (~350 MB) ==="

# GATK and Picard (steps 03a, 20, 29, chip-to-vcf) need the reference's .dict.
if [ -f "$REF_DICT" ]; then
  echo "[OK] Sequence dictionary already present."
else
  echo "Creating the sequence dictionary ${REF_DICT}..."
  DICT_PART="${REF_DICT%.dict}.part.dict"
  rm -f "$DICT_PART"
  # The dictionary goes next to the FASTA, so its directory is writable here.
  if run_in --rw "$REFDIR" "$GATK_IMAGE" \
       gatk CreateSequenceDictionary -R "$REF_FASTA_C" -O "$(cpath "$DICT_PART")" \
     && [ -s "$DICT_PART" ]; then
    mv -f "$DICT_PART" "$REF_DICT"
    echo "[OK] Sequence dictionary created."
  else
    rm -f "$DICT_PART"
    echo "[WARN] Could not create ${REF_DICT}. Steps 03a, 20, 29 and chip-to-vcf need it; re-run setup.sh."
    note_missing "the sequence dictionary ${REF_DICT##*/} (steps 03a, 20, 29 and chip-to-vcf)" "${SETUP_CMD} \"${GENOME_DIR}\""
  fi
fi

# Each from a fixed commit or release and checked before it is stored
# (scripts/lib/common.sh, install_data_file). A failed download does not stop
# setup: the step that needs the file says so, and setup.sh fetches it on its
# next run.
for name in $DATA_FILES; do
  case "$name" in
    delly_exclude) what="Delly exclude map (step 19)" ;;
    cytoband) what="GRCh38 chromosome bands (step 10)" ;;
    hla_dat) what="IPD-IMGT/HLA ${HLA_DB_RELEASE} (step 08)" ;;
    gencode_genes) what="GENCODE ${GENCODE_RELEASE} gene coordinates (step 08)" ;;
    gencode_exons) what="GENCODE ${GENCODE_RELEASE} exons (steps 06 and 17)" ;;
  esac
  if dest=$(data_file "$name"); then
    echo "[OK] ${what} already present."
  elif install_data_file "$name"; then
    echo "[OK] ${what}: ${dest}"
  else
    echo "[WARN] Could not install the ${what}. Re-run setup.sh to try again."
    note_missing "the ${what}" "${SETUP_CMD} \"${GENOME_DIR}\""
  fi
done
# gnomAD v3.1's chrM sites (step 20 annotates the calls with them), checked
# against GNOMAD_MITO_SHA256. A failed download does not stop setup either.
GNOMAD_MITO="${GENOME_DIR}/gnomad_mito/${GNOMAD_MITO_URL##*/}"
if [ -s "$GNOMAD_MITO" ] && [ "$(_digest sha256 "$GNOMAD_MITO")" = "$GNOMAD_MITO_SHA256" ]; then
  echo "[OK] gnomAD v3.1 mitochondrial sites (step 20) already present (sha256 matches GNOMAD_MITO_SHA256)."
else
  if [ -s "$GNOMAD_MITO" ]; then
    echo "[WARN] ${GNOMAD_MITO} is not the pinned file (sha256 differs from GNOMAD_MITO_SHA256): downloading it again."
  fi
  if fetch "$GNOMAD_MITO_URL" "$GNOMAD_MITO" sha256 "$GNOMAD_MITO_SHA256"; then
    echo "[OK] gnomAD v3.1 mitochondrial sites (step 20): ${GNOMAD_MITO}"
  else
    echo "[WARN] Could not install gnomAD's mitochondrial sites. Step 20 runs without its gnomAD columns; re-run setup.sh to try again."
    note_missing "gnomAD v3.1 mitochondrial sites (step 20)" "${SETUP_CMD} \"${GENOME_DIR}\""
  fi
fi
# Step 11's allele-frequency panel (ROH_AF_PANEL_URL: gnomAD genomes' common
# biallelic SNVs, built by og dev build roh-af), checked against
# ROH_AF_PANEL_SHA256 and indexed here with the tabix of GATK_IMAGE, as
# REVEL is (scripts/ci/check-container-helpers.sh keeps bgzip and tabix out
# of BCFTOOLS_IMAGE). Without it step 11 gives every site AF 0.4, so a
# failed download warns and does not stop setup.
ROH_AF_PANEL="${GENOME_DIR}/roh_af/${ROH_AF_PANEL_URL##*/}"
# Download only a missing panel or one whose sha256 is not the pinned one;
# then index the checked file, so a failed index is retried without
# downloading the panel again.
ROH_AF_PANEL_OK=0
if [ -s "$ROH_AF_PANEL" ] && [ "$(_digest sha256 "$ROH_AF_PANEL")" = "$ROH_AF_PANEL_SHA256" ]; then
  ROH_AF_PANEL_OK=1
else
  rm -f "${ROH_AF_PANEL}.tbi"
  if fetch "$ROH_AF_PANEL_URL" "$ROH_AF_PANEL" sha256 "$ROH_AF_PANEL_SHA256"; then
    ROH_AF_PANEL_OK=1
  else
    echo "[WARN] Could not install step 11's allele-frequency panel: step 11 gives every site AF 0.4 without it. Re-run setup.sh to try again."
  fi
fi
if [ "$ROH_AF_PANEL_OK" = 1 ]; then
  if [ -s "${ROH_AF_PANEL}.tbi" ]; then
    echo "[OK] ROH allele-frequency panel (step 11) already present (sha256 matches ROH_AF_PANEL_SHA256)."
  elif run_in --rw "${GENOME_DIR}/roh_af" "$GATK_IMAGE" tabix -f -s1 -b2 -e2 "$(cpath "$ROH_AF_PANEL")"; then
    echo "[OK] ROH allele-frequency panel (step 11): ${ROH_AF_PANEL}"
  else
    rm -f "${ROH_AF_PANEL}.tbi"
    echo "[WARN] Could not index step 11's allele-frequency panel (${ROH_AF_PANEL}): step 11 gives every site AF 0.4 until it is. Re-run setup.sh to try again."
  fi
fi
# A failed download does not stop setup here either: it goes on the list of
# what is missing, which the summary prints with the command that retries it.
if [ "${SETUP_EVERYTHING_SCOPE:-}" != base ]; then
  install_sample_qc_data \
    || note_missing "somalier's sites and VerifyBamID2's panel (step 33)" "${SETUP_CMD} --sample-qc-data \"${GENOME_DIR}\""
  install_strchive_data \
    || note_missing "STRchive ${STRCHIVE_VERSION} catalog (steps 09 and 9b)" "${SETUP_CMD} --strchive-data \"${GENOME_DIR}\""
fi
install_tandem_repeats \
  || note_missing "GIAB tandem-repeat map (step 04c)" "${SETUP_CMD} \"${GENOME_DIR}\""
# A plain setup warns when a PRS download fails. Everything lists it as missing.
if ! install_prs; then
  note_missing "PGS Catalog scores or pgsc_calc (step 25)" "${SETUP_CMD} --everything \"${GENOME_DIR}\""
fi

###############################################################################
if [ "${SETUP_EVERYTHING_SCOPE:-}" != base ]; then
  echo ""
  echo "=== Phase 5: AnnotSV annotation data (~5.3 GB + ~1.7 GB download, ~20 GB unpacked) ==="

  install_annotsv_data \
    || note_missing "AnnotSV data: the bundle, gnomAD-SV ${GNOMAD_SV_VERSION} and MANE ${MANE_RELEASE} (step 05)" "${SETUP_CMD} --annotsv-data \"${GENOME_DIR}\""

  ###############################################################################
  # Phase 6: Clinical reference data (default)
  ###############################################################################
  echo ""
  echo "=== Phase 6: Clinical reference data ==="
  if $CLINICAL_DATA; then
    # Each item it could not install is on the list already (note_missing).
    install_clinical_data || true
  else
    echo "[SKIP] Clinical data (--no-clinical-data): steps 13, 30, 23, 31 and 17 need the caches or bundle."
    echo "  Install later: ./scripts/setup.sh --clinical-data \"${GENOME_DIR}\""
  fi
fi

###############################################################################
# Summary
###############################################################################
# "Setup complete!" only when nothing is missing. Otherwise each missing item
# with the command that retries it, and the exit status SETUP_PARTIAL_EXIT (3),
# so a caller can tell a partial setup from a full one.
echo ""
echo "============================================"
if [ "${#MISSING[@]}" -eq 0 ]; then
  if [ "${SETUP_EVERYTHING_SCOPE:-}" = base ]; then echo "  Base setup complete!"; else echo "  Setup complete!"; fi
else
  n=${#MISSING[@]}
  echo "  Setup finished: core ready, ${n} item$([ "$n" -eq 1 ] || echo s) missing"
fi
echo ""
echo "  Reference genome:  ${REFDIR}/"
echo "  ClinVar database:  ${CLINVARDIR}/"
echo "  Docker images:     ${PULLED} pulled, ${SKIPPED} cached"
echo "============================================"
if [ "${#MISSING[@]}" -gt 0 ]; then
  print_missing
fi
if [ "${SETUP_EVERYTHING_SCOPE:-}" = base ]; then
  if [ "${#MISSING[@]}" -gt 0 ]; then exit "$SETUP_PARTIAL_EXIT"; fi
  exit 0
fi

echo ""
echo "Next steps:"
echo "  1. Place your FASTQ/BAM/VCF in: ${GENOME_DIR}/<sample_name>/"
echo "  2. Run: export GENOME_DIR=${GENOME_DIR}"
echo "  3. Run: ./scripts/validate-setup.sh <sample_name>"
echo "  4. Run: ./scripts/run-all.sh <sample_name> <male|female>"
if [ "${#MISSING[@]}" -gt 0 ]; then
  echo ""
  echo "Setup exits with status ${SETUP_PARTIAL_EXIT}: a partial setup. Run the retry commands above, or this setup again."
  exit "$SETUP_PARTIAL_EXIT"
fi
