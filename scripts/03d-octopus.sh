#!/usr/bin/env bash
# Octopus — Haplotype-aware Bayesian variant caller
# Input: sorted BAM + GRCh38 reference FASTA
# Output: VCF in $GENOME_DIR/<sample>/vcf_octopus/
#
# Octopus uses a haplotype-aware, reference-panel-free Bayesian genotype model.
# Designed as a 5th caller for benchmarking alongside DeepVariant, GATK, FreeBayes,
# and Strelka2, or as a potential FreeBayes replacement (faster, lower memory).
#
# Usage: ./scripts/03d-octopus.sh <sample_name> [male|female]
# Sex: with `male`, chrX and chrY are called with ploidy 1 outside the
# pseudoautosomal regions (assets/par_grch38.bed), as step 03 does, so a male
# sample gets no impossible heterozygous call there. Octopus sets ploidy per
# contig (--contig-ploidies), not per region, so those parts are called in a
# run of their own and joined with the rest. `female` calls every contig as
# Octopus does by default: diploid, except chrY and chrM, which it calls
# haploid. With no sex, the sex run-all.sh recorded for the sample in
# logs/run_status.tsv counts; with none there either, Octopus' default holds
# and a note says so.
#
# Rerun: a finished VCF with its index, called with the same INTERVALS, sex
# (male, or not), BAM (the same file content), reference and image, is kept
# (run-all.sh runs this step on every invocation with EXTRA_CALLERS=octopus);
# delete it to call again.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> [male|female]}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SEX=$(sample_sex "${2:-}")
ALIGN_DIR=${ALIGN_DIR:-aligned}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
OUTPUT_DIR="${SAMPLE_DIR}/vcf_octopus"

echo "=== Octopus: ${SAMPLE} ==="
echo "BAM: ${BAM}"
case "$SEX" in
  male) ploidy_line male ;;
  female) echo "Ploidy: female, every contig diploid but chrY and chrM (Octopus' default)" ;;
  *) echo "Ploidy: no sex given, chrX is called diploid (chrY and chrM haploid, Octopus' default). Pass male or female as the second argument." ;;
esac

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai" "$PAR_BED"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

# A finished VCF is reused only when it is complete (a whole BGZF file with a
# VCF header), its index is there, and it was called the way this run would
# call it: RUN_FILE, written last, holds INTERVALS, sex=male for a male
# sample, the BAM's path, the reference and the image on its first line, then
# the BAM's input_stamp (size, modification time and the sha256 of its first
# and last MiB). A file cut
# short by a killed run, a VCF without its .tbi, a VCF of other INTERVALS and
# one called from another BAM, even under the same name, are called again. A
# female sample and one of no sex are called the same way, so they share a
# record.
OUT="${OUTPUT_DIR}/${SAMPLE}.vcf.gz"
RUN_FILE="${OUTPUT_DIR}/${SAMPLE}.run"
RUN_DESC="INTERVALS=${INTERVALS:-}$([ "$SEX" != male ] || echo ' sex=male') bam=${BAM} reference=${REF_FASTA} image=${OCTOPUS_IMAGE}"
BAM_STAMP=$(input_stamp "$BAM")
RUN_KEY="${RUN_DESC}
${BAM_STAMP}"
if have_output "$OUT" "${OUT}.tbi" && [ -f "$RUN_FILE" ] && [ "$(cat "$RUN_FILE")" = "$RUN_KEY" ]; then
  echo "Octopus output already exists: ${OUT} (${RUN_DESC}, the same BAM), skipping."
  echo "Delete to re-run: rm -rf ${OUTPUT_DIR}"
  exit 0
fi
if [ -e "$OUT" ]; then
  echo "Calling again: ${OUT} is unfinished, has no index or was not called with ${RUN_DESC} on this BAM."
fi

mkdir -p "$OUTPUT_DIR"
# The index of an older VCF must never sit beside the new one.
rm -f "$RUN_FILE" "${OUT}.tbi"

# run_octopus OUT_C ARG...: Octopus in germline mode (default) on the BAM, into
# the container path OUT_C. Flags: -R reference FASTA, -I input BAM, -o output
# VCF, --threads worker threads, then ARG... (regions, ploidies).
# -w: Octopus makes its octopus-temp directory in the working directory, and
# the image sets none, so it would try / as an unprivileged user.
run_octopus() {
  local out=$1
  shift
  run_in -w "/genome/${SAMPLE}/vcf_octopus" --cpus "${THREADS}" --memory 16g \
    "${OCTOPUS_IMAGE}" \
    octopus \
      -R "${REF_FASTA_C}" \
      -I "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
      -o "$out" \
      --threads "${THREADS}" \
      "$@"
}

echo "Running Octopus (this takes 2-4 hours for 30X WGS)..."
if [ "$SEX" != male ]; then
  # Restrict to specific regions if INTERVALS is set (e.g., INTERVALS=chr22 for testing).
  # Expanded as ${REGION_ARGS[@]+...}: bash before 4.4 calls an empty array
  # unbound under set -u.
  REGION_ARGS=()
  if [ -n "${INTERVALS:-}" ]; then
    echo "Restricting to region: ${INTERVALS}"
    # One word per region: Octopus takes --regions as a list of words and
    # refuses "chr1:1-100 chr2" as one region.
    read -r -a IV <<< "${INTERVALS}"
    REGION_ARGS=(--regions "${IV[@]}")
  fi
  run_octopus "/genome/${SAMPLE}/vcf_octopus/${SAMPLE}.vcf.gz" ${REGION_ARGS[@]+"${REGION_ARGS[@]}"}
else
  # Two runs: the regions to call (INTERVALS, or every contig of the
  # reference; scatter_beds checks them) cut at the PARs by ploidy_units. The
  # haploid parts of chrX and chrY run with chrX=1 chrY=1; the rest with
  # MT=1 chrM=1, Octopus' default without its chrY=1, so chrY's PARs are
  # diploid like chrX's, as in step 03. The region files hold chrom:begin-end
  # with the BED's zero-based begin, which is how Octopus reads them.
  UNITS_DIR="${OUTPUT_DIR}/ploidy"
  mapfile -t UNITS < <(SCATTER=false scatter_beds "$UNITS_DIR" "${INTERVALS:-}")
  [ "${#UNITS[@]}" -gt 0 ] || { echo "ERROR: no calling regions (INTERVALS='${INTERVALS:-}')" >&2; exit 1; }
  SPLIT=$(ploidy_units male "${UNITS[@]}") || { echo "ERROR: could not cut the calling regions at ${PAR_BED}." >&2; exit 1; }
  [ -n "$SPLIT" ] || { echo "ERROR: no calling regions after the ploidy split" >&2; exit 1; }
  mapfile -t UNITS <<< "$SPLIT"
  PARTS=()
  for u in "${UNITS[@]}"; do
    awk -F'\t' '{ print $1 ":" $2 "-" $3 }' "$u" > "${u%.bed}.regions"
    if [[ "$u" == *.hap.bed ]]; then
      echo "  haploid run: $(wc -l < "$u" | tr -d ' ') region(s) of chrX and chrY outside the PARs"
      PLOIDIES=(--contig-ploidies chrX=1 chrY=1)
    else
      echo "  diploid run: $(wc -l < "$u" | tr -d ' ') region(s)"
      PLOIDIES=(--contig-ploidies MT=1 chrM=1)
    fi
    run_octopus "$(cpath "${u%.bed}.vcf.gz")" --regions-file "$(cpath "${u%.bed}.regions")" "${PLOIDIES[@]}"
    if ! have_output "${u%.bed}.vcf.gz"; then
      echo "ERROR: Octopus finished without a complete ${u%.bed}.vcf.gz." >&2
      exit 1
    fi
    [ -s "${u%.bed}.vcf.gz.tbi" ] || run_in --cpus 1 --memory 1g "$BCFTOOLS_IMAGE" \
      bcftools index -f -t "$(cpath "${u%.bed}.vcf.gz")"
    PARTS+=("$(cpath "${u%.bed}.vcf.gz")")
  done
  if [ "${#PARTS[@]}" -eq 1 ]; then
    mv -f "${UNITS[0]%.bed}.vcf.gz" "$OUT"
  else
    echo "  Joining the two runs (bcftools concat)..."
    run_in --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" \
      bcftools concat -a -D -Oz -o "$(cpath "${OUT}.tmp")" "${PARTS[@]}"
    mv -f "${OUT}.tmp" "$OUT"
  fi
  rm -rf "$UNITS_DIR"
fi

# Octopus indexes a .vcf.gz it writes; index it here when it did not.
if have_output "$OUT" && [ ! -s "${OUT}.tbi" ]; then
  run_in --cpus 1 --memory 1g "$BCFTOOLS_IMAGE" \
    bcftools index -f -t "/genome/${SAMPLE}/vcf_octopus/${SAMPLE}.vcf.gz"
fi
if ! have_output "$OUT" "${OUT}.tbi"; then
  echo "ERROR: Octopus finished without a complete ${OUT} and its .tbi." >&2
  exit 1
fi
printf '%s\n' "$RUN_KEY" > "$RUN_FILE"

echo "=== Octopus complete ==="
echo "VCF: ${OUTPUT_DIR}/${SAMPLE}.vcf.gz"
ls -lh "${OUTPUT_DIR}/${SAMPLE}.vcf.gz" 2>/dev/null || true
