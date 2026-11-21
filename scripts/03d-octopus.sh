#!/usr/bin/env bash
# Octopus — Haplotype-aware Bayesian variant caller
# Input: sorted BAM + GRCh38 reference FASTA
# Output: VCF in $GENOME_DIR/<sample>/vcf_octopus/
#
# Octopus uses a haplotype-aware, reference-panel-free Bayesian genotype model.
# Designed as a 5th caller for benchmarking alongside DeepVariant, GATK, FreeBayes,
# and Strelka2, or as a potential FreeBayes replacement (faster, lower memory).
#
# Rerun: a finished VCF with its index, called with the same INTERVALS, BAM
# (the same file content), reference and image, is kept (run-all.sh runs this
# step on every invocation with EXTRA_CALLERS=octopus); delete it to call
# again.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
ALIGN_DIR=${ALIGN_DIR:-aligned}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
OUTPUT_DIR="${SAMPLE_DIR}/vcf_octopus"

echo "=== Octopus: ${SAMPLE} ==="
echo "BAM: ${BAM}"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

# A finished VCF is reused only when it is complete (a whole BGZF file with a
# VCF header), its index is there, and it was called the way this run would
# call it: RUN_FILE, written last, holds INTERVALS, the BAM's path, the
# reference and the image on its first line, then the BAM's input_stamp (size,
# modification time and the sha256 of its first and last MiB). A file cut
# short by a killed run, a VCF without its .tbi, a VCF of other INTERVALS and
# one called from another BAM, even under the same name, are called again.
OUT="${OUTPUT_DIR}/${SAMPLE}.vcf.gz"
RUN_FILE="${OUTPUT_DIR}/${SAMPLE}.run"
RUN_DESC="INTERVALS=${INTERVALS:-} bam=${BAM} reference=${REF_FASTA} image=${OCTOPUS_IMAGE}"
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

# Restrict to specific regions if INTERVALS is set (e.g., INTERVALS=chr22 for testing).
# Expanded as ${REGION_ARGS[@]+...}: bash before 4.4 calls an empty array
# unbound under set -u.
REGION_ARGS=()
if [ -n "${INTERVALS:-}" ]; then
  echo "Restricting to region: ${INTERVALS}"
  REGION_ARGS=(--regions "${INTERVALS}")
fi

# Octopus in germline mode (default)
# Flags:
#   -R           Reference FASTA
#   -I           Input BAM
#   -o           Output VCF
#   --threads    Worker threads
#   --regions    Restrict to regions (optional, for testing)
echo "Running Octopus (this takes 2-4 hours for 30X WGS)..."
# -w: Octopus makes its octopus-temp directory in the working directory, and
# the image sets none, so it would try / as an unprivileged user.
run_in -w "/genome/${SAMPLE}/vcf_octopus" --cpus "${THREADS}" --memory 16g \
  "${OCTOPUS_IMAGE}" \
  octopus \
    -R "${REF_FASTA_C}" \
    -I "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
    -o "/genome/${SAMPLE}/vcf_octopus/${SAMPLE}.vcf.gz" \
    --threads "${THREADS}" \
    ${REGION_ARGS[@]+"${REGION_ARGS[@]}"}

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
