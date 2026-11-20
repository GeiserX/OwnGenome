#!/usr/bin/env bash
# TelomereHunter2: relative telomere content from a WGS BAM
# Output: tel_content (intratelomeric reads per million reads with 48-52% GC).
# A higher value means more telomeric sequence. It is a relative research
# measure for comparing samples from one platform, not an age.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}   # common.sh defaults to 8
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
ALIGN_DIR=${ALIGN_DIR:-aligned}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
OUTPUT_DIR="${SAMPLE_DIR}/telomere/${SAMPLE}"

echo "=== TelomereHunter: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Output: ${OUTPUT_DIR}"
echo "WARNING: This reads the entire BAM (~30-40GB). Takes 30-60 minutes."

for f in "$BAM" "${BAM}.bai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

# TelomereHunter2 sorts telomeric reads into intratelomeric, subtelomeric and
# junction classes by chromosome band. The BAM is GRCh38, so it always gets
# UCSC's GRCh38 bands: the copy setup.sh installs, or, when that is missing,
# the one the repository ships (assets/reference), installed here first.
if ! BANDS=$(data_file cytoband); then
  echo "The GRCh38 chromosome bands are not installed (${BANDS}); installing the copy in assets/reference."
  if ! install_data_file cytoband; then
    echo "ERROR: could not install the GRCh38 chromosome bands to ${BANDS}" >&2
    exit 1
  fi
fi
echo "Chromosome bands: ${BANDS}"

# -c: TelomereHunter2's worker count, THREADS as the container's CPU cap.
# --ulimit: it writes one BAM per contig and merges them at once. A BAM with
# decoy contigs (the HG002 test BAM has 2,580) needs more than the 1,024 open
# files a container gets by default ("Too many open files").
# --plotNone: nothing in the pipeline reads its plots.
run_in \
  --cpus "$THREADS" --memory 4g --ulimit nofile=65536:65536 \
  "${TELOMEREHUNTER_IMAGE}" \
  telomerehunter2 \
    -ibt "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
    -o "/genome/${SAMPLE}/telomere/${SAMPLE}" \
    -p "$SAMPLE" \
    -c "$THREADS" \
    --plotNone \
    -b "$(cpath "$BANDS")"

echo "=== TelomereHunter complete ==="
SUMMARY="${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_summary.tsv"
if [ -f "$SUMMARY" ]; then
  echo "Summary: $SUMMARY"
  TEL_CONTENT=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "tel_content") c = i; next} NR == 2 && c {print $c}' "$SUMMARY")
  echo "Telomere content: ${TEL_CONTENT}"
fi
