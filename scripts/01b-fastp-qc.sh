#!/usr/bin/env bash
# fastp QC + adapter trimming — pre-alignment quality control
# Input: paired-end FASTQ files from $GENOME_DIR/<sample>/fastq/
# Output: trimmed FASTQs + QC reports in $GENOME_DIR/<sample>/fastq_trimmed/
#
# Performs: adapter auto-detection (Illumina + BGI/MGI built-in), quality
# trimming (Q20 sliding window), polyG tail removal, length filtering (>=36bp),
# and per-read overlap-based PE adapter detection.
#
# Skip with: SKIP_TRIM=true
# The alignment step (02-alignment.sh) auto-detects trimmed FASTQs.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
R1="${SAMPLE_DIR}/fastq/${SAMPLE}_R1.fastq.gz"
R2="${SAMPLE_DIR}/fastq/${SAMPLE}_R2.fastq.gz"
OUTPUT_DIR="${SAMPLE_DIR}/fastq_trimmed"

# Skip if requested
if [ "${SKIP_TRIM:-false}" = "true" ]; then
  echo "=== fastp QC: SKIPPED (SKIP_TRIM=true) ==="
  echo "Alignment will use raw FASTQs."
  exit 0
fi

echo "=== fastp QC: ${SAMPLE} ==="
echo "R1: ${R1}"
echo "R2: ${R2}"

# Validate inputs
for f in "$R1" "$R2"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

# Skip if the trimmed pair already exists and was trimmed from the raw reads
# that are in fastq/ now: <sample>_fastp_inputs.tsv holds the input_stamp of
# R1 and R2 as fastp read them (scripts/lib/common.sh). A replaced raw FASTQ,
# or a pair trimmed before that record existed, is trimmed again. fastp writes
# into fastq_trimmed.part/, renamed to fastq_trimmed/ only when it finished,
# so a run that was killed leaves no half-written FASTQ for step 02 to align.
STAMP="${OUTPUT_DIR}/${SAMPLE}_fastp_inputs.tsv"
if [ -s "${OUTPUT_DIR}/${SAMPLE}_R1.fastq.gz" ] && [ -s "${OUTPUT_DIR}/${SAMPLE}_R2.fastq.gz" ] \
   && [ -s "${OUTPUT_DIR}/${SAMPLE}_fastp.json" ]; then
  if trimmed_reads_current "$SAMPLE_DIR" "$SAMPLE"; then
    echo "Trimmed FASTQs of these raw reads already exist in ${OUTPUT_DIR}/, skipping."
    echo "Delete them to re-run: rm -rf ${OUTPUT_DIR}"
    exit 0
  fi
  if [ -s "$STAMP" ]; then
    echo "The raw reads changed since ${OUTPUT_DIR}/ was trimmed from them: trimming again."
  else
    echo "${OUTPUT_DIR}/ has no record of the raw reads it was trimmed from (${STAMP##*/}): trimming again."
  fi
fi

PART_DIR="${OUTPUT_DIR}.part"
rm -rf "$PART_DIR"
mkdir -p "$PART_DIR"
# The stamp is taken before fastp reads the files: a file replaced while fastp
# runs then differs from it, and the next run trims again.
input_stamp "$R1" "$R2" > "${PART_DIR}/${STAMP##*/}"

# fastp: adapter trimming + QC
# Flags:
#   --detect_adapter_for_pe  Enable PE overlap-based adapter detection (off by default for PE)
#   --qualified_quality_phred 20  Mark bases below Q20 as unqualified
#   --cut_front / --cut_tail  Sliding window quality trimming from both ends
#   --cut_mean_quality 20    Window mean quality threshold
#   --length_required 36     Discard reads shorter than 36bp after trimming
#   -g / --trim_poly_g       Enable polyG tail trimming (NovaSeq/NextSeq two-color chemistry)
#   -R                       Report title (used by MultiQC for sample naming)
#   -w                       Worker threads (default 3, max 16 effective for I/O-bound work)
echo "Running fastp (adapter trimming + quality filtering)..."
run_in --cpus "${THREADS}" --memory 4g \
  "${FASTP_IMAGE}" \
  fastp \
    -i "/genome/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz" \
    -I "/genome/${SAMPLE}/fastq/${SAMPLE}_R2.fastq.gz" \
    -o "/genome/${SAMPLE}/fastq_trimmed.part/${SAMPLE}_R1.fastq.gz" \
    -O "/genome/${SAMPLE}/fastq_trimmed.part/${SAMPLE}_R2.fastq.gz" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --cut_front \
    --cut_tail \
    --cut_mean_quality 20 \
    --length_required 36 \
    -g \
    -R "${SAMPLE}" \
    -j "/genome/${SAMPLE}/fastq_trimmed.part/${SAMPLE}_fastp.json" \
    -h "/genome/${SAMPLE}/fastq_trimmed.part/${SAMPLE}_fastp.html" \
    -w "${THREADS}"

rm -rf "$OUTPUT_DIR"
mv "$PART_DIR" "$OUTPUT_DIR"

echo "=== fastp QC complete ==="
echo "Trimmed R1:   ${OUTPUT_DIR}/${SAMPLE}_R1.fastq.gz"
echo "Trimmed R2:   ${OUTPUT_DIR}/${SAMPLE}_R2.fastq.gz"
echo "JSON report:  ${OUTPUT_DIR}/${SAMPLE}_fastp.json"
echo "HTML report:  ${OUTPUT_DIR}/${SAMPLE}_fastp.html"
echo "Inputs:       ${STAMP}"
ls -lh "${OUTPUT_DIR}/${SAMPLE}_R1.fastq.gz" 2>/dev/null || true
