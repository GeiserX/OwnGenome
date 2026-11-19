#!/usr/bin/env bash
# ORA to FASTQ — Convert Illumina ORA-compressed files to standard FASTQ
# Input: .ora file(s) from Illumina DRAGEN sequencing
# Output: .fastq.gz in $GENOME_DIR/<sample>/fastq/
# NOTE: orad is a native binary (not Docker), typically at /opt/orad/bin/orad
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> <ora_reference_dir> <ora_file>}
ORA_REF=${2:?Usage: $0 <sample_name> <ora_reference_dir> <ora_file>}
ORA_FILE=${3:?Usage: $0 <sample_name> <ora_reference_dir> <ora_file>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
ORAD=${ORAD:-/opt/orad/bin/orad}
OUTPUT_DIR="${GENOME_DIR}/${SAMPLE}/fastq"

echo "=== ORA to FASTQ: ${SAMPLE} ==="
echo "ORA file: ${ORA_FILE}"
echo "ORA reference: ${ORA_REF}"
echo "Output: ${OUTPUT_DIR}/"

# Validate inputs
if [ ! -f "$ORA_FILE" ]; then
  echo "ERROR: ORA file not found: ${ORA_FILE}" >&2
  exit 1
fi

if [ ! -d "$ORA_REF" ]; then
  echo "ERROR: ORA reference directory not found: ${ORA_REF}" >&2
  exit 1
fi

if ! command -v "$ORAD" &>/dev/null; then
  echo "ERROR: orad binary not found at ${ORAD}" >&2
  echo "Download from Illumina and set ORAD=/path/to/orad" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# orad writes into a directory of its own, so this call's output can be told
# from files an earlier call left (step 01 runs once for R1 and once for R2).
# The step fails unless orad wrote at least one .fastq.gz and each one it wrote
# starts with a FASTQ record; only then are they moved into fastq/. An orad that exits 0 having
# written nothing, or a run that was killed, leaves nothing behind.
PART_DIR=$(mktemp -d "${OUTPUT_DIR}/.orad.XXXXXX")
trap 'rm -rf "$PART_DIR"' EXIT

"$ORAD" \
  --ora-reference "$ORA_REF" \
  --path "$PART_DIR" \
  "$ORA_FILE"

WROTE=()
for f in "$PART_DIR"/*.fastq.gz; do
  [ -e "$f" ] || continue
  first=$(gzip -cd "$f" 2>/dev/null | head -c 1) || true
  if [ "$first" != "@" ]; then
    echo "ERROR: orad exited 0 but its ${f##*/} holds no FASTQ record (empty or not gzip)." >&2
    echo "  Nothing was written to ${OUTPUT_DIR}/." >&2
    exit 1
  fi
  WROTE+=("$f")
done
if [ "${#WROTE[@]}" -eq 0 ]; then
  echo "ERROR: orad exited 0 but wrote no .fastq.gz for ${ORA_FILE}." >&2
  echo "  Check the ORA file, and that the ORA reference (${ORA_REF}) is the one for its run." >&2
  exit 1
fi
for f in "${WROTE[@]}"; do
  mv -f "$f" "${OUTPUT_DIR}/"
done

echo "=== ORA to FASTQ complete ==="
echo "FASTQ files: ${OUTPUT_DIR}/"
for f in "${WROTE[@]}"; do
  ls -lh "${OUTPUT_DIR}/${f##*/}"
done
