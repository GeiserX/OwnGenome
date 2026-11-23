#!/usr/bin/env bash
# mosdepth — fast per-base and per-region coverage statistics from BAM
# Input: sorted BAM + BAI from $GENOME_DIR/<sample>/aligned/
# Output: coverage distributions, thresholds, summary in $GENOME_DIR/<sample>/mosdepth/
#
# Complements indexcov (step 16): indexcov reads only the BAM index for a quick
# estimate; mosdepth reads actual alignments for precise per-base depth.
#
# A second pass measures how well each gene of assets/callability_genes.tsv
# (CALLABILITY_GENES names another table) could be read: for each gene on
# chr1-22, X and Y, the bases with 10 or more reads of mapping quality 20 or
# more, in <sample>_callability.bed.gz, and the same for each coding exon of
# those genes in assets/callability_exons.tsv (CALLABILITY_EXONS names another
# table), in <sample>_callability_exons.bed.gz. The report lists the genes
# where fewer than 90% of the bases, or of the bases of one coding exon, reach
# that (docs/16b-mosdepth.md).
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
OUTPUT_DIR="${SAMPLE_DIR}/mosdepth"
CALLABILITY_GENES=${CALLABILITY_GENES:-${PGP_ROOT}/assets/callability_genes.tsv}
CALLABILITY_EXONS=${CALLABILITY_EXONS:-${PGP_ROOT}/assets/callability_exons.tsv}
CALLABLE="${OUTPUT_DIR}/${SAMPLE}_callability.bed.gz"
CALLABLE_EXONS="${OUTPUT_DIR}/${SAMPLE}_callability_exons.bed.gz"

echo "=== mosdepth Coverage: ${SAMPLE} ==="
echo "BAM: ${BAM}"

# Validate inputs
if [ ! -f "$BAM" ]; then
  echo "ERROR: BAM not found: ${BAM}" >&2
  exit 1
fi
if [ ! -f "${BAM}.bai" ]; then
  echo "ERROR: BAM index not found: ${BAM}.bai" >&2
  echo "Create it with: samtools index ${BAM}" >&2
  exit 1
fi

if [ ! -s "$CALLABILITY_GENES" ]; then
  echo "ERROR: gene table not found or empty: ${CALLABILITY_GENES}" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# gene_bed and exon_bed print the BED rows (chrom, start, end, name) of the genes
# of the gene table with a position on chr1-22, X or Y (not the ALT-only ones;
# the report reads them from the table), and of the coding exons of those genes
# in the exon table (none when it is missing or names none of them).
gene_bed() {
  awk -F'\t' -v OFS='\t' '!/^#/ && $1 != "gene" && $4 ~ /^[0-9]+$/ {print $3, $4, $5, $1}' "$CALLABILITY_GENES"
}
exon_bed() {
  [ -s "$CALLABILITY_EXONS" ] || return 0
  gene_bed | awk -F'\t' -v OFS='\t' 'NR == FNR {g[$4] = 1; next}
    !/^#/ && $1 != "exon" && ($2 in g) {print $4, $5, $6, $1}' - "$CALLABILITY_EXONS"
}

# callability_pass: mosdepth over the genes of the table, reads of MAPQ 20 or
# more only, the bases at 10x or more per gene (the same command as the
# MOSDEPTH module). The coding exons of those genes ride in the same pass, so
# the BAM is read once: each row is measured on its own, and the rows named
# GENE:exonN go to the exon table, the others to the gene table. An exon row
# in the gene table would read as a gene, so the two stay apart.
# No --fast-mode here: without it mosdepth counts the bases where the two mates
# of a pair overlap once, not twice, and does not count the bases of a deletion
# in the CIGAR. With it, a short-insert pair or a deletion lifts bases to 10x
# that fewer than 10 fragments read, and a gene reads as read when it was not.
callability_pass() {
  local work="${OUTPUT_DIR}/callability.part"
  rm -rf "$work"
  mkdir -p "$work"
  gene_bed > "${work}/genes.bed"
  exon_bed > "${work}/exons.bed"
  echo "Gene callability (MAPQ >= 20, 10x) over $(grep -c . "${work}/genes.bed") genes of ${CALLABILITY_GENES}" \
    "and $(grep -c . "${work}/exons.bed") of their coding exons..."
  cat "${work}/genes.bed" "${work}/exons.bed" > "${work}/regions.bed"
  run_in --cpus "${THREADS}" --memory 4g \
    "${MOSDEPTH_IMAGE}" \
    mosdepth \
      --by "/genome/${SAMPLE}/mosdepth/callability.part/regions.bed" \
      --mapq 20 \
      --no-per-base \
      --threads "${THREADS}" \
      --thresholds 10 \
      "/genome/${SAMPLE}/mosdepth/callability.part/c" \
      "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
  local t="${work}/c.thresholds.bed.gz"
  if [ ! -f "$t" ]; then
    echo "WARNING: mosdepth wrote no gene table; the report says the gene callability was not measured." >&2
  elif [ -s "${work}/exons.bed" ]; then
    # The header goes to both tables.
    gzip -dc "$t" | awk -F'\t' -v g="${work}/g" -v e="${work}/e" \
      'NR == 1 {print > g; print > e; next} $4 ~ /:exon[0-9]+$/ {print > e; next} {print > g}'
    gzip -c "${work}/g" > "${work}/g.gz"
    gzip -c "${work}/e" > "${work}/e.gz"
    mv "${work}/e.gz" "$CALLABLE_EXONS"
    mv "${work}/g.gz" "$CALLABLE"
  else
    mv "$t" "$CALLABLE"
  fi
  rm -rf "$work"
}

# Skip if output already exists. mosdepth writes the summary last, so a
# non-empty summary means the run finished. The gene pass runs on its own
# when only its tables are missing (a folder from before one existed): the
# gene table, or the exon table when the exon table names a gene of it.
if have_output "${OUTPUT_DIR}/${SAMPLE}.mosdepth.summary.txt"; then
  echo "mosdepth output already exists in ${OUTPUT_DIR}/, skipping."
  if ! have_output "$CALLABLE" || { ! have_output "$CALLABLE_EXONS" && [ -n "$(exon_bed | head -1)" ]; }; then
    callability_pass
  fi
  echo "Delete to re-run: rm -rf ${OUTPUT_DIR}"
  exit 0
fi

# mosdepth flags:
#   --by 500        Window size for per-region coverage (500bp bins; WGS default)
#   --by <BED>      Capture BED for on-target coverage (WES; set CAPTURE_BED env var)
#   --fast-mode     Skip the CIGAR and the mate-overlap correction: faster, and fine
#                   for these depth summaries (the gene pass above runs without it)
#   --no-per-base   Omit the large per-base BED (saves ~3-5 GB per WGS sample)
#   --threads       Decompression threads (mosdepth uses 1 main + N decompression)
#   --thresholds    Report fraction of bases at these coverage thresholds

# Use capture BED for WES on-target coverage, or 500bp bins for WGS
if [ -n "${CAPTURE_BED:-}" ]; then
  if [ ! -f "${CAPTURE_BED}" ]; then
    echo "ERROR: CAPTURE_BED not found: ${CAPTURE_BED}" >&2
    exit 1
  fi
  # Compute container-relative path for the BED file
  CAPTURE_BED_REL="${CAPTURE_BED#"${GENOME_DIR}/"}"
  BY_FLAG="/genome/${CAPTURE_BED_REL}"
  echo "Using capture BED for WES on-target coverage: ${CAPTURE_BED}"
else
  BY_FLAG="500"
fi

echo "Computing coverage statistics..."
run_in --cpus "${THREADS}" --memory 4g \
  "${MOSDEPTH_IMAGE}" \
  mosdepth \
    --by "${BY_FLAG}" \
    --fast-mode \
    --no-per-base \
    --threads "${THREADS}" \
    --thresholds 1,5,10,15,20,30,50 \
    "/genome/${SAMPLE}/mosdepth/${SAMPLE}" \
    "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"

rm -f "$CALLABLE" "$CALLABLE_EXONS"
callability_pass

echo "=== mosdepth complete ==="
echo "Summary:        ${OUTPUT_DIR}/${SAMPLE}.mosdepth.summary.txt"
echo "Distribution:   ${OUTPUT_DIR}/${SAMPLE}.mosdepth.global.dist.txt"
echo "Region coverage: ${OUTPUT_DIR}/${SAMPLE}.regions.bed.gz"
echo "Thresholds:     ${OUTPUT_DIR}/${SAMPLE}.thresholds.bed.gz"
echo "Gene callability: ${CALLABLE}"
[ ! -f "$CALLABLE_EXONS" ] || echo "Exon callability: ${CALLABLE_EXONS}"

# Print summary
if [ -f "${OUTPUT_DIR}/${SAMPLE}.mosdepth.summary.txt" ]; then
  echo ""
  echo "Coverage summary:"
  head -5 "${OUTPUT_DIR}/${SAMPLE}.mosdepth.summary.txt"
fi
