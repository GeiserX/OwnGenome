#!/usr/bin/env bash
# FreeBayes — Alternative variant caller (SNPs + indels)
# Alternative to step 03 (DeepVariant). Outputs to vcf_freebayes/ to avoid conflicts.
# Input: sorted BAM + GRCh38 reference (.fasta + .fai)
# Output: VCF.gz in $GENOME_DIR/<sample>/vcf_freebayes/
# Usage: ./scripts/03b-freebayes.sh <sample_name> [male|female]
# Sex: with `male`, chrX and chrY are called with ploidy 1 (--ploidy 1)
# outside the pseudoautosomal regions (assets/par_grch38.bed), as step 03
# does, so a male sample gets no impossible heterozygous call there: their
# units are cut at the PARs and the parts outside are called on their own.
# `female` calls every contig diploid. With no sex, the sex run-all.sh
# recorded for the sample in logs/run_status.tsv counts; with none there
# either, every contig is called diploid and a note says so.
# Runtime: ~9 hours single-threaded for 30X WGS in one process; scattered
#   (below) about that divided by SCATTER_JOBS
# Memory: peaks at ~13 GB for full genome; needs 32 GB allocation for safety margin
#
# Scatter: FreeBayes runs on one thread, so the genome is split into units
# (chr1-22, X, Y and M one each, the other contigs together; or each region of
# INTERVALS="chr20 chr22") and SCATTER_JOBS of them (default THREADS/2, fewer
# when their 8 GB each do not fit in the container engine's memory) run at
# once, 1 CPU and 8 GB each. Their records are joined into the raw VCF and
# sorted as before. SCATTER=false runs one process over everything (for a
# male sample, two in turn: the haploid parts of chrX and chrY, and the
# rest). A region
# off its contig, or two regions that overlap, stop the step before FreeBayes
# starts.
#
# Rerun: a finished VCF called with the same INTERVALS, SCATTER, sex (male,
# or not), BAM (the same file content), reference and image is kept (run-all.sh runs this step
# on every invocation with EXTRA_CALLERS=freebayes); delete it to call again.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> [male|female]}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SEX=$(sample_sex "${2:-}")
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
ALIGN_DIR=${ALIGN_DIR:-aligned}
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
OUTPUT_DIR="${SAMPLE_DIR}/vcf_freebayes"
INTERVALS=${INTERVALS:-""}

echo "=== FreeBayes: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Reference: ${REF}"
echo "Output: ${OUTPUT_DIR}/${SAMPLE}.vcf.gz"
if [ -n "$INTERVALS" ]; then
  echo "Region: ${INTERVALS}"
fi
ploidy_line "$SEX"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai" "$PAR_BED"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

# A finished VCF is reused only when it was called the way this run would
# call it: RUN_FILE, written last, holds INTERVALS, SCATTER, sex=male for a
# male sample, the BAM's path, the reference and the image on its first line, then the BAM's input_stamp
# (size, modification time and the sha256 of its first and last MiB), so a
# realigned or replaced BAM is called again even under the same name, size
# and time. A subset VCF never stands in for a whole-genome request, and
# SCATTER=false after a scattered run calls again in one process. A female
# sample and one of no sex are called the same way, so they share a record.
OUT="${OUTPUT_DIR}/${SAMPLE}.vcf.gz"
RUN_FILE="${OUTPUT_DIR}/${SAMPLE}.run"
RUN_DESC="INTERVALS=${INTERVALS} scatter=${SCATTER:-true}$([ "$SEX" != male ] || echo ' sex=male') bam=${BAM} reference=${REF_FASTA} image=${FREEBAYES_IMAGE}"
BAM_STAMP=$(input_stamp "$BAM")
RUN_KEY="${RUN_DESC}
${BAM_STAMP}"
if have_output "$OUT" "${OUT}.tbi" && [ -f "$RUN_FILE" ] && [ "$(cat "$RUN_FILE")" = "$RUN_KEY" ]; then
  echo "Output already exists: ${OUT} (${RUN_DESC}, the same BAM)"
  echo "Skipping. Delete the file to re-run."
  exit 0
fi
if [ -e "$OUT" ]; then
  echo "Calling again: ${OUT} is unfinished or was not called with ${RUN_DESC} on this BAM."
fi
rm -f "$RUN_FILE"

# Step 1: Run FreeBayes over each unit (one thread each, unsorted VCF)
UNITS_DIR="${OUTPUT_DIR}/scatter"
mapfile -t UNITS < <(scatter_beds "$UNITS_DIR" "$INTERVALS")
[ "${#UNITS[@]}" -gt 0 ] || { echo "ERROR: no calling units (INTERVALS='${INTERVALS}')" >&2; exit 1; }
# One unit (SCATTER=false, or one region) runs as one process with 32 GB,
# also when a male sample's ploidy split makes it two: they run in turn.
ONE=false
[ "${#UNITS[@]}" -gt 1 ] || ONE=true
# The cut is kept only when it succeeds: one that fails part way must not
# leave the units before the failure to be called and cached as the result.
SPLIT=$(ploidy_units "$SEX" "${UNITS[@]}") || { echo "ERROR: could not cut the calling units at ${PAR_BED}." >&2; exit 1; }
[ -n "$SPLIT" ] || { echo "ERROR: no calling units after the ploidy split" >&2; exit 1; }
mapfile -t UNITS <<< "$SPLIT"
SCATTER_JOBS=$(scatter_jobs 8)
[ "$ONE" = false ] || SCATTER_JOBS=1
echo "Running FreeBayes: ${#UNITS[@]} unit(s), ${SCATTER_JOBS} at a time (this may take several hours for 30X WGS)..."

# call_unit BED: FreeBayes over one unit, into scatter/<unit>.vcf. Through a
# temporary name: a FreeBayes that fails leaves no partial VCF behind. A
# haploid unit (*.hap.bed, see ploidy_units) runs with --ploidy 1: every row
# of it is haploid, so the one ploidy holds for the whole run (a --cnv-map
# would say the same per region, but its rows must name the read group's
# sample).
call_unit() {
  local bed=$1 mem=8g
  local -a ploidy=()
  [ "$ONE" = false ] || mem=32g
  [[ "$bed" != *.hap.bed ]] || ploidy=(--ploidy 1)
  atomic_out "${bed%.bed}.vcf" run_in \
    --cpus 1 --memory "$mem" \
    "${FREEBAYES_IMAGE}" \
    freebayes -f "${REF_FASTA_C}" --targets "$(cpath "$bed")" ${ploidy[@]+"${ploidy[@]}"} \
      "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
}
run_parallel "$SCATTER_JOBS" call_unit "${UNITS[@]}"
# The raw VCF: the first unit's header, then every unit's records.
{
  grep '^#' "${UNITS[0]%.bed}.vcf"
  for u in "${UNITS[@]}"; do grep -v '^#' "${u%.bed}.vcf" || true; done
} > "${OUTPUT_DIR}/${SAMPLE}_raw.vcf.tmp"
mv -f "${OUTPUT_DIR}/${SAMPLE}_raw.vcf.tmp" "${OUTPUT_DIR}/${SAMPLE}_raw.vcf"
rm -rf "$UNITS_DIR"

# Step 2: Sort and compress in one bcftools call (no pipe whose first half can
# fail unseen), with its temporary files in the output directory, then index.
# The sorted VCF is written under a .tmp name and renamed when the sort
# succeeded, so a sort that dies half way leaves no truncated ${SAMPLE}.vcf.gz.
# The raw VCF is removed only after both succeeded.
echo "Sorting and compressing VCF..."
SORTED="$OUT"
if ! run_in \
  --cpus 4 --memory 4g \
  "${BCFTOOLS_IMAGE}" \
  bcftools sort -Oz -o "/genome/${SAMPLE}/vcf_freebayes/${SAMPLE}.vcf.gz.tmp" \
    -T "/genome/${SAMPLE}/vcf_freebayes/sort-tmp" \
    "/genome/${SAMPLE}/vcf_freebayes/${SAMPLE}_raw.vcf"; then
  rm -f "${SORTED}.tmp"
  echo "ERROR: bcftools sort failed; the raw VCF is kept: ${OUTPUT_DIR}/${SAMPLE}_raw.vcf" >&2
  exit 1
fi
# The index of an older VCF goes first: it must never sit beside the new one.
rm -f "${SORTED}.tbi"
mv -f "${SORTED}.tmp" "$SORTED"

echo "Indexing VCF..."
run_in \
  --cpus 1 --memory 1g \
  "${BCFTOOLS_IMAGE}" \
  bcftools index -f -t "/genome/${SAMPLE}/vcf_freebayes/${SAMPLE}.vcf.gz"

# Clean up raw unsorted VCF
rm -f "${OUTPUT_DIR}/${SAMPLE}_raw.vcf"
printf '%s\n' "$RUN_KEY" > "$RUN_FILE"

echo "=== FreeBayes complete ==="
echo "VCF: ${OUTPUT_DIR}/${SAMPLE}.vcf.gz"
echo ""
echo "Quick stats:"
echo "  Total variants: $(run_in "${BCFTOOLS_IMAGE}" bcftools stats "/genome/${SAMPLE}/vcf_freebayes/${SAMPLE}.vcf.gz" | grep '^SN' | grep 'number of records' | awk '{print $NF}' 2>/dev/null || echo 'run bcftools stats manually')"
echo ""
echo "NOTE: FreeBayes tends to call more variants than DeepVariant (higher sensitivity, more false positives)."
echo "Consider running bcftools filter or vcffilter for quality filtering."
