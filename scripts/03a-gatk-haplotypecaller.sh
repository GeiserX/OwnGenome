#!/usr/bin/env bash
# GATK HaplotypeCaller — Alternative variant caller (SNPs + indels)
# Alternative to step 03 (DeepVariant). Outputs to vcf_gatk/ to avoid conflicts.
# Input: sorted BAM + GRCh38 reference (with .dict and .fai)
# Output: VCF.gz in $GENOME_DIR/<sample>/vcf_gatk/
# Usage: ./scripts/03a-gatk-haplotypecaller.sh <sample_name> [male|female]
# Sex: with `male`, chrX and chrY are called with ploidy 1 (--sample-ploidy 1)
# outside the pseudoautosomal regions (assets/par_grch38.bed), as step 03
# does, so a male sample gets no impossible heterozygous call there: their
# units are cut at the PARs and the parts outside are called on their own.
# `female` calls every contig diploid. With no sex, the sex run-all.sh
# recorded for the sample in logs/run_status.tsv counts; with none there
# either, every contig is called diploid and a note says so.
# Reads flagged as duplicates by step 02 are skipped (GATK's default read
# filter).
#
# Scatter: HaplotypeCaller multithreads only its PairHMM, so the genome is
# split into units (chr1-22, X, Y and M one each, the other contigs together;
# or each region of INTERVALS="chr20 chr22") and SCATTER_JOBS of them (default
# THREADS/2, fewer when their 8 GB each do not fit in the container engine's
# memory) run at once, 2 CPUs and 8 GB each, then bcftools concat joins
# them in reference order. SCATTER=false runs one process over everything,
# with the same calls (for a male sample, two in turn: the haploid parts of
# chrX and chrY, and the rest). A region off its contig, or two regions that overlap,
# stop the step before GATK starts.
#
# Rerun: a finished VCF called with the same INTERVALS, SCATTER, sex (male,
# or not), BAM (the same file content), reference and image is kept (run-all.sh runs this step
# on every invocation with EXTRA_CALLERS=gatk); delete it to call again.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> [male|female]}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SEX=$(sample_sex "${2:-}")
INTERVALS=${INTERVALS:-""}

SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
ALIGN_DIR=${ALIGN_DIR:-aligned}
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
OUTPUT_DIR="${SAMPLE_DIR}/vcf_gatk"


echo "=== GATK HaplotypeCaller: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Reference: ${REF}"
echo "Threads: ${THREADS}"
if [ -n "$INTERVALS" ]; then
  echo "Intervals: ${INTERVALS}"
fi
ploidy_line "$SEX"
echo "Output: ${OUTPUT_DIR}/${SAMPLE}.vcf.gz"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai" "$REF_DICT" "$PAR_BED"; do
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
RUN_DESC="INTERVALS=${INTERVALS} scatter=${SCATTER:-true}$([ "$SEX" != male ] || echo ' sex=male') bam=${BAM} reference=${REF_FASTA} image=${GATK_IMAGE}"
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

UNITS_DIR="${OUTPUT_DIR}/scatter"
mapfile -t UNITS < <(scatter_beds "$UNITS_DIR" "$INTERVALS")
[ "${#UNITS[@]}" -gt 0 ] || { echo "ERROR: no calling units (INTERVALS='${INTERVALS}')" >&2; exit 1; }
# One unit (SCATTER=false, or one region) runs as one process with every CPU,
# also when a male sample's ploidy split makes it two: they run in turn.
ONE=false
[ "${#UNITS[@]}" -gt 1 ] || ONE=true
mapfile -t UNITS < <(ploidy_units "$SEX" "${UNITS[@]}")
[ "${#UNITS[@]}" -gt 0 ] || { echo "ERROR: no calling units after the ploidy split" >&2; exit 1; }
SCATTER_JOBS=$(scatter_jobs 8)
[ "$ONE" = false ] || SCATTER_JOBS=1

# call_unit BED: HaplotypeCaller over one unit, into scatter/<unit>.vcf.gz;
# a haploid unit (*.hap.bed, see ploidy_units) with --sample-ploidy 1.
call_unit() {
  local bed=$1 cpus=2 mem=8g
  local -a ploidy=()
  [ "$ONE" = false ] || { cpus=$THREADS mem=32g; }
  [[ "$bed" != *.hap.bed ]] || ploidy=(--sample-ploidy 1)
  run_in --cpus "$cpus" --memory "$mem" \
    "$GATK_IMAGE" \
    gatk HaplotypeCaller \
      -R "${REF_FASTA_C}" \
      -I "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
      -L "$(cpath "$bed")" \
      -O "$(cpath "${bed%.bed}.vcf.gz")" \
      --native-pair-hmm-threads "$cpus" \
      ${ploidy[@]+"${ploidy[@]}"} \
      -ERC NONE > "${bed%.bed}.log" 2>&1 \
    || { echo "ERROR: HaplotypeCaller failed on ${bed##*/}; its log:" >&2; tail -n 20 "${bed%.bed}.log" >&2; return 1; }
}

echo "=== [1/3] Running GATK HaplotypeCaller: ${#UNITS[@]} unit(s), ${SCATTER_JOBS} at a time ==="
run_parallel "$SCATTER_JOBS" call_unit "${UNITS[@]}"
PARTS=()
for u in "${UNITS[@]}"; do PARTS+=("$(cpath "${u%.bed}.vcf.gz")"); done
rm -f "${OUT}.tbi"
if [ "${#PARTS[@]}" -eq 1 ]; then
  mv -f "${UNITS[0]%.bed}.vcf.gz" "$OUT"
else
  echo "  Joining ${#PARTS[@]} units (bcftools concat)..."
  run_in --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" \
    bcftools concat -a -D -Oz -o "$(cpath "${OUT}.tmp")" "${PARTS[@]}"
  mv -f "${OUT}.tmp" "$OUT"
fi
rm -rf "$UNITS_DIR"

echo "=== [2/3] Indexing VCF with bcftools ==="
run_in --cpus 2 --memory 2g \
  "$BCFTOOLS_IMAGE" \
  bcftools index -ft "/genome/${SAMPLE}/vcf_gatk/${SAMPLE}.vcf.gz"
printf '%s\n' "$RUN_KEY" > "$RUN_FILE"

echo "=== [3/3] Variant statistics ==="
echo "VCF: ${OUTPUT_DIR}/${SAMPLE}.vcf.gz"
echo ""
echo "Quick stats:"
echo "  Total variants: $(run_in "$BCFTOOLS_IMAGE" bcftools stats "/genome/${SAMPLE}/vcf_gatk/${SAMPLE}.vcf.gz" | grep '^SN' | grep 'number of records' | awk '{print $NF}' 2>/dev/null || echo 'run bcftools stats manually')"

echo "=== GATK HaplotypeCaller complete ==="
