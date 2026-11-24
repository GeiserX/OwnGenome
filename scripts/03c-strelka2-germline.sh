#!/usr/bin/env bash
# Strelka2 — Alternative germline small variant caller (SNVs + indels)
# Alternative to step 03 (DeepVariant). Outputs to vcf_strelka2/ to avoid conflicts.
# Strelka2 is a SMALL VARIANT caller (SNVs + indels up to ~49 bp), NOT a structural
# variant caller. It complements Manta, which handles SVs (>50 bp).
# Input: sorted BAM + GRCh38 reference
# Output: germline variants VCF in $GENOME_DIR/<sample>/vcf_strelka2/
# Runtime: ~1-2 hours per 30X genome
# Usage: ./scripts/03c-strelka2-germline.sh <sample_name> [male|female]
# Sex: with `male`, chrX and chrY are called with ploidy 1 outside the
# pseudoautosomal regions (assets/par_grch38.bed), as step 03 does, so a male
# sample gets no impossible heterozygous call there. Strelka2 takes the
# ploidy as a VCF (--ploidy): <sample>.ploidy.vcf.gz in the output directory,
# one CN=1 record per haploid stretch of chrX and chrY, its column named
# after the BAM's read group sample (@RG SM). `female` calls every contig
# diploid. With no sex, the sex run-all.sh recorded for the sample in
# logs/run_status.tsv counts; with none there either, every contig is called
# diploid and a note says so.
#
# NOTE: Strelka2's scoring model is optimized for BWA-MEM alignments.
# When used with minimap2-aligned BAMs, SNP precision may be reduced due to:
#   - Missing XS (suboptimal alignment score) tags that BWA-MEM produces
#   - Different AS (alignment score) scaling between minimap2 and BWA-MEM
# For best results, use BWA-MEM2 alignments (scripts/02a-alignment-bwamem2.sh).
#
# Reference: Kim et al. Strelka2: fast and accurate calling of germline and
# somatic variants. Nature Methods (2018). https://doi.org/10.1038/s41592-018-0051-x
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
OUTPUT_DIR="${SAMPLE_DIR}/vcf_strelka2"

echo "=== Strelka2 Germline Calling: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Reference: ${REF}"
echo "Output: ${OUTPUT_DIR}"
echo "Threads: ${THREADS}"
ploidy_line "$SEX"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai" "$PAR_BED"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done


# A run is reused only when it was configured on this BAM: RUN_FILE, written
# right after the configure script, holds sex=male for a male sample, the
# BAM's path, the reference and the image on its first line, then the BAM's
# input_stamp (size, modification time and the sha256 of its first and last
# MiB). With it a finished run is not repeated (the configure script refuses
# a runDir that already holds a workflow, so a second run used to fail) and
# an interrupted one is resumed.
# Without it, or with another BAM (realigned or replaced, even under the same
# name), the old run is removed and the BAM is called again; so is a run of
# the other ploidy. A female sample and one of no sex are called the same
# way, so they share a record.
RUN_FILE="${OUTPUT_DIR}/${SAMPLE}.run"
RUN_DESC="$([ "$SEX" != male ] || echo 'sex=male ')bam=${BAM} reference=${REF_FASTA} image=${STRELKA_IMAGE}"
BAM_STAMP=$(input_stamp "$BAM")
RUN_KEY="${RUN_DESC}
${BAM_STAMP}"
SAME_RUN=false
if [ -f "$RUN_FILE" ] && [ "$(cat "$RUN_FILE")" = "$RUN_KEY" ]; then SAME_RUN=true; fi
if [ "$SAME_RUN" = true ] && [ -f "${OUTPUT_DIR}/results/variants/variants.vcf.gz" ] && [ -f "${OUTPUT_DIR}/results/variants/variants.vcf.gz.tbi" ]; then
  echo "Strelka2 already done for ${SAMPLE}: ${OUTPUT_DIR}/results/variants/variants.vcf.gz"
  echo "To run it again, delete ${OUTPUT_DIR}/ first."
  exit 0
fi

# Step 1: Configure Strelka2 germline workflow (skipped when an interrupted
# run of this BAM can be resumed)
if [ "$SAME_RUN" = true ] && [ -f "${OUTPUT_DIR}/runWorkflow.py" ]; then
  echo "[1/2] Found an unfinished Strelka2 run in ${OUTPUT_DIR}/; resuming it."
else
  if [ -d "$OUTPUT_DIR" ]; then
    # --root: a directory left by an earlier version of this script belongs to root.
    echo "Removing ${OUTPUT_DIR}/: it holds no run configured on ${RUN_DESC} with this BAM..."
    run_in --root "$STRELKA_IMAGE" \
      rm -rf "/genome/${SAMPLE}/vcf_strelka2"
  fi
  mkdir -p "$OUTPUT_DIR"
  PLOIDY_ARGS=()
  if [ "$SEX" = male ]; then
    # The haploid stretches: chrX and chrY of the reference, cut at the PARs.
    awk -F'\t' -v OFS='\t' '$1 == "chrX" || $1 == "chrY" { print $1, 0, $2 }' "${REF}.fai" > "${OUTPUT_DIR}/xy.bed"
    HAP=$(ploidy_units male "${OUTPUT_DIR}/xy.bed" | grep '\.hap\.bed$' || true)
    if [ -z "$HAP" ]; then
      echo "NOTE: the reference has no chrX or chrY outside the PARs; every contig is called diploid."
    else
      # Strelka2 matches the ploidy VCF's sample column to the BAM's read
      # group sample, so the column takes that name.
      SM=$(run_in --cpus 1 --memory 1g "$SAMTOOLS_IMAGE" \
        samtools view -H "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
        | awk -F'\t' '$1 == "@RG" { for (i = 2; i <= NF; i++) if ($i ~ /^SM:/) print substr($i, 4) }' | sort -u)
      if [ -z "$SM" ] || [ "$(printf '%s\n' "$SM" | wc -l | tr -d ' ')" -ne 1 ]; then
        echo "ERROR: the BAM's read groups name $(printf '%s\n' "$SM" | grep -c . || true) samples (@RG SM), not one: Strelka2's ploidy VCF needs that one name." >&2
        exit 1
      fi
      PLOIDY="${OUTPUT_DIR}/${SAMPLE}.ploidy.vcf"
      # A record spans POS+1 to END (Strelka2's rule), so a BED row s-e is
      # POS s; a row from base 1 starts at POS 1 instead, leaving base 1
      # diploid (it is in the N run at the start of GRCh38's chrX and chrY).
      awk -F'\t' -v OFS='\t' -v sm="$SM" '
        FNR == NR { if ($1 == "chrX" || $1 == "chrY") len[$1] = $2; next }
        FNR == 1 {
          print "##fileformat=VCFv4.2"
          if ("chrX" in len) print "##contig=<ID=chrX,length=" len["chrX"] ">"
          if ("chrY" in len) print "##contig=<ID=chrY,length=" len["chrY"] ">"
          print "##ALT=<ID=CNV,Description=\"Copy number variable region\">"
          print "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the region\">"
          print "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number\">"
          print "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sm }
        { print $1, ($2 > 0 ? $2 : 1), ".", "N", "<CNV>", ".", "PASS", "END=" $3, "CN", 1 }' \
        "${REF}.fai" "$HAP" > "$PLOIDY"
      run_in --cpus 1 --memory 1g "$BCFTOOLS_IMAGE" \
        bcftools sort -Oz -o "$(cpath "${PLOIDY}.gz")" "$(cpath "$PLOIDY")"
      run_in --cpus 1 --memory 1g "$BCFTOOLS_IMAGE" \
        bcftools index -f -t "$(cpath "${PLOIDY}.gz")"
      PLOIDY_ARGS=(--ploidy "$(cpath "${PLOIDY}.gz")")
      echo "Ploidy VCF: ${PLOIDY}.gz ($(grep -vc '^#' "$PLOIDY") haploid stretches, sample ${SM})"
    fi
    rm -f "${OUTPUT_DIR}/xy.bed" "${OUTPUT_DIR}/xy.hap.bed"
  fi
  echo "[1/2] Configuring Strelka2 germline workflow..."
  run_in --cpus "$THREADS" --memory 16g \
    "$STRELKA_IMAGE" \
    configureStrelkaGermlineWorkflow.py \
      --bam "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
      --referenceFasta "${REF_FASTA_C}" \
      ${PLOIDY_ARGS[@]+"${PLOIDY_ARGS[@]}"} \
      --runDir "/genome/${SAMPLE}/vcf_strelka2"
  printf '%s\n' "$RUN_KEY" > "$RUN_FILE"
fi

# Step 2: Run the workflow
echo "[2/2] Running Strelka2 (this takes 1-2 hours for 30X WGS)..."
run_in --cpus "$THREADS" --memory 16g \
  "$STRELKA_IMAGE" \
  "/genome/${SAMPLE}/vcf_strelka2/runWorkflow.py" \
    -m local \
    -j "$THREADS"

if [ ! -f "${OUTPUT_DIR}/results/variants/variants.vcf.gz" ] || [ ! -f "${OUTPUT_DIR}/results/variants/variants.vcf.gz.tbi" ]; then
  echo "ERROR: Strelka2 finished without ${OUTPUT_DIR}/results/variants/variants.vcf.gz(.tbi)." >&2
  echo "  Delete ${OUTPUT_DIR}/ and run this step again." >&2
  exit 1
fi

VARIANT_COUNT=$(run_in \
  "$BCFTOOLS_IMAGE" \
  bcftools stats "/genome/${SAMPLE}/vcf_strelka2/results/variants/variants.vcf.gz" \
  | grep '^SN' | grep 'number of records' | awk '{print $NF}')
VARIANT_COUNT=${VARIANT_COUNT:-unknown}

echo "=== Strelka2 complete ==="
echo "Total variants called: ${VARIANT_COUNT}"
echo "Results: ${OUTPUT_DIR}/results/variants/variants.vcf.gz"
echo ""
echo "View PASS variants only:"
echo "  bcftools view -f PASS ${OUTPUT_DIR}/results/variants/variants.vcf.gz"
echo ""
echo "Count variant types:"
echo "  bcftools stats ${OUTPUT_DIR}/results/variants/variants.vcf.gz | grep 'number of records'"
