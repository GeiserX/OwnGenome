#!/usr/bin/env bash
# 37-y-haplogroup.sh — Y-chromosome haplogroup (paternal line) with Yleaf
# Usage: ./scripts/37-y-haplogroup.sh <sample_name>
#
# Opt-in. Yleaf (YLEAF_IMAGE) reads the BAM's pileup at its Y-chromosome
# markers and predicts the haplogroup they support on its own tree, YFull
# v14.01. The marker positions come from setup.sh --yleaf-data; the pileup
# is made in the same SAMTOOLS_IMAGE run as the chrY BAM (the current tree's
# pileup in a second one), and Yleaf's own command line reads it (-rp, reuse
# the pileup). Runs only when step 16 (indexcov) infers male from
# the reads: the sex comes from indexcov/indexcov-indexcov.ped, not from an
# argument. A female sample, or one with no Y to read, is skipped with one line.
#
# Output: y_haplogroup/<sample>_y_haplogroup.txt, Yleaf's prediction table
# (Hg, Hg_marker, Total_reads, Valid_markers, QC-score). Hg NA means too few
# markers had reads to call one (the step says so). Also
# y_haplogroup/<sample>_chrY.bam (+ .bai): the BAM's chrY reads only, the file
# YFull and the Y-DNA projects read (docs/37-y-haplogroup.md). Nothing is uploaded.
# The step also places the sample on the YFull YTree release versions.env
# pins (YTREE_VERSION) (setup.sh --yleaf-data installs it in <data>/ytree) with
# og y tree-place: y_haplogroup/<sample>_y_tree.txt (the terminal branch,
# with the tree it comes from) and <sample>_y_tree_path.tsv (derived,
# ancestral and no-call counts per level). Without that folder the step says
# so and writes Yleaf's result only.
# With step 03's VCF in place it also lists the chrY SNVs that no public SNP
# table names, the ones that would place this line below its branch
# (og y private-snps): y_haplogroup/<sample>_y_private_snps.tsv. Without
# the VCF, or without the tree folder, that part is skipped with one line.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image YLEAF_IMAGE SAMTOOLS_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
BAM="${S}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam"
PED="${S}/indexcov/indexcov-indexcov.ped"
OUTDIR="${S}/y_haplogroup"
OUT="${OUTDIR}/${SAMPLE}_y_haplogroup.txt"
CHRY="${OUTDIR}/${SAMPLE}_chrY.bam"
YDATA="${GENOME_DIR}/reference/yleaf-${YLEAF_DATA_VERSION}/data"
YTREE="${YDATA}/ytree"
VCF="${S}/vcf/${SAMPLE}.vcf.gz"

echo "=== Y haplogroup (Yleaf): ${SAMPLE} ==="
if [ ! -s "${YDATA}/hg38/new_positions.txt" ]; then
  echo "ERROR: Yleaf's marker tables are not installed (${YDATA}): run scripts/setup.sh --yleaf-data ${GENOME_DIR}" >&2
  exit 1
fi
for f in "$BAM" "${BAM}.bai" "$REF_FASTA"; do
  [ -f "$f" ] || { echo "ERROR: File not found: ${f}" >&2; exit 1; }
done
if [ ! -s "$PED" ]; then
  echo "ERROR: no sex check for ${SAMPLE} (${PED}): run scripts/16-indexcov.sh ${SAMPLE} first; this step reads the sex it infers." >&2
  exit 1
fi
# goleft's .ped: sex in PED coding (1 male, 2 female) under the header name 'sex'.
SEX=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) { h = $i; sub(/^#/, "", h); if (h == "sex") c = i }; next }
  c { s = $c } END { print (s == 1 ? "male" : s == 2 ? "female" : "unknown") }' "$PED")
if [ "$SEX" != male ]; then
  echo "Step 37 skipped: indexcov infers ${SEX} from the reads (${PED#"$S"/}); a Y haplogroup needs a Y chromosome."
  exit 0
fi

# Yleaf's --ref-fasta takes only these names (a compressed FASTA stops it
# with a traceback), as the Nextflow pipeline's --reference does.
case "$REF_FASTA" in
  *.fa | *.fasta | *.fna) ;;
  *) echo "ERROR: Yleaf needs a plain reference FASTA ending in .fa, .fasta or .fna (REF_FASTA is ${REF_FASTA})" >&2; exit 1 ;;
esac
mkdir -p "$OUTDIR"
rm -rf "${OUTDIR}/yleaf" "$OUT" "${OUTDIR}/${SAMPLE}_y_tree.txt" "${OUTDIR}/${SAMPLE}_y_tree_path.tsv" \
  "${OUTDIR}/${SAMPLE}_y_private_snps.tsv" "$CHRY" "${CHRY}.bai"
TREE=false
[ -s "${YTREE}/version.txt" ] && TREE=true
# og y yleaf-positions writes Yleaf's marker positions and SAMTOOLS_IMAGE
# makes the pileup Yleaf would make (its flags at its default quality 20), in
# the run that also writes the chrY BAM. Yleaf then reads that pileup where it
# keeps its own (-rp), offline: --ref-fasta names the reference so it
# downloads none (it never reads a BAM's sequence), --no-update-check skips
# GitHub. It stops on an empty pileup, so with no read at any marker
# og y yleaf-no-reads writes the table Yleaf gives for too few markers.
# The current tree's pileup is ours, with its flags written out
# (-A -Q20 -q20): fused as -AQ20q1, samtools reads "20q1" as the -Q value
# and keeps reads of mapping quality 0.
O=$(cpath "$OUTDIR")
echo "[1/3] Yleaf's Y marker positions..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 2g "$YLEAF_IMAGE" \
  /pgp-bin/og y yleaf-positions --data "$(cpath "$YDATA")" "${O}/positions.txt"
if $TREE; then
  run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 2g "$YLEAF_IMAGE" \
    /pgp-bin/og y tree-positions --tree-dir "$(cpath "$YTREE")" "${O}/tree_positions.txt"
fi
echo "[2/3] Pileup at those positions, and the chrY reads as their own BAM (samtools)..."
run_in --cpus 1 --memory 2g "$SAMTOOLS_IMAGE" bash -euo pipefail -c \
  'samtools idxstats "$1" > "$2/idxstats.txt" && samtools mpileup -l "$2/positions.txt" -ABQ20q1 "$1" > "$2/pileup.txt"
   samtools view -b -o "$3" "$1" chrY && samtools index "$3"' \
  _ "$(cpath "$BAM")" "$O" "$(cpath "$CHRY")"
if $TREE; then
  run_in --cpus 1 --memory 2g "$SAMTOOLS_IMAGE" bash -euo pipefail -c \
    'samtools mpileup -l "$2/tree_positions.txt" -A -Q20 -q20 "$1" > "$2/tree_pileup.txt"' _ "$(cpath "$BAM")" "$O"
fi
echo "[3/3] Yleaf..."
echo "Reference: ${REF_FASTA}"
# Yleaf names the sample, and its folder, after the BAM without its extension.
STEM=$(basename "$BAM"); STEM=${STEM%.*}
if [ -s "${OUTDIR}/pileup.txt" ]; then
  mkdir -p "${OUTDIR}/yleaf/${STEM}"
  mv "${OUTDIR}/pileup.txt" "${OUTDIR}/yleaf/${STEM}/temp_haplogroup_pileup.pu"
  run_in --cpus 1 --memory 4g "$YLEAF_IMAGE" \
    Yleaf -bam "$(cpath "$BAM")" -o "${O}/yleaf" -rg hg38 --ref-fasta "$REF_FASTA_C" -tree yfull \
      -q 20 -t 1 -rp -force --no-update-check
  rm -f "${OUTDIR}/yleaf/${STEM}/temp_haplogroup_pileup.pu"
else
  mkdir -p "${OUTDIR}/yleaf"
  run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 2g "$YLEAF_IMAGE" \
    /pgp-bin/og y yleaf-no-reads --bam "$(cpath "$BAM")" --idxstats "${O}/idxstats.txt" "${O}/yleaf/hg_prediction.hg"
fi
rm -f "${OUTDIR}/idxstats.txt" "${OUTDIR}/pileup.txt"
PRED="${OUTDIR}/yleaf/hg_prediction.hg"
if [ ! -s "$PRED" ] || [ "$(grep -c . "$PRED")" -lt 2 ]; then
  # Yleaf logs a failed sample as "Skipping" and goes on: its log says why.
  echo "ERROR: Yleaf wrote no prediction (${PRED}); its log is ${OUTDIR}/yleaf/run.log" >&2
  exit 1
fi
TREE_LINE="Current tree: not installed (${YTREE}); run scripts/setup.sh --yleaf-data ${GENOME_DIR} to add it"
if $TREE; then
  TREE_LINE=$(run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 4g "$YLEAF_IMAGE" \
    /pgp-bin/og y tree-place --tree-dir "$(cpath "$YTREE")" --pileup "${O}/tree_pileup.txt" \
      --sample "$SAMPLE" --out-prefix "${O}/${SAMPLE}")
  rm -f "${OUTDIR}/tree_positions.txt" "${OUTDIR}/tree_pileup.txt"
fi

# The SNVs past the deepest branch the public trees have: the sample's own
# chrY calls minus every position the installed SNP tables name
# (og y private-snps), inside the chrY callable mask. Needs step 03's VCF
# and the tree folder's tables and mask.
PRIVATE_LINE="Private SNPs: skipped (the current tree's SNP tables are not installed)"
if $TREE; then
  if [ ! -s "${YTREE}/callable_mask_chrY.fasta.gz" ]; then
    PRIVATE_LINE="Private SNPs: skipped (${YTREE} has no chrY callable mask: run scripts/setup.sh --yleaf-data again)"
  elif [ -s "$VCF" ]; then
    PRIVATE_LINE=$(run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 2g "$YLEAF_IMAGE" \
      /pgp-bin/og y private-snps --vcf "$(cpath "$VCF")" --tree-dir "$(cpath "$YTREE")" \
        --sample "$SAMPLE" --out-prefix "${O}/${SAMPLE}")
  else
    PRIVATE_LINE="Private SNPs: skipped (no VCF at ${VCF}: run scripts/03-deepvariant.sh ${SAMPLE} first)"
  fi
fi

cp "$PRED" "$OUT"
read -r HG MARKERS QC < <(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
  NR == 2 { print $c["Hg"], $c["Valid_markers"], $c["QC-score"] }' "$OUT")
echo ""
if [ "$HG" = NA ]; then
  echo "Y haplogroup: insufficient markers (Yleaf found ${MARKERS} Y markers with haplogroup information; none passes its quality threshold)"
else
  echo "Y haplogroup: ${HG} (${MARKERS} markers, QC-score ${QC})"
fi
echo "  (Yleaf ${YLEAF_DATA_VERSION}, on its own tree, YFull v14.01)"
echo "$TREE_LINE"
echo "Results: ${OUT}"
if $TREE; then echo "Current tree: ${OUTDIR}/${SAMPLE}_y_tree.txt, per level: ${OUTDIR}/${SAMPLE}_y_tree_path.tsv"; fi
echo "$PRIVATE_LINE"
if [ -s "${OUTDIR}/${SAMPLE}_y_private_snps.tsv" ]; then
  echo "Private SNPs: ${OUTDIR}/${SAMPLE}_y_private_snps.tsv"
fi
echo "Y reads for YFull or a Y-DNA project: ${CHRY} (nothing is uploaded; docs/37-y-haplogroup.md lists the routes)"
