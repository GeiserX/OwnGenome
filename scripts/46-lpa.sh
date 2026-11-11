#!/usr/bin/env bash
# 46-lpa.sh — Lp(a): the LPA KIV-2 repeat count estimated from read depth
# Usage: ./scripts/46-lpa.sh <sample_name>
#
# Opt-in. The number of kringle IV type 2 (KIV-2) repeats in LPA sets the size
# of apo(a); a low total goes with a higher Lp(a) level, a mostly inherited
# cardiovascular risk factor. The step counts the repeats from the depth over
# the KIV-2 exons against the depth over single-copy LPA exons
# (og cn call on the LPA_KIV2 row of bin/cn_loci.tsv; the same code as the
# LPA Nextflow processes), with a 95% interval. Beside the count it reads the KIV-2 4925G>A splice variant
# (Coassin et al. 2017), which lowers the Lp(a) made from the chromosome that
# carries it: samtools mpileup at the site's seven GRCh38 copies (the six
# KIV-2 units and KIV-1), read as a fraction of copies in the same og cn
# call, into its own table.
# The two common risk variants, rs10455872 and rs3798220, are rows
# of step 42's traits panel; the report's Lp(a) section reads both. A research
# estimate, not a diagnostic test: the actionable test is a blood Lp(a)
# measurement.
#
# Input:  aligned/<sample>_sorted.bam with its .bai (ALIGN_DIR as for the other
#         BAM steps), from an aligner that reports reads with several equally
#         good places (minimap2, BWA-MEM): with no read of MAPQ 0 over KIV-2
#         the count is not assessed.
# Output: lpa/<sample>_lpa_kiv2.tsv          metric<TAB>value: status, the count,
#                                            its 95% interval, the depths, message
#         lpa/<sample>_lpa_kiv2_depth.tsv    samtools bedcov -j over the regions (all MAPQ)
#         lpa/<sample>_lpa_kiv2_depth_q1.tsv the same, reads of MAPQ 1 or more
#         lpa/<sample>_lpa_kiv2_4925.tsv     metric<TAB>value: the 4925G>A reading
#                                            (status, reads, fraction, copies, message);
#                                            'not assessed' when the check fails
#         lpa/<sample>_lpa_kiv2_4925_pileup.tsv samtools mpileup at the seven
#                                            copies of the 4925G>A site (every MAPQ)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image PYTHON_IMAGE SAMTOOLS_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
BAM="${S}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam"
OUTDIR="${S}/lpa"
TABLE="${OUTDIR}/${SAMPLE}_lpa_kiv2.tsv"
DEPTH="${OUTDIR}/${SAMPLE}_lpa_kiv2_depth.tsv"
DEPTH_Q1="${OUTDIR}/${SAMPLE}_lpa_kiv2_depth_q1.tsv"
PILEUP="${OUTDIR}/${SAMPLE}_lpa_kiv2_4925_pileup.tsv"
NULL="${OUTDIR}/${SAMPLE}_lpa_kiv2_4925.tsv"
BED="${OUTDIR}/regions.bed"
SITES="${OUTDIR}/sites.txt"

echo "=== Lp(a): LPA KIV-2 repeat count from read depth: ${SAMPLE} ==="
[ -f "$BAM" ] || { echo "ERROR: File not found: ${BAM}" >&2; exit 1; }
[ -f "${BAM}.bai" ] || { echo "ERROR: File not found: ${BAM}.bai" >&2; exit 1; }
echo "BAM: ${BAM}"
mkdir -p "$OUTDIR"
rm -f "$TABLE" "$DEPTH" "$DEPTH_Q1" "$PILEUP" "$NULL" "$BED" "$SITES"

echo "[1/3] The regions (og cn bed --locus LPA_KIV2) and the 4925G>A sites..."
# shellcheck disable=SC2016  # $1 and $2 belong to the inner shell
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 1g "$PYTHON_IMAGE" \
  sh -c '/pgp-bin/og cn bed --locus LPA_KIV2 > "$1" &&
    /pgp-bin/og cn marker-regions LPA_KIV2_4925 > "$2"' _ "$(cpath "$BED")" "$(cpath "$SITES")"

echo "[2/3] Read depth over the KIV-2 and single-copy LPA exons (samtools bedcov -j, all MAPQ and MAPQ >= 1)"
echo "      and the reads at the 4925G>A sites (samtools mpileup, every MAPQ)..."
# bedcov stops on a contig the BAM has no @SQ line for: without chr6 the
# tables stay empty and the count is not assessed. mpileup runs once per site
# (-r), every MAPQ, duplicates left out, bases of quality 13 or more. A
# failed mpileup leaves the count to run: the pileup then holds one line that
# is not mpileup output, so 4925G>A reads 'not assessed'. The LPA_DEPTH
# process runs the same commands.
# shellcheck disable=SC2016  # $1 to $7 belong to the inner bash, the program to awk
run_in --cpus 1 --memory 2g "$SAMTOOLS_IMAGE" bash -euo pipefail -c '
    if samtools view -H "$2" | awk -F "\t" "$5"; then
      samtools bedcov -j "$1" "$2" > "$3"
      samtools bedcov -j -Q 1 "$1" "$2" > "$4"
      ( while read -r r; do samtools mpileup -A -B -q 0 -Q 13 -d 100000 -r "$r" "$2" || exit 1; done < "$7" ) > "$6" \
        || echo "# samtools mpileup failed: 4925G>A not assessed" > "$6"
    else : > "$3"; : > "$4"; : > "$6"; fi' \
  _ "$(cpath "$BED")" "$(cpath "$BAM")" "$(cpath "$DEPTH")" "$(cpath "$DEPTH_Q1")" \
  '$1 == "@SQ" {for (i = 2; i <= NF; i++) if ($i == "SN:chr6") f = 1} END {exit !f}' \
  "$(cpath "$PILEUP")" "$(cpath "$SITES")"
rm -f "$BED" "$SITES"

echo "[3/3] The count, its interval and KIV-2 4925G>A beside it (og cn call --locus LPA_KIV2 --marker LPA_KIV2_4925)..."
# One run reads the count and the marker, so the marker is judged on the
# count's own MAPQ 0 share. A failure of the marker never costs the count (as
# step 35's SMN1 marker): when the count table is there, the 4925G>A table
# then says not assessed, and the step goes on. The LPA process writes the
# same table.
rc=0
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 1g "$PYTHON_IMAGE" \
  /pgp-bin/og cn call --locus LPA_KIV2 --bedcov "$(cpath "$DEPTH")" \
    --bedcov-mapq1 "$(cpath "$DEPTH_Q1")" --table "$(cpath "$TABLE")" \
    || rc=$?
if [ "$rc" -ne 0 ]; then
  # No count table: the count itself failed (exit 3 for a broken control).
  [ -s "$TABLE" ] || exit "$rc"
  echo "WARNING: the KIV-2 4925G>A check failed; ${NULL##*/} says not assessed." >&2
  printf 'metric\tvalue\nmarker\tLPA KIV-2 4925G>A\nrsid\trs1484919544 (first KIV-2 copy)\nstatus\tnot assessed\nalt_reads\t.\ndepth\t.\nfraction\t.\ncopies\t.\nkiv2_copies\t.\nmessage\t%s\n' \
    'LPA KIV-2 4925G>A: not assessed, the check failed; see the log of step 46.' > "$NULL"
fi

[ -s "$TABLE" ] || { echo "ERROR: no table written (${TABLE})" >&2; exit 1; }
echo ""
echo "Table: ${TABLE}"
echo "KIV-2 4925G>A: ${NULL}"
