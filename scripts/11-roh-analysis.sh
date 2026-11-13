#!/usr/bin/env bash
# ROH Analysis — Runs of homozygosity (consanguinity screening)
# Detects long stretches of homozygous DNA that indicate shared ancestry.
# The summary keeps chr1 to chr22 only and leaves out the part of each segment in a centromere
# or heterochromatin band (og roh summary's MASKED_BANDS), where few markers map.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
OUTPUT="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}_roh.txt"

echo "=== ROH Analysis: ${SAMPLE} ==="

for f in "$VCF" "${VCF}.tbi"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

# Auto-detect chip data: if FORMAT/PL is absent, use -G30 (genotype-only mode)
HAS_PL=$(run_in \
  -v "${GENOME_DIR}/${SAMPLE}/vcf:/data" \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -h "/data/${SAMPLE}.vcf.gz" | grep -c '##FORMAT=<ID=PL' || true)

ROH_FLAGS=(--AF-dflt 0.4)
if [ "${HAS_PL}" -eq 0 ]; then
  echo "  No FORMAT/PL tag detected (chip data) — using -G30 genotype-only mode"
  ROH_FLAGS+=(-G30)
fi

# Only PASS calls (and records with no filter, as chip VCFs have): DeepVariant's
# RefCall and other filtered records are not genotypes to count on. The same
# input as the Nextflow ROH module. pipefail: a failed view fails the step.
# (The flags and the validated sample name are plain words, safe inline.)
run_in \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}/${SAMPLE}/vcf:/data" \
  "${BCFTOOLS_IMAGE}" \
  bash -euo pipefail -c "bcftools view -f PASS,. -Ou /data/${SAMPLE}.vcf.gz | bcftools roh ${ROH_FLAGS[*]} -o /data/${SAMPLE}_roh.txt -"

# Summary: autosomal segments with 5 Mb or more outside GRCh38's centromere
# and heterochromatin bands (the acen, gvar and stalk rows of UCSC's
# cytoBand.hg38.txt, merged; MASKED_BANDS in og's step 11 code, the one copy),
# the threshold of docs/11-roh-analysis.md and of the reports; length_bp and
# length_mb are of that part. The Nextflow module runs the same command.
SUMMARY="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}_roh_summary.txt"
run_in \
  --cpus 1 --memory 1g \
  -v "${GENOME_DIR}/${SAMPLE}/vcf:/data" \
  -v "${PGP_ROOT}/bin:/pgp-bin:ro" \
  "${BCFTOOLS_IMAGE}" \
  /pgp-bin/og roh summary --sample "$SAMPLE" --roh "/data/${SAMPLE}_roh.txt" > "${SUMMARY}.tmp"
mv -f "${SUMMARY}.tmp" "$SUMMARY"

echo "=== ROH complete ==="
echo "Results: ${OUTPUT}"
echo "Summary: ${SUMMARY}"
echo ""
echo "Autosomal ROH of 5 Mb or more outside centromere and heterochromatin bands (potential consanguinity signal):"
awk -F'\t' 'NR > 3 {printf "%s:%s-%s  %sMB\n", $1, $2, $3, $5}' "$SUMMARY"
echo ""
echo "NOTE: no population allele frequencies were used. bcftools roh takes INFO/AC and INFO/AN where a record"
echo "      has them (in a single-sample VCF, the sample's own counts) and 0.4 elsewhere (--AF-dflt), so read"
echo "      the total and the short segments as rough; the long segments are the signal."
