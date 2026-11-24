#!/usr/bin/env bash
# ROH Analysis — Runs of homozygosity (consanguinity screening)
# Detects long stretches of homozygous DNA that indicate shared ancestry.
# The summary keeps chr1 to chr22 only and leaves out the part of each segment in a centromere
# or heterochromatin band (og roh summary's MASKED_BANDS), where few markers map.
# With the allele-frequency panel setup.sh installs (roh_af/, ROH_AF_PANEL_URL in versions.env)
# bcftools roh reads each site's population frequency from it and uses the panel's sites only;
# without it every site gets a flat 0.4 (--AF-dflt). docs/11-roh-analysis.md says which and why.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
OUTPUT="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}_roh.txt"
PANEL="${GENOME_DIR}/roh_af/${ROH_AF_PANEL_URL##*/}"
ANCESTRY="${GENOME_DIR}/${SAMPLE}/ancestry/${SAMPLE}_ancestry.tsv"
GROUP_AF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}_roh_af.tsv.gz"

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

# The allele frequencies. With the panel: its global AF (column 4), or the AF
# of one continental group when step 26 gives that group a probability of
# 0.8 or more (EUR reads gnomAD's nfe). bcftools roh then uses the panel's
# sites only, with matching alleles (--AF-file restricts the markers by
# itself; it refuses -T). A group run gets its own 4-column file, with the
# group's AF of 0.01 to 0.99, the global panel's range. The Nextflow ROH
# module picks and cuts the same way. Without the panel: 0.4 at every site,
# whatever the record's INFO says.
AF_DESC="0.4 at every site (--AF-dflt 0.4)"
AF_LINE=""
if [ -s "$PANEL" ] && [ -s "${PANEL}.tbi" ]; then
  GROUP=$( [ ! -s "$ANCESTRY" ] || awk -F'\t' '$1 ~ /^probability_(AFR|AMR|EAS|EUR|SAS)$/ && $2 + 0 >= 0.8 {print substr($1, 13), $2}' "$ANCESTRY")
  if [ -z "$GROUP" ]; then
    ROH_FLAGS=(--AF-file "/genome/roh_af/${PANEL##*/}")
    AF_DESC="global (panel ${PANEL##*/}); no step 26 group at 0.8 or more"
  else
    read -r POP PROB <<< "$GROUP"
    case "$POP" in AFR) COL=5 ;; AMR) COL=6 ;; EAS) COL=7 ;; EUR) COL=8 ;; SAS) COL=9 ;; esac
    GNOMAD_GROUP=$(echo "afr amr eas nfe sas" | cut -d' ' -f$((COL - 4)))
    rm -f "$GROUP_AF" "${GROUP_AF}.tbi"
    run_in \
      --cpus 2 --memory 1g \
      -v "${GENOME_DIR}/${SAMPLE}/vcf:/data" \
      "${BCFTOOLS_IMAGE}" \
      bash -euo pipefail -c "bgzip -dc '/genome/roh_af/${PANEL##*/}' \
        | awk -F'\t' -v OFS='\t' -v c=${COL} '\$1 !~ /^#/ && \$c != \".\" && \$c + 0 >= 0.01 && \$c + 0 <= 0.99 {print \$1, \$2, \$3, \$c}' \
        | bgzip > '/data/${SAMPLE}_roh_af.tsv.gz' && tabix -s1 -b2 -e2 '/data/${SAMPLE}_roh_af.tsv.gz'"
    ROH_FLAGS=(--AF-file "/data/${SAMPLE}_roh_af.tsv.gz")
    AF_DESC="${GNOMAD_GROUP} (panel ${PANEL##*/}); step 26 gives ${POP} ${PROB}"
  fi
  AF_LINE="# AF used: ${AF_DESC}"
else
  ROH_FLAGS=(--AF-dflt 0.4)
  echo "  No allele-frequency panel at ${PANEL} (scripts/setup.sh installs it): every site gets AF 0.4"
fi
echo "  Allele frequencies: ${AF_DESC}"
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
rm -f "$GROUP_AF" "${GROUP_AF}.tbi"
# The table names the frequencies it was found with, when they came from the
# panel (bcftools roh's own header gives its command line).
if [ -n "$AF_LINE" ]; then
  { printf '%s\n' "$AF_LINE"; cat "$OUTPUT"; } > "${OUTPUT}.tmp"
  mv -f "${OUTPUT}.tmp" "$OUTPUT"
fi

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
if [ -n "$AF_LINE" ]; then
  echo "NOTE: allele frequencies: ${AF_DESC}. Only the panel's"
  echo "      sites were used. The long segments are the signal; segment edges and the short segments are estimates."
else
  echo "NOTE: no population allele frequencies were used. With --AF-dflt 0.4 bcftools roh gives every site an"
  echo "      alternate allele frequency of 0.4, also where the record carries INFO/AC and INFO/AN. That inflates"
  echo "      the short segments, so read the total and the short segments as rough; the long segments are the signal."
fi
