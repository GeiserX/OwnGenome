#!/usr/bin/env bash
# Step 20's shifted-reference pass on the HG002 fixture, from what base case 36
# wrote: the control region was called again on chrM rotated by 8000 bp, and
# those calls were mapped back into the filtered calls. A pass that failed
# would leave step 20 green, with a WARNING and the linear calls alone, so the
# traces it leaves are checked here.
. "$(dirname "$0")/lib.sh"

M="${SAMPLE}/mito"
VCF="${M}/${SAMPLE}_chrM_filtered.vcf.gz"
check "step 20's filtered calls are there (case 36)" vcf_ok "$VCF"
check "they carry the shifted-pass header line" \
  has '^##og_shifted_pass=chrM:16025-16569 and chrM:1-575 from Mutect2 on chrM rotated by 8000 bp' "$(bcf view -h "$VCF" 2>/dev/null)"
check "the rotated pass's calls are kept" vcf_ok "${M}/shifted/shifted_mutect2.vcf.gz"
check "the merged calls are kept" vcf_ok "${M}/shifted/merged_mutect2.vcf.gz"
check "the rotated calls are on the rotated chrM's control region only (8025-9144)" \
  test "$(bcf query -f '%POS\n' "${M}/shifted/shifted_mutect2.vcf.gz" 2>/dev/null | awk '$1 < 8025 || $1 > 9144' | wc -l)" -eq 0
check_ge "control-region calls in the filtered calls" "$(vcf_count -r chrM:1-575,chrM:16025-16569 "$VCF")" 1
LEFT=$(find "${GENOME_DIR}/${M}/shifted" -name '*.bam' -o -name '*.sam' -o -name '*.img' 2>/dev/null | wc -l | tr -d ' ')
check_eq "no alignment is left in mito/shifted" "$LEFT" 0

finish
