#!/usr/bin/env bash
# Check the notices in setup's usage output and in the docs for each resource.
# A call without a data directory prints usage before touching Docker or the network.
set -euo pipefail

repo=$(cd "$(dirname "$0")/.." && pwd)
rc=0
output=$(env -u GENOME_DIR bash "$repo/scripts/setup.sh" 2>&1) || rc=$?
if [ "$rc" -ne 1 ] || ! grep -q '^Usage:' <<< "$output"; then
  printf 'Expected setup usage (exit 1), got exit %s:\n%s\n' "$rc" "$output" >&2
  exit 1
fi

while IFS='|' read -r notice pages; do
  if ! grep -Fxq -- "  $notice" <<< "$output"; then
    printf 'Missing setup notice: %s\n' "$notice" >&2
    exit 1
  fi
  read -r -a docs <<< "$pages"
  for page in licences.md "${docs[@]}"; do
    if ! grep -Fq -- "$notice" "$repo/docs/$page"; then
      printf 'Missing notice in docs/%s: %s\n' "$page" "$notice" >&2
      exit 1
    fi
  done
done <<'NOTICES'
AnnotSV's annotation bundle includes third-party data with their own terms.|00-reference-setup.md 05-annotsv.md
PCGR/CPSR reference data includes third-party databases with their own terms.|00-reference-setup.md 17-cpsr.md
Exomiser's hg38 and phenotype bundles include third-party data with their own terms.|48-exomiser.md
The REVEL authors make the scores freely available for non-commercial use.|00-reference-setup.md 30-vcfanno.md hardware-requirements.md
UCSC allows free non-commercial use of its hg19ToHg38 chain by independent researchers and nonprofit organisations. Other uses require a commercial licence.|43-archaic.md 44-ancient.md chip-data-guide.md
Illumina's DRAGEN ORA Decompression Software requires acceptance of its end-user licence.|01-ora-to-fastq.md
The PacBio Software License Agreement permits Paraphase use only with data from PacBio instruments.|35-paralogs.md
NOTICES

echo "Licence notices match setup output and documentation."
