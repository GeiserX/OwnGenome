#!/usr/bin/env bash
# synthetic-fasta.sh: writes synthetic.fa and synthetic.bed in the current
# directory, the same every time. synthetic.fa holds two 125 Mb sequences in
# 60-base lines, about 254 MB; synthetic.bed holds 300 rows on chr1, none
# overlapping, so the masked file and the summary line do not depend on how
# overlaps are counted. scripts/ci/oracle-diff.sh runs it for the case
# setup-mask-fasta synthetic-254mb.
set -euo pipefail
export LC_ALL=C
LEN=${SYNTHETIC_LEN:-125000000}
for name in chr1 chr2; do
  echo ">${name} synthetic"
  # yes and tr stop on SIGPIPE when head has enough: not a failure.
  (set +o pipefail; yes ACGTTGCAacgtNNGGCCTTAAGTCAacgtacgtGGGGCCCCAAAATTTTNacgtACGTAC 2> /dev/null | tr -d '\n' 2> /dev/null | head -c "$LEN" | fold -w 60)
  echo
done > synthetic.fa
awk 'BEGIN { for (i = 0; i < 300; i++) { s = i * 400000 + 1234 + i; printf "chr1\t%d\t%d\tr%d\n", s, s + 5000 + i, i } }' > synthetic.bed
