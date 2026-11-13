#!/usr/bin/env bash
# The Python that og cn replaced, as it was at the commit tests/parity/oracle.ref
# names for cn, for the cn oracle (inputs/cn/oracle.sh): copy_number.py imports
# hba_copy_number.py and lpa_kiv2.py, the markers lived in traits_report.py,
# and every case reads that commit's bin/cn_loci.tsv.
set -euo pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd)
sha=$(awk -F'\t' '$1 == "cn" {print $2}' "${root}/tests/parity/oracle.ref")
[ -n "$sha" ] || { echo "cn-oracle.sh: no cn row in oracle.ref" >&2; exit 1; }
git -C "$root" cat-file -e "${sha}^{commit}" 2> /dev/null || git -C "$root" fetch -q --depth=1 origin "$sha"
for f in hba_copy_number.py lpa_kiv2.py traits_report.py cn_loci.tsv; do
  git -C "$root" show "${sha}:bin/${f}" > "$f"
done
