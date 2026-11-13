#!/usr/bin/env bash
# The Python that og bloodtype replaced, as it was at the commit
# tests/parity/oracle.ref names for bloodtype, for its oracle
# (inputs/bloodtype/oracle.sh): bloodtype_report.py imports copy_number.py,
# and every case reads that commit's marker table and bin/cn_loci.tsv.
set -euo pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd)
sha=$(awk -F'\t' '$1 == "bloodtype" {print $2}' "${root}/tests/parity/oracle.ref")
[ -n "$sha" ] || { echo "bloodtype-oracle.sh: no bloodtype row in oracle.ref" >&2; exit 1; }
git -C "$root" cat-file -e "${sha}^{commit}" 2> /dev/null || git -C "$root" fetch -q --depth=1 origin "$sha"
for f in bin/copy_number.py bin/cn_loci.tsv assets/bloodtype_markers.tsv; do
  git -C "$root" show "${sha}:${f}" > "${f##*/}"
done
