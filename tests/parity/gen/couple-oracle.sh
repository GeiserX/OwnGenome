#!/usr/bin/env bash
# The oracle of og couple screen and og couple genes: bin/couple_report.py
# imports bin/collect_summary.py, bin/family_report.py, bin/paraphase_report.py
# and bin/cnvpytor_to_vcf.py, and oracle-diff.sh materialises one file per
# tool. This stages the four modules, the carrier list and the PAR file, at the
# commit tests/parity/oracle.ref names for couple-screen, and oracle.sh, which
# puts the program beside the modules under its own name and runs it.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
sha=$(awk -F'\t' '$1 == "couple-screen" { print $2 }' "${ROOT}/tests/parity/oracle.ref")
[ -n "$sha" ] || { echo "couple-oracle.sh: no couple-screen row in oracle.ref" >&2; exit 1; }
git -C "$ROOT" cat-file -e "${sha}^{commit}" 2> /dev/null || git -C "$ROOT" fetch -q --depth=1 origin "$sha"
for m in collect_summary family_report paraphase_report cnvpytor_to_vcf; do
  git -C "$ROOT" show "${sha}:bin/${m}.py" > "${m}.py"
done
git -C "$ROOT" show "${sha}:assets/carrier_genes.tsv" > carrier_genes.tsv
git -C "$ROOT" show "${sha}:assets/par_grch38.bed" > par.bed
cat > oracle.sh <<'SH'
# sh ./oracle.sh PROG ARGS: PROG is couple_report.py; the modules are staged here.
prog=$1
shift
d=$(mktemp -d)
cp ./collect_summary.py ./family_report.py ./paraphase_report.py ./cnvpytor_to_vcf.py "$d"/
cp "$prog" "${d}/couple_report.py"
rc=0
PYTHONDONTWRITEBYTECODE=1 python3 "${d}/couple_report.py" "$@" || rc=$?
rm -rf "$d"
exit "$rc"
SH
