#!/usr/bin/env bash
# The oracles of og sv add-callers and og sv postfilter: bin/cnvpytor_to_vcf.py
# imports bin/family_report.py (the PAR intervals) and bin/annotsv_postfilter.py
# imports both, and oracle-diff.sh materialises one file per tool. This stages
# the two modules, the rule table and the PAR file, at the commit
# tests/parity/oracle.ref names for sv-add-callers, and oracle.sh, which puts
# the program beside the modules under its own name and runs it.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
sha=$(awk -F'\t' '$1 == "sv-add-callers" { print $2 }' "${ROOT}/tests/parity/oracle.ref")
[ -n "$sha" ] || { echo "sv-oracle.sh: no sv-add-callers row in oracle.ref" >&2; exit 1; }
git -C "$ROOT" cat-file -e "${sha}^{commit}" 2> /dev/null || git -C "$ROOT" fetch -q --depth=1 origin "$sha"
git -C "$ROOT" show "${sha}:bin/family_report.py" > family_report.py
git -C "$ROOT" show "${sha}:bin/cnvpytor_to_vcf.py" > cnvpytor_to_vcf.py
git -C "$ROOT" show "${sha}:assets/sv_rule_regions.tsv" > rules.tsv
git -C "$ROOT" show "${sha}:assets/par_grch38.bed" > par.bed
cat > oracle.sh <<'SH'
# sh ./oracle.sh PROG ARGS: PROG is cnvpytor_to_vcf.py or annotsv_postfilter.py
# (the one that imports cnvpytor_to_vcf); the modules are staged here.
prog=$1
shift
name=cnvpytor_to_vcf.py
grep -q '^import cnvpytor_to_vcf' "$prog" && name=annotsv_postfilter.py
d=$(mktemp -d)
cp ./family_report.py ./cnvpytor_to_vcf.py "$d"/
cp "$prog" "${d}/${name}"
rc=0
PYTHONDONTWRITEBYTECODE=1 python3 "${d}/${name}" "$@" || rc=$?
rm -rf "$d"
exit "$rc"
SH
