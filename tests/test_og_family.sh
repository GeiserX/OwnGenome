#!/usr/bin/env bash
# og-leg: OG_REQUIRED
# `og family` through the og binary (OG, else bin/og), the way run-all.sh and
# the FAMILY workflow call it. The cases are cargo tests
# (rust/crates/og-steps/tests/family_report.rs) that call the library; this
# test checks that the binary still sends `og family` to that code.
# Checks:
#   1. plan --which prints the family of a sample of
#      tests/fixtures/families/two_families.tsv, and prints nothing (exit 0)
#      for a sample in no family;
#   2. plan --out writes the plan table, with its five columns;
#   3. prepare on a trio writes somalier's pedigree and says what it checks;
#   4. an unknown --family stops with exit 2 and names the file's families.
# Red control: check 1 on a sample of the other family gives another name.
# Without og the test is skipped with a note; with OG_REQUIRED set, the skip
# is a failure (the og-musl job of guard.yml sets it).
set -uo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
fails=0
ok() { echo "ok: $1"; }
bad() { echo "FAIL: $1"; fails=$((fails + 1)); }

OG=${OG:-${REPO}/bin/og}
if ! "$OG" version > /dev/null 2>&1; then
  if [ -n "${OG_REQUIRED:-}" ]; then
    echo "FAIL: og is not built (scripts/setup.sh --og), and OG_REQUIRED is set"
    exit 1
  fi
  echo "skip: og is not built for this machine (scripts/setup.sh --og)"
  exit 0
fi
two="${REPO}/tests/fixtures/families/two_families.tsv"

# 1. --which, and its control
got=$("$OG" family plan --families "$two" --sample partner_y --which 2> "${T}/err")
if [ "$got" = "family_b" ]; then
  ok "plan --which: partner_y is in family_b"
else
  bad "plan --which: partner_y gave '${got}': $(cat "${T}/err")"
fi
got=$("$OG" family plan --families "$two" --sample child_a --which 2> /dev/null)
if [ -n "$got" ] && [ "$got" != "family_b" ]; then
  ok "control: child_a gives another family (${got}), so the check above can fail"
else
  bad "control: child_a gave '${got}'; the check above cannot fail"
fi
rc=0
got=$("$OG" family plan --families "$two" --sample nobody --which 2> /dev/null) || rc=$?
if [ "$rc" -eq 0 ] && [ -z "$got" ]; then
  ok "plan --which: a sample in no family prints nothing and exits 0"
else
  bad "plan --which: nobody printed '${got}' (exit ${rc})"
fi

# 2. the plan table
if "$OG" family plan --families "$two" --family family_a --out "${T}/plan.tsv" > "${T}/log" 2>&1 \
  && [ "$(head -n 1 "${T}/plan.tsv")" = "$(printf 'family\tanalysis\tmembers\tdecision\treason')" ]; then
  ok "plan --out writes the plan table"
else
  bad "plan --out: $(cat "${T}/log"; head -n 1 "${T}/plan.tsv" 2> /dev/null)"
fi

# 3. prepare on a trio, with a grandparent named but not sequenced
printf 'sample\tsex\tfather\tmother\tpartner\nkid\tmale\tdad\tmum\t.\ndad\tmale\tgrandpa\t.\tmum\nmum\tfemale\t.\t.\tdad\n' > "${T}/ped.tsv"
args=(family prepare --pedigree "${T}/ped.tsv" --out-dir "${T}/out")
for s in kid dad mum; do
  printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\nchr1\t1\t.\tA\t<*>\t0\t.\tEND=10\tGT:GQ\t0/0:50\n' "$s" | gzip > "${T}/${s}.g.vcf.gz"
  args+=(--gvcf "${s}=${T}/${s}.g.vcf.gz")
done
want=$(printf 'family\tdad\t-9\t-9\t1\t-9\nfamily\tkid\tdad\tmum\t1\t-9\nfamily\tmum\t-9\t-9\t2\t-9')
if "$OG" "${args[@]}" > "${T}/log" 2>&1 && [ "$(sort "${T}/out/family.ped")" = "$want" ] \
  && grep -qF 'kid and mum: declared parent-child' "${T}/log"; then
  ok "prepare writes somalier's pedigree and names the pairs it checks"
else
  bad "prepare: $(cat "${T}/log" "${T}/out/family.ped" 2> /dev/null)"
fi

# 4. an unknown family
rc=0
"$OG" family plan --families "$two" --family family_c > "${T}/log" 2>&1 || rc=$?
if [ "$rc" -eq 2 ] && grep -qF "no family 'family_c'; its families are: family_a, family_b" "${T}/log"; then
  ok "plan: an unknown --family stops with exit 2 and lists the families"
else
  bad "plan --family family_c: exit ${rc}: $(cat "${T}/log")"
fi

echo "${fails} failure(s)"
[ "$fails" -eq 0 ]
