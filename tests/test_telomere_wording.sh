#!/usr/bin/env bash
# Step 10 measures a relative telomere content (intratelomeric reads per
# million reads with 48-52% GC). It is not an age, so no text that describes
# the step or its result may call it "biological age" (or ageing, aging).
#
# The files checked are the ones that describe step 10 or print its result.
# A listed file that a later change removes is skipped with a note; the step
# script and its Nextflow module must exist.
#
#   1. Control: the check run on a planted copy of the step script with the
#      phrase fails, so the check can fail.
#   2. The real files: none has the phrase.
set -uo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
FILES=(
  scripts/10-telomere-hunter.sh
  modules/local/telomere_hunter/main.nf
  docs/10-telomere-analysis.md
  docs/results/traits-and-somatic.md
  rust/crates/og-render/src/notes.rs
  rust/crates/og-render/src/sections/telomere.rs
  bin/render_report.py
)
REQUIRED="scripts/10-telomere-hunter.sh modules/local/telomere_hunter/main.nf"
PHRASE='biological[[:space:]-]+(age|ageing|aging)'
fails=0
ok()  { echo "PASS $*"; }
bad() { echo "FAIL $*"; fails=$((fails + 1)); }

# check ROOT: print each line of the listed files under ROOT that has the
# phrase; exit 1 when one does or a required file is missing.
check() {
  local root=$1 f rc=0
  for f in "${FILES[@]}"; do
    if [ ! -f "${root}/${f}" ]; then
      case " $REQUIRED " in
        *" $f "*) echo "missing: ${f}"; rc=1 ;;
        *) echo "note: ${f} is gone, skipped" ;;
      esac
      continue
    fi
    if grep -n -i -E "$PHRASE" "${root}/${f}" | sed "s|^|${f}:|"; then rc=1; fi
  done
  return "$rc"
}

T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
for f in "${FILES[@]}"; do
  [ -f "${REPO}/${f}" ] || continue
  mkdir -p "${T}/$(dirname "$f")"
  cp "${REPO}/${f}" "${T}/${f}"
done
printf '# Higher values = longer telomeres. Provides a Biological Age baseline.\n' >> "${T}/scripts/10-telomere-hunter.sh"
if out=$(check "$T"); then
  bad "control: a planted 'Biological Age' line in the step script passed the check"
elif grep -q '^scripts/10-telomere-hunter.sh:[0-9]*:.*Biological Age' <<< "$out"; then
  ok "control: the planted line fails the check and is named"
else
  bad "control: the check failed without naming the planted line: ${out}"
fi

if out=$(check "$REPO"); then
  ok "no file that describes step 10 calls its result a biological age"
else
  bad "step 10's wording:"
  printf '%s\n' "$out"
fi

[ "$fails" -eq 0 ] || { echo "${fails} check(s) failed"; exit 1; }
echo "all checks passed"
