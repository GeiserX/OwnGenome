#!/usr/bin/env bash
# report-corpus.sh: writes the summaries the report-render parity cases
# render, and the oracle that renders them with bin/render_report.py, in the
# current directory. scripts/ci/oracle-diff.sh runs it.
#
# The summaries:
#   <case>.json   one per summary-collect parity case: og summary collect
#                 run on that case's sample folder with its arguments, every
#                 section (no --only). A case og refuses has no summary.
#   m-<what>.json the summary of case "full" changed one way each with jq:
#                 the frame's sections in each state, new_since and the
#                 re-check, the manifest with and without a release, a
#                 sample name with characters HTML escapes and characters
#                 outside ASCII, counts missing or of another type.
#   a-<what>.json case c-full's summary changed one way each, for the
#                 ancestry chapter and the PRS card (see that block below).
# py/ holds the modules render_report.py imports, as they were at the
# oracle's commit (tests/parity/oracle.ref, report-render).
#
# oracle.sh renders with render_report.py, then keeps of its two reports
# the parts og report render writes so far: the frame and the sections of
# PORTED (each section's block in the text report, its frame and its line
# in the contents in the HTML report). A chapter left with no section goes,
# and so does the closing gene and rsID appendix, which og report index will
# write. Everything kept is compared byte for byte.
set -euo pipefail
export LC_ALL=C TZ=UTC
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "${here}/../../.." && pwd)
ORACLE_SHA=1ec630690e63796b4410d23ec65768200adb8dc0
og=${OG:-${root}/rust/target/release/og}

# --- the oracle's modules ------------------------------------------------------------
git -C "$root" cat-file -e "${ORACLE_SHA}^{commit}" 2> /dev/null \
  || git -C "$root" fetch -q --depth=1 origin "$ORACLE_SHA"
mkdir -p py
for m in collect_summary results_index genome_lookup report_charts live_results; do
  git -C "$root" show "${ORACLE_SHA}:bin/${m}.py" > "py/${m}.py"
done

# --- one summary per summary-collect case ---------------------------------------------
# Each generator writes its folders under folders/<generator>/, with the
# committed inputs it names beside them.
mkdir -p folders
awk -F'\t' '$1 == "summary-collect" { print $2 "\t" $3 "\t" $4 }' "${here}/../cases.tsv" > folders/cases.tsv
while IFS=$'\t' read -r case inputs args; do
  gen=${inputs%%,*}
  gen=${gen#gen:}
  dir="folders/${gen%.sh}"
  if [ ! -d "$dir" ]; then
    mkdir -p "$dir"
    IFS=, read -r -a items <<< "$inputs"
    for item in "${items[@]:1}"; do
      case "$item" in */oracle.sh) ;; *) cp -p "${here}/../inputs/${item}" "$dir"/ ;; esac
    done
    (cd "$dir" && bash "${here}/${gen}" > /dev/null)
  fi
  # The case's arguments, with every section and the summary written here.
  read -r -a argv <<< "$args"
  keep=()
  i=0
  while [ "$i" -lt "${#argv[@]}" ]; do
    case "${argv[$i]}" in
      --only | --out) i=$((i + 2)) ;;
      *) keep+=("${argv[$i]}"); i=$((i + 1)) ;;
    esac
  done
  (cd "$dir" && "$og" summary collect "${keep[@]}" --out "../../${case}.json" > /dev/null 2>&1) \
    || rm -f "${case}.json"
done < folders/cases.tsv
rm -rf folders

# --- the summary of case "full", changed one way each ---------------------------------
# Its callability section (130 kB of genes, not ported yet) is set to not run.
[ -s full.json ] || { echo "report-corpus.sh: og summary collect wrote no full.json" >&2; exit 1; }
base=base.tmp
jq -c '.sections.callability.state = "missing" | .sections.callability.data = {}' full.json > "$base"
mut() { jq -c "$2" "$base" > "m-$1.json"; }
for st in ok stale failed unreadable missing; do
  mut "qc-$st" "reduce (\"coverage\", \"sex_check\", \"sample_qc\", \"read_count\") as \$k (.;
    .sections[\$k].state = \"$st\" | .sections[\$k].note = \"note of \(\$k)\")"
  mut "variants-$st" ".sections.variants.state = \"$st\" | .sections.variants.note = \"why: <b>&'\\\"\""
done
mut qc-mixed '.sections.coverage.state = "stale" | .sections.coverage.note = "from 2020" | .sections.sample_qc.state = "failed"
  | .sections.read_count.state = "missing" | .sections.coverage.file_date = "2020-01-02 03:04:05"
  | .sections.sex_check.file_date = "2021-01-01 00:00:00"'
mut since-qc '.new_since = {"previous_utc": "2026-01-01T00:00:00Z", "sections": ["coverage"]}'
mut since-all '.new_since = {"previous_utc": null, "sections": [.sections | keys[]]}'
mut since-text '.new_since = {"previous_utc": "2026-01-01T00:00:00Z", "sections": "variants"}'
mut since-empty '.new_since = {}'
mut recheck-newer '.recheck = {"checked_utc": "2026-02-02T00:00:00Z", "items": [
  {"item": "clinvar_file_date", "label": "ClinVar", "class": "refresh-available"},
  {"item": "VEP_CACHE_RELEASE", "label": "", "class": "newer-release"},
  {"item": "PHARMCAT", "class": null}, {"item": "X", "label": "Y", "class": "no-newer-release"}]}'
mut recheck-none '.recheck = {"checked_utc": null, "items": []}'
mut release '.manifest = [["run", "git_commit", "0123456789abcdef"], ["run", "pipeline_version", "1.0.0"],
  ["image", "DEEPVARIANT_IMAGE", ("google/deepvariant" + ":1.9.0")], ["image", "PYTHON_IMAGE", ("python" + ":3.15")]]'
mut release-empty '.manifest = [["run", "git_commit", ""], ["run", "pipeline_version", ""]]'
mut no-manifest '.manifest = []'
mut started '.run.started_utc = "2026-03-03T03:03:03Z"'
mut sample-escape '.sample = "A&B <x> \"q\" '"'"'s"'
mut sample-unicode '.sample = "Muestra ñandú 𝔊 😀" | .generated_utc = "2026-10-01 ü"'
mut limits '.not_run = ["Step 06 <ClinVar> & co", "Step 07"] | .not_assessed = []'
mut sex '.run.declared_sex = "female" | .sections.sex_check.state = "ok" | .sections.sex_check.data.inferred_sex = "male"
  | .sections.sex_check.data.karyotype = "47,XXY" | .sections.sex_check.data.cn_chrX = 2
  | .sections.sample_qc.state = "ok" | .sections.sample_qc.data = {"inferred_sex": "male", "sex_check": "mismatch",
  "karyotype": "47,XXY", "x_depth_ratio": 0.98, "y_depth_ratio": 0.51, "contamination": "warn", "freemix": 0.031,
  "freemix_warn_above": 0.02, "verifybamid2_marker_check": "skipped", "same_person_as": "S2"}
  | .sections.read_count.state = "ok" | .sections.read_count.data = {"input_reads": 1000, "primary_records": 990,
  "missing_reads": 10, "read_count_check": "mismatch"}'
mut sex-agree '.run.declared_sex = "male" | .sections.sex_check.state = "ok" | .sections.sex_check.data.inferred_sex = "male"
  | .sections.sample_qc.state = "stale" | .sections.sample_qc.note = "from an earlier run"
  | .sections.sample_qc.data = {"sex_check": "not_checked", "sex_check_reason": "no <sites>", "contamination": "not_run"}
  | .sections.read_count.state = "ok" | .sections.read_count.data = {"input_reads": 5, "primary_records": 7,
  "missing_reads": -2, "read_count_check": "mismatch"}'
mut counts '.sections.variants.state = "ok" | .sections.variants.data = {"total": 4.5, "pass": true, "snps": "12", "indels": null}
  | .sections.coverage.state = "ok" | .sections.coverage.data.mean_depth = 31'
mut counts-float '.sections.coverage.state = "ok" | .sections.coverage.data.mean_depth = 30.05'
rm -f "$base"

# --- the ancestry chapter and the PRS card, from case c-full -----------------------
# a-<what>.json: c-full without its one PRS row that has no numbers (the
# Python's HTML report cannot render it), changed one way each: every state,
# a planted escape in the haplogroup links and names, a painting with its
# legend and the groups left out, the raw and the sex-specific scores, the Y
# placement's ancestral calls, new_since, and numbers where the summary
# writes text (a text-only case: the Python's HTML report cannot render it).
# A planted address that is not a web address (javascript:) is no case: og
# shows it as text where the Python made it a link; the cargo test
# ancestry_haplogroups.rs covers it. Nor is a legend colour that is not
# #rrggbb: og shows it grey where the Python put it in the CSS; the cargo
# test ancestry_fine.rs covers it, and a-escape plants a mixed-case one.
[ -s c-full.json ] || { echo "report-corpus.sh: og summary collect wrote no c-full.json" >&2; exit 1; }
abase=abase.tmp
jq -c '.sections.prs.data.scores |= map(select(.score != null))' c-full.json > "$abase"
amut() { jq -c "$2" "$abase" > "a-$1.json"; }
akeys='("ancestry_fine", "haplogroup", "y_haplogroup", "y_tree", "y_str", "archaic", "ancient", "prs")'
amut full '.'
for st in stale failed unreadable missing partial; do
  amut "$st" "reduce $akeys as \$k (.; .sections[\$k].state = \"$st\" | .sections[\$k].note = \"why: <b>&'\\\"\")"
done
amut escape '.sections.haplogroup.data.haplogroup = "H<b>&\"x'"'"'"
  | .sections.haplogroup.data.explore[0].url = "https://discover.familytreedna.com/mtdna/\"><script>alert(1)</script>"
  | .sections.haplogroup.data.explore[1].name = "<i>YFull</i> & co"
  | .sections.haplogroup.data.contamination_status = "YES" | .sections.haplogroup.data.contamination_level = "0.4<"
  | .sections.y_haplogroup.data.explore[1].url = "https://www.yfull.com/tree/\"x\"&y/"
  | .sections.y_tree.data.path[0].Branch = "<R>" | .sections.y_tree.data.path[0].Ancestral_SNPs = "M1&M2"
  | .sections.y_str.data.markers[0].marker = "DYS<393>"
  | .sections.ancient.data.nearest[0].group = "Spain_<BA>" | .sections.ancient.data.nearest[0].locality = "A & B"
  | .sections.ancestry_fine.data.continent[0].group = "West <Eurasia>"
  | .sections.ancestry_fine.data.legend[0].colour = "#AbCdEf" | .sections.ancestry_fine.data.legend[1].group = "<b>W</b>"
  | .sections.archaic.data.longest = ["chr1:1-2 <GENE>"]
  | .sections.prs.data.scores[0].condition = "C<A>D & \"x\"" | .sections.prs.data.scores[0].effect = "HR <1.5>"'
amut painting '.sections.ancestry_fine.data.painting = [{"code": "AFR", "group": "Africa", "fraction": 0.3},
    {"code": "WEA", "group": "West Eurasia", "fraction": 0.7}]
  | .sections.ancestry_fine.data.not_painted = []
  | .sections.ancestry_fine.data.left_out = [{"group": "Oceania (continental group)"}, {"group": "Oceania (regional group)"},
    {"group": "Papua (regional group)"}, {"group": "Odd (name) (x group)"}, {"group": "Plain"}]
  | .sections.ancestry_fine.data.min_references = 12
  | .sections.ancestry_fine.data.continent += [{"group": "Small", "fraction": 0.009, "low": 0.001, "high": 0.02}]
  | .sections.ancestry_fine.data.region += [{"group": "Wide", "fraction": 0.005, "low": 0.0, "high": 0.06},
    {"group": "Narrow", "fraction": 0.005, "low": 0.0, "high": 0.02}, {"group": "Unfit", "fraction": null, "low": null, "high": null}]'
amut no-painting '.sections.ancestry_fine.data.painting = [] | .sections.ancestry_fine.data.not_painted = []
  | .sections.ancestry_fine.data.left_out = [] | .sections.ancestry_fine.data.painting_png = ""'
amut prs-raw '.sections.prs.data.adjusted = false | del(.sections.prs.data.ancestry)
  | .sections.prs.data.scores |= map(.percentile = "" | .evaluated_in = null)'
amut prs-sex '.sections.prs.data.ancestry = {"population": "", "panel": "", "low_confidence": "True"}
  | .sections.prs.data.scores = [
    {"condition": "Prostate", "pgs_id": "PGS2", "score": "1", "matched": "5", "total": "6", "chrx": "0",
     "applies_to": "male", "sex_check": "not applicable to the declared sex", "percentile": "", "group": "",
     "evaluated_in": ["EUR"], "effect": ""},
    {"condition": "Breast", "pgs_id": "PGS3", "score": "2", "matched": "7", "total": "9", "chrx": "3 left out",
     "applies_to": "female", "sex_check": "no sex declared", "percentile": "40", "group": "AFR",
     "evaluated_in": [], "effect": ""},
    {"condition": "Height", "pgs_id": "PGS4", "score": "3", "matched": "1", "total": "2", "chrx": "",
     "applies_to": "", "sex_check": "", "percentile": "60", "group": "AFR", "evaluated_in": ["EUR", "AFR"],
     "effect": ""}]'
amut y-calls '.sections.y_tree.data.ancestral_in_terminal = "M1,M2" | .sections.y_tree.data.outvoted_ancestral = null
  | .sections.y_tree.data.path = [] | .sections.y_str.data.markers |= map(.count = "")
  | .sections.ancient.data.nearest = [] | .sections.ancient.data.plot_svg = ""
  | .sections.archaic.data.longest = []
  | .sections.haplogroup.data.explore = [] | del(.sections.haplogroup.data.contamination_status)
  | .sections.y_haplogroup.data.explore = []'
amut since '.new_since = {"previous_utc": "2026-01-01T00:00:00Z", "sections": ["y_str", "prs", "haplogroup"]}'
amut numbers '.sections.y_str.data.markers[0].count = 13 | .sections.y_str.data.markers[1].panel = true
  | .sections.ancient.data.nearest[0].rank = 1 | .sections.ancient.data.nearest[1].n = 5.5
  | .sections.y_tree.data.path[0].Level = 1 | .sections.y_tree.data.path[1].Derived = 2.0
  | .sections.archaic.data.neanderthal_bp = 51234567.0 | .sections.archaic.data.denisovan_bp = true
  | .sections.prs.data.scores[0].condition = 7 | .sections.prs.data.scores[1].score = 0.25
  | .sections.ancestry_fine.data.continent[1].fraction = 1'
rm -f "$abase"

# --- the oracle -----------------------------------------------------------------------
cat > oracle.sh << 'ORACLE'
#!/bin/sh
# The oracle of og report render (tests/parity/oracle.ref): render_report.py
# with the modules of py/, then the parts of its reports og writes so far.
#   sh ./oracle.sh PROG --summary F -o R.txt -o R.html
set -e
prog=$1
shift
d=$(mktemp -d)
cp py/*.py "$d"/
cp "$prog" "$d/render_report.py"
PYTHONPATH=$d python3 "$d/render_report.py" "$@"
rm -rf "$d"
# The sections og report render writes; the others are left out.
PORTED=" qc variants not-run not-assessed "
PORTED="${PORTED}ancestry-fine mt-haplogroup y-chromosome archaic ancient prs "
while [ $# -gt 0 ]; do
  case "$1" in
    -o)
      out=$2
      shift
      case "$out" in
        *.html) awk -v keep="$PORTED" -f html.awk "$out" > "$out.kept" ;;
        *.txt) awk -v keep="$PORTED" -f text.awk "$out" > "$out.kept" ;;
      esac
      mv "$out.kept" "$out" ;;
  esac
  shift
done
ORACLE

# html.awk: drops each section not in keep (its line in the contents and its
# frame) and every chart (og report charts will write them), then each
# chapter with no section left (its contents entry and its body), and the
# gene and rsID appendix.
cat > html.awk << 'AWK'
function wanted(line,   id) {
  if (!match(line, /id="[^"]*"/) && !match(line, /href="#[^"]*"/)) return 1
  id = substr(line, RSTART, RLENGTH)
  sub(/^(id="|href="#)/, "", id); sub(/"$/, "", id)
  return index(keep, " " id " ") > 0
}
{ lines[++n] = $0 }
END {
  # Pass 1: drop the sections not kept.
  m = 0
  for (i = 1; i <= n; i++) {
    l = lines[i]
    if (l ~ /^        <li><a href="#/ && !wanted(l)) continue
    if (l ~ /^    <figure class="chart">/) continue
    if (l ~ /^  <section class="sec[^"]*" id="/ && !wanted(l)) {
      while (lines[i] != "  </section>") i++
      continue
    }
    out[++m] = l
  }
  # Pass 2: drop the chapters with no section left, and the appendix.
  for (i = 1; i <= m; i++) {
    l = out[i]
    if (l ~ /^    <li><a href="#ch-terms">/) continue
    if (l ~ /^    <li><a href="#ch-/ && out[i + 1] == "      <ol>" && out[i + 2] == "      </ol>") { i += 3; continue }
    if (l ~ /^<section class="chapter" id="ch-/) {
      j = i + 4
      if (l ~ /id="ch-terms"/ || out[j] == "  </div>") {
        while (out[i] != "</section>") i++
        continue
      }
    }
    print l
  }
}
AWK

# text.awk: drops each section block not in keep, then each chapter left
# with no block. A block starts at its "## " title; the titles of the kept
# sections are listed here, as render_report.py writes them. A PRS row
# whose counts the summary has as null reads "not available", as og writes
# it; render_report.py writes "None" there (and stops its HTML report on
# that row).
cat > text.awk << 'AWK'
BEGIN {
  rule = sprintf("%80s", ""); gsub(/ /, "#", rule)
  title["## Quality control"] = "qc"; title["## Variant Calling (DeepVariant)"] = "variants"
  title["## Fine-Grained Ancestry (ADMIXTURE, RFMix)"] = "ancestry-fine"; title["## Mitochondrial Haplogroup"] = "mt-haplogroup"
  title["## Y-Chromosome Haplogroup (Yleaf)"] = "y-chromosome"; title["## Y Placement on the Current Public Tree"] = "y-chromosome"
  title["## Y-STR Markers (HipSTR)"] = "y-chromosome"; title["## Archaic Ancestry (IBDmix, step 43)"] = "archaic"
  title["## Ancient Populations (AADR PCA, step 44)"] = "ancient"; title["## Polygenic Risk Scores"] = "prs"
}
{ sub(/ None \(None\/None variants matched\) /, " not available (not available/not available variants matched) "); lines[++n] = $0 }
END {
  # The end part (steps not run onward) and the head (before the first
  # banner) stay; chapters are banner, title, banner, blank, then blocks.
  for (i = 1; i <= n; i++) if (lines[i] == "## Steps Not Run") { stop = i; break }
  for (i = 1; i < stop; i++) {
    if (lines[i] == rule && lines[i + 1] ~ /^#  [0-9]+\. /) {
      head = lines[i] "\n" lines[i + 1] "\n" lines[i + 2] "\n" lines[i + 3]
      i += 3; body = ""; on = 0
      while (i + 1 < stop && !(lines[i + 1] == rule && lines[i + 2] ~ /^#  [0-9]+\. /)) {
        i++
        if (lines[i] ~ /^## /) on = (lines[i] in title) && index(keep, " " title[lines[i]] " ") > 0
        if (on) body = body lines[i] "\n"
      }
      if (body != "") printf "%s\n%s", head, body
      continue
    }
    print lines[i]
  }
  for (i = stop; i <= n; i++) print lines[i]
}
AWK
