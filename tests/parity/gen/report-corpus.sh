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
# PGx, HLA and mitochondrial variants (port-report-pgx): each section in each
# state, the rows each writer prints, and texts HTML escapes.
pgx='.sections.pharmcat = {"title": "PharmCAT", "state": "ok", "source": "pharmcat/S.report.json", "file_date": null,
  "note": null, "data": {"version": "3.4.0", "html_report": true}}
  | .sections.cpic = {"title": "CPIC <recs>", "state": "ok", "source": "cpic/S_phenotypes.tsv", "file_date": null,
  "note": null, "data": {"genes": [{"gene": "CYP2D6", "diplotype": "*1/*4", "phenotype": "Intermediate <Metabolizer>",
  "status": "non-normal"}, {"gene": "CYP2C19", "diplotype": "*1/*17 & *2", "phenotype": "Rapid", "status": "ambiguous"},
  {"gene": "DPYD", "diplotype": null, "phenotype": null, "status": null}], "non_normal": 1, "ambiguous": 1,
  "not_called": 2, "unclassified": 1, "parse_failed": false, "warnings": ["two <results> & more"]}}
  | .sections.pypgx = {"title": "pypgx", "state": "ok", "source": "pypgx/S.tsv", "file_date": null, "note": null,
  "data": {"genes_total": 6, "genes_called": 3, "cyp2d6": "*1/*2", "cnv_call": {"CYP2D6": "2", "UGT2B17": "0 <del>"},
  "not_callable": ["GSTT1", "CYP2E1"], "low_gq": [{"gene": "UGT1A1", "diplotype": "*1/*28", "weakest": ""}],
  "depth_check": {"status": null, "gene_mapq1_share": "0.99", "body_mapq20_share": null},
  "comparison": {"conflicts": 2, "one_tool_only": 0}}}
  | .sections.cyrius = {"title": "Cyrius", "state": "ok", "source": "cyrius/S.tsv", "file_date": null, "note": null,
  "data": {"genotype": "*1/*2"}}
  | .cyp2d6 = {"calls": {"PharmCAT": "*1/*4", "pypgx": "*1/*2", "Cyrius": null}, "agree": false,
  "consensus": {"result": null, "passed_to_pharmcat": false, "reason": null}, "probable": {"diplotype": "*1/*2",
  "phenotype": "Normal Metabolizer", "label": "probable CYP2D6 call (ONE caller)"}}'
mut pgx-ok "$pgx"
mut pgx-stale "$pgx"' | .sections.cpic.state = "stale" | .sections.cpic.note = "from <2020>"
  | .sections.pypgx.state = "failed" | .sections.pypgx.note = "exit 1" | .sections.pypgx.data = {}
  | .sections.cyrius.state = "stale" | .sections.cyrius.note = "old" | .sections.pharmcat.data = {}
  | .cyp2d6.agree = null | .cyp2d6.consensus = {"result": "*1/*1", "passed_to_pharmcat": true, "reason": "agree"}
  | .cyp2d6.probable = null | .new_since = {"previous_utc": "2026-01-01", "sections": ["cyrius", "cpic"]}'
mut pgx-states "$pgx"' | .sections.pharmcat.state = "unreadable" | .sections.pharmcat.note = "bad json"
  | .sections.cpic.data.parse_failed = true | .sections.cpic.data.non_normal = 0 | .sections.cpic.data.ambiguous = 0
  | .sections.cpic.source = "S_phenotypes.tsv" | .sections.pypgx.state = "missing" | .sections.pypgx.data = {}
  | .cyp2d6.calls = {"PharmCAT": null, "pypgx": null, "Cyrius": null} | .cyp2d6.consensus = null | .cyp2d6.probable = null'
mut pgx-calls "$pgx"' | .sections.pharmcat.state = "failed" | .sections.pharmcat.note = "exit 2"
  | .sections.cpic.state = "missing" | .sections.pypgx.data.comparison.one_tool_only = 3
  | .sections.pypgx.data.depth_check = {"status": "low", "gene_mapq1_share": "0.5", "gene_mapq20_share": "0.4",
  "body_mapq1_share": "0.9", "body_mapq20_share": "0.8"} | .cyp2d6.agree = true
  | .cyp2d6.calls = {"PharmCAT": null, "pypgx": "*1/*2", "Cyrius": "*1/*2"}
  | .cyp2d6.consensus = {"result": "*1/*2", "passed_to_pharmcat": true, "reason": "pypgx & Cyrius agree"}'
hla='.databases.hla_database = "IPD-IMGT/HLA 3.60" | .sections.hla = {"title": "HLA", "state": "ok",
  "source": "hla_t1k/S_hla_genotype.tsv", "file_date": null, "note": null, "data": {"loci": [
  {"gene": "HLA-A", "alleles": ["A*02:01", "A*31:01"], "quality": ["40", "12"]},
  {"gene": "HLA-DRB5", "alleles": ["DRB5*01:01,DRB5*01:02,DRB5*01:03"], "quality": ["0"], "candidates": [12],
  "low_confidence": true, "withheld_from_pharmcat": true},
  {"gene": "HLA-B", "alleles": ["B*57:01", "B*08:01"], "quality": [0, "x"], "low_confidence": true},
  {"gene": "HLA-HFE", "alleles": ["HFE*001:01", "HFE*002:01"], "quality": ["40"], "label": "a <label>"},
  {"gene": "HLA-C", "alleles": []}],
  "derived": {"coeliac": {"status": "derived", "dq2_5": "present", "dq2_5_phase": "likely cis (DRB1*03)",
  "dq2_2": "absent", "dq8": "absent", "dqb1_02_dose": 1}, "alleles": [
  {"allele": "C*06:02", "state": "not_determined", "reason": "HLA-C not typed", "association": "psoriasis",
  "citation": "Nair 2006"}, {"allele": "B*57:01", "state": "present", "association": "abacavir & <x>",
  "citation": "Mallal 2008"}, {"allele": "B*27", "state": "uncertain", "association": "AS", "citation": "C"},
  {"allele": "A*31:01", "state": "absent", "association": "carbamazepine", "citation": "M"}]},
  "optitype": {"genes": [{"gene": "HLA-A", "agree": "yes", "t1k": "A*02:01", "optitype": "A*02:01", "note": "-"},
  {"gene": "HLA-B", "agree": "no", "t1k": "B*57:01", "optitype": "B*58:01", "note": "check <reads>"},
  {"gene": "HLA-C", "agree": "", "t1k": ".", "optitype": "C*07:01", "note": ""}], "reason": ""},
  "kir": {"reason": "", "typed": [{"gene": "KIR2DL1", "alleles": ["KIR2DL1*001", "KIR2DL1*002"], "quality": ["30", ""]}],
  "not_called": ["KIR2DS1", "KIR3DS1"]}}}'
mut hla-ok "$hla"
mut hla-none "$hla"' | .databases.hla_database = null | .sections.hla.state = "stale" | .sections.hla.note = "old"
  | .sections.hla.data.derived = {"coeliac": {"status": "not_determined", "reason": null, "dqb1_02_dose": null},
  "alleles": []} | .sections.hla.data.optitype = {"genes": [], "reason": "OptiType <failed>"}
  | .sections.hla.data.kir = {"reason": "KIR=false", "typed": [], "not_called": []}'
mut hla-failed "$hla"' | .sections.hla.state = "failed" | .sections.hla.note = "T1K <died>" | .sections.hla.data = {}'
mito='.sections.mito = {"title": "Mito", "state": "ok", "source": "mito/S.vcf.gz", "file_date": null, "note": null,
  "data": {"pass": 7, "heteroplasmic": 1, "polyc_heteroplasmic": [{"variant": "chrM:302 A>AC", "af": 0.3}],
  "polyc_tracts": ["m.303-315"], "annotation": {"source": "mito/S_annotated.tsv", "clinvar": "2026-10-01",
  "gnomad": null, "clinvar_plp": 1, "variants": [{"pos": "3243", "ref": "A", "alt": "G", "af": "0.40",
  "clinvar_significance": "Pathogenic", "clinvar_review": "criteria <x>", "clinvar_condition": "MELAS & co",
  "gnomad_af_hom": "0", "gnomad_af_het": "0", "polyc": false}, {"pos": "0", "ref": "A", "alt": "AC", "af": null,
  "clinvar_significance": ".", "gnomad_af_hom": ".", "polyc": true}, {"pos": 0, "ref": "C", "alt": "T",
  "clinvar_significance": "", "gnomad_af_hom": null}]},
  "confirmed": [{"hgvs": "m.3243A>G", "gene": "MT-TL1", "disease": "MELAS", "clinvar_id": "9589",
  "pmids": "2102678;123", "note": "Heteroplasmic.", "af": "0.40", "filter": null}],
  "deletions": {"total": 3, "list_floor_pct": 1.0, "listed": [{"bkp5": 8469, "bkp3": 13447, "length": 4977,
  "level_pct": 10.2, "known": {"disease": "The common deletion", "clinvar_id": "430680", "hgvs": "m.8483_13459del",
  "pmids": "2711184;2541333"}}, {"bkp5": 16000, "bkp3": 500, "length": 1068, "level_pct": 2.0},
  {"bkp5": 100.0, "bkp3": 200, "length": 99, "level_pct": 1.5}]},
  "copies": {"status": "called", "mtdna_copies": 200.0, "depth_ratio": 100.0}}}'
mut mito-ok "$mito"
mut mito-notes "$mito"' | .sections.mito.data = {"pass": 3, "heteroplasmic": 0, "annotation": {"note": "no <ClinVar>"},
  "confirmed_note": "assets/mito_confirmed.tsv is not in assets", "deletions_note": "EKLIPSE=false",
  "copies_note": "no idxstats"}'
mut mito-copies "$mito"' | .sections.mito.data.copies = {"status": "not assessed", "reason": "no_chrM_reads"}
  | .sections.mito.data.annotation.variants = [] | .sections.mito.data.confirmed = [] | .sections.mito.data.deletions = null
  | .sections.mito.state = "stale" | .sections.mito.note = "from <2020>"'
mut mito-failed "$mito"' | .sections.mito.state = "failed" | .sections.mito.note = "Mutect2 died" | .sections.mito.data = {}'
rm -f "$base"

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
PORTED="$PORTED pharmcat cpic cyp2d6 hla mito "
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
# frame) and every chart, then each chapter with no section left (its
# contents entry and its body), and the gene and rsID appendix.
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
    # The charts are og report charts' part, not ported yet.
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
# sections are listed here, as render_report.py writes them.
cat > text.awk << 'AWK'
BEGIN {
  rule = sprintf("%80s", ""); gsub(/ /, "#", rule)
  title["## Quality control"] = "qc"; title["## Variant Calling (DeepVariant)"] = "variants"
  title["## Pharmacogenomics (PharmCAT)"] = "pharmcat"; title["## CPIC Drug-Gene Recommendations"] = "cpic"
  title["## pypgx Pharmacogenomics"] = "cyp2d6"; title["## PGx copy number"] = "cyp2d6"
  title["## CYP2D6 across callers"] = "cyp2d6"; title["## HLA Typing (T1K)"] = "hla"
  title["## Mitochondrial Variants (Mutect2)"] = "mito"
}
{ lines[++n] = $0 }
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
