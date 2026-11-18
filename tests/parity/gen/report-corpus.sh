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
ORACLE_SHA=65095e975bea01c74ed92a4f1b8c1bd00d05caaa
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
# The sections of the traits, copy-number and repeat unit: every branch of
# each, in each state (the summary-collect folders never hold an
# ExpansionHunter result, and their blood groups render in text only).
sec() { printf '.sections.%s = {"title": "%s", "state": "%s", "note": %s, "file_date": "2021-03-04 05:06 UTC", "source": null, "data": %s}' "$1" "$1" "$2" "$3" "$4"; }
eh='{"records": 6, "stranger": true, "no_status": 2, "unsized": ["DMPK", "RFC1"], "intermediate": ["HTT"],
  "not_graded": {"MUC1": "reference >300 bp", "CEL": "<none>"},
  "flagged": [{"locus": "RFC1", "status": "pathogenic", "rfc1_motif": {"verdict": "expanded", "note": "AAGGG & <AAAGG>"}},
    {"locus": "HTT", "status": "premutation", "motif": {"verdict": "pure", "note": "CAG"}},
    {"locus": "ATXN1", "status": "full_mutation", "repeat_count": "60"},
    {"locus": "FXN", "status": "pathogenic", "repeat_count": "900/12", "allele_sizing": "in_repeat/spanning", "ci": "700-1000", "rfc1_motif": {"verdict": "x"}}],
  "loci": [{"repeat_count": "17/40", "locus": "HTT", "allele_sizing": "spanning/in_repeat", "ci": "17-17/38-45",
      "status": "premutation", "sizing": "in_repeat", "normal_max": 26, "zone": "reduced penetrance",
      "zone_label": "36 to 39 repeats", "zone_citation": "PMID 1", "motif": {"verdict": "pure", "note": "CAG"}},
    {"repeat_count": "11/11", "locus": "RFC1", "allele_sizing": "in_repeat/in_repeat", "ci": "x-1", "status": "normal",
      "sizing": "flanking_only", "normal_max": null},
    {"repeat_count": "5/7", "locus": "DMPK", "allele_sizing": "in_repeat", "ci": "12-20", "status": "normal",
      "sizing": "low_confidence"},
    {"locus": "MUC1", "repeat_count": ".", "status": "", "not_graded": "reference >300 bp", "sizing": "no_call", "ci": "."},
    {"locus": "AR", "repeat_count": "22", "status": "normal", "sizing": "sized", "ci": ""},
    {"repeat_count": "18/41", "locus": "HTT", "status": "normal", "sizing": "unknown"}],
  "key_loci": [{"locus": "C9ORF72"}, {"locus": "HTT"}, {"locus": "FMR1"}],
  "readings": [{"locus": "CFTR", "repeat": "TG", "genotype": "11/12", "status": "read", "reading": "a <b> reading",
      "source": "PMID 2", "five_t": true},
    {"locus": "UGT1A1", "repeat": "TA", "genotype": ".", "status": "not_read", "reading": "r2", "source": "s2",
      "agreement": "differs"},
    {"locus": "SLC6A4", "repeat": "5-HTTLPR", "genotype": "L/S", "status": "read", "reading": "r3", "source": "s3"}]}'
lpa='{"status": "too little depth", "message": "8.0x over the exons", "advice": "Measure Lp(a) once.",
  "flag": "none", "flag_text": "No genetic sign.",
  "variants": [{"rsid": "rs10455872", "genotype": "./.", "status": "no_call", "risk_alleles": null},
    {"rsid": "rs3798220", "genotype": "C/C", "status": "called", "risk_alleles": 2},
    {"rsid": "rs0", "genotype": "T/C", "status": "called", "risk_alleles": 1.0}]}'
lpa_null='{"status": "estimated", "kiv2_copies": "19.8", "ci95_low": "17.9", "ci95_high": "22.0", "class": "low",
  "message": "about 20", "advice": ".", "flag": "high", "flag_text": "Genetic signs.", "variants": [],
  "null_allele": {"rsid": "", "status": "uncertain", "alt_reads": "2", "depth": "150", "copies": "0.4", "message": "m"}}'
blood='{"abo": "AB & \"O\"", "rhd": "RhD négatif",
  "systems": [{"system": "ABO", "phenotype": "AB", "genotype": "A1/B", "confidence": "high"},
    {"system": "Kell<K>", "phenotype": "K-k+ é", "genotype": "KEL*02/KEL*02", "confidence": "medium"},
    {"system": "a-very-long-system", "phenotype": "a phenotype longer than thirty-four characters", "genotype": "KEL*01.07",
      "confidence": "low"}]}'
traits='{"called": 3, "not_covered": 2, "other_allele": 1,
  "rows": [{"rsid": "rs4988235", "gene": "MCM6", "trait": "Lactase <persistence>", "genotype": "A/G", "status": "called",
      "reading": "Likely persistent", "ancestry": "yes", "note": "European & others", "pmids": "1,2"},
    {"rsid": "rs1", "gene": "G1", "trait": "t1", "genotype": ".", "status": "not_covered", "reading": ".",
      "ancestry": "no", "note": ".", "pmids": "3"},
    {"rsid": "rs2", "gene": "G2", "trait": "t2", "genotype": "./.", "status": "not_covered", "reading": ".",
      "ancestry": "no", "note": ".", "pmids": "4"},
    {"rsid": "rs3", "gene": "G3", "trait": "t3", "genotype": "C/T", "status": "other_allele", "reading": "other",
      "ancestry": "no", "note": "n", "pmids": "5"}],
  "rules": [{"rule": "APOL1", "gene": "APOL1", "trait": "Kidney disease", "rsids": "rs73885319,rs60910145",
      "genotypes": "A/G,T/G", "phase": "read_backed", "result": "G0/G1", "status": "called", "reading": "one G1",
      "ancestry": "yes", "note": "African", "pmids": "6",
      "g2": {"result": "G1/G2", "text": "G2 read from the VCF.", "reading": "two risk variants <high>"}},
    {"rule": "HFE", "gene": "HFE", "trait": "Iron", "rsids": "rs1800562,rs1799945", "genotypes": "G/G,C/C",
      "phase": "not_needed", "result": "wt/wt", "status": "called", "reading": "no variant", "ancestry": "no",
      "note": ".", "pmids": "7"},
    {"rule": "GC", "gene": "GC", "trait": "Vitamin D", "rsids": "rs7041", "genotypes": ".", "phase": "unphased",
      "result": ".", "status": "not_covered", "reading": ".", "ancestry": "no", "note": ".", "pmids": "8",
      "g2": {"result": ".", "text": "not read", "reading": "r"}},
    {"rule": "MAPT", "gene": "MAPT", "trait": "t", "rsids": "rs1", "genotypes": "A/G", "phase": "other",
      "result": "H1/H2", "status": "called", "reading": "r", "ancestry": "no", "note": "n", "pmids": "9"}],
  "apoe": [{"rsid": "rs429358+rs7412", "gene": "APOE", "trait": "APOE e4", "genotype": "e3/e4", "status": "called",
      "reading": "one e4", "ancestry": "yes", "note": "risk factor", "pmids": "10"}]}'
paraphase='{"genes": [{"gene": "SMN1", "copy_number": "2", "status": "called", "detail": "."},
    {"gene": "PMS2", "copy_number": "no call", "status": "no_call", "detail": "one haplotype <x>"},
    {"gene": "HBA1", "copy_number": "1", "status": "called", "detail": ""}],
  "other_regions": 0, "other_regions_called": 0, "variants": {"smn1": 1.5, "hba": 2},
  "smn1": {"paraphase": null, "parascopy": {"text": "2 (quality 30)"}, "agree": null}}'
parascopy='{"smn1": {"text": "2 (quality 30, filter PASS)"}, "smn2": {}, "agcn": "4", "agcn_qual": "40",
  "agcn_filter": "PASS", "reading": "two SMN1 copies",
  "smn1_marker": {"rsid": null, "status": "present", "alt_reads": "10", "depth": "40", "copies": "1.0",
    "message": "a higher chance of 2+0"},
  "hba": {"status": "not assessed", "message": "not assessed: flank depth under 8x"},
  "cyp21a2": {"status": "called", "cyp21a2_copies": "2", "module_total": "4", "expected_total": "4",
    "finding": "conversion", "conversion_sites": "In2G", "reading": "a possible conversion"},
  "dosage": [{"locus": "HP", "status": "called", "copies": "3", "estimate": "3.02", "reading": "."},
    {"locus": "CCL3L1", "status": "called", "copies": "5", "estimate": "4.9", "reading": "more copies"},
    {"locus": "LCE3BC", "status": "not assessed", "copies": ".", "estimate": ".", "reason": "low_control_depth",
      "reading": "."}],
  "cn_intervals": {"rows": 4, "by_status": {"not assessed": 1, "called": 3}, "changes": []},
  "loci": [{"genes": "AMY1A,AMY1B,AMY1C", "text": "6 (quality 25)", "status": "called", "reading": "six copies"},
    {"genes": "C4A,C4B", "text": "not assessed (no read)", "status": "not assessed", "reading": "r"},
    {"genes": "CFH", "text": "4", "status": "called", "reading": "."}]}'
parascopy_b='{"smn1": {"text": "1"}, "cyp21a2": {"status": "not assessed", "reason": "psv_low_depth"},
  "dosage": [], "cn_intervals": {"rows": 2, "by_status": {}, "changes": null}}'
parascopy_c='{"cyp21a2": {"status": "not assessed", "message": "", "reason": ""}, "cn_intervals": {"rows": 3,
  "by_status": {"not assessed": 3}}, "hba": {"status": "called", "alpha_globin_genes": "4"}}'
mut genome-ok "$(sec expansions ok null "$eh") | $(sec lpa ok null "$lpa") | $(sec bloodtype ok null "$blood")
  | $(sec traits ok null "$traits") | $(sec paralogs ok null "$paraphase") | $(sec parascopy ok null "$parascopy")
  | .sections.paralogs.data.smn1.paraphase = \"2\" | .sections.paralogs.data.smn1.agree = \"agree\""
mut genome-stale "$(sec expansions stale '"from 2020"' "$eh") | $(sec lpa stale '"from <2020>"' "$lpa_null")
  | $(sec bloodtype stale '"old"' "$blood") | $(sec traits stale '"old"' "$traits")
  | $(sec paralogs stale '"old"' "$paraphase") | $(sec parascopy stale '"old"' "$parascopy_b")
  | .sections.paralogs.data.smn1.agree = \"disagree\"
  | .sample = \"A&B <x>\" | .new_since = {\"previous_utc\": \"2026-01-01T00:00:00Z\",
    \"sections\": [\"traits\", \"parascopy\", \"expansions\", \"lpa\", \"bloodtype\"]}"
mut genome-failed "$(sec expansions failed '"exit 1"' '{}') | $(sec lpa failed '"exit 2"' '{}')
  | $(sec bloodtype unreadable '"bad table"' '{}') | $(sec traits failed '"exit 3"' '{}')
  | $(sec paralogs failed '"exit 4"' '{}') | $(sec parascopy unreadable '"bad"' '{}')"
mut genome-other "$(sec expansions weird null '{}') | $(sec traits weird null '{}') | $(sec parascopy weird null "$parascopy_c")
  | $(sec paralogs missing null '{}')"
mut genome-stranger "$(sec expansions ok null '{"records": 0, "stranger": true, "flagged": [], "loci": []}')
  | $(sec parascopy failed '"exit 1"' '{}')"
mut genome-outdated "$(sec expansions ok null '{"records": 1, "stranger": true, "stranger_outdated": true,
  "loci": [{"locus": "AR", "repeat_count": "22"}]}') | $(sec parascopy ok null "$parascopy_c")
  | .sections.parascopy.data.cn_intervals.changes = [{\"name\": \"22q11.2 A-D\", \"locus\": \"x\",
    \"chrom\": \"chr22\", \"start\": \"18912230\", \"end\": \"21465659\", \"copies\": \"1\", \"expected\": \"2\",
    \"change\": \"deletion\"}, {\"name\": \"AZFc <b&r>\", \"locus\": \"y\", \"copies\": \"2\", \"expected\": \"1\",
    \"change\": \"DUPLICATION\"}]"
mut genome-no-stranger "$(sec expansions ok null '{"records": 2, "stranger": false, "no_status": 2,
  "not_graded": {}, "unsized": [], "intermediate": []}') | $(sec traits ok null '{"called": 0, "not_covered": 0,
  "other_allele": 0, "rows": [], "rules": [], "apoe": []}') | $(sec paralogs ok null "$paraphase")"
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
PORTED=" qc variants not-run not-assessed expansions lpa paralogs traits-panel bloodtype "
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
# frame), then each chapter with no section left (its contents entry and its
# body), and the gene and rsID appendix.
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
  title["## Repeat Expansions (ExpansionHunter)"] = "expansions"
  title["## Lp(a): LPA KIV-2 Repeat Count and Risk Variants"] = "lpa"
  title["## Paralog Genes (Paraphase, PacBio HiFi)"] = "paralogs"
  title["## SMN1/SMN2 Copy Number (Parascopy, short reads)"] = "paralogs"
  title["## Paralog Copy Number"] = "paralogs"
  title["## Traits Panel (step 42)"] = "traits-panel"
  title["## Blood Groups (predicted from the genome)"] = "bloodtype"
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
