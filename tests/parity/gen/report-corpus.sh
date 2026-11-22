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
# oracle.sh renders with render_report.py; both reports are compared byte
# for byte, whole.
set -euo pipefail
export LC_ALL=C TZ=UTC
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "${here}/../../.." && pwd)
ORACLE_SHA=dd9290f8fe342fcd010457258adb70082cecbd67
og=${OG:-${root}/rust/target/release/og}

# --- the oracle's modules ------------------------------------------------------------
git -C "$root" cat-file -e "${ORACLE_SHA}^{commit}" 2> /dev/null \
  || git -C "$root" fetch -q --depth=1 origin "$ORACLE_SHA"
mkdir -p py
for m in collect_summary results_index genome_lookup report_charts live_results; do
  git -C "$root" show "${ORACLE_SHA}:bin/${m}.py" > "py/${m}.py"
done
# og's index leaves out the acmg_sf section, whose card the clinical section
# already holds (rust/crates/og-render/src/index.rs), so a gene is not listed
# twice in the closing appendix; the oracle's index does the same.
sed 's/^        if sec.get("state") not in ("ok", "stale"):$/        if key == "acmg_sf" or sec.get("state") not in ("ok", "stale"):/' \
  py/results_index.py > py/results_index.py.new
mv py/results_index.py.new py/results_index.py
grep -q 'key == "acmg_sf"' py/results_index.py || { echo "report-corpus.sh: the acmg_sf skip did not apply" >&2; exit 1; }

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
  # The time of writing, fixed, so a report-render golden of a collected
  # summary (or of a change of one) stays true on every run.
  if [ -s "${case}.json" ]; then
    # sed, not jq: a summary may hold NaN or Infinity on purpose.
    sed -E 's/^( *"generated_utc": )"[^"]*"/\1"2026-01-02 03:04 UTC"/' "${case}.json" > "${case}.json.t"
    mv "${case}.json.t" "${case}.json"
  fi
done < folders/cases.tsv
rm -rf folders

# --- the summary of case "full", changed one way each ---------------------------------
# Its callability section (130 kB of genes) is set to not run, so each change stays small;
# the collected summaries render it.
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
# The opt-in step sections (port-report-optin): every one ran, with data
# down each branch of its lines, then the same in each other state, then
# the branches of empty and partial results.
optin='.sections.chip += {state: "ok", source: "chip/S1_chip.tsv", file_date: "2026-11-01 10:00:00", data: {
    calls: [{gene: "DNMT3A", protein_change: "R882H", vaf: "0.130", alt_reads: "6", depth: "46", note: "hotspot R882 <b>",
             status: "chip"}], other: 2, median_depth: 30.5, median_limit: 0.137, not_covered: ["JAK2", "TP53"]}}
  | .sections.viral += {state: "ok", source: "viral/S1_viral.tsv", data: {hhv6: "inherited integrated HHV-6B likely",
    human_depth: "30.00", genomes: 24, found: [
      {virus: "HHV-6B", accession: "NC_000898.1", reads: "400", breadth: "1.0000", relative_depth: "0.500",
       copies_per_cell: "1.0000", call: "inherited integrated HHV-6 likely", human_depth: "30.00", note: "."},
      {virus: "Human gammaherpesvirus 4", accession: "NC_007605.1", reads: "30", breadth: "0.9950",
       relative_depth: "0.010", copies_per_cell: ".", call: "viral DNA found", human_depth: "30.00", note: "x & y"}]}}
  | .sections.mei += {state: "ok", source: "mei/S1_mei_sites.tsv", data: {calls: 2, by_element: {Alu: 1, SVA: 1}, sites: [
      {gene: "FKTN", element: "SVA", condition: "Fukuyama congenital muscular dystrophy", inheritance: "autosomal recessive",
       status: "seen", mean_depth: "30.0", call: "chr9:105639660 SVA_F 1 copy PASS"},
      {gene: "TAF1", element: "SVA", condition: "X-linked dystonia-parkinsonism", inheritance: "X-linked",
       status: "not callable", mean_depth: ".", call: "."}]}}
  | .sections.mosaic += {state: "ok", source: "mosaic/S1_mocha_summary.tsv", data: {status: "assessed", message: ".",
    depth: "30.0", het_sites: "2090000", het_per_10mb: "3000", chromosomes: "chr1-chr22", floor_read: "0.020",
    floor_min: "0.05", floor: "0.050", calls_mocha: "3", calls_reported: "1", phasing: "input", not_reported: 2, calls: [
      {chrom: "chr20", start: "60000", end: "20000000", length_mb: "19.94", arms: "20p, to the p telomere",
       type: "CN-LOH", cell_fraction: "0.150", n_hets: "6600", lod_baf_phase: "900", reported: "yes", reason: "."},
      {chrom: "chr1", start: "500", end: "500", length_mb: "0.00", arms: "1p", type: "loss", cell_fraction: "0.2",
       n_hets: "1", lod_baf_phase: "9", reported: "yes", reason: "."}]}}
  | .sections.vntyper += {state: "ok", source: "vntyper/S1_vntyper.tsv", data: {status: "positive",
    confidence: "High_Precision*", variant: "59dupC", flag: "Not flagged", alt_depth: "416", region_depth: "7110",
    depth_score: "0.0585", vntr_mean_depth: "1075.32", vntr_uncovered_pct: "0.00", coverage_qc: "PASS",
    sensitivity: "At 1075x the 100x-or-more band: 160 of 161 runs (99%).", message: "A frameshift <dupC>.",
    advice: "Confirm with SNaPshot."}}
  | .sections.locityper += {state: "ok", source: "locityper/S1_locityper.tsv", data: {called: 2, passed: 1, loci: [
      {locus: "HLA-DQB1", region: "chr6:32652604-32669077", status: "called", genotype: "HG00733.1,NA19240.2",
       haplotype_1: "HG00733.1", haplotype_2: "NA19240.2", quality: "27.8", total_reads: "900",
       unexplained_reads: "12", weighted_distance: "0.200", warnings: ".", filter: "pass", reason: ".",
       primary: "T1K (step 08): DQB1*06:02 + DQB1*03:01"},
      {locus: "CYP2D6", region: "chr22:42122285-42152960", status: "called", genotype: "a,b", haplotype_1: "a",
       haplotype_2: "b", quality: "9", total_reads: "3000", unexplained_reads: "1200", weighted_distance: "50.00",
       warnings: ".", filter: "fail", reason: "weighted distance 50.00", primary: "pypgx (step 32): *1/*2"},
      {locus: "FCGR", region: "chr1:161506762-161681047", status: "not called", genotype: ".", haplotype_1: ".",
       haplotype_2: ".", quality: ".", total_reads: ".", unexplained_reads: ".", weighted_distance: ".", warnings: ".",
       filter: ".", reason: "Locityper wrote no result", primary: "no other step calls FCGR"},
      {locus: "KIR", region: "chr19:1-2", status: "called", genotype: "k,k", haplotype_1: "k", haplotype_2: "k",
       quality: "5", total_reads: "10", unexplained_reads: "1", weighted_distance: "1", warnings: ".", filter: "fail",
       reason: "."}]}}
  | .sections.exomiser += {state: "ok", source: "exomiser/S1_exomiser_top.tsv", data: {
    terms: [{id: "HP:0001156", label: "Brachydactyly"}, {id: "HP:0000001", label: null}], genes: [
      {rank: "1", gene: "GENEA", moi: "AD", exomiser_score: "0.998", phenotype_score: "0.800", variants: [
        {variant: "chr10:123256215 T>G", genotype: "0/1", consequence: "missense_variant",
         hgvs: "GENEA:ENST1:c.1688A>C:p.(Glu563Ala)", acmg: "LIKELY_PATHOGENIC", disease: "Synthetic syndrome one"}]},
      {rank: "2", gene: "GENEB", moi: "AR", exomiser_score: "0.9", phenotype_score: "0.8", variants: [
        {variant: "chr6:1000 G>A", genotype: "0/1", consequence: "splice_donor_variant", hgvs: "c.10+1G>A",
         acmg: "UNCERTAIN_SIGNIFICANCE", disease: "."},
        {variant: "chr6:2000 C>T + chr6:900 A>C", genotype: "1|0", consequence: "stop_gained", hgvs: "p.(Gln100*)",
         acmg: "PATHOGENIC", disease: "."}]},
      {rank: "3", gene: "GENEC", moi: "XR", exomiser_score: "0.5", phenotype_score: "0.4", variants: [
        {variant: "no locus", genotype: "1", consequence: "missense_variant", hgvs: "p.?", acmg: "VUS",
         disease: "D"}]}]}}
  | .sections.telomere += {state: "ok", source: "telomere/S1_summary.tsv",
    data: {tel_content: "1.5", intratelomeric_reads: 1234567, depth: 30.0, read_floor: 2000, low_reads: true}}
  | .sections.roh += {state: "ok", source: "vcf/S1_roh.txt", data: {segments: 3, total_mb: 37.5, largest_mb: 27,
    autosomal_over_5mb: [{region: "chr2:10000001-17000000", mb: 7.0},
                         {region: "chr15:23000001-50000000", mb: 27.0, masked_mb: 1.5}],
    masked_mb: 46.5, sex_chromosome_segments: 1, long_mb: 34.0, f_roh: 0.0128, relationship: "first cousins",
    upd_flag: {chrom: "chr15", mb: 27.0, share: 0.79, imprinted: [{region: "15q11-q13 (SNRPN, UBE3A)",
      maternal_upd: "Prader-Willi", paternal_upd: "Angelman", pmid: "22237428"}]}}}'
OPTIN_KEYS='"chip", "viral", "mei", "mosaic", "vntyper", "locityper", "exomiser", "telomere", "roh"'
mut optin-ok "$optin"
for st in stale failed unreadable odd; do
  mut "optin-$st" "$optin | reduce ($OPTIN_KEYS) as \$k (.; .sections[\$k].state = \"$st\" | .sections[\$k].note = \"why <\(\$k)>\")"
done
mut optin-since "$optin"' | .new_since = {"previous_utc": "2026-01-01T00:00:00Z", "sections": ["roh", "chip", "mosaic"]}'
mut optin-empty "$optin"' | .sections.chip.data = {calls: [], other: 0, median_depth: null, median_limit: null, not_covered: []}
  | .sections.viral.data = {hhv6: "no integrated HHV-6", human_depth: "30.00", found: [], genomes: 24}
  | .sections.mei.data = {calls: 0, by_element: {}, sites: []}
  | .sections.mosaic.data = {status: "not assessed", message: "phased BAF autocorrelation 0.05 over 0.03",
    calls: [], not_reported: 0}
  | .sections.vntyper.data = {status: "negative", variant: ".", coverage_qc: "PASS", sensitivity: "s", message: "m",
    advice: "a"}
  | .sections.locityper.data = {called: 0, passed: 0, loci: []}
  | .sections.exomiser.data = {terms: [], genes: []}
  | .sections.telomere.data = {tel_content: "334", read_check: "depth unknown: not checked"}
  | .sections.roh.data = {segments: 0, total_mb: 0.0, largest_mb: 0.0, autosomal_over_5mb: [], masked_mb: 0.0,
    sex_chromosome_segments: 0, long_mb: 0.0, f_roh: null, relationship: null, upd_flag: null}'
mut optin-partial "$optin"' | .sections.mosaic.data.calls = [] | .sections.mosaic.data.status = "assessed"
  | .sections.viral.data.hhv6 = "not assessed (no human depth)"
  | .sections.vntyper.data = {status: "possible", variant: "59dupC", confidence: "Low_Precision"}
  | .sections.telomere.data = {tel_content: 2, intratelomeric_reads: 40}
  | .sections.roh.data.upd_flag.imprinted = [] | .sections.roh.data.relationship = null | .sections.roh.data.f_roh = 1
  | .sections.chip.data.median_limit = 0.02 | .sections.chip.data.calls += [{gene: "TET2", protein_change: "Q1200fs",
    vaf: "0.2", alt_reads: "9", depth: "45", note: "truncating"}]'
# Follow-up 454: no unplaced unmapped reads, so no viral genome was looked at.
# og says so; render_report.py said no read was found (a golden, so the
# summary's clock is fixed and only the viral section changes).
mut optin-viral-na '.generated_utc = "2026-01-02 03:04 UTC" | .sections.viral += {state: "ok",
  source: "viral/S1_viral.tsv", data: {hhv6: "not assessed (the BAM has no unplaced unmapped reads)", found: [],
  genomes: 24}}'
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

# The clinical sections: case cb-full, which holds every one of them but the couple
# screen and callability, changed one way each. (all and all-history, the only
# collected summaries with a couple screen, run as text-only cases.)
[ -s cb-full.json ] || { echo "report-corpus.sh: og summary collect wrote no cb-full.json" >&2; exit 1; }
cmut() { jq -c "$2" cb-full.json > "m-$1.json"; }
couple='{"title": "Couple carrier screen", "step": "39", "state": "ok", "note": null, "file_date": "2021-03-04 05:06 UTC",
  "data": {"partner": "P<2>", "shared": [{"gene": "CFTR", "condition": "Cystic fibrosis", "risk": "1 in 4",
  "own": "c.1521_1523del (P)", "partner_findings": "c.350G>A (LP)"}], "x_linked": [{"gene": "DMD",
  "condition": "Duchenne", "risk": "1 in 4 sons", "own": "c.10G>T (P)", "partner_findings": "."}], "one_partner": 2}}'
cmut couple ".sections.couple = $couple"
cmut couple-empty '.sections.couple = {"title": "Couple carrier screen", "state": "ok", "note": null, "data": {}}'
for st in stale failed unreadable; do
  keys='"clinvar", "clinical", "cpsr", "slivar", "carrier", "couple"'
  [ "$st" = stale ] || keys="$keys, \"callability\""
  cmut "clinical-$st" ".sections.couple = $couple | reduce ($keys) as \$k (.;
    .sections[\$k].state = \"$st\" | .sections[\$k].note = \"note of \(\$k) <&>\")"
done
# shellcheck disable=SC2016  # jq's $h, not the shell's
cmut clinvar-many '.sections.clinvar.data.hits as $h | .sections.clinvar.data.hits = ($h + $h + $h + $h + $h + $h)
  | .sections.clinvar.data.hits[55].cdna = {"check": "flagged", "note": "clipped <x>", "alt_reads": 4, "alt_clipped": 3,
    "ref_reads": 10, "ref_clipped": 0, "construct": "pcDNA3"}
  | .sections.clinvar.data.hits[57].cdna = {"check": "clear", "note": "", "alt_reads": 4, "alt_clipped": 0,
    "ref_reads": 10, "ref_clipped": 0, "construct": "EGFP"}
  | .sections.clinvar.data.hits[58].cdna = {"check": "not_checked", "note": "", "alt_reads": null, "alt_clipped": null,
    "ref_reads": null, "ref_clipped": null, "construct": "lentiviral"}'
cmut clinvar-none '.sections.clinvar.data.hits = [] | .sections.clinvar.data.count = 0 | del(.sections.clinvar.data.deep_intronic)
  | .sections.clinvar.data.cdna_scan = {"note": "not run: no BAM", "loci": []} | del(.sections.clinvar.data.frequency)
  | del(.sections.clinvar.data.low_confidence) | del(.sections.clinvar.data.cdna_flagged)'
cmut clinvar-freq-note '.sections.clinvar.data.frequency = {"note": "no VEP output"} | .sections.clinvar.data.call_floor = null'
cmut cpsr-noown '.sections.cpsr.data.own_class_column = false | .sections.cpsr.data.classification = {}
  | .sections.cpsr.source = "cpsr/S.cpsr.grch38.html" | del(.sections.cpsr.data.cdna_flagged)'
cmut cpsr-none '.sections.cpsr.data = {"html_report": false}'
cmut slivar-phase '.sections.slivar.data.compound_het_phase.trans = 2.5 | .sections.slivar.data.compound_het_trans = []'
cmut slivar-unphased 'del(.sections.slivar.data.compound_het_phase)'
cmut sf-notes '.sections.clinical.data.acmg_sf.two_hit_note = "carrier list <unreadable>"
  | .sections.clinical.data.acmg_sf.structural_checked = true | .sections.clinical.data.acmg_sf.high_impact_weak = 2
  | .sections.clinical.data.acmg_sf.high_impact_frequency_unknown = 1 | del(.sections.clinical.data.acmg_sf.genes_on_list)'
cmut clinical-escape '.sections.clinvar.data.hits[0].gene = "A<&>'"'"'\"" | .sections.carrier.data.carriers[0].condition = "x <y> & z"
  | .sections.clinical.data.acmg_sf.clinvar_hits[0].gene = "<b>"'
cmut since-clinical '.new_since = {"previous_utc": "2026-01-01T00:00:00Z", "sections": ["clinvar", "carrier", "slivar"]}'
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
# prs-raw is a run without the panel, so no score has a percentile note: a
# withheld note (c-full's T2D row) would make og print the withheld line,
# which the oracle's commit predates (the cargo test prs_withheld.rs covers it).
amut prs-raw '.sections.prs.data.adjusted = false | del(.sections.prs.data.ancestry)
  | .sections.prs.data.scores |= map(.percentile = "" | .evaluated_in = null | .percentile_note = "")'
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
# The structural-variant, mosaic sex-loss and raw-data export sections:
# each call of step 33's sex-loss check, chrX with an allele fraction, the
# no-call reasons with tied counts, and the SV sections stale.
smut() { jq -c "$2" "$abase" > "s-$1.json"; }
sl='.sections.mosaic_sex_loss = {"title": "Mosaic loss of chromosome Y or X", "step": "33", "state": "ok",
  "source": "qc/S_sample_qc.tsv", "file_date": "2021-03-04 05:06 UTC", "note": null, "data": {'
smut sexloss-detected "$sl"'"chrom": "chrY", "call": "detected", "reason": "depth below the floor",
  "depth_ratio": 0.7123, "loss_depth": 0.2877, "loss_depth_low": 0.25, "loss_depth_high": 0.33, "loss_af": null}}'
smut sexloss-chrx "$sl"'"chrom": "chrX", "call": "detected", "reason": "allele fraction <off> & depth",
  "depth_ratio": null, "loss_depth": null, "loss_af": 0.125, "loss_af_low": 0.1, "loss_af_high": null, "x_het_sites": 40}}'
smut sexloss-odd "$sl"'"chrom": "", "call": "maybe <b>", "reason": null, "loss_depth": 0.05}}'
smut sexloss-not-estimated "$sl"'"chrom": "chrY", "call": "not_estimated", "reason": "too few reads on chrY"}}'
smut sexloss-gain-stale "$sl"'"chrom": "chrY", "call": "gain", "reason": "depth above 1"}} | .sections.mosaic_sex_loss.state = "stale"
  | .sections.mosaic_sex_loss.note = "from an older run <b>"'
smut sexloss-failed "$sl"'}} | .sections.mosaic_sex_loss.state = "failed" | .sections.mosaic_sex_loss.note = "step 33 failed"'
smut consumer-reasons '.sections.consumer_export.data.no_call_reasons = {"zeta": 5, "alpha": 5, "low <depth>": 9, "one": 1}
  | .sections.consumer_export.data.from_reference_blocks = 1234 | .sections.consumer_export.data.file_ancestrydna = "S_ancestry.txt"'
smut consumer-stale '.sections.consumer_export.state = "stale" | .sections.consumer_export.note = "older export"
  | .sections.consumer_export.data.no_call_reasons = {} | .sections.consumer_export.data.sites = 0'
# shellcheck disable=SC2016  # jq's $k, not the shell's
smut sv-stale 'reduce ("manta", "delly", "sv_consensus", "annotsv") as $k (.; .sections[$k].state = "stale"
  | .sections[$k].note = "why: <b>&'"'"'\"")'
smut sv-mixed '.sections.manta.state = "failed" | .sections.manta.note = "manta failed" | .sections.delly.state = "missing"
  | .sections.annotsv.state = "unreadable" | .sections.annotsv.note = "bad tsv"'
rm -f "$abase"

# --- the oracle -----------------------------------------------------------------------
cat > oracle.sh << 'ORACLE'
#!/bin/sh
# The oracle of og report render (tests/parity/oracle.ref): render_report.py
# with the modules of py/. Both reports are compared whole: every section,
# the charts in the HTML report and the closing gene and rsID appendix.
# A PRS row whose counts the summary has as null reads "not available", as
# og writes it; render_report.py writes "None" there (and stops its HTML
# report on that row), so the text report's "None" is read as og's words.
# og's frame for a phone and a keyboard, which render_report.py never had,
# is added to the Python's reports here, line for line as og writes it
# (rust/crates/og-render: html.rs, text.rs and READING in css.rs): the
# research-only line under the title of both reports, the skip link, the
# contents in a details element, main#main, each table in a scroll region
# named by its card's title (the manifest's "Run manifest"; in a card with
# several tables, plus the subheading above the table, else its number) and
# the CSS. The adapter only adds lines and wraps tables: every byte of each
# section the Python writes stays as it wrote it. An output with no golden
# is compared directly with this adapted oracle. An output with a golden
# (tests/parity/golden/report-render) is still compared with its golden,
# and each golden is still checked against this adapted oracle: a golden
# the oracle also writes is stale and fails.
#   sh ./oracle.sh PROG --summary F -o R.txt -o R.html
set -e
prog=$1
shift
d=$(mktemp -d)
cp py/*.py "$d"/
cp "$prog" "$d/render_report.py"
PYTHONPATH=$d python3 "$d/render_report.py" "$@"
rm -rf "$d"
frame_text() {
  awk '{ print }
$0 == "  GENOMICS ANALYSIS REPORT" {
  print "  Research software, not a clinical report: confirm any finding in a"
  print "  certified laboratory before acting on it."
}' "$1" > "$1.kept"
  mv "$1.kept" "$1"
}
frame_html() {
  awk 'function rep(s, a, b,   i, out) {
  out = ""
  while ((i = index(s, a)) > 0) { out = out substr(s, 1, i - 1) b; s = substr(s, i + length(a)) }
  return out s
}
NR == FNR {
  if ($0 ~ /^    <h2>.*<\/h2>$/ || index($0, "<summary>Run manifest (run_manifest.tsv)</summary>")) { c++; sh = "" }
  if ($0 ~ /^    <h3>.*<\/h3>$/) sh = substr($0, 9, length($0) - 13)
  t = $0
  k = gsub(/<table>/, "", t)
  if (k > 0) { if (sh != "") heads[c, sh]++; sh = ""; tables[c] += k }
  next
}
{ line = $0 }
line == "</style>" {
  print "  .skip { position: absolute; left: 8px; top: -60px; z-index: 10; padding: 8px 14px; border-radius: 6px;"
  print "           background: var(--card); color: var(--accent); font-weight: 600; }"
  print "  .skip:focus { top: 8px; }"
  print "  main:focus { outline: none; }"
  print "  .header .research-only { font-weight: 600; margin-bottom: 8px; }"
  print "  main { overflow-wrap: anywhere; }"
  print "  .sec, .card { min-width: 0; }"
  print "  .tscroll { overflow-x: auto; max-width: 100%; }"
  print "  .tscroll:focus-visible { outline: 2px solid var(--accent); }"
  print "  .toc-fold > summary { cursor: pointer; color: var(--accent); font-weight: 600; padding: 8px 0; }"
  print "  @media (max-width: 900px) { .toc-fold .toc-title { display: none; } }"
  print "  @media (min-width: 901px) { .toc-fold { position: sticky; top: 0; align-self: start; max-height: 100vh; overflow-y: auto; } }"
  print "  @supports selector(::details-content) {"
  print "    @media (min-width: 901px) {"
  print "      .toc-fold > summary { display: none; }"
  print "      .toc-fold::details-content { content-visibility: visible; display: block; }"
  print "    }"
  print "  }"
  print "  @media print { .toc-fold, .skip { display: none; } .tscroll { overflow: visible; } }"
}
line == "<body>" { print; print "<a class=\"skip\" href=\"#main\">Skip to results</a>"; next }
line == "<nav class=\"toc\" aria-label=\"Contents\">" { print "<details class=\"toc-fold\">"; print "<summary>Contents</summary>"; print; next }
line == "</nav>" { print; print "</details>"; next }
line == "<main>" { print "<main id=\"main\" tabindex=\"-1\">"; next }
line == "  <h1>Genomic Analysis Report</h1>" {
  print
  print "  <p class=\"research-only\">Research software, not a clinical report: confirm any finding in a certified laboratory before acting on it.</p>"
  next
}
line ~ /^    <h2>.*<\/h2>$/ { card++; n = 0; subh = ""; label = substr(line, 9, length(line) - 13) }
index(line, "<summary>Run manifest (run_manifest.tsv)</summary>") { card++; n = 0; subh = ""; label = "Run manifest" }
line ~ /^    <h3>.*<\/h3>$/ { subh = substr(line, 9, length(line) - 13) }
index(line, "table>") {
  out = ""
  s = line
  while ((i = index(s, "<table>")) > 0) {
    n++
    name = label
    if (tables[card] > 1) name = (subh != "" && heads[card, subh] == 1 ? label ": " subh : label ", table " n)
    subh = ""
    out = out substr(s, 1, i - 1) "<div class=\"tscroll\" role=\"region\" aria-label=\"" name "\" tabindex=\"0\"><table>"
    s = substr(s, i + 7)
  }
  line = rep(out s, "</table>", "</table></div>")
}
{ print line }' "$1" "$1" > "$1.kept"
  mv "$1.kept" "$1"
}
while [ $# -gt 0 ]; do
  case "$1" in
    -o)
      case "$2" in
        *.txt)
          sed 's| None (None/None variants matched) | not available (not available/not available variants matched) |' "$2" > "$2.kept"
          mv "$2.kept" "$2"
          frame_text "$2" ;;
        *.html) [ ! -e "$2" ] || frame_html "$2" ;;
      esac
      shift ;;
  esac
  shift
done
ORACLE
