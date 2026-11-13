#!/usr/bin/env python3
"""og pgx outside-calls: what reaches PharmCAT as an outside call.

CYP2D6 reaches PharmCAT only when pypgx and a second caller (Cyrius) give the
same diplotype and the depth check passed; HLA-A and HLA-B come from T1K.
Each case runs og on synthetic inputs (tests/fixtures/pgx/) and compares
both files it writes, byte for byte:

  1. disagree      pypgx and Cyrius differ: no outside call, 'indeterminate'
  2. one missing   pypgx only (Cyrius not run): no outside call
  3. agree         both give the same diplotype, depth ok: the call is passed
  4. depth         both agree, but the depth check found multi-mapped reads
  5. no call       Cyrius ran but made no call: one caller only
  6. HLA           T1K's HLA-A and HLA-B, truncated to two fields
  7. HLA edges     one allele (read as homozygous); an allele of quality 0

and what step 27 (og pgx cpic-report --consensus) makes of it:

  8. CYP2D6 held back: the CPIC report says so, names the drugs it affects,
     and drops the pypgx-only warning it prints without the table
  9. HLA passed on: the CPIC report lists it in the outside-call section and
     marks the gene PharmCAT reports with callSource OUTSIDE
 10. a stale report: PharmCAT's report has outside calls the consensus does
     not pass on (another HLA-B, a CYP2D6 now held back, or no readable
     table at all); they give no drug guidance and the report says to rerun
     step 07

and step 32's summary (og pgx pypgx-summary) and the probable CYP2D6 call
the reports show:

 11. (the star alleles of a diplotype and the defining-variant evidence of
     one call are og's own unit tests: cargo test -p og-steps pypgx_summary)
 13. pypgx alone with the depth check ok and PASS defining variants: the
     reports' probable call (bin/pgx_outside_calls.py, kept for the
     summary), while the outside-call file stays empty and step 36's row
     stays indeterminate; no probable call with a failed depth check, a
     usable Cyrius call, a low-GQ or non-PASS defining variant, or a
     summary without the columns
 14. the whole subcommand on pypgx results archives, pypgx's own tables
     (rows of pypgx 0.26.0's gene and allele tables, tests/fixtures/pgx/
     summary/) and an indexed VCF: GSTT1 not callable on this reference
     (not FAILED), a gene without results FAILED, pypgx's CNV call for a
     BAM gene only, the evidence columns per gene, byte for byte what the
     Python it replaced wrote with the real pypgx and pysam; the tables
     found beside the pypgx on PATH; both files or neither

Run: python3 tests/test_pgx_outside_calls.py (OG names og, else
tests/og_binary.py finds the build)
"""
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "tests"))
from og_binary import og_binary  # noqa: E402

OG = og_binary()
FX = os.path.join(REPO, "tests", "fixtures", "pgx")
FAILS = []

HEADER = "Gene\tResult\tOutside_call\tReason\tEvidence\n"
NO_HLA = ("HLA-A\tnot typed\tno\tHLA typing (step 08) did not run\t-\n"
          "HLA-B\tnot typed\tno\tHLA typing (step 08) did not run\t-\n")
MULTIMAPPED = ("CYP2D6 copy number unreliable: reads are multi-mapped "
               "(was this BAM aligned to a reference with ALT contigs?)")


def run(**inputs):
    """(calls file, consensus file) og writes for these inputs."""
    with tempfile.TemporaryDirectory() as tmp:
        calls, cons = os.path.join(tmp, "calls.tsv"), os.path.join(tmp, "consensus.tsv")
        cmd = [OG, "pgx", "outside-calls", "--calls", calls, "--consensus", cons]
        for opt, name in inputs.items():
            cmd += ["--" + opt.replace("_", "-"), os.path.join(FX, name)]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            return f"exit {r.returncode}: {r.stderr.strip()}", ""
        return open(calls).read(), open(cons).read()


def case(desc, want_calls, want_consensus, **inputs):
    calls, cons = run(**inputs)
    for what, got, want in (("outside calls", calls, want_calls), ("consensus", cons, want_consensus)):
        if got == want:
            print(f"PASS {desc}: {what}")
        else:
            print(f"FAIL {desc}: {what}\n--- got\n{got}--- want\n{want}---")
            FAILS.append(f"{desc}: {what}")


# 1. Written first: a naive "take pypgx" passes pypgx's call on and fails here.
case("1 disagree", "",
     HEADER + NO_HLA +
     "CYP2D6\tindeterminate\tno\tpypgx and Cyrius disagree\tpypgx: *1/*4; Cyrius: *1/*2; depth check: ok\n",
     pypgx="pypgx_cyp2d6_1_4.tsv", cyrius="cyrius_1_2.tsv", depth_check="depth_check_ok.tsv")

case("2 one missing", "",
     HEADER + NO_HLA +
     "CYP2D6\tindeterminate\tno\tone caller only: a second caller must agree (Cyrius: setup.sh --cyrius)"
     "\tpypgx: *1/*2; Cyrius: not run; depth check: ok\n",
     pypgx="pypgx_cyp2d6_1_2.tsv", depth_check="depth_check_ok.tsv")

case("3 agree", "CYP2D6\t*1/*2\n",
     HEADER + NO_HLA +
     "CYP2D6\t*1/*2\tyes\tpypgx and Cyrius agree; the depth check passed"
     "\tpypgx: *1/*2; Cyrius: *1/*2; depth check: ok\n",
     pypgx="pypgx_cyp2d6_1_2.tsv", cyrius="cyrius_1_2.tsv", depth_check="depth_check_ok.tsv")

case("4 depth unreliable", "",
     HEADER + NO_HLA +
     f"CYP2D6\tindeterminate\tno\t{MULTIMAPPED}\tpypgx: *1/*2; Cyrius: *1/*2; depth check: unreliable\n",
     pypgx="pypgx_cyp2d6_1_2.tsv", cyrius="cyrius_1_2.tsv", depth_check="depth_check_multimapped.tsv")

case("4b no depth check", "",
     HEADER + NO_HLA +
     "CYP2D6\tindeterminate\tno\tthe CYP2D6 depth check did not run"
     "\tpypgx: *1/*2; Cyrius: *1/*2; depth check: not run\n",
     pypgx="pypgx_cyp2d6_1_2.tsv", cyrius="cyrius_1_2.tsv")

case("5 Cyrius made no call", "",
     HEADER + NO_HLA +
     "CYP2D6\tindeterminate\tno\tone caller only: a second caller must agree (Cyrius: setup.sh --cyrius)"
     "\tpypgx: *1/*2; Cyrius: no call (None, Not_assigned_to_haplotypes); depth check: ok\n",
     pypgx="pypgx_cyp2d6_1_2.tsv", cyrius="cyrius_none.tsv", depth_check="depth_check_ok.tsv")

case("6 HLA from T1K", "HLA-A\t*02:01/*24:02\nHLA-B\t*07:02/*44:02\n",
     HEADER +
     "HLA-A\t*02:01/*24:02\tyes\tT1K (step 08)\tHLA-A*02:01:01 (quality 60); HLA-A*24:02:01 (quality 60)\n"
     "HLA-B\t*07:02/*44:02\tyes\tT1K (step 08)\tHLA-B*07:02:01 (quality 60); HLA-B*44:02:01 (quality 48)\n"
     "CYP2D6\tindeterminate\tno\tno caller made a call\tpypgx: not run; Cyrius: not run; depth check: not run\n",
     hla="t1k_genotype.tsv")

case("7 HLA one allele, and an allele of quality 0", "HLA-A\t*02:01/*02:01\n",
     HEADER +
     "HLA-A\t*02:01/*02:01\tyes\tT1K (step 08): one allele, read as homozygous\tHLA-A*02:01:01 (quality 60)\n"
     "HLA-B\tindeterminate\tno\tT1K gives an allele quality 0 or below (T1K says to ignore it)"
     "\tHLA-B*07:02:01 (quality 60); HLA-B*44:02:01 (quality 0)\n"
     "CYP2D6\tindeterminate\tno\tno caller made a call\tpypgx: not run; Cyrius: not run; depth check: not run\n",
     hla="t1k_genotype_homozygous_lowqual.tsv")

# --- step 27 with the consensus table ------------------------------------------
PARSE = [OG, "pgx"]
REPORTS = os.path.join(REPO, "tests", "fixtures", "pharmcat")


def cpic(report, **inputs):
    """The recommendations text of og pgx cpic-report, with a consensus
    table from og pgx outside-calls when consensus_from is given."""
    with tempfile.TemporaryDirectory() as tmp:
        cmd = [*PARSE, "cpic-report", "--sample", "T", "--outdir", tmp,
               "--report", os.path.join(REPORTS, report)]
        if "pypgx" in inputs:
            cmd += ["--pypgx", os.path.join(FX, inputs["pypgx"])]
        if "consensus_from" in inputs:
            cons = os.path.join(tmp, "consensus.tsv")
            gen = [OG, "pgx", "outside-calls", "--calls", os.path.join(tmp, "calls.tsv"), "--consensus", cons]
            for opt, name in inputs["consensus_from"].items():
                gen += ["--" + opt.replace("_", "-"), os.path.join(FX, name)]
            subprocess.run(gen, check=True, capture_output=True)
            cmd += ["--consensus", cons]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            return f"exit {r.returncode}: {r.stderr}"
        return open(os.path.join(tmp, "T_cpic_recommendations.txt")).read()


def report_has(desc, text, needle, present=True):
    if (needle in text) == present:
        print(f"PASS {desc}")
    else:
        print(f"FAIL {desc}: {'missing' if present else 'unexpected'} {needle!r}\n{text}")
        FAILS.append(desc)


# 8. The 3.4.0 report of the fixture: PharmCAT has no CYP2D6 result.
plain = cpic("report-3.4.0.json", pypgx="pypgx_cyp2d6_1_4.tsv")
report_has("8 control: without the table, pypgx's call alone draws a warning", plain,
           "WARNING: PharmCAT has no result for CYP2D6, but pypgx")
held = cpic("report-3.4.0.json", pypgx="pypgx_cyp2d6_1_4.tsv",
            consensus_from={"pypgx": "pypgx_cyp2d6_1_4.tsv", "depth_check": "depth_check_ok.tsv"})
report_has("8 the report has the outside-call section", held, "Calls From Other Tools (outside calls, step 36):")
report_has("8 CYP2D6 is indeterminate and why", held,
           "CYP2D6   indeterminate            not passed to PharmCAT: " + "one caller only")
report_has("8 it names the drugs CYP2D6 affects", held, "No drug guidance is given for CYP2D6 here. Drugs affected by CYP2D6: ")
report_has("8 no pypgx-only warning for the held-back gene", held, "PharmCAT has no result for CYP2D6", False)

# 9. PharmCAT's example report has HLA-B from an outside call (and CYP2D6 and
# MT-RNR1, which this consensus does not pass on: case 10).
passed = cpic("pharmcat-docs-example.json", consensus_from={"hla": "t1k_genotype_docs_example.tsv"})
report_has("9 HLA-B is listed as passed from T1K", passed,
           "HLA-B    *15:02/*57:01            passed to PharmCAT from T1K (step 08); PharmCAT reports it as an outside call.")
report_has("9 the gene line marks PharmCAT's outside call", passed, "[outside call]")
report_has("9 the confirmed HLA-B keeps its drug guidance", passed, "  HLA-B -- ")
report_has("9 no stale warning for the confirmed HLA-B", passed, "outside call HLA-B", False)


def medications(text):
    """The Affected Medications section of a recommendations text."""
    return text.split("Affected Medications:", 1)[-1].split("Calls From Other Tools", 1)[0]


# 10. Outside calls the consensus does not pass on.
plain = cpic("pharmcat-docs-example.json")
report_has("10 control: without the table, the example's CYP2D6 has drug guidance", medications(plain), "  CYP2D6 -- ")
report_has("10 control: without the table, the example's HLA-B has drug guidance", medications(plain), "  HLA-B -- ")
other = cpic("pharmcat-docs-example.json", consensus_from={"hla": "t1k_genotype.tsv"})
report_has("10 HLA-B of another call gets no drug guidance", medications(other), "  HLA-B -- ", False)
report_has("10 the report says the HLA-B outside call is stale", other,
           "WARNING: PharmCAT's report has the outside call HLA-B *15:02/*57:01, which the consensus does not pass on")
report_has("10 the consensus HLA-B is not shown as PharmCAT's", other,
           "HLA-B    *07:02/*44:02            passed to PharmCAT from T1K (step 08); PharmCAT's report does not list it")
held = cpic("pharmcat-docs-example.json", consensus_from={"pypgx": "pypgx_cyp2d6_1_4.tsv", "depth_check": "depth_check_ok.tsv"})
report_has("10 a held-back CYP2D6 gets no drug guidance from an older report", medications(held), "  CYP2D6 -- ", False)
report_has("10 the report says the CYP2D6 outside call is stale", held,
           "WARNING: PharmCAT's report has the outside call CYP2D6 *1/*3, which the consensus does not pass on")
report_has("10 and says how to fix it", held, "Rerun step 07, then this step. No drug guidance is given for CYP2D6 here.")
with tempfile.TemporaryDirectory() as tmp:
    r = subprocess.run([*PARSE, "cpic-report", "--sample", "T", "--outdir", tmp,
                        "--report", os.path.join(REPORTS, "pharmcat-docs-example.json"),
                        "--consensus", os.path.join(tmp, "missing.tsv")], capture_output=True, text=True)
    gone = open(os.path.join(tmp, "T_cpic_recommendations.txt")).read() if r.returncode == 0 else f"exit {r.returncode}"
report_has("10 with the table missing, no outside call gets drug guidance", medications(gone), "  HLA-B -- ", False)
report_has("10 with the table missing, the report says to run step 36 again", gone, "Run step 36 again, then steps 07 and 27.")

# --- step 32's summary and the probable CYP2D6 call -------------------------------
sys.path.insert(0, os.path.join(REPO, "bin"))
import pgx_outside_calls as poc  # noqa: E402

SUMMARY_HEADER = ["Gene", "Diplotype", "Phenotype", "CNV_call", "Source", "Status",
                  "Defining_variants", "Min_GQ", "Weakest_variant", "Defining_filter", "Quality"]


def same(desc, got, want):
    if got == want:
        print(f"PASS {desc}")
    else:
        print(f"FAIL {desc}: got {got!r}, want {want!r}")
        FAILS.append(desc)


# 13.
ROW = {"Gene": "CYP2D6", "Diplotype": "*1/*2", "Phenotype": "Normal Metabolizer", "Defining_filter": "PASS",
       "Quality": "ok"}
same("13 pypgx alone, depth ok, PASS defining variants: probable",
     poc.probable_cyp2d6(ROW, None, "ok"),
     {"diplotype": "*1/*2", "phenotype": "Normal Metabolizer", "label": "probable (one caller: pypgx)"})
same("13 Cyrius ran and made no call: still one caller, probable",
     (poc.probable_cyp2d6(ROW, ("None", "Not_assigned_to_haplotypes"), "ok") or {}).get("label"), poc.PROBABLE)
same("13 red control: a failed depth check gives no probable call", poc.probable_cyp2d6(ROW, None, "unreliable"), None)
same("13 a usable Cyrius call: the two callers decide, no probable call",
     poc.probable_cyp2d6(ROW, ("*1/*2", "PASS"), "ok"), None)
same("13 a low-GQ defining variant: no probable call",
     poc.probable_cyp2d6(dict(ROW, Quality=poc.LOW_GQ), None, "ok"), None)
same("13 a defining variant not PASS: no probable call",
     poc.probable_cyp2d6(dict(ROW, Defining_filter="not PASS: chr22:42128945 C>T RefCall"), None, "ok"), None)
same("13 an allele pypgx's table lacks: no probable call",
     poc.probable_cyp2d6(dict(ROW, Defining_filter="not in pypgx's table: *99"), None, "ok"), None)
same("13 a defining variant with no GQ: no probable call",
     poc.probable_cyp2d6(dict(ROW, Quality="GQ not in the VCF: check the reads"), None, "ok"), None)
same("13 a summary without the columns: no probable call",
     poc.probable_cyp2d6({"Diplotype": "*1/*2", "Phenotype": "NM"}, None, "ok"), None)
with tempfile.TemporaryDirectory() as tmp:
    summ = os.path.join(tmp, "pypgx.tsv")
    with open(summ, "w") as f:
        f.write("\t".join(SUMMARY_HEADER) + "\n"
                "CYP2D6\t*1/*2\tNormal Metabolizer\tNormal\tBAM\tcalled\t1/1 in the VCF\t50"
                "\tchr22:42126611 C>G GQ 50 AD 20,22\tPASS\tok\n")
    calls, cons = os.path.join(tmp, "c.tsv"), os.path.join(tmp, "k.tsv")
    subprocess.run([OG, "pgx", "outside-calls", "--calls", calls, "--consensus", cons, "--pypgx", summ,
                    "--depth-check", os.path.join(FX, "depth_check_ok.tsv")], check=True, capture_output=True)
    same("13 the outside-call file stays empty for a probable call", open(calls).read(), "")
    same("13 step 36's row stays indeterminate, one caller only", open(cons).read().splitlines()[-1].split("\t")[:4],
         ["CYP2D6", "indeterminate", "no", "one caller only: a second caller must agree (Cyrius: setup.sh --cyrius)"])

# 14. The subcommand on pypgx's results archives and tables and an indexed VCF.
SUMMARY_FX = os.path.join(FX, "summary")
# sample.vcf.gz (bgzip and tabix of bcftools 1.24) holds these records.
SUMMARY_VCF = """##fileformat=VCFv4.2
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr22,length=50818468>
##FILTER=<ID=RefCall,Description="ref">
##FORMAT=<ID=GT,Number=1,Type=String,Description="g">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="q">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="a">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
chr10\t94761900\t.\tC\tT\t50\tPASS\t.\tGT:GQ:AD\t0/1:60:0,31
chr10\t94775360\t.\tAGGGGGGGG\tA\t50\tPASS\t.\tGT:GQ:AD\t0/1:40:10,10
chr10\t94775367\t.\tA\tG\t50\tPASS\t.\tGT:GQ:AD\t0/1:48:20,19
chr10\t94781859\t.\tG\tA,C\t50\tPASS\t.\tGT:GQ:AD\t0/1:12:26,3,.
chr22\t42126611\t.\tC\tG\t50\tRefCall\t.\tGT:GQ:AD\t0/0:50:20,22
"""
# What bin/pgx_outside_calls.py pypgx-summary wrote for these inputs in the
# pypgx 0.26.0 image, with the real pypgx and pysam (the parent of the
# commit that deleted it).
SUMMARY_WANT = [
    SUMMARY_HEADER,
    ["GSTM1", "*0/*0", "no pypgx phenotype rule", "WholeDel2", "BAM", "called", "none", "N/A", "N/A", "none", "N/A"],
    ["GSTT1", "N/A", "N/A", "N/A", "BAM", "not callable on this reference", "N/A", "N/A", "N/A", "N/A", "N/A"],
    ["CYP2D6", "*1/*2x2", "Normal Metabolizer", "Duplication", "BAM", "called", "1/2 in the VCF", "50",
     "chr22:42126611 C>G GQ 50 AD 20,22", "not PASS: chr22:42126611 C>G RefCall; not in the VCF: chr22:42127941 G>A",
     "ok"],
    ["CYP2B6", "FAILED", "N/A", "N/A", "BAM", "failed", "N/A", "N/A", "N/A", "N/A", "N/A"],
    ["CYP2C19", "*2/*17", "Intermediate Metabolizer", "N/A", "VCF", "called", "3/3 in the VCF", "12",
     "chr10:94781859 G>A GQ 12 AD 26,3,.", "PASS", "low-confidence: check the reads"],
    ["CYP2C9", "Indeterminate", "Indeterminate", "N/A", "VCF", "no call", "N/A", "N/A", "N/A", "N/A", "N/A"],
    ["NUDT15", "*1/*99", "N/A", "N/A", "VCF", "called", "unknown", "N/A", "N/A", "not in pypgx's table: *99", "N/A"],
]
SUMMARY_LOG = ("Genes called: 4/7\nWARNING: star alleles resting on a variant below GQ 20, check the reads: "
               "CYP2C19 *2/*17 (chr10:94781859 G>A GQ 12 AD 26,3,.)\n")


def results_zip(path, genotype, phenotype, cnv):
    """A pypgx results archive: metadata.txt and the sample table data.tsv
    in one folder, as pypgx's Archive.to_file writes it."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("abc123/metadata.txt", "Gene=X\nAssembly=GRCh38\nSemanticType=SampleTable[Results]\n")
        z.writestr("abc123/data.tsv", "\tVariantData\tCNV\tCandidates\tGenotype\tPhenotype\n"
                   f"S\tx\t{cnv}\t*1\t{genotype}\t{phenotype}\n")


def summary(tmp, *extra, env=None):
    outbase, out = os.path.join(tmp, "out"), os.path.join(tmp, "summary.tsv")
    for gene, g, p, c in (("CYP2C19", "*2/*17", "Intermediate Metabolizer", "Normal"),
                          ("GSTM1", "*0/*0", "Indeterminate", "WholeDel2"),
                          ("CYP2D6", "*1/*2x2", "Normal Metabolizer", "Duplication"),
                          ("CYP2C9", "Indeterminate", "Indeterminate", "NA"),
                          ("NUDT15", "*1/*99", "", "Normal")):
        results_zip(os.path.join(outbase, gene, "results.zip"), g, p, c)
    r = subprocess.run([OG, "pgx", "pypgx-summary", "--outbase", outbase,
                        "--bam-genes", "GSTM1 GSTT1 CYP2D6 CYP2B6", "--vcf-genes", "CYP2C19 CYP2C9 NUDT15",
                        "--vcf", os.path.join(SUMMARY_FX, "sample.vcf.gz"), "--not-callable", "GSTT1",
                        "--out", out, *extra], capture_output=True, text=True, env=env)
    rows = [line.rstrip("\n").split("\t") for line in open(out)] if os.path.exists(out) else []
    return r, rows


same("14 the fixture VCF holds the records the test names",
     gzip.decompress(open(os.path.join(SUMMARY_FX, "sample.vcf.gz"), "rb").read()).decode(), SUMMARY_VCF)
with tempfile.TemporaryDirectory() as tmp:
    r, rows = summary(tmp, "--pypgx-data", SUMMARY_FX)
    same("14 the subcommand exits 0, without a warning", (r.returncode, r.stderr), (0, ""))
    same("14 the summary, row by row, as the Python wrote it", rows, SUMMARY_WANT)
    same("14 the log names the low-confidence call", r.stdout, SUMMARY_LOG)
with tempfile.TemporaryDirectory() as tmp:
    # pypgx's tables found beside the pypgx on PATH, for the Python its #!
    # line names (the image also has a python3.1 link to python3.12).
    data = os.path.join(tmp, "prefix", "lib", "python3.12", "site-packages", "pypgx", "api", "data")
    os.makedirs(data)
    for name in ("gene-table.csv", "allele-table.csv"):
        shutil.copy(os.path.join(SUMMARY_FX, name), data)
    os.symlink("python3.12", os.path.join(tmp, "prefix", "lib", "python3.1"))
    os.makedirs(os.path.join(tmp, "prefix", "bin"))
    with open(os.path.join(tmp, "prefix", "bin", "pypgx"), "w") as f:
        f.write("#!/usr/local/bin/python3.12\n")
    os.chmod(os.path.join(tmp, "prefix", "bin", "pypgx"), 0o755)
    env = dict(os.environ, PATH=os.path.join(tmp, "prefix", "bin") + os.pathsep + os.environ["PATH"])
    r, rows = summary(tmp, env=env)
    same("14 pypgx's tables are found beside the pypgx on PATH", (r.returncode, rows), (0, SUMMARY_WANT))
with tempfile.TemporaryDirectory() as tmp:
    env = dict(os.environ, PATH="/usr/bin:/bin")
    r, rows = summary(tmp, env=env)
    same("14 without pypgx's tables the evidence says not looked up, and why",
         (r.returncode, rows[3][6:] if len(rows) > 3 else None, "WARNING: no pypgx tables" in r.stderr),
         (0, ["not looked up", "N/A", "N/A", "N/A", "N/A"], True))

# 15. Both files or neither: a consensus table that cannot be written leaves
# no outside-call file behind.
with tempfile.TemporaryDirectory() as tmp:
    calls = os.path.join(tmp, "calls.tsv")
    r = subprocess.run([OG, "pgx", "outside-calls", "--calls", calls,
                        "--consensus", os.path.join(tmp, "no-such-dir", "consensus.tsv"),
                        "--hla", os.path.join(FX, "t1k_genotype.tsv")], capture_output=True, text=True)
    same("15 the consensus cannot be written: exit 1, and no outside-call file",
         (r.returncode, os.path.exists(calls), sorted(os.listdir(tmp))), (1, False, []))

if FAILS:
    print(f"{len(FAILS)} check(s) failed")
    sys.exit(1)
print("all checks passed")
