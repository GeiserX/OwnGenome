#!/usr/bin/env python3
# remnant: kept for the reports, the search index and the live site; deleted in lane R6b
"""collect_summary.py: what other Python still reads of the summary collector.

The collector itself is `og summary collect` (rust/crates/og-summary): it
reads one sample folder's step outputs into summary.json, the one file both
reports, the live site and the run delta read. collect() below runs it, for
bin/render_report.py (steps 24, HTML_REPORT and generate-report.sh) and for
the live site's fallback while a run has no summary.json yet. The rest are
the names other Python imports until its own port: the section keys, titles
and steps (render_report.py), the ClinVar review stars and the gene tables'
allele cells (the tests), the CSQ rule and the table readers
(results_index.py), and the HLA quality floor (render_report.py). `og dev
remnants` (tests/parity/remnants.tsv) checks that nothing else comes back.
Standard library only.
"""
import csv
import gzip
import json
import os
import shutil
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
# The repository's assets/ folder beside bin/: og summary collect reads the
# ACMG SF rules, the carrier list and the other gene tables there. The
# containers of step 24 and HTML_REPORT mount bin/ alone, so render_report.py
# --assets names the folder instead.
ASSETS_DIR = os.path.join(os.path.dirname(HERE), "assets")

# The summary's sections in report order: key, title, step. og summary
# collect writes the same keys, titles and steps (rust/crates/og-summary).
SECTIONS = [
    ("coverage", "Coverage (mosdepth)", "16b"),
    ("callability", "Gene callability (mosdepth)", "16b"),
    ("sex_check", "Sex check (indexcov)", "16"),
    ("sample_qc", "Sample identity and contamination (somalier, VerifyBamID2)", "33"),
    ("mosaic_sex_loss", "Mosaic loss of chromosome Y or X", "33"),
    ("read_count", "Reads kept by alignment", "02"),
    ("variants", "Variant calling", "03"),
    ("clinvar", "ClinVar screen", "06"),
    ("couple", "Couple carrier screen", "39"),
    ("carrier", "Carrier list", "06"),
    ("pharmcat", "PharmCAT", "07"),
    ("cpic", "CPIC drug-gene recommendations", "27"),
    ("pypgx", "pypgx", "32"),
    ("cyrius", "Cyrius (CYP2D6)", "21"),
    ("parascopy", "SMN1/SMN2 copy number (Parascopy, short reads)", "35"),
    ("paralogs", "Paralog genes (Paraphase, PacBio HiFi)", "35b"),
    ("hla", "HLA typing (T1K)", "08"),
    ("prs", "Polygenic risk scores", "25"),
    ("ancestry_fine", "Fine-grained ancestry (ADMIXTURE, RFMix)", "26b"),
    ("manta", "Structural variants (Manta)", "04"),
    ("delly", "Structural variants (Delly)", "19"),
    ("cnvpytor", "Copy number (CNVpytor)", "18"),
    ("sv_consensus", "SV consensus merge", "22"),
    ("annotsv", "Structural variants, class 4 and 5 (AnnotSV)", "05"),
    ("expansions", "Repeat expansions (ExpansionHunter)", "09"),
    ("telomere", "Telomere content (TelomereHunter)", "10"),
    ("roh", "Runs of homozygosity", "11"),
    ("haplogroup", "Mitochondrial haplogroup", "12"),
    ("y_haplogroup", "Y-chromosome haplogroup (Yleaf)", "37"),
    ("y_tree", "Y placement on the current public tree", "37"),
    ("y_str", "Y-STR markers (HipSTR)", "45"),
    ("consumer_export", "Raw-data export (23andMe v5 layout)", "40"),
    ("bloodtype", "Blood groups (predicted)", "41"),
    ("lpa", "Lp(a): LPA KIV-2 repeat count", "46"),
    ("vntyper", "MUC1 VNTR frameshifts (VNtyper)", "52"),
    ("locityper", "Locityper second opinion at hard loci", "53"),
    ("mito", "Mitochondrial variants", "20"),
    ("traits", "Traits panel", "42"),
    ("archaic", "Archaic ancestry (IBDmix)", "43"),
    ("ancient", "Ancient populations (AADR PCA)", "44"),
    ("viral", "Viral DNA in the unmapped reads", "47"),
    ("exomiser", "Symptom-driven ranking (Exomiser)", "48"),
    ("chip", "Clonal haematopoiesis (CHIP)", "49"),
    ("mei", "Mobile element insertions (MEGAnE)", "50"),
    ("mosaic", "Mosaic chromosomal alterations (MoChA)", "51"),
    ("cpsr", "Cancer predisposition (CPSR)", "17"),
    ("clinical", "Clinical variant filter", "23"),
    ("slivar", "Variant prioritization (slivar)", "31"),
]


def og():
    """The og binary: OG when set (the tests), else bin/og beside this file."""
    return os.environ.get("OG") or os.path.join(HERE, "og")


def collect(sample, sample_dir, declared_sex=None, out=None):
    """One sample folder's summary, as og summary collect writes it. With OUT
    the summary is also kept in that file; else it goes through a temporary
    file, removed after. og's messages go to stderr. RuntimeError when og
    writes no summary (the folder cannot be read, or the file cannot be
    written)."""
    tmp = None
    if out is None:
        # A private folder: og writes a sibling file and renames it over OUT.
        tmp = tempfile.mkdtemp(prefix="summary.")
        out = os.path.join(tmp, "summary.json")
    cmd = [og(), "summary", "collect", "--sample", sample, "--sample-dir", sample_dir,
           "--assets", ASSETS_DIR, "--out", out]
    if declared_sex:
        cmd += ["--declared-sex", declared_sex]
    try:
        try:
            # og reads bin/hla_alleles.tsv and bin/str_intermediate_ranges.tsv
            # beside this file, wherever the og binary sits.
            env = dict(os.environ, OG_BIN_TABLES=os.environ.get("OG_BIN_TABLES", HERE))
            rc = subprocess.run(cmd, stdout=subprocess.DEVNULL, env=env, check=False).returncode
        except OSError as e:
            raise RuntimeError(f"og summary collect did not start ({cmd[0]}): {e}") from e
        if rc != 0:
            raise RuntimeError(f"og summary collect wrote no summary of {sample} (exit {rc})")
        with open(out, encoding="utf-8") as f:
            return json.load(f)
    finally:
        if tmp:
            shutil.rmtree(tmp, ignore_errors=True)


# ClinVar review status -> stars. og clinvar hits (rust/crates/og-clinvar) holds the same table
# (tests/test_collect_summary.py checks that the two agree).
STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
}


# gnomAD allele frequency above which a P/LP hit is marked common.
COMMON_AF = 0.01


# Steps run-all.sh runs and records in logs/run_status.tsv: the rows of
# assets/steps.tsv with recorded=yes (tests/test_run_state.py keeps the two
# equal; the set stays written here because this file also runs where only
# bin/ is mounted). A section whose step is not here (variant calling, step
# 03, which run-all reuses on purpose) is never marked stale.
RUN_ALL_STEPS = {"03a", "03b", "03c", "03d", "04", "04a", "04b", "05", "06", "07", "08", "09", "09b", "10", "11",
                 "12", "13", "14", "15", "16", "16b", "17", "18", "19", "20", "21", "22", "23", "25", "26",
                 "26b", "27", "28", "29", "30", "31", "32", "33", "35", "35b", "37", "38", "39", "40", "41",
                 "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53"}


def open_text(path):
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", errors="replace")
    return open(path, errors="replace")


def info_map(field):
    out = {}
    for kv in field.split(";"):
        k, _, v = kv.partition("=")
        out[k] = v
    return out


def read_tsv(path):
    with open_text(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def hla_quality(q):
    """T1K's quality as a number; one that is not a number counts as -1, as in step 36."""
    try:
        return float(q)
    except ValueError:
        return -1.0


HLA_MAX_CANDIDATES = 3   # a quality-0 allele with more candidates is "no confident call"


# Pages that show more of a called branch, free and with no upload. The URL
# patterns were checked live on 2026-10-25 (docs/37-y-haplogroup.md,
# "Explore your branch further"): FamilyTreeDNA Discover answers its own
# HaplogroupSearch for R-M269 and H5a7, and YFull's tree pages answer 200 for
# a branch it has and 404 for one it does not.
EXPLORE_DOCS = "https://github.com/GeiserX/OwnGenome/blob/main/docs/"


def csq_rank(fmt, csq):
    """The rule of bin/csq_rank.awk on one CSQ value: (class, entries, freq).
    class is MANE (entries of the MANE Select transcript), else CANONICAL,
    else ANY (every entry); freq is rare, common, unknown (no value) or none
    (the CSQ has no frequency field): MAX_AF, else gnomADe_AF and gnomADg_AF,
    the highest over the entries."""
    ix = {k: i for i, k in enumerate(fmt)}
    entries = [v for v in (t.split("|") for t in csq.split(",") if t) if len(v) == len(fmt)]

    def get(v, k):
        return v[ix[k]] if k in ix else ""
    mane = [v for v in entries if get(v, "MANE_SELECT")]
    canon = [v for v in entries if get(v, "CANONICAL") == "YES"]
    cls, chosen = ("MANE", mane) if mane else ("CANONICAL", canon) if canon else ("ANY", entries)
    fields = ["MAX_AF"] if "MAX_AF" in ix else [f for f in ("gnomADe_AF", "gnomADg_AF") if f in ix]
    afs = []
    for v in entries:
        for f in fields:
            for a in get(v, f).split("&"):
                try:
                    afs.append(float(a))
                except ValueError:
                    pass
    freq = "none" if not fields else "unknown" if not afs else "rare" if max(afs) < COMMON_AF else "common"
    return cls, [{"impact": get(v, "IMPACT"), "gene": get(v, "SYMBOL"),
                  "consequence": get(v, "Consequence").replace("&", ",")} for v in chosen], freq
