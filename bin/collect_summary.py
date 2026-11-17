#!/usr/bin/env python3
# remnant: kept for the reports, the couple screen, the search index and the live site; deleted in lane R6b
"""collect_summary.py: what other Python still reads of the summary collector.

The collector itself is `og summary collect` (rust/crates/og-summary): it
reads one sample folder's step outputs into summary.json, the one file both
reports, the live site and the run delta read. collect() below runs it, for
bin/render_report.py (steps 24, HTML_REPORT and generate-report.sh) and for
the live site's fallback while a run has no summary.json yet. The rest are
the names other Python imports until its own port: the section keys, titles
and steps (render_report.py), ClinVar hit rows and the gene tables' allele
cells (couple_report.py), the CSQ rule and the table readers
(results_index.py), and the HLA quality floor (render_report.py). `og dev
remnants` (tests/parity/remnants.tsv) checks that nothing else comes back.
Standard library only.
"""
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import tempfile
import zlib

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


# The call floor of step 6 (scripts/06-clinvar-screen.sh writes a floor other
# than this default into the hits VCF as ##PGP_CALL_FLOOR=GQ=..,DP=..): a P/LP
# hit below either is a low-confidence call, listed but counted apart. A hits
# file without the line used this default, or was written before the floor.
CALL_FLOOR = {"GQ": 30, "DP": 10}


# gnomAD allele frequency above which a P/LP hit is marked common.
COMMON_AF = 0.01


# Steps run-all.sh runs and records in logs/run_status.tsv: the rows of
# assets/steps.tsv with recorded=yes (tests/test_run_state.py keeps the two
# equal; the set stays written here because this file also runs where only
# bin/ is mounted). A section whose step is not here (variant calling, step
# 03, which run-all reuses on purpose) is never marked stale.
RUN_ALL_STEPS = {"03a", "03b", "03c", "03d", "04", "04b", "05", "06", "07", "08", "09", "09b", "10", "11",
                 "12", "13", "14", "15", "16", "16b", "17", "18", "19", "20", "21", "22", "23", "25", "26",
                 "26b", "27", "28", "29", "30", "31", "32", "33", "35", "35b", "37", "38", "39", "40", "41",
                 "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53"}


SEVERITY_KINDS = ("hypomorph", "pseudodeficiency")


class Unreadable(Exception):
    pass


def open_text(path):
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", errors="replace")
    return open(path, errors="replace")


def vcf_records(path):
    """Yield the split data lines of a VCF (plain or gzip)."""
    try:
        with open_text(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                yield line.rstrip("\n").split("\t")
    except (OSError, EOFError, gzip.BadGzipFile, zlib.error) as e:
        raise Unreadable(f"{os.path.basename(path)}: {e}") from e


def info_map(field):
    out = {}
    for kv in field.split(";"):
        k, _, v = kv.partition("=")
        out[k] = v
    return out


def read_tsv(path):
    with open_text(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def clinvar_variants(v, sep):
    """'chr1:10:C:T|chr1:15:G:A' -> 'chr1:10 C>T<sep>chr1:15 G>A'."""
    out = []
    for x in v.split("|"):
        f = (x.split(":") + ["", "", "", ""])[:4]
        out.append(f"{f[0]}:{f[1]} {f[2]}>{f[3]}")
    return sep.join(out)


def clinvar_combined_note(info, gt=""):
    """Step 6's note on a hit that is a combined allele of two neighbouring
    records, or part of one, from the tags og clinvar join-adjacent
    writes; og clinvar hits words it the same way. "." when none.
    On an unphased het call (0/1) the genotype cannot say the two parts share
    a haplotype, so the note says they may not."""
    a = gt.split("/")
    phase = ("; unphased, so the parts may sit on different haplotypes"
             if len(a) == 2 and all(x.isdigit() for x in a) and a[0] != a[1] else "")
    if info.get("PGP_JOINED_FROM"):
        return "combined allele of " + clinvar_variants(info["PGP_JOINED_FROM"], " + ") + phase
    if info.get("PGP_PART_OF"):
        return "part of a combined allele: " + clinvar_variants(info["PGP_PART_OF"], ", ") + phase
    if info.get("PGP_ONE_ALLELE"):
        return ("possibly one allele with " + clinvar_variants(info["PGP_ONE_ALLELE"], ", ")
                + ": check ClinVar for a combined record")
    return "."


def low_confidence(dp, gq, floor):
    """True when the call's GQ or DP is known and below the floor."""
    return (gq is not None and gq < floor["GQ"]) or (dp is not None and dp < floor["DP"])


def half_call(gt):
    """True for a GT of one alt allele and one missing one: 1/., ./1, 1|."""
    return bool(re.fullmatch(r"[1-9][/|]\.|\.[/|][1-9]", gt or ""))


def clinvar_zygosity(gt):
    """A ClinVar hit's genotype as og clinvar hits words it. One alt
    allele and one missing one (1/., ./1) is the half of a multiallelic call
    step 06 split with the other alt written as missing (a 1/2 at a poly-T
    tract), or a half call: never het, which would say the other copy is ref."""
    if half_call(gt):
        return f"one copy, other copy not called ref ({gt})"
    return {"0/1": "het", "1/0": "het", "0|1": "het", "1|0": "het", "1/1": "hom", "1|1": "hom"}.get(gt, gt)


def clinvar_row(r, floor=None):
    info = info_map(r[7]) if len(r) > 7 else {}
    rev = info.get("CLNREVSTAT", "")
    genes = [g.split(":")[0] for g in info.get("GENEINFO", "").split("|") if g]
    gt = r[9].split(":")[0] if len(r) > 9 else ""
    zyg = clinvar_zygosity(gt)
    fmt = dict(zip(r[8].split(":"), r[9].split(":"))) if len(r) > 9 else {}
    dp, gq = (int(fmt[k]) if fmt.get(k, "").isdigit() else None for k in ("DP", "GQ"))
    return {
        "stars": STARS.get(rev, 0),
        "chrom": r[0], "pos": int(r[1]), "ref": r[3], "alt": r[4],
        "genotype": zyg,
        "gene": ",".join(genes) or ".",
        "significance": info.get("CLNSIG", "").replace("_", " ") or ".",
        "review_status": rev.replace("_", " ") or ".",
        "combined_allele": clinvar_combined_note(info, gt),
        "dp": dp, "gq": gq,
        "low_confidence": low_confidence(dp, gq, floor or CALL_FLOOR),
    }


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


def allele_key(locus):
    """'chr6:26092913:G:A' -> ('chr6', 26092913, 'G', 'A')."""
    c, pos, ref, alt = locus.split(":")
    return (c, int(pos), ref, alt)


def parse_alleles(text, kinds=False):
    """The variants or severity cell of assets/acmg_sf.tsv and
    assets/carrier_genes.tsv: '.' or entries separated by ';', each
    '[kind ]chrom:pos:ref:alt name'. {key: {'kind', 'name'}}."""
    out = {}
    for e in (text or ".").split(";"):
        e = e.strip()
        if e in ("", "."):
            continue
        parts = e.split(" ", 2 if kinds else 1)
        if kinds and len(parts) < 2:
            raise ValueError(f"severity entry needs a kind and a chrom:pos:ref:alt: {e}")
        kind = parts.pop(0) if kinds else None
        if kinds and kind not in SEVERITY_KINDS:
            raise ValueError(f"severity '{kind}' is not one of {', '.join(SEVERITY_KINDS)}: {e}")
        out[allele_key(parts[0])] = {"kind": kind, "name": parts[1] if len(parts) > 1 else parts[0]}
    return out
