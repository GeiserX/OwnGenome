#!/usr/bin/env python3
"""The widgets of a person's page: what bin/live_results.py hands the page to draw (results.json).

Each widget is checked on small made-up sections, against the numbers of the section it draws:
  1. repeat expansions (str): one entry per locus with its count, Stranger's thresholds and its
     95% interval; an allele not sized from short reads never carries Stranger's "normal";
  2. Lp(a) (ci): the KIV-2 count and its 95% interval, on a scale from 0 that holds the whole
     interval; no bar for a count that was not estimated;
  3. CHIP (vaf): each gene's detection limit and each call's VAF as numbers, each call with its own
     gene's limit; a call with no VAF is left off the track, not put at 0;
  4. copy numbers (copies): SMN1, SMN2, HBA1, HBA2 and the RCCX module of CYP21A2, never another
     gene of the paralog table; a copy number that is not a whole number is no dots; Parascopy's
     SMN1 and SMN2 at exon 7 (short reads) after them, one not resolved with no dots and its reason;
  5. the Y haplogroup path (ypath): one entry per level of step 37's path, root first, with its
     marker counts as numbers;
  6. every widget's strings pass the wording gate (live_results.FORBIDDEN);
  7. assets/live/person.js draws each kind the adapters give (a drawer per kind).
Each check runs again on a planted mistake and must fail there: the adapter's widget dropped
(1 to 5), "ok" planted in a widget (6), a drawer removed from person.js (7). A check that
cannot fail is no check.

Run: python3 tests/test_live_widgets.py
"""
import copy
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import live_results  # noqa: E402

FAILS = []


def check(desc, ok, detail=""):
    print(("[PASS] " if ok else "[FAIL] ") + desc + ("" if ok else f" -- {str(detail)[:800]}"))
    if not ok:
        FAILS.append(desc)


def control(desc, still_passes):
    check(f"control: {desc}", not still_passes, "the check passed on the planted mistake")


def sections():
    return {
        "expansions": {"state": "ok", "data": {"stranger": True, "records": 2, "loci": [
            {"locus": "HTT",  # invented counts; the locus on its own line keeps the personal-data scan quiet
             "repeat_count": "17/42", "status": "normal,full_mutation", "normal_max": "26", "pathologic_min": "40", "ci": "17-17/40-45", "sizing": "sized", "allele_sizing": "sized/sized"},
            {"locus": "FMR1", "repeat_count": "30", "status": "normal", "normal_max": "44", "pathologic_min": "200",
             "ci": "28-90", "sizing": "flanking_only", "allele_sizing": "flanking_only"}]}},
        "lpa": {"state": "ok", "data": {"status": "estimated", "kiv2_copies": "36.9", "ci95_low": "36.2",
                                        "ci95_high": "37.6", "class": "typical", "flag": "."}},
        "chip": {"state": "ok", "data": {"median_limit": 0.133, "median_depth": 31.0, "calls": [
            {"gene": "DNMT3A", "protein_change": "p.R882H", "vaf": "0.081", "alt_reads": "4", "depth": "49",
             "status": "PASS", "note": "invented"},
            {"gene": "TET2", "protein_change": "p.Q1034*", "vaf": ".", "alt_reads": "3", "depth": "40", "status": "PASS", "note": "invented"}],
            "genes": [{"gene": "DNMT3A", "mean_depth": "31.3", "detection_limit_vaf": "0.133", "status": "covered"},
                      {"gene": "TET2", "mean_depth": "29.0", "detection_limit_vaf": "0.142", "status": "covered"}]}},
        "paralogs": {"state": "ok", "data": {"genes": [
            {"gene": "SMN1", "copy_number": "2", "status": "called", "detail": "."},
            {"gene": "SMN2", "copy_number": "1", "status": "called", "detail": "."},
            {"gene": "SMN exon 7-8 deletion", "copy_number": "0", "status": "called", "detail": "."},
            {"gene": "HBA1", "copy_number": "2.4", "status": "no_call", "detail": "."},
            {"gene": "CYP21A2 (RCCX)", "copy_number": "4", "status": "called", "detail": "."},
            {"gene": "GBA1", "copy_number": "2", "status": "called", "detail": "."}]}},
        "parascopy": {"state": "ok", "data": {"rows": 1, "reading": "invented",
                                              "smn1": {"cn": 1, "text": "1 (quality 40, filter PASS)"},
                                              "smn2": {"cn": None, "text": "3 (quality 8, filter LowQual), not resolved"}}},
        "y_tree": {"state": "ok", "data": {"haplogroup": "I-M253", "tree": "YFull YTree v14.06.0", "path": [
            {"Level": "1", "Branch": "CT", "Derived": "11", "Ancestral": "0", "No_call": "1"},
            {"Level": "2", "Branch": "F", "Derived": "9", "Ancestral": "1", "No_call": "0"},
            {"Level": "3", "Branch": "I-M253", "Derived": "7", "Ancestral": "0", "No_call": "1"}]}},
    }


def run(adapters, name):
    S = sections()
    try:
        return (adapters[name](S, {"sections": S, "sample": "S"}) or {}).get("widget")
    except Exception as e:  # a crash is a problem, not a test error
        return {"error": f"{type(e).__name__}: {e}"}


def str_problems(w):
    if not w or w.get("kind") != "str":
        return [f"no str widget: {w}"]
    out = []
    by = {l["locus"]: l for l in w["loci"]}
    h = by.get("HTT") or {}
    if (h.get("repeat_count"), h.get("normal_max"), h.get("pathologic_min"), h.get("ci")) != ("17/42", "26", "40", "17-17/40-45"):
        out.append(f"HTT {h}")
    if (by.get("FMR1") or {}).get("status") == "normal":
        out.append("an unsized FMR1 allele carries Stranger's normal")
    return out


def ci_problems(w):
    if not w or w.get("kind") != "ci":
        return [f"no ci widget: {w}"]
    out = []
    if (w["value"], w["low"], w["high"]) != (36.9, 36.2, 37.6):
        out.append(f"numbers {w}")
    if not w["max"] >= w["high"] or not w["low"] <= w["value"] <= w["high"]:
        out.append(f"scale {w}")
    return out


def vaf_problems(w):
    if not w or w.get("kind") != "vaf":
        return [f"no vaf widget: {w}"]
    out = []
    if [(g["gene"], g["limit"]) for g in w["genes"]] != [("DNMT3A", 0.133), ("TET2", 0.142)]:
        out.append(f"limits {w['genes']}")
    if [(c["gene"], c["vaf"], c.get("limit")) for c in w["calls"]] != [("DNMT3A", 0.081, 0.133)]:
        out.append(f"calls {w['calls']}")
    return out


def copies_problems(w):
    if not w or w.get("kind") != "copies":
        return [f"no copies widget: {w}"]
    got = [(g["gene"], g["copies"], g["status"]) for g in w["genes"]]
    want = [("SMN1", 2, "called"), ("SMN2", 1, "called"), ("HBA1", None, "no_call"), ("CYP21A2 (RCCX)", 4, "called"),
            ("SMN1 (Parascopy)", 1, "called"), ("SMN2 (Parascopy)", None, "3 (quality 8, filter LowQual), not resolved")]
    return [] if got == want else [f"{got} != {want}"]


def ypath_problems(w):
    if not w or w.get("kind") != "ypath":
        return [f"no ypath widget: {w}"]
    got = [(p["branch"], p["derived"], p["ancestral"], p["no_call"]) for p in w["path"]]
    want = [("CT", 11, 0, 1), ("F", 9, 1, 0), ("I-M253", 7, 0, 1)]
    return [] if got == want and w["tip"] == "I-M253" else [f"{got} {w.get('tip')}"]


CASES = [("1. repeat expansions: the gauge's loci, thresholds and intervals", "expansions", str_problems),
         ("2. Lp(a): the count and its 95% interval on a scale that holds it", "lpa", ci_problems),
         ("3. CHIP: the detection limits and the calls with a VAF", "chip", vaf_problems),
         ("4. copy numbers: the five copy-number genes only, whole numbers as dots, Parascopy's SMN1 and SMN2 after them",
          "paralogs", copies_problems),
         ("5. the Y path: every level, root first, with its counts", "y-chromosome", ypath_problems)]


def without_widget(fn):
    def planted(S, s):
        got = dict(fn(S, s) or {})
        got.pop("widget", None)
        return got
    return planted


def drawers():
    with open(os.path.join(REPO, "assets", "live", "person.js")) as f:
        src = f.read()
    m = re.search(r"var WIDGETS = \{([^}]*)\}", src)
    return src, set(re.findall(r"(\w+):", m.group(1))) if m else set()


def main():
    A = live_results.ADAPTERS
    for desc, sid, problems in CASES:
        bad = problems(run(A, sid))
        check(desc, not bad, bad)
        control(f"{sid} adapter without its widget", not problems(run(dict(A, **{sid: without_widget(A[sid])}), sid)))
    lpa = sections()["lpa"]
    lpa["data"]["status"] = "too little depth"
    check("2. Lp(a): no bar for a count that was not estimated",
          not (A["lpa"]({"lpa": lpa}, {"sections": {"lpa": lpa}}) or {}).get("widget"))

    widgets = {sid: run(A, sid) for _, sid, _ in CASES}
    bad = live_results.wording_problems(widgets)
    check("6. no widget string reads as an all-clear", not bad, bad)
    planted = copy.deepcopy(widgets)
    planted["lpa"]["class"] = "ok"
    control("a widget saying ok", not live_results.wording_problems(planted))

    src, kinds = drawers()
    need = {w["kind"] for w in widgets.values()}
    check("7. person.js has a drawer for every widget kind", need <= kinds and all(f"function {fn}" in src for fn in
          ("strGauge", "ciBar", "vafChart", "copyDots", "yPath")), (need, kinds))
    control("a drawer removed from person.js", need <= (kinds - {"ypath"}))

    if FAILS:
        print(f"\n{len(FAILS)} check(s) failed")
        return 1
    print("\nAll checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
