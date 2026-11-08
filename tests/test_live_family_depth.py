#!/usr/bin/env python3
"""The family page's depth (bin/owngenome_live_family.py family.json), on the demo's invented trio
(tests/fixtures/live/demo/DEMO-FAMILY) and on small made-up summaries.

Checks:
  1. each child's stretches are its <child>_segments.tsv row for row, and the chromosome lengths
     the painting scales by are GRCh38's chr1 to chr22 and chrX;
  2. the compound-het pairs come in trans first, and every de novo row is there;
  3. a family without the trio files says they are not there (trio_files) and lists nothing;
  4. the members side by side: one row per gene or system any member's summary.json has (the CYP2D6
     callers and CPIC's CYP2D6 phenotype each a row of their own), one
     column per member, '.' where a member's summary lacks the row, "no summary yet" for a member
     without one; CPIC's phenotype quoted, never bare;
  5. a member's summary.json that links out of its run folder is never read;
  6. family.json passes the wording gate (live_results.FORBIDDEN).
Each check runs again on a planted mistake and must fail there: the segments read from the wrong
child (1), the pairs left in file order (2), a trio file planted (3), a member's summary dropped (4),
the inside check bypassed (5), "ok" planted in a row (6). A check that cannot fail is no check.

Run: python3 tests/test_live_family_depth.py
"""
import copy
import json
import os
import shutil
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import live_results  # noqa: E402
import owngenome_live  # noqa: E402
import owngenome_live_family as fam  # noqa: E402

DEMO = os.path.join(REPO, "tests", "fixtures", "live", "demo", "DEMO-FAMILY")
FAILS = []


def check(desc, ok, detail=""):
    print(("[PASS] " if ok else "[FAIL] ") + desc + ("" if ok else f" -- {str(detail)[:800]}"))
    if not ok:
        FAILS.append(desc)


def control(desc, still_passes):
    check(f"control: {desc}", not still_passes, "the check passed on the planted mistake")


def segment_problems(page, d):
    rows = fam.read_tsv(os.path.join(d, "DEMO-003_segments.tsv"))
    got = page["segments"].get("DEMO-003") or []
    out = []
    if [(r["side"], r["chrom"], r["start"], r["end"]) for r in got] != \
            [(r["side"], r["chrom"], r["start"], r["end"]) for r in rows] or not rows:
        out.append(f"{len(got)} stretches against {len(rows)} rows")
    if [c for c, _ in page["chrom_lengths"]] != [f"chr{i}" for i in range(1, 23)] + ["chrX"]:
        out.append(f"lengths {page['chrom_lengths'][:3]}")
    return out


def pair_problems(page):
    phases = [r["phase"] for r in page["compound_het"]]
    out = []
    if not phases or phases[0] != "in trans" or phases != sorted(phases, key=lambda p: fam.PHASE_ORDER.get(p, 1)):
        out.append(f"order {phases}")
    if len(page["de_novo"]) != 2:
        out.append(f"{len(page['de_novo'])} de novo rows")
    return out


def summary(cyp, cpic, hap, y=None):
    S = {"cpic": {"state": "ok", "data": {"genes": [{"gene": g, "diplotype": dip, "phenotype": ph} for g, dip, ph in cpic]}},
         "haplogroup": {"state": "ok", "data": {"haplogroup": hap}},
         "bloodtype": {"state": "ok", "data": {"abo": "A1", "rhd": "D+"}}}
    if y:
        S["y_tree"] = {"state": "ok", "data": {"haplogroup": y}}
    return {"cyp2d6": {"calls": {"pypgx": cyp, "Cyrius": cyp}, "agree": True}, "sections": S}


def side_problems(sb):
    by = {(r["group"], r["label"]): r["values"] for r in sb["rows"]}
    want = {("Medicines", "CYP2D6 (callers)"): ["pypgx *1/*4, Cyrius *1/*4 (callers agree)", "pypgx *1/*1, Cyrius *1/*1 (callers agree)",
                                      "no summary yet"],
            ("Medicines", "CYP2C19"): ['*1/*2 "Intermediate Metabolizer"', ".", "no summary yet"],
            ("Medicines", "CYP2D6"): ['*1/*4 "Intermediate Metabolizer"', ".", "no summary yet"],
            ("Medicines", "DPYD"): ['Reference/Reference "Normal Metabolizer"'] * 2 + ["no summary yet"],
            ("Haplogroups", "Y chromosome"): ["I-M253", ".", "no summary yet"],
            ("Blood groups", "ABO"): ["A1", "A1", "no summary yet"]}
    out = [f"{k}: {by.get(k)} != {v}" for k, v in want.items() if by.get(k) != v]
    if sb["members"] != ["A", "B", "C"]:
        out.append(f"members {sb['members']}")
    groups = [r["group"] for r in sb["rows"]]
    if groups != sorted(groups, key=["Medicines", "HLA", "Haplogroups", "Blood groups"].index):
        out.append(f"groups out of order {groups}")
    return out


def main():
    work = tempfile.mkdtemp(prefix="test_live_family_depth.")
    try:
        page = fam.family_page("DEMO-FAMILY", DEMO)
        bad = segment_problems(page, DEMO)
        check("1. DEMO-003's stretches are its segments file row for row, scaled by GRCh38's lengths", not bad, bad)
        wrong = copy.deepcopy(page)
        wrong["segments"]["DEMO-003"] = wrong["segments"]["DEMO-003"][::-1]
        control("the stretches in another order", not segment_problems(wrong, DEMO))

        bad = pair_problems(page)
        check("2. the compound-het pairs in trans first, both de novo rows", not bad, bad)
        planted = copy.deepcopy(page)
        planted["compound_het"] = fam.read_tsv(os.path.join(DEMO, "compound_het.tsv"))
        control("the pairs left in file order", not pair_problems(planted))

        bare = os.path.join(work, "bare")
        os.makedirs(bare)
        shutil.copy(os.path.join(DEMO, "members.tsv"), bare)
        b = fam.family_page("bare", bare)
        ok = lambda p: (p["trio_files"] == {"de_novo.tsv": False, "compound_het.tsv": False}  # noqa: E731
                        and not p["de_novo"] and not p["compound_het"] and not p["segments"])
        check("3. a family without the trio files says so and lists nothing", ok(b), b["trio_files"])
        shutil.copy(os.path.join(DEMO, "de_novo.tsv"), bare)
        control("a planted de_novo.tsv", ok(fam.family_page("bare", bare)))

        sums = {"A": summary("*1/*4", [("CYP2D6", "*1/*4", "Intermediate Metabolizer"),
                                       ("CYP2C19", "*1/*2", "Intermediate Metabolizer"),
                                       ("DPYD", "Reference/Reference", "Normal Metabolizer")], "H1", "I-M253"),
                "B": summary("*1/*1", [("DPYD", "Reference/Reference", "Normal Metabolizer")], "U5"),
                "C": None}
        sb = fam.side_by_side(["A", "B", "C"], sums)
        bad = side_problems(sb)
        check("4. side by side: a row per gene or system, '.' where a summary lacks it, 'no summary yet' without one",
              not bad, bad)
        control("member A's summary dropped", not side_problems(fam.side_by_side(["A", "B", "C"], dict(sums, A=None))))

        run = os.path.join(work, "root", "A")
        os.makedirs(run)
        outside = os.path.join(work, "outside.json")
        with open(outside, "w") as f:
            json.dump(sums["A"], f)
        os.symlink(outside, os.path.join(run, "summary.json"))

        class Live:
            root = os.path.realpath(os.path.join(work, "root"))
            inside = staticmethod(owngenome_live.Live.inside)
        check("5. a summary.json linking out of its run is not read", fam.member_summary(Live, run) is None)

        class Open(Live):
            inside = staticmethod(lambda p, top: True)
        control("the inside check bypassed", fam.member_summary(Open, run) is None)

        page["side_by_side"] = sb
        bad = live_results.wording_problems(page)
        check("6. family.json has no word that reads as an all-clear", not bad, bad)
        planted = copy.deepcopy(page)
        planted["side_by_side"]["rows"][0]["values"][0] = "ok"
        control("an ok in a side-by-side row", not live_results.wording_problems(planted))
    finally:
        shutil.rmtree(work, ignore_errors=True)
    if FAILS:
        print(f"\n{len(FAILS)} check(s) failed")
        return 1
    print("\nAll checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
