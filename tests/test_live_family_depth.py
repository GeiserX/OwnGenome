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
  6. family.json passes the wording gate (live_results.FORBIDDEN);
  7. a value from a stale section is flagged stale in its own cell, and only that one;
  8. every blood group system step 41 typed gets a row (ABO and RhD first, then the rest in order);
  9. each child's trio status: assessed with both parents sequenced and the file here; not assessed
     with the reason when a parent is not sequenced or not in the joint call; a missing file and a
     file linking out of the folder told apart, the link's target never written;
 10. compound-het pairs come only from the child's step 23 candidates: assessed, not, or not known;
 11. a member's summary.json of a shape the report never writes (sections a list, CYP2D6 calls a
     string, a CPIC gene named by a list) costs that member's column only: it says "summary not
     readable" there, and every other member's column is as before.
Each check runs again on a planted mistake and must fail there: the segments read from the wrong
child (1), the pairs left in file order (2), a trio file planted (3), a member's summary dropped (4),
the inside check bypassed (5), "ok" planted in a row (6), the stale section made complete (7), the
systems list emptied (8), the link replaced by a plain file (9), a bad summary read
silently as empty (11). A check that cannot fail is no check.

Run: python3 tests/test_live_family_depth.py
"""
import copy
import json
import os
import shutil
import subprocess
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


def stale_summary(state):
    S = summary("*1/*1", [("DPYD", "Reference/Reference", "Normal Metabolizer")], "H1")
    S["sections"]["bloodtype"] = {"state": state, "data": {
        "abo": "A1", "rhd": "D+", "systems": [{"system": "ABO", "phenotype": "A1"}, {"system": "RhD", "phenotype": "D+"},
                                              {"system": "Kell", "phenotype": "K-k+"},
                                              {"system": "Duffy", "phenotype": "Fy(a+b-)"}]}}
    return S


def source_summary(key, state):
    """Every side-by-side source present and complete, with PharmCAT's CYP2D6 call, the Y haplogroup from
    y_haplogroup alone (the fallback), and one section, key, in the given state."""
    S = summary("*1/*1", [("DPYD", "Reference/Reference", "Normal Metabolizer")], "H1")
    S["cyp2d6"]["calls"]["PharmCAT"] = "*1/*1"
    S["sections"]["hla"] = {"state": "ok", "data": {"loci": [{"gene": "HLA-A", "alleles": ["A*01:01", "A*02:01"]}]}}
    S["sections"]["y_haplogroup"] = {"state": "ok", "data": {"haplogroup": "R-M269"}}
    S["sections"][key]["state"] = state
    return S


def stale_problems(sb):
    """Member A's blood groups come from a stale section; B's from a complete one; C has no summary."""
    out = []
    for r in sb["rows"]:
        want = [r["group"] == "Blood groups", False, False]
        if r.get("stale") != want:
            out.append(f"{r['label']}: stale {r.get('stale')} != {want}")
    return out or ([] if sb["rows"] else ["no rows"])


def blood_problems(sb):
    labels = [r["label"] for r in sb["rows"] if r["group"] == "Blood groups"]
    want = ["ABO", "RhD", "Kell", "Duffy"]
    vals = {r["label"]: r["values"][0] for r in sb["rows"] if r["group"] == "Blood groups"}
    out = [] if labels == want else [f"blood rows {labels} != {want}"]
    if vals.get("Kell") != "K-k+" or vals.get("Duffy") != "Fy(a+b-)":
        out.append(f"values {vals}")
    return out


def trio_problems(st, want):
    out = []
    for child, key, assessed, words in want:
        got = (st.get(child) or {}).get(key) or {}
        if got.get("assessed") is not assessed or words not in (got.get("why") or ""):
            out.append(f"{child} {key}: {got} (want assessed={assessed}, '{words}')")
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

        sums7 = {"A": stale_summary("stale"), "B": summary("*1/*1", [], "U5"), "C": None}
        sb7 = fam.side_by_side(["A", "B", "C"], sums7)
        bad = stale_problems(sb7)
        check("7. a stale section's values are flagged stale in their own cells, and no other value is", not bad, bad)
        control("the stale section made complete",
                not stale_problems(fam.side_by_side(["A", "B", "C"], dict(sums7, A=stale_summary("ok")))))

        bad = blood_problems(sb7)
        check("8. every blood group system typed gets a row: ABO, RhD, then Kell and Duffy in order", not bad, bad)
        emptied = stale_summary("stale")
        emptied["sections"]["bloodtype"]["data"]["systems"] = []
        control("the systems list emptied", not blood_problems(fam.side_by_side(["A", "B", "C"], dict(sums7, A=emptied))))

        flags = {"cpic": {"CYP2D6 (callers)", "DPYD"}, "hla": {"HLA-A"}, "haplogroup": {"Mitochondrial"},
                 "y_haplogroup": {"Y chromosome"}}
        bad = []
        for key in flags:
            sb = fam.side_by_side(["A"], {"A": source_summary(key, "stale")})
            got = {r["label"] for r in sb["rows"] if r["stale"][0]}
            if got != flags[key]:
                bad.append(f"{key} stale flags {sorted(got)} != {sorted(flags[key])}")
        check("7. each stale source flags its own cells only: CPIC (and PharmCAT's CYP2D6 call), HLA, the "
              "haplogroups (the Y fallback too)", not bad, bad)
        sb = fam.side_by_side(["A"], {"A": source_summary("cpic", "ok")})
        control("the CPIC section made complete", {r["label"] for r in sb["rows"] if r["stale"][0]} == flags["cpic"])

        short = stale_summary("ok")
        short["sections"]["bloodtype"]["data"]["systems"] = [x for x in short["sections"]["bloodtype"]["data"]["systems"]
                                                             if x["system"] != "Kell"]
        merged = lambda a, b_: [r["label"] for r in fam.side_by_side(["A", "B"], {"A": a, "B": b_})["rows"]  # noqa: E731
                                if r["group"] == "Blood groups"]
        got = merged(short, stale_summary("ok"))
        kidd = lambda s, drop: [x for x in s + [{"system": "Kidd", "phenotype": "Jk(a+b+)"}]  # noqa: E731
                                if x["system"] != drop]
        a_, b2 = stale_summary("ok"), stale_summary("ok")
        a_["sections"]["bloodtype"]["data"]["systems"] = kidd(a_["sections"]["bloodtype"]["data"]["systems"], "Duffy")
        b2["sections"]["bloodtype"]["data"]["systems"] = kidd(b2["sections"]["bloodtype"]["data"]["systems"], "Kell")
        got2 = merged(a_, b2)
        check("8. a system only a later member has keeps step 41's order, also when each lacks one",
              got == ["ABO", "RhD", "Kell", "Duffy"] and got2 == ["ABO", "RhD", "Kell", "Duffy", "Kidd"], (got, got2))
        saved, fam.STEP41_SYSTEMS = fam.STEP41_SYSTEMS, ()
        try:
            control("step 41's order taken away", merged(a_, b2) == ["ABO", "RhD", "Kell", "Duffy", "Kidd"])
        finally:
            fam.STEP41_SYSTEMS = saved
        sys.path.insert(0, os.path.join(REPO, "tests"))
        from og_binary import og_binary
        step41 = subprocess.run([og_binary(), "bloodtype", "systems"], capture_output=True, text=True, check=True)
        check("8. the page's order of the systems is step 41's (og bloodtype systems)",
              list(fam.STEP41_SYSTEMS) == step41.stdout.split(), (fam.STEP41_SYSTEMS, step41.stdout))
        control("a system moved in the page's order",
                list(fam.STEP41_SYSTEMS[1:2] + fam.STEP41_SYSTEMS[:1] + fam.STEP41_SYSTEMS[2:]) == step41.stdout.split())

        tri = os.path.join(work, "trio")
        shutil.copytree(DEMO, tri)
        with open(os.path.join(tri, "members.tsv"), "a") as f:
            f.write("DEMO-004\tfemale\tDEMO-001\tDEMO-005\t\tyes\tyes\nDEMO-005\tfemale\t\t\t\tyes\tno\n")
        plan_p = os.path.join(tri, "family_plan.tsv")
        with open(plan_p) as f:
            plan_text = f.read()
        with open(plan_p, "w") as f:
            f.write(plan_text.replace("DEMO-001,DEMO-002,DEMO-003\t", "DEMO-001,DEMO-002,DEMO-003,DEMO-004\t"))
        secret = os.path.join(work, "elsewhere-private-name.tsv")
        shutil.copy(os.path.join(DEMO, "compound_het.tsv"), secret)
        os.remove(os.path.join(tri, "compound_het.tsv"))
        os.symlink(secret, os.path.join(tri, "compound_het.tsv"))
        want = [("DEMO-003", "de_novo", True, "both parents"),
                ("DEMO-004", "de_novo", False, "DEMO-005 (the mother) is not"),
                ("DEMO-003", "compound_het", False, "link out of the family folder"),
                ("DEMO-004", "compound_het", False, "link out of the family folder")]
        tp = fam.family_page("trio", tri)
        bad = trio_problems(tp["trio_status"], want)
        if tp["trio_file_state"] != {"de_novo.tsv": "here", "compound_het.tsv": "rejected"}:
            bad.append(f"file state {tp['trio_file_state']}")
        if "elsewhere-private-name" in json.dumps(tp):
            bad.append("the link's target is in family.json")
        check("9. trio status per child: assessed, a parent not sequenced, a file linking out refused", not bad, bad)
        with open(plan_p, "w") as f:
            f.write(plan_text)
        out_of_call = trio_problems(fam.family_page("trio", tri)["trio_status"],
                                    [("DEMO-004", "de_novo", False, "DEMO-004 is not in the joint call")])
        check("9. a child left out of the joint call is not assessed, and the page says so", not out_of_call, out_of_call)
        os.remove(os.path.join(tri, "de_novo.tsv"))
        miss = trio_problems(fam.family_page("trio", tri)["trio_status"],
                             [("DEMO-003", "de_novo", False, "de_novo.tsv is not here")])
        check("9. a missing de_novo.tsv is not assessed, said apart from a refused one", not miss, miss)
        os.remove(os.path.join(tri, "compound_het.tsv"))
        shutil.copy(secret, os.path.join(tri, "compound_het.tsv"))
        control("the link replaced by a plain file", not trio_problems(fam.family_page("trio", tri)["trio_status"], want[2:]))
        shutil.copy(os.path.join(DEMO, "de_novo.tsv"), tri)
        with open(plan_p, "w") as f:
            f.write(plan_text.replace("DEMO-001,DEMO-002,DEMO-003\t", "DEMO-002,DEMO-003,DEMO-004\t"))
        st9 = fam.family_page("trio", tri)["trio_status"]
        bad = trio_problems(st9, [("DEMO-003", "de_novo", False, "DEMO-001 (the father) is not"),
                                  ("DEMO-003", "compound_het", True, "a parent is sequenced"),
                                  ("DEMO-004", "de_novo", False, "DEMO-001 (the father) is not"),
                                  ("DEMO-004", "compound_het", False, "need a parent sequenced and in the joint call")])
        check("9. a sequenced parent left out of the joint call: no de novo, and no pairs when it was the only one",
              not bad, bad)
        with open(plan_p, "w") as f:
            f.write(plan_text)
        control("the parent back in the joint call", not trio_problems(
            fam.family_page("trio", tri)["trio_status"], [("DEMO-003", "de_novo", False, "DEMO-001 (the father) is not")]))

        st10 = fam.with_candidates(copy.deepcopy(st9), {"DEMO-003": "absent", "DEMO-004": "here"})
        bad = trio_problems(st10, [("DEMO-003", "compound_het", False, "wrote no candidates for DEMO-003"),
                                   ("DEMO-004", "compound_het", False, "need a parent")])
        st10 = fam.with_candidates(copy.deepcopy(st9), {"DEMO-003": "here"})
        bad += trio_problems(st10, [("DEMO-003", "compound_het", True, "candidates (step 23) are here")])
        for g, words in (("newer", "newer than compound_het.tsv"), (None, "run is not on this site")):
            st10 = fam.with_candidates(copy.deepcopy(st9), {"DEMO-003": g})
            bad += trio_problems(st10, [("DEMO-003", "compound_het", None, words)])
        check("10. compound-het pairs by the child's step 23 candidates: here assessed, absent not, newer or no run "
              "not known", not bad, bad)
        control("the candidates not known read as absent",
                not trio_problems(fam.with_candidates(copy.deepcopy(st9), {"DEMO-003": None}),
                                  [("DEMO-003", "compound_het", False, "wrote no candidates")]))
        run3 = os.path.join(work, "root", "DEMO-003")
        os.makedirs(os.path.join(run3, "clinical"))
        vcf = os.path.join(run3, "clinical", "DEMO-003_clinical.vcf.gz")
        none_yet = fam.candidates_given(Live, run3, "DEMO-003")
        open(vcf, "w").close()
        os.utime(vcf, (1000, 1000))
        here, newer = fam.candidates_given(Live, run3, "DEMO-003", 2000), fam.candidates_given(Live, run3, "DEMO-003", 500)
        os.remove(vcf)
        os.symlink(outside, vcf)
        got = (none_yet, here, newer, fam.candidates_given(Live, run3, "DEMO-003"), fam.candidates_given(Live, None, "DEMO-003"))
        check("10. candidates: none in the run, here, newer than compound_het.tsv, a link out of the run, no run",
              got == ("absent", "here", "newer", "absent", None), got)
        control("the inside check bypassed for the candidates", fam.candidates_given(Open, run3, "DEMO-003") == "absent")

        def malformed_problems():
            out = []
            for bad in ({"sections": ["not", "a", "mapping"]}, {"cyp2d6": "CYP2D6 *1/*4"},
                        {"sections": {"cpic": {"state": "ok", "data": {"genes": [{"gene": ["CYP2D6"]}]}}}}):
                try:
                    got = fam.side_by_side(["A", "B", "C"], dict(sums, B=bad))
                except Exception as e:  # noqa: BLE001 -- the very thing under test
                    out.append(f"{bad}: raised {type(e).__name__}")
                    continue
                cols = {r["label"]: r["values"] for r in got["rows"]}
                want = {r["label"]: r["values"] for r in fam.side_by_side(["A", "B", "C"], dict(sums, B=None))["rows"]}
                for label, v in want.items():
                    if cols.get(label, [None] * 3)[0] != v[0] or cols.get(label, [None] * 3)[1] != "summary not readable":
                        out.append(f"{bad}: {label} {cols.get(label)}")
            return out
        bad = malformed_problems()
        check("11. a member's summary of a shape the report never writes costs that member's column only", not bad, bad)
        real = fam.side_values

        def quiet(summary):  # the planted mistake: a bad summary read as one with no rows, never said
            try:
                return real(summary)
            except (AttributeError, TypeError):
                return []
        fam.side_values = quiet
        try:
            control("a bad summary shown as '.' in every row, as if it had none of them", not malformed_problems())
        finally:
            fam.side_values = real
    finally:
        shutil.rmtree(work, ignore_errors=True)
    if FAILS:
        print(f"\n{len(FAILS)} check(s) failed")
        return 1
    print("\nAll checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
