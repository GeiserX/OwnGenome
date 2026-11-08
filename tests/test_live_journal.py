#!/usr/bin/env python3
"""The Journal layer of a person's results (bin/live_results.py) and the finding page's route
(bin/owngenome_live_results.py).

Checks:
  1. the level mapping, on a small made-up summary: CPSR's pathogenic row is to act on, a ClinVar hit
     outside the carrier list is worth a look, a carrier-list row is for family planning (two changed
     copies of a carrier gene, and a full-mutation repeat, are to act on), a CPIC gene
     with changed function is for the prescriber, a contamination warning is data quality; each card
     takes its highest finding's level, the quality chapter is data quality, a failed card has none;
     every finding id is a safe path part (live_results.FINDING_ID) and sits in its level's section;
  2. the lead paragraph: with a finding to act on it says how many go to a doctor and names each one
     (marked, so a child's page can leave the names out), a data-quality problem comes first, and with no finding it says none needs a doctor and points at
     sections 8 and 9, never an all-clear (live_results.FORBIDDEN); with no clinical screen run it says the
     clinical screens have no result yet;
  3. through the real server: /runs/<run>/finding/<id>/ sends a browser with no session to the login and
     back, answers the finding page
     with the run and the id filled in and the server's Content-Security-Policy for an id of the run's
     findings, and a 404 for an id that is not one, for a path-like id and for a run with no such page;
  4. /runs.json of a running run counts from the summary its page last read only while that read is
     current (same steps done, same summary.json), and says whether a clinical screen has a result.
Each check runs again on a planted mistake and must fail there: 1 with the CPIC rule's level swapped
to review, 2 with the act findings left out of the lead and with their names unmarked and with the
old lead on a run with no clinical screen, 3 with the route answering any id, 4 with a current read. A check
that cannot fail is no check.

Run: python3 tests/test_live_journal.py
"""
import contextlib
import copy
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import urllib.parse

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.join(REPO, "tests"))
import collect_summary  # noqa: E402
import live_results  # noqa: E402
import owngenome_live  # noqa: E402

SAMPLE = "DEMO-001"
PASSWORD = "synthetic test password 9"
FAILS = []


def check(desc, ok, detail=""):
    print(("[PASS] " if ok else "[FAIL] ") + desc + ("" if ok else f" -- {str(detail)[:800]}"))
    if not ok:
        FAILS.append(desc)


def control(desc, still_passes):
    check(f"control: {desc}", not still_passes, "the check passed on the planted mistake")


def ok(data):
    return {"state": "ok", "data": data}


# A made-up summary: invented genes and positions, one row per rule under test.
SUMMARY = {"sections": {
    "cpsr": ok({"plp": [{"gene": "GENEA", "variant": "c.100del", "final": "Pathogenic", "clinvar": "Pathogenic"}]}),
    "clinvar": ok({"hits": [{"gene": "GENEB", "chrom": "chr2", "pos": 1000, "ref": "C", "alt": "T",
                             "significance": "Pathogenic", "stars": 2, "genotype": "het", "af": 0.0001}]}),
    "carrier": ok({"carriers": [{"gene": "GENEC", "variant": "c.5A>G", "condition": "Invented condition",
                                 "significance": "Pathogenic", "stars": 2, "inheritance": "AR", "genotype": "het"},
                                {"gene": "GENED", "variant": "c.9C>T", "condition": "Another invented condition",
                                 "significance": "Pathogenic", "stars": 2, "inheritance": "AR", "genotype": "hom"}]}),
    "expansions": ok({"loci": [{"locus": "FMR1", "repeat_count": "/".join(("30", "250")), "status": "pre_mutation,full_mutation",
                                "ci": "29-31/240-260", "allele_sizing": "spanning/flanking"}]}),
    "cpic": ok({"genes": [{"gene": "CYP2C19", "diplotype": "*2/*2", "phenotype": "Poor Metabolizer",
                           "status": "non-normal"},
                          {"gene": "CYP2C9", "diplotype": "*1/*1", "phenotype": "Normal Metabolizer", "status": "normal"}]}),
    "sample_qc": ok({"contamination": "warn", "freemix": 0.05}),
    "coverage": ok({"mean_depth": 31.2}),
    "lpa": {"state": "failed", "note": "invented failure", "data": {}},
}}

WANT = {"cpsr": "act", "clinvar": "review", "carrier-genec": "carrier", "carrier-gened": "act", "str-fmr1": "act",
        "pgx": "medicine", "contamination": "trust"}


def level_problems(res):
    out = []
    found = res["findings"]
    for prefix, lv in WANT.items():
        got = [f["level"] for f in found if f["id"].startswith(prefix)]
        if got != [lv]:
            out.append(f"{prefix}: {got}, want [{lv}]")
    fmr1 = [f for f in found if f["id"] == "str-fmr1"]
    if fmr1 and "full-mutation" not in fmr1[0]["title"]:
        out.append(f"a full-mutation repeat titled {fmr1[0]['title']!r}")
    if any(f["id"].startswith("pgx-cyp2c9") for f in found):
        out.append("a gene with a normal phenotype became a finding")
    for f in found:
        if not live_results.FINDING_ID.match(f["id"]):
            out.append(f"unsafe id {f['id']!r}")
        if f["section"] != live_results.LEVEL_SECTION[f["level"]]:
            out.append(f"{f['id']} in section {f['section']}")
    cards = {c["id"]: c for ch in res["chapters"] for c in ch["sections"]}
    for cid, lv in (("cpsr", "act"), ("cpic", "medicine"), ("carrier-list", "act"), ("clinvar", "review"),
                    ("qc", "trust"), ("lpa", None)):
        if cards[cid].get("level") != lv:
            out.append(f"card {cid}: {cards[cid].get('level')}, want {lv}")
    return out


def lead_text(segs):
    return "".join(x[0] for x in segs)


def lead_problems(res, screened=True):
    out = []
    text = lead_text(res["lead"])
    acts = [f for f in res["findings"] if f["level"] == "act"]
    trust = [f for f in res["findings"] if f["level"] == "trust"]
    if acts:
        n = len(acts)
        said = re.search(r"\b(\w+) findings? should go to a doctor", text)
        words = {"one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6, "seven": 7, "eight": 8, "nine": 9}
        if not said or words.get(said.group(1).lower(), int(said.group(1)) if said.group(1).isdigit() else None) != n:
            out.append(f"no count of findings for a doctor: {text[:300]}")
        out += [f"lead does not name {f['title']!r}" for f in acts if f["title"] not in text]
    elif screened and ("needs a doctor" not in text or "sections 8 and 9" not in text):
        out.append(f"no-finding lead: {text[:300]}")
    elif not screened and ("needs a doctor" in text or "clinical screens have no result yet" not in text):
        out.append(f"no clinical screen ran, yet the lead reads: {text[:300]}")
    if trust and not text.startswith("Read the data-quality section first"):
        out.append("a data-quality problem is not first")
    if live_results.FORBIDDEN.search(text):
        out.append(f"all-clear wording: {text[:300]}")
    # A child's tile: the lead without the names of the findings to act on.
    closed = "".join(x[0] if x[2:] != ["act"] else "" for x in res["lead"])
    out += [f"the closed lead still names {f['title']!r}" for f in acts if f["title"] in closed]
    return out


def request(port, method, path, body=None, headers=None):
    import http.client
    c = http.client.HTTPConnection("127.0.0.1", port, timeout=30)
    h = {"Host": f"127.0.0.1:{port}"}
    h.update(headers or {})
    c.request(method, path, body=body, headers=h)
    r = c.getresponse()
    data = r.read()
    hdr = {k.lower(): v for k, v in r.getheaders()}
    c.close()
    return r.status, hdr, data


def serve(root, state, pwfile):
    app = owngenome_live.App()
    with contextlib.redirect_stderr(io.StringIO()):
        live, web, ingest = owngenome_live.start(["--root", root, "--state-dir", state, "--port", "0", "--ingest-port", "0",
                                                  "--password-file", pwfile, "--quiet"], app)
    threading.Thread(target=web.serve_forever, daemon=True).start()
    return live, web, ingest


def stop(srv):
    live, web, ingest = srv
    live.stopping = True
    web.shutdown()
    web.server_close()
    ingest.server_close()


def route_problems(port, fid):
    out = []
    path = f"/runs/{SAMPLE}/finding/{urllib.parse.quote(fid)}/"
    code, h, _ = request(port, "GET", path)
    if code != 303 or h.get("location") != "/login?next=" + urllib.parse.quote(path):
        out.append(f"finding page without a session: {code} {h.get('location')}")
    code, h, _ = request(port, "POST", "/login", urllib.parse.urlencode({"password": PASSWORD, "next": "/"}),
                         {"Content-Type": "application/x-www-form-urlencoded"})
    m = re.match(r"owngenome_live=([^;]+)", h.get("set-cookie", ""))
    if not m:
        return out + [f"login {code}"]
    ck = {"Cookie": f"owngenome_live={m.group(1)}"}
    code, h, body = request(port, "GET", path, headers=ck)
    page = body.decode("utf-8", "replace")
    if code != 200 or f'data-finding="{fid}"' not in page or f'data-run="{SAMPLE}"' not in page \
            or "{{" in page or "script-src 'self'" not in h.get("content-security-policy", ""):
        out.append(f"the finding page: {code} {h.get('content-security-policy')} {page[:200]}")
    for bad in ("not-a-finding-of-this-run", "..%2F..%2Fsummary.json", "UPPER"):
        code, _, _ = request(port, "GET", f"/runs/{SAMPLE}/finding/{bad}/", headers=ck)
        if code != 404:
            out.append(f"finding {bad}: {code}, want 404")
    code, _, _ = request(port, "GET", f"/runs/NO-SUCH-RUN/finding/{urllib.parse.quote(fid)}/", headers=ck)
    if code != 404:
        out.append(f"a run that does not exist: {code}")
    return out


def main():
    work = tempfile.mkdtemp(prefix="test_live_journal.")
    srvs = []
    try:
        # 1.
        res = live_results.build("MADE-UP", copy.deepcopy(SUMMARY))
        bad = level_problems(res)
        check("1. each rule gives its level, each card its highest finding's level, each id is safe", not bad, bad)
        rows = live_results.finding_rows

        def swapped(summary):
            return [dict(f, level="review") if f["id"].startswith("pgx") else f for f in rows(summary)]
        live_results.finding_rows = swapped
        try:
            control("the CPIC rule's level swapped to review",
                    not level_problems(live_results.build("MADE-UP", copy.deepcopy(SUMMARY))))
        finally:
            live_results.finding_rows = rows

        # 2.
        bad = lead_problems(res)
        check("2. the lead counts and names the findings for a doctor, data quality first", not bad, bad)
        calm = copy.deepcopy(SUMMARY)
        for k in ("cpsr", "sample_qc", "carrier", "cpic", "expansions"):
            calm["sections"].pop(k)
        calm["sections"]["clinvar"] = ok({"hits": []})
        res2 = live_results.build("MADE-UP", calm)
        bad = lead_problems(res2)
        check("2. with nothing to act on the lead says none needs a doctor and points at sections 8 and 9", not bad,
              (bad, lead_text(res2["lead"])))
        # Only the medicine and quality steps finished: a medicine finding, no clinical screen. The lead must
        # not say nothing needs a doctor.
        bare = {"sections": {k: copy.deepcopy(SUMMARY["sections"][k]) for k in ("cpic", "sample_qc", "coverage")
                             if k in SUMMARY["sections"]}}
        res3 = live_results.build("MADE-UP", bare)
        bad = lead_problems(res3, screened=False)
        check("2. with no clinical screen run the lead says so, never that nothing needs a doctor", not bad,
              (bad, lead_text(res3["lead"])))
        old_lead = [["No finding in what these steps assessed needs a doctor; ", True],
                    ["sections 8 and 9 say what they do not see.", False]]
        control("a lead that calls a run with no clinical screen clear of doctor findings",
                not lead_problems(dict(res3, lead=old_lead), screened=False))
        planted = dict(res, lead=live_results.lead([f for f in res["findings"] if f["level"] != "act"], res["levels"],
                                                   SUMMARY))
        control("the act findings left out of the lead", not lead_problems(planted))
        untagged = dict(res, lead=[x[:2] for x in res["lead"]])
        control("the names of the findings to act on not marked for a child's page", not lead_problems(untagged))

        # 3.
        root = os.path.join(work, "root")
        os.makedirs(root)
        subprocess.run([sys.executable, os.path.join(REPO, "tests", "demo", "make_demo_sample.py"), "--genome-dir", root,
                        "--records", "2000"], check=True, stdout=subprocess.DEVNULL)
        d = os.path.join(root, SAMPLE)
        summary = collect_summary.collect(SAMPLE, d)
        with open(os.path.join(d, "summary.json"), "w") as f:
            json.dump(summary, f)
        found = live_results.build(SAMPLE, summary)["findings"]
        check("3. the demo sample has a finding to open", bool(found), found)
        pwfile = os.path.join(work, "password")
        with open(pwfile, "w") as f:
            f.write(owngenome_live.hash_password(PASSWORD) + "\n")
        srv = serve(root, os.path.join(work, "state"), pwfile)
        srvs.append(srv)
        fid = found[0]["id"] if found else "none"
        bad = route_problems(srv[1].server_address[1], fid)
        check("3. /runs/<run>/finding/<id>/ behind the login, the page for a finding of the run, 404 for any other",
              not bad, bad)
        real = live_results.FINDING_ID
        build = live_results.build
        live_results.FINDING_ID = re.compile(r".*")
        live_results.build = lambda *a, **k: dict(build(*a, **k), findings=[{"id": x} for x in (
            fid, "not-a-finding-of-this-run", "../../summary.json", "UPPER")])
        try:
            control("the route answering any id", not route_problems(srv[1].server_address[1], fid))
        finally:
            live_results.FINDING_ID = real
            live_results.build = build
        # 4. the home page's overview of a running run: the summary its page last read, only while current.
        live = srv[0]
        run = live.runs()[SAMPLE]
        real_state = live.base_state

        def running(r, *a, **k):
            st = real_state(r, *a, **k)
            return dict(st, run=dict(st["run"], state="running"))
        live.base_state = running
        try:
            st = running(run)
            done = tuple(x["step"] for x in st["steps"] if x["state"] == "done")
            stamp = live.file_stamp(os.path.join(d, "summary.json"))
            extra = copy.deepcopy(summary)
            extra.setdefault("sections", {})["cpsr"] = SUMMARY["sections"]["cpsr"]
            want = live_results.overview(summary)["levels"].get("act")

            def act_count(key):
                live.headline_cache[SAMPLE] = (key, None, extra)
                live.live_overview.pop(SAMPLE, None)
                item = [x for x in live.runs_index()["runs"] if x["name"] == SAMPLE][0]
                return (item.get("levels") or {}).get("act")
            got = act_count((SAMPLE, done, None, (1, 1)))
            check("4. the home page drops the page's summary once summary.json has changed since that read",
                  got == want, (got, want))
            control("a page summary read before summary.json changed is still used",
                    act_count((SAMPLE, done, None, stamp)) == want)
            got = act_count((SAMPLE, done + ("99",), None, stamp))
            check("4. ...and once the steps done have changed", got == want, (got, want))
            o = live_results.overview(summary)
            bare = {"sections": {k: SUMMARY["sections"][k] for k in ("cpic", "coverage")}}
            check("4. runs.json says whether a clinical screen has a result, so the home page claims nothing without one",
                  live_results.overview(bare)["screened"] is False and o["screened"] is True
                  and "screened" in [x for x in live.runs_index()["runs"] if x["name"] == SAMPLE][0],
                  (live_results.overview(bare)["screened"], o["screened"]))
        finally:
            live.base_state = real_state
            live.headline_cache.pop(SAMPLE, None)
    finally:
        for s in srvs:
            stop(s)
        shutil.rmtree(work, ignore_errors=True)
    if FAILS:
        print(f"\n{len(FAILS)} check(s) failed")
        return 1
    print("\nAll checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
