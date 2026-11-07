#!/usr/bin/env python3
"""A person's results and a family's page on the live site: bin/live_results.py, the results plugin
(bin/owngenome_live_results.py), /runs.json and the family plugin's family.json.

Checks:
  1. render_report.SECTION_NOTES names a caveat for every section id of CHAPTERS,
     and each name is a text constant of bin/render_report.py;
  2. results.json of DEMO-001 (tests/demo/make_demo_sample.py, collect_summary):
     the seven chapters with every section id once, and each card's caveat is
     its SECTION_NOTES constants, word for word, never empty;
  3. the wording gate: no string of results.json, of family.json
     (tests/fixtures/families_live family_a) or of the page assets
     (assets/live: the text of the pages, the string literals of the scripts)
     says a result is fine (live_results.FORBIDDEN), a quoted tool term and a
     caveat aside;
  4. ClinVar: with hits the card lists them in a table; with none it has no
     table and its empty text says ClinVar lists only part of what causes
     disease (never an all-clear);
  5. through the real server: /runs/<run>/results.json is behind the login
     (401 without a session), answers the cards with one, and writes one access
     log row per request; summary.json stays 404; /runs.json lists the run with
     its progress and attention count, behind the login too;
  6. APOE: the Traits card's hidden text, behind its warning, is each APOE row
     as the report prints it (trait, rsID, genotype, reading), never a dict;
  7. the adapters carry what the report says of a section, on small made-up
     sections: the ancient PCA plot (base64 in the summary) as an image, the
     Slivar phase lines, every Exomiser variant of a gene, a not-assessed
     secondary-findings list, a carrier list not made and an unparsed PharmCAT
     report as sentences (never zero counts), an unsized "normal" locus as not
     sized, AnnotSV's exon-edge check, a genotype rule and a deep intronic call.
Each check runs again on a planted mistake and must fail there: 1 with VNTYPER's
entry dropped, 2 with one card's caveat blanked, 3 with "ok" planted in a card
and "All clear" in a script, 4 with the hits emptied, 5 with the results route
registered without the login, 6 with a row printed as its Python dict, 7 with the Slivar adapter returning nothing. A check
that cannot fail is no check.

Run: python3 tests/test_live_results.py
"""
import contextlib
import base64
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
import owngenome_live_family  # noqa: E402
import render_report as rr  # noqa: E402

SAMPLE = "DEMO-001"
ASSETS = os.path.join(REPO, "assets", "live")
FAMILIES = os.path.join(REPO, "tests", "fixtures", "families_live")
PASSWORD = "synthetic test password 7"
FAILS = []


def check(desc, ok, detail=""):
    print(("[PASS] " if ok else "[FAIL] ") + desc + ("" if ok else f" -- {str(detail)[:800]}"))
    if not ok:
        FAILS.append(desc)


def control(desc, still_passes):
    check(f"control: {desc}", not still_passes, "the check passed on the planted mistake")


# --- 1. the caveat of every section ---------------------------------------------------------

def notes_problems(chapters, notes, module):
    ids = [sid for _, _, _, secs in chapters for sid, _, _ in secs]
    out = [f"{sid}: no caveat" for sid in ids if not notes.get(sid)]
    out += [f"{sid}: {n} is not a text constant" for sid, names in notes.items() for n in names
            if not isinstance(getattr(module, n, None), str) or not getattr(module, n).strip()]
    out += [f"{sid}: no such section" for sid in notes if sid not in ids]
    return out


# --- 2. one card per section, its caveat word for word --------------------------------------

def card_problems(res):
    out = []
    want = [sid for _, _, _, secs in rr.CHAPTERS for sid, _, _ in secs]
    got = [c["id"] for ch in res["chapters"] for c in ch["sections"]]
    if got != want:
        out.append(f"sections {got} != {want}")
    if [ch["key"] for ch in res["chapters"]] != [k for k, *_ in rr.CHAPTERS]:
        out.append("chapters out of the report's order")
    for ch in res["chapters"]:
        for c in ch["sections"]:
            exp = [getattr(rr, n) for n in rr.SECTION_NOTES.get(c["id"], ())]
            if not c["caveat"] or c["caveat"] != exp:
                out.append(f"{c['id']}: caveat {c['caveat']!r:.120}")
    return out


# --- 3. the wording gate --------------------------------------------------------------------

def asset_text(path):
    """What a page or script of assets/live puts on screen: the text of an HTML page, the string
    literals of a script (its comments and code left out)."""
    with open(path) as f:
        t = f.read()
    if path.endswith(".html"):
        t = re.sub(r"<(script|style)\b.*?</\1>", " ", t, flags=re.S | re.I)
        return [x.strip() for x in re.split(r"<[^>]*>", t) if x.strip()]
    t = re.sub(r"/\*.*?\*/", " ", t, flags=re.S)
    t = re.sub(r"(?m)^\s*//[^\n]*", " ", t)
    return [m.group(2) for m in re.finditer(r"(['\"])((?:\\.|(?!\1).)*)\1", t)]


def asset_problems(texts):
    out = []
    for name, items in texts.items():
        for s in items:
            if live_results.FORBIDDEN.search(re.sub(r'"[^"]*"', "", s)):
                out.append(f"{name}: {s[:120]}")
    return out


# --- 4. ClinVar with and without hits -------------------------------------------------------

def clinvar_card(summary):
    res = live_results.build(SAMPLE, summary)
    return next(c for ch in res["chapters"] for c in ch["sections"] if c["id"] == "clinvar")


def clinvar_problems(with_hits, without):
    out = []
    n = len(with_hits["table"]["rows"]) if with_hits.get("table") else 0
    if not n:
        out.append("hits but no table")
    if without.get("table"):
        out.append("no hit but a table")
    if "lists only part" not in (without.get("empty_text") or ""):
        out.append(f"empty text {without.get('empty_text')!r}")
    return out


# --- 6. APOE behind its warning -------------------------------------------------------------

APOE_ROW = {"rsid": "rs429358", "gene": "APOE", "trait": "APOE e4 marker", "genotype": "T/C",
            "status": "called", "reading": "one copy of the e4 allele", "ancestry": "", "note": "", "pmids": "1"}


def apoe_problems(hidden):
    out = []
    text = (hidden or {}).get("text") or ""
    if "{" in text or "'rsid'" in text:
        out.append(f"a dict in the text {text!r}")
    wants = [APOE_ROW["trait"], APOE_ROW["rsid"], APOE_ROW["genotype"], APOE_ROW["reading"]]
    for want in wants:
        if want not in text:
            out.append(f"{want!r} missing from {text!r}")
    return out


# --- 7. the adapters on made-up sections -----------------------------------------------------

def adapter_sections():
    sv = {"type": "DEL", "sv": "chr1:100-5000", "gt": "0/1", "length": 4900, "acmg_class": 5, "genes": ["GENE1"],
          "exon_edge": "yes", "exon_edge_detail": "GENE1 exon 2", "population_af": None, "population_source": "",
          "common": False, "support": None}
    return {
        "ancient": {"state": "ok", "data": {"nearest": [], "sites_used": 1, "panel_sites": 2,
                                            "plot_svg": base64.b64encode(b"<svg></svg>").decode()}},
        "slivar": {"state": "ok", "data": {"prioritized": 2, "compound_het_variants": 2, "compound_het_genes": 1,
                                           "compound_het_pairs": 1, "compound_het_phase": {"trans": 1},
                                           "compound_het_trans": [{"gene": "GENE2", "variants": "v1 v2", "reads": 3}]}},
        "exomiser": {"state": "ok", "data": {"genes": [{"rank": 1, "gene": "GENE3", "variants": [
            {"hgvs": "c.1A>G"}, {"hgvs": "c.2C>T"}]}]}},
        "clinical": {"state": "ok", "data": {"acmg_sf": {"unavailable": "no VEP output"}}},
        "carrier": {"state": "ok", "data": {"unavailable": "no ClinVar screen"}},
        "cpic": {"state": "ok", "data": {"parse_failed": True, "genes": [], "warnings": []}},
        "expansions": {"state": "ok", "data": {"stranger": True, "loci": [
            {"locus": "LOCUS1", "repeat_count": "20/30", "status": "normal", "sizing": "flanking_only"}]}},
        "annotsv": {"state": "ok", "data": {"count": 1, "rows": [sv]}},
        "traits": {"state": "ok", "data": {"rows": [], "rules": [
            {"trait": "Rule trait", "gene": "GENE4", "rsids": "rs1,rs2", "status": "called", "result": "AB",
             "phase": "read_backed", "reading": "rule reading"}]}},
        "clinvar": {"state": "ok", "data": {"hits": [], "deep_intronic": [
            {"gene": "GENE5", "hgvs": "c.100+500A>G", "chrom": "chr1", "pos": 1, "ref": "A", "alt": "G",
             "call": "het", "detail": "GQ 40"}]}},
    }


def adapter_problems(adapters):
    S, out = adapter_sections(), []
    def run(name, sec=None):
        try:
            return adapters[name](S, {"sections": S}) or {}
        except Exception as e:  # a crash is a problem, not a test error
            out.append(f"{name}: {type(e).__name__} {e}")
            return {}
    text = lambda got: json.dumps(got)  # noqa: E731
    if not (run("ancient").get("image") or {}).get("src", "").endswith(S["ancient"]["data"]["plot_svg"]):
        out.append("ancient: no plot image")
    if "in trans on 3 reads" not in text(run("slivar")):
        out.append("slivar: no phase lines")
    if len((run("exomiser").get("table") or {}).get("rows") or []) != 2:
        out.append("exomiser: not every variant")
    got = run("clinical_filter")
    if "not assessed" not in text(got) or "ClinVar P/LP in its genes" in text(got):
        out.append(f"clinical: not-assessed not said {got}")
    got = run("carrier")
    if "not made" not in text(got) or got.get("facts"):
        out.append(f"carrier: not-made not said {got}")
    got = run("cpic")
    if "could not be parsed" not in text(got) or got.get("facts"):
        out.append(f"cpic: parse failure not said {got}")
    rows = (run("expansions").get("table") or {}).get("rows") or [[]]
    if "not sized" not in rows[0]:
        out.append(f"expansions: unsized locus shown as {rows[0]}")
    if "on exon edges" not in text(run("structural")):
        out.append("structural: no exon-edge check")
    if "Rule trait" not in text(run("traits")):
        out.append("traits: no genotype rule")
    if "GENE5" not in text(run("clinvar")):
        out.append("clinvar: no deep intronic call")
    return out


# --- 5. the server --------------------------------------------------------------------------

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


def serve(root, state, pwfile, open_results=False):
    app = owngenome_live.App()
    if open_results:  # the planted mistake: the results route without the login
        route = app.route
        app.route = lambda m, p, fn, auth=True: route(m, p, fn, auth and "results" not in p)
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


def server_problems(port, state):
    out = []
    path = f"/runs/{SAMPLE}/results.json"
    code, _, _ = request(port, "GET", path)
    if code != 401:
        out.append(f"results.json without a session: {code}")
    code, _, _ = request(port, "GET", "/runs.json")
    if code != 401:
        out.append(f"runs.json without a session: {code}")
    code, h, _ = request(port, "POST", "/login", urllib.parse.urlencode({"password": PASSWORD, "next": "/"}),
                         {"Content-Type": "application/x-www-form-urlencoded"})
    m = re.match(r"owngenome_live=([^;]+)", h.get("set-cookie", ""))
    if not m:
        return out + [f"login {code}"]
    ck = {"Cookie": f"owngenome_live={m.group(1)}"}
    log = os.path.join(state, "access.log")
    before = open(log).read().count("\tresults\t") if os.path.isfile(log) else 0
    code, h, body = request(port, "GET", path, headers=ck)
    res = json.loads(body) if code == 200 else {}
    if code != 200 or len(res.get("chapters") or []) != 7 or res.get("sample") != SAMPLE:
        out.append(f"results.json with a session: {code} {str(res)[:200]}")
    after = open(log).read().count("\tresults\t") if os.path.isfile(log) else 0
    if after - before != 1:
        out.append(f"{after - before} access log rows for one results.json")
    code, _, _ = request(port, "GET", f"/runs/{SAMPLE}/summary.json", headers=ck)
    if code != 404:
        out.append(f"summary.json served: {code}")
    code, _, body = request(port, "GET", "/runs.json", headers=ck)
    idx = json.loads(body) if code == 200 else {}
    run = next((r for r in idx.get("runs") or [] if r["name"] == SAMPLE), None)
    if not run or run["url"] != f"/runs/{SAMPLE}/" or not run.get("results") or run.get("attention") != res.get("attention"):
        out.append(f"runs.json: {code} {run}")
    return out


def main():
    work = tempfile.mkdtemp(prefix="test_live_results.")
    srvs = []
    try:
        # 1.
        bad = notes_problems(rr.CHAPTERS, rr.SECTION_NOTES, rr)
        check(f"1. every section id has its caveat ({len(rr.SECTION_NOTES)} ids)", not bad, bad)
        planted = {k: v for k, v in rr.SECTION_NOTES.items() if k != "vntyper"}
        control("VNTYPER's caveat entry dropped", not notes_problems(rr.CHAPTERS, planted, rr))

        # 2.
        root = os.path.join(work, "root")
        os.makedirs(root)
        subprocess.run([sys.executable, os.path.join(REPO, "tests", "demo", "make_demo_sample.py"), "--genome-dir", root,
                        "--records", "2000"], check=True, stdout=subprocess.DEVNULL)
        d = os.path.join(root, SAMPLE)
        summary = collect_summary.collect(SAMPLE, d)
        with open(os.path.join(d, "summary.json"), "w") as f:
            json.dump(summary, f)
        res = live_results.build(SAMPLE, summary)
        bad = card_problems(res)
        check("2. results.json: the report's chapters, one card per section, each with its caveats word for word",
              not bad, bad)
        blank = copy.deepcopy(res)
        blank["chapters"][3]["sections"][0]["caveat"] = []
        control("one card's caveat blanked", not card_problems(blank))

        # 3.
        bad = live_results.wording_problems(res)
        check("3. no string of results.json reads as an all-clear", not bad, bad)
        fam = owngenome_live_family.family_page("family_a", os.path.join(FAMILIES, "family_a"))
        bad = live_results.wording_problems(fam)
        check("3. no string of family.json reads as an all-clear", not bad and fam.get("layout"), bad)
        planted = copy.deepcopy(res)
        planted["chapters"][0]["sections"][0]["headline"] = "Coverage: ok"
        control("'ok' planted in a card headline", not live_results.wording_problems(planted))
        texts = {}
        for p, _, fs in os.walk(ASSETS):
            for f in fs:
                if f.endswith((".html", ".js")):
                    texts[os.path.relpath(os.path.join(p, f), ASSETS)] = asset_text(os.path.join(p, f))
        bad = asset_problems(texts)
        check(f"3. no page or script text of assets/live reads as an all-clear ({len(texts)} files)",
              not bad and len(texts) >= 10, bad)
        planted = dict(texts, **{"person.js": texts["person.js"] + ["All clear: nothing to see"]})
        control("'All clear' planted in person.js", not asset_problems(planted))

        # 4.
        none = copy.deepcopy(summary)
        cv = none["sections"]["clinvar"]["data"]
        cv["hits"], cv["count"], cv["by_stars"] = [], 0, {}
        with_hits, without = clinvar_card(summary), clinvar_card(none)
        bad = clinvar_problems(with_hits, without)
        check("4. ClinVar: a table with hits; without, the sentence that ClinVar lists only part", not bad, bad)
        control("the demo's hits emptied on both sides", not clinvar_problems(without, without))

        # 6.
        S = {"traits": {"data": {"rows": [], "called": 0, "not_covered": 0, "apoe": [APOE_ROW]}}}
        hidden = live_results.a_traits(S, SAMPLE).get("hidden")
        bad = apoe_problems(hidden)
        check("6. APOE: the hidden text reads each row as the report does, never a dict", not bad, bad)
        control("the APOE row printed as its Python dict", not apoe_problems({"text": str(APOE_ROW)}))

        # 7.
        names = {"ancient": "a_ancient", "slivar": "a_slivar", "exomiser": "a_exomiser",
                 "clinical_filter": "a_clinical_filter", "carrier": "a_carrier", "cpic": "a_cpic",
                 "expansions": "a_expansions", "structural": "a_structural", "traits": "a_traits",
                 "clinvar": "a_clinvar"}
        adapters = {k: getattr(live_results, v) for k, v in names.items()}
        bad = adapter_problems(adapters)
        check("7. the adapters carry the report's plot, phase, variants, rules and not-assessed sentences", not bad, bad)
        control("the Slivar adapter returning nothing", not adapter_problems(dict(adapters, slivar=lambda S, s: {})))

        # 5.
        pwfile = os.path.join(work, "password")
        with open(pwfile, "w") as f:
            f.write(owngenome_live.hash_password(PASSWORD) + "\n")
        state = os.path.join(work, "state")
        srv = serve(root, state, pwfile)
        srvs.append(srv)
        bad = server_problems(srv[1].server_address[1], state)
        check("5. results.json and runs.json behind the login, one log row, summary.json 404", not bad, bad)
        state2 = os.path.join(work, "state-open")
        srv = serve(root, state2, pwfile, open_results=True)
        srvs.append(srv)
        control("the results route registered without the login",
                not server_problems(srv[1].server_address[1], state2))
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
