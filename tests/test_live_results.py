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
     (assets/live: the text of the pages, the string literals of the scripts,
     map.js among them) says a result is fine (live_results.FORBIDDEN), a quoted
     tool term and a caveat aside; nor does any string of the genome map
     documents `og map build` writes for DEMO-001 and DEMO-003 (Figure 1 shows
     their sentences), a layer's state code aside (skipped without og or
     cargo, except in CI);
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
     sized, AnnotSV's exon-edge check, a genotype rule and a deep intronic call;
  8. the veil's tags on DEMO-003 (make_demo_sample.py --sample DEMO-003, with its
     planted secondary finding and class 5 deletion): a tag is the list of every
     closed result a place names; the ClinVar row and viewer link of the planted
     hit, its genome-chart tick and its review-star bar carry the secondary
     finding's id, the deletion's span its own, every search item the fixed-panel
     tag (a hit quotes the summary's rows, prose), and nothing else is tagged;
     DEMO-001, with no result to act on, has no tag; every known deep intronic
     site's line carries the fixed-panel tag, never the id of a closed result in
     its gene (the card prints each site whether hit or not, so one missing line
     would name the gene); prose (a line, a fact, a headline) is tagged by the
     fixed-panel tag whatever it names, unless its card is declared as printing
     no gene there (LINES_NO_GENE, FACTS_NO_GENE, HEADLINES_NO_GENE), so no word
     is ever read as a gene; a gene cell is read as its list of symbols; a
     ClinVar row closed by its gene alone closes its "See the reads at" link with
     it; every other card's row that names the gene (the son's symptom ranking
     ranks it) carries one of the ids and no other row does, the ranking's
     finding carries the secondary finding's id, and so does a segment of the
     lead that prints the ranking, and a segment that prints two closed results'
     findings carries both ids; the home page's view of the run (overview) lists
     no such finding among its first ones and its closed lead names none; a CHIP
     call, a copy-number dot and the carrier list's two-hit fact and line carry
     the id of the closed result in their gene (the CHIP chart's panel of genes
     never does), and results.json carries the closed results as veil.spans (the
     family page closes a child's places while one is closed); a closed
     structural variant over two genes closes the rows and the finding that name
     only its second gene, a place over two closed results carries both ids, and
     veil.panels lists the cards of a fixed panel (live_results.FIXED_PANELS).
Each check runs again on a planted mistake and must fail there: 1 with VNTYPER's
entry dropped, 2 with one card's caveat blanked, 3 with "ok" planted in a card
and "All clear" in a script, 4 with the hits emptied, 5 with the results route
registered without the login, 6 with a row printed as its Python dict, 7 with the Slivar adapter returning nothing,
8 with a veil rule that tags nothing (rows, lines, facts and drawings), one that tags every row, one that tags
only a structural variant's first gene, one that keeps only the first id of a place, one that reads prose word by
word, and each entry of the declared lists dropped. A check
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


def map_documents(root, people):
    """[(run, document)] of `og map build` for each (run, results) of PEOPLE, as the map plugin builds them;
    None without an og (tests/test_live_demo.py og_binary)."""
    import gzip
    from test_live_demo import og_binary
    og = og_binary()
    if not og:
        return None
    ref = os.path.join(root, "reference")
    if not os.path.isdir(ref):
        shutil.copytree(os.path.join(REPO, "tests", "fixtures", "map", "gencode-slice", "reference"), ref)
    out = []
    for run, results in people:
        body = {k: v for k, v in results.items() if not k.startswith("_")}
        p = subprocess.run([og, "map", "build", "--root", root, "--run", run, "--assets", os.path.join(REPO, "assets"),
                            "--cache-dir", os.path.join(root, ".map-cache")], input=json.dumps({"results": body}).encode(),
                           capture_output=True, check=True)
        out.append((run, json.loads(gzip.decompress(p.stdout))))
    return out


def map_wording(doc):
    """The strings of a map document that read as an all-clear; a layer's state is a code, not words."""
    return live_results.wording_problems(doc, allow=("state",))


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


def veil_problems(res, idx):
    """What is wrong with the veil's tags of DEMO-003's results.json and search index."""
    import make_demo_sample
    gene = make_demo_sample.PLANTED["gene"]
    c, p, _, _ = make_demo_sample.PLANTED["hit"]
    sv_c, a, b = make_demo_sample.PLANTED["sv"]
    act = [f["id"] for f in res["findings"] if f["level"] == "act"]
    sf = next((i for i in act if i.startswith("sf-")), None)
    sv = next((i for i in act if i.startswith("sv-")), None)
    if not sf or not sv:
        return [f"no planted secondary finding and deletion among the results to act on: {act}"]
    v = res.get("veil") or {}
    out = []
    card = next(x for ch in res["chapters"] for x in ch["sections"] if x["id"] == "clinvar")
    # The planted hit sits inside the planted deletion and names its gene: its places name both closed results.
    both = sorted([sf, sv])
    rows = (v.get("rows") or {}).get("clinvar") or []
    want = [both if r[0] == gene and r[1].startswith(f"{c}:{p} ") else None for r in card["table"]["rows"]]
    if [sorted(t) if t else None for t in rows] != want:
        out.append(f"ClinVar rows tagged {rows}, want {want}")
    view = (v.get("viewer") or {}).get("clinvar") or []
    if [sorted(t) if t else None for t in view] != [both if x == f"{c}:{p}" else None for x in card.get("viewer") or []]:
        out.append(f"ClinVar viewer links tagged {view}")
    marks = {k: {m: sorted(t) for m, t in x.items()} for k, x in (v.get("marks") or {}).items()}
    specs = {sp["id"]: sp for sp in rr.chart_specs(res["_summary"])}
    hit = [m["id"] for m in specs["genome"]["marks"] if m.get("kind") == "hit" and m["pos"] == p and m["chrom"] == c]
    span = [m["id"] for m in specs["genome"]["marks"] if m.get("kind") == "sv" and (m["chrom"], m["pos"], m["end"]) == (sv_c, a, b)]
    bar = [m["id"] for m in specs["clinvar-stars"]["marks"] if gene in (m.get("detail") or "")]
    want = {"genome": {**{h: both for h in hit}, **{s_: both for s_ in span}}, "clinvar-stars": {x: both for x in bar}}
    if not hit or not span or not bar or marks != want:
        out.append(f"chart marks tagged {marks}, want {want}")
    # Every other card's rows: a row that names the gene carries one of the two ids, a row that does not carries
    # none; the symptom ranking's row of the gene carries the secondary finding's, and so does its finding.
    for cid, rows_ in ((x["id"], (x.get("table") or {}).get("rows") or []) for ch in res["chapters"] for x in ch["sections"]):
        if cid == "clinvar" or not rows_:
            continue
        got = (v.get("rows") or {}).get(cid) or [None] * len(rows_)
        bad = [(r, t) for r, t in zip(rows_, got) if bool(set(t or ()) & {sf, sv}) != (gene in " ".join(map(str, r)))]
        if bad or len(got) != len(rows_):
            out.append(f"{cid} rows tagged wrong: {bad or got}")
    exo = (v.get("rows") or {}).get("exomiser") or []
    ranked = [f["id"] for f in res["findings"] if f["level"] != "act" and gene in f["title"]]
    # A fixed-panel card's own findings (a medicine gene) carry the panel tag, whatever the closed result.
    named = {k: t for k, t in (v.get("findings") or {}).items() if t != [live_results.PANEL]}
    if not any(sf in (t or ()) for t in exo) or not ranked or set(named) != set(ranked) \
            or any(sf not in named[i] for i in ranked):
        out.append(f"the symptom ranking of {gene}: rows tagged {exo}, findings tagged {v.get('findings')}, want {ranked} -> {sf}")
    # A search item quotes the summary's rows: prose, so every one carries the fixed-panel tag.
    items = [it for items in idx["terms"].values() for it in items]
    if not items or any(it.get("veil") != [live_results.PANEL] for it in items):
        out.append(f"search items tagged {sorted({str(it.get('veil')) for it in items})}, want every one {[live_results.PANEL]}")
    return out


def veil_lines_problems():
    """What is wrong with the veil's tags of the ClinVar card's lines: a known deep intronic site is printed
    whether hit or not, a fixed panel, so every such line carries the panel tag (live_results.PANEL), the one in
    a closed secondary finding's gene and the others alike: a line closed alone would name the gene."""
    S = {"clinvar": {"state": "ok", "data": {"hits": [], "deep_intronic": [
        {"gene": "GLA", "hgvs": "NM_000169.3:c.640-801G>A", "chrom": "chrX", "pos": 101399747, "ref": "C",
         "alt": "T", "call": "hom", "detail": "GQ 40"},
        {"gene": "GENE5", "hgvs": "c.100+500A>G", "chrom": "chr1", "pos": 1, "ref": "A", "alt": "G",
         "call": "het", "detail": "GQ 40"}]}}}
    lines = live_results.a_clinvar(S, {"sections": S}).get("lines") or []
    found = [{"id": "sf-gla-test", "level": "act", "gene": "GLA", "detail": "GLA"}]
    v = live_results.veil({"sections": S}, found, [{"sections": [{"id": "clinvar", "lines": lines}]}])
    panel = [live_results.PANEL]
    want = [panel if ln.startswith(live_results.FIXED_LINES["clinvar"]) else None for ln in lines]
    got = (v.get("lines") or {}).get("clinvar")
    out = [] if want.count(panel) == 2 and got == want else [f"lines {lines} tagged {got}, want {want}"]
    # A segment of the lead is tagged by the findings it prints (lead()'s names), beside the act names the page
    # closes itself; a segment printing two closed results' findings carries both ids.
    rank = {"id": "exomiser-gla", "level": "review", "card": "exomiser", "gene": "GLA", "title": "Symptom ranking: GLA ranks 1"}
    two = found + [{"id": "sf-tp53-test", "level": "act", "gene": "TP53", "detail": "TP53"},
                   {"id": "exomiser-tp53", "level": "review", "card": "exomiser", "gene": "TP53",
                    "title": "Symptom ranking: TP53 ranks 2"}]
    lead = [[": Symptom ranking: GLA ranks 1 and Symptom ranking: TP53 ranks 2. ", False],
            [": Secondary finding in GLA. ", False, "act"], ["Mean depth", False]]
    got = live_results.veil({"sections": S}, two + [rank], [], lead, [[rank, two[2]], [found[0]], []]).get("lead")
    if got != [["sf-gla-test", "sf-tp53-test"], None, None]:
        out.append(f"the lead's segments tagged {got}, want both ids on the segment that prints both rankings")
    # A finding known by its gene alone (no position): its ClinVar row is closed by the gene, and the row's
    # "See the reads at" link, which shows only the position, with it.
    S["clinvar"]["data"]["hits"] = [{"gene": "GLA", "chrom": "chrX", "pos": 101398000, "ref": "G", "alt": "A"},
                                    {"gene": "GENE5", "chrom": "chr1", "pos": 5, "ref": "A", "alt": "G"}]
    card = dict(live_results.a_clinvar(S, {"sections": S}), id="clinvar")
    gene_only = [{"id": "carrier-gla-test", "level": "act", "gene": "GLA"}]
    v = live_results.veil({"sections": S}, gene_only, [{"sections": [card]}])
    rows, view = (v.get("rows") or {}).get("clinvar"), (v.get("viewer") or {}).get("clinvar")
    if rows != [["carrier-gla-test"], None] or view != rows:
        out.append(f"a gene-only finding: ClinVar rows tagged {rows}, viewer links {view}, want both "
                   "['carrier-gla-test', None]")
    return out


def widget_problems():
    """What is wrong with the veil's tags of a card's facts and drawing: a CHIP call and a copy-number dot, and a
    fact, that name a closed result's gene carry its id; the CHIP chart's list of genes, the same panel for
    everyone, is never tagged; nothing else is."""
    found = [{"id": "cpsr-tp53-test", "level": "act", "gene": "TP53", "detail": "TP53 c.524G>A"},
             {"id": "two-hits-smn1-test", "level": "act", "gene": "SMN1"}]
    S = {"chip": {"state": "ok", "data": {"median_limit": 0.05, "median_depth": 72, "calls": [
        {"gene": "DNMT3A", "protein_change": "p.R882H", "vaf": 0.12, "alt_reads": 9, "depth": 75, "status": "pass", "note": "a clone"},
        {"gene": "TP53", "protein_change": "p.R175H", "vaf": 0.08, "alt_reads": 6, "depth": 75, "status": "pass", "note": "a clone"}],
        "genes": [{"gene": "DNMT3A", "detection_limit_vaf": 0.05, "mean_depth": 75, "status": "ok"},
                  {"gene": "TP53", "detection_limit_vaf": 0.06, "mean_depth": 70, "status": "ok"}]}},
         "paralogs": {"state": "ok", "data": {"genes": [
             {"gene": "SMN1", "copy_number": 1, "status": "called", "detail": "."},
             {"gene": "SMN2", "copy_number": 2, "status": "called", "detail": "."}]}}}
    summ = {"sections": S}
    chip = dict(live_results.a_chip(S, summ), id="chip")
    par = dict(live_results.a_paralogs(S, summ), id="paralogs")
    par["facts"] = [["SMN1 copies", "1"], ["Genes read", "2"]]
    v = live_results.veil(summ, found, [{"sections": [chip, par]}])
    out = []
    want = {"chip": [None, ["cpsr-tp53-test"]], "paralogs": [["two-hits-smn1-test"], None]}
    if v.get("widget") != want:
        out.append(f"drawings tagged {v.get('widget')}, want {want}")
    # A fixed-panel card's facts are prose no list declares: the panel tag, whatever they name.
    if (v.get("facts") or {}) != {c: [[live_results.PANEL], [live_results.PANEL]] for c in ("chip", "paralogs")}:
        out.append(f"facts tagged {v.get('facts')}")
    if v.get("spans") != [list(x) for x in live_results.closable(found)]:
        out.append(f"spans {v.get('spans')}")
    # The carrier list's two-hit fact and line are tagged by the two-hit row's gene field, never by their words: the
    # line prints the inheritance code AR, and a closed result in AR leaves both open.
    th = {"genes": [{"gene": "CFTR", "verdict": "possible", "condition": "Cystic fibrosis", "deletion": "DEL chr7:117480025-117580025",
                     "variant": "c.1521_1523del", "genotype": "het", "why": "a deletion and a variant"}],
          "deletions_read": ["AnnotSV"], "deletions_not_read": []}
    S2 = {"carrier": {"state": "ok", "data": {"genes_on_list": 2, "carriers": [], "two_hits": th}}}
    car = dict(live_results.a_carrier(S2, {"sections": S2}), id="carrier-list")
    for g, want in (("CFTR", ["two-hits-cftr"]), ("AR", None)):
        found2 = [{"id": "two-hits-cftr", "level": "act", "gene": "CFTR"}] if g == "CFTR" else \
            [{"id": "cpsr-ar", "level": "act", "gene": "AR", "detail": "AR c.1A>G"}]
        v2 = live_results.veil({"sections": S2}, found2, [{"sections": [car]}])
        facts, lines = (v2.get("facts") or {}).get("carrier-list") or [], (v2.get("lines") or {}).get("carrier-list") or []
        if (facts[2:3] or [None])[0] != want or (lines[:1] or [None])[0] != want or (lines[1:2] or [None])[0] != [live_results.PANEL]:
            out.append(f"a closed result in {g}: the two-hit fact tagged {facts}, lines {lines}, want {want} and the panel tag")
    return out


def sv_genes_problems():
    """What is wrong with the veil's tags of a closed structural variant over two genes: a ClinVar row and a carrier
    row that name only its second gene (the hit sits elsewhere, so only the gene can tie it) carry its id, and so
    does the carrier finding; veil.spans lists it once per gene; veil.panels lists the cards of a fixed panel."""
    S = {"annotsv": {"state": "ok", "data": {"rows": [
        {"sv": "chr17:43044295-43170245", "type": "deletion", "length": -125950, "acmg_class": "5",
         "genes": ["BRCA1", "NBR2"], "gt": "0/1", "support": 2, "support_callers": ["Manta", "Delly"]}]}},
         "clinvar": {"state": "ok", "data": {"hits": [
             {"gene": "NBR2", "chrom": "chr1", "pos": 5, "ref": "A", "alt": "G", "significance": "Pathogenic"},
             {"gene": "GENE5", "chrom": "chr1", "pos": 9, "ref": "A", "alt": "G", "significance": "Pathogenic"}]}},
         "carrier": {"state": "ok", "data": {"genes_on_list": 2, "carriers": [
             {"gene": "NBR2", "variant": "c.1A>G", "genotype": "het", "significance": "Pathogenic"}]}},
         "paralogs": {"state": "ok", "data": {"genes": [{"gene": "SMN1", "copy_number": 2, "status": "called"}]}}}
    res = live_results.build("T", {"sections": S})
    v = res.get("veil") or {}
    sv = [f["id"] for f in res["findings"] if f["card"] == "structural-variants"]
    car = [f["id"] for f in res["findings"] if f["card"] == "carrier-list"]
    out = []
    if len(sv) != 1 or len(car) != 1:
        return [f"the planted findings: {[f['id'] for f in res['findings']]}"]
    if (v.get("rows") or {}).get("clinvar") != [[sv[0]], None]:
        out.append(f"ClinVar rows tagged {(v.get('rows') or {}).get('clinvar')}, want the deletion's id on the NBR2 row")
    if (v.get("rows") or {}).get("carrier-list") != [[sv[0]]]:
        out.append(f"carrier rows tagged {(v.get('rows') or {}).get('carrier-list')}")
    if (v.get("findings") or {}).get(car[0]) != [sv[0]]:
        out.append(f"the NBR2 carrier finding tagged {(v.get('findings') or {}).get(car[0])}")
    if sorted(s[1] for s in v.get("spans") or [] if s[0] == sv[0]) != ["BRCA1", "NBR2"]:
        out.append(f"spans {v.get('spans')}")
    if v.get("panels") != ["paralogs"]:
        out.append(f"panels {v.get('panels')}, want the copy-number card, the only fixed panel here")
    # A place that names two closed results keeps both ids: a class 4 duplication over BRCA1 and BRCA2 with a closed
    # CPSR result in each, so opening one never brings back a place that still names the other.
    S["annotsv"]["data"]["rows"] = [{"sv": "chr5:1000000-1300000", "type": "DUP", "length": 300000, "acmg_class": "4",
                                     "genes": ["BRCA1", "BRCA2"], "gt": "0/1", "support": 1, "caller": "Manta",
                                     "support_callers": ["Manta"]}]
    S["cpsr"] = {"state": "ok", "data": {"plp": [{"gene": g, "variant": "c.68_69del", "final": "Pathogenic"}
                                                 for g in ("BRCA1", "BRCA2")]}}
    res = live_results.build("T", {"sections": S})
    v = res.get("veil") or {}
    acts = [f["id"] for f in res["findings"] if f["card"] == "cpsr"]
    dup = [f["id"] for f in res["findings"] if f["card"] == "structural-variants"]
    want = sorted(acts)
    row = ((v.get("rows") or {}).get("structural-variants") or [None])[0]
    if len(acts) != 2 or len(dup) != 1 or sorted((v.get("findings") or {}).get(dup[0]) or []) != want \
            or sorted(row or []) != want:
        out.append(f"a duplication over two closed results: finding {(v.get('findings') or {}).get(dup[0] if dup else '')}, "
                   f"row {row}, want {want}")
    return out


def gene_and_position_problems():
    """Chart genes without position overlap, and a gene-free call with position overlap, need independent tags."""
    summary = {"sections": {
        "roh": {"state": "ok", "data": {}},
        "annotsv": {"state": "ok", "data": {"rows": [
            {"sv": "chr17:43044295-43170245", "type": "DEL", "length": -125950, "acmg_class": "5",
             "genes": ["BRCA1", "NBR2"], "gt": "0/1", "support": 1, "caller": "Manta"},
            {"sv": "chr5:1000000-1300000", "type": "DUP", "length": 300000, "acmg_class": "4",
             "genes": ["NKD2", "NBR2"], "gt": "0/1", "support": 1, "caller": "Manta"}]}},
        "clinvar": {"state": "ok", "data": {"hits": [
            {"gene": "NBR2", "chrom": "chr1", "pos": 5, "ref": "A", "alt": "G", "significance": "Pathogenic"}]}},
        "mosaic": {"state": "ok", "data": {"status": "assessed", "calls": [
            {"chrom": "chr17", "start": 43050000, "end": 43100000, "length_mb": 0.05,
             "type": "Loss", "arms": "17q", "cell_fraction": 0.2}]}}}}
    res = live_results.build("T", summary)
    v, out = res["veil"], []
    closed = "sv-chr17-43044295-43170245"
    for mark in ("hit-0", "sv-1"):
        got = v.get("marks", {}).get("genome", {}).get(mark)
        if got != [closed]:
            out.append(f"gene-only chart {mark}: {got}, want {[closed]}")
    if v.get("findings", {}).get("mca-chr17-43050000") != [closed]:
        out.append("coordinate-only mosaic finding not tagged with the deletion")
    if v.get("rows", {}).get("mosaic-mca") != [[closed]]:
        out.append("coordinate-only mosaic table row not tagged with the deletion")
    card = next(c for ch in res["chapters"] for c in ch["sections"] if c["id"] == "mosaic-mca")
    lines = [(ln, tag) for ln, tag in zip(card["lines"], v.get("lines", {}).get("mosaic-mca", []))
             if "chr17:43050000-43100000" in ln]
    if not lines or any(tag != [closed] for _, tag in lines):
        out.append(f"coordinate-only mosaic lines: {lines}")
    return out


# The veil's tables, pinned here so that an entry dropped from one turns this file red (each is checked by what
# veil() does with it, not by reading the table back): the fixed-panel cards, the fixed lines of a card (a sample
# line or row naming a gene, its card, and whether it is a row), the cards that name no gene, the cards whose lines,
# facts or headline name no gene, and prose that names a gene or a word spelled like one.
WANT_PANELS = ("mei", "expansions", "vntyper", "lpa", "chip", "paralogs", "cpic", "cyp2d6", "hla", "locityper",
               "traits-panel", "bloodtype")
WANT_LINES = (("clinvar", "Known deep intronic site (always reported): GLA c.640-801G>A chrX:1 C>T: hom", False),
              ("callability", ["GSTT1", "not on this reference", ".", "ALT-only"], True),
              ("callability", ["NBPF1", "not read", "69%", "segdup, collapsed"], True),
              ("roh", "Imprinted 6q24 (PLAGL1, HYMAI): from the mother: -; from the father: Transient neonatal diabetes "
                      "mellitus (6q24-related) (PMID 20301706)", False))
WANT_NO_GENE = (("prs", ["Asthma", "PGS002727", "68.3 (EUR)", "2605/2872", "HR 1.57 (EUR, PPM014751)"], "HR"),)
WANT_LINES_NO_GENE = ("qc", "clinvar", "mosaic", "mosaic-mca", "roh", "telomere", "consumer-export")
WANT_FACTS_NO_GENE = ("qc", "callability", "variants", "clinvar", "clinical-filter", "cpsr", "slivar", "exomiser",
                      "structural-variants", "mito", "mosaic", "mosaic-mca", "carrier-list", "couple", "pharmcat",
                      "ancestry-fine", "mt-haplogroup", "y-chromosome", "archaic", "ancient", "roh", "telomere",
                      "viral", "consumer-export")
WANT_FACT_PROSE = (("clinical-filter", "High-impact rare variants in its genes"),)
WANT_HEADLINES_NO_GENE = ("clinvar", "carrier", "roh", "pypgx", "annotsv", "viral", "archaic", "lpa", "cpsr", "mito")
# Prose, each with the gene it names or the gene its word is spelled like: the tag must be the panel tag either way.
# The carrier list's summary of the two-hit callers (it names SMN1, HBA1 and HBA2 as the callers it read) is prose of
# a card with no gene-free lines: the panel tag, as every other line of that card but a two-hit gene's own.
WANT_PROSE_LINES = (("carrier-list", "Two hits across callers (a deletion plus a variant on the other copy): none found; "
                                     "deletions read from Parascopy's SMN1 count (step 35)"),)
WANT_PROSE = (("KCNH2", "A change in KCNH2."), ("C9ORF72", "C9orf72 repeat"), ("HLA-A", "HLA-A*02:01"),
              ("SET", "the same set of genes"), ("MB", "0.0 MB in runs of 5 MB or more"), ("AR", "CFTR (Cystic fibrosis; AR)"),
              ("HP", "Seizure (HP:0001250)"), ("GC", "reads per million with 48-52% GC"))
WANT_GENE_CELLS = (("NKD2, SLC12A7 and 2 more", ["NKD2", "SLC12A7"]), ("SMN1 (Parascopy)", ["SMN1"]),
                   ("HLA-A", ["HLA-A"]), ("BRCA1, BRCA2", ["BRCA1", "BRCA2"]), (".", []))


def gene_of(text):
    """The first upper-case word of TEXT (a list is a table row) that looks like a gene symbol."""
    t = " ".join(map(str, text)) if isinstance(text, list) else text
    return next(w for w in re.findall(r"[A-Z][A-Z0-9]{2,}", t.split(":", 1)[-1] if not isinstance(text, list) else t)
                if w not in ("ALT", "SMN"))


def table_problems():
    """What is wrong with veil()'s tables (WANT_*): a fixed-panel card missing from veil.panels, a fixed line or row
    tagged by the gene it names instead of the panel tag, a card that names no gene tagged, a line, fact or headline
    of a card declared gene-free there tagged by the panel tag (or one of a card not declared not), prose read word by
    word, or a gene cell not read as its list of symbols."""
    out = []
    panel = [live_results.PANEL]
    found = lambda g: [{"id": "sf-" + g.lower(), "level": "act", "gene": g, "detail": g}]
    chapters = [{"sections": [{"id": c, "state": "ok"} for c in WANT_PANELS]}]
    got = live_results.veil({"sections": {}}, found("BRCA2"), chapters).get("panels") or []
    out += [f"card {c} is not in veil.panels" for c in WANT_PANELS if c not in got]
    for card, item, is_row in WANT_LINES:
        g = gene_of(item)
        sec = {"id": card, "table": {"head": ["a", "b", "c", "d"], "rows": [item]}} if is_row else {"id": card, "lines": [item]}
        v = live_results.veil({"sections": {}}, found(g), [{"sections": [sec]}])
        tag = ((v.get("rows" if is_row else "lines") or {}).get(card) or [None])[0]
        if tag != panel:
            out.append(f"{card}: {str(item)[:60]} tagged {tag} with a closed result in {g}, want the panel tag")
    for card, line in WANT_PROSE_LINES:
        v = live_results.veil({"sections": {}}, found("SMN1"), [{"sections": [{"id": card, "lines": [line]}]}])
        if ((v.get("lines") or {}).get(card) or [None])[0] != panel:
            out.append(f"{card}: {line[:60]} tagged {(v.get('lines') or {}).get(card)}, want the panel tag")
    for card, row, g in WANT_NO_GENE:
        sec = {"id": card, "table": {"head": ["a"] * len(row), "rows": [row]}, "lines": [" ".join(row)], "facts": [row[:2]]}
        v = live_results.veil({"sections": {}}, found(g), [{"sections": [sec]}])
        if any((v.get(k) or {}).get(card) for k in ("rows", "lines", "facts")):
            out.append(f"{card} names no gene, yet a closed result in {g} tags it: {v}")
    # A line or fact of a card declared gene-free is untagged even when it spells the closed gene; one of a card not
    # declared carries the panel tag whatever it says. A headline likewise, by the summary key it comes from.
    for card in sorted(set(WANT_LINES_NO_GENE) | {"cpsr", "slivar", "archaic"}):
        v = live_results.veil({"sections": {}}, found("BRCA2"), [{"sections": [{"id": card, "lines": ["a line, BRCA2"]}]}])
        tag = ((v.get("lines") or {}).get(card) or [None])[0]
        if tag != (None if card in WANT_LINES_NO_GENE else panel):
            out.append(f"the line of {card} tagged {tag}")
    for card in sorted(set(WANT_FACTS_NO_GENE) | set(WANT_PANELS)):
        v = live_results.veil({"sections": {}}, found("BRCA2"), [{"sections": [{"id": card, "facts": [["A count", "BRCA2"]]}]}])
        tag = ((v.get("facts") or {}).get(card) or [None])[0]
        if tag != (None if card in WANT_FACTS_NO_GENE else panel):
            out.append(f"the fact of {card} tagged {tag}")
    for card, label in WANT_FACT_PROSE:
        v = live_results.veil({"sections": {}}, found("BRCA2"), [{"sections": [{"id": card, "facts": [[label, "x"]]}]}])
        if ((v.get("facts") or {}).get(card) or [None])[0] != panel:
            out.append(f"the fact {label!r} of {card} is prose, yet not tagged by the panel tag")
    for key in sorted(set(WANT_HEADLINES_NO_GENE) | {"couple"}):
        sec = {"id": "c", "keys": [key], "headline": "BRCA2 and more"}
        v = live_results.veil({"sections": {}}, found("BRCA2"), [{"sections": [sec]}])
        if ((v.get("headline") or {}).get("c")) != (None if key in WANT_HEADLINES_NO_GENE else panel):
            out.append(f"the headline from {key} tagged {(v.get('headline') or {}).get('c')}")
    # Prose read with no gene field: the panel tag whatever it names, so the tag never tells the words apart.
    for g, text in WANT_PROSE:
        for spans in (live_results.closable(found(g)), live_results.closable(found("BRCA2"))):
            if live_results.veil_id(spans, text) != panel:
                out.append(f"{text!r} tagged {live_results.veil_id(spans, text)}, want the panel tag")
    for cell, want in WANT_GENE_CELLS:
        if live_results.gene_cell(cell) != want:
            out.append(f"gene cell {cell!r} read as {live_results.gene_cell(cell)}, want {want}")
    return out


def table_controls():
    """Each entry of the veil's tables dropped in turn, and the word rule back: table_problems() must fail every time."""
    names = ("FIXED_PANELS", "FIXED_LINES", "NAMES_NO_GENE", "LINES_NO_GENE", "FACTS_NO_GENE", "FACT_NAMES",
             "HEADLINES_NO_GENE")
    keep = {n: getattr(live_results, n) for n in names}
    first_prose = live_results.prose_tag
    try:
        for c in WANT_PANELS:
            live_results.FIXED_PANELS = {k: v for k, v in keep["FIXED_PANELS"].items() if k != c}
            control(f"FIXED_PANELS without {c}", not table_problems())
        live_results.FIXED_PANELS = keep["FIXED_PANELS"]
        for card, words in keep["FIXED_LINES"].items():
            for w in words:
                live_results.FIXED_LINES = dict(keep["FIXED_LINES"], **{card: tuple(x for x in words if x != w)})
                control(f"FIXED_LINES[{card}] without {w!r}", not table_problems())
        live_results.FIXED_LINES = keep["FIXED_LINES"]
        live_results.NAMES_NO_GENE = {}
        control("NAMES_NO_GENE empty", not table_problems())
        live_results.NAMES_NO_GENE = keep["NAMES_NO_GENE"]
        for n in ("LINES_NO_GENE", "FACTS_NO_GENE", "FACT_NAMES", "HEADLINES_NO_GENE"):
            for k in keep[n]:
                setattr(live_results, n, {x: y for x, y in keep[n].items() if x != k})
                control(f"{n} without {k}", not (table_problems() or widget_problems()))
            setattr(live_results, n, keep[n])
        live_results.LINES_NO_GENE = dict(keep["LINES_NO_GENE"], cpsr="its lines")
        control("LINES_NO_GENE with the CPSR card, whose lines name genes", not table_problems())
        live_results.LINES_NO_GENE = keep["LINES_NO_GENE"]
        live_results.prose_tag = lambda spans, text: live_results.veil_tags(
            spans, re.findall(r"[A-Za-z0-9][A-Za-z0-9-]*", str(text).upper()), text)
        control("a veil that reads prose word by word, the word rule", not table_problems())
        control("a veil that reads prose word by word, the two-hit line's AR", not widget_problems())
    finally:
        for n in names:
            setattr(live_results, n, keep[n])
        live_results.prose_tag = first_prose


# Every adapter is a fixed panel (live_results.FIXED_PANELS) or listed here as one whose rows say something about this
# genome alone: a new or renamed adapter turns this file red until someone decides which it is.
WANT_NOT_PANELS = ("qc", "callability", "variants", "clinvar", "clinical-filter", "cpsr", "slivar", "exomiser",
                   "structural-variants", "mito", "mosaic", "mosaic-mca", "carrier-list", "couple", "pharmcat",
                   "ancestry-fine", "mt-haplogroup", "y-chromosome", "archaic", "ancient", "roh", "prs", "telomere",
                   "viral", "consumer-export")


def table_columns():
    """Every column of every table() call in bin/live_results.py, read from its source."""
    import ast
    with open(os.path.join(REPO, "bin", "live_results.py")) as fh:
        tree = ast.parse(fh.read())
    return {e.value for n in ast.walk(tree) if isinstance(n, ast.Call) and getattr(n.func, "id", None) == "table"
            and n.args and isinstance(n.args[0], ast.List) for e in n.args[0].elts if isinstance(e, ast.Constant)}


def declared_problems(adapters=None):
    """What the veil's tables leave undeclared: a table column not in TABLE_COLUMNS, a fixed-panel card that is no
    adapter, an adapter in neither FIXED_PANELS nor WANT_NOT_PANELS."""
    adapters = live_results.ADAPTERS if adapters is None else adapters
    out = [f"column {c!r} is not in TABLE_COLUMNS" for c in sorted(table_columns() - set(live_results.TABLE_COLUMNS))]
    out += [f"FIXED_PANELS card {c} is no adapter" for c in live_results.FIXED_PANELS if c not in adapters]
    out += [f"adapter {c} is neither a fixed panel nor in WANT_NOT_PANELS" for c in adapters
            if c not in live_results.FIXED_PANELS and c not in WANT_NOT_PANELS]
    return out


def vocabulary_summary(gene):
    """A summary with two recessive carrier rows (inheritance AR), an insertion over NKD2 (type INS) and the score of
    'Triglycerides (TG)', plus a closed CPSR result in GENE."""
    car = [{"gene": g, "variant": v, "genotype": "het", "significance": "Pathogenic", "stars": 2, "inheritance": "AR",
            "condition": c} for g, v, c in (("CFTR", "c.1521_1523del", "Cystic fibrosis"),
                                            ("HBB", "c.20A>T", "Sickle cell disease"))]
    return {"sections": {
        "carrier": {"state": "ok", "data": {"genes_on_list": 2, "carriers": car}},
        "annotsv": {"state": "ok", "data": {"count": 1, "rows": [
            {"sv": "chr5:1000000-1000300", "type": "INS", "gt": "0/1", "length": 300, "acmg_class": "4", "caller": "Manta",
             "genes": ["NKD2"], "support": 1, "support_callers": ["Manta"]}]}},
        "exomiser": {"state": "ok", "data": {"terms": [], "genes": [
            {"rank": 1, "gene": "GJB2", "moi": "AR", "exomiser_score": 0.9, "variants": [{"hgvs": "c.35del"}]}]}},
        "roh": {"state": "ok", "data": {"autosomal_over_5mb": []}},
        "prs": {"state": "ok", "data": {"adjusted": True, "scores": [
            {"pgs_id": "PGS000066", "condition": "Triglycerides (TG)", "percentile": "63.3"}]}},
        "cpsr": {"state": "ok", "data": {"plp": [{"gene": gene, "variant": "c.68_69del", "final": "Pathogenic"}]}}}}


def named_tags(v):
    """The veil's places (rows, lines, facts, marks, ...) whose tag names a result, never the panel tag alone: the
    same nesting, emptied of the rest."""
    if v is None or (isinstance(v, list) and all(isinstance(x, str) for x in v)):  # a tag
        return v if v and v != [live_results.PANEL] else None
    if isinstance(v, dict):
        got = {k: named_tags(x) for k, x in v.items()}
        return {k: x for k, x in got.items() if x} or None
    got = [named_tags(x) for x in v]
    return got if any(got) else None


def vocabulary_problems():
    """What is wrong with the veil's tags of the report's own vocabulary: a closed result in AR (an inheritance code),
    INS (a structural-variant type) or TG (a trait's abbreviation in a score's name) tags nothing, while one in CFTR
    tags the CFTR carrier row, one in GJB2 the symptom ranking's row and one in NKD2 the insertion's row, finding and
    chart mark."""
    marks = {(sp["id"], m["id"]) for sp in rr.chart_specs(vocabulary_summary("AR")) for m in sp["marks"]}
    out = [f"chart mark {m} not drawn" for m in (("genome", "sv-0"), ("prs", "pgs-PGS000066")) if m not in marks]
    for g in ("AR", "INS", "TG"):
        res = live_results.build("T", vocabulary_summary(g))
        v = {k: x for k, x in (res.get("veil") or {}).items() if k not in ("spans", "panels")}
        v = named_tags(v) or {}
        if any(v.values()):
            out.append(f"a closed result in {g} tags {v}")
    for g, where in (("CFTR", ("rows", "carrier-list")), ("GJB2", ("rows", "exomiser")),
                     ("NKD2", ("rows", "structural-variants")), ("NKD2", ("marks", "genome"))):
        v = live_results.build("T", vocabulary_summary(g)).get("veil") or {}
        if not any((v.get(where[0]) or {}).get(where[1]) or []):
            out.append(f"a closed result in {g} leaves {where} untagged: {v.get(where[0])}")
    v = live_results.build("T", vocabulary_summary("NKD2")).get("veil") or {}
    sv = [k for k in (v.get("findings") or {}) if k.startswith("sv-")]
    if not sv:
        out.append(f"a closed result in NKD2 leaves the insertion's finding untagged: {v.get('findings')}")
    return out


def vocabulary_controls():
    """Each guard of vocabulary_problems() dropped in turn: the check must fail every time."""
    cols, names, marks = live_results.TABLE_COLUMNS, live_results.finding_names, live_results.mark_names
    try:
        live_results.TABLE_COLUMNS = dict(cols, Inheritance="gene")
        control("TABLE_COLUMNS that reads the Inheritance column as naming a gene", not vocabulary_problems())
        live_results.TABLE_COLUMNS = dict(cols, Type="gene")
        control("TABLE_COLUMNS that reads the SV Type column as naming a gene", not vocabulary_problems())
        live_results.TABLE_COLUMNS = cols
        live_results.finding_names = lambda f: (re.findall(r"[A-Z0-9]+", f"{f.get('title') or ''}"), f.get("detail") or "")
        control("a veil that reads a finding's title word by word", not vocabulary_problems())
        live_results.finding_names = names
        live_results.mark_names = lambda S, sp, m: (re.findall(r"[A-Z0-9]+", f"{m.get('label') or ''} {m.get('detail') or ''}"), "")
        control("a veil that reads a chart mark's label word by word", not vocabulary_problems())
    finally:
        live_results.TABLE_COLUMNS, live_results.finding_names, live_results.mark_names = cols, names, marks
    try:
        live_results.TABLE_COLUMNS = {k: x for k, x in cols.items() if k != "Inheritance"}
        control("TABLE_COLUMNS without the Inheritance column", not declared_problems())
    finally:
        live_results.TABLE_COLUMNS = cols
    renamed = {("traits" if k == "traits-panel" else k): f for k, f in live_results.ADAPTERS.items()}
    control("the traits-panel adapter renamed traits", not declared_problems(renamed))
    control("a new adapter nobody classified", not declared_problems(dict(live_results.ADAPTERS, newcard=None)))


def overview_problems():
    """What is wrong with the home page's view of a child's run (overview): a finding that names a closed result
    (a symptom ranking of its gene) is listed among the first findings, or the closed lead still names it."""
    r = {"attention": 0, "levels": {}, "chapters": [],
         "findings": [{"id": "sf-gla-test", "level": "act", "title": "Secondary finding in GLA"},
                      {"id": "exomiser-gla", "level": "review", "title": "Symptom ranking: GLA ranks 1"},
                      {"id": "pgx-demo", "level": "medicine", "title": "Medicine panel finding"},
                      {"id": "chip-tp53", "level": "review", "title": "A blood clone in TP53"}],
         "lead": [["Of everything, ", False], [": Secondary finding in GLA. ", False, "act"],
                  ["One result is worth a look", True], [": Symptom ranking: GLA ranks 1.", False]],
         "veil": {"findings": {"exomiser-gla": ["sf-gla-test"], "pgx-demo": [live_results.PANEL],
                               "chip-tp53": [live_results.PANEL]}, "lead": [None, None, None, ["sf-gla-test"]]}}
    keep = live_results.build
    try:
        live_results.build = lambda *a, **k: r
        o = live_results.overview({})
    finally:
        live_results.build = keep
    out = []
    if [f["id"] for f in o["top"]] != ["sf-gla-test", "pgx-demo", "chip-tp53"]:
        out.append(f"the home page's first findings: {o['top']}")
    for f in o["top"]:
        if f["id"] in ("pgx-demo", "chip-tp53") and f.get("veil") != [live_results.PANEL]:
            out.append(f"home panel finding has no panel tag: {f}")
    if "GLA" in o["lead_closed"] or "worth a look" not in o["lead_closed"]:
        out.append(f"the closed lead: {o['lead_closed']}")
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

        # 8.
        import results_index
        sys.path.insert(0, os.path.join(REPO, "tests", "demo"))
        check("8. DEMO-001, with no result to act on, has no veil tag", res.get("veil") == {}, res.get("veil"))
        subprocess.run([sys.executable, os.path.join(REPO, "tests", "demo", "make_demo_sample.py"), "--genome-dir", root,
                        "--records", "2000", "--sample", "DEMO-003", "--sex", "male", "--seed", "3"], check=True,
                       stdout=subprocess.DEVNULL)
        d3 = os.path.join(root, "DEMO-003")

        def tagged(sample_dir):
            s3 = collect_summary.collect("DEMO-003", sample_dir)
            r3 = live_results.build("DEMO-003", s3)
            r3["_summary"] = s3
            return r3, results_index.build(s3, sample_dir, rr.SECTION_OF)
        bad = veil_lines_problems()
        check("8. every known deep intronic site's line carries the fixed-panel tag, not a closed result's id, a ClinVar "
              "row closed by its gene closes its reads link, and a lead segment carries the ids of every closed result "
              "its findings name", not bad, bad)
        fixed = live_results.FIXED_LINES
        try:  # the card still prints its lines; veil() no longer knows they are a panel
            live_results.FIXED_LINES = type("Unknown", (dict,), {"get": lambda self, k, d=None: None})(fixed)
            control("a veil that tags a known deep intronic site's line by its gene", not veil_lines_problems())
        finally:
            live_results.FIXED_LINES = fixed
        bad = table_problems()
        check("8. every fixed-panel card, fixed line and gene-free card of the veil's tables closes as a panel or never, "
              "prose is the panel tag whatever it names unless declared gene-free, and a gene cell is its symbols", not bad, bad)
        table_controls()
        bad = declared_problems()
        check("8. every table column is declared as naming a gene or not (TABLE_COLUMNS), and every adapter as a fixed "
              "panel or not", not bad, bad)
        bad = vocabulary_problems()
        check("8. the report's own codes (AR as an inheritance, INS as an SV type, TG in a score's name) tag nothing; a "
              "row, finding or mark is tagged by its gene and place only", not bad, bad)
        vocabulary_controls()
        bad = overview_problems()
        check("8. a child's home tile lists no finding that names a closed result, and its closed lead names none",
              not bad, bad)
        overview = live_results.overview
        try:
            def no_panel(summary):
                r = overview(summary)
                for f in r["top"]:
                    f.pop("veil", None)
                return r
            live_results.overview = no_panel
            control("overview drops panel tags from medicine and review findings", not overview_problems())
        finally:
            live_results.overview = overview
        bad = widget_problems()
        check("8. a CHIP call, a copy-number dot and the two-hit fact and line carry the id of a closed result in their "
              "gene field, never by their words (AR), the CHIP chart's gene panel none, a fixed-panel card's facts the "
              "panel tag, and results.json carries the closed results as veil.spans", not bad, bad)
        bad = sv_genes_problems()
        check("8. a closed structural variant over two genes closes the rows and the finding that name only its second "
              "gene, veil.spans lists it per gene, veil.panels lists the fixed-panel cards, and a duplication over two "
              "closed results carries both ids", not bad, bad)
        bad = gene_and_position_problems()
        check("8. gene-only ClinVar and second-gene SV chart marks, and coordinate-only mosaic findings, rows and lines "
              "carry the closed result's tag", not bad, bad)
        names, spans = live_results.mark_names, live_results.closable
        try:
            for kind in ("hit", "sv"):
                live_results.mark_names = lambda S, sp, m: (([], "") if kind == "hit" else (names(S, sp, m)[0][:1], "")) \
                    if m.get("kind") == kind else names(S, sp, m)
                control(f"{kind} chart mark without every gene", not gene_and_position_problems())
            live_results.mark_names = names
            live_results.closable = lambda found: [(fid, g, None, None, None) for fid, g, *_ in spans(found)]
            control("closable without positions", not gene_and_position_problems())
        finally:
            live_results.mark_names, live_results.closable = names, spans
        first = live_results.closable
        try:
            live_results.closable = lambda found: first([{k: v for k, v in f.items() if k != "genes"} for f in found])
            control("a veil that tags only a structural variant's first gene", not sv_genes_problems())
        finally:
            live_results.closable = first
        bad = veil_problems(*tagged(d3))
        check("8. the planted ClinVar hit's row, viewer link, tick and star bar carry its secondary finding's id, the "
              "deletion's span its own, every other card's row of the gene and the symptom ranking's finding one of them, "
              "and every search item of the gene one of them", not bad, bad)
        # The planted mistakes: a rule that tags nothing, and one that tags every row with the first result.
        keep = live_results.veil_tags
        try:
            live_results.veil_tags = lambda spans, *a, **k: None
            control("a veil rule that tags nothing", not veil_problems(*tagged(d3)))
            control("a veil rule that tags no line", not veil_lines_problems())
            control("a veil rule that tags no fact or drawing", not widget_problems())
            live_results.veil_tags = lambda spans, *a, **k: [spans[0][0]] if spans else None
            control("a veil rule that tags every row", not veil_problems(*tagged(d3)))
            live_results.veil_tags = lambda spans, *a, **k: (keep(spans, *a, **k) or [None])[:1] if keep(spans, *a, **k) else None
            control("a veil rule that keeps only the first id of a place", not veil_problems(*tagged(d3)))
            control("a veil rule that keeps only the first id of a place, two closed results", not sv_genes_problems())
        finally:
            live_results.veil_tags = keep
        first_union = live_results.union
        try:
            live_results.union = lambda *t: next((x[:1] for x in t if x), None)
            control("a veil rule that keeps only the first id of a lead segment", not veil_lines_problems())
        finally:
            live_results.union = first_union

        # 3, the map documents: og map build of DEMO-001 and DEMO-003, with the results each page gets.
        map_docs = map_documents(root, [(SAMPLE, res), ("DEMO-003", tagged(d3)[0])])
        if map_docs is None:
            ci = os.environ.get("CI") == "true"
            print(("[FAIL] " if ci else "[SKIP] ") + "3. no og (and no cargo to build one) for the map documents")
            if ci:
                FAILS.append("no og in CI")
        else:
            bad = [(run, x) for run, doc in map_docs for x in map_wording(doc)]
            check(f"3. no string of the genome map documents reads as an all-clear ({len(map_docs)} documents)",
                  not bad and len(map_docs) == 2, bad)
            planted = copy.deepcopy(map_docs[0][1])
            planted["text"]["genome"] = [{"t": "All clear on the map.", "b": False}]
            control("'All clear' planted in a map document's sentence", not map_wording(planted))

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
