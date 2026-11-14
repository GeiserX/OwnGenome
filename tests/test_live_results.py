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
     planted secondary finding and class 5 deletion): the ClinVar row and viewer
     link of the planted hit, its genome-chart tick and its review-star bar carry
     the secondary finding's id, the deletion's span carries its own, the search
     index's items that name the gene carry one of them, and nothing else is
     tagged; DEMO-001, with no result to act on, has no tag; and a known deep
     intronic site's line in the gene of a closed secondary finding carries
     that finding's id (the card prints each site whether hit or not), a gene
     named before a full stop is that gene, and a ClinVar row closed by its gene
     alone closes its "See the reads at" link with it; every other card's row
     that names the gene (the son's symptom ranking ranks it) carries one of
     the ids and no other row does, the ranking's finding carries the secondary
     finding's id, and so does a segment of the lead that names the gene; the
     home page's view of the run (overview) lists no such finding among its
     first ones and its closed lead names none; a CHIP call, a copy-number dot
     and a fact that name a closed result's gene carry its id (the CHIP chart's
     panel of genes never does), and results.json carries the closed results as
     veil.spans (the family page closes a child's places while one is closed).
Each check runs again on a planted mistake and must fail there: 1 with VNTYPER's
entry dropped, 2 with one card's caveat blanked, 3 with "ok" planted in a card
and "All clear" in a script, 4 with the hits emptied, 5 with the results route
registered without the login, 6 with a row printed as its Python dict, 7 with the Slivar adapter returning nothing,
8 with a veil rule that tags nothing (rows, lines, facts and drawings) and one that tags every row. A check
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
    rows = (v.get("rows") or {}).get("clinvar") or []
    want = [sf if r[0] == gene and r[1].startswith(f"{c}:{p} ") else None for r in card["table"]["rows"]]
    if rows != want:
        out.append(f"ClinVar rows tagged {rows}, want {want}")
    view = (v.get("viewer") or {}).get("clinvar") or []
    if view != [sf if x == f"{c}:{p}" else None for x in card.get("viewer") or []]:
        out.append(f"ClinVar viewer links tagged {view}")
    marks = v.get("marks") or {}
    specs = {sp["id"]: sp for sp in rr.chart_specs(res["_summary"])}
    hit = [m["id"] for m in specs["genome"]["marks"] if m.get("kind") == "hit" and m["pos"] == p and m["chrom"] == c]
    span = [m["id"] for m in specs["genome"]["marks"] if m.get("kind") == "sv" and (m["chrom"], m["pos"], m["end"]) == (sv_c, a, b)]
    bar = [m["id"] for m in specs["clinvar-stars"]["marks"] if gene in (m.get("detail") or "")]
    want = {"genome": {**{h: sf for h in hit}, **{s_: sv for s_ in span}}, "clinvar-stars": {x: sf for x in bar}}
    if not hit or not span or not bar or marks != want:
        out.append(f"chart marks tagged {marks}, want {want}")
    # Every other card's rows: a row that names the gene carries one of the two ids, a row that does not carries
    # none; the symptom ranking's row of the gene carries the secondary finding's, and so does its finding.
    for cid, rows_ in ((x["id"], (x.get("table") or {}).get("rows") or []) for ch in res["chapters"] for x in ch["sections"]):
        if cid == "clinvar" or not rows_:
            continue
        got = (v.get("rows") or {}).get(cid) or [None] * len(rows_)
        bad = [(r, t) for r, t in zip(rows_, got) if (t in (sf, sv)) != (gene in " ".join(map(str, r)))]
        if bad or len(got) != len(rows_):
            out.append(f"{cid} rows tagged wrong: {bad or got}")
    exo = (v.get("rows") or {}).get("exomiser") or []
    ranked = [f["id"] for f in res["findings"] if f["level"] != "act" and gene in f["title"]]
    if sf not in exo or not ranked or v.get("findings") != {i: sf for i in ranked}:
        out.append(f"the symptom ranking of {gene}: rows tagged {exo}, findings tagged {v.get('findings')}, want {ranked} -> {sf}")
    tagged = {(t, it["section"]): it.get("veil") for t, items in idx["terms"].items() for it in items if it.get("veil")}
    named = {(gene, it["section"]) for it in idx["terms"].get(gene) or []}
    if set(tagged) != named or not named or any(x not in (sf, sv) for x in tagged.values()):
        out.append(f"search items tagged {tagged}, want every item of {gene}: {sorted(named)}")
    return out


def veil_lines_problems():
    """What is wrong with the veil's tags of the ClinVar card's lines: a known deep intronic site is printed
    whether hit or not, so its line names the gene of a closed secondary finding in that gene."""
    S = {"clinvar": {"state": "ok", "data": {"hits": [], "deep_intronic": [
        {"gene": "GLA", "hgvs": "NM_000169.3:c.640-801G>A", "chrom": "chrX", "pos": 101399747, "ref": "C",
         "alt": "T", "call": "hom", "detail": "GQ 40"},
        {"gene": "GENE5", "hgvs": "c.100+500A>G", "chrom": "chr1", "pos": 1, "ref": "A", "alt": "G",
         "call": "het", "detail": "GQ 40"}]}}}
    lines = live_results.a_clinvar(S, {"sections": S}).get("lines") or []
    found = [{"id": "sf-gla-test", "level": "act", "gene": "GLA", "detail": "GLA"}]
    v = live_results.veil({"sections": S}, found, [{"sections": [{"id": "clinvar", "lines": lines}]}])
    want = ["sf-gla-test" if "GLA " in ln else None for ln in lines]
    got = (v.get("lines") or {}).get("clinvar")
    out = [] if "sf-gla-test" in want and None in want and got == want else [f"lines {lines} tagged {got}, want {want}"]
    # A gene named at the end of a sentence ("in GLA.") is that gene.
    if live_results.veil_id(live_results.closable(found), "A change in GLA.") != "sf-gla-test":
        out.append("'A change in GLA.' carries no tag")
    # A segment of the lead that names the gene, beside the act names the page closes itself.
    lead = [[": Symptom ranking: GLA ranks 1. ", False], [": Secondary finding in GLA. ", False, "act"], ["Mean depth", False]]
    got = live_results.veil({"sections": S}, found, [], lead).get("lead")
    if got != ["sf-gla-test", None, None]:
        out.append(f"the lead's segments tagged {got}, want the names that are not the act names")
    # A finding known by its gene alone (no position): its ClinVar row is closed by the gene, and the row's
    # "See the reads at" link, which shows only the position, with it.
    S["clinvar"]["data"]["hits"] = [{"gene": "GLA", "chrom": "chrX", "pos": 101398000, "ref": "G", "alt": "A"},
                                    {"gene": "GENE5", "chrom": "chr1", "pos": 5, "ref": "A", "alt": "G"}]
    card = dict(live_results.a_clinvar(S, {"sections": S}), id="clinvar")
    gene_only = [{"id": "carrier-gla-test", "level": "act", "gene": "GLA"}]
    v = live_results.veil({"sections": S}, gene_only, [{"sections": [card]}])
    rows, view = (v.get("rows") or {}).get("clinvar"), (v.get("viewer") or {}).get("clinvar")
    if rows != ["carrier-gla-test", None] or view != rows:
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
    want = {"chip": [None, "cpsr-tp53-test"], "paralogs": ["two-hits-smn1-test", None]}
    if v.get("widget") != want:
        out.append(f"drawings tagged {v.get('widget')}, want {want}")
    if (v.get("facts") or {}) != {"paralogs": ["two-hits-smn1-test", None]}:
        out.append(f"facts tagged {v.get('facts')}")
    if v.get("spans") != [list(x) for x in live_results.closable(found)]:
        out.append(f"spans {v.get('spans')}")
    return out


def overview_problems():
    """What is wrong with the home page's view of a child's run (overview): a finding that names a closed result
    (a symptom ranking of its gene) is listed among the first findings, or the closed lead still names it."""
    r = {"attention": 0, "levels": {}, "chapters": [],
         "findings": [{"id": "sf-gla-test", "level": "act", "title": "Secondary finding in GLA"},
                      {"id": "exomiser-gla", "level": "review", "title": "Symptom ranking: GLA ranks 1"}],
         "lead": [["Of everything, ", False], [": Secondary finding in GLA. ", False, "act"],
                  ["One result is worth a look", True], [": Symptom ranking: GLA ranks 1.", False]],
         "veil": {"findings": {"exomiser-gla": "sf-gla-test"}, "lead": [None, None, None, "sf-gla-test"]}}
    keep = live_results.build
    try:
        live_results.build = lambda *a, **k: r
        o = live_results.overview({})
    finally:
        live_results.build = keep
    out = []
    if [f["id"] for f in o["top"]] != ["sf-gla-test"]:
        out.append(f"the home page's first findings: {o['top']}")
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
        check("8. a known deep intronic site's line in a closed secondary finding's gene carries its id, a gene "
              "before a full stop is that gene, a ClinVar row closed by its gene closes its reads link, and a lead "
              "segment naming the gene carries its id", not bad, bad)
        bad = overview_problems()
        check("8. a child's home tile lists no finding that names a closed result, and its closed lead names none",
              not bad, bad)
        bad = widget_problems()
        check("8. a CHIP call, a copy-number dot and a fact that name a closed result's gene carry its id, the CHIP "
              "chart's gene panel none, and results.json carries the closed results as veil.spans", not bad, bad)
        bad = veil_problems(*tagged(d3))
        check("8. the planted ClinVar hit's row, viewer link, tick and star bar carry its secondary finding's id, the "
              "deletion's span its own, every other card's row of the gene and the symptom ranking's finding one of them, "
              "and every search item of the gene one of them", not bad, bad)
        # The planted mistakes: a rule that tags nothing, and one that tags every row with the first result.
        keep = live_results.veil_id
        try:
            live_results.veil_id = lambda spans, *a, **k: None
            control("a veil rule that tags nothing", not veil_problems(*tagged(d3)))
            control("a veil rule that tags no line", not veil_lines_problems())
            control("a veil rule that tags no fact or drawing", not widget_problems())
            live_results.veil_id = lambda spans, *a, **k: spans[0][0] if spans else None
            control("a veil rule that tags every row", not veil_problems(*tagged(d3)))
        finally:
            live_results.veil_id = keep

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
