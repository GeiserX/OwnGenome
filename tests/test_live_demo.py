#!/usr/bin/env python3
"""The live site's public demo (tests/demo/build_live_demo.py), built and checked.

Checks:
  1. tests/fixtures/live/demo-001.weblog.jsonl is what tests/demo/make_demo_events.py
     writes from assets/steps.tsv now;
  2. every state hook the plugins bin/owngenome_live_*.py register has a value
     in the last frame of DEMO-001's page (the family page for the family
     plugin's own view), none of them an error;
  3. the frames follow the run: running at the first, finished at the last with
     every replayed step done; the charts and the search index appear only
     once the report step has ended, the re-check only after the run;
  4. the run-health card measures the invented machine, never the machine that
     builds the demo, and no frame holds a path of that machine;
  5. every link a frame gives is relative and reaches a file of the demo, and
     so is every link of runs.json, results.json and family.json;
  6. every page the build writes carries the server's Content-Security-Policy
     as a meta tag, and no page or script of the demo loads anything from
     outside;
  7. in a headless Chrome (skipped without Chrome, except in CI), served as a
     static host serves it (no policy header): the home page's people and
     family, the replay to its last frame, the charts, a search, the re-check
     banner, the run-health card, the family link and its pedigree, the viewer
     page, and no policy violation; the person page's ten sections, its
     contents, lead and ledger, #s/clinvar focusing the ClinVar card with its
     caveat, chart and flag, a finding's own page with its plan and evidence;
     the family page with a
     clickable pedigree and no step map or report link; the widgets of the
     partner's and the son's pages (the repeat gauge, the Lp(a) interval, the
     CHIP VAF, the copy-number dots, the Y path) and the AnnotSV and MoChA marks
     on their genome charts; the family page's relatedness grid, the son's
     painted stretches (one mark per row of his segments file), his de novo
     and compound-het tables and the members side by side; and at a phone's
     width (390 px) no page scrolls sideways (every section open, the finding
     page too) and the person page has its bottom tab bar; at 390 px the
     partner's page is under 8,000 CSS px with every finding to act on, for a
     prescriber, for the family or worth a look open, the rest folded, and a
     link to a card two folds deep (#s/roh) or to a folded section (#s7) opens
     it, while a wider page folds nothing; a section 5 card flagged with no
     finding of its own stays open on a phone, the folds a reader opened (a
     card's other lines too) stay open after the page turns wide and back, and
     a Ctrl+click on a link to the current place opens nothing; on the son's page a
     two planted CPSR results to act on name their genes and variants nowhere
     (their card included) until "Show the result" is clicked, and each click
     opens one result, there and on the finding's own page; they stay closed on
     both pages with family.json unavailable, with a state cut short (no
     plugins) and with no state at all (frames.json unavailable), and the son's
     search lists no hit and suggests no term from their closed card while a
     term of an open card is still found; the son's planted secondary finding
     and class 5 deletion (make_demo_sample.py PLANTED) name their gene and
     positions in no card, chart, chart title, link, symptom ranking or search
     suggestion until each is opened, a row, mark or ranking that names both
     comes back only once both are opened (the search too), and a reload
     closes them again; a place that names two closed results (a lead segment,
     a duplication over two genes) stays closed until both are opened, and the
     page after a click is the same whichever closed result shares the opened
     one's place; the home page and
     every finding page of the son (an open finding's rail and related list
     too) name them nowhere, and the symptom ranking's own page shows the gene
     once opened on purpose; on the family page, while the son has a planted
     result closed (the fixture's compound-het pair in KCNH2, a de novo
     candidate and a stretch painted from a grandparent planted in family.json),
     his painting's grandparents, his de novo and compound-het rows and his
     medicines and HLA cells of the side by side close, one box each, and no
     text or attribute names its gene, positions or that grandparent; each click
     opens one result, every place comes back once both are opened, and a reload
     closes them again; the page is the same, byte for byte, whichever result is
     closed (KCNH2, RYR1, HLA-B, a gene alone, a position, a copy-number span,
     a CHIP call) and wherever a grandparent's stretch lies, each view settled
     on knowing his results; while the page cannot read his results (no
     results.json, an unreadable one, one cut short, act findings without
     veil.spans or with malformed ones) every place stays closed and each box
     opens only its own place; with nothing of his to act on every place is
     open with no box; with a second child holding the same results, opening
     his leaves hers closed, and results opened close again when his run moves;
     a couple screen whose partner is a child with a closed result closes whole
     (counts, table, boxes, report link, the lead's words), the same page
     whether that result is in a gene of the screen or not, and with the
     partner's results unreadable its box and the cells' box each open only
     their own place; his map as an image stays closed until each result is
     opened from its box, the same page whichever result is closed, and with
     his results unreadable "Show the rows" leaves it closed and "Show the map"
     opens it and no table; on
     his page the copy-number and CHIP cards, a fixed panel, close whole while
     a result of his is closed (on a phone their line too), their box opens one
     result per click and every dot, call and fact is back once all are open;
     two of his pages that differ only in which result is closed (two SNV
     genes, two deletions of two genes each, SMN1 against HBA1, a result with a
     place on the map against one without, one ClinVar row against two, a link
     to each) are the same before any click, HTML, attributes, address and map
     canvas, and a closed deletion's second gene is named nowhere; DEMO-001's page
     mid-replay and at the end leads with the report (a short run band above
     the lead while it runs, none once finished), ends with the step timeline
     folded, keeps the wire to at most eight one-sentence lines with the plain
     finishes folded into one, and puts the HLA card's caveat in its margin
     with its absent and undetermined lines folded.
Checks 1, 2, 3, 5, the widgets (person.js without its widget line), the family depth (family.json
with the son's stretches emptied), the phone width and the son's closed result (person.js whose cards ignore
the veil, person.js and finding.js that read a missing state as a standalone adult, and search.js
without its filter; for the planted results, results.json without its veil tags, person.js whose table
rows ignore them, person.js whose findings ignore them, charts.js that keeps every mark, search.js that ignores a hit's tag
and finding.js whose findings ignore their tags; for the family page, de novo rows, compound-het rows, a painting,
a side by side's cells, a couple screen and a map image that ignore the veil, a closed painting that keeps each
stretch's grandparent, a box that opens every result with one click, de novo rows closed only when a result lies
under one, a painting closed only when a grandparent's stretch lies over a closed result, a side by side that closes
only the closed gene's cells or stays open for a wide span, a page that never reads the son's results, one opening
shared by every box or 'Show the grandparents' opening the side by side while the results cannot be read, a failed
request read as nothing closed, malformed veil.spans read as given, results without veil.spans read as nothing
closed or as closed with nothing to act on, an opened result kept per run only or per child only, a couple screen
closed only for a result in its genes or opened by 'Show the cells', a map image closed only for a result with a
position, and de novo rows that 'Show the map' opens; for the son's page, person.js that draws a fixed-panel
card's open rows or builds its line from its facts, that puts a closed card's id back or counts the closed rows,
and a veil that tags only a deletion's first gene), the phone page (person.js that never folds, and a route that opens only the nearest fold), the
folds (a flagged section 5 card folded, a redraw that forgets the folds, a Ctrl+click that opens one) and
the person page's layout (a wire showing every line or every plain finish, a card showing
every line, a step timeline left open) also run against a planted mistake and must
fail there (2: a build with one plugin's state hook dropped): a check that cannot fail is no check.

The family veil group described above runs from tests/test_owngenome_live_family.py in another Guard shard.
Run both files for every browser check.

Run: python3 tests/test_live_demo.py
"""
import base64
import copy
import difflib
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.join(REPO, "tests", "demo"))
sys.path.insert(0, os.path.join(REPO, "tests"))
import build_live_demo  # noqa: E402
import live_results  # noqa: E402
import make_demo_events  # noqa: E402
import owngenome_live  # noqa: E402
import run_state  # noqa: E402
from test_owngenome_live import Chrome, asset_problems, find_chrome  # noqa: E402

FAILS = []
S, FAM = build_live_demo.SAMPLE, build_live_demo.FAMILY
OTHERS = [n for n, _, _ in build_live_demo.OTHERS]


def check(desc, ok, detail=""):
    print(("[PASS] " if ok else "[FAIL] ") + desc + ("" if ok else f" -- {str(detail)[:800]}"))
    if not ok:
        FAILS.append(desc)


def control(desc, still_passes):
    """A planted mistake: the check must fail on it."""
    check(f"control: {desc}", not still_passes, "the check passed on the planted mistake")


def registered_hooks():
    app = owngenome_live.App()
    owngenome_live.configure(["--root", tempfile.gettempdir()], app)
    return [n for n, _ in app.state_hooks]


def missing_hooks(built, hooks):
    """The hooks with no value (or an error) on the demo's pages."""
    last = built["frames"][-1]["state"]["plugins"]
    fam = built["family"]["state"]["plugins"]
    bad = []
    for h in hooks:
        v = last.get(h)
        if h == "family":
            v = v if (v or {}).get("kind") == "member" and (fam.get(h) or {}).get("kind") == "family" else None
        if not v or (isinstance(v, dict) and "error" in v):
            bad.append(h)
    return bad


def replayed_steps():
    procs = set()
    with open(build_live_demo.WEBLOG) as f:
        for line in f:
            tr = json.loads(line).get("trace") or {}
            if tr.get("name"):
                procs.add(run_state.task_name(tr["name"])[0])
    return sorted(s["step"] for s in run_state.read_steps() if set(s["processes"]) & procs)


def follow_problems(frames, steps):
    """What in the frames does not follow the run."""
    out = []
    first, last = frames[0]["state"], frames[-1]["state"]
    if first["run"]["state"] != "running":
        out.append(f"first frame {first['run']['state']}")
    if last["run"]["state"] != "finished":
        out.append(f"last frame {last['run']['state']}")
    st = {s["step"]: s["state"] for s in last["steps"]}
    out += [f"step {n} {st.get(n)} at the end" for n in steps if st.get(n) != "done"]
    for fr in frames:
        p, s24 = fr["state"]["plugins"], {s["step"]: s["state"] for s in fr["state"]["steps"]}.get("24")
        for k in ("charts", "search"):
            if bool((p.get(k) or {}).get("file")) != (s24 == "done"):
                out.append(f"{k} at {fr['t']} with step 24 {s24}")
        if p.get("recheck") and fr is not frames[-1]:
            out.append(f"re-check at {fr['t']}, before the run ended")
    if not (last["plugins"].get("recheck") or {}).get("steps"):
        out.append("no re-check badge at the end")
    return out


def links(v):
    if isinstance(v, dict):
        for k, x in v.items():
            # a relatedness row's "segments" is a count, not a file
            if k in ("url", "pedigree", "segments", "report", "table") and isinstance(x, str) and not x.isdigit():
                yield x
            else:
                yield from links(x)
    elif isinstance(v, list):
        for x in v:
            yield from links(x)


def link_problems(out, built):
    """Links of the frames that are absolute or reach no file of the demo."""
    bad = []
    pages = [(S, fr["state"]) for fr in built["frames"]] + [(FAM, built["family"]["state"])]
    for run, st in pages:
        base = os.path.join(out, "runs", run)
        for u in links(st["plugins"]):
            if u.startswith("/"):
                bad.append(f"{run}: absolute {u}")
                continue
            p = os.path.normpath(os.path.join(base, u))
            if not (os.path.isfile(p) or os.path.isfile(os.path.join(p, "index.html"))):
                bad.append(f"{run}: {u} reaches nothing")
        for k, v in st["plugins"].items():
            if isinstance(v, dict) and isinstance(v.get("file"), str) and not os.path.isfile(os.path.join(base, v["file"])):
                bad.append(f"{run}: {k} file {v['file']} not in the demo")
    return bad


def json_link_problems(out):
    """Links of runs.json, results.json and family.json that are absolute or reach no file of the demo."""
    bad = []
    for rel, base in [("runs.json", out), (os.path.join("runs", FAM, "family.json"), os.path.join(out, "runs", FAM))] + \
            [(os.path.join("runs", n, "results.json"), os.path.join(out, "runs", n)) for n in [S] + OTHERS]:
        with open(os.path.join(out, rel)) as f:
            data = json.load(f)
        found = list(links(data))
        if rel.endswith("results.json"):
            found += [x for ch in data["chapters"] for c in ch["sections"] for x in c.get("links") or []]
        if not found:
            bad.append(f"{rel}: no link at all")
        for u in found:
            if u.startswith("https://"):  # an explore link of the report (a haplogroup tree), opened by the reader
                continue
            p = os.path.normpath(os.path.join(base, u))
            if u.startswith("/") or not (os.path.isfile(p) or os.path.isfile(os.path.join(p, "index.html"))):
                bad.append(f"{rel}: {u}")
    return bad


def located_term(out):
    """A term of the demo's search index with a mention at a locus (its hit links the viewer)."""
    with open(os.path.join(out, "runs", S, f"{S}_search_index.json")) as f:
        terms = json.load(f)["terms"]
    return next((t for t, items in sorted(terms.items())
                 if any(m.get("locus") for it in items for m in it.get("mentions") or [])), None)


class Static(SimpleHTTPRequestHandler):
    """A static host: files only, no Content-Security-Policy header. Nothing is cached: the checks rewrite
    results.json, map.json and the scripts between page loads, and a heuristically fresh copy in the browser's cache
    (10% of its age since Last-Modified) would load the old file and let a control pass on the planted mistake."""

    def log_message(self, *a):
        pass

    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        super().end_headers()


def browser_checks(out):
    chrome = find_chrome()
    if not chrome:
        ci = os.environ.get("CI") == "true"
        print(("[FAIL] " if ci else "[SKIP] ") + "no Chrome or Chromium to load the demo in")
        if ci:
            FAILS.append("no Chrome in CI")
        return
    handler = type("H", (Static,), {"__init__": lambda self, *a, **k: SimpleHTTPRequestHandler.__init__(
        self, *a, directory=out, **k)})
    srv = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{srv.server_address[1]}"
    b = Chrome(chrome, os.path.join(os.path.dirname(out), "chrome"))
    try:
        b.go(f"{base}/")
        people = b.wait("document.documentElement.dataset.pgpStatus === 'ready' && document.documentElement.dataset.pgpPeople")
        runs = b.eval("Array.from(document.querySelectorAll('.person-card a.person-link')).map(function (a) { "
                      "return a.getAttribute('href'); }).join(' ') + ' | ' + Array.from(document.querySelectorAll("
                      "'.family-head a')).map(function (a) { return a.getAttribute('href'); }).join(' ')")
        want = " ".join(f"runs/{n}/" for n in [S] + OTHERS) + f" | runs/{FAM}/ runs/{FAM}/"
        check("browser: the home page shows the three people of the family, each with a run linked, and the family",
              people == "3" and runs == want, (people, runs))
        b.go(f"{base}/runs/{S}/")
        first = b.wait("document.documentElement.dataset.pgpFrame && document.documentElement.dataset.pgpStatus === 'ready' "
                       "&& document.querySelectorAll('#map g.box').length")
        b.eval("var s = document.getElementById('replay-slider'); s.value = s.max; "
               "s.dispatchEvent(new Event('input')); true")
        n = b.wait("document.documentElement.dataset.pgpFrame === document.getElementById('replay-slider').max "
                   "&& document.documentElement.dataset.pgpCharts")
        boxes = b.eval("document.querySelectorAll('#map g.box').length")
        check("browser: the replay draws the step map, which grows as steps are planned",
              isinstance(first, int) and first >= 1 and isinstance(boxes, int) and boxes > 30, (first, boxes))
        check("browser: at the last frame the charts are drawn from the demo's charts file", n and int(n) > 0, n)
        state = b.eval("document.getElementById('run-state').textContent")
        check("browser: the last frame shows the run finished", state == "finished", state)
        term = located_term(out)
        ready = b.wait("var i = document.querySelector('input[type=search]'); i && !i.disabled && 'enabled'")
        check("browser: the search box opens on the demo's index", ready == "enabled",
              b.eval("Array.from(document.querySelectorAll('section.plugin')).map(function (s) { "
                     "return s.dataset.plugin + ': ' + s.textContent.slice(0, 120); }).join(' | ')"))
        b.eval("var i = document.querySelector('input[type=search]'); i.value = %s; "
               "i.dispatchEvent(new Event('input')); true" % json.dumps(term))
        hits = b.wait("!document.querySelector('input[type=search]').disabled && document.documentElement.dataset.pgpSearch")
        check(f"browser: a search for {term} finds it in the demo's index", hits and int(hits) > 0, (term, hits))
        cards = b.eval("Array.from(document.querySelectorAll('section.plugin')).map(function (s) { "
                       "return s.dataset.plugin + '=' + (s.hidden ? '' : s.textContent.length); }).join(' ')")
        banner = b.eval("(document.querySelector('section.plugin[data-plugin=recheck]') || {}).textContent || ''")
        check("browser: the re-check card shows the newer ClinVar", "newer ClinVar" in banner, banner[:300])
        ops = b.eval("(document.querySelector('section.plugin[data-plugin=ops]') || {}).textContent || ''")
        check("browser: the run-health card says the run ended and names the retried task",
              "retried" in ops, ops[:400])
        famlink = b.eval("var a = document.querySelector('section.plugin[data-plugin=family] a'); a && a.getAttribute('href')")
        check("browser: the family card links the family page", famlink == f"../../runs/{FAM}/", (famlink, cards))
        viewer = b.eval("var a = document.querySelector('a.viewer'); a && a.href")
        if viewer:
            b.go(viewer)
            text = b.wait("document.body && document.body.textContent.indexOf('not in the demo') >= 0 && 'ok'")
            check("browser: a viewer link opens the page that says the viewer is not in the demo", text == "ok", viewer)
        else:
            check("browser: the charts or the search give a viewer link", False, cards)
        b.go(f"{base}/runs/{FAM}/")
        ped = b.wait("var i = document.querySelector('img.pedigree'); i && i.complete && i.naturalWidth > 0 && 'drawn'")
        check("browser: the family page shows the pedigree", ped == "drawn",
              b.eval("var i = document.querySelector('img.pedigree'); i && i.src"))
        blocked = b.eval("document.documentElement.dataset.pgpBlocked || ''")
        check("browser: the policy blocked nothing on the family page", not blocked, blocked)
        b.go(f"{base}/runs/{S}/")
        b.wait("document.documentElement.dataset.pgpStatus === 'ready'")
        b.eval("var s = document.getElementById('replay-slider'); s.value = s.max; s.dispatchEvent(new Event('input')); true")
        b.wait("document.documentElement.dataset.pgpCharts")
        time.sleep(1)
        blocked = b.eval("document.documentElement.dataset.pgpBlocked || ''")
        check("browser: the policy blocked nothing on the sample page", not blocked, blocked)
        csp = b.eval("var m = document.querySelector('meta[http-equiv=Content-Security-Policy]'); m && m.content")
        check("browser: the page runs under the meta policy", csp == build_live_demo.CSP, csp)
        person_checks(b, base)
        layout_checks(b, base, out)
        widget_checks(b, base, out)
        family_depth_checks(b, base, out)
        phone_checks(b, base)
        phone_report_checks(b, base, out)
        fold_checks(b, base, out)
        veil_checks(b, base, out)
        planted_checks(b, base, out)
        map_checks(b, base, out)
        # The family veil group runs from test_owngenome_live_family.py in a different Guard shard.
        widget_veil_checks(b, base, out)
        infer_checks(b, base, out)
    finally:
        b.close()
        srv.shutdown()
        srv.server_close()


def person_checks(b, base):
    """The person page at the end of the run: the ten sections, the lead, the contents, a section link, a finding
    page; the family page's own view."""
    got = b.wait("document.documentElement.dataset.pgpResults && JSON.stringify(["
                 "document.querySelectorAll('main section.sec[data-journal]').length, "
                 "document.querySelectorAll('#rail [data-toc]').length, "
                 "(document.querySelector('#lead-slot .lede') || {textContent: ''}).textContent.length > 40, "
                 "document.querySelectorAll('#lead-slot .ledger a[data-level]').length])")
    check("browser: the person page has its ten sections, the contents linking each, a lead paragraph and the ledger",
          got == "[10,10,true,5]", got)
    b.eval("location.hash = '#s/clinvar'; true")
    focus = b.wait("document.activeElement && document.activeElement.id === 'card-clinvar' && 'card-clinvar'")
    card = b.eval("var c = document.getElementById('card-clinvar'); c && [!!c.querySelector('.caveat .caveat-lead'), "
                  "!!c.querySelector('figure[data-chart]'), c.classList.contains('attn'), "
                  "c.getBoundingClientRect().height > 0].join(' ')")
    check("browser: #s/clinvar focuses the ClinVar card, shown, with its caveat, its chart and its flag",
          focus == "card-clinvar" and card == "true true true true", (focus, card))
    link = b.eval("var a = document.querySelector('[data-finding] h3 a[href^=\"finding/\"]'); a && a.getAttribute('href')")
    check("browser: a finding links its own page", bool(link), link)
    if link:
        b.go(f"{base}/runs/{S}/{link}")
        page = b.wait("document.documentElement.dataset.pgpStatus === 'ready' && JSON.stringify(["
                      "!!document.querySelector('h1.h-display'), document.querySelectorAll('#plan li').length > 0, "
                      "!!document.getElementById('evidence'), document.documentElement.dataset.pgpBlocked || ''])")
        check("browser: the finding page shows its title, its plan and its evidence, and the policy blocks nothing",
              page == '[true,true,true,""]', (link, page))
    b.go(f"{base}/runs/{FAM}/")
    fam = b.wait("document.documentElement.dataset.pgpStatus === 'ready' && [!document.getElementById('map'), "
                 "!document.getElementById('report-link'), (document.querySelector('svg.pedigree-svg a') || {"
                 "getAttribute: function () { return null; } }).getAttribute('href')].join(' ')")
    check("browser: the family page has a clickable pedigree, no step map and no report link",
          fam == f"true true ../../runs/{S}/", fam)


# The widgets: (person, card, kind, the marks it must draw) on the demo's other people.
WIDGETS = [("DEMO-002", "expansions", "str", "[data-widget=str] .w-row", 8),
           ("DEMO-002", "lpa", "ci", "[data-widget=ci] .w-tick", 1),
           ("DEMO-002", "chip", "vaf", "[data-widget=vaf] .w-tick.w-call", 1),
           # five copy-number genes from Paraphase, then Parascopy's SMN1 and SMN2, on the partner's page: the son's keeps
           # the copy-number card closed while a result of his is closed (a fixed panel)
           ("DEMO-002", "paralogs", "copies", "[data-widget=copies] li", 7),
           ("DEMO-003", "y-chromosome", "ypath", "[data-widget=ypath] li", 5),
           ("DEMO-002", "clinvar", "MoChA bar", "figure[data-chart=genome] rect.mca", 1),
           ("DEMO-003", "clinvar", "AnnotSV span", "figure[data-chart=genome] rect.sv", 1)]


def widget_counts(b, base):
    got = {}
    for person, cid, kind, sel, _ in WIDGETS:
        b.go("about:blank")  # a new document each time, so the wait never reads the page before
        b.wait("location.href === 'about:blank'")
        b.go(f"{base}/runs/{person}/#s/{cid}")
        b.wait("document.documentElement.dataset.pgpResults && !!document.getElementById(%s) && 'drawn'"
               % json.dumps("card-" + cid))
        if kind == "ypath":
            acts = json.loads(b.eval("JSON.stringify((function () { var x = new XMLHttpRequest(); x.open('GET', 'results.json', false); x.send(); return JSON.parse(x.responseText).findings.filter(function (f) { return f.level === 'act'; }).map(function (f) { return f.id; }); }()))"))
            for fid in acts:
                open_closed(b, fid)
        if kind in ("MoChA bar", "AnnotSV span"):
            b.wait("document.documentElement.dataset.pgpCharts")
            time.sleep(0.5)
        got[kind] = b.eval("document.querySelectorAll(%s).length" % json.dumps(sel))
    return got


def widget_problems(got):
    return [f"{kind}: {got.get(kind)} marks, want {n}" for _, _, kind, _, n in WIDGETS if got.get(kind) != n]


def widget_checks(b, base, out):
    """Every widget kind drawn with its marks on the partner's and the son's pages."""
    b.call("Network.enable", session=b.session)  # the control below changes person.js: never serve it from cache
    b.call("Network.setCacheDisabled", {"cacheDisabled": True}, b.session)
    got = widget_counts(b, base)
    bad = widget_problems(got)
    check("browser: every widget and genome-chart mark is drawn with its marks (" + ", ".join(k for _, _, k, _, _ in WIDGETS)
          + ")", not bad, bad)
    js = os.path.join(out, "assets", "person.js")
    with open(js) as f:
        keep = f.read()
    try:
        with open(js, "w") as f:
            f.write(keep.replace("add(text, widget(of[1]));", ""))
        control("person.js without its widget line", not widget_problems(widget_counts(b, base)))
    finally:
        with open(js, "w") as f:
            f.write(keep)


def family_depth(b, base):
    b.go("about:blank")
    b.wait("location.href === 'about:blank'")
    b.go(f"{base}/runs/{FAM}/")
    b.wait("document.documentElement.dataset.pgpStatus === 'ready' && "
           "(document.documentElement.dataset.pgpKids || '').indexOf('DEMO-003:known') >= 0")
    # The son's results to act on close his tables here: open each on purpose (one per click) to count them.
    for _ in range(4):
        b.eval("var v = document.querySelector('.veil[data-veil-child=\"DEMO-003\"] button'); v && v.click(); true")
        time.sleep(0.2)
    return b.wait("JSON.stringify(["
                  "document.querySelectorAll('table.rel-grid td:not(.self):not(.none)').length, "
                  "document.querySelectorAll('[data-painting=\"DEMO-003\"] .paint-seg').length, "
                  "document.querySelectorAll('details[data-child=\"DEMO-003\"] table.rows').length, "
                  "document.querySelectorAll('section.side-by-side table.rows thead th').length])")


def family_depth_problems(text, rows):
    got = json.loads(text) if text else [0, 0, 0, 0]
    want = [6, rows, 3, 2 + 3]  # 3 pairs both ways; the stretches; painting, de novo, pairs; label columns + 3 people
    return [] if got == want else [f"{got} != {want}"]


def family_depth_checks(b, base, out):
    with open(os.path.join(build_live_demo.EXTRA, FAM, "DEMO-003_segments.tsv")) as f:
        rows = sum(1 for line in f if line.strip()) - 1
    bad = family_depth_problems(family_depth(b, base), rows)
    check("browser: the family page draws the relatedness grid, the son's stretches row for row, his de novo and "
          "compound-het tables and the members side by side", not bad, bad)
    fj = os.path.join(out, "runs", FAM, "family.json")
    with open(fj) as f:
        keep = f.read()
    try:
        planted = json.loads(keep)
        planted["segments"] = {}
        with open(fj, "w") as f:
            json.dump(planted, f)
        control("family.json with the son's stretches emptied", not family_depth_problems(family_depth(b, base), rows))
    finally:
        with open(fj, "w") as f:
            f.write(keep)


# Two results to act on planted on the son's page: invented genes and variants in two CPSR rows.
VEIL = (("GENEQ", "c.77del"), ("GENER", "c.88dup"))


def veil_text(b, base):
    """The son's pages: [the person page's text before any click, after one click on the CPSR card's 'Show the result',
    the first finding's own page after its 'Show the result' click]. textContent reads closed details too."""
    b.go("about:blank")
    b.wait("location.href === 'about:blank'")
    b.go(f"{base}/runs/{OTHERS[1]}/")
    b.wait("document.documentElement.dataset.pgpResults && document.documentElement.dataset.pgpStatus === 'ready' "
           "&& document.querySelector('.veil') && 'drawn'")
    before = b.eval("document.body.textContent")
    # The card's own button: a card holding two results opens one per click, as each finding's button does.
    b.eval("document.querySelector('#card-cpsr .veil button').click(); true")
    time.sleep(0.5)
    after = b.wait("document.body.textContent")
    fid = b.eval("var a = document.querySelector('[data-finding] h3 a[href*=\"finding/cpsr-\"]'); a && a.getAttribute('href')")
    page = ""
    if fid:
        b.go(f"{base}/runs/{OTHERS[1]}/{fid}")
        b.wait("document.documentElement.dataset.pgpStatus === 'ready' && document.querySelector('.veil') && 'closed'")
        b.eval("document.querySelector('.veil button').click(); true")
        page = b.wait("!document.querySelector('.veil') && document.body.textContent")
    return before, after, page


def veil_problems(before, after, page):
    (g1, v1), (g2, v2) = VEIL
    out = [f"{w} shown before any click" for w in (g1, v1, g2, v2) if w in (before or "")]
    shown = [w for w in (g1, v1) if w in (after or "")] or [w for w in (g2, v2) if w in (after or "")]
    if len(shown) != 2:
        out.append(f"the first click did not open one result: {shown}")
    # The other result stays closed after one click, on the person page and on the opened finding's page.
    other = (g2, v2) if shown and shown[0] == g1 else (g1, v1)
    out += [f"{w} opened by the click on another result" for w in other if w in (after or "")]
    if not page:
        out.append("no finding page of a planted result")
    else:
        mine = [w for pair in VEIL for w in pair if w in page]
        if len(mine) != 2:
            out.append(f"the finding page opened shows {mine}, want only its own gene and variant")
    return out


def veil_checks(b, base, out):
    """On a child's page each result to act on stays closed everywhere, its own section's card too, until it
    is opened; opening one opens only that one, on the person page and on its own page."""
    rj = os.path.join(out, "runs", OTHERS[1], "results.json")
    js = os.path.join(out, "assets", "person.js")
    with open(rj) as f:
        keep = f.read()
    with open(js) as f:
        keep_js = f.read()
    made = []
    try:
        with open(os.path.join(os.path.dirname(out), "root", OTHERS[1], "summary.json")) as f:
            summary = json.load(f)
        summary.setdefault("sections", {})["cpsr"] = {"state": "ok", "data": {"call_floor_dp": 10, "low_confidence": 0, "plp": [
            {"gene": g, "variant": v, "final": "Pathogenic", "clinvar": "Pathogenic"} for g, v in VEIL]}}
        res = live_results.build(OTHERS[1], summary)
        cpsr = [c for ch in res["chapters"] for c in ch["sections"] if c["id"] == "cpsr"]
        check("browser: the planted CPSR card has its table, so the closed check below has rows to hide",
              bool(cpsr and cpsr[0]["table"] and cpsr[0]["section"] == 1), cpsr)
        with open(rj, "w") as f:
            json.dump(res, f)
        # Their finding pages, as build_live_demo.person_site writes one per finding: an existing page, its id swapped.
        fdirs = os.path.join(out, "runs", OTHERS[1], "finding")
        model = sorted(os.listdir(fdirs))[0]
        with open(os.path.join(fdirs, model, "index.html")) as f:
            tpl = f.read()
        for fd in res["findings"]:
            if fd["card"] == "cpsr":
                os.makedirs(os.path.join(fdirs, fd["id"]))
                made.append(os.path.join(fdirs, fd["id"]))
                with open(os.path.join(fdirs, fd["id"], "index.html"), "w") as f:
                    f.write(tpl.replace(f'data-finding="{model}"', f'data-finding="{fd["id"]}"'))
        bad = veil_problems(*veil_text(b, base))
        check("browser: on the son's page two planted results to act on name their genes and variants nowhere "
              "until 'Show the result' is clicked, and each click opens one, on the person page and the finding page",
              not bad, bad)
        with open(js, "w") as f:
            f.write(keep_js.replace("} else if (closedFs.length) {", "} else if (false) {"))
        control("person.js whose cards ignore the veil", not veil_problems(*veil_text(b, base)))
        with open(js, "w") as f:
            f.write(keep_js)
        fid = [fd["id"] for fd in res["findings"] if fd["card"] == "cpsr"][0]
        unknown_checks(b, base, out, fid)
    finally:
        with open(rj, "w") as f:
            f.write(keep)
        with open(js, "w") as f:
            f.write(keep_js)
        for d in made:
            shutil.rmtree(d, ignore_errors=True)


# What the son's page shows, to scan for the planted results: its text (SVG titles and closed folds included),
# every link and every title attribute.
SHOWN_JS = ("[document.body.textContent].concat(Array.prototype.map.call(document.querySelectorAll('[href],[title]'), "
            "function (e) { return (e.getAttribute('href') || '') + ' ' + (e.getAttribute('title') || ''); })).join(' ')")


def planted_words():
    """(the gene, the hit's words, the deletion's words) of make_demo_sample.PLANTED, as the page could write them."""
    import make_demo_sample
    pl = make_demo_sample.PLANTED
    c, p, _, _ = pl["hit"]
    sc, a, z = pl["sv"]
    return pl["gene"], [f"{c}:{p}", f"{p:,}"], [f"{sc}:{a}", f"{a:,}", f"{z:,}"]


def son_page(b, base):
    b.go("about:blank")
    b.wait("location.href === 'about:blank'")
    b.go(f"{base}/runs/{OTHERS[1]}/")
    b.wait("var d = document.documentElement.dataset, i = document.querySelector('input[type=search]'); "
           "d.pgpResults && d.pgpCharts && d.pgpStatus === 'ready' && document.querySelector('.veil') && i && !i.disabled "
           "&& 'ready'", 15)
    time.sleep(0.5)


def son_search(b, q):
    b.eval("var i = document.querySelector('input[type=search]'); i.value = %s; i.dispatchEvent(new Event('input')); true"
           % json.dumps(q))
    hits = b.wait("document.documentElement.dataset.pgpSearch", 5)
    return hits, b.eval("(document.querySelector('.search-results') || {}).textContent || ''") or ""


def planted_closed(b, base):
    """[problems] of the son's page before any click: the planted gene and positions shown anywhere, or found by
    the search."""
    gene, hit, sv = planted_words()
    son_page(b, base)
    shown = b.eval(SHOWN_JS) or ""
    out = [f"{w} shown while closed" for w in [gene] + hit + sv if w in shown]
    hits, text = son_search(b, gene)
    if hits != "0" or any(w in text for w in hit + sv):
        out.append(f"the search for {gene} while closed: {hits} hits, {text[:200]}")
    hits, text = son_search(b, gene[:3])
    if gene in text:
        out.append(f"the search for {gene[:3]} suggests {gene} while closed: {text[:200]}")
    return out


# The handle of a closed finding, as the page counts it (person.js closedOrder()): the closed findings in the order
# of the sections, then any other; a result to act on is closed until opened, a finding tagged with one while that
# result is, a finding tagged PANEL while any result to act on is, a finding tagged with no finding of the page always.
HANDLE_JS = r"""fetch('results.json', {cache: 'no-store'}).then(function (r) { return r.json(); }).then(function (r) {
  var fid = %s, byId = {}, tags = (r.veil && r.veil.findings) || {};
  r.findings.forEach(function (f) { byId[f.id] = f; });
  function open(id) { return !!document.getElementById('f-' + id); }
  var acts = r.findings.filter(function (f) { return f.level === 'act'; });
  var anyShut = acts.some(function (f) { return !open(f.id); });
  function shut(f) {
    if (f.level === 'act') { return !open(f.id); }
    var t = tags[f.id];
    if (t === undefined) { t = ['*panel*']; }
    t = (Array.isArray(t) ? t : t ? [t] : []).filter(function (x) { return x !== f.id; });
    return t.some(function (x) { return x === '*panel*' ? anyShut : !byId[x] || !open(x); });
  }
  var order = [];
  (r.journal || []).forEach(function (s) { (s.findings || []).forEach(function (id) {
    if (byId[id] && shut(byId[id]) && order.indexOf(id) < 0) { order.push(id); } }); });
  r.findings.forEach(function (f) { if (shut(f) && order.indexOf(f.id) < 0) { order.push(f.id); } });
  return order.indexOf(fid) + 1;
})"""


def open_closed(b, fid):
    """Opens the closed result FID as a reader does who knows its place among the page's closed findings: its handle
    (#f/closed-<n>), never its id, which the page does not use to find a closed card; then its "Show the result".
    False when the handle names no closed card."""
    r = b.call("Runtime.evaluate", {"expression": HANDLE_JS % json.dumps(fid), "awaitPromise": True,
                                     "returnByValue": True}, b.session)
    n = (r.get("result") or {}).get("value") or 0
    if not n:
        return False
    b.eval("location.hash = %s; true" % json.dumps(f"#f/closed-{n}"))
    if not b.wait("location.hash === %s && document.getElementById('f-closed-%d') && location.hash"
                  % (json.dumps(f"#f/closed-{n}"), n), 5):
        return False
    b.eval("document.getElementById('f-closed-%d').querySelector('.veil button').click(); true" % n)
    return True


def planted_problems(b, base, out_dir):
    """[problems] of the planted results on the son's page: closed, then opened one at a time, then reloaded."""
    gene, hit, sv = planted_words()
    out = planted_closed(b, base)
    with open(os.path.join(out_dir, "runs", OTHERS[1], "results.json")) as f:
        ids = " ".join(x["id"] for x in json.load(f)["findings"] if x["level"] == "act")
    sf = next((i for i in ids.split() if i.startswith("sf-")), None)
    import make_demo_sample
    start = str(make_demo_sample.PLANTED["sv"][1])
    del_ = next((i for i in ids.split() if i.startswith("sv-") and start in i), None)
    if not sf or not del_:
        return out + [f"the planted results are not among the son's results to act on: {ids}"]
    rows_js = ("Array.prototype.filter.call(document.querySelectorAll('#card-clinvar tbody tr'), function (r) { "
               "return r.textContent.indexOf(%s) >= 0; }).length" % json.dumps(gene))
    tick_js = ("Array.prototype.filter.call(document.querySelectorAll('figure[data-chart=genome] a[data-mark^=hit-] "
               "title'), function (t) { return t.textContent.indexOf(%s) >= 0; }).length" % json.dumps(gene))
    ranked_js = ("Array.prototype.filter.call(document.querySelectorAll('[data-finding^=exomiser-]'), function (e) { "
                 "return e.textContent.indexOf(%s) >= 0; }).length" % json.dumps(gene))
    if not open_closed(b, sf):
        return out + [f"a link to {sf} lands on no closed card"]
    time.sleep(0.8)
    shown = b.eval(SHOWN_JS) or ""
    # The planted hit lies inside the planted deletion, in its gene: its ClinVar row, its tick and the symptom
    # ranking of the gene name both results, so they stay closed until the deletion is opened too.
    rows, tick, ranked = b.eval(rows_js), b.eval(tick_js), b.eval(ranked_js)
    if not all(w in shown for w in [gene, hit[0]]):
        out.append(f"opening the secondary finding did not bring back its own card ({gene}, {hit[0]})")
    if rows or tick or ranked:
        out.append(f"opening the secondary finding alone brought back places that name the closed deletion too: "
                   f"{rows} ClinVar rows, {tick} ticks, {ranked} rankings")
    out += [f"{w} shown after opening only the secondary finding" for w in sv if w in shown]
    # A search hit quotes the summary's rows (prose): it stays closed while any result is, the same whichever.
    hits, _ = son_search(b, gene)
    if hits not in ("0", None):
        out.append(f"the search for {gene} finds {hits} hits while the deletion is still closed")
    if not open_closed(b, del_):
        return out + [f"a link to {del_} lands on no closed card"]
    time.sleep(0.8)
    shown = b.eval(SHOWN_JS) or ""
    spans = b.eval("document.querySelectorAll('figure[data-chart=genome] a[data-mark^=sv-]').length")
    if sv[0] not in shown or spans != 2:
        out.append(f"opening the deletion did not bring back its span on the genome chart: {spans} spans")
    hits, _ = son_search(b, gene)
    if hits in (None, "0"):
        out.append(f"the search for {gene} finds nothing once every result is opened")
    rows, tick, ranked = b.eval(rows_js), b.eval(tick_js), b.eval(ranked_js)
    if rows != 1 or not tick or not ranked:
        out.append(f"opening both results did not bring back the hit's ClinVar row, tick and ranking: {rows} rows, "
                   f"{tick} ticks, {ranked} rankings")
    out += [f"after a reload: {x}" for x in planted_closed(b, base)]
    return out


def planted_pages(b, base, out):
    """[problems] of the son's other pages before any click: the home page and every finding page of his (an open
    finding's rail and related list included) showing the planted gene or a position; then the page of the
    finding that names the closed result (the symptom ranking), opened on purpose, shows the gene."""
    gene, hit, sv = planted_words()
    with open(os.path.join(out, "runs", OTHERS[1], "results.json")) as f:
        res = json.load(f)
    ids = [f["id"] for f in res.get("findings") or []]
    shown_js = "document.title + ' ' + " + SHOWN_JS
    out_ = []
    b.go("about:blank")
    b.wait("location.href === 'about:blank'")
    b.go(f"{base}/")
    b.wait("document.documentElement.dataset.pgpStatus === 'ready' && document.documentElement.dataset.pgpPeople", 15)
    shown = b.eval(shown_js) or ""
    out_ += [f"{w} on the home page" for w in [gene] + hit + sv if w in shown]
    for fid in ids:
        b.go("about:blank")
        b.wait("location.href === 'about:blank'")
        b.go(f"{base}/runs/{OTHERS[1]}/finding/{fid}/")
        b.wait("document.documentElement.dataset.pgpStatus === 'ready' && 'drawn'", 10)
        shown = b.eval(shown_js) or ""
        out_ += [f"{w} on the finding page {fid}" for w in [gene] + hit + sv if w in shown]
    named = [k for k, t in ((res.get("veil") or {}).get("findings") or {}).items() if t and k.startswith("exomiser-")]
    if len(ids) < 10 or not named:
        return out_ + [f"the son's findings are not the planted ones: {ids}, ranking {named}"]
    b.go("about:blank")
    b.wait("location.href === 'about:blank'")
    b.go(f"{base}/runs/{OTHERS[1]}/finding/{named[0]}/")
    if b.wait("document.documentElement.dataset.pgpStatus === 'ready' && document.querySelector('.veil') && 'closed'", 10):
        b.eval("document.querySelector('.veil button').click(); true")
        shown = b.wait("!document.querySelector('.veil') && document.body.textContent", 5) or ""
        if gene not in shown:
            out_.append(f"opening the page of {named[0]} did not show {gene}")
        out_ += [f"{w} shown on {named[0]} after opening it" for w in sv if w in shown]
    return out_


def planted_checks(b, base, out):
    """The son's planted secondary finding and class 5 deletion stay closed everywhere on his page: the ClinVar
    card's table and viewer link, the genome chart and the review-star chart (wherever it sits), their titles,
    the search; each opens on purpose with its own rows only, and a reload closes them again."""
    sd = os.path.join(out, "runs", OTHERS[1])
    rj = os.path.join(sd, "results.json")
    pj, cj, sj, fj = (os.path.join(out, "assets", x) for x in ("person.js", os.path.join("plugins", "charts.js"),
                                                                os.path.join("plugins", "search.js"), "finding.js"))
    keep = {}
    for p in (rj, pj, cj, sj, fj):
        with open(p) as f:
            keep[p] = f.read()

    def write(p, text):
        with open(p, "w") as f:
            f.write(text)
    try:
        bad = planted_problems(b, base, out)
        check("browser: the son's planted secondary finding and class 5 deletion name their gene and positions nowhere "
              "(cards, charts, titles, links, search) until each is opened, each opens only its own rows and marks (a "
              "place that names both once both are), and a reload closes them", not bad, bad)
        res = json.loads(keep[rj])
        res["veil"] = {}
        write(rj, json.dumps(res))
        bad = planted_closed(b, base)
        check("browser: results.json without veil metadata still keeps secondary places closed", not bad, bad)
        res = json.loads(keep[rj])
        def unsafe_tags(value):
            if isinstance(value, dict):
                return {k: unsafe_tags(v) for k, v in value.items()}
            if isinstance(value, list):
                return None if all(isinstance(v, str) for v in value) else [unsafe_tags(v) for v in value]
            return value
        for kind in ("rows", "viewer", "lines", "facts", "headline", "marks", "findings"):
            if kind in res["veil"]:
                res["veil"][kind] = unsafe_tags(res["veil"][kind])
        write(rj, json.dumps(res))
        control("results.json falsely exempts tagged secondary places", not planted_closed(b, base))
        write(rj, keep[rj])
        for p, a, z, desc in (
                (pj, "return !placeClosed(tags[i]);", "return true;", "person.js whose ClinVar rows ignore their tags"),
                (pj, "return !!(t && tagClosed(t));", "return false;", "person.js whose findings ignore their tags"),
                (cj, "return ask ? (sp.marks", "return false ? (sp.marks", "charts.js that keeps every mark"),
                (sj, "it.veil || null", "null", "search.js that ignores a hit's tag")):
            if a not in keep[p]:
                check(f"control: {desc} (the planted line is in {os.path.basename(p)})", False, a)
                continue
            write(p, keep[p].replace(a, z))
            control(desc, not planted_closed(b, base))
            write(p, keep[p])
        bad = planted_pages(b, base, out)
        check("browser: the home page and every finding page of the son name the planted gene and positions nowhere "
              "(the rail's next link and the related list included) until opened, and the symptom ranking's page "
              "opened on purpose shows it", not bad, bad)
        a = "return (tagOf(x) || []).some(function (t) { return t === PANEL ? !opened[x.id] && actShut() : !opened[t]; });"
        if a not in keep[fj]:
            check("control: finding.js whose findings ignore their tags (the planted line is in finding.js)", False, a)
        else:
            write(fj, keep[fj].replace(a, "return false;"))
            control("finding.js whose findings ignore their tags", not planted_pages(b, base, out))
            write(fj, keep[fj])
    finally:
        for p, text in keep.items():
            write(p, text)


# The family page: what it shows of the son while a result of his to act on is closed.
# Every text of the page and every attribute of every element: a gene, a position or a grandparent hidden in a
# title, a link, a label, a data attribute or a class is as shown as one in a cell.
FAMILY_SHOWN_JS = ("[document.body.textContent].concat(Array.prototype.map.call(document.querySelectorAll('*'), function (e) { "
                   "return Array.prototype.map.call(e.attributes, function (a) { return a.value; }).join(' '); })).join(' ')")
GRANDPARENT = "DEMO-GF"
PAIR_B = "150962400"  # the second variant of the fixture's compound-het pair in KCNH2
SON = OTHERS[1]
# The places that close for the son while a result of his is closed, one box each, in page order.
SON_BOXES = ["paint", "denovo", "het", "side"]


def no_cache(b):
    """Chrome keeps a script it fetched a moment ago without asking again: the controls below swap a script, so
    every request goes to the server from here on."""
    b.call("Network.enable", {}, b.session)
    b.call("Network.setCacheDisabled", {"cacheDisabled": True}, b.session)


def family_page(b, base):
    """The family page once it knows the son's results (or after 10 s, for a page that never learns them)."""
    b.go("about:blank")
    b.wait("location.href === 'about:blank'")
    b.go(f"{base}/runs/{FAM}/")
    b.wait("document.documentElement.dataset.pgpStatus === 'ready' && "
           "/DEMO-003:(known|failed)/.test(document.documentElement.dataset.pgpKids || '') && 'ready'", 10)
    time.sleep(0.3)


def boxes_of(b, child=SON):
    """The places closed for CHILD (data-veil-box of each of its closed boxes), in page order."""
    got = b.eval("JSON.stringify(Array.prototype.map.call(document.querySelectorAll('.veil[data-veil-child=\"" + child +
                 "\"]'), function (x) { return x.getAttribute('data-veil-box'); }))")
    return json.loads(got) if got else None


def open_box(b, box, child=SON):
    """Click the button of CHILD's closed box BOX; whether there was one."""
    ok = b.eval("var v = document.querySelector('.veil[data-veil-child=\"" + child + "\"][data-veil-box=\"" + box + "\"] button'); "
                "if (v) { v.click(); } !!v")
    time.sleep(0.4)
    return ok


def son_tables(b, child=SON):
    """[rows of each table] in CHILD's fold of the family page."""
    rows = b.eval("JSON.stringify(Array.prototype.map.call(document.querySelectorAll('details[data-child=\"" + child + "\"] "
                  "table.rows'), function (t) { return t.querySelectorAll('tbody tr').length; }))")
    return json.loads(rows) if rows else None


def side_rows(b):
    """The side by side as [[group, label, value of each member]], the group carried down."""
    got = b.eval("JSON.stringify(Array.prototype.map.call(document.querySelectorAll('#ch-side-by-side table.rows tbody tr'), "
                 "function (r) { return Array.prototype.map.call(r.cells, function (c) { return c.textContent; }); }))")
    rows, group = [], ""
    for r in json.loads(got) if got else []:
        group = r[0] or group
        rows.append([group] + r[1:])
    return rows


def side_problems_of(b, closed):
    """[problems] of the side by side: the son's medicines and HLA cells all "closed" when CLOSED, none otherwise; a
    parent's cell never closed, the haplogroups and blood groups never closed."""
    rows = side_rows(b)
    members = b.eval("JSON.stringify(Array.prototype.map.call(document.querySelectorAll('#ch-side-by-side table.rows thead th'), "
                     "function (x) { return x.textContent; }).slice(2))")
    members = json.loads(members) if members else []
    if SON not in members or not rows:
        return [f"no side by side with the son: {members}"]
    out = []
    for r in rows:
        for m, v in zip(members, r[2:]):
            want = closed and m == SON and r[0] in ("Medicines", "HLA")
            if (v == "closed") != want:
                out.append(f"{m}'s {r[1]} reads {v!r}")
    return out[:6]


def family_closed(b, base):
    """[problems] of the family page before any click: the son's planted gene, its positions or the grandparent of
    the stretch over it named anywhere, or a place of his open: the painting's grandparents, his de novo and
    compound-het rows and his medicines and HLA cells close, one box each."""
    gene, hit, sv = planted_words()
    family_page(b, base)
    shown = b.eval(FAMILY_SHOWN_JS) or ""
    out = [f"{w} shown on the family page while closed" for w in [gene, GRANDPARENT, PAIR_B] + hit + sv[:1] if w in shown]
    boxes, tables = boxes_of(b), son_tables(b)
    segs = b.eval("document.querySelectorAll('[data-painting=\"DEMO-003\"] .paint-seg').length")
    marked = b.eval("document.querySelectorAll('[data-painting=\"DEMO-003\"] .paint-seg.closed, "
                    "[data-painting=\"DEMO-003\"] .paint-seg.other').length")
    # Only the painting's table of stretches stays, every stretch drawn and listed from its parent, none marked.
    if boxes != SON_BOXES or not tables or len(tables) != 1 or tables[0] != segs or marked:
        out.append(f"closed places {boxes}, tables {tables}, {segs} stretches drawn, {marked} marked")
    out += side_problems_of(b, True)
    return out


def family_problems(b, base):
    """[problems] of the family page: closed; one click on the compound-het box opens one of the son's two planted
    results and leaves every place closed; a click on the de novo box opens the other and every place; a reload
    closes them again."""
    gene, hit, _ = planted_words()
    out = family_closed(b, base)
    if out:
        return out
    if not open_box(b, "het"):
        return ["no compound-het box"]
    shown = b.eval(FAMILY_SHOWN_JS) or ""
    if any(w in shown for w in [gene, GRANDPARENT, PAIR_B] + hit) or boxes_of(b) != SON_BOXES:
        out.append(f"one click opened a place while the other result is still closed: {boxes_of(b)}")
    open_box(b, "denovo")
    shown = b.eval(FAMILY_SHOWN_JS) or ""
    tables = son_tables(b)
    if not all(w in shown for w in [gene, GRANDPARENT, hit[0], PAIR_B]) or not tables or tables[1:] != [3, 3] or boxes_of(b):
        out.append(f"opening both results did not bring back every place: tables {tables}, closed {boxes_of(b)}")
    out += side_problems_of(b, False)
    out += [f"after a reload: {x}" for x in family_closed(b, base)]
    return out


def family_view(b, base, child=SON):
    """The family page's HTML before any click, or None when the page did not settle on knowing CHILD's results."""
    family_page(b, base)
    kids = b.eval("document.documentElement.dataset.pgpKids || ''") or ""
    if not re.search(r"\b" + re.escape(child) + r":known\b", kids):
        return None
    return b.eval("document.getElementById('family').innerHTML") or ""


# Closed results of different kinds and places, as veil.spans holds them.
CLOSED_CASES = (("RYR1, a medicines gene", [["sf-ryr1", "RYR1", "chr19", 38440000, 38440000]]),
                ("HLA-B", [["sf-hlab", "HLA-B", "chr6", 31353000, 31357000]]),
                ("a gene alone", [["two-hits-g0999", "G0999", None, None, None]]),
                ("a position under a de novo row", [["sf-pos", None, "chr3", 48211904, 48211904]]),
                ("a copy-number span in SMN1", [["sv-chr5-70925000-70953000", "SMN1", "chr5", 70925000, 70953000]]),
                ("a CHIP call in TP53", [["chip-tp53", "TP53", "chr17", 7675088, 7675088]]))


def invariance_problems(b, base, rj, res, fj, fams, cases=None, child=SON):
    """[problems]: the family page before any click differs between closed results of different genes and places
    (FAMS: name -> family.json text; CASES: name -> spans set in RJ, CHILD's results, from RES, its findings to act
    on renamed to match, None for RES as it is), so it would tell which is closed; or it never knew CHILD's results."""
    views = {}
    for name, spans in cases or (("the planted pair in KCNH2", None),) + CLOSED_CASES:
        r = copy.deepcopy(res)
        if spans is not None:
            act = [x for x in r.get("findings", []) if x.get("level") == "act"]
            r["veil"] = dict(r.get("veil") or {}, spans=spans)
            if act:
                r["findings"] = [x for x in r["findings"] if x.get("level") != "act"] + [dict(act[0], id=s_[0]) for s_ in spans]
        with open(rj, "w") as f:
            json.dump(r, f)
        for fname, text in fams.items():
            with open(fj, "w") as f:
                f.write(text)
            views[f"{name}, {fname}"] = family_view(b, base, child)
    with open(rj, "w") as f:
        json.dump(res, f)
    out = [f"the page never knew {child}'s results with {k!r}" for k, v in views.items() if v is None]
    first = next(iter(views))
    return out + [f"the page differs between {first!r} and {k!r}" for k, v in views.items() if v != views[first]]


def family_unknown(b, base, opens=True):
    """[problems] of the family page while it cannot read the son's results: a place of his open or a word of his
    planted results shown; then (OPENS) each box's own button, from a fresh page each, opens its own place and no
    other."""
    gene, hit, sv = planted_words()
    family_page(b, base)
    # A page still waiting for his results keeps everything closed anyway: the check is on the state it settles in.
    kids = b.eval("document.documentElement.dataset.pgpKids || ''") or ""
    shown = b.eval(FAMILY_SHOWN_JS) or ""
    out = [] if re.search(r"\bDEMO-003:(known|failed)\b", kids) else [f"the son's results never settled: {kids!r}"]
    out += [f"{w} shown with the son's results unreadable" for w in [gene, GRANDPARENT, PAIR_B] + hit + sv[:1] if w in shown]
    if boxes_of(b) != SON_BOXES or son_tables(b) is None or len(son_tables(b)) != 1:
        out.append(f"closed places {boxes_of(b)}, tables {son_tables(b)}")
    words = b.eval("JSON.stringify(Array.prototype.map.call(document.querySelectorAll('.veil[data-veil-child=\"DEMO-003\"] button'), "
                   "function (x) { return x.textContent; }))")
    if json.loads(words or "[]") != ["Show the grandparents", "Show the rows", "Show the rows", "Show the cells"]:
        out.append(f"the buttons read {words}")
    out += side_problems_of(b, True)
    if not opens or out:
        return out
    # What each place shows once open: the painting its grandparent, the de novo rows the planted hit (not the
    # pair's second variant), the compound-het rows the pair's second variant, the side by side his cells.
    for box in SON_BOXES:
        family_page(b, base)
        if not open_box(b, box):
            out.append(f"no {box} box")
            continue
        shown = b.eval(FAMILY_SHOWN_JS) or ""
        rest, tables = [x for x in SON_BOXES if x != box], son_tables(b) or []
        want = {"paint": GRANDPARENT in shown and PAIR_B not in shown and len(tables) == 1,
                "denovo": GRANDPARENT not in shown and PAIR_B not in shown and tables[1:] == [3],
                "het": GRANDPARENT not in shown and PAIR_B in shown and tables[1:] == [3],
                "side": GRANDPARENT not in shown and PAIR_B not in shown and len(tables) == 1}[box]
        if boxes_of(b) != rest or not want:
            out.append(f"the {box} box opened another place, or not its own: closed {boxes_of(b)}, tables {tables}")
        out += [f"after the {box} box: {x}" for x in side_problems_of(b, box != "side")]
    return out


def unknown_map_problems(b, base):
    """[problems] of the family page whose son has his map as an image, with his results unreadable: 'Show the
    rows' of a table opening the image, 'Show the map' not opening it, or (from a fresh page) 'Show the map' opening
    a table."""
    family_page(b, base)
    img = "!!document.querySelector('details[data-child=\"DEMO-003\"] img.segments')"
    out = []
    if not open_box(b, "denovo"):
        out.append("no de novo box")
    if b.eval(img):
        out.append("'Show the rows' opened the map too")
    if not open_box(b, "map") or not b.eval(img):
        out.append("'Show the map' did not open the map")
    family_page(b, base)
    if not open_box(b, "map") or not b.eval(img) or boxes_of(b) != ["denovo", "het", "side"] or son_tables(b):
        out.append(f"'Show the map' opened another place, or not the map: closed {boxes_of(b)}, tables {son_tables(b)}")
    return out


def map_problems(b, base):
    """[problems] of the family page whose son has his map as an image and two closed results: the image shown
    before both are opened from its box (one per click), or not shown after."""
    family_page(b, base)
    img = "!!document.querySelector('details[data-child=\"DEMO-003\"] img.segments')"
    out = []
    if b.eval(img) or boxes_of(b) != ["map", "denovo", "het", "side"]:
        out.append(f"the map shown while a result is closed, or a place open: {boxes_of(b)}")
    open_box(b, "map")
    if b.eval(img) or boxes_of(b) != ["map", "denovo", "het", "side"]:
        out.append(f"one click opened a place while the other result is closed: {boxes_of(b)}")
    open_box(b, "map")
    if not b.eval(img) or boxes_of(b):
        out.append(f"the map not back once both results are opened: closed {boxes_of(b)}")
    return out


def two_children_problems(b, base):
    """[problems] of a family with a second child, DEMO-004, whose run has the son's results under the same ids:
    opening every result of the son leaves each of DEMO-004's places closed."""
    family_page(b, base)
    out = []
    if boxes_of(b, "DEMO-004") != SON_BOXES:
        out.append(f"DEMO-004's closed places before: {boxes_of(b, 'DEMO-004')}")
    for _ in range(2):
        open_box(b, "denovo")
    if boxes_of(b):
        out.append(f"the son's places not opened: {boxes_of(b)}")
    other = b.eval("document.querySelector('details[data-child=\"DEMO-004\"]').textContent") or ""
    if boxes_of(b, "DEMO-004") != SON_BOXES or GRANDPARENT in other or planted_words()[1][0] in other:
        out.append(f"opening the son's results opened DEMO-004's: {boxes_of(b, 'DEMO-004')}")
    return out


def couple_problems(b, base):
    """[problems] of a family where DEMO-001, a partner of the couple screen, is also a child (of an unsequenced
    father) with a closed result in PRF1, a gene of the screen: the screen names it (its line's counts, its
    table, its boxes, its report link, the lead) before the result is opened from its box, or not after."""
    family_page(b, base)
    out = []
    sect = b.eval("var s = document.getElementById('ch-couples'); s ? s.textContent + ' ' + s.innerHTML : ''") or ""
    lead = b.eval("(document.getElementById('family-lead') || {}).textContent || ''") or ""
    links = b.eval("Array.prototype.map.call(document.querySelectorAll('a[href]'), function (a) { return a.getAttribute('href'); }).join(' ')") or ""
    bx = boxes_of(b, "DEMO-001") or []
    if ("PRF1" in sect or "70601223" in sect or "CPT2" in sect or "couple_DEMO" in links
            or "carrier screen stays closed in chapter III" not in lead):
        out.append("the screen of a couple whose partner has a closed result named before it is opened")
    if bx != ["couple:DEMO-001 DEMO-002", "side"]:
        out.append(f"DEMO-001's closed places: {bx}")
    if not open_box(b, "couple:DEMO-001 DEMO-002", "DEMO-001"):
        return out + ["no closed box under the screen"]
    sect = b.eval("document.getElementById('ch-couples').textContent") or ""
    links = b.eval("Array.prototype.map.call(document.querySelectorAll('a[href]'), function (a) { return a.getAttribute('href'); }).join(' ')") or ""
    if "PRF1" not in sect or "couple_DEMO" not in links or boxes_of(b, "DEMO-001"):
        out.append(f"the screen not back once the result is opened: closed {boxes_of(b, 'DEMO-001')}")
    return out


def nothing_closed_problems(b, base):
    """[problems] of the family page whose son has nothing to act on (no veil): a box of his, or a place of his
    closed: the painting's grandparent, his de novo and compound-het rows, his cells."""
    family_page(b, base)
    kids = b.eval("document.documentElement.dataset.pgpKids || ''") or ""
    shown = b.eval(FAMILY_SHOWN_JS) or ""
    tables = son_tables(b) or []
    out = [] if re.search(r"\bDEMO-003:known\b", kids) else [f"the son's results never known: {kids!r}"]
    if boxes_of(b) or GRANDPARENT not in shown or len(tables) != 3 or tables[1:] != [3, 3]:
        out.append(f"a place closed with nothing to act on: closed {boxes_of(b)}, tables {tables}")
    return out + side_problems_of(b, False)


def new_run_problems(b, base, fj, before, after):
    """[problems]: with family.json BEFORE, the son's results opened from the page; then family.json AFTER (his run
    in another folder) is fetched by the page itself, and his places must close again. The page asks every 30 s:
    here every 0.5 s."""
    b.call("Page.enable", {}, b.session)  # the script below runs only with the Page domain on
    fast = b.call("Page.addScriptToEvaluateOnNewDocument", {"source": "(function () { var s = window.setInterval; "
                  "window.setInterval = function (f, t) { return s(f, t === 30000 ? 500 : t); }; }());"}, b.session)
    try:
        with open(fj, "w") as f:
            f.write(before)
        family_page(b, base)
        out = []
        for _ in range(2):
            open_box(b, "denovo")
        if boxes_of(b):
            out.append(f"the son's places not opened: {boxes_of(b)}")
        with open(fj, "w") as f:
            f.write(after)
        moved = b.wait("var r = Array.prototype.some.call(document.querySelectorAll('a[href]'), function (a) { "
                       "return a.getAttribute('href').indexOf('-again/') >= 0; }); r && "
                       "/DEMO-003:known/.test(document.documentElement.dataset.pgpKids || '') && 'moved'", 10)
        time.sleep(1.5)
        if moved != "moved":
            out.append("the page never fetched the moved run")
        elif boxes_of(b) != SON_BOXES:
            out.append(f"opened results stayed open for the son's new run: closed {boxes_of(b)}")
        return out
    finally:
        b.call("Page.removeScriptToEvaluateOnNewDocument", {"identifier": (fast or {}).get("identifier", "")}, b.session)
        b.call("Page.disable", {}, b.session)


def gene_cells(b, child):
    """The values of CHILD's medicines and HLA cells in the side by side."""
    rows = side_rows(b)
    members = b.eval("JSON.stringify(Array.prototype.map.call(document.querySelectorAll('#ch-side-by-side table.rows thead th'), "
                     "function (x) { return x.textContent; }).slice(2))")
    members = json.loads(members) if members else []
    if child not in members:
        return []
    i = members.index(child)
    return [r[2 + i] for r in rows if r[0] in ("Medicines", "HLA")]


def unknown_couple_problems(b, base):
    """[problems] of the couple fixture while the page cannot read DEMO-001's results: from a fresh page each, 'Show
    the screen' opens the screen and leaves DEMO-001's cells closed, and 'Show the cells' opens the cells and leaves
    the screen closed."""
    out, box = [], "couple:DEMO-001 DEMO-002"
    screen = "!!document.querySelector('#ch-couples table.rows')"
    for which in (box, "side"):
        family_page(b, base)
        if boxes_of(b, "DEMO-001") != [box, "side"] or b.eval(screen) or set(gene_cells(b, "DEMO-001")) != {"closed"}:
            out.append(f"DEMO-001's places not all closed: {boxes_of(b, 'DEMO-001')}")
            continue
        if not open_box(b, which, "DEMO-001"):
            out.append(f"no {which} box")
            continue
        cells = set(gene_cells(b, "DEMO-001"))
        rest = [x for x in (box, "side") if x != which]
        ok = (b.eval(screen) and cells == {"closed"}) if which == box else (not b.eval(screen) and cells and "closed" not in cells)
        if boxes_of(b, "DEMO-001") != rest or not ok:
            out.append(f"the {which} box opened another place, or not its own: closed {boxes_of(b, 'DEMO-001')}, cells {sorted(cells)[:3]}")
    return out


def family_veil_checks(b, base, out):
    """On the family page, while a child has a result to act on still closed, or the page cannot read the child's
    results, every place that could carry one closes for that child, the same whichever result it is; each known
    result opens with its own click, and each place opens on its own while the results cannot be read."""
    no_cache(b)
    import make_demo_sample
    fj = os.path.join(out, "runs", FAM, "family.json")
    rj = os.path.join(out, "runs", SON, "results.json")
    rj1 = os.path.join(out, "runs", S, "results.json")
    fp = os.path.join(out, "assets", "family-page.js")
    keep = {}
    for p in (fj, fp, rj, rj1):
        with open(p) as f:
            keep[p] = f.read()

    def write(p, text):
        with open(p, "w") as f:
            f.write(text)

    def planted_control(desc, a, z, problems):
        """Run PROBLEMS against family-page.js with A replaced by Z: it must find a problem."""
        if a not in keep[fp]:
            check(f"control: {desc} (the planted line is in family-page.js)", False, a)
            return
        write(fp, keep[fp].replace(a, z))
        try:
            control(desc, not problems())
        finally:
            write(fp, keep[fp])
    try:
        # A de novo candidate at the planted hit, and the father's copy over it painted from a grandparent.
        fam = json.loads(keep[fj])
        c, pos, ref, alt = make_demo_sample.PLANTED["hit"]
        fam["de_novo"].append(dict(fam["de_novo"][0], chrom=c, pos=str(pos), ref=ref, alt=alt))
        segs = fam["segments"][SON]
        over = [s_ for s_ in segs if s_["chrom"] == c and s_["side"] == "paternal" and int(s_["start"]) <= pos <= int(s_["end"])]
        away = [s_ for s_ in segs if s_["chrom"] != c and s_["side"] == "paternal"]
        if len(over) != 1 or not away:
            check("browser: the son's stretches hold one paternal stretch over the planted hit, and others", False, over)
            return
        over[0]["ancestor"] = GRANDPARENT
        planted = json.dumps(fam)
        over[0]["ancestor"] = "."
        away[0]["ancestor"] = GRANDPARENT
        elsewhere = json.dumps(fam)
        write(fj, planted)
        bad = family_problems(b, base)
        check("browser: on the family page, while the son has a result to act on closed, his painting's grandparents, "
              "his de novo and compound-het rows and his medicines and HLA cells close, one box each, naming no gene, "
              "position or grandparent; each click opens one result, every place comes back once both are opened, and "
              "a reload closes them", not bad, bad)
        for desc, a, z in (
                ("de novo rows that ignore the veil", "var dshut = shutOf(c.child, 'denovo');", "var dshut = [];"),
                ("compound-het rows that ignore the veil", "var hshut = shutOf(c.child, 'het');", "var hshut = [];"),
                ("a painting that ignores the veil", "var shut = anc ? shutOf(child, 'paint') : [];", "var shut = [];"),
                ("a closed painting that keeps each stretch's grandparent", "ancestor: '.', basis: 'parent' };",
                 "ancestor: r.ancestor, basis: r.basis };"),
                ("a side by side whose cells ignore the veil", "closedCells[m] = true; return 'closed';", "closedCells[m] = true; return v;")):
            planted_control(desc, a, z, lambda: family_closed(b, base))
        planted_control("a box that opens every result with one click",
                        "function () { opened[openKey(child, shut[0])] = true; redraw(); });",
                        "function () { shut.forEach(function (t) { opened[openKey(child, t)] = true; }); redraw(); });",
                        lambda: family_problems(b, base))

        # The same page whichever result is closed and wherever it lies.
        res = json.loads(keep[rj])
        fams = {"a grandparent over the hit": planted, "a grandparent elsewhere": elsewhere}
        bad = invariance_problems(b, base, rj, res, fj, fams)
        check("browser: the family page is the same, byte for byte, whichever of the son's results is closed (a pair "
              "in KCNH2, RYR1, HLA-B, a gene alone, a position, a copy-number span in SMN1, a CHIP call in TP53) and "
              "wherever a grandparent's stretch lies, each view settled on knowing his results", not bad, bad)
        for desc, a, z in (
                ("de novo rows closed only when a result lies under one", "var dshut = shutOf(c.child, 'denovo');",
                 "var dshut = shutOf(c.child, 'denovo').filter(function (id) { return dn.some(function (r) { return "
                 "kids[c.child].spans.some(function (s) { return s[0] === id && s[2] === r.chrom && +r.pos >= s[3] && "
                 "+r.pos <= s[4]; }); }); });"),
                ("a painting closed only when a grandparent's stretch lies over a closed result",
                 "var anc = all.some(function (r) { return r.ancestor && r.ancestor !== '.'; });",
                 "var anc = all.some(function (r) { return r.ancestor && r.ancestor !== '.' && ((kids[child] || {}).spans || [])"
                 ".some(function (s) { return !s[2] || (s[2] === r.chrom && s[3] <= +r.end && s[4] >= +r.start); }); });"),
                ("a side by side that closes only the closed gene's cells",
                 "if (GENE_GROUPS[r.group] && shutBy[m] && shutBy[m].length) {",
                 "if (GENE_GROUPS[r.group] && shutBy[m] && shutBy[m].length && (kids[m].spans || []).some(function (s) { "
                 "return s[1] === r.label; })) {"),
                ("a side by side whose cells stay open for a span wider than 10 kb (a copy-number region)",
                 "if (GENE_GROUPS[r.group] && shutBy[m] && shutBy[m].length) {",
                 "if (GENE_GROUPS[r.group] && shutBy[m] && shutBy[m].length && (kids[m].spans || []).every(function (s) { "
                 "return !(s[4] - s[3] > 10000); })) {"),
                ("a page that never reads the son's results (every view the same, all closed as unreadable)",
                 "U.getJSON(url + 'results.json')", "U.getJSON(url + 'missing.json')")):
            planted_control(desc, a, z, lambda: invariance_problems(b, base, rj, res, fj, fams))
        write(fj, planted)

        # The page cannot read the son's results: no answer, an unreadable answer, act findings without veil.spans.
        os.replace(rj, rj + ".away")
        try:
            bad = family_unknown(b, base)
            check("browser: with the son's results.json unavailable every place of his stays closed, and each box's "
                  "own button (the painting's, the de novo, the compound-het, the side by side's) opens its own place "
                  "and no other", not bad, bad)
            planted_control("one opening shared by every box while the results cannot be read",
                            "return opened[openKey(child, 'box:' + box)] ? [] : [UNKNOWN];",
                            "return opened[openKey(child, 'box:denovo')] ? [] : [UNKNOWN];",
                            lambda: family_unknown(b, base))
            planted_control("'Show the grandparents' that opens the side by side while the results cannot be read",
                            "return opened[openKey(child, 'box:' + box)] ? [] : [UNKNOWN];",
                            "return opened[openKey(child, box === 'side' ? 'box:paint' : 'box:' + box)] ? [] : [UNKNOWN];",
                            lambda: family_unknown(b, base))
        finally:
            os.replace(rj + ".away", rj)
        write(rj, json.dumps({"schema_version": 1, "sample": SON, "chapters": [],
                              "note": "The results could not be read (KeyError)."}))
        bad = family_unknown(b, base, opens=False)
        check("browser: with the son's results answered as unreadable (a note, no findings) every place of his stays "
              "closed", not bad, bad)
        write(rj, keep[rj][:5000])
        bad = family_unknown(b, base, opens=False)
        check("browser: with the son's results.json cut short (a request that fails) every place of his stays closed",
              not bad, bad)
        planted_control("a failed request read as nothing closed", "cur.status = 'failed';",
                        "cur.status = 'known'; cur.spans = null; cur.acts = [];", lambda: family_unknown(b, base, opens=False))
        for desc, spans in (("[null]", [None]), ("[] while he has findings to act on", [])):
            bad_spans = json.loads(keep[rj])
            bad_spans["veil"]["spans"] = spans
            write(rj, json.dumps(bad_spans))
            bad = family_unknown(b, base, opens=False)
            check(f"browser: with the son's veil.spans malformed ({desc}) every place of his stays closed", not bad, bad)
            planted_control(f"malformed veil.spans ({desc}) read as given", "spans: spansOk ? v.spans : null",
                            "spans: Array.isArray(v.spans) ? v.spans : null", lambda: family_unknown(b, base, opens=False))
        clean = json.loads(keep[rj])
        clean["findings"] = [x for x in clean["findings"] if x.get("level") != "act"]
        clean.pop("veil", None)
        write(rj, json.dumps(clean))
        bad = nothing_closed_problems(b, base)
        check("browser: with nothing of the son's to act on (no veil) every place of his is open, with no box", not bad, bad)
        planted_control("results without veil.spans read as closed even with nothing to act on",
                        "if (!k || k.status !== 'known' || (!k.spans && k.acts.length)) {",
                        "if (!k || k.status !== 'known' || !k.spans) {", lambda: nothing_closed_problems(b, base))
        nospans = json.loads(keep[rj])
        nospans.pop("veil", None)
        write(rj, json.dumps(nospans))
        bad = family_unknown(b, base, opens=False)
        check("browser: with the son's results holding findings to act on but no veil.spans every place of his stays "
              "closed", not bad, bad)
        planted_control("results without veil.spans read as nothing closed",
                        "if (!k || k.status !== 'known' || (!k.spans && k.acts.length)) {",
                        "if (!k || k.status !== 'known') {", lambda: family_unknown(b, base, opens=False))
        write(rj, keep[rj])

        # A second child whose run holds the same results under the same ids: opening is per child.
        fam2 = json.loads(planted)
        fam2["state"]["children"].append(dict(fam2["state"]["children"][0], child="DEMO-004"))
        fam2["runs"]["DEMO-004"] = fam2["runs"][SON]
        fam2["segments"]["DEMO-004"] = fam2["segments"][SON]
        fam2["de_novo"] += [dict(r, child="DEMO-004") for r in fam2["de_novo"] if r["child"] == SON]
        fam2["compound_het"] += [dict(r, child="DEMO-004") for r in fam2["compound_het"] if r["child"] == SON]
        fam2["side_by_side"]["members"].append("DEMO-004")
        for r in fam2["side_by_side"]["rows"]:
            r["values"].append(r["values"][fam2["side_by_side"]["members"].index(SON)])
            if r.get("stale"):
                r["stale"].append(False)
        write(fj, json.dumps(fam2))
        bad = two_children_problems(b, base)
        check("browser: with two children holding the same results, opening each of the son's leaves the other "
              "child's places closed", not bad, bad)
        planted_control("an opened result kept per run only, not per child",
                        "function openKey(child, id) { return (kids[child] && kids[child].url) + '|' + child + '|' + id; }",
                        "function openKey(child, id) { return (kids[child] && kids[child].url) + '|' + id; }",
                        lambda: two_children_problems(b, base))
        write(fj, planted)

        # The son's run moves to another folder (same results) while the page is open: what was opened closes again.
        again = os.path.join(out, "runs", SON + "-again")
        shutil.copytree(os.path.join(out, "runs", SON), again)
        moved = json.loads(planted)
        moved["runs"][SON]["url"] = moved["runs"][SON]["url"].replace(f"/{SON}/", f"/{SON}-again/")
        try:
            bad = new_run_problems(b, base, fj, planted, json.dumps(moved))
            check("browser: results the reader opened close again when the son's run changes while the page is open",
                  not bad, bad)
            planted_control("an opened result kept per child only, not per run",
                            "function openKey(child, id) { return (kids[child] && kids[child].url) + '|' + child + '|' + id; }",
                            "function openKey(child, id) { return child + '|' + id; }",
                            lambda: new_run_problems(b, base, fj, planted, json.dumps(moved)))
        finally:
            shutil.rmtree(again, ignore_errors=True)
            write(fj, planted)

        # A partner of the couple screen who is also a child, with a closed result in a gene of the screen.
        famc = json.loads(keep[fj])
        famc["state"]["members"].insert(0, {"sample": GRANDPARENT, "sex": "male", "father": None, "mother": None,
                                            "partner": None, "sequenced": False})
        for m in famc["state"]["members"]:
            if m["sample"] == S:
                m["father"] = GRANDPARENT
        famc["state"]["children"].insert(0, {"child": S, "father": GRANDPARENT, "mother": None, "segments": None})
        rep = f"couple_{S}_{OTHERS[0]}.md"
        famc["state"]["couples"][0]["report"] = rep
        famc["state"]["reports"] = famc["state"].get("reports", []) + [{"name": rep, "url": rep}]
        write(fj, json.dumps(famc))
        r1 = json.loads(keep[rj1])
        r1.setdefault("veil", {})["spans"] = [["two-hits-prf1", "PRF1", None, None, None]]
        write(rj1, json.dumps(r1))
        bad = couple_problems(b, base)
        check("browser: a couple screen whose partner is a child with a closed result stays closed whole (counts, "
              "table, boxes, report link, the lead's words) until that result is opened from its box", not bad, bad)
        planted_control("a couple screen that ignores the veil", "var cs = coupleShut(f, c);", "var cs = [];",
                        lambda: couple_problems(b, base))
        couple_cases = (("PRF1, a gene of the screen", [["two-hits-prf1", "PRF1", None, None, None]]),
                        ("RYR1, a gene the screen does not name", [["sf-ryr1", "RYR1", "chr19", 38440000, 38440000]]))
        bad = invariance_problems(b, base, rj1, json.loads(keep[rj1]), fj, {"the couple": json.dumps(famc)},
                                  couple_cases, S)
        check("browser: the family page is the same, byte for byte, whether the partner's closed result is in a gene "
              "of the screen or not", not bad, bad)
        planted_control("a couple screen closed only when the closed result is in one of its genes",
                        "var cs = coupleShut(f, c);",
                        "var cs = coupleShut(f, c).filter(function (x) { return (kids[x[0]].spans || []).some(function (s) { "
                        "return ((f.couples || {})[c.a + ' ' + c.b] || []).some(function (r) { return r.gene === s[1]; }); }); });",
                        lambda: invariance_problems(b, base, rj1, json.loads(keep[rj1]), fj, {"the couple": json.dumps(famc)},
                                                    couple_cases, S))
        os.replace(rj1, rj1 + ".away")
        try:
            bad = unknown_couple_problems(b, base)
            check("browser: with the partner's results unavailable, 'Show the screen' opens the screen and leaves the "
                  "partner's cells closed, and 'Show the cells' opens the cells and leaves the screen closed", not bad, bad)
            planted_control("a couple screen that 'Show the cells' opens while the results cannot be read",
                            "return opened[openKey(child, 'box:' + box)] ? [] : [UNKNOWN];",
                            "return opened[openKey(child, box.indexOf('couple:') === 0 ? 'box:side' : 'box:' + box)] ? [] : [UNKNOWN];",
                            lambda: unknown_couple_problems(b, base))
        finally:
            os.replace(rj1 + ".away", rj1)
        write(rj1, keep[rj1])

        # The drawing as an image (no stretches in family.json).
        fam = json.loads(keep[fj])
        fam["segments"] = {}
        fam["state"]["children"][0]["segments"] = fam["state"]["pedigree"]
        write(fj, json.dumps(fam))
        bad = map_problems(b, base)
        check("browser: the son's map as an image stays closed while a result is, and its box opens one result per "
              "click, the map once both are opened", not bad, bad)
        planted_control("a map image that ignores the veil", "var shut = shutOf(c.child, 'map');", "var shut = [];",
                        lambda: map_problems(b, base))
        mapfam = {"the map as an image": json.dumps(fam)}
        bad = invariance_problems(b, base, rj, json.loads(keep[rj]), fj, mapfam)
        check("browser: with the son's map as an image the family page is the same, byte for byte, whichever of his "
              "results is closed", not bad, bad)
        planted_control("a map image closed only when the closed result has a position",
                        "var shut = shutOf(c.child, 'map');",
                        "var shut = shutOf(c.child, 'map').filter(function (id) { return (kids[c.child].spans || []).some("
                        "function (s) { return s[0] === id && s[2]; }); });",
                        lambda: invariance_problems(b, base, rj, json.loads(keep[rj]), fj, mapfam))
        write(fj, json.dumps(fam))
        os.replace(rj, rj + ".away")
        try:
            bad = unknown_map_problems(b, base)
            check("browser: with the son's results unavailable, 'Show the rows' leaves his map as an image closed, and "
                  "'Show the map' opens it", not bad, bad)
            planted_control("a map image that 'Show the rows' opens",
                            "return opened[openKey(child, 'box:' + box)] ? [] : [UNKNOWN];",
                            "return opened[openKey(child, 'box:denovo')] ? [] : [UNKNOWN];",
                            lambda: unknown_map_problems(b, base))
            planted_control("de novo rows that 'Show the map' opens",
                            "return opened[openKey(child, 'box:' + box)] ? [] : [UNKNOWN];",
                            "return opened[openKey(child, box === 'denovo' ? 'box:map' : 'box:' + box)] ? [] : [UNKNOWN];",
                            lambda: unknown_map_problems(b, base))
        finally:
            os.replace(rj + ".away", rj)
    finally:
        for p, text in keep.items():
            write(p, text)


# What a reader of the son's page gets before any click: the whole document (every text and every attribute, the
# map's timing attributes aside), the address, where the page is scrolled to and what has the focus, and a hash of
# the map's canvas.
PAGE_JS = r"""JSON.stringify((function () {
  var c = document.documentElement.cloneNode(true);
  Array.prototype.forEach.call(c.querySelectorAll('[data-pgp-map-parse-ms],[data-pgp-map-paint-ms]'), function (e) {
    e.removeAttribute('data-pgp-map-parse-ms'); e.removeAttribute('data-pgp-map-paint-ms'); });
  // The root's data attributes come from scripts that answer in any order: sort them by name.
  var at = Array.prototype.map.call(c.attributes, function (x) { return [x.name, x.value]; }).sort();
  at.forEach(function (x) { c.removeAttribute(x[0]); });
  at.forEach(function (x) { c.setAttribute(x[0], x[1]); });
  var cv = document.querySelector('.gm-canvas'), a = document.activeElement;
  return { html: c.outerHTML, href: location.href + ' scroll ' + Math.round(window.scrollY) + ' focus ' +
    (a ? a.tagName + '#' + a.id : ''), canvas: cv ? cv.toDataURL() : '' };
}()))"""


def page_now(b, base, hash_="", typed=False):
    """The son's page before any click, once its results, its charts and its map have drawn: (lines of its HTML,
    its address with its scroll position and focus, its canvas hash). TYPED: the page loads without HASH_, which is
    then set as a reader types it into the address bar (a hashchange); its canvas is left out then (a hashchange
    can move the map's camera)."""
    b.go("about:blank")
    b.wait("location.href === 'about:blank'")
    b.go(f"{base}/runs/{OTHERS[1]}/{'' if typed else hash_}")
    # The son's family answers on its own request: wait for its name in the kicker too, or a pair can differ by it.
    b.wait("var d = document.documentElement.dataset; d.pgpResults && d.pgpCharts && d.pgpStatus === 'ready' && (" + MAP_READY
           + ") && document.getElementById('kicker').textContent.indexOf('Family') >= 0 && 'ready'", 15)
    if typed:
        b.eval(f"location.hash = {json.dumps(hash_)}; true")
    time.sleep(0.8)  # the margin follows the camera 140 ms after it settles
    got = json.loads(b.eval(PAGE_JS) or "{}")
    return (got.get("html") or "").replace("><", ">\n<").split("\n"), got.get("href") or "", \
        "" if typed else hashlib.sha256((got.get("canvas") or "").encode()).hexdigest()[:12]


def plant_summary(out, og, change):
    """The son's results.json built again by live_results from his summary.json with CHANGE(sections) applied, and
    his map.json and rendered charts built again from that same summary. Returns what to write back afterwards."""
    son = OTHERS[1]
    rj, mj, cj = (os.path.join(out, "runs", son, n) for n in ("results.json", "map.json", f"{son}_charts.json"))
    keep = {}
    for p in (rj, mj, cj):
        with open(p, "rb") as fh:
            keep[p] = fh.read()
    root = os.path.join(os.path.dirname(out), "root")
    with open(os.path.join(root, son, "summary.json")) as fh:
        summary = json.load(fh)
    change(summary.setdefault("sections", {}))
    res = live_results.build(son, summary)
    with open(rj, "w") as fh:
        json.dump(res, fh)
    with open(cj, "w") as fh:
        json.dump(live_results.rr.charts_json(summary), fh)
    doc = build_live_demo.og_map(og, ["build", "--root", root, "--run", son, "--assets", os.path.join(REPO, "assets"),
                                      "--cache-dir", os.path.join(root, ".owngenome-live", "map")],
                                 json.dumps({"results": res}).encode())
    with open(mj, "wb") as fh:
        fh.write(doc)
    return keep


def cpsr(gene):
    """A CPSR result to act on in GENE, the same variant whatever the gene: a closed result known by its gene."""
    def change(S):
        S["cpsr"] = {"state": "ok", "data": {"call_floor_dp": 10, "low_confidence": 0, "plp": [
            {"gene": gene, "variant": "c.68_69del", "final": "Pathogenic", "clinvar": "Pathogenic"}]}}
    return change


def deep_sites(S):
    """The ClinVar card's line for each known deep intronic site of assets/deep_intronic_sites.tsv, reference calls:
    a fixed panel the son's demo summary leaves out."""
    with open(os.path.join(REPO, "assets", "deep_intronic_sites.tsv")) as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip() and not ln.startswith("#")]
    S["clinvar"]["data"]["deep_intronic"] = [dict(zip(rows[0], r), call="ref", detail="GQ 50") for r in rows[1:]]
    return len(rows) - 1


def with_sites(gene):
    """cpsr(GENE) with the known deep intronic sites' lines (deep_sites), the same on both pages."""
    def change(S):
        cpsr(gene)(S)
        deep_sites(S)
    return change


def sv_over(locus, genes, carrier=True):
    """A class 5 deletion over GENES, and (CARRIER) a carrier row that names only its second gene."""
    def change(S):
        rows = S["annotsv"]["data"]["rows"]
        S["annotsv"]["data"]["rows"] = rows + [{"sv": locus, "type": "DEL", "gt": "0/1", "length": "-100000",
                                                 "acmg_class": "5", "genes": genes, "caller": "Manta", "support": 1,
                                                 "support_callers": ["Manta"]}]
        if not carrier:
            return
        S["carrier"]["data"]["carriers"] = [{"gene": genes[1], "variant": "c.100C>T", "genotype": "het",
                                             "significance": "Pathogenic", "stars": 2, "inheritance": "AR"}]
    return change


def recessive(change):
    """CHANGE, plus two recessive carrier rows (CFTR and HBB, inheritance AR) on the carrier card."""
    def go(S):
        change(S)
        S["carrier"]["data"]["carriers"] = S["carrier"]["data"]["carriers"] + [
            {"gene": g, "variant": v, "genotype": "het", "significance": "Pathogenic", "stars": 2, "inheritance": "AR",
             "condition": c} for g, v, c in (("CFTR", "c.1521_1523del", "Cystic fibrosis"),
                                             ("HBB", "c.20A>T", "Sickle cell disease"))]
    return go


def insertion(change):
    """CHANGE, plus a class 4 insertion (type INS) over NKD2: a review finding with a place on the map."""
    def go(S):
        change(S)
        S["annotsv"]["data"]["rows"] = S["annotsv"]["data"]["rows"] + [
            {"sv": "chr5:1000000-1000300", "type": "INS", "gt": "0/1", "length": 300, "acmg_class": "4",
             "genes": ["NKD2"], "caller": "Manta", "support": 1, "support_callers": ["Manta"]}]
    return go


def second_hit(S):
    """A second ClinVar row in the gene of the son's planted secondary finding, the card's counts kept."""
    S["clinvar"]["data"]["hits"] = S["clinvar"]["data"]["hits"] + [dict(
        [h for h in S["clinvar"]["data"]["hits"] if h.get("gene") == "KCNH2"][0], pos=150970000)]


# Pairs of son's pages that differ only in which result is closed (or in what only a closed result names): the
# page before any click must be the same, byte for byte. The two sides of each pair differ in the thing it guards:
# a gene that a fixed-panel card's own finding names (the medicine card's CYP2D6 and NAT2) against one no card names,
# a gene of the carrier list's two-hit line (SMN1) against one it does not name, a gene of the callability card's
# list of genes not on this reference (GSTT1) or of hard genes not covered (NBPF1) against one not on it, a gene
# spelled like a word or unit the report prints (HR, MB) or like a code it prints (AR as an inheritance, INS as a
# structural-variant type, TG in a score's name) against one that is not, and one guessed address on two
# pages where it is a closed result's id on one side only, loaded with the page and typed after it.
GUESS = "#f/cpsr-brca1-c-68-69del"
INFER = (("two SNV results in different genes (BRCA1, BRCA2)", cpsr("BRCA1"), cpsr("BRCA2")),
         ("a deletion over BRCA1 and NBR2 against one over CFTR and CTTNBP2",
          sv_over("chr17:43044295-43144295", ["BRCA1", "NBR2"]), sv_over("chr7:117480025-117580025", ["CFTR", "CTTNBP2"])),
         ("SMN1 against HBA1, the same copy numbers", cpsr("SMN1"), cpsr("HBA1")),
         ("a result with a place on the map (BRCA2) against one without (MLH1)", cpsr("BRCA2"), cpsr("MLH1")),
         ("a closed result with one ClinVar row against two", lambda S: None, second_hit),
         ("a closed result in a gene of a known deep intronic site (ATM) against one in none (BRCA2)", with_sites("ATM"),
          with_sites("BRCA2")),
         ("a link to each of two closed results (#f/<id>)", cpsr("BRCA1"), cpsr("BRCA2"),
          ("#f/cpsr-brca1-c-68-69del", "#f/cpsr-brca2-c-68-69del")),
         ("a gene of the medicine card's own finding (CYP2D6) against the other one (NAT2)", cpsr("CYP2D6"), cpsr("NAT2")),
         ("a gene of the medicine card's own finding (CYP2D6) against one no card names (BRCA2)", cpsr("CYP2D6"),
          cpsr("BRCA2")),
         ("a gene the carrier list's two-hit line names (SMN1) against one it does not (BRCA2)", cpsr("SMN1"),
          cpsr("BRCA2")),
         ("a gene not on this reference (GSTT1) against one on it (BRCA2)", cpsr("GSTT1"), cpsr("BRCA2")),
         ("a hard gene the reads did not cover (NBPF1) against one they did (BRCA2)", cpsr("NBPF1"), cpsr("BRCA2")),
         ("a gene spelled like a word of the report (HR) against one that is not (BRCA2)", cpsr("HR"),
          cpsr("BRCA2")),
         ("a gene spelled like a unit of the report (MB) against one that is not (BRCA2)", cpsr("MB"), cpsr("BRCA2")),
         ("an inheritance code (AR) against a gene (BRCA2), with recessive carrier rows", recessive(cpsr("AR")),
          recessive(cpsr("BRCA2"))),
         ("a deletion over AR and OPHN1 against one over BRCA1 and NBR2, with recessive carrier rows",
          recessive(sv_over("chrX:67500000-67700000", ["AR", "OPHN1"], False)),
          recessive(sv_over("chr17:43044295-43144295", ["BRCA1", "NBR2"], False))),
         ("a deletion over AR and OPHN1 against one over BRCA1 and NBR2, each with a carrier row in its second gene "
          "(OPHN1 has no place on the map, NBR2 has one)", recessive(sv_over("chrX:67500000-67700000", ["AR", "OPHN1"])),
          recessive(sv_over("chr17:43044295-43144295", ["BRCA1", "NBR2"]))),
         ("a structural-variant type (INS) against a gene (BRCA2), with an insertion over NKD2", insertion(cpsr("INS")),
          insertion(cpsr("BRCA2"))),
         ("a trait's abbreviation in a score's name (TG) against a gene (BRCA2)", cpsr("TG"), cpsr("BRCA2")),
         ("one guessed address, a closed result's id on one page only (BRCA1 against BRCA2 closed)", cpsr("BRCA1"),
          cpsr("BRCA2"), (GUESS, GUESS)),
         ("the same guessed address typed after the page drew (BRCA1 against BRCA2 closed)", cpsr("BRCA1"),
          cpsr("BRCA2"), (GUESS, GUESS), True))


def infer_problems(b, base, out, og, cases=INFER):
    """[problems]: for each pair of CASES, the son's page before any click differs (its HTML, an attribute, the
    address or the map's canvas), with the number of differing lines and the first ones."""
    bad = []
    for what, one, two, *more in cases:
        got = []
        for change, h in zip((one, two), (more or [("", "")])[0]):
            keep = plant_summary(out, og, change)
            try:
                got.append(page_now(b, base, h, typed=more[1:] == [True]))
            finally:
                unplant(keep)
        (h1, a1, c1), (h2, a2, c2) = got
        diff = [x for x in difflib.unified_diff(h1, h2, lineterm="", n=0) if x[:1] in "+-" and x[:3] not in ("+++", "---")]
        if diff or a1 != a2 or c1 != c2:
            bad.append(f"{what}: {len(diff)} lines differ, address {a1 == a2}, canvas {c1 == c2}: {diff[:6]}")
    return bad


DEEP_JS = """(function () {
  var el = document.getElementById('f-closed-1'), r = el && el.getBoundingClientRect();
  var sec = el && el.closest('section.sec'), a = sec && document.querySelector('[data-toc="' + sec.id.slice(1) + '"]');
  return JSON.stringify({ found: !!el, y: Math.round(window.scrollY), inview: !!r && r.bottom > 0 && r.top < innerHeight,
    current: a ? a.getAttribute('aria-current') : null });
}())"""


def deep_link_problems(b, base):
    """[problems]: on the son's page, #f/closed-1 and a guessed #f/<id> each scroll to the first closed result's card,
    put it in view and mark its section in the contents."""
    bad = []
    for h in ("#f/closed-1", GUESS):
        page_now(b, base, h)
        got = json.loads(b.eval(DEEP_JS) or "{}")
        if not (got.get("found") and got.get("y", 0) > 0 and got.get("inview") and got.get("current") == "true"):
            bad.append(f"{h}: {got}")
    return bad


def sv_gene_problems(b, base, out, og):
    """[problems]: with a closed deletion over BRCA1 and NBR2, the page names NBR2 nowhere (a carrier row and its
    finding name only NBR2) before any click."""
    keep = plant_summary(out, og, sv_over("chr17:43044295-43144295", ["BRCA1", "NBR2"]))
    try:
        html = "\n".join(page_now(b, base)[0])
    finally:
        unplant(keep)
    return [f"NBR2 on the page: ...{html[max(0, html.find('NBR2') - 120):html.find('NBR2') + 40]}..."] if "NBR2" in html else []


def panel_words_problems(b, base, out, og):
    """[problems]: with a closed result in CYP2D6, a gene the medicine card's own finding names, the medicine card shows
    the fixed-panel box's words, never 'This section holds a result'."""
    keep = plant_summary(out, og, cpsr("CYP2D6"))
    try:
        page_now(b, base)
        text = b.eval("(document.querySelector('#card-cpic .card-text') || {}).textContent || ''")
    finally:
        unplant(keep)
    ok = "This section lists the same genes for everyone" in text and "This section holds" not in text
    return [] if ok else [f"the medicine card says {text[:160]!r}"]


ADDRESS_JS = ("window.__seen = window.__seen || []; if (!window.__heard) { window.__heard = 1; ['hashchange', 'popstate']"
              ".forEach(function (e) { window.addEventListener(e, function () { window.__seen.push(e + ' ' + location.hash); }); "
              "}); } true")


def address_problems(b, base, out):
    """[problems]: on the son's page, a result opened from its card and the reader's moves after it never put the
    id of a result his page can close into the address or a history entry: #f/closed-1, Show the result, a link to
    section 3, Back."""
    with open(os.path.join(out, "runs", OTHERS[1], "results.json")) as fh:
        res = json.load(fh)
    tags = (res.get("veil") or {}).get("findings") or {}
    ids = [f["id"] for f in res["findings"] if f["level"] == "act" or f["id"] in tags]
    seen = []
    b.go("about:blank")
    b.wait("location.href === 'about:blank'")
    b.go(f"{base}/runs/{OTHERS[1]}/#f/closed-1")
    b.wait("document.documentElement.dataset.pgpResults && document.documentElement.dataset.pgpStatus === 'ready' && "
           "document.querySelector('.finding[data-level=act] .veil button') && 'ready'", 15)
    b.eval(ADDRESS_JS)
    seen.append(b.eval("location.hash"))
    b.eval("document.querySelector('.finding[data-level=act] .veil button').click(); true")
    time.sleep(0.3)
    seen.append(b.eval("location.hash"))
    b.eval("location.hash = '#s3'; true")
    time.sleep(0.3)
    b.eval("history.back(); true")
    time.sleep(0.5)
    seen += [b.eval("location.hash")] + json.loads(b.eval("JSON.stringify(window.__seen)") or "[]")
    bad = [h for h in seen if any(i in (h or "") for i in ids)]
    if not str(seen[1]).startswith("#f/opened-"):
        bad.append(f"the opened result's address: {seen[1]!r}, want #f/opened-<n>")
    return [f"addresses {seen}: {bad}"] if bad else []


def home_now(b, base, out, change):
    """The home page with the son's view (runs.json, live_results.overview) built from his summary with CHANGE applied:
    the whole document."""
    p = os.path.join(out, "runs.json")
    with open(p, "rb") as fh:
        keep = fh.read()
    with open(os.path.join(os.path.dirname(out), "root", OTHERS[1], "summary.json")) as fh:
        summary = json.load(fh)
    change(summary.setdefault("sections", {}))
    idx = json.loads(keep)
    idx["runs"] = [dict(r, **live_results.overview(summary)) if r["name"] == OTHERS[1] else r for r in idx["runs"]]
    try:
        with open(p, "w") as fh:
            json.dump(idx, fh)
        b.go("about:blank")
        b.wait("location.href === 'about:blank'")
        b.go(f"{base}/")
        b.wait("document.documentElement.dataset.pgpStatus === 'ready' && document.documentElement.dataset.pgpPeople "
               "&& document.querySelector('#attention li') && 'ready'", 15)
        time.sleep(0.3)
        return b.eval("document.documentElement.outerHTML").replace("><", ">\n<").split("\n")
    finally:
        with open(p, "wb") as fh:
            fh.write(keep)


def renamed(gene):
    """The son's summary with his planted gene (KCNH2) named GENE everywhere: his first result to act on, the one the
    home page lists, is then in GENE."""
    def change(S):
        S.update(json.loads(json.dumps(S).replace("KCNH2", gene)))
    return change


HOME = (("his first result in BRCA1 against BRCA2", renamed("BRCA1"), renamed("BRCA2")),
        ("a gene of the medicine card's own finding (CYP2D6) against BRCA2", cpsr("CYP2D6"), cpsr("BRCA2")))


def home_problems(b, base, out, cases=HOME):
    """[problems]: the home page is the same whichever of the son's results is closed (his tile and the list of what
    needs attention: no link to a closed result's own page, whose address names its gene)."""
    bad = []
    for what, one, two in cases:
        h1, h2 = home_now(b, base, out, one), home_now(b, base, out, two)
        links = json.loads(b.eval("JSON.stringify(Array.from(document.querySelectorAll('#attention a')).filter(function (a) { var who = a.querySelector('.who .pc'); return who && who.textContent === " + json.dumps(SON) + "; }).map(function (a) { return a.getAttribute('href'); }))") or "[]")
        if any("finding/" in link for link in links):
            bad.append(f"home, {what}: a closed summary links to a finding page")
        diff = [x for x in difflib.unified_diff(h1, h2, lineterm="", n=0) if x[:1] in "+-" and x[:3] not in ("+++", "---")]
        if diff:
            bad.append(f"home, {what}: {len(diff)} lines differ: {diff[:4]}")
    return bad


def two_hits(change):
    """CHANGE, plus a two-hit gene in CFTR on the carrier list (a result to act on, two-hits-cftr), whose line prints
    the inheritance code AR."""
    def go(S):
        change(S)
        th = S["carrier"]["data"].setdefault("two_hits", {"deletions_read": [], "deletions_not_read": []})
        th["genes"] = [{"gene": "CFTR", "verdict": "possible", "condition": "Cystic fibrosis",
                        "deletion": "DEL chr7:117480025-117580025", "variant": "c.1521_1523del", "genotype": "het",
                        "why": "a deletion and a variant on the other copy"}]
    return go


def cpsrs(*genes):
    """A closed CPSR result in each of GENES, the same variant in each."""
    def change(S):
        S["cpsr"] = {"state": "ok", "data": {"call_floor_dp": 10, "low_confidence": 0, "plp": [
            {"gene": g, "variant": "c.68_69del", "final": "Pathogenic", "clinvar": "Pathogenic"} for g in genes]}}
    return change


def dup_over(change, locus, genes):
    """CHANGE, plus a class 4 duplication over GENES (a finding worth a look, with a place on the map)."""
    def go(S):
        change(S)
        S["annotsv"]["data"]["rows"] = S["annotsv"]["data"]["rows"] + [
            {"sv": locus, "type": "DUP", "gt": "0/1", "length": 300000, "acmg_class": "4", "genes": genes,
             "caller": "Manta", "support": 1, "support_callers": ["Manta"]}]
    return go


def cid(gene):
    """The id of cpsr(GENE)'s result."""
    return f"cpsr-{gene.lower()}-c-68-69del"


# Pairs of the son's pages that differ only in which result is closed, compared after the same click: the result
# each side opens through its handle (its place among the closed cards, as a reader who never saw its id opens it).
# A place tied to a closed result must open with the results it names and with nothing else, so the page after a
# click is again a function of the open results and the count of closed ones. The two-hit line prints the
# inheritance code AR: a closed result in AR (a CPSR result, or a deletion over AR) must not hold it closed.
CLICKS = (("a two-hit result in CFTR opened, the other closed result in AR against BRCA2", two_hits(cpsr("AR")),
           two_hits(cpsr("BRCA2")), ("two-hits-cftr", "two-hits-cftr")),
          ("a two-hit result in CFTR opened, a closed deletion over AR and OPHN1 against one over BRCA1 and NBR2",
           two_hits(sv_over("chrX:67500000-67700000", ["AR", "OPHN1"], False)),
           two_hits(sv_over("chr17:43044295-43144295", ["BRCA1", "NBR2"], False)), ("two-hits-cftr", "two-hits-cftr")))


def opened_page(b, base, fid):
    """The son's page after FID is opened through its handle: (lines of its HTML, its address, its canvas hash)."""
    page_now(b, base)
    if not open_closed(b, fid):
        return None
    time.sleep(0.8)
    got = json.loads(b.eval(PAGE_JS) or "{}")
    return (got.get("html") or "").replace("><", ">\n<").split("\n"), got.get("href") or "", \
        hashlib.sha256((got.get("canvas") or "").encode()).hexdigest()[:12]


def click_problems(b, base, out, og, cases=CLICKS):
    """[problems]: for each pair of CASES, the son's page after each side opens its result differs (its HTML, an
    attribute, the address or the map's canvas), with the number of differing lines and the first ones."""
    bad = []
    for what, one, two, fids in cases:
        got = []
        for change, fid in zip((one, two), fids):
            keep = plant_summary(out, og, change)
            try:
                got.append(opened_page(b, base, fid))
            finally:
                unplant(keep)
        if None in got:
            bad.append(f"{what}: a side has no closed card for {fids}")
            continue
        (h1, a1, c1), (h2, a2, c2) = got
        diff = [x for x in difflib.unified_diff(h1, h2, lineterm="", n=0) if x[:1] in "+-" and x[:3] not in ("+++", "---")]
        if diff or a1 != a2 or c1 != c2:
            bad.append(f"{what}: {len(diff)} lines differ, address {a1 == a2}, canvas {c1 == c2}: {diff[:6]}")
    return bad


# Two closed results, one opened (each in turn): the other's gene must be on the page no more often than before the
# click (the search box's example names BRCA2 on every page). A place that names both (the lead's segment of the
# findings worth a look, which prints the ClinVar findings in DEMOC and DEMOA; a duplication over BRCA1 and BRCA2)
# stays closed while either is.
TWO_CLOSED = (("two CPSR results named in one segment of the lead, DEMOC opened", cpsrs("DEMOC", "DEMOA"), cid("DEMOC"),
               "DEMOA"),
              ("the same, DEMOA opened", cpsrs("DEMOC", "DEMOA"), cid("DEMOA"), "DEMOC"),
              ("a class 4 duplication over BRCA1 and BRCA2, both closed, BRCA1 opened",
               dup_over(cpsrs("BRCA1", "BRCA2"), "chr5:1000000-1300000", ["BRCA1", "BRCA2"]), cid("BRCA1"), "BRCA2"),
              ("the same, BRCA2 opened", dup_over(cpsrs("BRCA1", "BRCA2"), "chr5:1000000-1300000", ["BRCA1", "BRCA2"]),
               cid("BRCA2"), "BRCA1"))


def two_closed_problems(b, base, out, og, cases=TWO_CLOSED):
    """[problems]: for each of CASES, the still closed result's gene is on the page more often after the click."""
    bad = []
    for what, change, fid, still in cases:
        keep = plant_summary(out, og, change)
        try:
            before = "\n".join(page_now(b, base)[0]).count(still)
            after = opened_page(b, base, fid)
        finally:
            unplant(keep)
        if after is None:
            bad.append(f"{what}: no closed card for {fid}")
            continue
        html = "\n".join(after[0])
        if html.count(still) > before:
            at = html.find(still, 0)
            bad.append(f"{what}: {still} {before} times before, {html.count(still)} after: "
                       f"...{html[max(0, at - 160):at + 40]}...")
    return bad


def apoe_pair_problems(b, base, out, og):
    """Opening the map's own gate cannot say which other result is still closed."""
    got, bad = [], []
    for gene in ("APOE", "BRCA2"):
        keep = plant_summary(out, og, cpsr(gene))
        try:
            stages = [page_now(b, base)]
            if not b.eval("(function () { var c = document.querySelector('.gm-chip[data-closed=gate]'); "
                          "if (!c) { return false; } c.querySelector('button').click(); "
                          "document.querySelector('.gm-chip[data-closed=gate] .gm-yes').click(); return true; }())"):
                return [f"{gene}: no APOE gate"]
            time.sleep(0.2)
            stages.append(page_snapshot(b))
            if not b.eval("!!document.querySelector('.gm-marks [data-id=apoe]')"):
                bad.append(f"{gene}: opening APOE did not draw its mark")
            if not b.eval("(function () { var c = Array.from(document.querySelectorAll('.gm-chip.opened')).find("
                          "function (e) { return e.textContent.indexOf('Opened: APOE') === 0; }); "
                          "if (!c) { return false; } c.querySelector('button').click(); return true; }())"):
                return bad + [f"{gene}: no close button for APOE"]
            time.sleep(0.2)
            stages.append(page_snapshot(b))
            if b.eval("!!document.querySelector('.gm-marks [data-id=apoe]')"):
                bad.append(f"{gene}: closing APOE kept its mark")
            got.append(stages)
        finally:
            unplant(keep)
    for stage, a, z in zip(("before", "opened", "closed again"), *got):
        if a != z:
            bad.append(f"APOE against BRCA2 {stage}: HTML {a[0] == z[0]}, address {a[1] == z[1]}, canvas {a[2] == z[2]}")
    return bad


def page_snapshot(b):
    """The same complete document, address and canvas comparison used by page_now, after a click."""
    got = json.loads(b.eval(PAGE_JS) or "{}")
    return (got.get("html") or "").replace("><", ">\n<").split("\n"), got.get("href") or "", \
        hashlib.sha256((got.get("canvas") or "").encode()).hexdigest()[:12]


def related_set_problems(b, base, out, og):
    """On each finding's own page, opening one of two results cannot open a related finding naming both."""
    made, bad = [], []
    folder = os.path.join(out, "runs", SON, "finding")
    # The demo server is static: give the planted findings the same page template as its original findings.
    original = "sf-kcnh2-chr7-150960000-c-t"
    with open(os.path.join(folder, original, "index.html")) as fh:
        template = fh.read()
    try:
        for genes in (["BRCA1", "BRCA2"], ["BRCA2", "BRCA1"]):
            # Related links use the finding's first gene. Reverse that as well as the CPSR order, so both
            # directions can actually list the duplication and expose a first-id reader.
            gene, fid = genes[0], cid(genes[0])
            keep = plant_summary(out, og, dup_over(cpsrs(*genes), "chr5:1000000-1300000", genes))
            try:
                with open(os.path.join(out, "runs", SON, "results.json")) as fh:
                    finding = next(f for f in json.load(fh)["findings"] if f["id"] == fid)
                if "sv-chr5-1000000-1300000" not in finding.get("related", []):
                    bad.append(f"{gene}: fixture never offers the duplication as a related result")
                dest = os.path.join(folder, fid)
                os.makedirs(dest)
                made.append(dest)
                with open(os.path.join(dest, "index.html"), "w") as fh:
                    fh.write(template.replace(original, fid))
                b.go("about:blank")
                b.wait("location.href === 'about:blank'")
                b.go(f"{base}/runs/{SON}/finding/{fid}/")
                if not b.wait("document.documentElement.dataset.pgpStatus === 'ready' && document.querySelector('.veil button') "
                              "&& 'ready'", 10):
                    bad.append(f"{gene}: finding never closed")
                    continue
                b.eval("document.querySelector('.veil button').click(); true")
                if not b.wait("!document.querySelector('.veil') && document.querySelector('#related') && 'opened'", 5):
                    bad.append(f"{gene}: finding never opened")
                related = b.eval("document.querySelector('#related').outerHTML") or ""
                if "sv-chr5-1000000-1300000" in related or "Large DUP" in related:
                    bad.append(f"{gene} opened alone: the duplication names the other closed result in Related results")
                if not b.eval("document.querySelector('#what').textContent.indexOf(" + json.dumps(gene) + ") >= 0"):
                    bad.append(f"{gene}: opening its finding did not show it")
            finally:
                unplant(keep)
    finally:
        for dest in made:
            shutil.rmtree(dest)
    return bad


def chart_position_problems(b, base, out, og):
    """Gene-only chart tags and coordinate-only mosaic tags, including reopening with two overlapping acts."""
    def change(S):
        dup_over(cpsr("KCNQ1"), "chr5:1000000-1300000", ["NKD2", "KCNQ1"])(S)
        second_hit(S)  # KCNH2, chr7:150970000: outside both closed results' positions, tied by gene alone.
        S["mosaic"] = {"state": "ok", "data": {"status": "assessed", "floor": 0.05, "depth": 30,
            "floor_read": 0.05, "floor_min": 0.05, "het_sites": 1000, "chromosomes": 22, "calls_reported": 1,
            "calls": [{"chrom": "chr7", "start": 150900000, "end": 151000000, "length_mb": 0.1,
                       "type": "Loss", "arms": "7q", "cell_fraction": 0.2}]}}
    keep = plant_summary(out, og, change)
    bad, states = [], []
    try:
        with open(os.path.join(out, "runs", SON, f"{SON}_charts.json")) as fh:
            marks = next(sp["marks"] for sp in json.load(fh)["charts"] if sp["id"] == "genome")
        hit = next(m["id"] for m in marks if m.get("kind") == "hit" and m.get("pos") == 150970000)
        sv = next(m["id"] for m in marks if m.get("kind") == "sv" and m.get("pos") == 1000000)
        # Inspect each surface separately. A link/title leak counts even when its visible row is absent.
        expression = """JSON.stringify((function () {
          var c = document.querySelector('#card-mosaic-mca');
          return [!!document.querySelector('figure[data-chart=genome] [data-mark=%s]'),
            !!document.querySelector('figure[data-chart=genome] [data-mark=%s]'),
            !!document.querySelector('[data-finding=mca-chr7-150900000]'),
            !!(c && Array.from(c.querySelectorAll('tbody tr')).some(function (r) { return r.textContent.includes('150900000'); })),
            !!document.querySelector('a[href*="150900000"]'),
            !!(c && Array.from(c.querySelectorAll('li')).some(function (e) { return e.textContent.includes('150900000'); }))];
        }()))""" % (hit, sv)
        page_now(b, base)
        states.append(json.loads(b.eval(expression)))
        for fid in ("sf-kcnh2-chr7-150960000-c-t", "sv-chr7-150955000-150965000", cid("KCNQ1")):
            if not open_closed(b, fid):
                return [f"no closed card for {fid}"]
            states.append(json.loads(b.eval(expression)))
        # Positive first: every planted mark and surface must reach the browser when opened. Otherwise a missing
        # rendered charts file could make the hidden-state assertions vacuous.
        if states[-1] != [True] * 6:
            bad.append(f"all opened: planted hit, SV, mosaic finding, row, viewer link and line must exist: {states[-1]}")
        for stage, got, want in zip(("all closed", "SNV alone opened", "SNV and deletion opened"), states,
                                   ([False] * 6, [False] * 6, [True, False, True, True, True, True])):
            if got != want:
                bad.append(f"{stage}: hit, SV, mosaic finding/row/link/line {got}, want {want}")
        # Reverse the overlapping results' opening order too. The SNV has no chip of its own on this demo map.
        page_now(b, base)
        if not open_closed(b, "sv-chr7-150955000-150965000"):
            bad.append("no closed deletion after reloading")
        elif json.loads(b.eval(expression)) != [False] * 6:
            bad.append(f"deletion alone opened: {b.eval(expression)}")
    finally:
        unplant(keep)
    return bad


def home_panel_problems(b, base, out, og):
    """A panel finding's title and id stay off a closed child's home view, but remain for adults and no-act runs."""
    p = os.path.join(out, "runs.json")
    with open(p, "rb") as fh:
        original = fh.read()
    idx = json.loads(original)
    idx["runs"] = [dict(r, state="running") for r in idx["runs"] if r["name"] == SON]
    bad = []

    def change(S):
        S.clear()
        cpsr("TP53")(S)
        S["chip"] = {"state": "ok", "data": {"calls": [{"gene": "TP53", "protein_change": "p.Demo",
            "vaf": 0.1, "alt_reads": 10, "depth": 100, "status": "reported"}]}}

    try:
        for case in ("child", "adult", "no act"):
            view = copy.deepcopy(idx)
            if case == "adult":
                view["families"] = []
            with open(p, "w") as fh:
                json.dump(view, fh)

            def plant(S):
                change(S)
                if case == "no act":
                    S["cpsr"]["data"]["plp"] = []

            html = "\n".join(home_now(b, base, out, plant))
            if case == "child":
                if "TP53" in html or "chip-tp53-p-demo" in html:
                    bad.append("closed child's home view exposes a panel finding's title or id")
            elif "A blood clone in TP53" not in html or "finding/chip-tp53-p-demo/" not in html:
                bad.append(f"{case}: the open panel title or finding link is missing")
    finally:
        with open(p, "wb") as fh:
            fh.write(original)
    return bad


def gate_regression_checks(b, base, out, og):
    """Each guard has a browser control that removes exactly the reviewed rule."""
    t0 = time.time()
    for desc, probe in (("APOE's gate before opening, opened and closed again gives the same HTML, address and canvas "
                         "with APOE or BRCA2 closed", apoe_pair_problems),
                        ("a finding page opened for either BRCA1 or BRCA2 keeps their duplication out of Related results",
                         related_set_problems),
                        ("a child's home view hides panel finding titles and links, while adults and no-act runs keep them",
                         home_panel_problems),
                        ("planted charts reach the browser and close by gene alone, and a gene-free mosaic finding, row, "
                         "viewer link and line close until both overlapping results open", chart_position_problems)):
        bad = probe(b, base, out, og)
        check("browser: " + desc, not bad, bad)
    for file, a, z, desc, probe in (
            ("map.js", "var tag = !!r.gate ||", "var tag =", "map-owned APOE row following overlapping act tags", apoe_pair_problems),
            ("finding.js", "return t.length ? t : null;", "return t.length ? t.slice(0, 1) : null;",
             "finding.js tagOf keeping only the first id", related_set_problems),
            ("home.js", "if (!closedSummary(r)) { return r.top || []; }", "if (true) { return r.top || []; }",
             "home.js ignoring count-based child summaries", lambda *args: fail_closed_problems(*args, cases=["home"]))):
        p = os.path.join(out, "assets", file)
        with open(p) as fh:
            source = fh.read()
        try:
            if source.count(a) != 1:
                check("control: " + desc, False, "missing mutation anchor")
                continue
            with open(p, "w") as fh:
                fh.write(source.replace(a, z))
            control(desc, not probe(b, base, out, og))
        finally:
            with open(p, "w") as fh:
                fh.write(source)
    names, spans = live_results.mark_names, live_results.closable
    try:
        for kind, desc in (("hit", "ClinVar chart mark without its gene"), ("sv", "SV chart mark with only its first gene")):
            live_results.mark_names = lambda S, sp, m: (([], "") if kind == "hit" else (names(S, sp, m)[0][:1], "")) \
                if m.get("kind") == kind else names(S, sp, m)
            control(desc, not chart_position_problems(b, base, out, og))
        live_results.mark_names = names
        live_results.closable = lambda found: [(fid, g, None, None, None) for fid, g, *_ in spans(found)]
        control("closable without coordinates", not chart_position_problems(b, base, out, og))
    finally:
        live_results.mark_names, live_results.closable = names, spans
    print(f"   gate regression checks: {time.time() - t0:.2f} s")


def fail_closed_problems(b, base, out, og, cases=None, producers=None):
    """Producer plantings: home rankings, ontology facts, chart geometry and structured mitochondrial sites."""
    bad = []
    cases = cases or ("home", "hpo", "contig", "mito", "producers", "metadata")
    if "home" in cases:
        for n in (1, 2, 3):
            def ranked(g):
                def change(S):
                    S.clear()
                    cpsrs(g, *["MLH1", "MSH2"][:n - 1])(S)
                    S["exomiser"] = {"state": "ok", "data": {"terms": [], "genes": [
                        {"gene": "BRCA1", "rank": 1, "moi": "AD", "exomiser_score": 0.9, "variants": []}]}}
                return change
            a, z = (home_now(b, base, out, ranked(g)) for g in ("BRCA1", "BRCA2"))
            if a != z:
                bad.append(f"home: {n} closed results change the whole document with the same BRCA1 ranking")
    if "hpo" in cases:
        import collect_summary
        with tempfile.TemporaryDirectory() as td:
            folder = os.path.join(td, "exomiser")
            os.mkdir(folder)
            with open(os.path.join(folder, f"{SON}_exomiser_top.tsv"), "w") as fh:
                fh.write("rank\tgene\tmoi\texomiser_score\tphenotype_score\tvariant\tgenotype\tconsequence\thgvs\tacmg\tdisease\n"
                         "1\tCFTR\tAR\t0.9\t0.8\tchr7:117559590\thet\tmissense\tc.1521_1523del\tPathogenic\tDemo condition\n")
            with open(os.path.join(folder, f"{SON}_exomiser_terms.tsv"), "w") as fh:
                fh.write("hpo_id\tlabel\nHP:0020055\tDiminished MLH1 expression\n")
            _, data = collect_summary.sec_exomiser(td, SON)
        def hpo(S):
            cpsrs("MLH1", "MSH2")(S)
            S["exomiser"] = {"state": "ok", "data": data}
        keep = plant_summary(out, og, hpo)
        try:
            page_now(b, base)
            for stage in ("closed", "partly opened", "all opened"):
                if stage == "partly opened":
                    assert open_closed(b, cid("MLH1"))
                if stage == "all opened":
                    with open(os.path.join(out, "runs", SON, "results.json")) as fh:
                        acts = [f["id"] for f in json.load(fh)["findings"] if f["level"] == "act"]
                    for fid in acts:
                        open_closed(b, fid)
                visible = b.eval("document.querySelector('#card-exomiser').textContent.includes('Diminished MLH1 expression')")
                if visible != (stage == "all opened"):
                    bad.append(f"hpo: label visibility {visible} while {stage}")
        finally:
            unplant(keep)
    if "contig" in cases:
        views = []
        for locus in ("chrUn_KI270442v1:100-200", "chr13:32300000-32300100"):
            keep = plant_summary(out, og, sv_over(locus, ["BRCA1"], False))
            try:
                views.append(page_now(b, base))
                with open(os.path.join(out, "runs", SON, "results.json")) as fh:
                    acts = [f["id"] for f in json.load(fh)["findings"] if f["level"] == "act"]
                for fid in acts:
                    assert open_closed(b, fid), fid
                note = b.eval("document.querySelector('figure[data-chart=genome]').textContent.includes('SVs on other contigs')")
                if note != locus.startswith("chrUn"):
                    bad.append("contig: the full chart note did not return after every act opened")
                assert b.eval("(function () { var b = document.querySelector('.gm-chip.opened button'); if (!b) return false; b.click(); return true; }())")
                if b.eval("document.querySelector('figure[data-chart=genome]').textContent.includes('SVs on other contigs')"):
                    bad.append("contig: re-closing an act leaves the chart note visible")
            finally:
                unplant(keep)
        if views[0] != views[1]:
            bad.append("contig: closed SV changes text, attributes, address, canvas or SVG geometry")
    if "mito" in cases:
        def mito(S):
            cpsrs("MT-TL1", "MT-TL2")(S)
            for row in S["cpsr"]["data"]["plp"]:
                row["variant"] = "chrM:3243 A>G"
            S["mito"] = {"state": "ok", "data": {"pass": 1, "heteroplasmic": 1,
                "heteroplasmic_sites": [{"pos": 3243, "ref": "A", "alt": "G", "af": 0.25}]}}
        keep = plant_summary(out, og, mito)
        try:
            page_now(b, base)
            expr = "JSON.stringify([!!document.querySelector('#card-mito tbody tr'), !!document.querySelector('figure[data-chart=mito] [data-mark=m3243-G]')])"
            for stage, fid in enumerate((None, "cpsr-mt-tl1-chrm-3243-a-g", "cpsr-mt-tl2-chrm-3243-a-g")):
                if fid and not open_closed(b, fid):
                    bad.append(f"mito: no closed result {fid}")
                    break
                got = json.loads(b.eval(expr))
                if got != [stage == 2] * 2:
                    bad.append(f"mito: row/mark after {stage} matching openings {got}")
        finally:
            unplant(keep)
    if "report-metadata" in cases:
        phone_at(b, 390)
        try:
            views = []
            for gene in ("MLH1", "BRCA2"):
                build = live_results.build
                def metadata(sample, summary, files=()):
                    token = "MLH1_METADATA_NOTE"
                    records = [token] * (2 if gene == "MLH1" else 1)
                    planted = dict(summary, not_assessed=records, not_run=records,
                                   databases=dict(summary.get("databases") or {}, test_metadata=token))
                    return build(sample, planted, files)
                live_results.build = metadata
                try:
                    keep = plant_summary(out, og, cpsr(gene))
                finally:
                    live_results.build = build
                try:
                    views.append(page_now(b, base))
                    if "MLH1_METADATA_NOTE" in "\n".join(views[-1][0]):
                        bad.append("report metadata: unclassified limitation, missing-step or version text reaches the DOM")
                    with open(os.path.join(out, "runs", SON, "results.json")) as fh:
                        acts = [f["id"] for f in json.load(fh)["findings"] if f["level"] == "act"]
                    for fid in acts:
                        assert open_closed(b, fid), fid
                    if "MLH1_METADATA_NOTE" not in (b.eval("document.documentElement.outerHTML") or ""):
                        bad.append("report metadata: text does not return after every act opens")
                finally:
                    unplant(keep)
            if views[0] != views[1]:
                bad.append("report metadata: extra records change closed text, attributes, address, canvas or SVG sizes")
        finally:
            b.call("Emulation.clearDeviceMetricsOverride", {}, b.session)
    if "trust" in cases:
        def qc_note(S):
            cpsr("MLH1")(S)
            S["sample_qc"] = {"state": "ok", "data": {"sex_check": "not_checked", "inferred_sex": "male",
                "sex_check_reason": "MLH1_TRUST_NOTE", "contamination": "not_checked"}}
        keep = plant_summary(out, og, qc_note)
        try:
            page_now(b, base)
            if "MLH1_TRUST_NOTE" in (b.eval("document.documentElement.outerHTML") or ""):
                bad.append("trust: a PANEL QC fact survives in the trust-line text or title")
            assert open_closed(b, cid("MLH1"))
            if "MLH1_TRUST_NOTE" in (b.eval("document.documentElement.outerHTML") or ""):
                bad.append("trust: opening one act exposes a QC fact while other acts remain closed")
            with open(os.path.join(out, "runs", SON, "results.json")) as fh:
                acts = [f["id"] for f in json.load(fh)["findings"] if f["level"] == "act"]
            for fid in acts:
                open_closed(b, fid)
            if "MLH1_TRUST_NOTE" not in (b.eval("document.querySelector('.trustline').outerHTML") or ""):
                bad.append("trust: the QC fact does not return after every act opens")
        finally:
            unplant(keep)
    if "no-act" in cases:
        def no_act(S):
            S.clear()
            S["exomiser"] = {"state": "ok", "data": {"terms": [], "genes": [
                {"gene": "BRCA1", "rank": 1, "moi": "AD", "exomiser_score": 0.9, "variants": []}]}}
        keep = plant_summary(out, og, no_act)
        try:
            page_now(b, base, "#f/exomiser-brca1")
            if b.eval("location.hash") != "#f/exomiser-brca1":
                bad.append("no-act: missing metadata changes an open finding's address")
        finally:
            unplant(keep)
    if "metadata" in cases:
        keep = plant_summary(out, og, cpsr("MLH1"))
        try:
            path = os.path.join(out, "runs", SON, "results.json")
            with open(path) as fh:
                result = json.load(fh)
            result["veil"] = {k: v for k, v in result["veil"].items() if k in ("spans", "panels")}
            card = next(c for ch in result["chapters"] for c in ch["sections"] if c["id"] == "qc")
            token = "MLH1_UNCLASSIFIED_NOTE"
            card.update(headline=token, facts=[["Unknown", token]], lines=[token], notes=[token],
                        table={"columns": ["Detail"], "rows": [[token]]})
            with open(path, "w") as fh:
                json.dump(result, fh)
            if token in "\n".join(page_now(b, base)[0]):
                bad.append("metadata: a missing tag exposes unclassified content")
            for finding in result["findings"]:
                if finding["level"] == "act":
                    assert open_closed(b, finding["id"]), finding["id"]
            if token not in (b.eval("document.documentElement.outerHTML") or ""):
                bad.append("metadata: content does not return after every act opens")
        finally:
            unplant(keep)
    if "producers" in cases:
        # Each adapter receives an unknown free-text field and an extra unstructured row, in turn. Neither
        # the token's spelling nor its record count earns an exemption from the default PANEL rule.
        for producer, adapter in list(live_results.ADAPTERS.items()):
            if producers and producer not in producers:
                continue
            def injected(S, summary):
                got = adapter(S, summary)
                got.update(facts=[["Mean depth", "MLH1 sample note"]], lines=["MLH1 sample note"],
                           table={"columns": ["Detail"], "rows": [["MLH1"], ["extra MLH1 record"]]})
                return got
            live_results.ADAPTERS[producer] = injected
            try:
                # Supply each adapter's real demo section. Missing sections are tested by the unit counterpart.
                views = []
                for gene in ("MLH1", "BRCA2"):
                    keep = plant_summary(out, og, cpsr(gene))
                    try:
                        views.append(page_now(b, base))
                        if any(token in "\n".join(views[-1][0]) for token in ("MLH1 sample note", "extra MLH1 record")):
                            bad.append(f"producer {producer}: unclassified token reached the closed DOM")
                    finally:
                        unplant(keep)
                if views[0] != views[1]:
                    bad.append(f"producer {producer}: unclassified planting changes the closed document")
            finally:
                live_results.ADAPTERS[producer] = adapter
    return bad


def panel_finding_problems(b, base, out, og, cases=("child", "missing tags", "adult", "no acts")):
    """A PANEL finding opens one act per click and reveals its content only after the last act opens."""
    folder = os.path.join(out, "runs", SON, "finding")
    original, fid = "sf-kcnh2-chr7-150960000-c-t", "chip-tp53-p-demo"
    dest = os.path.join(folder, fid)
    with open(os.path.join(folder, original, "index.html")) as fh:
        template = fh.read()
    fj = os.path.join(out, "runs", FAM, "family.json")
    with open(fj) as fh:
        family = fh.read()
    os.makedirs(dest)
    bad = []
    try:
        with open(os.path.join(dest, "index.html"), "w") as fh:
            fh.write(template.replace(original, fid))
        for case in cases:
            def change(S):
                S.clear()
                cpsrs("TP53", "BRCA2")(S)
                if case == "no acts":
                    S["cpsr"]["data"]["plp"] = []
                S["chip"] = {"state": "ok", "data": {"calls": [{"gene": "TP53", "protein_change": "p.Demo",
                    "vaf": 0.1, "alt_reads": 10, "depth": 100, "status": "reported"}]}}
            keep = plant_summary(out, og, change)
            try:
                rj = os.path.join(out, "runs", SON, "results.json")
                with open(rj) as fh:
                    res = json.load(fh)
                tags = (res.get("veil", {}).get("findings") or {}).get(fid)
                if case != "no acts" and tags != ["*panel*"]:
                    bad.append(f"{case}: producer did not classify the CHIP finding as PANEL: {tags}")
                if case == "missing tags":
                    del res["veil"]["findings"][fid]
                    with open(rj, "w") as fh:
                        json.dump(res, fh)
                fam = json.loads(family)
                if case == "adult":
                    for member in fam["state"]["members"]:
                        if member["sample"] == SON:
                            member["father"] = member["mother"] = None
                with open(fj, "w") as fh:
                    json.dump(fam, fh)
                b.go("about:blank")
                b.wait("location.href === 'about:blank'")
                b.go(f"{base}/runs/{SON}/finding/{fid}/")
                if not b.wait("document.documentElement.dataset.pgpStatus === 'ready' && 'ready'", 10):
                    bad.append(f"{case}: finding never ready")
                    continue
                for clicks in range(3 if case in ("child", "missing tags") else 1):
                    should_open = case in ("adult", "no acts") or clicks == 2
                    visible = b.eval("document.querySelector('#what').textContent.includes('A blood clone in TP53') "
                                     "&& document.title.includes('TP53') && !!document.querySelector('#evidence')")
                    shut = b.eval("!!document.querySelector('.veil button')")
                    if bool(visible) != should_open or bool(shut) == should_open:
                        bad.append(f"{case}: after {clicks} clicks PANEL finding visible={visible}, closed={shut}; "
                                   f"expected visible={should_open}")
                    if not should_open:
                        b.eval("var button = document.querySelector('.veil button'); button && button.click(); true")
                if case == "child":
                    b.go(f"{base}/runs/{SON}/finding/{fid}/")
                    if not b.wait("document.documentElement.dataset.pgpStatus === 'ready' && document.querySelector('.veil button') "
                                  "&& !document.querySelector('#evidence') && 'closed'", 10):
                        bad.append("child: a reload did not close the PANEL finding again")
            finally:
                unplant(keep)
    finally:
        with open(fj, "w") as fh:
            fh.write(family)
        shutil.rmtree(dest)
    return bad


def fail_closed_checks(b, base, out, og):
    started = time.monotonic()
    bad = panel_finding_problems(b, base, out, og)
    check("browser: a PANEL finding stays closed until one click per act, with missing tags, adult and no-act controls",
          not bad, bad)
    finding_path = os.path.join(out, "assets", "finding.js")
    with open(finding_path) as fh:
        finding_reader = fh.read()
    for old, new, desc in (
            ("return t === PANEL ? actShut() : !opened[t];", "return t === PANEL ? false : !opened[t];",
             "a finding page that ignores PANEL"),
            ("if (nextAct) { opened[nextAct.id] = true; }",
             "if (nextAct) { (r.findings || []).forEach(function (x) { if (x.level === 'act') { opened[x.id] = true; } }); }",
             "a PANEL finding button that opens every act in one click")):
        try:
            assert finding_reader.count(old) == 1
            with open(finding_path, "w") as fh:
                fh.write(finding_reader.replace(old, new))
            control(desc, not panel_finding_problems(b, base, out, og, ["child"]))
        finally:
            with open(finding_path, "w") as fh:
                fh.write(finding_reader)
    for case in ("home", "hpo", "contig", "mito", "producers", "metadata", "no-act", "trust", "report-metadata"):
        bad = fail_closed_problems(b, base, out, og, [case])
        check("browser: fail-closed producers " + case, not bad, bad)
    prose, row, veil = live_results.prose_tag, live_results.row_tag, live_results.veil
    try:
        live_results.FACTS_NO_GENE["exomiser"] = "unsafe HPO labels"
        control("an HPO fact declared gene-free", not fail_closed_problems(b, base, out, og, ["hpo"]))
        live_results.FACTS_NO_GENE.pop("exomiser", None)
        def no_mito(*args, **kwargs):
            got = veil(*args, **kwargs)
            got.get("rows", {}).pop("mito", None)
            return got
        live_results.veil = no_mito
        control("mitochondrial rows without structured positions", not fail_closed_problems(b, base, out, og, ["mito"]))
        live_results.veil = veil
        live_results.prose_tag = lambda spans, text: None
        control("unclassified producer prose left open", not fail_closed_problems(b, base, out, og, ["hpo"]))
        live_results.prose_tag = prose
        live_results.row_tag = lambda spans, columns, values, genes=(): \
            live_results.veil_tags(spans, live_results.gene_cell(values[0])) if columns == ["Detail"] \
            else row(spans, columns, values, genes)
        control("unclassified rows read as gene cells", not fail_closed_problems(b, base, out, og, ["producers"], ["qc"]))
    finally:
        live_results.FACTS_NO_GENE.pop("exomiser", None)
        live_results.prose_tag, live_results.row_tag, live_results.veil = prose, row, veil
    path = os.path.join(out, "assets", "person.js")
    with open(path) as fh:
        reader = fh.read()
    try:
        old = "parseFromString(sp.svg_closed,"
        assert reader.count(old) == 1
        with open(path, "w") as fh:
            fh.write(reader.replace(old, "parseFromString(sp.svg,"))
        control("chart notes and geometry outside PANEL", not fail_closed_problems(b, base, out, og, ["contig"]))
    finally:
        with open(path, "w") as fh:
            fh.write(reader)
    try:
        old = "tag === undefined ? tagClosed(PANEL)"
        assert reader.count(old) == 1
        with open(path, "w") as fh:
            fh.write(reader.replace(old, "tag === undefined ? false"))
        control("missing tag metadata treated as open", not fail_closed_problems(b, base, out, og, ["metadata"]))
    finally:
        with open(path, "w") as fh:
            fh.write(reader)
    try:
        old = "if (t === undefined) { t = ((results && results.findings) || []).some(function (x) { return x.level === 'act'; }) ? [PANEL] : null; }"
        assert reader.count(old) == 1
        with open(path, "w") as fh:
            fh.write(reader.replace(old, "if (t === undefined) { t = [PANEL]; }"))
        control("a no-act page treating every unclassified finding as closable", not fail_closed_problems(b, base, out, og, ["no-act"]))
    finally:
        with open(path, "w") as fh:
            fh.write(reader)
    try:
        old = "openFacts(qc)[0].filter(function (f) { return /depth|contamination|sex from/i.test(f[0]); })"
        assert reader.count(old) == 1
        with open(path, "w") as fh:
            fh.write(reader.replace(old, "(qc.facts || []).filter(function (f) { return /depth|contamination|sex from/i.test(f[0]); })"))
        control("trust line reading QC facts outside their veil", not fail_closed_problems(b, base, out, og, ["trust"]))
    finally:
        with open(path, "w") as fh:
            fh.write(reader)
    for old, new, desc in (
            ("(sec.n === 9 || sec.n === 10) && tagClosed(PANEL)", "false",
             "report metadata read outside PANEL"),
            ("sec.n === 9 && tagClosed(PANEL)", "false",
             "unclassified report metadata counted in the contents")):
        try:
            assert reader.count(old) == 1
            with open(path, "w") as fh:
                fh.write(reader.replace(old, new))
            control(desc, not fail_closed_problems(b, base, out, og, ["report-metadata"]))
        finally:
            with open(path, "w") as fh:
                fh.write(reader)
    print(f"   fail-closed producer checks: {time.monotonic() - started:.2f} s")


def infer_checks(b, base, out):
    """The son's page is a function of his open results and the count of his closed ones: two pages that differ only
    in which result is closed are the same before any click (INFER), and a closed deletion's second gene is named
    nowhere. Controls: person.js that puts a closed card's id back, that draws a fixed-panel card's open rows, or
    that counts the closed rows; a veil that tags only a deletion's first gene."""
    og = og_binary() or build_live_demo.OG
    if not os.path.isfile(os.path.join(out, "runs", OTHERS[1], "map.json")):
        print("[SKIP] infer checks: the demo has no map.json")
        return
    t0 = time.time()
    gate_regression_checks(b, base, out, og)
    bad = infer_problems(b, base, out, og)
    check("browser: two pages of the son that differ only in which result is closed are the same before any click, "
          "HTML, attributes, address and map canvas (" + "; ".join(x[0] for x in INFER) + ")", not bad, bad)
    bad = click_problems(b, base, out, og)
    check("browser: two pages of the son that differ only in which result is closed are the same after the same "
          "result is opened through its handle (" + "; ".join(x[0] for x in CLICKS) + ")", not bad, bad)
    bad = two_closed_problems(b, base, out, og)
    check("browser: a place that names two closed results stays closed while either is: opening one never shows the "
          "other's gene (" + "; ".join(x[0] for x in TWO_CLOSED) + ")", not bad, bad)
    bad = deep_link_problems(b, base)
    check("browser: #f/closed-1 and a guessed #f/<id> scroll to the first closed result, in view, its section marked in "
          "the contents", not bad, bad)
    bad = sv_gene_problems(b, base, out, og)
    check("browser: a closed deletion over two genes closes the rows that name only its second gene", not bad, bad)
    bad = panel_words_problems(b, base, out, og)
    check("browser: a fixed-panel card whose own finding names the closed result's gene shows the panel's words", not bad,
          bad)
    bad = address_problems(b, base, out)
    check("browser: opening a closed result and moving on never puts its id in the address or a history entry", not bad,
          bad)
    bad = home_problems(b, base, out)
    check("browser: the home page is the same whichever of the son's results is closed (" + "; ".join(x[0] for x in HOME)
          + ")", not bad, bad)
    home = os.path.join(out, "assets", "home.js")
    with open(home) as fh:
        keep_home = fh.read()
    a = "x.r.url + (homeShut(x.r, x.f) ? '' : 'finding/' + encodeURIComponent(x.f.id) + '/')"
    try:
        if a not in keep_home:
            check("control: home.js that links a closed result's own page (the planted line is in home.js)", False, a)
        else:
            with open(home, "w") as fh:
                fh.write(keep_home.replace(a, "x.r.url + 'finding/' + encodeURIComponent(x.f.id) + '/'"))
            control("home.js that links a closed result's own page", not home_problems(b, base, out, HOME[:1]))
    finally:
        with open(home, "w") as fh:
            fh.write(keep_home)
    first = live_results.panel_finding
    try:
        live_results.panel_finding = lambda f: False
        control("a veil that tags a fixed-panel card's own finding by the closed result's gene",
                not infer_problems(b, base, out, og, [x for x in INFER if "(NAT2)" in x[0]]))
    finally:
        live_results.panel_finding = first
    prose, two_names = live_results.prose_tag, live_results.two_hit_names
    try:  # the word rule back: prose tagged by its words read as gene symbols
        live_results.prose_tag = lambda spans, text: live_results.veil_tags(
            spans, re.findall(r"[A-Za-z0-9][A-Za-z0-9-]*", str(text).upper()), text)
        control("a veil that tags the carrier list's two-hit callers line by its words",
                not infer_problems(b, base, out, og, [x for x in INFER if "two-hit" in x[0]]))
        live_results.two_hit_names = lambda S: {}
        control("a veil that tags the two-hit gene lines by their words (AR)", not click_problems(b, base, out, og, CLICKS[:1]))
    finally:
        live_results.prose_tag = prose
        live_results.two_hit_names = two_names
    tags_, union_ = live_results.veil_tags, live_results.union
    try:  # one id per place, as before: the first closed result a place names
        live_results.veil_tags = lambda *a, **k: (tags_(*a, **k) or [None])[:1] if tags_(*a, **k) else None
        live_results.union = lambda *t: next((x[:1] for x in t if x), None)
        control("a veil that keeps one id per place (the lead's segment)", not two_closed_problems(b, base, out, og, TWO_CLOSED[:1]))
        control("a veil that keeps one id per place (the duplication)", not two_closed_problems(b, base, out, og, TWO_CLOSED[2:3]))
    finally:
        live_results.veil_tags, live_results.union = tags_, union_
    js = os.path.join(out, "assets", "person.js")
    with open(js) as fh:
        keep_js = fh.read()
    try:
        for desc, a, z, case in (
                # The guard is the pair of lines: either alone lands every guess the same way.
                ("person.js that picks a closed card by a guessed id",
                 "    if (f && !(veiled(f) && !m)) { return f; }\n    return closedOrder().length ? findingById(closedOrder()[0]) : null;",
                 "    return f || null;", [x for x in INFER if x[0].startswith("the same guessed")][0]),
                ("person.js that picks a closed card by a guessed id, at load",
                 "    if (f && !(veiled(f) && !m)) { return f; }\n    return closedOrder().length ? findingById(closedOrder()[0]) : null;",
                 "    return f || null;", [x for x in INFER if x[0].startswith("one guessed")][0]),
                ("person.js that puts a closed card's id back",
                 "function handle(f) { return veiled(f) ? 'closed-' + (closedOrder().indexOf(f.id) + 1) : f.id; }",
                 "function handle(f) { return f.id; }", INFER[0]),
                ("person.js that draws a fixed-panel card's open rows", "} else if (panelShut(c)) {", "} else if (false) {",
                 INFER[2]),
                ("person.js that counts the closed rows",
                 "'Some of its rows or numbers stay closed with a result to act on, until that result is opened.'",
                 "(ot[1] + (c.lines || []).length - lines.length) + ' rows stay closed with a result to act on.'", INFER[4])):
            if a not in keep_js:
                check(f"control: {desc} (the planted line is in person.js)", False, a)
                continue
            with open(js, "w") as fh:
                fh.write(keep_js.replace(a, z))
            control(desc, not infer_problems(b, base, out, og, [case]))
        for desc, a, z, still in (
                ("person.js whose fixed-panel card says 'This section holds a result' first", "} else if (panelShut(c)) {",
                 "} else if (panelShut(c) && !closedFs.length) {", lambda: not panel_words_problems(b, base, out, og)),
                ("person.js that reads only the first id of a tag", "if (Array.isArray(id)) { return id.some(tagClosed); }",
                 "if (Array.isArray(id)) { return tagClosed(id[0]); }",
                 lambda: not two_closed_problems(b, base, out, og, TWO_CLOSED[2:3])),
                ("person.js that writes an opened result's id into the address",
                 "function addr(f) { return !closable(f) ? f.id : veiled(f) ? handle(f) : 'opened-' + (openedOrder().indexOf(f.id) + 1); }",
                 "function addr(f) { return veiled(f) ? handle(f) : f.id; }", lambda: not address_problems(b, base, out))):
            if a not in keep_js:
                check(f"control: {desc} (the planted line is in person.js)", False, a)
                continue
            with open(js, "w") as fh:
                fh.write(keep_js.replace(a, z))
            control(desc, still())
    finally:
        with open(js, "w") as fh:
            fh.write(keep_js)
    first = live_results.closable
    try:
        live_results.closable = lambda found: first([{k: v for k, v in f.items() if k != "genes"} for f in found])
        control("a veil that tags only a deletion's first gene", not sv_gene_problems(b, base, out, og))
    finally:
        live_results.closable = first
    cols, names, marks = live_results.TABLE_COLUMNS, live_results.finding_names, live_results.mark_names
    try:
        live_results.TABLE_COLUMNS = dict(cols, Inheritance="gene")
        control("a veil that reads a table's Inheritance column as naming a gene",
                not infer_problems(b, base, out, og, [x for x in INFER if x[0].startswith("an inheritance code")]))
        live_results.TABLE_COLUMNS = cols
        live_results.finding_names = lambda f: (re.findall(r"[A-Z0-9]+", f"{f.get('title') or ''}"), f.get("detail") or "")
        control("a veil that reads a finding's title word by word",
                not infer_problems(b, base, out, og, [x for x in INFER if "(INS)" in x[0]]))
        live_results.finding_names = names
        live_results.mark_names = lambda S, sp, m: (re.findall(r"[A-Z0-9]+", f"{m.get('label') or ''} {m.get('detail') or ''}"), "")
        control("a veil that reads a chart mark's label word by word",
                not infer_problems(b, base, out, og, [x for x in INFER if "(INS)" in x[0]]))
    finally:
        live_results.TABLE_COLUMNS, live_results.finding_names, live_results.mark_names = cols, names, marks
    og_map = build_live_demo.og_map

    def blind(og_, args, body):
        """og map without results.veil.findings: it no longer knows which open findings a closed result ties."""
        d = json.loads(body)
        (d["results"].get("veil") or {}).pop("findings", None)
        return og_map(og_, args, json.dumps(d).encode())
    try:
        build_live_demo.og_map = blind
        control("og map that counts a finding a closed result ties among the open findings with no place",
                not infer_problems(b, base, out, og, [x for x in INFER if "(OPHN1 has no place" in x[0]]))
    finally:
        build_live_demo.og_map = og_map
    with open(js) as fh:
        keep_js = fh.read()
    route = "{ var hf = hashFinding(); target = document.getElementById("
    try:
        if route not in keep_js:
            check("control: person.js whose #f/closed-<n> lands nowhere (the planted line is in person.js)", False, route)
        else:
            with open(js, "w") as fh:
                fh.write(keep_js.replace(route, "{ var hf = null; h = 'f/nowhere'; target = document.getElementById("))
            control("person.js whose #f/closed-<n> lands nowhere", not deep_link_problems(b, base))
    finally:
        with open(js, "w") as fh:
            fh.write(keep_js)
    print(f"   infer checks: {time.time() - t0:.0f} s")


def widget_veil_checks(b, base, out):
    """On the son's page the copy-number card and the CHIP card, whose rows are a fixed panel (the same genes for
    everyone), close whole while any result of his to act on is closed: no dot, no call, no fact of theirs shows,
    and their line on a phone is the same words whatever the result; so do the ClinVar card's lines for the known
    deep intronic sites, a fixed panel inside it. Each "Show the result" of their box opens one result; once every
    result is opened every dot, call, fact and deep intronic line is back."""
    no_cache(b)
    rj = os.path.join(out, "runs", OTHERS[1], "results.json")
    js = os.path.join(out, "assets", "person.js")
    with open(rj) as f:
        keep = f.read()
    with open(js) as f:
        keep_js = f.read()

    def drawn():
        b.go("about:blank")
        b.wait("location.href === 'about:blank'")
        b.go(f"{base}/runs/{OTHERS[1]}/")
        b.wait("document.documentElement.dataset.pgpResults && document.documentElement.dataset.pgpStatus === 'ready' "
               "&& document.querySelector('#card-paralogs') && document.querySelector('.veil') && 'drawn'", 15)
        return shown()

    def phone_line():
        """The copy-number card's one line on a phone (its fold's summary, after the card's own title)."""
        b.call("Emulation.setDeviceMetricsOverride", {"width": 390, "height": 844, "deviceScaleFactor": 2, "mobile": False},
               b.session)
        try:
            b.go("about:blank")
            b.wait("location.href === 'about:blank'")
            b.go(f"{base}/runs/{OTHERS[1]}/")
            return b.wait("document.documentElement.dataset.pgpStatus === 'ready' && document.documentElement.dataset.pgpResults "
                          "&& (document.querySelector('[data-card=fold-paralogs] summary .what') || {}).textContent", 15) or ""
        finally:
            b.call("Emulation.setDeviceMetricsOverride", {"width": 1440, "height": 900, "deviceScaleFactor": 1, "mobile": False},
                   b.session)

    def phone_problems():
        line = phone_line()
        want = "The same genes for everyone, closed while a result about a child is closed."
        return [] if line == want else [f"the copy-number card's line on a phone: {line!r}, want {want!r}"]

    def shown():
        return b.eval("JSON.stringify([" + SHOWN_JS + ", (document.querySelector('[data-widget=copies]') || {}).textContent || '', "
                      "document.querySelectorAll('.finding[data-level=act] .veil').length, "
                      "document.querySelectorAll('#card-paralogs tbody tr,#card-chip tbody tr').length])")

    def problems(text, opened):
        page, cops, _, rows = json.loads(text) if text else ["", "", 0, 0]
        out_ = []
        deep = page.count("Known deep intronic site (always reported):")
        if opened and not ("SMN1" in cops and "SMN2" in cops and "p.R175H" in page and "SMN1 copies" in page
                       and deep == sites[0]):
            out_.append(f"every result opened, yet a dot, the call, the fact or a deep intronic line ({deep} of {sites[0]}) is "
                        f"missing: {cops[:200]}")
        if not opened and (rows or cops or "p.R175H" in page or "SMN1 copies" in page or deep):
            out_.append(f"a closed result, yet the copy-number or CHIP card shows its rows or the ClinVar card "
                        f"{deep} deep intronic lines: {cops[:200]}")
        return out_
    try:
        with open(os.path.join(os.path.dirname(out), "root", OTHERS[1], "summary.json")) as f:
            summary = json.load(f)
        S = summary.setdefault("sections", {})
        S["cpsr"] = {"state": "ok", "data": {"call_floor_dp": 10, "low_confidence": 0, "plp": [
            {"gene": "TP53", "variant": "c.524G>A", "final": "Pathogenic", "clinvar": "Pathogenic"},
            {"gene": "SMN1", "variant": "c.1A>G", "final": "Pathogenic", "clinvar": "Pathogenic"}]}}
        S["chip"] = {"state": "ok", "data": {"median_limit": 0.05, "median_depth": 72, "calls": [
            {"gene": g, "protein_change": ch, "vaf": v, "alt_reads": 6, "depth": 75, "status": "pass", "note": "a clone"}
            for g, ch, v in (("DNMT3A", "p.R882H", 0.12), ("TP53", "p.R175H", 0.08))],
            "genes": [{"gene": g, "detection_limit_vaf": 0.05, "mean_depth": 75, "status": "ok"} for g in ("DNMT3A", "TP53")]}}
        sites = [deep_sites(S)]
        S["paralogs"] = {"state": "ok", "data": {"genes": [{"gene": g, "copy_number": n, "status": "called", "detail": "."}
                                                           for g, n in (("SMN1", 2), ("SMN2", 2))]}}
        res = live_results.build(OTHERS[1], summary)
        with open(rj, "w") as f:
            json.dump(res, f)
        bad = problems(drawn(), False) + phone_problems()
        if not bad:
            # The copy-number card's box opens one result per click, until none is left closed.
            n = json.loads(shown())[2]
            for k in range(n):
                b.eval("var v = document.querySelector('#card-paralogs .veil button'); v && v.click(); true")
                time.sleep(0.3)
                left = json.loads(shown())[2]
                if left != n - k - 1:
                    bad.append(f"click {k + 1} on the copy-number card's box left {left} of {n} results closed")
                    break
            bad += problems(shown(), True)
        check("browser: on the son's page the copy-number and CHIP cards, a fixed panel, close whole while a result is "
              "closed (their line on a phone too), and so do the ClinVar card's deep intronic lines; the box opens one "
              "result per click, and every dot, call, fact and line is back once every result is opened", not bad, bad)
        for desc, a, z in (
                ("person.js that draws a fixed-panel card's open rows", "} else if (panelShut(c)) {", "} else if (false) {"),
                ("person.js whose card line reads a fixed-panel card's facts", ": g === 'show' && panelShut(c) ? PANEL_LINE",
                 ": false ? PANEL_LINE"),
                ("person.js that draws a fixed panel's lines", "if (id === PANEL) { return closedPage() && closedActs().length > 0; }",
                 "if (id === PANEL) { return false; }")):
            if a not in keep_js:
                check(f"control: {desc} (the planted line is in person.js)", False, a)
                continue
            with open(js, "w") as f:
                f.write(keep_js.replace(a, z))
            control(desc, not (problems(drawn(), False) + phone_problems()))
            with open(js, "w") as f:
                f.write(keep_js)
    finally:
        with open(rj, "w") as f:
            f.write(keep)
        with open(js, "w") as f:
            f.write(keep_js)


# finding.js as it was: a missing state or one with no plugins read as a standalone adult.
OLD_FINDING = ("if (!st || !st.plugins) { return Promise.resolve(UNKNOWN); }\n    var fam = st.plugins.family;",
               "var fam = st && st.plugins && st.plugins.family;")


def closed_pages(b, base, fid):
    """The son's person page and one planted finding's page, before any click: their text once each has drawn
    (or after 10 s, for a page whose state never comes)."""
    b.go("about:blank")
    b.wait("location.href === 'about:blank'")
    b.go(f"{base}/runs/{OTHERS[1]}/")
    b.wait("document.documentElement.dataset.pgpResults && document.documentElement.dataset.pgpStatus === 'ready' "
           "&& 'drawn'", 10)
    time.sleep(0.5)
    person = b.eval("document.body.textContent") or ""
    b.go(f"{base}/runs/{OTHERS[1]}/finding/{fid}/")
    b.wait("document.documentElement.dataset.pgpStatus === 'ready' && 'drawn'", 10)
    page = b.eval("document.body.textContent") or ""
    return person, page


def closed_problems(person, page):
    out = [f"{w} on the person page" for pair in VEIL for w in pair if w in person]
    out += [f"{w} on the finding page" for pair in VEIL for w in pair if w in page]
    if "A result to act on" not in page:
        out.append("the finding page does not show its closed title: " + page[:200])
    return out


def search_text(b, base, q):
    """[hits, the search panel's text] for q typed on the son's page."""
    b.go("about:blank")
    b.wait("location.href === 'about:blank'")
    b.go(f"{base}/runs/{OTHERS[1]}/")
    b.wait("document.documentElement.dataset.pgpResults && document.querySelector('.veil') && "
           "!document.querySelector('input[type=search]').disabled && 'ready'", 15)
    b.eval("var i = document.querySelector('input[type=search]'); i.value = %s; i.dispatchEvent(new Event('input')); true"
           % json.dumps(q))
    hits = b.wait("document.documentElement.dataset.pgpSearch", 5)
    return hits, b.eval("(document.querySelector('.search-results') || {}).textContent || ''") or ""


def search_problems(b, base):
    out = []
    hits, text = search_text(b, base, "GENEQ")
    if hits != "0" or "c.77del" in text:
        out.append(f"GENEQ, a gene of a closed card: {hits} hits, {text[:200]}")
    hits, text = search_text(b, base, "GENE")
    if "GENEQ" in text or "GENER" in text or "GENEOPEN" not in text:
        out.append(f"the suggestions for GENE: {text[:300]}")
    hits, text = search_text(b, base, "GENEOPEN")
    if hits != "1":
        out.append(f"GENEOPEN, a gene of an open card, not found: {hits} {text[:200]}")
    return out


def unknown_checks(b, base, out, fid):
    """The son's results to act on stay closed when the page cannot learn his family role, on the person page
    and on the finding's own page, and his search names nothing his closed card holds."""
    sd = os.path.join(out, "runs", OTHERS[1])
    fj = os.path.join(out, "runs", FAM, "family.json")
    frames_path = os.path.join(sd, "frames.json")
    idx = os.path.join(sd, f"{OTHERS[1]}_search_index.json")
    pj, fjs, sj = (os.path.join(out, "assets", x) for x in ("person.js", "finding.js", os.path.join("plugins", "search.js")))
    keep = {}
    for p in (frames_path, pj, fjs, sj) + ((idx,) if os.path.isfile(idx) else ()):
        with open(p) as f:
            keep[p] = f.read()
    frames = json.loads(keep[frames_path])

    def write(p, text):
        with open(p, "w") as f:
            f.write(text)

    def item(sec, text):
        return {"section": sec, "anchor": sec, "title": sec, "step": "24", "state": "ok", "source": "summary.json",
                "count": 1, "mentions": [{"text": text, "locus": None}]}
    try:
        # The search: an index naming the two planted genes in the closed CPSR card, and one gene of an open card.
        write(idx, json.dumps({"schema_version": 1, "sample": OTHERS[1], "terms": {
            "GENEQ": [item("cpsr", "GENEQ c.77del")], "GENER": [item("cpsr", "GENER c.88dup")],
            "GENEOPEN": [item("roh", "GENEOPEN in a run of homozygosity")]}}))
        with_search = copy.deepcopy(frames)
        for fr in with_search["frames"]:
            fr["state"].setdefault("plugins", {})["search"] = {"file": os.path.basename(idx), "mtime": 1}
        write(frames_path, json.dumps(with_search))
        bad = search_problems(b, base)
        check("browser: the son's search lists no hit and suggests no term of his closed card, and still finds an open one",
              not bad, bad)
        write(sj, keep[sj].replace("function hidden(it) { return !!(", "function hidden(it) { return false && !!("))
        control("search.js without its filter", not search_problems(b, base))
        write(sj, keep[sj])
        write(frames_path, keep[frames_path])

        # family.json unavailable
        os.replace(fj, fj + ".away")
        try:
            bad = closed_problems(*closed_pages(b, base, fid))
        finally:
            os.replace(fj + ".away", fj)
        check("browser: with family.json unavailable the son's results to act on stay closed on both pages", not bad, bad)

        # a state cut short: no plugins at all
        cut = copy.deepcopy(frames)
        for fr in cut["frames"]:
            fr["state"].pop("plugins", None)
        write(frames_path, json.dumps(cut))
        bad = closed_problems(*closed_pages(b, base, fid))
        check("browser: with a state that has no plugins the son's results to act on stay closed on both pages", not bad, bad)
        write(pj, keep[pj].replace("return !state || !state.plugins || !!(fam && fam.unknown);", "return false;"))
        write(fjs, keep[fjs].replace(OLD_FINDING[0], OLD_FINDING[1]))
        control("person.js and finding.js that read a state with no plugins as a standalone adult",
                not closed_problems(*closed_pages(b, base, fid)))
        write(pj, keep[pj])
        write(fjs, keep[fjs])

        # no state at all
        os.remove(frames_path)
        bad = closed_problems(*closed_pages(b, base, fid))
        check("browser: with no state at all (frames.json unavailable) the son's results to act on stay closed on both pages",
              not bad, bad)
        write(fjs, keep[fjs].replace(OLD_FINDING[0], OLD_FINDING[1]))
        control("finding.js that reads a missing state as a standalone adult", not closed_problems(*closed_pages(b, base, fid)))
    finally:
        for p, text in keep.items():
            write(p, text)
        if idx not in keep and os.path.isfile(idx):
            os.remove(idx)


# The person page reads as a report: what (frame) -> the page's layout facts, read in the browser.
LAYOUT_JS = """JSON.stringify((function () {
  var main = document.getElementById('content'), fold = document.getElementById('progress-fold');
  var band = document.getElementById('runhero'), lead = document.getElementById('lead-slot');
  var wire = Array.prototype.map.call(document.querySelectorAll('#found > li[role=button]'), function (li) {
    return (li.childNodes[1] || {}).textContent || ''; });
  var hla = document.getElementById('card-hla');
  var shownLines = hla ? hla.querySelectorAll('.card-text .lines-box > ul.lines > li').length : -1;
  return {
    last: main.lastElementChild && main.lastElementChild.id, folded: !!fold && !fold.open,
    band: !band.hidden, bandFirst: !!(band.compareDocumentPosition(lead) & Node.DOCUMENT_POSITION_FOLLOWING),
    wire: wire, hla: shownLines, hlaFold: !!(hla && hla.querySelector('details.lines-more')),
    margin: !!(hla && hla.querySelector('.card-margin .caveat .caveat-lead'))
  };
}()))"""


def person_layout(b, base, frame):
    """The layout facts of DEMO-001's page at a replay frame ('mid' or 'end')."""
    b.go("about:blank")
    b.wait("location.href === 'about:blank'")
    b.go(f"{base}/runs/{S}/")
    b.wait("document.documentElement.dataset.pgpStatus === 'ready' && document.documentElement.dataset.pgpFrame")
    b.eval("var s = document.getElementById('replay-slider'); s.value = %s; s.dispatchEvent(new Event('input')); true"
           % ("Math.floor(s.max / 2)" if frame == "mid" else "s.max"))
    b.wait("document.documentElement.dataset.pgpResults && document.documentElement.dataset.pgpFrame === %s"
           % ("String(Math.floor(document.getElementById('replay-slider').max / 2))" if frame == "mid"
              else "document.getElementById('replay-slider').max"))
    time.sleep(0.5)
    return json.loads(b.eval(LAYOUT_JS) or "{}")


def layout_problems(mid, end):
    out = []
    for name, got in (("mid", mid), ("end", end)):
        if got.get("last") != "progress" or not got.get("folded"):
            out.append(f"{name}: the step timeline is not the last block, folded: {got.get('last')}, folded {got.get('folded')}")
        wire = got.get("wire") or []
        if not wire or len(wire) > 8:
            out.append(f"{name}: the wire shows {len(wire)} lines, want 1 to 8")
        plain = [w for w in wire if "finished, see the report" in w]
        if len(plain) > 1:
            out.append(f"{name}: {len(plain)} plain 'finished' lines on the wire, want them folded into one")
        long = [w for w in wire if len(w) > 141]
        if long:
            out.append(f"{name}: wire lines longer than one sentence: {long[:2]}")
    if not mid.get("band") or not mid.get("bandFirst"):
        out.append(f"mid: the run band is not shown above the lead: {mid.get('band')}, {mid.get('bandFirst')}")
    if end.get("band"):
        out.append("end: the run band still shows on a finished run")
    if not (0 <= end.get("hla", -1) <= 6) or not end.get("hlaFold"):
        out.append(f"end: the HLA card shows {end.get('hla')} lines, the absent ones not folded ({end.get('hlaFold')})")
    if not end.get("margin"):
        out.append("end: the HLA card's caveat is not in its margin")
    return out


def layout_checks(b, base, out):
    """The person page answers 'is there anything serious?' first: a short band while the run goes, the step
    timeline last and folded, a wire of at most eight one-sentence lines with the plain finishes folded into one,
    a card's caveat in its margin and its absent or undetermined lines folded."""
    bad = layout_problems(person_layout(b, base, "mid"), person_layout(b, base, "end"))
    check("browser: the person page leads with the report, a short run band while it runs, the step timeline last "
          "and folded, a short wire, and each card's caveat in its margin with its quiet lines folded", not bad, bad)
    b.call("Network.enable", session=b.session)  # the controls below change live.js and person.js: never cached
    b.call("Network.setCacheDisabled", {"cacheDisabled": True}, b.session)
    files = {n: os.path.join(out, "assets", n) for n in ("live.js", "person.js")}
    keep = {n: open(f).read() for n, f in files.items()}
    planted = (("live.js", "said.slice(0, WIRE_SHOWN)", "said", "a wire that shows every line"),
               ("live.js", "if (m && hl[k].tone !== 'attention')", "if (false)", "a wire that lists every plain finish"),
               ("person.js", "var LINES_SHOWN = 6;", "var LINES_SHOWN = 99;", "a card that shows every line"),
               ("person.js", "folded = true; fold.open = ", "folded = true; fold.open = true || ", "a step timeline left open"))
    try:
        for name, a, z, desc in planted:
            if a not in keep[name]:
                check(f"control: {desc} (the planted line is in {name})", False, a)
                continue
            with open(files[name], "w") as f:
                f.write(keep[name].replace(a, z))
            control(desc, not layout_problems(person_layout(b, base, "mid"), person_layout(b, base, "end")))
            with open(files[name], "w") as f:
                f.write(keep[name])
    finally:
        for n, f in files.items():
            with open(f, "w") as fh:
                fh.write(keep[n])


def too_wide(b):
    return b.eval("document.documentElement.scrollWidth > window.innerWidth ? "
                  "document.documentElement.scrollWidth + ' > ' + window.innerWidth : ''")


def phone_checks(b, base):
    # mobile False: a mobile viewport zooms out to fit a page that overflows, which would hide the overflow.
    b.call("Emulation.setDeviceMetricsOverride", {"width": 390, "height": 844, "deviceScaleFactor": 2, "mobile": False},
           b.session)
    try:
        wide = {}
        for name, path, ready in (("home", "/", "document.documentElement.dataset.pgpStatus === 'ready'"),
                                  ("person", f"/runs/{S}/", "document.documentElement.dataset.pgpResults"),
                                  ("family", f"/runs/{FAM}/", "document.documentElement.dataset.pgpStatus === 'ready'"),
                                  ("partner", f"/runs/{OTHERS[0]}/", "document.documentElement.dataset.pgpResults"),
                                  ("finding", f"/runs/{OTHERS[0]}/finding/str-fmr1/",
                                   "document.documentElement.dataset.pgpStatus === 'ready'")):
            b.go(base + path)
            b.wait(ready)
            if name == "person":
                b.eval("var s = document.getElementById('replay-slider'); s.value = s.max; "
                       "s.dispatchEvent(new Event('input')); true")
                b.wait("document.documentElement.dataset.pgpFrame === document.getElementById('replay-slider').max")
                time.sleep(1)
                # Every section open: the cards folded into section 8 too.
                b.eval("Array.prototype.forEach.call(document.querySelectorAll('main details'), "
                       "function (d) { d.open = true; }); true")
                time.sleep(0.3)
                tabs = b.eval("getComputedStyle(document.getElementById('tabbar')).display")
                wide["person tab bar"] = "" if tabs == "grid" else f"tab bar display {tabs!r}"
            wide[name] = too_wide(b)
        check("browser: at 390 px no page scrolls sideways", not any(wide.values()), wide)
        b.eval("var d = document.createElement('div'); d.style.width = '900px'; d.style.height = '4px'; "
               "document.querySelector('main').appendChild(d); true")
        control("a 900 px wide block on the family page", not too_wide(b))
    finally:
        b.call("Emulation.clearDeviceMetricsOverride", {}, b.session)


PHONE_JS = r"""JSON.stringify((function () {
  // The content of a closed details keeps a layout box (content-visibility: hidden): ask whether it shows.
  var shown = function (e) { return !!e && e.checkVisibility() && e.getBoundingClientRect().height > 0; };
  var loud = Array.prototype.filter.call(document.querySelectorAll('#results [data-finding]'), function (e) {
    return /^(act|medicine|carrier|review)$/.test(e.dataset.level); });
  return { height: document.documentElement.scrollHeight, loud: loud.length,
    hidden: loud.filter(function (e) { return !shown(e); }).map(function (e) { return e.dataset.finding; }),
    deep: shown(document.getElementById('card-' + %s)), // section 7 open: its cards' own one-line folds show
    sec: shown(document.querySelector('#s7 details.secfold details.pfold > summary')),
    folds: document.querySelectorAll('#results details.pfold').length };
}()))"""
PHONE_MAX = 8000  # CSS px: the Journal design's phone page is about 4,600; every card open was about 37,000
PHONE_DEEP = "roh"  # a card of section 6: on a phone it sits in its own fold inside the section's fold


def phone_report(b, base):
    """DEMO-002's page at 390 px: its height, its findings to act on, a deep link to a card two folds deep and
    one to a folded section, each before and after the link."""
    b.call("Emulation.setDeviceMetricsOverride", {"width": 390, "height": 844, "deviceScaleFactor": 2, "mobile": False},
           b.session)
    try:
        b.go("about:blank")
        b.wait("location.href === 'about:blank'")
        b.go(f"{base}/runs/{OTHERS[0]}/")
        b.wait("document.documentElement.dataset.pgpResults")
        time.sleep(1)
        js = PHONE_JS % json.dumps(PHONE_DEEP)
        got = {"before": json.loads(b.eval(js) or "{}")}
        b.eval("location.hash = '#s/%s'; true" % PHONE_DEEP)
        time.sleep(0.5)
        got["card"] = json.loads(b.eval(js) or "{}")
        b.eval("location.hash = '#s7'; true")
        time.sleep(0.5)
        got["section"] = json.loads(b.eval(js) or "{}")
        return got
    finally:
        b.call("Emulation.clearDeviceMetricsOverride", {}, b.session)


def phone_report_problems(got):
    out = []
    pre = got.get("before") or {}
    if not pre.get("height") or pre["height"] >= PHONE_MAX:
        out.append(f"the page is {pre.get('height')} CSS px tall at 390 px, want under {PHONE_MAX}")
    if not pre.get("loud") or pre.get("hidden"):
        out.append(f"findings to act on, for a prescriber, for the family or worth a look closed behind a tap: "
                   f"{pre.get('hidden')} of {pre.get('loud')}")
    if pre.get("deep") or pre.get("sec"):
        out.append(f"card {PHONE_DEEP} or section 7 starts open on a phone: {pre.get('deep')}, {pre.get('sec')}")
    if not (got.get("card") or {}).get("deep"):
        out.append(f"#s/{PHONE_DEEP} leaves its card closed")
    if not (got.get("section") or {}).get("sec"):
        out.append("#s7 leaves section 7 closed")
    return out


def phone_report_checks(b, base, out):
    """On a phone the report is short: the findings that matter open, the rest one tap away, and a link opens
    the folds around its target. Wider, the page has no fold at all."""
    bad = phone_report_problems(phone_report(b, base))
    check(f"browser: at 390 px {OTHERS[0]}'s page is under {PHONE_MAX} CSS px with every finding to act on, for a "
          "prescriber, for the family or worth a look open, and a link to a folded card or section opens it", not bad, bad)
    b.go(f"{base}/runs/{OTHERS[0]}/")
    b.wait("document.documentElement.dataset.pgpResults")
    time.sleep(0.5)
    folds = b.eval("document.querySelectorAll('#results details.pfold').length")
    check("browser: wider than a phone the person page folds no card and no section", folds == 0, folds)
    b.call("Network.enable", session=b.session)
    b.call("Network.setCacheDisabled", {"cacheDisabled": True}, b.session)
    f = os.path.join(out, "assets", "person.js")
    keep = open(f).read()
    planted = (("var phoneQuery = window.matchMedia ?", "var phoneQuery = false ?", "a phone page with every card open"),
               ("d = d.parentElement && d.parentElement.closest('details')", "d = null", "a link that opens only the nearest fold"))
    try:
        for a, z, desc in planted:
            if a not in keep:
                check(f"control: {desc} (the planted line is in person.js)", False, a)
                continue
            with open(f, "w") as fh:
                fh.write(keep.replace(a, z))
            control(desc, not phone_report_problems(phone_report(b, base)))
    finally:
        with open(f, "w") as fh:
            fh.write(keep)


FLAGGED = "callability"  # a card of section 5 with no finding of its own, flagged for these checks


def phone_at(b, width):
    b.call("Emulation.setDeviceMetricsOverride", {"width": width, "height": 844, "deviceScaleFactor": 2, "mobile": False},
           b.session)


def fold_facts(b, base):
    """DEMO-002's page at 390 px with its section 5 card FLAGGED flagged: whether that card shows; whether a
    modified click (Ctrl) and a plain one on a link to the place the address names open its fold; and which
    folds the reader opened are still open after the page turns wide and back."""
    got = {}
    phone_at(b, 390)
    try:
        b.go("about:blank")
        b.wait("location.href === 'about:blank'")
        b.go(f"{base}/runs/{OTHERS[0]}/")
        b.wait("document.documentElement.dataset.pgpResults && document.querySelector('#results details.pfold')")
        time.sleep(0.5)
        got["flagged"] = b.eval("var e = document.getElementById('card-%s'); !!e && e.checkVisibility() && "
                                "e.getBoundingClientRect().height > 0" % FLAGGED)
        # The reader opens section 6's fold, the ROH card's fold in it and a card's other lines.
        got["opened"] = b.eval("""(function () {
          var s6 = document.querySelector('#s6 details.secfold'), roh = document.querySelector('details[data-card="fold-roh"]');
          var more = document.querySelector('#results article.result details.lines-more');
          if (!s6 || !roh || !more) { return ''; }
          s6.open = true; roh.open = true;
          for (var d = more; d; d = d.parentElement && d.parentElement.closest('details')) { d.open = true; }
          return more.closest('article').id; }())""")
        phone_at(b, 1280)
        b.wait("!document.querySelector('#results details.pfold') && 'wide'", 10)
        time.sleep(0.5)
        got["wide"] = b.eval("var a = document.getElementById(%s); !!a && !!a.querySelector('details.lines-more[open]')"
                             % json.dumps(got["opened"] or ""))
        phone_at(b, 390)
        b.wait("document.querySelector('#results details.pfold') && 'phone'", 10)
        time.sleep(0.5)
        got["back"] = b.eval("""JSON.stringify([!!document.querySelector('#s6 details.secfold[open]'),
          !!document.querySelector('details[data-card="fold-roh"][open]'),
          !!document.querySelector('#' + %s + ' details.lines-more[open]')])""" % json.dumps(got["opened"] or "x"))
        # A link to the place the address already names: Ctrl+click leaves its fold closed, a plain click opens it.
        b.eval("location.hash = '#s8'; true")
        b.wait("document.querySelector('#s8 details.secfold[open]') && 'open'", 5)
        b.eval("""window.addEventListener('click', function (e) { e.preventDefault(); });
          document.querySelector('#s8 details.secfold').open = false;
          document.querySelector('#tabbar a[href="#s8"]').dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true, ctrlKey: true}));
          true""")
        time.sleep(0.3)
        got["ctrl"] = b.eval("document.querySelector('#s8 details.secfold').open")
        b.eval("""document.querySelector('#tabbar a[href="#s8"]').dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true})); true""")
        time.sleep(0.3)
        got["plain"] = b.eval("document.querySelector('#s8 details.secfold').open")
        return got
    finally:
        b.call("Emulation.clearDeviceMetricsOverride", {}, b.session)


def fold_problems(got):
    out = []
    if not got.get("flagged"):
        out.append(f"the flagged section 5 card {FLAGGED} with no finding of its own is folded on a phone")
    if not got.get("opened"):
        out.append("no section 6 fold, ROH fold or card with other lines to open")
    if not got.get("wide"):
        out.append(f"the other lines opened in {got.get('opened')} closed when the page turned wide")
    if json.loads(got.get("back") or "null") != [True, True, True]:
        out.append(f"section 6, the ROH card and the other lines opened on the phone, open after wide and back: {got.get('back')}")
    if got.get("ctrl") is not False:
        out.append(f"a Ctrl+click on a link to the current place opened its fold: {got.get('ctrl')}")
    if got.get("plain") is not True:
        out.append(f"a plain click on a link to the current place left its fold closed: {got.get('plain')}")
    return out


def fold_checks(b, base, out):
    """On a phone a section 5 card the report flags with no finding of its own stays open with its section; a
    fold the reader opened stays open when the phone turns wide and back, inner folds too; a Ctrl+click on a
    link to the place the address names leaves the page as it is."""
    rj = os.path.join(out, "runs", OTHERS[0], "results.json")
    js = os.path.join(out, "assets", "person.js")
    keep, keep_js = open(rj).read(), open(js).read()
    b.call("Network.enable", session=b.session)
    b.call("Network.setCacheDisabled", {"cacheDisabled": True}, b.session)
    try:
        res = json.loads(keep)
        card = [c for ch in res["chapters"] for c in ch["sections"] if c["id"] == FLAGGED]
        sec5 = [x for x in res["journal"] if x["n"] == 5]
        check(f"browser: {FLAGGED} is a section 5 card with no finding of its own, so it can be flagged alone",
              bool(card and not card[0].get("findings") and sec5 and FLAGGED in sec5[0]["cards"]), (card, sec5))
        for c in card:
            c["attention"] = True
        with open(rj, "w") as f:
            json.dump(res, f)
        bad = fold_problems(fold_facts(b, base))
        check("browser: at 390 px a flagged section 5 card stays open, folds opened on a phone stay open after the "
              "page turns wide and back (other lines too), and a Ctrl+click leaves a fold as it is", not bad, bad)
        planted = (("function flaggedAlone(c) { return ", "function flaggedAlone(c) { return false && ",
                    "a flagged section 5 card folded with its section"),
                   ("function restoreFolds() { foldKeys()", "function restoreFolds() { return; foldKeys()",
                    "a redraw that forgets the folds"),
                   ("if (e.button !== 0 || e.ctrlKey", "if (false && e.ctrlKey", "a Ctrl+click that opens the fold"))
        for a, z, desc in planted:
            if a not in keep_js:
                check(f"control: {desc} (the planted line is in person.js)", False, a)
                continue
            with open(js, "w") as f:
                f.write(keep_js.replace(a, z))
            control(desc, not fold_problems(fold_facts(b, base)))
    finally:
        with open(rj, "w") as f:
            f.write(keep)
        with open(js, "w") as f:
            f.write(keep_js)


# --- the genome map (Figure 1, assets/live/map.js) ------------------------------------------

MARKS_TSV = os.path.join(REPO, "tests", "fixtures", "map", "karyogram-marks.tsv")


def og_binary():
    """An og that runs here: OG when set (Guard's unit-tests job and the docs builds set it to the og they
    built, as tests/og_binary.py reads it), else bin/og when scripts/setup.sh --og built it, else one cargo
    builds now from rust/ without the BAM feature (the map reads no BAM), kept in CARGO_TARGET_DIR or a folder
    of the temp dir. None when none of them works."""
    named = os.environ.get("OG")
    if named:
        if os.access(named, os.X_OK):
            return named
        print(f"OG={named} is not an executable og binary", file=sys.stderr)
        return None
    og = os.path.join(REPO, "bin", "og")
    if subprocess.run([og, "map", "--help"], capture_output=True, check=False).returncode == 0:
        return og
    cargo = shutil.which("cargo") or os.path.expanduser("~/.cargo/bin/cargo")
    if not os.path.isfile(cargo):
        return None
    target = os.environ.get("CARGO_TARGET_DIR") or os.path.join(tempfile.gettempdir(), "owngenome-test-og")
    p = subprocess.run([cargo, "build", "--locked", "--quiet", "-p", "og", "--no-default-features"], capture_output=True,
                       text=True, cwd=os.path.join(REPO, "rust"), env=dict(os.environ, CARGO_TARGET_DIR=target), check=False)
    exe = os.path.join(target, "debug", "og")
    if p.returncode != 0 or not os.path.isfile(exe):
        print(f"cargo could not build og ({p.returncode}): {p.stderr[-2000:]}", file=sys.stderr)
        return None
    return exe


def karyogram_marks():
    """(run, id, chrom, mb, level) of each mark the karyogram drew, recorded before the map replaced it."""
    rows = []
    with open(MARKS_TSV) as f:
        for line in f:
            if line.startswith(("#", "run\t")) or not line.strip():
                continue
            run, fid, c, m, lv = line.rstrip("\n").split("\t")
            rows.append((run, fid, "chr" + c, float(m), lv))
    return rows


def read_map(path):
    with open(path, "rb") as f:
        return json.loads(gzip.decompress(f.read()))


def write_map(path, doc):
    with open(path, "wb") as f:
        f.write(gzip.compress(json.dumps(doc, separators=(",", ":")).encode(), mtime=0))


# What the map shows: its state, the findings it draws (its list of marks), its timings, level and place,
# its chips, its margin's words and every text and title it holds, and the address.
MAP_JS = r"""JSON.stringify((function () {
  var f = document.getElementById('genome-map');
  if (!f) { return { state: 'absent', marks: [] }; }
  var t = document.querySelector('.gm-tip');
  return { state: f.dataset.pgpMap || '', level: f.dataset.level || '', chrom: f.dataset.chrom || '', view: f.dataset.view || '',
    parse: Number(f.dataset.pgpMapParseMs), paint: Number(f.dataset.pgpMapPaintMs),
    marks: Array.prototype.map.call(f.querySelectorAll('.gm-marks a'), function (a) { return { id: a.dataset.id,
      level: a.dataset.level, chrom: a.dataset.chrom, s: Number(a.dataset.s), e: Number(a.dataset.e) }; }),
    chips: Array.prototype.map.call(f.querySelectorAll('.gm-chip'), function (c) { return c.textContent; }),
    title: (document.getElementById('gm-title') || {}).textContent || '',
    layers: Array.prototype.map.call(f.querySelectorAll('.gm-layer'), function (l) { return [l.dataset.layer, l.dataset.state,
      l.querySelector('dd').textContent]; }),
    toast: (f.querySelector('.gm-toast') || { hidden: true }).hidden ? '' : f.querySelector('.gm-toast').textContent,
    tip: t && !t.hidden ? t.textContent : '',
    text: [f.textContent].concat(Array.prototype.map.call(f.querySelectorAll('[href],[title],[aria-label]'), function (e) {
      return (e.getAttribute('href') || '') + ' ' + (e.getAttribute('title') || '') + ' ' + (e.getAttribute('aria-label') || ''); })).join(' '),
    hash: location.hash };
}()))"""
MAP_READY = ("document.documentElement.dataset.pgpResults && (function () { var f = document.getElementById('genome-map'); "
             "return f && /^(drawn|failed|nodecompress)$/.test(f.dataset.pgpMap || '') "
             "&& (f.dataset.pgpMap !== 'drawn' || !!f.dataset.level); }())")


def map_page(b, base, run, hash_=""):
    """A fresh load of RUN's page with HASH_; the map's facts once it settled (None when it never drew)."""
    b.go("about:blank")
    b.wait("location.href === 'about:blank'")
    b.go(f"{base}/runs/{run}/{hash_}")
    if not b.wait(MAP_READY):
        return json.loads(b.eval(MAP_JS))
    time.sleep(0.6)  # the margin follows the camera 140 ms after it settles
    return json.loads(b.eval(MAP_JS))


def map_now(b):
    time.sleep(1.5)  # the camera eases for about 0.6 s, then the margin follows 140 ms after it settles
    return json.loads(b.eval(MAP_JS))


# The page's closed results to act on and the map's chips, read at once (a chip's button rebuilds the chips at
# once: nothing to wait for): the act cards in page order (ids: each one's id, its handle while it is closed) and
# the places of the closed ones among them (closed: a closed card names no id, so a result is known by its place), the act chip's count, whether its confirm step
# shows, where the keyboard is (first: the chips' first button; canvas: the map), and the opened chips.
CHIP_JS = r"""JSON.stringify((function () {
  var ids = [], closed = [];
  Array.prototype.forEach.call(document.querySelectorAll('.finding[data-level=act]'), function (a, i) {
    ids.push(a.dataset.finding);
    if (a.querySelector('.veil')) { closed.push(i); } });
  var c = document.querySelector('.gm-chip[data-closed=act]'), f = document.activeElement;
  var first = document.querySelector('.gm-chips button');
  return { closed: closed, ids: ids, buttons: !!first, chip: c ? c.firstChild.textContent : '', confirm: !!(c && c.querySelector('.gm-confirm')),
    focus: f && f.classList.contains('gm-yes') ? 'gm-yes' : f === first ? 'first' : f === document.querySelector('.gm-canvas') ? 'canvas'
      : (f ? f.tagName + '.' + f.className : ''),
    opened: Array.prototype.map.call(document.querySelectorAll('.gm-chip.opened'), function (o) { var a = o.querySelector('a');
      var raw = a ? a.getAttribute('href') : '';
      return { kind: o.dataset.opened, text: o.textContent, raw: raw, href: opened(raw) }; }) };
  // The finding an "Opened in its section" link names: #f/opened-<n> is the n-th finding the page can close and
  // shows opened, in the order of the sections (person.js openedOrder()).
  function opened(raw) {
    var m = /^#f\/opened-(\d+)$/.exec(raw);
    if (!m) { return raw; }
    if (!window.__pgpResults) {
      var x = new XMLHttpRequest(); x.open('GET', 'results.json', false); x.send(); window.__pgpResults = JSON.parse(x.responseText);
    }
    var r = window.__pgpResults, tags = (r.veil && r.veil.findings) || {}, byId = {}, order = [];
    r.findings.forEach(function (f) { byId[f.id] = f; });
    function ok(f) { return f && (f.level === 'act' || (tags[f.id] && tags[f.id] !== f.id)) && !!document.getElementById('f-' + f.id); }
    (r.journal || []).forEach(function (s) { (s.findings || []).forEach(function (id) { if (ok(byId[id]) && order.indexOf(id) < 0) { order.push(id); } }); });
    r.findings.forEach(function (f) { if (ok(f) && order.indexOf(f.id) < 0) { order.push(f.id); } });
    return '#f/' + (order[Number(m[1]) - 1] || raw);
  }
}()))"""


def chips_now(b):
    time.sleep(0.05)  # a chip's button draws the page and the chips again before its click returns
    return json.loads(b.eval(CHIP_JS))


def closed_words(n):
    return "" if n == 0 else "1 result closed" if n == 1 else f"{n} results closed"


def chip_open(b):
    """The act chip's Open, then its Show the result, as a keyboard reader does it. Returns (the id it opened, its
    place among the act cards, problems): one Open opens exactly one of the page's closed results, the count falls by one, the keyboard
    goes to Show the result and then to the chips' first button (the map when no chip is left), and an
    "Opened in its section" chip never names a result that is still closed."""
    before = chips_now(b)
    b.eval("var o = document.querySelector('.gm-chip[data-closed=act] button'); o.focus(); o.click(); true")
    mid = chips_now(b)
    bad = []
    if not mid["confirm"] or mid["closed"] != before["closed"]:
        bad.append(f"Open: confirm step {mid['confirm']}, closed {before['closed']} -> {mid['closed']}")
    if mid["focus"] != "gm-yes":
        bad.append(f"Open: the keyboard is on {mid['focus']!r}, not Show the result")
    b.eval("document.querySelector('.gm-chip[data-closed=act] .gm-yes').click(); true")
    after = chips_now(b)
    gone = [x for x in before["closed"] if x not in after["closed"]]
    n = len(before["closed"]) - 1
    if len(gone) != 1 or len(after["closed"]) != n:
        bad.append(f"Show the result opened {gone}, not one: closed {before['closed']} -> {after['closed']}")
    if not after["chip"].startswith(closed_words(n)) or bool(after["chip"]) != (n > 0):
        bad.append(f"Show the result: the chip reads {after['chip']!r}, want {closed_words(n)!r}")
    if after["focus"] != ("first" if after["buttons"] else "canvas"):
        bad.append(f"Show the result: the keyboard is on {after['focus']!r}")
    named = [o["href"] for o in after["opened"] if o["kind"] == "section" and o["href"][3:] not in after["ids"]]
    if named:
        bad.append(f"an Opened in its section chip names a result that is not open: {named}")
    raw = [o["raw"] for o in after["opened"] if o["kind"] == "section" and not re.match(r"^#f/opened-\d+$", o["raw"])]
    if raw:
        bad.append(f"an Opened in its section chip links a result by its id: {raw}")
    leaked = [after["ids"][i] for i in after["closed"] if not re.match(r"^closed-\d+$", after["ids"][i] or "")]
    if leaked:
        bad.append(f"a closed card carries an id, not a handle: {leaked}")
    if len(gone) != 1:
        return None, None, bad
    return after["ids"][gone[0]], gone[0], bad


def parity_problems(b, base, rows):
    """The karyogram's marks the map does not place with the same level and chromosome within 0.01 Mb."""
    bad = []
    for run in sorted({r[0] for r in rows}):
        got = map_page(b, base, run)
        if got["state"] != "drawn":
            bad.append(f"{run}: map {got['state']}")
            continue
        by = {m["id"]: m for m in got["marks"]}
        for _, fid, c, m, lv in (r for r in rows if r[0] == run):
            g = by.get(fid)
            if not g:
                bad.append(f"{run}: {fid} not on the map")
            elif g["level"] != lv or g["chrom"] != c or abs(g["s"] / 1e6 - m) > 0.01:
                bad.append(f"{run}: {fid} at {g['chrom']} {g['s'] / 1e6:.3f} Mb {g['level']}, want {c} {m} {lv}")
    return bad


def click_canvas(b, x, y, move_only=False):
    """A mouse at (x, y) of the map's canvas, in CSS px from its corner: moved there, and clicked unless MOVE_ONLY."""
    b.eval("document.querySelector('.gm-canvas').scrollIntoView({block: 'center'}); true")
    time.sleep(0.3)
    r = json.loads(b.eval("JSON.stringify(document.querySelector('.gm-canvas').getBoundingClientRect())"))
    px, py = r["left"] + x, r["top"] + y
    b.call("Input.dispatchMouseEvent", {"type": "mouseMoved", "x": px, "y": py}, b.session)
    if not move_only:
        for t in ("mousePressed", "mouseReleased"):
            b.call("Input.dispatchMouseEvent", {"type": t, "x": px, "y": py, "button": "left", "clickCount": 1}, b.session)
    time.sleep(0.6)


def genome_xy(b, doc, chrom, bp):
    """Where BP of CHROM sits on the genome level of the canvas (map.js gx and rowY, at the desktop width)."""
    cw = b.eval("document.querySelector('.gm-canvas').clientWidth")
    rows = [c for c in doc["ref"]["chroms"] if not (c["n"] == "chrY" and doc["people"][0].get("has_y") is False)]
    i = [c["n"] for c in rows].index(chrom)
    ch = b.eval("document.querySelector('.gm-canvas').clientHeight")
    rh = max(22, min(44, (ch - 18) / len(rows)))
    return 40 + bp * (cw - 40 - 22) / rows[0]["len"], 9 + i * rh + rh / 2 - 1


def state_problems(b, base):
    """The genome and chromosome levels, the links that open each, a click on a row, Esc, and a place that
    cannot be read."""
    bad = []
    g = map_page(b, base, S)
    if not (g["level"] == "genome" and "all 24 chromosomes" in g["title"] and g["view"] == "genome"):
        bad.append(f"genome: {g['level']} {g['title']!r} {g['view']}")
    click_canvas(b, 200, 9 + 2 * 22 + 11)
    c = map_now(b)
    if not (c["level"] == "chromosome" and c["chrom"] == "chr3" and c["hash"] == "#chr3"):
        bad.append(f"a click on chromosome 3's row: {c['level']} {c['chrom']} {c['hash']}")
    b.eval("document.querySelector('.gm-canvas').dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true})); true")
    c = map_now(b)
    if c["level"] != "genome" or c["hash"] != "#genome":
        bad.append(f"Esc: {c['level']} {c['hash']}")
    for link, chrom, view in (("#chr16", "chr16", "chr16"), ("#chr16:10.0-40.0M", "chr16", "chr16:10.0-40.0M"),
                              ("#view=chromosome&chr=11", "chr11", "chr11"), ("#chr2:20000001-60000000", "chr2", "chr2:20.0-60.0M")):
        c = map_page(b, base, S, link)
        if not (c["level"] == "chromosome" and c["chrom"] == chrom and c["view"] == view):
            bad.append(f"{link}: {c['level']} {c['chrom']} {c['view']}")
        elif any(m["chrom"] != chrom for m in c["marks"]):
            bad.append(f"{link}: marks of another chromosome listed: {c['marks']}")
    b.eval("location.hash = '#genome'; true")  # after the last link's place
    c = map_now(b)
    if c["level"] != "genome":
        bad.append(f"#genome after a place: {c['level']}")
    b.eval("location.hash = '#chrQ'; true")
    c = map_now(b)
    if "Could not read the place" not in c["toast"] or c["level"] != "genome":
        bad.append(f"#chrQ: toast {c['toast']!r}, level {c['level']}")
    return bad


def map_veil_problems(b, base, out):
    """On the son's closed page: his closed result to act on (the planted deletion) is in no mark, tooltip,
    chip, sentence, link or address of the map until opened; the chip's Open (after its confirm) opens one
    result per click in the page's order (chip_open), so his result with no place opens first, in its section;
    then Open and the card's "Show the result" each bring back exactly the deletion's mark, on the map and in
    the card; Close it again and a reload close it. On DEMO-001, an adult's page, APOE stays closed until
    opened, from the chip or its card's fold."""
    bad = []
    son = OTHERS[1]
    doc = read_map(os.path.join(out, "runs", son, "map.json"))
    f = doc["per"][son]["findings"]
    acts = [i for i, lv in enumerate(f["level"]) if lv == "act" and not f["gate"][i]]
    if not acts:
        return ["the son's map has no result to act on to close"]
    i = acts[0]
    fid, gene = f["id"][i], f["gene"][i]
    words = [gene, f"{f['chrom'][i]}:{f['s'][i] + 1}", f"{f['s'][i] + 1:,}", fid]
    x, y = genome_xy(b, doc, f["chrom"][i], (f["s"][i] + f["e"][i]) / 2)

    def leaks(c, but=""):
        return [w for w in words if w not in but and (w in c["text"] or w in c["tip"] or w in c["hash"])] + \
            [m["id"] for m in c["marks"] if m["id"] == fid]
    c = map_page(b, base, son)
    click_canvas(b, x, y, move_only=True)
    c = map_now(b)
    if c["state"] != "drawn":
        return [f"the son's map: {c['state']}"]
    if leaks(c):
        bad.append(f"closed: the map shows {leaks(c)}")
    # The chip counts every result to act on the page keeps closed, with a place on the map or not.
    want = b.eval("new Set(Array.prototype.filter.call(document.querySelectorAll('.finding[data-level=act]'), "
                  "function (a) { return a.querySelector('.veil'); }).map(function (a) { return a.dataset.finding; })).size")
    label = "1 result closed" if want == 1 else f"{want} results closed"
    if not want or not any(t.startswith(label + " ") for t in c["chips"]):
        bad.append(f"closed: no chip counting the page's {want} closed results: {c['chips']}")
    # Open opens the page's closed results one per click, in the page's order, so the son's own result with no
    # place (his first closed card) opens before the deletion: the map then names it only in its section chip.
    order = chips_now(b)["closed"]
    if not want:
        return bad + ["the son's page shows no closed result to open"]
    k, at = 0, None
    for k in range(want):
        got, at, more = chip_open(b)
        if more:
            return bad + [f"Open {k + 1}: {x}" for x in more]
        if at != order[k]:
            return bad + [f"Open {k + 1} opened {got}, the act card at {at}, not the page's next closed one {order[k]} ({order})"]
        if got == fid:
            break
        section = [o for o in chips_now(b)["opened"] if o["kind"] == "section"]
        if [o["href"] for o in section] != [f"#f/{got}"] or "Opened in its section: " not in section[0]["text"]:
            bad.append(f"Open {k + 1} opened {got}, which has no place: section chips {section}")
        c = map_now(b)
        if leaks(c, " ".join(o["text"] for o in section)):  # the opened result names its own gene
            bad.append(f"Open {k + 1} ({got}): the map shows {leaks(c)}")
    c = map_now(b)
    card = b.eval("!!document.querySelector('#f-%s h3') && !document.querySelector('#f-%s .veil')" % (fid, fid))
    if [m["id"] for m in c["marks"] if m["id"] == fid] != [fid] or not card:
        bad.append(f"Open: mark {[m['id'] for m in c['marks']]}, card shown {card}")
    if not any(t.startswith("Opened: " + f["title"][i]) for t in c["chips"]):
        bad.append(f"Open: chips {c['chips']}")
    if fid in c["hash"] or gene in c["hash"]:
        bad.append(f"Open: the address names it: {c['hash']}")
    click_canvas(b, x, y, move_only=True)
    c = map_now(b)
    if gene not in c["tip"]:
        bad.append(f"Open: its tooltip at its place does not name it: {c['tip']!r}")
    # The results opened before it in their section close again from their chips, and those chips go.
    b.eval("Array.prototype.forEach.call(document.querySelectorAll('.gm-chip.opened[data-opened=section]'), "
           "function () { document.querySelector('.gm-chip.opened[data-opened=section] button').click(); }); true")
    if [o for o in chips_now(b)["opened"] if o["kind"] == "section"] or chips_now(b)["closed"] != order[:k]:
        bad.append(f"Close it again on the section chips: chips {chips_now(b)['opened']}, closed {chips_now(b)['closed']}")
    b.eval("document.querySelector('.gm-chip.opened[data-opened=\\'1\\'] button').click(); true")
    c = map_now(b)
    card = at in chips_now(b)["closed"]
    if leaks(c) or not card:
        bad.append(f"Close it again: map shows {leaks(c)}, card closed {card}")
    b.eval("document.querySelectorAll('.finding[data-level=act]')[%d].querySelector('.veil button').click(); true" % at)
    c = map_now(b)
    if [m["id"] for m in c["marks"] if m["id"] == fid] != [fid]:
        bad.append("the card's Show the result does not open it on the map")
    c = map_page(b, base, son)
    if leaks(c):
        bad.append(f"a reload does not close it: {leaks(c)}")
    # APOE on an adult's page.
    c = map_page(b, base, S)
    if any(m["id"] == "apoe" for m in c["marks"]) or not any(t.startswith("APOE closed") for t in c["chips"]):
        bad.append(f"APOE on DEMO-001: marks {[m['id'] for m in c['marks']]}, chips {c['chips']}")
    b.eval("document.querySelector('.gm-chip[data-closed=gate] button').click(); true")
    time.sleep(0.3)
    b.eval("document.querySelector('.gm-chip[data-closed=gate] .gm-yes').click(); true")
    c = map_now(b)
    fold = b.eval("(document.querySelector('details[data-gate=apoe]') || {}).open")
    if not any(m["id"] == "apoe" for m in c["marks"]) or fold is not True:
        bad.append(f"Open APOE: marks {[m['id'] for m in c['marks']]}, its card's fold open {fold}")
    b.eval("document.querySelector('details[data-gate=apoe] summary').click(); true")
    c = map_now(b)
    if any(m["id"] == "apoe" for m in c["marks"]):
        bad.append("closing APOE's card fold leaves it on the map")
    b.eval("document.querySelector('details[data-gate=apoe] summary').click(); true")
    c = map_now(b)
    if not any(m["id"] == "apoe" for m in c["marks"]):
        bad.append("opening APOE's card fold does not open it on the map")
    c = map_page(b, base, S)
    if any(m["id"] == "apoe" for m in c["marks"]):
        bad.append("a reload does not close APOE")
    return bad


def overlap_problems(b, base, out):
    """On the son's closed page, a second result to act on over the first carries only the first one's tag
    (rust/crates/og-map/README.md). Opening the first leaves the second closed: in no mark, sentence, link,
    title or address of the map."""
    son = OTHERS[1]
    path = os.path.join(out, "runs", son, "map.json")
    with open(path, "rb") as fh:
        keep = fh.read()
    try:
        doc = read_map(path)
        f = doc["per"][son]["findings"]
        acts = [i for i, lv in enumerate(f["level"]) if lv == "act" and not f["gate"][i]]
        if not acts:
            return ["the son's map has no result to act on to close"]
        i, n = acts[0], len(f["id"])
        fid = f["id"][i]
        for k in f:
            if isinstance(f[k], list) and len(f[k]) == n:
                f[k].append(f[k][i])
        plant = {"id": "planted-overlap", "gene": "PLANTEDGENE", "title": "Planted overlapping result",
                 "href": "finding/planted-overlap/", "phrase": "Planted overlap phrase."}
        for k, v in plant.items():
            f[k][-1] = v
        write_map(path, doc)
        c = map_page(b, base, son, "#" + f["chrom"][i])
        if c["state"] != "drawn":
            return [f"the planted map: {c['state']}"]
        if not open_closed(b, fid):
            return [f"a link to {fid} lands on no closed card"]
        c = map_now(b)
        bad = []
        if fid not in [m["id"] for m in c["marks"]]:
            bad.append(f"the first result did not open: marks {[m['id'] for m in c['marks']]}")
        leaked = [w for w in plant.values() if w in c["text"] or w in c["tip"] or w in c["hash"]]
        leaked += [m["id"] for m in c["marks"] if m["id"] == "planted-overlap"]
        if leaked:
            bad.append(f"opening the first shows the second: {leaked}")
        return bad
    finally:
        with open(path, "wb") as fh:
            fh.write(keep)


# The figure as a reader gets it before any click: its HTML (the chips aside, which count the closed results
# on purpose, and the timing attributes), its marks, margin and layers, and a hash of the canvas's pixels.
FIGURE_JS = r"""JSON.stringify((function () {
  var f = document.getElementById('genome-map');
  if (!f) { return null; }
  var c = f.cloneNode(true);
  Array.prototype.forEach.call(c.querySelectorAll('.gm-chips'), function (x) { x.textContent = ''; });
  c.removeAttribute('data-pgp-map-parse-ms');
  c.removeAttribute('data-pgp-map-paint-ms');
  return { html: c.outerHTML, canvas: f.querySelector('.gm-canvas').toDataURL() };
}()))"""


def plant_son(out, extra, og):
    """The son's results.json with the findings EXTRA added, its findings' veil tags (results.veil.findings)
    computed again by live_results.veil as the server computes them, and his map.json built again from it by og
    map build, as the map plugin builds it: og tags each row itself. Returns what to write back afterwards."""
    son = OTHERS[1]
    rj, mj = (os.path.join(out, "runs", son, n) for n in ("results.json", "map.json"))
    keep = {}
    for p in (rj, mj):
        with open(p, "rb") as fh:
            keep[p] = fh.read()
    res = json.loads(keep[rj])
    res["findings"] = res["findings"] + [dict(f) for f in extra]
    sec = {"act": 1, "medicine": 2, "carrier": 3, "review": 4}
    for f in extra:  # each in its level's section, as live_results lists them: its card shows on the page
        for j in res["journal"]:
            if j["n"] == sec[f["level"]]:
                j["findings"] = j["findings"] + [f["id"]]
    tags = live_results.veil({}, res["findings"], []).get("findings") or {}
    res["veil"] = dict(res.get("veil") or {}, findings=dict((res.get("veil") or {}).get("findings") or {}, **tags))
    with open(rj, "w") as fh:
        json.dump(res, fh)
    root = os.path.join(os.path.dirname(out), "root")
    doc = build_live_demo.og_map(og, ["build", "--root", root, "--run", son, "--assets", os.path.join(REPO, "assets"),
                                      "--cache-dir", os.path.join(root, ".owngenome-live", "map")],
                                 json.dumps({"results": res}).encode())
    with open(mj, "wb") as fh:
        fh.write(doc)
    return keep


def unplant(keep):
    for p, data in keep.items():
        with open(p, "wb") as fh:
            fh.write(data)


def row_tag(out, fid):
    """The id whose tag og put on FID's row of the son's map, or None."""
    son = OTHERS[1]
    doc = read_map(os.path.join(out, "runs", son, "map.json"))
    f = doc["per"][son]["findings"]
    if fid not in f["id"]:
        return None
    v = f["veil"][f["id"].index(fid)]
    return doc["dict"]["veil"][v] if v >= 0 else None


def son_figure(b, base, hash_=""):
    """The son's map before any click: (the map's facts, the figure with its canvas hashed)."""
    c = map_page(b, base, OTHERS[1], hash_)
    fig = json.loads(b.eval(FIGURE_JS) or "null") or {}
    fig["canvas"] = hashlib.sha256((fig.get("canvas") or "").encode()).hexdigest()[:12]
    return c, fig


def act(fid, gene, viewer, title):
    return {"id": fid, "level": "act", "card": "clinical-filter", "gene": gene, "viewer": viewer, "title": title,
            "detail": f"{gene} {viewer}", "text": f"{title}, planted."}


# An open finding over a closed result whose card names another gene: og tags the open finding's row (the
# whole gene span) with the closed result, while the page shows its card. A map that hides that row only
# while the result is closed places the closed result inside that gene. Each planting: (what, the closed
# results, the open findings in both pages, the open finding to watch, chromosome).
UNDER = (("a gene-only review finding (CFTR) over a closed deletion naming another gene",
          [act("pl-iota", "PLIOTA", "chr7:117400000-117700000", "Planted iota deletion")],
          [{"id": "pl-kappa", "level": "review", "card": "exomiser", "gene": "CFTR", "title": "Planted kappa review",
            "detail": "CFTR", "text": "A planted review finding that names CFTR only."}], "pl-kappa", "chr7"),
         # The page's own CYP2D6 medicine finding closes with the medicine card while any result is closed (a fixed
         # panel), so a carrier finding stands in for a gene-only finding of a card that is no panel.
         ("a gene-only carrier finding (CYP2D6) over a closed variant naming another gene",
          [act("pl-theta", "PLTHETA", "chr22:42125000", "Planted theta result")],
          [{"id": "pl-lambda", "level": "carrier", "card": "carrier-list", "gene": "CYP2D6",
            "title": "Planted lambda carrier", "detail": "CYP2D6", "text": "A planted carrier finding that names CYP2D6 only."}],
          "pl-lambda", "chr22"))


def under_problems(b, base, out, og):
    """For each planting of UNDER: the son's map before any click is the same, figure HTML, canvas, marks,
    margin and layers, with the closed result and without it, at the genome level and on its chromosome; and
    the open finding is drawn in both."""
    bad = []
    for what, closed, opened, watch, chrom in UNDER:
        got = {}
        for name, extra in (("without", opened), ("with", opened + closed)):
            keep = plant_son(out, extra, og)
            try:
                if name == "with" and row_tag(out, watch) != closed[0]["id"]:
                    bad.append(f"{what}: og did not tag {watch}'s row with {closed[0]['id']} ({row_tag(out, watch)}): "
                               "the planting tests nothing")
                got[name] = [son_figure(b, base, h) for h in ("", "#" + chrom)]
                with open(os.path.join(out, "runs", OTHERS[1], "results.json")) as fh:
                    title = [f["title"] for f in json.load(fh)["findings"] if f["id"] == watch][0]
                if not b.eval("document.body.innerText.indexOf(%s) >= 0" % json.dumps(title)):
                    bad.append(f"{what}: the page does not show {watch}'s card ({title}): the planting tests nothing")
            finally:
                unplant(keep)
        for (c0, f0), (c1, f1), lv in zip(got["without"], got["with"], ("genome", chrom)):
            if watch not in [m["id"] for m in c0["marks"]]:
                bad.append(f"{what}, {lv}: {watch} is not on the map even with no closed result")
            if f0 != f1:
                diff = [k for k in ("html", "canvas") if f0.get(k) != f1.get(k)]
                bad.append(f"{what}, {lv}: the figure differs with the closed result ({diff}): marks "
                           f"{[m['id'] for m in c0['marks']]} vs {[m['id'] for m in c1['marks']]}")
    return bad


def second_open_problems(b, base, out, og):
    """Two overlapping closed results on the son's page (og tags the second with the first): opening the
    second by its own card draws exactly it on the map, never the first."""
    two = [act("pl-alpha", "PLALPHA", "chr17:43050000-43060000", "Planted alpha result"),
           act("pl-beta", "PLBETA", "chr17:43055000", "Planted beta result")]
    keep = plant_son(out, two, og)
    try:
        if row_tag(out, "pl-beta") != "pl-alpha":
            return [f"og did not tag pl-beta's row with pl-alpha ({row_tag(out, 'pl-beta')}): the planting tests nothing"]
        c = map_page(b, base, OTHERS[1], "#chr17")
        ids = [m["id"] for m in c["marks"]]
        if "pl-alpha" in ids or "pl-beta" in ids:
            return [f"closed: marks {ids}"]
        if not open_closed(b, "pl-beta"):
            return ["a link to pl-beta lands on no closed card"]
        c = map_now(b)
        ids = [m["id"] for m in c["marks"]]
        bad = []
        if "pl-beta" not in ids:
            bad.append(f"opened by its card, pl-beta is not on the map: marks {ids}, chips {c['chips']}")
        if "pl-alpha" in ids or "PLALPHA" in c["text"]:
            bad.append(f"opening pl-beta shows pl-alpha: marks {ids}")
        return bad
    finally:
        unplant(keep)


# A finding the server tags with a closed result because it names its gene in its genes (results.veil.findings:
# an open structural variant over two genes), once with og's tag on its row too (CFTR, inside the deletion's span)
# and once with none (BRCA1).
NAMED = (act("pl-x", "PLXGENE", "chr7:117400000-117700000", "Planted x deletion"),
         {"id": "pl-q", "level": "review", "card": "exomiser", "gene": "CFTR", "genes": ["CFTR", "PLXGENE"],
          "title": "Symptom ranking: PLXGENE ranks 1", "detail": "CFTR",
          "text": "A planted ranking that names the deletion's gene."},
         {"id": "pl-r", "level": "review", "card": "exomiser", "gene": "BRCA1", "genes": ["BRCA1", "PLXGENE"],
          "title": "Symptom ranking: PLXGENE ranks 2", "detail": "BRCA1",
          "text": "A planted ranking that names the deletion's gene."})
# A row of the map with no finding in results.json (a newer finding the map holds before the page has it),
# on the closed deletion's place and with its tag.
GHOST = {"id": "pl-ghost", "level": "review", "gene": "GHOSTGENE", "title": "Ghost review", "href": "finding/pl-ghost/",
         "phrase": "Ghost review phrase."}


def named_problems(b, base, out, og):
    """On the son's closed page, with the planted deletion pl-x closed: the findings the server closes with it
    (pl-q, pl-r) and a tagged map row the page has no finding for (pl-ghost) are in no mark, sentence, link,
    title or address, and the figure is the one of the page without them; opening pl-x by its card draws pl-x,
    pl-q and pl-r."""
    son = OTHERS[1]
    _, before = son_figure(b, base)
    keep = plant_son(out, NAMED, og)
    try:
        with open(os.path.join(out, "runs", son, "results.json")) as fh:
            tags = (json.load(fh).get("veil") or {}).get("findings") or {}
        if tags.get("pl-q") != ["pl-x"] or tags.get("pl-r") != ["pl-x"] or row_tag(out, "pl-r"):
            return [f"the server's tags {tags}, og's tag on pl-r {row_tag(out, 'pl-r')}: the planting tests nothing"]
        path = os.path.join(out, "runs", son, "map.json")
        doc = read_map(path)
        f = doc["per"][son]["findings"]
        if "pl-x" not in f["id"] or "pl-x" not in doc["dict"]["veil"]:
            return ["og did not place or tag pl-x: the planting tests nothing"]
        i, n = f["id"].index("pl-x"), len(f["id"])
        for k in f:
            if isinstance(f[k], list) and len(f[k]) == n:
                f[k].append(f[k][i])
        for k, v in GHOST.items():
            f[k][-1] = v
        f["veil"][-1] = doc["dict"]["veil"].index("pl-x")
        f["veils"][-1] = [f["veil"][-1]]
        write_map(path, doc)
        c, fig = son_figure(b, base)
        bad = []
        ids = [m["id"] for m in c["marks"]]
        shown = [x for x in ("pl-x", "pl-q", "pl-r", "pl-ghost") if x in ids]
        shown += [w for w in ("PLXGENE", "Ghost", "GHOSTGENE") if w in c["text"] or w in c["tip"] or w in c["hash"]]
        if shown:
            bad.append(f"closed: the map shows {shown}")
        if fig != before:
            bad.append(f"closed: the figure differs from the page without them "
                       f"({[k for k in ('html', 'canvas') if fig.get(k) != before.get(k)]})")
        if not open_closed(b, "pl-x"):
            return bad + ["a link to pl-x lands on no closed card"]
        c = map_now(b)
        ids = [m["id"] for m in c["marks"]]
        if not {"pl-x", "pl-q", "pl-r"} <= set(ids):
            bad.append(f"opened by its card: marks {ids}")
        return bad
    finally:
        unplant(keep)


def chip_problems(b, base, out, og):
    """The son's map before any click is the same, its chips included, whether his planted closed result has a
    place on the map (chr13) or none (no viewer, a gene not in GENCODE): the chip counts the results the page
    keeps closed, never only the placed ones."""
    got, bad = {}, []
    for name, fid, viewer in (("no place", "pl-nowhere", ""), ("a place", "pl-somewhere", "chr13:32340000")):
        keep = plant_son(out, [act(fid, "PL" + fid[3:].upper(), viewer, "Planted closed result")], og)
        try:
            placed = fid in read_map(os.path.join(out, "runs", OTHERS[1], "map.json"))["per"][OTHERS[1]]["findings"]["id"]
            if placed != bool(viewer):
                bad.append(f"{name}: og placed {fid} {placed}: the planting tests nothing")
            c, fig = son_figure(b, base)
            got[name] = (c["chips"], fig)
        finally:
            unplant(keep)
    if got["no place"] != got["a place"]:
        bad.append(f"the chips {got['no place'][0]} with no place, {got['a place'][0]} with one"
                   + ("" if got["no place"][1] == got["a place"][1] else ", and the figure differs"))
    return bad


# Two planted closed results, one with a place on the map and one without, in both orders of the page.
OPEN_ORDER = (("unplaced pl-a, placed pl-b", [act("pl-a", "PLAGENE", "", "Planted a result"),
                                            act("pl-b", "PLBGENE", "chr13:32340000", "Planted b result")]),
              ("placed pl-a, unplaced pl-b", [act("pl-a", "PLAGENE", "chr13:32340000", "Planted a result"),
                                            act("pl-b", "PLBGENE", "", "Planted b result")]))


def open_order_problems(b, base, out, og):
    """For each planting of OPEN_ORDER, the chip's Open clicked until nothing is closed (chip_open checks each
    click): the ids it opens, click by click, are the same in both, so the order never tells which closed result
    has a place. A result with no place opens with an "Opened in its section" chip that links its card; its
    Close it again closes it on the page and takes the chip away. Before that, Open then Keep it closed opens
    nothing and gives the keyboard back to the chips."""
    seqs, bad = {}, []
    for k, (what, extra) in enumerate(OPEN_ORDER):
        keep = plant_son(out, extra, og)
        try:
            placed = read_map(os.path.join(out, "runs", OTHERS[1], "map.json"))["per"][OTHERS[1]]["findings"]["id"]
            if [x["id"] in placed for x in extra] != [bool(x["viewer"]) for x in extra]:
                return [f"{what}: og placed {[x['id'] for x in extra if x['id'] in placed]}: the planting tests nothing"]
            map_page(b, base, OTHERS[1])
            first = chips_now(b)
            if k == 0:
                b.eval("var o = document.querySelector('.gm-chip[data-closed=act] button'); o.focus(); o.click(); true")
                b.eval("document.querySelector('.gm-chip[data-closed=act] .gm-no').click(); true")
                kept = chips_now(b)
                if kept["closed"] != first["closed"] or kept["confirm"] or kept["focus"] != "first":
                    return [f"{what}: Keep it closed: closed {kept['closed']}, confirm {kept['confirm']}, "
                            f"keyboard on {kept['focus']!r}"]
            seq = []
            for _ in first["closed"]:
                got, _, more = chip_open(b)
                if more:
                    return bad + [f"{what}, Open {len(seq) + 1}: {x}" for x in more]
                seq.append(got)
                opened = chips_now(b)["opened"]
                section = [o["href"] for o in opened if o["kind"] == "section"]
                if (got in placed) == (f"#f/{got}" in section):
                    return [f"{what}: {got} (placed {got in placed}) opened with section chips {section}"]
            seqs[what] = seq
            if k == 0:
                gone = [x["id"] for x in extra if not x["viewer"]][0]
                at = chips_now(b)["ids"].index(gone) if gone in chips_now(b)["ids"] else -1
                was = [o for o in chips_now(b)["opened"] if o["kind"] == "section"]
                raw = [o["raw"] for o in was if o["href"] == f"#f/{gone}"]
                b.eval("var o = document.querySelector('.gm-chip.opened[data-opened=section] a[href=\"%s\"]'); "
                       "o.parentNode.parentNode.querySelector('button').click(); true" % (raw[0] if raw else ""))
                after = chips_now(b)
                left = [o for o in after["opened"] if o["kind"] == "section"]
                if after["closed"] != [at] or len(left) != len(was) - 1:
                    bad.append(f"{what}: Close it again on {gone}'s section chip: closed {after['closed']}, "
                               f"opened chips {after['opened']}")
        finally:
            unplant(keep)
    if len(set(map(tuple, seqs.values()))) != 1:
        bad.append(f"the order Open opens them tells which has a place: {seqs}")
    return bad


def runs_sentence(b, base):
    """Whether DEMO-001's planted runs sentence on chr1 is in the margin at the first frame of the replay and
    at the last: (first, last, the level and chromosome at the first frame)."""
    map_page(b, base, S, "#chr1")
    seen = []
    for to in ("0", "s.max"):
        b.eval("var s = document.getElementById('replay-slider'); s.value = %s; s.dispatchEvent(new Event('input')); true" % to)
        c = map_now(b)
        seen.append(("PLANTEDRUNS" in c["text"], c["level"] + " " + c["chrom"]))
    return seen[0][0], seen[1][0], seen[0][1]


def runs_sentence_problems(b, base, out):
    """A replay before the runs-of-homozygosity step ended has no sentence about them on an open chromosome;
    at its end it has it."""
    path = os.path.join(out, "runs", S, "map.json")
    with open(path, "rb") as fh:
        keep = fh.read()
    try:
        doc = read_map(path)
        doc.setdefault("text", {}).setdefault("chrom", {})["chr1"] = [{"t": "PLANTEDRUNS sentence.", "b": False}]
        write_map(path, doc)
        first, last, where = runs_sentence(b, base)
        if where != "chromosome chr1":
            return [f"the first frame is at {where}: the check tests nothing"]
        return (["the first frame says the runs sentence"] if first else []) + ([] if last else ["the last frame lacks it"])
    finally:
        with open(path, "wb") as fh:
            fh.write(keep)


# A second figure on DEMO-001's page with a veil that keeps nothing closed: drawn from map.json, then reloaded
# from a file that is not there.
SECOND_JS = r"""(function () {
  var o = { url: 'map.json', veil: function () { return 'open'; }, now: function () { return true; },
    known: function () { return true; }, stepDone: function () { return true; }, open: function () {},
    close: function () {}, fallback: function () { return []; } };
  var g = window.PGPMap.create(o);
  g.node.id = 'gm-second';
  document.body.appendChild(g.node);
  window.pgpSecond = { g: g, o: o };
  g.refresh('k1');
  return true;
}())"""
SECOND_FACTS = r"""JSON.stringify((function () {
  var f = document.getElementById('gm-second');
  return { state: f.dataset.pgpMap, marks: f.querySelectorAll('.gm-marks a').length,
    narr: f.querySelector('.gm-narration').textContent, layers: f.querySelector('.gm-layers').textContent,
    tip: f.querySelector('.gm-tip').hidden ? '' : f.querySelector('.gm-tip').textContent };
}()))"""


def failed_reload_problems(b, base):
    """A figure whose reload fails draws nothing more of the document it had: no mark link, no margin sentence,
    no layer line, no tooltip."""
    map_page(b, base, S)
    b.eval(SECOND_JS)
    if not b.wait("document.getElementById('gm-second').dataset.pgpMap === 'drawn'"):
        return ["the second figure never drew"]
    time.sleep(0.6)
    drawn = json.loads(b.eval(SECOND_FACTS))
    if not drawn["marks"] or not drawn["narr"]:
        return [f"the second figure drew no mark or sentence: the check tests nothing ({drawn})"]
    b.eval("pgpSecond.o.url = 'no-such-map.json'; pgpSecond.g.refresh('k2'); true")
    if not b.wait("document.getElementById('gm-second').dataset.pgpMap === 'failed'"):
        return ["the second figure's reload did not fail"]
    time.sleep(1.0)
    c = json.loads(b.eval(SECOND_FACTS))
    return [f"after the failed reload: {c}"] if c["marks"] or c["narr"] or c["layers"] or c["tip"] else []


def replay_apoe_problems(b, base):
    """DEMO-001 at the first frame of its replay, before the traits step: the map has no APOE chip, no
    switch to open it and no APOE mark, as its card has no fold yet."""
    map_page(b, base, S)
    b.eval("var s = document.getElementById('replay-slider'); s.value = 0; s.dispatchEvent(new Event('input')); true")
    c = map_now(b)
    bad = []
    if b.eval("!!document.querySelector('details[data-gate=apoe]')"):
        bad.append("the traits card has its APOE fold at the first frame: the check tests nothing")
    if any("APOE" in t for t in c["chips"]) or any(m["id"] == "apoe" for m in c["marks"]):
        bad.append(f"before the traits step: chips {c['chips']}, marks {[m['id'] for m in c['marks']]}")
    return bad


def missing_layer_problems(b, base, out):
    """A map whose ROH layer is missing says so in words, in the margin, never a blank."""
    path = os.path.join(out, "runs", S, "map.json")
    with open(path, "rb") as fh:
        keep = fh.read()
    try:
        doc = read_map(path)
        doc["layers"]["roh"] = dict(doc["layers"]["roh"], state="missing", count=0, note="planted: no ROH text")
        doc["per"][S].pop("roh", None)
        write_map(path, doc)
        c = map_page(b, base, S, "#chr1")
        row = [x for x in c["layers"] if x[0] == "roh"]
        if not row or row[0][1] != "missing" or "Not on the map: planted: no ROH text" not in row[0][2]:
            return [f"the missing ROH layer reads {row}"]
        return []
    finally:
        with open(path, "wb") as fh:
            fh.write(keep)


def fallback_problems(b, base, out):
    """A browser with no DecompressionStream: the sentence and the findings with a place, each through the veil."""
    son = OTHERS[1]
    doc = read_map(os.path.join(out, "runs", son, "map.json"))
    f = doc["per"][son]["findings"]
    gene = [f["gene"][i] for i, lv in enumerate(f["level"]) if lv == "act"][0]
    b.call("Page.enable", session=b.session)  # a script for each new document needs the Page domain on
    sid = b.call("Page.addScriptToEvaluateOnNewDocument", {"source": "delete window.DecompressionStream;"}, b.session)["identifier"]
    try:
        c = map_page(b, base, son)
        rows = b.eval("document.querySelectorAll('.gm-fallback table tbody tr').length")
    finally:
        b.call("Page.removeScriptToEvaluateOnNewDocument", {"identifier": sid}, b.session)
    bad = []
    if c["state"] != "nodecompress" or "This browser cannot open the map file" not in c["text"]:
        bad.append(f"no DecompressionStream: {c['state']}")
    if not rows:
        bad.append("no DecompressionStream: no findings table")
    if gene in c["text"]:
        bad.append(f"no DecompressionStream: the table names the closed {gene}")
    return bad


def map_screens(b, base, dest):
    """Screenshots of the map, dark, at 1440 x 900 and 390 x 844: genome and chromosome, a closed page."""
    os.makedirs(dest, exist_ok=True)
    taken = []
    try:
        for w, h in ((1440, 900), (390, 844)):
            b.call("Emulation.setDeviceMetricsOverride", {"width": w, "height": h, "deviceScaleFactor": 1, "mobile": False},
                   b.session)
            for run, link, name in ((S, "", "genome"), (S, "#chr16", "chromosome"), (OTHERS[1], "", "closed")):
                got = map_page(b, base, run, link)
                want = "chromosome" if link else "genome"
                if got["state"] != "drawn" or got["level"] != want:
                    print(f"ERROR: {run}{link}: the map is {got['state']} at {got['level'] or 'no'} level, "
                          f"want drawn at {want}: no screenshot")
                    continue
                b.eval("document.getElementById('genome-map').scrollIntoView({block: 'start'}); true")
                time.sleep(0.5)
                png = b.call("Page.captureScreenshot", {"format": "png"}, b.session)["data"]
                p = os.path.join(dest, f"map-{name}-{w}x{h}.png")
                with open(p, "wb") as fh:
                    fh.write(base64.b64decode(png))
                taken.append(p)
    finally:
        b.call("Emulation.clearDeviceMetricsOverride", {}, b.session)
    return taken


def map_checks(b, base, out):
    """Figure 1, the genome map: drawn from map.json with every mark the karyogram drew, in time; its levels,
    links and a place it cannot read; the veil and opening on purpose; a missing layer; no DecompressionStream;
    no sideways scroll at 390 px."""
    if not os.path.isfile(os.path.join(out, "runs", S, "map.json")):
        ci = os.environ.get("CI") == "true"
        print(("[FAIL] " if ci else "[SKIP] ") + "the demo has no map.json: no og to build it")
        if ci:
            FAILS.append("no map.json in CI")
        return
    b.call("Network.enable", session=b.session)  # the controls below change files: never serve them from cache
    b.call("Network.setCacheDisabled", {"cacheDisabled": True}, b.session)
    # A desktop window, so the whole canvas fits in view for the clicks.
    b.call("Emulation.setDeviceMetricsOverride", {"width": 1440, "height": 900, "deviceScaleFactor": 1, "mobile": False},
           b.session)
    rows = karyogram_marks()
    bad = parity_problems(b, base, rows)
    check(f"map: the DEMO-001 to DEMO-003 person pages draw the map from map.json and every karyogram mark "
          f"({len(rows)} marks: same level and chromosome, within 0.01 Mb)", not bad and len(rows) >= 9, bad)
    path = os.path.join(out, "runs", S, "map.json")
    with open(path, "rb") as fh:
        keep = fh.read()
    try:
        doc = read_map(path)
        f = doc["per"][S]["findings"]
        drop = f["id"].index(next(r[1] for r in rows if r[0] == S))
        for k in f:
            f[k] = [v for j, v in enumerate(f[k]) if j != drop]
        write_map(path, doc)
        control("a map.json without one of the karyogram's findings", not parity_problems(b, base, rows))
    finally:
        with open(path, "wb") as fh:
            fh.write(keep)
    g = map_page(b, base, S)
    check(f"map: first paint within 150 ms of the document arriving, parsed within 30 ms ({g.get('paint')} ms, "
          f"{g.get('parse')} ms)", g.get("paint", 999) <= 150 and g.get("parse", 999) <= 30, g)
    bad = state_problems(b, base)
    check("map: genome and chromosome levels; #chr16, #chr16:10.0-40.0M, chr=11 and a 1-based span open their view; "
          "a click on a row opens it, Esc closes it; an unknown place gives a toast", not bad, bad)
    js = os.path.join(out, "assets", "person.js")
    with open(js) as fh:
        keep_js = fh.read()
    bad = map_veil_problems(b, base, out)
    check("map: the son's closed result is in no mark, tooltip, chip, sentence, link or address until opened; Open "
          "(after its confirm) opens one result per click, in the page's order, and the count falls by one; Open and "
          "the card's Show the result each bring back exactly its mark there and in the card; "
          "Close it again and a reload close it; APOE stays closed on an adult's page until opened, from the chip "
          "or its card's fold", not bad, bad)
    mjs = os.path.join(out, "assets", "map.js")
    with open(mjs) as fh:
        keep_mjs = fh.read()
    one = "opts.open(id);  // "
    every = "f.level === 'act' && veiled(f); })"
    try:
        with open(js, "w") as fh:
            fh.write(keep_js.replace("function mapVeil(id, gated) {", "function mapVeil(id, gated) { return 'open';"))
        control("person.js whose map veil keeps nothing closed", not map_veil_problems(b, base, out))
        with open(js, "w") as fh:
            fh.write(keep_js.replace(every, "f.level === 'act'; })"))
        control("person.js whose map chip counts every result to act on, closed or open", keep_js.count(every) != 1 or
                not map_veil_problems(b, base, out))
        with open(js, "w") as fh:
            fh.write(keep_js)
        with open(mjs, "w") as fh:
            fh.write(keep_mjs.replace(one, "(gated ? [id] : closedActs()).forEach(function (x) { opts.open(x); });  // "))
        control("map.js whose Open opens every closed result at once", keep_mjs.count(one) != 1 or
                not map_veil_problems(b, base, out))
    finally:
        with open(js, "w") as fh:
            fh.write(keep_js)
        with open(mjs, "w") as fh:
            fh.write(keep_mjs)
    bad = overlap_problems(b, base, out)
    check("map: a second closed result to act on over the first stays closed when the first is opened", not bad, bad)
    own = "var own = M.closable.indexOf(r.id) < 0 || stateOf(r.id) !== 'closed';"
    try:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs.replace(own, "var own = true;"))
        control("map.js whose rows ask only their tag, never their own result", keep_mjs.count(own) != 1 or
                not overlap_problems(b, base, out))
    finally:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs)
    og = og_binary() or build_live_demo.OG
    bad = under_problems(b, base, out, og)
    check("map: an open finding over a closed result naming another gene stays on the map; before any click the "
          "figure is the same with the closed result and without it", not bad, bad)
    bad = second_open_problems(b, base, out, og)
    check("map: the second of two overlapping closed results, opened by its own card, is drawn, and the first "
          "is not", not bad, bad)
    known = "(opts.known && opts.known(r.id)) || "
    try:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs.replace(known, ""))
        control("map.js whose page findings follow another result's tag", keep_mjs.count(known) != 1 or
                not under_problems(b, base, out, og))
        control("map.js whose page findings follow another result's tag (the second result opened)",
                keep_mjs.count(known) != 1 or not second_open_problems(b, base, out, og))
    finally:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs)
    bad = named_problems(b, base, out, og)
    check("map: a finding the server closes with a closed result, and a tagged map row the page has no finding "
          "for, stay off the map until the result is opened", not bad, bad)
    now = "&& (!opts.now || opts.now(r.id));"
    try:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs.replace(now, "&& true;"))
        control("map.js whose rows ignore what the page shows now", keep_mjs.count(now) != 1 or
                not named_problems(b, base, out, og))
        with open(mjs, "w") as fh:
            fh.write(keep_mjs.replace(known, "true || "))
        control("map.js that takes every row for a finding of the page", keep_mjs.count(known) != 1 or
                not named_problems(b, base, out, og))
    finally:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs)
    bad = chip_problems(b, base, out, og)
    check("map: the closed-results chip counts the results the page keeps closed, with a place on the map or "
          "without", not bad, bad)
    try:
        with open(js, "w") as fh:
            fh.write(keep_js.replace("closed: mapClosed, ", ""))
        control("person.js that does not tell the map which results it keeps closed",
                keep_js.count("closed: mapClosed, ") != 1 or not chip_problems(b, base, out, og))
    finally:
        with open(js, "w") as fh:
            fh.write(keep_js)
    bad = open_order_problems(b, base, out, og)
    check("map: the chip's Open opens one closed result per click in the page's order, whether it has a place or "
          "not; one with no place gets an Opened in its section chip until closed again; the keyboard follows",
          not bad, bad)
    concat = "if (page) { return page.filter(function (id) { return !M.gated[id]; }); }"
    placed_first = ("if (page) { var placed = M.closable.filter(function (id) { return !M.gated[id] && stateOf(id) === "
                    "'closed' && page.indexOf(id) >= 0; }); return placed.concat(page.filter(function (id) { return "
                    "placed.indexOf(id) < 0; })); }")
    alive = "stateOf(id) !== 'closed' && opts.title && opts.title(id);"
    pushed = "inSection.push(id);"
    focus = "(b || canvas).focus({ preventScroll: true });"
    try:
        for desc, a, z in (("map.js whose Open opens the closed results with a place first", concat, placed_first),
                           ("map.js whose section chip stays for a closed result", alive, "true;"),
                           ("map.js with no section chip for a result with no place", pushed, ""),
                           ("map.js that leaves the keyboard where the chips were rebuilt", focus, "void b;")):
            with open(mjs, "w") as fh:
                fh.write(keep_mjs.replace(a, z))
            control(desc, keep_mjs.count(a) != 1 or not open_order_problems(b, base, out, og))
    finally:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs)
    bad = replay_apoe_problems(b, base)
    check("map: before the traits step a replay has no APOE chip, switch or mark", not bad, bad)
    gated = "if (M.gated[id] && opts.now && !opts.now(id)) { return; }"
    try:
        with open(js, "w") as fh:
            fh.write(keep_js.replace("if (!f) { var t = cards()", "if (!f) { return true; var t = cards()"))
        with open(mjs, "w") as fh:
            fh.write(keep_mjs.replace(gated, ""))
        control("person.js and map.js that show APOE's chip and row before its card", keep_mjs.count(gated) != 1 or
                not replay_apoe_problems(b, base))
    finally:
        with open(js, "w") as fh:
            fh.write(keep_js)
        with open(mjs, "w") as fh:
            fh.write(keep_mjs)
    bad = missing_layer_problems(b, base, out)
    check("map: a missing layer says why in words, never a blank", not bad, bad)
    bad = runs_sentence_problems(b, base, out)
    check("map: a replay says an open chromosome's runs sentence only once the runs' step ended", not bad, bad)
    waits = "&& rohLayer() === M.layers.roh"
    try:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs.replace(waits, ""))
        control("map.js whose runs sentence does not wait for its step", keep_mjs.count(waits) != 1 or
                not runs_sentence_problems(b, base, out))
    finally:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs)
    bad = failed_reload_problems(b, base)
    check("map: after a failed reload the figure draws nothing more of its old document", not bad, bad)
    undone = [("if (!M || failure) { return; }", "if (!M) { return; }", 2), ("if (!raf && !failure) {", "if (!raf) {", 1),
              ("window.cancelAnimationFrame(raf);\n      raf = 0;\n      window.clearTimeout(settleT);\n", "", 1)]
    try:
        planted_mjs = keep_mjs
        for a, z, _ in undone:
            planted_mjs = planted_mjs.replace(a, z)
        with open(mjs, "w") as fh:
            fh.write(planted_mjs)
        control("map.js that keeps drawing after a failed reload",
                any(keep_mjs.count(a) != n for a, _, n in undone) or not failed_reload_problems(b, base))
    finally:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs)
    bad = fallback_problems(b, base, out)
    check("map: a browser with no DecompressionStream gets the sentence and a findings table through the veil", not bad, bad)
    b.call("Emulation.setDeviceMetricsOverride", {"width": 390, "height": 844, "deviceScaleFactor": 2, "mobile": False}, b.session)
    try:
        wide = {}
        for link in ("", "#chr16"):
            map_page(b, base, S, link)
            wide[link or "genome"] = too_wide(b)
    finally:
        b.call("Emulation.clearDeviceMetricsOverride", {}, b.session)
    check("map: at 390 px the map, open or not, scrolls nothing sideways", not any(wide.values()), wide)
    if os.environ.get("MAP_SCREENS"):
        print("map screenshots: " + ", ".join(map_screens(b, base, os.environ["MAP_SCREENS"])))


def screens_main(dest, site):
    """--map-screens DEST SITE: the map's screenshots of an already built demo (the docs workflow's artifact)."""
    chrome = find_chrome()
    if not chrome:
        print("ERROR: no Chrome or Chromium to take the screenshots")
        return 1
    handler = type("H", (Static,), {"__init__": lambda self, *a, **k: SimpleHTTPRequestHandler.__init__(
        self, *a, directory=site, **k)})
    srv = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    work = tempfile.mkdtemp(prefix="map-screens.")
    b = Chrome(chrome, os.path.join(work, "chrome"))
    try:
        taken = map_screens(b, f"http://127.0.0.1:{srv.server_address[1]}", dest)
    finally:
        b.close()
        srv.shutdown()
        srv.server_close()
        shutil.rmtree(work, ignore_errors=True)
    print("\n".join(taken))
    return 0 if len(taken) == 6 else 1


def main():
    work = tempfile.mkdtemp(prefix="test_live_demo.")
    try:
        # 1. the weblog fixture
        check("1. the weblog fixture is what make_demo_events.py writes now", make_demo_events.main(["--check"]) == 0)
        changed = os.path.join(work, "changed.jsonl")
        with open(build_live_demo.WEBLOG) as f:
            lines = f.readlines()
        with open(changed, "w") as f:
            f.writelines(lines[:-1])
        control("a fixture missing its last message", make_demo_events.main(["--check", "--out", changed]) == 0)

        out = os.path.join(work, "demo")
        og = og_binary() or build_live_demo.OG
        built = build_live_demo.build(out, work=os.path.join(work, "root"), og=og)
        frames = built["frames"]
        print(f"built: {len(frames)} frames, hooks {built['hooks']}")

        # 2. every plugin's state hook
        hooks = registered_hooks()
        check("2. the build calls every registered state hook", sorted(built["hooks"]) == sorted(hooks), (built["hooks"], hooks))
        bad = missing_hooks(built, hooks)
        check(f"2. every state hook ({', '.join(hooks)}) has a value on the demo's pages", not bad and len(hooks) >= 5, bad)
        dropped = build_live_demo.build(os.path.join(work, "demo-dropped"), work=os.path.join(work, "root-dropped"),
                                        drop_hook=hooks[0], og=og)
        control(f"a build that drops the {hooks[0]} hook", not missing_hooks(dropped, hooks))

        # 3. the frames follow the run
        steps = replayed_steps()
        bad = follow_problems(frames, steps)
        check(f"3. the frames follow the run ({len(frames)} frames, {len(steps)} steps)", not bad and len(steps) > 30, bad)
        early = copy.deepcopy(frames)
        early[0]["state"]["plugins"]["charts"] = early[-1]["state"]["plugins"]["charts"]
        control("charts shown before the report step ended", not follow_problems(early, steps))

        # 4. the invented machine, no path of this one
        ops = frames[-1]["state"]["plugins"]["ops"]
        disk = ops.get("disk") or [{}]
        check("4. the run-health card measures the invented machine",
              ops["memory"]["host_bytes"] == build_live_demo.DEMO_RAM
              and disk[0].get("size_bytes") == build_live_demo.DEMO_DISK.f_blocks * build_live_demo.DEMO_DISK.f_frsize,
              (ops["memory"], disk))
        with open(os.path.join(out, "runs", S, "frames.json")) as f:
            text = f.read()
        leaks = [p for p in (work, os.path.realpath(work), REPO, os.path.expanduser("~")) if p in text]
        check("4. no frame holds a path of the machine that built it", not leaks, leaks)

        # 5. links
        bad = link_problems(out, built)
        check("5. every link of the frames is relative and reaches a file of the demo", not bad, bad)
        planted = copy.deepcopy(built)
        planted["family"]["state"]["plugins"]["family"]["pedigree"] = f"/family/{FAM}/pedigree.svg"
        control("an absolute pedigree link", not link_problems(out, planted))
        bad = json_link_problems(out)
        check("5. every link of runs.json, results.json and family.json is relative and reaches a file", not bad, bad)
        with open(os.path.join(out, "runs.json")) as f:
            idx = json.load(f)
        with open(os.path.join(out, "runs.json"), "w") as f:
            json.dump(dict(idx, runs=[dict(r, url="/" + r["url"]) for r in idx["runs"]]), f)
        control("an absolute link in runs.json", not json_link_problems(out))
        with open(os.path.join(out, "runs.json"), "w") as f:
            json.dump(idx, f)

        # 6. the policy on every page, nothing from outside
        pages, texts = [], {}
        for p, _, fs in os.walk(out):
            for f in fs:
                rel = os.path.relpath(os.path.join(p, f), out)
                if f.endswith((".html", ".js", ".css")) and rel not in [os.path.join("runs", n, "report.html") for n in [S] + OTHERS]:
                    with open(os.path.join(p, f)) as fh:
                        texts[rel] = fh.read()
                    if f.endswith(".html"):
                        pages.append(rel)
        meta = f'<meta http-equiv="Content-Security-Policy" content="{build_live_demo.CSP}">'
        nometa = [p for p in pages if meta not in texts[p]]
        check(f"6. every page the build writes carries the policy ({len(pages)} pages)", not nometa and len(pages) >= 4, nometa)
        check("6. the policy is the server's, frame-ancestors aside",
              build_live_demo.CSP + "; frame-ancestors 'none'" == owngenome_live.CSP, build_live_demo.CSP)
        bad = asset_problems(texts)
        check("6. no page or script of the demo loads anything from outside", not bad, bad)

        # 7. in a browser
        browser_checks(out)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    if FAILS:
        print(f"\n{len(FAILS)} check(s) failed")
        return 1
    print("\nAll checks passed")
    return 0


if __name__ == "__main__":
    if sys.argv[1:2] == ["--map-screens"] and len(sys.argv) == 4:
        sys.exit(screens_main(sys.argv[2], sys.argv[3]))
    sys.exit(main())
