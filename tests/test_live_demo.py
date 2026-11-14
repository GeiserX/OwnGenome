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
     suggestion until each is opened, opening one brings back its own rows,
     marks and ranking only, and a reload closes them again; the home page and
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
     his page a CHIP call, a copy-number dot and a card's fact in the gene of a
     closed result stay out of sight (on a phone the card's line too) until it
     is opened; DEMO-001's page
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
position, and de novo rows that 'Show the map' opens; for the son's page, person.js that draws
every item, shows every fact or builds a card's line from every fact), the phone page (person.js that never folds, and a route that opens only the nearest fold), the
folds (a flagged section 5 card folded, a redraw that forgets the folds, a Ctrl+click that opens one) and
the person page's layout (a wire showing every line or every plain finish, a card showing
every line, a step timeline left open) also run against a planted mistake and must
fail there (2: a build with one plugin's state hook dropped): a check that cannot fail is no check.

The genome map of each person page (map_checks) is checked the same way, each check with its planted mistake: its
levels from the genome down to the reads and the addresses that open them, Esc and the zoom keys, the search box,
Nearest and the finding's card, Paint the copies, the veil and opening on purpose, and the frame budget on a
genome-scale document (perf_checks, with a 20 ms busy loop in draw() as its control).

Run: python3 tests/test_live_demo.py
"""
import base64
import copy
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
    """A static host: files only, no Content-Security-Policy header."""

    def log_message(self, *a):
        pass


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
        family_veil_checks(b, base, out)
        widget_veil_checks(b, base, out)
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
           # five copy-number genes from Paraphase, then Parascopy's SMN1 and SMN2
           ("DEMO-003", "paralogs", "copies", "[data-widget=copies] li", 7),
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


def planted_problems(b, base):
    """[problems] of the planted results on the son's page: closed, then opened one at a time, then reloaded."""
    gene, hit, sv = planted_words()
    out = planted_closed(b, base)
    ids = b.eval("Array.prototype.map.call(document.querySelectorAll('[data-finding][data-level=act]'), "
                 "function (e) { return e.getAttribute('data-finding'); }).join(' ')") or ""
    sf = next((i for i in ids.split() if i.startswith("sf-")), None)
    import make_demo_sample
    start = str(make_demo_sample.PLANTED["sv"][1])
    del_ = next((i for i in ids.split() if i.startswith("sv-") and start in i), None)
    if not sf or not del_:
        return out + [f"the planted results are not among the son's results to act on: {ids}"]
    b.eval("document.querySelector('#f-%s .veil button').click(); true" % sf)
    time.sleep(0.8)
    shown = b.eval(SHOWN_JS) or ""
    rows = b.eval("Array.prototype.filter.call(document.querySelectorAll('#card-clinvar tbody tr'), function (r) { "
                  "return r.textContent.indexOf(%s) >= 0; }).length" % json.dumps(gene))
    tick = b.eval("Array.prototype.filter.call(document.querySelectorAll('figure[data-chart=genome] a[data-mark^=hit-] "
                  "title'), function (t) { return t.textContent.indexOf(%s) >= 0; }).length" % json.dumps(gene))
    if rows != 1 or not tick or not all(w in shown for w in [gene, hit[0]]):
        out.append(f"opening the secondary finding did not bring back its ClinVar row and tick: {rows} rows, {tick} ticks")
    ranked = b.eval("Array.prototype.filter.call(document.querySelectorAll('[data-finding^=exomiser-]'), function (e) { "
                    "return e.textContent.indexOf(%s) >= 0; }).length" % json.dumps(gene))
    if not ranked:
        out.append(f"opening the secondary finding did not bring back the symptom ranking of {gene}")
    out += [f"{w} shown after opening only the secondary finding" for w in sv if w in shown]
    hits, _ = son_search(b, gene)
    if hits in (None, "0"):
        out.append(f"the search for {gene} finds nothing once its result is opened")
    b.eval("document.querySelector('#f-%s .veil button').click(); true" % del_)
    time.sleep(0.8)
    shown = b.eval(SHOWN_JS) or ""
    spans = b.eval("document.querySelectorAll('figure[data-chart=genome] a[data-mark^=sv-]').length")
    if sv[0] not in shown or spans != 2:
        out.append(f"opening the deletion did not bring back its span on the genome chart: {spans} spans")
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
    named = [k for k in ((res.get("veil") or {}).get("findings") or {}) if k.startswith("exomiser-")]
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
        bad = planted_problems(b, base)
        check("browser: the son's planted secondary finding and class 5 deletion name their gene and positions nowhere "
              "(cards, charts, titles, links, search) until each is opened, each opens only its own rows and marks, "
              "and a reload closes them", not bad, bad)
        res = json.loads(keep[rj])
        res["veil"] = {}
        write(rj, json.dumps(res))
        control("results.json without its veil tags", not planted_closed(b, base))
        write(rj, keep[rj])
        for p, a, z, desc in (
                (pj, "return !(tags[i] && tagClosed(tags[i]));", "return true;", "person.js whose ClinVar rows ignore their tags"),
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
        a = "!!tagOf(x) && !opened[tagOf(x)]"
        if a not in keep[fj]:
            check("control: finding.js whose findings ignore their tags (the planted line is in finding.js)", False, a)
        else:
            write(fj, keep[fj].replace(a, "false"))
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


def widget_veil_checks(b, base, out):
    """On the son's page a CHIP call and a copy-number dot in the gene of a closed result stay out of sight until that
    result is opened; the other dot stays drawn. A CHIP call is also a finding of its card that names its gene, so
    the CHIP card closes whole with the result; SMN1 at two copies is no finding, so only its dot's tag closes it."""
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
               "&& document.querySelector('[data-widget=copies]') && document.querySelector('.veil') && 'drawn'", 15)
        return shown()

    def phone_line():
        """The copy-number card's one line on a phone (its fold's summary, after the card's own title, the same for
        everyone), drawn from its first open fact."""
        b.call("Emulation.setDeviceMetricsOverride", {"width": 390, "height": 844, "deviceScaleFactor": 2, "mobile": False},
               b.session)
        try:
            b.go("about:blank")
            b.wait("location.href === 'about:blank'")
            b.go(f"{base}/runs/{OTHERS[1]}/")
            return b.wait("document.documentElement.dataset.pgpStatus === 'ready' && document.documentElement.dataset.pgpResults "
                          "&& (document.querySelector('[data-card=fold-paralogs] summary .what') || {}).textContent", 15) or ""
        finally:
            b.call("Emulation.clearDeviceMetricsOverride", {}, b.session)

    def phone_problems():
        line = phone_line()
        return [] if line and "SMN1" not in line else [f"the copy-number card's line on a phone: {line!r}"]

    def shown():
        return b.eval("JSON.stringify([" + SHOWN_JS + ", (document.querySelector('[data-widget=copies]') || {}).textContent || ''])")

    def problems(text, opened):
        page, cops = json.loads(text) if text else ["", ""]
        out_ = []
        if "SMN2" not in cops:
            out_.append(f"the open dot is not drawn: {cops[:120]}")
        if ("p.R175H" in page) != opened or ("SMN1" in cops) != opened:
            out_.append(f"the CHIP call or the SMN1 dot of the closed results {'not ' if opened else ''}shown: "
                        f"{'p.R175H' in page} | {cops[:200]}")
        # The copy-number card's facts name SMN1: closed in its margin and in its line of section 8.
        if ("SMN1 copies" in page) != opened:
            out_.append(f"the card's fact on SMN1 {'not ' if opened else ''}shown")
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
        S["paralogs"] = {"state": "ok", "data": {"genes": [{"gene": g, "copy_number": n, "status": "called", "detail": "."}
                                                           for g, n in (("SMN1", 2), ("SMN2", 2))]}}
        res = live_results.build(OTHERS[1], summary)
        with open(rj, "w") as f:
            json.dump(res, f)
        bad = problems(drawn(), False) + phone_problems()
        if not bad:
            for _ in range(2):  # the CPSR card holds both results: one click opens one
                b.eval("document.querySelector('#card-cpsr .veil button').click(); true")
                time.sleep(0.5)
            bad = problems(shown(), True)
        check("browser: on the son's page a CHIP call, a copy-number dot and a card's fact in the gene of a closed "
              "result stay out of sight (the fact in the card's line on a phone too) until it is opened", not bad, bad)
        for desc, a, z in (
                ("person.js that draws every item", "add(text, widget(of[1]));", "add(text, widget(c.widget));"),
                ("person.js that shows every fact", "facts = (c.facts || []).filter(function (f, i) { return !(ft[i] && tagClosed(ft[i])); });",
                 "facts = c.facts || [];"),
                ("person.js whose card line reads every fact",
                 "(openFacts(c)[0][0] ? openFacts(c)[0][0][0] + ': ' + openFacts(c)[0][0][1] : NOT_ALL_CLEAR)",
                 "(c.facts && c.facts[0] ? c.facts[0][0] + ': ' + c.facts[0][1] : NOT_ALL_CLEAR)")):
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
# once: nothing to wait for): the closed act cards in page order, the act chip's count, whether its confirm step
# shows, where the keyboard is (first: the chips' first button; canvas: the map), and the opened chips.
CHIP_JS = r"""JSON.stringify((function () {
  var ids = [];
  Array.prototype.forEach.call(document.querySelectorAll('.finding[data-level=act]'), function (a) {
    if (a.querySelector('.veil') && ids.indexOf(a.dataset.finding) < 0) { ids.push(a.dataset.finding); } });
  var c = document.querySelector('.gm-chip[data-closed=act]'), f = document.activeElement;
  var first = document.querySelector('.gm-chips button');
  return { closed: ids, buttons: !!first, chip: c ? c.firstChild.textContent : '', confirm: !!(c && c.querySelector('.gm-confirm')),
    focus: f && f.classList.contains('gm-yes') ? 'gm-yes' : f === first ? 'first' : f === document.querySelector('.gm-canvas') ? 'canvas'
      : (f ? f.tagName + '.' + f.className : ''),
    opened: Array.prototype.map.call(document.querySelectorAll('.gm-chip.opened'), function (o) { var a = o.querySelector('a');
      return { kind: o.dataset.opened, text: o.textContent, href: a ? a.getAttribute('href') : '' }; }) };
}()))"""


def chips_now(b):
    time.sleep(0.05)  # a chip's button draws the page and the chips again before its click returns
    return json.loads(b.eval(CHIP_JS))


def closed_words(n):
    return "" if n == 0 else "1 result closed" if n == 1 else f"{n} results closed"


def chip_open(b):
    """The act chip's Open, then its Show the result, as a keyboard reader does it. Returns (the id it opened,
    problems): one Open opens exactly one of the page's closed results, the count falls by one, the keyboard
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
    named = [o["href"] for o in after["opened"] if o["kind"] == "section" and o["href"][3:] in after["closed"]]
    if named:
        bad.append(f"an Opened in its section chip names a closed result: {named}")
    return (gone[0] if len(gone) == 1 else None), bad


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
    for k in range(want):
        got, more = chip_open(b)
        if more:
            return bad + [f"Open {k + 1}: {x}" for x in more]
        if got != order[k]:
            return bad + [f"Open {k + 1} opened {got}, not the page's next closed result {order[k]} ({order})"]
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
    card = b.eval("!!document.querySelector('#f-%s .veil')" % fid)
    if leaks(c) or not card:
        bad.append(f"Close it again: map shows {leaks(c)}, card closed {card}")
    b.eval("document.querySelector('#f-%s .veil button').click(); true" % fid)
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
        b.eval("document.querySelector('#f-%s .veil button').click(); true" % fid)
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
         ("the page's own CYP2D6 medicine finding over a closed variant naming another gene",
          [act("pl-theta", "PLTHETA", "chr22:42125000", "Planted theta result")], [], "pgx-cyp2d6", "chr22"))


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
        b.eval("document.querySelector('#f-pl-beta .veil button').click(); true")
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


# A finding the server tags with a closed result because it names its gene (results.veil.findings: a symptom
# ranking), once with og's tag on its row too (CFTR, inside the deletion's span) and once with none (BRCA1).
NAMED = (act("pl-x", "PLXGENE", "chr7:117400000-117700000", "Planted x deletion"),
         {"id": "pl-q", "level": "review", "card": "exomiser", "gene": "CFTR", "title": "Symptom ranking: PLXGENE ranks 1",
          "detail": "CFTR", "text": "A planted ranking that names the deletion's gene."},
         {"id": "pl-r", "level": "review", "card": "exomiser", "gene": "BRCA1", "title": "Symptom ranking: PLXGENE ranks 2",
          "detail": "BRCA1", "text": "A planted ranking that names the deletion's gene."})
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
        if tags.get("pl-q") != "pl-x" or tags.get("pl-r") != "pl-x" or row_tag(out, "pl-r"):
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
        b.eval("document.querySelector('#f-pl-x .veil button').click(); true")
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
                got, more = chip_open(b)
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
                b.eval("var o = document.querySelector('.gm-chip.opened[data-opened=section] a[href=\"#f/%s\"]'); "
                       "o.parentNode.parentNode.querySelector('button').click(); true" % gone)
                after = chips_now(b)
                if after["closed"] != [gone] or any(o["kind"] == "section" and gone in o["href"] for o in after["opened"]):
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


# The closer levels: (address, level, view the address keeps, words the title starts or ends with, genes drawn).
CLOSER = (("#chr11:4.5-6.5M", "region", "chr11:4.500-6.500M", "Chromosome 11, band 11p15", 10),
          ("#chr11:5.150-5.350M", "gene", "chr11:5.150-5.350M", "and its neighbours", 10),
          ("#chr11:5225001-5235000", "reads", "chr11:5225001-5235000", "and its neighbours", 1))
# The card, the Nearest button, the genes drawn and the search box's lines, read at once.
NEAR_JS = r"""JSON.stringify((function () {
  var f = document.getElementById('genome-map'), n = f.querySelector('.gm-nearest'), c = f.querySelector('.gm-card');
  return { near: n.hidden ? '' : n.textContent, card: c.hidden ? '' : c.textContent, genes: f.dataset.genes || '',
    crumb: f.querySelector('.gm-crumb').textContent,
    found: Array.prototype.map.call(f.querySelectorAll('.gm-find-list li'), function (li) {
      return { kind: li.dataset.kind || '', text: li.textContent }; }) };
}()))"""


def map_key(b, key, target="canvas"):
    """A key pressed on the map's canvas, or on the page's body."""
    sel = "document.querySelector('.gm-canvas')" if target == "canvas" else "document.body"
    b.eval("%s.dispatchEvent(new KeyboardEvent('keydown', {key: %s, bubbles: true})); true" % (sel, json.dumps(key)))


def closer_problems(b, base):
    """The region, gene and reads levels from their addresses (the title, the path, the genes lane, the address
    kept); Esc steps out one level at a time; the gauge's Region and Gene stops land on their level."""
    bad = []
    for link, lv, view, words, genes in CLOSER:
        c = map_page(b, base, S, link)
        time.sleep(0.8)  # the chromosome's gene file arrives after the first frame
        c, n = map_now(b), json.loads(b.eval(NEAR_JS))
        if c["level"] != lv or c["view"] != view:
            bad.append(f"{link}: level {c['level']}, view {c['view']}, want {lv} {view}")
        if not (c["title"].startswith(words) or c["title"].endswith(words)):
            bad.append(f"{link}: title {c['title']!r}")
        if not n["genes"].isdigit() or int(n["genes"]) < genes:
            bad.append(f"{link}: {n['genes']!r} genes drawn, want {genes} or more")
        if "11p15" not in n["crumb"] or not n["crumb"].startswith("Genome › chr11"):
            bad.append(f"{link}: path {n['crumb']!r}")
    c = map_page(b, base, S, "#chr11:5.150-5.350M")
    steps = []
    for _ in range(3):
        map_key(b, "Escape")
        steps.append(map_now(b)["level"])
    if steps != ["region", "chromosome", "genome"]:
        bad.append(f"Esc from the gene level steps through {steps}")
    for stop, lv in (("region", "region"), ("gene", "gene")):
        map_page(b, base, S)
        b.eval("document.querySelector('.gm-gauge [data-stop=%s]').click(); true" % stop)
        c = map_now(b)
        if c["level"] != lv:
            bad.append(f"the {stop} stop from the genome: level {c['level']}")
    c = map_page(b, base, S, "#chr11:100.0-120.0M")
    map_key(b, "+")
    z = map_now(b)
    map_key(b, "ArrowRight")
    r = map_now(b)
    if z["view"] != "chr11:105.0-115.0M" or r["view"] != "chr11:106.0-116.0M":
        bad.append(f"+ then → from chr11:100.0-120.0M: {z['view']}, then {r['view']}")
    return bad


def map_find(b, q):
    """Q typed in the map's search box: its lines once they settled (the gene names arrive with the first key)."""
    b.eval("var i = document.getElementById('gm-find'); i.focus(); i.value = %s; "
           "i.dispatchEvent(new Event('input')); true" % json.dumps(q))
    time.sleep(0.8)
    b.eval("document.getElementById('gm-find').dispatchEvent(new Event('input')); true")
    time.sleep(0.2)
    return json.loads(b.eval(NEAR_JS))["found"]


def find_problems(b, base, out):
    """The search box: / puts the cursor in it; a gene, a band and a finding are found and Enter flies there (a
    finding with its card); on the son's closed page his closed result is never a line, by its gene, its id or
    its title, while his open finding still is."""
    bad = []
    map_page(b, base, S)
    map_key(b, "/", target="body")
    if b.eval("document.activeElement && document.activeElement.id") != "gm-find":
        bad.append("/ does not put the cursor in the map's search box")
    for q, kind, words, chrom, lv in (("HBB", "gene", "HBB · gene · chr11", "chr11", ("gene", "reads")),
                                      ("11p15.4", "band", "11p15.4 · band", "chr11", ("chromosome", "region")),
                                      ("DEMOC", "finding", "DEMOC · ", "chr16", ("gene", "reads"))):
        map_page(b, base, S)
        found = map_find(b, q)
        if not found or found[0]["kind"] != kind or not found[0]["text"].startswith(words):
            bad.append(f"{q}: lines {found}")
            continue
        b.eval("document.getElementById('gm-find').dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', bubbles: true})); true")
        c, n = map_now(b), json.loads(b.eval(NEAR_JS))
        if c["chrom"] != chrom or c["level"] not in lv:
            bad.append(f"{q}: Enter flies to {c['chrom']} at {c['level']}, want {chrom} at {lv}")
        if kind == "finding" and "DEMOC" not in n["card"]:
            bad.append(f"{q}: Enter opens no card for it: {n['card']!r}")
    son = OTHERS[1]
    doc = read_map(os.path.join(out, "runs", son, "map.json"))
    f = doc["per"][son]["findings"]
    i = [j for j, lv in enumerate(f["level"]) if lv == "act" and not f["gate"][j]][0]
    fid, gene, title = f["id"][i], f["gene"][i], f["title"][i]
    map_page(b, base, son)
    for q in (gene, fid, title):
        found = map_find(b, q)
        named = [x for x in found if (x["kind"] == "finding" and gene in x["text"]) or fid in x["text"] or title in x["text"]]
        if named:
            bad.append(f"the son's closed result is a line for {q!r}: {named}")
    found = map_find(b, "DEMOC")
    if not any(x["kind"] == "finding" for x in found):
        bad.append(f"the son's open finding DEMOC is not found: {found}")
    return bad


def nearest_problems(b, base, out):
    """Nearest: with no finding of the open chromosome in view it names the closest one (else the first on the
    next chromosome that has one) and flies there with its card; Esc closes the card; on the son's closed page
    it never leads to his closed result."""
    bad = []
    for link, words, chrom in (("#chr11:100.0-101.0M", "Nearest: DEMOB at 108.3 Mb ›", "chr11"),
                               ("#chr3", "Nearest: DEMOB on chr11 ›", "chr11")):
        map_page(b, base, S, link)
        # The margin, Nearest with it, follows the settled camera; on a busy machine that takes longer.
        b.wait("!document.querySelector('#genome-map .gm-nearest').hidden", 5)
        n = json.loads(b.eval(NEAR_JS))
        if n["near"] != words:
            bad.append(f"{link}: Nearest reads {n['near']!r}, want {words!r}")
            continue
        b.eval("document.querySelector('.gm-nearest').click(); true")
        c = map_now(b)
        for _ in range(10):  # a busy machine eases the camera slower
            if any("demob" in m["id"] for m in c["marks"]):
                break
            c = map_now(b)
        n = json.loads(b.eval(NEAR_JS))
        if c["chrom"] != chrom or not any("demob" in m["id"] for m in c["marks"]) or "DEMOB" not in n["card"]:
            bad.append(f"{link}: Nearest flew to {c['chrom']} {c['view']}, marks {[m['id'] for m in c['marks']]}, "
                       f"card {n['card']!r}")
        map_key(b, "Escape")
        time.sleep(0.2)
        if json.loads(b.eval(NEAR_JS))["card"] or map_now(b)["level"] == "chromosome":
            bad.append(f"{link}: Esc leaves the card open or steps out with a card open")
    son = OTHERS[1]
    doc = read_map(os.path.join(out, "runs", son, "map.json"))
    f = doc["per"][son]["findings"]
    i = [j for j, lv in enumerate(f["level"]) if lv == "act" and not f["gate"][j]][0]
    others = [f["s"][j] for j in range(len(f["id"])) if f["chrom"][j] == f["chrom"][i] and j != i]
    # A view next to his closed result, with every other finding of that chromosome farther away.
    a = max(0, f["s"][i] - 3_000_000)
    if not any(abs(x - f["s"][i]) < 6_000_000 for x in others):
        n = map_page(b, base, son, f"#{f['chrom'][i]}:{a + 1}-{a + 2_000_000}") and json.loads(b.eval(NEAR_JS))
        if f["gene"][i] in n["near"] or f["gene"][i] in n["card"]:
            bad.append(f"the son's Nearest leads to his closed result: {n['near']!r}")
    else:
        bad.append("the son's closed result has another finding within 6 Mb: no view to check Nearest on")
    return bad


def origin_problems(b, base, out):
    """Paint the copies: Origin is off with its reason, on a run in a family and on one in none (DEMO-001's
    frames with the family plugin taken out)."""
    bad = []
    js = "JSON.stringify([document.querySelector('.gm-paint [data-paint=origin]').getAttribute('aria-disabled'), " \
         "document.querySelector('.gm-paint-why').textContent])"
    c = map_page(b, base, S)
    got = json.loads(b.eval(js))
    if got != ["true", "Origin: which parent each stretch came from is painted on the family page."]:
        bad.append(f"in a family: {got}")
    path = os.path.join(out, "runs", S, "frames.json")
    with open(path, "rb") as fh:
        keep = fh.read()
    try:
        frames = json.loads(keep)
        for fr in frames["frames"] if isinstance(frames, dict) else frames:
            ((fr.get("state") or {}).get("plugins") or {}).pop("family", None)
        with open(path, "w") as fh:
            json.dump(frames, fh)
        c = map_page(b, base, S)
        got = json.loads(b.eval(js))
        if got != ["true", "Origin needs a family run with a sequenced parent: this run is in none."]:
            bad.append(f"in no family: {got}")
    finally:
        with open(path, "wb") as fh:
            fh.write(keep)
    return bad


def genome_scale(out):
    """DEMO-001's map document and chr1's gene file made genome-scale, as a 30x genome's are: about 4,000 runs
    of homozygosity over the 24 chromosomes (every 50th one 5 Mb or longer) and 5,000 genes with 40,000 exon
    blocks on chr1. Invented, made by arithmetic, so the same every time. Returns what to write back."""
    mp, gp = os.path.join(out, "runs", S, "map.json"), os.path.join(out, "map", "genes", "chr1.json")
    keep = {}
    for p in (mp, gp):
        with open(p, "rb") as fh:
            keep[p] = fh.read()
    doc = read_map(mp)
    chroms = doc["ref"]["chroms"]
    total = sum(c["len"] for c in chroms)
    roh, bins, k = {}, {"size": 1_000_000}, 0
    for c in chroms:
        n = round(4000 * c["len"] / total)
        step = c["len"] // (n + 1)
        rows = {"ds": [], "l": [], "markers": [], "q": [], "long": [], "art": []}
        cov = [0] * (c["len"] // 1_000_000 + 1)
        last = 0
        for j in range(n):
            k += 1
            s0 = (j + 1) * step - step // 2
            ln = 6_000_000 if k % 50 == 0 else 100_000 + (k * 7919) % 1_300_000
            ln = min(ln, c["len"] - s0 - 1)
            rows["ds"].append(s0 - last)
            last = s0
            rows["l"].append(ln)
            rows["markers"].append(ln // 1400 + 1)
            rows["q"].append(50.0)
            rows["long"].append(1 if ln >= 5_000_000 else 0)
            rows["art"].append(0)
            for bi in range(s0 // 1_000_000, (s0 + ln - 1) // 1_000_000 + 1):
                cov[bi] += min(s0 + ln, (bi + 1) * 1_000_000) - max(s0, bi * 1_000_000)
        roh[c["n"]] = rows
        bins[c["n"]] = cov
    doc["per"][S]["roh"], doc["per"][S]["roh_bins"] = roh, bins
    doc["layers"]["roh"] = dict(doc["layers"]["roh"], state="ok", count=k)
    write_map(mp, doc)
    length = chroms[0]["len"]
    g = {"ds": [], "l": [], "strand": [], "name": [], "type": []}
    blocks, last = [], 0
    for j in range(5000):
        s0 = 10_000 + j * (length - 200_000) // 5000
        ln = 20_000 + (j * 104_729) % 100_000
        g["ds"].append(s0 - last)
        last = s0
        g["l"].append(ln)
        g["strand"].append("+" if j % 2 else "-")
        g["name"].append(f"SYNG{j + 1}")
        g["type"].append("protein_coding" if j % 3 else "lncRNA")
        blocks += [(s0 + x * ln // 8, s0 + x * ln // 8 + 300 + (j * 31 + x * 97) % 900) for x in range(8)]
    blocks.sort()
    merged = []
    for a, z in blocks:
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], z)
        else:
            merged.append([a, z])
    x, last = {"ds": [], "l": []}, 0
    for a, z in merged:
        x["ds"].append(a - last)
        x["l"].append(z - a)
        last = a
    write_map(gp, {"v": 1, "gencode": "v50", "chrom": "chr1", "g": g, "x": x})
    return keep


# Inside the page, before its scripts: the time each animation frame's callbacks take (the map's draw()), plus
# the time the pointer handlers took since the last frame. A pointer event is timed from a listener on window
# that runs before the canvas's (capture) to one that runs after it (bubble).
FRAME_HOOK = r"""(function () {
  var raf = window.requestAnimationFrame.bind(window), pending = 0, t0 = 0;
  window.__mapFrames = [];
  window.addEventListener('pointermove', function () { t0 = performance.now(); }, true);
  window.addEventListener('pointermove', function () { if (t0) { pending += performance.now() - t0; t0 = 0; } }, false);
  window.requestAnimationFrame = function (cb) {
    return raf(function (t) { var a = performance.now(); cb(t); window.__mapFrames.push(performance.now() - a + pending); pending = 0; });
  };
}());"""
PAN_PLACES = ("#chr1:100.0-160.0M", "#chr1:100.0-105.0M", "#chr1:100.000-100.200M")


def map_pan(b, base):
    """A scripted pan of 300 steps on DEMO-001's map, 100 at each of the chromosome, region and gene levels:
    the work of each frame in ms (pointer handlers plus draw()), the first paint and parse times, the JS heap in
    MB and the levels and genes the pan drew at."""
    b.call("Page.enable", session=b.session)
    sid = b.call("Page.addScriptToEvaluateOnNewDocument", {"source": FRAME_HOOK}, b.session)["identifier"]
    work, seen = [], []
    try:
        first = map_page(b, base, S)
        for place in PAN_PLACES:
            b.eval(f"location.hash = '{place}'; true")
            time.sleep(1.5)  # settled, and the gene file of chr1 in
            r = json.loads(b.eval("document.querySelector('.gm-canvas').scrollIntoView({block: 'center'}); "
                                  "JSON.stringify(document.querySelector('.gm-canvas').getBoundingClientRect())"))
            time.sleep(0.3)
            x, y = r["left"] + r["width"] / 2, r["top"] + r["height"] - 30
            b.eval("window.__mapFrames = []; true")
            b.call("Input.dispatchMouseEvent", {"type": "mousePressed", "x": x, "y": y, "button": "left", "buttons": 1,
                                                "clickCount": 1}, b.session)
            for i in range(100):
                x += 4 if (i // 25) % 2 == 0 else -4
                b.call("Input.dispatchMouseEvent", {"type": "mouseMoved", "x": x, "y": y, "button": "left", "buttons": 1},
                       b.session)
                time.sleep(0.02)
            b.call("Input.dispatchMouseEvent", {"type": "mouseReleased", "x": x, "y": y, "button": "left", "buttons": 0,
                                                "clickCount": 1}, b.session)
            time.sleep(0.3)
            work += json.loads(b.eval("JSON.stringify(window.__mapFrames)"))
            f = json.loads(b.eval("JSON.stringify([document.getElementById('genome-map').dataset.level, "
                                  "document.getElementById('genome-map').dataset.genes])"))
            seen.append(f)
        heap = b.eval("performance.memory ? performance.memory.usedJSHeapSize : -1")
    finally:
        b.call("Page.removeScriptToEvaluateOnNewDocument", {"identifier": sid}, b.session)
    return {"work": work, "paint": first.get("paint"), "parse": first.get("parse"), "heap": heap / 1e6, "seen": seen}


def pct(xs, q):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(q * len(xs)))] if xs else float("inf")


def perf_problems(got):
    """The budget: at least 200 frames drawn, through the chromosome, region and gene levels with genes drawn
    below the chromosome level; work per frame median 8 ms and 95th percentile 16 ms at most; first paint within
    150 ms of the document arriving, parsed within 30 ms; JS heap 60 MB at most."""
    bad, w = [], got["work"]
    if len(w) < 200:
        bad.append(f"the pan drew {len(w)} frames, want 200 or more")
    if [s[0] for s in got["seen"]] != ["chromosome", "region", "gene"] or not all(int(s[1] or 0) > 0 for s in got["seen"][1:]):
        bad.append(f"the pan ran at {got['seen']}")
    if pct(w, 0.5) > 8 or pct(w, 0.95) > 16:
        bad.append(f"work per frame: median {pct(w, 0.5):.2f} ms, p95 {pct(w, 0.95):.2f} ms")
    if not (got["paint"] <= 150 and got["parse"] <= 30):
        bad.append(f"first paint {got['paint']} ms, parse {got['parse']} ms")
    if not 0 < got["heap"] <= 60:
        bad.append(f"JS heap {got['heap']:.1f} MB")
    return bad


def perf_checks(b, base, out):
    """The frame budget on a genome-scale document, and the same check with a 20 ms busy loop in draw()."""
    keep = genome_scale(out)
    mjs = os.path.join(out, "assets", "map.js")
    with open(mjs) as fh:
        keep_mjs = fh.read()
    head = "function draw() {\n      raf = 0;"
    try:
        got = map_pan(b, base)
        w = got["work"]
        print(f"map pan: {len(w)} frames, work median {pct(w, 0.5):.2f} ms, p95 {pct(w, 0.95):.2f} ms, max "
              f"{max(w or [0]):.2f} ms; first paint {got['paint']} ms, parse {got['parse']} ms, heap {got['heap']:.1f} MB")
        bad = window_problems(b, base)
        check("map: past 200 rows the In view table's Show all is a scrolled window that shows every mark drawn",
              not bad, bad)
        win = "var a = clamp(Math.floor(inWrap.scrollTop / ROW_PX) - 10, 0, cap)"
        with open(mjs, "w") as fh:
            fh.write(keep_mjs.replace(win, "var a = 0"))
        control("map.js whose scrolled window always starts at the first row", keep_mjs.count(win) != 1 or
                not window_problems(b, base))
        with open(mjs, "w") as fh:
            fh.write(keep_mjs)
        bad = perf_problems(got)
        check("map: on a genome-scale document, a 300-step pan through the chromosome, region and gene levels works "
              "each frame in 8 ms at the median and 16 ms at the 95th percentile, first paint within 150 ms (parse "
              "30 ms), JS heap 60 MB at most", not bad, bad)
        with open(mjs, "w") as fh:
            fh.write(keep_mjs.replace(head, "function draw() {\n      var bz = performance.now(); "
                                            "while (performance.now() - bz < 20) { void 0; }\n      raf = 0;"))
        slow = map_pan(b, base)
        print(f"map pan, busy loop: work median {pct(slow['work'], 0.5):.2f} ms, p95 {pct(slow['work'], 0.95):.2f} ms")
        control("map.js with a 20 ms busy loop in draw()", keep_mjs.count(head) != 1 or not perf_problems(slow))
    finally:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs)
        unplant(keep)


# --- the whole-genome strip and the In view table (map.js drawsAt(), strip, table) -------------------

# drawsAt() of the page's map.js asked for every row of each run's document and every gene of the demo's gene
# files, as og map build --predicate-table writes its own answers: one tab-separated line per row.
PRED_JS = r"""(async function (runs) {
  var LV = window.PGPMap.LEVELS, on = window.PGPMap.drawsAt, out = {}, genes = [];
  async function gz(u) {
    var r = await fetch(u, { cache: 'no-store' });
    if (!r.ok) { throw new Error(u + ': ' + r.status); }
    return JSON.parse(await new Response(r.body.pipeThrough(new DecompressionStream('gzip'))).text());
  }
  function line(lines, layer, c, s, e, id, row) {
    lines.push([layer, c, s, e, id].concat(LV.map(function (lv) { return on(layer, row, lv) ? 1 : 0; })).join('\t'));
  }
  function each(o, fn) { var at = 0; (o.ds || []).forEach(function (d, i) { at += d; fn(at, at + o.l[i], i); }); }
  var chroms = [];
  for (var i = 1; i <= 22; i++) { chroms.push('chr' + i); }
  chroms.push('chrX', 'chrY');
  for (var c = 0; c < chroms.length; c++) {
    var g = (await gz('/map/genes/' + chroms[c] + '.json')).g || {};
    each(g, function (s, e, j) { line(genes, 'genes', chroms[c], s, e, g.name[j], { name: g.name[j] }); });
  }
  for (var k = 0; k < runs.length; k++) {
    var d = await gz('/runs/' + runs[k] + '/map.json'), per = d.per[runs[k]], f = per.findings, lines = [];
    f.id.forEach(function (id, j) {
      line(lines, 'findings', f.chrom[j], f.s[j], f.e[j], id, { id: id, level: f.level[j], s: f.s[j], e: f.e[j] });
    });
    Object.keys(per.roh || {}).forEach(function (cn) {
      var r = per.roh[cn];
      each(r, function (s, e, j) { line(lines, 'roh', cn, s, e, '', { long: r.long[j] }); });
    });
    out[runs[k]] = lines.concat(genes);
  }
  return JSON.stringify(out);
})(%s)"""

# A run of homozygosity of 7 Mb planted in DEMO-001's step 11 text for the predicate check: the demo's own runs
# are all shorter than 5 Mb, and a long one is the row whose answer differs at the genome level.
LONG_RG = "RG\t{run}\tchr2\t60000001\t67000000\t7000000\t5000\t80.0\n"


def predicate_lines(og, out, run):
    """(exit code, lines) of og map build --predicate-table for RUN of the demo, with the results its page gets."""
    root = os.path.join(os.path.dirname(out), "root")
    with open(os.path.join(out, "runs", run, "results.json"), "rb") as fh:
        body = b'{"results": ' + fh.read() + b"}"
    p = subprocess.run([og, "map", "build", "--root", root, "--run", run, "--assets", os.path.join(REPO, "assets"),
                        "--cache-dir", os.path.join(root, ".owngenome-live", "map"), "--predicate-table"],
                       input=body, capture_output=True, check=False)
    return p.returncode, p.stdout.decode().splitlines()


def plant_long_run(out, og):
    """DEMO-001's step 11 text with LONG_RG added and its map.json built again from it. Returns what to write back."""
    root = os.path.join(os.path.dirname(out), "root")
    txt, mj = os.path.join(root, S, "vcf", f"{S}_roh.txt"), os.path.join(out, "runs", S, "map.json")
    keep = {}
    for p in (txt, mj):
        with open(p, "rb") as fh:
            keep[p] = fh.read()
    with open(txt, "ab") as fh:
        fh.write(LONG_RG.format(run=S).encode())
    with open(os.path.join(out, "runs", S, "results.json"), "rb") as fh:
        body = b'{"results": ' + fh.read() + b"}"
    doc = build_live_demo.og_map(og, ["build", "--root", root, "--run", S, "--assets", os.path.join(REPO, "assets"),
                                      "--cache-dir", os.path.join(root, ".owngenome-live", "map")], body)
    with open(mj, "wb") as fh:
        fh.write(doc)
    return keep


def predicate_problems(b, base, out, og):
    """drawsAt() of map.js against og map build --predicate-table, row by row, on the three demo runs (DEMO-001
    with a planted long run): the same rows and the same answer at each level; and the comparison covers a finding,
    a long run, a shorter run and a gene, so it cannot pass on rows that agree by chance."""
    bad, runs = [], [S] + OTHERS
    map_page(b, base, S)
    r = b.call("Runtime.evaluate", {"expression": PRED_JS % json.dumps(runs), "awaitPromise": True, "returnByValue": True},
               b.session)
    js = json.loads(r.get("result", {}).get("value") or "{}")
    seen = set()
    for run in runs:
        code, lines = predicate_lines(og, out, run)
        if code != 0 or not lines or lines[0] != "layer\tchrom\ts\te\tid\tgenome\tchromosome\tregion\tgene\treads":
            bad.append(f"{run}: og map build --predicate-table exit {code}, header {lines[:1]}")
            continue
        want, got = sorted(lines[1:]), sorted(js.get(run) or [])
        if want != got:
            bad.append(f"{run}: og and drawsAt() differ: og only {sorted(set(want) - set(got))[:4]}, "
                       f"drawsAt() only {sorted(set(got) - set(want))[:4]}")
        for x in want:
            c = x.split("\t")
            seen.add((c[0], c[5]))
    if not {("findings", "1"), ("roh", "1"), ("roh", "0"), ("genes", "0")} <= seen:
        bad.append(f"the rows compared miss a finding, a long run, a shorter run or a gene: {sorted(seen)}")
    return bad


# The In view table, the canvas and the strip, read at once: the table's rows (Show all first), what the last
# frame drew one by one (canvas.pgpDrawn), the strip's ticks (null while it is hidden), the count, the line it
# announces, the density lines, the sort, each row's cells and the current row.
INVIEW_JS = r"""JSON.stringify((function () {
  var f = document.getElementById('genome-map'), cv = f.querySelector('.gm-canvas'), st = f.querySelector('.gm-strip');
  var all = f.querySelector('.gm-showall');
  if (all && all.textContent === 'Show all') { all.click(); }
  var trs = f.querySelectorAll('.gm-table tbody tr[data-key]');
  return { rows: Array.prototype.map.call(trs, function (r) { return r.dataset.key; }), drawn: (cv.pgpDrawn || []).slice(),
    ticks: st.hidden ? null : (st.pgpDrawn || []).slice(), stripH: st.hidden ? 0 : st.clientHeight,
    count: (f.querySelector('.gm-inview-n') || {}).textContent || '', live: f.querySelector('.gm-inview-live').textContent,
    dense: f.querySelector('.gm-inview-dense').textContent, open: f.querySelector('.gm-inview').open,
    summary: f.querySelector('.gm-inview > summary').textContent,
    sort: Array.prototype.map.call(f.querySelectorAll('.gm-table th'), function (t) { return t.getAttribute('aria-sort') || ''; }),
    cells: Array.prototype.map.call(trs, function (r) { return Array.prototype.map.call(r.cells, function (c) { return c.textContent; }); }),
    current: Array.prototype.map.call(f.querySelectorAll('.gm-table tr[aria-current=true]'), function (r) { return r.dataset.key; }),
    card: (f.querySelector('.gm-card') || {}).textContent || '', text: f.querySelector('.gm-inview').textContent };
}()))"""


def inview_now(b):
    return json.loads(b.eval(INVIEW_JS))


def inview_page(b, base, run, link=""):
    """RUN's page at LINK, settled with the gene file in when the view is close enough for genes: (map, table)."""
    c = map_page(b, base, run, link)
    if c.get("level") in ("region", "gene", "reads"):
        time.sleep(0.8)  # the chromosome's gene file arrives after the first frame
        c = map_now(b)
    return c, inview_now(b)


SETTLED_JS = ("JSON.stringify((function () { var f = document.getElementById('genome-map'); return [f.dataset.view || '', "
              "f.dataset.level || '', (f.querySelector('.gm-inview-n') || {}).textContent || '', "
              "(f.querySelector('.gm-canvas').pgpDrawn || []).length, f.dataset.genes || '']; }()))")


def settled_at(b, view):
    """Waits until the map has settled on VIEW (its address as the map writes it): the margin and the table follow
    the camera, then the gene file of the view if it needs one. (map, table)."""
    end, last = time.time() + 8, None
    while time.time() < end:
        time.sleep(0.25)
        now = json.loads(b.eval(SETTLED_JS))
        if now[0] == view and now == last:
            break
        last = now
    return json.loads(b.eval(MAP_JS)), inview_now(b)


def go_view(b, link):
    """The open page moved to LINK (an address as the map writes it) without a reload, settled."""
    b.eval("location.hash = %s; true" % json.dumps(link))
    return settled_at(b, link[1:])


def count_words(n):
    return "1 feature" if n == 1 else f"{n:,} features"


# Views at every level, several per level: a page load for each run's genome level, then the rest by address.
INVIEW_VIEWS = ((S, "", "genome"), (S, "#chr11", "chromosome"), (S, "#chr22", "chromosome"),
                (S, "#chr22:40.5-43.5M", "region"), (S, "#chr11:4.0-7.0M", "region"), (S, "#chr19:43.5-46.5M", "region"),
                (S, "#chr11:5.150-5.350M", "gene"), (S, "#chr22:42.000-42.200M", "gene"),
                (S, "#chr11:5225001-5250000", "reads"), (S, "#chr22:42110001-42140000", "reads"),
                (OTHERS[1], "", "genome"), (OTHERS[1], "#chr16", "chromosome"))


def inview_problems(b, base, views=INVIEW_VIEWS):
    """On views at every level: the table's rows are exactly what the canvas drew one by one, each once, its count
    says how many, and the shorter runs drawn as a band have a line of their own at the genome level; the strip is
    hidden at the genome level and over an open chromosome ticks exactly the findings the genome level draws; T
    folds the table and opens it again."""
    bad, genome, page = [], {}, None
    for run, link, lv in views:
        if run != page:
            c, t = inview_page(b, base, run, link)
            page = run
        else:
            c, t = go_view(b, link)
        where = f"{run}{link or ' genome'}"
        if c.get("level") != lv:
            bad.append(f"{where}: level {c.get('level')}, want {lv}")
            continue
        if sorted(t["rows"]) != sorted(t["drawn"]) or len(set(t["rows"])) != len(t["rows"]):
            bad.append(f"{where}: table only {sorted(set(t['rows']) - set(t['drawn']))[:3]}, "
                       f"canvas only {sorted(set(t['drawn']) - set(t['rows']))[:3]}")
        if t["count"] != count_words(len(t["rows"])):
            bad.append(f"{where}: the count says {t['count']!r} for {len(t['rows'])} rows")
        if lv == "genome" and not any(k.startswith("findings\t") for k in t["rows"]):
            bad.append(f"{where}: no finding in the table: {t['rows'][:3]}")
        if lv in ("region", "gene", "reads") and not any(k.startswith("genes\t") for k in t["rows"]):
            bad.append(f"{where}: no gene in the table")
        if lv == "genome":
            genome[run] = sorted(k for k in t["rows"] if k.startswith("findings\t"))
            if t["ticks"] is not None or t["stripH"]:
                bad.append(f"{where}: the strip shows at the genome level")
            if "shorter runs of homozygosity, drawn as a faint band" not in t["dense"]:
                bad.append(f"{where}: no line for the shorter runs: {t['dense']!r}")
        elif t["ticks"] is None or t["stripH"] != 34 or (run in genome and sorted(t["ticks"]) != genome[run]):
            bad.append(f"{where}: strip {t['stripH']} px, ticks {t['ticks']}, want {genome.get(run)}")
    if views is INVIEW_VIEWS:
        folds = []
        for _ in range(2):
            map_key(b, "t", target="body")
            folds.append(inview_now(b)["open"])
        if folds != [False, True]:
            bad.append(f"T folds and opens the table: {folds}")
    return bad


def inview_veil_problems(b, base, adult=True):
    """A closed result is never a row of the table or a tick of the strip: APOE on DEMO-001's page (unless not
    ADULT), the son's closed result to act on on his; once the son's is opened from the chip it is both, so the
    check sees it."""
    bad = []
    son = OTHERS[1]
    if adult:
        bad += inview_apoe_problems(b, base)
    c, t = inview_page(b, base, son, "#chr7")
    closed = [k for k in t["rows"] + (t["ticks"] or []) if k.startswith("findings\tchr7\t")]
    if closed or "KCNH2" in t["text"]:
        bad.append(f"{son}: his closed result shows at chr7: {closed}")
    # The chip opens his closed results one per click, in the page's order: open them until the one at chr7 is.
    for _ in range(4):
        _, cbad = chip_open(b)
        bad += cbad
        time.sleep(0.5)
        t = inview_now(b)
        if any(k.startswith("findings\tchr7\t") for k in t["rows"]) or cbad:
            break
    if not any(k.startswith("findings\tchr7\t") for k in t["rows"]) or \
            not any(k.startswith("findings\tchr7\t") for k in (t["ticks"] or [])):
        bad.append(f"{son}: opened, his result is not a row and a tick at chr7: {t['rows']}, {t['ticks']}")
    return bad


def inview_apoe_problems(b, base):
    """APOE, closed on every page until opened, is no row and no tick on DEMO-001's page."""
    bad = []
    c, t = inview_page(b, base, S)
    if any("\tapoe" in k for k in t["rows"]):
        bad.append(f"{S}: APOE, closed, is a row: {t['rows']}")
    c, t = go_view(b, "#chr19")
    if any("\tapoe" in k for k in (t["ticks"] or []) + t["rows"]):
        bad.append(f"{S}: APOE, closed, is a strip tick or a row at chr19")
    return bad


def size_bp(words):
    n, unit = words.replace(",", "").split()
    return float(n) * {"bp": 1, "kb": 1e3, "Mb": 1e6}[unit]


def sort_problems(b, base):
    """The table's order: by place with the findings first, aria-sort ascending on Where; Where again reverses it
    (descending); Size sorts by size, smallest first, with aria-sort moved there."""
    bad = []
    c, t = inview_page(b, base, S, "#chr22:40.5-43.5M")
    keys = t["rows"]
    place = sorted(keys, key=lambda k: (int(k.split("\t")[2]), int(k.split("\t")[3]), k))
    first = sorted([k for k in place if k.startswith("findings\t")], key=place.index)
    if keys != first + [k for k in place if k not in first] or place[0] in first:
        bad.append(f"default order {keys[:4]}: want the findings first, then by place (and a gene before the finding)")
    if t["sort"] != ["", "", "ascending", "", "", "", ""]:
        bad.append(f"default aria-sort {t['sort']}")
    b.eval("document.querySelector('.gm-table th button[data-col=where]').click(); true")
    r = inview_now(b)
    if r["rows"] != keys[::-1] or r["sort"][2] != "descending":
        bad.append(f"Where again: {r['rows'][:3]}, aria-sort {r['sort']}")
    b.eval("document.querySelector('.gm-table th button[data-col=size]').click(); true")
    r = inview_now(b)
    sizes = [size_bp(x[3]) for x in r["cells"]]
    if sizes != sorted(sizes) or r["sort"] != ["", "", "", "ascending", "", "", ""]:
        bad.append(f"Size: {sizes[:5]}, aria-sort {r['sort']}")
    return bad


ANNOUNCE_JS = ("window.__said = []; new MutationObserver(function () { window.__said.push(document.querySelector("
               "'.gm-inview-live').textContent); }).observe(document.querySelector('.gm-inview-live'), "
               "{ childList: true, characterData: true, subtree: true }); true")


def announce_problems(b, base):
    """The table announces its count once when it changes and never when it does not: three pans along a stretch
    with nothing in view say nothing; then a view with four marks says "4 features in view", once."""
    bad = []
    c, t = inview_page(b, base, S, "#chr3:100.0-101.0M")
    b.eval(ANNOUNCE_JS)
    for _ in range(3):
        map_key(b, "ArrowRight")
        time.sleep(0.4)
    settled_at(b, "chr3:100.300-101.300M")
    said = b.eval("JSON.stringify(window.__said)")
    if t["count"] != "0 features" or said != "[]":
        bad.append(f"panning an empty stretch ({t['count']}) announced {said}")
    c, t = go_view(b, "#chr11")
    said = json.loads(b.eval("JSON.stringify(window.__said)"))
    if said != [f"{t['count']} in view"] or t["count"] != "4 features":
        bad.append(f"#chr11 ({t['count']}) announced {said}")
    return bad


def fly_problems(b, base):
    """A row's button flies the map to its feature: a finding's row from the genome level to its place with its
    card, the row then current; a run's row to a view that holds the run, with no card."""
    bad = []
    c, t = inview_page(b, base, S)
    key = next((k for k in t["rows"] if k.startswith("findings\tchr16\t")), None)
    if not key:
        return [f"no chr16 finding row at the genome level: {t['rows']}"]
    b.eval("document.querySelector('.gm-table tr[data-key=\"%s\"] .gm-go').click(); true" % key.replace("\t", "\\t"))
    for _ in range(6):  # a busy machine eases the camera slower
        c = map_now(b)
        if c["chrom"] == "chr16":
            break
    t = inview_now(b)
    s, e = int(key.split("\t")[2]), int(key.split("\t")[3])
    if c["chrom"] != "chr16" or key not in t["drawn"] or "DEMOC" not in t["card"] or t["current"] != [key]:
        bad.append(f"the DEMOC row flew to {c['chrom']} {c['view']}, drawn {key in t['drawn']}, card {t['card'][:40]!r}, "
                   f"current {t['current']}")
    c, t = go_view(b, "#chr11")
    key = next((k for k in t["rows"] if k.startswith("roh\t")), None)
    if not key:
        return bad + [f"no run of homozygosity row at chr11: {t['rows']}"]
    b.eval("document.querySelector('.gm-table tr[data-key=\"%s\"] .gm-go').click(); true" % key.replace("\t", "\\t"))
    c, t = map_now(b), inview_now(b)
    s, e = int(key.split("\t")[2]), int(key.split("\t")[3])
    m = re.match(r"chr11:([0-9.]+)-([0-9.]+)M$", c["view"])
    if not m or not float(m.group(1)) * 1e6 <= s + 1e5 or not float(m.group(2)) * 1e6 >= e - 1e5 or t["card"] or \
            t["current"] != [key]:
        bad.append(f"the run's row {key!r} flew to {c['view']}, card {t['card'][:40]!r}, current {t['current']}")
    return bad


def strip_x(b, doc, chrom, frac, x0=128, pad=22):
    """Where FRAC of CHROM sits on the strip, in page px (map.js stripGeo() and sx(); X0 and PAD are a phone's
    64 and 10)."""
    r = json.loads(b.eval("JSON.stringify(document.querySelector('.gm-strip').getBoundingClientRect())"))
    rows = [c for c in doc["ref"]["chroms"] if not (c["n"] == "chrY" and doc["people"][0].get("has_y") is False)]
    total = sum(c["len"] for c in rows) + 5e6 * (len(rows) - 1)
    off = 0
    for c in rows:
        if c["n"] == chrom:
            return (r["left"] + x0 + (off + frac * c["len"]) * (r["width"] - pad - x0) / total,
                    r["top"] + r["height"] / 2 - (2 if r["height"] > 20 else 0), c["len"])
        off += c["len"] + 5e6
    raise ValueError(chrom)


def drag(b, x0, y, x1):
    b.call("Input.dispatchMouseEvent", {"type": "mouseMoved", "x": x0, "y": y}, b.session)
    b.call("Input.dispatchMouseEvent", {"type": "mousePressed", "x": x0, "y": y, "button": "left", "buttons": 1,
                                        "clickCount": 1}, b.session)
    for i in range(1, 9):
        b.call("Input.dispatchMouseEvent", {"type": "mouseMoved", "x": x0 + (x1 - x0) * i / 8, "y": y, "button": "left",
                                            "buttons": 1}, b.session)
        time.sleep(0.02)
    b.call("Input.dispatchMouseEvent", {"type": "mouseReleased", "x": x1, "y": y, "button": "left", "buttons": 0,
                                        "clickCount": 1}, b.session)


def hash_now(b, wait=0.35):
    """The address shortly after an action: a click or a brush on the strip pushes its place at once, a pan
    replaces it 250 ms after the last move."""
    time.sleep(wait)
    return b.eval("location.hash")


def view_mb(h):
    m = re.match(r"#chr[0-9XY]+:([0-9.]+)-([0-9.]+)M$", h or "")
    return (float(m.group(1)), float(m.group(2))) if m else None


def strip_problems(b, base, out, phone_only=False):
    """The strip: a click opens a chromosome whole, a drag elsewhere brushes a new view on the chromosome where it
    started, a drag on the brush pans; at 390 px it shows 12 px high and a click there does nothing (a picture
    only), and the table starts folded with its count in its summary."""
    bad = []
    doc = read_map(os.path.join(out, "runs", S, "map.json"))
    if not phone_only:
        map_page(b, base, S, "#chr11")
        b.eval("document.querySelector('.gm-strip').scrollIntoView({block: 'center'}); true")
        time.sleep(0.3)
        x, y, _ = strip_x(b, doc, "chr16", 0.5)
        drag(b, x, y, x)
        h = hash_now(b)
        if h != "#chr16":
            bad.append(f"a click on chr16 of the strip opens {h}")
        time.sleep(1.2)  # the map opens chr16 before the next drag
        x0, y, n = strip_x(b, doc, "chr2", 0.2)
        x1 = strip_x(b, doc, "chr2", 0.4)[0]
        drag(b, x0, y, x1)
        h = hash_now(b)
        v = view_mb(h)
        if not h.startswith("#chr2:") or not v or abs(v[0] - 0.2 * n / 1e6) > 4 or abs(v[1] - 0.4 * n / 1e6) > 4:
            bad.append(f"a brush over chr2 20-40% opens {h}")
        if v:
            time.sleep(1.2)  # the map opens the brushed view before the brush is dragged
            mid = strip_x(b, doc, "chr2", (v[0] + v[1]) / 2 * 1e6 / n)[0]
            drag(b, mid, y, mid + 12)
            p = view_mb(hash_now(b, 0.6))
            if not p or not (p[0] - v[0]) > 5:
                bad.append(f"a drag on the brush from {v} moved the view to {p}")
    b.call("Emulation.setDeviceMetricsOverride", {"width": 390, "height": 844, "deviceScaleFactor": 2, "mobile": False},
           b.session)
    try:
        map_page(b, base, S, "#chr11")
        t = inview_now(b)
        if t["stripH"] != 12 or t["open"] or not t["summary"].startswith("In view · ") or not t["summary"].endswith("features"):
            bad.append(f"at 390 px: strip {t['stripH']} px, table open {t['open']}, summary {t['summary']!r}")
        b.eval("document.querySelector('.gm-strip').scrollIntoView({block: 'center'}); true")
        time.sleep(0.3)
        x, y, _ = strip_x(b, doc, "chr16", 0.5, 64, 10)
        drag(b, x, y, x)
        h = hash_now(b)
        if h != "#chr11":
            bad.append(f"at 390 px a click on chr16 of the strip moved the map to {h}")
    finally:
        b.call("Emulation.setDeviceMetricsOverride", {"width": 1440, "height": 900, "deviceScaleFactor": 1, "mobile": False},
               b.session)
    return bad


WINDOW_JS = r"""(function () {
  var w = document.querySelector('.gm-inview-wrap'), seen = {}, steps = 0;
  return new Promise(function (done) {
    function look() {
      Array.prototype.forEach.call(w.querySelectorAll('tbody tr[data-key]'), function (r) { seen[r.dataset.key] = 1; });
      if (w.scrollTop + w.clientHeight >= w.scrollHeight - 1 || steps++ > 400) { done(JSON.stringify(Object.keys(seen))); return; }
      w.scrollTop += 240;
      window.requestAnimationFrame(function () { window.requestAnimationFrame(look); });
    }
    look();
  });
}())"""


def window_problems(b, base):
    """On the genome-scale document, chromosome 1 open: 12 rows and the words for the rest, then Show all lists
    every run in a scrolled window (fewer rows in the table at once than in view), and scrolling it through shows
    each mark the canvas drew, once."""
    bad = []
    map_page(b, base, S, "#chr1")
    first = json.loads(b.eval("JSON.stringify([document.querySelectorAll('.gm-table tbody tr[data-key]').length, "
                              "document.querySelector('.gm-inview-foot').textContent])"))
    t = inview_now(b)  # Show all
    n = len(t["drawn"])
    if n <= 200 or first[0] != 12 or not first[1].startswith(f"and {n - 12:,} more in view · sorted by where · findings first"):
        bad.append(f"before Show all: {first} for {n} marks")
    virtual = b.eval("document.querySelector('.gm-inview-wrap').classList.contains('gm-virtual')")
    if not virtual or not 12 < len(t["rows"]) < n or t["count"] != count_words(n):
        bad.append(f"Show all: window {virtual}, {len(t['rows'])} rows in the table of {n} ({t['count']})")
    r = b.call("Runtime.evaluate", {"expression": WINDOW_JS, "awaitPromise": True, "returnByValue": True}, b.session)
    seen = json.loads(r.get("result", {}).get("value") or "[]")
    if sorted(seen) != sorted(t["drawn"]):
        bad.append(f"scrolled through, the window showed {len(seen)} of the {n} marks drawn")
    return bad


def inview_checks(b, base, out, og):
    """The whole-genome strip, the In view table and drawsAt(), each check with its planted mistake."""
    mjs = os.path.join(out, "assets", "map.js")
    with open(mjs) as fh:
        keep_mjs = fh.read()

    def planted(desc, a, z, fn, more=()):
        """FN's problems with map.js changed from A to Z (and each pair of MORE): the check must fail there."""
        pairs = [(a, z)] + list(more)
        try:
            text = keep_mjs
            for x, y in pairs:
                text = text.replace(x, y)
            with open(mjs, "w") as fh:
                fh.write(text)
            control(desc, any(keep_mjs.count(x) != 1 for x, _ in pairs) or not fn())
        finally:
            with open(mjs, "w") as fh:
                fh.write(keep_mjs)

    keep = plant_long_run(out, og)
    try:
        bad = predicate_problems(b, base, out, og)
        check("map: drawsAt() in map.js gives og map build --predicate-table's answer for every row of the three demo "
              "runs and every gene, at every level", not bad, bad)
        planted("map.js whose drawsAt() drops the long runs at the genome level",
                "return i > 0 || !!(row && row.long);", "return i > 0;", lambda: predicate_problems(b, base, out, og))
    finally:
        unplant(keep)
    bad = inview_problems(b, base)
    check("map: on views at every level the In view table lists exactly what the canvas draws one by one, with its "
          "count and a line for the runs drawn as a band; the strip shows only over an open chromosome and ticks the "
          "findings the genome level draws; T folds and opens the table", not bad, bad)
    planted("map.js whose table lists every run of homozygosity, drawn one by one or not",
            "if (drawsAt('roh', { long: R.o.long[i] }, lv)) { rows.push(", "if (true) { rows.push(",
            lambda: inview_problems(b, base, INVIEW_VIEWS[:1]))
    planted("map.js whose strip ticks the findings of the open chromosome only",
            "        if (!drawsAt('findings', r, 'genome')) { return; }",
            "        if (!drawsAt('findings', r, 'genome') || r.chrom !== chrom().n) { return; }",
            lambda: inview_problems(b, base, INVIEW_VIEWS[:2]))
    bad = inview_veil_problems(b, base)
    check("map: a closed result is no row of the In view table and no tick of the strip (APOE on an adult's page, "
          "the son's result to act on), until opened", not bad, bad)
    planted("map.js whose table lists every finding, closed or not",
            "drawn().forEach(function (r) { if (meets(", "M.findings.forEach(function (r) { if (meets(",
            lambda: inview_veil_problems(b, base, adult=False))
    planted("map.js whose strip ticks every finding, closed or not",
            "      drawn().forEach(function (r) {\n        if (!drawsAt('findings', r, 'genome'))",
            "      M.findings.forEach(function (r) {\n        if (!drawsAt('findings', r, 'genome'))",
            lambda: inview_veil_problems(b, base, adult=False))
    bad = sort_problems(b, base)
    check("map: the In view table sorts by place with the findings first (aria-sort ascending on Where), Where again "
          "reverses it, and Size sorts by size", not bad, bad)
    planted("map.js whose table sorts by place alone",
            "return (a.layer === 'findings' ? 0 : 1) - (b.layer === 'findings' ? 0 : 1) || place(a, b);",
            "return place(a, b);", lambda: sort_problems(b, base))
    bad = announce_problems(b, base)
    check("map: the In view table announces its count once when it changes, never while it stays the same", not bad, bad)
    planted("map.js whose table announces its count on every settle",
            "if (n !== announced) { announced = n; inLive", "if (true) { announced = n; inLive",
            lambda: announce_problems(b, base))
    bad = fly_problems(b, base)
    check("map: a row of the In view table flies the map to its feature (a finding with its card) and becomes the "
          "current row", not bad, bad)
    planted("map.js whose rows fly without the finding's card",
            "if (x.layer === 'findings') { pick(x.r.id); } else { pick(null); }", "pick(null);",
            lambda: fly_problems(b, base))
    bad = strip_problems(b, base, out)
    check("map: the strip opens a chromosome on a click, brushes a new view and pans on a drag; at 390 px it is a "
          "12 px picture and the table starts folded", not bad, bad)
    planted("map.js whose strip ignores a click",
            "if (i >= 0) { openChrom(i, undefined, undefined, true); }\n      } else if (b0.mode",
            "void i;\n      } else if (b0.mode", lambda: strip_problems(b, base, out))
    planted("map.js whose strip takes clicks on a phone",
            "if (!M || phone || !b0) { return; }", "if (!M || !b0) { return; }",
            lambda: strip_problems(b, base, out, phone_only=True),
            more=[("if (!M || phone || e.button !== 0 || level() === 'genome') { return; }",
                   "if (!M || e.button !== 0 || level() === 'genome') { return; }")])


def map_screens(b, base, dest):
    """Screenshots of the map, dark, at 1440 x 900 and 390 x 844: genome, chromosome (with the whole-genome strip),
    region and gene, a closed page, and the In view table of a region (opened on a phone, where it starts folded)."""
    os.makedirs(dest, exist_ok=True)
    taken = []
    try:
        for w, h in ((1440, 900), (390, 844)):
            b.call("Emulation.setDeviceMetricsOverride", {"width": w, "height": h, "deviceScaleFactor": 1, "mobile": False},
                   b.session)
            for run, link, name, want in ((S, "", "genome", "genome"), (S, "#chr16", "chromosome", "chromosome"),
                                          (S, "#chr11:4.5-6.5M", "region", "region"),
                                          (S, "#chr11:5.150-5.350M", "gene", "gene"), (OTHERS[1], "", "closed", "genome"),
                                          (S, "#chr22:40.5-43.5M", "inview", "region")):
                got = map_page(b, base, run, link)
                if got["level"] in ("region", "gene"):
                    time.sleep(0.8)  # the chromosome's gene file
                    got = map_now(b)
                if got["state"] != "drawn" or got["level"] != want:
                    print(f"ERROR: {run}{link}: the map is {got['state']} at {got['level'] or 'no'} level, "
                          f"want drawn at {want}: no screenshot")
                    continue
                if name == "inview":
                    b.eval("var d = document.querySelector('.gm-inview'); d.open = true; d.scrollIntoView({block: 'start'}); window.scrollBy(0, -80); true")
                else:
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
    bad = closer_problems(b, base)
    check("map: the region, gene and reads levels open from their addresses with their titles, paths and genes; "
          "Esc steps out one level at a time; the Region and Gene stops land on their level; + and → zoom and pan",
          not bad, bad)
    lvl = "if (bpp > CHROM_BPP) { return 'chromosome'; }"
    try:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs.replace(lvl, "return 'chromosome';"))
        control("map.js whose levels stop at the chromosome", keep_mjs.count(lvl) != 1 or not closer_problems(b, base))
    finally:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs)
    bad = find_problems(b, base, out)
    check("map: / opens the search; a gene, a band and a finding are found and Enter flies there (a finding with "
          "its card); the son's closed result is never a line, by its gene, id or title", not bad, bad)
    drawn_ = "var fs = drawn();\n      fs.filter("
    try:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs.replace(drawn_, "var fs = M.findings;\n      fs.filter("))
        control("map.js whose search looks through every finding, closed or not", keep_mjs.count(drawn_) != 1 or
                not find_problems(b, base, out))
    finally:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs)
    bad = nearest_problems(b, base, out)
    check("map: Nearest names the closest finding (else the first on the next chromosome) and flies there with its "
          "card; Esc closes the card first; it never leads to the son's closed result", not bad, bad)
    near_ = "var ch = chrom(), v = view(), fs = drawn(),"
    try:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs.replace(near_, "var ch = chrom(), v = view(), fs = M.findings,"))
        control("map.js whose Nearest leads to every finding, closed or not", keep_mjs.count(near_) != 1 or
                not nearest_problems(b, base, out))
    finally:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs)
    bad = origin_problems(b, base, out)
    check("map: Paint the copies keeps Origin off with its reason, on a run in a family and on one in none", not bad, bad)
    off = "'aria-pressed': 'false', 'aria-disabled': 'true',"
    try:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs.replace(off, "'aria-pressed': 'false',"))
        control("map.js whose Origin button is not off", keep_mjs.count(off) != 1 or not origin_problems(b, base, out))
    finally:
        with open(mjs, "w") as fh:
            fh.write(keep_mjs)
    inview_checks(b, base, out, og)
    perf_checks(b, base, out)
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
    return 0 if len(taken) == 12 else 1


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
