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
     term of an open card is still found; DEMO-001's page
     mid-replay and at the end leads with the report (a short run band above
     the lead while it runs, none once finished), ends with the step timeline
     folded, keeps the wire to at most eight one-sentence lines with the plain
     finishes folded into one, and puts the HLA card's caveat in its margin
     with its absent and undetermined lines folded.
Checks 1, 2, 3, 5, the widgets (person.js without its widget line), the family depth (family.json
with the son's stretches emptied), the phone width and the son's closed result (person.js whose cards ignore
the veil, person.js and finding.js that read a missing state as a standalone adult, and search.js
without its filter), the phone page (person.js that never folds, and a route that opens only the nearest fold), the
folds (a flagged section 5 card folded, a redraw that forgets the folds, a Ctrl+click that opens one) and
the person page's layout (a wire showing every line or every plain finish, a card showing
every line, a step timeline left open) also run against a planted mistake and must
fail there (2: a build with one plugin's state hook dropped): a check that cannot fail is no check.

Run: python3 tests/test_live_demo.py
"""
import copy
import json
import os
import shutil
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
            f.write(keep.replace("add(text, widget(c.widget));", ""))
        control("person.js without its widget line", not widget_problems(widget_counts(b, base)))
    finally:
        with open(js, "w") as f:
            f.write(keep)


def family_depth(b, base):
    b.go("about:blank")
    b.wait("location.href === 'about:blank'")
    b.go(f"{base}/runs/{FAM}/")
    return b.wait("document.documentElement.dataset.pgpStatus === 'ready' && JSON.stringify(["
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
        built = build_live_demo.build(out, work=os.path.join(work, "root"))
        frames = built["frames"]
        print(f"built: {len(frames)} frames, hooks {built['hooks']}")

        # 2. every plugin's state hook
        hooks = registered_hooks()
        check("2. the build calls every registered state hook", sorted(built["hooks"]) == sorted(hooks), (built["hooks"], hooks))
        bad = missing_hooks(built, hooks)
        check(f"2. every state hook ({', '.join(hooks)}) has a value on the demo's pages", not bad and len(hooks) >= 5, bad)
        dropped = build_live_demo.build(os.path.join(work, "demo-dropped"), work=os.path.join(work, "root-dropped"),
                                        drop_hook=hooks[0])
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
    sys.exit(main())
