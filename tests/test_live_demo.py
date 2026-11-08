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
     page too) and the person page has its bottom tab bar; on the son's page a
     two planted CPSR results to act on name their genes and variants nowhere
     (their card included) until "Show the result" is clicked, and each click
     opens one result, there and on the finding's own page.
Checks 1, 2, 3, 5, the widgets (person.js without its widget line), the family depth (family.json
with the son's stretches emptied), the phone width and the son's closed result (person.js whose cards ignore
the veil) also run against a planted mistake and must
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
        widget_checks(b, base, out)
        family_depth_checks(b, base, out)
        phone_checks(b, base)
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
           ("DEMO-003", "paralogs", "copies", "[data-widget=copies] li", 5),
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
            f.write(keep.replace("add(art, widget(c.widget));", ""))
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
    finally:
        with open(rj, "w") as f:
            f.write(keep)
        with open(js, "w") as f:
            f.write(keep_js)
        for d in made:
            shutil.rmtree(d, ignore_errors=True)


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
