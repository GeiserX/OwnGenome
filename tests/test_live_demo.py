#!/usr/bin/env python3
"""Build the public live demo and check its reports, map, phone layout and report gate.

The gate property walks all report payloads, including unknown fields. Immutable act records,
family membership, transport routes and other people's records are its controls. It compares
whole documents, accessibility text, addresses, pixels, storage and resource requests.
"""
import base64
import copy
import gzip
import json
import os
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

    # A path answered after a wait, once (the gate's slow-reload check): {path: seconds}.
    DELAY = {}

    def do_GET(self):
        wait = Static.DELAY.pop(self.path.split("?")[0], 0)
        if wait:
            time.sleep(wait)
        super().do_GET()


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
        map_checks(b, base, out)
        # The family veil group runs from test_owngenome_live_family.py in a different Guard shard.
        gate_checks(b, base, out)
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
           ("DEMO-003", "clinvar", "AnnotSV span", "figure[data-chart=genome] rect.sv", 2)]


def open_gate(b):
    assert b.wait("document.documentElement.dataset.pgpStatus === 'ready'"), 'page did not settle: '+str(b.eval("JSON.stringify([location.pathname,document.documentElement.dataset,document.body.textContent.slice(0,200)])"))
    for _ in range(10):
        closed = b.eval("document.querySelector('[data-opening^=closed]') !== null")
        if not closed:
            return
        b.eval("document.querySelector('[data-opening^=closed]').click(); true")
        if not b.eval("document.querySelector('[data-opening^=closed]') !== null"):
            assert b.wait("!document.documentElement.dataset.pgpGated"), 'gate did not release'
            return
    raise AssertionError('gate did not release in ten clicks')


def widget_counts(b, base):
    got = {}
    for person, cid, kind, sel, _ in WIDGETS:
        b.go("about:blank")  # a new document each time, so the wait never reads the page before
        b.wait("location.href === 'about:blank'")
        b.go(f"{base}/runs/{person}/#s/{cid}")
        open_gate(b)
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
    assert b.wait("document.documentElement.dataset.pgpStatus === 'ready'"), 'family never ready'
    for _ in range(10):
        n = b.eval("document.querySelectorAll('[data-opening^=closed]').length")
        if not n:
            break
        b.eval("document.querySelector('[data-opening^=closed]').click(); true")
        assert b.wait("document.querySelectorAll('[data-opening^=closed]').length < %d" % n), 'family opening did not settle'
    assert b.wait("!document.querySelector('.report-gate')"), 'full family projection never settled'
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


def parity_problems(b, base, rows):
    """The karyogram's marks the map does not place with the same level and chromosome within 0.01 Mb."""
    bad = []
    for run in sorted({r[0] for r in rows if r[0] != OTHERS[1]}):
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
    """A browser with no DecompressionStream: the sentence and the findings with a place."""
    b.call("Page.enable", session=b.session)  # a script for each new document needs the Page domain on
    sid = b.call("Page.addScriptToEvaluateOnNewDocument", {"source": "delete window.DecompressionStream;"}, b.session)["identifier"]
    try:
        c = map_page(b, base, S)
        rows = b.eval("document.querySelectorAll('.gm-fallback table tbody tr').length")
    finally:
        b.call("Page.removeScriptToEvaluateOnNewDocument", {"identifier": sid}, b.session)
    bad = []
    if c["state"] != "nodecompress" or "This browser cannot open the map file" not in c["text"]:
        bad.append(f"no DecompressionStream: {c['state']}")
    if not rows:
        bad.append("no DecompressionStream: no findings table")
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
                if name == 'closed':
                    b.go('about:blank'); b.wait("location.href === 'about:blank'")
                    b.go(base+'/runs/'+run+'/')
                    if not b.wait("document.documentElement.dataset.pgpGated === 'true'"):
                        print(f"ERROR: {run}: no closed report gate: no screenshot")
                        continue
                    target = 'report-gate'
                else:
                    got = map_page(b, base, run, link)
                    want = "chromosome" if link else "genome"
                    if got["state"] != "drawn" or got["level"] != want:
                        print(f"ERROR: {run}{link}: the map is {got['state']} at {got['level'] or 'no'} level, "
                              f"want drawn at {want}: no screenshot")
                        continue
                    target = 'genome-map'
                b.eval("document.getElementById(%s).scrollIntoView({block: 'start'}); true" % json.dumps(target))
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
    mjs = os.path.join(out, 'assets', 'map.js')
    with open(mjs) as fh:
        keep_mjs = fh.read()
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


# Draws the same numbers in every load, so the order of the closed controls is the same in two loads that are
# compared; the order property (gate_order) runs without it.
SEEDED = "crypto.getRandomValues = function (a) { for (var i = 0; i < a.length; i++) { a[i] = %s; } return a; };"


def open_act(b, title, scope=''):
    """Open the act whose own record has TITLE. A closed control opens an act the reader cannot choose, so each is
    tried in turn and a wrong one is closed again. Returns the number of clicks."""
    closed, clicks = json.dumps(scope + '[data-opening^=closed]'), 0
    for k in range(10):
        n = b.eval("document.querySelectorAll(%s).length" % closed)
        assert n > k, 'no closed control opens ' + title
        m = b.eval("document.querySelectorAll(%s).length" % json.dumps(scope + '[data-opening^=opened]')) + 1
        b.eval("document.querySelectorAll(%s)[%d].click(); true" % (closed, k))
        clicks += 1
        if n == 1:
            return clicks
        got = json.dumps(scope + '[data-finding=opened-%d] .opened-act h2' % m)
        assert b.wait("!!document.querySelector(%s)" % got), 'an opening never showed'
        if b.eval("document.querySelector(%s).textContent" % got) == title:
            return clicks
        b.eval("document.querySelector(%s).click(); true" % json.dumps(scope + '[data-opening=opened-%d]' % m))
        clicks += 1
        assert b.wait("document.querySelectorAll(%s).length === %d" % (closed, n)), 'closing again never settled'
    raise AssertionError('no closed control opens ' + title)


# How far a run is, never what it found: the home and family pages show it for a gated child too (core.js PROGRESS).
PROGRESS = ['name', 'kind', 'url', 'state', 'started_utc', 'last_heard_utc', 'finished_utc', 'planned', 'done', 'running', 'failed']


def progress(run):
    return {k: run[k] for k in PROGRESS if k in run}


def resources(b, name):
    return b.eval("performance.getEntriesByType('resource').filter(function (r) { return r.name.endsWith(%s); }).length"
                  % json.dumps(name))


def gate_checks(b, base, out, family_only=False):
    """A single ownership property over every payload, including fields added after this test."""
    import itertools
    start = time.monotonic()
    son = OTHERS[1]
    helper = open(os.path.join(REPO, 'tests', 'demo', 'gate-property.js')).read()
    paths = [os.path.join(out, 'runs', son, n) for n in os.listdir(os.path.join(out, 'runs', son))
             if n.endswith('.json')]
    paths += [os.path.join(out, 'runs.json'), os.path.join(out, 'runs', FAM, 'family.json')]
    saved = {p: open(p, 'rb').read() for p in paths}
    result_path = os.path.join(out, 'runs', son, 'results.json')
    original = json.loads(saved[result_path])
    acts = [f for f in original['findings'] if f['level'] == 'act']
    # Placed, unplaced and multi-gene acts with stable identities and their own titles. Other records are not controls.
    acts = (acts + [dict(acts[0], id='gate-third', gene='HBA1', genes=['HBA1', 'HBA2'], viewer=None,
                         title='Gate fixture result in HBA1', detail='HBA1 and HBA2')])[:3]
    assert len({f['title'] for f in acts}) == 3, 'the acts need their own titles'

    def write(p, d):
        if p.endswith('/map.json'):
            write_map(p, d)
        else:
            with open(p, 'w') as f:
                json.dump(d, f)

    def walk(d):
        b.eval(helper)
        return b.eval('GateProperty.walk(' + json.dumps(d) + ', function () { return false; })')

    payload_cache = {}

    def payloads(count, changed=False, keep_first=False):
        cache_key = (count, changed, keep_first if changed else False)
        if cache_key in payload_cache:
            for p, raw in payload_cache[cache_key].items():
                with open(p, 'wb') as f:
                    f.write(raw)
            return
        for p, raw in saved.items():
            d = json.loads(gzip.decompress(raw)) if p.endswith('/map.json') else json.loads(raw)
            if p == result_path:
                d['findings'] = [f for f in d['findings'] if f['level'] != 'act'] + acts[:count]
                d.setdefault('veil', {})['acts'] = [f['id'] for f in acts[:count]]
                if changed:
                    nonacts = [f for f in original['findings'] if f['level'] != 'act']
                    d = walk(d)
                    changed_acts = d['findings'][len(nonacts):len(nonacts)+count]
                    d['findings'] = d['findings'][:len(nonacts)] + [d['findings'][-1]] + [
                        f if keep_first and i == 0 else dict(changed_acts[i], id=f['id'], level='act')
                        for i, f in enumerate(acts[:count])] + [None]
                    d['veil']['acts'] = [f['id'] for f in acts[:count]]
            elif changed and '/'+son+'/' in p:
                if p.endswith('/frames.json'):
                    # Family membership and URLs are routing controls, not report content.
                    family_plugin = json.loads(raw)['frames'][-1]['state']['plugins']['family']
                    d = walk(d)
                    for fr in d['frames']:
                        fr['state']['plugins']['family'] = family_plugin
                else:
                    d = walk(d)
            elif changed and p.endswith('/runs.json'):
                for i, run in enumerate(d['runs']):
                    if run['name'] == son:
                        d['runs'][i] = dict(walk(run), **dict(progress(run), name=son, url=run['url']))
                        # How many results to act on is what the closed page may say.
                        d['runs'][i]['levels']['act'] = run['levels']['act']
            elif changed and p.endswith('/family.json'):
                # One generic walk, then restore the other open owners and immutable routing controls.
                changed_family = walk(d)
                for j, child in enumerate(d['state']['children']):
                    if child['child'] == son:
                        d['state']['children'][j] = dict(changed_family['state']['children'][j], child=son)
                for key in ('relationships', 'ibd_pairs'):
                    for j, row in enumerate(d.get(key, [])):
                        if son in (row.get('sample_a'), row.get('sample_b')):
                            d[key][j] = dict(changed_family[key][j], sample_a=row['sample_a'], sample_b=row['sample_b'])
                for key in ('trio_files', 'trio_file_state'):
                    if key in d:
                        d[key] = changed_family[key]
                # Plan rows about the parents alone are theirs; a row naming the child is the child's.
                d['state']['plan'] = [dict(changed_family['state']['plan'][j], members=r['members'])
                                      if son in str(r.get('members', '')).split(',') else r
                                      for j, r in enumerate(d['state']['plan'])]
                d['state']['plan'].append(dict(analysis='future_child_analysis', members=son, decision='runs', reason='new'))
                d['state']['reports'] = changed_family['state']['reports']
                for key in ('segments', 'trio_status'):
                    if son in d.get(key, {}):
                        d[key][son] = changed_family[key][son]
                for key in ('de_novo', 'compound_het'):
                    d[key] = [dict(changed_family[key][i], child=son) if r.get('child') == son else r
                              for i, r in enumerate(d.get(key, []))]
                    d[key] += [dict(child=son, future_child_field='new child record')]
                sb = d['side_by_side']; i = sb['members'].index(son)
                for j, row in enumerate(sb['rows']):
                    changed_row = changed_family['side_by_side']['rows'][j]
                    row['values'][i] = changed_row['values'][i]
                    row['stale'][i] = changed_row['stale'][i]
                    row['future_child_note'] = changed_row['future_payload_field']
                sb['rows'].append(dict(group='Future child group', label='Future child field',
                                       values=['no summary yet' if m != son else 'new child value' for m in sb['members']],
                                       stale=[False]*len(sb['members'])))
                sb['rows'].append(copy.deepcopy(sb['rows'][0]))
                sb['rows'].reverse()
                if son in d.get('runs', {}):
                    d['runs'][son] = dict(changed_family['runs'][son], **dict(progress(d['runs'][son]), url=d['runs'][son]['url']))
                d['future_child_payload'] = changed_family['future_payload_field']
                d = {k: d[k] for k in reversed(d)}
            write(p, d)
        payload_cache[cache_key] = {p: open(p, 'rb').read() for p in saved}

    def visit(route):
        b.go('about:blank')
        assert b.wait("location.href === 'about:blank'"), 'navigation did not start'
        b.go(base + route)
        assert b.wait("document.readyState === 'complete' && document.documentElement.dataset.pgpStatus === 'ready'"), route + ' never ready'
        # The gate is synchronous once its inputs settle; no arbitrary multi-second settling sleep.
        if route.startswith('/runs/'+son+'/'):
            assert b.wait("document.documentElement.dataset.pgpGated === 'true'"), route + ' not gated'
        b.eval(helper)

    def snap():
        b.eval(helper)
        s = b.eval('GateProperty.snapshot()')
        ax = b.call('Accessibility.getFullAXTree', {}, b.session)
        s['accessibility'] = [(n.get('role', {}).get('value'), n.get('name', {}).get('value'))
                              for n in ax.get('nodes', []) if not n.get('ignored')]
        return s

    routes = ['/runs/'+son+'/', '/', '/runs/'+FAM+'/',
              '/runs/'+son+'/finding/'+acts[0]['id']+'/', '/runs/'+son+'/#chr7']
    if family_only:
        routes = ['/runs/'+FAM+'/']
    # Full Chrome asks for a fallback favicon after load, outside page readiness. Give the static
    # fixture a fixed inline icon. Every request the page makes still stays in the snapshot.
    html_paths = {os.path.join(out, route.split('#')[0].lstrip('/'), 'index.html') for route in routes}
    html_saved = {p: open(p, 'rb').read() for p in html_paths}
    for p, raw in html_saved.items():
        assert b'<head>' in raw, 'fixture has no head: '+p
        with open(p, 'wb') as f:
            f.write(raw.replace(b'<head>', b'<head><link rel="icon" href="data:,">', 1))
    b.call('Network.enable', {}, b.session)
    b.call('Network.setCacheDisabled', {'cacheDisabled': True}, b.session)
    b.call('Page.enable', session=b.session)  # scripts for each new document need the Page domain
    seed = b.call('Page.addScriptToEvaluateOnNewDocument', {'source': SEEDED % '0'}, b.session)['identifier']
    bad = []
    try:
        b.eval(helper)
        future = walk({'future_field': {'s': 'new', 'n': 2, 'b': True, 'a': []}})
        check('gate walker changes a new field, primitives, array length and key order',
              future != {'future_field': {'s': 'new', 'n': 2, 'b': True, 'a': []}}
              and future['future_field']['n'] != 2 and len(future['future_field']['a']) == 1, future)
        for width, media in ((1440, 'screen'), (390, 'screen'), (1440, 'print')):
            b.call('Emulation.setDeviceMetricsOverride', dict(width=width, height=900, deviceScaleFactor=1, mobile=False), b.session)
            b.call('Emulation.setEmulatedMedia', {'media': media, 'features': [{'name':'prefers-color-scheme', 'value':'dark'}]}, b.session)
            for route in routes:
                snapshots, partials = [], []
                for changed in (False, True):
                    payloads(3, changed); visit(route); snapshots.append(snap())
                    if route != '/' and '#' not in route:
                        payloads(3, changed, keep_first=True); visit(route)
                        scope = '[data-gate-child=%s] ' % son if route == '/runs/'+FAM+'/' else ''
                        completed = resources(b, '/'+son+'/results.json')
                        clicks = open_act(b, acts[0]['title'], scope)
                        if scope:
                            # Each click asks for the family again, and the gate for the child's results.
                            assert b.wait("performance.getEntriesByType('resource').filter(function(r){return r.name.endsWith('/"+son+"/results.json');}).length >= "+str(completed + clicks)), 'family projection refresh never settled'
                        partials.append(snap())
                for what, pair in (('closed', snapshots), ('partially opened', partials)):
                    if pair and pair[0] != pair[1]:
                        keys = [k for k in pair[0] if pair[0][k] != pair[1][k]]
                        x, y = pair[0]['dom'], pair[1]['dom']
                        at = next((i for i, (m, n) in enumerate(zip(x, y)) if m != n), min(len(x), len(y)))
                        bad.append((width, media, route, what, keys, x[max(0, at-60):at+60], y[max(0, at-60):at+60],
                                    [v['requests'] for v in pair] if 'requests' in keys else []))
        print('   gate matrix time: %.2f s' % (time.monotonic()-start))
        check('gate property: every mutable payload produces identical DOM, AX text, address, pixels, SVG, storage and requests', not bad, bad)
        if not family_only:
            b.call('Emulation.setEmulatedMedia', {'media': 'screen'}, b.session)
            b.call('Emulation.setDeviceMetricsOverride', dict(width=1440, height=900, deviceScaleFactor=1, mobile=False), b.session)
            transitions = []
            for count in (1, 2, 3):
                payloads(count)
                for order in itertools.permutations(range(count)):
                    visit('/runs/'+son+'/')
                    for stage, i in enumerate(order):
                        if stage < count - 1:
                            before = snap()
                            open_act(b, acts[i]['title'])
                            assert b.wait("document.documentElement.dataset.pgpGated === 'true'"), 'released early'
                            b.eval("document.querySelector('[data-opening=opened-%d]').click(); true" % (stage + 1))
                            if snap() != before:
                                transitions.append((count, order, stage, 'close-again differs'))
                            open_act(b, acts[i]['title'])
                        else:
                            b.eval("document.querySelector('[data-opening^=closed]').click(); true")
                            assert b.wait("!document.documentElement.dataset.pgpGated && !!document.querySelector('#results .finding')"), 'last opening never released'
                            # Search something while open; closing must replace this whole document and its caches.
                            b.wait("document.querySelector('input[type=search]') && !document.querySelector('input[type=search]').disabled")
                            b.eval("var s=document.querySelector('input[type=search]'); s.value='KCNH2'; s.dispatchEvent(new Event('input')); true")
                            b.eval("PGPUI.closeAct(%s); true" % json.dumps(acts[i]['id']))
                            assert b.wait("document.documentElement.dataset.pgpGated === 'true'"), 'close did not clean up'
                            if b.eval("!!document.querySelector('input[type=search],.search-results,.search-status,canvas,svg,img')"):
                                transitions.append((count, order, stage, 'old view survives'))
                            b.call('Page.reload', {}, b.session)
                            assert b.wait("document.querySelectorAll('[data-opening^=closed]').length === %d" % count), 'reload kept opening'
            print('   gate matrix and transitions time: %.2f s' % (time.monotonic()-start))
            check('gate transitions: 1/2/3 acts, every order, close-again at every stage, search then close and reload', not transitions, transitions)
            payloads(3)
            guesses = []
            for h in ('#f/'+acts[0]['id'], '#f/closed-99', '#f/opened-1', '#f/closed-1', '#s/exomiser', '?chr=17&start=43050000#f/closed-99'):
                visit('/runs/'+son+'/'+h)
                guesses.append(b.eval('location.hash'))
                assert b.eval('location.search') == '', 'a gated address kept a data query'
            check('gate handles: typed and guessed ids share the first opaque closed handle', len(set(guesses)) == 1 and guesses[0] == '#f/closed-1', guesses)
            visit('/runs/'+son+'/')
            b.eval("location.hash='f/closed-2'; true")
            assert b.wait("location.hash === '#f/closed-2'")
            b.eval("location.hash='f/closed-3'; true")
            assert b.wait("location.hash === '#f/closed-3'")
            b.eval("history.back(); true")
            assert b.wait("location.hash === '#f/closed-2'")
            b.eval("history.forward(); true")
            assert b.wait("location.hash === '#f/closed-3'")
            b.eval("document.querySelector('[data-opening=closed-1]').click(); true")
            check('gate handles follow their act after another act opens', b.eval("location.hash") == '#f/closed-2')
            gate_addresses(b, base, out, acts, visit, snap)
            gate_slow_release(b, base, acts, visit, payloads)
    finally:
        b.call('Page.removeScriptToEvaluateOnNewDocument', {'identifier': seed}, b.session)
        for p, raw in dict(saved, **html_saved).items():
            with open(p, 'wb') as f:
                f.write(raw)
        b.call('Emulation.clearDeviceMetricsOverride', {}, b.session)
        b.call('Emulation.setEmulatedMedia', {'media': ''}, b.session)
    if family_only:
        gate_family_independence(b, base, out)
    if not family_only:
        gate_order(b, base, out)
        gate_controls(b, base, out)
        gate_transparency(b, base, out)
    print('gate browser time: %.2f s' % (time.monotonic()-start))


def gate_addresses(b, base, out, acts, visit, snap):
    """A gated child's finding address answers the same page for any id, and keeps no id in the address; once the
    report opens, its addresses are the full report's (#f/<id>)."""
    son = OTHERS[1]
    src = os.path.join(out, 'runs', son, 'finding', acts[0]['id'], 'index.html')
    page = open(src).read()
    made, seen = [], []
    try:
        # The live server answers every well-formed id of a gated child with this page (bin/owngenome_live_results.py).
        for fid in ('two-hits-cftr', 'zz-not-a-finding'):
            d = os.path.join(out, 'runs', son, 'finding', fid)
            os.makedirs(d)
            made.append(d)
            with open(os.path.join(d, 'index.html'), 'w') as f:
                f.write(page.replace(acts[0]['id'], fid))
        for fid in (acts[0]['id'], 'two-hits-cftr', 'zz-not-a-finding'):
            visit('/runs/'+son+'/finding/'+fid+'/')
            s = snap()
            seen.append({k: s[k] for k in ('dom', 'text', 'title', 'address', 'accessibility')})
    finally:
        for d in made:
            shutil.rmtree(d)
    check('gate addresses: a gated finding page is the same for an act, a guessed act and nonsense, with no id in the address',
          all(x == seen[0] for x in seen) and acts[0]['id'] not in seen[0]['address'], [x['address'] for x in seen])
    visit('/runs/'+son+'/')
    open_gate(b)
    assert b.wait("!!document.querySelector('#results .finding')"), 'full report never loaded'
    fid = acts[0]['id']
    got = b.eval("[!!document.querySelector('a[href=\"#f/%s\"]'), !!document.querySelector('a[href^=\"#f/opened-\"],a[href^=\"#f/closed-\"]')]" % fid)
    b.eval("location.hash = 'f/%s'; true" % fid)
    kept = b.wait("location.hash === '#f/%s' && !!document.getElementById('f-%s')" % (fid, fid), timeout=5)
    check('gate addresses: the opened report links and keeps #f/<id>, as the full report does', got == [True, False] and kept, (got, kept))


def gate_slow_release(b, base, acts, visit, payloads):
    """Opening the last act releases the report however slowly the reloaded page's results arrive."""
    son = OTHERS[1]
    payloads(2)
    visit('/runs/'+son+'/')
    open_act(b, acts[0]['title'])
    Static.DELAY['/runs/'+son+'/results.json'] = 6
    try:
        b.eval("document.querySelector('[data-opening^=closed]').click(); true")
        released = b.wait("!document.documentElement.dataset.pgpGated && !!document.querySelector('#results .finding')", timeout=25)
    finally:
        Static.DELAY.clear()
    check('gate release: a reload whose results.json takes 6 s still releases the report', released, released)


def gate_order(b, base, out):
    """Only an unopened act changes, its rule and gene and so its place in the producer's order; the opened act is the
    same. The pages must not differ, and which act a closed control opens comes from a draw of this page load."""
    son = OTHERS[1]
    rj = os.path.join(out, 'runs', son, 'results.json')
    raw = open(rj).read()
    data = json.loads(raw)
    acts = [f for f in data['findings'] if f['level'] == 'act']
    others = [f for f in data['findings'] if f['level'] != 'act']
    opened, closed = acts[0], acts[1]
    tp53 = dict(closed, id='cpsr-tp53-c-743g-a', gene='TP53', genes=['TP53'], title='Cancer predisposition (TP53)',
                detail='TP53 c.743G>A', viewer='chr17:7674220')
    # The producer lists findings by rule and gene: the KCNH2 act before a later act, a TP53 CPSR act before it.
    fixtures = {'later': [opened, closed], 'earlier': [tp53, opened]}
    pages = {}
    try:
        for name, order in fixtures.items():
            with open(rj, 'w') as f:
                json.dump(dict(data, findings=order + others, veil=dict(data.get('veil') or {}, acts=[x['id'] for x in order])), f)
            for route, scope in (('/runs/'+son+'/', ''), ('/runs/'+FAM+'/', '[data-gate-child=%s] ' % son)):
                b.go('about:blank'); assert b.wait("location.href === 'about:blank'")
                b.go(base + route)
                assert b.wait("!!document.querySelector(%s)" % json.dumps(scope + '[data-opening^=closed]')), route + ' not gated'
                completed = resources(b, '/'+son+'/results.json')
                clicks = open_act(b, opened['title'], scope)
                assert b.wait("!!document.querySelector(%s)" % json.dumps(scope + '.opened-act')), 'opening never settled'
                if scope:
                    assert b.wait("performance.getEntriesByType('resource').filter(function(r){return r.name.endsWith('/"+son+"/results.json');}).length >= "+str(completed + clicks)), 'family refresh never settled'
                time.sleep(0.3)
                pages[(name, route)] = b.eval("[document.documentElement.outerHTML, document.body.innerText, location.href, document.title]")
        bad = [r for r in ('/runs/'+son+'/', '/runs/'+FAM+'/') if pages[('later', r)] != pages[('earlier', r)]]
        check('gate order: an unopened act earlier or later in the report leaves the page with an opened act unchanged', not bad, bad)

        def first_opens(draw):
            seed = b.call('Page.addScriptToEvaluateOnNewDocument', {'source': SEEDED % draw}, b.session)['identifier']
            try:
                b.go('about:blank'); assert b.wait("location.href === 'about:blank'")
                b.go(base + '/runs/'+son+'/')
                assert b.wait("!!document.querySelector('[data-opening=closed-1]')"), 'not gated'
                b.eval("document.querySelector('[data-opening=closed-1]').click(); true")
                assert b.wait("!!document.querySelector('[data-finding=opened-1] .opened-act h2')")
                return b.eval("document.querySelector('[data-finding=opened-1] .opened-act h2').textContent")
            finally:
                b.call('Page.removeScriptToEvaluateOnNewDocument', {'identifier': seed}, b.session)

        def draws():
            return {first_opens('0'), first_opens('4294967295')}
        got = draws()
        check('gate order: which act the first closed control opens comes from the page load, not the report', len(got) == 2, got)
        js = os.path.join(out, 'assets', 'core.js')
        source = open(js).read()
        anchor = ': shuffled(ids);'
        assert source.count(anchor) == 1, 'mutation anchor moved: shuffled'
        try:
            with open(js, 'w') as f:
                f.write(source.replace(anchor, ': ids.slice();'))
            control('gate order: closed controls in the producer order', len(draws()) == 2)
        finally:
            with open(js, 'w') as f:
                f.write(source)
    finally:
        with open(rj, 'w') as f:
            f.write(raw)


def gate_controls(b, base, out):
    """Each gate guard must fail when its production implementation is deliberately broken."""
    son = OTHERS[1]
    js = os.path.join(out, 'assets', 'core.js')
    rj = os.path.join(out, 'runs', son, 'results.json')
    source = open(js).read(); raw = open(rj).read()
    data = json.loads(raw); data['future_payload_field'] = 'FUTURE_GATE_SENTINEL'
    with open(rj, 'w') as f:
        json.dump(data, f)
    seen = {}

    def visit(route=None):
        b.go('about:blank'); assert b.wait("location.href === 'about:blank'")
        b.go(base+(route or '/runs/'+son+'/'))
        if not route:
            assert b.wait("document.documentElement.dataset.pgpGated === 'true'"), 'fixture never gated'

    def partial():
        visit()
        b.eval("document.querySelector('[data-opening^=closed]').click(); true")
        ok = bool(b.wait("document.documentElement.dataset.pgpGated === 'true' && document.querySelectorAll('[data-opening^=closed]').length === 1", timeout=2))
        seen['released'] = bool(b.wait("!!document.documentElement.dataset.pgpResults && !document.documentElement.dataset.pgpGated", timeout=10)) if not ok else False
        return ok

    def raw_view():
        visit()
        # Both requests have settled when the optional deliberately broken view signals it is ready.
        b.wait("document.readyState === 'complete'", timeout=2)
        if 'data-gate-raw-read' in open(js).read():
            assert b.wait("document.documentElement.dataset.gateRawRead === 'done'", timeout=2), 'raw probe did not settle'
        return not b.eval("document.documentElement.outerHTML.includes('FUTURE_GATE_SENTINEL')")

    def handoff():
        # A view that asks a closed run for its files gets no answer at all, whatever it reads the body with.
        visit()
        b.eval("window.gateRaw = null; fetch('results.json').then(function (r) { return r.text(); }).then(function (t) { window.gateRaw = t.length; }); true")
        got = b.wait("window.gateRaw", timeout=2)
        seen['handed'] = bool(got)
        return not got

    def family_rows():
        # The family view never receives a gated child's trio rows, even ones it would not draw.
        visit('/runs/'+FAM+'/')
        assert b.wait("!!document.querySelector('[data-gate-child=%s]')" % son), 'family never gated'
        b.eval("window.gateRows = null; PGPUI.getJSON('family.json').then(function (d) { window.gateRows = "
               "(d.de_novo || []).concat(d.compound_het || []).filter(function (r) { return r.child === %s; }).length; }); true"
               % json.dumps(son))
        assert b.wait("window.gateRows !== null"), 'family read never settled'
        n = b.eval("window.gateRows")
        seen['rows'] = n
        return n == 0

    def cleanup():
        visit(); open_gate(b)
        assert b.wait("!!document.querySelector('#results .finding')"), 'full report never loaded'
        assert b.wait("document.querySelector('input[type=search]') && !document.querySelector('input[type=search]').disabled"), 'search never loaded'
        b.eval("var s=document.querySelector('input[type=search]'); s.value='KCNH2'; s.dispatchEvent(new Event('input')); true")
        fid = next(f['id'] for f in data['findings'] if f['level']=='act')
        # Read in the same turn as the close: the document replacement that follows would hide a kept query.
        got = b.eval('PGPUI.closeAct('+json.dumps(fid)+"); var s = document.querySelector('input[type=search]');"
                     "[PGPUI.gateOwnsPage() && !document.querySelector('input[type=search],.search-results,.search-status,canvas,svg,img'), s ? s.value : '']")
        seen['query'] = got[1]
        return got[0]

    def address():
        b.go('about:blank'); assert b.wait("location.href === 'about:blank'")
        b.go(base+'/runs/'+son+'/?chr=17&start=43050000#f/closed-99')
        assert b.wait("document.documentElement.dataset.pgpGated === 'true'")
        return b.eval("location.search === '' && location.hash === '#f/closed-1'")

    def stale_control():
        visit('/runs/'+FAM+'/')
        assert b.wait("!!document.querySelector('[data-gate-child=%s] [data-opening^=closed]')" % son)
        b.eval("window.gateOldControl=document.querySelector('[data-gate-child=%s] [data-opening^=closed]'); true" % son)
        changed = copy.deepcopy(data)
        for f in changed['findings']:
            if f['level'] == 'act':
                f['text'] += ' REPLACEMENT_ACT_RECORD'
        try:
            with open(rj, 'w') as f:
                json.dump(changed, f)
            b.eval("window.gateRefreshDone=false; PGPUI.getJSON('family.json').then(function(){window.gateRefreshDone=true;}); true")
            assert b.wait('window.gateRefreshDone'), 'family refresh never settled'
            b.eval("gateOldControl.click(); true")
            # Read the outcome while the replacement is still what the server answers: the redraw the click starts
            # reads results.json again.
            opened = b.wait("!!document.querySelector('[data-gate-child=%s] .opened-act')" % son, timeout=1.5)
            seen['replacement'] = bool(opened) and b.eval(
                "document.querySelector('[data-gate-child=%s] .opened-act').textContent.includes('REPLACEMENT_ACT_RECORD')" % son)
            return not opened
        finally:
            with open(rj, 'w') as f:
                json.dump(data, f)

    def refresh_cleanup():
        visit(); open_gate(b)
        assert b.wait("!!document.querySelector('#results .finding')")
        # Hold document replacement so the guard must remove the old report synchronously.
        b.eval("window.gateTransferStarted=false; crypto.subtle.digest=function(){window.gateTransferStarted=true; return new Promise(function(){});}; true")
        changed=copy.deepcopy(data)
        next(f for f in changed['findings'] if f['level']=='act')['text'] += ' REFRESHED_ACT_RECORD'
        try:
            with open(rj, 'w') as f:
                json.dump(changed, f)
            b.eval("PGPUI.getJSON('results.json'); true")
            assert b.wait('window.gateTransferStarted'), 'refresh never reached document replacement'
            ok = b.eval("PGPUI.gateOwnsPage() && !document.querySelector('#results,canvas')")
            seen['report'] = b.eval("!!document.querySelector('#results .finding')")
            return ok
        finally:
            with open(rj, 'w') as f:
                json.dump(data, f)

    gate_open = "if (ownsPage || closedRun(requestURL)) { return new Promise(function () {}); }"
    gate_late = "return got.then(function (response) {\n        if (ownsPage) { return new Promise(function () {}); }"
    mutations = [
        ('stale control opens a replacement act', ['if (!same(byId(active, id), byId(rec, id))) { changed(); return; }'], [''], stale_control,
         lambda: seen.get('replacement')),
        ('refresh leaves the old report visible', ['if (shut(current)) { ownsPage = true; eventSources.forEach(function (es) { es.close(); }); cleanPage(); transfer();'],
         ['if (shut(current)) { transfer();'], refresh_cleanup, lambda: seen.get('report')),
        ('keep a map query in the closed address', ['location.pathname + to'], ['location.pathname + location.search + to'], address, lambda: True),
        ('release one act early', ["rec.acts.some(function (f) { return !rec.opened[f.id]; })"],
         ["rec.acts.filter(function (f) { return !rec.opened[f.id]; }).length > 1"], partial, lambda: seen.get('released')),
        ('one view reads the raw payload', ['add(document.body, wrap);'],
         ["add(document.body, wrap); jsonAt(current.url + 'results.json').then(function (r) { add(wrap, el('p', {}, JSON.stringify(r))); document.documentElement.setAttribute('data-gate-raw-read', 'done'); });"],
         raw_view, lambda: True),
        ('a view receives a closed run file', [gate_open, gate_late], ['', 'return got.then(function (response) {'], handoff,
         lambda: seen.get('handed')),
        ('the family view receives a gated child\'s trio rows', ['copy((raw[key] || []).filter(function (r) { return !gated(r.child); }))'],
         ['copy(raw[key] || [])'], family_rows, lambda: seen.get('rows')),
        ('skip cleanup after closing', ['ownsPage = true; eventSources.forEach(function (es) { es.close(); }); cleanPage(); transfer(); return true;'],
         ['return true;'], cleanup, lambda: seen.get('query') == 'KCNH2')]
    try:
        for desc, befores, afters, probe, broke in mutations:
            for before in befores:
                assert source.count(before) == 1, 'mutation anchor moved: '+desc
            seen.clear()
            check('gate guard: '+desc, probe())
            planted = source
            for before, after in zip(befores, afters):
                planted = planted.replace(before, after)
            with open(js, 'w') as f:
                f.write(planted)
            seen.clear()
            still_passes = probe()
            control('gate: '+desc, still_passes)
            assert still_passes or broke(), 'control failed without the broken outcome: '+desc
            print('   observed guard failure: '+desc)
            with open(js, 'w') as f:
                f.write(source)
    finally:
        with open(js, 'w') as f:
            f.write(source)
        with open(rj, 'w') as f:
            f.write(raw)


# A summary that gives an adult a result to act on (BRCA1) and two symptom rankings, one of the same gene: the
# home page's top entries and the person page come from a fresh producer run, not from the demo's files.
SYNTH = {"sections": {
    "cpsr": {"state": "ok", "data": {"plp": [{"gene": "BRCA1", "variant": "c.68_69del", "final": "Pathogenic",
                                               "clinvar": "Pathogenic"}]}},
    "exomiser": {"state": "ok", "data": {"genes": [{"gene": "BRCA1", "rank": 1, "exomiser_score": 0.9, "moi": "AD"},
                                                     {"gene": "TTN", "rank": 2, "exomiser_score": 0.8, "moi": "AR"}]}},
    "clinvar": {"state": "ok", "data": {"hits": []}},
    "coverage": {"state": "ok", "data": {"mean_depth": 31.2}}}}


def gate_transparency(b, base, out):
    """The gate leaves adult and no-act pages equal to the renderer without it: DOM, attributes, address, canvas and
    the requests, in order and counted."""
    son = OTHERS[1]
    js = os.path.join(out, 'assets', 'core.js')
    files = [os.path.join(out, 'runs', S, 'results.json'), os.path.join(out, 'runs', son, 'results.json'),
             os.path.join(out, 'runs.json')]
    source = open(js).read()
    saved = {p: open(p).read() for p in files}
    # The same current renderer with the loader absent is the transparent reference.
    marker = '  // --- the report gate'
    assert source.count(marker) == 1, 'gate marker moved'
    bare = source[:source.index(marker)] + '''
  PGPUI.ready = Promise.resolve();
  PGPUI.gateOwnsPage = function () { return false; };
  PGPUI.openedActs = function () { return {}; };
  PGPUI.closeAct = function () { return false; };
}());
'''

    def capture(route, apoe=False, timeout=20):
        b.go('about:blank'); assert b.wait("location.href === 'about:blank'")
        b.go(base + route)
        if route == '/':
            assert b.wait("document.documentElement.dataset.pgpStatus === 'ready' && !!document.documentElement.dataset.pgpPeople", timeout=timeout), 'home never drawn'
        elif '/finding/' in route:
            assert b.wait("document.documentElement.dataset.pgpStatus === 'ready' && !!document.querySelector('h1')", timeout=timeout), 'finding never drawn'
        else:
            assert b.wait("document.documentElement.dataset.pgpResults && !!document.querySelector('[data-pgp-map=drawn]')", timeout=timeout), 'open report not drawn'
            assert b.wait("!!document.querySelector('#rail a[href*=DEMO-FAMILY]') && !!document.documentElement.dataset.pgpCharts"), 'family or charts never settled'
        if apoe:
            # A redraw while the address names a place: the contents rail must come back as it did.
            b.eval("document.querySelector('details[data-gate=apoe] summary').click(); true")
        time.sleep(0.5)
        got = b.eval("""(function () {
          var clone=document.documentElement.cloneNode(true);
          clone.querySelectorAll('[data-pgp-map-paint-ms],[data-pgp-map-parse-ms]').forEach(function (n) {
            n.removeAttribute('data-pgp-map-paint-ms'); n.removeAttribute('data-pgp-map-parse-ms');
          });
          return [clone.outerHTML, location.hash, Array.from(document.querySelectorAll('canvas')).map(function(c){return c.toDataURL();}),
            // The browser lists a request when it ends: sort by when each one started.
            performance.getEntriesByType('resource').slice().sort(function (a, z) { return a.startTime - z.startTime; })
              .map(function (r) { return new URL(r.name).pathname; })];
        }())""")
        # Scripts and the views' own files load in parallel, in any order, with or without the gate: count every
        # request, and keep the order of the files the gate reads first.
        got[3] = [sorted(got[3]), [p for p in got[3] if p.endswith(('/state.json', '/frames.json', '/family.json', '/results.json'))]]
        return got

    def write(p, d):
        with open(p, 'w') as f:
            json.dump(d, f)

    def no_act():
        d = json.loads(saved[files[1]]); d['findings'] = [f for f in d['findings'] if f['level'] != 'act']
        d['levels']['act'] = 0
        write(files[1], d)
        idx = json.loads(saved[files[2]])
        for r in idx['runs']:
            if r['name'] == son:
                r['levels']['act'] = 0
                r['top'] = [t for t in r['top'] if t['level'] != 'act']
                r['lead'] = r['lead_closed'] = 'No finding in what these steps assessed needs a doctor.'
        write(files[2], idx)

    def synthetic():
        no_act()
        write(files[0], live_results.build(S, SYNTH))
        o = live_results.overview(SYNTH)
        idx = json.load(open(files[2]))
        for r in idx['runs']:
            if r['name'] == S:
                r.update(attention=o['attention'], levels=o['levels'], top=o['top'], lead=o['lead'],
                         lead_closed=o['lead_closed'], screened=o['screened'])
        write(files[2], idx)

    first = next(f['id'] for f in json.loads(saved[files[0]])['findings'] if f['level'] != 'act')
    cases = [('adult', None, '/runs/'+S+'/', False), ('adult, a redraw at a section address', None, '/runs/'+S+'/#s3', True),
             ('adult finding page', None, '/runs/'+S+'/finding/'+first+'/', False),
             ('child with no act', no_act, '/runs/'+son+'/', False), ('home with no gated child', no_act, '/', False),
             ('fresh adult with an act and a same-gene ranking', synthetic, '/runs/'+S+'/', False)]
    try:
        for label, prep, route, apoe in cases:
            for p in files:
                with open(p, 'w') as f:
                    f.write(saved[p])
            if prep:
                prep()
            pair = []
            for text in (bare, source):
                with open(js, 'w') as f:
                    f.write(text)
                pair.append(capture(route, apoe))
            detail = [('dom', 'hash', 'canvas', 'requests')[i] for i in range(4) if pair[0][i] != pair[1][i]]
            if pair[0][3] != pair[1][3]:
                detail += [pair[0][3], pair[1][3]]
            if pair[0][0] != pair[1][0]:
                at = next((i for i, (a, z) in enumerate(zip(pair[0][0], pair[1][0])) if a != z), 0)
                detail += [pair[0][0][max(0, at-80):at+120], pair[1][0][max(0, at-80):at+120]]
            check('gate transparency: '+label+' has identical DOM, attributes, address, canvas and requests', pair[0] == pair[1], detail)
        for p in files:
            with open(p, 'w') as f:
                f.write(saved[p])
        no_act()
        reference = capture('/runs/'+son+'/')
        for desc, before, after in (('closing an adult or no-act report', 'ownsPage = shut(current);', 'ownsPage = true;'),
                                    ('reading a file twice', 'var kept = shared[requestURL];', 'var kept = null;')):
            assert source.count(before) == 1, 'mutation anchor moved: '+desc
            with open(js, 'w') as f:
                f.write(source.replace(before, after))
            try:
                still_equal = capture('/runs/'+son+'/', timeout=4) == reference
            except AssertionError:
                still_equal = False
            control('gate transparency: '+desc, still_equal)
    finally:
        with open(js, 'w') as f:
            f.write(source)
        for p in files:
            with open(p, 'w') as f:
                f.write(saved[p])


def gate_family_independence(b, base, out):
    """Two children may have the same act IDs; opening one must never open the other."""
    son, other = OTHERS[1], 'DEMO-004'
    fj = os.path.join(out, 'runs', FAM, 'family.json')
    raw = open(fj).read(); family=json.loads(raw)
    member = next(m for m in family['state']['members'] if m['sample']==son)
    family['state']['members'].append(dict(member, sample=other))
    child = next(c for c in family['state']['children'] if c['child']==son)
    family['state']['children'].append(dict(child, child=other))
    family['runs'][other]=dict(family['runs'][son], url='../../runs/'+other+'/')
    sb=family['side_by_side']; i=sb['members'].index(son); sb['members'].append(other)
    for row in sb['rows']:
        row['values'].append(row['values'][i]); row['stale'].append(row['stale'][i])
    target=os.path.join(out,'runs',other)
    assert not os.path.exists(target), 'independent-child fixture already exists'
    os.makedirs(target)
    results=json.load(open(os.path.join(out,'runs',son,'results.json')))
    results['sample']=other; results['person']['sample']=other
    with open(os.path.join(target,'results.json'),'w') as f:
        json.dump(results,f)
    try:
        with open(fj,'w') as f:
            json.dump(family,f)
        b.go('about:blank'); assert b.wait("location.href === 'about:blank'")
        b.go(base+'/runs/'+FAM+'/')
        assert b.wait("!!document.querySelector('[data-gate-child=DEMO-004]')"), 'second child never gated'
        before=b.eval("document.querySelector('[data-gate-child=DEMO-004]').outerHTML")
        for remaining in range(sum(f['level']=='act' for f in results['findings']),0,-1):
            b.eval("document.querySelector('[data-gate-child=DEMO-003] [data-opening^=closed]').click(); true")
            assert b.wait("document.querySelectorAll('[data-gate-child=DEMO-003] [data-opening^=closed]').length === %d" % (remaining-1)), 'opening did not settle'
        after=b.eval("document.querySelector('[data-gate-child=DEMO-004]').outerHTML")
        check('family gate openings are independent for children sharing act IDs', before==after)
    finally:
        with open(fj,'w') as f:
            f.write(raw)
        shutil.rmtree(target)


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
