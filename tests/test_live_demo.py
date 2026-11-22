#!/usr/bin/env python3
"""Build the public live demo and check its reports, map, phone layout and report gate.

The gate property walks all report payloads, including unknown fields. Immutable act records,
family membership, transport routes and other people's records are its controls. It compares
whole documents, accessibility text, addresses, pixels, storage and resource requests.

The genome map of each person page (map_checks) is checked the same way, each check with its planted mistake: its
levels from the genome down to the reads and the addresses that open them, Esc and the zoom keys, the search box,
Nearest and the finding's card, Paint the copies, the veil and opening on purpose, and the frame budget on a
genome-scale document (perf_checks, with a 20 ms busy loop in draw() as its control).

Run: python3 tests/test_live_demo.py
"""
import base64
import copy
import gzip
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
sys.path.insert(0, os.path.join(REPO, "tests", "browser"))
import build_live_demo  # noqa: E402
import live_results  # noqa: E402
import make_demo_events  # noqa: E402
import owngenome_live  # noqa: E402
import run_state  # noqa: E402
from cdp import Chrome, asset_problems, find_chrome  # noqa: E402

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
        if self.path.split("?")[0] not in Static.BFCACHE:
            self.send_header("Cache-Control", "no-store")
        super().end_headers()

    # Paths answered without no-store, so Chrome may keep their page in the back-forward cache (gate_controls).
    BFCACHE = set()

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
        # The last opening loads the page again: a read that meets that navigation means the gate released.
        try:
            left = b.eval("document.querySelector('[data-opening^=closed]') !== null")
        except RuntimeError as e:
            if 'navigated' not in str(e):
                raise
            left = False
        if not left:
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
    finding with its card); APOE, closed on every page until opened, is never a line, by its gene, its id or its
    title, while an open finding still is. (The son's page is closed whole by the report gate: it has no map.)"""
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
    doc = read_map(os.path.join(out, "runs", S, "map.json"))
    f = doc["per"][S]["findings"]
    i = [j for j in range(len(f["id"])) if f["gate"][j]][0]
    fid, gene, title = f["id"][i], f["gene"][i], f["title"][i]
    map_page(b, base, S)
    for q in (gene, fid, title):
        found = map_find(b, q)
        # The gene itself stays findable (a line of kind "gene", from the gene names); its finding never is.
        named = [x for x in found if x["kind"] == "finding" and (gene in x["text"] or fid in x["text"] or title in x["text"])]
        if named:
            bad.append(f"{gene}, closed, is a line for {q!r}: {named}")
    found = map_find(b, "DEMOC")
    if not any(x["kind"] == "finding" for x in found):
        bad.append(f"the open finding DEMOC is not found: {found}")
    return bad


def nearest_problems(b, base, out):
    """Nearest: with no finding of the open chromosome in view it names the closest one (else the first on the
    next chromosome that has one) and flies there with its card; Esc closes the card; it never leads to APOE,
    closed on every page until opened."""
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
    doc = read_map(os.path.join(out, "runs", S, "map.json"))
    f = doc["per"][S]["findings"]
    i = [j for j in range(len(f["id"])) if f["gate"][j]][0]
    others = [f["s"][j] for j in range(len(f["id"])) if f["chrom"][j] == f["chrom"][i] and j != i]
    # A view next to the closed result, with every other finding of that chromosome farther away.
    a = max(0, f["s"][i] - 3_000_000)
    if not any(abs(x - f["s"][i]) < 6_000_000 for x in others):
        n = map_page(b, base, S, f"#{f['chrom'][i]}:{a + 1}-{a + 2_000_000}") and json.loads(b.eval(NEAR_JS))
        if f["gene"][i] in n["near"] or f["gene"][i] in n["card"]:
            bad.append(f"Nearest leads to {f['gene'][i]}, closed: {n['near']!r}")
    else:
        bad.append(f"{f['gene'][i]}, closed, has another finding within 6 Mb: no view to check Nearest on")
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


def unplant(keep):
    for p, data in keep.items():
        with open(p, "wb") as fh:
            fh.write(data)


def predicate_problems(b, base, out, og, runs=None):
    """drawsAt() of map.js against og map build --predicate-table, row by row, on the three demo runs (DEMO-001
    with a planted long run; RUNS only when given): the same rows and the same answer at each level; and the
    comparison covers a finding, a long run, a shorter run and a gene, so it cannot pass on rows that agree by
    chance."""
    bad, runs = [], list(runs or [S] + OTHERS)
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
                (OTHERS[0], "", "genome"), (OTHERS[0], "#chr16", "chromosome"))


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


def inview_veil_problems(b, base):
    """A closed result is never a row of the table or a tick of the strip: APOE on DEMO-001's page. The son's page,
    closed whole by the report gate while a result to act on is unopened, has no map, table or strip at all."""
    bad = inview_apoe_problems(b, base)
    son = OTHERS[1]
    b.go("about:blank")
    b.wait("location.href === 'about:blank'")
    b.go(f"{base}/runs/{son}/#chr7")
    if not b.wait("document.documentElement.dataset.pgpGated === 'true'"):
        bad.append(f"{son}: his page is not closed by the report gate")
    elif b.eval("!!document.querySelector('#genome-map, .gm-inview, .gm-strip, canvas')"):
        bad.append(f"{son}: his closed page has a map, a table or a strip")
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


def sort_problems(b, base, default_only=False):
    """The table's order: by place with the findings first, aria-sort ascending on Where; Where again reverses it
    (descending); Size sorts by size, smallest first, with aria-sort moved there. DEFAULT_ONLY: the first order
    alone."""
    bad = []
    c, t = inview_page(b, base, S, "#chr22:40.5-43.5M")
    keys = t["rows"]
    place = sorted(keys, key=lambda k: (int(k.split("\t")[2]), int(k.split("\t")[3]), k))
    first = sorted([k for k in place if k.startswith("findings\t")], key=place.index)
    if keys != first + [k for k in place if k not in first] or place[0] in first:
        bad.append(f"default order {keys[:4]}: want the findings first, then by place (and a gene before the finding)")
    if t["sort"] != ["", "", "ascending", "", "", "", ""]:
        bad.append(f"default aria-sort {t['sort']}")
    if default_only:
        return bad
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


def announce_problems(b, base, pans_only=False):
    """The table announces its count once when it changes and never when it does not: three pans along a stretch
    with nothing in view say nothing; then a view with four marks says "4 features in view", once. PANS_ONLY: the
    pans alone."""
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
    if pans_only:
        return bad
    c, t = go_view(b, "#chr11")
    said = json.loads(b.eval("JSON.stringify(window.__said)"))
    if said != [f"{t['count']} in view"] or t["count"] != "4 features":
        bad.append(f"#chr11 ({t['count']}) announced {said}")
    return bad


SHOWALL_JS = r"""(function () {
  var f = document.getElementById('genome-map'), a = f.querySelector('.gm-showall'), out = [];
  if (!a) { return JSON.stringify(null); }
  for (var i = 0; i < 2; i++) {
    a = f.querySelector('.gm-showall');
    a.focus();
    var was = a.textContent;
    a.click();
    var now = document.activeElement;
    out.push(was + ' > ' + (now && now.classList.contains('gm-showall') ? now.textContent : (now ? now.tagName : '')));
  }
  return JSON.stringify(out);
}())"""


def showall_problems(b, base):
    """The keyboard stays on Show all when the table is built again (Show fewer, then Show all again), and the
    count's live line sits outside the foldable table, so a screen reader hears it while the table is folded."""
    bad = []
    inview_page(b, base, S, "#chr22:40.5-43.5M")  # 21 rows with its genes
    if b.eval("(function () { var l = document.querySelector('#genome-map .gm-inview-live'); "
              "return !l || !!l.closest('details'); }())"):
        bad.append("the count's live line is missing or inside the foldable table")
    said = json.loads(b.eval(SHOWALL_JS))
    if said is None:
        return bad + ["no Show all at chr22:40.5-43.5M: the check tests nothing"]
    if sorted(said) != ["Show all > Show fewer", "Show fewer > Show all"]:
        bad.append(f"after Show all and Show fewer the keyboard is on {said}")
    return bad


def fly_problems(b, base, finding_only=False):
    """A row's button flies the map to its feature: a finding's row from the genome level to its place with its
    card, the row then current; a run's row to a view that holds the run, with no card. FINDING_ONLY: the
    finding's row alone."""
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
    if finding_only:
        return bad
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


def strip_problems(b, base, out, phone_only=False, click_only=False):
    """The strip: a click opens a chromosome whole, a drag elsewhere brushes a new view on the chromosome where it
    started, a drag on the brush pans; at 390 px it shows 12 px high and a click there does nothing (a picture
    only), and the table starts folded with its count in its summary. PHONE_ONLY: the 390 px half alone;
    CLICK_ONLY: the first click alone."""
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
        if click_only:
            return bad
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


PLACES = (("", "chr7", "chr7", "chromosome"), ("#chr7", "#chr11:4.5-6.5M", "chr11", "region"),
          ("#chr7", "genome", "", "genome"))


def place_find_problems(b, base, places=PLACES):
    """The search box takes a place as the address writes it: chr7 opens the chromosome, #chr11:4.5-6.5M flies to
    the region, genome steps back out to the whole genome (PLACES only when given)."""
    bad = []
    for link, q, chrom, lv in places:
        map_page(b, base, S, link)
        found = map_find(b, q)
        if not found or found[0]["kind"] != "place":
            bad.append(f"{q}: lines {found}")
            continue
        b.eval("document.getElementById('gm-find').dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', bubbles: true})); true")
        for _ in range(5):  # a busy machine eases the camera slower, and the margin follows the settled camera
            c = map_now(b)
            if c["level"] == lv and (not chrom or c["chrom"] == chrom):
                break
        if c["level"] != lv or (chrom and c["chrom"] != chrom):
            bad.append(f"{q}: Enter flies to {c['chrom']} at {c['level']}, want {chrom or 'the genome'} at {lv}")
    return bad


def card_focus_problems(b, base):
    """Tab from a finding's mark to its card's link: the view settles again and the keyboard stays on the link."""
    map_page(b, base, S, "#chr11")
    if not b.eval("(function () { var a = document.querySelector('#genome-map .gm-marks a[data-id]'); "
                  "if (!a) { return false; } a.focus(); return true; }())"):
        return ["no finding's mark on chr11: the check tests nothing"]
    time.sleep(0.6)  # the focus picks the finding; its card is drawn and the view settles
    if not b.eval("(function () { var l = document.querySelector('#genome-map .gm-card:not([hidden]) a'); "
                  "if (!l) { return false; } l.focus(); return document.activeElement === l; }())"):
        return ["the finding's card has no link to focus: the check tests nothing"]
    time.sleep(1.0)  # the mark's blur asks for a frame; the margin settles 140 ms after it
    if not b.eval("!!(document.activeElement && document.activeElement.closest && "
                  "document.activeElement.closest('#genome-map .gm-card'))"):
        return ["the keyboard left the card's link once the view settled: "
                + str(b.eval("document.activeElement ? document.activeElement.tagName : ''"))]
    return []


ROWS_JS = r"""JSON.stringify(Array.prototype.map.call(
  document.querySelectorAll('#gm-second .gm-inview tbody tr[data-layer=findings] .gm-go'),
  function (x) { return x.textContent; }))"""


def reload_table_problems(b, base, out):
    """A reload whose document gives a finding new words under the same place and id: the In view table shows the
    new words, not the ones it had."""
    path = os.path.join(out, "runs", S, "map.json")
    other = os.path.join(out, "runs", S, "map-reload.json")
    try:
        map_page(b, base, S)
        b.eval(SECOND_JS)
        if not b.wait("document.getElementById('gm-second').dataset.pgpMap === 'drawn'"):
            return ["the second figure never drew"]
        time.sleep(0.6)
        before = json.loads(b.eval(ROWS_JS))
        doc = read_map(path)
        titles = doc["per"][S]["findings"]["title"]
        if not before or titles.count(before[0]) != 1:
            return [f"no finding row to rename: the check tests nothing ({before})"]
        titles[titles.index(before[0])] = "Planted new words for this finding"
        write_map(other, doc)
        b.eval("pgpSecond.o.url = 'map-reload.json'; pgpSecond.g.refresh('k2'); true")
        if b.wait("%s.indexOf('Planted new words for this finding') >= 0" % ROWS_JS.replace("JSON.stringify", ""), 5):
            return []
        return [f"after the reload the table still reads {json.loads(b.eval(ROWS_JS))[:3]}"]
    finally:
        if os.path.exists(other):
            os.remove(other)


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
                "return i > 0 || !!(row && row.long);", "return i > 0;",
                lambda: predicate_problems(b, base, out, og, runs=[S]))
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
    check("map: a closed result is no row of the In view table and no tick of the strip (APOE on an adult's page); "
          "the son's page, closed by the report gate, has no map", not bad, bad)
    planted("map.js whose table lists every finding, closed or not",
            "drawn().forEach(function (r) { if (meets(", "M.findings.forEach(function (r) { if (meets(",
            lambda: inview_apoe_problems(b, base))
    planted("map.js whose strip ticks every finding, closed or not",
            "      drawn().forEach(function (r) {\n        if (!drawsAt('findings', r, 'genome'))",
            "      M.findings.forEach(function (r) {\n        if (!drawsAt('findings', r, 'genome'))",
            lambda: inview_apoe_problems(b, base))
    bad = sort_problems(b, base)
    check("map: the In view table sorts by place with the findings first (aria-sort ascending on Where), Where again "
          "reverses it, and Size sorts by size", not bad, bad)
    planted("map.js whose table sorts by place alone",
            "return (a.layer === 'findings' ? 0 : 1) - (b.layer === 'findings' ? 0 : 1) || place(a, b);",
            "return place(a, b);", lambda: sort_problems(b, base, default_only=True))
    bad = announce_problems(b, base)
    check("map: the In view table announces its count once when it changes, never while it stays the same", not bad, bad)
    planted("map.js whose table announces its count on every settle",
            "if (n !== announced) { announced = n; inLive", "if (true) { announced = n; inLive",
            lambda: announce_problems(b, base, pans_only=True))
    bad = showall_problems(b, base)
    check("map: the keyboard stays on Show all and Show fewer, and the count is said outside the foldable table",
          not bad, bad)
    planted("map.js whose Show all drops the keyboard", "if (more) { more.focus(", "if (false) { more.focus(",
            lambda: showall_problems(b, base))
    planted("map.js whose count's live line folds with the table", "main.insertBefore(inLive, inview.nextSibling);",
            "inview.insertBefore(inLive, inWrap);", lambda: showall_problems(b, base))
    bad = fly_problems(b, base)
    check("map: a row of the In view table flies the map to its feature (a finding with its card) and becomes the "
          "current row", not bad, bad)
    planted("map.js whose rows fly without the finding's card",
            "if (x.layer === 'findings') { pick(x.r.id); } else { pick(null); }", "pick(null);",
            lambda: fly_problems(b, base, finding_only=True))
    bad = strip_problems(b, base, out)
    check("map: the strip opens a chromosome on a click, brushes a new view and pans on a drag; at 390 px it is a "
          "12 px picture and the table starts folded", not bad, bad)
    planted("map.js whose strip ignores a click",
            "if (i >= 0) { openChrom(i, undefined, undefined, true); }\n      } else if (b0.mode",
            "void i;\n      } else if (b0.mode", lambda: strip_problems(b, base, out, click_only=True))
    planted("map.js whose strip takes clicks on a phone",
            "if (!M || phone || !b0) { return; }", "if (!M || !b0) { return; }",
            lambda: strip_problems(b, base, out, phone_only=True),
            more=[("if (!M || phone || e.button !== 0 || level() === 'genome') { return; }",
                   "if (!M || e.button !== 0 || level() === 'genome') { return; }")])
    bad = place_find_problems(b, base)
    check("map: the search box takes a place as the address writes it (chr7, #chr11:4.5-6.5M, genome) and Enter "
          "flies there", not bad, bad)
    planted("map.js whose search takes a place only with a colon",
            "var p = /^(genome|chr[0-9XYMT]+(:.*)?)$/i.test(pq) ? parsePlace(pq) : null;",
            "var p = /^chr/i.test(pq) && /[:]/.test(pq) ? parsePlace(pq) : null;", lambda: place_find_problems(b, base, PLACES[:1]))
    bad = card_focus_problems(b, base)
    check("map: the keyboard stays on the card's link when the view settles again", not bad, bad)
    planted("map.js whose card is built again on every settle",
            "if (sig === cardSig && !card.hidden) { return; }", "", lambda: card_focus_problems(b, base))
    bad = reload_table_problems(b, base, out)
    check("map: a reloaded document with new words for a finding under the same key shows them in the In view table",
          not bad, bad)
    planted("map.js whose table keeps its rows across a reload",
            "        tSig = '';  // a new document may", "        void 0;  // a new document may",
            lambda: reload_table_problems(b, base, out))


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
                                          (S, "#chr11:5.150-5.350M", "gene", "gene"), (OTHERS[1], "", "closed", "gate"),
                                          (S, "#chr22:40.5-43.5M", "inview", "region")):
                if name == "closed":
                    b.go('about:blank'); b.wait("location.href === 'about:blank'")
                    b.go(base+'/runs/'+run+'/')
                    if not b.wait("document.documentElement.dataset.pgpGated === 'true'"):
                        print(f"ERROR: {run}: no closed report gate: no screenshot")
                        continue
                    b.eval("document.getElementById('report-gate').scrollIntoView({block: 'start'}); true")
                    time.sleep(0.5)
                    png = b.call("Page.captureScreenshot", {"format": "png"}, b.session)["data"]
                    p = os.path.join(dest, f"map-{name}-{w}x{h}.png")
                    with open(p, "wb") as fh:
                        fh.write(base64.b64decode(png))
                    taken.append(p)
                    continue
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
          "its card); APOE, closed, is never a line, by its gene, id or title", not bad, bad)
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
          "card; Esc closes the card first; it never leads to APOE, closed", not bad, bad)
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
    og = og_binary() or build_live_demo.OG
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
                    # Whether the son is a child (child, the server's one rule) is a routing control, as his url is.
                    d['runs'][son] = dict(changed_family['runs'][son], **dict(progress(d['runs'][son]), url=d['runs'][son]['url'],
                                                                              child=d['runs'][son].get('child')))
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
                            said = type_search(b, 'KCNH2')
                            if SEARCH_CLOSED not in said:
                                # A child's opened report searches as before the gate: it says what it leaves out.
                                transitions.append((count, order, stage, 'search lost its closed note'))
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
            gate_slow_release(b, base, out, acts, visit, payloads)
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
        # Opened on a guessed address, the report opens on the run's page, as main answers no such finding.
        landed = {}
        for fid in (acts[0]['id'], 'zz-not-a-finding'):
            visit('/runs/'+son+'/finding/'+fid+'/')
            open_gate(b)
            b.wait("document.readyState === 'complete' && document.documentElement.dataset.pgpStatus === 'ready' && !document.documentElement.dataset.pgpGated")
            landed[fid] = b.eval('location.pathname')
        check('gate addresses: opened on a guessed finding address, the report opens on the run page',
              landed == {acts[0]['id']: '/runs/'+son+'/finding/'+acts[0]['id']+'/', 'zz-not-a-finding': '/runs/'+son+'/'}, landed)
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


def gate_slow_release(b, base, out, acts, visit, payloads):
    """Opening the last act releases the report however slowly the reloaded page and its results arrive: the opening
    never depends on the time. A gate that accepts the opening only within 1.5 s must fail here."""
    son = OTHERS[1]
    payloads(2)

    def release():
        visit('/runs/'+son+'/')
        open_act(b, acts[0]['title'])
        Static.DELAY['/runs/'+son+'/'] = 4
        Static.DELAY['/runs/'+son+'/results.json'] = 6
        try:
            b.eval("document.querySelector('[data-opening^=closed]').click(); true")
            return b.wait("!document.documentElement.dataset.pgpGated && !!document.querySelector('#results .finding')", timeout=25)
        finally:
            Static.DELAY.clear()
    released = release()
    check('gate release: a reload whose page takes 4 s and whose results.json takes 6 s still releases the report',
          released, released)
    js = os.path.join(out, 'assets', 'core.js')
    source = open(js).read()
    pairs = (("sessionStorage.setItem(TRANSFER, JSON.stringify({ path: target,", "sessionStorage.setItem(TRANSFER, JSON.stringify({ at: Date.now(), path: target,"),
             ("t.path === location.pathname) { ticket = t; }", "t.path === location.pathname && Date.now() - t.at < 1500) { ticket = t; }"))
    try:
        planted = source
        for x, y in pairs:
            assert source.count(x) == 1, 'mutation anchor moved: ' + x
            planted = planted.replace(x, y)
        with open(js, 'w') as f:
            f.write(planted)
        control('gate release: an opening accepted only within 1.5 s', bool(release()))
    finally:
        with open(js, 'w') as f:
            f.write(source)


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
        # A release's own location.replace can land after about:blank; go there again until it stays.
        assert any(b.go('about:blank') or b.wait("location.href === 'about:blank'") for _ in range(3)), 'about:blank never stayed'
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

    def ticket():
        # The opening ticket crosses one load only: after closing an opened report again, the closed page keeps no
        # ticket (it names the opened results' places in the report and a digest of every result's record).
        visit(); open_gate(b)
        assert b.wait("!!document.querySelector('#results .finding')"), 'full report never loaded'
        fid = next(f['id'] for f in data['findings'] if f['level'] == 'act')
        b.eval('PGPUI.closeAct('+json.dumps(fid)+'); true')
        assert b.wait("document.documentElement.dataset.pgpGated === 'true' && document.readyState === 'complete'"), 'never closed'
        left = b.eval("[sessionStorage.getItem('og-report-opening'), window.name]")
        seen['ticket'] = left[0]
        return left[0] is None and not left[1].startswith('og-report-opening:')

    # What <html> may keep once a child's open report is cleared: the gate's own two attributes.
    LEFT = ("Array.prototype.map.call(document.documentElement.attributes, function (a) { return a.name; }).filter(function (n) {"
            " return (n.indexOf('data-pgp-') === 0 && n !== 'data-pgp-status' && n !== 'data-pgp-gated') || n === 'data-finding'; })")

    def opened():
        visit(); open_gate(b)
        assert b.wait("!!document.querySelector('#results .finding') && !!document.documentElement.dataset.pgpResults"), 'full report never loaded'

    def close_attrs():
        # Read in the same turn as the close, before the document is replaced.
        opened()
        fid = next(f['id'] for f in data['findings'] if f['level'] == 'act')
        seen['attrs'] = b.eval('PGPUI.closeAct(' + json.dumps(fid) + '); ' + LEFT)
        return not seen['attrs']

    def pagehide_attrs():
        opened()
        seen['attrs'] = b.eval("window.dispatchEvent(new PageTransitionEvent('pagehide', { persisted: false })); " + LEFT)
        return not seen['attrs'] and b.eval("!document.querySelector('#results')")

    def back():
        # Leave the open report for another page, then Back: from the back-forward cache the page loads afresh, closed.
        # The page is served without no-store here, so Chrome keeps it in that cache (a no-store page never enters it).
        Static.BFCACHE.add('/runs/%s/' % son)
        try:
            opened()
            b.go(base + '/')
            assert b.wait("document.documentElement.dataset.pgpStatus === 'ready' && !!document.documentElement.dataset.pgpPeople"), 'home never drawn'
            b.eval('history.back(); true')
            b.wait("location.pathname === '/runs/%s/' && document.readyState === 'complete'" % son)
            got = b.wait("document.documentElement.dataset.pgpGated === 'true' && performance.getEntriesByType('navigation')[0].type", timeout=10)
            seen['back'] = got or b.eval("[performance.getEntriesByType('navigation')[0].type, document.body.textContent.length]")
        finally:
            Static.BFCACHE.clear()
        return got == 'reload'

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
         ['return true;'], cleanup, lambda: seen.get('query') == 'KCNH2'),
        ('keep the opening ticket after reading it', ['    sessionStorage.removeItem(TRANSFER);\n'], [''], ticket,
         lambda: seen.get('ticket') is not None),
        ('keep the page attributes after closing', ["if (a.name.indexOf('data-pgp-') === 0 || a.name === 'data-finding') { root.removeAttribute(a.name); }"],
         [''], close_attrs, lambda: bool(seen.get('attrs'))),
        ('keep the open report on pagehide', ['ownsPage = true; cleaned = true; eventSources.forEach(function (es) { es.close(); }); cleanPage(); }'],
         ['ownsPage = true; cleaned = true; eventSources.forEach(function (es) { es.close(); }); }'], pagehide_attrs,
         lambda: bool(seen.get('attrs'))),
        ('show the cleared report again from the back-forward cache', ["if (e.persisted && cleaned) { location.reload(); }"], [''],
         back, lambda: seen.get('back') != 'reload')]
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


# Each network read a page script starts, in the order it starts it, below every wrapper of the page's own: the
# path, and for fetch its cache mode and credentials. Installed before any script of the page runs.
REQUEST_LOG = r"""(function () {
  var log = window.ogRequests = [], f = window.fetch, E = window.EventSource;
  window.fetch = function (u, o) {
    log.push(['fetch', new URL(u && u.url || u, location.href).pathname, (o && o.cache) || 'default', (o && o.credentials) || 'same-origin']);
    return f.apply(this, arguments);
  };
  if (E) {
    var W = function (u, o) { log.push(['events', new URL(u, location.href).pathname]); return new E(u, o); };
    W.prototype = E.prototype; W.CONNECTING = E.CONNECTING; W.OPEN = E.OPEN; W.CLOSED = E.CLOSED;
    window.EventSource = W;
  }
}());"""


def bare_core(source):
    """core.js with the report gate cut out: the transparent reference a page without the gate is compared to."""
    marker = '  // --- the report gate'
    assert source.count(marker) == 1, 'gate marker moved'
    return source[:source.index(marker)] + '''
  PGPUI.ready = Promise.resolve();
  PGPUI.gateOwnsPage = function () { return false; };
  PGPUI.openedActs = function () { return {}; };
  PGPUI.closeAct = function () { return false; };
}());
'''


SEARCH_CLOSED = 'Results this page keeps closed are left out of the search until they are opened.'


# The gate transparency check compares two renders' DOM as text. A browser keeps an element's attributes in the order
# they were set, and the page's scripts set the <html> data-pgp-* attributes as their files arrive, in any order. So
# each element of the clone gets its attributes back sorted by name before it is serialized: the order no longer
# counts, while a changed value, an added or a missing attribute still makes the two texts differ.
SORT_ATTRS = """[clone].concat(Array.from(clone.querySelectorAll('*'))).forEach(function (n) {
            var as=Array.from(n.attributes).sort(function (a, z) { return a.name < z.name ? -1 : a.name > z.name ? 1 : 0; });
            as.forEach(function (a) { n.removeAttributeNode(a); });
            as.forEach(function (a) { n.setAttributeNodeNS(a); });
          });"""


def type_search(b, term):
    """Type TERM in the header search once its index has loaded; the text of what it shows."""
    assert b.wait("document.querySelector('input[type=search]') && !document.querySelector('input[type=search]').disabled"), 'search never loaded'
    b.eval("var s=document.querySelector('input[type=search]'); s.value=%s; s.dispatchEvent(new Event('input')); true" % json.dumps(term))
    assert b.wait("document.documentElement.dataset.pgpSearch !== undefined"), 'search never answered'
    return b.eval("document.querySelector('.search-results').textContent")


def gate_transparency(b, base, out):
    """The gate leaves adult and no-act pages equal to the renderer without it: DOM, attributes, address, canvas and
    the requests, counted, each with its fetch options, and the files the gate reads first in their order."""
    son = OTHERS[1]
    js = os.path.join(out, 'assets', 'core.js')
    fam = os.path.join(out, 'runs', FAM, 'family.json')
    files = [os.path.join(out, 'runs', S, 'results.json'), os.path.join(out, 'runs', son, 'results.json'),
             os.path.join(out, 'runs.json'), fam]
    source = open(js).read()
    saved = {p: open(p).read() for p in files}
    bare = bare_core(source)

    def capture(route, apoe=False, timeout=20, search=None, family=True):
        b.go('about:blank'); assert b.wait("location.href === 'about:blank'")
        b.go(base + route)
        if route == '/':
            assert b.wait("document.documentElement.dataset.pgpStatus === 'ready' && !!document.documentElement.dataset.pgpPeople", timeout=timeout), 'home never drawn'
        elif '/finding/' in route:
            assert b.wait("document.documentElement.dataset.pgpStatus === 'ready' && !!document.querySelector('h1')", timeout=timeout), 'finding never drawn'
        else:
            assert b.wait("document.documentElement.dataset.pgpResults && !!document.querySelector('[data-pgp-map=drawn]')", timeout=timeout), 'open report not drawn'
            assert b.wait("(%s || !!document.querySelector('#rail a[href*=DEMO-FAMILY]')) && !!document.documentElement.dataset.pgpCharts"
                          % ('false' if family else 'true')), 'family or charts never settled'
        if apoe:
            # A redraw while the address names a place: the contents rail must come back as it did.
            b.eval("document.querySelector('details[data-gate=apoe] summary').click(); true")
        if search:
            type_search(b, search)
        time.sleep(0.5)
        got = b.eval("""(function () {
          var clone=document.documentElement.cloneNode(true);
          clone.querySelectorAll('[data-pgp-map-paint-ms],[data-pgp-map-parse-ms]').forEach(function (n) {
            n.removeAttribute('data-pgp-map-paint-ms'); n.removeAttribute('data-pgp-map-parse-ms');
          });
          """ + SORT_ATTRS + """
          return [clone.outerHTML, location.hash, Array.from(document.querySelectorAll('canvas')).map(function(c){return c.toDataURL();}),
            performance.getEntriesByType('resource').map(function (r) { return new URL(r.name).pathname; }), window.ogRequests || []];
        }())""")
        # Scripts and the views' own files load in parallel, in any order, with or without the gate: count every
        # request with its options, and keep the order of the files the gate reads first.
        reads = got.pop()
        first = [r for r in reads if r[1].endswith(('/state.json', '/frames.json', '/family.json', '/results.json', '/events'))]
        # Chrome asks for a fallback favicon after load whenever it likes; no script of the page asks for it.
        got[3] = [sorted(p for p in got[3] if p != '/favicon.ico'), sorted(json.dumps(r) for r in reads), first]
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

    def no_family():
        # A family that cannot be read (its family.json answers 404): the page asks for it once, as without the gate.
        os.remove(fam)

    def no_act_no_family():
        no_act()
        no_family()

    term = located_term(out)
    first = next(f['id'] for f in json.loads(saved[files[0]])['findings'] if f['level'] != 'act')
    cases = [('adult', None, '/runs/'+S+'/', {}), ('adult, a redraw at a section address', None, '/runs/'+S+'/#s3', {'apoe': True}),
             ('adult finding page', None, '/runs/'+S+'/finding/'+first+'/', {}),
             ('child with no act, a search typed', no_act, '/runs/'+son+'/', {'search': term}),
             ('home with no gated child', no_act, '/', {}),
             ('fresh adult with an act and a same-gene ranking', synthetic, '/runs/'+S+'/', {}),
             ('adult whose family cannot be read', no_family, '/runs/'+S+'/', {'family': False}),
             ('adult finding page whose family cannot be read', no_family, '/runs/'+S+'/finding/'+first+'/', {}),
             ('child with no act whose family cannot be read', no_act_no_family, '/runs/'+son+'/', {'family': False})]
    log = b.call('Page.addScriptToEvaluateOnNewDocument', {'source': REQUEST_LOG}, b.session)['identifier']
    try:
        for label, prep, route, how in cases:
            for p in files:
                with open(p, 'w') as f:
                    f.write(saved[p])
            if prep:
                prep()
            pair = []
            for text in (bare, source):
                with open(js, 'w') as f:
                    f.write(text)
                pair.append(capture(route, **how))
            detail = [('dom', 'hash', 'canvas', 'requests')[i] for i in range(4) if pair[0][i] != pair[1][i]]
            for k in range(3):
                if pair[0][3][k] != pair[1][3][k]:
                    x, y = pair[0][3][k], pair[1][3][k]
                    detail += [('without the gate only', [r for r in x if r not in y]), ('with it only', [r for r in y if r not in x]), x, y]
            if pair[0][0] != pair[1][0]:
                at = next((i for i, (a, z) in enumerate(zip(pair[0][0], pair[1][0])) if a != z), 0)
                detail += [pair[0][0][max(0, at-80):at+120], pair[1][0][max(0, at-80):at+120]]
            check('gate transparency: '+label+' has identical DOM, attributes, address, canvas and requests', pair[0] == pair[1], detail)
            if how.get('search'):
                # A child's page keeps results closed until opened, with or without results to act on, and says so
                # under the hits, as before the gate.
                check('gate transparency: a search on a child with no act says that closed results are left out',
                      SEARCH_CLOSED in pair[1][0], term)
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
        b.call('Page.removeScriptToEvaluateOnNewDocument', {'identifier': log}, b.session)
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
