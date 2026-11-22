#!/usr/bin/env python3
"""The report's charts in the HTML report and on the live page, on the invented
sections of tests/fixtures/report_charts/ and on DEMO-001, the invented demo
sample.

The cargo tests in rust/crates/og-render/tests/charts.rs check og report
charts and hold it to the output of bin/report_charts.py frozen at 65095e97.
What stays here checks what bin/render_report.py still ships:

  7. the HTML report: no <script>, every href="#..." of the page resolves to
     an id, each chart sits in its section below that section's caveat, and
     the hits table carries the id the charts link to. It stays until
     og report render embeds the figures (charts::figure in html.rs); that
     lane ports it to its cargo tests;
  8. two renders of the demo give the same charts.json apart from
     generated_utc, with its schema version and the five charts;
 10. the live page (bin/owngenome_live.py with bin/owngenome_live_charts.py and
     assets/live/plugins/charts.js): the plugin serves <run>_charts.json and
     names it in state.json; in a headless Chrome, signed in, the five charts
     appear, a click on a ClinVar tick opens the panel with its caveat and a
     viewer/#<locus> link, the page fetched nothing off its origin and the
     policy blocked nothing. Red control: the same page without charts.js
     fails the check. Without bin/owngenome_live.py on the branch this part says
     SKIP; without Chrome it says SKIP, and FAIL in CI.

The demo's charts file it serves comes from bin/render_report.py --charts-json.
"""
import contextlib
import copy
import io
import json
import os
import re
import shutil
import sys
import tempfile
import threading

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import render_report  # noqa: E402

FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS += 1


def fixture():
    with open(os.path.join(REPO, "tests", "fixtures", "report_charts", "sections.json")) as f:
        return {"sample": "S", "generated_utc": "2026-01-01 00:00 UTC", "sections": json.load(f)["sections"]}


def href_problems(page):
    ids = set(re.findall(r'\sid="([^"]+)"', page))
    return sorted({h for h in re.findall(r'href="#([^"]*)"', page) if h not in ids})


def demo_checks(work):
    sys.path.insert(0, os.path.join(REPO, "tests", "demo"))
    import make_demo_sample
    g = f"{work}/demo"
    os.makedirs(g)
    make_demo_sample.traits = lambda d: None
    with contextlib.redirect_stdout(io.StringIO()):
        make_demo_sample.main(["--genome-dir", g, "--records", "20000"])
    sd = f"{g}/DEMO-001"
    outs = []
    for i in (1, 2):
        with contextlib.redirect_stdout(io.StringIO()):
            rc = render_report.main(["--sample", "DEMO-001", "--sample-dir", sd, "--json", f"{work}/s{i}.json",
                                     "-o", f"{work}/r{i}.html", "--charts-json", f"{work}/c{i}.json"])
        with open(f"{work}/c{i}.json") as f:
            outs.append(json.load(f))
        check(f"8. render {i} of the demo exits 0", rc == 0)
    a, b = (dict(x, generated_utc=None) for x in outs)
    check("8. charts.json of two renders is the same apart from generated_utc", a == b and outs[0]["charts"])
    b2 = copy.deepcopy(b)
    b2["charts"][0]["marks"][0]["value"] = -1
    check("8. red control: one changed value is caught", a != b2)
    check("8. charts.json carries its schema version and the five charts of the demo",
          outs[0]["schema_version"] == 1 and [c["id"] for c in outs[0]["charts"]]
          == ["clinvar-stars", "genome", "mito", "ancestry-fine", "prs"], [c["id"] for c in outs[0]["charts"]])

    # 7. the offline page, with the fixture's sections on the demo (53 hits, ROH, poly-C)
    with open(f"{work}/s1.json") as f:
        s = json.load(f)
    s["sections"].update(fixture()["sections"])
    page = render_report.html_report(s)
    check("7. the HTML report has no <script", "<script" not in page.lower())
    check("7. every href=\"#...\" resolves to an id of the page", href_problems(page) == [], href_problems(page))
    check("7. red control: a link to a missing id is caught",
          href_problems(page.replace('id="clinvar-hits"', 'id="gone"')) == ["clinvar-hits"])
    for sid, key in (("clinvar", "CLINVAR_FREQ_NOTE"), ("roh", "ROH_NOTE"), ("mito", "MITO_NOTE"),
                     ("ancestry-fine", "ANCESTRY_FINE_NOTE"), ("prs", "PRS_ANCESTRY_NOTE")):
        sec = re.search(rf'<section class="sec[^"]*" id="{sid}">(.*?)\n  </section>', page, re.S).group(1)
        cav = sec.find(render_report.E(getattr(render_report, key)))
        fig = sec.find('<figure class="chart">')
        check(f"7. section {sid}: one chart, below its caveat",
              sec.count('<figure class="chart">') == 1 and -1 < cav < fig, (cav, fig))
    check("7. the hits table carries the id the charts link to, and its rows keep their shape",
          page.count('<a id="clinvar-hits"></a>') == 1 and len(re.findall(r"<tr><td>chr", page)) >= 50)


# The Python live server and its genome folder for check 10, moved from
# tests/test_owngenome_live.py (its checks are cargo tests of og live serve now).
FIX = os.path.join(REPO, "tests", "fixtures", "live")
RUN = "DEMO-001"
PASSWORD = "synthetic test password 7"


def make_root(work):
    root = os.path.join(work, "genome")
    d = os.path.join(root, RUN)
    shutil.copytree(os.path.join(FIX, RUN), d)
    os.makedirs(os.path.join(d, "clinvar"))
    os.makedirs(os.path.join(d, "extra"))
    with open(os.path.join(d, "clinvar", "hits.tsv"), "w") as f:
        f.write("gene\tvariant\nGENEA\tchr1:100 A>G\n")
    with open(os.path.join(d, "clinvar", "unlinked.tsv"), "w") as f:
        f.write("not linked\n")
    with open(os.path.join(d, "extra", "note.txt"), "w") as f:
        f.write("plugin file\n")
    outside = os.path.join(work, "outside.txt")
    with open(outside, "w") as f:
        f.write("outside the root\n")
    os.symlink(outside, os.path.join(d, "clinvar", "escape.tsv"))
    # The report links CPSR's folder, which holds its HTML report and its table
    os.makedirs(os.path.join(d, "cpsr"))
    with open(os.path.join(d, "cpsr", f"{RUN}.cpsr.grch38.html"), "w") as f:
        f.write("<!doctype html><title>CPSR</title>\n")
    with open(os.path.join(d, "cpsr", f"{RUN}.cpsr.grch38.classification.tsv.gz"), "w") as f:
        f.write("not linked\n")
    os.makedirs(os.path.join(d, "twopages"))
    for n in ("a.html", "b.html"):
        with open(os.path.join(d, "twopages", n), "w") as f:
            f.write("<!doctype html>\n")
    with open(os.path.join(d, f"{RUN}_report.html"), "w") as f:
        f.write('<!doctype html><html><head><style>body{color:#eee}</style></head><body>'
                '<a href="clinvar/hits.tsv">hits</a> <a href="clinvar/escape.tsv">escape</a> '
                '<a href="cpsr/">the cpsr folder</a> <a href="twopages/">two pages</a> '
                '<a href="https://example.org/x">outside</a></body></html>\n')
    return root


def serve(argv, app=None):
    import owngenome_live
    live, web, ingest = owngenome_live.start(argv, app)
    for s in (web, ingest):
        threading.Thread(target=s.serve_forever, daemon=True).start()
    return live, web, ingest


def stop(*servers):
    for live, web, ingest in servers:
        live.stopping = True
        web.shutdown()
        ingest.shutdown()
        web.server_close()
        ingest.server_close()


def page_problems(b, base, run):
    """What is wrong with the live page's charts after a sign-in, from a click on a ClinVar tick."""
    b.go(f"{base}/runs/{run}/")
    if not b.wait("location.pathname === '/login' && !!document.getElementById('password')"):
        return ["no login form"]
    b.eval(f"document.getElementById('password').value = {json.dumps(PASSWORD)}; "
           "document.querySelector('form').submit(); true")
    if not b.wait(f"location.pathname === '/runs/{run}/' && document.documentElement.dataset.pgpStatus === 'ready'"):
        return ["the run page did not load"]
    out = []
    n = b.wait("document.documentElement.dataset.pgpCharts", 15)
    if n != "5":
        out.append(f"{n} charts shown, not 5")
    else:
        b.eval("document.querySelector('figure[data-chart=\"genome\"] a[data-mark^=\"hit-\"]')"
               ".dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true})); true")
        side = b.wait("(function () { var d = document.querySelector('.chart-detail'); if (!d || d.hidden) { "
                      "return null; } var v = d.querySelector('a.viewer'); return JSON.stringify({text: d.textContent, "
                      "viewer: v ? v.getAttribute('href') : '', at: location.pathname}); })()", 10)
        side = json.loads(side) if side else {}
        if render_report.CLINVAR_FREQ_NOTE not in side.get("text", ""):
            out.append(f"the panel has no caveat: {side}")
        if not re.fullmatch(r"viewer/#chr[0-9XY]+:[0-9]+", side.get("viewer", "")):
            out.append(f"no viewer/#<locus> link: {side}")
        if side.get("at") != f"/runs/{run}/":
            out.append(f"the click left the page: {side}")
    off = b.eval("performance.getEntriesByType('resource').map(function (e) { return e.name; })"
                 f".filter(function (u) {{ return u.indexOf({json.dumps(base + '/')}) !== 0; }})")
    if off:
        out.append(f"fetched off the origin: {off}")
    blocked = b.eval("document.documentElement.dataset.pgpBlocked || ''")
    if blocked:
        out.append(f"the policy blocked: {blocked}")
    return out


def live_checks(work, charts_file):
    import owngenome_live_charts

    class FakeApp:
        globs, hooks = [], []

        def run_file(self, g):
            self.globs.append(g)

        def state_hook(self, name, fn):
            self.hooks.append((name, fn))
    fa = FakeApp()
    owngenome_live_charts.register(fa)
    d = os.path.join(work, "plugin", "RUNX")
    os.makedirs(d)
    before = owngenome_live_charts.charts_state(d, {}, [])
    shutil.copy(charts_file, os.path.join(d, "RUNX_charts.json"))
    after = owngenome_live_charts.charts_state(d, {}, [])
    check("10. the plugin serves *_charts.json and names the run's file in state.json, once it exists",
          fa.globs == ["*_charts.json"] and [n for n, _ in fa.hooks] == ["charts"] and before == {}
          and after.get("file") == "RUNX_charts.json" and isinstance(after.get("mtime"), int), (fa.globs, before, after))
    os.utime(os.path.join(d, "RUNX_charts.json"), ns=(after["mtime"] + 1000, after["mtime"] + 1000))
    check("10. a rewrite within the same second gives a new mtime key",
          owngenome_live_charts.charts_state(d, {}, [])["mtime"] == after["mtime"] + 1000)

    if not os.path.isfile(os.path.join(REPO, "bin", "owngenome_live.py")):
        print("[SKIP] 10. the live page in Chrome: bin/owngenome_live.py is not on this branch")
        return
    sys.path.insert(0, os.path.join(REPO, "tests", "browser"))
    import cdp
    import owngenome_live
    chrome = cdp.find_chrome()
    if not chrome:
        ci = os.environ.get("CI") == "true"
        print(("[FAIL] " if ci else "[SKIP] ") + "10. no Chrome or Chromium to load the live page in")
        if ci:
            check("10. Chrome in CI", False)
        return
    root = make_root(work)
    shutil.copy(charts_file, os.path.join(root, RUN, f"{RUN}_charts.json"))
    pw = os.path.join(work, "password")
    with open(pw, "w") as f:
        f.write(owngenome_live.hash_password(PASSWORD) + "\n")
    bare = os.path.join(work, "assets-without-charts")
    shutil.copytree(os.path.join(REPO, "assets", "live"), bare)
    os.remove(os.path.join(bare, "plugins", "charts.js"))
    for tag, assets in (("live", os.path.join(REPO, "assets", "live")), ("control", bare)):
        argv = ["--root", root, "--state-dir", os.path.join(work, f"state-{tag}"), "--port", "0", "--ingest-port",
                "0", "--password-file", pw, "--assets", assets, "--quiet"]
        with contextlib.redirect_stderr(io.StringIO()):
            srv = serve(argv, owngenome_live.App())
        b = cdp.Chrome(chrome, os.path.join(work, f"chrome-{tag}"))
        try:
            problems = page_problems(b, f"http://127.0.0.1:{srv[1].server_address[1]}", RUN)
        finally:
            b.close()
            stop(srv)
        if tag == "live":
            check("10. live page: five charts; a click on a ClinVar tick opens the panel with its caveat and a "
                  "viewer/# link; nothing fetched off the origin; nothing blocked", problems == [], problems)
        else:
            check("10. red control: the same page without charts.js is caught", problems != [], problems)


def main():
    with tempfile.TemporaryDirectory() as work:
        demo_checks(work)
        live_checks(work, f"{work}/c1.json")
    print(f"\n{FAILS} failed" if FAILS else "\nAll checks passed")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
