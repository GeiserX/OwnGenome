#!/usr/bin/env python3
"""The report's charts on the live page, on DEMO-001, the invented demo sample.

The charts' own checks (bin/report_charts.py, og report charts) are cargo
tests now: rust/crates/og-render/tests/charts.rs. What stays here is the
browser check of the live page, which needs the live site's Python server:

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
import io
import json
import os
import re
import shutil
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import render_report  # noqa: E402

FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS += 1




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


def page_problems(b, base, run):
    """What is wrong with the live page's charts after a sign-in, from a click on a ClinVar tick."""
    import test_owngenome_live as tl
    b.go(f"{base}/runs/{run}/")
    if not b.wait("location.pathname === '/login' && !!document.getElementById('password')"):
        return ["no login form"]
    b.eval(f"document.getElementById('password').value = {json.dumps(tl.PASSWORD)}; "
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
    sys.path.insert(0, os.path.join(REPO, "tests"))
    import owngenome_live
    import test_owngenome_live as tl
    chrome = tl.find_chrome()
    if not chrome:
        ci = os.environ.get("CI") == "true"
        print(("[FAIL] " if ci else "[SKIP] ") + "10. no Chrome or Chromium to load the live page in")
        if ci:
            check("10. Chrome in CI", False)
        return
    root = tl.make_root(work)
    shutil.copy(charts_file, os.path.join(root, tl.RUN, f"{tl.RUN}_charts.json"))
    pw = os.path.join(work, "password")
    with open(pw, "w") as f:
        f.write(owngenome_live.hash_password(tl.PASSWORD) + "\n")
    bare = os.path.join(work, "assets-without-charts")
    shutil.copytree(os.path.join(REPO, "assets", "live"), bare)
    os.remove(os.path.join(bare, "plugins", "charts.js"))
    for tag, assets in (("live", os.path.join(REPO, "assets", "live")), ("control", bare)):
        argv = ["--root", root, "--state-dir", os.path.join(work, f"state-{tag}"), "--port", "0", "--ingest-port",
                "0", "--password-file", pw, "--assets", assets, "--quiet"]
        with contextlib.redirect_stderr(io.StringIO()):
            srv = tl.serve(argv, owngenome_live.App())
        b = tl.Chrome(chrome, os.path.join(work, f"chrome-{tag}"))
        try:
            problems = page_problems(b, f"http://127.0.0.1:{srv[1].server_address[1]}", tl.RUN)
        finally:
            b.close()
            tl.stop(srv)
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
