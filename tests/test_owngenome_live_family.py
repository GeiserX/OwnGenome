#!/usr/bin/env python3
"""bin/owngenome_live_family.py and assets/live/plugins/family.js on tests/fixtures/families_live.

Checks:
  1. the family folders (family_a with members.tsv, family_b with pedigree.tsv
     only) are listed as runs of kind "family"; the sample folder child_a is
     not;
  2. family_a's state: its members, each child with its parents and its
     chromosome map (child_a has one, child_b none), the relationship check,
     the shared segments, the couple screen counted from step 39's table, the
     plan, and exactly the result files that exist;
  3. child_a's state links to family_a; a sample in no family gets null; while a
     family's pedigree cannot be read, a sample no readable family lists gets
     "unknown" (it may be that family's child), and child_a still its family;
  4. GET /family/<family>/pedigree.svg is the drawing of the declared parents;
     a listed file is served as written; notes.txt (in the folder, written by
     no step), a name with '..', a symlink out of the folder, a sample's
     folder and an unknown family are 404;
  5. family.js loads nothing from outside, sets no innerHTML and has no eval;
     node --check reads it when node is installed;
  6. with bin/owngenome_live.py present (the live-server change), the same through
     the real server: the run list, state.json's plugins.family, the route
     behind the login, and a family page that loads the script; in a headless
     Chrome (skipped without Chrome, except in CI) the card, the pedigree and a
     child's map load and the page's policy blocks nothing; then a 1-in-4 shared
     row draws the four children and a 1-in-2 row gives its chance in words with no
     drawing, a missing trio file and one linking out of the folder are worded
     apart (the link's target never shown), and the side by side lists every
     blood group system with a stale value marked in its cell;
  7. the demo family's report-gate property and two children opening on their own (test_live_demo.py's
     gate_checks on the family page), here so that browser time is not added to test_live_demo.py's;
  8. the real server answers every finding id of a child with a result to act on with the same page, and keeps
     the 404 for an adult's id that is not a finding.
  9. the report gate on the real server, in live mode: one rule of who is a child; the home page counts a
     child's results to act on from the summary its results.json reads (a stale summary.json too); family.json
     has no flagged-card count for a child whose report is closed, so a parent's page is the same whichever card
     holds the child's results; and an adult's and a no-act child's pages make the requests of the page without
     the gate, in its order, the event stream included, and are the same page (the genome map drawn with OG); an
     adult's page at a section address is the page without the gate whichever of the family and the results
     answers first; each member run of family.json says whether it is a child by the one rule; a no-act child
     with a stale summary.json keeps its numbers on the home page; a family file that cannot be read gives a
     fixed 500 and the server then reads that family as unreadable; on the family page a gated child reads
     "Report closed" and, once opened there, the real count; a child of one family is closed on another family's
     page where it is a partner.
Checks 2 and 4 are run again on a copy of the fixture with two parents swapped
and must fail there: a check that cannot fail is no check.

Run: python3 tests/test_owngenome_live_family.py
"""
import http.client
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import types
import xml.etree.ElementTree as ET

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import owngenome_live_family  # noqa: E402

FIX = os.path.join(REPO, "tests", "fixtures", "families_live")
JS = os.path.join(REPO, "assets", "live", "plugins", "family.js")
NS = "{http://www.w3.org/2000/svg}"
FAILS = []

CHILDREN = [{"child": "child_a", "father": "father_a", "mother": "mother_a",
             "segments": "/family/family_a/child_a_segments.svg"},
            {"child": "child_b", "father": "father_a", "mother": "mother_a", "segments": None},
            {"child": "father_a", "father": "gf_a", "mother": "gm_a", "segments": None},
            {"child": "mother_a", "father": None, "mother": "gm_m", "segments": None}]
FILES = ["child_a_segments.svg", "couple_father_a_mother_a.md", "couple_father_a_mother_a.tsv", "family_plan.tsv",
         "family_report.md", "ibd_pairs.tsv", "ibd_report.md", "members.tsv", "relationships.tsv"]
DESCENT = {"child_a": ("father_a", "mother_a"), "child_b": ("father_a", "mother_a"), "father_a": ("gf_a", "gm_a"),
           "mother_a": (".", "gm_m")}


def check(name, ok, detail=""):
    print(("[ OK ] " if ok else "[FAIL] ") + name + ("" if ok else "  %s" % (detail,)))
    if not ok:
        FAILS.append(name)


def control(name, passed):
    print(("[ OK ] control: " if not passed else "[FAIL] control: ") + name +
          (" was caught" if not passed else " was NOT caught"))
    if passed:
        FAILS.append("control: " + name)


class FakeApp:
    """The plugin interface of bin/owngenome_live.py's App, as docs/live-site.md describes it."""

    def __init__(self, root):
        self.config = types.SimpleNamespace(root=root)
        self.routes, self.hooks, self.kinds, self.logs = [], [], [], []

    def route(self, method, pattern, fn, auth=True):
        self.routes.append((method, re.compile(pattern), fn, auth))

    def state_hook(self, name, fn):
        self.hooks.append((name, fn))

    def run_kind(self, fn):
        self.kinds.append(fn)

    def run_file(self, glob):
        raise AssertionError("the family plugin serves its files through its own route")

    def log(self, action, run="-", section="-", status="-"):
        self.logs.append((action, run, section, status))

    def runs(self):
        out = {}
        for fn in self.kinds:
            for r in fn(self):
                out[r["name"]] = r
        return out

    def get(self, path):
        req = FakeReq(self)
        for m, rx, fn, auth in self.routes:
            mt = rx.fullmatch(path)
            if mt and m == "GET":
                fn(req, mt)
                return req.code, req.ctype, req.body, auth
        return 404, None, b"", None


class FakeReq:
    def __init__(self, app):
        self.app, self.code, self.ctype, self.body = app, None, None, b""

    def run_dir(self, name):
        r = self.app.runs().get(name)
        if r:
            return r["dir"]
        d = os.path.join(self.app.config.root, name)
        return d if os.path.isdir(d) else None  # a sample run, as bin/owngenome_live.py lists one

    def send_body(self, code, ctype, body, headers=()):
        self.code, self.ctype, self.body = code, ctype, body.encode() if isinstance(body, str) else body

    def send_json(self, obj, code=200):
        self.send_body(code, "application/json", json.dumps(obj))


def descents(svg):
    out = {}
    for g in ET.fromstring(svg).iter(NS + "g"):
        if g.get("class") == "descent":
            for c in g.get("data-children").split():
                out[c] = (g.get("data-father"), g.get("data-mother"))
    return out


def state_checks(app, root):
    hook = dict(app.hooks)["family"]
    s = hook(os.path.join(root, "family_a"), {}, [])
    res = [("family_a's children with their parents and maps", s.get("children") == CHILDREN, s.get("children")),
           ("family_a's members: seven, gm_m not sequenced",
            [(m["sample"], m["sequenced"]) for m in s.get("members", [])] ==
            [("child_a", True), ("child_b", True), ("father_a", True), ("mother_a", True), ("gf_a", True),
             ("gm_a", True), ("gm_m", False)], s.get("members"))]
    return s, res


def route_checks(app):
    code, ctype, body, auth = app.get("/family/family_a/pedigree.svg")
    ok = code == 200 and ctype == "image/svg+xml" and auth is True
    return [("the pedigree is drawn from the declared parents, behind the login",
             ok and descents(body) == DESCENT, (code, ctype, auth, ok and descents(body)))]


def family_browser_checks(work):
    """Run the demo's family veil group here so its browser time is in another Guard shard."""
    import test_live_demo as demo
    chrome = demo.find_chrome()
    if not chrome:
        if os.environ.get("CI") == "true":
            check("family browser group has Chrome in CI", False)
        else:
            print("[SKIP] no Chrome or Chromium for the family browser group")
        return
    og = demo.og_binary() or demo.build_live_demo.OG
    out = os.path.join(work, "family-browser", "demo")
    demo.build_live_demo.build(out, work=os.path.join(work, "family-browser", "root"), og=og)
    handler = type("H", (demo.Static,), {"__init__": lambda self, *a, **k: demo.SimpleHTTPRequestHandler.__init__(
        self, *a, directory=out, **k)})
    srv = demo.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    demo.threading.Thread(target=srv.serve_forever, daemon=True).start()
    b = demo.Chrome(chrome, os.path.join(work, "family-browser", "chrome"))
    before, t0 = len(demo.FAILS), demo.time.time()
    try:
        demo.gate_checks(b, f"http://127.0.0.1:{srv.server_address[1]}", out, family_only=True)
        FAILS.extend(demo.FAILS[before:])
    finally:
        b.close()
        srv.shutdown()
        srv.server_close()
    print(f"   family browser group: {demo.time.time() - t0:.2f} s")


def main():
    work = tempfile.mkdtemp(prefix="pgp-live-family-", dir=os.environ.get("TMPDIR"))
    try:
        root = os.path.join(work, "genome")
        shutil.copytree(FIX, root, symlinks=True)
        os.makedirs(os.path.join(root, "loner", "logs"))
        app = FakeApp(root)
        owngenome_live_family.register(app)

        # 1. the runs
        runs = app.runs()
        check("family_a and family_b are runs of kind family; child_a is not",
              sorted(runs) == ["family_a", "family_b"] and {r["kind"] for r in runs.values()} == {"family"}, runs)

        # 2. the family's state
        s, res = state_checks(app, root)
        for n, ok, d in res:
            check(n, ok, d)
        check("the relationship check as step 38 wrote it",
              [(r["a"], r["b"], r["verdict"]) for r in s["relationships"]] ==
              [("child_a", "father_a", "ok"), ("child_a", "child_b", "ok")], s["relationships"])
        check("the shared segments as step 38 wrote them",
              s["ibd_pairs"] == [{"a": "child_a", "b": "child_b", "declared": "full siblings", "shared_cm": "2550.0",
                                  "degree": "first (full siblings)"}], s["ibd_pairs"])
        check("the couple screen counted from step 39's table (1 shared, 1 X-linked), the grandparents not screened",
              s["couples"] == [{"a": "father_a", "b": "mother_a",
                                "report": "/family/family_a/couple_father_a_mother_a.md",
                                "table": "/family/family_a/couple_father_a_mother_a.tsv", "screened": True,
                                "shared": 1, "x_linked": 1},
                               {"a": "gf_a", "b": "gm_a", "report": None, "table": None, "screened": False,
                                "shared": 0, "x_linked": 0}], s["couples"])
        check("the plan rows", [r["analysis"] for r in s["plan"]] == ["joint_call", "couple_carrier"], s["plan"])
        check("exactly the result files that exist are listed", s["files"] == FILES, s["files"])
        check("family_b (pedigree.tsv only): a couple not screened yet",
              dict(app.hooks)["family"](os.path.join(root, "family_b"), {}, [])["couples"] ==
              [{"a": "partner_x", "b": "partner_y", "report": None, "table": None, "screened": False, "shared": 0,
                "x_linked": 0}])

        # 3. a sample's page
        hook = dict(app.hooks)["family"]
        check("child_a's page links to family_a", hook(os.path.join(root, "child_a"), {}, []) ==
              {"kind": "member", "families": [{"name": "family_a", "url": "/runs/family_a/"}]})
        check("a sample in no family gets null", hook(os.path.join(root, "loner"), {}, []) is None)

        # 4. the route
        for n, ok, d in route_checks(app):
            check(n, ok, d)
        code, ctype, body, _ = app.get("/family/family_a/couple_father_a_mother_a.md")
        with open(os.path.join(FIX, "family_a", "couple_father_a_mother_a.md"), "rb") as f:
            check("a listed file is served as written", code == 200 and body == f.read() and ctype.startswith("text/plain"),
                  (code, ctype))
        code, ctype, body, _ = app.get("/family/family_a/child_a_segments.svg")
        check("a child's map is served as an SVG", code == 200 and ctype == "image/svg+xml", (code, ctype))
        for path, why in (("/family/family_a/notes.txt", "a file no step writes"),
                          ("/family/family_a/..", "a name of dots"),
                          ("/family/child_a/pedigree.svg", "a sample's folder"),
                          ("/family/nobody/pedigree.svg", "an unknown family")):
            code = app.get(path)[0]
            check("404 for %s" % why, code == 404, code)
        rep = os.path.join(root, "family_a", "ibd_report.md")
        os.replace(rep, rep + ".moved")
        os.symlink(os.path.join(REPO, "README.md"), rep)
        code = app.get("/family/family_a/ibd_report.md")[0]
        listed = "ibd_report.md" in hook(os.path.join(root, "family_a"), {}, [])["files"]
        check("404 for a result file that is a symlink out of the folder, and it is not listed",
              code == 404 and not listed, (code, listed))
        os.remove(rep)
        os.replace(rep + ".moved", rep)
        check("every read is in the access log", ("family", "family_a", "pedigree.svg", 200) in app.logs and
              ("family", "family_a", "notes.txt", 404) in app.logs, app.logs)
        empty = os.path.join(root, "family_empty")
        os.makedirs(empty)
        with open(os.path.join(FIX, "family_a", "members.tsv")) as f, \
                open(os.path.join(empty, "members.tsv"), "w") as g:
            g.write(f.readline())
        es = hook(empty, {}, [])
        code = app.get("/family/family_empty/pedigree.svg")[0]
        check("a members file with a header and no rows: an error on the page, 404 for the pedigree",
              "error" in (es or {}) and code == 404, (es, code))
        got = (hook(os.path.join(root, "loner"), {}, []), hook(os.path.join(root, "child_a"), {}, []))
        check("while a family cannot be read, a sample no readable family lists is unknown, never in no family",
              got == ({"kind": "member", "families": [], "unknown": True},
                      {"kind": "member", "families": [{"name": "family_a", "url": "/runs/family_a/"}]}), got)
        shutil.rmtree(empty)
        check("...and once it reads again, in no family", hook(os.path.join(root, "loner"), {}, []) is None)

        # Planted mistake: child_a's father swapped for its grandfather in the copy.
        bad = os.path.join(work, "swapped")
        shutil.copytree(FIX, bad)
        p = os.path.join(bad, "family_a", "members.tsv")
        with open(p) as f:
            text = f.read()
        with open(p, "w") as f:
            f.write(text.replace("child_a\tfemale\tfather_a\tmother_a", "child_a\tfemale\tgf_a\tmother_a"))
        app2 = FakeApp(bad)
        owngenome_live_family.register(app2)
        control("a swapped parent in the state", state_checks(app2, bad)[1][0][1])
        control("a swapped parent in the drawing", route_checks(app2)[0][1])

        # 5. the page script
        with open(JS) as f:
            js = f.read()
        urls = re.findall(r"(?:https?:)?//[A-Za-z0-9.-]+\.[a-z]{2,}", js)
        check("family.js loads nothing from outside, sets no innerHTML, has no eval",
              not urls and "innerHTML" not in js and not re.search(r"\beval\s*\(|new\s+Function\s*\(", js), urls)
        node = shutil.which("node")
        if node:
            r = subprocess.run([node, "--check", JS], capture_output=True, text=True, timeout=60)
            check("node --check reads family.js", r.returncode == 0, r.stderr)
        elif os.environ.get("CI"):
            check("node is on the CI runner", False)
        else:
            print("[SKIP] no node to syntax-check family.js")

        # 6. the real server
        real_server(root)
        family_browser_checks(work)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    if FAILS:
        print("%d check(s) failed: %s" % (len(FAILS), ", ".join(FAILS)))
        return 1
    print("all checks passed")
    return 0


def real_server(root):
    import owngenome_live  # on main: no skip, so a renamed server fails here instead of passing quietly
    state = os.path.join(os.path.dirname(root), "state")
    live, web, ingest = owngenome_live.start(["--root", root, "--state-dir", state, "--port", "0", "--ingest-port", "0",
                                              "--quiet"])
    threading.Thread(target=web.serve_forever, daemon=True).start()
    port = web.server_address[1]

    def get(path):
        c = http.client.HTTPConnection("127.0.0.1", port, timeout=30)
        c.request("GET", path)
        r = c.getresponse()
        return r.status, r.getheader("Content-Type"), r.read()

    try:
        code, _, body = get("/")
        check("the run list shows family_a as a family", code == 200 and b"family_a" in body and b"family" in body, code)
        code, _, body = get("/runs/family_a/state.json")
        fam = json.loads(body).get("plugins", {}).get("family") if code == 200 else None
        check("state.json carries plugins.family for family_a",
              fam and fam["kind"] == "family" and fam["children"] == CHILDREN, (code, fam))
        code, _, body = get("/runs/child_a/state.json")
        fam = json.loads(body).get("plugins", {}).get("family") if code == 200 else None
        check("state.json of child_a links to family_a", fam and fam["families"][0]["name"] == "family_a", (code, fam))
        code, ctype, body = get("/family/family_a/pedigree.svg")
        check("the real server draws the pedigree", code == 200 and ctype == "image/svg+xml"
              and descents(body) == DESCENT, (code, ctype))
        code, _, _ = get("/family/family_a/notes.txt")
        check("the real server does not serve notes.txt", code == 404, code)
        code, _, body = get("/runs/family_a/")
        check("the family page loads family.js", code == 200 and b"/assets/plugins/family.js" in body, code)
        finding_route_checks(live, get, root)
        browser_checks(port, os.path.dirname(root))
        gate_server_checks(live, get, root, port, os.path.dirname(root))
    finally:
        web.shutdown()
        ingest.server_close()
        web.server_close()


# A summary with one result to act on (BRCA1): what child_a and loner (an adult in no family) both get below.
ACT_SUMMARY = {"sections": {"cpsr": {"state": "ok", "data": {"plp": [
    {"gene": "BRCA1", "variant": "c.68_69del", "final": "Pathogenic", "clinvar": "Pathogenic"}]}},
    "clinvar": {"state": "ok", "data": {"hits": []}}}}


def finding_route_checks(live, get, root):
    """A child with a result to act on: every well-formed finding id answers the same page, so a guessed address
    cannot tell whether the child has that result. An adult keeps the plain 404 for an id that is not a finding."""
    import live_results
    status = os.path.join(root, "loner", "logs", "run_status.tsv")
    with open(status, "w") as f:
        f.write("meta\tstarted_epoch\t1767347940\nstep\t06\tplanned\n")
    real = live.summary
    live.summary = lambda name, d, st: ACT_SUMMARY if name in ("child_a", "loner") else real(name, d, st)
    act = next(f["id"] for f in live_results.build("child_a", ACT_SUMMARY)["findings"] if f["level"] == "act")

    def answers(run):
        # The page names the id it was asked for (data-finding), as the address does: the rest is the same.
        return [(lambda r: (r[0], r[2].replace(x.encode(), b"ID")))(get("/runs/%s/finding/%s/" % (run, x)))
                for x in (act, "two-hits-cftr", "zz-not-a-finding")]
    try:
        child, adult = answers("child_a"), answers("loner")
        check("a child's finding route answers the same page for its act, a guessed act and nonsense",
              [c[0] for c in child] == [200] * 3 and len({c[1] for c in child}) == 1, [c[0] for c in child])
        check("an adult's finding route: the page for its finding, 404 for any other id",
              [c[0] for c in adult] == [200, 404, 404], [c[0] for c in adult])
        # The server loads its plugins itself: patch every loaded copy of the module.
        mods = [m for m in list(sys.modules.values()) if getattr(m, "child_run", None) is not None and hasattr(m, "register")]
        keep = [m.child_run for m in mods]
        for m in mods:
            m.child_run = lambda *a: False
        try:
            control("a server that tells a child's guessed finding id apart",
                    bool(mods) and [c[0] for c in answers("child_a")] == [200] * 3)
        finally:
            for m, k in zip(mods, keep):
                m.child_run = k
    finally:
        live.summary = real
        os.remove(status)


def browser_checks(port, work):
    sys.path.insert(0, os.path.join(REPO, "tests"))
    import test_owngenome_live as tpl
    chrome = tpl.find_chrome()
    if not chrome:
        ci = os.environ.get("CI") == "true"
        print(("[FAIL] " if ci else "[SKIP] ") + "no Chrome or Chromium to load the family page in")
        if ci:
            FAILS.append("no Chrome in CI")
        return
    b = tpl.Chrome(chrome, os.path.join(work, "chrome"))
    try:
        base = "http://127.0.0.1:%d" % port
        b.go(base + "/runs/family_a/")
        card = b.wait("var c = document.querySelector('section[data-plugin=family]'); "
                      "c && !c.hidden && c.querySelector('img.pedigree') && c.textContent")
        check("browser: the family card shows on family_a's page", card and "Couple carrier screen" in card, card)
        loaded = b.wait("var i = document.querySelector('img.pedigree'); i && i.complete && i.naturalWidth > 0 && "
                        "i.naturalWidth", 10)
        check("browser: the pedigree image loads", loaded, loaded)
        # child_a has a run whose results cannot be read: closed whole, chromosome image included. Its parents have
        # no run here, so nothing of theirs is closed: their couple screen shows.
        got = b.wait("var g = document.querySelector('[data-gate-child=child_a]'); "
                     "[!!g, !!document.querySelector('img.segments'), !!document.querySelector('#ch-couples table')].join()", 10)
        check("browser: child_a's unreadable report stays closed whole; its parents' couple screen shows",
              got == "true,false,true", got)
        blocked = b.eval("document.documentElement.dataset.pgpBlocked || ''")
        check("browser: the policy blocked nothing on the family page", not blocked, blocked)
        depth_checks(b, base, os.path.join(work, "genome"))
        b.go(base + "/runs/child_a/")
        closed = b.wait("!!document.querySelector('#report-gate')")
        check("browser: an unreadable child's page shows the fixed gate", bool(closed), closed)

    finally:
        b.close()


def plp(gene, variant):
    return {"gene": gene, "variant": variant, "final": "Pathogenic", "clinvar": "Pathogenic"}


# Two reports with the same two results to act on, held by different cards: two CPSR rows, or one CPSR row and a
# class 5 deletion (which flags the structural-variant card). A parent's page must not tell them apart.
SV = {"state": "ok", "data": {"count": 1, "consensus": 1, "single_caller": 0, "exon_edge_checked": False, "exon_edge": 0,
                              "common": 0, "rows": [{"sv": "chr17:43044295-43125483", "type": "DEL", "gt": "0/1",
                                                     "length": -81188, "acmg_class": 5, "genes": ["BRCA1"]}]}}
TWO_IN_ONE_CARD = {"sections": {"cpsr": {"state": "ok", "data": {"plp": [plp("BRCA1", "c.68_69del"), plp("TP53", "c.743G>A")]}}}}
TWO_CARDS = {"sections": {"cpsr": {"state": "ok", "data": {"plp": [plp("BRCA1", "c.68_69del")]}}, "annotsv": SV}}
NO_ACT = {"sections": {"clinvar": {"state": "ok", "data": {"hits": []}}, "annotsv": dict(SV, data=dict(SV["data"], count=0, rows=[]))}}


def gate_server_checks(live, get, root, port, work):
    """The report gate's server side and its live-mode requests (docstring, check 9)."""
    import owngenome_live
    import live_results
    sys.path.insert(0, os.path.join(REPO, "tests"))
    fams = [{"name": "b", "members": [{"sample": "x", "father": None, "mother": None}]},
            {"name": "a", "members": [{"sample": "x", "father": "p", "mother": None}]},
            {"name": "c", "members": []}]
    child_of = getattr(owngenome_live, "child_of", lambda *a: None)
    got = [child_of(fams, "x"), child_of(fams[:1], "x"), child_of(fams, "y"), child_of(fams[:2], "y")]
    check("one rule of who is a child: the first family by name decides; a run no readable family lists is a child "
          "while a family cannot be read", got == [True, False, True, False], got)
    check("the finding route asks the same rule", __import__("owngenome_live_results").child_run(live, "father_a")
          and not __import__("owngenome_live_results").child_run(live, "gf_a"))
    made = []

    def run(name, summary=None):
        d = os.path.join(root, name)
        os.makedirs(os.path.join(d, "logs"), exist_ok=True)
        if not os.path.exists(os.path.join(d, "logs", "run_status.tsv")):
            with open(os.path.join(d, "logs", "run_status.tsv"), "w") as f:
                f.write("meta\tstarted_epoch\t1767347940\nstep\t06\tplanned\n")
            made.append(d)
        p = os.path.join(d, "summary.json")
        if summary is None:
            if os.path.exists(p):
                os.remove(p)
        else:
            with open(p, "w") as f:
                json.dump(summary, f)

    def family():
        return json.loads(get("/runs/family_a/family.json")[2])

    def index():
        return {r["name"]: r for r in json.loads(get("/runs.json")[2])["runs"]}

    real = live.summary
    try:
        # father_a is a child of gf_a and gm_a; gf_a an adult with a run here.
        run("gf_a")
        # A stale summary.json with no result to act on, while the summary results.json reads (summary(): the step
        # files read again) has one. The home page must count it, so it closes the child as the child's page does.
        run("child_a", NO_ACT)
        live.summary = lambda name, d, st: TWO_IN_ONE_CARD if name == "child_a" else real(name, d, st)
        acts = (index().get("child_a") or {}).get("levels", {}).get("act")
        check("home: a child's results to act on are counted from the summary its results.json reads, not a stale "
              "summary.json", acts == 2, acts)
        check("family.json: no flagged-card count for that child", family()["runs"]["child_a"]["attention"] is None,
              family()["runs"]["child_a"])
        live.summary = real
        run("child_a")
        # The same two results to act on in one card or in two: the flagged-card count differs (0 and 1).
        counts = [live_results.overview(x)["attention"] for x in (TWO_IN_ONE_CARD, TWO_CARDS)]
        assert counts[0] != counts[1], counts
        seen = []
        for summary in (TWO_IN_ONE_CARD, TWO_CARDS):
            run("father_a", summary)
            seen.append(family()["runs"].get("father_a"))
        check("family.json: a child's entry is the same whichever card holds its results to act on, with no count",
              seen[0] == seen[1] and seen[0] and seen[0]["attention"] is None, seen)
        run("father_a", NO_ACT)
        n = family()["runs"]["father_a"]["attention"]
        check("family.json: a child with no result to act on keeps its count", n == live_results.overview(NO_ACT)["attention"], n)
        flags = {k: v.get("child") for k, v in family()["runs"].items()}
        check("family.json: each member run says whether its pages treat it as a child, by the one rule",
              flags.get("father_a") is True and flags.get("gf_a") is False, flags)
        # A no-act child while summary.json is older than a finished step: the home page keeps the numbers it had
        # without the gate (summary.json's), as long as neither summary has a result to act on.
        stale, now = live_results.overview(STALE_SUMMARY), live_results.overview(NO_ACT)
        assert (stale["levels"], stale["top"], stale["lead"]) != (now["levels"], now["top"], now["lead"])
        run("father_a", STALE_SUMMARY)
        live.summary = lambda name, d, st: NO_ACT if name == "father_a" else real(name, d, st)
        got = {k: (index().get("father_a") or {}).get(k) for k in ("attention", "levels", "top", "lead")}
        live.summary = real
        check("home: a no-act child with a stale summary.json keeps the numbers it has without the gate",
              got == {k: stale[k] for k in got}, got)
        run("father_a", NO_ACT)
        family_failure_checks(live, get, root)
        live_browser_checks(live, port, work, run, real)
    finally:
        live.summary = real
        for d in made:
            shutil.rmtree(d)
        run("child_a")


def family_failure_checks(live, get, root):
    """A family file that cannot be read: family.json answers a fixed error, never a dropped connection, and the server
    then reads that family as unreadable, as the page does when family.json fails (one rule of who is a child)."""
    rel = os.path.join(root, "family_a", "relationships.tsv")
    mode = os.stat(rel).st_mode
    os.chmod(rel, 0)
    try:
        try:
            code, _, body = get("/runs/family_a/family.json")
        except (OSError, http.client.HTTPException) as e:
            code, body = type(e).__name__, b""
        check("family.json with an unreadable file: a fixed error", code == 500 and body == b"family not readable\n",
              (code, body[:80]))
        check("...and the server then reads that family as unreadable: an adult it lists is a child until it reads again",
              live.child("gf_a"), live.child("gf_a"))
    finally:
        os.chmod(rel, mode)
    code = get("/runs/family_a/family.json")[0]
    check("...and once the file reads again, family.json answers and the adult is an adult", code == 200 and not live.child("gf_a"),
          (code, live.child("gf_a")))


# The map's per-load timings, the live free-disk figure and the run's elapsed time change between two loads of the
# same code: blanked before two pages are compared.
PAGE_JS = """(function () {
  var clone = document.documentElement.cloneNode(true);
  clone.querySelectorAll('[data-pgp-map-paint-ms],[data-pgp-map-parse-ms]').forEach(function (n) {
    n.removeAttribute('data-pgp-map-paint-ms'); n.removeAttribute('data-pgp-map-parse-ms'); });
  clone.querySelectorAll('.ops-disk, #timing .v, #timing .r').forEach(function (n) { n.textContent = ''; });
  return clone.outerHTML;
}())"""


def differences(x, y, n=4):
    """Up to N places where two pages differ, a little of each side."""
    import difflib
    a, b = x.split("><"), y.split("><")
    out = []
    for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(None, a, b, autojunk=False).get_opcodes():
        if tag != "equal" and len(out) < n:
            out.append(("><".join(a[i1:i2])[:160], "><".join(b[j1:j2])[:160]))
    return out


def slow_routes(live, delays):
    """The server's routes, each answer whose path ends with a key of DELAYS held that many seconds: fixes which of two
    answers arrives first. Returns the routes to put back."""
    import urllib.parse
    keep = list(live.app.routes)

    def wrap(fn):
        def run_(req, m):
            path = urllib.parse.urlsplit(req.path).path
            time.sleep(sum(v for k, v in delays.items() if path.endswith(k)))
            return fn(req, m)
        return run_
    live.app.routes[:] = [(a, rx, wrap(fn), auth) for a, rx, fn, auth in keep]
    return keep


def live_browser_checks(live, port, work, run, real):
    """In live mode (state.json, then the event stream): an adult's and a no-act child's page make the requests of
    the page without the gate, in its order; a parent's page is the same whichever card holds a child's results."""
    import test_live_demo as demo
    chrome = demo.find_chrome()
    if not chrome:
        if os.environ.get("CI") == "true":
            check("live gate browser checks have Chrome in CI", False)
        else:
            print("[SKIP] no Chrome or Chromium for the live gate checks")
        return
    assets = os.path.join(work, "gate-assets")
    shutil.copytree(os.path.join(REPO, "assets", "live"), assets)
    js = os.path.join(assets, "core.js")
    source = open(js).read()
    keep = live.cfg.assets
    live.cfg.assets = assets
    # The map plugin runs og from bin/; CI's og is OG (scripts/setup.sh --og puts one in bin/.og), so the pages
    # compared below have a map.
    # The server loads its plugins itself: the map plugin is the module behind its routes.
    maps = [fn.__globals__ for _, _, fn, _ in live.app.routes if "NO_OG" in getattr(fn, "__globals__", {})][:1]
    og_keep = [m["OG"] for m in maps]
    if os.environ.get("OG") and os.access(os.environ["OG"], os.X_OK):
        for m in maps:
            m["OG"] = os.environ["OG"]
    # gf_a's own report, so its page draws the contents and the family beside them.
    adult = live.summary = lambda name, d, st: NO_ACT if name == "gf_a" else real(name, d, st)
    b = demo.Chrome(chrome, os.path.join(work, "gate-chrome"))
    base = "http://127.0.0.1:%d" % port
    t0 = demo.time.time()
    try:
        b.call("Page.enable", session=b.session)  # scripts for each new document need the Page domain
        log = b.call("Page.addScriptToEvaluateOnNewDocument", {"source": demo.REQUEST_LOG}, b.session)["identifier"]

        def capture(route, text):
            with open(js, "w") as f:
                f.write(text)
            b.go("about:blank")
            assert b.wait("location.href === 'about:blank'")
            b.go(base + route)
            assert b.wait("document.documentElement.dataset.pgpStatus === 'ready' && "
                          "!!document.querySelector('#rail a[href*=family_a]')"), route + " never drawn"
            b.wait("/^(drawn|failed|nodecompress|unavailable)$/.test((document.getElementById('genome-map') || {dataset: {}}).dataset.pgpMap || '')", 15)
            demo.time.sleep(1.0)
            return b.eval("[%s, location.href, (window.ogRequests || []).map(function (r) { "
                          "return r.join(' '); }), (document.getElementById('genome-map') || {dataset: {}}).dataset.pgpMap || '']" % PAGE_JS)
        pages = {}
        for label, route, summary in (("an adult", "/runs/gf_a/", None), ("a no-act child", "/runs/father_a/", NO_ACT)):
            # The child's results come from its summary as results.json reads it (a step file read now).
            live.summary = lambda name, d, st: summary if summary and name == "father_a" else adult(name, d, st)
            try:
                pair = [capture(route, demo.bare_core(source)), capture(route, source)]
            finally:
                live.summary = adult
            pages[label] = pair
            check("live mode: the page of %s makes the requests of the page without the gate, in its order, and is the "
                  "same page" % label, pair[0][:3] == pair[1][:3],
                  [p[2] for p in pair] if pair[0][2] != pair[1][2] else differences(pair[0][0], pair[1][0]))
        drawn = pages["an adult"][1][3]
        if os.environ.get("CI") == "true" or os.environ.get("OG"):
            check("live mode: the adult's page compared has its genome map drawn", drawn == "drawn", drawn)
        else:
            print("[SKIP] no og (OG) here: the live pages are compared without a map (map %r)" % drawn)
        first = [r.split(" ")[1].rsplit("/", 1)[-1] for r in pages["an adult"][0][2]][:3]
        print("   live mode, an adult's reads without the gate: " + ", ".join(pages["an adult"][0][2]))
        check("live mode: without the gate the page asks for state.json, then the event stream, then family.json",
              first == ["state.json", "events", "family.json"], first)
        # A section address on an adult's page: the marks it sets (the contents' aria-current, the card's .target)
        # depend on whether the family or the results answer first, without the gate as with it. Each order is forced
        # by holding the other answer back, twice each; the map's document, asked for once the results are in, comes
        # last, as it does when nothing is held.
        marks = []
        for first_, held in (("results.json", "/family.json"), ("family.json", "/results.json")):
            keep_routes = slow_routes(live, {held: 0.6, "/map.json": 1.5})
            try:
                for route in ("/runs/gf_a/#s3", "/runs/gf_a/#s/clinvar"):
                    for _ in range(2):
                        pair = [capture(route, demo.bare_core(source)), capture(route, source)]
                        marks.append((first_, route, pair[0][0].count('aria-current="true"'), pair[0][:3] == pair[1][:3]))
                        if pair[0][:3] != pair[1][:3]:
                            print("   %s first, %s: %r" % (first_, route, differences(pair[0][0], pair[1][0]) or
                                                          [p[1:3] for p in pair]))
            finally:
                live.app.routes[:] = keep_routes
        check("live mode: an adult's page at a section address is the page without the gate, whichever of the family "
              "and the results answers first", all(m[3] for m in marks), marks)
        check("...and the two orders give the page two outcomes, so the check sees the order",
              len({(m[0], m[2] > 0) for m in marks}) == 2 and len({m[2] > 0 for m in marks}) == 2, marks)
        # A parent's page, whichever card holds the child's two results to act on.
        seen = []
        for summary in (TWO_IN_ONE_CARD, TWO_CARDS):
            run("father_a", summary)
            seen.append(capture("/runs/gf_a/", source)[0])
        check("live mode: a parent's page is the same whichever card holds a child's results to act on",
              seen[0] == seen[1], "the pages differ")
        b.call("Page.removeScriptToEvaluateOnNewDocument", {"identifier": log}, b.session)
        family_gate_checks(b, base, live, run)
    finally:
        live.summary = real
        live.cfg.assets = keep
        for m, o in zip(maps, og_keep):
            m["OG"] = o
        b.close()
        shutil.rmtree(assets, ignore_errors=True)
    print("   live gate browser checks: %.2f s" % (demo.time.time() - t0))


def family_gate_checks(b, base, live, run):
    """The family page: a gated child's member line is the fixed closed line, never "Nothing flagged"; once the child's
    results are opened there, the line and the pedigree's mark give the real count, as without the gate. A person who
    is a child in one family and a partner in another is closed on both family pages (one rule of who is a child)."""
    import live_results
    run("father_a", TWO_CARDS)
    prev = live.summary
    live.summary = lambda name, d, st: TWO_CARDS if name == "father_a" else prev(name, d, st)
    n = live_results.overview(TWO_CARDS)["attention"]
    assert n, "the released child needs a count to show"
    z = os.path.join(live.root, "family_z")
    os.makedirs(z)
    with open(os.path.join(z, "members.tsv"), "w") as f:
        f.write("sample\tsex\tfather\tmother\tpartner\tin_pedigree\tsequenced\n"
                "father_a\tmale\t\t\tzara\tyes\tyes\nzara\tfemale\t\t\tfather_a\tyes\tno\n")
    line = ("Array.prototype.filter.call(document.querySelectorAll('.person-card'), function (c) { "
            "return c.querySelector('.person-name').textContent === 'father_a'; }).map(function (c) { "
            "var l = c.querySelector('.person-look'); return l ? l.textContent : ''; })[0]")
    flag = ("Array.prototype.filter.call(document.querySelectorAll('.ped-person'), function (g) { "
            "return g.querySelector('.ped-name').textContent === 'father_a'; }).map(function (g) { "
            "var f = g.querySelector('.ped-flag'); return f ? f.textContent : ''; })[0]")
    try:
        b.go("about:blank")
        assert b.wait("location.href === 'about:blank'")
        b.go(base + "/runs/family_a/")
        assert b.wait("!!document.querySelector('[data-gate-child=father_a] [data-opening^=closed]')"), "father_a never gated"
        closed = [b.eval(line), b.eval(flag)]
        for _ in range(4):
            if not b.eval("document.querySelectorAll('[data-gate-child=father_a] [data-opening^=closed]').length"):
                break
            b.eval("document.querySelector('[data-gate-child=father_a] [data-opening^=closed]').click(); true")
            b.wait("!document.querySelector('[data-gate-child=father_a]') || "
                   "!!document.querySelector('[data-gate-child=father_a] [data-opening^=opened]')", 10)
        b.wait("!document.querySelector('[data-gate-child=father_a]') && /worth a look|Nothing flagged/.test(%s || '')" % line, 15)
        opened = [b.eval(line), b.eval(flag)]
        check("family page: a gated child's member line is the fixed closed line; once opened, the real count shows on "
              "the line and the pedigree", closed == ["Report closed", ""] and opened == ["%d worth a look" % n, "▲ %d" % n],
              (closed, opened))
        b.go("about:blank")
        assert b.wait("location.href === 'about:blank'")
        b.go(base + "/runs/family_z/")
        got = b.wait("document.documentElement.dataset.pgpStatus === 'ready' && !!document.querySelector('.person-card') && "
                     "(document.querySelector('[data-gate-child=father_a]') ? 'closed' : 'open')", 15)
        check("family page: a person who is a child in another family is closed on this family's page too", got == "closed", got)
    finally:
        live.summary = prev
        shutil.rmtree(z, ignore_errors=True)
        run("father_a", NO_ACT)


# A stale blood-group section with two systems past ABO and RhD, for child_a's summary.json.
STALE_SUMMARY = {"sections": {"bloodtype": {"state": "stale", "note": "older than its input", "data": {
    "abo": "A1", "rhd": "D+", "systems": [{"system": "ABO", "phenotype": "A1"}, {"system": "RhD", "phenotype": "D+"},
                                          {"system": "Kell", "phenotype": "K-k+"}, {"system": "Duffy", "phenotype": "Fy(a+b-)"}]}},
    "haplogroup": {"state": "ok", "data": {"haplogroup": "H1"}}}}


def depth_checks(b, base, root):
    """The couple screen's risk drawing, a child's trio files (one missing, one linking out of the folder) and the
    side by side (a stale value, every blood group system), on the family page in the browser."""
    fam_dir = os.path.join(root, "family_a")
    tsv = os.path.join(fam_dir, "couple_father_a_mother_a.tsv")
    with open(tsv) as f:
        text = f.read()
    text = text.replace("\tshared\t25%\t", "\tshared\t1 in 4 per pregnancy\t")
    text += ("GENE4\tAR\tACMG\tinvented condition four\tshared\t1 in 2 per pregnancy\tfather_a\tinvented finding"
             "\tmother_a\tinvented finding\tinvented\n")
    with open(tsv, "w") as f:
        f.write(text)
    secret = os.path.join(os.path.dirname(root), "outside-folder-name.tsv")
    with open(secret, "w") as f:
        f.write("child\tgene\tvariant_a\tvariant_b\ta_from\tb_from\tphase\n")
    os.symlink(secret, os.path.join(fam_dir, "compound_het.tsv"))
    with open(os.path.join(root, "child_a", "summary.json"), "w") as f:
        json.dump(STALE_SUMMARY, f)
    # Give the formatting fixtures readable no-act results. Unknown results deliberately cannot open.
    log = os.path.join(root, "child_a", "logs", "run_status.tsv")
    with open(log) as f:
        old_log = f.read()
    with open(log, "a") as f:
        f.write("step\t06\tok\n")
    parent_dirs = []
    for parent in ("father_a", "mother_a"):
        d = os.path.join(root, parent)
        os.makedirs(os.path.join(d, "logs"), exist_ok=True)
        with open(os.path.join(d, "logs", "run_status.tsv"), "w") as f:
            f.write("step\t06\tok\n")
        parent_dirs.append(d)
    clin = os.path.join(root, "child_a", "clinical")
    try:
        b.go(base + "/runs/family_a/")
        ready = b.wait("document.documentElement.dataset.pgpStatus === 'ready' && "
                       "!!document.querySelector('div[data-shared=GENE4]')")
        four = b.eval("var k = document.querySelector('div[data-shared=GENE1] .kids4 .grid4'); k ? k.children.length : -1")
        check("browser: a 1-in-4 shared row draws the four children", ready and four == 4, (ready, four))
        two = b.eval("var d = document.querySelector('div[data-shared=GENE4]'); d ? "
                     "[!!d.querySelector('.kids4'), d.textContent].join('|') : ''")
        check("browser: a 1-in-2 shared row gives the chance in words, without the drawing",
              two.startswith("false|") and "1 in 2 per pregnancy" in two, two)
        dn = b.eval("var p = document.querySelector('details[data-child=child_a] p[data-trio=de_novo]'); p && p.textContent")
        ch = b.eval("var p = document.querySelector('details[data-child=child_a] p[data-trio=compound_het]'); "
                    "p && p.textContent")
        page = b.eval("document.body.textContent") or ""
        check("browser: a missing trio file and one linking out of the folder are told apart, the link's target unshown",
              dn and "Not assessed" in dn and "de_novo.tsv is not here" in dn and ch and "link out of the family folder"
              in ch and "outside-folder-name" not in page, (dn, ch))
        rows = b.eval("Array.from(document.querySelectorAll('#ch-side-by-side tbody tr')).map(function (r) { "
                      "return Array.from(r.cells).map(function (c) { return c.textContent; }).join('|'); }).join('\\n')") or ""
        kell = next((r for r in rows.split("\n") if "|Kell|" in r), "")
        hap = next((r for r in rows.split("\n") if "|Mitochondrial|" in r), "")
        legend = b.eval("var p = document.querySelector('#ch-side-by-side p[data-stale]'); p && p.dataset.stale")
        check("browser: side by side has a row per blood group system, a stale value marked in its cell and no other",
              "|Duffy|" in rows and kell.split("|")[2] == "K-k+ (older than its input)" and hap.split("|")[2] == "H1"
              and legend == "4", (kell, hap, legend))
        os.remove(os.path.join(fam_dir, "compound_het.tsv"))
        shutil.copy(secret, os.path.join(fam_dir, "compound_het.tsv"))

        def comp_line():
            b.go(base + "/runs/family_a/")
            b.wait("document.documentElement.dataset.pgpStatus === 'ready' && "
                   "!!document.querySelector('details[data-child=child_a] p[data-trio=compound_het]')")
            return b.eval("var p = document.querySelector('details[data-child=child_a] p[data-trio=compound_het]'); "
                          "p && p.textContent") or ""
        ch = comp_line()
        check("browser: a child whose run has no step 23 candidates is not assessed for compound-het pairs",
              ch.startswith("Not assessed") and "wrote no candidates for child_a" in ch, ch)
        os.makedirs(clin)
        vcf = os.path.join(clin, "child_a_clinical.vcf.gz")
        open(vcf, "w").close()
        ch = comp_line()
        check("browser: candidates newer than compound_het.tsv leave the compound-het line not known",
              ch.startswith("Not known whether assessed") and "newer than compound_het.tsv" in ch, ch)
        old = os.path.getmtime(os.path.join(fam_dir, "compound_het.tsv")) - 60
        os.utime(vcf, (old, old))
        ch = comp_line()
        check("browser: an assessed child's empty compound-het line says pairs come only from step 23's candidates",
              ch.startswith("Assessed (") and "step 23" in ch, ch)
    finally:
        os.remove(os.path.join(fam_dir, "compound_het.tsv"))
        os.remove(os.path.join(root, "child_a", "summary.json"))
        shutil.rmtree(clin, ignore_errors=True)
        with open(log, "w") as f:
            f.write(old_log)
        for d in parent_dirs:
            shutil.rmtree(d)


if __name__ == "__main__":
    sys.exit(main())
