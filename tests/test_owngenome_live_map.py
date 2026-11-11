#!/usr/bin/env python3
"""bin/owngenome_live_map.py in a running bin/owngenome_live.py, with a fake og beside it.

The plugin runs bin/og from its own folder, so the test copies it (and the results plugin) into a
scratch bin/ with a fake og that records its arguments and stdin and answers as told.

Checks:
  1. without a session every route is 401 and og is not run;
  2. GET /runs/<run>/map.json answers og's bytes unchanged (sha256 equal to a direct og call with
     the same arguments and stdin), as application/gzip with the site's headers (CSP, nosniff,
     no-referrer, no-store, DENY);
  3. og gets `map build --root ROOT --run RUN --assets <repo>/assets --cache-dir STATE/map`, and its
     stdin is {"results": the body GET /runs/<run>/results.json sends}; live_results.build is a stub
     whose answer carries its arguments, so a plugin that calls it unlike the results plugin differs;
  4. /map/genes/<chrom>.json and /map/genes/names.json answer og's bytes of `og map genes`;
  5. og exit 3 is 404 with og's reason; exit 127 or no og is 503 with the setup line, and og's
     first stderr line on the server's stderr; any other exit is 500 with og's first stderr line
     on the server's stderr and none of og's bytes sent;
  6. a family's run or an unknown run is 404 and og is not run;
  7. each request writes one access-log row 'map' whose run and section are only the run or '-'
     and 'genes': no finding id, no gene and no chromosome; the server has no POST route for it.
Each check also runs against planted mistakes (mutants of the plugin) and must fail there: a check
that cannot fail is no check.

Run: python3 tests/test_owngenome_live_map.py
"""
import contextlib
import hashlib
import http.client
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
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)
import live_results  # noqa: E402
import owngenome_live  # noqa: E402

FIX = os.path.join(REPO, "tests", "fixtures", "live")
RUN = "DEMO-001"
PASSWORD = "synthetic test password 7"
# A results body with a finding, so a plugin that wrote a finding or a gene to the log is caught.
RESULTS = {"schema_version": 1, "sample": RUN, "chapters": [], "note": None,
           "findings": [{"id": "clinvar-PLANTGENE1-1", "gene": "PLANTGENE1", "title": "planted finding"}]}
FAKE_OG = r'''#!PYTHON
import hashlib, json, os, sys
args = sys.argv[1:]
data = sys.stdin.buffer.read()
with open(os.environ["FAKE_OG_LOG"], "a") as f:
    f.write(json.dumps({"args": args, "stdin": data.decode()}) + "\n")
mode = os.environ.get("FAKE_OG_MODE", "ok")
if mode == "nomap":
    sys.stderr.write("og map build: no map for this run: the planted reason\n")
    sys.exit(3)
if mode == "notbuilt":
    sys.stderr.write("ERROR: og is not built for x86_64; run scripts/setup.sh --og\n")
    sys.exit(127)
if mode == "crash":
    sys.stdout.buffer.write(b"half a document")
    sys.stderr.write("og map build: the planted first line\nand a second line\n")
    sys.exit(1)
# Every byte value, and enough of them to fill a pipe, from the arguments and stdin.
seed = hashlib.sha256(json.dumps(args).encode() + data).digest()
sys.stdout.buffer.write(b"\x1f\x8b" + bytes(range(256)) * 1200 + seed + b"\r\n\x00end")
'''
FAILS = []


def check(desc, ok, detail="", quiet=False):
    """Records one check; QUIET runs (the mutants) print nothing."""
    if not quiet:
        print(("[PASS] " if ok else "[FAIL] ") + desc + ("" if ok else f" -- {str(detail)[:600]}"))
    if not ok:
        FAILS.append(desc)


def request(port, path, headers=None, method="GET"):
    """(status, lower-cased headers, body) of one request to the browser port."""
    c = http.client.HTTPConnection("127.0.0.1", port, timeout=30)
    h = {"Host": f"127.0.0.1:{port}"}
    h.update(headers or {})
    c.request(method, path, headers=h)
    r = c.getresponse()
    data = r.read()
    hdr = {k.lower(): v for k, v in r.getheaders()}
    c.close()
    return r.status, hdr, data


def login(port):
    """The Cookie header of a fresh session."""
    c = http.client.HTTPConnection("127.0.0.1", port, timeout=30)
    c.request("POST", "/login", urllib.parse.urlencode({"password": PASSWORD, "next": "/"}),
              {"Host": f"127.0.0.1:{port}", "Content-Type": "application/x-www-form-urlencoded"})
    r = c.getresponse()
    r.read()
    m = re.match(r"owngenome_live=([^;]+)", r.getheader("Set-Cookie") or "")
    c.close()
    return {"Cookie": f"owngenome_live={m.group(1)}"} if m else {}


def og_calls(work):
    """The fake og's calls since the last read (arguments and stdin), then forgets them."""
    p = os.path.join(work, "og.log")
    calls = [json.loads(x) for x in open(p)] if os.path.exists(p) else []
    if os.path.exists(p):
        os.remove(p)
    return calls


def scenario(plugin_src, quiet=False):
    """Runs every check against a server with PLUGIN_SRC as its map plugin; returns the failures."""
    start = len(FAILS)
    work = os.path.realpath(tempfile.mkdtemp(prefix="pgp-live-map-"))
    srv, err = None, io.StringIO()
    keep_build, keep_env = live_results.build, dict(os.environ)
    try:
        root = os.path.join(work, "genome")
        shutil.copytree(os.path.join(FIX, RUN), os.path.join(root, RUN))
        with open(os.path.join(root, RUN, f"{RUN}_report.html"), "w") as f:  # gives the run files to build()
            f.write('<a href="tables/linked.tsv">a linked table</a>')
        os.makedirs(os.path.join(root, "FAM-A"))  # a family's run, from the run_kind below
        state = os.path.join(work, "state")
        plugins = os.path.join(work, "bin")
        os.makedirs(plugins)
        with open(os.path.join(plugins, "owngenome_live_map.py"), "w") as f:
            f.write(plugin_src)
        shutil.copy(os.path.join(BIN, "owngenome_live_results.py"), plugins)
        og = os.path.join(plugins, "og")
        with open(og, "w") as f:
            f.write(FAKE_OG.replace("PYTHON", sys.executable))
        os.chmod(og, 0o755)
        pw = os.path.join(work, "password")
        with open(pw, "w") as f:
            f.write(owngenome_live.hash_password(PASSWORD, iterations=1000) + "\n")
        os.environ.update(FAKE_OG_LOG=os.path.join(work, "og.log"), FAKE_OG_MODE="ok")
        # The stub's answer depends on every argument, so a plugin that calls build() unlike the
        # results plugin sends og a different body.
        live_results.build = lambda name, summary, files=None: json.loads(json.dumps(
            dict(RESULTS, sample=name, summary=summary, files=sorted(files or []))))

        app = owngenome_live.App()
        app.run_kind(lambda a: [{"name": "FAM-A", "dir": os.path.join(root, "FAM-A"), "kind": "family"}])
        with contextlib.redirect_stderr(err):
            live, web, ingest = owngenome_live.start(
                ["--root", root, "--state-dir", state, "--port", "0", "--ingest-port", "0", "--password-file", pw,
                 "--plugin-dir", plugins, "--quiet"], app)
        srv = (live, web, ingest)
        live.summary = lambda name, d, st: {"planted": True}
        for s in (web, ingest):
            threading.Thread(target=s.serve_forever, daemon=True).start()
        port = web.server_address[1]
        routes = [f"/runs/{RUN}/map.json", "/map/genes/chr7.json", "/map/genes/names.json"]
        log_path = os.path.join(state, "access.log")

        def log_rows():
            """The access log's rows, without the header."""
            return [r.split("\t") for r in open(log_path).read().splitlines()[1:]]

        # 1. the login
        codes = [request(port, r)[0] for r in routes]
        check("1. without a session every map route is 401 and og is not run",
              codes == [401] * 3 and not og_calls(work), codes, quiet)
        auth = login(port)
        n_rows = len(log_rows())

        # 2, 3. the person's document
        code, hd, body = request(port, routes[0], auth)
        calls = og_calls(work)
        _, _, results = request(port, f"/runs/{RUN}/results.json", auth)
        direct = subprocess.run([og] + calls[0]["args"], input=calls[0]["stdin"].encode(), capture_output=True,
                                check=False) \
            if calls else None
        og_calls(work)
        check("2. map.json is og's bytes unchanged (sha256 equal to a direct og call)",
              code == 200 and direct is not None and len(body) > 300000
              and hashlib.sha256(body).hexdigest() == hashlib.sha256(direct.stdout).hexdigest(),
              (code, len(body), direct and len(direct.stdout)), quiet)
        want = {"content-type": "application/gzip", "x-content-type-options": "nosniff",
                "referrer-policy": "no-referrer", "cache-control": "no-store", "x-frame-options": "DENY",
                "content-security-policy": owngenome_live.CSP}
        check("2. map.json has the site's headers", all(hd.get(k) == v for k, v in want.items()), hd, quiet)
        check("3. og runs `map build` with the root, the run, the repo's assets and STATE/map",
              len(calls) == 1 and calls[0]["args"] == ["map", "build", "--root", os.path.realpath(root), "--run", RUN,
                                                       "--assets", os.path.join(work, "assets"),
                                                       "--cache-dir", os.path.join(state, "map")],
              calls, quiet)
        check("3. og's stdin is {\"results\": the body of results.json}",
              len(calls) == 1 and json.loads(calls[0]["stdin"]) == {"results": json.loads(results)}
              and json.loads(results)["findings"] == RESULTS["findings"] and json.loads(results)["files"],
              (calls, results[:200]), quiet)

        # 4. the gene files
        ok = True
        for path, want_args in ((routes[1], ["--chrom", "chr7"]), (routes[2], ["--names"])):
            code, hd, body = request(port, path, auth)
            calls = og_calls(work)
            direct = subprocess.run([og] + calls[0]["args"], input=b"", capture_output=True, check=False) \
                if calls else None
            og_calls(work)
            ok = ok and code == 200 and hd.get("content-type") == "application/gzip" \
                and hd.get("cache-control") == "no-store" and direct is not None and body == direct.stdout \
                and calls[0]["args"] == ["map", "genes", "--root", os.path.realpath(root),
                                         "--cache-dir", os.path.join(state, "map")] + want_args
        check("4. the gene files are og's bytes of `og map genes --chrom` and `--names`", ok, (code, calls), quiet)

        # 5. og's exit codes
        os.environ["FAKE_OG_MODE"] = "nomap"
        code, _, body = request(port, routes[0], auth)
        check("5. og exit 3 is 404 with og's reason", code == 404 and b"the planted reason" in body, (code, body), quiet)
        os.environ["FAKE_OG_MODE"] = "notbuilt"
        with contextlib.redirect_stderr(err):
            code, _, body = request(port, routes[1], auth)
        os.rename(og, og + ".away")
        code2, _, body2 = request(port, routes[2], auth)
        os.rename(og + ".away", og)
        setup = b"The map needs og on this host: run scripts/setup.sh --og"
        check("5. og exit 127, or no og at all, is 503 with the setup line",
              code == 503 and setup in body and code2 == 503 and setup in body2, (code, body, code2, body2), quiet)
        check("5. og exit 127 puts og's first stderr line on the server's stderr",
              "og is not built for x86_64" in err.getvalue(), err.getvalue()[-300:], quiet)
        os.environ["FAKE_OG_MODE"] = "crash"
        with contextlib.redirect_stderr(err):
            code, _, body = request(port, routes[0], auth)
        check("5. any other og exit is 500, logs og's first stderr line, and sends none of og's bytes",
              code == 500 and b"half a document" not in body and "the planted first line" in err.getvalue()
              and "and a second line" not in err.getvalue(), (code, body, err.getvalue()[-300:]), quiet)
        os.environ["FAKE_OG_MODE"] = "ok"
        og_calls(work)

        # 6. runs with no person's map
        codes = [request(port, f"/runs/{r}/map.json", auth)[0] for r in ("FAM-A", "NO-SUCH-RUN")]
        check("6. a family's run and an unknown run are 404, and og is not run",
              codes == [404, 404] and not og_calls(work), codes, quiet)

        # 7. the access log
        code, _, _ = request(port, routes[0], auth, method="POST")
        rows = log_rows()[n_rows:]
        rows = [r for r in rows if r[3] not in ("results",)]
        maps = [r for r in rows if r[3] == "map"]
        text = open(log_path).read()
        check("7. one 'map' row per request, naming only the run or 'genes'",
              len(maps) == 9 and all((r[4], r[5]) in ((RUN, "-"), ("-", "genes"), ("FAM-A", "-"),
                                                         ("NO-SUCH-RUN", "-")) for r in maps),
              (len(maps), maps), quiet)
        check("7. the access log names no finding, no gene and no chromosome",
              not any(w in text for w in ("PLANTGENE1", "clinvar-PLANTGENE1-1", "chr7", "planted")), text, quiet)
        check("7. there is no POST route to the map", code == 404, code, quiet)
    finally:
        live_results.build = keep_build
        os.environ.clear()
        os.environ.update(keep_env)
        if srv:
            live, web, ingest = srv
            live.stopping = True
            for s in (web, ingest):
                s.shutdown()
                s.server_close()
        shutil.rmtree(work, ignore_errors=True)
    return len(FAILS) - start


# Planted mistakes: each must make a check fail.
MUTANTS = [
    ("a route without the login", 'app.route("GET", r"/runs/(%s)/map\\.json" % NAME, build)',
     'app.route("GET", r"/runs/(%s)/map\\.json" % NAME, build, auth=False)'),
    ("the log row names a finding", 'stdin, run["name"], "-")',
     'stdin, run["name"], body["findings"][0]["id"])'),
    ("the log row names the chromosome", 'b"", "-",\n                     "genes")',
     'b"", "-",\n                     m.group(1))'),
    ("og's bytes changed", "return req.send_body(200, GZIP, p.stdout)",
     "return req.send_body(200, GZIP, p.stdout[:-1])"),
    ("exit 3 read as a failure", "if code == 3:", "if code == 33:"),
    ("og missing read as a failure", "if code == 127:", "if code == 172:"),
    ("og's stdout sent on a failure", 'return req.send_body(500, TEXT, "The map could not be built: the server\'s log '
     'says why.\\n")\n\n\ndef register', "return req.send_body(500, TEXT, p.stdout)\n\n\ndef register"),
    ("stdin without the results", 'json.dumps({"results": body})', 'json.dumps({"results": {}})'),
    ("og's 127 reason dropped", 'if err:  # e.g.', 'if False:  # e.g.'),
    ("build() given other files", 'live.run_files(run) if run["dir"] else set())', 'set())'),
    ("build() given another name", 'return live_results.build(name, summary,',
     'return live_results.build(name + "x", summary,'),
    ("a family's run built", 'if not run or run.get("kind") == "family":', "if not run:"),
    ("the assets of the wrong folder", 'ASSETS = os.path.join(os.path.dirname(HERE), "assets")',
     'ASSETS = os.path.join(HERE, "assets")'),
]


def main():
    """Runs the checks on the plugin, then on each planted mistake."""
    with open(os.path.join(BIN, "owngenome_live_map.py")) as f:
        src = f.read()
    scenario(src)
    real = len(FAILS)
    for name, old, new in MUTANTS:
        if old not in src:
            check(f"mutant '{name}': its line is in owngenome_live_map.py", False, old)
            continue
        keep = len(FAILS)
        n = scenario(src.replace(old, new, 1), quiet=True)
        del FAILS[keep:]
        check(f"mutant '{name}' fails a check ({n} failed)", n > 0)
    print(f"\n{'FAILED: ' + str(len(FAILS)) + ' check(s)' if FAILS else 'All checks passed.'}")
    return 1 if FAILS or real else 0


if __name__ == "__main__":
    sys.exit(main())
