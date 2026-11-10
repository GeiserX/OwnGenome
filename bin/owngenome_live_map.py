"""owngenome_live_map.py: the genome map's files on the live site, a plugin of bin/owngenome_live.py.

  GET /runs/<run>/map.json        `og map build` of a person's run: stdin is {"results": the body
                                  /runs/<run>/results.json sends}, the answer is og's stdout (gzip)
  GET /map/genes/<chrom>.json     `og map genes --chrom <chrom>`: one chromosome's genes (gzip)
  GET /map/genes/names.json       `og map genes --names`: every gene's name and place (gzip)

Every route needs the login. The plugin decides nothing about the map: og builds each file and
caches it under <state-dir>/map, and the plugin sends og's bytes unchanged, with the site's
headers. og writes a file only once it is complete and its exit code sets the status, so the
plugin waits for og to end before it answers: 0 is 200; 3 (no map for this run, no gene files)
is 404 with og's reason; 127 (og not built on this host) is 503; anything else is 500, with og's
first stderr line on the server's stderr. There is no POST route: opening a closed result is the
browser's (the document tags each row a closed result can hide). Each request writes one access-log
row 'map', with the run (or 'genes' as the section for a gene file) and never a finding, a gene or a
chromosome.

This is Python only while the site is: the Rust port of the live server replaces this file.
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OG = os.path.join(HERE, "og")
ASSETS = os.path.join(os.path.dirname(HERE), "assets")
NAME = r"[A-Za-z0-9][A-Za-z0-9._-]*"
OG_SECONDS = 300
GZIP = "application/gzip"
TEXT = "text/plain; charset=utf-8"
NO_OG = "The map needs og on this host: run scripts/setup.sh --og\n"


def results_body(live, run):
    """What GET /runs/<run>/results.json sends (bin/owngenome_live_results.py), built the same way."""
    import live_results
    name = run["name"]
    d = run["dir"] or os.path.join(live.root, name)
    summary = live.summary(name, d, live.base_state(run))
    if not summary:
        return {"schema_version": 1, "sample": name, "chapters": [],
                "note": "No step has finished yet: the results appear here as the steps end."}
    return live_results.build(name, summary, live.run_files(run) if run["dir"] else set())


def relay(app, req, args, stdin, run, section):
    """Runs og with ARGS and STDIN and sends its answer; one access-log row."""
    try:
        p = subprocess.run([OG, "map"] + args, input=stdin, capture_output=True, timeout=OG_SECONDS)
        code = p.returncode
    except FileNotFoundError:
        p, code = None, 127
    except (OSError, subprocess.TimeoutExpired) as e:
        p, code = None, type(e).__name__
    err = (p.stderr.decode(errors="replace").strip().splitlines() or [""])[0] if p else ""
    if code == 0:
        app.log("map", run, section, 200)
        return req.send_body(200, GZIP, p.stdout)
    if code == 3:
        app.log("map", run, section, 404)
        return req.send_body(404, TEXT, err + "\n")
    if code == 127:
        app.log("map", run, section, 503)
        return req.send_body(503, TEXT, NO_OG)
    print(f"owngenome_live: og map {args[0]} failed ({code}): {err}", file=sys.stderr)
    app.log("map", run, section, 500)
    return req.send_body(500, TEXT, "The map could not be built: the server's log says why.\n")


def register(app):
    def cache_dir():
        return os.path.join(app.config.state_dir, "map")

    def build(req, m):
        live = req.live
        run = live.run(m.group(1))
        if not run or run.get("kind") == "family":
            app.log("map", m.group(1), "-", 404)
            return req.send_body(404, TEXT, "not found\n")
        try:
            body = results_body(live, run)
        except Exception as e:  # a reader bug: no map, never a page down
            print(f"owngenome_live: results of {run['name']} not built for the map: {type(e).__name__}",
                  file=sys.stderr)
            app.log("map", run["name"], "-", 500)
            return req.send_body(500, TEXT, "The map could not be built: the server's log says why.\n")
        stdin = json.dumps({"results": body}).encode()
        return relay(app, req, ["build", "--root", live.root, "--run", run["name"], "--assets", ASSETS,
                                "--cache-dir", cache_dir()], stdin, run["name"], "-")

    def genes(req, m):
        which = ["--names"] if m.group(1) == "names" else ["--chrom", m.group(1)]
        return relay(app, req, ["genes", "--root", req.live.root, "--cache-dir", cache_dir()] + which, b"", "-",
                     "genes")

    app.route("GET", r"/runs/(%s)/map\.json" % NAME, build)
    app.route("GET", r"/map/genes/(chr[0-9A-Za-z]{1,2}|names)\.json", genes)
