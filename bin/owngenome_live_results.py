"""owngenome_live_results.py: a person's results on the live site, a plugin of bin/owngenome_live.py.

  GET /runs/<run>/results.json   bin/live_results.py's build() of the run's summary: the
                                 report's seven chapters, one card per section with its
                                 facts, its rows, its caveats and its attention flag, and
                                 the Journal's findings, levels, lead and ten sections
  GET /runs/<run>/finding/<id>/  one finding's page (assets/live/finding.html): 404 unless
                                 <id> is a finding of that run's results.json now; for a
                                 child (child_run) with a result to act on, the page for
                                 any well-formed id

The summary is the one the server already reads for the step headlines
(Live.summary): the run's summary.json when step 24 wrote it after every done
step, else collect_summary.collect() of the folder, read again only when a step
ends. A reader error gives no cards and a note, never a 500. The raw
summary.json is never served: no run_file glob is added, and results.json holds
what the report prints. Every request writes one access-log row.
"""
import html
import os
import sys
import urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

NAME = r"[A-Za-z0-9][A-Za-z0-9._-]*"


def register(app):
    def results(req, m):
        live = req.live
        run = live.run(m.group(1))
        if not run or run.get("kind") == "family":
            app.log("results", m.group(1), "-", 404)
            return req.send_body(404, "text/plain; charset=utf-8", "not found\n")
        name = run["name"]
        d = run["dir"] or os.path.join(live.root, name)
        out = {"schema_version": 1, "sample": name, "chapters": [], "note": None}
        try:
            import live_results
            st = live.base_state(run)
            summary = live.summary(name, d, st)
            if summary:
                files = live.run_files(run) if run["dir"] else set()
                out = live_results.build(name, summary, files)
            else:
                out["note"] = "No step has finished yet: the results appear here as the steps end."
        except Exception as e:  # a reader bug must not take the page down
            print(f"owngenome_live: results of {name} not built: {type(e).__name__}", file=sys.stderr)
            out = {"schema_version": 1, "sample": name, "chapters": [],
                   "note": f"The results could not be read ({type(e).__name__}); the report and the Progress tab "
                           "still work."}
        app.log("results", name, "-", 200)
        return req.send_json(out)

    app.route("GET", r"/runs/(%s)/results\.json" % NAME, results)

    def finding(req, m):
        """One finding's page. The page itself reads ../../results.json; this only checks that the id is
        a finding of the run now, so a typed or stale address is a plain 404, never an empty page. A child with a
        result to act on gets the page for any well-formed id: the page keeps the child's report closed
        (assets/live/core.js) until every such result is opened, and a 404 would tell a guessed id apart from a
        result the child has."""
        live = req.live
        # A page: without a session it sends the browser to the login and back, as the viewer does.
        if live.password_hash and req.sid is None:
            app.log("finding", m.group(1), "-", 303)
            return req.send_body(303, "text/plain", "", [("Location", "/login?next=" + urllib.parse.quote(
                urllib.parse.urlsplit(req.path).path))])
        run = live.run(m.group(1))
        fid = m.group(2)
        ok = False
        if run and run.get("kind") != "family" and live_results_mod().FINDING_ID.match(fid):
            try:
                d = run["dir"] or os.path.join(live.root, run["name"])
                summary = live.summary(run["name"], d, live.base_state(run))
                found = live_results_mod().build(run["name"], summary)["findings"] if summary else []
                ok = any(f["id"] == fid for f in found) or (
                    any(f["level"] == "act" for f in found) and child_run(live, run["name"]))
            except Exception as e:  # a reader bug is a 404 here, never a 500
                print(f"owngenome_live: finding of {run['name']} not checked: {type(e).__name__}", file=sys.stderr)
        if not ok:
            app.log("finding", m.group(1), "-", 404)
            return req.send_body(404, "text/plain; charset=utf-8", "not found\n")
        with open(os.path.join(live.cfg.assets, "finding.html")) as f:
            page = f.read()
        page = (page.replace("{{assets}}", "/assets").replace("{{run}}", html.escape(run["name"]))
                .replace("{{finding}}", html.escape(fid)).replace("{{mode}}", "live"))
        app.log("finding", run["name"], fid, 200)
        return req.send_body(200, "text/html; charset=utf-8", page)

    # auth=False: finding() asks for the login itself, so a bookmarked finding leads to the login form.
    app.route("GET", r"/runs/(%s)/finding/([^/]{1,200})/" % NAME, finding, auth=False)


def child_run(live, name):
    """Whether the person page treats NAME as a child's run: a family lists it with a father or a mother, or a
    family's members cannot be read and no readable family lists NAME."""
    fams = [r for r in live.runs().values() if r.get("kind") == "family"]
    listed = False
    for f in fams:
        for m in f.get("members") or []:
            if m.get("sample") == name:
                listed = True
                if m.get("father") or m.get("mother"):
                    return True
    return not listed and any(not (f.get("members") or []) for f in fams)


def live_results_mod():
    import live_results
    return live_results

