"""owngenome_live_results.py: a person's results on the live site, a plugin of bin/owngenome_live.py.

  GET /runs/<run>/results.json   bin/live_results.py's build() of the run's summary: the
                                 report's seven chapters, one card per section with its
                                 facts, its rows, its caveats and its attention flag

The summary is the one the server already reads for the step headlines
(Live.summary): the run's summary.json when step 24 wrote it after every done
step, else collect_summary.collect() of the folder, read again only when a step
ends. A reader error gives no cards and a note, never a 500. The raw
summary.json is never served: no run_file glob is added, and results.json holds
what the report prints. Every request writes one access-log row.
"""
import os
import sys

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

