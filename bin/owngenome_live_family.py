"""owngenome_live_family.py: the family page of the live site (a bin/owngenome_live.py plugin).

A family folder under GENOME_DIR (one holding members.tsv or pedigree.tsv, as
steps 38 and 39 and the Nextflow family workflow leave it) is listed as a run
of kind "family". Its page shows the pedigree, drawn from that file at request
time by `og live family svg`, each child's chromosome map
(<child>_segments.svg), the relationship check, the shared segments and the
couple carrier screen, as the steps wrote them. A sample's page links to the
family it is in.

  GET /family/<family>/<file>   pedigree.svg (drawn), or one of the family's
                                result files listed in state.json; any other
                                name is 404
  GET /runs/<family>/family.json  the family page's data (assets/live/family.html):
                                the state below, the pedigree's layout (so each
                                person is a link), the relationship and shared-
                                segment rows with their measured numbers, each
                                couple's gene rows, each child's inherited
                                stretches (<child>_segments.tsv), the de novo and
                                compound-het candidates (de_novo.tsv,
                                compound_het.tsv) with whether each child was
                                assessed and why not (trio_status; a file linking
                                out of the folder is "rejected" in trio_file_state),
                                each member's run progress, attention count and
                                whether their pages treat them as a child (child),
                                and the members side by side: CYP2D6, the CPIC
                                genes, HLA, the haplogroups and every blood group
                                system, each value with its stale flag (from each
                                summary.json only); 500 with a fixed line when a
                                file of the family cannot be read, and the family
                                then lists no members until it can

state.json gets plugins.family: {"kind": "family", ...} on a family's page,
{"kind": "member", "families": [...]} on the page of a sample in a family,
{"kind": "member", "families": [], "unknown": true} when no readable family
lists the sample but a family's pedigree could not be read, null otherwise.

og builds every answer (`og live family`, rust/crates/og-live/src/family.rs):
state.json's part, family.json, the family list and the drawing. og reads the
files the steps wrote, counts their rows and draws the pedigree the file
declares; it computes nothing about genetics. The plugin adds what only the
site knows (each member run's progress, attention count and child flag) and
sends og's answer. OG names og when set, else bin/og. While og cannot list the
families, every family folder counts as unreadable, so every run is treated
as a child (fail closed). This is Python only while the site is: the Rust
port of the live server replaces this file.
"""
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# OG names og when set (tests/og_binary.py reads it the same way), else bin/og.
OG = os.environ.get("OG") or os.path.join(HERE, "og")
OG_SECONDS = 120

NAME = r"[A-Za-z0-9][A-Za-z0-9._-]*"
NAME_RE = re.compile("^%s$" % NAME)
MARKERS = ("members.tsv", "pedigree.tsv")
TYPES = {".svg": "image/svg+xml", ".md": "text/plain; charset=utf-8", ".tsv": "text/plain; charset=utf-8"}


def inside(path, top):
    real, top = os.path.realpath(path), os.path.realpath(top)
    return real == top or real.startswith(top.rstrip(os.sep) + os.sep)


def is_family(d):
    return bool(d) and os.path.isdir(d) and any(os.path.isfile(os.path.join(d, m)) for m in MARKERS)


def family_dirs(root):
    """[(name, dir)] of the family folders under root, sorted."""
    try:
        names = sorted(os.listdir(root))
    except OSError:
        return []
    out = []
    for n in names:
        d = os.path.join(root, n)
        if not n.startswith(".") and NAME_RE.match(n) and inside(d, root) and is_family(d):
            out.append((n, d))
    return out


def og(args, stdin=None, og_bin=None):
    """og live family ARGS: (exit code, stdout, its first stderr line); 127 when og is not on this host."""
    try:
        p = subprocess.run([og_bin or OG, "live", "family"] + args, input=stdin, capture_output=True,
                           timeout=OG_SECONDS, check=False)
    except FileNotFoundError:
        return 127, b"", "og is not on this host: run scripts/setup.sh --og"
    except (OSError, subprocess.TimeoutExpired) as e:
        return 1, b"", type(e).__name__
    return p.returncode, p.stdout, (p.stderr.decode(errors="replace").strip().splitlines() or [""])[0]


def og_json(args, stdin=None, og_bin=None):
    code, out, err = og(args, stdin, og_bin)
    if code != 0:
        raise OSError(err or "og live family %s exited %s" % (args[0], code))
    return json.loads(out)


def family_state(name, d):
    """state.json's plugins.family on the family's page (og live family state)."""
    return og_json(["state", d, "--name", name])


def family_page(name, d, site=None, root=None, og_bin=None):
    """family.json (og live family page). SITE holds what only the site knows: "runs" (each member run's entry),
    "candidates" and "summaries" (each child's and member's run folder, read only inside ROOT)."""
    args = ["page", d, "--name", name]
    if site is None:
        return og_json(args, None, og_bin)
    args += ["--root", root or os.path.dirname(os.path.abspath(d)), "--input", "-"]
    return og_json(args, json.dumps(site).encode(), og_bin)


def draw(d, title, og_bin=None):
    """The pedigree's SVG (og live family svg), or None when the folder holds no pedigree that names anyone."""
    code, out, err = og(["svg", d, "--title", title], None, og_bin)
    if code == 2:
        return None
    if code != 0:
        raise OSError(err or "og live family svg exited %s" % code)
    return out


def families(root):
    """Each family folder under root as /runs.json lists it: {"name", "dir", "kind", "members", "unreadable"}.
    members is [] while the pedigree cannot be read; unreadable while family.json cannot be built from the folder."""
    if not family_dirs(root):
        return []
    try:
        return og_json(["families", root])
    except (OSError, ValueError) as e:
        print(f"owngenome_live: the families could not be read: {e}", file=sys.stderr)
        return [{"name": n, "dir": d, "kind": "family", "members": [], "unreadable": True} for n, d in family_dirs(root)]


def member_state(root, sample):
    """The families that list the sample, and child: whether its pages treat it as a child, by the one rule over
    every family (bin/owngenome_live.py child_of). A family that cannot be read may be the one that lists it as a
    child: with no readable family listing it, the answer is unknown, never "in no family" ({"kind": "member",
    "families": [], "unknown": true}), and its pages keep results to act on closed, as the home page does."""
    from owngenome_live import child_of
    fams = families(root)
    child = child_of(fams, sample)
    mine = [{"name": f["name"], "url": "/runs/%s/" % f["name"]} for f in fams
            if any(m["sample"] == sample for m in f["members"])]
    if mine:
        return {"kind": "member", "families": mine, "child": child}
    return {"kind": "member", "families": [], "unknown": True} if child else None


def register(app):
    def root():
        return os.path.realpath(app.config.root)

    def runs(app_):
        return families(root())

    def page_data(req, m):
        name = m.group(1)
        d = req.run_dir(name)
        if not d or not is_family(d) or not inside(d, root()):
            app.log("family", name, "family.json", 404)
            return req.send_body(404, "text/plain; charset=utf-8", "not found\n")
        try:
            out = page_json(req, name, d)
        except Exception as e:  # an unreadable file is a fixed error, never a dropped connection
            print(f"owngenome_live: family.json of {name} not built: {type(e).__name__}", file=sys.stderr)
            app.log("family", name, "family.json", 500)
            return req.send_body(500, "text/plain; charset=utf-8", "family not readable\n")
        app.log("family", name, "family.json", 200)
        return req.send_json(out)

    def page_json(req, name, d):
        live = getattr(req, "live", None)
        runs, candidates, summaries = {}, {}, {}
        listed = live.runs() if live else {}
        for p in family_state(name, d).get("members") or []:
            r = listed.get(p["sample"]) if live else None
            if not r:
                continue
            rd = r["dir"] or os.path.join(live.root, r["name"])
            summaries[p["sample"]] = rd
            if r.get("kind") == "family":
                continue
            candidates[p["sample"]] = rd
            st = live.base_state(r)
            planned = [s for s in st["steps"] if s["state"] != "not planned"]
            n, mtime = live.attention(rd)
            # A child whose own page keeps the report closed has no count of flagged cards here: the count changes
            # with which cards the unopened results flag, and a parent's page would show it.
            if live.gated(p["sample"], rd, st, listed):
                n = None
            # child: whether the person's pages treat them as a child, by the one rule over every family (child_of),
            # so the family page closes the same people as their own pages do.
            runs[p["sample"]] = {"url": "/runs/%s/" % p["sample"], "state": st["run"]["state"],
                                 "planned": len(planned), "done": sum(1 for s in planned if s["state"] == "done"),
                                 "attention": n, "results": mtime is not None, "child": live.child(p["sample"], listed)}
        # og tells each child's compound-het line by its run's step 23 candidates, and puts the members side by
        # side from each sequenced member's summary.json alone (as the home page's attention count reads it:
        # never a file linking out of its run).
        return family_page(name, d, {"runs": runs, "candidates": candidates, "summaries": summaries},
                           live.root if live else None)

    def hook(run_dir, state, events):
        if not run_dir:
            return None
        name = os.path.basename(os.path.normpath(run_dir))
        if is_family(run_dir) and inside(run_dir, root()):
            return family_state(name, run_dir)
        return member_state(root(), name)

    def serve(req, m):
        name, f = m.group(1), m.group(2)
        d = req.run_dir(name)
        if not d or not is_family(d) or not inside(d, root()):
            app.log("family", name, f, 404)
            return req.send_body(404, "text/plain; charset=utf-8", "not found\n")
        if f == "pedigree.svg":
            try:
                svg = draw(d, name)
            except OSError as e:
                print(f"owngenome_live: pedigree.svg of {name} not drawn: {e}", file=sys.stderr)
                app.log("family", name, f, 503)
                return req.send_body(503, "text/plain; charset=utf-8", "The family page needs og on this host\n")
            if svg is None:
                app.log("family", name, f, 404)
                return req.send_body(404, "text/plain; charset=utf-8", "not found\n")
            app.log("family", name, f, 200)
            return req.send_body(200, TYPES[".svg"], svg)
        code, out, _ = og(["files", d])
        files = json.loads(out) if code == 0 else []
        p = os.path.join(d, f)
        if f not in files or not inside(p, d):
            app.log("family", name, f, 404)
            return req.send_body(404, "text/plain; charset=utf-8", "not found\n")
        with open(p, "rb") as fh:
            body = fh.read()
        app.log("family", name, f, 200)
        return req.send_body(200, TYPES.get(os.path.splitext(f)[1], "text/plain; charset=utf-8"), body)

    app.run_kind(runs)
    app.state_hook("family", hook)
    app.route("GET", r"/family/(%s)/(%s)" % (NAME, NAME), serve)
    app.route("GET", r"/runs/(%s)/family\.json" % NAME, page_data)
