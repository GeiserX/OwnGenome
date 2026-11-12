"""owngenome_live_family.py: the family page of the live site (a bin/owngenome_live.py plugin).

A family folder under GENOME_DIR (one holding members.tsv or pedigree.tsv, as
steps 38 and 39 and the Nextflow family workflow leave it) is listed as a run
of kind "family". Its page shows the pedigree, drawn from that file at request
time by bin/pedigree_svg.py, each child's chromosome map
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
                                each member's run progress and attention count,
                                and the members side by side: CYP2D6, the CPIC
                                genes, HLA, the haplogroups and every blood group
                                system, each value with its stale flag (from each
                                summary.json only)

state.json gets plugins.family: {"kind": "family", ...} on a family's page,
{"kind": "member", "families": [...]} on the page of a sample in a family,
{"kind": "member", "families": [], "unknown": true} when no readable family
lists the sample but a family's pedigree could not be read, null otherwise.

The plugin computes nothing about genetics: it reads the files the steps
wrote, counts their rows and draws the pedigree the file declares.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pedigree_svg  # noqa: E402

NAME = r"[A-Za-z0-9][A-Za-z0-9._-]*"
NAME_RE = re.compile("^%s$" % NAME)
MARKERS = ("members.tsv", "pedigree.tsv")
# The family's result files the page may open, besides the drawn pedigree and the
# per-child and per-couple files (named from the members).
FIXED = ("family_report.md", "ibd_report.md", "relationships.tsv", "ibd_pairs.tsv", "family_plan.tsv",
         "members.tsv", "pedigree.tsv")
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


def read_tsv(path):
    """Rows of a tab-separated file with a header line, as dicts; [] when it is not there or is a
    link out of its own folder (family.json serves only what the family folder itself holds)."""
    if not os.path.isfile(path) or not inside(path, os.path.dirname(os.path.abspath(path))):
        return []
    rows = []
    with open(path, errors="replace") as f:
        head = f.readline().rstrip("\r\n").lstrip("#").split("\t")
        for line in f:
            line = line.rstrip("\r\n")
            if line.strip() and not line.startswith("#"):
                vals = line.split("\t")
                rows.append(dict(zip(head, vals + [""] * (len(head) - len(vals)))))
    return rows


def members_of(d):
    order, members, _ = pedigree_svg.load(d)
    if not order:
        raise ValueError("the pedigree names no one")
    return order, members


def couple_names(order, members):
    """couple_<a>_<b> of each pair of partners, a and b sorted by name as step 39 names them."""
    out = []
    for a, b in pedigree_svg.couples(order, members):
        if members[a]["partner"] == b or members[b]["partner"] == a:
            x, y = sorted((a, b))
            out.append((x, y, "couple_%s_%s" % (x, y)))
    return out


def served_files(d, order, members):
    """The file names of the family folder the page may open: those that exist."""
    names = list(FIXED) + ["%s_segments.svg" % s for s in order]
    for _, _, base in couple_names(order, members):
        names += [base + ".md", base + ".tsv"]
    return sorted(n for n in set(names) if os.path.isfile(os.path.join(d, n)) and inside(os.path.join(d, n), d))


def family_state(name, d):
    try:
        order, members = members_of(d)
    except (OSError, ValueError) as e:
        return {"kind": "family", "error": "the pedigree could not be read: %s" % e}
    files = served_files(d, order, members)
    url = "/family/%s/" % name

    def link(f):
        return url + f if f in files else None

    people = [{"sample": s, "sex": members[s]["sex"], "father": members[s]["father"],
               "mother": members[s]["mother"], "partner": members[s]["partner"],
               "sequenced": members[s]["sequenced"]} for s in order]
    children = [{"child": p["sample"], "father": p["father"], "mother": p["mother"],
                 "segments": link("%s_segments.svg" % p["sample"])}
                for p in people if p["father"] or p["mother"]]
    rel = [{"a": r.get("sample_a"), "b": r.get("sample_b"), "declared": r.get("declared"),
            "verdict": r.get("verdict"), "reason": r.get("reason")}
           for r in read_tsv(os.path.join(d, "relationships.tsv"))]
    ibd = [{"a": r.get("sample_a"), "b": r.get("sample_b"), "declared": r.get("declared"),
            "shared_cm": r.get("shared_cm"), "degree": r.get("degree")}
           for r in read_tsv(os.path.join(d, "ibd_pairs.tsv"))]
    couples = []
    for a, b, base in couple_names(order, members):
        rows = read_tsv(os.path.join(d, base + ".tsv"))
        results = [r.get("result") for r in rows]
        couples.append({"a": a, "b": b, "report": link(base + ".md"), "table": link(base + ".tsv"),
                        "screened": bool(rows) or base + ".md" in files,
                        "shared": results.count("shared"), "x_linked": results.count("x_linked")})
    plan = [{"analysis": r.get("analysis"), "members": r.get("members"), "decision": r.get("decision"),
             "reason": r.get("reason")} for r in read_tsv(os.path.join(d, "family_plan.tsv"))]
    return {"kind": "family", "name": name, "pedigree": url + "pedigree.svg", "members": people,
            "children": children, "relationships": rel, "ibd_pairs": ibd, "couples": couples, "plan": plan,
            "reports": [{"name": f, "url": url + f} for f in files if f.endswith(".md")],
            "files": files}


def layout_of(order, members):
    """The pedigree drawing's geometry (bin/pedigree_svg.py layout) as data, so the page draws each
    person as a link: people with their centre, partner lines and the descent lines of each sibship."""
    pos, gen, pairs, widest, ngen = pedigree_svg.layout(order, members)
    h = pedigree_svg.SIZE / 2
    people = [{"sample": s, "x": round(pos[s][0], 1), "y": round(pos[s][1], 1), "gen": gen[s], "sex": members[s]["sex"],
               "sequenced": bool(members[s]["sequenced"])} for s in order]
    lines = []
    for a, b in pairs:
        (xa, ya), (xb, yb) = sorted((pos[a], pos[b]))
        lines.append([xa + h, ya, xb - h, yb])
    sibships = {}
    for s in order:
        f, mo = members[s]["father"], members[s]["mother"]
        f, mo = (f if f in pos else None), (mo if mo in pos else None)
        if f or mo:
            sibships.setdefault((f, mo), []).append(s)
    for (f, mo), kids in sibships.items():
        if f and mo:
            ox, oy = (pos[f][0] + pos[mo][0]) / 2, (pos[f][1] + pos[mo][1]) / 2
        else:
            ox, oy = pos[f or mo][0], pos[f or mo][1] + h
        bar = min(pos[k][1] for k in kids) - h - 22
        xs = [pos[k][0] for k in kids] + [ox]
        lines += [[ox, oy, ox, bar], [min(xs), bar, max(xs), bar]]
        lines += [[pos[k][0], bar, pos[k][0], pos[k][1] - h] for k in kids]
    width = 2 * pedigree_svg.MARGIN + pedigree_svg.SIZE + (widest - 1) * pedigree_svg.STEP_X
    height = pedigree_svg.TOP + pedigree_svg.SIZE + (ngen - 1) * pedigree_svg.STEP_Y + 30
    return {"width": width, "height": height, "size": pedigree_svg.SIZE, "people": people,
            "lines": [[round(v, 1) for v in ln] for ln in lines]}


def family_page(name, d):
    """family.json without the member runs: the state, the layout, the full pair rows and couple rows."""
    st = family_state(name, d)
    import render_report  # the couple screen's caveat, word for word
    out = {"state": st, "layout": None, "relationships": [], "ibd_pairs": [], "couples": {},
           "couple_note": render_report.COUPLE_NOTE, "segments": {}, "de_novo": [], "compound_het": [],
           "trio_files": {}, "side_by_side": None}
    if st.get("error"):
        return out
    order, members = members_of(d)
    out["layout"] = layout_of(order, members)
    keep = ("sample_a", "sample_b", "declared", "expected_relatedness", "relatedness", "ibs0", "ibs2", "sites",
            "verdict", "reason")
    out["relationships"] = [{k: r.get(k) for k in keep} for r in read_tsv(os.path.join(d, "relationships.tsv"))]
    keep = ("sample_a", "sample_b", "declared", "segments", "shared_cm", "ibd2_cm", "longest_cm", "relatedness",
            "degree")
    out["ibd_pairs"] = [{k: r.get(k) for k in keep} for r in read_tsv(os.path.join(d, "ibd_pairs.tsv"))]
    order_of = {"shared": 0, "x_linked": 1, "one_partner": 2}
    for a, b, base in couple_names(order, members):
        rows = read_tsv(os.path.join(d, base + ".tsv"))
        keep = ("gene", "inheritance", "condition", "result", "offspring_risk", "partner_a", "findings_a",
                "partner_b", "findings_b", "explanation")
        rows = sorted(({k: r.get(k) for k in keep} for r in rows), key=lambda r: order_of.get(r.get("result"), 3))
        out["couples"]["%s %s" % (a, b)] = rows
    out.update(child_depth(d, order, members))
    import report_charts  # GRCh38's lengths, for the painting of each child's stretches
    out["chrom_lengths"] = [[c, n] for c, n in report_charts.GRCH38 if c != "chrY"]
    # Step 38 writes the verdict "ok" for a pair that matches its declared relation. The page says that
    # in words: no word on the page reads as an all-clear (state.json keeps the file's own value).
    out["state"] = dict(st, relationships=[dict(r, verdict=verdict_words(r.get("verdict")))
                                           for r in st.get("relationships") or []])
    for r in out["relationships"]:
        r["verdict"] = verdict_words(r.get("verdict"))
    return out


SEGMENT_KEYS = ("side", "parent", "chrom", "start", "end", "length_bp", "sites", "discordant", "ancestor", "basis")
DE_NOVO_KEYS = ("child", "chrom", "pos", "ref", "alt", "child_gt", "father_gt", "mother_gt", "child_dp", "child_gq",
                "child_ab", "quality")
COMP_KEYS = ("child", "gene", "variant_a", "variant_b", "a_from", "b_from", "phase")
PHASE_ORDER = {"in trans": 0, "unknown": 1, "in cis": 2}


def child_depth(d, order, members):
    """Each child's stretches (step 38's <child>_segments.tsv) and the trio tables, as the steps wrote them:
    {segments: {child: [rows]}, de_novo: [rows], compound_het: [rows, in trans first]}."""
    segs = {}
    for s in order:
        if not (members[s]["father"] or members[s]["mother"]):
            continue
        rows = read_tsv(os.path.join(d, "%s_segments.tsv" % s))
        if rows:
            segs[s] = [{k: r.get(k) for k in SEGMENT_KEYS} for r in rows]
    files = {f: file_state(d, f) for f in TRIO_FILES}
    comp = [{k: r.get(k) for k in COMP_KEYS} for r in read_tsv(os.path.join(d, "compound_het.tsv"))]
    comp.sort(key=lambda r: PHASE_ORDER.get(r.get("phase"), 1))
    return {"segments": segs,
            "de_novo": [{k: r.get(k) for k in DE_NOVO_KEYS} for r in read_tsv(os.path.join(d, "de_novo.tsv"))],
            "compound_het": comp,
            # read_tsv's own test: a link out of the family folder counts as not there
            "trio_files": {f: files[f] == "here" for f in TRIO_FILES},
            "trio_file_state": files,
            "trio_status": trio_status(d, order, members, files)}


TRIO_FILES = ("de_novo.tsv", "compound_het.tsv")
# Step 41's own order of the blood group systems, as `og bloodtype systems`
# prints it (tests/test_live_family_depth.py checks the two agree).
STEP41_SYSTEMS = ("ABO", "RhD", "RhCE", "Kell", "Duffy", "Kidd", "MNS", "Lutheran", "Diego", "Dombrock", "Colton",
                  "Vel", "Langereis", "Junior", "H", "Secretor", *(f"HPA-{n}" for n in range(1, 36)),
                  "HNA-1", "HNA-3", "HNA-4", "HNA-5")


def file_state(d, f):
    """'here', 'rejected' (a link out of the family folder: read_tsv does not read it, and the page never
    shows where it points) or 'missing'."""
    p = os.path.join(d, f)
    if not os.path.lexists(p):
        return "missing"
    if not inside(p, os.path.abspath(d)):
        return "rejected"
    return "here" if os.path.isfile(p) else "missing"


def trio_status(d, order, members, files):
    """Per child, whether step 38 assessed it for de novo candidates and compound-het pairs, and why not:
    {child: {"de_novo": {"assessed": bool, "why": text}, "compound_het": {...}}}. The rule is step 38's
    (og family report): de novo candidates need both parents sequenced, compound-het pairs a sequenced
    parent; a file that is not here, or is a link out of the folder, means no child was assessed."""
    plan = {r.get("analysis"): r for r in read_tsv(os.path.join(d, "family_plan.tsv"))}
    trio = plan.get("trio_phasing") or {}

    def absent(f):
        if files[f] == "rejected":
            return "%s is a link out of the family folder, so this page does not read it" % f
        why = "%s is not here" % f
        if trio.get("decision"):
            why += ": trio phasing %s%s" % (trio["decision"], " (%s)" % trio["reason"] if trio.get("reason") else "")
        return why

    joint = plan.get("joint_call") or {}
    called = set((joint.get("members") or "").split(",")) if joint.get("decision") == "runs" else None

    def seq(p):
        return bool(p) and p in members and bool(members[p]["sequenced"]) and (called is None or p in called)

    need = "sequenced" if called is None else "sequenced and in the joint call"

    def not_called(p):
        if not members[p]["sequenced"]:
            return "%s is not sequenced" % p
        return "%s is not in the joint call (family_plan.tsv)" % p

    def gaps(m):
        return "; ".join(("no %s is named" % w) if not m[w] else "%s (the %s) is not" % (m[w], w)
                         for w in ("father", "mother") if not seq(m[w]))

    out = {}
    for s in order:
        m = members[s]
        if not (m["father"] or m["mother"]):
            continue
        st = {}
        if files["de_novo.tsv"] != "here":
            st["de_novo"] = {"assessed": False, "why": absent("de_novo.tsv")}
        elif not seq(s):
            st["de_novo"] = {"assessed": False, "why": not_called(s)}
        elif seq(m["father"]) and seq(m["mother"]):
            st["de_novo"] = {"assessed": True, "why": "both parents are sequenced"}
        else:
            st["de_novo"] = {"assessed": False, "why": "de novo candidates need both parents %s: %s" % (need, gaps(m))}
        if files["compound_het.tsv"] != "here":
            st["compound_het"] = {"assessed": False, "why": absent("compound_het.tsv")}
        elif not seq(s):
            st["compound_het"] = {"assessed": False, "why": not_called(s)}
        elif seq(m["father"]) or seq(m["mother"]):
            st["compound_het"] = {"assessed": True, "why": "a parent is sequenced; pairs come only from the child's "
                                  "clinical-filter candidates (step 23), so with none given no pair is listed"}
        else:
            st["compound_het"] = {"assessed": False, "why": "compound-het pairs need a parent %s" % need}
        out[s] = st
    return out


def candidates_given(live, d, sample, since=None):
    """Whether the clinical filter (step 23) wrote <sample>'s candidates, the only variants step 38 pairs:
    "here", "absent" (the run is here without them; a file linking out of the run counts as absent),
    "newer" (written after since, compound_het.tsv's mtime, so step 38 may not have had them) or None
    (the run is not on this site: the page cannot tell)."""
    if not d:
        return None
    p = os.path.join(d, "clinical", "%s_clinical.vcf.gz" % sample)
    if not (os.path.isfile(p) and live.inside(p, os.path.realpath(d)) and live.inside(p, live.root)):
        return "absent"
    return "newer" if since is not None and os.path.getmtime(p) > since else "here"


def with_candidates(status, given):
    """trio_status with the compound-het line told by the child's step 23 candidates ({child: one of
    candidates_given's answers}). assessed is None where the page cannot tell."""
    for child, st in status.items():
        c = st.get("compound_het") or {}
        if not c.get("assessed"):
            continue
        g = given.get(child)
        if g == "here":
            c["why"] = ("a parent is sequenced and the child's clinical-filter candidates (step 23) are here; "
                        "pairs come only from them")
        elif g == "absent":
            st["compound_het"] = {"assessed": False, "why": "the clinical filter (step 23) wrote no candidates for "
                                  "%s (clinical/%s_clinical.vcf.gz), so step 38 had no variant to pair" % (child, child)}
        elif g == "newer":
            st["compound_het"] = {"assessed": None, "why": "%s's clinical-filter candidates (step 23) are newer than "
                                  "compound_het.tsv, so step 38 may not have paired them" % child}
        else:
            st["compound_het"] = {"assessed": None, "why": "pairs come only from the child's clinical-filter "
                                  "candidates (step 23), and %s's run is not on this site to tell whether step 23 "
                                  "wrote them" % child}
    return status


def _sec(summary, key):
    sec = ((summary or {}).get("sections") or {}).get(key) or {}
    return sec.get("data") or {} if sec.get("state") in ("ok", "stale") else None


def _stale(summary, *keys):
    secs = (summary or {}).get("sections") or {}
    return any(((secs.get(k) or {}).get("state")) == "stale" for k in keys)


def side_values(summary):
    """[(group, label, value, stale)] of one member's summary.json, as its report prints them; nothing
    computed. stale is true when the value comes from a section older than its input (state "stale")."""
    out = []
    c = (summary or {}).get("cyp2d6") or {}
    calls = ", ".join("%s %s" % (k, v) for k, v in (c.get("calls") or {}).items() if v)
    if calls:
        agree = c.get("agree")
        # summary.json's cyp2d6 block reads pypgx and Cyrius only when complete; its PharmCAT call is cpic's
        out.append(("Medicines", "CYP2D6 (callers)", calls + (" (callers agree)" if agree else " (callers differ)"
                                                             if agree is False else ""),
                    bool((c.get("calls") or {}).get("PharmCAT")) and _stale(summary, "cpic")))
    cpic = _sec(summary, "cpic")
    for g in (cpic or {}).get("genes") or []:
        ph = g.get("phenotype")
        out.append(("Medicines", g.get("gene"), "%s%s" % (g.get("diplotype") or ".", ' "%s"' % ph if ph else ""),
                    _stale(summary, "cpic")))
    for l in (_sec(summary, "hla") or {}).get("loci") or []:
        out.append(("HLA", l.get("gene"), " ".join(l.get("alleles") or []) or "not typed", _stale(summary, "hla")))
    mt = _sec(summary, "haplogroup")
    if mt is not None:
        out.append(("Haplogroups", "Mitochondrial", mt.get("haplogroup") or ".", _stale(summary, "haplogroup")))
    ykey = "y_tree" if _sec(summary, "y_tree") else "y_haplogroup"
    y = _sec(summary, "y_tree") or _sec(summary, "y_haplogroup")
    if y is not None:
        out.append(("Haplogroups", "Y chromosome", y.get("haplogroup") or ".", _stale(summary, ykey)))
    bt = _sec(summary, "bloodtype")
    if bt is not None:
        old = _stale(summary, "bloodtype")
        systems = [x for x in bt.get("systems") or [] if x.get("system")]
        by = {x["system"]: x.get("phenotype") or "." for x in systems}
        out.append(("Blood groups", "ABO", bt.get("abo") or by.get("ABO", "."), old))
        out.append(("Blood groups", "RhD", bt.get("rhd") or by.get("RhD", "."), old))
        # every other system step 41 typed, in its order (Kell, Duffy, Kidd, MNS...)
        for x in systems:
            if x["system"] not in ("ABO", "RhD"):
                out.append(("Blood groups", x["system"], by[x["system"]], old))
    return out


def side_by_side(names, summaries):
    """The members' results in one table: {members, rows: [{group, label, values, stale}]}. A member with
    no summary.json yet says so in its column; a row a member's summary does not have is '.'. stale holds,
    per member, whether that value comes from a section older than its input. A summary whose fields are
    not of the shape the report writes costs its own column only: it says so there, and the family page
    still draws."""
    per, bad = {}, set()
    for n in names:
        try:
            got = side_values(summaries.get(n))
            for g, label, _, _ in got:
                hash((g, label))  # a label that is a list or a dict cannot key the table's rows
            per[n] = got
        except (AttributeError, KeyError, TypeError, ValueError, IndexError):
            per[n], bad = [], bad | {n}
    # Each member's rows keep their own order (step 41's for the blood groups): a row first seen in a later
    # member goes right after the row before it in that member's list, not at the end.
    keys = []
    for n in names:
        prev = None
        for g, label, _, _ in per[n]:
            if (g, label) not in keys:
                keys.insert(keys.index(prev) + 1 if prev else len(keys), (g, label))
            prev = (g, label)
    groups = ["Medicines", "HLA", "Haplogroups", "Blood groups"]
    step41 = {s: i for i, s in enumerate(STEP41_SYSTEMS)}   # for systems two members each lack one of
    pos = {k: i for i, k in enumerate(keys)}
    keys.sort(key=lambda k: (groups.index(k[0]) if k[0] in groups else 9,
                             step41.get(k[1], len(step41)) if k[0] == "Blood groups" else 0, pos[k]))
    rows = [{"group": g, "label": label, "values": [], "stale": []} for g, label in keys]
    for r in rows:
        for n in names:
            if n in bad:
                r["values"].append("summary not readable")
                r["stale"].append(False)
            elif summaries.get(n) is None:
                r["values"].append("no summary yet")
                r["stale"].append(False)
            else:
                hit = next((x for x in per[n] if (x[0], x[1]) == (r["group"], r["label"])), None)
                r["values"].append(hit[2] if hit else ".")
                r["stale"].append(bool(hit and hit[3]))
    return {"members": list(names), "rows": rows}


def verdict_words(v):
    return "matches the declared relation" if (v or "").strip().lower() == "ok" else v


def member_summary(live, d):
    """A member run's summary.json as a dict, or None (none yet, unreadable, or a link out of the run)."""
    import json
    p = os.path.join(d, "summary.json") if d else None
    if not p or not os.path.isfile(p) or not live.inside(p, os.path.realpath(d)) or not live.inside(p, live.root):
        return None
    try:
        with open(p) as f:
            got = json.load(f)
    except (OSError, ValueError):
        return None
    return got if isinstance(got, dict) else None


def member_state(root, sample):
    """The families that list the sample. A family whose pedigree cannot be read may be the one that
    lists it as a child: with no readable family listing it, the answer is unknown, never "in no
    family" ({"kind": "member", "families": [], "unknown": true}), and its pages keep results to act
    on closed, as the home page does for a run no readable family lists."""
    fams, unread = [], False
    for name, d in family_dirs(root):
        try:
            order, _ = members_of(d)
        except (OSError, ValueError):
            unread = True
            continue
        if sample in order:
            fams.append({"name": name, "url": "/runs/%s/" % name})
    if fams:
        return {"kind": "member", "families": fams}
    return {"kind": "member", "families": [], "unknown": True} if unread else None


def register(app):
    def root():
        return os.path.realpath(app.config.root)

    def runs(app_):
        out = []
        for n, d in family_dirs(root()):
            try:
                order, members = members_of(d)
                people = [{"sample": s, "sex": members[s]["sex"], "father": members[s]["father"],
                           "mother": members[s]["mother"], "partner": members[s]["partner"],
                           "sequenced": members[s]["sequenced"]} for s in order]
            except (OSError, ValueError):
                people = []
            out.append({"name": n, "dir": d, "kind": "family", "members": people})
        return out

    def page_data(req, m):
        name = m.group(1)
        d = req.run_dir(name)
        if not d or not is_family(d) or not inside(d, root()):
            app.log("family", name, "family.json", 404)
            return req.send_body(404, "text/plain; charset=utf-8", "not found\n")
        out = family_page(name, d)
        live = getattr(req, "live", None)
        runs = {}
        for p in out["state"].get("members") or []:
            r = live.run(p["sample"]) if live else None
            if not r or r.get("kind") == "family":
                continue
            st = live.base_state(r)
            planned = [s for s in st["steps"] if s["state"] != "not planned"]
            n, mtime = live.attention(r["dir"] or os.path.join(live.root, r["name"]))
            runs[p["sample"]] = {"url": "/runs/%s/" % p["sample"], "state": st["run"]["state"],
                                 "planned": len(planned), "done": sum(1 for s in planned if s["state"] == "done"),
                                 "attention": n, "results": mtime is not None}
        out["runs"] = runs
        given = {}
        comp_p = os.path.join(d, "compound_het.tsv")
        since = os.path.getmtime(comp_p) if os.path.isfile(comp_p) else None
        for c in out.get("trio_status") or {}:
            r = live.run(c) if live else None
            given[c] = (candidates_given(live, r["dir"] or os.path.join(live.root, r["name"]), c, since)
                        if r and r.get("kind") != "family" else None)
        with_candidates(out.get("trio_status") or {}, given)
        # The members side by side, from each sequenced member's summary.json alone (read as the home
        # page's attention count reads it: never a file linking out of its run).
        seq = [p["sample"] for p in out["state"].get("members") or [] if p.get("sequenced")]
        sums = {}
        for s_ in seq:
            r = live.run(s_) if live else None
            sums[s_] = member_summary(live, r["dir"] or os.path.join(live.root, r["name"])) if r else None
        out["side_by_side"] = side_by_side(seq, sums) if seq else None
        app.log("family", name, "family.json", 200)
        return req.send_json(out)

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
        try:
            order, members = members_of(d)
        except (OSError, ValueError):
            app.log("family", name, f, 404)
            return req.send_body(404, "text/plain; charset=utf-8", "not found\n")
        if f == "pedigree.svg":
            app.log("family", name, f, 200)
            return req.send_body(200, TYPES[".svg"], pedigree_svg.draw(order, members, name))
        p = os.path.join(d, f)
        if f not in served_files(d, order, members) or not inside(p, d):
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
