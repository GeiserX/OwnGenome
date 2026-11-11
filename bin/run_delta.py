#!/usr/bin/env python3
"""run_delta.py: what changed since the previous run of the same sample.

The command is `og run delta` (rust/crates/og-runlog, docs/delta-report.md):
scripts/generate-report.sh runs it after the text report, before it keeps a
copy of this run's summary.json and run_manifest.tsv under
<sample>/history/<UTC time>/, and it writes <out>.txt, <out>.html and
<out>.json. This module is the same comparison for bin/collect_summary.py,
whose "New since" badges need it. Standard library only.

What it compares:
  inputs   the run manifest (bin/write_manifest.sh): versions.env values, the
           image each *_IMAGE names and the digest it resolved to, the data
           rows (ClinVar file date, VEP cache, CPSR bundle, pypgx bundle, HLA
           database, ancestry panel, PGS headers), the pipeline commit and the
           declared sex. When and by which script it was written are left out.
  sections a report section that ran in one run and not in the other.
  ClinVar  hits added or dropped (by chrom:pos ref>alt), a hit whose
           significance changed (reclassified), whose review status or stars
           changed (re-starred), or whose genotype changed.
  CPSR     the count of each classification.
  PGx      each CPIC gene's diplotype and phenotype, pypgx's CYP2D6 call and
           the CYP2D6 consensus.
  PRS      scores added or dropped, and a percentile that moved by
           PRS_BAND points or more, or appeared or went away.
For each result change it lists the changed inputs that feed that section, or
says that the inputs it tracks are the same and the call itself differs.
changed_sections() names the summary sections the changes touch, and
read_recheck() the last database re-check (recheck/state.json beside the
summary, bin/owngenome_recheck.py).
"""
import json
import os

PRS_BAND = 5.0  # percentile points

# Changed inputs that can explain a change in a section: an input whose name
# contains one of these (case-insensitive). The pipeline commit explains any.
WHY_TOKENS = {
    "clinvar": ["clinvar", "DEEPVARIANT", "BCFTOOLS", "MINIMAP2", "BWA", "GATK", "FASTP"],
    "cpsr": ["PCGR", "VEP", "DEEPVARIANT"],
    "cpic": ["PHARMCAT", "PYPGX", "CYRIUS", "DEEPVARIANT", "BCFTOOLS"],
    "pypgx": ["PYPGX", "pypgx_bundle"],
    "cyp2d6": ["PYPGX", "pypgx_bundle", "CYRIUS", "PHARMCAT"],
    "prs": ["pgs:", "PGSC", "PLINK2", "ancestry_panel", "DEEPVARIANT"],
}
IGNORED_RUN_KEYS = ("written_utc", "written_by")
# The summary section of each result area (diff_states names its own).
SECTION_OF_AREA = {"ClinVar": "clinvar", "CPSR": "cpsr", "PGx (CPIC)": "cpic", "PGx (pypgx)": "pypgx",
                   "PGx (CYP2D6 consensus)": "pypgx", "PRS": "prs"}


def read_manifest(path):
    """{(section, key): value} of a run_manifest.tsv; {} when there is none."""
    rows = {}
    if not path or not os.path.isfile(path):
        return rows
    with open(path, errors="replace") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 3 and not c[0].startswith("#"):
                rows[(c[0], c[1])] = "\t".join(c[2:])
    return rows


def manifest_items(rows):
    """The manifest as {item: value}, one entry per thing that can change.
    An image is keyed by its *_IMAGE name, its value is 'ref (digest)', so a
    new tag and a new digest of the same tag are each one change."""
    out = {}
    for (sec, key), val in rows.items():
        if sec == "run":
            if key not in IGNORED_RUN_KEYS:
                out[key] = val
        elif sec == "version":
            if key.endswith("_IMAGE"):
                out[key] = f"{val} ({rows.get(('image', val), 'not in the manifest')})"
            else:
                out[key] = val
        elif sec == "data":
            out[key] = val
    return out


def change(area, kind, item, before, after, why=None, section=None):
    c = {"area": area, "kind": kind, "item": item, "before": before, "after": after}
    if why is not None:
        c["why"] = why
    if section is not None:
        c["section"] = section  # the summary.json section it changes: the report marks it "new since"
    return c


def diff_inputs(prev, cur):
    a, b = manifest_items(prev), manifest_items(cur)
    out = []
    for k in sorted(set(a) | set(b)):
        if a.get(k) != b.get(k):
            out.append(change("inputs", "input", k, a.get(k, "absent"), b.get(k, "absent")))
    return out


def why_for(area, inputs):
    tokens = [t.lower() for t in WHY_TOKENS.get(area, [])] + ["git_commit"]
    hit = [c for c in inputs if any(t in c["item"].lower() for t in tokens)]
    if not hit:
        return "the inputs this delta tracks for it are the same; the call itself differs"
    return "changed inputs: " + "; ".join(f"{c['item']} {c['before']} -> {c['after']}" for c in hit)


def section(s, key):
    return ((s.get("sections") or {}).get(key)) or {}


def has_data(sec):
    return sec.get("state") in ("ok", "stale")


def diff_states(p, c):
    out = []
    keys = sorted(set(p.get("sections") or {}) | set(c.get("sections") or {}))
    for k in keys:
        a, b = section(p, k), section(c, k)
        # A section only one summary knows (a newer version added it) and that did not run reads as
        # missing on both sides, so a pipeline upgrade lists no noise rows.
        sa, sb = a.get("state", "missing"), b.get("state", "missing")
        if sa != sb:
            out.append(change(b.get("title") or a.get("title") or k, "section", k, sa, sb, section=k))
    return out


def hit_key(h):
    return f"{h.get('chrom')}:{h.get('pos')} {h.get('ref')}>{h.get('alt')}"


def diff_clinvar(p, c):
    a, b = section(p, "clinvar"), section(c, "clinvar")
    if not (has_data(a) and has_data(b)):
        return []
    ha = {hit_key(h): h for h in a.get("data", {}).get("hits", [])}
    hb = {hit_key(h): h for h in b.get("data", {}).get("hits", [])}
    out = []
    for k in sorted(set(ha) | set(hb)):
        x, y = ha.get(k), hb.get(k)
        gene = (y or x).get("gene", ".")
        item = f"{k} ({gene})"
        if x is None:
            out.append(change("ClinVar", "added", item, "absent", f"{y.get('significance')}, {y.get('stars')} stars"))
        elif y is None:
            out.append(change("ClinVar", "dropped", item, f"{x.get('significance')}, {x.get('stars')} stars", "absent"))
        else:
            if x.get("significance") != y.get("significance"):
                out.append(change("ClinVar", "reclassified", item, x.get("significance"), y.get("significance")))
            if x.get("stars") != y.get("stars") or x.get("review_status") != y.get("review_status"):
                out.append(change("ClinVar", "re-starred", item,
                                  f"{x.get('stars')} stars ({x.get('review_status')})",
                                  f"{y.get('stars')} stars ({y.get('review_status')})"))
            if x.get("genotype") != y.get("genotype"):
                out.append(change("ClinVar", "genotype", item, x.get("genotype"), y.get("genotype")))
    return out


def diff_cpsr(p, c):
    a, b = section(p, "cpsr"), section(c, "cpsr")
    if not (has_data(a) and has_data(b)):
        return []
    ca, cb = a.get("data", {}).get("classification") or {}, b.get("data", {}).get("classification") or {}
    return [change("CPSR", "classification count", k, ca.get(k, 0), cb.get(k, 0))
            for k in sorted(set(ca) | set(cb)) if ca.get(k, 0) != cb.get(k, 0)]


def diff_cpic(p, c):
    a, b = section(p, "cpic"), section(c, "cpic")
    if not (has_data(a) and has_data(b)):
        return []
    ga = {g.get("gene"): g for g in a.get("data", {}).get("genes", [])}
    gb = {g.get("gene"): g for g in b.get("data", {}).get("genes", [])}
    out = []
    for k in sorted(set(ga) | set(gb), key=str):
        x, y = ga.get(k), gb.get(k)
        fx = f"{x.get('diplotype')} ({x.get('phenotype')})" if x else "absent"
        fy = f"{y.get('diplotype')} ({y.get('phenotype')})" if y else "absent"
        if fx != fy:
            out.append(change("PGx (CPIC)", "call", k, fx, fy))
    return out


def diff_cyp2d6(p, c):
    out = []
    a, b = section(p, "pypgx"), section(c, "pypgx")
    if has_data(a) and has_data(b):
        x, y = a.get("data", {}).get("cyp2d6"), b.get("data", {}).get("cyp2d6")
        if x != y:
            out.append(change("PGx (pypgx)", "call", "CYP2D6", x, y))
    x = ((p.get("cyp2d6") or {}).get("consensus") or {}).get("result")
    y = ((c.get("cyp2d6") or {}).get("consensus") or {}).get("result")
    if x != y:
        out.append(change("PGx (CYP2D6 consensus)", "call", "CYP2D6", x or "none", y or "none"))
    return out


def pct(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def shown(v):
    return "none" if v is None else v


def diff_prs(p, c):
    a, b = section(p, "prs"), section(c, "prs")
    if not (has_data(a) and has_data(b)):
        return []
    sa = {r.get("pgs_id"): r for r in a.get("data", {}).get("scores", [])}
    sb = {r.get("pgs_id"): r for r in b.get("data", {}).get("scores", [])}
    out = []
    for k in sorted(set(sa) | set(sb), key=str):
        x, y = sa.get(k), sb.get(k)
        item = f"{k} ({(y or x).get('condition', '')})"
        if x is None or y is None:
            out.append(change("PRS", "added" if x is None else "dropped", item,
                              "absent" if x is None else "scored", "scored" if x is None else "absent"))
            continue
        px, py = pct(x.get("percentile")), pct(y.get("percentile"))
        if (px is None) != (py is None) or (px is not None and abs(px - py) >= PRS_BAND):
            # A percentile of 0 is a value: only a missing one reads "none".
            out.append(change("PRS", "percentile", item, shown(x.get("percentile")), shown(y.get("percentile"))))
    return out


def compare(prev_summary, prev_manifest, summary, manifest):
    inputs = diff_inputs(prev_manifest, manifest)
    results = []
    for area, fn in (("clinvar", diff_clinvar), ("cpsr", diff_cpsr), ("cpic", diff_cpic),
                     ("cyp2d6", diff_cyp2d6), ("prs", diff_prs)):
        for ch in fn(prev_summary, summary):
            ch["why"] = why_for(area, inputs)
            ch["section"] = SECTION_OF_AREA.get(ch["area"], area)
            results.append(ch)
    return inputs + diff_states(prev_summary, summary) + results


def changed_sections(changes):
    """The summary sections a delta changes, in order: what the report marks 'New since'."""
    return list(dict.fromkeys(c["section"] for c in changes if c.get("section")))


def newest_history(sample_dir):
    """The newest history/<stamp>/ folder with a summary.json, or None."""
    h = os.path.join(sample_dir, "history")
    try:
        names = sorted(n for n in os.listdir(h) if os.path.isfile(os.path.join(h, n, "summary.json")))
    except OSError:
        return None
    return os.path.join(h, names[-1]) if names else None


def read_recheck(path):
    """<sample>/recheck/state.json (bin/owngenome_recheck.py) as a short block, or None."""
    try:
        with open(path) as f:
            st = json.load(f)
    except (OSError, ValueError):
        return None
    if not isinstance(st, dict) or not isinstance(st.get("items"), list):
        return None
    # A check of an earlier run: <sample>/run_manifest.tsv was written after it, so its "newer than
    # the run used" may already be what this run used. Not shown; the next check replaces it.
    man = read_manifest(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(path))), "run_manifest.tsv"))
    if utc14(man.get(("run", "written_utc"))) > utc14(st.get("checked_utc")):
        return None
    return {"checked_utc": st.get("checked_utc"),
            "items": [{k: i.get(k) for k in ("item", "label", "class", "ran_with", "upstream", "action")}
                      for i in st["items"] if isinstance(i, dict)],
            "affected_steps": st.get("affected_steps") or [],
            "rerun": st.get("rerun") if isinstance(st.get("rerun"), dict) else None}


def utc14(text):
    """'2026-09-02 12:00:00 UTC' or '2026-09-02T12:00:00Z' as '20260902120000'; '' without one."""
    return "".join(c for c in str(text or "") if c.isdigit())[:14]


if __name__ == "__main__":
    raise SystemExit("bin/run_delta.py has no command line any more: run bin/og run delta (docs/delta-report.md)")
