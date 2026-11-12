#!/usr/bin/env python3
"""live_results.py: one person's results as the live site draws them (results.json).

  build(sample, summary, files=()) -> dict

Pure functions, standard library only. The live site's results plugin
(bin/owngenome_live_results.py) and the public demo (tests/demo/build_live_demo.py)
call build() on a summary.json (bin/collect_summary.py). The answer follows the
report's own chapters (bin/render_report.py CHAPTERS): seven questions, one card per
section id, the same anchors as report.html. Each card carries:

  facts       [[label, value]] the numbers the report prints
  table       {"columns": [...], "rows": [[...]]} the rows the report lists, or null
  lines       the report's own sentences (render_report's *_lines helpers), plain text
  caveat      the section's caveats word for word (render_report.SECTION_NOTES)
  attention   true only for a flag the report or the live page already raises: a
              headline tone (bin/owngenome_live.py HEADLINES), a failed, unreadable or
              stale result, a repeat outside Stranger's "normal" or in an intermediate
              zone, an AnnotSV class 4 or 5 row, a CHIP or MoChA call, a two-hit gene
              or a couple's shared or X-linked gene. Nothing says a section is fine:
              a card without attention is "for information".

Beside the cards, "veil" tags the rows that name a result to act on outside its own card (every card's
table rows and lines, the viewer links, the charts' marks, the other findings and the lead's names) with
that result's id, so a child's page closes them with it (veil() below).

Nothing here reads a file or computes anything about genetics; the adapters pick
fields the summary already holds. A section id with no adapter falls back to its
report headline and a link to the report.
"""
import hashlib
import re

import render_report as rr

SCHEMA = 1
# Words that would tell a reader a result is fine. tests/test_live_results.py scans every string of
# results.json and family.json, and the page assets, for them. "Fine-grained" is a method's name.
FORBIDDEN = re.compile(r"\b(ok|(?<![-/])fine(?!-grained)|healthy|all.clear|nothing found|normal range)\b", re.I)
EMPTY = {
    # The live page's no-hit sentence (render_report's ClinVar card shows a gray 0, never green).
    "clinvar": "No P/LP match in this ClinVar file at the call floor. ClinVar lists only part of what can cause "
               "disease, so this is information, and the caveats below still apply.",
}
STATE_WORDS = {"failed": "Failed", "unreadable": "Unreadable", "stale": "Older than its input"}


def plain(text):
    """A report line as plain text: markdown bold dropped, white space folded, and a check's bare
    "ok" (pypgx's depth check) said as "passed", so no line reads as an all-clear."""
    text = re.sub(r"\s+", " ", str(text).replace("**", "")).strip()
    return re.sub(r"(?<=: )ok\b", "passed", text)


def cell(v):
    if v is None or v == "":
        return "."
    if isinstance(v, float):
        return f"{v:g}"
    if isinstance(v, (list, tuple)):
        return ", ".join(cell(x) for x in v) or "."
    return str(v)


def num(v):
    """A summary value as a float, or None ('.', '', 'NA', a range)."""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if x == x and x not in (float("inf"), float("-inf")) else None


def table(columns, rows, limit=200):
    rows = [[cell(c) for c in r] for r in rows]
    if not rows:
        return None
    out = {"columns": list(columns), "rows": rows[:limit]}
    if len(rows) > limit:
        out["more"] = len(rows) - limit
    return out


def _headline_tone(key, sec, summary):
    """(text, tone) from the live site's headline table, or (None, 'info')."""
    try:
        import owngenome_live
        h = owngenome_live.headline(key, sec, summary)
    except Exception:  # the page must not go down on one reader
        h = None
    return (plain(h[0]), h[1]) if h else (None, "info")


# --- the adapters: one per section id --------------------------------------------------------
# Each takes (S, summary) where S is summary["sections"], and returns a dict with any of
# facts, table, lines, attention, links (relative file patterns), extra.

def a_qc(S, s):
    cov, sexc, q, rdc = (S.get(k) or {} for k in ("coverage", "sex_check", "sample_qc", "read_count"))
    facts = []
    md = (cov.get("data") or {}).get("mean_depth")
    if isinstance(md, (int, float)):
        facts.append(["Mean depth", f"{md:.1f}x"])
    if rdc.get("state") in ("ok", "stale"):
        d = rdc["data"]
        facts.append(["Reads kept by alignment", rr.read_count_text(d)
                      + ("" if d.get("read_count_check") != "mismatch" else " (reads lost: see the report)")])
    inferred = (sexc.get("data") or {}).get("inferred_sex")
    declared = (s.get("run") or {}).get("declared_sex")
    if inferred:
        facts.append(["Inferred sex (indexcov)", inferred])
    facts.append(["Declared sex", declared or "not declared"])
    attention = False
    if q.get("state") in ("ok", "stale"):
        d = q["data"]
        check = d.get("sex_check")
        word = {"ok": "Declared and measured sex agree", "mismatch": "Declared and measured sex differ",
                "not_checked": f"Not checked: {d.get('sex_check_reason') or 'no reason given'}"}.get(check, check or ".")
        facts.append(["Sex from the reads (somalier)", f"{d.get('inferred_sex') or '.'}. {word}"])
        attention = check == "mismatch"
        cont = d.get("contamination")
        if cont == "ok":
            facts.append(["Contamination (VerifyBamID2)", f"Below the contamination flag (FREEMIX {d.get('freemix')}, "
                                                          f"flag above {d.get('freemix_warn_above')})"])
        else:
            facts.append(["Contamination (VerifyBamID2)", rr.freemix_text(d)])
            attention = attention or cont == "warn"
        if d.get("same_person_as"):
            facts.append(["The same person as", d["same_person_as"]])
            attention = True
    lines = [x for x in rr.karyotype_lines(S)]
    if lines[:-1]:
        attention = True
    return {"facts": facts, "lines": lines, "attention": attention,
            "links": [r"(^|/)multiqc[^/]*\.html$", r"(^|/)fastp[^/]*\.html$", r"_fastp\.html$"]}


def a_callability(S, s):
    d = S["callability"]["data"]
    rows = [[g["gene"], g.get("status"), "." if g.get("share") is None else f"{100 * g['share']:.0f}%",
             ", ".join(g.get("reasons") or [])] for g in d.get("not_read") or []]
    return {"facts": [["Hard-to-read genes read", f"{d.get('read')} of {d.get('catalog_genes')}"],
                      ["Not read, no other step calls them", len(d.get("not_read") or [])],
                      ["Not read here, called by a step of their own", len(d.get("not_read_owned") or [])]],
            "table": table(["Gene", "Status", "Bases read", "Why"], rows)}


def a_variants(S, s):
    d = S["variants"]["data"]
    return {"facts": [["Total variants", d.get("total")], ["PASS variants", d.get("pass")], ["SNPs", d.get("snps")],
                      ["Indels", d.get("indels")]], "links": [r"(^|/)vep/.*summary\.html$"]}


def clinvar_call(h):
    """A hit's Call cell: the low-confidence mark, the cDNA read check and the combined-allele note,
    each as the report prints it."""
    parts = [rr.call_text(h) or ("low confidence" if h.get("low_confidence") else "")]
    if "cdna" in h:
        parts.append("reads: " + rr.cdna_text(h))
    parts.append(rr.combined_cell(h))
    return "; ".join(x for x in parts if x) or "."


def a_clinvar(S, s):
    d = S["clinvar"]["data"]
    by = d.get("by_stars") or {}
    facts = [["P/LP matches", d.get("count", 0)],
             ["By review stars (4/3/2/1/0)", " / ".join(str(by.get(k, 0)) for k in ("4", "3", "2", "1", "0"))],
             ["ClinVar file date", (s.get("databases") or {}).get("clinvar_release") or "unknown"]]
    if "low_confidence" in d:
        f = d.get("call_floor") or {}
        facts.append([f"Low-confidence calls (GQ below {f.get('GQ')} or DP below {f.get('DP')}), not counted",
                      d["low_confidence"]])
    rows = []
    for h in (d.get("hits") or [])[:50]:
        rows.append([h.get("gene"), f"{h.get('chrom')}:{h.get('pos')} {h.get('ref')}>{h.get('alt')}", h.get("genotype"),
                     h.get("significance"), "★" * int(h.get("stars") or 0) or "0",
                     rr.af_text(h) or "frequency unknown", clinvar_call(h)])
    lines = [x for x in (rr.cdna_scan_lines(d) or [])]
    for x in d.get("deep_intronic") or []:
        # Always reported by the report, hit or not: each site's call stays visible on the card.
        lines.append(f"Known deep intronic site (always reported): {x.get('gene')} {x.get('hgvs')} "
                     f"{x.get('chrom')}:{x.get('pos')} {x.get('ref')}>{x.get('alt')}: {rr.deep_intronic_call(x)}")
    freq = rr.clinvar_freq_line(d)
    if freq:
        lines.insert(0, freq)
    viewer = [f"{h.get('chrom')}:{h.get('pos')}" for h in (d.get("hits") or [])[:50]]
    return {"facts": facts, "lines": lines,
            "table": table(["Gene", "Variant", "Genotype", "Significance", "Stars", "Frequency", "Call"], rows),
            "viewer": viewer}


def a_clinical_filter(S, s):
    d = S["clinical"]["data"]
    sf = d.get("acmg_sf") or {}
    bi = d.get("by_impact") or {}
    facts = [["Variants after the clinical filter", d.get("variants")],
             ["By impact (high / moderate / low)", f"{bi.get('HIGH', 0)} / {bi.get('MODERATE', 0)} / {bi.get('LOW', 0)}"],
             ["Genes", d.get("genes")]]
    if sf.get("unavailable"):
        # Not assessed is never zero findings: the report's sentence instead of the counts.
        return {"facts": facts, "lines": rr.sf_apart_lines(sf), "attention": True}
    rows = [[r[0], r[1], r[2], r[3]] for r in rr.sf_rows(sf)] if sf else []
    high = f"{rr.high_rare(sf)}{rr.high_apart(sf)}" if sf.get("high_impact_checked") else "not checked"
    facts += [[f"{sf.get('version') or 'ACMG SF'}: ClinVar P/LP in its genes", len(sf.get("clinvar_hits") or [])],
              ["High-impact rare variants in its genes", high]]
    return {"facts": facts,
            "table": table(["Gene", "Variant", "Kind", "Rule"], rows),
            "attention": bool(sf.get("clinvar_hits") or sf.get("high_impact") or sf.get("structural"))}


def a_cpsr(S, s):
    d = S["cpsr"]["data"]
    facts = [["CPSR report", "Ran" if d.get("html_report") else "No HTML"]]
    facts += [[k.replace("_", " "), v] for k, v in (d.get("classification") or {}).items()]
    rows = [[x.get("gene"), x.get("variant"), x.get("final"), x.get("clinvar") or ".",
             "low confidence" if x.get("low_confidence") else "."] for x in d.get("plp") or []]
    return {"facts": facts, "lines": rr.cpsr_lines(d, s.get("databases") or {}),
            "table": table(["Gene", "Variant", "Final class", "ClinVar", "Call"], rows), "links": [r"^cpsr/"]}


def a_slivar(S, s):
    d = S["slivar"]["data"]
    return {"facts": [["Prioritized variants", d.get("prioritized")],
                      ["Compound-het candidate variants", d.get("compound_het_variants")],
                      ["Compound-het candidate genes", d.get("compound_het_genes")]],
            "lines": rr.comphet_phase_text(d) if "compound_het_phase" in d else []}


def a_exomiser(S, s):
    d = S["exomiser"]["data"]
    rows = []
    for g in d.get("genes") or []:
        for v in g.get("variants") or [{}]:
            rows.append([g.get("rank"), g.get("gene"), g.get("moi"), g.get("exomiser_score"),
                         v.get("hgvs") or v.get("variant"), v.get("genotype"), v.get("acmg"), v.get("disease")])
    return {"facts": [["Symptoms searched", ", ".join(f"{t.get('label')} ({t.get('id')})" for t in d.get("terms") or [])
                       or "none"]],
            "table": table(["Rank", "Gene", "Inheritance", "Score", "Variant", "Genotype", "ACMG", "Disease"], rows)}


def a_structural(S, s):
    m, dl, cn, sv = ((S.get(k) or {}).get("data") or {} for k in ("manta", "delly", "cnvpytor", "sv_consensus"))
    av = (S.get("annotsv") or {}).get("data") or {}
    facts = [["Manta SVs (total / PASS)", f"{m.get('total', '.')} / {m.get('pass', '.')}"],
             ["Delly SVs (total / PASS)", f"{dl.get('total', '.')} / {dl.get('pass', '.')}"],
             ["CNVpytor CNVs (significant)", f"{cn.get('total', '.')} ({cn.get('significant', '.')})"],
             ["Consensus SVs (2+ callers)", sv.get("consensus", ".")],
             ["AnnotSV class 4 and 5", av.get("count", ".")]]
    rows = [[x.get("type"), x.get("sv"), x.get("gt"), x.get("length"), x.get("acmg_class"),
             rr.annotsv_genes_text(x) if x.get("genes") is not None else ".", rr.annotsv_support_text(x),
             "; ".join(rr.annotsv_flags(x)) if "exon_edge" in x else "."]
            for x in av.get("rows") or []]
    return {"facts": facts, "table": table(["Type", "Where", "Genotype", "Length", "Class", "Genes", "Support",
                                            "Checks"], rows),
            "attention": bool(av.get("count"))}


def a_mei(S, s):
    d = S["mei"]["data"]
    rows = [[x.get("gene"), x.get("element"), x.get("condition"), x.get("inheritance"), x.get("status"),
             x.get("mean_depth")] for x in d.get("sites") or []]
    return {"facts": [["Insertions called", d.get("calls")],
                      ["By element", ", ".join(f"{k} {v}" for k, v in (d.get("by_element") or {}).items()) or "none"]],
            "table": table(["Gene", "Element", "Condition", "Inheritance", "Known site", "Depth"], rows),
            "attention": any(x.get("status") not in (None, "", "not seen") for x in d.get("sites") or [])}


def a_expansions(S, s):
    d = S["expansions"]["data"]
    rows, flagged, outside = [], False, 0
    for l in d.get("loci") or []:
        status = l.get("status") or ""
        out = any(x and x != "normal" for x in status.split(","))
        outside += out
        flagged = flagged or out or bool(l.get("zone"))
        shown = rr.str_status(l)
        count = rr.str_count(l) if isinstance(l.get("repeat_count"), str) else l.get("repeat_count")
        rows.append([l.get("locus"), count,
                     f'"{shown}"' if shown not in ("", "not graded", "not sized") else (shown or "not graded"),
                     l.get("normal_max"), l.get("pathologic_min"), l.get("ci"), l.get("sizing"),
                     f'{l.get("zone_label")} ({l.get("zone")})' if l.get("zone") else "."])
    facts = [["Loci tested", d.get("records")],
             ["Outside Stranger's \"normal\"", outside if d.get("stranger") else "Stranger not run"],
             ["In an intermediate zone", ", ".join(d.get("intermediate") or []) or "none"],
             ["Not sized from short reads", ", ".join(d.get("unsized") or []) or "none"]]
    return {"facts": facts, "attention": flagged,
            "table": table(["Locus", "Repeats", "Stranger status", "Normal up to", "Pathologic from", "95% interval",
                            "Sizing", "Intermediate zone"], rows),
            # The gauge: each allele's count against Stranger's thresholds, its 95% interval as a whisker.
            # status is render_report's str_status: never Stranger's "normal" for an allele not sized.
            "widget": {"kind": "str", "loci": [dict({k: l.get(k) for k in ("locus", "repeat_count", "normal_max",
                                                                               "pathologic_min", "ci", "allele_sizing")},
                                                    status=rr.str_status(l))
                                               for l in d.get("loci") or []]}}


def a_vntyper(S, s):
    d = S["vntyper"]["data"]
    return {"facts": [["Result", d.get("status")], ["Confidence", d.get("confidence")],
                      ["Mean depth over the repeat", d.get("vntr_mean_depth")]],
            "lines": rr.vntyper_lines(d), "attention": d.get("status") in ("positive", "possible")}


def a_mito(S, s):
    d = S["mito"]["data"]
    a = d.get("annotation") or {}
    facts = [["PASS variants on chrM", d.get("pass")], ["Heteroplasmic sites", d.get("heteroplasmic")]]
    if "clinvar_plp" in a:
        facts.append(["ClinVar pathogenic or likely pathogenic", a.get("clinvar_plp")])
    rows = [[f"m.{x.get('pos')}{x.get('ref')}>{x.get('alt')}", f"{100 * x['af']:.1f}%" if isinstance(x.get("af"), (int, float))
             else "."] for x in d.get("heteroplasmic_sites") or []]
    return {"facts": facts, "lines": rr.mito_lines(d), "table": table(["Site", "Allele fraction"], rows),
            "attention": (a.get("clinvar_plp") or 0) > 0 or bool(d.get("confirmed"))}


def a_lpa(S, s):
    d = S["lpa"]["data"]
    facts = [["KIV-2 repeats", f"{d.get('kiv2_copies')} (95% interval {d.get('ci95_low')} to {d.get('ci95_high')})"
              if d.get("status") == "estimated" else d.get("status")], ["Class", d.get("class")]]
    if d.get("null_allele"):
        facts.append(["KIV-2 4925G>A null allele", rr.marker_state(d["null_allele"])])
    out = {"facts": facts, "lines": rr.lpa_lines(d), "attention": d.get("flag") == "high"}
    v, lo, hi = num(d.get("kiv2_copies")), num(d.get("ci95_low")), num(d.get("ci95_high"))
    if d.get("status") == "estimated" and v is not None and lo is not None and hi is not None:
        # The count and its 95% interval on a scale from 0: the page draws the bar, the facts say it in words.
        out["widget"] = {"kind": "ci", "label": "KIV-2 repeats, both chromosomes", "value": v, "low": lo, "high": hi,
                         "unit": "repeats", "max": max(60, 10 * int(hi * 1.25 / 10 + 1)), "class": d.get("class")}
    return out


def a_chip(S, s):
    d = S["chip"]["data"]
    rows = [[c.get("gene"), c.get("protein_change"), c.get("vaf"), c.get("alt_reads"), c.get("depth"), c.get("status")]
            for c in d.get("calls") or []]
    genes = [{"gene": g.get("gene"), "limit": num(g.get("detection_limit_vaf")), "depth": num(g.get("mean_depth")),
              "status": g.get("status")} for g in d.get("genes") or []]
    # Each call carries its own gene's limit: the genes' limits differ, and a call is read against its own.
    limit = {g["gene"]: g["limit"] for g in genes}
    calls = [{"gene": c.get("gene"), "change": c.get("protein_change"), "vaf": num(c.get("vaf")),
              "limit": limit.get(c.get("gene"))} for c in d.get("calls") or []]
    out = {"facts": [["CHIP calls", len(d.get("calls") or [])],
                     ["Typical detection limit at this depth (VAF)", d.get("median_limit")]],
           "table": table(["Gene", "Change", "VAF", "Alt reads", "Depth", "Status"], rows),
           "lines": rr.chip_lines(d, s.get("sample")), "attention": bool(d.get("calls"))}
    if genes or calls:
        # Each call's VAF against the smallest VAF this depth can see in each gene (chip_report's own limit):
        # under the limit nothing can be seen, so an empty track is not an all-clear.
        out["widget"] = {"kind": "vaf", "max": 0.5, "genes": genes, "calls": [c for c in calls if c["vaf"] is not None]}
    return out


def a_mosaic(S, s):
    q = S["mosaic_sex_loss"]["data"]
    call = rr.MOSAIC_CALL.get(q.get("call"), (q.get("call"), ""))[0]
    return {"facts": [[q.get("chrom") or "chrY or chrX", call]], "lines": [q.get("reason") or ""] + rr.mosaic_lines(q),
            "attention": q.get("call") in ("detected", "gain")}


def a_mosaic_mca(S, s):
    d = S["mosaic"]["data"]
    rows = [[c.get("chrom"), c.get("start"), c.get("end"), c.get("length_mb"), c.get("type"), c.get("cell_fraction")]
            for c in d.get("calls") or []]
    return {"facts": [["Status", d.get("status")], ["Calls reported", d.get("calls_reported", len(rows))],
                      ["Smallest cell fraction it can see", d.get("floor")]],
            "table": table(["Chromosome", "Start", "End", "Length (Mb)", "Type", "Cell fraction"], rows),
            "lines": rr.mosaic_mca_lines(d), "attention": bool(d.get("calls"))}


def a_carrier(S, s):
    d = S["carrier"]["data"]
    if d.get("unavailable"):
        return {"lines": rr.carrier_lines(d), "attention": True}
    rows = [[c.get("gene"), c.get("variant"), c.get("genotype"), c.get("significance"), c.get("stars"),
             c.get("inheritance"), c.get("condition")] for c in d.get("carriers") or []]
    two = (d.get("two_hits") or {}).get("genes") or []
    return {"facts": [["Genes on the carrier list", d.get("genes_on_list")], ["Carrier findings", len(rows)],
                      ["Genes with two hits", ", ".join(g if isinstance(g, str) else g.get("gene", ".") for g in two)
                       or "none"]],
            "table": table(["Gene", "Variant", "Genotype", "Significance", "Stars", "Inheritance", "Condition"], rows),
            "lines": rr.two_hit_lines(d.get("two_hits")), "attention": bool(two)}


def a_couple(S, s):
    d = S["couple"]["data"]
    rows = [[lab, x.get("gene"), x.get("condition"), x.get("risk"), x.get("own"), x.get("partner_findings")]
            for lab, key in (("Both carry", "shared"), ("X-linked", "x_linked")) for x in d.get(key) or []]
    return {"facts": [["Partner", d.get("partner")], ["Genes where both carry a change", len(d.get("shared") or [])],
                      ["X-linked findings", len(d.get("x_linked") or [])],
                      ["Genes where only one partner carries a change",
                       f"{d.get('one_partner', 0)}, unless the other partner carries something the screen does not see"]],
            "table": table(["Result", "Gene", "Condition", "Chance for a child", "This person", "Partner"], rows)}


# The genes the copy-number dots draw (the carrier screen's copy-number genes): SMN1 and SMN2, the two
# alpha-globin genes and the RCCX module of CYP21A2.
COPY_GENES = ("SMN1", "SMN2", "HBA1", "HBA2", "CYP21A2")


def a_paralogs(S, s):
    lines, facts = [], []
    if (S.get("parascopy") or {}).get("state") in ("ok", "stale"):
        d = S["parascopy"]["data"]
        lines += rr.parascopy_lines(d) + rr.parascopy_locus_lines(d) + rr.cn_interval_lines(d.get("cn_intervals"))
        if d.get("smn1_marker"):
            # the silent-carrier marker beside the SMN1 psCN
            facts = [["SMN1 copies", (d.get("smn1") or {}).get("text", ".")],
                     ["SMN1 g.27134T>G silent-carrier marker", rr.marker_state(d["smn1_marker"])]]
    rows, dots = [], []
    if (S.get("paralogs") or {}).get("state") in ("ok", "stale"):
        for g in S["paralogs"]["data"].get("genes") or []:
            rows.append([g.get("gene"), g.get("copy_number"), g.get("status"), g.get("detail")])
            if str(g.get("gene") or "").split(" ")[0] in COPY_GENES:
                n = num(g.get("copy_number"))
                dots.append({"gene": g.get("gene"), "copies": int(n) if n is not None and n == int(n) and 0 <= n <= 20
                             else None, "status": g.get("status") or "."})
    if (S.get("parascopy") or {}).get("state") in ("ok", "stale"):
        # Parascopy's SMN1 and SMN2 copies at exon 7 (short reads) as dots too; its count is None unless
        # psCN passed its filter, and then the dots say why there is none.
        d = S["parascopy"]["data"]
        for g in ("smn1", "smn2") if d.get("rows") else ():  # rows: step 35's SMN table had rows
            x = d.get(g) or {}
            if not x:
                continue
            n = num(x.get("cn"))
            whole = n is not None and n == int(n) and 0 <= n <= 20
            dots.append({"gene": f"{g.upper()} (Parascopy)", "copies": int(n) if whole else None,
                         "status": "called" if whole else (x.get("text") or "not called")})
    out = {"facts": facts, "lines": lines, "table": table(["Gene", "Copies", "Status", "Detail"], rows)}
    if dots:
        out["widget"] = {"kind": "copies", "genes": dots}
    return out


def a_pharmcat(S, s):
    d = S["pharmcat"]["data"]
    return {"facts": [["PharmCAT report", "Ran"], ["PharmCAT version", d.get("version") or "unknown"],
                      ["Its own HTML report", "beside it" if d.get("html_report") else "not written"]],
            "links": [r"(^|/)pharmcat/[^/]+\.html$"]}


def a_cpic(S, s):
    d = S["cpic"]["data"]
    if d.get("parse_failed"):
        return {"lines": ["The PharmCAT report could not be parsed: see the CPIC report."]
                + [f"WARNING: {w}" for w in d.get("warnings") or []], "attention": True}
    rows = [[g.get("gene"), g.get("diplotype"), f'"{g.get("phenotype")}"' if g.get("phenotype") else ".", g.get("status")]
            for g in d.get("genes") or []]
    return {"facts": [["Genes with a changed-function phenotype", d.get("non_normal", 0)],
                      ["Genes with more than one possible result", d.get("ambiguous", 0)],
                      ["Genes not called", d.get("not_called", 0)]],
            "table": table(["Gene", "Diplotype", "Phenotype (CPIC's words)", "Status"], rows),
            "lines": [f"WARNING: {w}" for w in d.get("warnings") or []]}


def a_cyp2d6(S, s):
    c = s.get("cyp2d6") or {}
    rows = [[rr.CYP2D6_LABELS.get(k, k), v or "not run"] for k, v in (c.get("calls") or {}).items()]
    agree = c.get("agree")
    facts = [["pypgx and Cyrius agree", "yes" if agree else "no" if agree is False else "fewer than two usable calls"]]
    v = rr.cyp2d6_verdict(c)
    if v:
        facts.append(["Step 36", plain(v)])
    lines = rr.pypgx_lines(S["pypgx"]["data"]) if (S.get("pypgx") or {}).get("state") in ("ok", "stale") else []
    return {"facts": facts, "table": table(["Caller", "CYP2D6"], rows), "lines": lines, "attention": agree is False}


def a_hla(S, s):
    d = S["hla"]["data"]
    rows = [[l.get("gene"), rr.hla_text(l)] for l in d.get("loci") or []]
    return {"facts": [["HLA database", (s.get("databases") or {}).get("hla_database") or "unknown"]],
            "table": table(["Gene", "Alleles (T1K quality)"], rows), "lines": rr.hla_derived_lines(d.get("derived"))}


def a_locityper(S, s):
    d = S["locityper"]["data"]
    rows = [[l.get("locus"), l.get("genotype"), l.get("quality"), l.get("filter"), l.get("primary")]
            for l in d.get("loci") or []]
    return {"facts": [["Loci called", f"{d.get('called', 0)} of {len(d.get('loci') or [])}"],
                      ["Passing the authors' filter", d.get("passed", 0)]],
            "table": table(["Locus", "Panel haplotypes", "Quality", "Filter", "Primary caller"], rows)}


def a_ancestry_fine(S, s):
    d = S["ancestry_fine"]["data"]
    rows = [[lab, r.get("group"), rr.share(r.get("fraction")), f"{rr.share(r.get('low'))} to {rr.share(r.get('high'))}"]
            for level, lab in rr.ANCESTRY_FINE_LEVELS for r in rr.ancestry_fine_rows(d, level)]
    return {"table": table(["Level", "Group", "Share", "Bootstrap range"], rows),
            "facts": [["Chromosome painting", rr.ancestry_fine_painted(d)]]}


def a_mt_haplogroup(S, s):
    d = S["haplogroup"]["data"]
    return {"facts": [["Haplogroup", d.get("haplogroup")], ["Two mtDNA haplogroups in the reads",
                                                            rr.contamination_text(d)]],
            "explore": [{"name": x.get("name"), "url": x.get("url")} for x in d.get("explore") or []]}


def a_y(S, s):
    facts, explore, rows = [], [], []
    y = (S.get("y_haplogroup") or {}).get("data") or {}
    t = (S.get("y_tree") or {}).get("data") or {}
    if y:
        facts += [["Haplogroup (Yleaf)", y.get("haplogroup")], ["QC score", y.get("qc_score")]]
        explore = [{"name": x.get("name"), "url": x.get("url")} for x in y.get("explore") or []]
    if t:
        facts.append([f"Placement on {t.get('tree')}", t.get("haplogroup")])
        rows = [[p.get("Level"), p.get("Branch"), p.get("Derived"), p.get("Ancestral"), p.get("No_call")]
                for p in t.get("path") or []]
    st = (S.get("y_str") or {}).get("data") or {}
    if st.get("markers"):
        facts.append(["Y-STR markers read", sum(1 for m in st["markers"] if m.get("count"))])
    out = {"facts": facts, "explore": explore,
           "table": table(["Level", "Branch", "Derived", "Ancestral", "No call"], rows)}
    if t and t.get("path"):
        # The path from the root of the tree to the placement, one branch per level, with its marker counts.
        out["widget"] = {"kind": "ypath", "tree": t.get("tree"), "tip": t.get("haplogroup"),
                         "path": [{"level": p.get("Level"), "branch": p.get("Branch"), "derived": num(p.get("Derived")),
                                   "ancestral": num(p.get("Ancestral")), "no_call": num(p.get("No_call"))}
                                  for p in t["path"]]}
    return out


def a_archaic(S, s):
    d = S["archaic"]["data"]
    return {"facts": [["Neanderthal", f"{d.get('neanderthal_percent')}%"], ["Denisovan", f"{d.get('denisovan_percent')}%"]],
            "lines": rr.archaic_lines(d)}


def a_ancient(S, s):
    d = S["ancient"]["data"]
    rows = [[r.get("rank"), r.get("group"), r.get("dates"), r.get("country"), r.get("distance")]
            for r in (d.get("nearest") or [])[:10]]
    out = {"table": table(["Rank", "Group", "Dates", "Country", "Distance"], rows),
           "facts": [["Sites used", f"{d.get('sites_used')} of {d.get('panel_sites')}"]]}
    svg = d.get("plot_svg")
    if isinstance(svg, str) and svg.strip():
        # collect_summary stores it base64-encoded (the report puts it in a data URL as is); raw SVG is encoded here.
        import base64
        if svg.lstrip().startswith("<"):
            svg = base64.b64encode(svg.encode()).decode()
        out["image"] = {"alt": "PCA of ancient genomes with this genome projected", "src": "data:image/svg+xml;base64," + svg}
    return out


def a_roh(S, s):
    d = S["roh"]["data"]
    return {"facts": [["Segments", d.get("segments")], ["Total", f"{d.get('total_mb')} Mb"],
                      ["Largest", f"{d.get('largest_mb')} Mb"], ["F_ROH", d.get("f_roh")],
                      ["Nearest band", d.get("relationship")]],
            "lines": rr.roh_lines(d), "attention": bool(d.get("upd_flag"))}


def a_traits(S, s):
    d = S["traits"]["data"]
    rows = [[r.get("gene"), r.get("trait"), (r.get("genotype") or ".") + (" (ancestry)" if r.get("ancestry") == "yes" else ""),
             (r.get("reading") or ".") + (f". {r['note']}" if r.get("note") else "")]
            for r in d.get("rows") or [] if r.get("status") == "called" and (r.get("gene") or "").upper() != "APOE"]
    lines = []
    for r in d.get("rules") or []:
        # The report's genotype rules: several sites read together, with phase and APOL1's G2.
        lines.append(f"{r.get('trait')} ({r.get('gene')} {str(r.get('rsids') or '').replace(',', ', ')}): "
                     f"{rr.rule_result(r)}. {r.get('reading')}")
        if r.get("g2"):
            lines.append(f"With G2 (rs71785313): {rr.g2_result(r)}. {r['g2'].get('text')} {r['g2'].get('reading')}")
    out = {"facts": [["Rows read", d.get("called")], ["No genotype at the position", d.get("not_covered")]],
           "table": table(["Gene", "Trait", "Genotype", "Reading"], rows), "lines": lines}
    apoe = d.get("apoe")
    if apoe:
        # Behind its warning, as the report shows it: the page puts it in a closed <details>.
        # collect_summary keeps one dict per APOE row; each reads as the report's text line does.
        rows = apoe if isinstance(apoe, list) else [apoe]
        text = "; ".join(f"{r.get('trait')} ({r.get('rsid')}): {r.get('genotype')}. {r.get('reading')}"
                         if isinstance(r, dict) else cell(r) for r in rows)
        out["hidden"] = {"summary": rr.APOE_WARNING, "text": plain(text)}
    return out


def a_prs(S, s):
    d = S["prs"]["data"]
    if not d.get("adjusted"):
        # The report shows a percentile only when it is ancestry-adjusted.
        return {"facts": [["Scores", len(d.get("scores") or [])]],
                "lines": ["Raw scores only: without the ancestry panel no percentile is shown (docs/25-prs.md)."]}
    rows = [[r.get("condition"), r.get("pgs_id"), rr.prs_percentile(r) or ".", f"{r.get('matched')}/{r.get('total')}",
             rr.prs_effect(r) or "."] for r in d.get("scores") or []]
    return {"table": table(["Condition", "Score", "Percentile", "Variants matched", "Effect"], rows),
            "lines": [rr.prs_note(d)]}


def a_bloodtype(S, s):
    d = S["bloodtype"]["data"]
    rows = [[x.get("system"), x.get("phenotype"), x.get("genotype"), x.get("confidence")] for x in d.get("systems") or []]
    return {"facts": [["ABO", d.get("abo")], ["RhD", d.get("rhd")]],
            "table": table(["System", "Predicted phenotype", "Genotype", "Confidence"], rows)}


def a_telomere(S, s):
    d = S["telomere"]["data"]
    return {"facts": [["Telomere content", f"{d.get('tel_content')} ({rr.TELOMERE_UNIT})"]],
            "lines": rr.telomere_lines(d), "attention": bool(d.get("low_reads"))}


def a_viral(S, s):
    d = S["viral"]["data"]
    rows = [[x.get("virus"), x.get("reads"), x.get("breadth"), x.get("copies_per_cell"), x.get("call")]
            for x in d.get("found") or []]
    return {"facts": [["HHV-6", d.get("hhv6")], ["Genomes searched", d.get("genomes")]],
            "table": table(["Virus", "Reads", "Breadth", "Copies per cell", "Call"], rows)}


def a_consumer(S, s):
    d = S["consumer_export"]["data"]
    return {"facts": [["Sites", d.get("sites")], ["Called", d.get("called")], ["No call", d.get("no_call")]],
            "lines": rr.consumer_export_lines(d)}


ADAPTERS = {
    "qc": a_qc, "callability": a_callability, "variants": a_variants, "clinvar": a_clinvar,
    "clinical-filter": a_clinical_filter, "cpsr": a_cpsr, "slivar": a_slivar, "exomiser": a_exomiser,
    "structural-variants": a_structural, "mei": a_mei, "expansions": a_expansions, "vntyper": a_vntyper,
    "mito": a_mito, "lpa": a_lpa, "chip": a_chip, "mosaic": a_mosaic, "mosaic-mca": a_mosaic_mca,
    "carrier-list": a_carrier, "couple": a_couple, "paralogs": a_paralogs, "pharmcat": a_pharmcat, "cpic": a_cpic,
    "cyp2d6": a_cyp2d6, "hla": a_hla, "locityper": a_locityper, "ancestry-fine": a_ancestry_fine,
    "mt-haplogroup": a_mt_haplogroup, "y-chromosome": a_y, "archaic": a_archaic, "ancient": a_ancient, "roh": a_roh,
    "traits-panel": a_traits, "prs": a_prs, "bloodtype": a_bloodtype, "telomere": a_telomere, "viral": a_viral,
    "consumer-export": a_consumer,
}
# Report charts (bin/report_charts.py) drawn inside a card: chart id -> section id.
CHART_SECTION = {"clinvar": "clinvar", "genome": "clinvar", "mito": "mito", "ancestry": "ancestry-fine", "prs": "prs"}


def caveats(sid):
    return [getattr(rr, n) for n in rr.SECTION_NOTES.get(sid, ()) if isinstance(getattr(rr, n, None), str)]


def section_state(secs):
    st = {x.get("state") or "missing" for x in secs}
    if st & {"failed", "unreadable"}:
        return "failed" if not st & {"ok", "stale"} else "part"
    if "stale" in st:
        return "stale"
    if "ok" in st:
        return "incomplete" if "missing" in st else "ran"
    return "missing"


def _links(patterns, files):
    out = []
    for p in patterns or []:
        rx = re.compile(p)
        out += sorted(f for f in files if rx.search(f) and f not in out)
    return out[:4]


def card(sid, title, keys, summary, files):
    S = summary.get("sections") or {}
    secs = [S.get(k) or {"state": "missing"} for k in keys]
    state = section_state(secs)
    steps = []
    for k, sec in zip(keys, secs):
        st = sec.get("step") or rr.STEP_OF.get(k)
        if st and st not in steps:
            steps.append(st)
    cav = caveats(sid)
    out = {"id": sid, "title": title, "keys": list(keys), "steps": steps, "state": state, "anchor": sid,
           "attention": False, "headline": None, "facts": [], "table": None, "lines": [], "caveat": cav,
           "empty_text": EMPTY.get(sid), "links": [], "notes": [], "flagged": False}
    for k, sec in zip(keys, secs):
        if sec.get("state") in ("failed", "unreadable", "stale"):
            out["notes"].append(f"{sec.get('title') or k}: {plain(rr.stale_note(sec))}")
            out["attention"] = True
    if state == "missing":
        return out
    tones = []
    for k, sec in zip(keys, secs):
        if sec.get("state") in ("ok", "stale"):
            text, tone = _headline_tone(k, sec, summary)
            if text and not out["headline"] and not text.endswith("finished, see the report"):
                out["headline"] = text
            tones.append(tone)
    if "attention" in tones:
        out["attention"] = True
        out["flagged"] = True
    if out["headline"]:
        # The step headline carries its caveat (the Progress tab shows no other); the card prints it below.
        for c in sorted({plain(c) for c in cav}, key=len, reverse=True):
            if out["headline"].endswith(" " + c):
                out["headline"] = out["headline"][:-len(c) - 1].rstrip()
                break
    fn = ADAPTERS.get(sid)
    if fn and state not in ("failed",):
        try:
            got = fn(S, summary) or {}
        except (KeyError, TypeError, ValueError, AttributeError, IndexError, ZeroDivisionError) as e:
            got = {"lines": [f"This card could not be drawn ({type(e).__name__}); the report shows the section."]}
        out["facts"] = [[plain(a), plain(cell(b))] for a, b in got.get("facts") or []]
        out["table"] = got.get("table")
        notes = {plain(c) for c in cav} | {plain(n) for n in rr.caveat_notes()}
        out["lines"] = [plain(x) for x in got.get("lines") or [] if plain(x) and plain(x) not in notes
                        and plain(x) != out["headline"]]
        out["attention"] = out["attention"] or bool(got.get("attention"))
        out["flagged"] = out.get("flagged", False) or bool(got.get("attention"))
        out["links"] = _links(got.get("links"), files)
        for k in ("widget", "explore", "image", "hidden", "viewer"):
            if got.get(k):
                out[k] = got[k]
    return out


# --- the Journal: one level per result, the findings, the lead ---------------------------------
# The live site reads each person's results as a report in ten numbered sections (docs/design/journal/).
# A finding is one row the summary already holds (a ClinVar hit, a carrier gene, a changed-function
# medicine gene) with its level on the one scale every page uses. Nothing here grades anything the
# report does not: each rule below names the field it reads, and a card the report flags that no rule
# turns into a finding still becomes one, at "review", worded by its own headline.

LEVELS = ("act", "medicine", "review", "carrier", "trust", "info")
LEVEL_WORDS = {"act": "Take to a doctor", "medicine": "Tell your prescriber", "review": "Worth a look",
               "carrier": "For family planning", "trust": "Data quality", "info": "For information"}
# (number, key, title, question). The live run uses the same numbers, so a section being written sits
# where it will sit in the finished report.
JOURNAL = (
    (1, "act", "To act on", "Is there anything to take to a doctor?"),
    (2, "medicines", "Medicines", "Which medicines should a prescriber know about?"),
    (3, "family", "For your family", "What could I pass on to a child?"),
    (4, "look", "Worth a look", "What deserves a conversation, but is not a diagnosis?"),
    (5, "trust", "Can I trust this?", "This is about the data, not about the person."),
    (6, "origins", "Origins", "Where do my ancestors come from?"),
    (7, "traits", "Traits and scores", "What do studies say about traits and risk scores?"),
    (8, "assessed", "Assessed, no finding", "What each step looked at, and what it does not see."),
    (9, "cannot", "What it cannot see", "What no step of this run looks at."),
    (10, "sources", "Sources and files", "The databases, the report and the files behind this page."),
)
LEVEL_SECTION = {"act": 1, "medicine": 2, "carrier": 3, "review": 4, "trust": 5}
CHAPTER_SECTION = {"quality": 5, "clinical": 8, "carrier": 3, "medicines": 2, "ancestry": 6, "traits": 7, "other": 10}
# What to do, by level: the same few steps for every finding of a level, never advice about one gene.
PLAN = {
    "act": ["Bring this page to a doctor or a genetic counsellor.",
            "Ask for a clinical test in a certified laboratory to confirm it: this pipeline is research software.",
            "Ask what it means for parents, brothers, sisters and children: close relatives may carry it too."],
    "medicine": ["Show this to whoever prescribes for you, and to your pharmacist.",
                 "Do not change a medicine on your own: the guidance is written for the prescriber.",
                 "Ask for a clinical pharmacogenetic test before a prescription changes."],
    "review": ["Read the evidence and its caveat below.",
               "Talk it over with a doctor or a genetic counsellor if it bears on your health or a pregnancy.",
               "A measurement or a clinical test settles it: this result alone is not a diagnosis."],
    "carrier": ["A carrier is not affected by the condition.",
                "It matters when a partner carries a change in the same gene: compare with the couple screen.",
                "Relatives planning a pregnancy may carry it too."],
    "trust": ["Read nothing else in this report until this is understood.",
              "Check the sample's label and the run's inputs, then run it again."],
}
STR_WORDS = {"pre_mutation": "premutation range", "full_mutation": "full-mutation range"}
CARRIER_NOTE = "A carrier has one changed copy of a gene where it takes two to cause the condition."


def slug(*parts):
    """A finding's id: lower case, digits and dashes, at most 80 characters (the route's own pattern)."""
    raw = "-".join(str(p) for p in parts if p not in (None, ""))
    t = re.sub(r"[^a-z0-9]+", "-", raw.lower()).strip("-")
    if len(t) > 80:  # a long allele: cut, and keep two ids apart by a hash of the whole key
        t = t[:71].rstrip("-") + "-" + hashlib.sha1(raw.encode()).hexdigest()[:8]
    return t or "finding"


FINDING_ID = re.compile(r"^[a-z0-9][a-z0-9-]{0,79}$")


def _f(level, card_id, key, title, text, gene=None, detail=None, evidence=(), viewer=None, af=None):
    return {"id": slug(*key), "level": level, "card": card_id, "title": plain(title), "text": plain(text),
            "gene": gene, "detail": plain(detail) if detail else None,
            "evidence": [[plain(a), plain(cell(b))] for a, b in evidence if b not in (None, "", ".")],
            "viewer": viewer, "af": af}


def _data(S, k):
    sec = S.get(k) or {}
    if not isinstance(sec, dict) or sec.get("state") not in ("ok", "stale"):
        return {}
    d = sec.get("data") or {}
    # Data that is not a mapping gives no rows here; the section's card says it could not be drawn (card()).
    return d if isinstance(d, dict) else {}


def _items(v):
    """A summary's list of rows; a lone row that is not in a list is read as one row (and fails as one)."""
    return v if isinstance(v, list) else [v] if v else []


class _Row:
    """with _Row(out, card[, level]): one rule turning one row of the summary into a finding. A row the rule
    cannot read (a field of the wrong type, a stars value that is not a number) costs that row only: it
    becomes a finding on its card that says the row could not be read, and every other row and card still
    renders. Its words name no gene or variant: the row is unreadable. Its level is the highest the rule
    can give (level "act" for a rule that can give one): a child's page veils a card that holds a finding
    to act on, and the card's table may still draw the row the rule could not read."""

    def __init__(self, out, card_id, level="review"):
        self.out, self.card_id, self.level = out, card_id, level

    def __enter__(self):
        return self

    def __exit__(self, kind, err, tb):
        if kind is None or not issubclass(kind, (KeyError, TypeError, ValueError, AttributeError, IndexError,
                                                 ZeroDivisionError)):
            return False
        self.out.append(_f(self.level, self.card_id, ("unreadable", self.card_id, len(self.out)),
                           "A row of this section could not be read",
                           f"One row of this section holds a value of a kind the page cannot read ({kind.__name__}). "
                           f"The other rows show as usual; the report and its files have the row as the step wrote it."))
        return True


def _copies(n):
    return "two copies" if n == "hom" else "one copy" if n == "het" else (n or "genotype not given")


def finding_rows(summary):
    """[finding] from the rows of the summary, each with its level, in the order of the rules."""
    S = summary.get("sections") or {}
    out = []
    # The genes the later rules leave to an earlier one. A row that is not a mapping names no gene here;
    # its own rule below turns it into the one unreadable finding.
    sf = _data(S, "clinical").get("acmg_sf") or {}
    sf = sf if isinstance(sf, dict) else {}
    def gene(x):
        g = x if isinstance(x, str) else x.get("gene") if isinstance(x, dict) else None
        return g if isinstance(g, str) else None
    sf_genes = {gene(x) for x in _items(sf.get("clinvar_hits")) if isinstance(x, dict)}
    car = _data(S, "carrier")
    carrier_genes = {gene(c) for c in _items(car.get("carriers")) if isinstance(c, dict)}
    th = car.get("two_hits") if isinstance(car.get("two_hits"), dict) else {}
    two_hit = {gene(g) for g in _items(th.get("genes"))} - {None}
    # 1. ACMG secondary findings: the gene's own reporting rule already said "report". The report's own
    # rows (rr.sf_rows): one it cannot read costs this rule, never the rest.
    sf_list = []
    with _Row(out, "clinical-filter", "act"):
        sf_list = rr.sf_rows(sf) if sf and not sf.get("unavailable") else []
    for r in sf_list:
        with _Row(out, "clinical-filter", "act"):
            out.append(_f("act", "clinical-filter", ("sf", r[0], r[1]), f"Secondary finding in {r[0]}",
                          f"A change in {r[0]} that the ACMG secondary-findings list reports: {r[3]}.",
                          gene=r[0], detail=f"{r[0]} {r[1]} {r[2] if r[2] != '.' else ''}".strip(),
                          evidence=[["Variant", r[1]], ["Genotype", r[2]], ["Why it is reported", r[3]]]))
    # 2. Cancer predisposition, CPSR's pathogenic and likely pathogenic rows.
    for x in _items(_data(S, "cpsr").get("plp")):
        with _Row(out, "cpsr", "act"):
            out.append(_f("act", "cpsr", ("cpsr", x.get("gene"), x.get("variant")),
                          f"Cancer predisposition ({x.get('gene')})",
                          f"CPSR grades this change in {x.get('gene')} as {x.get('final')}: it can raise the risk of some "
                          f"cancers. It is not a diagnosis of cancer.", gene=x.get("gene"),
                          detail=f"{x.get('gene')} {x.get('variant')}",
                          evidence=[["Class", x.get("final")], ["ClinVar", x.get("clinvar")],
                                    ["Call", "low confidence" if x.get("low_confidence") else None]]))
    # 3. Two hits in one carrier-list gene: possibly both copies changed.
    for g in sorted(x for x in two_hit if x):
        with _Row(out, "carrier-list", "act"):
            out.append(_f("act", "carrier-list", ("two-hits", g), f"Two changes in {g}",
                          f"Two changes were found in {g}, a gene where it takes two to cause the condition. Whether they "
                          f"sit on different copies decides whether it matters.", gene=g))
    # 4. ClinVar P/LP hits not already a secondary finding or a carrier row.
    for h in _items(_data(S, "clinvar").get("hits"))[:50]:
        with _Row(out, "clinvar"):
            g = h.get("gene")
            if g in sf_genes or g in carrier_genes or g in two_hit:
                continue
            var = f"{h.get('chrom')}:{h.get('pos')} {h.get('ref')}>{h.get('alt')}"
            out.append(_f("review", "clinvar", ("clinvar", g, h.get("chrom"), h.get("pos"), h.get("ref"), h.get("alt")),
                          f"{h.get('significance') or 'Pathogenic'} change in {g}",
                          f"ClinVar lists this change in {g} as {str(h.get('significance') or '').lower()}; this genome "
                          f"has {_copies(h.get('genotype'))}. Whether it matters depends on how the gene is inherited: "
                          f"a genetic counsellor reads it with the rest of this report.",
                          gene=g, detail=f"{g} {var} {h.get('genotype') or ''}".strip(),
                          evidence=[["Significance", h.get("significance")], ["ClinVar review", "★" * int(h.get("stars") or 0) or "0 stars"],
                                    ["Genotype", h.get("genotype")], ["Frequency", rr.af_text(h) or "frequency unknown"],
                                    ["Call", clinvar_call(h)]],
                          viewer=f"{h.get('chrom')}:{h.get('pos')}", af=h.get("af")))
    # 5. Large deletions and duplications AnnotSV ranks class 4 or 5.
    for x in _items(_data(S, "annotsv").get("rows")):
        with _Row(out, "structural-variants", "act"):
            genes = ", ".join(x.get("genes") or []) or "no gene"
            out.append(_f("act" if str(x.get("acmg_class")) == "5" else "review", "structural-variants",
                          ("sv", x.get("sv")), f"Large {x.get('type') or 'structural change'} over {genes}",
                          f"A {x.get('type') or 'structural change'} of {abs(int(num(x.get('length')) or 0)):,} bases that "
                          f"AnnotSV ranks class {x.get('acmg_class')}. Short reads can misplace large changes: confirm it.",
                          gene=(x.get("genes") or [None])[0], detail=x.get("sv"),
                          evidence=[["Class", x.get("acmg_class")], ["Callers", rr.annotsv_support_text(x)],
                                    ["Genotype", x.get("gt")]], viewer=str(x.get("sv") or "").replace(" ", "") or None))
            # Every gene it covers, not only the first: a closed one closes the rows that name any of them (closable()).
            gl = x.get("genes")
            out[-1]["genes"] = [str(g) for g in gl if g] if isinstance(gl, list) else []
    # 6. Repeat expansions outside Stranger's "normal", or in an intermediate zone.
    for l in _items(_data(S, "expansions").get("loci")):
        with _Row(out, "expansions", "act"):
            st = [x for x in (l.get("status") or "").split(",") if x and x != "normal"]
            if not st and not l.get("zone"):
                continue
            full = any("full" in x for x in st)
            lead_st = next((x for x in st if "full" in x), st[0]) if st else None
            what = (STR_WORDS.get(lead_st, lead_st.replace("_", " ")) if lead_st
                    else f"{l.get('zone_label') or 'intermediate'} range")
            count = rr.str_count(l) if l.get("repeat_count") else "not given"
            out.append(_f("act" if full else "review", "expansions", ("str", l.get("locus")),
                          f"{l.get('locus')} repeat in the {what}",
                          f"{l.get('locus')} has {count} repeats: Stranger reads it as \"{what}\". "
                          f"{'See a clinical geneticist.' if full else 'It may grow when passed on; worth a word with a genetic counsellor.'}",
                          gene=l.get("locus"), detail=f"{l.get('locus')} {count} repeats",
                          evidence=[["Repeats", count], ["95% interval", l.get("ci")],
                                    ["Pathologic from", l.get("pathologic_min")],
                                    ["Sizing", l.get("allele_sizing") or l.get("sizing")]]))
    with _Row(out, "vntyper", "act"):
        v = _data(S, "vntyper")
        if v.get("status") in ("positive", "possible"):
            out.append(_f("act" if v["status"] == "positive" else "review", "vntyper", ("muc1", v["status"]),
                          f"MUC1 kidney gene: {v['status']}", f"VNtyper reads the MUC1 repeat as {v['status']}; it needs a "
                          f"clinical confirmation.", gene="MUC1", evidence=[["Confidence", v.get("confidence")]]))
    with _Row(out, "lpa"):
        lp = _data(S, "lpa")
        if lp.get("flag") == "high":
            out.append(_f("review", "lpa", ("lpa",), "Lp(a) likely high",
                          "The KIV-2 count suggests a high Lp(a); one blood measurement settles it.", gene="LPA",
                          evidence=[["KIV-2 repeats", lp.get("kiv2_copies")], ["Class", lp.get("class")]]))
    for c in _items(_data(S, "chip").get("calls")):
        with _Row(out, "chip"):
            out.append(_f("review", "chip", ("chip", c.get("gene"), c.get("protein_change")),
                          f"A blood clone in {c.get('gene')} (CHIP)",
                          f"{c.get('gene')} {c.get('protein_change') or ''} is carried by a share of the blood cells (VAF "
                          f"{c.get('vaf')}). It becomes common with age; most people with it never get a blood cancer.",
                          gene=c.get("gene"), evidence=[["VAF", c.get("vaf")], ["Reads", f"{c.get('alt_reads')} of {c.get('depth')}"]]))
    for c in _items(_data(S, "mosaic").get("calls")):
        with _Row(out, "mosaic-mca"):
            out.append(_f("review", "mosaic-mca", ("mca", c.get("chrom"), c.get("start")),
                          f"A {c.get('type') or 'change'} on {c.get('chrom')} in some blood cells",
                          f"MoChA sees a {c.get('type') or 'change'} of {c.get('length_mb')} Mb on {c.get('chrom')} in a share "
                          f"of the blood cells ({c.get('cell_fraction')}).",
                          evidence=[["Where", f"{c.get('chrom')}:{c.get('start')}-{c.get('end')}"],
                                    ["Cell fraction", c.get("cell_fraction")]]))
    for x in _items(_data(S, "exomiser").get("genes"))[:3]:
        with _Row(out, "exomiser"):
            out.append(_f("review", "exomiser", ("exomiser", x.get("gene")), f"Symptom ranking: {x.get('gene')} ranks {x.get('rank')}",
                          "A ranking to bring to a clinical geneticist, not a diagnosis.", gene=x.get("gene"),
                          evidence=[["Score", x.get("exomiser_score")], ["Inheritance", x.get("moi")]]))
    # 7. Carriers: one changed copy in a gene of the carrier list, and the couple screen.
    for c in _items(car.get("carriers")):
        with _Row(out, "carrier-list", "act"):
            if c.get("gene") in two_hit:
                continue
            gt = str(c.get("genotype") or "")
            if gt in ("hom", "hemizygous"):
                # Both copies changed, or the only copy (X in a male, or inside a deletion): not a carrier.
                out.append(_f("act", "carrier-list", ("carrier", c.get("gene"), c.get("variant")),
                              f"{c.get('condition') or c.get('gene')}: {'two changed copies' if gt == 'hom' else 'the only copy changed'}",
                              f"This genome has {'two changed copies' if gt == 'hom' else 'no unchanged copy'} of "
                              f"{c.get('gene')}, a gene of {c.get('condition') or 'a condition on the carrier list'} "
                              f"({c.get('inheritance') or 'inheritance not given'}). That is not the pattern of a carrier: "
                              f"a doctor or genetic counsellor reads it.",
                              gene=c.get("gene"), detail=f"{c.get('gene')} {c.get('variant') or ''} {gt}".strip(),
                              evidence=[["Genotype", gt], ["Significance", c.get("significance")], ["Stars", c.get("stars")],
                                        ["Inheritance", c.get("inheritance")]]))
                continue
            out.append(_f("carrier", "carrier-list", ("carrier", c.get("gene"), c.get("variant")),
                          f"{c.get('condition') or c.get('gene')} carrier" if c.get("condition") else f"Carrier in {c.get('gene')}",
                          f"{CARRIER_NOTE} This genome has one changed copy of {c.get('gene')}: a carrier, not affected.",
                          gene=c.get("gene"), detail=f"{c.get('gene')} {c.get('variant') or ''}".strip(),
                          evidence=[["Significance", c.get("significance")], ["Stars", c.get("stars")],
                                    ["Inheritance", c.get("inheritance")]]))
    cp = _data(S, "couple")
    for lab, key in (("Both partners carry", "shared"), ("X-linked", "x_linked")):
        for x in _items(cp.get(key)):
            with _Row(out, "couple"):
                out.append(_f("carrier", "couple", ("couple", key, x.get("gene")), f"{lab} {x.get('gene')}"
                              if key == "shared" else f"X-linked finding in {x.get('gene')}",
                              f"{x.get('condition') or x.get('gene')}: {x.get('risk') or 'see the couple screen'} for each "
                              f"pregnancy, with {cp.get('partner')}.", gene=x.get("gene"),
                              evidence=[["Chance for a child", x.get("risk")], ["Partner", cp.get("partner")]]))
    with _Row(out, "paralogs"):
        smn1 = [g for g in _data(S, "paralogs").get("genes") or [] if g.get("gene") == "SMN1" and g.get("status") == "called"]
        one = smn1 and str(smn1[0].get("copy_number")) == "1"
        by = "Paraphase"
        if not smn1:  # no Paraphase call: Parascopy's resolved count from short reads
            cn = ((_data(S, "parascopy").get("smn1") or {}).get("cn"))
            one, by = cn == 1, "Parascopy"
        if one:
            out.append(_f("carrier", "paralogs", ("smn1",), "Spinal muscular atrophy carrier (SMN1)",
                          "One copy of SMN1 where most people have two: a carrier of spinal muscular atrophy, not "
                          "affected.", gene="SMN1", evidence=[["SMN1 copies", "1"], ["Called by", by]]))
    # 8. Medicines: each gene whose called phenotype changes function.
    for g in _items(_data(S, "cpic").get("genes")):
        with _Row(out, "cpic"):
            if g.get("status") not in ("non-normal", "ambiguous"):
                continue
            amb = g["status"] == "ambiguous"
            out.append(_f("medicine", "cpic", ("pgx", g.get("gene")),
                          f"{g.get('gene')}: {g.get('phenotype')}" if not amb else f"{g.get('gene')}: more than one possible result",
                          f"{g.get('gene')} {g.get('diplotype')} reads as \"{g.get('phenotype')}\" in CPIC's words. The drugs "
                          f"it bears on are in the CPIC report; tell whoever prescribes.", gene=g.get("gene"),
                          detail=f"{g.get('gene')} {g.get('diplotype')}",
                          evidence=[["Diplotype", g.get("diplotype")], ["Phenotype", g.get("phenotype")], ["Status", g.get("status")]]))
    # 9. Trust: a problem with the data itself.
    with _Row(out, "qc"):
        q = _data(S, "sample_qc")
        if q.get("sex_check") == "mismatch":
            out.append(_f("trust", "qc", ("sex",), "Declared and measured sex differ",
                          "The sex the reads show differs from the one declared: a swapped or mislabelled sample, or a real "
                          "sex-chromosome difference.", evidence=[["Declared", q.get("declared_sex")], ["From the reads", q.get("inferred_sex")]]))
        if q.get("contamination") == "warn":
            out.append(_f("trust", "qc", ("contamination",), "Another person's DNA may be mixed in",
                          rr.freemix_text(q), evidence=[["FREEMIX", q.get("freemix")]]))
        if q.get("same_person_as"):
            out.append(_f("trust", "qc", ("same-person",), f"The same person as {q['same_person_as']}",
                          "Two runs of this server hold the same genome.", evidence=[["Matches", q["same_person_as"]]]))
    seen, uniq = set(), []
    for x in out:
        if x["id"] not in seen:
            seen.add(x["id"])
            uniq.append(x)
    return uniq


def level_of(card, found):
    """A card's level: its highest finding's, else data quality for the quality chapter, else review
    when the report flags it, else information. A failed or missing card has no level (None)."""
    if card.get("state") in ("missing", "failed"):
        return None
    mine = [f["level"] for f in found if f["card"] == card["id"]]
    for lv in LEVELS:
        if lv in mine:
            return lv
    if card.get("chapter") == "quality":
        return "trust"
    return "review" if card.get("flagged") else "info"


def journal(chapters, summary):
    """findings, levels and the ten sections of one person's report, from the cards build() made."""
    found = finding_rows(summary)
    cards = []
    for ch in chapters:
        for c in ch["sections"]:
            c["chapter"] = ch["key"]
            cards.append(c)
    by_id = {c["id"]: c for c in cards}
    have = {f["card"] for f in found}
    # A card the report flags (its headline or its own rule), for which no row rule found a finding,
    # still becomes one: worded by its headline, at review.
    for c in cards:
        if c.get("flagged") and c["id"] not in have and c.get("state") not in ("missing", "failed"):
            lv = "trust" if c["chapter"] == "quality" else "review"
            found.append(_f(lv, c["id"], ("card", c["id"]), c["title"], c.get("headline") or
                            f"The report flags this section: {c['title']}.", evidence=(c.get("facts") or [])[:4]))
    for f in found:
        c = by_id.get(f["card"]) or {}
        f["section"] = LEVEL_SECTION[f["level"]]
        f["steps"] = c.get("steps") or []
        f["caveat"] = c.get("caveat") or []
        f["anchor"] = c.get("anchor") or f["card"]
        f["plan"] = PLAN[f["level"]]
        f["level_words"] = LEVEL_WORDS[f["level"]]
    for c in cards:
        c["level"] = level_of(c, found)
        lv = [f["level"] for f in found if f["card"] == c["id"]]
        top = next((x for x in LEVELS if x in lv), None)
        home = CHAPTER_SECTION.get(c["chapter"], 10)
        if c["chapter"] == "clinical" and top in LEVEL_SECTION:
            home = LEVEL_SECTION[top]
        c["section"] = home
        c["findings"] = [f["id"] for f in found if f["card"] == c["id"]]
    for f in found:
        f["related"] = [g["id"] for g in found if g["id"] != f["id"] and
                        ((f["gene"] and g["gene"] == f["gene"]) or g["card"] == f["card"])][:6]
    sections = []
    for n, key, title, question in JOURNAL:
        sections.append({"n": n, "key": key, "title": title, "question": question,
                         "findings": [f["id"] for f in found if f["section"] == n],
                         "cards": [c["id"] for c in cards if c["section"] == n]})
    levels = {lv: sum(1 for f in found if f["level"] == lv) for lv in LEVELS if lv != "info"}
    levels["info"] = sum(1 for c in cards if c["level"] == "info")
    return found, levels, sections


def _names(found, n=3):
    t = [f["title"] for f in found[:n]]
    more = len(found) - len(t)
    s = ", ".join(t[:-1]) + (" and " if len(t) > 1 else "") + t[-1] if t else ""
    return s + (f", and {more} more" if more > 0 else "")


def _count(n, one, many):
    words = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine"}
    return f"{words.get(n, n)} {one if n == 1 else many}"


def screened(summary):
    """True when a clinical screen (ClinVar, the clinical filter or CPSR) has a result: only then can the
    lead or the home page say that no finding needs a doctor."""
    S = summary.get("sections") or {}
    return any((S.get(k) or {}).get("state") in ("ok", "stale") for k in ("clinvar", "clinical", "cpsr"))


def lead(found, levels, summary):
    """The report's lead paragraph: [[text, bold]] segments (a third item "act" marks the names of the
    findings to act on), plain sentences that name the findings by level. It never says a section is fine: with no finding of a level it says what was assessed."""
    S = summary.get("sections") or {}
    by = {lv: [f for f in found if f["level"] == lv] for lv in LEVELS}
    out = []
    if by["trust"]:
        out += [["Read the data-quality section first: ", True], [_names(by["trust"]) + ". ", False]]
    if by["act"]:
        n = len(by["act"])
        out += [["Of everything this genome was screened for, ", False],
                [_count(n, "finding should", "findings should") + " go to a doctor", True],
                [": " + _names(by["act"]) + ". ", False, "act"]]  # "act": the names a child's page keeps closed
    elif screened(summary):
        out += [["No finding in what these steps assessed needs a doctor; ", True],
                ["sections 8 and 9 say what they do not see. ", False]]
    else:  # no clinical screen has a result: saying nothing needs a doctor would read as an all-clear
        out += [["The clinical screens have no result yet; ", True],
                ["section 9 lists the steps that did not run. ", False]]
    if by["medicine"]:
        genes = sorted({f["gene"] for f in by["medicine"] if f["gene"]})
        out += [[_count(len(by["medicine"]), "medicine result", "medicine results").capitalize(), True],
                [f" bear{'s' if len(by['medicine']) == 1 else ''} on how drugs work ({', '.join(genes)}): tell whoever prescribes. ", False]]
    if by["carrier"]:
        out += [[f"A carrier of {_names(by['carrier'], 2)}", True],
                [": carriers are not affected; it matters when a partner carries the same gene. ", False]]
    if by["review"]:
        out += [[_count(len(by["review"]), "result is", "results are").capitalize() + " worth a look", True],
                [": " + _names(by["review"]) + ". ", False]]
    cov = _data(S, "coverage").get("mean_depth")
    q = _data(S, "sample_qc")
    checked = q.get("contamination") == "ok" and q.get("sex_check") == "ok"  # both ran and raised nothing
    if not by["trust"] and isinstance(cov, (int, float)):
        out.append([f"Mean depth {cov:.1f}×" + ("; the identity and contamination checks raised no flag." if checked
                    else "; the identity and contamination checks have no result here: section 5 says why."), False])
    if out:
        out[-1][0] = out[-1][0].rstrip()
    return out


# --- the veil's tags -----------------------------------------------------------------------------
# A child's page keeps each result to act on closed (assets/live/person.js veiled()). The page closes a
# finding's own card with it, but other places can name the same gene or variant: any card's table and
# lines, the charts' marks and their titles, another finding (a symptom ranking), the lead, the search index. Each such row carries the id of the
# result it belongs to, so the page closes it with that result and opens it with it. The rule: a row
# belongs to a closed result when it names the result's gene, or its position falls inside the result's
# position or span; a row that is the result's own position or span belongs to that result first, else
# the first such result in finding order. The veil lives in the browser: the files still hold every row.
# The cards whose rows are a fixed panel: the same genes or loci for everyone, known before the run. A row or
# dot of one that closes with a result would tell which gene the result is in by what is missing, so a child's
# page closes each of these cards whole while any result to act on is closed, behind one box with the same words
# whatever the result (results.json veil.panels; assets/live/person.js panelShut()). The one table: a card of a
# fixed panel is listed here and nowhere else.
FIXED_PANELS = {
    "expansions": "the repeat loci of the catalog",
    "vntyper": "MUC1",
    "lpa": "LPA",
    "chip": "the CHIP genes and their detection limits",
    "paralogs": "SMN1, SMN2, HBA1, HBA2, the RCCX module and the other paralog genes",
    "cpic": "the CPIC genes",
    "cyp2d6": "CYP2D6",
    "hla": "the HLA genes",
    "locityper": "the Locityper loci",
    "traits-panel": "the trait panel's genes",
    "bloodtype": "the blood group systems",
}
LOCUS_RE = re.compile(r"(?<![A-Za-z0-9])(?:chr)?([0-9]{1,2}|X|Y|M|MT):([0-9][0-9,]*)(?:-([0-9][0-9,]*))?")
WORD_RE = re.compile(r"[A-Z0-9][A-Z0-9._-]*")


def _loci(text):
    """[(chrom, start, end)] of every chrom:pos or chrom:start-end in TEXT."""
    out = []
    for m in LOCUS_RE.finditer(str(text or "")):
        a = int(m.group(2).replace(",", ""))
        b = int(m.group(3).replace(",", "")) if m.group(3) else a
        out.append(("chr" + m.group(1), min(a, b), max(a, b)))
    return out


def closable(found):
    """[(id, gene, chrom, start, end)] of the findings a child's page keeps closed (level act), one per gene the
    finding names (its genes, a structural variant's every gene, else its gene; None for none): the position or
    span from the finding's viewer locus or its detail, else None (its genes alone). A row that names any of a
    result's genes belongs to that result."""
    out = []
    for f in found:
        if f.get("level") != "act":
            continue
        loc = (_loci(f.get("viewer")) or _loci(f.get("detail")) or [(None, None, None)])[0]
        genes = []
        for g in list(f.get("genes") or []) + [f.get("gene")]:
            g = str(g or "").upper()
            if g and g not in genes:
                genes.append(g)
        for g in genes or [None]:
            out.append((f["id"], g) + loc)
    return out


def closable_of(summary):
    """closable() of the findings build() gives SUMMARY, without building the cards."""
    return closable(finding_rows(summary or {}))


def veil_id(spans, text="", chrom=None, pos=None, end=None):
    """The id of the closed result (spans: closable()) a row belongs to, or None. TEXT is what the row
    shows (its gene names, its loci are read too); chrom, pos and end its own position when it has one."""
    if not spans:
        return None
    locs = []
    if chrom and pos is not None:
        try:
            p = int(pos)
            e = int(end) if end not in (None, "") else p
            c = str(chrom)
            locs.append((c if c.startswith("chr") else "chr" + c, min(p, e), max(p, e)))
        except (TypeError, ValueError):
            pass
    locs += _loci(text)
    found = WORD_RE.findall(str(text or "").upper())
    words = set(found) | {w.rstrip("._-") for w in found}  # "in KCNH2." names KCNH2; HLA-A stays whole
    for fid, _gene, c, a, b in spans:
        if c and (c, a, b) in locs:
            return fid
    for fid, gene, c, a, b in spans:
        if gene and gene in words:
            return fid
        if c and any(lc == c and ls <= b and le >= a for lc, ls, le in locs):
            return fid
    return None


def veil(summary, found, chapters, lead_parts=()):
    """The rows a child's page closes with a result to act on, by the result's id: each card's table rows
    (rows: card -> one id or None per row; the ClinVar card's by each hit's own position), each card's viewer
    links (viewer: card -> one per link), each card's lines (lines: card -> one per line; a known deep
    intronic site is printed whether hit or not), the marks of the report's charts (marks: chart -> {mark id:
    id}), every other finding that names a closed result (findings: finding id -> id; a symptom ranking can
    rank the gene) and the lead's segments that name one (lead: one per segment of LEAD_PARTS, the "act"
    names left out: the page closes those itself), each card's facts (facts: card -> one per fact) and the items
    of its drawing (widget: card -> one per CHIP call or copy-number dot), and the cards of a fixed panel
    (panels: FIXED_PANELS, which close whole while any result is closed). spans is closable() itself: the family
    page (assets/live/family-page.js) reads which results are closed from it, and closes every place of the child
    that could carry one while one is. Empty when nothing can be closed."""
    spans = closable(found)
    if not spans:
        return {}
    out = {"rows": {}, "viewer": {}, "lines": {}, "facts": {}, "widget": {}, "marks": {}, "findings": {}, "lead": [],
           "spans": [list(s) for s in spans],
           "panels": [c["id"] for ch in chapters for c in ch["sections"] if c["id"] in FIXED_PANELS and c.get("state") != "missing"]}
    S = summary.get("sections") or {}
    hits = [h for h in _items(_data(S, "clinvar").get("hits"))[:50]]
    tags = [veil_id(spans, str(h.get("gene") or ""), h.get("chrom"), h.get("pos")) if isinstance(h, dict) else None
            for h in hits]
    if any(tags):
        out["rows"]["clinvar"] = tags
    for ch in chapters:
        for c in ch["sections"]:
            rows = (c.get("table") or {}).get("rows") or []
            if c["id"] != "clinvar" and rows:
                tags = [veil_id(spans, " ".join(str(x) for x in r) if isinstance(r, (list, tuple)) else str(r))
                        for r in rows]
                if any(tags):
                    out["rows"][c["id"]] = tags
            tags = [veil_id(spans, v) for v in c.get("viewer") or []]
            if c["id"] == "clinvar":  # one link per row: a row closed by its gene alone closes its link too
                rt = out["rows"].get("clinvar") or []
                tags = [(rt[i] if i < len(rt) else None) or t for i, t in enumerate(tags)]
            if any(tags):
                out["viewer"][c["id"]] = tags
            tags = [veil_id(spans, ln) for ln in c.get("lines") or []]
            if any(tags):
                out["lines"][c["id"]] = tags
            tags = [veil_id(spans, " ".join(str(x) for x in f)) for f in c.get("facts") or []]
            if any(tags):
                out["facts"][c["id"]] = tags
            # A drawing's items: the CHIP calls (gene and change) and the copy-number dots (gene). The CHIP
            # chart's list of genes with their detection limits is the same panel for everyone: never tagged.
            w = c.get("widget") or {}
            items = w.get("calls") if w.get("kind") == "vaf" else w.get("genes") if w.get("kind") == "copies" else None
            tags = [veil_id(spans, f"{x.get('gene') or ''} {x.get('change') or ''}") for x in items or []]
            if any(tags):
                out["widget"][c["id"]] = tags
    try:
        specs = rr.chart_specs(summary)
    except (KeyError, TypeError, ValueError, AttributeError, IndexError, ZeroDivisionError):
        specs = []
    for sp in specs:
        for m in sp.get("marks") or []:
            if m.get("kind") == "roh":  # a run of homozygosity is no result's source
                continue
            fid = veil_id(spans, f"{m.get('label') or ''} {m.get('detail') or ''}", m.get("chrom"), m.get("pos"),
                          m.get("end"))
            if fid:
                out["marks"].setdefault(sp["id"], {})[m["id"]] = fid
    for f in found:
        if f.get("level") != "act":
            fid = veil_id(spans, f"{f.get('gene') or ''} {f.get('detail') or ''} {f.get('title') or ''} "
                                 f"{f.get('viewer') or ''}")
            if fid:
                out["findings"][f["id"]] = fid
    tags = [None if len(x) > 2 else veil_id(spans, x[0]) for x in lead_parts or ()]
    out["lead"] = tags if any(tags) else []
    return {k: v for k, v in out.items() if v}


def build(sample, summary, files=()):
    """results.json of one person from a summary.json; files: the run's paths the report links (for links)."""
    summary = summary or {}
    files = set(files or ())
    S = summary.get("sections") or {}
    chapters = []
    for cid, ctitle, question, secs in rr.CHAPTERS:
        cards = [card(sid, title, keys, summary, files) for sid, title, keys in secs]
        chapters.append({"key": cid, "title": ctitle, "question": question, "sections": cards,
                         "attention": sum(1 for c in cards if c["attention"])})
    sexc = (S.get("sex_check") or {}).get("data") or {}
    person = {"sample": sample, "declared_sex": (summary.get("run") or {}).get("declared_sex"),
              "inferred_sex": sexc.get("inferred_sex"), "karyotype": sexc.get("karyotype")}
    found, levels, sections = journal(chapters, summary)
    lead_parts = lead(found, levels, summary)
    return {"schema_version": SCHEMA, "sample": sample, "generated_utc": summary.get("generated_utc"),
            "person": person, "chapters": chapters,
            "findings": found, "levels": levels, "journal": sections, "lead": lead_parts,
            "level_words": LEVEL_WORDS, "databases": {k: plain(cell(v)) for k, v in (summary.get("databases") or {}).items()
                                                      if isinstance(v, (str, int, float))},
            "attention": sum(c["attention"] for c in chapters),
            "not_assessed": [plain(x) for x in summary.get("not_assessed") or []],
            "not_run": [plain(x) for x in summary.get("not_run") or []],
            "new_since": summary.get("new_since") or None,
            "charts": CHART_SECTION, "veil": veil(summary, found, chapters, lead_parts)}


def attention_count(summary):
    """The number of cards worth a look, for the home page (cheap: no file read)."""
    return build("", summary)["attention"]


def overview(summary):
    """What the home page shows of one person: the attention count, the counts by level, the findings to
    act on first (at most four, by level) and the lead as plain text. No file read."""
    r = build("", summary)
    order = {lv: i for i, lv in enumerate(LEVELS)}
    order["trust"] = -1  # a problem with the data comes before any finding read from it
    # A finding that names a result to act on (veil: a symptom ranking of its gene) is not listed: a child's tile
    # closes only the results to act on themselves, and the result it names already comes first.
    named = (r.get("veil") or {}).get("findings") or {}
    top = sorted((f for f in r["findings"] if f["level"] != "info" and f["id"] not in named),
                 key=lambda f: order[f["level"]])[:4]
    tags = (r.get("veil") or {}).get("lead") or []
    return {"attention": r["attention"], "levels": r["levels"], "screened": screened(summary),
            "top": [{"id": f["id"], "level": f["level"], "title": f["title"]} for f in top],
            "lead": "".join(x[0] for x in r["lead"]),
            # The same lead with the names of the findings to act on left out, and every other segment that names
            # one of them: a child's tile shows this one.
            "lead_closed": "".join(": closed until opened on purpose. " if x[2:] == ["act"]
                                   else closed_lead(x[0]) if i < len(tags) and tags[i] else x[0]
                                   for i, x in enumerate(r["lead"]))}


def closed_lead(t):
    """A segment of the lead that names a result to act on, without the names (person.js closedLead()): the
    list after a colon, the genes in brackets, or what follows "of"."""
    w = "closed until opened on purpose"
    if t.lstrip().startswith(":"):
        return f": {w}. "
    if re.search(r"\(.*\)", t):
        return re.sub(r"\(.*\)", f"({w})", t)
    if " of " in t:
        return re.sub(r" of .*$", f" result, {w}", t)
    return f"{w}. "


def strings(obj):
    """Every string in a JSON value, with its key path: the wording gate's input."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            for p, s in strings(v):
                yield (k,) + p, s
    elif isinstance(obj, list):
        for v in obj:
            yield from strings(v)
    elif isinstance(obj, str):
        yield (), obj


def wording_problems(obj, allow=("src", "url", "caveat", "summary")):
    """[(path, text)] of the strings that would read as an all-clear. Allowed: data URLs and links, the
    caveats (render_report's constants word for word, checked against them separately), and a quoted
    tool term (CPIC's "Normal Metabolizer", Stranger's "normal")."""
    out = []
    for path, s in strings(obj):
        if path and path[-1] in allow:
            continue
        t = re.sub(r'"[^"]*"', "", s)
        if FORBIDDEN.search(t):
            out.append((".".join(map(str, path)), s[:160]))
    return out
