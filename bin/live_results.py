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

Nothing here reads a file or computes anything about genetics; the adapters pick
fields the summary already holds. A section id with no adapter falls back to its
report headline and a link to the report.
"""
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
            "widget": {"kind": "str", "loci": [{k: l.get(k) for k in ("locus", "repeat_count", "status", "normal_max",
                                                                          "pathologic_min", "ci")}
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
    return {"facts": facts, "lines": rr.lpa_lines(d), "attention": d.get("flag") == "high"}


def a_chip(S, s):
    d = S["chip"]["data"]
    rows = [[c.get("gene"), c.get("protein_change"), c.get("vaf"), c.get("alt_reads"), c.get("depth"), c.get("status")]
            for c in d.get("calls") or []]
    return {"facts": [["CHIP calls", len(d.get("calls") or [])],
                      ["Typical detection limit at this depth (VAF)", d.get("median_limit")]],
            "table": table(["Gene", "Change", "VAF", "Alt reads", "Depth", "Status"], rows),
            "lines": rr.chip_lines(d, s.get("sample")), "attention": bool(d.get("calls"))}


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


def a_paralogs(S, s):
    lines = []
    if (S.get("parascopy") or {}).get("state") in ("ok", "stale"):
        d = S["parascopy"]["data"]
        lines += rr.parascopy_lines(d) + rr.parascopy_locus_lines(d) + rr.cn_interval_lines(d.get("cn_intervals"))
    rows = []
    if (S.get("paralogs") or {}).get("state") in ("ok", "stale"):
        rows = [[g.get("gene"), g.get("copy_number"), g.get("status"), g.get("detail")]
                for g in S["paralogs"]["data"].get("genes") or []]
    return {"lines": lines, "table": table(["Gene", "Copies", "Status", "Detail"], rows)}


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
    return {"facts": facts, "explore": explore,
            "table": table(["Level", "Branch", "Derived", "Ancestral", "No call"], rows)}


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
           "empty_text": EMPTY.get(sid), "links": [], "notes": []}
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
        out["links"] = _links(got.get("links"), files)
        for k in ("widget", "explore", "image", "hidden", "viewer"):
            if got.get(k):
                out[k] = got[k]
    return out


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
    return {"schema_version": SCHEMA, "sample": sample, "generated_utc": summary.get("generated_utc"),
            "person": person, "chapters": chapters,
            "attention": sum(c["attention"] for c in chapters),
            "not_assessed": [plain(x) for x in summary.get("not_assessed") or []],
            "not_run": [plain(x) for x in summary.get("not_run") or []],
            "new_since": summary.get("new_since") or None,
            "charts": CHART_SECTION}


def attention_count(summary):
    """The number of cards worth a look, for the home page (cheap: no file read)."""
    return build("", summary)["attention"]


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
