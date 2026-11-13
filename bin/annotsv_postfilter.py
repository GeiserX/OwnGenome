#!/usr/bin/env python3
# remnant: kept for collect_summary and the docs demo; deleted in lane R3
"""annotsv_postfilter.py: what is left of step 05's class 4 and 5 checks once
they moved to `og sv postfilter` (rust/crates/og-steps/src/sv). Only the
names the reports and the docs demo still import are kept: the table's
column names, the FORMAT reader of its sample column, and the join of step
22's consensus support that the reports make when they are written.
`og dev remnants` (tests/parity/remnants.tsv) checks that nothing else
comes back. Standard library only.
"""
import gzip
import os

CONSENSUS_BP = 1000
KEEP = ["AnnotSV_ID", "ID", "SV_chrom", "SV_start", "SV_end", "SV_length", "SV_type", "Gene_name", "ACMG_class",
        "AnnotSV_ranking_score"]
ADDED = ["Caller", "Exon_edge", "Exon_edge_detail", "Population_AF", "Population_source", "Common", "AnnotSV_class",
         "Depth_support", "Rule", "Rule_detail"]


def open_text(path):
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", errors="replace")
    return open(path, errors="replace")


def chrom_key(c):
    """AnnotSV writes chromosomes without 'chr'; step 22 writes them with it."""
    return c[3:] if c.lower().startswith("chr") else c


def sv_type(row):
    t = (row.get("SV_type") or "").upper()
    return "DUP" if t.startswith("DUP") else t


def info_value(info, key):
    for kv in (info or "").split(";"):
        k, _, v = kv.partition("=")
        if k == key:
            return v
    return ""


def sample_format(row):
    """{FORMAT key: value} of the first sample column, the one after FORMAT
    in AnnotSV's header; {} when the row has none."""
    cols = sample_columns(row)
    if not cols:
        return {}
    return dict(zip((row.get("FORMAT") or "").split(":"), (row.get(cols[1]) or "").split(":")))


def sample_columns(row):
    """["FORMAT", <the sample column>] of AnnotSV's header, the two columns
    sample_format reads; [] when the row has no sample column."""
    keys = list(row)
    if "FORMAT" not in keys:
        return []
    i = keys.index("FORMAT") + 1
    return [] if i >= len(keys) or keys[i] == "Annotation_mode" else keys[i - 1:i + 1]


def read_consensus(path):
    """(caller names, records) of step 22's consensus VCF; a record is
    (chrom, pos, end, type, supp, supp_vec)."""
    names, recs = [], []
    with open_text(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            c = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                names = c[9:]
                continue
            if len(c) < 8:
                continue
            t = info_value(c[7], "SVTYPE").upper()
            try:
                pos = int(c[1])
                end = int(info_value(c[7], "END") or pos)
                supp = int(info_value(c[7], "SUPP") or 0)
            except ValueError:
                continue
            recs.append((chrom_key(c[0]), pos, end, "DUP" if t.startswith("DUP") else t, supp,
                         info_value(c[7], "SUPP_VEC")))
    return names, recs


def consensus_support(rows, path):
    """Support and Support_callers for each checked row from the consensus
    at path; None for both when there is no consensus file."""
    if not path or not os.path.isfile(path):
        return [{"Support": None, "Support_callers": None} for _ in rows]
    names, recs = read_consensus(path)
    by = {}
    for rec in recs:
        by.setdefault((rec[0], rec[3]), []).append(rec)
    out = []
    for r in rows:
        best = None
        try:
            start, end = int(r["SV_start"]), int(r["SV_end"])
        except (KeyError, ValueError):
            start = end = None
        if start is not None:
            for rec in by.get((chrom_key(r.get("SV_chrom", "")), sv_type(r)), []):
                d = abs(rec[1] - start) + abs(rec[2] - end)
                if abs(rec[1] - start) <= CONSENSUS_BP and abs(rec[2] - end) <= CONSENSUS_BP and (best is None or d < best[0]):
                    best = (d, rec)
        if best:
            vec = best[1][5]
            callers = [names[i] if i < len(names) else f"caller {i + 1}" for i, x in enumerate(vec) if x == "1"]
            out.append({"Support": best[1][4], "Support_callers": ",".join(callers)})
        else:
            out.append({"Support": 1, "Support_callers": r.get("Caller") or "."})
    return out
