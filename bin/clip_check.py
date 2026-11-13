#!/usr/bin/env python3
"""clip_check.py: the names bin/collect_summary.py reads of the cDNA read check.

The cDNA read check itself is og cdna check (rust/crates/og-bam/src/cdna.rs,
docs/06-clinvar-screen.md): steps 06 and 17 run it, and it writes
<sample>_clinvar_cdna_check.tsv and <sample>_cpsr_cdna_check.tsv. What is
left here is what the report reads beside those tables: the table's columns,
the note of a flagged hit, CPSR's classes in one spelling and the variant of
a CPSR row. og's own copies of these are in cdna.rs; the parity cases of
og-bam keep the two the same.
"""
import re

FLAG_NOTE = "possible cDNA/transgene contamination: check the reads"
COLUMNS = ["chrom", "pos", "ref", "alt", "alt_reads", "alt_clipped", "ref_reads", "ref_clipped",
           "cdna_check", "cdna_note", "construct_site"]


def norm_class(v):
    """A CPSR or ClinVar class in one spelling: 'Likely_Pathogenic' and 'LP'
    -> 'likely pathogenic'; '' for an empty or NA cell."""
    v = (v or "").replace("_", " ").strip().lower()
    return {"p": "pathogenic", "lp": "likely pathogenic", "b": "benign", "lb": "likely benign",
            "na": "", ".": ""}.get(v, v)


PLP = ("pathogenic", "likely pathogenic")


def cpsr_variant(row):
    """(chrom, pos, ref, alt) of a CPSR row from VAR_ID or GENOMIC_CHANGE;
    None when neither can be read."""
    f = (row.get("VAR_ID") or "").split("_")
    if len(f) >= 4 and f[-3].isdigit():
        return "_".join(f[:-3]), int(f[-3]), f[-2], f[-1]
    m = re.match(r"^(\S+):g\.(\d+)([A-Za-z]+)>([A-Za-z]+)$", row.get("GENOMIC_CHANGE") or "")
    if m:
        return m.group(1), int(m.group(2)), m.group(3), m.group(4)
    return None
