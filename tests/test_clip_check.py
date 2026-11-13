#!/usr/bin/env python3
"""bin/clip_check.py: what is left of it after og cdna check took the read
check (rust/crates/og-bam, its unit tests and parity cases): the names
bin/collect_summary.py reads beside step 06's and step 17's tables.

  - norm_class spells every CPSR and ClinVar class one way;
  - cpsr_variant reads VAR_ID, else GENOMIC_CHANGE, else nothing;
  - COLUMNS and FLAG_NOTE are what og cdna check writes (the header and the
    note of a flagged row of an og-bam parity golden).

Run: python3 tests/test_clip_check.py
"""
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import clip_check as cc  # noqa: E402

FAILS = []


def expect(desc, got, want):
    if got == want:
        print(f"PASS {desc}")
    else:
        print(f"FAIL {desc}: got {got!r}, want {want!r}")
        FAILS.append(desc)


def main():
    expect("norm_class: underscores, case and spaces", cc.norm_class(" Likely_Pathogenic "), "likely pathogenic")
    expect("norm_class: LP", cc.norm_class("LP"), "likely pathogenic")
    expect("norm_class: NA and None are empty", (cc.norm_class("NA"), cc.norm_class(None)), ("", ""))
    expect("norm_class: a class it does not know stays", cc.norm_class("VUS"), "vus")
    expect("PLP", cc.PLP, ("pathogenic", "likely pathogenic"))
    expect("cpsr_variant: VAR_ID", cc.cpsr_variant({"VAR_ID": "13_32340300_G_A"}), ("13", 32340300, "G", "A"))
    expect("cpsr_variant: a contig with an underscore", cc.cpsr_variant({"VAR_ID": "HLA_A_5_C_T"}),
           ("HLA_A", 5, "C", "T"))
    expect("cpsr_variant: GENOMIC_CHANGE when VAR_ID does not read",
           cc.cpsr_variant({"VAR_ID": ".", "GENOMIC_CHANGE": "2:g.47800000A>G"}), ("2", 47800000, "A", "G"))
    expect("cpsr_variant: neither", cc.cpsr_variant({"GENOMIC_CHANGE": "2:g.A>G"}), None)
    golden = os.path.join(REPO, "rust/crates/og-bam/tests/data/reads/want/check-spliced/files/out.tsv")
    with open(golden) as fh:
        lines = [line.rstrip("\n").split("\t") for line in fh]
    expect("COLUMNS: the header og cdna check writes", lines[0], cc.COLUMNS)
    expect("FLAG_NOTE: the note of its flagged rows", {r[9] for r in lines[1:] if r[8] == "flagged"}, {cc.FLAG_NOTE})
    if FAILS:
        print(f"{len(FAILS)} check(s) failed")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
