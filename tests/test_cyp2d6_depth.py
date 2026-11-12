#!/usr/bin/env python3
"""The CYP2D6 depth check (og pgx cyp2d6-bed, -measure, -check), by meaning.

`cargo test -p og-bam` (tests/parity.rs) holds og to the files under
rust/crates/og-bam/tests/data/parity/, byte for byte: what the Python tool
og replaced wrote on each case (want/<case>/), or what og writes where it
fixes that tool on purpose (og/<case>/). This test says what those files
must mean, so a golden rewritten with another verdict fails here even when
og agrees with it.

Five depth tables (tests/fixtures/pgx/depth/<case>.q0.bed and .q1.bed, the
mosdepth region files of all reads and of MAPQ >= 1):

  clean_no_alt          the fixture's CYP2D slice mapped to the no-ALT
                        reference (ALT depth A/B run 37177862776): ok
  multimapped_with_alt  the same reads mapped to the Broad hg38 FASTA with
                        ALT contigs (same run): unreliable, multi-mapped
  deleted_both          CYP2D6 deleted on both copies (*5/*5): almost no reads
                        at the gene, normal flanks: ok, the deletion is real
  deleted_one           one copy deleted: half the flank depth: ok
  no_flank_reads        no reads at all: unreliable

Also: the check file has the status, the message, the four depths and the
MAPQ shares; `cyp2d6-bed` prints CYP2D6 between its two flanks; an input
without a CYP2D6 row, with regions other than those (a zero-length one
included), a region listed twice with two depths, or an output that cannot
be written exits 2 and writes nothing. The MAPQ shares beside the verdict
(`cyp2d6-measure`, then `cyp2d6-check --measure`): on synthetic reads at
CYP2D6, reads in REP6 at MAPQ 0 lower the share over Cyrius's region but not
over the gene body; without the measure file every share it gives reads 'not
measured'; the shares never change the status; a file that is not
measure's exits 2.

Run: python3 tests/test_cyp2d6_depth.py
"""
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARITY = os.path.join(REPO, "rust", "crates", "og-bam", "tests", "data", "parity")
FAILS = []

OK = "CYP2D6 depth is fit for a copy-number call"
UNRELIABLE = ("CYP2D6 copy number unreliable: reads are multi-mapped "
              "(was this BAM aligned to a reference with ALT contigs?)")
NO_DEPTH = ("CYP2D6 copy number unreliable: no reads in the flanks of CYP2D6 "
            "(chr22:42.05-42.25 Mb), so its depth cannot be compared")
ROWS = ["metric", "status", "message", "gene_depth_all", "gene_depth_mapq1", "flank_depth_all",
        "flank_depth_mapq1", "gene_mapq1_share", "gene_mapq20_share", "body_depth_all",
        "body_mapq1_share", "body_mapq20_share"]


def expect(desc, got, want):
    if got == want:
        print(f"PASS {desc}")
    else:
        print(f"FAIL {desc}: got {got!r}, want {want!r}")
        FAILS.append(desc)


def golden(case):
    """og's golden of a case when it has one (a deliberate fix), else the oracle's."""
    d = os.path.join(PARITY, "og", case)
    return d if os.path.isdir(d) else os.path.join(PARITY, "want", case)


def exit_code(case):
    return int(open(os.path.join(golden(case), "exit")).read().strip())


def written(case):
    """The files the case writes, by name."""
    d = os.path.join(golden(case), "files")
    return sorted(os.listdir(d)) if os.path.isdir(d) else []


def table(case, name):
    path = os.path.join(golden(case), "files", name)
    rows = [line.rstrip("\n").split("\t", 1) for line in open(path)]
    return {k: v for k, v in rows}, [k for k, _ in rows]


CASES = [
    ("clean_no_alt", "ok", OK),
    ("multimapped_with_alt", "unreliable", UNRELIABLE),
    ("deleted_both", "ok", OK),
    ("deleted_one", "ok", OK),
    ("no_flank_reads", "unreliable", NO_DEPTH),
]

for case, status, message in CASES:
    for suffix, measured in (("", False), ("-measured", True)):
        c = f"check-{case}{suffix}"
        expect(f"{c}: exit code", exit_code(c), 0)
        rows, order = table(c, "check.tsv")
        expect(f"{c}: status", rows.get("status"), status)
        expect(f"{c}: message", rows.get("message"), message)
        expect(f"{c}: header and the eleven rows", order, ROWS)
        shares = [rows.get(k) for k in ("gene_mapq20_share", "body_depth_all", "body_mapq1_share", "body_mapq20_share")]
        if measured:
            expect(f"{c}: the measure file gives every share", "not measured" in shares, False)
        else:
            expect(f"{c}: without the measure file the other shares are not measured", shares, ["not measured"] * 4)

rows, _ = table("check-clean_no_alt", "check.tsv")
expect("clean_no_alt: the verdict's MAPQ >= 1 share over Cyrius's region (16.79 / 18.53)",
       rows.get("gene_mapq1_share"), "0.91")
expect("clean_no_alt: depths, length-weighted",
       [rows.get(k) for k in ("gene_depth_all", "gene_depth_mapq1", "flank_depth_all", "flank_depth_mapq1")],
       ["18.53", "16.79", "29.68", "29.67"])

bed = open(os.path.join(golden("cyp2d6-bed"), "stdout")).read().splitlines()
expect("cyp2d6-bed: CYP2D6 between two flanks", [line.split("\t")[3] for line in bed], ["flank", "CYP2D6", "flank"])

for c, why in (("check-no-gene", "an input without a CYP2D6 row"),
               ("check-other-region", "regions other than the ones cyp2d6-bed prints"),
               ("check-zero-length", "a zero-length region"),
               ("check-duplicate-row", "a region listed twice with two depths"),
               ("check-bad-measure", "a measure file that is not measure's"),
               ("check-out-unwritable", "an output that cannot be written")):
    expect(f"{why}: exit 2 and nothing written", (exit_code(c), written(c)), (2, []))

# The synthetic reads of make_parity.py: Cyrius's region tiled, MAPQ 0 in
# REP6 and around the body, the body at MAPQ 60 except some at MAPQ 10.
rows = [line.rstrip("\n").split("\t") for line in open(os.path.join(golden("measure-cyp2d6_synthetic"), "files", "shares.tsv"))]
m = {(r, int(q)): float(d) for r, q, d in rows[1:]}
expect("measure: the MAPQ 0 reads outside the body leave the body's MAPQ >= 1 depth whole",
       m[("body", 1)], m[("body", 0)])
expect("measure: the MAPQ 10 reads count at MAPQ >= 1, not at >= 20",
       [m[(k, 20)] < m[(k, 1)] for k in ("gene", "body")], [True, True])
expect("measure: over Cyrius's region the MAPQ 0 reads lower the MAPQ >= 1 depth",
       m[("gene", 1)] < m[("gene", 0)], True)

if FAILS:
    print(f"{len(FAILS)} check(s) failed")
    sys.exit(1)
print("all checks passed")
