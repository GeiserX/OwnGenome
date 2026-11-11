#!/usr/bin/env python3
"""chip_report.py: what the reports still read of step 49 (clonal haematopoiesis) in Python.

The step itself runs in Rust: `og chip prepare`, `og chip report` and
`og chip check-table` (rust/crates/og-chip, docs/49-chip.md). This module
keeps only the names bin/collect_summary.py and tests/demo/make_demo_sample.py
import: the two tables' columns, the gene table reader and the detection
limit, the smallest VAF a clone needs to be found with probability
DETECT_POWER, taken as at least MIN_ALT_READS reads with the change out of the
mean depth (binomial). og chip computes the same limit the same way.

Standard library only (it runs in the plain python image).
"""
import math
import sys

TABLE_COLUMNS = ["gene", "truncating", "missense", "note", "pmids"]
OUT_COLUMNS = ["gene", "chrom", "pos", "ref", "alt", "consequence", "protein_change", "exon", "vaf",
               "alt_reads", "depth", "filter", "status", "class", "note"]
DEPTH_COLUMNS = ["gene", "transcript", "chrom", "coding_bases", "mean_depth", "detection_limit_vaf", "status"]
MIN_ALT_READS = 3     # reads with the change a call needs, and the detection limit assumes
DETECT_POWER = 0.8    # the detection limit: the VAF found this often at the mean depth


def die(msg):
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


def read_table(path):
    """The rows of assets/chip_genes.tsv as dicts (og chip check-table checks the rest of its schema)."""
    rows = []
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if not ln.startswith("#") and ln.strip()]
    if not lines or lines[0].split("\t") != TABLE_COLUMNS:
        die(f"{path}: the header must be {' '.join(TABLE_COLUMNS)} (tab separated)")
    for ln in lines[1:]:
        f = ln.split("\t")
        if len(f) != len(TABLE_COLUMNS):
            die(f"{path}: '{ln[:60]}' has {len(f)} columns, not {len(TABLE_COLUMNS)}")
        rows.append(dict(zip(TABLE_COLUMNS, f)))
    return rows


def p_at_least(k, n, f):
    """P(X >= k) for X ~ Binomial(n, f)."""
    return 1.0 - sum(math.comb(n, i) * f ** i * (1 - f) ** (n - i) for i in range(k))


def detection_limit(depth, k=MIN_ALT_READS, power=DETECT_POWER):
    """The smallest VAF found with probability `power` at `depth` (k reads needed); None below k reads."""
    n = int(depth)
    if n < k:
        return None
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if p_at_least(k, n, mid) >= power:
            hi = mid
        else:
            lo = mid
    return hi
