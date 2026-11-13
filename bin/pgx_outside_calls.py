#!/usr/bin/env python3
# remnant: kept for collect_summary (and its test); deleted in lane R3
"""pgx_outside_calls.py: what is left of step 36's outside calls once they
moved to `og pgx outside-calls`, `og pgx hla-concordance` and `og pgx
pypgx-summary` (rust/crates/og-steps/src/pgx_rules.rs and pypgx_summary.rs).
Only the rules the summary still imports are kept: what counts as a call,
the probable CYP2D6 call of pypgx alone, and the names of step 32's no-call
statuses. `og dev remnants` (tests/parity/remnants.tsv) checks that nothing
else comes back. Standard library only.
"""

HLA_GENES = ("HLA-A", "HLA-B")
NO_CALL = {"", "none", "n/a", "na", "failed", "indeterminate", "not called", "unknown", "."}
PROBABLE = "probable (one caller: pypgx)"
# Step 32's summary: the status of a gene the reference cannot hold, and the
# quality flag of a call resting on a defining variant below the GQ floor.
NOT_CALLABLE = "not callable on this reference"
LOW_GQ = "low-confidence: check the reads"


def is_call(diplotype):
    d = (diplotype or "").strip()
    return d.lower() not in NO_CALL and "/" in d and ";" not in d


def same_diplotype(a, b):
    """Equal up to the order of the two alleles."""
    def key(d):
        return sorted(x.strip() for x in d.split("/"))
    return key(a) == key(b)


def probable_cyp2d6(row, cyrius, depth_status):
    """{'diplotype', 'phenotype', 'label'} for pypgx's CYP2D6 call when it is
    the only usable call, else None. row: pypgx's CYP2D6 summary row (a
    dict); cyrius: (genotype, filter), or None when Cyrius did not run;
    depth_status: the depth check's status, or None without its file (step 32
    turns the row into Indeterminate when the check fails, so a call in the
    row passed it).

    Not probable: no call, a usable Cyrius call (then the two callers decide),
    a depth check that did not pass, a defining variant that is missing from
    the VCF or not PASS (a summary written before these columns has none), a
    star allele pypgx's table lacks, or a defining variant with a GQ below the
    floor or no GQ."""
    if not row or not is_call(row.get("Diplotype")):
        return None
    if cyrius is not None and is_call(cyrius[0]) and cyrius[1] == "PASS":
        return None
    if depth_status is not None and depth_status != "ok":
        return None
    if (row.get("Defining_filter") or "") not in ("PASS", "none") or row.get("Quality") not in ("ok", "N/A"):
        return None
    return {"diplotype": row["Diplotype"].strip(), "phenotype": (row.get("Phenotype") or "").strip(),
            "label": PROBABLE}
