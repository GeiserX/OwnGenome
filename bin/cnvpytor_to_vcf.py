#!/usr/bin/env python3
# remnant: kept for couple_report; deleted in lane R4
"""cnvpytor_to_vcf.py: what is left of step 05's input builder once it moved
to `og sv add-callers` (rust/crates/og-steps/src/sv). Only the read-depth
rule for long Delly calls is kept, because bin/couple_report.py (step 39)
still imports it; `og sv add-callers` and `og sv postfilter` apply the same
rule (og-steps::sv::depth_contradicts). `og dev remnants`
(tests/parity/remnants.tsv) checks that nothing else comes back. Standard
library only.
"""

# A Delly DEL or DUP longer than this needs read depth to agree: its FORMAT
# RDCN (Delly's read-depth copy number) must be below (DEL) or above (DUP) the
# copies a normal region has. 1,000 bp is Delly's own line between small and
# large SVs in its read counts (src/delly.h, indelsize). Short-read Delly can
# mark PASS a long SV whose RDCN shows no change; a short one rests on split
# reads alone. og sv add-callers (step 05) uses the same rule.
DEPTH_MIN_SIZE = 1000


def delly_depth_contradicts(kind, chrom, pos, end, fmt, vals, sex, in_par):
    """True when a Delly DEL or DUP longer than DEPTH_MIN_SIZE has a FORMAT
    RDCN that does not show the change: not below the normal copies for a
    DEL, not above them for a DUP. The normal is 1 on a male's chrX or chrY
    outside the PARs and 2 elsewhere; with no sex, either 1 or 2 there, so
    only an RDCN at or past both is held against the call. A missing or
    negative RDCN (no reads in the flanks) is no evidence of a change. A
    record without RDCN (not Delly 2.x) is never contradicted."""
    if end - pos <= DEPTH_MIN_SIZE or "RDCN" not in fmt:
        return False
    i = fmt.index("RDCN")
    cn = vals[i] if i < len(vals) else "."
    if not cn.isdigit():
        return True
    sex_chrom = chrom in ("chrX", "chrY") and not in_par
    if sex_chrom and sex == "male":
        low = high = 1
    elif sex_chrom and sex != "female":
        low, high = 1, 2
    else:
        low = high = 2
    if kind == "DEL":
        return int(cn) >= high
    return int(cn) <= low
