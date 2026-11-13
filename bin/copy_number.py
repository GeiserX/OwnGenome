#!/usr/bin/env python3
"""copy_number.py: the catalog reader and the depth-ratio engine of
bin/cn_loci.tsv, as the catalog tests (tests/test_copy_number.py,
test_hba_copy_number.py, test_lpa_kiv2.py) import them. No step runs it.

The commands are og cn's (rust/crates/og-cn, docs/cn-loci.md): `og cn bed`
prints a locus's segments, `og cn call` turns `samtools bedcov` over them
into the step's table, `og cn intervals` reads the interval rows and `og cn
check-bam` refuses a BAM aligned with ALT contigs, and `og bloodtype
report` reads step 41's RHD row with og cn's engine. This file keeps
load_catalog, assess, cut_copies, BrokenControl and CatalogError for those
tests until they move to og too.

One method, one catalog row per locus: the depth over the locus's target
segments against the depth over single-copy control segments, as
`samtools bedcov` sums it. With T the target depth and C the control depth
(each combined as the row's target_depth and control_depth say),
    estimate = scale x T / C - offset
then guarded by a minimum control depth, flank agreement and the share of the
target depth from reads of mapping quality 0, and called as the row's `call`
column says. docs/cn-loci.md describes the columns and how to add a locus.
Standard library only.
"""
import csv
import math
import os

CATALOG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cn_loci.tsv")
COLUMNS = ["locus", "kind", "chrom", "targets", "controls", "details", "target_depth", "control_depth",
           "scale", "offset", "bedcov", "mapq0_min", "mapq0_max", "min_control_depth", "max_control_ratio",
           "call", "interval", "readings", "input", "owner", "cross_check", "pmids"]
KINDS = ("gene", "interval", "organelle")
# GRCh38 names of the primary chromosomes (the ClinGen BED spells chrX 'chrx':
# scripts/ci/build-cn-intervals.py normalises it, and a row with any other
# spelling is refused here).
CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]}
INPUTS = ("bedcov", "cnvpytor-level")
Z95 = 1.96
# The findings a psv row needs a reading for.
PSV_FINDINGS = ("none", "deletion", "conversion", "pseudogene_deletion")
BASES = "ACGT"


class CatalogError(ValueError):
    pass


class BrokenControl(ValueError):
    """A control segment with no read while the other controls are covered."""


def _num(row, col, default=None):
    v = row[col]
    if v == ".":
        return default
    try:
        x = float(v)
    except ValueError:
        x = math.nan
    if not math.isfinite(x):
        raise CatalogError(f"{row['locus']}: {col} is {v!r}, not a number or '.'")
    return x


def _segments(row, col):
    """[(name, start, end, label)] of a targets, controls or details cell."""
    if row[col] == ".":
        return []
    out = []
    for part in row[col].split(";"):
        bits = part.split(":", 2)
        try:
            name, (s, e) = bits[0], bits[1].split("-")
            start, end = int(s), int(e)
        except (IndexError, ValueError):
            raise CatalogError(f"{row['locus']}: {col} segment {part!r} is not NAME:start-end[:label]")
        if not name or end <= start or start < 0:
            raise CatalogError(f"{row['locus']}: {col} segment {part!r} has no name or an empty span")
        out.append((name, start, end, bits[2] if len(bits) > 2 else ""))
    return out


def _readings(row):
    out = {}
    if row["readings"] == ".":
        return out
    for part in row["readings"].split("|"):
        k, sep, text = part.partition("=")
        if not sep or not k or not text:
            raise CatalogError(f"{row['locus']}: reading {part!r} is not KEY=text")
        out[k] = text
    return out


def _call(row):
    mode, sep, args = row["call"].partition(":")
    try:
        vals = [float(x) for x in args.split(",")]
    except ValueError:
        vals = []
    if not all(math.isfinite(v) for v in vals):
        vals = []
    if not sep or mode not in ("round", "cut", "class", "psv") or not vals \
            or (mode in ("round", "psv") and len(vals) != 1) or (mode == "class" and len(vals) != 2) \
            or vals != sorted(vals):
        raise CatalogError(f"{row['locus']}: call {row['call']!r} is not round:D, cut:A,B..., class:LOW,HIGH or psv:D")
    return mode, vals


def _psv_sites(row, details):
    """[(name, gene 0-based position, pseudogene 0-based position, gene
    bases, pseudogene bases, text)] of a psv row: every details segment is one
    base, labelled 'PSEUDO_POS GENE>PSEUDO text' (bases among ACGT; In2G's
    gene base is A or C)."""
    out = []
    for name, s, e, label in details:
        bits = label.split(" ", 2)
        try:
            pos = int(bits[0])
            gene, _, pseudo = bits[1].partition(">")
        except (IndexError, ValueError):
            raise CatalogError(f"{row['locus']}: PSV {name} label {label!r} is not 'PSEUDO_POS GENE>PSEUDO text'")
        if e - s != 1 or pos < 0 or not gene or not pseudo or set(gene + pseudo) - set(BASES) \
                or set(gene) & set(pseudo) or len(bits) < 3:
            raise CatalogError(f"{row['locus']}: PSV {name} must be one base with distinct gene and pseudogene "
                               f"bases and a text (got {s}-{e}, {label!r})")
        out.append((name, s, pos, set(gene), set(pseudo), bits[2]))
    return out


def parse_row(row):
    """The catalog row checked, with its parsed fields under '_' keys."""
    loc = row["locus"]
    if row["kind"] not in KINDS:
        raise CatalogError(f"{loc}: kind {row['kind']!r} (one of {', '.join(KINDS)})")
    if row["input"] not in INPUTS:
        raise CatalogError(f"{loc}: input {row['input']!r} (one of {', '.join(INPUTS)})")
    if row["chrom"] not in CHROMS:
        raise CatalogError(f"{loc}: chrom {row['chrom']!r} (GRCh38 names, chr1 to chrM)")
    p = dict(row)
    p["_targets"], p["_controls"], p["_details"] = (_segments(row, c) for c in ("targets", "controls", "details"))
    p["_call"] = _call(row)
    p["_readings"] = _readings(row)
    p["_psv"] = _psv_sites(row, p["_details"]) if p["_call"][0] == "psv" else []
    if p["_call"][0] == "psv" and (not p["_psv"] or row["target_depth"] != "sum"
                                   or not set(PSV_FINDINGS) <= set(p["_readings"])):
        raise CatalogError(f"{loc}: call psv needs PSV sites in details, target_depth sum (gene and pseudogene) "
                           f"and a reading for each of {', '.join(PSV_FINDINGS)}")
    if not p["_targets"]:
        raise CatalogError(f"{loc}: no target segment")
    if row["input"] == "bedcov":
        if not p["_controls"]:
            raise CatalogError(f"{loc}: a bedcov row needs control segments")
        if row["target_depth"] not in ("sum", "pooled") or row["control_depth"] not in ("flanks", "pooled"):
            raise CatalogError(f"{loc}: target_depth sum|pooled and control_depth flanks|pooled")
        p["_scale"], p["_offset"] = _num(row, "scale"), _num(row, "offset", 0.0)
        if not p["_scale"] or p["_scale"] <= 0:
            raise CatalogError(f"{loc}: scale must be a positive number")
        p["_min_control"] = _num(row, "min_control_depth", 0.0)
        p["_max_ratio"] = _num(row, "max_control_ratio")
        p["_mapq0_min"], p["_mapq0_max"] = _num(row, "mapq0_min"), _num(row, "mapq0_max")
        p["_passes"] = row["bedcov"].split(",")
        if (p["_mapq0_min"] is not None or p["_mapq0_max"] is not None) and len(p["_passes"]) != 2:
            raise CatalogError(f"{loc}: a MAPQ 0 guard needs two bedcov passes (all reads, MAPQ >= 1)")
        if p["_call"][0] == "class" and row["interval"] != "cv":
            raise CatalogError(f"{loc}: call class needs interval cv")
        if row["interval"] not in ("cv", "."):
            raise CatalogError(f"{loc}: interval {row['interval']!r} (cv or '.')")
    names = [n for n, *_ in p["_targets"] + p["_controls"] + p["_details"]]
    if len(names) != len(set(names)):
        raise CatalogError(f"{loc}: a segment name is used twice")
    return p


def load_catalog(path=None):
    """{locus: parsed row} of the catalog, in file order."""
    path = path or CATALOG
    with open(path) as f:
        lines = [ln for ln in f if ln.strip() and not ln.startswith("#")]
    rows = list(csv.DictReader(lines, delimiter="\t", quoting=csv.QUOTE_NONE))
    if not rows or list(rows[0].keys()) != COLUMNS:
        raise CatalogError(f"{path}: the header must be {' '.join(COLUMNS)}")
    out, seen = {}, set()
    for r in rows:
        if None in r or None in r.values():
            raise CatalogError(f"{path}: row {r.get('locus')!r} has the wrong number of columns")
        if r["locus"] in out:
            raise CatalogError(f"{path}: locus {r['locus']} is listed twice")
        p = parse_row(r)
        for n, *_ in p["_targets"] + p["_controls"] + p["_details"]:
            if n in seen:
                raise CatalogError(f"{path}: segment name {n} is used by two loci")
            seen.add(n)
        out[r["locus"]] = p
    return out


def regions(row):
    """[(chrom, start, end, name)] of a row's segments, sorted by position."""
    segs = row["_controls"] + row["_targets"] + row["_details"]
    return sorted(((row["chrom"], s, e, n) for n, s, e, _ in segs), key=lambda r: (r[1], r[2]))


def _combine(mode, segs, m):
    if mode == "sum":
        return sum(m[n] for n, *_ in segs)
    if mode == "flanks":
        return sum(m[n] for n, *_ in segs) / len(segs)
    return sum(m[n] * (e - s) for n, s, e, _ in segs) / sum(e - s for _, s, e, _ in segs)


def _cv(values):
    """Coefficient of variation (sample standard deviation over the mean)."""
    mu = sum(values) / len(values)
    if mu <= 0 or len(values) < 2:
        return 0.0
    return math.sqrt(sum((v - mu) ** 2 for v in values) / (len(values) - 1)) / mu


def cut_copies(value, cuts):
    """The number of cut points a value reaches (cut:0.5,1.5 gives 0, 1 or 2)."""
    return sum(value >= c for c in cuts)


def reading(row, key):
    r = row["_readings"]
    if str(key) in r:
        return r[str(key)]
    if isinstance(key, int) and "more" in r and key > max(int(k) for k in r if k.isdigit()):
        return r["more"].format(n=key)
    return ""


def ploidy(row, sex):
    """The copies of the row's chromosome a sample has: 2, or 1 on a male's
    chrX or chrY; None where it cannot be said (a chrY row on a female or on
    a sample of no known sex, a chrX row without a sex)."""
    if row["chrom"] not in ("chrX", "chrY"):
        return 2
    if sex == "male":
        return 1
    if sex == "female" and row["chrom"] == "chrX":
        return 2
    return None


def assess(row, m, q1=None, sex=None):
    """The locus's block for mean depths m (every read) and q1 (reads of
    MAPQ >= 1, or None for a row without a MAPQ 0 guard; mapq0_share is then
    None).

    status is 'called' or 'not assessed'; reason names the guard that
    stopped it: no_data, low_control_depth, controls_disagree,
    low_mapq0_share, high_mapq0_share, below_floor, between_counts,
    sex_unknown; status 'not applicable' (reason no_y) for a chrY row on a
    female. sex matters only for a row on chrX or chrY: there the controls
    sit on the same chromosome, so the scale follows the sample's count
    there (ploidy). expected is the count with no change."""
    guarded = q1 is not None
    if not guarded and (row["_mapq0_min"] is not None or row["_mapq0_max"] is not None):
        raise ValueError(f"{row['locus']} has a MAPQ 0 guard: the MAPQ >= 1 depths are needed")
    q1 = q1 if guarded else {}
    tg, ct = row["_targets"], row["_controls"]
    b = {"locus": row["locus"], "status": "not assessed", "reason": None, "estimate": None, "total": None,
         "copies": None, "ci95": None, "class": None, "target_depth": None, "control_depth": None,
         "depths": dict(m), "mapq0_share": None, "details": {}, "reading": "", "expected": None}
    # expected: the copies with no change, the ploidy times the reference
    # units (each target segment is one with target_depth sum: AZFc's 4 DAZ)
    units = len(row["_targets"]) if row.get("target_depth") == "sum" else 1
    pl = ploidy(row, sex)
    b["expected"] = None if pl is None else pl * units
    if row["chrom"] == "chrY" and sex == "female":
        b.update(status="not applicable", reason="no_y")
        return b
    if pl is None:
        b["reason"] = "sex_unknown"
        return b
    if not m or any(n not in m for n, *_ in tg + ct):
        b["reason"] = "no_data"
        return b
    for n, *_ in ct:
        rest = [c for c in ct if c[0] != n]
        if m[n] == 0 and rest and _combine(row["control_depth"], rest, m) >= row["_min_control"] > 0:
            raise BrokenControl(f"BROKEN_CONTROL: control segment {n} of {row['locus']} has no read while the "
                                f"other controls have {_combine(row['control_depth'], rest, m):.1f}x; check its "
                                "coordinates in bin/cn_loci.tsv (or the sample lacks that stretch)")
    T, C = _combine(row["target_depth"], tg, m), _combine(row["control_depth"], ct, m)
    tq = _combine(row["target_depth"], tg, {n: q1.get(n, 0.0) for n, *_ in tg})
    # The MAPQ 0 share is known only when the MAPQ >= 1 pass was given.
    share = (1 - tq / T if T > 0 else 0.0) if guarded else None
    b.update(target_depth=T, control_depth=C, mapq0_share=share)
    low = min(m[n] for n, *_ in ct) if row["control_depth"] == "flanks" else C
    if low < row["_min_control"] or C <= 0:
        b["reason"] = "low_control_depth"
        return b
    if row["control_depth"] == "flanks" and row["_max_ratio"] is not None and len(ct) > 1:
        vals = [m[n] for n, *_ in ct]
        if min(vals) <= 0 or max(vals) / min(vals) > row["_max_ratio"]:
            b["reason"] = "controls_disagree"
            return b
    total = row["_scale"] * pl / 2 * T / C
    est = total - row["_offset"]
    b.update(total=total, estimate=est,
             details={n: row["_scale"] * pl / 2 * m[n] / C for n, *_ in row["_details"] if n in m})
    # With no copy left there is no read to judge the aligner by (a row with
    # an offset has a floor of its own below, after this guard).
    no_copy = row["_offset"] == 0 and round(est) == 0
    if row["_mapq0_min"] is not None and not no_copy and share < row["_mapq0_min"]:
        b["reason"] = "low_mapq0_share"
        return b
    if row["_mapq0_max"] is not None and share > row["_mapq0_max"]:
        b["reason"] = "high_mapq0_share"
        return b
    if row["_offset"] > 0 and est <= 0:
        b["reason"] = "below_floor"
        return b
    mode, vals = row["_call"]
    if row["interval"] == "cv":
        rel = math.sqrt(_cv([m[n] for n, *_ in tg]) ** 2 / len(tg) + _cv([m[n] for n, *_ in ct]) ** 2 / len(ct))
        b["ci95"] = (max(0.0, total * math.exp(-Z95 * rel) - row["_offset"]),
                     total * math.exp(Z95 * rel) - row["_offset"])
    if mode == "psv":
        # the module total (gene plus pseudogene copies); og cn reads the sites
        n = round(est)
        if abs(est - n) > vals[0]:
            b["reason"] = "between_counts"
            return b
        b["module_total"] = n
    elif mode == "round":
        n = round(est)
        if abs(est - n) > vals[0]:
            b["reason"] = "between_counts"
            return b
        b["copies"] = n
        b["reading"] = reading(row, n)
    elif mode == "cut":
        b["copies"] = cut_copies(est, vals)
        b["reading"] = reading(row, b["copies"])
    else:
        lo, hi = b["ci95"]
        b["class"] = "low" if hi < vals[0] else "high" if lo > vals[1] else "typical"
        b["reading"] = reading(row, b["class"])
    b["status"] = "called"
    return b
