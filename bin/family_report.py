#!/usr/bin/env python3
# remnant: kept for couple_report, ibd_report, pedigree_svg, annotsv_postfilter, cnvpytor_to_vcf; deleted in lane R4
"""family_report.py: what is left of the family analysis of step 38 once it
moved to `og family` (rust/crates/og-steps/src/family). Only the pedigree,
relatedness and PAR rules that other Python tools still import are kept
here; `og dev remnants` (tests/parity/remnants.tsv) checks that nothing
else comes back. Standard library only.
"""
import re
import sys

COLUMNS = ("family", "sample", "sex", "father", "mother", "partner")
MISSING = {"", ".", "0", "-9"}
SEX = {"male": "male", "m": "male", "1": "male",
       "female": "female", "f": "female", "2": "female",
       "unknown": None, "": None, ".": None, "0": None, "-9": None}
NAME = re.compile(r"^[A-Za-z0-9._-]+$")

# Relatedness (twice the kinship coefficient) bins, as KING draws them: the
# geometric middle between two degrees.
BINS = ((0.707, "same person"), (0.354, "first degree"),
        (0.177, "second degree"), (0.0884, "third degree"))


def die(lines, code=2):
    for line in lines:
        print("ERROR: " + line, file=sys.stderr)
    sys.exit(code)


def missing(v):
    return v is None or v.strip() in MISSING


def read_pedigree(path):
    """Rows of the pedigree TSV as dicts, and the problems found. The first
    non-blank line is the header (a leading '#' is dropped); later lines
    starting with '#' are comments."""
    rows, errors, header = [], [], None
    with open(path) as f:
        for n, raw in enumerate(f, 1):
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            if header is None:
                header = [c.strip().lower() for c in line.lstrip("#").split("\t")]
                if "sample" not in header:
                    return [], ["%s: the first line must be a header with a 'sample' column "
                                "(columns: %s, tab separated); it is '%s'" % (path, ", ".join(COLUMNS), line)]
                bad = [c for c in header if c not in COLUMNS]
                if bad:
                    errors.append("%s: unknown column%s %s; the columns are %s"
                                  % (path, "s" if len(bad) > 1 else "", ", ".join(bad), ", ".join(COLUMNS)))
                continue
            if line.startswith("#"):
                continue
            vals = [v.strip() for v in line.split("\t")]
            if len(vals) > len(header):
                errors.append("%s line %d: %d fields, the header has %d (fields are tab separated)"
                              % (path, n, len(vals), len(header)))
                continue
            vals += [""] * (len(header) - len(vals))
            rec = dict(zip(header, vals))
            rec["line"] = n
            rows.append(rec)
    if header is None:
        errors.append("%s is empty" % path)
    return rows, errors


def build_members(rows, declared, path, source="the samplesheet"):
    """Every person the pedigree names: its rows, and the parents and
    partners it names without a row of their own (not sequenced, no
    parents known). Returns (members, errors); a member is a dict with
    sex, father, mother, partner (None when not given) and row (bool)."""
    errors, members = [], {}

    for r in rows:
        where = "%s line %d" % (path, r["line"])
        s = r.get("sample", "").strip()
        if missing(s):
            errors.append("%s: no sample name" % where)
            continue
        if bad_name(s):
            errors.append("%s: sample name '%s' may hold only letters, digits, '.', '_' and '-' "
                          "and may not start with '-'" % (where, s))
            continue
        if s in members:
            errors.append("%s: sample '%s' has a second row; each person has one" % (where, s))
            continue
        sex_raw = r.get("sex", "").strip().lower()
        if sex_raw not in SEX:
            errors.append("%s: sex '%s' of '%s' is not male, female or unknown" % (where, r.get("sex"), s))
            sex_raw = ""
        m = {"sex": SEX[sex_raw], "row": True}
        for k in ("father", "mother", "partner"):
            v = r.get(k, "").strip()
            if missing(v):
                m[k] = None
            elif bad_name(v):
                errors.append("%s: %s '%s' of '%s' is not a valid sample name" % (where, k, v, s))
                m[k] = None
            elif v == s:
                errors.append("%s: '%s' is named as its own %s" % (where, s, k))
                m[k] = None
            else:
                m[k] = v
        if m["father"] and m["father"] == m["mother"]:
            errors.append("%s: '%s' has the same person, '%s', as father and mother" % (where, s, m["father"]))
        d = declared.get(s)
        if d:
            if m["sex"] and m["sex"] != d:
                errors.append("%s: the pedigree says '%s' is %s, %s says %s"
                              % (where, s, m["sex"], source, d))
            m["sex"] = m["sex"] or d
        members[s] = m

    # People named but without a row: parents and partners not in the file.
    for s, m in list(members.items()):
        for k, sex in (("father", "male"), ("mother", "female"), ("partner", None)):
            p = m.get(k)
            if p and p not in members:
                members[p] = {"sex": sex, "father": None, "mother": None, "partner": None, "row": False}

    for s, m in members.items():
        if m["father"] and members[m["father"]]["sex"] == "female":
            errors.append("'%s' is the father of '%s' but its sex is female" % (m["father"], s))
        if m["mother"] and members[m["mother"]]["sex"] == "male":
            errors.append("'%s' is the mother of '%s' but its sex is male" % (m["mother"], s))
    fathers = {m["father"] for m in members.values() if m["father"]}
    mothers = {m["mother"] for m in members.values() if m["mother"]}
    for p in sorted(fathers & mothers):
        errors.append("'%s' is named as a father and as a mother" % p)

    for s, m in members.items():
        p = m["partner"]
        if p and members[p]["row"] and members[p]["partner"] not in (None, s):
            errors.append("'%s' names '%s' as partner, but '%s' names '%s'"
                          % (s, p, p, members[p]["partner"]))

    # No one is their own ancestor.
    for s in members:
        seen, todo = set(), [s]
        while todo:
            x = todo.pop()
            for k in ("father", "mother"):
                p = members.get(x, {}).get(k)
                if not p:
                    continue
                if p == s:
                    errors.append("'%s' is among its own ancestors: the parent columns form a loop" % s)
                    todo = []
                    break
                if p not in seen:
                    seen.add(p)
                    todo.append(p)
    return members, errors


def bad_name(v):
    return not NAME.match(v) or v.startswith("-")


def kinship_table(members):
    """phi(a, b): the kinship coefficient of two members, from the pedigree.
    Founders (no parent known) are unrelated; an unknown parent is a founder."""
    depth_memo, phi_memo = {}, {}

    def depth(x):
        if x not in depth_memo:
            m = members[x]
            ps = [p for p in (m["father"], m["mother"]) if p]
            depth_memo[x] = 1 + max(depth(p) for p in ps) if ps else 0
        return depth_memo[x]

    def phi(a, b):
        if a == b:
            m = members[a]
            if m["father"] and m["mother"]:
                return 0.5 * (1 + phi(m["father"], m["mother"]))
            return 0.5
        key = (a, b) if a < b else (b, a)
        if key not in phi_memo:
            if depth(a) < depth(b):
                a, b = b, a
            m = members[a]
            if not m["father"] and not m["mother"]:
                phi_memo[key] = 0.0
            else:
                phi_memo[key] = 0.5 * sum(phi(p, b) for p in (m["father"], m["mother"]) if p)
        return phi_memo[key]

    return phi


def declared_relationship(a, b, members, expected):
    ma, mb = members[a], members[b]
    if a in (mb["father"], mb["mother"]) or b in (ma["father"], ma["mother"]):
        return "parent-child"
    shared = {p for p in (ma["father"], ma["mother"]) if p} & {p for p in (mb["father"], mb["mother"]) if p}
    if len(shared) == 2:
        return "full siblings"
    if len(shared) == 1:
        return "half siblings"
    if ma["partner"] == b or mb["partner"] == a:
        return "partners" if expected == 0 else "partners, %s" % degree(expected)
    return degree(expected) if expected > 0 else "unrelated"


def degree(r):
    for lo, name in BINS:
        if r >= lo:
            return name
    return "unrelated"


def read_par(path):
    par = {}
    if path:
        with open(path) as f:
            for line in f:
                c = line.split("\t")
                if len(c) >= 3:
                    par.setdefault(c[0], []).append((int(c[1]), int(c[2])))
    return par


def in_par(par, chrom, pos):
    # BED: 0-based start, end exclusive; pos is 1-based
    return any(s < pos <= e for s, e in par.get(chrom, ()))
