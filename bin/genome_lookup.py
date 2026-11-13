#!/usr/bin/env python3
# remnant: kept for results_index (names_of, walk_mentions, flat, parse_query, RSID_RE); deleted in lane R5c
"""genome_lookup.py: what is left of the genotype lookup once it moved to
`og lookup` (rust/crates/og-lookup; scripts/owngenome-lookup.sh runs
`og lookup query`). Only the query reading and the summary.json names rule
that bin/results_index.py builds the search index with are kept here, so the
search and a lookup find the same mentions; `og dev remnants`
(tests/parity/remnants.tsv) checks that nothing else comes back. Standard
library only.
"""
import re

# Keys of a summary.json entry that name a gene. "locus" is how the repeat
# expansions, the CFTR/UGT1A1/SLC6A4 readings and the Locityper rows name
# theirs (HTT, CFTR, HLA-DRB1); their regions sit under other keys (region,
# chrom), so a locus value is a name.
GENE_KEYS = ("gene", "genes", "Gene", "SYMBOL", "symbol", "locus")


def norm_chrom(c):
    """chr7 / 7 -> 7; chrM / M / MT -> MT: contig names compared across files."""
    c = c[3:] if c.lower().startswith("chr") else c
    return "MT" if c.upper() in ("M", "MT") else c.upper()


def parse_query(q):
    """('rsid', 'rs123') | ('site', ('7', 117559590)) | ('gene', 'CFTR'), or raise ValueError."""
    q = q.strip()
    if re.fullmatch(r"(?i)rs\d+", q):
        return "rsid", q.lower()
    m = re.fullmatch(r"(?i)(?:chr)?([0-9]{1,2}|X|Y|M|MT)[:_-]([0-9][0-9,]*)", q)
    if m:
        pos = int(m.group(2).replace(",", ""))
        if pos < 1:
            raise ValueError(f"position must be 1 or more: {q}")
        return "site", (norm_chrom(m.group(1)), pos)
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", q):
        return "gene", q
    raise ValueError(f"not an rsID (rs123), a position (chr7:117559590) or a gene symbol (CFTR): {q!r}")


RSID_RE = re.compile(r"(?i)^rs\d+$")


def names_of(obj):
    """(genes, rsids) a summary.json dict names itself: the genes of its gene
    keys (GENE_KEYS; a string is split on , ; & / ( ) and spaces and keeps
    the words that hold an upper-case letter, so "VDBP (GC)" gives VDBP and
    GC and the prose of "SMN exon 7-8 deletion" gives SMN only; a list gives
    its strings, a dict its keys), upper-cased, and the rsIDs of its rsid and
    rsids keys, lower-cased. '.' and '' are no name. bin/results_index.py builds the
    search index with this same function, so a lookup and the search find
    the same mentions."""
    genes, rsids = set(), set()
    if not isinstance(obj, dict):
        return genes, rsids
    for k in GENE_KEYS:
        v = obj.get(k)
        if isinstance(v, str):
            names = [n for n in re.split(r"[,;&/ ()]+", v) if n != n.lower()]
        elif isinstance(v, list):
            names = [x for x in v if isinstance(x, str)]
        elif isinstance(v, dict):
            names = list(v)
        else:
            continue
        genes.update(n.upper() for n in names if n and n != ".")
    for k in ("rsid", "rsids"):
        v = obj.get(k)
        for x in ([v] if isinstance(v, str) else v if isinstance(v, list) else []):
            if isinstance(x, str) and RSID_RE.match(x):
                rsids.add(x.lower())
    return genes, rsids


def walk_mentions(obj, want_genes, want_rsid, path=()):
    """Every dict in OBJ that names one of WANT_GENES (upper-case, by a gene
    key) or WANT_RSID (by rsid or rsids): the first dict on each branch that
    names it, not the dicts inside it."""
    if isinstance(obj, dict):
        genes, rsids = names_of(obj)
        if genes & want_genes or (want_rsid and want_rsid in rsids):
            yield obj
            return
        for v in obj.values():
            yield from walk_mentions(v, want_genes, want_rsid)
    elif isinstance(obj, list):
        for v in obj:
            yield from walk_mentions(v, want_genes, want_rsid)


def flat(d, limit=300):
    parts = []
    for k, v in d.items():
        if isinstance(v, (str, int, float, bool)) and v not in ("", None):
            parts.append(f"{k}={v}")
        elif isinstance(v, list) and v and all(isinstance(x, (str, int, float)) for x in v):
            parts.append(f"{k}={','.join(map(str, v))}")
        elif isinstance(v, dict) and v and all(isinstance(x, (str, int, float)) for x in v.values()):
            parts.append(f"{k}={'; '.join(f'{a}: {b}' for a, b in v.items())}")
    s = ", ".join(parts)
    return s if len(s) <= limit else s[:limit - 3] + "..."
