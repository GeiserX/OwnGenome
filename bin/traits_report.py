# remnant: kept for collect_summary, demo and catalog tests; deleted in lane R3.
"""Step 42 metadata still read by the Python report collectors.

The executable implementation is og traits table, phase-sites and check.
Copy-aware marker readers are og cn marker-read.
"""
COLUMNS = ["chrom", "pos", "rsid", "gene", "ref", "alt", "trait", "hom_ref", "het", "hom_alt",
           "ancestry", "note", "pmids"]
OUT_COLUMNS = ["rsid", "gene", "chrom", "pos", "trait", "genotype", "status", "reading", "ancestry",
               "note", "pmids"]
RULE_OUT_COLUMNS = ["rule", "gene", "trait", "rsids", "chrom", "pos", "genotypes", "gq", "phase", "result",
                    "status", "reading", "ancestry", "note", "pmids"]

TRAITS_TABLE_RULES = ("APOE",)

class TableError(Exception):
    pass


def read_table(path):
    """The rows of the traits table as dicts. Lines starting with # are comments;
    the first other line is the header, which must be COLUMNS."""
    rows, header = [], None
    with open(path, encoding="utf-8") as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cells = line.split("\t")
            if header is None:
                if cells != COLUMNS:
                    raise TableError(f"{path}:{n}: header is {cells}, expected {COLUMNS}")
                header = cells
                continue
            if len(cells) != len(COLUMNS):
                raise TableError(f"{path}:{n}: {len(cells)} columns, expected {len(COLUMNS)}")
            row = dict(zip(COLUMNS, cells))
            row["line"] = n
            rows.append(row)
    if header is None:
        raise TableError(f"{path}: no header line")
    return rows


MARKERS = {
    # KIV-2 4925G>A: the last base of the second KIV-2 exon (its donor splice
    # site), G>A on LPA's minus strand, so C>T on the plus strand. GRCh38 holds
    # it at the end of each of the six KIV-2 units (exons 5, 7, ..., 15 of
    # ENST00000316300.10; NM_005577.4:c.733G>A at the first, p.Ala245Thr) and
    # at the near-identical last exon of KIV-1 (exon 3, c.391): seven copies,
    # the block og cn call --locus LPA_KIV2 counts as N + 2 copies (N KIV-2 repeats and
    # one KIV-1 per chromosome). dbSNP names the first KIV-2 copy only.
    "LPA_KIV2_4925": {
        "name": "LPA KIV-2 4925G>A", "chrom": "chr6", "ref": "C", "alt": "T",
        "dbsnp": ("rs1484919544", "chr6", 160640667, "C", "T"),
        "copies": (160646214, 160640667, 160635123, 160629574, 160624030, 160618484, 160612932),
        "context": "GACTCCTTACCTTGTTCGGAA",
        "min_depth": 100, "min_alt": 5, "min_fraction": 0.008, "max_absent": 2,
    },
    # SMN1 g.27134T>G is NM_000344.4:c.*3+80T>G in intron 7 (plus strand).
    # SMN2 holds the same 21 bases 875,420 bases upstream, as at exon 7.
    "SMN1_27134": {
        "name": "SMN1 g.27134T>G", "chrom": "chr5", "ref": "T", "alt": "G",
        "dbsnp": ("rs143838139", "chr5", 70952074, "T", "G"),
        "copies": (70952074, 70076654),
        "context": "TTTTTGAACATTTAAAAAGTT",
        "min_depth": 10, "min_alt": 3, "min_fraction": 0.05, "max_absent": 1,
    },
}
MARKER_STATUSES = ("present", "not seen", "uncertain", "too little depth", "not assessed")
