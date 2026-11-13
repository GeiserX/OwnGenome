"""Writes the parity goldens of `og cdna check`, `og cdna scan`,
`og comphet phase` and `og bam contigs`.

Run with pysam, from the repository root of a commit that still has
bin/clip_check.py (with its read check), bin/cdna_scan.py and
bin/comphet_phase.py, for example the parent of the pull request that
deleted them:

    python3 -I rust/crates/og-bam/tests/make_reads_parity.py

It writes tests/data/reads/: the input BAMs, tables and VCFs it makes from
made-up reads (the fixtures of the unit tests those tools had, and
tests/fixtures/comphet/make_fixture.py), cases.tsv (one row per case:
name, tool, arguments) and, for each case, the Python's exit code, its
stdout and every file it wrote. `{d}` in an argument is tests/data and
`{o}` the case's output directory. tests/parity.rs runs og on the same rows.
The `og/` folders are written by hand, and this script leaves them as they
are.
"""

import contextlib
import gzip
import io
import os
import shutil
import sys

import pysam

REPO = os.getcwd()
DATA = os.path.join(REPO, "rust/crates/og-bam/tests/data")
OUT = os.path.join(DATA, "reads")
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.join(REPO, "tests/fixtures/comphet"))
import cdna_scan as cs  # noqa: E402
import clip_check as cc  # noqa: E402
import comphet_phase as cp  # noqa: E402
import make_fixture  # noqa: E402

CASES = []


def case(name, tool, *args):
    CASES.append((name, tool, list(args)))


def put(name, text, gz=False):
    path = os.path.join(OUT, "in", name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        # mtime=0, so a gzip input is the same bytes on every run.
        fh.write(gzip.compress(text.encode(), mtime=0) if gz else text.encode())
    return "{d}/reads/in/" + name


def bam(name, contigs, reads):
    """A sorted, indexed BAM of reads (contig, name, flag, start, cigar, seq, sa)."""
    path = os.path.join(OUT, "in", name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": c, "LN": n} for c, n in contigs]}
    tid = {c: i for i, (c, _) in enumerate(contigs)}
    unsorted = path + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for contig, qname, flag, start, cigar, seq, sa in reads:
            a = pysam.AlignedSegment()
            a.query_name, a.flag = qname, flag
            a.reference_id, a.reference_start, a.mapping_quality = tid[contig], start, 60
            a.cigarstring = cigar
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            a.next_reference_id, a.next_reference_start, a.template_length = -1, -1, 0
            if sa:
                a.set_tag("SA", sa)
            out.write(a)
    pysam.sort("-o", path, unsorted)
    os.remove(unsorted)
    pysam.index(path)
    return "{d}/reads/in/" + name


def sam_bam(name, sam_text):
    """A sorted, indexed BAM of SAM text with its header."""
    path = os.path.join(OUT, "in", name)
    sam = path + ".sam"
    with open(sam, "w") as fh:
        fh.write(sam_text)
    unsorted = path + ".unsorted.bam"
    with pysam.AlignmentFile(sam, "r") as src, pysam.AlignmentFile(unsorted, "wb", template=src) as out:
        for r in src:
            out.write(r)
    pysam.sort("-o", path, unsorted)
    os.remove(unsorted)
    os.remove(sam)
    pysam.index(path)
    return "{d}/reads/in/" + name


# --- the read check: test_clip_check.py's fixture ------------------------------------

def genome(n=10000):
    """A made-up contig: the same bases every run."""
    s, out = 7, []
    for _ in range(n):
        s = (s * 75 + 74) % 65537
        out.append("ACGT"[s % 4])
    return "".join(out)


G = genome()
EXON_A, EXON_B = (1000, 1200), (3000, 3200)
POS0 = 1150
REF = G[POS0]
ALT = "ACGT"[("ACGT".index(REF) + 1) % 4]


def with_alt(seq, start):
    i = POS0 - start
    return seq[:i] + ALT + seq[i + 1:]


def clip_reads(kind):
    reads = []
    for i in range(6):
        start = 1112 + 2 * i
        if kind == "both":
            aligned = EXON_A[1] - start
            reads.append(("chrT", f"ref{i}", 0, start, f"{aligned}M{100 - aligned}S",
                          G[start:EXON_A[1]] + G[EXON_B[0]:EXON_B[0] + 100 - aligned],
                          f"chrT,{EXON_B[0] + 1},+,{aligned}S{100 - aligned}M,60,0;"))
        else:
            reads.append(("chrT", f"ref{i}", 0, start, "100M", G[start:start + 100], ""))
    for i in range(8):
        start = 1110 + 2 * i
        if kind == "genomic":
            reads.append(("chrT", f"alt{i}", 0, start, "100M", with_alt(G[start:start + 100], start), ""))
        else:
            aligned = EXON_A[1] - start
            seq = with_alt(G[start:EXON_A[1]], start) + G[EXON_B[0]:EXON_B[0] + 100 - aligned]
            reads.append(("chrT", f"alt{i}", 0, start, f"{aligned}M{100 - aligned}S", seq,
                          f"chrT,{EXON_B[0] + 1},+,{aligned}S{100 - aligned}M,60,0;"))
    if kind != "genomic":
        reads.append(("chrT", "suppl", 2048, 1140, "60M40H", with_alt(G[1140:EXON_A[1]], 1140),
                      f"chrT,{EXON_B[0] + 1},+,60S40M,60,0;"))
        # Not counted: a duplicate, a secondary, a QC failure and a repeat of alt0.
        seq = with_alt(G[1110:1210], 1110)
        reads += [("chrT", "dup", 1024, 1110, "100M", seq, ""), ("chrT", "sec", 256, 1110, "100M", seq, ""),
                  ("chrT", "qc", 512, 1110, "100M", seq, ""), ("chrT", "alt0", 0, 1110, "100M", seq, "")]
    return reads


def indel_reads():
    """Reads with a 3 bp deletion and a TT insertion after POS0, and reads without."""
    dref = G[POS0:POS0 + 4]
    reads = []
    for i in range(4):
        reads.append(("chrT", f"del{i}", 0, 1100 + i, f"{51 - i}M3D{46 + i}M",
                      G[1100 + i:POS0 + 1] + G[POS0 + 4:1200 + i], ""))
        reads.append(("chrT", f"ins{i}", 0, 1100 + i, f"{51 - i}M2I{49 + i}M",
                      G[1100 + i:POS0 + 1] + "TT" + G[POS0 + 1:1200 + i], ""))
        reads.append(("chrT", f"plain{i}", 0, 1100 + i, "100M", G[1100 + i:1200 + i], ""))
    # Neither: a read that stops inside the deletion's gap, one with only
    # 2 of the 3 bases deleted, one with GA inserted.
    reads.append(("chrT", "stop", 0, 1100, "52M", G[1100:1152], ""))
    reads.append(("chrT", "part", 0, 1100, "51M2D47M", G[1100:POS0 + 1] + G[POS0 + 3:1200], ""))
    reads.append(("chrT", "other", 0, 1100, "51M2I49M", G[1100:POS0 + 1] + "GA" + G[POS0 + 1:1200], ""))
    return reads, dref


def check_cases():
    exons = put("exons.bed", f"chrT\t{EXON_A[0]}\t{EXON_A[1]}\nchrT\t{EXON_B[0]}\t{EXON_B[1]}\n"
                "# a comment line\nchrT\tstart\tend\n")
    head = "chrom\tpos\tref\talt\tgenotype\tclinvar_id\tgeneinfo\tclnsig\tclnrevstat\n"
    hits = put("hits.tsv", head
               + f"chrT\t{POS0 + 1}\t{REF}\t{ALT}\t0/1\t1\tGENET:1\tPathogenic\tx\n"
               + "chrQ\t500\tA\tG\t0/1\t2\tGENEQ:2\tPathogenic\tx\n"
               + f"chrT\t{POS0 + 1}\t{REF}AC\tGT\t0/1\t3\tGENET:1\tPathogenic\tx\n"
               + f"T\t{POS0 + 1}\t{REF.lower()}\t{ALT.lower()}\t0/1\t4\tGENET:1\tPathogenic\tx\n"
               + f"chrT\t{POS0 + 1}\t{REF}\t<DEL>\t0/1\t5\tGENET:1\tPathogenic\tx\n")
    catalog = put("construct.tsv", "# a made-up construct site\nchrom\tpos\tref\talt\tgene\tname\n"
                  f"T\t{POS0 + 1}\t{REF}\t{ALT}\tGENET\tX1Y\nchrQ\t500\tA\tG\tGENEQ\tQ2R\n")
    for kind in ("spliced", "genomic", "both"):
        b = bam(f"clip_{kind}.bam", [("chrT", len(G))], clip_reads(kind))
        case(f"check-{kind}", "cdna-check", "--hits", hits, "--bam", b, "--exons", exons, "--out", "{o}/out.tsv")
        case(f"check-{kind}-construct", "cdna-check", "--hits", hits, "--bam", b, "--exons", exons,
             "--construct", catalog, "--out", "{o}/out.tsv")
    spliced = "{d}/reads/in/clip_spliced.bam"
    case("check-skip", "cdna-check", "--hits", hits, "--skip", "no BAM", "--out", "{o}/out.tsv")
    case("check-skip-construct", "cdna-check", "--hits", hits, "--skip", "no BAM", "--construct", catalog,
         "--out", "{o}/out.tsv")
    case("check-shipped-construct", "cdna-check", "--hits", hits, "--skip", "no BAM", "--construct",
         "{r}/assets/construct_variants.tsv", "--out", "{o}/out.tsv")

    reads, dref = indel_reads()
    b = bam("indels.bam", [("chrT", len(G))], reads)
    ih = put("indel_hits.tsv", "pos\tchrom\tref\talt\n"
             f"{POS0 + 1}\tchrT\t{dref}\t{dref[0]}\n{POS0 + 1}\tchrT\t{REF}\t{REF}TT\n"
             f"{POS0 + 1}\tchrT\t{G[POS0:POS0 + 2]}\t{ALT}{G[POS0 + 1]}\n")
    case("check-indels", "cdna-check", "--hits", ih, "--bam", b, "--exons", exons, "--out", "{o}/out.tsv")

    # CPSR: P/LP rows by either class column, VAR_ID or GENOMIC_CHANGE, without chr.
    cls = put("cpsr.classification.tsv.gz",
              "VAR_ID\tGENOMIC_CHANGE\tCLASSIFICATION\tCPSR_CLASSIFICATION\tSYMBOL\n"
              f"T_{POS0 + 1}_{REF}_{ALT}\t.\tPathogenic\tVUS\tGENET\n"
              f".\tT:g.{POS0 + 1}{REF}>{ALT}\tVUS\tLikely_Pathogenic\tGENET\n"
              f"T_{POS0 + 1}_{REF}_{ALT}\t.\tVUS\tVUS\tGENET\n"
              f"T_{POS0 + 1}_{REF}_{ALT}\t.\tLP\tNA\tGENET\n"
              ".\tbad\tPathogenic\tPathogenic\tGENEX\n", gz=True)
    case("check-cpsr", "cdna-check", "--cpsr", cls, "--bam", spliced, "--exons", exons, "--construct", catalog,
         "--out", "{o}/out.tsv")
    case("check-cpsr-skip", "cdna-check", "--cpsr", cls, "--skip", "no BAM", "--out", "{o}/out.tsv")

    # Refused: a missing exon file, a hits table without a column, a catalog
    # without one, no reads and no reason, both sources.
    case("check-no-exons", "cdna-check", "--hits", hits, "--bam", spliced, "--exons", "{d}/reads/in/none.bed",
         "--out", "{o}/out.tsv")
    bad = put("bad_hits.tsv", "chrom\tpos\tref\nchrT\t1\tA\n")
    case("check-hits-without-alt", "cdna-check", "--hits", bad, "--skip", "x", "--out", "{o}/out.tsv")
    nocol = put("bad_construct.tsv", "chrom\tpos\tref\talt\tgene\nchrT\t1\tA\tC\tG\n")
    case("check-construct-without-name", "cdna-check", "--hits", hits, "--skip", "x", "--construct", nocol,
         "--out", "{o}/out.tsv")
    badpos = put("badpos_construct.tsv", "chrom\tpos\tref\talt\tgene\tname\nchrT\tx\tA\tC\tG\tN\n")
    case("check-construct-bad-pos", "cdna-check", "--hits", hits, "--skip", "x", "--construct", badpos,
         "--out", "{o}/out.tsv")
    case("check-no-reads", "cdna-check", "--hits", hits, "--bam", spliced, "--out", "{o}/out.tsv")
    case("check-two-sources", "cdna-check", "--hits", hits, "--cpsr", cls, "--skip", "x", "--out", "{o}/out.tsv")


# --- the scan: test_cdna_scan.py's fixture -------------------------------------------

E1, E2, E3 = (1000, 1200), (3000, 3150), (5000, 5200)
LENGTH = 20000


def across(contig, name, left, right, aligned):
    clip = 100 - aligned
    return (contig, name, 0, left[1] - aligned, f"{aligned}M{clip}S", "A" * 100,
            f"{contig},{right[0] + 1},+,{aligned}S{clip}M,60,0;")


def scan_reads(kind, contig="chrT"):
    if kind == "genomic":
        return [(contig, f"g{i}", 0, e[1] - 50 + 7 * i, "100M", "A" * 100, "") for e in (E1, E2, E3) for i in range(4)] \
            + [(contig, f"s{i}", 0, e[0] - 50 + 7 * i, "100M", "A" * 100, "") for e in (E2, E3) for i in range(4)]
    if kind == "one":
        return [across(contig, f"a{i}", E1, E2, 60 + 5 * i) for i in range(4)]
    if kind == "thin":
        return [across(contig, "a0", E1, E2, 60), across(contig, "b0", E2, E3, 60), across(contig, "c0", E1, E3, 60)]
    reads = ([across(contig, f"a{i}", E1, E2, 60 + 5 * i) for i in range(3)]
             + [(contig, "n0", 0, E1[1] - 50, f"50M{E2[0] - E1[1]}N50M", "A" * 100, "")]
             + [across(contig, f"b{i}", E2, E3, 60 + 5 * i) for i in range(3)]
             + [across(contig, f"c{i}", E1, E3, 60 + 5 * i) for i in range(2)])
    if kind == "busy":
        a = across(contig, "a0", E1, E2, 60)
        reads += [
            # A duplicate, a secondary, the supplementary record of a0 (counted once).
            (contig, "d", 1024) + a[3:], (contig, "s", 256) + a[3:],
            (contig, "a0", 2048, E2[0], "60S40M", "A" * 100, f"{contig},{a[3] + 1},+,60M40S,60,0;"),
            # Across two junctions with N gaps: counts on each, once in the locus.
            (contig, "n1", 0, E1[1] - 50, f"50M{E2[0] - E1[1]}N150M{E3[0] - E2[1]}N50M", "A" * 250, ""),
            # Pieces that do not pair: the other strand, out of read order, too far apart.
            (contig, "x1", 0, E1[1] - 60, "60M40S", "A" * 100, f"{contig},{E2[0] + 1},-,60S40M,60,0;"),
            (contig, "x2", 0, E1[1] - 60, "60M40S", "A" * 100, f"{contig},{E2[0] + 1},+,10S40M50S,60,0;"),
            (contig, "x3", 0, E1[1] - 60, "60M40S", "A" * 100, f"{contig},{E2[0] + 1},+,64S36M,60,0;"),
            # A reverse read whose SA piece is on its own strand pairs.
            (contig, "r1", 16, E1[1] - 60, "60M40S", "A" * 100, f"{contig},{E2[0] + 1},-,60S40M,60,0;"),
            (contig, "r2", 16, E1[1] - 60, "60M40S", "A" * 100, f"{contig},{E2[0] + 1},-,57S43M,60,0;"),
            # Too short an intron.
            (contig, "i1", 0, E1[1] - 50, "50M20N50M", "A" * 100, ""),
        ]
    return reads


def vector_reads():
    """A 150 bp read: bases 0-50 on a vector contig (the primary), 50-100 at
    E1's end and 100-150 at E2's start (both supplementary), twice."""
    out = []
    for q in ("v1", "v2"):
        out += [
            ("chrV", q, 0, 0, "50M100S", "A" * 150,
             f"chrT,{E1[1] - 49},+,50S50M50S,60,0;chrT,{E2[0] + 1},+,100S50M,60,0;"),
            ("chrT", q, 2048, E1[1] - 50, "50S50M50S", "A" * 150,
             f"chrV,1,+,50M100S,60,0;chrT,{E2[0] + 1},+,100S50M,60,0;"),
            ("chrT", q, 2048, E2[0], "100S50M", "A" * 150,
             f"chrV,1,+,50M100S,60,0;chrT,{E1[1] - 49},+,50S50M50S,60,0;"),
        ]
    return out


def scan_cases():
    exons = put("scan_exons.bed", "".join(f"chrT\t{a}\t{b}\n" for a, b in (E1, E2, E3)))
    for kind in ("spliced", "genomic", "one", "thin", "busy"):
        b = bam(f"scan_{kind}.bam", [("chrT", LENGTH)], scan_reads(kind))
        case(f"scan-{kind}", "cdna-scan", "--bam", b, "--exons", exons, "--out", "{o}/scan.tsv")
    b = bam("scan_vector.bam", [("chrV", 500), ("chrT", LENGTH)], vector_reads()
            + [across("chrT", f"a{i}", E1, E3, 60 + 5 * i) for i in range(2)])
    case("scan-vector", "cdna-scan", "--bam", b, "--exons", exons, "--out", "{o}/scan.tsv")
    # Two contigs in header order chr2, chr10, each with the spliced gene:
    # the Python lists chr10 first (contig names sorted as text); og lists
    # them in the BAM's order (og/scan-genome-order).
    multi = put("scan_exons_multi.bed", "".join(f"{c}\t{a}\t{b}\n" for c in ("chr10", "chr2", "chrU")
                                                for a, b in (E1, E2, E3)))
    b = bam("scan_two_contigs.bam", [("chr2", LENGTH), ("chr10", LENGTH), ("chrZ", LENGTH)],
            scan_reads("spliced", "chr2") + scan_reads("spliced", "chr10") + scan_reads("spliced", "chrZ"))
    case("scan-genome-order", "cdna-scan", "--bam", b, "--exons", multi, "--out", "{o}/scan.tsv")
    # Exons 10 bp apart are two regions: a read across both (w) is fetched
    # in each, and still counts once on its junction.
    close = put("scan_exons_close.bed", "chrT\t1000\t1200\nchrT\t1210\t1300\nchrT\t3000\t3150\n")
    b = bam("scan_close.bam", [("chrT", LENGTH)],
            [("chrT", f"w{i}", 0, 1150, f"50M10D90M{3000 - 1300}N50M", "A" * 190, "") for i in range(3)]
            + [("chrT", f"z{i}", 0, 1150, f"50M{3000 - 1200}N50M", "A" * 100, "") for i in range(2)])
    case("scan-close-regions", "cdna-scan", "--bam", b, "--exons", close, "--out", "{o}/scan.tsv")
    case("scan-skip", "cdna-scan", "--skip", "no BAM", "--out", "{o}/scan.tsv")
    spliced = "{d}/reads/in/scan_spliced.bam"
    case("scan-no-exons", "cdna-scan", "--bam", spliced, "--exons", "{d}/reads/in/none.bed", "--out", "{o}/scan.tsv")
    other = put("scan_exons_other.bed", "chrQ\t1\t100\n")
    case("scan-no-shared-contig", "cdna-scan", "--bam", spliced, "--exons", other, "--out", "{o}/scan.tsv")
    empty = put("scan_exons_empty.bed", "# nothing\n")
    case("scan-empty-exons", "cdna-scan", "--bam", spliced, "--exons", empty, "--out", "{o}/scan.tsv")
    case("scan-no-reads", "cdna-scan", "--bam", spliced, "--out", "{o}/scan.tsv")


# --- the phase: tests/fixtures/comphet ------------------------------------------------

def phase_cases():
    for swap in (False, True):
        tag = "swap" if swap else "fixture"
        seq = make_fixture.reference()
        vs = make_fixture.variants(seq, swap)
        b = sam_bam(f"comphet_{tag}.bam", make_fixture.sam(seq, vs))
        vcf_text = make_fixture.comphet_vcf(vs)
        vcf = put(f"comphet_{tag}.vcf", vcf_text)
        rows = ["GENE\tCHROM\tPOS\tREF\tALT\tIMPACT\tConsequence\tGT\tMAX_AF"]
        for line in vcf_text.splitlines():
            if not line.startswith("#"):
                f = line.split("\t")
                rows.append(f"G\t{f[0]}\t{f[1]}\t{f[3]}\t{f[4]}\tHIGH\tstop_gained\t0/1\t0.0001")
        variants = put(f"comphet_{tag}_variants.tsv", "\n".join(rows + ["G\tchrP\t1\tA\tC\tLOW\tx\t0/1\t."]) + "\n")
        case(f"phase-{tag}", "comphet-phase", "--vcf", vcf, "--bam", b, "--out", "{o}/pairs.tsv",
             "--variants", variants, "--variants-out", "{o}/variants.tsv")
    vcf_text = make_fixture.comphet_vcf(make_fixture.variants(make_fixture.reference()))
    gz = put("comphet_fixture.vcf.gz", vcf_text, gz=True)
    b = "{d}/reads/in/comphet_fixture.bam"
    case("phase-gz", "comphet-phase", "--vcf", gz, "--bam", b, "--out", "{o}/pairs.tsv")
    case("phase-skip", "comphet-phase", "--vcf", gz, "--skip", "no BAM", "--out", "{o}/pairs.tsv")
    nogt = put("comphet_variants_nogt.tsv", "CHROM\tPOS\tREF\tALT\nchrP\t1000\tx\ty\n")
    case("phase-skip-no-gt", "comphet-phase", "--vcf", gz, "--skip", "no BAM", "--out", "{o}/pairs.tsv",
         "--variants", nogt, "--variants-out", "{o}/variants.tsv")
    # Pairs the reads do not phase: another contig, a contig the BAM lacks,
    # a complex allele, a pair id naming three variants, a gene with a slash.
    head = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    rec = lambda c, p, r, a, e: f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\tslivar_comphet={e}\n"  # noqa: E731
    odd = put("comphet_odd.vcf", head
              + rec("chrP", 1000, "A", "C", "S/GENEX/1/chrQ/500/G/T")
              + rec("chrQ", 500, "G", "T", "S/GENEX/1/chrP/1000/A/C")
              + rec("chrZ", 10, "A", "C", "S/GENEZ/2/chrZ/20/A/T")
              + rec("chrP", 5000, "AC", "GTT", "S/HLA/A/3/chrP/5400/A/C")
              + rec("chrP", 7000, "A", "C", "S/GENEY/4/chrP/7001/A/C,S/GENEY/4/chrP/7002/A/C")
              + "chrP\t8000\t.\tA\tC\t.\tPASS\tOTHER=1\n")
    case("phase-odd", "comphet-phase", "--vcf", odd, "--bam", b, "--out", "{o}/pairs.tsv")
    nochr = put("comphet_nochr.vcf", head
                + rec("P", 5000, *make_fixture.variants(make_fixture.reference())[2][2:4],
                      "S/GENET/9/P/5400/{}/{}".format(*make_fixture.variants(make_fixture.reference())[3][2:4])))
    case("phase-nochr", "comphet-phase", "--vcf", nochr, "--bam", b, "--out", "{o}/pairs.tsv")
    bad = put("comphet_bad.vcf", head + rec("chrP", 1000, "A", "C", "S/GENE/1/chrP"))
    case("phase-bad-entry", "comphet-phase", "--vcf", bad, "--bam", b, "--out", "{o}/pairs.tsv")
    case("phase-no-vcf", "comphet-phase", "--vcf", "{d}/reads/in/none.vcf", "--bam", b, "--out", "{o}/pairs.tsv")
    case("phase-bam-and-skip", "comphet-phase", "--vcf", gz, "--bam", b, "--skip", "x", "--out", "{o}/pairs.tsv")
    case("phase-variants-alone", "comphet-phase", "--vcf", gz, "--skip", "x", "--out", "{o}/pairs.tsv",
         "--variants", nogt)


def contigs_cases():
    case("contigs-one", "bam-contigs", "{d}/reads/in/scan_spliced.bam")
    case("contigs-three", "bam-contigs", "{d}/reads/in/scan_two_contigs.bam")
    case("contigs-missing", "bam-contigs", "{d}/reads/in/none.bam")


def bam_contigs(argv):
    """The one-liners of step 32 and the PYPGX module it replaces."""
    try:
        print(chr(10).join(pysam.AlignmentFile(argv[0]).references))
    except (OSError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    return 0


MAINS = {"cdna-check": cc.main, "cdna-scan": cs.main, "comphet-phase": cp.main, "bam-contigs": bam_contigs}


def run(name, tool, args):
    d = os.path.join(OUT, "want", name)
    shutil.rmtree(d, ignore_errors=True)
    files = os.path.join(d, "files")
    os.makedirs(files)
    real = [a.replace("{d}", DATA).replace("{r}", REPO).replace("{o}", files) for a in args]
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
        try:
            rc = MAINS[tool](real)
        except SystemExit as e:
            rc = e.code if isinstance(e.code, int) else 1
    with open(os.path.join(d, "exit"), "w") as fh:
        fh.write(f"{rc}\n")
    with open(os.path.join(d, "stdout"), "w") as fh:
        fh.write(out.getvalue())


def main():
    for sub in ("in", "want"):
        shutil.rmtree(os.path.join(OUT, sub), ignore_errors=True)
    check_cases()
    scan_cases()
    phase_cases()
    contigs_cases()
    with open(os.path.join(OUT, "cases.tsv"), "w") as fh:
        fh.write("# case<TAB>tool<TAB>arguments...: made by ../../make_reads_parity.py; "
                 "{d} tests/data, {r} the repository, {o} the output directory\n")
        for name, tool, args in CASES:
            fh.write("\t".join([name, tool] + args) + "\n")
    for name, tool, args in CASES:
        run(name, tool, args)
    print(f"{len(CASES)} cases written to {OUT}")


if __name__ == "__main__":
    main()
