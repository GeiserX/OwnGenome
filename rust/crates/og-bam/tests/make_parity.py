"""Writes the parity goldens of `og pgx cyp2d6-*` and `og str motif`.

Run in PYPGX_IMAGE (versions.env), which has the pysam the two Python
tools read BAMs with, from the repository root, with the commit that
still has bin/cyp2d6_depth_check.py and bin/rfc1_motif.py:

    docker run --rm --network none -u "$(id -u):$(id -g)" -e HOME=/tmp \\
      -v "$PWD:/repo" -w /repo "$PYPGX_IMAGE" \\
      python3 -I rust/crates/og-bam/tests/make_parity.py

It writes tests/data/parity/: the input BAMs, VCFs and catalogs it makes,
cases.tsv (one row per case: name, tool, arguments), and for each case a
directory with the Python's exit code, its stdout and every file it wrote.
`{d}` in an argument is tests/data, `{r}` the repository root and `{o}`
the case's output directory. tests/parity.rs runs og on the same rows.
"""

import contextlib
import gzip
import io
import json
import os
import random
import re
import shutil
import sys
import tempfile

import pysam

REPO = os.getcwd()
DATA = os.path.join(REPO, "rust/crates/og-bam/tests/data")
OUT = os.path.join(DATA, "parity")
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.join(REPO, "tests"))
import cyp2d6_depth_check as dc  # noqa: E402
import rfc1_motif as rm  # noqa: E402
import test_rfc1_motif as t  # noqa: E402

DEPTH = "{r}/tests/fixtures/pgx/depth"
CASES = []


def case(name, tool, *args):
    CASES.append((name, tool, list(args)))


def write(name, text, gz=False):
    path = os.path.join(OUT, "in", name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        # mtime=0, so a gzip input is the same bytes on every run.
        fh.write(gzip.compress(text.encode(), mtime=0) if gz else text.encode())
    return "{d}/parity/in/" + name


def sites_only(vcf):
    """The VCF without its FORMAT and sample columns: the motif count reads
    columns 1 to 8 only, and a repeat-count genotype in a tracked file trips
    the personal-data check."""
    return "".join("\t".join(line.rstrip("\n").split("\t")[:8]) + "\n" for line in vcf.splitlines(keepends=True))


def bam_from(reads, path, header):
    unsorted = path + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for a in reads:
            out.write(a)
    pysam.sort("-o", path, unsorted)
    pysam.index(path)
    os.remove(unsorted)


# --- CYP2D6 -------------------------------------------------------------------
def cyp2d6_inputs():
    src = os.path.join(DATA, "cyp2d6.bam")
    with pysam.AlignmentFile(src) as b:
        reads = list(b.fetch())
        hdr = b.header.to_dict()
    # The same reads on a contig named without chr, and on another contig.
    for name, sn in (("cyp2d6_nochr.bam", "22"), ("cyp2d6_other.bam", "chrZ")):
        h = json.loads(json.dumps(hdr))
        for sq in h["SQ"]:
            if sq["SN"] == "chr22":
                sq["SN"] = sn
        new = pysam.AlignmentHeader.from_dict(h)
        out = []
        for r in reads:
            d = r.to_dict()
            for k in ("ref_name", "next_ref_name"):
                if d[k] == "chr22":
                    d[k] = sn
            out.append(pysam.AlignedSegment.from_dict(d, new))
        bam_from(out, os.path.join(OUT, "in", name), h)
    # Reads that test each floor and flag: MAPQ 0 in REP6, MAPQ 10 in the
    # body, a duplicate, a deleted stretch, a pair whose mates overlap.
    h = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr22", "LN": 50818468}]}
    head = pysam.AlignmentHeader.from_dict(h)
    rng = random.Random(11)
    reads = []

    def read(name, pos, cigar, mapq, flag=0):
        a = pysam.AlignedSegment(head)
        a.query_name, a.flag, a.reference_id, a.reference_start = name, flag, 0, pos
        a.cigarstring, a.mapping_quality = cigar, mapq
        qlen = sum(int(n) for n, op in re.findall(r"(\d+)([MIS=X])", cigar))
        a.query_sequence = "".join(rng.choice("ACGT") for _ in range(qlen))
        a.query_qualities = pysam.qualitystring_to_array("I" * qlen)
        a.next_reference_id, a.next_reference_start = -1, -1
        reads.append(a)

    for i in range(40):
        read(f"rep6_{i}", 42131000 + i * 20, "150M", 0)
    for i in range(40):
        read(f"body_{i}", 42127000 + i * 50, "100M", 10 if i % 3 == 0 else 60)
    read("dup", 42128000, "150M", 60, 0x400)
    read("sec", 42128000, "150M", 60, 0x100)
    read("del", 42126400, "50M200D50M", 60)
    read("clip", 42123100, "20S130M", 60)
    read("pair", 42129000, "150M", 60, 0x1 | 0x2 | 0x40)
    read("pair", 42129050, "150M", 60, 0x1 | 0x2 | 0x80 | 0x10)
    bam_from(reads, os.path.join(OUT, "in", "cyp2d6_synthetic.bam"), h)


def cyp2d6_cases():
    case("cyp2d6-bed", "cyp2d6-bed")
    case("cyp2d6-bed-stray", "cyp2d6-bed", "stray")
    for bam in ("{d}/cyp2d6.bam", "{d}/parity/in/cyp2d6_nochr.bam", "{d}/parity/in/cyp2d6_synthetic.bam"):
        name = os.path.basename(bam)[:-4]
        case(f"measure-{name}", "cyp2d6-measure", "--bam", bam, "--out", "{o}/shares.tsv")
    case("measure-no-chr22", "cyp2d6-measure", "--bam", "{d}/parity/in/cyp2d6_other.bam", "--out", "{o}/shares.tsv")
    case("measure-missing-bam", "cyp2d6-measure", "--bam", "{d}/parity/in/none.bam", "--out", "{o}/shares.tsv")
    case("measure-no-out", "cyp2d6-measure", "--bam", "{d}/cyp2d6.bam")
    shares = write("shares.tsv", "")  # filled by measure, below
    for c in ("clean_no_alt", "multimapped_with_alt", "deleted_both", "deleted_one", "no_flank_reads"):
        args = ["--all", f"{DEPTH}/{c}.q0.bed", "--mapq1", f"{DEPTH}/{c}.q1.bed"]
        case(f"check-{c}", "cyp2d6-check", *args, "--out", "{o}/check.tsv")
        case(f"check-{c}-measured", "cyp2d6-check", *args, "--measure", shares, "--out", "{o}/check.tsv")
    clean0 = open(os.path.join(REPO, "tests/fixtures/pgx/depth/clean_no_alt.q0.bed")).read()
    q1 = f"{DEPTH}/clean_no_alt.q1.bed"

    def variant(name, text, gz=False):
        case(f"check-{name}", "cyp2d6-check", "--all", write(f"{name}.q0.bed" + (".gz" if gz else ""), text, gz),
             "--mapq1", q1, "--out", "{o}/check.tsv")

    variant("gz", clean0, gz=True)
    variant("comments-and-blank-lines", "# mosdepth\n\n" + clean0 + "   \n")
    variant("crlf", clean0.replace("\n", "\r\n"))
    variant("extra-named-rows", clean0 + "chr22\t1\t2\tother\t5\n")
    lines = clean0.splitlines(keepends=True)
    # A region listed twice: the same depth, then two depths (the Python
    # averages them; og refuses the file, a golden under og/).
    variant("duplicate-row-same", lines[0] + lines[1] + lines[1] + lines[2])
    variant("duplicate-row", lines[0] + lines[1] + lines[1].replace("\t18.53", "\t2.00") + lines[2])
    variant("flank-names", clean0.replace("\tflank\t", "\tFlank_left\t", 1))
    variant("no-gene", lines[0] + lines[2])
    variant("no-flank", lines[1])
    variant("other-region", clean0.replace("42050000", "42050001"))
    variant("zero-length", lines[0] + lines[1] + lines[2] + "chr22\t42123192\t42123192\tCYP2D6\t0\n")
    variant("short-row", clean0 + "chr22\t1\t2\tCYP2D6\n")
    variant("bad-number", clean0.replace("18.53", "x"))
    variant("nan-mean", clean0.replace("18.53", "nan"))
    variant("empty", "")
    case("check-bad-measure", "cyp2d6-check", "--all", f"{DEPTH}/clean_no_alt.q0.bed", "--mapq1", q1,
         "--measure", write("bad_measure.tsv", "a\tb\tc\n"), "--out", "{o}/check.tsv")
    case("check-measure-zero", "cyp2d6-check", "--all", f"{DEPTH}/clean_no_alt.q0.bed", "--mapq1", q1, "--measure",
         write("zero_measure.tsv", "region\tmin_mapq\tmean_depth\ngene\t0\t0.0000\ngene\t20\t0.0000\n"
               "body\t0\t0.0000\nbody\t1\t0.0000\n\n"), "--out", "{o}/check.tsv")
    case("check-missing-input", "cyp2d6-check", "--all", "{d}/parity/in/none.bed", "--mapq1", q1, "--out", "{o}/check.tsv")
    case("check-no-mapq1", "cyp2d6-check", "--all", q1, "--out", "{o}/check.tsv")
    case("check-out-unwritable", "cyp2d6-check", "--all", q1, "--mapq1", q1, "--out", "{o}/no/such/dir/check.tsv")


# --- RFC1 and the other motif loci ----------------------------------------------
S2, E2 = 39350999, 39351024  # a made-up FAME-like locus on chr4: (TTTTA)5


def motif_inputs():
    for kind in ("benign", "pathogenic", "mixed", "acagg", "unknown", "few"):
        t.write_bam(pysam, os.path.join(OUT, "in", f"rfc1_{kind}.bam"), kind)
    vcf = sites_only(t.VCF)
    write("rfc1.vcf", vcf)
    write("rfc1.vcf.gz", vcf, gz=True)
    write("no_rfc1.vcf", vcf.replace("REPID=RFC1", "REPID=HTT").replace("VARID=RFC1", "VARID=HTT"))
    # Reads at the second locus: TTTCA expanded on most, TTTTA on some, the
    # reference allele on others, each spanning with its flanks.
    rng = random.Random(5)
    left = "".join(rng.choice("ACGT") for _ in range(300))
    right = "".join(rng.choice("ACGT") for _ in range(300))
    reads = []
    h = {"HD": {"VN": "1.6", "SO": "coordinate"},
         "SQ": [{"SN": "chr4", "LN": 190214555}, {"SN": "chr1", "LN": 248956422}]}
    head = pysam.AlignmentHeader.from_dict(h)
    with pysam.AlignmentFile(os.path.join(OUT, "in", "rfc1_pathogenic.bam")) as b:
        for r in b.fetch(until_eof=True):
            reads.append(pysam.AlignedSegment.from_dict(r.to_dict(), head))

    def spanning(name, unit, units):
        rep = unit * units
        ins = len(rep) - (E2 - S2)
        rlen = 150 - 40 - len(rep)
        a = pysam.AlignedSegment(head)
        a.query_name, a.flag, a.reference_id, a.reference_start = name, 0, 0, S2 - 40
        a.cigarstring = f"40M{E2 - S2}M" + (f"{ins}I" if ins else "") + f"{rlen}M"
        a.query_sequence = left[-40:] + rep + right[:rlen]
        a.query_qualities = pysam.qualitystring_to_array("I" * 150)
        a.mapping_quality = 60
        a.next_reference_id, a.next_reference_start = -1, -1
        return a

    reads += [spanning(f"fame{i}", "TTTCA", 16) for i in range(5)]
    reads += [spanning(f"famb{i}", "TTTTA", 14) for i in range(2)]
    reads += [spanning(f"famr{i}", "TTTTA", 5) for i in range(4)]
    bam_from(reads, os.path.join(OUT, "in", "two_loci.bam"), h)
    vcf_head = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    rec = lambda chrom, pos, end, rid, status: (  # noqa: E731
        f"{chrom}\t{pos}\t.\tA\t<STR40>\t.\tPASS\tEND={end};VARID={rid};REPID={rid};STR_STATUS={status}\n")
    rfc1_line = vcf.splitlines(keepends=True)[-1]
    write("two_loci.vcf", vcf_head + rfc1_line + rec("chr4", S2 + 1, E2, "FAMEX", "normal,pre_mutation")
          + rec("chr4", 1000, 1010, "QUIET", "normal") + rec("chr4", 2000, 2010, "NOMOTIF", "full_mutation"))
    write("none_flagged.vcf", vcf_head + rfc1_line + rec("chr4", S2 + 1, E2, "FAMEX", "normal"))
    cat = [
        {"LocusId": "RFC1", "ReferenceRegion": [f"chr4:{t.START}-{t.START + 25}", f"chr4:{t.START + 25}-{t.END}"],
         "BenignMotifs": ["AAAAG", "AAAGG"], "PathogenicMotifs": ["AAGGG", "ACAGG"], "Disease": "CANVAS",
         "InheritanceMode": "AR"},
        {"LocusId": "FAMEX", "ReferenceRegion": f"chr4:{S2}-{E2}", "BenignMotifs": ["TTTTA"],
         "PathogenicMotifs": ["TTTCA"], "Disease": "FAME", "InheritanceMode": "AD"},
        {"LocusId": "QUIET", "ReferenceRegion": "chr4:999-1010", "PathogenicMotifs": ["CAG"]},
        {"LocusId": "NOMOTIF", "ReferenceRegion": "chr4:1999-2010", "PathogenicMotifs": []},
        {"LocusId": "ABSENT", "ReferenceRegion": "chr9:1-10", "PathogenicMotifs": ["GAA"], "Disease": None},
    ]
    write("catalog.json", json.dumps(cat, indent=2))
    # The catalog's RFC1 as a recessive locus with no Disease, and an
    # unknown contig.
    cat2 = [dict(cat[0], Disease=None), dict(cat[1], ReferenceRegion="chrUn:5-30")]
    del cat2[0]["Disease"]
    write("catalog_nodisease.json", json.dumps(cat2))
    write("catalog_disease_null.json", json.dumps([dict(cat[0], Disease=None)]))
    write("catalog_bad_region.json", json.dumps([dict(cat[0], ReferenceRegion="chr4-1-2")]))
    write("catalog_no_id.json", json.dumps([{k: v for k, v in cat[0].items() if k != "LocusId"}]))
    write("catalog_broken.json", "[{")
    # A BAM without chr4.
    h1 = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 248956422}]}
    bam_from([], os.path.join(OUT, "in", "no_chr4.bam"), h1)


def motif_cases():
    vcf = "{d}/parity/in/rfc1.vcf"
    out = ["--out", "{o}/motif.tsv"]
    for kind in ("benign", "pathogenic", "mixed", "acagg", "unknown", "few"):
        case(f"motif-{kind}", "str-motif", "--bam", f"{{d}}/parity/in/rfc1_{kind}.bam", "--vcf", vcf, *out)
    p = "{d}/parity/in/rfc1_pathogenic.bam"
    case("motif-pathogenic-no-softclip", "str-motif", "--bam", p, "--vcf", vcf, "--no-softclip", *out)
    case("motif-pathogenic-no-mates", "str-motif", "--bam", p, "--vcf", vcf, "--no-mates", *out)
    case("motif-pathogenic-neither", "str-motif", "--bam", p, "--vcf", vcf, "--no-mates", "--no-softclip", *out)
    case("motif-vcf-gz", "str-motif", "--bam", p, "--vcf", "{d}/parity/in/rfc1.vcf.gz", *out)
    case("motif-skip", "str-motif", "--vcf", vcf, "--skip", "no BAM", *out)
    case("motif-skip-and-bam", "str-motif", "--bam", p, "--vcf", vcf, "--skip", "a reason", *out)
    case("motif-no-rfc1-record", "str-motif", "--bam", p, "--vcf", "{d}/parity/in/no_rfc1.vcf", *out)
    case("motif-no-contig", "str-motif", "--bam", "{d}/parity/in/no_chr4.bam", "--vcf", vcf, *out)
    case("motif-missing-bam", "str-motif", "--bam", "{d}/parity/in/none.bam", "--vcf", vcf, *out)
    case("motif-missing-vcf", "str-motif", "--bam", p, "--vcf", "{d}/parity/in/none.vcf", *out)
    case("motif-no-bam-no-skip", "str-motif", "--vcf", vcf, *out)
    case("motif-empty-skip", "str-motif", "--vcf", vcf, "--skip", "", *out)
    case("motif-flagged-no-catalog", "str-motif", "--bam", p, "--vcf", vcf, "--flagged", *out)
    two, tv = "{d}/parity/in/two_loci.bam", "{d}/parity/in/two_loci.vcf"
    cat = "{d}/parity/in/catalog.json"
    case("motif-catalog-rfc1", "str-motif", "--bam", two, "--vcf", tv, "--catalog", cat, *out)
    case("motif-catalog-flagged", "str-motif", "--bam", two, "--vcf", tv, "--catalog", cat, "--flagged", *out)
    case("motif-catalog-flagged-skip", "str-motif", "--vcf", tv, "--catalog", cat, "--flagged", "--skip", "no BAM", *out)
    case("motif-catalog-none-flagged", "str-motif", "--bam", two, "--vcf", "{d}/parity/in/none_flagged.vcf",
         "--catalog", cat, "--flagged", *out)
    case("motif-catalog-disease-null", "str-motif", "--bam", two, "--vcf", tv,
         "--catalog", "{d}/parity/in/catalog_disease_null.json", *out)
    case("motif-catalog-nodisease", "str-motif", "--bam", two, "--vcf", tv,
         "--catalog", "{d}/parity/in/catalog_nodisease.json", *out)
    case("motif-catalog-unknown-contig", "str-motif", "--bam", two, "--vcf", tv,
         "--catalog", "{d}/parity/in/catalog_nodisease.json", "--flagged", *out)
    for bad in ("bad_region", "no_id", "broken"):
        case(f"motif-catalog-{bad.replace('_', '-')}", "str-motif", "--bam", two, "--vcf", tv,
             "--catalog", f"{{d}}/parity/in/catalog_{bad}.json", *out)


def run(tool, argv):
    buf = io.StringIO()
    rc = 0
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
        try:
            if tool == "str-motif":
                rc = rm.main(argv)
            else:
                rc = dc.main([tool.split("-", 1)[1]] + argv)
        except SystemExit as e:  # argparse
            rc = e.code
    return rc, buf.getvalue()


def main():
    # og/ (what og writes where it differs on purpose) is written by hand.
    for d in ("in", "want"):
        shutil.rmtree(os.path.join(OUT, d), ignore_errors=True)
    os.makedirs(os.path.join(OUT, "in"))
    cyp2d6_inputs()
    motif_inputs()
    cyp2d6_cases()
    motif_cases()
    # The measure file the check cases read: measure on the fixture slice.
    with tempfile.TemporaryDirectory() as tmp:
        assert dc.main(["measure", "--bam", os.path.join(DATA, "cyp2d6.bam"), "--out", f"{tmp}/s.tsv"]) == 0
        shutil.copy(f"{tmp}/s.tsv", os.path.join(OUT, "in", "shares.tsv"))
    rows = []
    for name, tool, args in CASES:
        with tempfile.TemporaryDirectory() as tmp:
            argv = [a.replace("{d}", DATA).replace("{r}", REPO).replace("{o}", tmp) for a in args]
            rc, stdout = run(tool, argv)
            d = os.path.join(OUT, "want", name)
            os.makedirs(d)
            with open(os.path.join(d, "exit"), "w") as fh:
                fh.write(f"{rc}\n")
            with open(os.path.join(d, "stdout"), "w") as fh:
                fh.write(stdout)
            for root, _, files in os.walk(tmp):
                for f in files:
                    rel = os.path.relpath(os.path.join(root, f), tmp)
                    dest = os.path.join(d, "files", rel)
                    os.makedirs(os.path.dirname(dest), exist_ok=True)
                    shutil.copy(os.path.join(root, f), dest)
        rows.append("\t".join([name, tool] + args))
    with open(os.path.join(OUT, "cases.tsv"), "w") as fh:
        fh.write("# case<TAB>tool<TAB>arguments...: made by ../../make_parity.py; {d} tests/data, {r} the repository, "
                 "{o} the output directory\n")
        fh.write("".join(r + "\n" for r in rows))
    print(f"{len(rows)} cases")


if __name__ == "__main__":
    main()
