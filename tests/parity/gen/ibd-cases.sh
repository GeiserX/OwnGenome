#!/usr/bin/env bash
# Inputs of og family ibd's parity cases, one folder each, as
# tests/test_ibd_report.py builds them: Beagle's phased VCFs
# (phased/<chrom>.vcf.gz), hap-ibd's segments (ibd/<chrom>.ibd.gz), PLINK
# maps (maps/<chrom>.map), chroms.txt and pedigree.tsv.
#   one/      four members on chr1: a parent and child, full siblings, a
#             short segment in hap-ibd's other order, a member without a
#             pedigree row
#   two/      one/ plus chr2
#   stitch/   stretches joined across one opposite homozygote, not two; the
#             chromosome ends reached the same way
#   none/     no chromosome phased
#   7col/ otherchrom/ nosample/ mapbp/   damaged inputs (exit 2)
#   nomap/    chroms.txt names chr7, which has no map (exit 2)
#   short/    a phased record with fewer sample columns than the header
#   twice/    chroms.txt names chr1 twice
#   bench/    22 chromosomes of 50,000 markers, five members (the speed case)
# gzip files are written with mtime 0, so a run is reproducible.
set -euo pipefail
python3 - <<'PY'
import gzip
import os
import random

MB = 1000000
GT = {0: "0|0", 1: "0|1", 2: "1|1"}
PEDIGREE = "sample\tsex\tfather\tmother\tpartner\nA\tmale\tB\tX\t.\nC\tfemale\tB\tX\t.\nB\tmale\t.\t.\tX\n"


def gz(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as f:
        f.write(text.encode())


def put(path, text):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as f:
        f.write(text)


def chrom(d, c, samples, first, last, rows, map_end, dos=lambda s, mb: 1, step=MB):
    put(f"{d}/maps/{c}.map", "".join(f"{c}\t.\t{i:.4f}\t{max(1, i * MB)}\n" for i in range(map_end + 1)))
    head = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"
    body = "".join(f"{c}\t{p}\t.\tA\tG\t.\tPASS\t.\tGT\t" + "\t".join(GT[dos(s, p)] for s in samples) + "\n"
                   for p in range(first, last + 1, step))
    gz(f"{d}/phased/{c}.vcf.gz", head + body)
    gz(f"{d}/ibd/{c}.ibd.gz", "".join("\t".join(str(x) for x in r) + "\n" for r in rows))


def case(d, chroms, pedigree=PEDIGREE):
    put(f"{d}/chroms.txt", "".join(c + "\n" for c in chroms))
    put(f"{d}/pedigree.tsv", pedigree)


S4 = ["A", "B", "C", "D"]
DOS1 = lambda s, p: (1 if 10 * MB <= p <= 13 * MB else 2) if s == "D" else 0  # noqa: E731
CHR1 = [("A", 1, "B", 1, "chr1", 1 * MB, 50 * MB, 49.0), ("A", 1, "B", 2, "chr1", 50 * MB, 99 * MB, 49.0),
        ("A", 1, "C", 1, "chr1", 1 * MB, 99 * MB, 98.0), ("A", 2, "C", 2, "chr1", 20 * MB, 44500000, 24.5),
        ("D", 2, "B", 1, "chr1", 10 * MB, 13 * MB, 3.0)]
for d in ("one", "two"):
    chrom(d, "chr1", S4, MB, 99 * MB, CHR1, 100, DOS1)
case("one", ["chr1"])
chrom("two", "chr2", S4, 2 * MB, 49 * MB, [("B", 2, "A", 2, "chr2", 2 * MB, 49 * MB, 47.0)], 50, DOS1)
case("two", ["chr1", "chr2"])


def dos9(s, p):
    mb = p // MB
    return {"Q": 2 if mb in (62, 66) else 0, "R": 2 if mb <= 3 or mb >= 51 else 0, "S": 2 if mb == 99 else 0}.get(s, 0)


chrom("stitch", "chr3", ["P", "Q", "R", "S"], MB, 99 * MB,
      [("P", 1, "Q", 1, "chr3", 1 * MB, 30 * MB, 29.0), ("P", 1, "Q", 2, "chr3", 35 * MB, 60 * MB, 25.0),
       ("P", 1, "Q", 1, "chr3", 70 * MB, 99 * MB, 29.0), ("P", 2, "R", 1, "chr3", 5 * MB, 50 * MB, 45.0),
       ("P", 2, "S", 2, "chr3", 3 * MB, 97 * MB, 94.0)], 100, dos9)
case("stitch", ["chr3"], "sample\tsex\tfather\tmother\tpartner\nP\tmale\t.\t.\t.\n")
case("none", [])

for d, rows in (("7col", [("A", 1, "B", 1, "chr1", MB, 5 * MB)]),
                ("otherchrom", [("A", 1, "B", 1, "chr9", MB, 5 * MB, 4.0)]),
                ("nosample", [("A", 1, "Z", 1, "chr1", MB, 5 * MB, 4.0)]),
                ("mapbp", CHR1), ("nomap", CHR1), ("short", CHR1), ("twice", CHR1)):
    chrom(d, "chr1", S4, MB, 99 * MB, rows, 100)
    case(d, ["chr1"])
with open("mapbp/maps/chr1.map", "a") as f:
    f.write("chr1\t.\t3.0\n")
case("nomap", ["chr1", "chr7"])
case("twice", ["chr1", "chr1"])
vcf = gzip.open("short/phased/chr1.vcf.gz", "rt").read().split("\n")
vcf[5] = vcf[5].rsplit("\t", 1)[0]
gz("short/phased/chr1.vcf.gz", "\n".join(vcf))

# The speed case: a parent and child, two siblings and an unrelated member,
# genotypes drawn per marker, hap-ibd segments broken every 20 cM or so.
rnd = random.Random(38)
S5 = ["DAD", "KID", "SIB", "MUM", "OTHER"]
chroms = ["chr%d" % i for i in range(1, 23)]
for c in chroms:
    first, step, n = MB, 1000, 50000
    pos = range(first, first + step * n, step)
    rows = []
    for a, b, ha, hb in (("DAD", "KID", 1, 1), ("DAD", "SIB", 2, 1), ("KID", "SIB", 1, 1), ("KID", "SIB", 2, 2),
                         ("MUM", "KID", 1, 2)):
        s = first
        while s < first + step * n - 5 * MB:
            e = min(s + rnd.randint(8, 30) * MB, first + step * (n - 1))
            rows.append((a, ha, b, hb, c, s, e, (e - s) / MB))
            s = e + rnd.randint(0, 3) * MB
    g = {(s, p): rnd.choice((0, 1, 1, 2)) for s in S5 for p in pos}
    chrom("bench", c, S5, first, first + step * (n - 1), rows, 52, lambda s, p: g[s, p], step)
case("bench", chroms, "sample\tsex\tfather\tmother\tpartner\nKID\tmale\tDAD\tMUM\t.\nSIB\tfemale\tDAD\tMUM\t.\n")
PY
