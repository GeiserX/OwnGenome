#!/usr/bin/env bash
# sex-loss.sh: writes the mosdepth windows and VCFs of the qc-sex-loss
# parity cases in the current directory, the same every time: the synthetic
# genomes of the old unit test of bin/mosaic_sex_loss.py (a man, 30% of his
# Y lost, 47,XYY, drifting depth, a sliced BAM, contigs without chr; a woman
# with and without one X lost in 30% of her cells, shifted allele fractions
# with no depth change, a noisy library), a truncated gzip and a VCF with a
# read count past 32 bits. scripts/ci/oracle-diff.sh runs it.
set -euo pipefail
export LC_ALL=C
python3 - <<'PY'
import gzip
import random


def windows(f, chrom, start, end, depth, rnd, prefix="chr"):
    for s in range(start - start % 500, end, 500):
        d = max(0.0, rnd.gauss(depth, depth * 0.1)) if depth else 0.0
        f.write(f"{prefix}{chrom}\t{s}\t{s + 500}\t{d:.2f}\n")


def regions(path, sex, x_depth=None, y_depth=None, seed=1, prefix="chr", slice_only=False):
    rnd = random.Random(seed)
    with gzip.open(path, "wt") as f:
        windows(f, "1", 1_000_000, 4_000_000, 30, rnd, prefix)
        windows(f, "1", 4_000_000, 4_100_000, 0, rnd, prefix)
        windows(f, "2", 5_000_000, 6_000_000, 30, rnd, prefix)
        windows(f, "2", 6_000_000, 6_050_000, 400, rnd, prefix)
        if slice_only:
            windows(f, "Y", 2_790_000, 2_800_000, 15, rnd, prefix)
            return
        xd = x_depth if x_depth is not None else (15 if sex == "male" else 30)
        windows(f, "X", 1_000_000, 2_000_000, 30, rnd, prefix)
        windows(f, "X", 73_200_000, 74_500_000, xd, rnd, prefix)
        windows(f, "X", 89_500_000, 90_500_000, 30 if sex == "male" else 0, rnd, prefix)
        if sex == "male":
            yd = y_depth if y_depth is not None else 15
            windows(f, "Y", 0, 2_781_479, 0, rnd, prefix)
            windows(f, "Y", 2_781_479, 3_050_044, yd, rnd, prefix)
            windows(f, "Y", 3_050_044, 6_235_111, 15, rnd, prefix)
            windows(f, "Y", 6_748_713, 7_574_481, yd, rnd, prefix)
            windows(f, "Y", 7_574_481, 9_046_914, 60, rnd, prefix)
            windows(f, "Y", 10_316_944, 10_544_039, 900, rnd, prefix)
            windows(f, "Y", 26_673_214, 56_673_214, 0, rnd, prefix)


def drift_regions(path, y_level, seed=5):
    rnd = random.Random(seed)
    with gzip.open(path, "wt") as f:
        for chrom, start in (("1", 1_000_000), ("2", 5_000_000)):
            for s in range(start, start + 10_000_000, 200_000):
                windows(f, chrom, s, s + 200_000, 30 * rnd.gauss(1, 0.08), rnd)
        windows(f, "Y", 2_781_479, 3_000_000, 15 * y_level, rnd)


def vcf(path, d=0.0, spread=0.0, seed=2, extra=()):
    rnd = random.Random(seed)
    lines = ["##fileformat=VCFv4.2\n", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"]

    def site(chrom, pos, p):
        n = rnd.randint(20, 40)
        p = min(0.95, max(0.05, p + rnd.gauss(0, spread))) if spread else p
        a = sum(rnd.random() < p for _ in range(n))
        lines.append(f"{chrom}\t{pos}\t.\tA\tG\t40\tPASS\t.\tGT:AD:DP\t0/1:{n - a},{a}:{n}\n")
    for i in range(3000):
        site("chr1", 1_000_000 + 1000 * i, 0.5)
    for i in range(600):
        site("chrX", 73_200_000 + 2000 * i, 0.5 + (d if rnd.random() < 0.5 else -d))
    lines += ["chrX\t73100000\t.\tA\tG\t40\tPASS\t.\tGT:AD\t1/1:0,30\n",
              "chrX\t73100100\t.\tA\tG\t40\tRefCall\t.\tGT:AD\t0/1:5,25\n",
              "chrX\t73100200\t.\tAT\tA\t40\tPASS\t.\tGT:AD\t0/1:5,25\n",
              "chrX\t1500000\t.\tA\tG\t40\tPASS\t.\tGT:AD\t0/1:2,28\n",
              "chrX\t73100300\t.\tA\tG\t40\tPASS\t.\tGT:AD\t0/1:1,5\n"]
    lines += list(extra)
    with gzip.open(path, "wt") as f:
        f.writelines(lines)


regions("m.bed.gz", "male")
regions("m30.bed.gz", "male", y_depth=10.5)
regions("xyy.bed.gz", "male", y_depth=30)
regions("f.bed.gz", "female")
regions("f30.bed.gz", "female", x_depth=30 * 0.85)
regions("slice.bed.gz", "male", slice_only=True)
regions("nochr.bed.gz", "male", y_depth=10.5, prefix="")
drift_regions("drift.bed.gz", 1.13)
d30 = 0.3 / (2 * (2 - 0.3))
vcf("f.vcf.gz")
vcf("f30.vcf.gz", d=d30)
vcf("fcn.vcf.gz", d=d30)
vcf("noisy.vcf.gz", spread=0.06)
vcf("huge.vcf.gz", extra=["chr1\t999\t.\tA\tG\t40\tPASS\t.\tGT:AD\t0/1:0,4294967296\n"])
with open("m.bed.gz", "rb") as f:
    data = f.read()
with open("trunc.bed.gz", "wb") as f:
    f.write(data[: len(data) // 2])
PY
