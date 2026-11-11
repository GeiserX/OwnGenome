#!/usr/bin/env bash
# clinvar-synthetic.sh: writes synthetic.fa (+ .fai) and synthetic.vcf in the
# current directory, the same every time (a fixed seed), for the og clinvar
# join-adjacent cases build-synthetic and prepare-synthetic: three 300,000-base
# contigs with soft-masked and N stretches, and about 120,000 records a few
# bases apart (SNVs, deletions and insertions, about one in twenty with a REF
# that disagrees with the reference), with unphased, phased, half and
# homozygous genotypes in a few phase sets. Many neighbours can be one
# allele, so the joiner's buffer, its order and its reference checks all run.
set -euo pipefail
export LC_ALL=C
python3 - <<'PY'
import random

R = random.Random(20261111)
LEN = 300000
names = ["s1", "s2", "s3"]
seqs = {}
for n in names:
    s = [R.choice("ACGT") for _ in range(LEN)]
    for _ in range(40):
        a = R.randrange(LEN - 2000)
        b = a + R.randrange(50, 1500)
        for i in range(a, b):
            s[i] = s[i].lower()
    for _ in range(5):
        a = R.randrange(LEN - 500)
        for i in range(a, a + R.randrange(10, 300)):
            s[i] = "N"
    seqs[n] = "".join(s)
with open("synthetic.fa", "w") as f, open("synthetic.fa.fai", "w") as fai:
    off = 0
    for n in names:
        head = f">{n} synthetic\n"
        off += len(head)
        fai.write(f"{n}\t{LEN}\t{off}\t70\t71\n")
        body = "".join(seqs[n][i:i + 70] + "\n" for i in range(0, LEN, 70))
        f.write(head + body)
        off += len(body)

GTS = ["0/1"] * 6 + ["1/0", "1/1", "./1", "1/.", "0/0", "0|1", "1|0", "0|1", "1|0", "1|.", ".|1", "1|1", "2/0"]
with open("synthetic.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
            + "".join(f"##contig=<ID={n},length={LEN}>\n" for n in names)
            + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
            '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG\n")
    for n in names:
        s = seqs[n].upper()
        pos = 100
        ps = 100
        while pos < LEN - 100:
            pos += R.choice([1, 1, 2, 3, 4, 5, 6, 8, 10, 11, 12, 15, 20])
            if R.random() < 0.02:
                ps = pos
            kind = R.random()
            if kind < 0.6:
                ref = s[pos - 1]
                alt = R.choice([b for b in "ACGT" if b != ref] or ["A"])
            elif kind < 0.8:
                k = R.randrange(2, 8)
                ref = s[pos - 1:pos - 1 + k]
                alt = ref[0]
            else:
                ref = s[pos - 1]
                alt = ref + "".join(R.choice("ACGT") for _ in range(R.randrange(1, 6)))
            if R.random() < 0.05:
                ref = R.choice("ACGT") + ref[1:]
            if R.random() < 0.03:
                ref = ref.lower()
            gt = R.choice(GTS)
            if "|" in gt and R.random() < 0.8:
                fmt, smp = "GT:GQ:DP:PS", f"{gt}:{R.randrange(60)}:{R.randrange(40)}:{ps}"
            else:
                fmt, smp = "GT:GQ:DP", f"{gt}:{R.randrange(60)}:{R.randrange(40)}"
            f.write(f"{n}\t{pos}\t.\t{ref}\t{alt}\t{R.randrange(1, 60)}\tPASS\t.\t{fmt}\t{smp}\n")
PY
