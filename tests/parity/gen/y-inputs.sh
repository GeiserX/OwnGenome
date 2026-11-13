#!/usr/bin/env bash
# The inputs of the Y cases of step 37 (y-tree-positions, y-tree-place,
# y-private-snps, y-yleaf-positions), all synthetic:
#   t/        the tree and tables of tests/test_y_tree_place.py, with a
#             pseudoautosomal SNP, and one pileup per case of that test
#             (t-<case>.pileup), plus the read-mark and contig cases
#   big/      a random tree of 4,000 branches with aliases, names on both
#             tables, invalid table rows and PAR positions; pileups along a
#             few lineages with no-calls, mixes and read marks (big-<n>.pileup)
#   odd/      a tree with branches without an id, ties broken by id, escapes,
#             extra keys and falsy values
#   p/        the private-SNP tables, a gzipped and a plain mask, a mask with
#             two records, and the VCFs (the test's records, and 3,000 random)
#   yd/       Yleaf's data folder: hg38/new_positions.txt with a repeated
#             position, and hg_prediction_tables/tree.json; yd-notree/ lacks
#             the tree
set -euo pipefail
python3 - <<'PY'
import gzip
import json
import os
import random


def put(path, text, gz=False):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with (gzip.open(path, "wt") if gz else open(path, "w")) as f:
        f.write(text)


def table(rows):
    return "".join(f"chry\t{n}\tR-{n}\t{p}\t{a}->{d}\t{a}\t{d}\n" for n, (p, a, d) in rows.items())


# --- t/: the unit test's tree ---------------------------------------------------------
TREE = {"id": "", "snps": "", "children": [
    {"id": "A", "snps": "a1, a2/a2alt, par", "children": [
        {"id": "B", "snps": "b1, b2, b3, bnot", "children": [
            {"id": "C", "snps": "c1, c2", "children": []}]},
        {"id": "X", "snps": "x1", "children": [
            {"id": "Y", "snps": "y1", "children": [
                {"id": "Z", "snps": "z1", "children": [
                    {"id": "W", "snps": "w1", "children": []}]}]}]},
        {"id": "P", "snps": "p1", "children": [
            {"id": "Q", "snps": "q1, q1b", "children": []}]}]}]}
FTDNA = {"a1": (3000100, "C", "T"), "a2alt": (3000200, "G", "A"), "b1": (3000300, "A", "G"),
         "b2": (3000400, "T", "C"), "b3": (3000500, "C", "G"), "c1": (3000600, "G", "T"),
         "x1": (3000800, "A", "C"), "y1": (3000900, "T", "A"), "z1": (3001000, "G", "C"),
         "w1": (3001100, "C", "A"), "p1": (3001200, "T", "G"), "q1": (3001300, "A", "T"),
         "q1b": (3001300, "A", "T"), "par": (2000000, "A", "G")}
YFULL = {"a1": (3999999, "C", "T"), "c2": (3000700, "A", "T")}
put("t/current_tree.json", json.dumps(TREE))
put("t/positions_ftdna_hg38.txt", table(FTDNA))
put("t/positions_yfull_hg38.txt", table(YFULL))
put("t/version.txt", "tree\tYFull YTree test\npositions\ttest tables\n")


def pos_of(name):
    return FTDNA.get(name) or YFULL[name]


def pileup(path, states, depth=12):
    out = []
    for name, st in states.items():
        dp = depth
        if isinstance(st, tuple):
            st, dp = st
        pos, anc, der = pos_of(name)
        if st == "mix":
            bases = anc * (dp // 2) + der * (dp - dp // 2)
        else:
            bases = (der if st == "D" else anc) * dp
        out.append(f"chrY\t{pos}\tN\t{dp}\t{bases}\t{'I' * dp}\n")
    put(path, "".join(out))


A = {"a1": "D", "a2alt": "D"}
B = {"b1": "D", "b2": "D", "b3": "D"}
for case, states in {
    "deep": {**A, **B, "c1": "D", "c2": "D", "x1": "A"},
    "split": {**A, "b1": "D", "b2": "D", "b3": "A", "x1": "A"},
    "nocall": {**A, **B, "c1": ("D", 9), "c2": "mix", "x1": "A"},
    "blocked": {**A, **B, "x1": "A", "y1": "D"},
    "outvoted": {**A, **B, "x1": "A", "y1": "D", "z1": "D"},
    "unblocked": {**A, **B, "y1": "D", "z1": "D"},
    "twoanc": {"a1": "A", "a2alt": "A", **B, "c1": "D", "c2": "D"},
    "twosingles": {**A, **B, "x1": "A", "y1": "A", "z1": "D", "w1": "D"},
    "onesite": {**A, **B, "p1": "A", "q1": "D"},
    "none": {"a1": "A", "a2alt": "A"},
}.items():
    pileup(f"t-{case}.pileup", states)
# Read marks, indels, lower case, '.' and ',' against the reference column,
# a repeated position (the last line wins), and lines a reader skips.
marks = [f"chrY\t{pos_of(n)[0]}\tN\t12\t{pos_of(n)[2] * 12}\t{'I' * 12}\n" for n in ("a1", "a2alt", "b1", "b2", "b3")]
marks += [f"chrY\t{pos_of('c1')[0]}\tt\t5\t^I..,$+2ACt-1Ga*\tIIIII\n",
          f"chrY\t{pos_of('c2')[0]}\tt\t5\t^I.$,+2GG,t-1A.\tIIIII\n",
          f"chrY\t{pos_of('x1')[0]}\ta\t12\tCCCCCCCCCCCC\tI\n",
          f"chrY\t{pos_of('x1')[0]}\ta\t12\t,,,,,.......\tI\n",
          "chrY\tx3000800\tA\t1\tA\tI\n", "short\tline\n", "\n"]
put("t-marks.pileup", "".join(marks))
# The pileup of blocked (x1 ancestral blocks Y: B), with a chr1 line at x1's
# position that says derived: read as chrY it would unblock Y.
put("t-contigs.pileup", open("t-blocked.pileup").read()
    + f"chr1\t{pos_of('x1')[0]}\tN\t12\t{pos_of('x1')[2] * 12}\t{'I' * 12}\n")

# --- big/: a random tree -----------------------------------------------------------------
rng = random.Random(37)
names = [f"S{i}" for i in range(12000)]
rng.shuffle(names)
ni = iter(names)
nodes = [{"id": "", "snps": "", "children": []}]
parents = [0]
for i in range(1, 4000):
    parent = rng.choice(parents[-60:] if rng.random() < 0.8 else parents)
    groups = []
    for _ in range(rng.choice([0, 1, 1, 2, 3, 5])):
        g = [next(ni)]
        if rng.random() < 0.2:
            g.append(next(ni))
        groups.append("/".join(g))
    node = {"id": f"N{rng.randint(0, 999):03d}-{i}" if rng.random() < 0.97 else "", "snps": ", ".join(groups),
            "tmrca": rng.randint(0, 9000), "children": []}
    nodes[parent]["children"].append(node)
    nodes.append(node)
    parents.append(i)
put("big/current_tree.json", json.dumps(nodes[0]))
ft, yf = [], []
for n in names:
    r = rng.random()
    pos = rng.randint(2_700_000, 57_000_000)
    a, d = rng.sample("ACGT", 2)
    row = f"chry\t{n}\tR\t{pos}\t{a}->{d}\t{a}\t{d}\n"
    if r < 0.55:
        ft.append(row)
    elif r < 0.8:
        yf.append(row)
    elif r < 0.85:
        ft.append(row)
        b, e = rng.sample("ACGT", 2)
        yf.append(f"chry\t{n}\tR\t{rng.randint(3_000_000, 50_000_000)}\t{b}->{e}\t{b}\t{e}\n")
    elif r < 0.88:
        ft.append(f"chry\t{n}\tR\t{pos}\t{a}->{a}\t{a}\t{a}\n")
        yf.append(row)
    elif r < 0.9:
        ft.append(f"chry\t{n}\tR\t{pos}x\t{a}->{d}\t{a}\t{d}\n")
    elif r < 0.92:
        ft.append(f"chry\t{n}\tR\t{pos}\t{a}->N\t{a}\tN\n")
    elif r < 0.93:
        ft.append(f"chry\t{n}\tR\t{pos}\n")
    elif r < 0.95:
        yf.append(f"chry\t{n}\tR\t{pos}\t{a}->{d}\t{a.lower()}\t{d}\n")
    # the rest are on no table
rng.shuffle(ft)
put("big/positions_ftdna_hg38.txt", "".join(ft) + "chry\tS1\tR\t2781479\tA->G\tA\tG\n")
put("big/positions_yfull_hg38.txt", "".join(yf) + "chry\tS2\tR\t56887903\tA->G\tA\tG\r\n")
put("big/version.txt", "tree\tsynthetic tree\r\npositions\tsynthetic tables\n")
place = {}
for line in ft + yf:
    c = line.rstrip("\n").split("\t")
    if len(c) >= 7 and c[1] not in place and c[3].isdigit() and c[5] in "ACGT" and c[6] in "ACGT" \
            and len(c[5]) == len(c[6]) == 1 and c[5] != c[6]:
        place[c[1]] = (int(c[3]), c[5], c[6])
for k in range(4):
    leaf = rng.randrange(1, len(nodes))
    on = set()
    # parents[] holds no parent links; find the path by walking the tree.
    stack = [(nodes[0], [0])]
    path = []
    index = {id(nd): i for i, nd in enumerate(nodes)}
    while stack:
        nd, p = stack.pop()
        if index[id(nd)] == leaf:
            path = p
            break
        for ch in nd["children"]:
            stack.append((ch, p + [index[id(ch)]]))
    on = set(path)
    out = {}
    for i, nd in enumerate(nodes):
        for g in nd["snps"].split(","):
            for n in (x.strip() for x in g.split("/")):
                if n in place and place[n][0] not in out:
                    pos, a, d = place[n]
                    st = "D" if i in on else ("A" if rng.random() < 0.96 else "D")
                    r = rng.random()
                    dp = rng.randint(8, 30)
                    if r < 0.04:
                        bases = a * (dp // 2) + d * (dp - dp // 2)
                    elif r < 0.06:
                        dp = 9
                        bases = (d if st == "D" else a) * dp
                    elif r < 0.08:
                        continue
                    else:
                        b = d if st == "D" else a
                        bases = "".join(b.lower() if rng.random() < 0.5 else b for _ in range(dp))
                        if rng.random() < 0.1:
                            bases = "^]" + bases + "$-1A+3ACG"
                    out[pos] = f"chrY\t{pos}\t{rng.choice('ACGTN')}\t{dp}\t{bases}\t{'I' * dp}\n"
    put(f"big-{k}.pileup", "".join(out[p] for p in sorted(out)))

# --- odd/: ids, ties, escapes and falsy values ---------------------------------------------
odd = ('{"formed": "-", "tmrca": 1.5e3, "id": null, "snps": 0, "children": ['
       '{"snps": "o1", "x": [1, {"y": [null, true]}], "children": ['
       '{"id": "M\\u00e9-b", "snps": "o2,\\u006f3"}, {"id": "M\\u00e9-a", "snps": "o4/o4b, o5"}]},'
       '{"id": "K", "snps": "", "children": null},'
       '{"id": false, "snps": "o6", "children": [{"id": "L", "snps": "o7"}]}]}')
put("odd/current_tree.json", odd)
put("odd/positions_ftdna_hg38.txt", table({"o1": (3100000, "A", "G"), "o2": (3100100, "C", "T"), "o3": (3100200, "G", "A"),
                                           "o4b": (3100300, "T", "C"), "o5": (3100400, "A", "C"),
                                           "o6": (3100500, "C", "G"), "o7": (3100600, "G", "C")}))
put("odd/positions_yfull_hg38.txt", "")
odd_states = {3100000: "G", 3100100: "T", 3100200: "A", 3100300: "C", 3100400: "C", 3100500: "G", 3100600: "C"}
put("odd.pileup", "".join(f"chrY\t{p}\tN\t10\t{b * 10}\t{'I' * 10}\n" for p, b in odd_states.items()))
put("odd-anc.pileup", "".join(f"chrY\t{p}\tN\t10\t{'A' * 10}\t{'I' * 10}\n" for p in odd_states))

# --- p/: private SNPs -----------------------------------------------------------------------
put("p/positions_ftdna_hg38.txt", "chrY\tpub1\tR\t3100100\tG->A\tG\tA\nchrY\tbad\tR\tx\tG->A\tG\tA\nshort\n")
put("p/positions_yfull_hg38.txt", "chrY\tpub2\tR\t3100200\tG->A\tG\tA\n")
MASK_LEN, Z = 3102000, (3101250, 3101350)
seq = "P" * (Z[0] - 1) + "Z" * (Z[1] - Z[0] + 1) + "P" * (MASK_LEN - Z[1])
lines = "".join(seq[i:i + 60] + "\n" for i in range(0, len(seq), 60))
put("p/callable_mask_chrY.fasta.gz", ">chrY\n" + lines, gz=True)
put("mask-plain.fa", ">chrY\r\n" + lines.replace("\n", "  \r\n"))
put("mask-two.fa.gz", ">chrY\n" + lines + ">chrX\nPPPP\n", gz=True)
HEAD = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
REC = [
    (3100100, "G", "A", "PASS", "1:0,30:30:50"), (3100200, "G", "A", "PASS", "1:0,30:30:50"),
    (3100300, "C", "T", "PASS", "1:0,30:30:50"), (3100400, "c", "t", "PASS", "1/1:1,29:30:50"),
    (2781479, "C", "T", "PASS", "1:0,30:30:50"), (2781480, "C", "T", "PASS", "1:0,30:30:50"),
    (56887903, "C", "T", "PASS", "1:0,30:30:50"), (3100500, "C", "T", "RefCall", "1:0,30:30:50"),
    (3100600, "C", "T", "PASS", "0/1:15,15:30:50"), (3100700, "C", "CTT", "PASS", "1:0,30:30:50"),
    (3100800, "C", "T,G", "PASS", "1:0,20,10:30:50"), (3100900, "C", "T", "PASS", "1:0,9:9:50"),
    (3101000, "C", "T", "PASS", "1:0,30:30:5"), (3101100, "C", "T", "PASS", "1:6,24:30:50"),
    (3101200, "C", "T", "PASS", "1:0,300:300:50"), (3101300, "C", "T", "PASS", "1:0,30:30:50"),
    (3101500, "C", "T", ".", "1|1:0,31:31:20"), (3101600, "C", "T", "", "1:0,33:33:99"),
    (3200000, "C", "T", "PASS", "1:0,30:30:50"),
    # the mask's edges: the last P before the Z run, its first and last Z,
    # the P after it, the mask's last base and the one past its end
    (3101249, "A", "G", "PASS", "1:0,30:30:50"), (3101250, "A", "G", "PASS", "1:0,30:30:50"),
    (3101350, "A", "G", "PASS", "1:0,30:30:50"), (3101351, "A", "G", "PASS", "1:0,30:30:50"),
    (3102000, "A", "G", "PASS", "1:0,30:30:50"), (3102001, "A", "G", "PASS", "1:0,30:30:50"),
]
body = "".join(f"chrY\t{p}\t.\t{r}\t{a}\t50\t{f}\t.\tGT:AD:DP:GQ\t{s}\n" for p, r, a, f, s in REC)
put("p.vcf", HEAD + body)
put("p.vcf.gz", HEAD + body, gz=True)
# even.vcf: four kept calls (a float median); dpvaf.vcf: no AD, DP and VAF.
put("even.vcf", HEAD + "".join(f"Y\t{p}\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP:GQ\t1:0,{d}:{d}:60\n"
                               for p, d in ((3100300, 30), (3100310, 31), (3100320, 40), (3100330, 25))))
put("dpvaf.vcf", HEAD + "chrY\t3100300\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ:VAF\t1:20:40:0.95\n"
    "chrY\t3100310\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ:VAF\t1:20:40:nan\n"
    "chrY\t3100320\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ:VAF\t1:20:40:.\n"
    "chrY\t3100330\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP:GQ\t1:x,3:20:40\n"
    "chrY\t3100340\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP:GQ\t1:0,0:0:40\n"
    "chrY\t3100350\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP\t1:0,20:20\n"
    "chrY\t3100360\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP:GQ\t1:0,20\n")
rng = random.Random(38)
rows = []
for i in range(3000):
    pos = rng.choice([3100100, 3100200]) if rng.random() < 0.05 else rng.randint(2_700_000, 3_110_000)
    ref = rng.choice("ACGTacgt") if rng.random() < 0.97 else "CA"
    alt = rng.choice(["A", "C", "G", "T", "t", "T,G", "<NON_REF>", "AT"])
    filt = rng.choice(["PASS"] * 8 + [".", "", "RefCall"])
    gt = rng.choice(["1", "1/1", "1|1", "0/1", ".", "1/2", "./1"])
    dp = rng.choice([rng.randint(5, 60), rng.randint(60, 300), 0])
    a = int(dp * rng.choice([1, 1, 0.95, 0.8, 0.5]))
    fmt = rng.choice(["GT:AD:DP:GQ", "GT:DP:GQ:VAF", "GT:AD:GQ", "GT:GQ"])
    vals = {"GT": gt, "AD": f"{dp - a},{a}" if rng.random() < 0.95 else rng.choice(["", ".", "3"]),
            "DP": str(dp), "GQ": rng.choice([str(rng.randint(0, 99)), "."]),
            "VAF": rng.choice([f"{a / dp if dp else 0:.4f}", "nan", "."])}
    chrom = rng.choice(["chrY"] * 20 + ["Y", "chrX"])
    rows.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\t{filt}\t.\t{fmt}\t{':'.join(vals[k] for k in fmt.split(':'))}\n")
put("rand.vcf.gz", HEAD + "".join(rows), gz=True)

# --- yd/: Yleaf's data folder ---------------------------------------------------------------
pos = "".join(f"chry\tM{i}\tR-M{i}\t{p}\tC->T\tC\tT\n" for i, p in enumerate([14549643, 10182081, 13753719, 10182081, 2655180]))
put("yd/hg38/new_positions.txt", pos + "chry\tbad\tR\tx1\tC->T\tC\tT\nshort\tline\n")
put("yd/hg_prediction_tables/tree.json", '{"ROOT (Y-Chromosome \'Adam\')": []}\n')
put("yd-notree/hg38/new_positions.txt", pos)
PY
