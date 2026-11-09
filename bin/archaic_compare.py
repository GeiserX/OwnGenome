#!/usr/bin/env python3
"""Compare step 43's Neanderthal totals with Chen et al. 2020's released calls.

Used by the archaic jobs of public-evidence.yml only (not part of the
pipeline). Standard library only.

  chen     the released per-person calls of Chen et al. (the file "Neanderthal
           sequence in 1000 genome.50kb.txt": chr, start, end, LOD, maxLOD,
           size, pop, anc, ID, tab separated), merged per person and
           chromosome: one person's base pairs and the spread of a group.
  compare  the same person run through step 43 on one or more inputs (a
           segments TSV per chromosome), beside the released figure.

Run: python3 bin/archaic_compare.py compare --chen-file F --sample NA12878 \
         --group EUR --chroms "21 22" --mine label=SEGMENTS.tsv [--mine label=...]
"""
import argparse
import collections
import statistics
import sys


def norm(c):
    return c[3:] if c.startswith("chr") else c


def merged_bp(ivs):
    """Total length of the union of (start, end) intervals."""
    total, cur = 0, None
    for s, e in sorted(ivs):
        if cur and s <= cur[1]:
            cur[1] = max(cur[1], e)
        else:
            if cur:
                total += cur[1] - cur[0]
            cur = [s, e]
    return total + (cur[1] - cur[0] if cur else 0)


def read_chen(path):
    """({person: {chrom: [(start, end)]}}, {person: group})."""
    calls = collections.defaultdict(lambda: collections.defaultdict(list))
    group = {}
    with open(path) as f:
        for line in f:
            c = line.split()
            if len(c) < 9 or c[0] == "chr":
                continue
            calls[c[8]][norm(c[0])].append((int(c[1]), int(c[2])))
            group[c[8]] = c[7]
    return calls, group


def read_mine(path):
    """{chrom: [(start, end)]} of the Neanderthal rows of a segments TSV."""
    out = collections.defaultdict(list)
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        if head[:5] != ["chrom", "start", "end", "length_bp", "ancestry"]:
            sys.exit(f"{path}: not a step 43 segments file")
        for line in f:
            c = line.rstrip("\n").split("\t")
            if c[4] == "Neanderthal":
                out[norm(c[0])].append((int(c[1]), int(c[2])))
    return out


def chrom_order(c):
    return (0, int(c)) if c.isdigit() else (1, c)


def cmd_compare(a):
    calls, group = read_chen(a.chen_file)
    if a.sample not in calls:
        sys.exit(f"{a.sample} is not in {a.chen_file}")
    peers = [p for p, g in group.items() if g == a.group]
    ivs = collections.defaultdict(lambda: collections.defaultdict(list))
    labels = []
    for item in a.mine:
        label, _, path = item.partition("=")
        if label not in labels:
            labels.append(label)
        for c, v in read_mine(path).items():
            ivs[label][c].extend(v)
    mine = {l: {c: merged_bp(v) for c, v in ivs[l].items()} for l in labels}
    chroms = sorted({norm(c) for c in a.chroms.replace(",", " ").split()}, key=chrom_order)
    if not chroms:
        sys.exit("--chroms names no chromosome")

    def chen_bp(p, c):
        return merged_bp(calls[p].get(c, []))

    print(f"| chromosome | Chen {a.sample} (bp) | Chen {a.group} median | Chen {a.group} max | "
          + " | ".join(f"{l} (bp)" for l in labels) + " |")
    print("|---|---|---|---|" + "---|" * len(labels))
    tot_me, tot_med, tot_max = 0, 0.0, 0
    tot_mine = dict.fromkeys(labels, 0)
    for c in chroms:
        vals = [chen_bp(p, c) for p in peers]
        me, med, mx = chen_bp(a.sample, c), statistics.median(vals), max(vals)
        tot_me, tot_med, tot_max = tot_me + me, tot_med + med, tot_max + mx
        cells = []
        for l in labels:
            v = mine[l].get(c, 0)
            tot_mine[l] += v
            cells.append(str(v))
        print(f"| {c} | {me} | {med:.0f} | {mx} | " + " | ".join(cells) + " |")
    print(f"| all | {tot_me} | {tot_med:.0f} | {tot_max} | " + " | ".join(str(tot_mine[l]) for l in labels) + " |")
    print()
    # The person-level spread over the same chromosomes.
    sums = sorted(sum(chen_bp(p, c) for c in chroms) for p in peers)
    print(f"Chen et al., {len(peers)} {a.group} people, these {len(chroms)} chromosomes: "
          f"lowest {sums[0]}, median {statistics.median(sums):.0f}, highest {sums[-1]} bp; "
          f"{a.sample} {tot_me} bp.")
    for l in labels:
        ratio = tot_mine[l] / tot_me if tot_me else float("nan")
        print(f"{l}: {tot_mine[l]} bp, {ratio:.2f} times the released {a.sample} figure.")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compare")
    c.add_argument("--chen-file", required=True)
    c.add_argument("--sample", required=True)
    c.add_argument("--group", default="EUR")
    c.add_argument("--chroms", required=True, help="the chromosomes analysed, separated by spaces or commas")
    c.add_argument("--mine", action="append", required=True, metavar="LABEL=SEGMENTS.tsv")
    c.set_defaults(fn=cmd_compare)
    a = p.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
