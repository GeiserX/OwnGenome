#!/usr/bin/env python3
"""Compare step 43's Neanderthal totals with Chen et al. 2020's released calls.

Used by the archaic jobs of public-evidence.yml only (not part of the
pipeline). Standard library only.

The released calls are Chen et al.'s file "Neanderthal sequence in 1000
genome.50kb.txt" (chr, start, end, LOD, maxLOD, size, pop, anc, ID, tab
separated), merged per person and chromosome.

  compare  the same person run through step 43 on one or more inputs (a
           segments TSV), beside the released figure, and how much of each
           lies under the released calls.
  genomes  the raw IBDmix output of one step 43 run, per archaic genome
           (Altai, Vindija, Denisova) at LOD 4 and 50 kb, and the union, each
           beside the released calls.

Run: python3 bin/archaic_compare.py compare --chen-file F --sample NA12878 \
         --group EUR --chroms "21 22" --mine label=SEGMENTS.tsv [--mine label=...]
"""
import argparse
import collections
import os
import statistics
import sys

ARCHAICS = ("Altai", "Vindija", "Denisova")


def ibdmix_segments(path, min_lod=4.0, min_length=50000):
    """IBDmix's segments for SAMPLE in one of its output files, kept as step 43
    keeps them (og archaic report: LOD 4 and 50 kb or longer), as (start, end)."""
    out = []
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        if head[:5] != ["ID", "chrom", "start", "end", "slod"]:
            sys.exit(f"ERROR: {path} is not IBDmix output (header {' '.join(head[:5])})")
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 5 or c[0] != "SAMPLE":
                continue
            s, e, lod = int(c[2]), int(c[3]), float(c[4])
            if lod >= min_lod and e - s >= min_length:
                out.append((s, e))
    return out


def norm(c):
    """The chromosome name without a leading chr."""
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


def intersect_bp(a, b):
    """Base pairs under both interval lists (each merged first)."""
    def merged(ivs):
        """The intervals with overlaps joined."""
        out = []
        for s, e in sorted(ivs):
            if out and s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        return out
    a, b, i, j, total = merged(a), merged(b), 0, 0, 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi > lo:
            total += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total


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
    """Sort key: numbered chromosomes first, in order."""
    return (0, int(c)) if c.isdigit() else (1, c)


def cmd_compare(a):
    """Print the per-chromosome table and the totals of the compare command."""
    calls, group = read_chen(a.chen_file)
    if a.sample not in calls:
        sys.exit(f"{a.sample} is not in {a.chen_file}")
    peers = [p for p, g in group.items() if g == a.group]
    if not peers:
        sys.exit(f"no person of group {a.group} in {a.chen_file}")
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
        """One person's merged base pairs on one chromosome."""
        return merged_bp(calls[p].get(c, []))

    print(f"| chromosome | Chen {a.sample} (bp) | Chen {a.group} median | Chen {a.group} max | "
          + " | ".join(f"{l} (bp)" for l in labels) + " |")
    print("|---|---|---|---|" + "---|" * len(labels))
    tot_me = 0
    # the "all" median and maximum come from person-level totals, not from sums of
    # per-chromosome medians or maxima, which no person need have
    sums = sorted(sum(chen_bp(p, c) for c in chroms) for p in peers)
    tot_med, tot_max = statistics.median(sums), sums[-1]
    tot_mine = dict.fromkeys(labels, 0)
    for c in chroms:
        vals = [chen_bp(p, c) for p in peers]
        me, med, mx = chen_bp(a.sample, c), statistics.median(vals), max(vals)
        tot_me += me
        cells = []
        for l in labels:
            v = mine[l].get(c, 0)
            tot_mine[l] += v
            cells.append(str(v))
        print(f"| {c} | {me} | {med:.0f} | {mx} | " + " | ".join(cells) + " |")
    print(f"| all | {tot_me} | {tot_med:.0f} | {tot_max} | " + " | ".join(str(tot_mine[l]) for l in labels) + " |")
    print()
    # The person-level spread over the same chromosomes.
    print(f"Chen et al., {len(peers)} {a.group} people, these {len(chroms)} chromosomes: "
          f"lowest {sums[0]}, median {statistics.median(sums):.0f}, highest {sums[-1]} bp; "
          f"{a.sample} {tot_me} bp.")
    for l in labels:
        ratio = tot_mine[l] / tot_me if tot_me else float("nan")
        both = sum(intersect_bp(ivs[l].get(c, []), calls[a.sample].get(c, [])) for c in chroms)
        print(f"{l}: {tot_mine[l]} bp, {ratio:.2f} times the released {a.sample} figure; "
              f"{both} bp under the released calls, {tot_mine[l] - both} bp only here, {tot_me - both} bp only released.")


def cmd_genomes(a):
    """Print the per-archaic-genome table of the genomes command."""
    calls, _ = read_chen(a.chen_file)
    if a.sample not in calls:
        sys.exit(f"{a.sample} is not in {a.chen_file}")
    chroms = []
    with open(a.spans) as f:
        next(f)
        chroms = list(dict.fromkeys(norm(line.split("\t")[0]) for line in f if line.strip()))
    per = collections.defaultdict(lambda: collections.defaultdict(list))
    for c in chroms:
        for n in ARCHAICS:
            per[n][c] += ibdmix_segments(os.path.join(a.ibdmix_dir, f"chr{c}.{n}.txt"))
    per["Altai+Vindija"] = {c: per["Altai"][c] + per["Vindija"][c] for c in chroms}
    released = sum(merged_bp(calls[a.sample].get(c, [])) for c in chroms)
    print(f"| IBDmix output | bp | under the released {a.sample} calls | only here | only released |")
    print("|---|---|---|---|---|")
    for n in ("Altai", "Vindija", "Denisova", "Altai+Vindija"):
        mine = sum(merged_bp(per[n][c]) for c in chroms)
        both = sum(intersect_bp(per[n][c], calls[a.sample].get(c, [])) for c in chroms)
        print(f"| {n} | {mine} | {both} | {mine - both} | {released - both} |")
    print(f"\nReleased {a.sample}: {released} bp on these {len(chroms)} chromosomes.")


def main(argv=None):
    """Parse the command line and run the chosen command."""
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compare")
    c.add_argument("--chen-file", required=True)
    c.add_argument("--sample", required=True)
    c.add_argument("--group", default="EUR")
    c.add_argument("--chroms", required=True, help="the chromosomes analysed, separated by spaces or commas")
    c.add_argument("--mine", action="append", required=True, metavar="LABEL=SEGMENTS.tsv")
    c.set_defaults(fn=cmd_compare)
    g = sub.add_parser("genomes")
    g.add_argument("--chen-file", required=True)
    g.add_argument("--sample", required=True)
    g.add_argument("--ibdmix-dir", required=True)
    g.add_argument("--spans", required=True)
    g.set_defaults(fn=cmd_genomes)
    a = p.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
