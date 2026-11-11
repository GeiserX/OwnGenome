"""ancestry_fine_report.py: what other Python still reads of step 26b.

Step 26b's commands (panel-map, check-map, order, sample-bed, chroms and
report) run in the Rust binary as og ancestry-fine (rust/crates/og-ancestry).
This module keeps only the names two Python programs import:

  bin/collect_summary.py          label, COLOURS and MIN_REFERENCES, for the
                                  report's group names and painting legend
  tests/demo/make_demo_sample.py  label, CHROM_LEN, painting_shares,
                                  write_rows and painting_png, for the docs
                                  demo, which is built without og

The group names, the colours and MIN_REFERENCES must stay the same as
rust/crates/og-ancestry/src/lib.rs; tests/test_ancestry_fine_report.py checks
them against what og writes.

A group name is ours: a set of reference populations we put together and
named, not a nationality. docs/26b-ancestry-fine.md lists them.
"""
import os
import struct
import zlib

LABELS = {
    "continent": {
        "AFR": "Africa",
        "WEA": "West Eurasia (Europe, Middle East, North Africa)",
        "CSA": "Central and South Asia",
        "EAS": "East Asia",
        "AMR": "Americas (Native American)",
        "OCE": "Oceania",
    },
    "region": {
        "WAF": "West Africa", "EAF": "East Africa", "CAF": "Central African rainforest", "SAF": "Southern Africa",
        "NWE": "Northwestern Europe", "NEE": "Northeastern Europe", "IBE": "Iberia", "ITA": "Italy",
        "SAR": "Sardinia", "CAU": "Caucasus", "LEV": "Levant", "BDN": "Bedouin", "NAF": "North Africa",
        "CSA": "Central and South Asia", "EAS": "East Asia", "SIB": "Siberia",
        "AMR": "Americas (Native American)", "OCE": "Oceania",
    },
}

# Painting colours of the continental groups, and grey where nothing is painted.
COLOURS = {"AFR": (230, 159, 0), "WEA": (0, 114, 178), "CSA": (0, 158, 115), "EAS": (204, 121, 167),
           "AMR": (213, 94, 0), "OCE": (86, 180, 233)}
UNPAINTED = (225, 225, 225)

# GRCh38 lengths of chr1 to chr22, for the painting picture.
CHROM_LEN = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
             138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
             83257441, 80373285, 58617616, 64444167, 46709983, 50818468]

# The fewest unrelated reference people a group needs to be modelled
# (rust/crates/og-ancestry/src/lib.rs says why).
MIN_REFERENCES = 10


def label(level, code):
    return LABELS.get(level, {}).get(code, code)


def write_rows(path, header, rows):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    os.replace(tmp, path)


def chrom_key(c):
    c = c[3:] if c.startswith("chr") else c
    return (0, int(c)) if c.isdigit() else (1, c)


def painting_shares(rows):
    total = sum(r[2] - r[1] for r in rows)
    by = {}
    for r in rows:
        by[r[4]] = by.get(r[4], 0) + r[2] - r[1]
    return total, {k: v / total for k, v in by.items()} if total else {}


def png_bytes(width, height, pixels):
    """A plain RGB PNG of rows of (r, g, b) tuples."""
    raw = b"".join(b"\x00" + bytes(v for px in row for v in px) for row in pixels)

    def chunk(kind, data):
        return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data) & 0xffffffff)

    return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b""))


def painting_png(rows, width=1000, bar=6, gap=6):
    """Chromosomes 1 to 22 as two bars each (the two haplotypes), drawn to
    scale, each painted segment in its group's colour, the rest grey."""
    scale = (width - 20) / CHROM_LEN[0]
    height = gap + 22 * (2 * bar + 1 + gap)
    img = [[(255, 255, 255)] * width for _ in range(height)]
    for i, n in enumerate(CHROM_LEN):
        y0 = gap + i * (2 * bar + 1 + gap)
        for hap in (0, 1):
            y = y0 + hap * (bar + 1)
            for yy in range(y, y + bar):
                for x in range(10, 10 + int(n * scale)):
                    img[yy][x] = UNPAINTED
    for chrom, start, end, hap, code in rows:
        key = chrom_key(chrom)
        if key[0] != 0 or not 1 <= key[1] <= 22:
            continue
        i = key[1] - 1
        y = gap + i * (2 * bar + 1 + gap) + (hap - 1) * (bar + 1)
        x0 = 10 + int(start * scale)
        x1 = max(x0 + 1, 10 + int(end * scale))
        colour = COLOURS.get(code, (0, 0, 0))
        for yy in range(y, y + bar):
            for x in range(x0, min(x1, width - 10)):
                img[yy][x] = colour
    return png_bytes(width, height, img)
