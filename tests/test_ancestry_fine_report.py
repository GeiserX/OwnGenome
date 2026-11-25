#!/usr/bin/env python3
# og-leg: OG_REQUIRED
"""og ancestry-fine, step 26b's commands, on made-up inputs:

  panel-map    a group of fewer than 10 unrelated people (one Oceanian, three
               San) is no reference and no painting reference, and the
               install log names it as left out
  check-map    panel-map's table passes; one Oceanian put back as reference
               and painting reference, a short row, a flag other than 0 or 1
               and a table with no reference person each stop it
  chroms       a chromosome is painted only with 100 SNVs spanning 5 cM
  report       fractions with the bootstrap range, merged painting rows, the
               painted shares, a PNG of the right size, the left-out group

The Rust unit tests (rust/crates/og-ancestry) and tests/parity check each
command in detail; this test runs them as the pipeline does, through bin/og.
Both pipeline reports of the table are checked by
rust/crates/og-render/tests/ancestry_fine.rs.

Run: OG=bin/og python3 tests/test_ancestry_fine_report.py
"""
import os
import struct
import subprocess
import sys
import tempfile
import zlib

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

OG = os.environ.get("OG", "")
FAILS = []


def expect(desc, got, want):
    if got == want:
        print(f"PASS {desc}")
    else:
        print(f"FAIL {desc}: got {got!r}, want {want!r}")
        FAILS.append(desc)


def og(*args, cwd=None):
    return subprocess.run([OG, "ancestry-fine", *args], capture_output=True, text=True, cwd=cwd)


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


def run_all(tmp):
    # --- panel-map and check-map --------------------------------------------------
    # 10 unrelated Yoruba, 3 unrelated San, one unrelated Papuan and two
    # related ones: in gnomAD's metadata 29 of the 30 Oceanians are related.
    small = ([(f"y{i}", "YRI", "false") for i in range(10)] + [(f"s{i}", "San", "false") for i in range(3)]
             + [("p0", "PapuanHighlands", "false"), ("p1", "Bougainville", "true"), ("p2", "PapuanSepik", "true")])
    write(os.path.join(tmp, "meta.tsv"), "s\tpopulation\tsample_filters.all_samples_related\n"
          + "".join("\t".join(m) + "\n" for m in small))
    write(os.path.join(tmp, "bcf.txt"), "".join(m[0] + "\n" for m in small))
    samples = os.path.join(tmp, "samples.tsv")
    r = og("panel-map", "--meta", os.path.join(tmp, "meta.tsv"), "--bcf-samples", os.path.join(tmp, "bcf.txt"),
           "--paint-per-group", "50", "--out", samples)
    expect("panel-map: exits 0", (r.returncode, r.stderr), (0, ""))
    expect("panel-map: the install log names the groups left out",
           "Left out of the model, fewer than 10 unrelated reference people: "
           "Oceania (continental and regional group), Southern Africa (regional group)" in r.stdout, True)
    rows = [line.split("\t") for line in open(samples).read().splitlines()]
    by = {row[0]: row for row in rows[1:]}
    expect("panel-map: the header", rows[0], ["sample", "population", "continent", "region", "reference", "paint"])
    expect("panel-map: the one unrelated Oceanian keeps its groups but is no reference", by["p0"][2:], ["OCE", "OCE", "0", "0"])
    expect("panel-map: a regional group of 3 (San) is no reference", sorted({by[f"s{i}"][4] for i in range(3)}), ["0"])
    expect("panel-map: the 10 West Africans are references and painting references",
           sorted({tuple(by[f"y{i}"][4:]) for i in range(10)}), [("1", "1")])
    expect("check-map: panel-map's table passes", og("check-map", "--map", samples).returncode, 0)
    one = rows[0:1] + [row[:4] + ["1", "1"] if row[0] == "p0" else row for row in rows[1:]]
    write(os.path.join(tmp, "one.tsv"), "".join("\t".join(row) + "\n" for row in one))
    r = og("check-map", "--map", os.path.join(tmp, "one.tsv"))
    expect("check-map: stops on one Oceanian put back, and says why",
           (r.returncode, "regional group OCE has 1 of the 10 reference people it needs" in r.stderr,
            "continental group OCE has 1 of the 2 painting references RFMix needs" in r.stderr), (1, True, True))
    header = "sample\tpopulation\tcontinent\tregion\treference\tpaint\n"
    for name, text in (("short", header + "y0\tYRI\tAFR\n"), ("flag", header + "y10\tYRI\tAFR\tWAF\ttrue\t0\n"),
                       ("empty", header)):
        write(os.path.join(tmp, f"{name}.tsv"), text)
        expect(f"check-map: stops on the {name} table", og("check-map", "--map", os.path.join(tmp, f"{name}.tsv")).returncode, 1)

    # --- chroms ---------------------------------------------------------------------
    maps = os.path.join(tmp, "maps")
    os.makedirs(maps)
    for c in ("chr20", "chr22", "chr19"):
        write(os.path.join(maps, f"{c}.map"), f"{c} . 0.0 1000000\n{c} . 10.0 11000000\n")   # 1 cM per Mb
    write(os.path.join(tmp, "counts.tsv"), "chr20\t500\t2000000\t8000000\nchr22\t500\t2000000\t4000000\n"
          "chr21\t999\t1\t40000000\nchr19\t50\t1000000\t11000000\nhomref\t7\n")
    r = og("chroms", "--counts", os.path.join(tmp, "counts.tsv"), "--maps", maps)
    expect("chroms: 6 cM painted; 2 cM, no map and 50 SNVs not", r.stdout, "chr20\n")

    # --- report -----------------------------------------------------------------------
    model = os.path.join(tmp, "model")
    os.makedirs(model)
    write(os.path.join(model, "levels.tsv"), "continent\t3\tEAS,WEA,AFR\nregion\t2\tIBE,ITA\n")
    # The panel's samples.tsv, beside model/: an Oceanian the model left out.
    write(os.path.join(tmp, "samples.tsv"), header + "i1\tIBS\tWEA\tIBE\t1\t1\nt1\tTSI\tWEA\tITA\t1\t1\n"
          "p0\tPapuanHighlands\tOCE\tOCE\t0\t0\n")
    radm = os.path.join(tmp, "admixture")
    os.makedirs(radm)
    write(os.path.join(radm, "S_continent.3.Q"), "0.05 0.9 0.05\n")
    write(os.path.join(radm, "S_continent.boot.Q"), "".join(f"0.05 {0.8 + i / 100:.2f} {0.15 - i / 100:.2f}\n" for i in range(11)))
    write(os.path.join(radm, "S_region.2.Q"), "0.6 0.4\n")
    paint = os.path.join(tmp, "paint")
    os.makedirs(paint)
    write(os.path.join(paint, "chr20.msp.tsv"),
          "#Subpopulation order/codes: AFR=0\tEAS=1\tWEA=2\n#chm\tspos\tepos\tsgpos\tegpos\tn snps\tS.0\tS.1\n"
          "chr20\t1000\t5000\t1\t2\t10\t2\t2\nchr20\t5000\t9000\t2\t3\t10\t2\t0\nchr20\t9000\t11000\t3\t4\t10\t2\t0\n")
    write(os.path.join(paint, "failed.txt"), "chr19\n")   # RFMix failed on chr19; chr20 was still painted
    rep = os.path.join(tmp, "out")
    r = og("report", "--sample", "S", "--model", model, "--admixture", radm, "--paint", paint, "--panel", "test",
           "--out", rep)
    expect("report: exits 0", (r.returncode, r.stderr), (0, ""))
    table = [line.rstrip("\n").split("\t") for line in open(os.path.join(rep, "S_ancestry_fine.tsv"))]
    expect("report: the header", table[0], ["level", "code", "group", "fraction", "se", "low", "high"])
    cont = [row for row in table if row[0] == "continent"]
    expect("report: continental groups, largest first", [row[1] for row in cont], ["WEA", "EAS", "AFR"])
    expect("report: West Eurasia's fraction and bootstrap range", cont[0][3:], ["0.9000", "0.0332", "0.8025", "0.8975"])
    expect("report: a level without bootstrap replicates has no range",
           [row[4:] for row in table if row[0] == "region"][0], [".", ".", "."])
    bed = [line.rstrip("\n").split("\t") for line in open(os.path.join(rep, "S_painting.bed"))][1:]
    expect("report: haplotype 1 one West Eurasian row, haplotype 2 split at 5000",
           bed, [["chr20", "999", "11000", "1", "WEA"], ["chr20", "999", "5000", "2", "WEA"],
                 ["chr20", "4999", "11000", "2", "AFR"]])
    expect("report: painted shares by length", {row[1]: row[3] for row in table if row[0] == "painting"},
           {"WEA": "0.7000", "AFR": "0.3000"})
    png = open(os.path.join(rep, "S_painting.png"), "rb").read()
    w, h = struct.unpack(">II", png[16:24])
    expect("report: a PNG of 1000 pixels by the 22 chromosomes", (png[:8], w, h), (b"\x89PNG\r\n\x1a\n", 1000, 6 + 22 * 19))
    idat = png[png.index(b"IDAT") + 4:png.index(b"IEND") - 8]
    expect("report: the PNG's pixel data is whole", len(zlib.decompress(idat)), h * (1 + 3 * w))
    text = open(os.path.join(rep, "S_ancestry_fine.txt")).read()
    expect("report: stdout is the text", r.stdout, text)
    expect("report: the text names the group, the fraction and its range",
           "West Eurasia (Europe, Middle East, North Africa)     90%  (80% to 90%)" in text, True)
    expect("report: the text says a name is ours, not a nationality", "not a nationality" in text, True)
    expect("report: the text says which group the model left out and why",
           "Left out of the model: Oceania (continental and regional group). The panel has fewer than 10 unrelated" in text, True)

    os.remove(os.path.join(paint, "failed.txt"))
    os.remove(os.path.join(paint, "chr20.msp.tsv"))
    r = og("report", "--sample", "S", "--model", model, "--admixture", radm, "--paint", paint, "--panel", "test",
           "--out", rep)
    expect("report: without a painting the text says none", "Chromosome painting: none." in r.stdout, True)


def main():
    if not OG:
        print("skip: OG is not set (OG=bin/og)")
        return 1 if os.environ.get("OG_REQUIRED") else 0
    with tempfile.TemporaryDirectory() as tmp:
        run_all(tmp)
    print("")
    if FAILS:
        print(f"RESULT: {len(FAILS)} FAILED")
        return 1
    print("RESULT: ALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
