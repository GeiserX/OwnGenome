#!/usr/bin/env python3
"""og ancestry-fine, step 26b's commands, on made-up inputs, and what is left
of bin/ancestry_fine_report.py:

  panel-map    a group of fewer than 10 unrelated people (one Oceanian, three
               San) is no reference and no painting reference, and the
               install log names it as left out
  check-map    panel-map's table passes; one Oceanian put back as reference
               and painting reference, a short row, a flag other than 0 or 1
               and a table with no reference person each stop it
  order        the column of each group in a supervised ADMIXTURE .Q
  sample-bed   the .bed code of each genotype (ALT is the .bim's A1), chr22 in
               the VCF meets 22 in the .bim, the .P.in and .ped rows of the
               kept SNVs, and the bootstrap draws only blocks that hold SNVs
  chroms       a chromosome is painted only with 100 SNVs spanning 5 cM
  report       fractions with the bootstrap range, merged painting rows, the
               painted shares, a PNG of the right size, the left-out group
               and the chromosome RFMix failed on, in the table, the text and
               both pipeline reports (collect_summary.py, render_report.py)
  the module   the group names, colours and MIN_REFERENCES collect_summary.py
               reads from it are the ones og writes

The Rust unit tests (rust/crates/og-ancestry) and tests/parity check each
command in detail; this test runs them as the pipeline does, through bin/og.

Run: OG=bin/og python3 tests/test_ancestry_fine_report.py
"""
import os
import struct
import subprocess
import sys
import tempfile
import zlib

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import ancestry_fine_report as af  # noqa: E402

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
           f"Left out of the model, fewer than {af.MIN_REFERENCES} unrelated reference people: "
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

    # --- order ----------------------------------------------------------------------
    write(os.path.join(tmp, "ref.pop"), "EUR\nEUR\nAFR\nEAS\nAFR\n")
    write(os.path.join(tmp, "ref.3.Q"), "1 0 0\n0.99 0.01 0\n0 0 1\n0 1 0\n0 0.02 0.98\n")
    r = og("order", "--pop", os.path.join(tmp, "ref.pop"), "--q", os.path.join(tmp, "ref.3.Q"))
    expect("order: each column is the group whose people fill it", r.stdout, "EUR,EAS,AFR\n")

    # --- sample-bed -----------------------------------------------------------------
    write(os.path.join(tmp, "ref.bim"), "".join(f"22\t22:{p}\t0\t{p}\t{alt}\t{ref}\n" for p, ref, alt in
                                                [(100, "A", "G"), (200, "C", "T"), (300, "G", "A"), (400, "T", "C"), (500, "A", "C")]))
    write(os.path.join(tmp, "g.vcf"), "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
          "chr22\t100\t.\tA\tG\t.\t.\t.\tGT\t1/1\n"     # hom ALT: code 0
          "chr22\t200\t.\tC\tT\t.\t.\t.\tGT\t0|1\n"     # het: code 2
          "chr22\t300\t.\tG\tA\t.\t.\t.\tGT\t0/0\n"     # hom REF: code 3
          "chr22\t400\t.\tT\tG\t.\t.\t.\tGT\t0/1\n")    # another ALT: left out; 500 missing: left out
    write(os.path.join(tmp, "continent.P"), "".join(f"0.{i} 0.{9 - i}\n" for i in range(5)))
    adm = os.path.join(tmp, "adm")
    r = og("sample-bed", "--bim", os.path.join(tmp, "ref.bim"), "--vcf", os.path.join(tmp, "g.vcf"), "--sample", "S",
           "--level", f"continent:2:{os.path.join(tmp, 'continent.P')}", "--bootstrap", "7", "--min-sites", "3",
           "--out-dir", adm)
    expect("sample-bed: exits 0 and counts the genotyped SNVs", (r.returncode, r.stdout),
           (0, "3 of the 5 SNVs ADMIXTURE uses genotyped\n"))
    expect("sample-bed: .bed magic, then hom ALT 0, het 2, hom REF 3",
           open(os.path.join(adm, "S_continent.bed"), "rb").read(), b"\x6c\x1b\x01\x00\x02\x03")
    expect("sample-bed: the .P.in rows of the kept SNVs only",
           open(os.path.join(adm, "S_continent.2.P.in")).read(), "0.0 0.9\n0.1 0.8\n0.2 0.7\n")
    expect("sample-bed: .ped allele pairs (2 is A1): hom ALT, het, hom REF",
           open(os.path.join(adm, "S.geno")).read(), "2 2\n1 2\n1 1\n")
    sh = open(os.path.join(adm, "run_admixture.sh")).read()
    expect("sample-bed: run_admixture.sh runs 7 replicates of sample S", ("-le 7 ]" in sh, "S_${level}.bed" in sh), (True, True))
    # The block bootstrap draws only blocks that hold SNVs: 4 SNVs in 3 blocks
    # round up to blocks of 2, so there are 2 blocks and every replicate holds
    # 3 x 2 SNVs.
    awk_prog = sh.split("-v k=\"$k\" '", 1)[1].split("' \"S.geno\"", 1)[0]
    boot = os.path.join(tmp, "boot")
    os.makedirs(boot)
    write(os.path.join(boot, "g"), "".join(f"1 {i}\n" for i in range(4)))
    write(os.path.join(boot, "m"), "".join(f"22 s{i} 0 {i}\n" for i in range(4)))
    write(os.path.join(boot, "p"), "".join(f"0.{i} 0.5\n" for i in range(4)))
    sizes = set()
    for seed in range(1, 21):
        subprocess.run(["awk", "-v", f"r={seed}", "-v", "nb=3", "-v", "k=2", awk_prog, "g", "m", "p"], cwd=boot, check=True)
        sizes.add(len(open(os.path.join(boot, "rep.map")).readlines()))
    expect("sample-bed: every bootstrap replicate draws 3 non-empty blocks of 2 SNVs (20 seeds)", sizes, {6})

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
    expect("report: the table names the group left out and the chromosome not painted",
           [row[:3] for row in table if row[0] in ("left_out", "not_painted")],
           [["left_out", "OCE", "Oceania (continental group)"], ["left_out", "OCE", "Oceania (regional group)"],
            ["not_painted", "chr19", "RFMix failed"]])
    expect("report: the text says which group the model left out and why",
           "Left out of the model: Oceania (continental and regional group). The panel has fewer than 10 unrelated" in text, True)
    expect("report: the text says chr19 was not painted, and chr20 still was",
           ("Not painted: RFMix failed on chr19" in text, "of chr20:" in text), (True, True))

    # --- the module and og agree on the names collect_summary.py reads ---------------
    expect("module: each group name in og's table is the module's label",
           [row[2] for row in table[1:] if row[0] in ("continent", "region", "painting")],
           [af.label("region" if row[0] == "region" else "continent", row[1])
            for row in table[1:] if row[0] in ("continent", "region", "painting")])
    expect("module: MIN_REFERENCES is 10, as og's", af.MIN_REFERENCES, 10)
    expect("module: a colour for each continental group", sorted(af.COLOURS), sorted(af.LABELS["continent"]))

    # Both pipeline reports (collect_summary.py, then render_report.py) carry them.
    import collect_summary  # noqa: E402
    import render_report  # noqa: E402
    sd = os.path.join(tmp, "S")
    os.makedirs(os.path.join(sd, "ancestry_fine"))
    for name in ("S_ancestry_fine.tsv", "S_painting.png"):
        with open(os.path.join(rep, name), "rb") as src, open(os.path.join(sd, "ancestry_fine", name), "wb") as dst:
            dst.write(src.read())
    summ = collect_summary.collect("S", sd)
    txt, html = render_report.text_report(summ), render_report.html_report(summ)
    for kind, doc in (("text", txt), ("HTML", html)):
        expect(f"report: the {kind} report names the group left out",
               "Left out of the model: Oceania (continental and regional group). The panel has fewer than 10" in doc, True)
        expect(f"report: the {kind} report names the chromosome not painted", "not painted: RFMix failed on chr19" in doc, True)
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
