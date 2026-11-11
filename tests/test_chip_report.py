#!/usr/bin/env python3
"""Step 49 (clonal haematopoiesis): og chip on synthetic bcftools output and on
assets/chip_genes.tsv, and what bin/chip_report.py still holds for the reports.

Without og (OG unset) only part A runs:
  A. bin/chip_report.py's detection limit: about 14% at 30x (binomial, 3
     reads, 8 in 10), lower at 60x, none under 3x; the summary reader and
     the report lines show the call and the limit (bin/collect_summary.py
     sec_chip, bin/render_report.py chip_lines), the median limit included.
With OG naming the og binary, also:
  1. og chip check-table passes assets/chip_genes.tsv; a malformed class
     token, a gene listed twice, a gene that counts no class and a bad
     PubMed id each fail it (exit 1);
  2. the CHIP classes: DNMT3A R882H, a TET2 frameshift, JAK2 V617F, a PPM1D
     stop in its last exon and a TP53 missense in the DNA-binding domain are
     chip; DNMT3A M880L, JAK2 V617L, a PPM1D stop in exon 5 and a TP53
     missense outside the domain are other_change, not counted;
  3. the germline exclusion: the same R882H at VAF 0.48 (an inherited
     heterozygote) and 0.97 is likely_germline, at 0.34 it is still chip;
     a non-PASS call is filtered and one with 2 reads is low_support;
  4. what is left out: a synonymous or intronic change, a consequence bcftools
     marks downstream of a stop, and a change on another transcript;
  5. the depth table reads a base with no pileup line as depth 0, and a gene
     off the reference as not_on_reference; the summary sentence says the
     limit, that CHIP starts at 2%, that a negative result does not rule a
     clone out, and that it is not a diagnosis; og's limit is the Python's;
  6. prepare picks the MANE Select transcript of a synthetic Ensembl GFF3 (not
     the other one), renames the contigs (2 -> chr2), pads the coding exons by
     10 bases, leaves out a contig the reference lacks, and stops on a gene
     with no MANE Select transcript; the report command writes both tables,
     and the reports read them.

Run: OG=rust/target/release/og python3 tests/test_chip_report.py
"""
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import chip_report as cr  # noqa: E402
import collect_summary as cs  # noqa: E402
import render_report as rr  # noqa: E402

TABLE = os.path.join(REPO, "assets", "chip_genes.tsv")
OG = os.environ.get("OG", "")
FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS += 1


def write(d, name, text):
    p = os.path.join(d, name)
    with open(p, "w") as fh:
        fh.write(text)
    return p


def og(*args):
    return subprocess.run([OG, "chip", *args], capture_output=True, text=True)


def chip_tables(d, sample, calls, depth):
    """A sample folder with the two tables, as og chip report writes them."""
    os.makedirs(os.path.join(d, sample, "chip"))
    write(os.path.join(d, sample, "chip"), f"{sample}_chip.tsv",
          "\t".join(cr.OUT_COLUMNS) + "\n" + "".join("\t".join(c) + "\n" for c in calls))
    write(os.path.join(d, sample, "chip"), f"{sample}_chip_depth.tsv",
          "\t".join(cr.DEPTH_COLUMNS) + "\n" + "".join("\t".join(r) + "\n" for r in depth))
    return os.path.join(d, sample)


# A. what Python keeps
lim30, lim60 = cr.detection_limit(30), cr.detection_limit(60)
check("the 30x detection limit is between 13% and 16%", 0.13 < lim30 < 0.16, lim30)
check("P(>=3 reads) at that VAF is 0.8", abs(cr.p_at_least(3, 30, lim30) - 0.8) < 1e-6, cr.p_at_least(3, 30, lim30))
check("the 60x limit is lower", lim60 < lim30, (lim60, lim30))
check("no limit below 3x", cr.detection_limit(2.9) is None)
with tempfile.TemporaryDirectory() as d:
    call = ["DNMT3A", "chr2", "1150", "C", "T", "missense", "R882H", "2", "0.130", "6", "46", "PASS", "chip",
            "missense", "hotspot R882"]
    s = chip_tables(d, "S", [call], [["DNMT3A", "T2", "chr2", "392", "30.0", f"{lim30:.3f}", "covered"],
                                     ["JAK2", "T3", "chr9", "0", "", "", "not_on_reference"]])
    _, sec = cs.sec_chip(s, "S")
    check("sec_chip reads the chip call and the median limit",
          [c["protein_change"] for c in sec.get("calls", [])] == ["R882H"] and abs(sec["median_limit"] - 0.137) < 0.001
          and sec["not_covered"] == ["JAK2"], sec)
    text = "\n".join(rr.chip_lines(sec, "S"))
    check("the report lines show the call, the limit and the uncovered gene",
          "DNMT3A R882H: VAF 0.130 (6 of 46 reads)" in text and "about 14% or more" in text and "JAK2" in text, text)
    # With an even count of covered genes the median depth is the mean of the two
    # middle depths, and the reported limit must be the one at that depth.
    e = chip_tables(d, "E", [], [[f"G{x}", "T", "chr1", "10", f"{x}.0", f"{cr.detection_limit(x):.4f}", "covered"]
                                 for x in (20, 30, 60, 80)])
    _, sec4 = cs.sec_chip(e, "E")
    check("with 4 covered genes (20, 30, 60, 80x) the limit is the one at the median 45x",
          sec4.get("median_depth") == 45 and abs(sec4["median_limit"] - cr.detection_limit(45)) < 1e-9, sec4)
    check("no call reads 'No CHIP call.'", "No CHIP call." in rr.chip_lines({"calls": [], "median_limit": None}, "S"))
    check("the table reader reads assets/chip_genes.tsv", "DNMT3A" in [r["gene"] for r in cr.read_table(TABLE)])

if not OG:
    print("skip: OG is not set (OG=rust/target/release/og): only part A ran")
    print(f"\n{'All checks passed' if FAILS == 0 else f'{FAILS} check(s) failed'}")
    sys.exit(1 if FAILS or os.environ.get("OG_REQUIRED") else 0)

# Synthetic MANE transcripts in prepare's GFF3 (exons as start, end, rank).
GENES = {
    "DNMT3A": ("ENST00000321117", "chr2", "-", [(25234000, 25234500, 23), (25240000, 25240300, 22)]),
    "TET2": ("ENST00000380013", "chr4", "+", [(105200000, 105210000, 3)]),
    "JAK2": ("ENST00000381652", "chr9", "+", [(5073700, 5073800, 14)]),
    "PPM1D": ("ENST00000305921", "chr17", "+", [(60600000, 60600200, 5), (60663000, 60666280, 6)]),
    "TP53": ("ENST00000269305", "chr17", "-", [(7673700, 7674000, 7)]),
}
GFF = "##gff-version 3\n" + "".join(
    f"{c}\tensembl\tgene\t{ex[0][0]}\t{ex[-1][1]}\t.\t{st}\t.\tID=gene:{g};Name={g}\n"
    f"{c}\tensembl\tmRNA\t{ex[0][0]}\t{ex[-1][1]}\t.\t{st}\t.\tID=transcript:{tx};Parent=gene:{g};transcript_id={tx}\n"
    + "".join(f"{c}\tensembl\texon\t{s}\t{e}\t.\t{st}\t.\tParent=transcript:{tx};rank={r}\n" for s, e, r in ex)
    + "###\n" for g, (tx, c, st, ex) in GENES.items())

WORK = tempfile.mkdtemp()
GFF_PATH = write(WORK, "chip.gff3", GFF)
BED_PATH = write(WORK, "chip.bed", "chr2\t25233990\t25234510\tDNMT3A\n")
EMPTY = write(WORK, "empty.tsv", "")


def call(chrom, pos, csq, gene, tx, aa, flt="PASS", ad="40,6", af="0.13", dp="46"):
    return [chrom, str(pos), "C", "T", flt, f"{csq}|{gene}|{tx}|protein_coding|-|{aa}|{pos}C>T", ad, af, dp]


def rows_of(calls):
    calls_path = write(WORK, "calls.tsv", "".join("\t".join(c) + "\n" for c in calls))
    out = os.path.join(WORK, "S_chip.tsv")
    p = og("report", "--genes", TABLE, "--gff", GFF_PATH, "--bed", BED_PATH, "--calls", calls_path, "--depth", EMPTY,
           "--sample", "S", "--out", out, "--depth-out", os.path.join(WORK, "S_chip_depth.tsv"))
    if p.returncode:
        return p.stderr
    lines = open(out).read().splitlines()
    return [dict(zip(lines[0].split("\t"), ln.split("\t"))) for ln in lines[1:]]


def status_of(c):
    rows = rows_of([c])
    if isinstance(rows, str):
        return rows
    return (rows[0]["status"], rows[0]["class"], rows[0]["protein_change"]) if rows else None


# 1. the table
p = og("check-table", TABLE)
check("og chip check-table passes assets/chip_genes.tsv", p.returncode == 0 and p.stderr == "", p.stderr)
with tempfile.TemporaryDirectory() as d:
    head = "gene\ttruncating\tmissense\tnote\tpmids\n"
    for desc, body in [
        ("a malformed missense token", "DNMT3A\tall\tR88X2\tn\t25426837\n"),
        ("a malformed truncating spec", "DNMT3A\texons:9-2\tnone\tn\t25426837\n"),
        ("a gene listed twice", "TET2\tall\tnone\tn\t25426837\nTET2\tall\tnone\tn\t25426837\n"),
        ("a gene that counts no class", "TET2\tnone\tnone\tn\t25426837\n"),
        ("a bad PubMed id", "TET2\tall\tnone\tn\tPMID1\n"),
    ]:
        p = og("check-table", write(d, "t.tsv", head + body))
        check(f"check-table exits 1 on {desc}", p.returncode == 1 and "t.tsv: " in p.stderr, (p.returncode, p.stderr))

# 2. classes
D, T = "missense", "ENST00000321117"
check("DNMT3A R882H is chip", status_of(call("chr2", 25234373, D, "DNMT3A", T, "882R>882H")) == ("chip", "missense", "R882H"),
      status_of(call("chr2", 25234373, D, "DNMT3A", T, "882R>882H")))
check("DNMT3A M880L (outside the hotspot) is other_change",
      status_of(call("chr2", 25234380, D, "DNMT3A", T, "880M>880L"))[0] == "other_change")
check("a TET2 frameshift is chip, labelled fs",
      status_of(call("chr4", 105205000, "frameshift", "TET2", "ENST00000380013", "1200QRS>1200P*"))
      == ("chip", "truncating", "Q1200fs"),
      status_of(call("chr4", 105205000, "frameshift", "TET2", "ENST00000380013", "1200QRS>1200P*")))
check("a DNMT3A splice acceptor change is chip (nearest exon)",
      status_of(call("chr2", 25234502, "splice_acceptor", "DNMT3A", T, ""))[0] == "chip")
check("JAK2 V617F is chip", status_of(call("chr9", 5073770, D, "JAK2", "ENST00000381652", "617V>617F"))[0] == "chip")
check("JAK2 V617L is other_change", status_of(call("chr9", 5073770, D, "JAK2", "ENST00000381652", "617V>617L"))[0] == "other_change")
check("a PPM1D stop in its last exon (6) is chip",
      status_of(call("chr17", 60664000, "stop_gained", "PPM1D", "ENST00000305921", "500R>500*"))[0] == "chip")
check("a PPM1D stop in exon 5 is other_change",
      status_of(call("chr17", 60600100, "stop_gained", "PPM1D", "ENST00000305921", "100R>100*"))[0] == "other_change")
check("a TP53 missense at codon 248 (DNA-binding domain) is chip",
      status_of(call("chr17", 7673800, D, "TP53", "ENST00000269305", "248R>248Q"))[0] == "chip")
check("a TP53 missense at codon 340 is other_change",
      status_of(call("chr17", 7673800, D, "TP53", "ENST00000269305", "340E>340K"))[0] == "other_change")

# 3. germline exclusion and support
R882 = dict(chrom="chr2", pos=25234373, csq=D, gene="DNMT3A", tx=T, aa="882R>882H")
check("R882H at VAF 0.48 is likely_germline", status_of(call(**R882, ad="20,19", af="0.48", dp="39"))[0] == "likely_germline")
check("R882H at VAF 0.97 is likely_germline", status_of(call(**R882, ad="1,35", af="0.97", dp="36"))[0] == "likely_germline")
check("R882H at VAF 0.34 is still chip", status_of(call(**R882, ad="25,13", af="0.34", dp="38"))[0] == "chip")
check("a non-PASS R882H is filtered", status_of(call(**R882, flt="weak_evidence"))[0] == "filtered")
check("R882H with 2 reads is low_support", status_of(call(**R882, ad="30,2", af="0.06", dp="32"))[0] == "low_support")
check("the VAF falls back on AD when AF is missing", rows_of([call(**R882, ad="30,10", af=".", dp="40")])[0]["vaf"] == "0.250")

# 4. left out
check("a synonymous change is left out", status_of(call("chr2", 25234373, "synonymous", "DNMT3A", T, "882R")) is None)
check("an intronic change is left out", status_of(call("chr2", 25236000, "intron", "DNMT3A", T, "")) is None)
check("a consequence downstream of a stop is left out", status_of(call("chr2", 25234373, "*missense", "DNMT3A", T, "882R>882H")) is None)
check("a change on another transcript is left out", status_of(call("chr2", 25234373, D, "DNMT3A", "ENST00000264709", "838R>838H")) is None)

# 5. depth and the summary
with tempfile.TemporaryDirectory() as d:
    bed = write(d, "b.bed", "chr2\t100\t110\tDNMT3A\nchr4\t200\t210\tTET2\n")
    dep = write(d, "d.tsv", "".join(f"chr2\t{p}\t30\n" for p in range(101, 106)))   # 5 of 10 bases at 30x
    out, dout = os.path.join(d, "S_chip.tsv"), os.path.join(d, "S_chip_depth.tsv")
    p = og("report", "--genes", TABLE, "--gff", GFF_PATH, "--bed", bed, "--calls", EMPTY, "--depth", dep,
           "--sample", "S", "--out", out, "--depth-out", dout)
    lines = open(dout).read().splitlines() if p.returncode == 0 else []
    rows = {ln.split("\t")[0]: dict(zip(cr.DEPTH_COLUMNS, ln.split("\t"))) for ln in lines[1:]}
    check("a base with no pileup line counts as depth 0 (mean 15)", rows.get("DNMT3A", {}).get("mean_depth") == "15.0",
          rows.get("DNMT3A"))
    check("a gene with no reads is too_shallow", rows.get("TET2", {}).get("status") == "too_shallow"
          and rows["TET2"]["mean_depth"] == "0.0", rows.get("TET2"))
    check("a gene off the reference is not_on_reference", rows.get("JAK2", {}).get("status") == "not_on_reference",
          rows.get("JAK2"))
    check("every table gene has a depth row", sorted(rows) == sorted(r["gene"] for r in cr.read_table(TABLE)))
    check("og's 15x limit is the Python's", rows.get("DNMT3A", {}).get("detection_limit_vaf") == f"{cr.detection_limit(15):.3f}",
          rows.get("DNMT3A"))
    s = p.stdout.splitlines()[-1] if p.stdout else ""
    for want in ("about 26%", "15x", "2%", "does not rule", "Not a diagnosis"):
        check(f"the summary says '{want}'", want in s, s)
    p = og("report", "--genes", TABLE, "--gff", GFF_PATH, "--bed", bed, "--calls", EMPTY, "--depth", EMPTY,
           "--sample", "S", "--out", out, "--depth-out", dout)
    check("the summary of no covered gene says no result either way", "no result either way" in p.stdout, p.stdout)

# 6. prepare and the report command
ENSEMBL = """##gff-version 3
2\tensembl\tgene\t1000\t5000\t.\t-\t.\tID=gene:G1;Name=DNMT3A;biotype=protein_coding
2\tensembl\tmRNA\t1000\t5000\t.\t-\t.\tID=transcript:T1;Parent=gene:G1;biotype=protein_coding;tag=basic;transcript_id=T1
2\tensembl\tCDS\t1100\t1200\t.\t-\t0\tID=CDS:P1;Parent=transcript:T1
2\tensembl\tmRNA\t1000\t5000\t.\t-\t.\tID=transcript:T2;Parent=gene:G1;biotype=protein_coding;tag=basic,MANE_Select;transcript_id=T2
2\tensembl\texon\t1000\t1300\t.\t-\t.\tParent=transcript:T2;rank=2
2\tensembl\tCDS\t1050\t1300\t.\t-\t0\tID=CDS:P2;Parent=transcript:T2
2\tensembl\texon\t4000\t5000\t.\t-\t.\tParent=transcript:T2;rank=1
2\tensembl\tCDS\t4000\t4100\t.\t-\t0\tID=CDS:P2;Parent=transcript:T2
9\tensembl\tgene\t100\t900\t.\t+\t.\tID=gene:G2;Name=JAK2;biotype=protein_coding
9\tensembl\tmRNA\t100\t900\t.\t+\t.\tID=transcript:T3;Parent=gene:G2;biotype=protein_coding;tag=MANE_Select;transcript_id=T3
9\tensembl\texon\t100\t900\t.\t+\t.\tParent=transcript:T3;rank=1
9\tensembl\tCDS\t200\t800\t.\t+\t0\tID=CDS:P3;Parent=transcript:T3
"""
with tempfile.TemporaryDirectory() as d:
    genes = write(d, "genes.tsv", "gene\ttruncating\tmissense\tnote\tpmids\nDNMT3A\tall\tR882\tn\t25426837\n"
                                  "JAK2\tnone\tV617F\tn\t15858187\n")
    gff = write(d, "in.gff3", ENSEMBL)
    fai = write(d, "ref.fai", "chr2\t10000\t6\t60\t61\n")
    out_gff, out_bed = os.path.join(d, "o.gff3"), os.path.join(d, "o.bed")
    p = og("prepare", "--gff", gff, "--genes", genes, "--fai", fai, "--out-gff", out_gff, "--out-bed", out_bed)
    check("prepare runs", p.returncode == 0, p.stderr)
    g = open(out_gff).read() if p.returncode == 0 else ""
    check("prepare keeps the MANE Select transcript and its children", "transcript:T2" in g and "rank=2" in g, g)
    check("prepare drops the other transcript", "T1" not in g, g)
    check("prepare names contigs chr2 and chr9", "\nchr2\t" in g and "\nchr9\t" in g and "\n2\t" not in g, g)
    bed = open(out_bed).read() if p.returncode == 0 else ""
    check("the BED holds the coding exons padded by 10, on the reference's contigs only",
          bed == "chr2\t1039\t1310\tDNMT3A\nchr2\t3989\t4110\tDNMT3A\n", repr(bed))
    bad = write(d, "bad.gff3", ENSEMBL.replace("tag=MANE_Select;", "tag=basic;"))
    p2 = og("prepare", "--gff", bad, "--genes", genes, "--fai", fai, "--out-gff", out_gff, "--out-bed", out_bed)
    check("prepare stops on a gene with no MANE Select transcript", p2.returncode == 1 and "JAK2" in p2.stderr, p2.stderr)

    calls = write(d, "calls.tsv", "chr2\t1150\tC\tT\tPASS\tmissense|DNMT3A|T2|protein_coding|-|882R>882H|1150C>T\t40,6\t0.130\t46\n")
    depth = write(d, "depth.tsv", "".join(f"chr2\t{p}\t30\n" for p in range(1040, 1311)) + "".join(f"chr2\t{p}\t30\n" for p in range(3990, 4111)))
    os.makedirs(os.path.join(d, "S", "chip"))
    out, dout = os.path.join(d, "S", "chip", "S_chip.tsv"), os.path.join(d, "S", "chip", "S_chip_depth.tsv")
    p3 = og("report", "--genes", genes, "--gff", out_gff, "--bed", out_bed, "--calls", calls, "--depth", depth,
            "--sample", "S", "--out", out, "--depth-out", dout)
    check("report runs and prints one CHIP call and the limit", p3.returncode == 0 and "1 CHIP call " in p3.stdout
          and "about 14%" in p3.stdout, p3.stdout + p3.stderr)
    lines = open(out).read().splitlines() if p3.returncode == 0 else []
    check("the calls table has the columns in order and one chip row",
          lines[:1] == ["\t".join(cr.OUT_COLUMNS)] and len(lines) == 2 and "\tchip\t" in lines[1], lines)
    dl = open(dout).read().splitlines() if p3.returncode == 0 else []
    check("the depth table has DNMT3A at 30x and JAK2 not on the reference",
          len(dl) == 3 and "\t30.0\t" in dl[1] and dl[2].endswith("not_on_reference"), dl)
    _, sec = cs.sec_chip(os.path.join(d, "S"), "S")
    check("the reports read og's tables", [c["protein_change"] for c in sec.get("calls", [])] == ["R882H"]
          and sec["not_covered"] == ["JAK2"], sec)
    p4 = og("report", "--genes", genes, "--gff", out_gff, "--bed", out_bed, "--calls", os.path.join(d, "none.tsv"),
            "--depth", depth, "--sample", "S", "--out", out, "--depth-out", dout)
    check("a missing calls file is an ERROR line and exit 1, not a traceback",
          p4.returncode == 1 and p4.stderr.startswith("ERROR: ") and "Traceback" not in p4.stderr, p4.stderr)

print(f"\n{'All checks passed' if FAILS == 0 else f'{FAILS} check(s) failed'}")
sys.exit(1 if FAILS else 0)
