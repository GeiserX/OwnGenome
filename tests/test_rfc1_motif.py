#!/usr/bin/env python3
"""og str motif: which motif do the reads at RFC1 carry? By meaning.

This file holds the made-up read sets at the GRCh38 RFC1 repeat
(chr4:39348425-39348479, the reference (AAAAG)11) with random flanks, which
the e2e cases (expansions-and-rfc1-1-motif, str-catalog-1-strchive,
nextflow-report-parity-1-sections) and rust/crates/og-bam/tests/make_parity.py
import, and a test of what og must write on them:
  benign      reads across the repeat with an expanded AAAAG allele, and
              reads whose AAAAG part is soft-clipped: benign_motif
  pathogenic  AAGGG reads of four kinds: soft-clipped next to the repeat,
              unmapped and placed at their mate's position, mates mapped to
              another contig (mate rescue), and reverse-strand CCCTT reads
              across the repeat; three reads of an expanded AAAAG allele;
              and eight reads across an unexpanded (AAAAG)11 allele, which
              are not counted: pathogenic_motif, with AAGGG and AAAAG shares
  mixed       AAAAG and AAAGG reads only: benign_motif, both shares listed
  acagg       ACAGG reads: pathogenic_motif
  unknown     AAGAG reads (not a known benign or pathogenic motif):
              uncertain_motif
  few         two repeat reads: not_resolved
and duplicates, secondary and QC-failed alignments, a homopolymer and a run
of three units are not counted.

`cargo test -p og-bam` (tests/parity.rs) runs og on BAMs of these reads and
holds it to what the Python tool it replaced wrote there, byte for byte
(rust/crates/og-bam/tests/data/parity/want/motif-*). This test says what
those files must mean, so a golden rewritten with another verdict fails here
even when og agrees with it. Controls, each a check that fails when the
rescue it names is off: the pathogenic case counts its 6 soft-clipped reads
(--no-softclip counts 6 fewer) and its 3 rescued mates (--no-mates counts 3
fewer).

Run: python3 tests/test_rfc1_motif.py
"""
import os
import random
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARITY = os.path.join(REPO, "rust", "crates", "og-bam", "tests", "data", "parity")

FAILS = []
CHROM, START, END = "chr4", 39348424, 39348479  # 0-based start, end exclusive: 55 bases, (AAAAG)11
OTHER = "chr1"
RNG = random.Random(7)
LEFT = "".join(RNG.choice("ACGT") for _ in range(300))
RIGHT = "".join(RNG.choice("ACGT") for _ in range(300))
VCF = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
       f"{CHROM}\t{START + 1}\t.\tA\t<STR100>\t.\tPASS\tEND={END};REF=11;RL=55;RU=AARRG;VARID=RFC1;REPID=RFC1;"
       "STR_STATUS=full_mutation\tGT:SO:REPCN\t" + ":".join(["/".join("11"), "INREPEAT/INREPEAT", "/".join(["100", "150"])]) + "\n")


def expect(desc, got, want):
    if got == want:
        print(f"PASS {desc}")
    else:
        print(f"FAIL {desc}: got {got!r}, want {want!r}")
        FAILS.append(desc)


def expect_true(desc, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'} {desc}" + ("" if ok else f": {detail}"))
    if not ok:
        FAILS.append(desc)


class Read:
    """A read: the attributes a SAM or BAM writer needs."""

    def __init__(self, name, contig, start, cigar, seq, read1=True, unmapped=False, mate=None, reverse=False,
                 dup=False, secondary=False, qcfail=False):
        self.query_name, self.contig, self.reference_start = name, contig, start
        self.reference_id = 0 if contig == CHROM else 1
        self.cigar_text = cigar
        self.query_sequence = seq
        self.is_read1, self.is_read2 = read1, not read1
        self.is_unmapped, self.is_reverse = unmapped, reverse
        self.is_duplicate, self.is_secondary, self.is_supplementary, self.is_qcfail = dup, secondary, False, qcfail
        self.is_paired = True
        mate = mate or (contig, start)
        self.mate_contig, self.next_reference_start = mate
        self.next_reference_id = 0 if mate[0] == CHROM else 1
        self.mate_is_unmapped = False


def spanning(name, unit, units, offset=40, read1=True, **kw):
    """A 150-base read aligned across the repeat: offset flank bases, the
    repeat (55 reference bases plus an insertion of the rest), flank."""
    rep = unit * units
    right = 150 - offset - len(rep)
    assert right > 0 and len(rep) >= END - START
    ins = len(rep) - (END - START)
    cig = f"{offset}M{END - START}M" + (f"{ins}I" if ins else "") + f"{right}M"
    return Read(name, CHROM, START - offset, cig, LEFT[-offset:] + rep + RIGHT[:right], read1=read1, **kw)


def clipped(name, unit, read1=True):
    """A read whose 60 aligned bases end at the repeat start; the 90 repeat bases are soft-clipped."""
    seq = LEFT[-60:] + (unit * 18)
    return Read(name, CHROM, START - 60, "60M90S", seq, read1=read1)


def placed_unmapped(name, unit):
    """An unmapped read placed at its mapped mate's position, all repeat."""
    return Read(name, CHROM, START - 100, "", (unit * 30), read1=False, unmapped=True)


def rescued(name, unit):
    """A flank read at the locus and its repeat mate mapped to another contig."""
    a = Read(name, CHROM, START - 120, "100M", LEFT[-120:-20], read1=True, mate=(OTHER, 5000))
    b = Read(name, OTHER, 5000, "150M", unit * 30, read1=False, mate=(CHROM, START - 120))
    return a, b


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def case(kind):
    """(reads at the locus, mates elsewhere by name) of one case."""
    at, elsewhere = [], {}
    if kind == "benign":
        at += [spanning(f"b{i}", "AAAAG", 20) for i in range(6)]
        at += [clipped(f"bc{i}", "AAAAG") for i in range(4)]
    elif kind == "pathogenic":
        at += [clipped(f"pc{i}", "AAGGG") for i in range(6)]
        at += [placed_unmapped(f"pu{i}", "AAGGG") for i in range(2)]
        for i in range(3):
            a, b = rescued(f"pm{i}", "AAGGG")
            at.append(a)
            elsewhere[a.query_name] = b
        at += [Read(f"pr{i}", CHROM, START - 40, "150M", revcomp(LEFT[-40:] + "AAGGG" * 22), reverse=True)
               for i in range(2)]
        at += [spanning(f"aaaag{i}", "AAAAG", 20) for i in range(3)]
        # Not counted: reads across an unexpanded (AAAAG)11 allele, and one
        # across 12 AAGGG units would be (more than the reference's 11).
        at += [spanning(f"ref{i}", "AAAAG", 11) for i in range(8)]
        # Not counted: a duplicate, a secondary and a QC-failed alignment.
        at.append(clipped("pc0", "AAGGG", read1=True))  # the same read again
        at.append(Read("dup", CHROM, START - 60, "60M90S", LEFT[-60:] + "AAGGG" * 18, dup=True))
        at.append(Read("sec", CHROM, START - 60, "60M90S", LEFT[-60:] + "AAGGG" * 18, secondary=True))
        at.append(Read("qc", CHROM, START - 60, "60M90S", LEFT[-60:] + "AAGGG" * 18, qcfail=True))
        # Not counted either: a homopolymer and a run of three units.
        at.append(Read("poly", CHROM, START - 60, "60M90S", LEFT[-60:] + "A" * 90))
        at.append(Read("short", CHROM, START - 60, "60M90S", LEFT[-60:] + "AAGGG" * 3 + RIGHT[:75]))
    elif kind == "mixed":
        at += [spanning(f"m{i}", "AAAAG", 20) for i in range(6)]
        at += [spanning(f"n{i}", "AAAGG", 20) for i in range(4)]
    elif kind == "acagg":
        at += [clipped(f"a{i}", "ACAGG") for i in range(5)]
    elif kind == "unknown":
        at += [clipped(f"u{i}", "AAGAG") for i in range(5)]
    elif kind == "few":
        at += [clipped(f"f{i}", "AAGGG") for i in range(2)]
    return at, elsewhere


def write_bam(pysam, path, kind):
    at, elsewhere = case(kind)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CHROM, "LN": 190214555}, {"SN": OTHER, "LN": 248956422}]}
    unsorted = path + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for r in at + list(elsewhere.values()):
            a = pysam.AlignedSegment()
            a.query_name = r.query_name
            a.query_sequence = r.query_sequence
            a.query_qualities = pysam.qualitystring_to_array("I" * len(r.query_sequence))
            a.flag = (1 | (64 if r.is_read1 else 128) | (4 if r.is_unmapped else 0) | (16 if r.is_reverse else 0)
                      | (256 if r.is_secondary else 0) | (512 if r.is_qcfail else 0) | (1024 if r.is_duplicate else 0))
            a.reference_id = r.reference_id
            a.reference_start = r.reference_start
            a.next_reference_id = r.next_reference_id
            a.next_reference_start = r.next_reference_start
            if not r.is_unmapped:
                a.cigarstring = r.cigar_text
                a.mapping_quality = 60
            out.write(a)
    pysam.sort("-o", path, unsorted)
    pysam.index(path)
    os.remove(unsorted)


def golden(case):
    """og's golden of a parity case when it has one (a deliberate fix), else the oracle's."""
    d = os.path.join(PARITY, "og", case)
    return d if os.path.isdir(d) else os.path.join(PARITY, "want", case)


def table(case):
    path = os.path.join(golden(case), "files", "motif.tsv")
    lines = [line.rstrip("\n").split("\t") for line in open(path)]
    return dict(zip(lines[0], lines[1]))


def exit_code(case):
    return int(open(os.path.join(golden(case), "exit")).read().strip())


def main():
    for kind, verdict, reads in (("benign", "benign_motif", "10"), ("pathogenic", "pathogenic_motif", "16"),
                                 ("mixed", "benign_motif", "10"), ("acagg", "pathogenic_motif", "5"),
                                 ("unknown", "uncertain_motif", "5"), ("few", "not_resolved", "2")):
        t = table(f"motif-{kind}")
        expect(f"{kind}: exit 0, {verdict}, {reads} repeat reads",
               (exit_code(f"motif-{kind}"), t["verdict"], t["repeat_reads"]), (0, verdict, reads))
    t = table("motif-pathogenic")
    expect("pathogenic: the locus from the VCF, AAGGG dominant, both shares; the 8 reference-allele reads not counted",
           (t["chrom"], t["start"], t["end"], t["dominant_motif"], t["motif_shares"]),
           (CHROM, str(START + 1), str(END), "AAGGG", "AAGGG=0.81,AAAAG=0.19"))
    expect_true("pathogenic: the note asks for repeat-primed PCR", "repeat-primed PCR" in t["note"], t["note"])
    expect_true("benign: the note says the flag is not supported",
                "the flag is not supported" in table("motif-benign")["note"], table("motif-benign")["note"])
    expect("mixed: AAAAG and AAAGG counted apart", table("motif-mixed")["motif_shares"], "AAAAG=0.60,AAAGG=0.40")
    expect("unknown: AAGAG named by its smallest form, class unknown",
           (table("motif-unknown")["dominant_motif"], table("motif-unknown")["dominant_class"]), ("AAGAG", "unknown"))
    expect("control: --no-softclip loses the 6 clipped reads", table("motif-pathogenic-no-softclip")["repeat_reads"], "10")
    expect("control: --no-mates loses the 3 rescued mates", table("motif-pathogenic-no-mates")["repeat_reads"], "13")
    t = table("motif-skip")
    expect("--skip: not_checked with the reason", (t["verdict"], t["note"]), ("not_checked", "no BAM"))
    expect("a VCF without an RFC1 record: not_checked", table("motif-no-rfc1-record")["verdict"], "not_checked")
    if FAILS:
        print(f"{len(FAILS)} check(s) failed")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
