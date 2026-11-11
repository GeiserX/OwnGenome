#!/usr/bin/env bash
# clinvar-cases.sh: writes the small inputs of the og clinvar cases in the
# current directory, the same every time:
#   ref.fa, ref.fa.fai      one 450-base contig per joining scenario (c01 to
#                           c32), each the same two stretches of GRCh38 chr22
#                           (from 30000000 and 30000200, here at 100 and 300)
#                           with N between; lc holds them in lower case
#   join.vcf                the scenarios of tests/test_clinvar_join_adjacent.py
#                           and more, one contig each, header without a PASS
#                           line and with old PGP_ INFO lines
#   join-edge.vcf           lines a joiner must pass or skip: a ## line between
#                           records, a # line, an empty line, short records,
#                           a second #CHROM line, no final line ending
#   filled.vcf              candidates as bcftools +fill-from-fasta fills them,
#                           each kind finish keeps or drops
#   tag.vcf                 a hits VCF for the tag pass
#   hits.vcf, hits-floor.vcf, hits-gene.vcf   hits VCFs for og clinvar hits
#   sites.tsv, sites-vcf.txt, sites-gvcf.txt, empty.txt   the always-report
#                           list and its calls
# scripts/ci/oracle-diff.sh runs it for the clinvar-* cases of
# tests/parity/cases.tsv.
set -euo pipefail
export LC_ALL=C
python3 - <<'PY'
S1 = "AGACTCTGTCGCAAAAAAAAAATTAATGTAACCATTTAGCCATTACCCAGTTCCTTAAAAT"
S2 = "GCCTGTAATCCCAGCACTTTGGGAGGCTGAG"
seq = ["N"] * 450
seq[99:99 + len(S1)] = S1
seq[299:299 + len(S2)] = S2
SEQ = "".join(seq)


def ref_at(pos, n):
    return SEQ[pos - 1:pos - 1 + n]


# Each scenario: (records as (pos, ref, alt, format, sample)), on its own contig.
U = "GT"
P = "GT:PS"
SCEN = [
    [(137, "AGC", "A", U, "0/1"), (142, "T", "TGGG", U, "0/1")],
    [(303, "T", "A", U, "1/0"), (308, "T", "C", U, "0/1")],
    [(303, "T", "A", P, "0|1:100"), (308, "T", "C", P, "0|1:100")],
    [(303, "T", "A", P, "0|1:100"), (308, "T", "C", P, "0|1:200")],
    [(303, "T", "A", P, "0|1:100"), (308, "T", "C", P, "1|0:100")],
    [(303, "T", "A", P, "0|1:100"), (308, "T", "C", U, "0/1")],
    [(303, "T", "A", U, "0/1"), (308, "T", "C", U, "1/1")],
    [(303, "T", "A", U, "./1"), (308, "T", "C", U, "./1")],
    [(303, "T", "A", P, "1|.:100"), (308, "T", "C", P, "1|0:100")],
    [(303, "T", "A", P, "1|0:100"), (308, "T", "C", P, "1|.:100")],
    [(303, "T", "A", P, ".|1:100"), (308, "T", "C", P, ".|1:100")],
    [(303, "T", "A", P, "1|.:100"), (308, "T", "C", P, "0|1:100")],
    [(303, "TGT", "T", U, "0/1"), (304, "G", "A", U, "0/1")],
    [(303, "T", "A", U, "0/1"), (313, "G", "T", U, "0/1")],
    [(303, "T", "A", U, "0/1"), (314, "C", "A", U, "0/1")],
    [(310, "C", "T", U, "0/1")],
    [(303, "C", "A", U, "0/1"), (308, "T", "C", U, "0/1")],
    [(301, "CCTGTAATCCCA", "C", U, "0/1"), (305, "T", "G", U, "0/1"),
     (308, "T", "C", U, "0/1"), (314, "C", "A", U, "0/1")],
    [(303, "T", "A", U, "0/1"), (306, "A", "G", U, "0/1"), (309, "C", "T", U, "0/1")],
    [(303, "t", "a", U, "0/1"), (308, "t", "c", U, "0/1")],
    [(303, "T", "<DEL>", U, "0/1"), (304, "G", "*", U, "0/1"), (305, "T", "A,G", U, "0/1"),
     (306, "A", ".", U, "0/1"), (308, "T", "C", U, "0/1")],
    [(303, "T", "T", U, "0/1"), (308, "T", "T", U, "0/1")],
    [(303, "T", "A", U, "1/1"), (308, "T", "C", U, "1/1")],
    [(303, "T", "A", U, "1"), (308, "T", "C", U, "1")],
    [(303, "T", "A", U, "0/10"), (308, "T", "C", U, "10/0")],
    [(303, "T", "A", "GT:GQ:DP", "0/1:40:30"), (308, "T", "C", "DP:GT", "30:1/0")],
    [(303, "T", "A", "GQ:PS:GT", "40:7:1|0"), (308, "T", "C", "GT:PS:PS", "1|0:9:7")],
    [(303, "T", "A", "GT:PS", "0|1:."), (308, "T", "C", "GT", "0|1")],
    [(303, "T", "A", U, "0/1"), (303, "T", "G", U, "0/1"), (305, "T", "C", U, "0/1"),
     (308, "T", "C", U, "0/1")],
    [(137, "AGC", "A", U, "0/1"), (142, "T", "TGGG", U, "0/1"), (303, "T", "A", U, "0/1")],
    [(303, "T", "A", U, "0/0"), (308, "T", "C", U, "0/0")],
    [(303, "T", "A", U, "0|1"), (308, "T", "C", U, "0|1")],
]
CONTIGS = [f"c{i + 1:02d}" for i in range(len(SCEN))] + ["lc"]

with open("ref.fa", "w") as f, open("ref.fa.fai", "w") as fai:
    off = 0
    for name in CONTIGS:
        s = SEQ.lower() if name == "lc" else SEQ
        head = f">{name}\n"
        off += len(head)
        fai.write(f"{name}\t{len(s)}\t{off}\t60\t61\n")
        body = "".join(s[i:i + 60] + "\n" for i in range(0, len(s), 60))
        f.write(head + body)
        off += len(body)

HEAD = ("##fileformat=VCFv4.2\n"
        + "".join(f"##contig=<ID={c},length=450>\n" for c in CONTIGS)
        + '##INFO=<ID=PGP_JOINED_FROM,Number=1,Type=String,Description="an older run">\n'
        + '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        + '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
        + '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        + '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n')
COLS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"


def rec(c, pos, ref, alt, fmt, smp):
    return f"{c}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\t{fmt}\t{smp}\n"


with open("join.vcf", "w") as f:
    f.write(HEAD + COLS)
    for c, recs in zip(CONTIGS, SCEN):
        for r in recs:
            f.write(rec(c, *r))
    for r in SCEN[19] + SCEN[0]:
        f.write(rec("lc", *r))

with open("join-edge.vcf", "w") as f:
    f.write(HEAD + COLS)
    f.write(rec("c01", 137, "AGC", "A", U, "0/1"))
    f.write("##a header line among the records\n")
    f.write("#a comment\n\n")
    f.write("c01\t140\t.\tC\tT\t50\tPASS\t.\n")
    f.write(rec("c01", 142, "T", "TGGG", U, "0/1"))
    f.write(COLS)
    f.write(rec("c02", 303, "T", "A", U, "0/1"))
    f.write(rec("c02", 308, "T", "C", U, "0/1").rstrip("\n"))


def filled(pos, ref, alt, seq, frm):
    return f"c01\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tPGP_JOIN_SEQ={seq};PGP_JOINED_FROM={frm}\tGT\t0/1\n"


with open("filled.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n"
            '##INFO=<ID=PGP_JOIN_SEQ,Number=1,Type=String,Description="Step 06, temporary">\n'
            '##INFO=<ID=PGP_JOINED_FROM,Number=1,Type=String,Description="Step 06">\n' + COLS)
    f.write(filled(303, "TNNNNT", "ANNNNC", ref_at(303, 6), "c01:303:T:A|c01:308:T:C"))
    f.write(filled(303, "TNNNNT", "ANNNNC", ref_at(303, 6).lower(), "c01:303:t:a|c01:308:t:c"))
    f.write(filled(303, "CNNNNT", "ANNNNC", ref_at(303, 6), "c01:303:C:A|c01:308:T:C"))
    f.write(filled(303, "TNNNNG", "ANNNNC", ref_at(303, 6), "c01:303:T:A|c01:308:G:C"))
    f.write(filled(303, "TNNNNT", "ANNNNC", ref_at(303, 5), "c01:303:T:A|c01:308:T:C"))
    f.write(filled(303, "TGT", "T", ref_at(303, 3), "c01:303:TGT:T|c01:304:G:A"))
    f.write(filled(303, "TNNNNT", "TNNNNT", ref_at(303, 6), "c01:303:T:T|c01:308:T:T"))
    f.write(filled(303, "TNNNNT", "ANNNNC", ref_at(303, 6), "c01:303:T:A|c01:308:T:C|c01:310:C:T"))
    f.write(filled(303, "TNNNNT", "ANNNNC", ref_at(303, 6), "c01:0303:T:A|c01:308:T:C"))
    f.write("c01\t137\t.\tAGCNNNNT\tANNNNTGGG\t.\tPASS\tPGP_JOINED_FROM=c01:137:AGC:A|c01:142:T:TGGG;"
            f"PGP_JOIN_SEQ={ref_at(137, 6)}\tGT\t0/1\n")
    f.write("c01\t137\t.\tAGC\tA\t.\tPASS\t.\tGT\t0/1\n")
    f.write("c01\t137\t.\tAGC\n")

CV = "GENEINFO={g};CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter"
TAG = [
    ("chr22", 37, "9900001", "AGC", "A", CV.format(g="SYNTH1:1"), "GT", "0/1"),
    ("chr22", 38, "9900003", "GCCAT", "CATGGG", "PGP_JOINED_FROM=chr22:37:AGC:A|chr22:42:T:TGGG;"
     + CV.format(g="SYNTH1:1"), "GT", "0/1"),
    ("chr22", 39, "9900013", "CCAT", "ATGGG", "PGP_JOINED_FROM=chr22:37:AGC:A|chr22:43:A:G;"
     + CV.format(g="SYNTH1:1"), "GT", "0/1"),
    ("chr22", 42, "9900002", "T", "TGGG", CV.format(g="SYNTH1:1"), "GT", "0/1"),
    ("chr22", 210, "9900004", "C", "T", CV.format(g="SYNTH2:2"), "GT", "0/1"),
    ("chr22", 300, "9900005", "A", "G", CV.format(g="SYNTH3:3"), "GT", "0/1"),
    ("chr22", 305, "9900006", "C", "T", CV.format(g="SYNTH3:3"), "GT", "0/1"),
    ("chr22", 308, "9900016", "CAT", "C", CV.format(g="SYNTH9:9|SYNTH3:3"), "GT", "1/0"),
    ("chr22", 400, "9900007", "A", "G", CV.format(g="SYNTH4:4"), "GT", "0/1"),
    ("chr22", 404, "9900008", "C", "T", CV.format(g="SYNTH5:5"), "GT", "0/1"),
    ("chr22", 500, "9900009", "A", "G", CV.format(g="SYNTH6:6"), "GT:PS", "0|1:7"),
    ("chr22", 504, "9900010", "C", "T", CV.format(g="SYNTH6:6"), "GT:PS", "1|0:7"),
    ("chr22", 600, "9900011", "A", "G", CV.format(g="SYNTH7:7"), "GT:PS", "1|.:8"),
    ("chr22", 605, "9900012", "C", "T", CV.format(g="SYNTH7:7"), "GT:PS", "1|0:8"),
    ("chr22", 700, "9900014", "A", "G", "GENEINFO=SYNTH8:8", "GT", "0/1"),
    ("chr22", 711, "9900015", "C", "T", "GENEINFO=SYNTH8:8", "GT", "0/1"),
    ("chr22", 720, "9900017", "A", "G", ".", "GT", "0/1"),
    ("chr22", 722, "9900018", "C", "T", "GENEINFO=", "GT", "0/1"),
    ("chr23", 300, "9900019", "A", "G", CV.format(g="SYNTH3:3"), "GT", "0/1"),
    ("chr23", 302, "9900020", "C", "T", CV.format(g="SYNTH3:3"), "GT", "0/1"),
    ("chr23", 304, "9900021", "C", "T", CV.format(g="SYNTH3:3"), "GT", "0/1"),
]
with open("tag.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n" + COLS)
    for t in TAG:
        f.write("\t".join(map(str, t[:5])) + "\t50\tPASS\t" + "\t".join(t[5:]) + "\n")
    f.write("chr23\t310\t.\tA\tG\t50\tPASS\t" + CV.format(g="SYNTH3:3") + "\n")

STATUSES = ["practice_guideline", "reviewed_by_expert_panel",
            "criteria_provided,_multiple_submitters,_no_conflicts", "criteria_provided,_single_submitter",
            "criteria_provided,_conflicting_classifications", "criteria_provided,_conflicting_interpretations",
            "no_assertion_criteria_provided", "no_classification_provided", "Practice_guideline", ""]
GTS = [("GT:GQ:DP", "0/1:50:30"), ("GT:DP:GQ", "1|0:9:50"), ("GT:GQ:DP", "1/1:29:10"), ("GT", "1/."),
       ("GT:GQ", ".|2:."), ("GT", "2/2"), ("GT:GQ:DP", "0/0:30:9"), ("GT:DP", "01/1:x"), ("GT:GQ:GQ", "0|1:5:40"),
       ("DP:GT", "12:./1")]
HHEAD = "##fileformat=VCFv4.2\n" + COLS
with open("hits.vcf", "w") as f:
    f.write(HHEAD)
    for i, st in enumerate(STATUSES):
        info = f"GENEINFO=G{i}:{i}|H{i}:9;CLNSIG=Likely_pathogenic;CLNREVSTAT={st}" if st else "CLNSIG=Pathogenic"
        fmt, smp = GTS[i]
        f.write(f"chr1\t{100 + i}\t{i}\tC\tT\t.\tPASS\t{info}\t{fmt}\t{smp}\n")
    f.write("chr1\t200\t.\tAG\tA\t.\tPASS\tPGP_JOINED_FROM=chr1:200:A:G|chr1:203:C:T;CLNREVSTAT=reviewed_by_expert_panel"
            "\tGT\t0/1\n")
    f.write("chr1\t201\t.\tG\tA\t.\tPASS\tPGP_PART_OF=chr1:200:AG:A|chr1:199:TG:T;CLNSIG=Pathogenic\tGT\t1|0\n")
    f.write("chr1\t202\t.\tG\tA\t.\tPASS\tPGP_ONE_ALLELE=chr1:209:C:A;CLNSIG=x;CLNSIG=Benign;GENEINFO=Z:1\tGT\t0/1\n")
    f.write("chr1\t203\t.\tG\tA\t.\tPASS\tPGP_PART_OF=chr1:200:AG:A;PGP_JOINED_FROM=chr1:1:A|x;flag\tGT\t0/2\n")
    f.write("chr1\t204\t.\tG\tA\t.\tPASS\tPGP_ONE_ALLELE=;GENEINFO=:1\n")
    f.write("chr1\t205\t.\tG\n")
with open("hits-floor.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##PGP_CALL_FLOOR=GQ=45,DP=40\n" + COLS)
    f.write("chr1\t100\t.\tC\tT\t.\tPASS\tCLNREVSTAT=practice_guideline\tGT:GQ:DP\t0/1:44:40\n")
    f.write("chr1\t101\t.\tC\tT\t.\tPASS\tCLNREVSTAT=practice_guideline\tGT:GQ:DP\t0/1:45:39\n")
    f.write("chr1\t102\t.\tC\tT\t.\tPASS\tCLNREVSTAT=practice_guideline\tGT:GQ:DP\t0/1:45:40\n")
    f.write("##PGP_CALL_FLOOR=GQ=x,DP=5,GQ=7\n")
    f.write("chr1\t103\t.\tC\tT\t.\tPASS\tCLNREVSTAT=practice_guideline\tGT:GQ:DP\t0/1:6:5\n")
    f.write("chr1\t104\t.\tC\tT\t.\tPASS\tCLNREVSTAT=practice_guideline\tGT:GQ:DP\t0/1:7:4\n")
with open("hits-gene.vcf", "w") as f:
    f.write(HHEAD)
    f.write("chr1\t100\t.\tC\tT\t.\tPASS\tGENEINFO=A:1||B:2\tGT\t0/1\n")
    f.write("chr1\t101\t.\tC\tT\t.\tPASS\tGENEINFO=|\tGT\t0/1\n")
    f.write("chr1\t102\t.\tC\tT\t.\tPASS\tGENEINFO=A:1|:2\tGT\t0/1\n")

with open("sites.tsv", "w") as f:
    f.write("# the always-report list of the test\nchrom\tpos\tref\talt\tgene\thgvs\tclinvar_id\n")
    for i, (c, p, r, a) in enumerate([
            ("chr1", 100, "A", "G"), ("chr1", 200, "C", "T"), ("chr1", 300, "G", "A"), ("chr1", 400, "T", "C"),
            ("chr1", 500, "A", "T"), ("chr1", 600, "A", "T"), ("chr1", 700, "A", "T"), ("chr1", 800, "A", "T"),
            ("chr1", 900, "A", "T"), ("chr2", "0100", "A", "G"), ("chr1", 1000, "A", "T"), ("chr1", 1100, "A", "T"),
            ("chr1", 1200, "A", "T")]):
        f.write(f"{c}\t{p}\t{r}\t{a}\tGENE{i}\tc.{i}A>T\t{i}\n")
    f.write("chr1\t1300\tA\tT\tshort\n")
with open("sites-vcf.txt", "w") as f:
    f.write("chr1\t100\t100\tA\tG\t0/1\t40\t30\n"
            "chr1\t200\t200\tC\tT\t0/0\t50\t30\n"
            "chr1\t300\t300\tG\tA\t0|0\t12\t30\n"
            "chr1\t400\t400\tT\tG\t0/1\t40\t30\n"
            "chr1\t400\t400\tT\tA\t1/1\t40\t30\n"
            "chr1\t900\t900\tA\tT\t./.\t.\t.\n"
            "chr2\t100\t100\tA\tG\t1|0\t.\t20\n"
            "chr1\t1000\t1000\tA\tC\t0/1\t40\t30\n"
            "chr1\t1000\t1000\tA\tT\t./1\t40\t30\n")
with open("sites-gvcf.txt", "w") as f:
    f.write("chr1\t490\t510\tA\t<*>\t0/0\t40\t.\n"
            "chr1\t600\t600\tA\t<*>\t0/0\t5\t.\n"
            "chr1\t690\t710\tA\t<*>\t./.\t0\t.\n"
            "chr1\t895\t905\tA\t<*>\t0/0\t.\t.\n"
            "chr1\t1100\t1100\tA\tT,<*>\t0/1\t30\t20\n"
            "chr1\t1150\t1250\tA\t<*>\t0|0\t30\t.\n")
open("empty.txt", "w").close()
PY
