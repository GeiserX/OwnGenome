#!/usr/bin/env bash
# summary-tools.sh: writes the inputs of the qc-sample, prs-table and mito-*
# parity cases in the current directory, the same every time: somalier and
# VerifyBamID2 tables (a pedigree sex left standing, X and Y disagreeing, a
# low-quality sample, a contaminated one, a count of inf), pgsc_calc output
# folders (no panel, a panel with its ancestry, a chrX row its matches did
# not keep, a below-threshold log) with the score folders og prs format
# writes, and step 20's chrM VCF, a ClinVar VCF with its tabix index,
# gnomAD's chrM sites, eKLIPse CSVs and samtools idxstats. The fixtures of
# the old unit tests of bin/collect_summary.py's step tools.
# scripts/ci/oracle-diff.sh runs it.
set -euo pipefail
export LC_ALL=C
python3 - <<'PY'
import gzip
import os
import struct
import zlib


def put(path, text, gz=False):
    if os.path.dirname(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    data = text.encode() if isinstance(text, str) else text
    if gz:
        data = gzip.compress(data, mtime=0)
    with open(path, "wb") as f:
        f.write(data)


# --- step 33: somalier relate and VerifyBamID2 -------------------------------------

HEAD = ("#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\toriginal_pedigree_sex\t"
        "gt_depth_mean\tn_hom_ref\tn_het\tn_hom_alt\tp_middling_ab\tX_depth_mean\tX_n\tX_hom_ref\t"
        "X_het\tX_hom_alt\tY_depth_mean\tY_n\n")


def row(sid="S", sex="1", ped="-9", x_n="39", x_het="0", x_hom="39", y_depth="14.0", mid="0.010",
        hom_ref="60"):
    return (f"fam\t{sid}\t-9\t-9\t{sex}\t-9\t{ped}\t30.0\t{hom_ref}\t40\t20\t{mid}\t15.0\t{x_n}\t0\t"
            f"{x_het}\t{x_hom}\t{y_depth}\t5\n")


put("male.samples.tsv", HEAD + row() + row("T", sex="2", x_n="40", x_het="20", x_hom="20"))
put("ped-kept.samples.tsv", HEAD + row(ped="male", x_n="2", x_hom="2"))
put("ped-female.samples.tsv", HEAD + row(ped="female"))
put("lowq.samples.tsv", HEAD + row(ped="male", mid="0.080") + row("U", ped="male", mid="nan"))
put("xy.samples.tsv", HEAD + row(sex="-2", ped="female", x_n="40", x_het="20", x_hom="20"))
put("xy-unknown.samples.tsv", HEAD + row(sex="-2", ped="female", x_n="4", x_hom="4"))
put("female.samples.tsv", HEAD + row(sex="2", x_n="40", x_het="20", x_hom="20", y_depth="0.1"))
put("inf.samples.tsv", HEAD + row(hom_ref="inf"))
put("crlf.samples.tsv", (HEAD + "\n" + row()).replace("\n", "\r\n"))
put("s.samples.tsv.gz", HEAD + row(), gz=True)
put("clean.selfSM", "#SEQ_ID\tRG\tCHIP_ID\t#SNPS\t#READS\tAVG_DP\tFREEMIX\tFREELK1\nS\tNA\tNA\t10000\t5000\t2.1\t0.0012\t1\n")
put("dirty.selfSM", "#SEQ_ID\tRG\tCHIP_ID\t#SNPS\t#READS\tAVG_DP\tFREEMIX\nS\tNA\tNA\t10000\t5000\t2.1\t0.0412\n")
put("nofreemix.selfSM", "#SEQ_ID\tFREEMIX\nS\tNA\n")
put("s.pairs.tsv", "#sample_a\tsample_b\trelatedness\tibs0\nS\tT\t0.95\t0\nU\tS\t0.991\t0\nS\tV\t0.5\t9\n"
    "S\tS\t1.0\t0\nW\tS\tnan\t0\n\nT\tU\t0.99\t0\n")

# --- steps 25 and 26: pgsc_calc's output and og prs format's folders -----------------


def score_folder(d, pid, trait, rows, chrx="1\tleft out", applies="", evals="", amb="0\tdropped"):
    put(f"{d}/{pid}.txt.gz", f"#pgs_id={pid}\n#pgs_name={pid}\n" + (f"#trait_reported={trait}\n" if trait else "")
        + "#genome_build=GRCh38\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n"
        + "1\t1100\tA\tG\t0.5\n" * rows, gz=True)
    with open(f"{d}/chrx_rows.tsv", "a") as f:
        if f.tell() == 0:
            f.write("pgs_id\tchrx_rows\tchrx\n")
        f.write(f"{pid}\t{chrx}\n")
    with open(f"{d}/score_notes.tsv", "a") as f:
        if f.tell() == 0:
            f.write("pgs_id\tapplies_to\tambiguous_rows\tambiguous\n")
        f.write(f"{pid}\t{applies}\t{amb}\n")
    with open(f"{d}/score_evaluations.tsv", "a") as f:
        if f.tell() == 0:
            f.write("pgs_id\tppm_id\tancestry\teffect\n")
        if evals:
            f.write(f"{pid}\t{evals}\n")


os.makedirs("pgs", exist_ok=True)
score_folder("pgs", "PGS000001", "Coronary artery disease", 3, evals="PPM000591;PPM000594\tEUR;AFR\tHR 1.53;HR 1.27")
score_folder("pgs", "PGS000662", "Prostate cancer", 4, chrx="0\tleft out", applies="male",
             evals="PPM022982;PPM015482\tEUR;AFR\tHR 1.72;OR 1.86", amb="2\tkept")
put("pgs/PGS000002.txt.gz", "#pgs_id=PGS000002\n#genome_build=GRCh38\n"
    "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n" + "1\t1\tA\t\t1\n" * 5, gz=True)
os.makedirs("pgs-plain", exist_ok=True)
score_folder("pgs-plain", "PGS000001", "Coronary artery disease", 3)
os.makedirs("pgs-old", exist_ok=True)
put("pgs-old/PGS000001.txt.gz", "#trait_reported=Old\nchr_name\tchr_position\teffect_allele\teffect_weight\n"
    "1\t1\tA\t1\n", gz=True)
os.makedirs("pgs-x", exist_ok=True)
score_folder("pgs-x", "PGS000001", "Coronary artery disease", 4, chrx="3\tscored")
os.makedirs("pgs-bad-evals", exist_ok=True)
score_folder("pgs-bad-evals", "PGS000001", "Coronary artery disease", 3, evals="PPM000591\tEuropean\tHR 1.53")

SUMMARY = ("dataset,accession,ambiguous,is_multiallelic,match_flipped,duplicate_best_match,duplicate_ID,"
           "match_IDs,match_status,count,score_pass,match_rate\n"
           "sample,PGS000001,false,false,false,false,false,NA,matched,2,true,0.67\n"
           "sample,PGS000001,false,false,false,false,false,NA,unmatched,1,true,0.67\n"
           "sample,PGS000662_hmPOS_GRCh38,false,false,false,false,false,NA,matched,3.0,true,0.75\n"
           "sample,PGS000002,false,false,false,false,false,NA,unmatched,5,false,0.0\n")
put("raw/sample/match/sample_summary.csv", SUMMARY)
put("raw/sample/score/aggregated_scores.txt.gz", "sampleset\tFID\tIID\tPGS\tSUM\tDENOM\n"
    "sample\tHG\tHG\tPGS000001\t62.0\t6\nsample\tHG\tHG\tPGS000662\t-0.000012345678901\t8\n", gz=True)
put("adj/sample/match/sample_summary.csv", SUMMARY)
put("adj/sample/score/sample_pgs.txt.gz", "sampleset\tFID\tIID\tPGS\tSUM\tZ_MostSimilarPop\tpercentile_MostSimilarPop\n"
    "sample\tHG\tHG\tPGS000001\t1.25\t1.6\t94.44\nsample\tHG\tHG\tPGS000662\t3\t0.2\t\n"
    "reference\tR1\tR1\tPGS000001\t10.0\t0.1\t12.3\n", gz=True)
put("adj/sample/score/sample_popsimilarity.txt.gz",
    "sampleset\tFID\tIID\tPC1\tPC2\tPC10\tSuperPop\tRF_P_AFR\tRF_P_EUR\tMostSimilarPop\t"
    "MostSimilarPop_LowConfidence\tREFERENCE\n"
    "reference\tR1\tR1\t-1\t2\t0\tAFR\t0.9\t0.1\tAFR\tFalse\tTrue\n"
    "sample\tHG\tHG\t0.012\t-0.034\t0.5\t\t0.05\t0.95\tEUR\tFalse\tFalse\n", gz=True)
LOG = "row_nr,accession,chr_name,chr_position,effect_allele,other_allele,match_status,dataset\n"
put("adj/sample/match/sample_log.csv.gz", LOG + "0,PGS000001,1,1100,A,G,matched,sample\n"
    "1,PGS000001_hmPOS_GRCh38,X,51000000,A,G,matched,sample\n2,PGS000001,X,51000100,A,G,excluded,sample\n"
    '3,PGS000001,23,51000200,A,G,"matched",sample\n', gz=True)
put("two-targets/sample/match/sample_summary.csv", SUMMARY)
put("two-targets/sample/score/sample_pgs.txt.gz", "sampleset\tFID\tIID\tPGS\tSUM\n"
    "sample\tHG\tHG\tPGS000001\t1.25\n", gz=True)
put("two-targets/sample/score/sample_popsimilarity.txt.gz", "sampleset\tPC1\tMostSimilarPop\tREFERENCE\n"
    "sample\t1\tEUR\tFalse\nsample\t2\tEUR\tFalse\n", gz=True)
put("below.log", "ERROR Score PGS000001 fails minimum matching threshold (33.33% variants match)\n"
    "ERROR Score PGS000662_hmPOS_GRCh38 fails minimum matching threshold (7% variants match)\n"
    "ERROR pgscatalog.core.lib.pgsexceptions.ZeroMatchesError: All scores fail to meet match threshold 0.75\n")
os.makedirs("empty", exist_ok=True)

# --- step 20: the chrM calls, ClinVar, gnomAD, eKLIPse and idxstats ------------------


def bgzf_block(data):
    c = zlib.compressobj(6, zlib.DEFLATED, -15)
    cdata = c.compress(data) + c.flush()
    return (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
            + struct.pack("<H", len(cdata) + 25) + cdata
            + struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data)))


def tabixed_vcf(path, header, blocks):
    """A bgzipped VCF of one BGZF block per item of blocks and its .tbi: bin 0
    of each sequence spans all its lines, and the linear index."""
    raw, lin, span, coff = b"", {}, {}, 0
    for i, lines in enumerate(blocks):
        data = (header if i == 0 else "").encode()
        for line in lines:
            chrom, pos = line.split("\t")[0], int(line.split("\t")[1])
            w, ix = (pos - 1) >> 14, lin.setdefault(chrom, [])
            while len(ix) <= w:
                ix.append(None)
            start = (coff << 16) | len(data)
            if ix[w] is None:
                ix[w] = start
            data += (line + "\n").encode()
            first = span.get(chrom, (start, 0))[0]
            span[chrom] = (first, (coff << 16) | len(data))
        blk = bgzf_block(data)
        raw += blk
        coff += len(blk)
    with open(path, "wb") as f:
        f.write(raw + bgzf_block(b""))
    names = list(lin)
    nm = b"".join(n.encode() + b"\0" for n in names)
    tbi = b"TBI\x01" + struct.pack("<8i", len(names), 2, 1, 2, 0, ord("#"), 0, len(nm)) + nm
    for n in names:
        ix, last = lin[n], 0
        for j, v in enumerate(ix):
            ix[j] = last = v if v is not None else last
        tbi += struct.pack("<iIi", 1, 0, 1) + struct.pack("<QQ", *span[n])
        tbi += struct.pack("<i", len(ix)) + struct.pack(f"<{len(ix)}Q", *ix)
    with open(path + ".tbi", "wb") as f:
        f.write(gzip.compress(tbi, mtime=0))


put("S_chrM_filtered.vcf.gz", "##fileformat=VCFv4.2\n##contig=<ID=chrM,length=16569>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    "chrM\t73\t.\tA\tG\t.\tPASS\t.\tGT:AF\t1:1.00\n"
    "chrM\t302\t.\tAC\tA,ACC\t.\tPASS\t.\tGT:AF\t0/1/2:0.3,0.4\n"
    "chrM\t3243\t.\tA\tG\t.\tPASS\t.\tGT:AF\t0/1:0.3\n"
    "chrM\t9000\t.\tC\tT\t.\tstrand_bias\t.\tGT:AF\t0/1:0.1\n"
    "chrM\t16189\t.\tT\tC\t.\tPASS\t.\tGT\t0/1\n"
    "chrM\t16190\t.\tCCA\tCGA\t.\tPASS\t.\tGT:DP\t0/1:5\n"
    "chrM\t16200\t.\tA\n", gz=True)
CV = "##fileformat=VCFv4.1\n##fileDate=2026-10-05\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
tabixed_vcf("clinvar_chr.vcf.gz", CV, [
    ["chr1\t100\t1\tC\tT\t.\t.\tCLNSIG=Benign"],
    ["chrM\t73\t9\tA\tG\t.\t.\tCLNSIG=not_provided;CLNREVSTAT=no_classification_provided;CLNDN=not_provided",
     "chrM\t302\t8\tA\tAC\t.\t.\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter;CLNDN=not_provided",
     "chrM\t303\t10\tC\tCC,T\t.\t.\tCLNSIG=Likely_benign;CLNDN=x",
     "chrM\t3243\t7\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;CLNDN=MELAS",
     "chrM\t16189\t6\tT\tA\t.\t.\tCLNSIG=Uncertain_significance;CLNREVSTAT=criteria_provided,_single_submitter;"
     "CLNDN=not_provided",
     "chrM\t16191\t5\tA\tG\t.\t.\tCLNSIG=Benign;FLAG"]])
tabixed_vcf("clinvar_mt.vcf.gz", CV.replace("##fileDate=2026-10-05\n", ""), [
    ["MT\t3243\t7\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;CLNDN=MELAS"]])
tabixed_vcf("clinvar_nochrm.vcf.gz", CV, [["chr1\t100\t1\tC\tT\t.\t.\tCLNSIG=Benign"]])
put("gnomad.chrM.vcf.bgz", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chrM\t302\t.\tA\tAC,ACC\t.\tPASS\tAN=56434;AF_hom=0.01,0.02;AF_het=0.2;max_hl=0.99\n"
    "chrM\t3243\t.\tA\tG\t.\tPASS\tAN=56434;AF_hom=0.00002;AF_het=0.0003;max_hl=0.9\n"
    "chrM\t16190\t.\tC\tG\t.\tPASS\tAN=1\nshort\tline\n", gz=True)
put("bad-pos.vcf", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrM\tx\t.\tA\tG\t.\tPASS\t.\n")

EK = ('"Title";"5\' breakpoint";"3\' breakpoint";"Freq";"Freq For";"Freq Rev";"5\' Blast";'
      '"3\' Blast";"5\' Depth";"3\' Depth";"Repetition"\n')


def ek_row(b5, b3, lvl, rep="none", depth="390.0"):
    v = str(lvl).replace(".", ",")
    return f'"S";"{b5}";"{b3}";"{v}";"{v}";"{v}";"12";"14";"{depth}";"402,7";"{rep}"\n'


put("eKLIPse_deletions.csv", EK + ek_row(8469, 13447, 10.2) + ek_row(8482, 13460, 2.5) + ek_row(9000, 14000, 3.1)
    + ek_row(5000, 7000, 0.4567, rep="") + ek_row(16000, 500, 2.0, rep="ACCTA"))
put("crlf_deletions.csv", (EK + ek_row(8469, 13447, 10.2)).replace("\n", "\r\n"))
put("empty_deletions.csv", EK)
put("notek.csv", "pos\tref\n")
put("short_deletions.csv", EK + '"S";"1";"2"\n')
put("inf_deletions.csv", EK + ek_row(1, 2, "inf"))
put("nan_depth_deletions.csv", EK + ek_row(1, 2, 1.0, depth="nan"))


def idx(chrm, per, m="chrM", prefix="chr"):
    lines = [f"{prefix}{i}\t50000000\t{per}\t0" for i in range(1, 23)]
    lines += [f"{prefix}X\t156040895\t{per}\t0", f"{m}\t16569\t{chrm}\t0", "chrUn_x\t1000\t7\t0", "*\t0\t0\t500"]
    return "\n".join(lines) + "\n"


put("wgs.idxstats", idx(331380, 10_000_000))
put("mt.idxstats", idx(331380, 10_000_000, m="MT", prefix=""))
put("chrmt.idxstats", idx(331380, 10_000_000, m="chrMT"))
put("slice.idxstats", idx(331380, 500))
put("nochrm.idxstats", idx(0, 10_000_000))
put("crlf.idxstats", idx(165690, 10_000_000).replace("\n", "\r\n"))
put("bad.idxstats", "chr1\t100\tmany\t0\n")
PY
