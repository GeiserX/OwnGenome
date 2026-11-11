#!/usr/bin/env bash
# synthetic-vep.sh: writes synthetic VEP VCFs in the current directory, the
# same every time (a fixed seed), for the og vcf csq-rank cases:
#   vep-edge.vcf     20,000 records built to hit each branch of the rule:
#                    one, two or no MANE and canonical entries with equal
#                    impacts (the first must win), unknown impacts, entries
#                    shorter than the Format, '&'-joined, signed, nan, inf,
#                    hex and exponent frequencies, '.' and old PGP_ items in
#                    INFO, GT, GQ, DP and AD in any order and any shape,
#                    records with no FORMAT or no sample, CRLF lines
#   vep-gnomad.vcf   the same kind of records with gnomADe_AF and
#                    gnomADg_AF and no MAX_AF
set -euo pipefail
export LC_ALL=C
python3 - <<'PY'
import random

R = random.Random(20261110)
FULL = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE",
        "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position",
        "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "SYMBOL_SOURCE",
        "HGNC_ID", "CANONICAL", "MANE_SELECT", "MANE_PLUS_CLINICAL", "SIFT", "PolyPhen", "MAX_AF",
        "gnomADe_AF", "gnomADg_AF", "CLIN_SIG"]
IMPACTS = ["HIGH", "MODERATE", "LOW", "MODIFIER", "", "high", "OTHER"]
FREQS = ["", "0.001", "0.5", "1e-3", "1E+2", ".5", "5.", "0.01", "0.0099", "0.010", "nan", "inf",
         "-0.1", "&", "0.2&0.004", "&0.3", "0.1&", "1e999", "00.5", "0x10", " 0.1", "0"]


def header(fields):
    return ("##fileformat=VCFv4.2\n"
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
            f'Format: {"|".join(fields)}">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")


def entry(fields, t):
    v = dict.fromkeys(fields, "")
    v["Allele"] = "G"
    v["Consequence"] = R.choice(["stop_gained", "missense_variant", "intron_variant", "", "a&b"])
    v["IMPACT"] = R.choice(IMPACTS)
    v["SYMBOL"] = R.choice(["GENE1", "GENE2", ""])
    v["Feature"] = f"ENST{t:011d}"
    v["CANONICAL"] = R.choice(["YES", "YES", "", "NO"])
    v["MANE_SELECT"] = R.choice(["NM_1.1", "", "", ""])
    for f in ("MAX_AF", "gnomADe_AF", "gnomADg_AF"):
        if f in v:
            v[f] = R.choice(FREQS)
    out = [v[f] for f in fields]
    if R.random() < 0.1:
        out = out[:R.randint(0, len(out))]
    return "|".join(out) if R.random() > 0.03 else ""


def info(fields):
    items = [R.choice([".", "DP=5", "DB", "", "PGP_IMPACT=HIGH", "PGP_CALL_WHY=GQ1", "PGP_CALL=ok",
                       "PGP_IMPACTX=1", "PGP_TX", "AF=0.5"]) for _ in range(R.randint(0, 3))]
    if R.random() < 0.9:
        csq = ",".join(entry(fields, t) for t in range(R.randint(0, 6)))
        items.insert(R.randint(0, len(items)), "CSQ=" + csq)
    return ";".join(items) if R.random() > 0.03 else "."


def sample():
    vals = {"GT": R.choice(["0/1", "1/1", "0|1", "1|1", "1/2", "0/0", "./.", "1", "0", "2/2", "1|0",
                            "01/1", "0/1/2", "", ".", "/1"]),
            "GQ": R.choice(["50", "12", "30", "29", "0", "", ".", "5.0", "007"]),
            "DP": R.choice(["30", "9", "10", "", "."]),
            "AD": R.choice(["15,15", "3,2", "12,5", "0,30", "0,0", "1", "5,", "", "10,0,5", "1,2,3",
                            "0,1", "40,10", "2,8"])}
    keys = ["GT", "GQ", "DP", "AD"]
    R.shuffle(keys)
    keys = keys[:R.randint(0, 4)]
    vs = [vals[k] for k in keys]
    if R.random() < 0.15:
        vs = vs[:R.randint(0, len(vs))]
    return ":".join(keys), ":".join(vs)


def edge(path, fields, n):
    with open(path, "w", newline="") as f:
        f.write(header(fields))
        for i in range(n):
            cols = ["chr1", str(100 + i), ".", "A", "G", "50", "PASS", info(fields)]
            k = R.random()
            if k < 0.05:
                pass
            elif k < 0.1:
                cols.append(sample()[0])
            else:
                cols += list(sample())
            line = "\t".join(cols)
            f.write(line + ("\r\n" if R.random() < 0.01 else "\n"))


edge("vep-edge.vcf", FULL, 20000)
edge("vep-gnomad.vcf", [x for x in FULL if x != "MAX_AF"], 5000)
PY
