#!/usr/bin/env bash
# synthetic-vep-500k.sh: writes vep-500k.vcf in the current directory, the
# same every time (a fixed seed): 500,000 VEP-like records with 1 to 8 CSQ
# entries of 32 fields, about 230 MB, with MANE and canonical entries of
# equal impact, for the og vcf csq-rank case synthetic-500k (ties at scale,
# and its speed).
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
with open("vep-500k.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n"
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
            f'Format: {"|".join(FULL)}">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
    imps = ["MODIFIER", "LOW", "MODERATE", "HIGH"]
    nf = len(FULL)
    for i in range(500000):
        trs = []
        mane = R.randint(-1, 7)
        for t in range(R.randint(1, 8)):
            v = [""] * nf
            v[0], v[1], v[2], v[3], v[6] = "G", "intron_variant", R.choice(imps), f"GENE{t}", f"ENST{t:011d}"
            v[23] = "YES" if R.random() < 0.4 else ""
            v[24] = "NM_1.1" if t == mane or R.random() < 0.05 else ""
            v[28] = f"{R.random() / 10:.5f}"
            v[29] = f"{R.random() / 10:.4g}"
            trs.append("|".join(v))
        gq, dp, a, b = R.randint(0, 60), R.randint(0, 50), R.randint(0, 30), R.randint(0, 30)
        gt = R.choice(["0/1", "1/1", "0|1", "1/2"])
        f.write(f"chr1\t{1000 + i * 10}\t.\tA\tG\t50\tPASS\tAC=1;AN=2;CSQ={','.join(trs)}"
                f"\tGT:GQ:DP:AD\t{gt}:{gq}:{dp}:{a},{b}\n")
PY
