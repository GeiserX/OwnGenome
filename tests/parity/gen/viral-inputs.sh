#!/usr/bin/env bash
# The inputs of the viral-report cases: a viral FASTA of five small genomes and,
# per case, samtools coverage of the viral and the human reads and samtools
# depth -aa of the viral ones, with chosen depths.
set -euo pipefail
python3 - <<'PY'
COV_HEAD = "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
GENOMES = [("NC_001664.4", "Human betaherpesvirus 6A, variant A DNA, complete virion genome, isolate U1102", 2000),
           ("NC_000898.1", "Human herpesvirus 6B, complete genome", 2000),
           ("NC_007605.1", "Human gammaherpesvirus 4, complete genome", 1500),
           ("NC_004102.1", "Hepatitis C virus genotype 1, complete genome", 900),
           ("NC_001422.1", "Escherichia phage phiX174, complete genome", 500)]
A, B, EBV, HCV, PHIX = (g[0] for g in GENOMES)


def put(name, text):
    with open(name, "w") as f:
        f.write(text)


def viral(tag, depths, reads, human=None):
    cov, dep = [COV_HEAD], []
    for acc, _, length in GENOMES:
        fn = depths.get(acc, lambda p: 0)
        d = [fn(p) for p in range(1, length + 1)]
        c = sum(1 for x in d if x)
        cov.append(f"{acc}\t1\t{length}\t{reads.get(acc, 0)}\t{c}\t{100 * c / length:.4g}\t{sum(d) / length:.6g}\t30\t60\n")
        dep += [f"{acc}\t{p}\t{x}\n" for p, x in enumerate(d, 1)]
    if human is None:
        human = [("chr1", 100000, 100000, 30.0), ("chr22", 50000, 50000, 30.0)]
    h = [COV_HEAD] + [f"{r}\t1\t{ln}\t1000\t{cb}\t{100 * cb / ln:.4g}\t{m:.6g}\t30\t60\n" for r, ln, cb, m in human]
    put(f"{tag}.vcov.tsv", "".join(cov))
    put(f"{tag}.vdepth.tsv", "".join(dep))
    put(f"{tag}.hcov.tsv", "".join(h))


flat = lambda x: (lambda p: x)  # noqa: E731
put("v.fa", "".join(f">{acc} {title}\n" + "A" * length + "\n\n" for acc, title, length in GENOMES))
put("crlf.fa", "".join(f">{acc} {title}\r\n" + "A" * length + "\r\n" for acc, title, length in GENOMES))
viral("ici", {B: flat(15)}, {B: 400})
viral("relative-30", {B: flat(9)}, {B: 400})
viral("pileup", {B: lambda p: 100 if p < 50 else 0}, {B: 30})
viral("low", {EBV: lambda p: 1 if p % 2 else 0, HCV: flat(2)}, {EBV: 40, HCV: 25})
viral("few", {EBV: lambda p: 1 if p <= 100 else 0}, {EBV: 12})
viral("phix", {PHIX: flat(31), B: flat(15), A: flat(17)}, {PHIX: 900, B: 400, A: 420})
viral("no-human", {B: flat(15)}, {B: 400}, human=[("chrX", 1000, 1000, 30.0)])
viral("human-mix", {}, {}, human=[("chr1", 1000, 500, 10.0), ("22", 1000, 1000, 40.0), ("chrX", 1000, 1000, 99.0),
                                  ("chrUn_KI270302v1", 1000, 1000, 99.0), ("chrEBV", 1000, 1000, 99.0)])
lines = open("ici.vdepth.tsv").readlines()
put("cut.vdepth.tsv", "".join(lines[:-10]))
put("gone.vdepth.tsv", "".join(x for x in lines if not x.startswith(HCV + "\t")))
put("bad.hcov.tsv", "chr1\t1\t2\n")
PY
