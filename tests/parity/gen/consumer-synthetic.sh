#!/usr/bin/env bash
# consumer-synthetic.sh: writes synthetic.g.vcf.gz (BGZF, about 2 million
# records) and synthetic-sites.tsv.gz (about 200,000 sites) in the current
# directory, the same every time. The gVCF mixes what DeepVariant writes:
# reference blocks, het and hom SNVs, deletions with the block that starts
# inside them, a variant and a block at one position, no-call blocks and
# filtered calls, on chr1, chr2, chrX, chrY and chrM; a few sites carry
# another reference base. scripts/ci/oracle-diff.sh runs it for the case
# consumer export-synthetic.
set -euo pipefail
python3 - <<'PY'
import gzip, struct, zlib

state = 12345
def rnd(n):
    global state
    state = (state * 6364136223846793005 + 1442695040888963407) % (1 << 64)
    return (state >> 33) % n

def base(contig, pos):
    return "ACGT"[(pos * 2654435761 + len(contig)) % 7 % 4]

out, sites = [], []
out.append("##fileformat=VCFv4.2\n")
for c in ("chr1", "chr2", "chrX", "chrY", "chrM"):
    out.append(f"##contig=<ID={c}>\n")
out.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
chroms = {"chr1": "1", "chr2": "2", "chrX": "X", "chrY": "Y", "chrM": "MT"}
sizes = {"chr1": 10800000, "chr2": 8400000, "chrX": 6000000, "chrY": 1800000, "chrM": 16569}
n = 0
for contig, size in sizes.items():
    pos = 1 + rnd(50)
    while pos < size:
        r = base(contig, pos)
        k = rnd(100)
        gt1 = "1" if contig in ("chrY", "chrM") and rnd(2) else None
        if k < 62:
            end = min(size, pos + rnd(40))
            gt = "0" if gt1 else "0/0"
            out.append(f"{contig}\t{pos}\t.\t{r}\t<*>\t0\t.\tEND={end}\tGT:GQ:MIN_DP\t{gt}:30:20\n")
            pos = end + 1
        elif k < 85:
            alt = "ACGT"[(("ACGT".index(r)) + 1 + rnd(3)) % 4]
            gt = gt1 or ("0/1", "1/1", "0|1")[rnd(3)]
            flt = "PASS" if rnd(20) else "LowQual"
            out.append(f"{contig}\t{pos}\t.\t{r}\t{alt},<*>\t40\t{flt}\t.\tGT:GQ\t{gt}:40\n")
            if rnd(5) == 0:
                out.append(f"{contig}\t{pos}\t.\t{r}\t<*>\t0\t.\tEND={pos + 3}\tGT:GQ\t0/0:20\n")
                pos += 3
            pos += 1
        elif k < 93:
            ln = 2 + rnd(4)
            ref = "".join(base(contig, p) for p in range(pos, pos + ln))
            gt = ("0/1", "1/1", "0/0")[rnd(3)]
            flt = "RefCall" if gt == "0/0" else "PASS"
            out.append(f"{contig}\t{pos}\t.\t{ref}\t{ref[0]},<*>\t30\t{flt}\t.\tGT:GQ\t{gt}:30\n")
            out.append(f"{contig}\t{pos + 1}\t.\t{base(contig, pos + 1)}\t<*>\t0\t.\tEND={pos + ln + 2}\tGT:GQ\t0/0:20\n")
            pos += ln + 3
        elif k < 97:
            end = min(size, pos + rnd(10))
            out.append(f"{contig}\t{pos}\t.\t{r}\t<*>\t0\t.\tEND={end}\tGT:GQ\t./.:0\n")
            pos = end + 1
        else:
            pos += 1 + rnd(30)
        n += 1
    for p in range(1 + rnd(5), size, 1):
        if rnd(150):
            continue
        ref38 = base(contig, p) if rnd(400) else "ACGT"[("ACGT".index(base(contig, p)) + 1) % 4]
        sites.append((f"rs{len(sites) + 1}", chroms[contig], p + 10000, ref38, p, ref38, "-" if rnd(4) == 0 else "+"))

data = "".join(out).encode()
with open("synthetic.g.vcf.gz", "wb") as f:
    for i in range(0, len(data), 65280):
        chunk = data[i:i + 65280]
        co = zlib.compressobj(6, zlib.DEFLATED, -15)
        body = co.compress(chunk) + co.flush()
        f.write(b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", len(body) + 25))
        f.write(body + struct.pack("<II", zlib.crc32(chunk), len(chunk)))
    f.write(bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000"))
head = (f"## template_sites: {len(sites) + 1000}\n## dropped: 1000\n#rsid\tchrom\tpos37\tref37\tpos38\tref38\tstrand\n")
with gzip.open("synthetic-sites.tsv.gz", "wt") as f:
    f.write(head + "".join("\t".join(str(x) for x in s) + "\n" for s in sites))
print(f"{n} records, {len(sites)} sites")
PY
