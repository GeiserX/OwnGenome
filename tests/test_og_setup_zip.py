#!/usr/bin/env python3
"""og setup unzip and og setup genetic-map against the Python they replaced.

scripts/setup.sh unpacked Exomiser's data zips with
zipfile.ZipFile(z).extractall(dest), and cut the PLINK genetic map of steps
26b and 51 out of its zip with two small inline programs (REF_26B below,
copied from setup.sh before og). Zips are built here with Python's zipfile;
og and the Python run on each, and the test checks:

  1. unzip: stored and deflated members, ZIP64 records, nested and empty
     directories, a member named twice (the last one wins), names with `.`
     and empty parts, a non-ASCII name: the same tree, file for file and
     byte for byte, as extractall writes;
  2. unzip refuses, with exit 2 and nothing written, a member whose name
     climbs out of the folder (`../x`, `/x`), where extractall silently
     writes it elsewhere; and a member with a flipped byte fails its CRC
     check (exit 2);
  3. genetic-map: CHROM.map is the member byte for byte and CHROM.rfmix.map
     the same as the Python wrote (fields 1, 4 and 3 of each non-blank line,
     whitespace runs, tabs, CRLF and a last line without a newline), for
     two chromosomes in one call; without --rfmix only the .map is written;
  4. genetic-map stops with exit 2 on a row with fewer than four fields and
     on a chromosome the zip lacks, and leaves the files of that chromosome
     as they were (the Python wrote .map and an empty .rfmix.map there);
  5. the control: a tree with one byte changed is told apart from the
     Python's, so the comparison can fail.

Run: OG=rust/target/release/og python3 tests/test_og_setup_zip.py
"""
import os
import subprocess
import sys
import tempfile
import warnings
import zipfile

OG = os.environ.get("OG", "")
FAILS = []

# scripts/setup.sh's step 26b map program before og (run as python3 -c).
REF_26B = '''
import sys, zipfile
z, c, out = sys.argv[1:]
with zipfile.ZipFile(z) as f:
    text = f.read(f"chr_in_chrom_field/plink.chr{c}.GRCh38.map").decode()
rows = [line.split() for line in text.splitlines() if line.strip()]
open(f"{out}/{c}.map", "w").write(text)
open(f"{out}/{c}.rfmix.map", "w").write("".join(f"{r[0]}\\t{r[3]}\\t{r[2]}\\n" for r in rows))
'''


def check(desc, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS.append(desc)


def og(*args):
    return subprocess.run([OG, *args], capture_output=True, text=True)


def tree(root):
    """{relative path: bytes, or None for a directory}."""
    out = {}
    for d, dirs, files in os.walk(root):
        for x in dirs:
            out[os.path.relpath(os.path.join(d, x), root)] = None
        for x in files:
            p = os.path.join(d, x)
            out[os.path.relpath(p, root)] = open(p, "rb").read()
    return out


def make_zip(path, members, zip64=False):
    with zipfile.ZipFile(path, "w") as z:
        for name, data, method in members:
            if name.endswith("/"):
                z.writestr(zipfile.ZipInfo(name), b"")
                continue
            info = zipfile.ZipInfo(name)
            info.compress_type = method
            with z.open(info, "w", force_zip64=zip64) as f:
                f.write(data)


def part_unzip(tmp):
    big = os.urandom(300_000) + b"A" * 2_000_000
    members = [
        ("2602_phenotype/", b"", None),
        ("2602_phenotype/phenotype.mv.db", big, zipfile.ZIP_DEFLATED),
        ("2602_phenotype/rw/stored.txt", b"stored as is\n", zipfile.ZIP_STORED),
        ("2602_phenotype/empty/", b"", None),
        ("2602_phenotype/twice.txt", b"first\n", zipfile.ZIP_DEFLATED),
        ("2602_phenotype/twice.txt", b"second\n", zipfile.ZIP_STORED),
        ("2602_phenotype/./odd//parts.txt", b"parts\n", zipfile.ZIP_DEFLATED),
        ("2602_phenotype/café.txt", b"utf-8 name\n", zipfile.ZIP_DEFLATED),
    ]
    for zip64 in (False, True):
        z = os.path.join(tmp, f"data{int(zip64)}.zip")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # the duplicate name, on purpose
            make_zip(z, members, zip64)
        py, ours = os.path.join(tmp, f"py{int(zip64)}"), os.path.join(tmp, f"og{int(zip64)}")
        zipfile.ZipFile(z).extractall(py)
        r = og("setup", "unzip", z, ours)
        a, b = tree(py), tree(ours)
        kind = "ZIP64" if zip64 else "plain"
        check(f"unzip ({kind}): exit 0", r.returncode == 0, r.stderr)
        check(f"unzip ({kind}): the same tree as extractall, byte for byte ({len(a)} entries)",
              a == b, sorted(set(a) ^ set(b)) or [k for k in a if a[k] != b.get(k)])
        if not zip64:
            changed = dict(b)
            changed["2602_phenotype/twice.txt"] = b"secone\n"
            check("control: a tree with one byte changed differs from extractall's", changed != a)

    for bad in ("../evil.txt", "top/../../evil.txt", "/abs/evil.txt"):
        z = os.path.join(tmp, "bad.zip")
        make_zip(z, [("fine.txt", b"x", zipfile.ZIP_STORED), (bad, b"y", zipfile.ZIP_STORED)])
        dest = os.path.join(tmp, "bad", "out")
        r = og("setup", "unzip", z, dest)
        check(f"unzip refuses the member {bad!r} (exit 2, nothing written)",
              r.returncode == 2 and not os.path.exists(dest) and not os.path.exists(os.path.join(tmp, "bad", "evil.txt")),
              (r.returncode, r.stderr))

    z = os.path.join(tmp, "crc.zip")
    make_zip(z, [("a.txt", b"some content here", zipfile.ZIP_STORED)])
    data = bytearray(open(z, "rb").read())
    at = data.find(b"some content here")
    data[at] ^= 1
    open(z, "wb").write(data)
    r = og("setup", "unzip", z, os.path.join(tmp, "crc"))
    check("unzip: a flipped byte fails the CRC check (exit 2)", r.returncode == 2 and "CRC" in r.stderr,
          (r.returncode, r.stderr))


MAP20 = "chr20 . 0.0 60000\nchr20\t.\t0.12\t70001\n\n   \nchr20  .  1.5   80000\n"
MAP21 = "chr21 . 0 1\r\nchr21 . 0.5 2\r\nchr21 . 1.25 3"


def part_genetic_map(tmp):
    z = os.path.join(tmp, "plink.GRCh38.map.zip")
    make_zip(z, [("chr_in_chrom_field/", b"", None),
                 ("chr_in_chrom_field/plink.chrchr20.GRCh38.map", MAP20.encode(), zipfile.ZIP_DEFLATED),
                 ("chr_in_chrom_field/plink.chrchr21.GRCh38.map", MAP21.encode(), zipfile.ZIP_DEFLATED),
                 ("chr_in_chrom_field/plink.chrchr22.GRCh38.map", b"chr22 . 0.0 1\nchr22 . 0.4\n", zipfile.ZIP_DEFLATED)])
    py, ours = os.path.join(tmp, "py-maps"), os.path.join(tmp, "og-maps")
    os.makedirs(py)
    os.makedirs(ours)
    for c in ("chr20", "chr21"):
        subprocess.run([sys.executable, "-c", REF_26B, z, c, py], check=True)
    r = og("setup", "genetic-map", "--zip", z, "--out", ours, "--rfmix", "chr20", "chr21")
    check("genetic-map --rfmix chr20 chr21: exit 0", r.returncode == 0, r.stderr)
    a, b = tree(py), tree(ours)
    check("genetic-map: the same .map and .rfmix.map files as the Python, byte for byte",
          a == b, sorted(set(a) ^ set(b)) or [k for k in a if a[k] != b.get(k)])
    check("control: the chr21 .rfmix.map holds three rows", a.get("chr21.rfmix.map", b"").count(b"\n") == 3, a)

    only = os.path.join(tmp, "only-map")
    os.makedirs(only)
    r = og("setup", "genetic-map", "--zip", z, "--out", only, "chr20")
    check("genetic-map without --rfmix writes only chr20.map",
          r.returncode == 0 and sorted(os.listdir(only)) == ["chr20.map"], (r.returncode, os.listdir(only)))

    keep = os.path.join(tmp, "keep")
    os.makedirs(keep)
    open(os.path.join(keep, "chr22.rfmix.map"), "w").write("old\n")
    r = og("setup", "genetic-map", "--zip", z, "--out", keep, "--rfmix", "chr22")
    check("genetic-map: a row with three fields stops it (exit 2) and names the line",
          r.returncode == 2 and "line 2" in r.stderr, (r.returncode, r.stderr))
    check("genetic-map: and leaves that chromosome's files as they were",
          sorted(os.listdir(keep)) == ["chr22.rfmix.map"]
          and open(os.path.join(keep, "chr22.rfmix.map")).read() == "old\n", os.listdir(keep))
    r = og("setup", "genetic-map", "--zip", z, "--out", keep, "chrX")
    check("genetic-map: a chromosome the zip lacks stops it (exit 2)", r.returncode == 2, (r.returncode, r.stderr))


def main():
    if not (OG and os.access(OG, os.X_OK)):
        print("skip: OG is not set (OG=rust/target/release/og)")
        return 1 if os.environ.get("OG_REQUIRED") else 0
    with tempfile.TemporaryDirectory() as tmp:
        part_unzip(tmp)
        part_genetic_map(tmp)
    if FAILS:
        print(f"{len(FAILS)} check(s) failed")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
