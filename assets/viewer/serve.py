#!/usr/bin/env python3
"""What is left of the local genome viewer's Python once it moved to
`og view serve` (rust/crates/og-live/src/viewer.rs, started by
scripts/owngenome-view.sh). Only the names the live site's viewer,
bin/owngenome_live_viewer.py, imports are kept: sample_files, build_config,
send_file and LIVE_CSP, and what they call. That plugin, and this file with
it, go when the live site's server moves to og.

Both viewers pick Manta's diploid SVs before step 22's consensus. The local
viewer (og) also refuses a reversed locus; the live site's viewer never takes
a locus from the command line.
"""

import os
import re
import urllib.parse

CSP = ("default-src 'none'; script-src 'self'; connect-src 'self'; img-src 'self' data: blob:; "
       "style-src 'self' 'unsafe-inline'; font-src 'self' data:; worker-src 'self' blob:; "
       "base-uri 'none'; form-action 'none'; frame-ancestors 'none'")
# Only the live site's viewer can be embedded by another page of that site.
LIVE_CSP = CSP.replace("frame-ancestors 'none'", "frame-ancestors 'self'")
LOCUS_RE = re.compile(r"^([A-Za-z0-9._-]+)(?::([0-9][0-9,]*)(?:-([0-9][0-9,]*))?)?$")


def read_fai(fai):
    """Contig names and lengths of a .fai, in order."""
    out = []
    with open(fai) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0]:
                out.append((p[0], int(p[1])))
    if not out:
        raise ValueError(f"{fai} names no contig")
    return out


def check_locus(locus, contigs):
    """The locus as igv.js reads it (chr, chr:pos or chr:start-end), checked
    against the reference; the first contig when none is given."""
    if not locus:
        return contigs[0][0]
    m = LOCUS_RE.match(locus)
    if not m:
        raise ValueError(f"'{locus}' is not a locus: give chr1, chr1:12345 or chr1:12000-13000")
    names = dict(contigs)
    if m.group(1) not in names:
        raise ValueError(f"'{m.group(1)}' is not a contig of the reference")
    for x in (m.group(2), m.group(3)):
        if x and not 1 <= int(x.replace(",", "")) <= names[m.group(1)]:
            raise ValueError(f"{x} is outside {m.group(1)} (1 to {names[m.group(1)]})")
    return locus


def index_of(names):
    """The first of these index files that exists."""
    return next((n for n in names if os.path.isfile(n)), None)


def sample_files(sample_dir, sample):
    """(reads, vcf, sv_vcf) of a sample folder, by the names the steps write,
    None for one that is not there: the default aligned/ and vcf/ folders, as
    owngenome-view.sh picks them with none of READS, VCF, SV_VCF, ALIGN_DIR or
    VCF_DIR set, and the ones the live site's viewer (bin/owngenome_live_viewer.py)
    serves."""
    def first(*names):
        return next((p for p in (os.path.join(sample_dir, n) for n in names) if os.path.isfile(p)), None)
    return (first(f"aligned/{sample}_sorted.bam", f"aligned/{sample}_sorted.cram"),
            first(f"vcf/{sample}.vcf.gz"),
            first("manta/results/variants/diploidSV.vcf.gz", f"sv_merged/{sample}_sv_consensus.vcf.gz"))


def build_config(sample, fasta, reads=None, vcf=None, sv_vcf=None, locus=None):
    """(config, files): the igv.js browser options and the list of files the
    server may send, /data/<index>/<basename> each."""
    files = []

    def url(path):
        files.append(os.path.abspath(path))
        return f"data/{len(files) - 1}/{urllib.parse.quote(os.path.basename(path))}"

    fai = fasta + ".fai"
    if not os.path.isfile(fasta) or not os.path.isfile(fai):
        raise ValueError(f"the reference needs its FASTA and .fai: {fasta}")
    contigs = read_fai(fai)
    cfg = {"sample": sample,
           "locus": check_locus(locus, contigs),
           "reference": {"id": "reference", "name": os.path.basename(fasta),
                         "fastaURL": url(fasta), "indexURL": url(fai)},
           "tracks": []}
    if reads:
        cram = reads.endswith(".cram")
        stem = os.path.splitext(reads)[0]
        idx = index_of([reads + ".crai", stem + ".crai"] if cram
                       else [reads + ".bai", stem + ".bai", reads + ".csi"])
        if not idx:
            raise ValueError(f"no index beside {reads}: run samtools index on it")
        cfg["tracks"].append({"name": f"{sample} reads", "type": "alignment", "format": "cram" if cram else "bam",
                              "url": url(reads), "indexURL": url(idx)})
    for path, name in ((vcf, f"{sample} small variants"), (sv_vcf, f"{sample} structural variants")):
        if not path:
            continue
        t = {"name": name, "type": "variant", "format": "vcf", "url": url(path)}
        idx = index_of([path + ".tbi", path + ".csi"])
        if idx:
            t["indexURL"] = url(idx)
        else:
            t["indexed"] = False
        cfg["tracks"].append(t)
    return cfg, files


def parse_range(header, size):
    """(start, end) inclusive for a single 'bytes=' range, None when there is
    none or it is not one we serve whole-file for, 'bad' when unsatisfiable."""
    m = re.match(r"^bytes=(\d*)-(\d*)$", (header or "").strip())
    if not m or (not m.group(1) and not m.group(2)):
        return None
    if m.group(1):
        start = int(m.group(1))
        end = int(m.group(2)) if m.group(2) else size - 1
    else:
        start, end = max(0, size - int(m.group(2))), size - 1
    if start >= size or end < start:
        return "bad"
    return start, min(end, size - 1)


def send_file(handler, head, ctype, src):
    """Send the file src, or the one byte range the request's Range header
    asks for (206; 416 past the end), 1 MB at a time. head(code, ctype,
    length, extra headers) writes the status line and headers."""
    size = os.path.getsize(src)
    rng = parse_range(handler.headers.get("Range"), size)
    if rng == "bad":
        return head(416, "text/plain; charset=utf-8", 0, [("Content-Range", f"bytes */{size}")])
    start, end = rng if rng else (0, size - 1)
    length = end - start + 1 if size else 0
    extra = [("Accept-Ranges", "bytes")] + ([("Content-Range", f"bytes {start}-{end}/{size}")] if rng else [])
    head(206 if rng else 200, ctype, length, extra)
    if handler.command == "HEAD" or not length:
        return
    with open(src, "rb") as f:
        f.seek(start)
        left = length
        while left > 0:
            chunk = f.read(min(1 << 20, left))
            if not chunk:
                break
            handler.wfile.write(chunk)
            left -= len(chunk)
