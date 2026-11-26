//! `og dev install-manifest`: what an install tree holds, file by file, and
//! the compare of two such lists. It is the check every slice of the Rust
//! port of `scripts/setup.sh` must pass: install one component with bash
//! and with og into two empty folders, and the two trees must hold the
//! same content (`tests/setup/ab.sh`).
//!
//! # The manifest
//!
//! `og dev install-manifest DIR [--exclude PATH]...` writes to stdout a
//! `#og-install-manifest 1` line, the column line
//! `path type size sha256 mode digests`, then one row per entry under DIR,
//! sorted by path bytes. `type` is `f` (a file: its size and sha256), `d`
//! (a folder), `l` (a symlink: `->` and its target in the sha256 column) or
//! `o` (anything else). `mode` is the permission bits in octal (`-` for a
//! symlink). Paths are relative to DIR, so trees at two places have equal
//! manifests. Left out: the folder `downloads` at the top (raw downloads,
//! kept or not by choice), every file whose name ends in `.log` (tool logs,
//! with times and hosts in them), and each `--exclude` path with all under
//! it. In a name, a tab, a newline, a carriage return, a backslash, any
//! other control byte and any byte that is not UTF-8 are written as `\t`,
//! `\n`, `\r`, `\\` and `\xHH`.
//!
//! `digests` holds, for a file of a kind a compare mode reads, the sha256 of
//! what that mode compares (`kind=hex`, comma separated, `-` for none). The
//! manifest carries them so a tree can be compared after it is gone (the
//! reference manifest of a full install):
//!
//! - `gzip`: the decompressed content of a gzip or BGZF file (every
//!   member), so the header's time and OS byte do not count. `bad` when it
//!   does not decompress.
//! - `vcf`: a VCF (plain or gzip) with the value of its `##fileDate=` line
//!   and the `Date=` part of each `##bcftools_<cmd>Command=` line masked.
//!   Everything else, the bcftools version and command lines included, is
//!   compared.
//! - `index`: a `.tbi` or `.csi`: its header (format, columns, the CSI shift
//!   and depth) and its sequence names, when the index is valid for the data
//!   file beside it; `invalid` when it is not. Valid means every chunk ends
//!   inside the data file and the chunk starts checked (up to 32, spread over
//!   the index) begin a BGZF block. Bins and offsets are left out: they move
//!   with any byte of the data file, which the data file's own row compares.
//! - `jellyfish`: a jellyfish file (9 digits, then its JSON header) without
//!   that header and its alignment padding: the header holds the host name,
//!   the folder and the time of the count.
//! - `numeric`: a text file of at most 16 MiB with a decimal number in it:
//!   its lines and whitespace-separated words, each word with a decimal
//!   point or an exponent rounded to 4 significant digits, every other word
//!   as it is. Two runs of a multithreaded float program agree on it unless a
//!   value sits on a rounding edge, which makes the compare red, never green.
//!
//! # The compare
//!
//! `og dev install-manifest --compare A.tsv B.tsv [--allow ALLOW.tsv]`
//! prints one line per difference, `DIFF<TAB>kind<TAB>path<TAB>detail`, with
//! kind `removed` (in A only), `added` (in B only), `type`, `content` (a
//! file's bytes, or what its compare mode reads, or a symlink's target) or
//! `mode`. A difference in bytes that the path's compare mode finds equal is
//! printed as `ALLOWED` with the mode and the row's reason. Lines starting
//! with `#` are comments in every input.
//!
//! The allowlist (`tests/setup/install-manifest-allow.tsv`) has the columns
//! `path mode component reason cannot_compare`. `path` is a pattern over the
//! whole relative path: `*` and `?` match within one segment, `**` across
//! segments; the first matching row counts. `mode` is one of `gzip`, `vcf`,
//! `index`, `jellyfish`, `numeric` or `skip` (content not compared), and may
//! add `perm` (`gzip,perm`, or `perm` alone): the permission bits may
//! differ. `component` and `reason` must not be empty; `cannot_compare` says
//! why no other mode can compare the file and is required for `skip`, `-`
//! otherwise. A stamp (`DONE`, `INSTALLED`, `RELEASE`, `version.txt`,
//! `*.stamp`, `*.installed`, `*.release`) is always compared byte for byte:
//! its bytes decide whether the next setup run reinstalls, so a row whose
//! last segment can match a stamp name is refused, and no row applies to a
//! stamp.
//!
//! # For tests/setup/ab.sh
//!
//! `--control B.tsv [--allow ALLOW.tsv]` prints the smallest non-empty file
//! compared byte for byte; `--flip DIR PATH` flips one bit of the byte in
//! the middle of DIR/PATH (twice restores it). `--compare-logs A B` compares
//! two setup logs with times, durations, speeds and download progress
//! masked, and the missing list (`Missing (...)`, `  - item`, `retry:`)
//! line by line.
//!
//! Exit 0 when nothing differs outside the allowlist, 1 when something
//! does, 2 when an argument, a manifest or the allowlist is wrong or a file
//! cannot be read.

use flate2::read::MultiGzDecoder;
use regex::Regex;
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, HashMap};
use std::fmt::Write as _;
use std::fs;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::os::unix::ffi::OsStrExt;
use std::os::unix::fs::PermissionsExt;
use std::path::Path;
use std::sync::LazyLock;

const USAGE: &str = "usage: og dev install-manifest DIR [--exclude PATH]...
       og dev install-manifest --compare A.tsv B.tsv [--allow ALLOW.tsv]
       og dev install-manifest --control B.tsv [--allow ALLOW.tsv]
       og dev install-manifest --flip DIR PATH
       og dev install-manifest --compare-logs A.log B.log
";
const MAGIC: &str = "#og-install-manifest 1";
const COLUMNS: &str = "path\ttype\tsize\tsha256\tmode\tdigests";
const ALLOW_COLUMNS: &str = "path\tmode\tcomponent\treason\tcannot_compare";
/// Files up to this size are read once into memory; larger ones are
/// streamed, a second time for their gzip digests.
const SMALL: u64 = 16 << 20;
/// Names a row may never cover: a sample of every stamp shape.
const STAMP_SAMPLES: [&str; 10] = [
    "DONE",
    "INSTALLED",
    "RELEASE",
    "version.txt",
    "catalog.stamp",
    "chr1.stamp",
    "counts.stamp",
    "gnomad-sv.installed",
    "mane_select_refseq.release",
    "x.stamp",
];

/// Runs `og dev install-manifest ARGS`; returns the exit code.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let result = match args.first().map(String::as_str) {
        Some("--compare") => compare_cli(&args[1..], out),
        Some("--control") => control_cli(&args[1..], out),
        Some("--flip") => flip_cli(&args[1..], out),
        Some("--compare-logs") => logs_cli(&args[1..], out),
        _ => manifest_cli(args, out, err),
    };
    match result {
        Ok(code) => code,
        Err(msg) => {
            if msg.is_empty() {
                let _ = write!(err, "{USAGE}");
            } else {
                let _ = writeln!(err, "og dev install-manifest: {msg}");
            }
            2
        }
    }
}

fn shown(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn io_err(p: &Path) -> impl Fn(std::io::Error) -> String + '_ {
    move |e| format!("{}: {e}", shown(p))
}

// ---------------------------------------------------------------- manifest

fn manifest_cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> Result<u8, String> {
    let mut dir = None;
    let mut exclude = Vec::new();
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--exclude" => {
                let p = it.next().ok_or("--exclude needs a path")?;
                exclude.push(p.trim_matches('/').as_bytes().to_vec());
            }
            "-h" | "--help" => return Err(String::new()),
            s if s.starts_with('-') => return Err(format!("unknown option {s}")),
            _ if dir.is_some() => return Err(String::new()),
            _ => dir = Some(a.clone()),
        }
    }
    let dir = dir.ok_or_else(String::new)?;
    let root = Path::new(&dir);
    if !fs::symlink_metadata(root).map_err(io_err(root))?.is_dir() {
        return Err(format!("{dir} is not a folder"));
    }
    let mut w = Walk {
        exclude,
        rows: Vec::new(),
        skipped: 0,
        bytes: 0,
    };
    w.walk(root, &[])?;
    w.rows.sort();
    let mut text = format!("{MAGIC}\n{COLUMNS}\n");
    for (path, rest) in &w.rows {
        let _ = writeln!(text, "{}\t{rest}", escape(path));
    }
    out.write_all(text.as_bytes())
        .map_err(|e| format!("stdout: {e}"))?;
    let _ = writeln!(
        err,
        "og dev install-manifest: {} entries, {} bytes hashed, {} left out",
        w.rows.len(),
        w.bytes,
        w.skipped
    );
    Ok(0)
}

struct Walk {
    exclude: Vec<Vec<u8>>,
    /// (relative path bytes, the row's columns after the path)
    rows: Vec<(Vec<u8>, String)>,
    skipped: usize,
    bytes: u64,
}

impl Walk {
    fn walk(&mut self, dir: &Path, rel: &[u8]) -> Result<(), String> {
        for entry in fs::read_dir(dir).map_err(io_err(dir))? {
            let entry = entry.map_err(io_err(dir))?;
            let name = entry.file_name();
            let mut path = rel.to_vec();
            if !path.is_empty() {
                path.push(b'/');
            }
            path.extend_from_slice(name.as_bytes());
            let full = entry.path();
            let meta = fs::symlink_metadata(&full).map_err(io_err(&full))?;
            let ft = meta.file_type();
            if self.exclude.contains(&path)
                || (rel.is_empty() && name.as_bytes() == b"downloads")
                || (!ft.is_dir() && name.as_bytes().ends_with(b".log"))
            {
                self.skipped += 1;
                continue;
            }
            let mode = format!("{:04o}", meta.permissions().mode() & 0o7777);
            let row = if ft.is_symlink() {
                let target = fs::read_link(&full).map_err(io_err(&full))?;
                format!("l\t-\t->{}\t-\t-", escape(target.as_os_str().as_bytes()))
            } else if ft.is_dir() {
                format!("d\t-\t-\t{mode}\t-")
            } else if ft.is_file() {
                let (sha, digests) = digest_file(&full, meta.len())?;
                self.bytes += meta.len();
                format!("f\t{}\t{sha}\t{mode}\t{digests}", meta.len())
            } else {
                format!("o\t-\t-\t{mode}\t-")
            };
            self.rows.push((path.clone(), row));
            if ft.is_dir() {
                self.walk(&full, &path)?;
            }
        }
        Ok(())
    }
}

fn hex(h: Sha256) -> String {
    h.finalize()
        .iter()
        .fold(String::with_capacity(64), |mut s, b| {
            let _ = write!(s, "{b:02x}");
            s
        })
}

fn sha_of(bytes: &[u8]) -> String {
    let mut h = Sha256::new();
    h.update(bytes);
    hex(h)
}

/// A full install is hundreds of GB: keep each file out of the page cache
/// once read, so a manifest taken on a busy machine does not push its
/// work into swap.
fn drop_cache(f: &fs::File) {
    #[cfg(target_os = "linux")]
    let _ = rustix::fs::fadvise(f, 0, None, rustix::fs::Advice::DontNeed);
    #[cfg(not(target_os = "linux"))]
    let _ = f;
}

/// The sha256 of a file and its `digests` column.
fn digest_file(path: &Path, len: u64) -> Result<(String, String), String> {
    let e = io_err(path);
    let mut f = fs::File::open(path).map_err(&e)?;
    let mut head = [0u8; 64];
    let n = read_full(&mut f, &mut head).map_err(&e)?;
    let head = &head[..n];
    f.seek(SeekFrom::Start(0)).map_err(&e)?;
    let gz = head.starts_with(&[0x1f, 0x8b]);
    let name = path.file_name().map(|n| n.as_bytes()).unwrap_or_default();
    let mut digests: Vec<(&str, String)> = Vec::new();
    let sha;
    if len <= SMALL {
        let mut buf = Vec::with_capacity(len as usize);
        f.read_to_end(&mut buf).map_err(&e)?;
        sha = sha_of(&buf);
        if gz {
            let mut plain = Vec::new();
            match MultiGzDecoder::new(&buf[..]).read_to_end(&mut plain) {
                Ok(_) => {
                    digests.push(("gzip", sha_of(&plain)));
                    if plain.starts_with(b"##fileformat=VCF") {
                        digests.push(("vcf", vcf_digest(&plain[..]).map_err(&e)?));
                    }
                }
                Err(_) => digests.push(("gzip", "bad".into())),
            }
        } else if buf.starts_with(b"##fileformat=VCF") {
            digests.push(("vcf", vcf_digest(&buf[..]).map_err(&e)?));
        } else if let Some(d) = numeric_digest(&buf) {
            digests.push(("numeric", d));
        }
        if let Some(d) = jellyfish_digest(&buf[..]).map_err(&e)? {
            digests.push(("jellyfish", d));
        }
    } else {
        let mut h = Sha256::new();
        let mut buf = vec![0u8; 1 << 20];
        loop {
            let n = f.read(&mut buf).map_err(&e)?;
            if n == 0 {
                break;
            }
            h.update(&buf[..n]);
        }
        sha = hex(h);
        if gz {
            f.seek(SeekFrom::Start(0)).map_err(&e)?;
            digests.extend(gzip_stream_digests(&f));
        } else if head.starts_with(b"##fileformat=VCF") {
            f.seek(SeekFrom::Start(0)).map_err(&e)?;
            digests.push(("vcf", vcf_digest(BufReader::new(&f)).map_err(&e)?));
        }
        if is_jellyfish(head) {
            f.seek(SeekFrom::Start(0)).map_err(&e)?;
            if let Some(d) = jellyfish_digest(BufReader::new(&f)).map_err(&e)? {
                digests.push(("jellyfish", d));
            }
        }
    }
    drop_cache(&f);
    if name.ends_with(b".tbi") || name.ends_with(b".csi") {
        digests.push(("index", index_digest(path)));
    }
    let col = if digests.is_empty() {
        "-".to_owned()
    } else {
        digests
            .iter()
            .map(|(k, v)| format!("{k}={v}"))
            .collect::<Vec<_>>()
            .join(",")
    };
    Ok((sha, col))
}

fn read_full(f: &mut impl Read, buf: &mut [u8]) -> std::io::Result<usize> {
    let mut n = 0;
    while n < buf.len() {
        match f.read(&mut buf[n..])? {
            0 => break,
            k => n += k,
        }
    }
    Ok(n)
}

/// The gzip and, for a VCF, the vcf digest of a large gzip file, in one
/// streaming pass.
fn gzip_stream_digests(f: &fs::File) -> Vec<(&'static str, String)> {
    let mut r = BufReader::with_capacity(
        1 << 20,
        MultiGzDecoder::new(BufReader::with_capacity(1 << 20, f)),
    );
    let mut raw = Sha256::new();
    let mut vcf = Sha256::new();
    let mut is_vcf = None;
    let mut in_header = true;
    let mut line = Vec::new();
    loop {
        line.clear();
        match r.read_until(b'\n', &mut line) {
            Ok(0) => break,
            Ok(_) => {}
            Err(_) => return vec![("gzip", "bad".into())],
        }
        raw.update(&line);
        let v = *is_vcf.get_or_insert_with(|| line.starts_with(b"##fileformat=VCF"));
        if v {
            if in_header && !line.starts_with(b"##") {
                in_header = false;
            }
            vcf.update(if in_header {
                vcf_mask(&line)
            } else {
                std::borrow::Cow::Borrowed(&line[..])
            });
        }
    }
    let mut d = vec![("gzip", hex(raw))];
    if is_vcf == Some(true) {
        d.push(("vcf", hex(vcf)));
    }
    d
}

static BCFTOOLS_DATE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"^(##bcftools_\w+Command=.*; Date=).*$").expect("a valid pattern")
});

/// One VCF header line with its volatile part masked.
fn vcf_mask(line: &[u8]) -> std::borrow::Cow<'_, [u8]> {
    if line.starts_with(b"##fileDate=") {
        return std::borrow::Cow::Borrowed(b"##fileDate=<masked>\n");
    }
    if line.starts_with(b"##bcftools_") {
        let s = String::from_utf8_lossy(line);
        let t = s.trim_end_matches(['\n', '\r']);
        if let Some(c) = BCFTOOLS_DATE.captures(t) {
            return std::borrow::Cow::Owned(format!("{}<masked>\n", &c[1]).into_bytes());
        }
    }
    std::borrow::Cow::Borrowed(line)
}

fn vcf_digest(mut r: impl BufRead) -> std::io::Result<String> {
    let mut h = Sha256::new();
    let mut header = true;
    let mut line = Vec::new();
    loop {
        line.clear();
        if r.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if header && !line.starts_with(b"##") {
            header = false;
        }
        h.update(if header {
            vcf_mask(&line)
        } else {
            std::borrow::Cow::Borrowed(&line[..])
        });
    }
    Ok(hex(h))
}

fn is_jellyfish(head: &[u8]) -> bool {
    head.len() > 10 && head[..9].iter().all(u8::is_ascii_digit) && head[9] == b'{'
}

/// The digest of a jellyfish file's data, after the 9-digit length, the
/// JSON header and the padding to its `alignment`; None when it is not one.
fn jellyfish_digest(mut r: impl Read) -> std::io::Result<Option<String>> {
    let mut len = [0u8; 10];
    if read_full(&mut r, &mut len)? < 10 || !is_jellyfish(&[&len[..], b"x"].concat()) {
        return Ok(None);
    }
    let n: usize = std::str::from_utf8(&len[..9])
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    if n == 0 {
        return Ok(None);
    }
    let mut json = vec![0u8; n - 1];
    if read_full(&mut r, &mut json)? < json.len() {
        return Ok(None);
    }
    let json = [b"{".as_slice(), &json].concat();
    let align = Regex::new(r#""alignment"\s*:\s*(\d+)"#)
        .expect("a valid pattern")
        .captures(&String::from_utf8_lossy(&json))
        .and_then(|c| c[1].parse::<usize>().ok())
        .filter(|&a| a > 0)
        .unwrap_or(8);
    let pad = (align - (9 + n) % align) % align;
    let mut skip = vec![0u8; pad];
    read_full(&mut r, &mut skip)?;
    let mut h = Sha256::new();
    let mut buf = vec![0u8; 1 << 20];
    loop {
        let k = r.read(&mut buf)?;
        if k == 0 {
            break;
        }
        h.update(&buf[..k]);
    }
    Ok(Some(hex(h)))
}

/// The numeric digest of a small text file with a decimal number in it.
fn numeric_digest(buf: &[u8]) -> Option<String> {
    let text = std::str::from_utf8(buf).ok()?;
    if text.contains('\0') {
        return None;
    }
    let mut out = String::with_capacity(text.len());
    let mut any = false;
    for line in text.lines() {
        let mut first = true;
        for w in line.split_whitespace() {
            if !first {
                out.push(' ');
            }
            first = false;
            let decimal = w.contains(['.', 'e', 'E']) && w.bytes().any(|b| b.is_ascii_digit());
            match w.parse::<f64>() {
                Ok(v) if decimal && v.is_finite() => {
                    any = true;
                    let v = if v == 0.0 { 0.0 } else { v };
                    let _ = write!(out, "{v:.3e}");
                }
                _ => out.push_str(w),
            }
        }
        out.push('\n');
    }
    any.then(|| sha_of(out.as_bytes()))
}

/// The `index` digest of a .tbi or .csi: its header and names, or
/// `invalid` when it does not fit the data file beside it.
fn index_digest(path: &Path) -> String {
    match read_index(path) {
        Ok(Some(d)) => d,
        _ => "invalid".into(),
    }
}

struct Cursor<'a>(&'a [u8]);

impl<'a> Cursor<'a> {
    fn take(&mut self, n: usize) -> Option<&'a [u8]> {
        let s: &'a [u8] = self.0;
        if s.len() < n {
            return None;
        }
        let (a, b) = s.split_at(n);
        self.0 = b;
        Some(a)
    }
    fn i32(&mut self) -> Option<i32> {
        Some(i32::from_le_bytes(self.take(4)?.try_into().ok()?))
    }
    fn u64(&mut self) -> Option<u64> {
        Some(u64::from_le_bytes(self.take(8)?.try_into().ok()?))
    }
    fn count(&mut self) -> Option<usize> {
        usize::try_from(self.i32()?).ok()
    }
}

fn read_index(path: &Path) -> std::io::Result<Option<String>> {
    let s = path.as_os_str().as_bytes();
    let data = Path::new(std::ffi::OsStr::from_bytes(&s[..s.len() - 4]));
    let raw = fs::read(path)?;
    let mut plain = Vec::new();
    MultiGzDecoder::new(&raw[..]).read_to_end(&mut plain)?;
    let mut c = Cursor(&plain);
    let mut h = Sha256::new();
    let mut chunks: Vec<(u64, u64)> = Vec::new();
    let Some(magic) = c.take(4) else {
        return Ok(None);
    };
    h.update(magic);
    let csi = match magic {
        b"TBI\x01" => false,
        b"CSI\x01" => true,
        _ => return Ok(None),
    };
    let n_ref;
    // htslib's pseudo-bin: its two "chunks" hold offsets and record counts.
    let mut meta_bin = 37450u32;
    if csi {
        let (Some(shift), Some(depth), Some(l_aux)) = (c.take(4), c.take(4), c.count()) else {
            return Ok(None);
        };
        h.update(shift);
        h.update(depth);
        let depth = u32::from_le_bytes(depth.try_into().unwrap_or_default());
        meta_bin = (((1u64 << (3 * depth.min(9) + 3)) - 1) / 7 + 1) as u32;
        let Some(aux) = c.take(l_aux) else {
            return Ok(None);
        };
        h.update(aux);
        let Some(n) = c.count() else { return Ok(None) };
        n_ref = n;
    } else {
        let Some(n) = c.count() else { return Ok(None) };
        n_ref = n;
        let Some(fields) = c.take(24) else {
            return Ok(None);
        };
        h.update(fields);
        let Some(l_nm) = c.count() else {
            return Ok(None);
        };
        let Some(names) = c.take(l_nm) else {
            return Ok(None);
        };
        h.update(names);
    }
    h.update((n_ref as u64).to_le_bytes());
    for _ in 0..n_ref {
        let Some(n_bin) = c.count() else {
            return Ok(None);
        };
        for _ in 0..n_bin {
            let Some(bin) = c
                .take(4)
                .map(|b| u32::from_le_bytes(b.try_into().unwrap_or_default()))
            else {
                return Ok(None);
            };
            if csi && c.u64().is_none() {
                return Ok(None);
            }
            let Some(n_chunk) = c.count() else {
                return Ok(None);
            };
            for _ in 0..n_chunk {
                let (Some(b), Some(e)) = (c.u64(), c.u64()) else {
                    return Ok(None);
                };
                if bin != meta_bin {
                    chunks.push((b, e));
                }
            }
        }
        if !csi {
            let Some(n_intv) = c.count() else {
                return Ok(None);
            };
            if c.take(8 * n_intv).is_none() {
                return Ok(None);
            }
        }
    }
    // Valid for the data file: every chunk ends inside it, and the starts
    // checked begin a BGZF block.
    let Ok(mut d) = fs::File::open(data) else {
        return Ok(None);
    };
    let size = d.metadata()?.len();
    if chunks
        .iter()
        .any(|&(b, e)| b > e || (b >> 16) > size || (e >> 16) > size)
    {
        return Ok(None);
    }
    let mut starts: Vec<u64> = chunks.iter().map(|&(b, _)| b >> 16).collect();
    starts.sort_unstable();
    starts.dedup();
    let step = starts.len().div_ceil(32).max(1);
    for &off in starts.iter().step_by(step) {
        let mut b = [0u8; 18];
        d.seek(SeekFrom::Start(off))?;
        if read_full(&mut d, &mut b)? < 18 || b[..4] != [0x1f, 0x8b, 8, 4] || b[12..14] != *b"BC" {
            return Ok(None);
        }
    }
    drop_cache(&d);
    Ok(Some(hex(h)))
}

/// A path's bytes as one TSV field: valid UTF-8 as it is, except `\`, tab,
/// newline, carriage return and other control bytes; any byte that is not
/// UTF-8 as `\xHH`.
fn escape(bytes: &[u8]) -> String {
    let mut s = String::with_capacity(bytes.len());
    for chunk in bytes.utf8_chunks() {
        for c in chunk.valid().chars() {
            match c {
                '\\' => s.push_str("\\\\"),
                '\t' => s.push_str("\\t"),
                '\n' => s.push_str("\\n"),
                '\r' => s.push_str("\\r"),
                c if (c as u32) < 0x20 || c == '\x7f' => {
                    let _ = write!(s, "\\x{:02x}", c as u32);
                }
                c => s.push(c),
            }
        }
        for b in chunk.invalid() {
            let _ = write!(s, "\\x{b:02x}");
        }
    }
    s
}

/// The bytes an escaped path stands for (the inverse of [`escape`]).
fn unescape(s: &str) -> Option<Vec<u8>> {
    let mut out = Vec::with_capacity(s.len());
    let b = s.as_bytes();
    let mut i = 0;
    while i < b.len() {
        if b[i] != b'\\' {
            out.push(b[i]);
            i += 1;
            continue;
        }
        match b.get(i + 1)? {
            b'\\' => out.push(b'\\'),
            b't' => out.push(b'\t'),
            b'n' => out.push(b'\n'),
            b'r' => out.push(b'\r'),
            b'x' => {
                out.push(u8::from_str_radix(s.get(i + 2..i + 4)?, 16).ok()?);
                i += 2;
            }
            _ => return None,
        }
        i += 2;
    }
    Some(out)
}

// ---------------------------------------------------------------- compare

#[derive(Clone, PartialEq, Eq)]
struct Entry {
    kind: String,
    size: String,
    sha: String,
    mode: String,
    digests: BTreeMap<String, String>,
}

fn read_manifest(path: &str) -> Result<BTreeMap<String, Entry>, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let mut lines = text.lines();
    if lines.next() != Some(MAGIC) {
        return Err(format!(
            "{path} is not an install manifest (no {MAGIC} line)"
        ));
    }
    let mut rows = lines.filter(|l| !l.starts_with('#'));
    if rows.next() != Some(COLUMNS) {
        return Err(format!("{path}: the column line is not {COLUMNS:?}"));
    }
    let mut map = BTreeMap::new();
    for (i, line) in rows.enumerate() {
        let f: Vec<&str> = line.split('\t').collect();
        let [p, kind, size, sha, mode, digests] = f[..] else {
            return Err(format!(
                "{path}: row {} has {} columns, not 6",
                i + 1,
                f.len()
            ));
        };
        if !matches!(kind, "f" | "d" | "l" | "o") {
            return Err(format!("{path}: {p} has type {kind:?}"));
        }
        let mut d = BTreeMap::new();
        if digests != "-" {
            for kv in digests.split(',') {
                let (k, v) = kv
                    .split_once('=')
                    .ok_or_else(|| format!("{path}: {p}: digest {kv:?}"))?;
                d.insert(k.to_owned(), v.to_owned());
            }
        }
        let e = Entry {
            kind: kind.to_owned(),
            size: size.to_owned(),
            sha: sha.to_owned(),
            mode: mode.to_owned(),
            digests: d,
        };
        if map.insert(p.to_owned(), e).is_some() {
            return Err(format!("{path}: {p} is listed twice"));
        }
    }
    Ok(map)
}

const MODES: [&str; 6] = ["gzip", "vcf", "index", "jellyfish", "numeric", "skip"];

struct Allow {
    pattern: String,
    /// The compare mode for content; None compares the bytes.
    content: Option<String>,
    perm: bool,
    reason: String,
}

fn is_stamp(path: &str) -> bool {
    let name = path.rsplit('/').next().unwrap_or(path);
    matches!(name, "DONE" | "INSTALLED" | "RELEASE" | "version.txt")
        || [".stamp", ".installed", ".release"]
            .iter()
            .any(|s| name.ends_with(s))
}

fn read_allow(path: &str) -> Result<Vec<Allow>, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let mut rows = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty());
    if rows.next() != Some(ALLOW_COLUMNS) {
        return Err(format!("{path}: the first row is not {ALLOW_COLUMNS:?}"));
    }
    let mut out = Vec::new();
    for line in rows {
        let f: Vec<&str> = line.split('\t').collect();
        let [pattern, mode, component, reason, cannot] = f[..] else {
            return Err(format!("{path}: {line:?} has {} columns, not 5", f.len()));
        };
        if pattern.is_empty() || component.trim().is_empty() || reason.trim().is_empty() {
            return Err(format!(
                "{path}: {pattern}: a row needs a path, a component and a reason"
            ));
        }
        let last = pattern.rsplit('/').next().unwrap_or(pattern);
        if let Some(s) = STAMP_SAMPLES
            .iter()
            .find(|s| glob(last.as_bytes(), s.as_bytes()))
        {
            return Err(format!(
                "{path}: {pattern} can match a stamp ({s}); stamps are compared byte for byte, never allowed"
            ));
        }
        let mut a = Allow {
            pattern: pattern.to_owned(),
            content: None,
            perm: false,
            reason: reason.to_owned(),
        };
        for k in mode.split(',') {
            match k {
                "perm" => a.perm = true,
                k if MODES.contains(&k) && a.content.is_none() => a.content = Some(k.to_owned()),
                _ => {
                    return Err(format!(
                        "{path}: {pattern}: mode {mode:?}: one of {} with an optional perm",
                        MODES.join(", ")
                    ));
                }
            }
        }
        let skip = a.content.as_deref() == Some("skip");
        if skip == (cannot.trim().is_empty() || cannot == "-") {
            return Err(format!(
                "{path}: {pattern}: cannot_compare must say why no mode can compare a skipped file, and be - for any other row"
            ));
        }
        out.push(a);
    }
    Ok(out)
}

/// `*` and `?` within one segment, `**` across segments; the whole path.
fn glob(p: &[u8], s: &[u8]) -> bool {
    match p {
        [] => s.is_empty(),
        [b'*', b'*', rest @ ..] => (0..=s.len()).any(|i| glob(rest, &s[i..])),
        [b'*', rest @ ..] => (0..=s.len())
            .take_while(|&i| i == 0 || s[i - 1] != b'/')
            .any(|i| glob(rest, &s[i..])),
        [b'?', rest @ ..] => matches!(s, [c, tail @ ..] if *c != b'/' && glob(rest, tail)),
        [c, rest @ ..] => matches!(s, [d, tail @ ..] if d == c && glob(rest, tail)),
    }
}

fn row_for<'a>(allow: &'a [Allow], path: &str) -> Option<&'a Allow> {
    if is_stamp(path) {
        return None;
    }
    allow
        .iter()
        .find(|r| glob(r.pattern.as_bytes(), path.as_bytes()))
}

enum Verdict {
    Same,
    Allowed(String),
    Differs(String),
}

/// Two files of one path under the path's row.
fn content(ea: &Entry, eb: &Entry, row: Option<&Allow>) -> Verdict {
    let short = |h: &str| h.chars().take(12).collect::<String>();
    if ea.kind == "l" {
        return if ea.sha == eb.sha {
            Verdict::Same
        } else {
            Verdict::Differs(format!("{} -> {}", ea.sha, eb.sha))
        };
    }
    if ea.size == eb.size && ea.sha == eb.sha {
        // Equal bytes, but an index must still fit its data file.
        if row.and_then(|r| r.content.as_deref()) == Some("index")
            && [ea, eb]
                .iter()
                .any(|e| e.digests.get("index").is_none_or(|v| v == "invalid"))
        {
            return Verdict::Differs("index: the index does not fit its data file".into());
        }
        return Verdict::Same;
    }
    let bytes = format!(
        "size {} -> {}, sha256 {} -> {}",
        ea.size,
        eb.size,
        short(&ea.sha),
        short(&eb.sha)
    );
    let Some(mode) = row.and_then(|r| r.content.as_deref()) else {
        return Verdict::Differs(bytes);
    };
    if mode == "skip" {
        return Verdict::Allowed(format!("skip: {bytes}"));
    }
    match (ea.digests.get(mode), eb.digests.get(mode)) {
        (Some(a), Some(b)) if a == b && !matches!(a.as_str(), "bad" | "invalid") => {
            Verdict::Allowed(format!("{mode} equal: {bytes}"))
        }
        (Some(a), Some(b)) => {
            Verdict::Differs(format!("{mode} {} -> {}: {bytes}", short(a), short(b)))
        }
        (a, b) => Verdict::Differs(format!(
            "{mode}: {} has no {mode} digest: {bytes}",
            if a.is_none() {
                "A"
            } else if b.is_none() {
                "B"
            } else {
                "?"
            }
        )),
    }
}

fn compare_cli(args: &[String], out: &mut dyn Write) -> Result<u8, String> {
    let (a, b, allow) = match args {
        [a, b] => (a, b, None),
        [a, b, flag, al] if flag == "--allow" => (a, b, Some(al)),
        _ => return Err(String::new()),
    };
    let (ma, mb) = (read_manifest(a)?, read_manifest(b)?);
    let allow = match allow {
        Some(p) => read_allow(p)?,
        None => Vec::new(),
    };
    // (path, kind, detail, Some(reason) when allowed)
    let mut lines: Vec<(&str, &str, String, Option<String>)> = Vec::new();
    for (p, ea) in &ma {
        let Some(eb) = mb.get(p) else {
            lines.push((p, "removed", format!("{} in A only", ea.kind), None));
            continue;
        };
        if ea.kind != eb.kind {
            lines.push((p, "type", format!("{} -> {}", ea.kind, eb.kind), None));
            continue;
        }
        let row = row_for(&allow, p);
        match content(ea, eb, row) {
            Verdict::Same => {}
            Verdict::Allowed(d) => lines.push((p, "content", d, row.map(|r| r.reason.clone()))),
            Verdict::Differs(d) => lines.push((p, "content", d, None)),
        }
        if ea.mode != eb.mode {
            let reason = row.filter(|r| r.perm).map(|r| r.reason.clone());
            lines.push((p, "mode", format!("{} -> {}", ea.mode, eb.mode), reason));
        }
    }
    for (p, eb) in &mb {
        if !ma.contains_key(p) {
            lines.push((p, "added", format!("{} in B only", eb.kind), None));
        }
    }
    lines.sort_by(|x, y| x.0.cmp(y.0).then(x.1.cmp(y.1)));
    let mut text = String::new();
    let (mut bad, mut allowed) = (0usize, 0usize);
    for (p, kind, detail, reason) in &lines {
        match reason {
            Some(r) => {
                allowed += 1;
                let _ = writeln!(text, "ALLOWED\t{kind}\t{p}\t{detail}\t{r}");
            }
            None => {
                bad += 1;
                let _ = writeln!(text, "DIFF\t{kind}\t{p}\t{detail}");
            }
        }
    }
    let _ = writeln!(
        text,
        "# {} entries in A, {} in B: {bad} differences, {allowed} allowed",
        ma.len(),
        mb.len()
    );
    out.write_all(text.as_bytes())
        .map_err(|e| format!("stdout: {e}"))?;
    Ok(u8::from(bad > 0))
}

fn control_cli(args: &[String], out: &mut dyn Write) -> Result<u8, String> {
    let (m, allow) = match args {
        [m] => (m, None),
        [m, flag, al] if flag == "--allow" => (m, Some(al)),
        _ => return Err(String::new()),
    };
    let m = read_manifest(m)?;
    let allow = match allow {
        Some(p) => read_allow(p)?,
        None => Vec::new(),
    };
    let pick = m
        .iter()
        .filter(|(p, e)| {
            e.kind == "f"
                && e.size != "0"
                && row_for(&allow, p)
                    .and_then(|r| r.content.as_ref())
                    .is_none()
        })
        .min_by_key(|(p, e)| (e.size.parse::<u64>().unwrap_or(u64::MAX), p.as_str()));
    match pick {
        Some((p, _)) => {
            writeln!(out, "{p}").map_err(|e| format!("stdout: {e}"))?;
            Ok(0)
        }
        None => Ok(1),
    }
}

fn flip_cli(args: &[String], out: &mut dyn Write) -> Result<u8, String> {
    let [dir, path] = args else {
        return Err(String::new());
    };
    let rel = unescape(path).ok_or_else(|| format!("{path} is not an escaped manifest path"))?;
    let full = Path::new(dir).join(std::ffi::OsStr::from_bytes(&rel));
    let e = io_err(&full);
    let mut f = fs::OpenOptions::new()
        .read(true)
        .write(true)
        .open(&full)
        .map_err(&e)?;
    let len = f.metadata().map_err(&e)?.len();
    if len == 0 {
        return Err(format!("{} is empty", shown(&full)));
    }
    let at = len / 2;
    let mut b = [0u8; 1];
    f.seek(SeekFrom::Start(at)).map_err(&e)?;
    f.read_exact(&mut b).map_err(&e)?;
    b[0] ^= 1;
    f.seek(SeekFrom::Start(at)).map_err(&e)?;
    f.write_all(&b).map_err(&e)?;
    writeln!(out, "flipped one bit at byte {at} of {path}").map_err(|e| format!("stdout: {e}"))?;
    Ok(0)
}

// ---------------------------------------------------------------- logs

static CURL_ROW: LazyLock<Regex> = LazyLock::new(|| {
    // A progress row of curl: a percentage, then sizes, speeds and times
    // only. curl 8 leaves a time it does not know blank, so the number of
    // fields varies.
    Regex::new(r"^\s*\d{1,3}(\s+(\d+(\.\d+)?[kMGTP]?|[\d-]*:[\d:-]*)){7,11}\s*$")
        .expect("a valid pattern")
});
static PULL: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(
        r"^(Unable to find image '[^']+' locally$|[0-9a-f]{12}: (Pulling fs layer|Waiting|Downloading|Verifying Checksum|Download complete|Extracting|Pull complete|Already exists)|Digest: sha256:|Status: (Downloaded newer image|Image is up to date) for |\S+: Pulling from )",
    )
    .expect("a valid pattern")
});
static MASKS: LazyLock<Vec<(Regex, &'static str)>> = LazyLock::new(|| {
    [
        (
            r"\d{4}-\d{2}-\d{2}([T ]\d{2}:\d{2}(:\d{2}(\.\d+)?)?(Z|[+-]\d{2}:?\d{2})?)?",
            "<time>",
        ),
        (r"\b\d{1,2}:\d{2}:\d{2}\b", "<time>"),
        (r"\b\d+(\.\d+)?\s?[kMGT]i?B/s\b", "<speed>"),
        (
            r"\b\d+(\.\d+)?\s?(ms|s|sec|secs|seconds|min|mins|minutes|h|hours)\b",
            "<duration>",
        ),
    ]
    .into_iter()
    .map(|(r, t)| (Regex::new(r).expect("a valid pattern"), t))
    .collect()
});

/// A setup log's lines with progress dropped and times masked.
fn masked_log(text: &str) -> Vec<String> {
    text.split(['\n', '\r'])
        .map(str::trim_end)
        .filter(|l| {
            !l.is_empty()
                && !CURL_ROW.is_match(l)
                && !l.trim_start().starts_with("% Total")
                && !l.trim_start().starts_with("Dload")
                && !PULL.is_match(l)
        })
        .map(|l| {
            MASKS.iter().fold(l.to_owned(), |s, (re, to)| {
                re.replace_all(&s, *to).into_owned()
            })
        })
        .collect()
}

/// The missing list of a setup run: its summary and each item with its
/// retry command, in order.
fn missing_list(lines: &[String]) -> Vec<&str> {
    lines
        .iter()
        .map(String::as_str)
        .filter(|l| {
            l.starts_with("Missing (")
                || l.starts_with("  - ")
                || l.trim_start().starts_with("retry: ")
                || l.contains("Setup complete!")
                || l.contains("Setup finished:")
                || l.contains("exits with status")
        })
        .collect()
}

fn logs_cli(args: &[String], out: &mut dyn Write) -> Result<u8, String> {
    let [a, b] = args else {
        return Err(String::new());
    };
    let read = |p: &str| {
        fs::read(p)
            .map(|v| String::from_utf8_lossy(&v).into_owned())
            .map_err(|e| format!("{p}: {e}"))
    };
    let (la, lb) = (masked_log(&read(a)?), masked_log(&read(b)?));
    let mut text = String::new();
    let mut bad = 0usize;
    let (ma, mb) = (missing_list(&la), missing_list(&lb));
    if ma != mb {
        bad += 1;
        let _ = writeln!(
            text,
            "DIFF\tmissing\tA: {}\tB: {}",
            ma.join(" | "),
            mb.join(" | ")
        );
    }
    // Lines in one log and not the other, counted, in the order they come.
    let mut count: HashMap<&str, i64> = HashMap::new();
    for l in &la {
        *count.entry(l).or_default() += 1;
    }
    for l in &lb {
        *count.entry(l).or_default() -= 1;
    }
    for (side, lines, sign) in [("A", &la, 1), ("B", &lb, -1)] {
        for l in lines.iter() {
            let c = count.get_mut(l.as_str()).expect("counted");
            if *c * sign > 0 {
                *c -= sign;
                bad += 1;
                let _ = writeln!(text, "DIFF\tlog\t{side} only\t{l}");
            }
        }
    }
    let _ = writeln!(
        text,
        "# {} lines in A, {} in B after masking: {bad} differences",
        la.len(),
        lb.len()
    );
    out.write_all(text.as_bytes())
        .map_err(|e| format!("stdout: {e}"))?;
    Ok(u8::from(bad > 0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn glob_segments() {
        assert!(glob(b"a/*.lock", b"a/x.lock"));
        assert!(!glob(b"a/*.lock", b"a/b/x.lock"));
        assert!(glob(b"a/**.lock", b"a/b/x.lock"));
        assert!(glob(
            b"**/__pycache__/*.pyc",
            b"tools/c/x/__pycache__/y.pyc"
        ));
        assert!(glob(b"a/?.txt", b"a/x.txt"));
        assert!(!glob(b"a/?.txt", b"a//.txt"));
        assert!(!glob(b"a/x", b"a/x/y"));
    }

    #[test]
    fn escapes_round_trip() {
        let raw = b"a\tb\nc\\d\x01\xffe\r";
        assert_eq!(escape(raw), "a\\tb\\nc\\\\d\\x01\\xffe\\r");
        assert_eq!(unescape(&escape(raw)).unwrap(), raw);
        assert_eq!(escape("é".as_bytes()), "é");
    }

    #[test]
    fn stamps() {
        for s in STAMP_SAMPLES {
            assert!(is_stamp(&format!("a/{s}")), "{s}");
        }
        assert!(!is_stamp("a/catalog.json"));
    }

    #[test]
    fn numeric_rounds_decimals_only() {
        let a = numeric_digest(b"id1 Pop 0.123456 -0.000001234\n12345678 x\n").unwrap();
        let b = numeric_digest(b"id1  Pop 0.1234561 -0.0000012341\n12345678 x\n").unwrap();
        assert_eq!(a, b);
        let c = numeric_digest(b"id1 Pop 0.123456 -0.000001234\n12345679 x\n").unwrap();
        assert_ne!(a, c, "an integer word is compared exactly");
        let d = numeric_digest(b"id1 Pop 0.1236 -0.000001234\n12345678 x\n").unwrap();
        assert_ne!(a, d, "a change in the fourth digit counts");
        assert!(numeric_digest(b"no numbers 12 here\n").is_none());
    }

    #[test]
    fn log_masks() {
        let l = masked_log(
            "2026-11-26T10:00:01Z start\r  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current\n                                 Dload  Upload   Total   Spent    Left  Speed\n100  1234  100  1234    0     0   5000      0 --:--:-- --:--:-- --:--:--  5012\ndone in 12 s at 3.4 MB/s\n",
        );
        assert_eq!(l, ["<time> start", "done in <duration> at <speed>"]);
        // curl 8 leaves a time it does not know blank; docker run says when
        // it pulls an image first, which only the first of two runs does.
        let l = masked_log(
            "100 192.3M 100 192.3M   0      0  9.51M      0   00:20   00:20          6.83M\n100 597.6k 100 597.6k   0      0 689.5k      0                              0\nUnable to find image 'staphb/bcftools:1.24' locally\nkept\n",
        );
        assert_eq!(l, ["kept"]);
        let l = masked_log("100 apples\n");
        assert_eq!(l, ["100 apples"]);
    }
}
