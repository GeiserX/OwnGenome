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
//! kept or not by choice) and each `--exclude` path with all under it. A
//! file whose name ends in `.log` (a tool log, with times and hosts in it)
//! is listed for its presence only: `f`, then `-` for its size, sha256,
//! mode and digests. In a name, a tab, a newline, a carriage return, a
//! backslash, any other control byte and any byte that is not UTF-8 are
//! written as `\t`, `\n`, `\r`, `\\` and `\xHH`.
//! A file that changes while it is read (its size, modification time or
//! identity) stops the manifest with an error: its row would describe
//! neither version.
//!
//! `digests` holds, for a file of a kind a compare mode reads, the sha256 of
//! what that mode compares (`kind=hex`, comma separated, `-` for none). The
//! manifest carries them so a tree can be compared after it is gone (the
//! reference manifest of a full install). The `gzip` and `vcf` digests
//! start with the container, `bgzf:`, `gzip:` or (a plain VCF) `plain:`:
//! tabix and every region read need BGZF, so plain gzip is not the same
//! file. BGZF means every gzip member is a BGZF block (htslib's header
//! check), each starting where the one before ends, up to the end, and
//! each member ending where its block header says and decoding to at most
//! 64 KiB.
//!
//! - `gzip`: the decompressed content of a gzip or BGZF file (every
//!   member), so the header's time and OS byte do not count. `bad` when it
//!   does not decompress.
//! - `vcf`: a VCF (plain or gzip) with the value of its `##fileDate=` line
//!   and the `Date=` part of each `##bcftools_<cmd>Command=` line masked.
//!   Everything else, the bcftools version and command lines included, is
//!   compared. A header line longer than 16 MiB (its newline included) is
//!   compared as it is.
//! - `index`: a `.tbi` or `.csi` whose data file beside it is a regular file,
//!   not a link: its header (format, columns, the CSI shift and depth), its
//!   sequence names, and per sequence what does not move with the data file's
//!   bytes: the record counts of htslib's pseudo-bin (mapped and unmapped)
//!   and, in a `.tbi`, the length of the linear index (the last 16 kb window
//!   holding a record); then the count of records with no position, when the
//!   index has it. `invalid` when it does not fit the data file, when bytes
//!   other than that count follow the last sequence, or when the index is
//!   over 64 MiB. Fit means every chunk ends inside the data file, and every
//!   chunk start and end points into a BGZF block that decodes whole, to at
//!   most 64 KiB, CRC and length matching: no further than its decoded
//!   length, or at offset 0 at the end of the file. Offsets and bin numbers
//!   are left out: they move with any byte of the data file, which the data
//!   file's own row compares. htslib merges a bin into its parent when its
//!   chunks span less than 64 KiB of the compressed file, so even the set of
//!   bins moves with the compressed sizes.
//! - `jellyfish`: a jellyfish file (9 digits, then its JSON header of at
//!   most 16 MiB): the header with the values of `hostname`, `pwd` and
//!   `time` masked (the host, the folder and the time of the count), then
//!   the counts. Jellyfish 2 pads its header with NUL bytes to its
//!   alignment and counts them in the 9-digit length, so the counts start
//!   right after that length; the trailing NULs are left out of the
//!   header, since a longer host name needs fewer of them. The rest of the
//!   header (the layout, the hash matrix, the command line) is compared:
//!   it decides how the counts are read.
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
//! printed as `ALLOWED` with the mode and the row's reason. A manifest has
//! no comments: after its first two lines every line is a row, and a path
//! may start with `#`. In the allowlist, lines starting with `#` are
//! comments.
//!
//! The allowlist (`tests/setup/install-manifest-allow.tsv`) has the columns
//! `path mode component reason cannot_compare`. `path` is a pattern over the
//! whole relative path: `*` and `?` match within one segment, `**` across
//! segments; the first matching row counts. `mode` is one of `gzip`, `vcf`,
//! `index`, `jellyfish`, `numeric` or `skip` (content not compared), and may
//! add `perm` (`gzip,perm`, or `perm` alone): the permission bits may
//! differ. `component` and `reason` must not be empty; `cannot_compare` says
//! why no other mode can compare the file: required for `skip`, and exactly
//! `-` for any other mode. A manifest row that names one digest kind twice
//! is refused. A stamp (`DONE`, `INSTALLED`, `RELEASE`, `version.txt`,
//! `pins.tsv`, `*.stamp`, `*.installed`, `*.release`, `*.sha256`) is
//! always compared byte for byte: its bytes decide whether the next setup
//! run reinstalls or reuses, so a row whose last segment can match any
//! stamp name is refused, and no row applies to a stamp.
//!
//! # For tests/setup/ab.sh
//!
//! `--control B.tsv [--allow ALLOW.tsv] [--not-in SEED.tsv]` prints the
//! smallest non-empty file compared byte for byte, leaving out a file whose
//! row SEED.tsv holds as it is (a file the seed put there, not the
//! component); `--flip DIR PATH` flips one bit of the byte in
//! the middle of DIR/PATH (twice restores it); PATH must name a regular file
//! inside DIR, reached through folders only. `--compare-logs A B` compares
//! two setup logs with times, durations, speeds and download progress
//! masked and download retries left out (fetch's `Download attempt i/n
//! failed: URL` and `Checksum download attempt i/n failed: URL`, and curl's
//! `curl: (N) ...`; nothing else), and the missing list (`Missing (...)`, `  - item`, `retry:`)
//! line by line, unmasked.
//!
//! Exit 0 when nothing differs outside the allowlist, 1 when something
//! does, 2 when an argument, a manifest or the allowlist is wrong or a file
//! cannot be read.

use flate2::bufread::GzDecoder;
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

pub(crate) const USAGE: &str = "usage: og dev install-manifest DIR [--exclude PATH]...
       og dev install-manifest --compare A.tsv B.tsv [--allow ALLOW.tsv]
       og dev install-manifest --control B.tsv [--allow ALLOW.tsv] [--not-in SEED.tsv]
       og dev install-manifest --flip DIR PATH
       og dev install-manifest --compare-logs A.log B.log
";
const MAGIC: &str = "#og-install-manifest 1";
const COLUMNS: &str = "path\ttype\tsize\tsha256\tmode\tdigests";
const ALLOW_COLUMNS: &str = "path\tmode\tcomponent\treason\tcannot_compare";
/// Files up to this size are read once into memory; larger ones are
/// streamed. Decompressed content is always streamed.
const SMALL: u64 = 16 << 20;
/// The longest VCF header line masked, and the longest jellyfish header
/// read: past it a line is compared as it is, and a file is no jellyfish
/// file, so memory stays bounded on any input.
const HEADER_MAX: usize = 16 << 20;
/// The largest index read, compressed and decompressed; the real ones are a
/// few MB. A larger one is `invalid`, which never compares equal.
const INDEX_MAX: u64 = 64 << 20;
/// The most a BGZF block decodes to.
const BGZF_MAX: u64 = 64 << 10;
/// Stamp names: whole names, and endings. `pins.tsv` (the fine panel's
/// pins) and `*.sha256` (a saved region's checksum) decide a reuse.
const STAMP_NAMES: [&str; 5] = ["DONE", "INSTALLED", "RELEASE", "version.txt", "pins.tsv"];
const STAMP_ENDINGS: [&str; 4] = [".stamp", ".installed", ".release", ".sha256"];

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
                || (rel.is_empty() && ft.is_dir() && name.as_bytes() == b"downloads")
            {
                self.skipped += 1;
                continue;
            }
            let mode = format!("{:04o}", meta.permissions().mode() & 0o7777);
            let row = if ft.is_file() && name.as_bytes().ends_with(b".log") {
                // A log's bytes hold times and hosts: its presence only.
                "f\t-\t-\t-\t-".to_owned()
            } else if ft.is_symlink() {
                let target = fs::read_link(&full).map_err(io_err(&full))?;
                format!("l\t-\t->{}\t-\t-", escape(target.as_os_str().as_bytes()))
            } else if ft.is_dir() {
                format!("d\t-\t-\t{mode}\t-")
            } else if ft.is_file() {
                let (sha, digests) = digest_file(&full, &meta)?;
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

/// Whether two observations of a file are the same file with the same size
/// and modification time.
fn same_file(a: &fs::Metadata, b: &fs::Metadata) -> bool {
    use std::os::unix::fs::MetadataExt;
    a.dev() == b.dev()
        && a.ino() == b.ino()
        && a.len() == b.len()
        && a.mtime() == b.mtime()
        && a.mtime_nsec() == b.mtime_nsec()
}

/// The sha256 of a file and its `digests` column. `meta` is what the walk
/// saw: a file that is not that file, or changes while it is read, is an
/// error, as its row would describe neither version.
fn digest_file(path: &Path, meta: &fs::Metadata) -> Result<(String, String), String> {
    let e = io_err(path);
    let len = meta.len();
    let mut f = fs::File::open(path).map_err(&e)?;
    let changed = || format!("{}: the file changed while it was read", shown(path));
    if !same_file(meta, &f.metadata().map_err(&e)?) {
        return Err(changed());
    }
    let read;
    let mut head = [0u8; 64];
    let n = read_full(&mut f, &mut head).map_err(&e)?;
    let head = &head[..n];
    f.seek(SeekFrom::Start(0)).map_err(&e)?;
    let gz = head.starts_with(&[0x1f, 0x8b]);
    let name = path.file_name().map(|n| n.as_bytes()).unwrap_or_default();
    let mut digests: Vec<(&str, String)> = Vec::new();
    let sha;
    if len <= SMALL {
        let buf = read_capped(&f, len).map_err(&e)?;
        read = buf.len() as u64;
        sha = sha_of(&buf);
        if gz {
            let mut hr = Hashing {
                inner: &buf[..],
                h: Sha256::new(),
                n: 0,
                walk: BgzfWalk::default(),
            };
            let d = gzip_stream_digests(Members::new(&mut hr));
            std::io::copy(&mut hr, &mut std::io::sink()).map_err(&e)?;
            digests.extend(with_container(hr.walk.kind(), d));
        } else if buf.starts_with(b"##fileformat=VCF") {
            let d = vcf_digest(&buf[..]).map_err(&e)?;
            digests.push(("vcf", format!("plain:{d}")));
        } else if let Some(d) = numeric_digest(&buf) {
            digests.push(("numeric", d));
        }
        if let Some(d) = jellyfish_digest(&buf[..]).map_err(&e)? {
            digests.push(("jellyfish", d));
        }
    } else {
        let (s, d, n) = stream_digests(BufReader::with_capacity(1 << 20, (&f).take(len + 1)), head)
            .map_err(&e)?;
        sha = s;
        read = n;
        digests.extend(d);
    }
    // An index is parsed from the same open file, so a file put in its
    // place meanwhile cannot lend it its digest; the check after it covers
    // both reads.
    if name.ends_with(b".tbi") || name.ends_with(b".csi") {
        digests.push(("index", index_digest(path, &f)));
    }
    if read != len || !same_file(meta, &f.metadata().map_err(&e)?) {
        return Err(changed());
    }
    drop_cache(&f);
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

/// A file of `len` bytes read whole, but never more than `len + 1` bytes:
/// one byte more shows the file grew, at once, and memory stays at its
/// size whatever a writer appends meanwhile.
fn read_capped(r: impl Read, len: u64) -> std::io::Result<Vec<u8>> {
    let mut buf = Vec::with_capacity(len as usize + 1);
    r.take(len + 1).read_to_end(&mut buf)?;
    Ok(buf)
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

/// A reader that hashes every byte it hands on, so the one read a decoder
/// makes of a file also gives the file's sha256 and its container.
struct Hashing<R> {
    inner: R,
    h: Sha256,
    /// Bytes read.
    n: u64,
    walk: BgzfWalk,
}

impl<R: Read> Read for Hashing<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let k = self.inner.read(buf)?;
        self.h.update(&buf[..k]);
        self.n += k as u64;
        self.walk.update(&buf[..k]);
        Ok(k)
    }
}

/// The BSIZE of a BGZF block header, by htslib's own check (gzip, deflate,
/// FEXTRA, XLEN 6 and the BC subfield of length 2), or None.
fn bgzf_bsize(h: &[u8]) -> Option<u64> {
    (h.len() >= 18
        && h[..3] == [0x1f, 0x8b, 8]
        && h[3] & 4 != 0
        && h[10..12] == [6, 0]
        && h[12..14] == *b"BC"
        && h[14..16] == [2, 0])
    .then(|| u64::from(u16::from_le_bytes([h[16], h[17]])))
}

/// Whether a gzip file is BGZF, from its raw bytes in chunks of any size
/// and the member ends the decoder reports: every member a BGZF block,
/// each one starting where the one before ends (BSIZE + 1 bytes on), each
/// member ending where its block header says, and the last one ending the
/// file.
#[derive(Default)]
struct BgzfWalk {
    /// Bytes seen.
    pos: u64,
    /// Where the next block starts.
    next: u64,
    /// The header of the block at `next`, while it comes in.
    head: Vec<u8>,
    /// Block ends the headers declared that no member has reached yet.
    ends: std::collections::VecDeque<u64>,
    broken: bool,
}

impl BgzfWalk {
    fn update(&mut self, mut data: &[u8]) {
        while !data.is_empty() && !self.broken {
            if self.pos < self.next {
                let k =
                    usize::try_from(self.next - self.pos).map_or(data.len(), |g| g.min(data.len()));
                self.pos += k as u64;
                data = &data[k..];
                continue;
            }
            let k = (18 - self.head.len()).min(data.len());
            self.head.extend_from_slice(&data[..k]);
            self.pos += k as u64;
            data = &data[k..];
            if self.head.len() == 18 {
                match bgzf_bsize(&self.head) {
                    // A block holds its 18-byte header and 8-byte trailer.
                    Some(b) if b + 1 >= 26 => {
                        self.next += b + 1;
                        self.ends.push_back(self.next);
                    }
                    _ => self.broken = true,
                }
                self.head.clear();
            }
        }
    }

    /// The decoder read a whole member of `decoded` bytes, ending at raw
    /// offset `at`. The walk has always read that member's header by then.
    /// A BGZF block decodes to at most 64 KiB.
    fn member_end(&mut self, at: u64, decoded: u64) {
        if self.ends.pop_front() != Some(at) || decoded > BGZF_MAX {
            self.broken = true;
        }
    }

    fn kind(&self) -> &'static str {
        if !self.broken
            && self.head.is_empty()
            && self.ends.is_empty()
            && self.next > 0
            && self.pos == self.next
        {
            "bgzf"
        } else {
            "gzip"
        }
    }
}

/// The digests of a gzip file with its container in front, `bad` as it is.
fn with_container(kind: &str, digests: Vec<(&'static str, String)>) -> Vec<(&'static str, String)> {
    digests
        .into_iter()
        .map(|(k, v)| (k, if v == "bad" { v } else { format!("{kind}:{v}") }))
        .collect()
}

/// A file's digests: each compare mode's kind and value.
type Digests = Vec<(&'static str, String)>;

/// The sha256, the mode digests and the length of a large file, from one
/// read of `r`: a full install is hundreds of GB, and a second read doubles
/// the time.
fn stream_digests(r: impl Read, head: &[u8]) -> std::io::Result<(String, Digests, u64)> {
    let mut hr = Hashing {
        inner: r,
        h: Sha256::new(),
        n: 0,
        walk: BgzfWalk::default(),
    };
    let mut digests = Vec::new();
    let gz = head.starts_with(&[0x1f, 0x8b]);
    if gz {
        digests.extend(gzip_stream_digests(Members::new(&mut hr)));
    } else if head.starts_with(b"##fileformat=VCF") {
        digests.push(("vcf", format!("plain:{}", vcf_digest(&mut hr)?)));
    } else if is_jellyfish(head)
        && let Some(d) = jellyfish_digest(BufReader::with_capacity(1 << 20, &mut hr))?
    {
        digests.push(("jellyfish", d));
    }
    // Whatever a decoder left (bytes after the gzip stream, the end of a
    // file it stopped reading) still counts in the sha256.
    std::io::copy(&mut hr, &mut std::io::sink())?;
    if gz {
        digests = with_container(hr.walk.kind(), digests);
    }
    Ok((hex(hr.h), digests, hr.n))
}

/// The members of a gzip stream decoded one after another, as
/// `MultiGzDecoder` reads them, telling the BGZF walk where each one ends.
struct Members<'a, R: Read> {
    dec: Option<GzDecoder<BufReader<&'a mut Hashing<R>>>>,
    /// Bytes the current member decoded to so far.
    decoded: u64,
}

impl<'a, R: Read> Members<'a, R> {
    fn new(hr: &'a mut Hashing<R>) -> Self {
        Members {
            dec: Some(GzDecoder::new(BufReader::with_capacity(1 << 20, hr))),
            decoded: 0,
        }
    }
}

impl<R: Read> Read for Members<'_, R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        while let Some(dec) = self.dec.as_mut() {
            let k = dec.read(buf)?;
            if k > 0 || buf.is_empty() {
                self.decoded += k as u64;
                return Ok(k);
            }
            // The member and its trailer are read: what the buffer still
            // holds comes after it. A broken walk stopped counting, and the
            // file is plain gzip whatever follows.
            let mut b = self.dec.take().expect("a decoder").into_inner();
            if !b.get_ref().walk.broken {
                let at = b.get_ref().walk.pos - b.buffer().len() as u64;
                b.get_mut().walk.member_end(at, self.decoded);
            }
            self.decoded = 0;
            if !b.fill_buf()?.is_empty() {
                self.dec = Some(GzDecoder::new(b));
            }
        }
        Ok(0)
    }
}

/// The gzip and, for a VCF, the vcf digest of a gzip file, in one
/// streaming pass over its content in 1 MiB chunks.
fn gzip_stream_digests(mut r: impl Read) -> Vec<(&'static str, String)> {
    let mut raw = Sha256::new();
    let mut vcf: Option<VcfHasher> = None;
    let mut buf = vec![0u8; 1 << 20];
    let mut first = true;
    loop {
        let k = match read_full(&mut r, &mut buf) {
            Ok(0) => break,
            Ok(k) => k,
            Err(_) => return vec![("gzip", "bad".into())],
        };
        let chunk = &buf[..k];
        if first {
            first = false;
            if chunk.starts_with(b"##fileformat=VCF") {
                vcf = Some(VcfHasher::default());
            }
        }
        raw.update(chunk);
        if let Some(v) = vcf.as_mut() {
            v.update(chunk);
        }
    }
    let mut d = vec![("gzip", hex(raw))];
    if let Some(v) = vcf {
        d.push(("vcf", v.finish()));
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

/// The vcf digest, fed in chunks of any size. Each `##` header line is
/// held until its end and hashed masked; one longer than `HEADER_MAX`, its
/// newline included, is hashed unmasked, however the chunks fall. From the
/// first line that does not start with `##`, the bytes are hashed as they
/// come.
#[derive(Default)]
struct VcfHasher {
    h: Sha256,
    body: bool,
    line: Vec<u8>,
    /// The line being read passed `HEADER_MAX`: its rest is hashed as it is.
    long: bool,
}

impl VcfHasher {
    fn update(&mut self, mut data: &[u8]) {
        while !data.is_empty() {
            if self.body {
                self.h.update(data);
                return;
            }
            let end = data.iter().position(|&b| b == b'\n').map(|i| i + 1);
            let (part, rest) = data.split_at(end.unwrap_or(data.len()));
            data = rest;
            if self.long {
                self.h.update(part);
            } else {
                self.line.extend_from_slice(part);
                // Past the cap the line is hashed as it is, whether or not
                // this chunk also holds its newline.
                if self.line.len() > HEADER_MAX {
                    self.body = !self.line.starts_with(b"##");
                    self.h.update(&self.line);
                    self.line.clear();
                    self.long = true;
                }
            }
            if end.is_some() {
                if !self.long {
                    self.end_line();
                }
                self.long = false;
            }
        }
    }

    fn end_line(&mut self) {
        if self.line.starts_with(b"##") {
            self.h.update(vcf_mask(&self.line));
        } else {
            self.body = true;
            self.h.update(&self.line);
        }
        self.line.clear();
    }

    fn finish(mut self) -> String {
        if !self.line.is_empty() {
            self.end_line();
        }
        hex(self.h)
    }
}

fn vcf_digest(mut r: impl Read) -> std::io::Result<String> {
    let mut v = VcfHasher::default();
    let mut buf = vec![0u8; 1 << 20];
    loop {
        let k = r.read(&mut buf)?;
        if k == 0 {
            break;
        }
        v.update(&buf[..k]);
    }
    Ok(v.finish())
}

fn is_jellyfish(head: &[u8]) -> bool {
    head.len() > 10 && head[..9].iter().all(u8::is_ascii_digit) && head[9] == b'{'
}

static JF_VOLATILE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r#""(hostname|pwd|time)"\s*:\s*"(?:[^"\\]|\\.)*""#).expect("a valid pattern")
});

/// The digest of a jellyfish file: its JSON header with the host, folder
/// and time masked, then its counts; None when it is not one, or its
/// header is longer than `HEADER_MAX`. Jellyfish 2 pads the header with
/// NUL bytes to its `alignment` inside the 9-digit length, so the counts
/// start right after that length and the NULs are not hashed: a longer
/// host name needs fewer of them.
fn jellyfish_digest(mut r: impl Read) -> std::io::Result<Option<String>> {
    let mut len = [0u8; 10];
    if read_full(&mut r, &mut len)? < 10 || !is_jellyfish(&[&len[..], b"x"].concat()) {
        return Ok(None);
    }
    let n: usize = std::str::from_utf8(&len[..9])
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    if n == 0 || n > HEADER_MAX {
        return Ok(None);
    }
    let mut json = vec![0u8; n - 1];
    if read_full(&mut r, &mut json)? < json.len() {
        return Ok(None);
    }
    json.truncate(json.iter().rposition(|&b| b != 0).map_or(0, |i| i + 1));
    let json = String::from_utf8_lossy(&[b"{".as_slice(), &json].concat()).into_owned();
    let header = JF_VOLATILE.replace_all(&json, r#""$1":"<masked>""#);
    let mut h = Sha256::new();
    h.update(format!("{}\n", header.len()));
    h.update(header.as_bytes());
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
fn index_digest(path: &Path, f: &fs::File) -> String {
    match read_index(path, f) {
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

/// The index digest of the open index file `f` at `path`, read from its
/// start; the data file is `path` less `.tbi` or `.csi`.
fn read_index(path: &Path, mut f: &fs::File) -> std::io::Result<Option<String>> {
    let s = path.as_os_str().as_bytes();
    let data = Path::new(std::ffi::OsStr::from_bytes(&s[..s.len() - 4]));
    if f.metadata()?.len() > INDEX_MAX {
        return Ok(None);
    }
    f.seek(SeekFrom::Start(0))?;
    let mut plain = Vec::new();
    MultiGzDecoder::new(BufReader::new(f))
        .take(INDEX_MAX + 1)
        .read_to_end(&mut plain)?;
    if plain.len() as u64 > INDEX_MAX {
        return Ok(None);
    }
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
        // The pseudo-bin's record counts (mapped, unmapped), when present.
        let mut counts: Option<(u64, u64)> = None;
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
            for i in 0..n_chunk {
                let (Some(b), Some(e)) = (c.u64(), c.u64()) else {
                    return Ok(None);
                };
                if bin != meta_bin {
                    chunks.push((b, e));
                } else if i == 1 {
                    counts = Some((b, e));
                }
            }
        }
        match counts {
            Some((mapped, unmapped)) => {
                h.update([1]);
                h.update(mapped.to_le_bytes());
                h.update(unmapped.to_le_bytes());
            }
            None => h.update([0]),
        }
        if !csi {
            let Some(n_intv) = c.count() else {
                return Ok(None);
            };
            // How many 16 kb windows the linear index spans: where the
            // last record ends, whatever the offsets.
            h.update((n_intv as u64).to_le_bytes());
            if c.take(8 * n_intv).is_none() {
                return Ok(None);
            }
        }
    }
    // The records with no position, which htslib writes after the last
    // sequence: the count, or nothing, and no other bytes.
    match c.0.len() {
        0 => h.update([0]),
        8 => {
            h.update([1]);
            h.update(c.0);
        }
        _ => return Ok(None),
    }
    // Valid for the data file: a regular file, not a link (the manifest
    // compares a link's text, never its target's bytes, and a FIFO would
    // block the read), every chunk ends inside it, at offset 0 when at its
    // end, and every chunk start and end points into a BGZF block.
    if !fs::symlink_metadata(data).is_ok_and(|m| m.is_file()) {
        return Ok(None);
    }
    let Ok(d) = fs::File::open(data) else {
        return Ok(None);
    };
    let size = d.metadata()?.len();
    let outside = |v: u64| (v >> 16) > size || ((v >> 16) == size && v & 0xffff != 0);
    if chunks
        .iter()
        .any(|&(b, e)| b > e || outside(b) || outside(e))
    {
        return Ok(None);
    }
    // Each block once, with the furthest offset into it, in file order.
    let mut blocks: BTreeMap<u64, u64> = BTreeMap::new();
    for v in chunks.iter().flat_map(|&(b, e)| [b, e]) {
        let w = blocks.entry(v >> 16).or_default();
        *w = (*w).max(v & 0xffff);
    }
    for (&at, &within) in &blocks {
        if !in_block(&d, size, at, within)? {
            return Ok(None);
        }
    }
    drop_cache(&d);
    Ok(Some(hex(h)))
}

/// Whether offset `within` of the BGZF block at file offset `at` exists:
/// up to the block's decoded length (its end counts), or 0 at the end of
/// the file. The block must decode as one whole gzip member, its CRC and
/// ISIZE matching what it decodes to.
fn in_block(d: &fs::File, size: u64, at: u64, within: u64) -> std::io::Result<bool> {
    use std::os::unix::fs::FileExt;
    let read = |buf: &mut [u8], off: u64| match d.read_exact_at(buf, off) {
        Ok(()) => Ok(true),
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(false),
        Err(e) => Err(e),
    };
    if at == size {
        return Ok(within == 0);
    }
    let mut b = [0u8; 18];
    if !read(&mut b, at)? {
        return Ok(false);
    }
    let Some(bsize) = bgzf_bsize(&b) else {
        return Ok(false);
    };
    // BSIZE is the block's length less one: the 18-byte header (XLEN 6),
    // the deflate data, then the CRC and ISIZE.
    let end = at + bsize + 1;
    if end > size || end < at + 26 {
        return Ok(false);
    }
    let mut block = vec![0u8; (bsize + 1) as usize];
    if !read(&mut block, at)? {
        return Ok(false);
    }
    // The decoder checks the CRC and ISIZE; the member must use the whole
    // block, no more and no less, and decode to at most 64 KiB.
    let mut dec = GzDecoder::new(&block[..]);
    let mut plain = Vec::new();
    if (&mut dec)
        .take(BGZF_MAX + 1)
        .read_to_end(&mut plain)
        .is_err()
        || plain.len() as u64 > BGZF_MAX
        || !dec.into_inner().is_empty()
    {
        return Ok(false);
    }
    Ok(within <= plain.len() as u64)
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

/// Whether `v` is a value the manifest writes for the digest kind `k`: a
/// sha256, with its container for `gzip` and `vcf`, or the kind's failure
/// mark.
fn digest_ok(k: &str, v: &str) -> bool {
    let sha = |s: &str| s.len() == 64 && s.bytes().all(|b| matches!(b, b'0'..=b'9' | b'a'..=b'f'));
    let with = |kinds: &[&str]| {
        v.split_once(':')
            .is_some_and(|(c, s)| kinds.contains(&c) && sha(s))
    };
    match k {
        "gzip" => v == "bad" || with(&["bgzf", "gzip"]),
        "vcf" => with(&["bgzf", "gzip", "plain"]),
        "index" => v == "invalid" || sha(v),
        "jellyfish" | "numeric" => sha(v),
        _ => false,
    }
}

fn read_manifest(path: &str) -> Result<BTreeMap<String, Entry>, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let mut lines = text.lines();
    if lines.next() != Some(MAGIC) {
        return Err(format!(
            "{path} is not an install manifest (no {MAGIC} line)"
        ));
    }
    // No comment lines: a path may start with `#`.
    let mut rows = lines;
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
                if !digest_ok(k, v) {
                    return Err(format!("{path}: {p}: {k} digest {v:?}"));
                }
                if d.insert(k.to_owned(), v.to_owned()).is_some() {
                    return Err(format!("{path}: {p} has two {k} digests"));
                }
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
    STAMP_NAMES.contains(&name) || STAMP_ENDINGS.iter().any(|s| name.ends_with(s))
}

/// The stamp name or ending the pattern of one path segment can match, if
/// any: `*` matches any bytes, `?` one byte.
fn stamp_match(seg: &str) -> Option<&'static str> {
    let seg = seg.as_bytes();
    if let Some(n) = STAMP_NAMES.iter().find(|n| glob(seg, n.as_bytes())) {
        return Some(n);
    }
    let fits = |p: &[u8], s: &[u8]| p.iter().zip(s).all(|(&a, &b)| a == b'?' || a == b);
    STAMP_ENDINGS.iter().copied().find(|end| {
        let end = end.as_bytes();
        match seg.iter().rposition(|&b| b == b'*') {
            // No `*`: the pattern's last bytes must fit the ending.
            None => seg.len() >= end.len() && fits(&seg[seg.len() - end.len()..], end),
            // A `*`, then a tail: the `*` takes any bytes the tail leaves.
            Some(i) => {
                let tail = &seg[i + 1..];
                if tail.len() >= end.len() {
                    fits(&tail[tail.len() - end.len()..], end)
                } else {
                    fits(tail, &end[end.len() - tail.len()..])
                }
            }
        }
    })
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
        if let Some(s) = stamp_match(last) {
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
        let cannot = cannot.trim();
        let fits = if skip {
            !cannot.is_empty() && cannot != "-"
        } else {
            cannot == "-"
        };
        if !fits {
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
    let Some((m, mut rest)) = args.split_first() else {
        return Err(String::new());
    };
    let (mut allow, mut seed) = (None, None);
    while let [flag, value, tail @ ..] = rest {
        match flag.as_str() {
            "--allow" if allow.is_none() => allow = Some(value),
            "--not-in" if seed.is_none() => seed = Some(value),
            _ => return Err(String::new()),
        }
        rest = tail;
    }
    if !rest.is_empty() {
        return Err(String::new());
    }
    let m = read_manifest(m)?;
    let allow = match allow {
        Some(p) => read_allow(p)?,
        None => Vec::new(),
    };
    let seed = match seed {
        Some(p) => read_manifest(p)?,
        None => BTreeMap::new(),
    };
    let pick = m
        .iter()
        .filter(|(p, e)| {
            e.kind == "f"
                && e.size != "0"
                && e.size != "-"
                // A file the seed put there as it is: not the component's.
                && seed.get(*p) != Some(*e)
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
    let rel = Path::new(std::ffi::OsStr::from_bytes(&rel));
    // Only a regular file inside DIR, reached through folders: never a
    // path that climbs out or starts at /, and never through a symlink.
    let mut full = Path::new(dir).to_path_buf();
    let mut parts = rel.components().peekable();
    if parts.peek().is_none() {
        return Err(format!("{path} names no file"));
    }
    while let Some(c) = parts.next() {
        let std::path::Component::Normal(c) = c else {
            return Err(format!("{path} is not a path inside {dir}"));
        };
        full.push(c);
        let t = fs::symlink_metadata(&full)
            .map_err(io_err(&full))?
            .file_type();
        let (ok, want) = match parts.peek() {
            Some(_) => (t.is_dir(), "folder"),
            None => (t.is_file(), "regular file"),
        };
        if !ok {
            return Err(format!("{path}: {} is not a {want}", shown(&full)));
        }
    }
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
static RETRY: LazyLock<Regex> = LazyLock::new(|| {
    // A download that failed and is tried again: fetch's notice for a file
    // or a checksum list, and curl's own message. A rerun meets another
    // server state, so these lines come and go between two good runs. A
    // check failure (`failed its check`), the last `ERROR: could not
    // download` and a pull retry stay compared.
    Regex::new(r"^(  (Checksum d|D)ownload attempt \d+/\d+ failed: \S|curl: \(\d+\) )")
        .expect("a valid pattern")
});

static PULL: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(
        r"^(Unable to find image '[^']+' locally$|[0-9a-f]{12}: (Pulling fs layer|Waiting|Downloading|Verifying Checksum|Download complete|Extracting|Pull complete|Already exists)|Digest: sha256:|Status: (Downloaded newer image|Image is up to date) for |\S+: Pulling from )",
    )
    .expect("a valid pattern")
});
static DATE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"\d{4}-\d{2}-\d{2}([T ]\d{2}:\d{2}(:\d{2}(\.\d+)?)?(Z|[+-]\d{2}:?\d{2})?)?")
        .expect("a valid pattern")
});

/// A line with each date or timestamp masked, but not a date inside a
/// name or a path, such as `release-2026-11-26.tar.gz` or
/// `releases/2026-11-26/x`: a letter, digit, `_`, `.`, `-` or `/` right
/// before it, or a letter, digit, `_` or `/` right after it, or a `.` or `-`
/// then a letter or digit, makes it part of the name. A `.` that ends a
/// sentence does not.
fn mask_dates(line: &str) -> String {
    let word = |c: char| c.is_alphanumeric() || c == '_';
    let mut out = String::with_capacity(line.len());
    let mut last = 0;
    for m in DATE.find_iter(line) {
        let before = line[..m.start()]
            .chars()
            .next_back()
            .is_some_and(|c| word(c) || matches!(c, '.' | '-' | '/'));
        let mut rest = line[m.end()..].chars();
        let after = match rest.next() {
            Some('.' | '-') => rest.next().is_some_and(word),
            c => c.is_some_and(|c| word(c) || c == '/'),
        };
        if before || after {
            continue;
        }
        out.push_str(&line[last..m.start()]);
        out.push_str("<time>");
        last = m.end();
    }
    out.push_str(&line[last..]);
    out
}

static MASKS: LazyLock<Vec<(Regex, &'static str)>> = LazyLock::new(|| {
    [
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
                && !RETRY.is_match(l)
        })
        .map(|l| {
            MASKS.iter().fold(mask_dates(l), |s, (re, to)| {
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
    let (ta, tb) = (read(a)?, read(b)?);
    let (la, lb) = (masked_log(&ta), masked_log(&tb));
    let mut text = String::new();
    let mut bad = 0usize;
    // The missing list as it is: a mask would hide a date in a retry path.
    let lines = |t: &str| -> Vec<String> {
        t.split(['\n', '\r'])
            .map(str::trim_end)
            .filter(|l| !l.is_empty())
            .map(str::to_owned)
            .collect()
    };
    let (ua, ub) = (lines(&ta), lines(&tb));
    let (ma, mb) = (missing_list(&ua), missing_list(&ub));
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
    fn a_small_file_read_stops_one_byte_past_its_size() {
        // A writer that keeps appending while the file is read: the read
        // stops one byte past the size the walk saw.
        let got = read_capped(std::io::repeat(7).take(1 << 20), 10).unwrap();
        assert_eq!(got.len(), 11);
    }

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
        for s in STAMP_NAMES {
            assert!(is_stamp(&format!("a/{s}")), "{s}");
        }
        for s in STAMP_ENDINGS {
            assert!(is_stamp(&format!("a/chr1{s}")), "{s}");
        }
        assert!(!is_stamp("a/catalog.json"));
        for p in [
            "DONE",
            "D?NE",
            "*",
            "**",
            "x.stamp",
            "*.stamp",
            "*mp",
            "*.st?mp",
            "foo*",
            "foo.*",
            "?.release",
            "*.inst*",
            "version.*",
            "*e",
            "pins.tsv",
            "*.tsv",
            "p?ns.*",
            "region.sha256",
            "*.sha2?6",
            "*56",
        ] {
            assert!(stamp_match(p).is_some(), "{p} can match a stamp");
        }
        for p in [
            "*.stamps",
            "DONE.txt",
            "x.release.gz",
            "?????",
            "*.vcf.gz",
            "*chr.vcf.gz.tbi",
            "stamp",
            "*.tbi",
            "DONE?",
            "pins.tsv.gz",
            "*.sha256sum",
            "*.csv",
        ] {
            assert_eq!(stamp_match(p), None, "{p} matches no stamp");
        }
        assert!(is_stamp("reference/fine/pins.tsv"));
        assert!(is_stamp("reference/locityper/region.bed.sha256"));
    }

    /// The vcf digest line by line, all in memory: what VcfHasher must give
    /// for any chunking.
    fn vcf_lines(bytes: &[u8]) -> String {
        let mut h = Sha256::new();
        let mut header = true;
        for line in bytes.split_inclusive(|&b| b == b'\n') {
            header = header && line.starts_with(b"##");
            h.update(if header { vcf_mask(line) } else { line.into() });
        }
        hex(h)
    }

    #[test]
    fn the_vcf_digest_is_the_same_in_any_chunks() {
        let cases: [&[u8]; 6] = [
            b"##fileformat=VCFv4.2\n##fileDate=20261126\n##bcftools_viewCommand=view x; Date=Thu Nov 26 2026\n#CHROM\tPOS\n1\t5\n##not header\n",
            b"##fileformat=VCFv4.2\r\n##fileDate=1\r\n#CHROM\r\n",
            b"##fileformat=VCFv4.2\n##fileDate=20261126",
            b"##fileformat=VCFv4.2\n",
            b"",
            b"##fileformat=VCFv4.2\n#CHROM\n1\t5",
        ];
        for bytes in cases {
            for size in [1, 2, 3, 7, 64, bytes.len().max(1)] {
                let mut v = VcfHasher::default();
                for c in bytes.chunks(size) {
                    v.update(c);
                }
                assert_eq!(v.finish(), vcf_lines(bytes), "{size}: {bytes:?}");
            }
            assert_eq!(vcf_digest(bytes).unwrap(), vcf_lines(bytes));
        }
    }

    #[test]
    fn a_header_line_past_the_limit_is_not_held() {
        let mut v = VcfHasher::default();
        let mut whole = Sha256::new();
        let chunk = vec![b'a'; 1 << 20];
        v.update(b"##fileDate=");
        whole.update(b"##fileDate=");
        for _ in 0..((2 * HEADER_MAX) >> 20) {
            v.update(&chunk);
            whole.update(&chunk);
            assert!(v.line.len() <= HEADER_MAX, "{} bytes held", v.line.len());
        }
        assert!(v.long);
        let rest = b"\n##fileDate=1\n#CHROM\n1\t5\n";
        v.update(rest);
        // The long line is hashed as it is; the lines after it are masked.
        whole.update(b"\n##fileDate=<masked>\n#CHROM\n1\t5\n");
        assert_eq!(v.finish(), hex(whole));
    }

    /// One header line `##fileDate=` padded to `len` bytes with its
    /// newline, then a record line.
    fn date_line(len: usize) -> Vec<u8> {
        let mut v = b"##fileDate=".to_vec();
        v.resize(len - 1, b'7');
        v.extend_from_slice(b"\n#CHROM\n1\t5\n");
        v
    }

    #[test]
    fn the_header_cap_holds_wherever_the_newline_falls() {
        for (len, masked) in [(HEADER_MAX, true), (HEADER_MAX + 1, false)] {
            let bytes = date_line(len);
            let mut want = Sha256::new();
            if masked {
                want.update(b"##fileDate=<masked>\n");
            } else {
                want.update(&bytes[..len]);
            }
            want.update(b"#CHROM\n1\t5\n");
            let want = hex(want);
            // The line whole in one chunk, then split just before its
            // newline, so the cap is crossed in the chunk with the newline
            // or the one before it.
            for cut in [bytes.len(), len - 1, len - 2] {
                let mut v = VcfHasher::default();
                v.update(&bytes[..cut]);
                v.update(&bytes[cut..]);
                assert_eq!(v.finish(), want, "line of {len} bytes, cut at {cut}");
            }
        }
    }

    #[test]
    fn bgzf_is_told_from_plain_gzip() {
        let body = b"a\tb\n".repeat(40_000);
        let mut w =
            og_io::bgzf::Writer::new(Vec::new(), std::num::NonZero::new(1).unwrap()).unwrap();
        w.write_all(&body).unwrap();
        let bgzf = w.finish().unwrap();
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gz.write_all(&body).unwrap();
        let gz = gz.finish().unwrap();
        /// Hands its bytes on `step` at a time.
        struct Steps<'a>(&'a [u8], usize);
        impl Read for Steps<'_> {
            fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
                let k = self.0.len().min(self.1).min(buf.len());
                buf[..k].copy_from_slice(&self.0[..k]);
                self.0 = &self.0[k..];
                Ok(k)
            }
        }
        let kind = |bytes: &[u8], step: usize| {
            let mut hr = Hashing {
                inner: Steps(bytes, step),
                h: Sha256::new(),
                n: 0,
                walk: BgzfWalk::default(),
            };
            let _ = gzip_stream_digests(Members::new(&mut hr));
            std::io::copy(&mut hr, &mut std::io::sink()).unwrap();
            hr.walk.kind()
        };
        for step in [1, 7, 4096, usize::MAX] {
            assert_eq!(kind(&bgzf, step), "bgzf", "{step}");
            assert_eq!(kind(&gz, step), "gzip", "{step}");
        }
        // BGZF with gzip after it, or cut short, is not BGZF.
        assert_eq!(kind(&[&bgzf[..], &gz[..]].concat(), 4096), "gzip");
        assert_eq!(kind(&bgzf[..bgzf.len() - 1], 4096), "gzip");
        assert_eq!(kind(b"", 1), "gzip");
        // The first block's BSIZE spans the second block too: the headers
        // still chain to the end, but the first member ends before.
        let first = bgzf_bsize(&bgzf).unwrap() as usize + 1;
        let second = bgzf_bsize(&bgzf[first..]).unwrap() as usize + 1;
        let mut forged = bgzf.clone();
        forged[16..18].copy_from_slice(&((first + second - 1) as u16).to_le_bytes());
        for step in [1, 4096, usize::MAX] {
            assert_eq!(kind(&forged, step), "gzip", "{step}");
        }
    }

    #[test]
    fn a_jellyfish_header_past_the_limit_is_not_read() {
        let input = b"999999999{\"alignment\":8}";
        let mut r = &input[..];
        assert_eq!(jellyfish_digest(&mut r).unwrap(), None);
        assert_eq!(
            r.len(),
            input.len() - 10,
            "only the length and the brace read"
        );
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

    /// The digests a small file of plain gzip, or a plain VCF, gets in
    /// memory.
    fn small(bytes: &[u8]) -> (String, Vec<(&'static str, String)>) {
        let mut d = Vec::new();
        if bytes.starts_with(&[0x1f, 0x8b]) {
            let mut plain = Vec::new();
            MultiGzDecoder::new(bytes).read_to_end(&mut plain).unwrap();
            d.push(("gzip", format!("gzip:{}", sha_of(&plain))));
            d.push(("vcf", format!("gzip:{}", vcf_digest(&plain[..]).unwrap())));
        } else {
            d.push(("vcf", format!("plain:{}", vcf_digest(bytes).unwrap())));
        }
        (sha_of(bytes), d)
    }

    #[test]
    fn a_large_file_is_read_once_for_all_its_digests() {
        // A reader that hands out each byte once and cannot go back.
        struct Once<'a>(&'a [u8], usize);
        impl Read for Once<'_> {
            fn read(&mut self, b: &mut [u8]) -> std::io::Result<usize> {
                let k = b.len().min(self.0.len() - self.1).min(7);
                b[..k].copy_from_slice(&self.0[self.1..self.1 + k]);
                self.1 += k;
                Ok(k)
            }
        }
        let vcf =
            b"##fileformat=VCFv4.2\n##fileDate=20261126\n#CHROM\tPOS\n1\t5\n1\t9\n".repeat(50);
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        std::io::Write::write_all(&mut gz, &vcf).unwrap();
        let gz = gz.finish().unwrap();
        for bytes in [&vcf[..], &gz[..]] {
            let mut r = Once(bytes, 0);
            let (sha, digests, n) = stream_digests(&mut r, &bytes[..64]).unwrap();
            assert_eq!((sha, digests), small(bytes));
            assert_eq!(n, bytes.len() as u64);
            assert_eq!(r.1, bytes.len());
        }
        // The same content as BGZF: the same digests, with its container.
        let mut w =
            og_io::bgzf::Writer::new(Vec::new(), std::num::NonZero::new(1).unwrap()).unwrap();
        w.write_all(&vcf).unwrap();
        let bg = w.finish().unwrap();
        let (_, got, _) = stream_digests(&mut Once(&bg, 0), &bg[..64]).unwrap();
        let (_, plain) = small(&gz);
        let want: Vec<_> = plain
            .into_iter()
            .map(|(k, v)| (k, v.replacen("gzip:", "bgzf:", 1)))
            .collect();
        assert_eq!(got, want);
        // Bytes after the gzip stream end the decoding; they still count in the sha.
        let mut bad = gz.clone();
        bad.extend_from_slice(b"not gzip at all");
        let mut r = Once(&bad, 0);
        let (sha, _, n) = stream_digests(&mut r, &bad[..64]).unwrap();
        assert_eq!(sha, sha_of(&bad));
        assert_eq!(n, bad.len() as u64);
        assert_eq!(r.1, bad.len());
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
