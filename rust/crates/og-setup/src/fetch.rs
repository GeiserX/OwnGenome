//! `og setup fetch`: the one way setup downloads a file, and the table of
//! every file it downloads (`assets/data-sources.tsv`).
//!
//! A download goes to `DEST.part` (an interrupted one is resumed), is
//! checked, and only then renamed to DEST, so a DEST that exists is whole.
//! The check is the checksum when one is given (a value, or a checksum file
//! at a URL: its only line, or the line naming the file), else a gzip check
//! for a `.gz` DEST, else non-empty. Each try downloads, then checks: a
//! broken transfer keeps the `.part` for the next try, a download that fails
//! its check is deleted and fetched anew. `FETCH_TRIES` tries (default 3),
//! `FETCH_WAIT` seconds apart (default 5, then 10, ...).
//!
//! Reuse is by checksum: a DEST already there that matches its checksum is
//! kept and no request is made for it (a checksum read from a URL costs
//! that one request). A DEST that does not match is never replaced
//! silently: a file in the downloads folder (a source of the table, or any
//! `URL DEST` under `OG_DOWNLOADS`; it may be one a user put there) is
//! refused and left as it is; a `URL DEST` call elsewhere, which setup makes
//! for a file it installs itself, says it downloads the file again. With no
//! checksum, a file in the downloads folder that passes the check a new
//! download would get (a whole gzip, or not empty) is used; elsewhere
//! nothing proves an old file is the right one, so it is downloaded again.
//!
//! `curl` does the transfer, and `wget` when there is no curl, with the
//! same options `scripts/lib/common.sh` gave them: a connection over 30 s
//! or a transfer slower than 1 KB/s for `FETCH_STALL` seconds (default 120)
//! fails the try instead of hanging. The URL is curl's last argument, and
//! the waits run the `sleep` command, so a test harness's own curl, wget
//! and sleep stand in for them.
//!
//! The table: one row per file (or per family of files, `{part}`), the
//! installer and steps that need it, where it is stored, its URL, its
//! checksum and its size. `${NAME}` in a URL or a value is the pin NAME of
//! `versions.env`; `${NAME:-default}` is the environment variable NAME,
//! or the default. A value `@ID` is the checksum list that row ID stores.
//! A path under `downloads/` is in the downloads folder: `OG_DOWNLOADS`,
//! or `$GENOME_DIR/downloads`; any other path is under `GENOME_DIR`.

use flate2::read::MultiGzDecoder;
use indexmap::IndexMap;
use md5::Md5;
use sha2::{Digest, Sha256};
use std::fmt;
use std::fs::{self, File};
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// The table, relative to the checkout.
pub const TABLE: &str = "assets/data-sources.tsv";
/// The table's columns, in order.
pub const COLUMNS: [&str; 10] = [
    "id",
    "installers",
    "steps",
    "path",
    "parts",
    "url",
    "check",
    "value",
    "bytes",
    "licence",
];

// --- checks ---------------------------------------------------------------

/// How a download is checked.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// No checksum: a whole gzip file for a `.gz` DEST, else non-empty.
    None,
    Md5,
    Sha256,
    /// The BSD checksum `sum` prints (Ensembl's CHECKSUMS files).
    Sum,
    /// `BYTES:SHA256-OF-FIRST-MIB:SHA256-OF-LAST-MIB`, and a whole gzip file.
    GzipEnds,
    /// Google Cloud Storage's CRC-32C, base64 as its `x-goog-hash` header.
    Crc32c,
}

impl Kind {
    pub fn parse(s: &str) -> Result<Kind, String> {
        Ok(match s {
            "" => Kind::None,
            "md5" => Kind::Md5,
            "sha256" => Kind::Sha256,
            "sum" => Kind::Sum,
            "gzip-ends" => Kind::GzipEnds,
            "crc32c" => Kind::Crc32c,
            other => return Err(format!("ERROR: unknown checksum type '{other}'")),
        })
    }

    pub fn name(self) -> &'static str {
        match self {
            Kind::None => "",
            Kind::Md5 => "md5",
            Kind::Sha256 => "sha256",
            Kind::Sum => "sum",
            Kind::GzipEnds => "gzip-ends",
            Kind::Crc32c => "crc32c",
        }
    }
}

/// Why a file failed its check.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Bad {
    /// The file is not the one: a new download may fix it.
    Wrong(String),
    /// An error no new download fixes (an unreadable file, a malformed pin).
    Fatal(String),
}

fn read_err(path: &Path, e: &io::Error) -> Bad {
    Bad::Fatal(format!("ERROR: cannot read {}: {e}", path.display()))
}

/// The digest of a whole file, as the shell tools print it.
pub fn digest(kind: Kind, path: &Path) -> io::Result<String> {
    let mut f = File::open(path)?;
    let mut buf = vec![0u8; 1 << 20];
    match kind {
        Kind::Md5 => {
            let mut h = Md5::new();
            loop_read(&mut f, &mut buf, |b| h.update(b))?;
            Ok(hex(&h.finalize()))
        }
        Kind::Sha256 => {
            let mut h = Sha256::new();
            loop_read(&mut f, &mut buf, |b| h.update(b))?;
            Ok(hex(&h.finalize()))
        }
        Kind::Sum => {
            let mut s: u32 = 0;
            loop_read(&mut f, &mut buf, |b| {
                for &x in b {
                    s = ((s >> 1) | ((s & 1) << 15)).wrapping_add(u32::from(x)) & 0xffff;
                }
            })?;
            Ok(s.to_string())
        }
        Kind::Crc32c => {
            let mut c = Crc32c::default();
            loop_read(&mut f, &mut buf, |b| c.update(b))?;
            Ok(c.base64())
        }
        Kind::None | Kind::GzipEnds => Err(io::Error::other(format!(
            "no digest for check '{}'",
            kind.name()
        ))),
    }
}

fn loop_read(f: &mut File, buf: &mut [u8], mut each: impl FnMut(&[u8])) -> io::Result<()> {
    loop {
        let n = f.read(buf)?;
        if n == 0 {
            return Ok(());
        }
        each(&buf[..n]);
    }
}

fn hex(b: &[u8]) -> String {
    b.iter().map(|x| format!("{x:02x}")).collect()
}

/// CRC-32C (Castagnoli), bit by bit over a 256-entry table.
#[derive(Debug)]
pub struct Crc32c {
    table: [u32; 256],
    crc: u32,
}

impl Default for Crc32c {
    fn default() -> Self {
        let mut table = [0u32; 256];
        for (i, t) in table.iter_mut().enumerate() {
            let mut c = i as u32;
            for _ in 0..8 {
                c = if c & 1 == 1 {
                    (c >> 1) ^ 0x82f6_3b78
                } else {
                    c >> 1
                };
            }
            *t = c;
        }
        Crc32c { table, crc: !0 }
    }
}

impl Crc32c {
    pub fn update(&mut self, b: &[u8]) {
        for &x in b {
            self.crc = self.table[((self.crc ^ u32::from(x)) & 0xff) as usize] ^ (self.crc >> 8);
        }
    }

    pub fn value(&self) -> u32 {
        !self.crc
    }

    /// Big-endian bytes in base64, as Google writes `crc32c=...`.
    pub fn base64(&self) -> String {
        base64(&self.value().to_be_bytes())
    }
}

fn base64(b: &[u8]) -> String {
    const A: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::new();
    for ch in b.chunks(3) {
        let n = ch.len();
        let v = (u32::from(ch[0]) << 16)
            | (u32::from(*ch.get(1).unwrap_or(&0)) << 8)
            | u32::from(*ch.get(2).unwrap_or(&0));
        for i in 0..4 {
            if i <= n {
                out.push(A[((v >> (18 - 6 * i)) & 63) as usize] as char);
            } else {
                out.push('=');
            }
        }
    }
    out
}

/// True when `path` is one or more whole gzip members and nothing else
/// (what `gzip -t` accepts).
pub fn gzip_ok(path: &Path) -> bool {
    let Ok(f) = File::open(path) else {
        return false;
    };
    let mut d = MultiGzDecoder::new(io::BufReader::with_capacity(1 << 20, f));
    let mut sink = vec![0u8; 1 << 20];
    loop {
        match d.read(&mut sink) {
            Ok(0) => break,
            Ok(_) => {}
            Err(_) => return false,
        }
    }
    // An empty file has no member.
    fs::metadata(path).is_ok_and(|m| m.len() > 0)
}

fn sha256_range(f: &mut File, from: SeekFrom, len: u64) -> io::Result<String> {
    f.seek(from)?;
    let mut h = Sha256::new();
    let mut buf = vec![0u8; 1 << 16];
    let mut left = len;
    while left > 0 {
        let want = buf.len().min(usize::try_from(left).unwrap_or(usize::MAX));
        let n = f.read(&mut buf[..want])?;
        if n == 0 {
            break;
        }
        h.update(&buf[..n]);
        left -= n as u64;
    }
    Ok(hex(&h.finalize()))
}

/// The check of one file: `Ok` when FILE is good. NAME is the download's
/// name in messages; DEST decides the gzip check when there is no checksum.
pub fn check(kind: Kind, want: &str, file: &Path, dest: &Path, name: &str) -> Result<(), Bad> {
    match kind {
        Kind::GzipEnds => {
            const MIB: u64 = 1 << 20;
            let mut it = want.splitn(3, ':');
            let (bytes, first, last) = (
                it.next().unwrap_or(""),
                it.next().unwrap_or(""),
                it.next().unwrap_or(""),
            );
            let size = fs::metadata(file).map_err(|e| read_err(file, &e))?.len();
            let mut f = File::open(file).map_err(|e| read_err(file, &e))?;
            let head = sha256_range(&mut f, SeekFrom::Start(0), MIB.min(size))
                .map_err(|e| read_err(file, &e))?;
            let tail = sha256_range(&mut f, SeekFrom::Start(size.saturating_sub(MIB)), MIB)
                .map_err(|e| read_err(file, &e))?;
            if size.to_string() != bytes || head != first || tail != last || !gzip_ok(file) {
                return Err(Bad::Wrong(format!(
                    "{name} failed its pinned size, SHA-256 range or gzip check"
                )));
            }
            Ok(())
        }
        Kind::None => {
            if dest.to_string_lossy().ends_with(".gz") {
                if !gzip_ok(file) {
                    return Err(Bad::Wrong(format!("{name} is not a complete gzip file")));
                }
            } else if fs::metadata(file).map_or(true, |m| m.len() == 0) {
                return Err(Bad::Wrong(format!("{name} downloaded empty")));
            }
            Ok(())
        }
        _ => {
            let mut want = want.to_owned();
            if kind == Kind::Sum {
                if want.is_empty() || !want.bytes().all(|b| b.is_ascii_digit()) {
                    return Err(Bad::Fatal(format!(
                        "ERROR: '{want}' is not a sum checksum for {name}"
                    )));
                }
                want = want.trim_start_matches('0').to_owned();
                if want.is_empty() {
                    want.push('0');
                }
            }
            let got = digest(kind, file).map_err(|e| read_err(file, &e))?;
            if got != want {
                return Err(Bad::Wrong(format!(
                    "{} checksum of {name} is {got}, expected {want}",
                    kind.name()
                )));
            }
            Ok(())
        }
    }
}

// --- transfers ------------------------------------------------------------

/// The environment knobs of a fetch, read once.
#[derive(Debug, Clone)]
pub struct Knobs {
    pub tries: u32,
    /// `None`: 5 s times the try's number.
    pub wait: Option<String>,
    pub stall: u64,
}

impl Knobs {
    /// From `FETCH_TRIES`, `FETCH_WAIT` and `FETCH_STALL`, with the shell's
    /// messages for a bad value.
    pub fn from_env() -> Result<Knobs, String> {
        Self::from(&|k| std::env::var(k).ok())
    }

    pub fn from(env: &dyn Fn(&str) -> Option<String>) -> Result<Knobs, String> {
        let raw = env("FETCH_TRIES").unwrap_or_default();
        let t = if raw.is_empty() {
            "3".to_owned()
        } else {
            raw.clone()
        };
        let tries = if t.bytes().all(|b| b.is_ascii_digit()) {
            t.trim_start_matches('0').parse::<u32>().unwrap_or(0)
        } else {
            0
        };
        if tries == 0 {
            return Err(format!(
                "ERROR: FETCH_TRIES must be a whole number above 0, got '{raw}'"
            ));
        }
        let wait = env("FETCH_WAIT").filter(|w| !w.is_empty());
        if let Some(w) = &wait
            && !w.bytes().all(|b| b.is_ascii_digit())
        {
            return Err(format!(
                "ERROR: FETCH_WAIT must be a whole number of seconds, got '{w}'"
            ));
        }
        let raw = env("FETCH_STALL").unwrap_or_default();
        let s = if raw.is_empty() {
            "120".to_owned()
        } else {
            raw.clone()
        };
        let stall = if s.bytes().all(|b| b.is_ascii_digit()) {
            s.trim_start_matches('0').parse::<u64>().unwrap_or(0)
        } else {
            0
        };
        if stall == 0 {
            return Err(format!(
                "ERROR: FETCH_STALL must be a whole number of seconds above 0, got '{raw}'"
            ));
        }
        Ok(Knobs { tries, wait, stall })
    }

    fn pause(&self, try_no: u32) {
        let secs = self
            .wait
            .clone()
            .unwrap_or_else(|| (u64::from(try_no) * 5).to_string());
        // The sleep command, as the shell ran it.
        let ran = Command::new("sleep")
            .arg(&secs)
            .stdin(Stdio::null())
            .status()
            .is_ok_and(|s| s.success());
        if !ran && let Ok(n) = secs.parse::<u64>() {
            std::thread::sleep(std::time::Duration::from_secs(n));
        }
    }
}

/// The first executable NAME on PATH (`command -v`).
pub fn on_path(name: &str) -> Option<PathBuf> {
    let path = std::env::var_os("PATH")?;
    std::env::split_paths(&path)
        .map(|d| d.join(name))
        .find(|p| is_executable(p))
}

fn is_executable(p: &Path) -> bool {
    use std::os::unix::fs::PermissionsExt;
    fs::metadata(p).is_ok_and(|m| m.is_file() && m.permissions().mode() & 0o111 != 0)
}

/// Where one transfer writes.
pub enum Out<'a> {
    File(&'a Path),
    Stdout,
}

/// One transfer of URL, resuming a partial OUT; `Err` with a message when
/// no transfer could start, `Ok(false)` when it failed. To stdout, the body
/// is returned.
pub fn get(url: &str, out: Out, knobs: &Knobs) -> Result<(bool, Vec<u8>), String> {
    let stall = knobs.stall.to_string();
    if on_path("curl").is_some() {
        let bound = [
            "--connect-timeout",
            "30",
            "--speed-limit",
            "1024",
            "--speed-time",
            stall.as_str(),
        ];
        match out {
            Out::Stdout => {
                let o = Command::new("curl")
                    .arg("-fsSL")
                    .args(bound)
                    .arg(url)
                    .stdin(Stdio::null())
                    .stderr(Stdio::inherit())
                    .output()
                    .map_err(|e| format!("ERROR: cannot run curl: {e}"))?;
                Ok((o.status.success(), o.stdout))
            }
            Out::File(p) => {
                let run = |resume: bool| -> Result<Option<i32>, String> {
                    let mut c = Command::new("curl");
                    c.arg("-fL").args(bound);
                    if resume {
                        c.args(["-C", "-"]);
                    }
                    c.arg("-o").arg(p).arg(url);
                    c.stdin(Stdio::null()).stdout(Stdio::null());
                    c.status()
                        .map(|s| s.code())
                        .map_err(|e| format!("ERROR: cannot run curl: {e}"))
                };
                let code = run(true)?;
                if code == Some(33) {
                    // The server cannot resume: one fresh download.
                    let _ = fs::remove_file(p);
                    return Ok((run(false)? == Some(0), vec![]));
                }
                Ok((code == Some(0), vec![]))
            }
        }
    } else if on_path("wget").is_some() {
        let bound = [
            "--connect-timeout=30".to_owned(),
            format!("--read-timeout={stall}"),
            "--tries=1".to_owned(),
        ];
        match out {
            Out::Stdout => {
                let o = Command::new("wget")
                    .arg("-q")
                    .args(&bound)
                    .args(["-O", "-"])
                    .arg(url)
                    .stdin(Stdio::null())
                    .stderr(Stdio::inherit())
                    .output()
                    .map_err(|e| format!("ERROR: cannot run wget: {e}"))?;
                Ok((o.status.success(), o.stdout))
            }
            Out::File(p) => {
                let s = Command::new("wget")
                    .arg("-c")
                    .args(&bound)
                    .arg("-O")
                    .arg(p)
                    .arg(url)
                    .stdin(Stdio::null())
                    .stdout(Stdio::null())
                    .status()
                    .map_err(|e| format!("ERROR: cannot run wget: {e}"))?;
                Ok((s.success(), vec![]))
            }
        }
    } else {
        Err("ERROR: neither curl nor wget is installed.".to_owned())
    }
}

/// The checksum a checksum file names for NAME: its only non-empty line's
/// first field, else the first field of the last line whose last field's
/// base name is NAME.
pub fn checksum_in(text: &str, name: &str) -> Option<String> {
    let (mut lines, mut first, mut hit) = (0, String::new(), None);
    for l in text.lines() {
        let f: Vec<&str> = l.split_whitespace().collect();
        let Some(last) = f.last() else { continue };
        lines += 1;
        first = f[0].to_owned();
        if last.rsplit('/').next() == Some(name) {
            hit = Some(f[0].to_owned());
        }
    }
    if lines == 1 { Some(first) } else { hit }
}

/// The file name a URL ends in, without its query or fragment.
fn url_name(url: &str) -> &str {
    let path = url.split(['?', '#']).next().unwrap_or(url);
    path.rsplit('/').next().unwrap_or(path)
}

/// The bytes of a `file://` URL: its path as written, else percent-decoded.
fn read_local(url: &str) -> io::Result<Vec<u8>> {
    let raw = url.strip_prefix("file://").unwrap_or(url);
    match fs::read(raw) {
        Err(e) if e.kind() == io::ErrorKind::NotFound => {
            let b = raw.as_bytes();
            let (mut out, mut i) = (Vec::with_capacity(b.len()), 0);
            while i < b.len() {
                let hex = |c: u8| (c as char).to_digit(16);
                match (
                    b[i],
                    b.get(i + 1).copied().and_then(hex),
                    b.get(i + 2).copied().and_then(hex),
                ) {
                    (b'%', Some(h), Some(l)) => {
                        out.push((h * 16 + l) as u8);
                        i += 3;
                    }
                    (c, _, _) => {
                        out.push(c);
                        i += 1;
                    }
                }
            }
            use std::os::unix::ffi::OsStrExt;
            fs::read(Path::new(std::ffi::OsStr::from_bytes(&out)))
        }
        x => x,
    }
}

fn is_url(s: &str) -> bool {
    ["http://", "https://", "ftp://", "file://"]
        .iter()
        .any(|p| s.starts_with(p))
}

/// What to do with a DEST that is there but fails its check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OnWrong {
    /// Download it again, and say so (a file setup installs itself).
    Replace,
    /// Leave it and fail (a file in the downloads folder).
    Refuse,
}

/// One fetch request.
#[derive(Debug, Clone)]
pub struct Fetch<'a> {
    pub url: &'a str,
    pub dest: &'a Path,
    pub kind: Kind,
    /// A checksum, or the URL of a checksum file.
    pub want: &'a str,
    pub on_wrong: OnWrong,
    /// The table row, named in messages.
    pub source: Option<&'a str>,
}

/// How a fetch ended well.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Done {
    Downloaded,
    Reused,
}

/// Fetches; progress and errors go to `err`.
/// The error carries nothing: every failure is written to .err. first.
#[allow(clippy::result_unit_err)]
pub fn fetch(r: &Fetch, knobs: &Knobs, err: &mut dyn Write) -> Result<Done, ()> {
    let name = url_name(r.url);
    let mut want = r.want.to_owned();
    if r.kind != Kind::None && want.trim().is_empty() {
        let _ = writeln!(err, "ERROR: no {} value given for {name}", r.kind.name());
        return Err(());
    }
    if r.kind != Kind::None && want.starts_with("file://") {
        // A local checksum list is read here: wget has no file:// scheme,
        // and a path with a space is not a valid URL for curl.
        match read_local(&want) {
            Ok(body) => match checksum_in(&String::from_utf8_lossy(&body), name) {
                Some(l) if !l.is_empty() => want = l,
                _ => {
                    let _ = writeln!(err, "ERROR: {want} names no checksum for {name}");
                    return Err(());
                }
            },
            Err(e) => {
                let _ = writeln!(err, "ERROR: cannot read the checksum list {want}: {e}");
                return Err(());
            }
        }
    } else if r.kind != Kind::None && is_url(&want) {
        let mut line = None;
        for i in 1..=knobs.tries {
            match get(&want, Out::Stdout, knobs) {
                Ok((true, body)) => {
                    // A read that worked but names no checksum is not retried.
                    line = Some(checksum_in(&String::from_utf8_lossy(&body), name));
                    break;
                }
                Ok(_) => {}
                Err(m) => {
                    let _ = writeln!(err, "{m}");
                }
            }
            let _ = writeln!(
                err,
                "  Checksum download attempt {i}/{} failed: {want}",
                knobs.tries
            );
            if i < knobs.tries {
                knobs.pause(i);
            }
        }
        match line.flatten() {
            Some(l) if !l.is_empty() => want = l,
            _ => {
                let _ = writeln!(
                    err,
                    "ERROR: could not read the checksum of {name} from {want}"
                );
                return Err(());
            }
        }
    }
    let label = r.source.unwrap_or(name);
    if let Some(d) = r.dest.parent().filter(|d| !d.as_os_str().is_empty())
        && let Err(e) = fs::create_dir_all(d)
    {
        let _ = writeln!(err, "ERROR: cannot create {}: {e}", d.display());
        return Err(());
    }
    // The downloads folder may be shared: one fetch of a file at a time, from
    // the check of what is there to the rename of what was checked.
    let _lock = if r.on_wrong == OnWrong::Refuse {
        match lock(r.dest, label, err) {
            Ok(l) => l,
            Err(()) => return Err(()),
        }
    } else {
        None
    };
    // In the downloads folder a file with no checksum is used when it passes
    // the check a new download would get (whole gzip, or non-empty).
    if r.dest.is_file() && (r.kind != Kind::None || r.on_wrong == OnWrong::Refuse) {
        match check(r.kind, &want, r.dest, r.dest, name) {
            Ok(()) => {
                let what = match r.kind {
                    Kind::None if is_gz(r.dest) => "gzip check",
                    Kind::None => "check (not empty)",
                    k => k.name(),
                };
                let _ = writeln!(
                    err,
                    "  {label}: {} is already here and matches its {what}; not downloaded.",
                    r.dest.display()
                );
                return Ok(Done::Reused);
            }
            Err(Bad::Fatal(m)) => {
                let _ = writeln!(err, "{m}");
                return Err(());
            }
            Err(Bad::Wrong(why)) => match r.on_wrong {
                OnWrong::Replace => {
                    let _ = writeln!(
                        err,
                        "  {} is not the pinned file ({why}): downloading it again.",
                        r.dest.display()
                    );
                }
                OnWrong::Refuse => {
                    let _ = writeln!(
                        err,
                        "ERROR: {} is there but is not the file of {label}: {why}. It was left as it is: \
                         move it away (or run og setup prune {label}) and run setup again.",
                        r.dest.display()
                    );
                    return Err(());
                }
            },
        }
    }
    let part = PathBuf::from(format!("{}.part", r.dest.display()));
    let mut bad = String::new();
    for i in 1..=knobs.tries {
        let ok = match get(r.url, Out::File(&part), knobs) {
            Ok((ok, _)) => ok,
            Err(m) => {
                let _ = writeln!(err, "{m}");
                false
            }
        };
        if ok {
            match check(r.kind, &want, &part, r.dest, name) {
                Ok(()) => {
                    let done = if r.on_wrong == OnWrong::Refuse {
                        publish_new(&part, r.dest)
                    } else {
                        fs::rename(&part, r.dest)
                    };
                    return match done {
                        Ok(()) => Ok(Done::Downloaded),
                        Err(e) if e.kind() == io::ErrorKind::AlreadyExists => {
                            let _ = writeln!(
                                err,
                                "ERROR: {} appeared while {label} was downloading; it was left as it is. \
                                 The checked download is {}: keep one of them and run setup again.",
                                r.dest.display(),
                                part.display()
                            );
                            Err(())
                        }
                        Err(e) => {
                            let _ = writeln!(
                                err,
                                "ERROR: cannot rename {} to {}: {e}",
                                part.display(),
                                r.dest.display()
                            );
                            Err(())
                        }
                    };
                }
                Err(Bad::Fatal(m)) => {
                    let _ = writeln!(err, "{m}");
                    return Err(());
                }
                Err(Bad::Wrong(why)) => {
                    if let Err(e) = fs::remove_file(&part) {
                        let _ = writeln!(
                            err,
                            "ERROR: {why}, and {} cannot be removed: {e}",
                            part.display()
                        );
                        return Err(());
                    }
                    let _ = writeln!(
                        err,
                        "  Download attempt {i}/{} failed its check: {why}. Removed it.",
                        knobs.tries
                    );
                    bad = why;
                }
            }
        } else {
            bad.clear();
            let _ = writeln!(
                err,
                "  Download attempt {i}/{} failed: {}",
                knobs.tries, r.url
            );
        }
        if i < knobs.tries {
            knobs.pause(i);
        }
    }
    if bad.is_empty() {
        let _ = writeln!(
            err,
            "ERROR: could not download {} (the partial file {} is kept and resumed next time)",
            r.url,
            part.display()
        );
    } else {
        let _ = writeln!(
            err,
            "ERROR: {bad}. Removed the download after {} tries; run again to fetch it anew.",
            knobs.tries
        );
    }
    Err(())
}

/// Holds DEST's lock: the hidden file `.NAME.lock` beside it, kept (it is
/// empty), as removing it would let a third fetch lock another file. In a
/// folder this user cannot write there is no lock (and nothing to publish).
fn lock(dest: &Path, label: &str, err: &mut dyn Write) -> Result<Option<File>, ()> {
    let name = dest
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let p = dest.with_file_name(format!(".{name}.lock"));
    let f = match File::options()
        .create(true)
        .truncate(false)
        .write(true)
        .open(&p)
    {
        Ok(f) => f,
        Err(e)
            if matches!(
                e.kind(),
                io::ErrorKind::PermissionDenied | io::ErrorKind::ReadOnlyFilesystem
            ) =>
        {
            return Ok(None);
        }
        Err(e) => {
            let _ = writeln!(err, "ERROR: cannot open the lock {}: {e}", p.display());
            return Err(());
        }
    };
    match f.try_lock() {
        Ok(()) => return Ok(Some(f)),
        Err(fs::TryLockError::WouldBlock) => {
            let _ = writeln!(
                err,
                "  {label}: another fetch of {} is running; waiting for it.",
                dest.display()
            );
        }
        Err(fs::TryLockError::Error(e)) => {
            let _ = writeln!(err, "ERROR: cannot lock {}: {e}", p.display());
            return Err(());
        }
    }
    f.lock().map(|()| Some(f)).map_err(|e| {
        let _ = writeln!(err, "ERROR: cannot lock {}: {e}", p.display());
    })
}

/// Moves PART to DEST unless DEST is there (`AlreadyExists`): a hard link,
/// which never replaces, or a rename where the filesystem has no links.
pub fn publish_new(part: &Path, dest: &Path) -> io::Result<()> {
    match fs::hard_link(part, dest) {
        Ok(()) => fs::remove_file(part),
        Err(e) if e.kind() == io::ErrorKind::AlreadyExists => Err(e),
        Err(_) if fs::symlink_metadata(dest).is_ok() => {
            Err(io::Error::from(io::ErrorKind::AlreadyExists))
        }
        Err(_) => fs::rename(part, dest),
    }
}

fn is_gz(p: &Path) -> bool {
    p.extension().is_some_and(|e| e == "gz")
}

/// True when DEST is inside the downloads folder of this run (OG_DOWNLOADS,
/// or GENOME_DIR/downloads), so a wrong file there is refused, never replaced.
pub fn in_downloads(dest: &Path, env: &dyn Fn(&str) -> Option<String>) -> bool {
    let get = |k: &str| env(k).filter(|v| !v.is_empty());
    let Some(d) = get("OG_DOWNLOADS")
        .map(PathBuf::from)
        .or_else(|| get("GENOME_DIR").map(|g| Path::new(&g).join("downloads")))
    else {
        return false;
    };
    resolved(dest).starts_with(resolved(&d))
}

/// P with its longest existing head resolved (symlinks, `.` and `..`) and
/// the rest, which does not exist yet, joined to it with `..` applied.
pub fn resolved(p: &Path) -> PathBuf {
    use std::path::Component;
    let abs = std::path::absolute(p).unwrap_or_else(|_| p.to_path_buf());
    let comps: Vec<Component> = abs.components().collect();
    for k in (1..=comps.len()).rev() {
        let head: PathBuf = comps[..k].iter().collect();
        if let Ok(mut out) = head.canonicalize() {
            for c in &comps[k..] {
                match c {
                    Component::ParentDir => {
                        out.pop();
                    }
                    Component::Normal(n) => out.push(n),
                    _ => {}
                }
            }
            return out;
        }
    }
    abs
}

// --- the table ------------------------------------------------------------

/// One row of the table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Source {
    pub id: String,
    pub installers: Vec<String>,
    pub steps: String,
    pub path: String,
    pub parts: Vec<String>,
    pub url: String,
    pub check: String,
    pub value: String,
    pub bytes: Option<u64>,
    pub licence: String,
    /// The table's line number, for messages.
    pub line: usize,
}

impl Source {
    /// True when the file lives in the downloads folder.
    pub fn is_download(&self) -> bool {
        self.path.starts_with("downloads/")
    }
}

/// `parts`: words, with `A..B` for the whole numbers A to B.
pub fn expand_parts(s: &str) -> Result<Vec<String>, String> {
    let mut out = vec![];
    for w in s.split_whitespace() {
        if let Some((a, b)) = w.split_once("..") {
            let (a, b): (u32, u32) = (
                a.parse().map_err(|_| format!("bad range '{w}'"))?,
                b.parse().map_err(|_| format!("bad range '{w}'"))?,
            );
            if a > b {
                return Err(format!("bad range '{w}'"));
            }
            out.extend((a..=b).map(|n| n.to_string()));
        } else {
            out.push(w.to_owned());
        }
    }
    Ok(out)
}

/// The table, keyed by id.
#[derive(Debug, Clone, Default)]
pub struct Sources {
    pub rows: IndexMap<String, Source>,
}

impl Sources {
    pub fn parse(text: &str) -> Result<Sources, String> {
        let mut s = Sources::default();
        let mut header = false;
        let mut paths: IndexMap<String, String> = IndexMap::new();
        for (n, l) in text.lines().enumerate() {
            let n = n + 1;
            if l.starts_with('#') || l.trim().is_empty() {
                continue;
            }
            let f: Vec<&str> = l.split('\t').collect();
            if !header {
                if f != COLUMNS {
                    return Err(format!(
                        "{TABLE}:{n}: the header must be {}",
                        COLUMNS.join(" ")
                    ));
                }
                header = true;
                continue;
            }
            if f.len() != COLUMNS.len() {
                return Err(format!(
                    "{TABLE}:{n}: {} fields, not {}",
                    f.len(),
                    COLUMNS.len()
                ));
            }
            let id = f[0].to_owned();
            if id.is_empty()
                || !id
                    .bytes()
                    .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'_')
            {
                return Err(format!("{TABLE}:{n}: id '{id}' is not [a-z0-9_]+"));
            }
            Kind::parse(f[6]).map_err(|e| format!("{TABLE}:{n}: {e}"))?;
            let parts = expand_parts(f[4]).map_err(|e| format!("{TABLE}:{n}: {e}"))?;
            if parts.is_empty() == (f[3].contains("{part}") || f[5].contains("{part}")) {
                return Err(format!(
                    "{TABLE}:{n}: {{part}} in the path or URL goes with a parts list, and only with one"
                ));
            }
            let bytes = if f[8].is_empty() {
                None
            } else {
                Some(
                    f[8].parse()
                        .map_err(|_| format!("{TABLE}:{n}: bytes '{}' is not a number", f[8]))?,
                )
            };
            if f[3].is_empty() || f[3].starts_with('/') || f[3].split('/').any(|p| p == "..") {
                return Err(format!(
                    "{TABLE}:{n}: path '{}' must be relative, inside GENOME_DIR",
                    f[3]
                ));
            }
            // Two families may share a path template with different parts.
            if let Some(other) = paths.insert(format!("{}\t{}", f[3], f[4]), id.clone()) {
                return Err(format!("{TABLE}:{n}: path {} is also {other}'s", f[3]));
            }
            let row = Source {
                id: id.clone(),
                installers: f[1]
                    .split(',')
                    .map(str::to_owned)
                    .filter(|x| !x.is_empty())
                    .collect(),
                steps: f[2].to_owned(),
                path: f[3].to_owned(),
                parts,
                url: f[5].to_owned(),
                check: f[6].to_owned(),
                value: f[7].to_owned(),
                bytes,
                licence: f[9].to_owned(),
                line: n,
            };
            if s.rows.insert(id.clone(), row).is_some() {
                return Err(format!("{TABLE}:{n}: id {id} is in the table twice"));
            }
        }
        if !header {
            return Err(format!("{TABLE}: no header line"));
        }
        for r in s.rows.values() {
            if let Some(o) = r.value.strip_prefix('@')
                && !s.rows.contains_key(o)
            {
                return Err(format!(
                    "{TABLE}:{}: {} names @{o}, which has no row",
                    r.line, r.id
                ));
            }
        }
        Ok(s)
    }

    pub fn read(root: &Path) -> Result<Sources, String> {
        let p = root.join(TABLE);
        let t = fs::read_to_string(&p)
            .map_err(|e| format!("ERROR: cannot read {}: {e}", p.display()))?;
        Sources::parse(&t).map_err(|e| format!("ERROR: {e}"))
    }

    pub fn get(&self, id: &str) -> Result<&Source, String> {
        self.rows
            .get(id)
            .ok_or_else(|| format!("ERROR: {id} is not a source of {TABLE}"))
    }
}

/// Where a run reads pins and settings.
/// How a run reads one environment variable.
pub type Env = Box<dyn Fn(&str) -> Option<String>>;

pub struct Ctx {
    pub root: PathBuf,
    pins: IndexMap<String, String>,
    env: Env,
}

impl fmt::Debug for Ctx {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Ctx")
            .field("root", &self.root)
            .finish_non_exhaustive()
    }
}

/// The `NAME="value"` lines of a `versions.env`.
pub fn versions_pins(text: &str) -> IndexMap<String, String> {
    let mut m = IndexMap::new();
    for l in text.lines() {
        let Some((k, v)) = l.split_once('=') else {
            continue;
        };
        if k.is_empty()
            || !k
                .bytes()
                .all(|b| b.is_ascii_uppercase() || b.is_ascii_digit() || b == b'_')
        {
            continue;
        }
        let v = v.trim_start();
        let v = if let Some(rest) = v.strip_prefix('"') {
            rest.split('"').next().unwrap_or("")
        } else {
            v.split_whitespace().next().unwrap_or("")
        };
        m.insert(k.to_owned(), v.to_owned());
    }
    m
}

impl Ctx {
    pub fn new(root: &Path, env: Env) -> Result<Ctx, String> {
        let p = root.join("versions.env");
        let t = fs::read_to_string(&p)
            .map_err(|e| format!("ERROR: cannot read {}: {e}", p.display()))?;
        Ok(Ctx {
            root: root.to_path_buf(),
            pins: versions_pins(&t),
            env,
        })
    }

    pub fn from_env(root: &Path) -> Result<Ctx, String> {
        Ctx::new(root, Box::new(|k| std::env::var(k).ok()))
    }

    fn var(&self, k: &str) -> Option<String> {
        (self.env)(k).filter(|v| !v.is_empty())
    }

    /// A setting from the environment (empty is unset).
    pub fn setting(&self, k: &str) -> Option<String> {
        self.var(k)
    }

    /// `${NAME}` from versions.env, `${NAME:-default}` from the environment.
    pub fn expand(&self, s: &str, part: Option<&str>) -> Result<String, String> {
        let mut out = String::new();
        let mut rest = s;
        while let Some(i) = rest.find("${") {
            out.push_str(&rest[..i]);
            let after = &rest[i + 2..];
            let j = after
                .find('}')
                .ok_or_else(|| format!("unclosed ${{ in '{s}'"))?;
            let inner = &after[..j];
            let v = if let Some((k, d)) = inner.split_once(":-") {
                self.var(k).unwrap_or_else(|| d.to_owned())
            } else {
                self.pins
                    .get(inner)
                    .cloned()
                    .ok_or_else(|| format!("{inner} is not a pin of versions.env (in '{s}')"))?
            };
            out.push_str(&v);
            rest = &after[j + 1..];
        }
        out.push_str(rest);
        Ok(match part {
            Some(p) => out.replace("{part}", p),
            None => out,
        })
    }

    /// The downloads folder: OG_DOWNLOADS, or GENOME_DIR/downloads.
    pub fn downloads(&self) -> Result<PathBuf, String> {
        if let Some(d) = self.var("OG_DOWNLOADS") {
            return Ok(PathBuf::from(d));
        }
        self.genome_dir().map(|g| g.join("downloads"))
    }

    pub fn genome_dir(&self) -> Result<PathBuf, String> {
        self.var("GENOME_DIR").map(PathBuf::from).ok_or_else(|| {
            "ERROR: GENOME_DIR is not set (the folder setup installs into)".to_owned()
        })
    }

    fn part_of<'a>(&self, r: &Source, part: Option<&'a str>) -> Result<Option<&'a str>, String> {
        match (r.parts.is_empty(), part) {
            (true, None) => Ok(None),
            (true, Some(_)) => Err(format!("ERROR: {} has no parts", r.id)),
            (false, None) => Err(format!(
                "ERROR: {} needs --part (one of: {})",
                r.id,
                r.parts.join(" ")
            )),
            (false, Some(p)) if r.parts.iter().any(|x| x == p) => Ok(Some(p)),
            (false, Some(p)) => Err(format!("ERROR: {p} is not a part of {}", r.id)),
        }
    }

    /// Where row R (part PART) is stored.
    pub fn path(&self, r: &Source, part: Option<&str>) -> Result<PathBuf, String> {
        let part = self.part_of(r, part)?;
        let p = self
            .expand(&r.path, part)
            .map_err(|e| format!("ERROR: {}: {e}", r.id))?;
        match p.strip_prefix("downloads/") {
            Some(rest) => Ok(self.downloads()?.join(rest)),
            None => Ok(self.genome_dir()?.join(p)),
        }
    }

    pub fn url(&self, r: &Source, part: Option<&str>) -> Result<String, String> {
        let part = self.part_of(r, part)?;
        self.expand(&r.url, part)
            .map_err(|e| format!("ERROR: {}: {e}", r.id))
    }

    /// The check value, `@ID` read from that row's stored file.
    pub fn value(&self, s: &Sources, r: &Source, part: Option<&str>) -> Result<String, String> {
        let part = self.part_of(r, part)?;
        if let Some(o) = r.value.strip_prefix('@') {
            let p = self.path(s.get(o)?, None)?;
            return Ok(format!("file://{}", p.display()));
        }
        self.expand(&r.value, part)
            .map_err(|e| format!("ERROR: {}: {e}", r.id))
    }
}

/// Every pin of the table written out in it (not `${...}` or `@...`): its
/// literal URLs (up to the first `${` or `{part}`) and checksums.
pub fn literal_pins(s: &Sources) -> Vec<(String, String)> {
    let mut out = vec![];
    for r in s.rows.values() {
        let url = r
            .url
            .split("${")
            .next()
            .unwrap_or("")
            .split("{part}")
            .next()
            .unwrap_or("");
        if url.len() > "https://".len() + 3 && r.url == url {
            out.push((r.id.clone(), url.to_owned()));
        }
        if !r.value.is_empty()
            && r.value != "-"
            && !r.value.contains("${")
            && !r.value.starts_with('@')
            && !is_url(&r.value)
        {
            out.push((r.id.clone(), r.value.clone()));
        }
    }
    out
}

// --- the command ----------------------------------------------------------

const USAGE: &str =
    "usage: og setup fetch URL DEST [md5|sha256|sum|gzip-ends|crc32c VALUE|CHECKSUM_URL]
       og setup fetch --get URL OUT|-
       og setup fetch --source ID [--part P] [--root DIR]
       og setup source ID [--part P] [--root DIR] path|url|check|value|bytes
       og setup path [--root DIR]
";

/// The checkout: --root, or the folder above `bin/` of this binary's
/// launcher (bin/.og/og-*), or the current directory's.
pub fn root_of(arg: Option<&str>) -> PathBuf {
    if let Some(r) = arg {
        return PathBuf::from(r);
    }
    if let Some(r) = std::env::var_os("PGP_ROOT").filter(|v| !v.is_empty()) {
        return PathBuf::from(r);
    }
    if let Ok(exe) = std::env::current_exe() {
        // bin/.og/og-linux-x86_64 -> the checkout.
        if let Some(r) = exe.parent().and_then(Path::parent).and_then(Path::parent)
            && r.join(TABLE).is_file()
        {
            return r.to_path_buf();
        }
    }
    PathBuf::from(".")
}

/// Takes `--name VALUE` out of ARGS.
pub fn take(args: &mut Vec<String>, name: &str) -> Result<Option<String>, String> {
    match args.iter().position(|a| a == name) {
        None => Ok(None),
        Some(i) if i + 1 < args.len() => {
            let v = args.remove(i + 1);
            args.remove(i);
            Ok(Some(v))
        }
        Some(_) => Err(format!("{name} needs a value")),
    }
}

/// Takes the flag NAME out of ARGS.
pub fn flag(args: &mut Vec<String>, name: &str) -> bool {
    match args.iter().position(|a| a == name) {
        Some(i) => {
            args.remove(i);
            true
        }
        None => false,
    }
}

/// `og setup fetch`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let mut a = args.to_vec();
    let knobs = match Knobs::from_env() {
        Ok(k) => k,
        Err(m) => {
            let _ = writeln!(err, "{m}");
            return 1;
        }
    };
    if flag(&mut a, "--get") {
        let [url, dest] = a.as_slice() else {
            let _ = write!(err, "{USAGE}");
            return 2;
        };
        let target = if dest == "-" {
            Out::Stdout
        } else {
            Out::File(Path::new(dest))
        };
        return match get(url, target, &knobs) {
            Ok((ok, body)) => match out.write_all(&body).and_then(|()| out.flush()) {
                Ok(()) => u8::from(!ok),
                Err(e) => {
                    let _ = writeln!(err, "ERROR: cannot write the download of {url}: {e}");
                    1
                }
            },
            Err(m) => {
                let _ = writeln!(err, "{m}");
                1
            }
        };
    }
    let source = match (
        take(&mut a, "--source"),
        take(&mut a, "--part"),
        take(&mut a, "--root"),
    ) {
        (Ok(s), Ok(p), Ok(r)) => s.map(|s| (s, p, r)),
        (Err(m), _, _) | (_, Err(m), _) | (_, _, Err(m)) => {
            let _ = writeln!(err, "og setup fetch: {m}");
            return 2;
        }
    };
    if let Some((id, part, root)) = source {
        if !a.is_empty() {
            let _ = write!(err, "{USAGE}");
            return 2;
        }
        let res = (|| -> Result<(PathBuf, String, Kind, String), String> {
            let ctx = Ctx::from_env(&root_of(root.as_deref()))?;
            let s = Sources::read(&ctx.root)?;
            let r = s.get(&id)?;
            let part = part.as_deref();
            Ok((
                ctx.path(r, part)?,
                ctx.url(r, part)?,
                Kind::parse(&r.check)?,
                ctx.value(&s, r, part)?,
            ))
        })();
        let (dest, url, kind, value) = match res {
            Ok(x) => x,
            Err(m) => {
                let _ = writeln!(err, "{m}");
                return 1;
            }
        };
        if url.is_empty() || url == "-" {
            let _ = writeln!(
                err,
                "ERROR: {id} has no URL of its own: its installer builds it"
            );
            return 1;
        }
        let label = match &part {
            Some(p) => format!("{id} {p}"),
            None => id.clone(),
        };
        let r = Fetch {
            url: &url,
            dest: &dest,
            kind,
            want: &value,
            on_wrong: OnWrong::Refuse,
            source: Some(&label),
        };
        return match fetch(&r, &knobs, err) {
            Ok(_) => {
                let _ = writeln!(out, "{}", dest.display());
                0
            }
            Err(()) => 1,
        };
    }
    let (url, dest, kind, want) = match a.as_slice() {
        [u, d] => (u, d, "", ""),
        [u, d, k, w] => (u, d, k.as_str(), w.as_str()),
        [u, d, k] => (u, d, k.as_str(), ""),
        _ => {
            let _ = write!(err, "{USAGE}");
            return 2;
        }
    };
    let kind = match Kind::parse(kind) {
        Ok(k) => k,
        Err(m) => {
            let _ = writeln!(err, "{m}");
            return 1;
        }
    };
    let dest = Path::new(dest);
    let on_wrong = if in_downloads(dest, &|k| std::env::var(k).ok()) {
        OnWrong::Refuse
    } else {
        OnWrong::Replace
    };
    let r = Fetch {
        url,
        dest,
        kind,
        want,
        on_wrong,
        source: None,
    };
    u8::from(fetch(&r, &knobs, err).is_err())
}

/// `og setup source ID [--part P] FIELD`: one field of a row, expanded.
pub fn source_cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let mut a = args.to_vec();
    let (part, root) = match (take(&mut a, "--part"), take(&mut a, "--root")) {
        (Ok(p), Ok(r)) => (p, r),
        (Err(m), _) | (_, Err(m)) => {
            let _ = writeln!(err, "og setup source: {m}");
            return 2;
        }
    };
    let [id, field] = a.as_slice() else {
        let _ = write!(err, "{USAGE}");
        return 2;
    };
    let res = (|| -> Result<String, String> {
        let ctx = Ctx::from_env(&root_of(root.as_deref()))?;
        let s = Sources::read(&ctx.root)?;
        let r = s.get(id)?;
        let part = part.as_deref();
        Ok(match field.as_str() {
            "path" => ctx.path(r, part)?.display().to_string(),
            "url" => ctx.url(r, part)?,
            "check" => r.check.clone(),
            "value" => ctx.value(&s, r, part)?,
            "bytes" => r.bytes.map(|b| b.to_string()).unwrap_or_default(),
            other => {
                return Err(format!(
                    "og setup source: no field '{other}' (path, url, check, value, bytes)"
                ));
            }
        })
    })();
    match res {
        Ok(v) => {
            let _ = writeln!(out, "{v}");
            0
        }
        Err(m) => {
            let _ = writeln!(err, "{m}");
            1
        }
    }
}

/// `og setup path`: the downloads folder.
pub fn path_cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let mut a = args.to_vec();
    let root = match take(&mut a, "--root") {
        Ok(r) => r,
        Err(m) => {
            let _ = writeln!(err, "og setup path: {m}");
            return 2;
        }
    };
    if !a.is_empty() {
        let _ = write!(err, "{USAGE}");
        return 2;
    }
    match Ctx::from_env(&root_of(root.as_deref())).and_then(|c| c.downloads()) {
        Ok(d) => {
            let _ = writeln!(out, "{}", d.display());
            0
        }
        Err(m) => {
            let _ = writeln!(err, "{m}");
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crc32c_matches_the_standard_check_value() {
        let mut c = Crc32c::default();
        c.update(b"123456789");
        assert_eq!(c.value(), 0xe306_9283);
        assert_eq!(base64(&[0xe3, 0x06, 0x92, 0x83]), "4waSgw==");
    }

    #[test]
    fn base64_pads() {
        assert_eq!(base64(b"f"), "Zg==");
        assert_eq!(base64(b"fo"), "Zm8=");
        assert_eq!(base64(b"foo"), "Zm9v");
    }

    #[test]
    fn checksum_files_are_read_as_the_shell_read_them() {
        assert_eq!(checksum_in("abc  x.gz\n", "other"), Some("abc".into()));
        let many = "1 2 a/x.gz\n3 4 b/y.gz\n\n";
        assert_eq!(checksum_in(many, "y.gz"), Some("3".into()));
        assert_eq!(checksum_in(many, "z.gz"), None);
    }

    #[test]
    fn knobs_refuse_what_the_shell_refused() {
        let env = |t: &'static str, w: &'static str, s: &'static str| {
            move |k: &str| match k {
                "FETCH_TRIES" => Some(t.to_owned()),
                "FETCH_WAIT" => Some(w.to_owned()),
                "FETCH_STALL" => Some(s.to_owned()),
                _ => None,
            }
        };
        assert!(
            Knobs::from(&env("x", "", ""))
                .unwrap_err()
                .contains("FETCH_TRIES must be")
        );
        assert!(Knobs::from(&env("0", "", "")).is_err());
        assert!(
            Knobs::from(&env("2", "1s", ""))
                .unwrap_err()
                .contains("FETCH_WAIT must be")
        );
        assert!(
            Knobs::from(&env("2", "", "00"))
                .unwrap_err()
                .contains("FETCH_STALL must be")
        );
        let k = Knobs::from(&env("08", "", "")).unwrap();
        assert_eq!((k.tries, k.stall, k.wait), (8, 120, None));
    }

    #[test]
    fn parts_expand() {
        assert_eq!(expand_parts("1..3 X").unwrap(), ["1", "2", "3", "X"]);
        assert!(expand_parts("3..1").is_err());
    }
}
