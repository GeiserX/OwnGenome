//! `og view serve`: the local genome viewer. It serves igv.js and one
//! sample's reads and calls on 127.0.0.1, and nothing else from the disk:
//!
//! ```text
//!   /                 index.html of --assets (assets/viewer)
//!   /viewer.js        the page's script, from --assets
//!   /igv.min.js       the pinned igv.js release owngenome-view.sh checked by sha256
//!   /config.json      the reference, the tracks and the start locus (build_config)
//!   /data/<n>/<name>  the n-th file of the config, its name percent-encoded in
//!                     the URL: the reference FASTA and its .fai, the BAM or CRAM
//!                     and its index, the VCF and the SV VCF and their .tbi
//! ```
//!
//! `scripts/owngenome-view.sh` finds the sample's files and runs it. The
//! server binds 127.0.0.1 only (there is no bind option), answers byte-range
//! requests (igv.js reads a BAM, CRAM, FASTA or indexed VCF in pieces), and
//! sends a Content-Security-Policy that keeps the page on this server, so
//! igv.js's own calls to outside services are blocked by the browser.
//!
//! It has no login and needs none on loopback, but every account on the
//! machine can reach 127.0.0.1. Two requests are refused with 403: a `Host`
//! header that is not 127.0.0.1, localhost or [::1] (a web page reaching it
//! through a DNS name pointed at 127.0.0.1), and a request a proxy forwarded
//! (a `Forwarded`, `X-Forwarded-*` or `X-Real-IP` header): a reverse proxy
//! would hand the sample to whoever reaches the proxy, and the viewer has no
//! password to stop them. An ssh port forward adds no such header, so the
//! documented way to reach a remote viewer keeps working.
//!
//! The live site's server reuses [`sample_files`], [`build_config`],
//! [`parse_range`] and the policies here for its per-run viewer.

use indexmap::IndexMap;
use pycompat::json::{self, Dump, Value};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::net::{Shutdown, TcpListener, TcpStream};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, mpsc};
use std::thread::{self, JoinHandle};
use std::time::Duration;

/// The only address the viewer listens on.
pub const BIND: &str = "127.0.0.1";
/// The port when none is given; the report's Viewer links use it.
pub const DEFAULT_PORT: u16 = 8765;
/// The policy of every answer: the page talks to this server only.
pub const CSP: &str = "default-src 'none'; script-src 'self'; connect-src 'self'; img-src 'self' data: blob:; \
     style-src 'self' 'unsafe-inline'; font-src 'self' data:; worker-src 'self' blob:; \
     base-uri 'none'; form-action 'none'; frame-ancestors 'none'";
/// The live site's viewer: the same policy, but a page of the same site may
/// frame it (the genome map's reads panel).
pub const LIVE_CSP: &str = "default-src 'none'; script-src 'self'; connect-src 'self'; img-src 'self' data: blob:; \
     style-src 'self' 'unsafe-inline'; font-src 'self' data:; worker-src 'self' blob:; \
     base-uri 'none'; form-action 'none'; frame-ancestors 'self'";
/// The `Host` names a request may carry (its port aside).
const ALLOWED_HOSTS: [&str; 3] = ["127.0.0.1", "localhost", "[::1]"];
/// Headers a proxy adds when it forwards a request (lower case).
pub const FORWARDED_HEADERS: [&str; 5] = [
    "forwarded",
    "x-forwarded-for",
    "x-forwarded-host",
    "x-forwarded-proto",
    "x-real-ip",
];
/// Why a forwarded request is refused.
const FORWARDED_TEXT: &str = "the viewer has no password, so it does not answer a request a proxy \
     forwarded: reach it over an ssh port forward instead (docs/viewer.md)";
/// The longest request head read, in bytes.
const MAX_HEAD: usize = 64 * 1024;
/// How much of a file one write sends.
const CHUNK: usize = 1 << 20;

const USAGE: &str = "usage: og view serve --sample NAME --fasta FASTA --igv-js FILE --assets DIR\n\
     \x20                    [--reads BAM|CRAM] [--vcf VCF] [--sv-vcf VCF] [--locus LOCUS] [--port N]\n";

/// The contig names and lengths of a `.fai`, in order.
pub fn read_fai(fai: &Path) -> Result<Vec<(String, u64)>, String> {
    let text =
        fs::read_to_string(fai).map_err(|e| format!("cannot read {}: {e}", fai.display()))?;
    let mut out = Vec::new();
    for line in text.lines() {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() >= 2 && !p[0].is_empty() {
            let len = p[1].trim().parse::<u64>().map_err(|_| {
                format!(
                    "{}: '{}' is not the length of {}",
                    fai.display(),
                    p[1],
                    p[0]
                )
            })?;
            out.push((p[0].to_owned(), len));
        }
    }
    if out.is_empty() {
        return Err(format!("{} names no contig", fai.display()));
    }
    Ok(out)
}

/// `chr`, `chr:pos` or `chr:start-end` split into its parts; `None` when
/// the text has another shape. A position is digits with optional commas,
/// starting with a digit.
fn split_locus(locus: &str) -> Option<(&str, Option<&str>, Option<&str>)> {
    let name_end = locus
        .find(|c: char| !(c.is_ascii_alphanumeric() || matches!(c, '.' | '_' | '-')))
        .unwrap_or(locus.len());
    let (name, rest) = locus.split_at(name_end);
    if name.is_empty() {
        return None;
    }
    let pos = |s: &str| {
        s.starts_with(|c: char| c.is_ascii_digit())
            && s.chars().all(|c| c.is_ascii_digit() || c == ',')
    };
    if rest.is_empty() {
        return Some((name, None, None));
    }
    let rest = rest.strip_prefix(':')?;
    let (start, end) = match rest.split_once('-') {
        Some((s, e)) => (s, Some(e)),
        None => (rest, None),
    };
    if !pos(start) || end.is_some_and(|e| !pos(e)) {
        return None;
    }
    Some((name, Some(start), end))
}

/// A position as written (commas allowed); `None` past `u64`.
fn position(text: &str) -> Option<u64> {
    text.replace(',', "").parse().ok()
}

/// The locus as igv.js reads it (`chr`, `chr:pos` or `chr:start-end`),
/// checked against the reference: the contig exists, each position lies on
/// it, and the end does not come before the start. The first contig when
/// none is given.
pub fn check_locus(locus: Option<&str>, contigs: &[(String, u64)]) -> Result<String, String> {
    let Some(locus) = locus.filter(|l| !l.is_empty()) else {
        return Ok(contigs[0].0.clone());
    };
    let Some((name, start, end)) = split_locus(locus) else {
        return Err(format!(
            "'{locus}' is not a locus: give chr1, chr1:12345 or chr1:12000-13000"
        ));
    };
    let Some(&(_, len)) = contigs.iter().find(|(n, _)| n == name) else {
        return Err(format!("'{name}' is not a contig of the reference"));
    };
    let mut at = Vec::new();
    for x in [start, end].into_iter().flatten() {
        match position(x) {
            Some(p) if (1..=len).contains(&p) => at.push(p),
            _ => return Err(format!("{x} is outside {name} (1 to {len})")),
        }
    }
    if let [s, e] = at[..]
        && e < s
    {
        return Err(format!(
            "'{locus}' ends before it starts: give the smaller position first, as in {name}:{e}-{s}"
        ));
    }
    Ok(locus.to_owned())
}

/// The first of these files that exists.
fn index_of(names: &[PathBuf]) -> Option<PathBuf> {
    names.iter().find(|n| n.is_file()).cloned()
}

/// `path` with `suffix` added to its name.
fn plus(path: &Path, suffix: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

/// The reads, the small-variant VCF and the SV VCF of a sample folder, by
/// the names the steps write; `None` for one that is not there.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SampleFiles {
    pub reads: Option<PathBuf>,
    pub vcf: Option<PathBuf>,
    pub sv_vcf: Option<PathBuf>,
}

/// The files of `sample_dir`, as `scripts/owngenome-view.sh` picks them
/// with none of `READS`, `VCF`, `SV_VCF`, `ALIGN_DIR` or `VCF_DIR` set: the
/// BAM, else the CRAM, of `aligned/`; the VCF of `vcf/`; and Manta's
/// diploid SVs (step 4), else step 22's consensus. Manta comes first because
/// the consensus keeps only the calls two callers agree on, and the reads
/// are where a call one caller made gets checked.
pub fn sample_files(sample_dir: &Path, sample: &str) -> SampleFiles {
    let first = |names: &[String]| {
        names
            .iter()
            .map(|n| sample_dir.join(n))
            .find(|p| p.is_file())
    };
    SampleFiles {
        reads: first(&[
            format!("aligned/{sample}_sorted.bam"),
            format!("aligned/{sample}_sorted.cram"),
        ]),
        vcf: first(&[format!("vcf/{sample}.vcf.gz")]),
        sv_vcf: first(&[
            "manta/results/variants/diploidSV.vcf.gz".to_owned(),
            format!("sv_merged/{sample}_sv_consensus.vcf.gz"),
        ]),
    }
}

/// `urllib.parse.quote` of a file name: every byte but letters, digits and
/// `_.-~/` as `%XX`.
pub fn quote(name: &str) -> String {
    let mut out = String::new();
    for b in name.bytes() {
        if b.is_ascii_alphanumeric() || matches!(b, b'_' | b'.' | b'-' | b'~' | b'/') {
            out.push(b as char);
        } else {
            out.push_str(&format!("%{b:02X}"));
        }
    }
    out
}

/// `urllib.parse.unquote`: each `%XX` as its byte, read as UTF-8 with
/// replacement; a `%` without two hex digits stays as it is.
pub fn unquote(text: &str) -> String {
    let b = text.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        let hex = |c: u8| (c as char).to_digit(16);
        if b[i] == b'%'
            && i + 2 < b.len()
            && let (Some(h), Some(l)) = (hex(b[i + 1]), hex(b[i + 2]))
        {
            out.push((h * 16 + l) as u8);
            i += 3;
        } else {
            out.push(b[i]);
            i += 1;
        }
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// What the viewer shows: the igv.js options (`config.json`) and the files
/// the server may send, `data/<index>/<name>` each.
#[derive(Debug, Clone)]
pub struct Config {
    pub json: Value,
    pub files: Vec<PathBuf>,
    pub locus: String,
    pub sample: String,
}

impl Config {
    /// `config.json` as the page reads it.
    pub fn to_json(&self) -> String {
        let opts = Dump {
            indent: Some(1),
            ..Dump::default()
        };
        json::dumps(&self.json, &opts).expect("the config holds no float")
    }
}

fn file_name(path: &Path) -> String {
    path.file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default()
}

fn obj(pairs: Vec<(&str, Value)>) -> Value {
    Value::Object(
        pairs
            .into_iter()
            .map(|(k, v)| (k.to_owned(), v))
            .collect::<IndexMap<_, _>>(),
    )
}

fn s(text: &str) -> Value {
    Value::Str(text.to_owned())
}

/// The config for `sample`: the reference (`fasta`, with its `.fai`
/// beside it), the reads (with their index, required), the VCF and the SV
/// VCF (read in pieces with a `.tbi` or `.csi`, whole without), and the
/// start locus ([`check_locus`]).
pub fn build_config(
    sample: &str,
    fasta: &Path,
    reads: Option<&Path>,
    vcf: Option<&Path>,
    sv_vcf: Option<&Path>,
    locus: Option<&str>,
) -> Result<Config, String> {
    let mut files: Vec<PathBuf> = Vec::new();
    let mut url = |path: &Path| {
        files.push(std::path::absolute(path).unwrap_or_else(|_| path.to_owned()));
        Value::Str(format!(
            "data/{}/{}",
            files.len() - 1,
            quote(&file_name(path))
        ))
    };
    let fai = plus(fasta, ".fai");
    if !fasta.is_file() || !fai.is_file() {
        return Err(format!(
            "the reference needs its FASTA and .fai: {}",
            fasta.display()
        ));
    }
    let contigs = read_fai(&fai)?;
    let locus = check_locus(locus, &contigs)?;
    let reference = obj(vec![
        ("id", s("reference")),
        ("name", s(&file_name(fasta))),
        ("fastaURL", url(fasta)),
        ("indexURL", url(&fai)),
    ]);
    let mut tracks = Vec::new();
    if let Some(reads) = reads {
        let cram = reads.to_string_lossy().ends_with(".cram");
        let stem = reads.with_extension("");
        let idx = index_of(&if cram {
            vec![plus(reads, ".crai"), plus(&stem, ".crai")]
        } else {
            vec![
                plus(reads, ".bai"),
                plus(&stem, ".bai"),
                plus(reads, ".csi"),
            ]
        })
        .ok_or_else(|| {
            format!(
                "no index beside {}: run samtools index on it",
                reads.display()
            )
        })?;
        tracks.push(obj(vec![
            ("name", Value::Str(format!("{sample} reads"))),
            ("type", s("alignment")),
            ("format", s(if cram { "cram" } else { "bam" })),
            ("url", url(reads)),
            ("indexURL", url(&idx)),
        ]));
    }
    for (path, name) in [(vcf, "small variants"), (sv_vcf, "structural variants")] {
        let Some(path) = path else { continue };
        let mut t = vec![
            ("name", Value::Str(format!("{sample} {name}"))),
            ("type", s("variant")),
            ("format", s("vcf")),
            ("url", url(path)),
        ];
        match index_of(&[plus(path, ".tbi"), plus(path, ".csi")]) {
            Some(idx) => t.push(("indexURL", url(&idx))),
            None => t.push(("indexed", Value::Bool(false))),
        }
        tracks.push(obj(t));
    }
    let json = obj(vec![
        ("sample", s(sample)),
        ("locus", s(&locus)),
        ("reference", reference),
        ("tracks", Value::Array(tracks)),
    ]);
    Ok(Config {
        json,
        files,
        locus,
        sample: sample.to_owned(),
    })
}

/// What a `Range` header asks of a file of `size` bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Range {
    /// No range, or one this server answers with the whole file.
    Whole,
    /// Bytes `start..=end`.
    Part(u64, u64),
    /// A range past the end: 416.
    Bad,
}

/// One `bytes=start-end`, `bytes=start-` or `bytes=-suffix` range. A list
/// of ranges, or anything else, is answered with the whole file.
pub fn parse_range(header: Option<&str>, size: u64) -> Range {
    let Some(spec) = header.map(str::trim).and_then(|h| h.strip_prefix("bytes=")) else {
        return Range::Whole;
    };
    let Some((a, b)) = spec.split_once('-') else {
        return Range::Whole;
    };
    let digits = |t: &str| t.chars().all(|c| c.is_ascii_digit());
    if !digits(a) || !digits(b) || (a.is_empty() && b.is_empty()) {
        return Range::Whole;
    }
    // Digits that overflow u64 ask for more than any file holds.
    let num = |t: &str| t.parse::<u64>().unwrap_or(u64::MAX);
    let (start, end) = if a.is_empty() {
        (size.saturating_sub(num(b)), size.saturating_sub(1))
    } else {
        let start = num(a);
        (
            start,
            if b.is_empty() {
                size.saturating_sub(1)
            } else {
                num(b)
            },
        )
    };
    if start >= size || end < start {
        return Range::Bad;
    }
    Range::Part(start, end.min(size - 1))
}

/// Where a route's answer comes from.
enum Source {
    Bytes(Vec<u8>),
    File(PathBuf),
}

/// One logged request: `"GET / HTTP/1.1" 200 -`, as Python's
/// http.server writes it.
pub type Logger = Box<dyn Fn(&str) + Send + Sync>;

struct Shared {
    routes: HashMap<String, (&'static str, Source)>,
    log: Option<Logger>,
}

/// The routes of `cfg`: the page and its script from `assets`, igv.js,
/// the config and its files.
fn routes(cfg: &Config, assets: &Path, igv_js: &Path) -> HashMap<String, (&'static str, Source)> {
    let html = "text/html; charset=utf-8";
    let js = "text/javascript; charset=utf-8";
    let mut r = HashMap::new();
    r.insert(
        "/".to_owned(),
        (html, Source::File(assets.join("index.html"))),
    );
    r.insert(
        "/index.html".to_owned(),
        (html, Source::File(assets.join("index.html"))),
    );
    r.insert(
        "/viewer.js".to_owned(),
        (js, Source::File(assets.join("viewer.js"))),
    );
    r.insert(
        "/igv.min.js".to_owned(),
        (js, Source::File(igv_js.to_owned())),
    );
    r.insert(
        "/config.json".to_owned(),
        (
            "application/json",
            Source::Bytes(cfg.to_json().into_bytes()),
        ),
    );
    for (i, f) in cfg.files.iter().enumerate() {
        r.insert(
            format!("/data/{i}/{}", file_name(f)),
            ("application/octet-stream", Source::File(f.clone())),
        );
    }
    r
}

/// A parsed request head.
struct Request {
    line: String,
    method: String,
    target: String,
    headers: Vec<(String, String)>,
}

impl Request {
    fn header(&self, name: &str) -> Option<&str> {
        self.headers
            .iter()
            .find(|(k, _)| k == name)
            .map(|(_, v)| v.as_str())
    }
}

/// Reads the request line and headers; `None` for anything that is not an
/// HTTP request head within [`MAX_HEAD`] bytes.
fn read_head(stream: &TcpStream) -> Option<Request> {
    let mut r = BufReader::new(stream.take(MAX_HEAD as u64));
    let mut line = String::new();
    r.read_line(&mut line).ok()?;
    let line = line.trim_end_matches(['\r', '\n']).to_owned();
    let mut parts = line.split(' ');
    let (method, target, version) = (parts.next()?, parts.next()?, parts.next()?);
    if parts.next().is_some() || !version.starts_with("HTTP/") || method.is_empty() {
        return None;
    }
    let mut headers = Vec::new();
    loop {
        let mut h = String::new();
        if r.read_line(&mut h).ok()? == 0 {
            return None;
        }
        let h = h.trim_end_matches(['\r', '\n']);
        if h.is_empty() {
            break;
        }
        let (k, v) = h.split_once(':')?;
        headers.push((k.trim().to_ascii_lowercase(), v.trim().to_owned()));
    }
    Some(Request {
        method: method.to_owned(),
        target: target.to_owned(),
        line,
        headers,
    })
}

/// The status line's reason.
fn reason(code: u16) -> &'static str {
    match code {
        200 => "OK",
        206 => "Partial Content",
        400 => "Bad Request",
        403 => "Forbidden",
        404 => "Not Found",
        405 => "Method Not Allowed",
        416 => "Range Not Satisfiable",
        _ => "",
    }
}

/// Writes the status line and headers: the viewer's policy and the
/// no-store, nosniff, no-referrer and same-origin headers on every answer.
fn head(
    out: &mut dyn Write,
    code: u16,
    ctype: &str,
    length: u64,
    extra: &[(&str, String)],
) -> io::Result<()> {
    let mut h = format!(
        "HTTP/1.1 {code} {}\r\nServer: owngenome-view\r\nContent-Type: {ctype}\r\nContent-Length: {length}\r\n\
         Cache-Control: no-store\r\nContent-Security-Policy: {CSP}\r\nX-Content-Type-Options: nosniff\r\n\
         Referrer-Policy: no-referrer\r\nCross-Origin-Resource-Policy: same-origin\r\nConnection: close\r\n",
        reason(code)
    );
    for (k, v) in extra {
        h.push_str(&format!("{k}: {v}\r\n"));
    }
    h.push_str("\r\n");
    out.write_all(h.as_bytes())
}

/// Sends `src`, or the one byte range `range` asks for (206, or 416 past
/// the end), [`CHUNK`] bytes at a time. Returns the status sent.
pub fn send_file(
    out: &mut dyn Write,
    is_head: bool,
    ctype: &str,
    src: &Path,
    range: Option<&str>,
) -> io::Result<u16> {
    let mut f = File::open(src)?;
    let size = f.metadata()?.len();
    let rng = parse_range(range, size);
    if rng == Range::Bad {
        head(
            out,
            416,
            "text/plain; charset=utf-8",
            0,
            &[("Content-Range", format!("bytes */{size}"))],
        )?;
        return Ok(416);
    }
    let (start, end) = match rng {
        Range::Part(a, b) => (a, b),
        _ => (0, size.saturating_sub(1)),
    };
    let length = if size == 0 { 0 } else { end - start + 1 };
    let mut extra = vec![("Accept-Ranges", "bytes".to_owned())];
    if rng != Range::Whole {
        extra.push(("Content-Range", format!("bytes {start}-{end}/{size}")));
    }
    let code = if rng == Range::Whole { 200 } else { 206 };
    head(out, code, ctype, length, &extra)?;
    if is_head || length == 0 {
        return Ok(code);
    }
    f.seek(SeekFrom::Start(start))?;
    let mut left = length;
    let mut buf = vec![0u8; CHUNK.min(length as usize)];
    while left > 0 {
        let want = buf.len().min(left as usize);
        let n = f.read(&mut buf[..want])?;
        if n == 0 {
            break;
        }
        out.write_all(&buf[..n])?;
        left -= n as u64;
    }
    Ok(code)
}

/// A short text answer.
fn fail(out: &mut dyn Write, is_head: bool, code: u16, text: &str) -> io::Result<u16> {
    let body = format!("{text}\n");
    let extra: Vec<(&str, String)> = if code == 405 {
        vec![("Allow", "GET, HEAD".to_owned())]
    } else {
        Vec::new()
    };
    head(
        out,
        code,
        "text/plain; charset=utf-8",
        body.len() as u64,
        &extra,
    )?;
    if !is_head {
        out.write_all(body.as_bytes())?;
    }
    Ok(code)
}

/// The `Host` header without its port, as the allow-list compares it.
fn host_name(host: &str) -> &str {
    match host.rfind(':') {
        Some(i) if i + 1 < host.len() && host[i + 1..].bytes().all(|b| b.is_ascii_digit()) => {
            &host[..i]
        }
        _ => host,
    }
}

/// Answers one request on `stream`.
fn answer(shared: &Shared, stream: &mut TcpStream) -> io::Result<()> {
    let Some(req) = read_head(stream) else {
        let _ = fail(stream, false, 400, "bad request");
        return Ok(());
    };
    let is_head = req.method == "HEAD";
    let code = if req.method != "GET" && !is_head {
        fail(stream, false, 405, "only GET and HEAD")
    } else if !ALLOWED_HOSTS.contains(&host_name(req.header("host").unwrap_or(""))) {
        fail(
            stream,
            is_head,
            403,
            "the Host header must be 127.0.0.1 or localhost",
        )
    } else if FORWARDED_HEADERS.iter().any(|h| req.header(h).is_some()) {
        fail(stream, is_head, 403, FORWARDED_TEXT)
    } else {
        // The config percent-encodes each file name (a space, '#' or '?' in
        // it), so the browser asks for the encoded form; the routes hold
        // the plain name.
        let path = req.target.split(['?', '#']).next().unwrap_or("");
        match shared.routes.get(&unquote(path)) {
            None => fail(stream, is_head, 404, "not found"),
            Some((ctype, Source::Bytes(data))) => head(stream, 200, ctype, data.len() as u64, &[])
                .and_then(|()| {
                    if is_head {
                        Ok(())
                    } else {
                        stream.write_all(data)
                    }
                })
                .map(|()| 200),
            Some((ctype, Source::File(src))) => {
                match send_file(stream, is_head, ctype, src, req.header("range")) {
                    Err(e) if e.kind() == io::ErrorKind::NotFound => {
                        fail(stream, is_head, 404, "not found")
                    }
                    other => other,
                }
            }
        }
    };
    if let (Some(log), Ok(code)) = (&shared.log, &code) {
        log(&format!("\"{}\" {code} -", req.line));
    }
    code.map(|_| ())
}

fn handle(shared: &Shared, mut stream: TcpStream) {
    let _ = stream.set_read_timeout(Some(Duration::from_secs(30)));
    let _ = stream.set_write_timeout(Some(Duration::from_secs(60)));
    // A client that goes away mid-answer (igv.js cancels ranges it no
    // longer needs) is not an error of the server.
    let _ = answer(shared, &mut stream);
    let _ = stream.flush();
    let _ = stream.shutdown(Shutdown::Both);
}

/// Accepts connections until `stop` is set, one thread each.
fn serve(listener: TcpListener, shared: Arc<Shared>, stop: Arc<AtomicBool>) {
    for conn in listener.incoming() {
        if stop.load(Ordering::SeqCst) {
            break;
        }
        let Ok(stream) = conn else { continue };
        let shared = Arc::clone(&shared);
        thread::spawn(move || handle(&shared, stream));
    }
}

/// The viewer running in a thread; dropping it stops it.
pub struct Server {
    port: u16,
    stop: Arc<AtomicBool>,
    thread: Option<JoinHandle<()>>,
}

impl Server {
    /// Starts the viewer for `cfg` on 127.0.0.1:`port` (0 picks a free
    /// port), with the page from `assets` and igv.js from `igv_js`; each
    /// request goes to `log` when one is given.
    pub fn start(
        cfg: &Config,
        assets: &Path,
        igv_js: &Path,
        port: u16,
        log: Option<Logger>,
    ) -> io::Result<Server> {
        let listener = TcpListener::bind((BIND, port))?;
        let port = listener.local_addr()?.port();
        let shared = Arc::new(Shared {
            routes: routes(cfg, assets, igv_js),
            log,
        });
        let stop = Arc::new(AtomicBool::new(false));
        let s = Arc::clone(&stop);
        let thread = thread::spawn(move || serve(listener, shared, s));
        Ok(Server {
            port,
            stop,
            thread: Some(thread),
        })
    }

    /// The port it listens on.
    pub fn port(&self) -> u16 {
        self.port
    }
}

impl Drop for Server {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::SeqCst);
        // Wake the accept loop so it sees the flag.
        let _ = TcpStream::connect((BIND, self.port));
        if let Some(t) = self.thread.take() {
            let _ = t.join();
        }
    }
}

/// The options of `og view serve`.
#[derive(Debug, Default)]
struct Args {
    sample: Option<String>,
    fasta: Option<PathBuf>,
    reads: Option<PathBuf>,
    vcf: Option<PathBuf>,
    sv_vcf: Option<PathBuf>,
    locus: Option<String>,
    igv_js: Option<PathBuf>,
    assets: Option<PathBuf>,
    port: Option<String>,
}

fn parse_args(args: &[String]) -> Result<Args, String> {
    let mut a = Args::default();
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        let (flag, inline) = match arg.split_once('=') {
            Some((f, v)) if f.starts_with("--") => (f, Some(v.to_owned())),
            _ => (arg.as_str(), None),
        };
        let mut value = || -> Result<String, String> {
            match &inline {
                Some(v) => Ok(v.clone()),
                None => it
                    .next()
                    .cloned()
                    .ok_or_else(|| format!("argument {flag}: expected one argument")),
            }
        };
        let path = |v: String| Some(PathBuf::from(v));
        match flag {
            "--sample" => a.sample = Some(value()?),
            "--fasta" => a.fasta = path(value()?),
            "--reads" => a.reads = path(value()?),
            "--vcf" => a.vcf = path(value()?),
            "--sv-vcf" => a.sv_vcf = path(value()?),
            "--locus" => a.locus = Some(value()?),
            "--igv-js" => a.igv_js = path(value()?),
            "--assets" => a.assets = path(value()?),
            "--port" => a.port = Some(value()?),
            _ => return Err(format!("unrecognized arguments: {arg}")),
        }
    }
    Ok(a)
}

/// `og view serve ARGS`: build the config, print the address and serve
/// until stopped. Exit 2 on a usage error or a config the reference
/// refuses, 1 when the port cannot be bound.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if args.iter().any(|a| a == "-h" || a == "--help") {
        let _ = write!(out, "{USAGE}");
        return 0;
    }
    let usage = |err: &mut dyn Write, why: &str| {
        let _ = writeln!(err, "{USAGE}og view serve: error: {why}");
        2
    };
    let a = match parse_args(args) {
        Ok(a) => a,
        Err(why) => return usage(err, &why),
    };
    let (Some(sample), Some(fasta), Some(igv_js), Some(assets)) =
        (&a.sample, &a.fasta, &a.igv_js, &a.assets)
    else {
        return usage(
            err,
            "the following arguments are required: --sample, --fasta, --igv-js, --assets",
        );
    };
    let port = match a.port.as_deref().map(str::parse::<u16>) {
        None => DEFAULT_PORT,
        Some(Ok(p)) => p,
        Some(Err(_)) => {
            return usage(
                err,
                &format!(
                    "argument --port: '{}' is not a port number",
                    a.port.as_deref().unwrap_or_default()
                ),
            );
        }
    };
    let cfg = match build_config(
        sample,
        fasta,
        a.reads.as_deref(),
        a.vcf.as_deref(),
        a.sv_vcf.as_deref(),
        a.locus.as_deref(),
    ) {
        Ok(c) => c,
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            return 2;
        }
    };
    for (what, p) in [
        ("igv.min.js", igv_js.clone()),
        ("the viewer page", assets.join("index.html")),
        ("the viewer script", assets.join("viewer.js")),
    ] {
        if !p.is_file() {
            let _ = writeln!(err, "ERROR: {what} not found: {}", p.display());
            return 2;
        }
    }
    let listener = match TcpListener::bind((BIND, port)) {
        Ok(l) => l,
        Err(e) => {
            let _ = writeln!(
                err,
                "ERROR: cannot listen on {BIND}:{port}: {e}. Set OWNGENOME_VIEW_PORT to a free port."
            );
            return 1;
        }
    };
    let bound = listener.local_addr().map_or(port, |a| a.port());
    let _ = writeln!(
        out,
        "Viewer for {}: http://{BIND}:{bound}/#{}",
        cfg.sample, cfg.locus
    );
    let _ = writeln!(out, "Press Ctrl-C to stop it.");
    let _ = out.flush();
    // The request log goes through `err`, which the caller may hold locked
    // (og's main holds stderr's lock while a verb runs). A handler thread
    // that wrote to stderr itself would wait on that lock forever, before it
    // closes its connection. So the handlers send their lines here, and
    // this thread, the one that owns `err`, writes them.
    let (tx, rx) = mpsc::channel::<String>();
    let shared = Arc::new(Shared {
        routes: routes(&cfg, assets, igv_js),
        log: Some(Box::new(move |line: &str| {
            let _ = tx.send(line.to_owned());
        })),
    });
    let server = thread::spawn(move || serve(listener, shared, Arc::new(AtomicBool::new(false))));
    for line in rx {
        let _ = writeln!(err, "{line}");
        let _ = err.flush();
    }
    let _ = server.join();
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn contigs() -> Vec<(String, u64)> {
        vec![("chr1".into(), 3000), ("chr2".into(), 1000)]
    }

    #[test]
    fn a_locus_is_checked_against_the_reference() {
        let c = contigs();
        assert_eq!(check_locus(None, &c).unwrap(), "chr1");
        assert_eq!(check_locus(Some(""), &c).unwrap(), "chr1");
        for ok in [
            "chr2",
            "chr1:451",
            "chr1:1,000-2,000",
            "chr1:5-5",
            "chr1:1-3000",
        ] {
            assert_eq!(check_locus(Some(ok), &c).unwrap(), ok);
        }
        for bad in [
            "chr9:100",
            "chr1:99999",
            "chr1;alert(1)",
            "BRCA1 x",
            "chr1:0",
            "chr1:",
            "chr1:-5",
            "chr1:5-",
            ":5",
            "chr1:99999999999999999999999",
        ] {
            assert!(check_locus(Some(bad), &c).is_err(), "{bad} was accepted");
        }
    }

    #[test]
    fn a_reversed_locus_is_refused() {
        let e = check_locus(Some("chr1:2,000-1000"), &contigs()).unwrap_err();
        assert!(
            e.contains("ends before it starts") && e.contains("chr1:1000-2000"),
            "{e}"
        );
    }

    #[test]
    fn ranges_as_http_reads_them() {
        assert_eq!(parse_range(None, 100), Range::Whole);
        assert_eq!(parse_range(Some("bytes=10-19"), 100), Range::Part(10, 19));
        assert_eq!(parse_range(Some(" bytes=90- "), 100), Range::Part(90, 99));
        assert_eq!(parse_range(Some("bytes=-28"), 100), Range::Part(72, 99));
        assert_eq!(parse_range(Some("bytes=-500"), 100), Range::Part(0, 99));
        assert_eq!(parse_range(Some("bytes=50-5000"), 100), Range::Part(50, 99));
        assert_eq!(parse_range(Some("bytes=100-"), 100), Range::Bad);
        assert_eq!(parse_range(Some("bytes=20-10"), 100), Range::Bad);
        assert_eq!(parse_range(Some("bytes=-"), 100), Range::Whole);
        assert_eq!(parse_range(Some("bytes=1-2,5-6"), 100), Range::Whole);
        assert_eq!(parse_range(Some("items=1-2"), 100), Range::Whole);
    }

    #[test]
    fn names_are_quoted_and_unquoted_as_urllib_does() {
        assert_eq!(quote("my ref#1?.fasta"), "my%20ref%231%3F.fasta");
        assert_eq!(quote("año.bam"), "a%C3%B1o.bam");
        assert_eq!(unquote("my%20ref%231%3F.fasta"), "my ref#1?.fasta");
        assert_eq!(unquote("100%"), "100%");
        assert_eq!(unquote("%zz%4"), "%zz%4");
        assert_eq!(unquote("..%2F..%2Fx"), "../../x");
    }

    #[test]
    fn the_host_port_is_stripped_before_the_allow_list() {
        assert_eq!(host_name("127.0.0.1:8765"), "127.0.0.1");
        assert_eq!(host_name("[::1]:8765"), "[::1]");
        assert_eq!(host_name("[::1]"), "[::1]");
        assert_eq!(host_name("localhost"), "localhost");
        assert_eq!(host_name("evil.example:80"), "evil.example");
    }
}
