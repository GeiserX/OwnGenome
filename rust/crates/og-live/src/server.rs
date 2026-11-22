//! `og live serve`: the live site's server. It shows each run under the
//! genome folder while it goes: the step map, the event stream, the report
//! and its files, the family pages, the genome viewer, the map and the
//! Prometheus metrics.
//!
//! ```text
//! GET /                        the people under --root, by family (assets/live/home.html)
//! GET /runs.json               each run's progress and each family's members
//! GET /assets/<file>, /assets/plugins/<file>
//! GET /runs/<run>/             a person's page (index.html), or a family's (family.html)
//! GET /runs/<run>/state.json   the run's state, the headline of each finished
//!                              section and what the plugins add
//! GET /runs/<run>/events       the stored messages as Server-Sent Events
//!                              (Last-Event-ID resumes)
//! GET /runs/<run>/report.html  the run's <run>_report.html
//! GET /runs/<run>/<relpath>    a file that report links to by a relative href, or
//!                              one a plugin names; nothing else from disk
//! ```
//!
//! and the plugins' routes ([`crate::plugins`]). The plugins are compiled
//! in: there is no plugin folder.
//!
//! Not here yet: the password, the sessions, TLS, the proxy options and the
//! ingest listener. Until they land, `og live serve` listens on loopback
//! only, takes no password, and refuses the flags that would make a site
//! less safe if they were ignored (exit 2). It reads the messages another
//! process stores under `STATE_DIR/runs/<run>/events.jsonl` (the Python
//! server's ingest listener today), so both can run beside each other on one
//! state folder. The access log (`STATE_DIR/access.log`, TSV) records who
//! opened what, never a query string, a request body or a result.

use crate::csp::{self, CSP, DENY, RUN_HTML_CSP};
use crate::family::{inside, realpath};
use crate::ops::OpsConfig;
use crate::plugins::{self, ResultsFns};
use crate::runs::{self, Run};
use crate::state::{self, Live, dumps, lock};
use crate::static_files as sf;
use crate::viewer::unquote;
use std::collections::BTreeMap;
use std::fs::{self, OpenOptions};
use std::io::{self, BufRead, BufReader, Read, Write};
use std::net::{IpAddr, Shutdown, SocketAddr, TcpListener, TcpStream, ToSocketAddrs};
use std::path::{Path, PathBuf};
use std::sync::atomic::Ordering;
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};
use std::time::Duration;

/// How long a client may take to send its request (an event stream writes
/// every [`HEARTBEAT`]).
pub const REQUEST_SECONDS: u64 = 60;
/// How often an open event stream writes when nothing happens.
pub const HEARTBEAT: Duration = Duration::from_secs(15);
/// The largest request head read.
const MAX_HEAD: usize = 64 * 1024;
/// The cookie the login sets (`og live serve` sets none yet; logout clears
/// it).
pub const COOKIE: &str = "owngenome_live";

const USAGE: &str =
    "usage: og live serve --root GENOME_DIR [--state-dir DIR] [--bind ADDR] [--port N]
                     [--max-streams N] [--access-log FILE] [--access-log-max-mb N]
                     [--config FILE] [--quiet]
                     [--viewer-fasta FILE --viewer-igv-js FILE]
                     [--disk-warn-percent N] [--mem-warn-percent N] [--stall-minutes N]
                     [--metrics-token-file FILE] [--eta-box NAME] [--eta-cores N]
";

/// What the site is started with.
#[derive(Clone)]
pub struct Config {
    /// GENOME_DIR: the folder holding one folder per run.
    pub root: PathBuf,
    /// Where the messages are read and the access log and map cache are
    /// written (`ROOT/.owngenome-live`).
    pub state_dir: PathBuf,
    pub bind: String,
    pub port: u16,
    /// `assets/live`.
    pub assets: PathBuf,
    /// `assets/`: the viewer's page, the map's and the summary's tables.
    pub repo_assets: PathBuf,
    /// `steps.tsv`.
    pub steps: PathBuf,
    /// Event streams open at once.
    pub max_streams: usize,
    /// `STATE_DIR/access.log`.
    pub access_log: PathBuf,
    /// Then it moves to `access.log.1`; 0 never.
    pub access_log_max_mb: f64,
    /// No request line on stderr.
    pub quiet: bool,
    pub ops: OpsConfig,
    /// The bearer token Prometheus sends to `/metrics` (with a password).
    pub metrics_token_file: Option<PathBuf>,
    /// The genome viewer's reference FASTA, with its `.fai` beside it.
    pub viewer_fasta: Option<PathBuf>,
    /// The checked igv.min.js.
    pub viewer_igv_js: Option<PathBuf>,
    /// The og the map runs (`og map build`, `og map genes`): this binary.
    pub og: PathBuf,
    /// The producer for results, findings and the map. The CLI configures it.
    pub results: Option<ResultsFns>,
    /// How often a quiet event stream writes.
    pub heartbeat: Duration,
}

/// The first ancestor of this binary that holds `rel`.
fn beside_exe(rel: &str) -> Option<PathBuf> {
    let exe = std::env::current_exe().ok()?;
    exe.ancestors()
        .skip(1)
        .find(|d| d.join(rel).exists())
        .map(Path::to_path_buf)
}

impl Config {
    /// The defaults for genome folder `root`: assets and steps.tsv of the
    /// checkout this binary sits in.
    pub fn new(root: impl Into<PathBuf>) -> Config {
        let root = root.into();
        let repo = beside_exe("assets/live/index.html").unwrap_or_default();
        let state_dir = root.join(".owngenome-live");
        Config {
            access_log: state_dir.join("access.log"),
            state_dir,
            root,
            bind: "127.0.0.1".into(),
            port: 8787,
            assets: repo.join("assets").join("live"),
            repo_assets: repo.join("assets"),
            steps: og_runlog::steps::default_path().unwrap_or_else(|| repo.join("steps.tsv")),
            max_streams: 32,
            access_log_max_mb: 10.0,
            quiet: false,
            ops: OpsConfig::default(),
            metrics_token_file: None,
            viewer_fasta: None,
            viewer_igv_js: None,
            og: std::env::current_exe().unwrap_or_else(|_| PathBuf::from("og")),
            results: None,
            heartbeat: HEARTBEAT,
        }
    }
}

// --- the access log -----------------------------------------------------------------------

/// `STATE_DIR/access.log`: one TSV row per request that opens something.
pub struct AccessLog {
    path: PathBuf,
    max_bytes: u64,
    lock: Mutex<()>,
}

/// The UTC time now as `2026-11-20T10:00:00Z`.
pub fn utc_now() -> String {
    jiff::Timestamp::now()
        .strftime("%Y-%m-%dT%H:%M:%SZ")
        .to_string()
}

impl AccessLog {
    pub const HEADER: &str = "utc\tclient\tsession8\taction\trun\tsection\tstatus\n";

    pub fn new(path: &Path, max_mb: f64) -> io::Result<AccessLog> {
        if let Some(d) = std::path::absolute(path)?.parent() {
            fs::create_dir_all(d)?;
        }
        Ok(AccessLog {
            path: path.to_owned(),
            max_bytes: (max_mb * 1024.0 * 1024.0) as u64,
            lock: Mutex::new(()),
        })
    }

    /// One row: the time, then each field with tabs and line breaks as
    /// spaces, cut at 200 characters, `-` for none.
    pub fn write(&self, fields: &[&str]) {
        let mut row = utc_now();
        for f in fields {
            let f = if f.is_empty() { "-" } else { f };
            let clean: String = f
                .chars()
                .map(|c| {
                    if matches!(c, '\t' | '\r' | '\n') {
                        ' '
                    } else {
                        c
                    }
                })
                .take(200)
                .collect();
            row.push('\t');
            row.push_str(&clean);
        }
        row.push('\n');
        let _g = lock(&self.lock);
        if self.max_bytes > 0 && fs::metadata(&self.path).is_ok_and(|m| m.len() > self.max_bytes) {
            let mut old = self.path.clone().into_os_string();
            old.push(".1");
            let _ = fs::rename(&self.path, old);
        }
        let new = !self.path.exists();
        if let Ok(mut f) = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&self.path)
        {
            let text = if new {
                format!("{}{row}", Self::HEADER)
            } else {
                row
            };
            let _ = f.write_all(text.as_bytes());
        }
    }
}

// --- one request --------------------------------------------------------------------------

/// A request and the connection its answer goes to.
pub struct Req<'a> {
    pub live: &'a Live,
    /// `GET` (also for a HEAD, then [`Req::head_only`]) or `POST`.
    pub method: String,
    pub head_only: bool,
    /// The path, still quoted, without its query.
    pub path: String,
    pub query: String,
    headers: Vec<(String, String)>,
    pub client: String,
    out: &'a mut dyn Write,
}

/// The status line's reason.
fn reason(code: u16) -> &'static str {
    match code {
        200 => "OK",
        206 => "Partial Content",
        301 => "Moved Permanently",
        302 => "Found",
        303 => "See Other",
        400 => "Bad Request",
        401 => "Unauthorized",
        403 => "Forbidden",
        404 => "Not Found",
        409 => "Conflict",
        413 => "Request Entity Too Large",
        416 => "Requested Range Not Satisfiable",
        500 => "Internal Server Error",
        501 => "Not Implemented",
        503 => "Service Unavailable",
        _ => "",
    }
}

/// The text of a refusal without its own.
fn deny_text(code: u16) -> &'static str {
    match code {
        401 => "login required",
        403 => "forbidden",
        404 => "not found",
        429 => "too many failed logins",
        _ => "error",
    }
}

/// `urllib.parse.quote(text)`: letters, digits and `_.-~/` as they are.
pub fn quote(text: &str) -> String {
    crate::viewer::quote(text)
}

impl Req<'_> {
    /// Header `name` (any case), the first one.
    pub fn header(&self, name: &str) -> Option<&str> {
        let n = name.to_ascii_lowercase();
        self.headers
            .iter()
            .find(|(k, _)| *k == n)
            .map(|(_, v)| v.as_str())
    }

    /// The status line and headers: the type and length, `extra`, then the
    /// policy `csp`, no sniffing, no referrer, no cache and frame option
    /// `frame`.
    pub fn head(
        &mut self,
        code: u16,
        ctype: &str,
        length: Option<u64>,
        extra: &[(&str, String)],
        csp: &str,
        frame: &str,
    ) -> io::Result<()> {
        let date = jiff::Timestamp::now().strftime("%a, %d %b %Y %H:%M:%S GMT");
        let mut h = format!(
            "HTTP/1.1 {code} {}\r\nServer: owngenome-live\r\nDate: {date}\r\nContent-Type: {ctype}\r\n",
            reason(code)
        );
        if let Some(n) = length {
            h.push_str(&format!("Content-Length: {n}\r\n"));
        }
        for (k, v) in extra {
            h.push_str(&format!("{k}: {v}\r\n"));
        }
        for (k, v) in csp::common(csp, frame) {
            h.push_str(&format!("{k}: {v}\r\n"));
        }
        h.push_str("Connection: close\r\n\r\n");
        self.out.write_all(h.as_bytes())
    }

    /// The bytes after the head, unless the request is a HEAD.
    pub fn body(&mut self, body: &[u8]) -> io::Result<()> {
        if !self.head_only {
            self.out.write_all(body)?;
        }
        self.out.flush()
    }

    /// A whole answer under policy `csp`.
    pub fn send(
        &mut self,
        code: u16,
        ctype: &str,
        body: &[u8],
        extra: &[(&str, String)],
        csp: &str,
    ) -> io::Result<()> {
        self.head(code, ctype, Some(body.len() as u64), extra, csp, DENY)?;
        self.body(body)
    }

    /// A whole answer under the site's policy.
    pub fn send_body(&mut self, code: u16, ctype: &str, body: impl AsRef<[u8]>) -> io::Result<()> {
        self.send(code, ctype, body.as_ref(), &[], CSP)
    }

    /// An empty redirect.
    pub fn redirect(&mut self, code: u16, location: &str) -> io::Result<()> {
        self.send(
            code,
            "text/plain",
            b"",
            &[("Location", location.to_owned())],
            CSP,
        )
    }

    /// `json.dumps(v)`.
    pub fn send_json(&mut self, v: &pycompat::json::Value) -> io::Result<()> {
        self.send_body(200, "application/json", dumps(v))
    }

    /// A refusal: `text`, else the code's own line.
    pub fn deny(&mut self, code: u16, text: Option<&str>) -> io::Result<()> {
        let t = format!("{}\n", text.unwrap_or_else(|| deny_text(code)));
        self.send_body(code, "text/plain; charset=utf-8", t)
    }

    /// One access-log row for this request.
    pub fn log(&self, action: &str, run: &str, section: &str, status: u16) {
        self.live
            .access
            .write(&[&self.client, "-", action, run, section, &status.to_string()]);
    }

    /// The folder of run `name`, when the site knows it and it has one.
    pub fn run_dir(&self, name: &str) -> Option<PathBuf> {
        self.live.run(name).and_then(|r| r.dir)
    }

    /// The writer for an answer streamed by hand (the event stream, a file
    /// by ranges).
    pub fn out(&mut self) -> &mut dyn Write {
        &mut *self.out
    }
}

// --- the dispatch -------------------------------------------------------------------------

/// Whether the Host header names a host the site answers to: loopback or
/// the bind address, with any port.
fn host_ok(req: &Req<'_>) -> bool {
    let h = req.header("host").unwrap_or("").trim().to_lowercase();
    let name = if let Some(r) = h.strip_prefix('[') {
        match r.split_once(']') {
            Some((inner, rest)) if port_ok(rest) => inner.to_owned(),
            _ => return false,
        }
    } else {
        let (n, rest) = match h.find(':') {
            Some(i) => (&h[..i], &h[i..]),
            None => (h.as_str(), ""),
        };
        if n.is_empty() || !port_ok(rest) {
            return false;
        }
        n.to_owned()
    };
    req.live
        .allowed
        .iter()
        .any(|a| a.trim_matches(|c| c == '[' || c == ']') == name)
}

/// `""` or `:` and digits.
fn port_ok(rest: &str) -> bool {
    rest.is_empty()
        || rest
            .strip_prefix(':')
            .is_some_and(|p| !p.is_empty() && p.bytes().all(|b| b.is_ascii_digit()))
}

/// A POST from another site is refused. Browsers say where a request comes
/// from in Sec-Fetch-Site; with Referrer-Policy no-referrer they send
/// `Origin: null` for the site's own forms.
fn origin_ok(req: &Req<'_>) -> bool {
    if let Some(site) = req.header("sec-fetch-site").filter(|s| !s.is_empty()) {
        return site == "same-origin" || site == "none";
    }
    let o = req.header("origin").unwrap_or("");
    if o.is_empty() || o == "null" {
        return true;
    }
    netloc(o).to_lowercase() == req.header("host").unwrap_or("").to_lowercase()
}

/// The network location of a URL (`urlsplit(url).netloc`).
fn netloc(url: &str) -> &str {
    let Some(i) = url.find("//") else { return "" };
    let scheme_ok = url[..i]
        .strip_suffix(':')
        .is_some_and(|s| s.is_empty() || s.chars().next().is_some_and(|c| c.is_ascii_alphabetic()))
        || i == 0;
    if !scheme_ok {
        return "";
    }
    let rest = &url[i + 2..];
    &rest[..rest.find(['/', '?', '#']).unwrap_or(rest.len())]
}

/// Answers one request.
pub fn handle(req: &mut Req<'_>) -> io::Result<()> {
    if !host_ok(req) {
        req.log("host-refused", "-", "-", 403);
        return req.deny(403, None);
    }
    if req.method == "POST" && !origin_ok(req) {
        req.log("origin-refused", "-", "-", 403);
        return req.deny(403, None);
    }
    let path = req.path.clone();
    if path == "/login" {
        return login(req);
    }
    if path == "/logout" && req.method == "POST" {
        req.log("logout", "-", "-", 303);
        let clear = format!("{COOKIE}=; Path=/; HttpOnly; SameSite=Strict; Max-Age=0");
        return req.send(
            303,
            "text/plain",
            b"",
            &[("Location", "/login".into()), ("Set-Cookie", clear)],
            CSP,
        );
    }
    if req.method == "GET"
        && let Some(name) = path.strip_prefix("/assets/")
    {
        return asset(req, name);
    }
    if plugins::route(req, &path)? {
        return Ok(());
    }
    if req.method != "GET" {
        return req.deny(404, None);
    }
    if path == "/" {
        return index(req);
    }
    if path == "/runs.json" {
        req.log("runs", "-", "-", 200);
        let v = req.live.runs_index();
        return req.send_json(&v);
    }
    let Some(tail) = path.strip_prefix("/runs/") else {
        return req.deny(404, None);
    };
    let (seg, rest) = match tail.find('/') {
        Some(i) => (&tail[..i], Some(&tail[i + 1..])),
        None => (tail, None),
    };
    if seg.is_empty() {
        return req.deny(404, None);
    }
    let Some(run) = req.live.run(&unquote(seg)) else {
        req.log("run", "-", "-", 404);
        return req.deny(404, None);
    };
    match rest {
        None => {
            let to = format!("/runs/{}/", run.name);
            req.redirect(301, &to)
        }
        Some(rest) => run_path(req, &run, rest),
    }
}

/// `/login`: the login page; a POST goes home (no password is set).
fn login(req: &mut Req<'_>) -> io::Result<()> {
    if req.method == "POST" {
        return req.redirect(303, "/");
    }
    let next = sf::query_value(&req.query, "next").unwrap_or_else(|| "/".into());
    match sf::login_page(&req.live.cfg.assets, "", &next) {
        Ok(page) => req.send_body(200, "text/html; charset=utf-8", page),
        Err(_) => req.deny(500, None),
    }
}

fn asset(req: &mut Req<'_>, name: &str) -> io::Result<()> {
    let top = realpath(&req.live.cfg.assets);
    match sf::asset(&top, name) {
        Some((ctype, body)) => req.send_body(200, ctype, body),
        None => req.deny(404, None),
    }
}

/// The home page: the people grouped by family, drawn by home.js from
/// `/runs.json`; the plain list in it is what a browser without scripts
/// shows.
fn index(req: &mut Req<'_>) -> io::Result<()> {
    let runs: Vec<(String, String)> = req
        .live
        .runs()
        .into_values()
        .map(|r| (r.name, r.kind))
        .collect();
    req.log("index", "-", "-", 200);
    match sf::home_page(&req.live.cfg.assets, &runs) {
        Ok(page) => req.send_body(200, "text/html; charset=utf-8", page),
        Err(_) => req.deny(500, None),
    }
}

fn run_path(req: &mut Req<'_>, run: &Run, rest: &str) -> io::Result<()> {
    let name = run.name.as_str();
    match rest {
        "" => {
            let page = sf::run_page(&req.live.cfg.assets, name, run.is_family());
            req.log("run", name, "-", 200);
            return match page {
                Ok(p) => req.send_body(200, "text/html; charset=utf-8", p),
                Err(_) => req.deny(500, None),
            };
        }
        "state.json" => {
            req.log("state", name, "-", 200);
            let v = req.live.state(run);
            return req.send_json(&v);
        }
        "events" => return events(req, run),
        _ => {}
    }
    let root = req.live.root.clone();
    let d = run.dir.as_deref();
    if rest == "report.html" {
        let p = d.map(|d| d.join(format!("{name}_report.html")));
        let ok = matches!((&p, d), (Some(p), Some(d)) if p.is_file() && inside(p, &realpath(d)));
        let page = p.filter(|_| ok).and_then(|p| fs::read(p).ok());
        let Some(page) = page else {
            req.log("report", name, "-", 404);
            return req.deny(404, None);
        };
        let page =
            plugins::report_filter(req.live, name, String::from_utf8_lossy(&page).into_owned());
        req.log("report", name, "-", 200);
        return req.send(
            200,
            "text/html; charset=utf-8",
            page.as_bytes(),
            &[],
            RUN_HTML_CSP,
        );
    }
    let rel = sf::normpath(&unquote(rest));
    let Some(d) = d.filter(|_| {
        !rel.starts_with("..") && !rel.starts_with('/') && runs::file_allowed(run, &root, &rel)
    }) else {
        req.log("file", name, &rel, 404);
        return req.deny(404, None);
    };
    let p = d.join(&rel);
    let top = realpath(d);
    if p.is_dir()
        && inside(&p, &top)
        && inside(&p, &root)
        && let Some(page) = runs::sole_html(run, &root, &rel, &runs::run_files(name, Some(d)))
    {
        req.log("file", name, &rel, 302);
        let to = format!("/runs/{}/{}", quote(name), quote(&format!("{rel}/{page}")));
        return req.redirect(302, &to);
    }
    if !p.is_file() || !inside(&p, &top) || !inside(&p, &root) {
        req.log("file", name, &rel, 404);
        return req.deny(404, None);
    }
    let Ok(body) = fs::read(&p) else {
        req.log("file", name, &rel, 404);
        return req.deny(404, None);
    };
    let ext = sf::ext(&rel).to_lowercase();
    req.log("file", name, &rel, 200);
    let csp = if [".html", ".htm", ".svg"].contains(&ext.as_str()) {
        RUN_HTML_CSP
    } else {
        CSP
    };
    req.send(200, sf::type_or_bytes(&ext), &body, &[], csp)
}

/// Gives back an event stream's place when it ends.
struct StreamSlot<'a>(&'a Live);

impl Drop for StreamSlot<'_> {
    fn drop(&mut self) {
        self.0.streams.fetch_sub(1, Ordering::SeqCst);
    }
}

/// `/runs/<run>/events`: the run's stored messages as Server-Sent Events,
/// from the one after `Last-Event-ID`, then each new one as it is stored; a
/// comment line every heartbeat while nothing comes.
fn events(req: &mut Req<'_>, run: &Run) -> io::Result<()> {
    let live = req.live;
    let name = run.name.as_str();
    let max = live.max_streams.load(Ordering::SeqCst);
    if live
        .streams
        .fetch_update(Ordering::SeqCst, Ordering::SeqCst, |n| {
            (n < max).then_some(n + 1)
        })
        .is_err()
    {
        req.log("events", name, "-", 503);
        return req.deny(503, Some("too many open event streams"));
    }
    let _slot = StreamSlot(live);
    let after: usize = req
        .header("last-event-id")
        .and_then(|v| v.trim().parse::<i64>().ok())
        .map_or(0, |n| n.max(0) as usize);
    req.log("events", name, "-", 200);
    req.head(
        200,
        "text/event-stream",
        None,
        &[("X-Accel-Buffering", "no".into())],
        CSP,
        DENY,
    )?;
    if req.head_only {
        return req.out().flush();
    }
    let out = req.out();
    out.write_all(b"retry: 3000\n\n")?;
    out.flush()?;
    let mut sent = after;
    while !live.stopping.load(Ordering::SeqCst) {
        let new = live
            .store
            .wait(name, sent, live.cfg.heartbeat, &live.stopping);
        if new.is_empty() {
            out.write_all(b": heartbeat\n\n")?;
        } else {
            let mut text = String::new();
            for (i, payload) in state::messages(&new, &live.steps, sent) {
                text.push_str(&format!(
                    "id: {i}\ndata: {}\n\n",
                    state::dumps_compact(&payload)
                ));
            }
            out.write_all(text.as_bytes())?;
            sent += new.len();
        }
        out.flush()?;
    }
    Ok(())
}

// --- the connection -----------------------------------------------------------------------

/// A request's headers: (name lowered, value), in order.
type Headers = Vec<(String, String)>;

/// The head of one request: method, target and headers (names lowered);
/// `None` for one that cannot be read.
fn read_head(r: &mut impl BufRead) -> Option<(String, String, Headers)> {
    let mut head = Vec::new();
    loop {
        let mut line = Vec::new();
        let n = r
            .by_ref()
            .take((MAX_HEAD - head.len()) as u64 + 1)
            .read_until(b'\n', &mut line)
            .ok()?;
        if n == 0 {
            if head.is_empty() {
                return None;
            }
            break;
        }
        head.extend_from_slice(&line);
        if head.len() > MAX_HEAD {
            return None;
        }
        if line == b"\r\n" || line == b"\n" {
            break;
        }
    }
    let text = String::from_utf8_lossy(&head).into_owned();
    let mut lines = text.split('\n').map(|l| l.trim_end_matches('\r'));
    let first = lines.next()?;
    let mut parts = first.split(' ');
    let (method, target, version) = (parts.next()?, parts.next()?, parts.next()?);
    if parts.next().is_some() || !version.starts_with("HTTP/") {
        return None;
    }
    let mut headers = Vec::new();
    for l in lines {
        if l.is_empty() {
            break;
        }
        let (k, v) = l.split_once(':')?;
        headers.push((k.trim().to_ascii_lowercase(), v.trim().to_owned()));
    }
    Some((method.to_owned(), target.to_owned(), headers))
}

fn serve_conn(live: &Live, stream: TcpStream) {
    let timeout = Some(Duration::from_secs(REQUEST_SECONDS));
    let _ = stream.set_read_timeout(timeout);
    let _ = stream.set_write_timeout(timeout);
    let client = stream
        .peer_addr()
        .map(|a| a.ip().to_string())
        .unwrap_or_else(|_| "-".into());
    let Ok(read_half) = stream.try_clone() else {
        return;
    };
    let mut reader = BufReader::new(read_half);
    let mut out = io::BufWriter::new(&stream);
    let Some((method, target, headers)) = read_head(&mut reader) else {
        let _ = out.write_all(
            b"HTTP/1.1 400 Bad Request\r\nContent-Type: text/plain; charset=utf-8\r\nContent-Length: 12\r\n\
              Connection: close\r\n\r\nbad request\n",
        );
        return;
    };
    let (path, query) = match target.split_once('?') {
        Some((p, q)) => (p, q),
        None => (target.as_str(), ""),
    };
    let path = path.split('#').next().unwrap_or("").to_owned();
    if !live.cfg.quiet {
        eprintln!("og live serve: {method} {path}");
    }
    let (m, head_only) = match method.as_str() {
        "GET" => ("GET", false),
        "HEAD" => ("GET", true),
        "POST" => ("POST", false),
        _ => {
            let _ = out.write_all(
                b"HTTP/1.1 501 Not Implemented\r\nContent-Type: text/plain; charset=utf-8\r\nContent-Length: 20\r\n\
                  Connection: close\r\n\r\nunsupported method\r\n",
            );
            return;
        }
    };
    let mut req = Req {
        live,
        method: m.into(),
        head_only,
        path,
        query: query.split('#').next().unwrap_or("").to_owned(),
        headers,
        client,
        out: &mut out,
    };
    let _ = handle(&mut req);
    let _ = out.flush();
    drop(out);
    let _ = stream.shutdown(Shutdown::Write);
    // No route reads a body: what the client still sends is read and
    // dropped once the answer is out, so closing does not reset the
    // connection before the client has read it.
    let _ = stream.set_read_timeout(Some(Duration::from_secs(1)));
    let _ = io::copy(&mut reader.take(MAX_HEAD as u64), &mut io::sink());
}

/// A running site: [`Server::start`] binds and serves on threads, dropping
/// it stops it.
pub struct Server {
    live: Arc<Live>,
    addr: SocketAddr,
    accept: Option<JoinHandle<()>>,
}

impl Server {
    /// Reads steps.tsv, opens the access log, binds `cfg.bind:cfg.port` (0:
    /// any free port) and serves.
    pub fn start(cfg: Config) -> Result<Server, String> {
        let steps = og_runlog::steps::read_steps(&cfg.steps).map_err(|e| format!("steps: {e}"))?;
        fs::create_dir_all(&cfg.state_dir)
            .map_err(|e| format!("--state-dir {}: {e}", cfg.state_dir.display()))?;
        let access = AccessLog::new(&cfg.access_log, cfg.access_log_max_mb)
            .map_err(|e| format!("--access-log {}: {e}", cfg.access_log.display()))?;
        let host = cfg.bind.trim_matches(|c| c == '[' || c == ']').to_owned();
        let listener = TcpListener::bind((host.as_str(), cfg.port))
            .map_err(|e| format!("cannot listen on {}:{}: {e}", cfg.bind, cfg.port))?;
        let addr = listener.local_addr().map_err(|e| e.to_string())?;
        let live = Arc::new(Live::new(cfg, steps, access));
        let l2 = Arc::clone(&live);
        let accept = thread::spawn(move || {
            for conn in listener.incoming() {
                if l2.stopping.load(Ordering::SeqCst) {
                    break;
                }
                let Ok(stream) = conn else { continue };
                let l3 = Arc::clone(&l2);
                thread::spawn(move || serve_conn(&l3, stream));
            }
        });
        Ok(Server {
            live,
            addr,
            accept: Some(accept),
        })
    }

    /// The port it listens on.
    pub fn port(&self) -> u16 {
        self.addr.port()
    }

    /// The state it serves.
    pub fn live(&self) -> &Live {
        &self.live
    }

    /// Serves until the process ends.
    pub fn wait(mut self) {
        if let Some(h) = self.accept.take() {
            let _ = h.join();
        }
    }
}

impl Drop for Server {
    fn drop(&mut self) {
        self.live.stopping.store(true, Ordering::SeqCst);
        // Wakes the accept loop so it sees the flag.
        let _ = TcpStream::connect_timeout(&self.addr, Duration::from_secs(1));
        if let Some(h) = self.accept.take() {
            let _ = h.join();
        }
    }
}

// --- the command line ---------------------------------------------------------------------

/// Whether `addr` is a loopback address, or a name that resolves only to
/// loopback.
pub fn is_loopback(addr: &str) -> bool {
    if addr == "localhost" {
        return true;
    }
    let bare = addr.trim_matches(|c| c == '[' || c == ']');
    if let Ok(ip) = bare.parse::<IpAddr>() {
        return ip.is_loopback();
    }
    match (bare, 0).to_socket_addrs() {
        Ok(it) => {
            let all: Vec<SocketAddr> = it.collect();
            !all.is_empty() && all.iter().all(|a| a.ip().is_loopback())
        }
        Err(_) => false,
    }
}

/// The flags that take no value.
const SWITCHES: [&str; 2] = ["quiet", "cookie_secure"];

/// The flags of a site with a password, TLS or a proxy: refused until they
/// are ported, because a site that ignored them would be less safe than
/// its owner asked for.
const NOT_YET: [&str; 6] = [
    "password_file",
    "allowed_host",
    "trust_proxy",
    "tls_cert",
    "tls_key",
    "cookie_secure",
];

/// Flags read and not used yet: the ingest listener's and the sessions'.
/// `og live serve` starts neither, so they change nothing.
const UNUSED: [&str; 6] = [
    "ingest_bind",
    "ingest_port",
    "max_body_kb",
    "session_hours",
    "session_idle_minutes",
    "login_max_failures",
];

/// Every flag, by its name with `_`.
const FLAGS: [&str; 26] = [
    "root",
    "state_dir",
    "bind",
    "port",
    "ingest_bind",
    "ingest_port",
    "password_file",
    "allowed_host",
    "trust_proxy",
    "tls_cert",
    "tls_key",
    "cookie_secure",
    "session_hours",
    "session_idle_minutes",
    "max_streams",
    "max_body_kb",
    "login_max_failures",
    "access_log",
    "access_log_max_mb",
    "assets",
    "steps",
    "quiet",
    "viewer_fasta",
    "viewer_igv_js",
    "metrics_token_file",
    "plugin_dir",
];

/// The ops plugin's flags.
const OPS_FLAGS: [&str; 5] = [
    "disk_warn_percent",
    "mem_warn_percent",
    "stall_minutes",
    "eta_box",
    "eta_cores",
];

fn known(k: &str) -> bool {
    FLAGS.contains(&k) || OPS_FLAGS.contains(&k)
}

/// The flags as given: each name (with `_`) and its values, in order.
type Flags = BTreeMap<String, Vec<String>>;

/// `KEY=VALUE` lines of `--config`: any flag, by its name with `_` or `-`.
fn read_config(path: &str) -> Result<Flags, String> {
    let text = fs::read(path).map_err(|e| format!("{path}: {e}"))?;
    let text = String::from_utf8_lossy(&text);
    let mut out = Flags::new();
    for (n, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let bad = || format!("{path}:{}: not a KEY=VALUE line of a known flag", n + 1);
        let (k, v) = line.split_once('=').ok_or_else(bad)?;
        let k = k.trim().replace('-', "_");
        if !known(&k) || k == "config" {
            return Err(bad());
        }
        let v = v.trim().to_owned();
        let v = if SWITCHES.contains(&k.as_str()) {
            if ["1", "true", "yes"].contains(&v.to_lowercase().as_str()) {
                "1".into()
            } else {
                String::new()
            }
        } else {
            v
        };
        if k == "allowed_host" || k == "trust_proxy" {
            out.entry(k).or_default().push(v);
        } else {
            out.insert(k, vec![v]);
        }
    }
    Ok(out)
}

/// The command line: (the flags of `--config`, then those given, which win)
/// or the error to print.
fn parse_args(args: &[String]) -> Result<Flags, String> {
    let mut given = Flags::new();
    let mut config = None;
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        i += 1;
        let Some(flag) = a.strip_prefix("--") else {
            return Err(format!("unrecognized arguments: {a}"));
        };
        let (name, inline) = match flag.split_once('=') {
            Some((n, v)) => (n, Some(v.to_owned())),
            None => (flag, None),
        };
        let k = name.replace('-', "_");
        if k != "config" && !known(&k) {
            return Err(format!("unrecognized arguments: {a}"));
        }
        if SWITCHES.contains(&k.as_str()) {
            if inline.is_some() {
                return Err(format!("argument --{name}: ignored explicit argument"));
            }
            given.insert(k, vec!["1".into()]);
            continue;
        }
        let v = match inline {
            Some(v) => v,
            None if i < args.len() && !args[i].starts_with("--") => {
                i += 1;
                args[i - 1].clone()
            }
            None => return Err(format!("argument --{name}: expected one argument")),
        };
        if k == "config" {
            config = Some(v);
        } else if k == "allowed_host" || k == "trust_proxy" {
            given.entry(k).or_default().push(v);
        } else {
            given.insert(k, vec![v]);
        }
    }
    let mut out = match config {
        Some(p) => read_config(&p)?,
        None => Flags::new(),
    };
    for (k, v) in given {
        if k == "allowed_host" || k == "trust_proxy" {
            out.entry(k).or_default().extend(v);
        } else {
            out.insert(k, v);
        }
    }
    Ok(out)
}

fn one<'a>(f: &'a Flags, k: &str) -> Option<&'a str> {
    f.get(k).and_then(|v| v.last()).map(String::as_str)
}

fn number<T: std::str::FromStr>(f: &Flags, k: &str, kind: &str) -> Result<Option<T>, String> {
    match one(f, k) {
        None => Ok(None),
        Some(v) => v.trim().parse().map(Some).map_err(|_| {
            format!(
                "argument --{}: invalid {kind} value: '{v}'",
                k.replace('_', "-")
            )
        }),
    }
}

/// The site's [`Config`] from the command line, or (exit code, message).
pub fn configure(args: &[String]) -> Result<Config, (u8, String)> {
    let usage = |m: String| (2u8, format!("{USAGE}og live serve: error: {m}"));
    let f = parse_args(args).map_err(usage)?;
    for k in NOT_YET {
        if f.get(k).is_some_and(|v| v.iter().any(|x| !x.is_empty())) {
            return Err((
                2,
                format!(
                    "og live serve: --{} is not ported yet: the password, TLS and proxy options \
                     land with the login (bin/owngenome_live.py has them until then)",
                    k.replace('_', "-")
                ),
            ));
        }
    }
    if f.contains_key("plugin_dir") {
        return Err((
            2,
            "og live serve: --plugin-dir: the plugins are compiled into og".into(),
        ));
    }
    if std::env::var("OWNGENOME_LIVE_PASSWORD_HASH").is_ok_and(|h| !h.trim().is_empty()) {
        return Err((
            2,
            "og live serve: OWNGENOME_LIVE_PASSWORD_HASH is set, and the password is not ported \
             yet (bin/owngenome_live.py has it until then)"
                .into(),
        ));
    }
    let Some(root) = one(&f, "root").filter(|r| !r.is_empty()) else {
        return Err(usage(
            "--root is required (on the command line or in --config)".into(),
        ));
    };
    let mut cfg = Config::new(root);
    cfg.results = Some(crate::results::FNS);
    if let Some(s) = one(&f, "state_dir").filter(|s| !s.is_empty()) {
        cfg.state_dir = PathBuf::from(s);
    }
    cfg.access_log = one(&f, "access_log")
        .filter(|s| !s.is_empty())
        .map_or_else(|| cfg.state_dir.join("access.log"), PathBuf::from);
    if let Some(b) = one(&f, "bind") {
        cfg.bind = b.to_owned();
    }
    if let Some(p) = number::<u16>(&f, "port", "int").map_err(usage)? {
        cfg.port = p;
    }
    if let Some(n) = number::<usize>(&f, "max_streams", "int").map_err(usage)? {
        cfg.max_streams = n;
    }
    if let Some(n) = number::<f64>(&f, "access_log_max_mb", "float").map_err(usage)? {
        cfg.access_log_max_mb = n;
    }
    for k in UNUSED {
        if k.ends_with("hours") || k.ends_with("minutes") {
            number::<f64>(&f, k, "float").map_err(usage)?;
        } else if k != "ingest_bind" {
            number::<i64>(&f, k, "int").map_err(usage)?;
        }
    }
    if let Some(a) = one(&f, "assets") {
        cfg.assets = PathBuf::from(a);
        if let Some(p) = cfg.assets.parent() {
            cfg.repo_assets = p.to_path_buf();
        }
    }
    if let Some(s) = one(&f, "steps") {
        cfg.steps = PathBuf::from(s);
    }
    cfg.quiet = one(&f, "quiet").is_some_and(|v| !v.is_empty());
    let path = |k: &str| one(&f, k).filter(|v| !v.is_empty()).map(PathBuf::from);
    cfg.viewer_fasta = path("viewer_fasta");
    cfg.viewer_igv_js = path("viewer_igv_js");
    cfg.metrics_token_file = path("metrics_token_file");
    cfg.ops = OpsConfig {
        root: Some(cfg.root.clone()),
        eta_box: one(&f, "eta_box").map(str::to_owned),
        eta_cores: number::<usize>(&f, "eta_cores", "int").map_err(usage)?,
        disk_warn_percent: one(&f, "disk_warn_percent").map(str::to_owned),
        mem_warn_percent: one(&f, "mem_warn_percent").map(str::to_owned),
        stall_minutes: number::<f64>(&f, "stall_minutes", "float").map_err(usage)?,
    };
    if !is_loopback(&cfg.bind) {
        return Err((
            2,
            format!(
                "ERROR: --bind {} is not loopback, and og live serve takes no password yet: \
                 it listens on loopback only until the login is ported",
                cfg.bind
            ),
        ));
    }
    if !cfg.root.is_dir() {
        return Err((
            1,
            format!("ERROR: --root {} is not a folder", cfg.root.display()),
        ));
    }
    Ok(cfg)
}

/// `og live serve ARGS`: serves until it is stopped. Exit 2 for flags it
/// cannot take, 1 when it cannot start.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if args.iter().any(|a| a == "-h" || a == "--help") {
        let _ = write!(out, "{USAGE}");
        return 0;
    }
    let cfg = match configure(args) {
        Ok(c) => c,
        Err((code, m)) => {
            let _ = writeln!(err, "{m}");
            return code;
        }
    };
    let bind = cfg.bind.clone();
    let srv = match Server::start(cfg) {
        Ok(s) => s,
        Err(e) => {
            let _ = writeln!(err, "og live serve: {e}");
            return 1;
        }
    };
    let _ = writeln!(
        err,
        "og live serve: no password hash set: anyone who can reach the loopback port sees the runs"
    );
    let _ = writeln!(err, "og live serve: browser http://{bind}:{}/", srv.port());
    let _ = err.flush();
    srv.wait();
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(a: &[&str]) -> Vec<String> {
        a.iter().map(|x| (*x).to_owned()).collect()
    }

    #[test]
    fn hosts_with_and_without_ports() {
        assert!(port_ok("") && port_ok(":8787") && !port_ok(":") && !port_ok(":x"));
        assert_eq!(netloc("http://127.0.0.1:8787/x"), "127.0.0.1:8787");
        assert_eq!(netloc("null"), "");
    }

    #[test]
    fn flags_and_config() {
        let dir = std::env::temp_dir().join(format!("og-live-cfg-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let c = dir.join("live.conf");
        fs::write(
            &c,
            format!(
                "# x\nroot={}\nport=9000\nquiet=yes\nviewer-fasta=/r.fa\n",
                dir.display()
            ),
        )
        .unwrap();
        let cfg = configure(&args(&["--config", c.to_str().unwrap(), "--port=9001"])).unwrap();
        assert_eq!(cfg.port, 9001);
        assert!(cfg.quiet);
        assert_eq!(cfg.viewer_fasta.as_deref(), Some(Path::new("/r.fa")));
        assert_eq!(cfg.state_dir, dir.join(".owngenome-live"));
        let refused = |a: &[&str]| configure(&args(a)).err().map(|e| e.0);
        let root = dir.to_str().unwrap();
        assert_eq!(refused(&["--root", root, "--password-file", "/x"]), Some(2));
        assert_eq!(refused(&["--root", root, "--bind", "0.0.0.0"]), Some(2));
        assert_eq!(refused(&["--root", root, "--tls-cert", "a"]), Some(2));
        assert_eq!(refused(&["--root", root, "--port", "x"]), Some(2));
        assert_eq!(refused(&["--root", root, "--nope"]), Some(2));
        assert_eq!(refused(&["--port", "1"]), Some(2));
        assert_eq!(refused(&["--root", root, "--ingest-port", "8788"]), None);
        fs::write(&c, "nope=1\n").unwrap();
        assert_eq!(refused(&["--config", c.to_str().unwrap()]), Some(2));
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn loopback() {
        assert!(is_loopback("127.0.0.1") && is_loopback("::1") && is_loopback("[::1]"));
        assert!(is_loopback("localhost") && !is_loopback("0.0.0.0") && !is_loopback("192.0.2.1"));
    }

    #[test]
    fn access_rows() {
        let dir = std::env::temp_dir().join(format!("og-live-log-{}", std::process::id()));
        let log = AccessLog::new(&dir.join("a.log"), 0.0).unwrap();
        log.write(&["127.0.0.1", "", "file", "S1", "a\tb\nc", "200"]);
        let text = fs::read_to_string(dir.join("a.log")).unwrap();
        let mut lines = text.lines();
        assert_eq!(lines.next(), Some(AccessLog::HEADER.trim_end()));
        let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
        assert_eq!(&row[1..], ["127.0.0.1", "-", "file", "S1", "a b c", "200"]);
        fs::remove_dir_all(&dir).unwrap();
    }
}
