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
//! With a password hash every route but `/login`, `/assets/` and
//! `/metrics` (which takes its own bearer token) needs a session
//! ([`crate::auth`]). The ingest listener ([`crate::ingest`]) stores the
//! run's messages under `STATE_DIR/runs/<run>/events.jsonl`, the file the
//! pages read. The access log (`STATE_DIR/access.log`, TSV) records who
//! opened what, never a query string, a request body, a result or the
//! password.
//!
//! Without a password hash the site serves loopback only, and it refuses
//! what would put it in front of other people: the proxy, Host and TLS
//! flags exit 2, and a request a proxy forwarded (`Forwarded`,
//! `X-Forwarded-*`, `X-Real-IP`) gets 403. The ingest port listens on
//! loopback only and takes a message only at the path that holds the
//! token the server keeps in `STATE_DIR/ingest-token`.

use crate::auth::{self, Auth, Net};
use crate::csp::{self, CSP, DENY, RUN_HTML_CSP};
use crate::family::{inside, realpath};
use crate::ops::OpsConfig;
use crate::plugins::{self, ResultsFns};
use crate::runs::{self, Run};
use crate::state::{self, Live, dumps, lock};
use crate::static_files as sf;
use crate::viewer::unquote;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::fs::{self, OpenOptions};
use std::io::{self, BufRead, BufReader, Read, Write};
use std::net::{IpAddr, Shutdown, SocketAddr, TcpListener, TcpStream, ToSocketAddrs};
use std::os::unix::fs::DirBuilderExt;
use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};
use std::time::{Duration, Instant};

/// How long a client may take to send its whole request, head and body
/// (after the TLS handshake); each write of the answer gets as long. An
/// event stream, which only writes once its request is read, lives on.
pub const REQUEST_SECONDS: u64 = 60;
/// The connections each listener serves at once; one more is closed at
/// once.
pub const MAX_CONNECTIONS: usize = 256;
/// How often an open event stream writes when nothing happens.
pub const HEARTBEAT: Duration = Duration::from_secs(15);
/// The largest request head read.
const MAX_HEAD: usize = 64 * 1024;
/// The cookie a login sets.
pub const COOKIE: &str = "owngenome_live";

const USAGE: &str =
    "usage: og live serve --root GENOME_DIR [--state-dir DIR] [--bind ADDR] [--port N]
                     [--ingest-bind ADDR] [--ingest-port N] [--password-file FILE]
                     [--allowed-host NAME ...] [--trust-proxy CIDR ...]
                     [--tls-cert FILE --tls-key FILE] [--cookie-secure]
                     [--session-hours H] [--session-idle-minutes M] [--login-max-failures N]
                     [--max-streams N] [--max-body-kb N] [--access-log FILE]
                     [--access-log-max-mb N] [--config FILE] [--quiet]
                     [--viewer-fasta FILE --viewer-igv-js FILE]
                     [--disk-warn-percent N] [--mem-warn-percent N] [--stall-minutes N]
                     [--metrics-token-file FILE] [--eta-box NAME] [--eta-cores N]
       og live serve hash-password
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
    /// Each person's results (`results.json`); `None` until `og live
    /// results` lands, and then the results, the findings and the map say
    /// they are not here yet.
    pub results: Option<ResultsFns>,
    /// How often a quiet event stream writes.
    pub heartbeat: Duration,
    /// The password and the sessions.
    pub auth: Arc<Auth>,
    /// Host names the browser may use besides loopback and the bind
    /// address.
    pub allowed_host: Vec<String>,
    /// The proxies whose `X-Forwarded-For` names the client.
    pub trust_proxy: Vec<Net>,
    /// The certificate chain and key of HTTPS.
    pub tls: Option<(PathBuf, PathBuf)>,
    /// Mark the cookie `Secure` (TLS ends at a proxy).
    pub cookie_secure: bool,
    /// Where the ingest listener listens; `None`: none is started.
    pub ingest: Option<(String, u16)>,
    /// The largest message the ingest listener takes.
    pub max_body_kb: usize,
    /// The token the ingest listener takes, made by [`Server::start`]
    /// (`STATE_DIR/ingest-token`); `None`: it takes no message.
    pub ingest_token: Option<String>,
    /// How long a client may take to send its request ([`REQUEST_SECONDS`]).
    pub request_timeout: Duration,
    /// The connections each listener serves at once ([`MAX_CONNECTIONS`]).
    pub max_connections: usize,
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
            auth: Arc::new(Auth::default()),
            allowed_host: Vec::new(),
            trust_proxy: Vec::new(),
            tls: None,
            cookie_secure: false,
            ingest: None,
            max_body_kb: 1024,
            ingest_token: None,
            request_timeout: Duration::from_secs(REQUEST_SECONDS),
            max_connections: MAX_CONNECTIONS,
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
    /// The peer, or the client a trusted proxy names.
    pub client: String,
    /// The session id (SHA-256 of the cookie) when it is live.
    pub sid: Option<String>,
    /// Its first 8 characters for the access log, else `-`.
    pub session8: String,
    input: &'a mut dyn BufRead,
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
        self.live.access.write(&[
            &self.client,
            &self.session8,
            action,
            run,
            section,
            &status.to_string(),
        ]);
    }

    /// The first `n` bytes of the request's body.
    pub fn read_body(&mut self, n: usize) -> io::Result<Vec<u8>> {
        let mut b = Vec::with_capacity(n);
        Read::take(&mut *self.input, n as u64).read_to_end(&mut b)?;
        Ok(b)
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

/// Whether a request says it came through a proxy: a `Forwarded`,
/// `X-Forwarded-*` or `X-Real-IP` header. Without a password the site is
/// for this machine only, and a proxy would hand it to anyone who reaches
/// the proxy.
fn forwarded(req: &Req<'_>) -> bool {
    req.headers
        .iter()
        .any(|(k, _)| k == "forwarded" || k == "x-real-ip" || k.starts_with("x-forwarded-"))
}

/// Answers one request.
pub fn handle(req: &mut Req<'_>) -> io::Result<()> {
    if req.live.cfg.auth.hash.is_none() && forwarded(req) {
        req.log("forwarded-refused", "-", "-", 403);
        return req.deny(
            403,
            Some("a request through a proxy needs a password (docs/live-site.md)"),
        );
    }
    if !host_ok(req) {
        req.log("host-refused", "-", "-", 403);
        return req.deny(403, None);
    }
    if req.method == "POST" && !origin_ok(req) {
        req.log("origin-refused", "-", "-", 403);
        return req.deny(403, None);
    }
    let authed = auth::begin(req);
    let path = req.path.clone();
    if path == "/login" {
        return auth::login(req);
    }
    if path == "/logout" && req.method == "POST" {
        return auth::logout(req);
    }
    if req.method == "GET"
        && let Some(name) = path.strip_prefix("/assets/")
    {
        return asset(req, name);
    }
    // /metrics takes a bearer token of its own (crate::plugins).
    if !authed && path != "/metrics" {
        if req.method != "GET" {
            return req.deny(404, None);
        }
        return auth::refuse(req);
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
    let token = auth::cookie_token(req.header("cookie"));
    let out = req.out();
    out.write_all(b"retry: 3000\n\n")?;
    out.flush()?;
    let mut sent = after;
    while !live.stopping.load(Ordering::SeqCst) {
        let new = live
            .store
            .wait(name, sent, live.cfg.heartbeat, &live.stopping);
        if live.cfg.auth.hash.is_some() && !live.cfg.auth.alive(token.as_deref()) {
            break; // logged out or expired: an open stream ends with its session
        }
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
pub(crate) fn read_head(r: &mut dyn BufRead) -> Option<(String, String, Headers)> {
    let mut head = Vec::new();
    loop {
        let mut line = Vec::new();
        let n = Read::take(&mut *r, (MAX_HEAD - head.len()) as u64 + 1)
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

/// A socket whose reads share one deadline, so a client that sends its
/// request a byte at a time cannot hold the connection past it. A write
/// keeps the socket's own write timeout.
pub(crate) struct ReadBy {
    pub(crate) tcp: TcpStream,
    pub(crate) end: Instant,
}

impl Read for ReadBy {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let left = self.end.saturating_duration_since(Instant::now());
        if left.is_zero() {
            return Err(io::Error::new(
                io::ErrorKind::TimedOut,
                "the request took too long",
            ));
        }
        self.tcp.set_read_timeout(Some(left))?;
        (&self.tcp).read(buf)
    }
}

impl Write for ReadBy {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        (&self.tcp).write(buf)
    }
    fn flush(&mut self) -> io::Result<()> {
        (&self.tcp).flush()
    }
}

/// One TLS connection, read and written by turns from the one thread that
/// serves it.
struct Shared(Rc<RefCell<rustls::StreamOwned<rustls::ServerConnection, ReadBy>>>);

impl Read for Shared {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.0.borrow_mut().read(buf)
    }
}

impl Write for Shared {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0.borrow_mut().write(buf)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.0.borrow_mut().flush()
    }
}

fn serve_conn(live: &Live, stream: TcpStream, tls: Option<&Arc<rustls::ServerConfig>>) {
    let peer = stream
        .peer_addr()
        .map(|a| a.ip().to_string())
        .unwrap_or_else(|_| "-".into());
    let timeout = live.cfg.request_timeout;
    let Some(cfg) = tls else {
        let _ = stream.set_write_timeout(Some(timeout));
        let Ok(read_half) = stream.try_clone() else {
            return;
        };
        let mut reader = BufReader::new(ReadBy {
            tcp: read_half,
            end: Instant::now() + timeout,
        });
        let mut out = io::BufWriter::new(&stream);
        serve_io(live, &peer, &mut reader, &mut out);
        drop(out);
        let _ = stream.shutdown(Shutdown::Write);
        // What the client still sends is read and dropped once the answer
        // is out, so closing does not reset the connection before the
        // client has read it.
        reader.get_mut().end = Instant::now() + Duration::from_secs(1);
        let _ = io::copy(&mut reader.take(MAX_HEAD as u64), &mut io::sink());
        return;
    };
    let Some(tls) = auth::accept_tls(cfg, stream) else {
        return;
    };
    // The request's deadline starts after the handshake, which has its own.
    let (conn, sock) = tls.into_parts();
    let _ = sock.set_write_timeout(Some(timeout));
    let sock = ReadBy {
        tcp: sock,
        end: Instant::now() + timeout,
    };
    let tls = Rc::new(RefCell::new(rustls::StreamOwned::new(conn, sock)));
    let mut reader = BufReader::new(Shared(Rc::clone(&tls)));
    let mut out = io::BufWriter::new(Shared(Rc::clone(&tls)));
    serve_io(live, &peer, &mut reader, &mut out);
    drop(out);
    let mut t = tls.borrow_mut();
    t.conn.send_close_notify();
    let _ = t.flush();
    let _ = t.sock.tcp.shutdown(Shutdown::Write);
}

/// One connection being served; dropping it frees its place.
struct Slot(Arc<AtomicUsize>);

impl Drop for Slot {
    fn drop(&mut self) {
        self.0.fetch_sub(1, Ordering::SeqCst);
    }
}

/// Serves each connection of `listener` with `serve`, on a thread of its
/// own, until `live` is stopping. At most `cfg.max_connections` are served
/// at once: one more is closed at once, so a flood of connections cannot
/// start a thread each.
pub(crate) fn accept_loop<F>(live: Arc<Live>, listener: TcpListener, serve: F) -> JoinHandle<()>
where
    F: Fn(&Live, TcpStream) + Clone + Send + 'static,
{
    thread::spawn(move || {
        let open = Arc::new(AtomicUsize::new(0));
        for conn in listener.incoming() {
            if live.stopping.load(Ordering::SeqCst) {
                break;
            }
            let Ok(stream) = conn else { continue };
            let before = open.fetch_add(1, Ordering::SeqCst);
            let slot = Slot(Arc::clone(&open));
            if before >= live.cfg.max_connections {
                continue; // drops the stream, which closes it, and the slot
            }
            let (l, serve) = (Arc::clone(&live), serve.clone());
            thread::spawn(move || {
                let _slot = slot;
                serve(&l, stream);
            });
        }
    })
}

/// Reads one request from `reader` and answers it on `out`.
fn serve_io(live: &Live, peer: &str, reader: &mut dyn BufRead, out: &mut dyn Write) {
    let Some((method, target, headers)) = read_head(reader) else {
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
    // Every X-Forwarded-For field in the order sent, as one list: a proxy
    // may add its entry as a field of its own.
    let forwarded_for: Vec<&str> = headers
        .iter()
        .filter(|(k, _)| k == "x-forwarded-for")
        .map(|(_, v)| v.as_str())
        .collect();
    let forwarded_for = (!forwarded_for.is_empty()).then(|| forwarded_for.join(","));
    let client = auth::client(peer, forwarded_for.as_deref(), &live.cfg.trust_proxy);
    let mut req = Req {
        live,
        method: m.into(),
        head_only,
        path,
        query: query.split('#').next().unwrap_or("").to_owned(),
        headers,
        client,
        sid: None,
        session8: "-".into(),
        input: reader,
        out,
    };
    let _ = handle(&mut req);
    let _ = req.out.flush();
}

/// A running site: [`Server::start`] binds and serves on threads, dropping
/// it stops it.
pub struct Server {
    live: Arc<Live>,
    addr: SocketAddr,
    ingest_addr: Option<SocketAddr>,
    accept: Option<JoinHandle<()>>,
    ingest: Option<JoinHandle<()>>,
}

/// Binds `bind:port` (a bracketed IPv6 address too).
fn listen(bind: &str, port: u16, what: &str) -> Result<TcpListener, String> {
    let host = bind.trim_matches(|c| c == '[' || c == ']');
    TcpListener::bind((host, port))
        .map_err(|e| format!("cannot listen{what} on {bind}:{port}: {e}"))
}

impl Server {
    /// Reads steps.tsv, opens the access log, binds `cfg.bind:cfg.port` (0:
    /// any free port) and the ingest listener of `cfg.ingest`, and serves.
    pub fn start(mut cfg: Config) -> Result<Server, String> {
        let steps = og_runlog::steps::read_steps(&cfg.steps).map_err(|e| format!("steps: {e}"))?;
        // Mode 755, whatever the umask: the ingest port needs a state folder
        // no other account can write to.
        fs::DirBuilder::new()
            .recursive(true)
            .mode(0o755)
            .create(&cfg.state_dir)
            .map_err(|e| format!("--state-dir {}: {e}", cfg.state_dir.display()))?;
        let access = AccessLog::new(&cfg.access_log, cfg.access_log_max_mb)
            .map_err(|e| format!("--access-log {}: {e}", cfg.access_log.display()))?;
        let tls = match &cfg.tls {
            Some((cert, key)) => Some(auth::tls_config(cert, key)?),
            None => None,
        };
        let listener = listen(&cfg.bind, cfg.port, "")?;
        let addr = listener.local_addr().map_err(|e| e.to_string())?;
        // The address that was bound, not the name: with no password the
        // site, and always the ingest port, stay on loopback.
        if cfg.auth.hash.is_none() && !addr.ip().is_loopback() {
            return Err(format!(
                "--bind {} listens on {addr}, not loopback, and no password hash is set",
                cfg.bind
            ));
        }
        let ingest_listener = match &cfg.ingest {
            Some((bind, port)) => {
                let l = listen(bind, *port, " for ingest")?;
                let a = l.local_addr().map_err(|e| e.to_string())?;
                if !a.ip().is_loopback() {
                    return Err(format!("--ingest-bind {bind} listens on {a}, not loopback"));
                }
                Some(l)
            }
            None => None,
        };
        if ingest_listener.is_some() {
            cfg.ingest_token = Some(crate::ingest::token(&cfg.state_dir)?);
        }
        let mut live = Live::new(cfg, steps, access);
        let more: Vec<String> = live
            .cfg
            .allowed_host
            .iter()
            .map(|h| h.to_lowercase())
            .collect();
        live.allowed.extend(more);
        let live = Arc::new(live);
        let (ingest_addr, ingest) = match ingest_listener {
            Some(l) => (
                Some(l.local_addr().map_err(|e| e.to_string())?),
                Some(crate::ingest::serve(Arc::clone(&live), l)),
            ),
            None => (None, None),
        };
        let accept = accept_loop(Arc::clone(&live), listener, move |l, s| {
            serve_conn(l, s, tls.as_ref())
        });
        Ok(Server {
            live,
            addr,
            ingest_addr,
            accept: Some(accept),
            ingest,
        })
    }

    /// The port it listens on.
    pub fn port(&self) -> u16 {
        self.addr.port()
    }

    /// The port of the ingest listener, when it has one.
    pub fn ingest_port(&self) -> Option<u16> {
        self.ingest_addr.map(|a| a.port())
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
        // Wakes each accept loop so it sees the flag.
        for (a, h) in [
            (Some(self.addr), self.accept.take()),
            (self.ingest_addr, self.ingest.take()),
        ] {
            if let (Some(a), Some(h)) = (a, h) {
                let _ = TcpStream::connect_timeout(&a, Duration::from_secs(1));
                let _ = h.join();
            }
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
/// `env_hash` is `OWNGENOME_LIVE_PASSWORD_HASH`, read only without
/// `--password-file`.
pub fn configure_with(args: &[String], env_hash: Option<&str>) -> Result<Config, (u8, String)> {
    let usage = |m: String| (2u8, format!("{USAGE}og live serve: error: {m}"));
    let f = parse_args(args).map_err(usage)?;
    if f.contains_key("plugin_dir") {
        return Err((
            2,
            "og live serve: --plugin-dir: the plugins are compiled into og".into(),
        ));
    }
    let Some(root) = one(&f, "root").filter(|r| !r.is_empty()) else {
        return Err(usage(
            "--root is required (on the command line or in --config)".into(),
        ));
    };
    let mut cfg = Config::new(root);
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
    let ingest_bind = one(&f, "ingest_bind").unwrap_or("127.0.0.1").to_owned();
    let ingest_port = number::<u16>(&f, "ingest_port", "int")
        .map_err(usage)?
        .unwrap_or(8788);
    cfg.ingest = Some((ingest_bind, ingest_port));
    if let Some(n) = number::<usize>(&f, "max_body_kb", "int").map_err(usage)? {
        // A Content-Length is read as a signed 64-bit number.
        if n.checked_mul(1024)
            .is_none_or(|b| i64::try_from(b).is_err())
        {
            return Err(usage(format!(
                "argument --max-body-kb: {n} is too large (at most {})",
                i64::MAX / 1024
            )));
        }
        cfg.max_body_kb = n;
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
    cfg.cookie_secure = one(&f, "cookie_secure").is_some_and(|v| !v.is_empty());
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
    let list = |k: &str| -> Vec<String> {
        f.get(k)
            .map(|v| v.iter().filter(|x| !x.is_empty()).cloned().collect())
            .unwrap_or_default()
    };
    cfg.allowed_host = list("allowed_host");
    cfg.trust_proxy = list("trust_proxy")
        .iter()
        .map(|c| Net::parse(c))
        .collect::<Result<_, _>>()
        .map_err(|e| (2, format!("ERROR: {e}")))?;
    let (cert, key) = (path("tls_cert"), path("tls_key"));
    // An empty --password-file (or password_file= in --config) is a
    // mistake, never a request for a site with no login.
    if one(&f, "password_file").is_some_and(|v| v.trim().is_empty()) {
        return Err(usage(
            "argument --password-file: empty: give the file 'og live serve hash-password' \
             wrote, or leave the flag out for a loopback site with no login"
                .into(),
        ));
    }
    let hash = auth::read_password_hash(path("password_file").as_deref(), env_hash)
        .map_err(|e| (1, format!("ERROR: password hash: {e}")))?;
    let mut a = Auth::new(hash);
    // NaN or inf would make a login last for ever.
    let limit = |k: &str| -> Result<Option<f64>, (u8, String)> {
        match number::<f64>(&f, k, "float").map_err(usage)? {
            Some(v) if !v.is_finite() => Err(usage(format!(
                "argument --{}: invalid float value: '{v}'",
                k.replace('_', "-")
            ))),
            v => Ok(v),
        }
    };
    if let Some(h) = limit("session_hours")? {
        a.session_hours = h;
    }
    if let Some(m) = limit("session_idle_minutes")? {
        a.session_idle_minutes = m;
    }
    if let Some(n) = number::<usize>(&f, "login_max_failures", "int").map_err(usage)? {
        a.login_max_failures = n;
    }
    let password = a.hash.is_some();
    cfg.auth = Arc::new(a);
    if !password && let Some(flag) = network_flag(&cfg, cert.is_some() || key.is_some()) {
        return Err((
            2,
            format!(
                "ERROR: {flag} is for a site reached through a proxy or over the network, and \
                 no password hash is set, so anyone who reaches it would see the runs: run \
                 'og live serve hash-password' and pass --password-file (or \
                 OWNGENOME_LIVE_PASSWORD_HASH), or leave {flag} out"
            ),
        ));
    }
    cfg.tls = match (cert, key) {
        (Some(c), Some(k)) => Some((c, k)),
        (None, None) => None,
        _ => return Err((1, "ERROR: --tls-cert and --tls-key go together".into())),
    };
    if let Some((bind, _)) = &cfg.ingest
        && !is_loopback(bind)
    {
        return Err((
            2,
            format!(
                "ERROR: --ingest-bind {bind} is not loopback: the ingest port takes the messages \
                 of runs on this machine only (in a container, give it the host's network)"
            ),
        ));
    }
    if !password && !is_loopback(&cfg.bind) {
        return Err((
            2,
            format!(
                "ERROR: --bind {} is not loopback and no password hash is set: run \
                 'og live serve hash-password' and pass --password-file (or \
                 OWNGENOME_LIVE_PASSWORD_HASH)",
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

/// The first flag of a site reached through a proxy or over the network
/// that `cfg` sets (`tls`: `--tls-cert` or `--tls-key` was given).
fn network_flag(cfg: &Config, tls: bool) -> Option<&'static str> {
    if !cfg.trust_proxy.is_empty() {
        Some("--trust-proxy")
    } else if !cfg.allowed_host.is_empty() {
        Some("--allowed-host")
    } else if tls {
        Some("--tls-cert and --tls-key")
    } else if cfg.cookie_secure {
        Some("--cookie-secure")
    } else {
        None
    }
}

/// [`configure_with`] the `OWNGENOME_LIVE_PASSWORD_HASH` of the process.
pub fn configure(args: &[String]) -> Result<Config, (u8, String)> {
    let env = std::env::var("OWNGENOME_LIVE_PASSWORD_HASH").ok();
    configure_with(args, env.as_deref())
}

/// The lines `og live serve` prints as it starts, before it listens: what
/// is open to whom.
pub fn start_notes(cfg: &Config) -> Vec<String> {
    let mut out = Vec::new();
    if cfg.auth.hash.is_none() {
        out.push(
            "og live serve: no password hash set: anyone who can reach the loopback port sees \
             the runs"
                .into(),
        );
    }
    if cfg.auth.hash.is_some() && cfg.tls.is_none() && !is_loopback(&cfg.bind) {
        out.push(format!(
            "og live serve: WARNING: --bind {} is not loopback and serves plain HTTP: the \
             password, the session cookie and every page cross the network unencrypted; add \
             --tls-cert and --tls-key, or bind 127.0.0.1 behind a proxy that does TLS",
            cfg.bind
        ));
    }
    out
}

/// `og live serve ARGS`: serves until it is stopped. Exit 2 for flags it
/// cannot take, 1 when it cannot start. `og live serve hash-password`
/// prints the hash of a password.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if args.first().is_some_and(|a| a == "hash-password") {
        use std::io::IsTerminal;
        let stdin = io::stdin();
        let tty = stdin.is_terminal();
        return auth::hash_password_cli(&args[1..], &mut stdin.lock(), tty, out, err);
    }
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
    serve_cli(cfg, err)
}

/// `bind` as the host of a URL: an IPv6 address in brackets.
fn url_host(bind: &str) -> String {
    if bind.contains(':') && !bind.starts_with('[') {
        format!("[{bind}]")
    } else {
        bind.to_owned()
    }
}

/// Prints the start notes, starts the site of `cfg` and serves until the
/// process ends; 1 when it cannot start.
pub fn serve_cli(cfg: Config, err: &mut dyn Write) -> u8 {
    for line in start_notes(&cfg) {
        let _ = writeln!(err, "{line}");
    }
    let bind = cfg.bind.clone();
    let scheme = if cfg.tls.is_some() { "https" } else { "http" };
    let srv = match Server::start(cfg) {
        Ok(s) => s,
        Err(e) => {
            let _ = writeln!(err, "og live serve: {e}");
            let _ = err.flush();
            return 1;
        }
    };
    let _ = write!(
        err,
        "og live serve: browser {scheme}://{}:{}/",
        url_host(&bind),
        srv.port()
    );
    if let (Some(p), Some((ib, _))) = (srv.ingest_port(), &srv.live().cfg.ingest) {
        let _ = write!(
            err,
            " ; ingest http://{}:{p}/<token>/ingest/<run> (the token is in {})",
            url_host(ib),
            srv.live().cfg.state_dir.join("ingest-token").display()
        );
    }
    let _ = writeln!(err);
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
        let refused = |a: &[&str]| configure_with(&args(a), None).err().map(|e| e.0);
        let root = dir.to_str().unwrap();
        assert_eq!(refused(&["--root", root, "--password-file", "/x"]), Some(1));
        assert_eq!(refused(&["--root", root, "--bind", "0.0.0.0"]), Some(2));
        assert_eq!(refused(&["--root", root, "--tls-cert", "a"]), Some(2));
        assert_eq!(refused(&["--root", root, "--port", "x"]), Some(2));
        assert_eq!(refused(&["--root", root, "--nope"]), Some(2));
        assert_eq!(refused(&["--port", "1"]), Some(2));
        assert_eq!(refused(&["--root", root, "--ingest-port", "8788"]), None);
        for v in ["nan", "inf", "-inf"] {
            assert_eq!(
                refused(&["--root", root, "--session-hours", v]),
                Some(2),
                "{v}"
            );
            assert_eq!(
                refused(&["--root", root, "--session-idle-minutes", v]),
                Some(2),
                "{v}"
            );
        }
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
