//! HTTP through `curl`, failing closed.
//!
//! Every answer that is not what the caller asked for is a
//! [`LookupError`]: an HTTP error, no answer, an empty body, a body that is
//! not JSON when JSON was asked for. 4xx answers other than 429 fail at
//! once; 5xx, 429, timeouts and connection errors are asked again after 5,
//! 10, 15 s. A request to one host waits for the one before it to that host
//! to end, so parallel lookups never hit one upstream twice at once.
//!
//! `curl` does the transfer (it is on every machine OwnGenome runs on, and
//! `scripts/lib/common.sh` downloads with it too), so `og` carries no TLS
//! stack. The URL, the headers and any token go to curl on its standard
//! input as a config file, never on its command line, where `ps` would show
//! them.

use crate::LookupError;
use pycompat::json::{Value, loads};
use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Duration;

pub const USER_AGENT: &str = "OwnGenome freshness check";
/// Seconds without a byte before a request is given up (the socket timeout
/// the Python tool had).
pub const TIMEOUT: u64 = 30;
pub const RETRIES: u32 = 3;

/// What a request wants back.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Want {
    /// A non-empty JSON document.
    Json,
    /// Non-empty text.
    Text,
    /// No body: the headers are the answer.
    Head,
}

/// One request.
#[derive(Debug, Clone)]
pub struct Request<'a> {
    pub url: &'a str,
    pub method: &'a str,
    pub headers: Vec<(String, String)>,
    pub data: Option<&'a [u8]>,
    pub want: Want,
    /// `None`: the client's own.
    pub retries: Option<u32>,
    /// Seconds without progress; `None`: [`TIMEOUT`].
    pub timeout: Option<u64>,
}

impl<'a> Request<'a> {
    pub fn get(url: &'a str, want: Want) -> Self {
        Request {
            url,
            method: if want == Want::Head { "HEAD" } else { "GET" },
            headers: Vec::new(),
            data: None,
            want,
            retries: None,
            timeout: None,
        }
    }

    pub fn header(mut self, name: &str, value: &str) -> Self {
        self.headers.push((name.to_owned(), value.to_owned()));
        self
    }
}

/// The answer's headers, from the last response of a redirect chain.
#[derive(Debug, Clone, Default)]
pub struct Headers(Vec<(String, String)>);

impl Headers {
    pub fn from_pairs(pairs: Vec<(String, String)>) -> Headers {
        Headers(pairs)
    }

    /// The first header of this name, ignoring case.
    pub fn get(&self, name: &str) -> Option<&str> {
        self.0
            .iter()
            .find(|(n, _)| n.eq_ignore_ascii_case(name))
            .map(|(_, v)| v.as_str())
    }

    /// Parses curl's `--dump-header` file: the header block of the last
    /// response.
    fn parse(raw: &str) -> Headers {
        let mut blocks: Vec<Vec<(String, String)>> = vec![];
        let mut cur: Option<Vec<(String, String)>> = None;
        for line in raw.lines() {
            let line = line.trim_end_matches('\r');
            if line.starts_with("HTTP/") {
                if let Some(b) = cur.take() {
                    blocks.push(b);
                }
                cur = Some(vec![]);
            } else if let (Some(b), Some((n, v))) = (cur.as_mut(), line.split_once(':')) {
                b.push((n.trim().to_owned(), v.trim().to_owned()));
            }
        }
        blocks.extend(cur);
        Headers(blocks.pop().unwrap_or_default())
    }
}

/// A body as the caller wanted it.
#[derive(Debug, Clone)]
pub enum Body {
    Json(Value),
    Text(String),
    None,
}

/// Where requests go and how hard they try. Shared between threads.
pub struct Client {
    pub retries: u32,
    /// Seconds per retry step: the n-th retry waits `n * backoff`.
    pub backoff: u64,
    /// Print each host on stderr the first time it is asked (recheck).
    pub announce: Option<&'static str>,
    pub github: GitHub,
    hosts: Mutex<Vec<String>>,
    lanes: Mutex<HashMap<String, Arc<Mutex<()>>>>,
}

impl Default for Client {
    fn default() -> Self {
        Client::new(RETRIES)
    }
}

/// The host and port of a URL (Python's `urlsplit(url).netloc`).
pub fn netloc(url: &str) -> &str {
    let rest = url.split_once("://").map_or(url, |(_, r)| r);
    let end = rest.find(['/', '?', '#']).unwrap_or(rest.len());
    let auth = &rest[..end];
    auth.rsplit_once('@').map_or(auth, |(_, h)| h)
}

/// The scheme and host of a URL, the part a registry token is bound to.
pub fn origin(url: &str) -> (&str, &str) {
    (url.split_once("://").map_or("", |(s, _)| s), netloc(url))
}

/// `next` resolved against `base` (an absolute URL, an absolute path, or a
/// path relative to `base`'s folder).
pub fn join(base: &str, next: &str) -> String {
    if next.contains("://") {
        return next.to_owned();
    }
    let (scheme, host) = origin(base);
    if next.starts_with("//") {
        return format!("{scheme}:{next}");
    }
    if next.starts_with('/') {
        return format!("{scheme}://{host}{next}");
    }
    let path_start = base.find("://").map_or(0, |i| i + 3) + host.len();
    let path = &base[path_start..];
    let path = &path[..path.find(['?', '#']).unwrap_or(path.len())];
    let dir = &path[..path.rfind('/').map_or(0, |i| i + 1)];
    format!("{scheme}://{host}{dir}{next}")
}

/// A value for curl's config file: in double quotes, with `\` and `"`
/// escaped and no line break (a line break would start a new option).
/// The error never names the value: it can be an `Authorization` header,
/// and lookup errors end up in the report and the public issue.
fn quoted(v: &str) -> Result<String, LookupError> {
    if v.contains(['\n', '\r']) {
        return Err(LookupError::Failed(
            "a request value holds a line break".into(),
        ));
    }
    Ok(format!(
        "\"{}\"",
        v.replace('\\', "\\\\").replace('"', "\\\"")
    ))
}

/// A scratch folder for one request's body and headers (and the
/// self-test's fixtures), removed on drop.
pub(crate) struct Scratch(PathBuf);

impl Scratch {
    /// The folder.
    pub(crate) fn path(&self) -> &std::path::Path {
        &self.0
    }

    /// A new folder only this user can read, under a name nobody can guess
    /// and created here, never one that was already there: on a shared
    /// temporary folder another user cannot plant it and read an answer.
    pub(crate) fn new() -> std::io::Result<Scratch> {
        use std::os::unix::fs::DirBuilderExt;
        static N: AtomicUsize = AtomicUsize::new(0);
        let mut last = None;
        for _ in 0..16 {
            let nanos = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map_or(0, |d| d.as_nanos());
            let salt = std::collections::hash_map::RandomState::new();
            let tag = std::hash::BuildHasher::hash_one(&salt, (nanos, N.load(Ordering::Relaxed)));
            let d = std::env::temp_dir().join(format!(
                "og-fresh-{}-{}-{tag:016x}",
                std::process::id(),
                N.fetch_add(1, Ordering::Relaxed)
            ));
            match std::fs::DirBuilder::new().mode(0o700).create(&d) {
                Ok(()) => return Ok(Scratch(d)),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => last = Some(e),
                Err(e) => return Err(e),
            }
        }
        Err(last.unwrap_or_else(|| std::io::Error::other("no scratch folder")))
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// One transfer: the status (0 when there was no answer) or curl's error.
enum Outcome {
    Status(u16, Vec<u8>, Headers),
    NoAnswer(String),
}

impl Client {
    pub fn new(retries: u32) -> Client {
        Client {
            retries,
            backoff: 5,
            announce: None,
            github: GitHub::from_env(),
            hosts: Mutex::new(Vec::new()),
            lanes: Mutex::new(HashMap::new()),
        }
    }

    /// The hosts asked so far, in the order they were first asked.
    pub fn hosts(&self) -> Vec<String> {
        self.hosts.lock().expect("hosts lock").clone()
    }

    fn lane(&self, host: &str) -> Arc<Mutex<()>> {
        let mut lanes = self.lanes.lock().expect("lanes lock");
        lanes.entry(host.to_owned()).or_default().clone()
    }

    fn seen(&self, host: &str) {
        let mut hosts = self.hosts.lock().expect("hosts lock");
        if !hosts.iter().any(|h| h == host) {
            hosts.push(host.to_owned());
            if let Some(who) = self.announce {
                eprintln!("{who}: asking {host} (release names only, nothing about the sample)");
            }
        }
    }

    fn transfer(&self, req: &Request<'_>) -> Result<Outcome, LookupError> {
        let io =
            |e: std::io::Error| LookupError::Failed(format!("{}: cannot run curl ({e})", req.url));
        let scratch = Scratch::new().map_err(io)?;
        let body = scratch.0.join("body");
        let heads = scratch.0.join("headers");
        let timeout = req.timeout.unwrap_or(TIMEOUT);
        let mut cfg = String::new();
        cfg += "silent\nshow-error\nlocation\nmax-redirs = 10\n";
        cfg += &format!("connect-timeout = {timeout}\nspeed-limit = 1\nspeed-time = {timeout}\n");
        // A host that trickles a byte a second would pass the speed limit
        // forever: every attempt also ends after ten timeouts in all.
        cfg += &format!("max-time = {}\n", timeout * 10);
        cfg += &format!("user-agent = {}\n", quoted(USER_AGENT)?);
        for (n, v) in &req.headers {
            cfg += &format!("header = {}\n", quoted(&format!("{n}: {v}"))?);
        }
        if req.want == Want::Head {
            cfg += "head\n";
        } else if req.method != "GET" {
            cfg += &format!("request = {}\n", quoted(req.method)?);
        }
        if let Some(data) = req.data {
            let p = scratch.0.join("data");
            std::fs::write(&p, data).map_err(io)?;
            cfg += &format!("data-binary = {}\n", quoted(&format!("@{}", p.display()))?);
        }
        cfg += &format!("output = {}\n", quoted(&body.to_string_lossy())?);
        cfg += &format!("dump-header = {}\n", quoted(&heads.to_string_lossy())?);
        cfg += "write-out = \"%{http_code}\"\n";
        cfg += &format!("url = {}\n", quoted(req.url)?);
        let mut child = Command::new("curl")
            .args(["--disable", "--config", "-"])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(io)?;
        child
            .stdin
            .take()
            .expect("piped stdin")
            .write_all(cfg.as_bytes())
            .map_err(io)?;
        let out = child.wait_with_output().map_err(io)?;
        if !out.status.success() {
            let msg = String::from_utf8_lossy(&out.stderr).trim().to_owned();
            let msg = msg.strip_prefix("curl: ").unwrap_or(&msg).to_owned();
            return Ok(Outcome::NoAnswer(if msg.is_empty() {
                format!("curl exited {}", out.status.code().unwrap_or(-1))
            } else {
                msg
            }));
        }
        let code: u16 = String::from_utf8_lossy(&out.stdout)
            .trim()
            .parse()
            .unwrap_or(0);
        let raw = std::fs::read(&body).unwrap_or_default();
        let h = Headers::parse(&std::fs::read_to_string(&heads).unwrap_or_default());
        Ok(Outcome::Status(code, raw, h))
    }
}

impl Fetch for Client {
    fn github_api(&self) -> &GitHub {
        &self.github
    }

    /// The body as wanted and the headers, or why not.
    fn request(&self, req: &Request<'_>) -> Result<(Body, Headers), LookupError> {
        let host = netloc(req.url).to_owned();
        self.seen(&host);
        let lane = self.lane(&host);
        let retries = req.retries.unwrap_or(self.retries);
        let mut last = String::new();
        for attempt in 0..=retries {
            if attempt > 0 {
                std::thread::sleep(Duration::from_secs(self.backoff * u64::from(attempt)));
            }
            let outcome = {
                let _turn = lane
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                self.transfer(req)?
            };
            let (code, raw, headers) = match outcome {
                Outcome::NoAnswer(why) => {
                    last = format!("no answer ({why})");
                    continue;
                }
                Outcome::Status(c, r, h) => (c, r, h),
            };
            if !(200..300).contains(&code) {
                last = format!("HTTP {code}");
                if (400..500).contains(&code) && code != 429 {
                    break;
                }
                continue;
            }
            if req.want == Want::Head {
                return Ok((Body::None, headers));
            }
            let text = String::from_utf8_lossy(&raw).into_owned();
            if text.trim().is_empty() {
                last = "empty answer".into();
                continue;
            }
            if req.want == Want::Text {
                return Ok((Body::Text(text), headers));
            }
            match loads(&text) {
                Ok(v) => return Ok((Body::Json(v), headers)),
                Err(_) => {
                    last = "answer is not JSON".into();
                }
            }
        }
        Err(LookupError::Failed(format!("{}: {last}", req.url)))
    }
}

/// Where lookups get their answers: [`Client`] on the network, a stub in
/// the tests.
pub trait Fetch: Sync {
    /// The body as wanted and the headers, or why not.
    fn request(&self, req: &Request<'_>) -> Result<(Body, Headers), LookupError>;

    fn json(&self, url: &str) -> Result<Value, LookupError> {
        self.json_with(Request::get(url, Want::Json))
    }

    fn json_with(&self, req: Request<'_>) -> Result<Value, LookupError> {
        match self.request(&req)?.0 {
            Body::Json(v) => Ok(v),
            _ => Err(LookupError::Shape(format!("{}: no JSON body", req.url))),
        }
    }

    fn text(&self, url: &str) -> Result<String, LookupError> {
        match self.request(&Request::get(url, Want::Text))?.0 {
            Body::Text(t) => Ok(t),
            _ => Err(LookupError::Shape(format!("{url}: no text body"))),
        }
    }

    fn head(&self, url: &str) -> Result<Headers, LookupError> {
        self.head_with(Request::get(url, Want::Head))
    }

    fn head_with(&self, req: Request<'_>) -> Result<Headers, LookupError> {
        Ok(self.request(&req)?.1)
    }

    /// Where this fetcher's GitHub API is.
    fn github_api(&self) -> &GitHub;

    /// The GitHub API at `path`.
    fn github(&self, path: &str) -> Result<Value, LookupError> {
        let gh = self.github_api();
        let url = gh.url(path);
        self.json_with(gh.request(&url, "GET", None))
    }
}

/// Where the GitHub API is, and the token it gets.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GitHub {
    pub api: String,
    pub token: Option<String>,
}

impl GitHub {
    /// `FRESHNESS_GITHUB_API` (a test server), else api.github.com, and
    /// `GITHUB_TOKEN` when it is set.
    pub fn from_env() -> GitHub {
        GitHub {
            api: std::env::var("FRESHNESS_GITHUB_API")
                .unwrap_or_else(|_| "https://api.github.com".into()),
            token: std::env::var("GITHUB_TOKEN")
                .ok()
                .map(|t| t.trim().to_owned())
                .filter(|t| !t.is_empty()),
        }
    }

    /// A request to `path` of the API, with the token when there is one.
    pub fn request<'a>(
        &self,
        url: &'a str,
        method: &'a str,
        data: Option<&'a [u8]>,
    ) -> Request<'a> {
        let mut r = Request::get(url, Want::Json).header("Accept", "application/vnd.github+json");
        r.method = method;
        r.data = data;
        if let Some(t) = &self.token {
            r = r.header("Authorization", &format!("Bearer {t}"));
        }
        if data.is_some() {
            r = r.header("Content-Type", "application/json");
        }
        r
    }

    pub fn url(&self, path: &str) -> String {
        format!("{}{path}", self.api)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn urls() {
        assert_eq!(netloc("https://quay.io/v2/x"), "quay.io");
        assert_eq!(netloc("http://127.0.0.1:8080?a"), "127.0.0.1:8080");
        assert_eq!(
            join("https://a.io/v2/x/tags/list?n=1", "/v2/x/tags/list?last=b"),
            "https://a.io/v2/x/tags/list?last=b"
        );
        assert_eq!(
            join("https://a.io/v2/x/list", "more"),
            "https://a.io/v2/x/more"
        );
        assert_eq!(join("https://a.io/v2", "https://b.io/x"), "https://b.io/x");
        assert_eq!(origin("https://a.io:1/x"), ("https", "a.io:1"));
    }

    #[test]
    fn last_header_block_wins() {
        let h = Headers::parse(
            "HTTP/1.1 302 Found\r\nLocation: /x\r\nLast-Modified: old\r\n\r\nHTTP/2 200\r\nlast-modified: new\r\n\r\n",
        );
        assert_eq!(h.get("Last-Modified"), Some("new"));
        assert_eq!(h.get("Location"), None);
    }

    #[test]
    fn scratch_folders_are_new_and_private() {
        use std::os::unix::fs::PermissionsExt;
        let (a, b) = (Scratch::new().unwrap(), Scratch::new().unwrap());
        assert_ne!(a.0, b.0);
        let mode = std::fs::metadata(&a.0).unwrap().permissions().mode();
        assert_eq!(mode & 0o777, 0o700);
        let path = a.0.clone();
        drop(a);
        assert!(!path.exists());
    }

    #[test]
    fn config_values_are_quoted() {
        assert_eq!(quoted(r#"a"b\c"#).unwrap(), r#""a\"b\\c""#);
        assert!(quoted("a\nurl = x").is_err());
        let secret = "Authorization: Bearer s3cr3t-value\n";
        let Err(LookupError::Failed(msg)) = quoted(secret) else {
            panic!("a line break must fail");
        };
        assert!(!msg.contains("s3cr3t"), "the error names the value: {msg}");
    }
}
