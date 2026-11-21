//! What the `port_live_server*` tests share: the synthetic DEMO-001 run of
//! `tests/fixtures/live` in a scratch genome folder, `og live serve` on it
//! in this process, and plain HTTP requests to it. Each test file includes
//! it with `#[path]`; on its own it holds no test.
#![allow(dead_code)]

use og_live::server::{Config, Server};
use std::fs;
use std::io::{Read, Write};
use std::net::TcpStream;
use std::path::{Path, PathBuf};
use std::time::Duration;

pub const RUN: &str = "DEMO-001";

/// The repository's root: four levels above this crate's tests.
pub fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(3)
        .expect("rust/crates/og-live sits three levels below the root")
        .to_path_buf()
}

/// A test's own folder, removed when it is dropped (also on a panic).
pub struct Work {
    pub dir: PathBuf,
}

impl Work {
    pub fn new(name: &str) -> Work {
        let d = Path::new(env!("CARGO_TARGET_TMPDIR"))
            .join(format!("og-live-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        Work {
            dir: fs::canonicalize(&d).unwrap(),
        }
    }

    pub fn root(&self) -> PathBuf {
        self.dir.join("genome")
    }

    pub fn state(&self) -> PathBuf {
        self.dir.join("state")
    }
}

impl Drop for Work {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.dir);
    }
}

fn copy_dir(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for e in fs::read_dir(from).unwrap() {
        let e = e.unwrap();
        let p = e.path();
        if p.is_dir() {
            copy_dir(&p, &to.join(e.file_name()));
        } else {
            fs::copy(&p, to.join(e.file_name())).unwrap();
        }
    }
}

pub fn write(p: &Path, text: &str) {
    fs::create_dir_all(p.parent().unwrap()).unwrap();
    fs::write(p, text).unwrap();
}

/// The genome folder of the checks: DEMO-001 with a report that links a
/// table, a symlink out of the root, CPSR's folder and a folder with two
/// pages, and files it does not link.
pub fn make_root(w: &Work) -> PathBuf {
    let root = w.root();
    let d = root.join(RUN);
    copy_dir(&repo().join("tests/fixtures/live").join(RUN), &d);
    write(
        &d.join("clinvar/hits.tsv"),
        "gene\tvariant\nGENEA\tchr1:100 A>G\n",
    );
    write(&d.join("clinvar/unlinked.tsv"), "not linked\n");
    write(&w.dir.join("outside.txt"), "outside the root\n");
    std::os::unix::fs::symlink(w.dir.join("outside.txt"), d.join("clinvar/escape.tsv")).unwrap();
    write(
        &d.join(format!("cpsr/{RUN}.cpsr.grch38.html")),
        "<!doctype html><title>CPSR</title>\n",
    );
    write(
        &d.join(format!("cpsr/{RUN}.cpsr.grch38.classification.tsv.gz")),
        "not linked\n",
    );
    write(&d.join("twopages/a.html"), "<!doctype html>\n");
    write(&d.join("twopages/b.html"), "<!doctype html>\n");
    write(
        &d.join(format!("{RUN}_report.html")),
        "<!doctype html><html><head><style>body{color:#eee}</style></head><body>\
         <a href=\"clinvar/hits.tsv\">hits</a> <a href=\"clinvar/escape.tsv\">escape</a> \
         <a href=\"cpsr/\">the cpsr folder</a> <a href=\"twopages/\">two pages</a> \
         <a href=\"https://example.org/x\">outside</a></body></html>\n",
    );
    root
}

/// The site's config for `w`: the repository's assets and steps.tsv, its
/// own state folder, port 0, no request lines.
pub fn config(w: &Work) -> Config {
    let r = repo();
    let mut c = Config::new(w.root());
    c.state_dir = w.state();
    c.access_log = w.state().join("access.log");
    c.port = 0;
    c.assets = r.join("assets/live");
    c.repo_assets = r.join("assets");
    c.steps = r.join("assets/steps.tsv");
    c.quiet = true;
    c
}

/// The events file of `run` under the state folder.
pub fn events_path(w: &Work, run: &str) -> PathBuf {
    w.state().join("runs").join(run).join("events.jsonl")
}

/// Appends messages to `run`'s events file, as the ingest listener stores
/// them.
pub fn store(w: &Work, run: &str, lines: &[&str]) {
    let p = events_path(w, run);
    fs::create_dir_all(p.parent().unwrap()).unwrap();
    let mut f = fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(p)
        .unwrap();
    for l in lines {
        f.write_all(l.trim().as_bytes()).unwrap();
        f.write_all(b"\n").unwrap();
    }
}

/// The fixture's messages (tests/fixtures/live/events.jsonl).
pub fn fixture_events() -> Vec<String> {
    fs::read_to_string(repo().join("tests/fixtures/live/events.jsonl"))
        .unwrap()
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(str::to_owned)
        .collect()
}

pub fn start(c: Config) -> Server {
    Server::start(c).unwrap()
}

/// One answer: status, headers (names lowered) and body.
#[derive(Debug, Clone)]
pub struct Resp {
    pub status: u16,
    pub headers: Vec<(String, String)>,
    pub body: Vec<u8>,
}

impl Resp {
    pub fn header(&self, name: &str) -> Option<&str> {
        self.headers
            .iter()
            .find(|(k, _)| k == name)
            .map(|(_, v)| v.as_str())
    }

    pub fn text(&self) -> String {
        String::from_utf8_lossy(&self.body).into_owned()
    }

    pub fn json(&self) -> pycompat::json::Value {
        pycompat::json::loads(&self.text())
            .unwrap_or_else(|e| panic!("not JSON ({e}): {}", self.text()))
    }
}

/// Sends `head` (with its blank line) and `body` as they are; reads the
/// whole answer.
pub fn raw(port: u16, head: &str, body: &[u8]) -> Resp {
    let mut s = TcpStream::connect(("127.0.0.1", port)).unwrap();
    s.set_read_timeout(Some(Duration::from_secs(30))).unwrap();
    s.write_all(head.as_bytes()).unwrap();
    s.write_all(body).unwrap();
    let mut all = Vec::new();
    let _ = s.read_to_end(&mut all);
    parse(&all)
}

pub fn parse(all: &[u8]) -> Resp {
    let split = all
        .windows(4)
        .position(|w| w == b"\r\n\r\n")
        .unwrap_or_else(|| panic!("no head in {:?}", String::from_utf8_lossy(all)));
    let head = String::from_utf8_lossy(&all[..split]).into_owned();
    let mut lines = head.split("\r\n");
    let status = lines
        .next()
        .and_then(|l| l.split(' ').nth(1))
        .and_then(|c| c.parse().ok())
        .unwrap_or(0);
    let headers = lines
        .filter_map(|l| l.split_once(':'))
        .map(|(k, v)| (k.trim().to_ascii_lowercase(), v.trim().to_owned()))
        .collect();
    Resp {
        status,
        headers,
        body: all[split + 4..].to_vec(),
    }
}

/// A request with `Host: 127.0.0.1:<port>` unless `headers` names one.
pub fn request(port: u16, method: &str, path: &str, headers: &[(&str, &str)]) -> Resp {
    let mut head = format!("{method} {path} HTTP/1.1\r\n");
    if !headers.iter().any(|(k, _)| k.eq_ignore_ascii_case("host")) {
        head.push_str(&format!("Host: 127.0.0.1:{port}\r\n"));
    }
    for (k, v) in headers {
        head.push_str(&format!("{k}: {v}\r\n"));
    }
    head.push_str("Connection: close\r\n\r\n");
    raw(port, &head, b"")
}

pub fn get(port: u16, path: &str) -> Resp {
    request(port, "GET", path, &[])
}

/// Whether an answer carries the site's policy: nothing from outside, no
/// framing, no sniffing, no cache.
pub fn policy_ok(r: &Resp) -> bool {
    let csp = r.header("content-security-policy").unwrap_or("");
    csp.contains("default-src 'none'")
        && csp.contains("frame-ancestors 'none'")
        && r.header("x-frame-options") == Some("DENY")
        && r.header("x-content-type-options") == Some("nosniff")
        && r.header("cache-control") == Some("no-store")
        && r.header("referrer-policy") == Some("no-referrer")
}

/// The access log's rows, without the header.
pub fn log_rows(w: &Work) -> Vec<Vec<String>> {
    fs::read_to_string(w.state().join("access.log"))
        .unwrap_or_default()
        .lines()
        .skip(1)
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect()
}
