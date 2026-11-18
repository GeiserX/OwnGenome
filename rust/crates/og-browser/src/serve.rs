//! A directory served on 127.0.0.1 for a page under test: GET only, one
//! request per connection, the headers the test asks for on every answer.
//! Only files inside the directory are served, after links are followed.

use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::net::{Shutdown, TcpListener, TcpStream};
use std::path::{Component, Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};
use std::time::Duration;

/// The open connections of a server: a clone of each socket, to close it,
/// and the thread that answers it.
type Open = Arc<Mutex<Vec<(TcpStream, JoinHandle<()>)>>>;

/// A static file server on a free port of 127.0.0.1. It stops when dropped.
pub struct StaticServer {
    port: u16,
    stop: Arc<AtomicBool>,
    accept: Option<JoinHandle<()>>,
    /// Each open connection (a clone of its socket) and its thread.
    conns: Open,
}

impl StaticServer {
    /// Serve the files under `dir`, with `headers` (name, value) added to
    /// every answer. A path ending in `/` serves its `index.html`.
    pub fn start(dir: &Path, headers: &[(&str, &str)]) -> std::io::Result<StaticServer> {
        let root = dir.canonicalize()?;
        let listener = TcpListener::bind("127.0.0.1:0")?;
        let port = listener.local_addr()?.port();
        let stop = Arc::new(AtomicBool::new(false));
        let conns: Open = Arc::new(Mutex::new(Vec::new()));
        let (root, extra) = (
            Arc::new(root),
            Arc::new(
                headers
                    .iter()
                    .map(|(k, v)| format!("{k}: {v}\r\n"))
                    .collect::<String>(),
            ),
        );
        let (flag, open) = (Arc::clone(&stop), Arc::clone(&conns));
        let accept = thread::spawn(move || {
            for conn in listener.incoming() {
                if flag.load(Ordering::SeqCst) {
                    break;
                }
                // A connection the server could not keep a handle to (to
                // close it on drop) is closed at once, never answered.
                if let Ok(c) = conn
                    && let Ok(peer) = c.try_clone()
                {
                    let (root, extra) = (Arc::clone(&root), Arc::clone(&extra));
                    let h = thread::spawn(move || answer(c, &root, &extra));
                    let mut open = open.lock().unwrap_or_else(|e| e.into_inner());
                    open.retain(|(_, h)| !h.is_finished());
                    open.push((peer, h));
                }
            }
        });
        Ok(StaticServer {
            port,
            stop,
            accept: Some(accept),
            conns,
        })
    }

    pub fn port(&self) -> u16 {
        self.port
    }

    /// `http://127.0.0.1:<port>/<path>`.
    pub fn url(&self, path: &str) -> String {
        format!(
            "http://127.0.0.1:{}/{}",
            self.port,
            path.trim_start_matches('/')
        )
    }
}

impl Drop for StaticServer {
    /// Stop accepting, close every open connection and wait for its thread,
    /// so nothing is served once the server is gone.
    fn drop(&mut self) {
        self.stop.store(true, Ordering::SeqCst);
        // Wake the accept loop so it sees the flag.
        let _ = TcpStream::connect(("127.0.0.1", self.port));
        if let Some(h) = self.accept.take() {
            let _ = h.join();
        }
        let open = std::mem::take(&mut *self.conns.lock().unwrap_or_else(|e| e.into_inner()));
        for (c, h) in open {
            let _ = c.shutdown(Shutdown::Both);
            let _ = h.join();
        }
    }
}

fn answer(mut c: TcpStream, root: &Path, extra: &str) {
    let _ = c.set_read_timeout(Some(Duration::from_secs(10)));
    let mut r = BufReader::new(&c);
    let mut first = String::new();
    if r.read_line(&mut first).is_err() {
        return;
    }
    // Read the rest of the head; a GET has no body.
    let mut line = String::new();
    while r.read_line(&mut line).is_ok_and(|n| n > 2) {
        line.clear();
    }
    let mut parts = first.split_whitespace();
    let (method, target) = (parts.next().unwrap_or(""), parts.next().unwrap_or(""));
    let (status, ctype, body) = if method != "GET" {
        (
            "405 Method Not Allowed",
            "text/plain",
            b"method not allowed\n".to_vec(),
        )
    } else {
        match file_of(root, target)
            .filter(|p| inside(root, p))
            .and_then(|p| fs::read(&p).ok().map(|b| (p, b)))
        {
            Some((p, b)) => ("200 OK", content_type(&p), b),
            None => ("404 Not Found", "text/plain", b"not found\n".to_vec()),
        }
    };
    let head = format!(
        "HTTP/1.1 {status}\r\nContent-Type: {ctype}\r\nContent-Length: {}\r\n{extra}Connection: close\r\n\r\n",
        body.len()
    );
    let _ = c
        .write_all(head.as_bytes())
        .and_then(|()| c.write_all(&body));
    let _ = c.shutdown(Shutdown::Write);
}

/// The file under `root` that the request target names, or `None` for a
/// path that leaves `root`.
fn file_of(root: &Path, target: &str) -> Option<PathBuf> {
    let path = target.split(['?', '#']).next()?;
    let mut rel = percent_decode(path.strip_prefix('/')?)?;
    if rel.is_empty() || rel.ends_with('/') {
        rel.push_str("index.html");
    }
    let rel = Path::new(&rel);
    if !rel.components().all(|c| matches!(c, Component::Normal(_))) {
        return None;
    }
    Some(root.join(rel))
}

/// Whether `p`, its links followed, is still under `root` (already
/// canonical). A missing file is not.
fn inside(root: &Path, p: &Path) -> bool {
    p.canonicalize().is_ok_and(|f| f.starts_with(root))
}

fn percent_decode(s: &str) -> Option<String> {
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'%' {
            let hex = std::str::from_utf8(b.get(i + 1..i + 3)?).ok()?;
            out.push(u8::from_str_radix(hex, 16).ok()?);
            i += 3;
        } else {
            out.push(b[i]);
            i += 1;
        }
    }
    String::from_utf8(out).ok()
}

fn content_type(p: &Path) -> &'static str {
    match p.extension().and_then(|e| e.to_str()) {
        Some("html") => "text/html; charset=utf-8",
        Some("js") => "text/javascript; charset=utf-8",
        Some("css") => "text/css; charset=utf-8",
        Some("json") => "application/json",
        Some("svg") => "image/svg+xml",
        Some("png") => "image/png",
        _ => "application/octet-stream",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_target_maps_into_the_root_only() {
        let root = Path::new("/srv");
        assert_eq!(file_of(root, "/"), Some(PathBuf::from("/srv/index.html")));
        assert_eq!(
            file_of(root, "/a/b.js?x=1"),
            Some(PathBuf::from("/srv/a/b.js"))
        );
        assert_eq!(
            file_of(root, "/my%20page.html"),
            Some(PathBuf::from("/srv/my page.html"))
        );
        assert_eq!(
            file_of(root, "/running/"),
            Some(PathBuf::from("/srv/running/index.html"))
        );
        assert_eq!(file_of(root, "/../etc/passwd"), None);
        assert_eq!(file_of(root, "/a/%2e%2e/%2e%2e/x"), None);
        assert_eq!(file_of(root, "/%zz"), None);
        assert_eq!(file_of(root, "no-slash"), None);
    }

    #[test]
    fn it_serves_a_file_with_the_headers_and_404_otherwise() {
        let d = std::env::temp_dir().join(format!("og-browser-serve-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        fs::write(d.join("a.json"), "{}").unwrap();
        let srv =
            StaticServer::start(&d, &[("Content-Security-Policy", "default-src 'none'")]).unwrap();
        let get = |path: &str| {
            let mut c = TcpStream::connect(("127.0.0.1", srv.port())).unwrap();
            write!(c, "GET {path} HTTP/1.1\r\nHost: x\r\n\r\n").unwrap();
            let mut s = String::new();
            std::io::Read::read_to_string(&mut c, &mut s).unwrap();
            s
        };
        let ok = get("/a.json");
        assert!(ok.starts_with("HTTP/1.1 200 OK\r\n"), "{ok}");
        assert!(
            ok.contains("\r\nContent-Type: application/json\r\n"),
            "{ok}"
        );
        assert!(
            ok.contains("\r\nContent-Security-Policy: default-src 'none'\r\n"),
            "{ok}"
        );
        assert!(ok.ends_with("\r\n\r\n{}"), "{ok}");
        let missing = get("/b.json");
        assert!(
            missing.starts_with("HTTP/1.1 404 Not Found\r\n"),
            "{missing}"
        );
        assert!(missing.contains("Content-Security-Policy"), "{missing}");
        assert_eq!(
            srv.url("/x/y.html"),
            format!("http://127.0.0.1:{}/x/y.html", srv.port())
        );
        drop(srv);
        fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn a_link_out_of_the_root_is_not_served() {
        let d = crate::Scratch::new(&std::env::temp_dir(), "serve-link");
        let (root, outside) = (d.path().join("site"), d.path().join("secret.txt"));
        fs::create_dir(&root).unwrap();
        fs::write(&outside, "secret").unwrap();
        fs::write(root.join("in.txt"), "fine").unwrap();
        std::os::unix::fs::symlink(&outside, root.join("out.txt")).unwrap();
        std::os::unix::fs::symlink(root.join("in.txt"), root.join("alias.txt")).unwrap();
        let root = root.canonicalize().unwrap();
        assert!(!inside(&root, &root.join("out.txt")));
        assert!(inside(&root, &root.join("alias.txt")));
        assert!(!inside(&root, &root.join("missing.txt")));
        let srv = StaticServer::start(&root, &[]).unwrap();
        let get = |path: &str| {
            let mut c = TcpStream::connect(("127.0.0.1", srv.port())).unwrap();
            write!(c, "GET {path} HTTP/1.1\r\nHost: x\r\n\r\n").unwrap();
            let mut s = String::new();
            std::io::Read::read_to_string(&mut c, &mut s).unwrap();
            s
        };
        assert!(get("/out.txt").starts_with("HTTP/1.1 404 "));
        assert!(get("/alias.txt").ends_with("\r\n\r\nfine"));
    }

    #[test]
    fn dropping_the_server_closes_an_open_connection() {
        let d = crate::Scratch::new(&std::env::temp_dir(), "serve-drop");
        let srv = StaticServer::start(d.path(), &[]).unwrap();
        // A client that sends half a request and waits: its thread blocks on read.
        let mut c = TcpStream::connect(("127.0.0.1", srv.port())).unwrap();
        write!(c, "GET /x HTTP/1.1\r\n").unwrap();
        // Wait until the server holds the connection and its thread runs.
        let held = |s: &StaticServer| {
            s.conns
                .lock()
                .unwrap()
                .iter()
                .any(|(_, h)| !h.is_finished())
        };
        let t = std::time::Instant::now();
        while !held(&srv) {
            assert!(t.elapsed() < Duration::from_secs(5), "no connection held");
            std::thread::sleep(Duration::from_millis(10));
        }
        let t = std::time::Instant::now();
        drop(srv);
        // Without the shutdown the drop would wait out the 10 s read timeout.
        assert!(t.elapsed() < Duration::from_secs(5), "{:?}", t.elapsed());
    }
}
