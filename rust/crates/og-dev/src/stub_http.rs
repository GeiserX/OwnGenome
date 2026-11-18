//! `og dev stub-http RULE...`: a local HTTP server with scripted answers, for
//! the shell tests of downloads (`tests/test_fetch_retry.sh`,
//! `tests/test_e2e_data_cache.sh`, `tests/test_og_check.sh`).
//!
//! It listens on 127.0.0.1 on a free port, prints that port and a newline on
//! stdout, and answers until it is killed. A test reads the port from the
//! file it sent stdout to.
//!
//! Each RULE is `PATH=ANSWER`. A request is answered by the rule whose PATH
//! is the request's path (the query dropped), else by the rule `*=ANSWER`,
//! else with 404. The answers:
//!
//! - `404`: not found.
//! - `dir:DIR`: the file at DIR plus the request's path, or 404 when there is
//!   none or the path holds a `..` component.
//! - `flat:DIR`: the file in DIR named as the last component of the path, or
//!   404.
//! - `drip:LEN:CHUNK:MS`: 200 with `Content-Length: LEN`, then CHUNK bytes of
//!   `x` every MS milliseconds until LEN bytes are sent or the client goes.
//!   `drip:N:N:0` is a whole file of N bytes; `drip:100000:1000:60000` sends
//!   1000 bytes and then holds the connection, a stalled transfer.
//!
//! Every answer closes the connection after it (HTTP/1.0), and a HEAD
//! request gets the headers without the body.

use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::path::{Component, Path, PathBuf};
use std::time::Duration;

const USAGE: &str = "usage: og dev stub-http PATH=ANSWER... (ANSWER: 404, dir:DIR, flat:DIR or drip:LEN:CHUNK:MS; PATH * for any path)";

/// What one rule answers.
#[derive(Clone, Debug, PartialEq)]
enum Answer {
    NotFound,
    Dir(PathBuf),
    Flat(PathBuf),
    Drip { len: u64, chunk: u64, ms: u64 },
}

/// The rules, in the order given: (path, answer).
type Rules = Vec<(String, Answer)>;

fn number(s: &str, what: &str, rule: &str) -> Result<u64, String> {
    s.parse()
        .map_err(|_| format!("{what} '{s}' of rule '{rule}' is not a whole number"))
}

fn parse(args: &[String]) -> Result<Rules, String> {
    if args.is_empty() {
        return Err("no rule given".into());
    }
    let mut rules = Vec::new();
    for rule in args {
        let (path, answer) = rule
            .split_once('=')
            .ok_or_else(|| format!("rule '{rule}' is not PATH=ANSWER"))?;
        if path != "*" && !path.starts_with('/') {
            return Err(format!("the path of rule '{rule}' is neither * nor /..."));
        }
        let answer = if answer == "404" {
            Answer::NotFound
        } else if let Some(d) = answer.strip_prefix("dir:") {
            Answer::Dir(PathBuf::from(d))
        } else if let Some(d) = answer.strip_prefix("flat:") {
            Answer::Flat(PathBuf::from(d))
        } else if let Some(spec) = answer.strip_prefix("drip:") {
            let parts: Vec<&str> = spec.split(':').collect();
            let [len, chunk, ms] = parts[..] else {
                return Err(format!("rule '{rule}': drip takes LEN:CHUNK:MS"));
            };
            let (len, chunk, ms) = (
                number(len, "LEN", rule)?,
                number(chunk, "CHUNK", rule)?,
                number(ms, "MS", rule)?,
            );
            if chunk == 0 && len > 0 {
                return Err(format!("rule '{rule}': CHUNK must be at least 1"));
            }
            Answer::Drip { len, chunk, ms }
        } else {
            return Err(format!("rule '{rule}': unknown answer '{answer}'"));
        };
        rules.push((path.to_owned(), answer));
    }
    Ok(rules)
}

/// The answer for a request path: its own rule, else `*`, else 404.
fn pick<'a>(rules: &'a Rules, path: &str) -> &'a Answer {
    rules
        .iter()
        .find(|(p, _)| p == path)
        .or_else(|| rules.iter().find(|(p, _)| p == "*"))
        .map_or(&Answer::NotFound, |(_, a)| a)
}

/// The file a `dir:` or `flat:` answer serves, if there is one.
fn file_for(answer: &Answer, path: &str) -> Option<PathBuf> {
    let rel = Path::new(path.trim_start_matches('/'));
    let f = match answer {
        Answer::Dir(d) => {
            if !rel.components().all(|c| matches!(c, Component::Normal(_))) {
                return None;
            }
            d.join(rel)
        }
        Answer::Flat(d) => d.join(rel.file_name()?),
        _ => return None,
    };
    f.is_file().then_some(f)
}

/// Reads the request head; the method and the path without its query.
fn request(stream: &mut TcpStream) -> Option<(String, String)> {
    let mut head = Vec::new();
    let mut buf = [0u8; 4096];
    while !head.windows(4).any(|w| w == b"\r\n\r\n") && !head.windows(2).any(|w| w == b"\n\n") {
        let n = stream.read(&mut buf).ok()?;
        if n == 0 || head.len() > 65536 {
            break;
        }
        head.extend_from_slice(&buf[..n]);
    }
    let line = String::from_utf8_lossy(&head);
    let mut words = line.lines().next()?.split_whitespace();
    let method = words.next()?.to_owned();
    let target = words.next()?;
    let path = target.split(['?', '#']).next().unwrap_or("/").to_owned();
    Some((method, path))
}

fn status(stream: &mut TcpStream, code: &str, len: u64) -> std::io::Result<()> {
    write!(
        stream,
        "HTTP/1.0 {code}\r\nContent-Length: {len}\r\nConnection: close\r\n\r\n"
    )?;
    stream.flush()
}

/// Answers one connection. Write errors (the client went) end it quietly.
fn handle(mut stream: TcpStream, rules: &Rules) -> std::io::Result<()> {
    let Some((method, path)) = request(&mut stream) else {
        return Ok(());
    };
    if method != "GET" && method != "HEAD" {
        return status(&mut stream, "501 Not Implemented", 0);
    }
    let body = method == "GET";
    let answer = pick(rules, &path);
    match answer {
        Answer::NotFound => status(&mut stream, "404 Not Found", 0),
        Answer::Dir(_) | Answer::Flat(_) => match file_for(answer, &path).map(std::fs::read) {
            Some(Ok(data)) => {
                status(&mut stream, "200 OK", data.len() as u64)?;
                if body {
                    stream.write_all(&data)?;
                }
                stream.flush()
            }
            _ => status(&mut stream, "404 Not Found", 0),
        },
        Answer::Drip { len, chunk, ms } => {
            status(&mut stream, "200 OK", *len)?;
            if !body {
                return Ok(());
            }
            let mut left = *len;
            while left > 0 {
                let n = left.min(*chunk);
                stream.write_all(&vec![b'x'; n as usize])?;
                stream.flush()?;
                left -= n;
                if left > 0 {
                    std::thread::sleep(Duration::from_millis(*ms));
                }
            }
            Ok(())
        }
    }
}

/// Answers every connection on LISTENER, each on its own thread, until the
/// process is killed or the listener fails.
fn serve(listener: TcpListener, rules: Rules) -> std::io::Error {
    let rules = std::sync::Arc::new(rules);
    loop {
        match listener.accept() {
            Ok((stream, _)) => {
                let rules = std::sync::Arc::clone(&rules);
                std::thread::spawn(move || {
                    let _ = handle(stream, &rules);
                });
            }
            Err(e) => return e,
        }
    }
}

/// `og dev stub-http RULE...`: prints the port, then serves until killed.
/// Exit 2 on a bad rule, 1 when the server cannot listen or stops.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let rules = match parse(args) {
        Ok(r) => r,
        Err(e) => {
            let _ = writeln!(err, "{USAGE}\nog dev stub-http: error: {e}");
            return 2;
        }
    };
    let listener = match TcpListener::bind(("127.0.0.1", 0)) {
        Ok(l) => l,
        Err(e) => {
            let _ = writeln!(err, "og dev stub-http: cannot listen on 127.0.0.1: {e}");
            return 1;
        }
    };
    let port = match listener.local_addr() {
        Ok(a) => a.port(),
        Err(e) => {
            let _ = writeln!(err, "og dev stub-http: {e}");
            return 1;
        }
    };
    if writeln!(out, "{port}").and_then(|()| out.flush()).is_err() {
        let _ = writeln!(err, "og dev stub-http: cannot write the port to stdout");
        return 1;
    }
    let e = serve(listener, rules);
    let _ = writeln!(err, "og dev stub-http: stopped: {e}");
    1
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Instant;

    fn rules(args: &[&str]) -> Rules {
        parse(&args.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>()).expect("rules parse")
    }

    /// A server on a free port with RULES; its port.
    fn start(args: &[&str]) -> u16 {
        let listener = TcpListener::bind(("127.0.0.1", 0)).unwrap();
        let port = listener.local_addr().unwrap().port();
        let r = rules(args);
        std::thread::spawn(move || serve(listener, r));
        port
    }

    /// Sends REQUEST and reads until the server closes or WAIT passes.
    fn get(port: u16, request: &str, wait: Duration) -> (Vec<u8>, bool) {
        let mut s = TcpStream::connect(("127.0.0.1", port)).unwrap();
        s.write_all(request.as_bytes()).unwrap();
        s.set_read_timeout(Some(Duration::from_millis(100)))
            .unwrap();
        let start = Instant::now();
        let mut got = Vec::new();
        let mut buf = [0u8; 65536];
        loop {
            match s.read(&mut buf) {
                Ok(0) => return (got, true),
                Ok(n) => got.extend_from_slice(&buf[..n]),
                Err(_) if start.elapsed() < wait => {}
                Err(_) => return (got, false),
            }
        }
    }

    fn split(resp: &[u8]) -> (String, Vec<u8>) {
        let at = resp
            .windows(4)
            .position(|w| w == b"\r\n\r\n")
            .expect("a header");
        (
            String::from_utf8_lossy(&resp[..at]).into_owned(),
            resp[at + 4..].to_vec(),
        )
    }

    const LONG: Duration = Duration::from_secs(10);

    #[test]
    fn rules_parse_and_bad_ones_are_refused() {
        assert_eq!(
            rules(&["/a=drip:10:5:0", "*=404"]),
            vec![
                (
                    "/a".to_owned(),
                    Answer::Drip {
                        len: 10,
                        chunk: 5,
                        ms: 0
                    }
                ),
                ("*".to_owned(), Answer::NotFound)
            ]
        );
        for bad in [
            "--og-scaffold-probe",
            "a=404",
            "/a=drip:1:2",
            "/a=drip:x:1:0",
            "/a=drip:5:0:0",
            "/a=teapot",
        ] {
            assert!(parse(&[bad.to_owned()]).is_err(), "{bad} was accepted");
        }
        assert!(parse(&[]).is_err());
        let mut out = Vec::new();
        let mut err = Vec::new();
        assert_eq!(cli(&["--og-scaffold-probe".into()], &mut out, &mut err), 2);
        assert!(out.is_empty());
    }

    #[test]
    fn a_path_without_a_rule_is_404_and_star_answers_the_rest() {
        let port = start(&["/whole.bin=drip:7:7:0"]);
        let (resp, closed) = get(port, "GET /other HTTP/1.1\r\nHost: x\r\n\r\n", LONG);
        assert!(closed && split(&resp).0.starts_with("HTTP/1.0 404 "));
        let port = start(&["/a=drip:3:3:0", "*=drip:5:5:0"]);
        let (resp, _) = get(port, "GET /b?q=1 HTTP/1.1\r\n\r\n", LONG);
        assert_eq!(split(&resp).1, b"xxxxx");
        let (resp, _) = get(port, "GET /a HTTP/1.1\r\n\r\n", LONG);
        assert_eq!(split(&resp).1, b"xxx");
    }

    #[test]
    fn a_whole_drip_is_sent_and_closed_and_head_has_no_body() {
        let port = start(&["*=drip:100000:100000:0"]);
        let (resp, closed) = get(port, "GET /whole.bin HTTP/1.1\r\n\r\n", LONG);
        let (head, body) = split(&resp);
        assert!(closed);
        assert!(head.starts_with("HTTP/1.0 200 OK") && head.contains("Content-Length: 100000"));
        assert_eq!(body.len(), 100000);
        let (resp, closed) = get(port, "HEAD /whole.bin HTTP/1.1\r\n\r\n", LONG);
        assert!(closed && split(&resp).1.is_empty());
    }

    #[test]
    fn a_stalled_drip_sends_its_first_chunk_and_then_holds_the_connection() {
        let port = start(&["/stall.bin=drip:100000:1000:60000"]);
        let (resp, closed) = get(
            port,
            "GET /stall.bin HTTP/1.1\r\n\r\n",
            Duration::from_millis(800),
        );
        let (head, body) = split(&resp);
        assert!(!closed, "the stalled transfer was closed");
        assert!(head.contains("Content-Length: 100000"));
        assert_eq!(body.len(), 1000);
    }

    #[test]
    fn dir_and_flat_serve_files_and_refuse_dot_dot() {
        let root = std::env::temp_dir().join(format!("og-stub-http-{}", std::process::id()));
        std::fs::create_dir_all(root.join("v1")).unwrap();
        std::fs::write(root.join("v1/SUMS"), b"sums\n").unwrap();
        std::fs::write(root.join("top"), b"top\n").unwrap();
        let d = format!("*=dir:{}", root.join("v1").display());
        let port = start(&[d.as_str()]);
        let (resp, _) = get(port, "GET /SUMS HTTP/1.1\r\n\r\n", LONG);
        assert_eq!(split(&resp).1, b"sums\n");
        let (resp, _) = get(port, "GET /../top HTTP/1.1\r\n\r\n", LONG);
        assert!(split(&resp).0.starts_with("HTTP/1.0 404 "));
        let (resp, _) = get(port, "GET /none HTTP/1.1\r\n\r\n", LONG);
        assert!(split(&resp).0.starts_with("HTTP/1.0 404 "));
        let f = format!("*=flat:{}", root.join("v1").display());
        let port = start(&[f.as_str()]);
        let (resp, _) = get(port, "GET /deep/path/SUMS HTTP/1.1\r\n\r\n", LONG);
        assert_eq!(split(&resp).1, b"sums\n");
        let _ = std::fs::remove_dir_all(&root);
    }
}
