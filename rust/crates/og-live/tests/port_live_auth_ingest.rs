//! `og live serve`'s login, sessions, HTTPS and ingest listener on the
//! synthetic DEMO-001 run (tests/fixtures/live). These are the checks of
//! the Python server's `tests/test_owngenome_live.py` that need the login
//! or the ingest port, ported case for case: the password hash, the login
//! and its limit on failed logins, the session cookie, the Host check with
//! `--allowed-host`, what the ingest port takes and stores, an event
//! stream that ends with its session, the plugin routes behind the login
//! and a silent client on the HTTPS port. The D12 refusals are
//! `port_live_auth_ingest_d12.rs`; the login form in Chrome is
//! `port_live_auth_ingest_browser.rs`.

#[path = "port_live_server_util.rs"]
mod util;

use og_live::auth::{self, Auth};
use og_live::server::{self, Config, Server};
use std::fs;
use std::io::{BufRead, BufReader, Cursor, Read, Write};
use std::net::TcpStream;
use std::path::Path;
use std::process::Command;
use std::sync::{Arc, Barrier};
use std::thread;
use std::time::{Duration, Instant};
use util::*;

const PASSWORD: &str = "synthetic test password 7";
/// The same password hashed by `bin/owngenome_live.py`'s format with
/// Python's `hashlib.pbkdf2_hmac` (salt `00112233445566778899aabbccddeeff`,
/// 1000 rounds): a site made by the Python server keeps its login.
const PY_HASH: &str = "pbkdf2_sha256$1000$00112233445566778899aabbccddeeff$\
                       2d0f8cdad5ceae90d47e704316f775d61bd586f3c81c2a18614c665f0cbc80b2";

/// The site of `w` with a password (a hash of 1000 rounds, so a debug
/// build checks it fast), `--allowed-host live.example.org`, an ingest
/// listener on a free port and a quick heartbeat.
fn locked_config(w: &Work) -> Config {
    let mut c = config(w);
    c.auth = Arc::new(Auth::new(Some(
        auth::hash_password(PASSWORD, 1000).unwrap(),
    )));
    c.allowed_host = vec!["live.example.org".into()];
    c.ingest = Some(("127.0.0.1".into(), 0));
    c.heartbeat = Duration::from_millis(200);
    c
}

fn locked_site(name: &str) -> (Work, Server) {
    let w = Work::new(name);
    make_root(&w);
    let srv = start(locked_config(&w));
    (w, srv)
}

/// `POST /login` with `password`: the answer and the session cookie it set.
fn login(port: u16, password: &str) -> (Resp, Option<String>) {
    let body = format!("password={}&next=/", password.replace(' ', "+"));
    let head = format!(
        "POST /login HTTP/1.1\r\nHost: 127.0.0.1:{port}\r\nContent-Type: \
         application/x-www-form-urlencoded\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
        body.len()
    );
    let r = raw(port, &head, body.as_bytes());
    let cookie = r
        .header("set-cookie")
        .and_then(|c| c.strip_prefix("owngenome_live="))
        .map(|c| c.split(';').next().unwrap_or("").to_owned())
        .filter(|c| !c.is_empty());
    (r, cookie)
}

fn with_cookie(port: u16, method: &str, path: &str, cookie: &str) -> Resp {
    let c = format!("owngenome_live={cookie}");
    request(port, method, path, &[("Cookie", c.as_str())])
}

/// The path a message of `run` is posted to: after the token the server
/// made (`STATE_DIR/ingest-token`).
fn ingest_path(w: &Work, run: &str) -> String {
    let token = fs::read_to_string(w.state().join("ingest-token")).expect("the ingest token");
    format!("/{}/ingest/{run}", token.trim())
}

/// Posts `body` to the ingest port at `path`.
fn post(port: u16, path: &str, body: &str, extra: &[(&str, &str)]) -> Resp {
    let mut head = format!(
        "POST {path} HTTP/1.1\r\nHost: 127.0.0.1\r\nContent-Type: application/json\r\n\
         Content-Length: {}\r\n",
        body.len()
    );
    for (k, v) in extra {
        head.push_str(&format!("{k}: {v}\r\n"));
    }
    head.push_str("Connection: close\r\n\r\n");
    raw(port, &head, body.as_bytes())
}

/// The first line of the answer to a hand-made request, or "" when the
/// server closed without one.
fn status_line(port: u16, head: &str, body: &[u8]) -> String {
    let mut s = TcpStream::connect(("127.0.0.1", port)).unwrap();
    s.set_read_timeout(Some(Duration::from_secs(10))).unwrap();
    let _ = s.write_all(head.as_bytes());
    let _ = s.write_all(body);
    let mut buf = [0u8; 200];
    let n = s.read(&mut buf).unwrap_or(0);
    String::from_utf8_lossy(&buf[..n])
        .split("\r\n")
        .next()
        .unwrap_or("")
        .to_owned()
}

#[test]
fn hash_password_writes_a_600000_round_hash_and_a_python_hash_still_logs_in() {
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let mut input = Cursor::new(format!("{PASSWORD}\n"));
    let code = auth::hash_password_cli(&[], &mut input, false, &mut out, &mut err);
    let h = String::from_utf8_lossy(&out).trim().to_owned();
    let all = format!("{h}{}", String::from_utf8_lossy(&err));
    assert_eq!(code, 0, "{all}");
    assert!(h.starts_with("pbkdf2_sha256$600000$"), "{h}");
    assert!(!all.contains(PASSWORD));
    assert!(auth::verify_password(PASSWORD, &h));
    assert!(!auth::verify_password(&format!("{PASSWORD}x"), &h));
    assert!(auth::verify_password(PASSWORD, PY_HASH));
    assert!(!auth::verify_password("synthetic test password 8", PY_HASH));
    let mut short = Cursor::new("seven77\n");
    let code = auth::hash_password_cli(&[], &mut short, false, &mut Vec::new(), &mut err);
    assert_eq!(code, 2);
    let extra = ["--rounds".to_owned()];
    let code = auth::hash_password_cli(
        &extra,
        &mut Cursor::new(""),
        false,
        &mut Vec::new(),
        &mut err,
    );
    assert_eq!(code, 2);
}

#[test]
fn the_server_reads_only_a_hash() {
    let w = Work::new("hash-read");
    make_root(&w);
    let root = w.root().to_string_lossy().into_owned();
    let args = |extra: &[&str]| -> Vec<String> {
        let mut a = vec!["--root".to_owned(), root.clone()];
        a.extend(extra.iter().map(|x| (*x).to_owned()));
        a
    };
    let plain = w.dir.join("password");
    write(&plain, &format!("{PASSWORD}\n"));
    let plain = plain.to_string_lossy().into_owned();
    match server::configure_with(&args(&["--password-file", &plain]), None) {
        Err((1, m)) => assert!(m.contains("password hash") && !m.contains(PASSWORD), "{m}"),
        other => panic!(
            "a plain password was taken: {:?}",
            other.map(|c| c.auth.hash.clone())
        ),
    }
    let file = w.dir.join("hash");
    write(&file, &format!("{PY_HASH}\n"));
    let c =
        server::configure_with(&args(&["--password-file", &file.to_string_lossy()]), None).unwrap();
    assert_eq!(c.auth.hash.as_deref(), Some(PY_HASH));
    let c = server::configure_with(&args(&[]), Some(PY_HASH)).unwrap();
    assert_eq!(c.auth.hash.as_deref(), Some(PY_HASH));
    // --password-file wins over the variable, which is read only without it.
    let c = server::configure_with(
        &args(&["--password-file", &file.to_string_lossy()]),
        Some("x"),
    )
    .unwrap();
    assert_eq!(c.auth.hash.as_deref(), Some(PY_HASH));
    let c = server::configure_with(
        &args(&["--session-hours", "2", "--login-max-failures", "3"]),
        None,
    )
    .unwrap();
    assert_eq!((c.auth.session_hours, c.auth.login_max_failures), (2.0, 3));
    let (code, _) = server::configure_with(&args(&["--bind", "0.0.0.0"]), None)
        .err()
        .unwrap();
    assert_eq!(code, 2, "a non-loopback --bind without a password hash");
    assert!(server::configure_with(&args(&["--bind", "0.0.0.0"]), Some(PY_HASH)).is_ok());
}

#[test]
fn login_cookie_and_the_access_log() {
    let (w, srv) = locked_site("login");
    let p = srv.port();
    let page = get(p, &format!("/runs/{RUN}/"));
    assert_eq!(page.status, 303);
    assert_eq!(page.header("location"), Some("/login?next=/runs/DEMO-001/"));
    assert_eq!(get(p, &format!("/runs/{RUN}/state.json")).status, 401);
    assert_eq!(get(p, "/runs.json").status, 401);
    let form = get(p, "/login?next=/runs/DEMO-001/");
    assert_eq!(form.status, 200);
    assert!(form.text().contains("/runs/DEMO-001/") && policy_ok(&form));
    let (bad, none) = login(p, "wrong password");
    assert_eq!((bad.status, none), (401, None));
    assert!(
        log_rows(&w)
            .iter()
            .any(|r| r[3] == "login-fail" && r[6] == "401")
    );
    let (ok, cookie) = login(p, PASSWORD);
    let cookie = cookie.expect("a session cookie");
    let set = ok.header("set-cookie").unwrap_or("");
    assert_eq!(ok.status, 303);
    assert_eq!(ok.header("location"), Some("/"));
    assert!(
        set.contains("HttpOnly") && set.contains("SameSite=Strict") && set.contains("Path=/"),
        "{set}"
    );
    assert!(
        set.contains("Max-Age=43200") && !set.contains("Secure"),
        "{set}"
    );
    assert!(
        cookie.len() == 43
            && cookie
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b == b'-' || b == b'_')
    );
    assert_eq!(
        with_cookie(p, "GET", &format!("/runs/{RUN}/state.json"), &cookie).status,
        200
    );
    assert_eq!(with_cookie(p, "GET", "/", "made-up").status, 303);
    // The log keeps 8 characters of the session's hash, never the cookie.
    let s8 = &auth::session_id(&cookie)[..8];
    let rows = log_rows(&w);
    assert!(
        rows.iter().any(|r| r[2] == s8 && r[3] == "login"),
        "{rows:?}"
    );
    assert!(
        rows.iter().any(|r| r[2] == s8 && r[3] == "state"),
        "{rows:?}"
    );
    let text = fs::read_to_string(w.state().join("access.log")).unwrap();
    assert!(!text.contains(&cookie) && !text.contains(PASSWORD));
    // A next that leaves the site goes home.
    let body = format!(
        "password={}&next=//evil.example.net/",
        PASSWORD.replace(' ', "+")
    );
    let head = format!(
        "POST /login HTTP/1.1\r\nHost: 127.0.0.1:{p}\r\nContent-Length: {}\r\n\r\n",
        body.len()
    );
    assert_eq!(raw(p, &head, body.as_bytes()).header("location"), Some("/"));
    // Logout ends the session and clears the cookie.
    let out = with_cookie(p, "POST", "/logout", &cookie);
    assert_eq!((out.status, out.header("location")), (303, Some("/login")));
    assert!(out.header("set-cookie").unwrap_or("").contains("Max-Age=0"));
    assert_eq!(
        with_cookie(p, "GET", &format!("/runs/{RUN}/state.json"), &cookie).status,
        401
    );
}

#[test]
fn a_session_ends_after_its_hours_or_idle_minutes() {
    let w = Work::new("expire");
    make_root(&w);
    let mut c = locked_config(&w);
    let mut a = Auth::new(c.auth.hash.clone());
    a.session_idle_minutes = 0.01; // 0.6 s
    c.auth = Arc::new(a);
    let srv = start(c);
    let p = srv.port();
    let (_, cookie) = login(p, PASSWORD);
    let cookie = cookie.unwrap();
    assert_eq!(with_cookie(p, "GET", "/runs.json", &cookie).status, 200);
    thread::sleep(Duration::from_millis(900));
    assert_eq!(with_cookie(p, "GET", "/runs.json", &cookie).status, 401);
    let w2 = Work::new("expire-hours");
    make_root(&w2);
    let mut c = locked_config(&w2);
    let mut a = Auth::new(c.auth.hash.clone());
    a.session_hours = 0.0002; // 0.72 s
    c.auth = Arc::new(a);
    let srv = start(c);
    let p = srv.port();
    let (_, cookie) = login(p, PASSWORD);
    let cookie = cookie.unwrap();
    for _ in 0..3 {
        assert_eq!(with_cookie(p, "GET", "/runs.json", &cookie).status, 200);
        thread::sleep(Duration::from_millis(150));
    }
    thread::sleep(Duration::from_millis(600));
    assert_eq!(
        with_cookie(p, "GET", "/runs.json", &cookie).status,
        401,
        "use does not extend the hours"
    );
}

#[test]
fn five_wrong_passwords_then_429_even_for_the_right_one() {
    let (w, srv) = locked_site("limit");
    let p = srv.port();
    let codes: Vec<u16> = (0..6)
        .map(|i| login(p, &format!("wrong {i}")).0.status)
        .collect();
    assert_eq!(codes, [401, 401, 401, 401, 401, 429]);
    let (r, cookie) = login(p, PASSWORD);
    assert_eq!((r.status, cookie), (429, None));
    assert_eq!(r.header("retry-after"), Some("900"));
    assert!(
        log_rows(&w)
            .iter()
            .any(|r| r[3] == "login-blocked" && r[6] == "429")
    );
}

#[test]
fn parallel_wrong_passwords_get_at_most_five_checks() {
    let (_w, srv) = locked_site("burst");
    let p = srv.port();
    let n = 12;
    let gate = Arc::new(Barrier::new(n));
    let hands: Vec<_> = (0..n)
        .map(|i| {
            let gate = Arc::clone(&gate);
            thread::spawn(move || {
                gate.wait();
                login(p, &format!("parallel wrong {i}")).0.status
            })
        })
        .collect();
    let got: Vec<u16> = hands.into_iter().map(|h| h.join().unwrap()).collect();
    assert!(got.iter().all(|c| *c == 401 || *c == 429), "{got:?}");
    assert!(got.iter().filter(|c| **c == 401).count() <= 5, "{got:?}");
}

#[test]
fn a_login_body_with_a_bad_length_is_413_before_it_is_read() {
    let (_w, srv) = locked_site("login-length");
    let p = srv.port();
    for n in ["-1", "abc", "4097"] {
        let head = format!(
            "POST /login HTTP/1.1\r\nHost: 127.0.0.1:{p}\r\nContent-Type: \
             application/x-www-form-urlencoded\r\nContent-Length: {n}\r\n\r\n"
        );
        let line = status_line(p, &head, b"password=x");
        assert!(line.contains(" 413 "), "Content-Length {n}: {line}");
    }
}

#[test]
fn host_check_with_an_allowed_host() {
    let (_w, srv) = locked_site("allowed-host");
    let p = srv.port();
    let (_, cookie) = login(p, PASSWORD);
    let c = format!("owngenome_live={}", cookie.unwrap());
    let evil = request(
        p,
        "GET",
        "/",
        &[("Cookie", &c), ("Host", "evil.example.net")],
    );
    let named = request(
        p,
        "GET",
        "/",
        &[("Cookie", &c), ("Host", "live.example.org:8443")],
    );
    assert_eq!((evil.status, named.status), (403, 200));
}

#[test]
fn the_ingest_port_takes_cut_messages_and_nothing_else() {
    let (w, srv) = locked_site("ingest");
    let p = srv.port();
    let ip = srv.ingest_port().expect("an ingest listener");
    let path = ingest_path(&w, RUN);
    for line in fixture_events() {
        assert_eq!(post(ip, &path, line.trim(), &[]).status, 200, "{line}");
    }
    let stored = fs::read_to_string(events_path(&w, RUN)).unwrap();
    assert_eq!(stored.lines().count(), fixture_events().len());
    let leaks = |text: &str| -> Vec<&str> {
        [
            "HP:0001250",
            "DEMO-001.vcf.gz",
            "bcftools roh",
            "commandLine",
        ]
        .into_iter()
        .filter(|x| text.contains(x))
        .collect()
    };
    assert!(leaks(&stored).is_empty(), "{:?}", leaks(&stored));
    assert!(stored.contains(r#""outdir":"/data/genome""#));
    let fixture = fs::read_to_string(repo().join("tests/fixtures/live/events.jsonl")).unwrap();
    assert_eq!(
        leaks(&fixture).len(),
        4,
        "control: the raw messages hold every leak"
    );
    // The browser port takes no message, even with a session.
    let (_, cookie) = login(p, PASSWORD);
    let c = format!("owngenome_live={}", cookie.unwrap());
    let head = format!(
        "POST /ingest/{RUN} HTTP/1.1\r\nHost: 127.0.0.1:{p}\r\nCookie: {c}\r\nContent-Length: 2\r\n\r\n"
    );
    assert_eq!(raw(p, &head, b"{}").status, 404);
    // Too large, by its Content-Length, before the body.
    let big = format!(
        "POST {path} HTTP/1.1\r\nHost: 127.0.0.1\r\nContent-Type: application/json\r\nContent-Length: {}\r\n\r\n",
        2 * 1024 * 1024
    );
    assert!(status_line(ip, &big, &[b'x'; 1024]).contains(" 413 "));
    let neg = format!("POST {path} HTTP/1.1\r\nHost: 127.0.0.1\r\nContent-Length: -1\r\n\r\n");
    assert!(status_line(ip, &neg, br#"{"event":"started"}"#).contains(" 413 "));
    let none = format!("POST {path} HTTP/1.1\r\nHost: 127.0.0.1\r\n\r\n");
    assert!(status_line(ip, &none, b"").contains(" 411 "));
    let nan = r#"{"event": "process_started", "trace": {"attempt": NaN, "exit": Infinity}}"#;
    assert_eq!(post(ip, &path, nan, &[]).status, 400);
    assert_eq!(post(ip, &path, "not json", &[]).status, 400);
    for bad in ["../x", "-x"] {
        assert_eq!(
            post(ip, &ingest_path(&w, bad), r#"{"event":"started"}"#, &[]).status,
            400,
            "{bad}"
        );
    }
    let origin = [("Origin", "https://evil.example.net")];
    assert_eq!(
        post(ip, &path, r#"{"event":"started"}"#, &origin).status,
        403
    );
    assert_eq!(request(ip, "GET", &path, &[]).status, 404);
    let n = fs::read_to_string(events_path(&w, RUN))
        .unwrap()
        .lines()
        .count();
    assert_eq!(n, fixture_events().len(), "a refused message was stored");
    // A trace with numbers where text goes is stored as text.
    let numeric =
        r#"{"event":"process_started","trace":{"name":123,"process":7,"tag":9,"status":1}}"#;
    assert_eq!(post(ip, &path, numeric, &[]).status, 200);
    let last = fs::read_to_string(events_path(&w, RUN)).unwrap();
    let last = last.lines().last().unwrap();
    assert!(
        last.contains(r#""name":"123""#) && last.contains(r#""status":"1""#),
        "{last}"
    );
    assert!(
        !fs::read_to_string(w.state().join("access.log"))
            .unwrap_or_default()
            .contains(PASSWORD)
    );
}

#[test]
fn an_open_stream_ends_with_its_session_and_sends_nothing_after() {
    let (w, srv) = locked_site("stream-end");
    let p = srv.port();
    let ip = srv.ingest_port().unwrap();
    let path = ingest_path(&w, RUN);
    for line in fixture_events() {
        assert_eq!(post(ip, &path, line.trim(), &[]).status, 200);
    }
    let n0 = fixture_events().len();
    let (_, cookie) = login(p, PASSWORD);
    let cookie = cookie.unwrap();
    let mut s = TcpStream::connect(("127.0.0.1", p)).unwrap();
    s.set_read_timeout(Some(Duration::from_secs(5))).unwrap();
    write!(
        s,
        "GET /runs/{RUN}/events HTTP/1.1\r\nHost: 127.0.0.1:{p}\r\nCookie: owngenome_live={cookie}\r\n\r\n"
    )
    .unwrap();
    let mut r = BufReader::new(s);
    let (mut ids, mut line) = (0, String::new());
    while ids < n0 {
        line.clear();
        assert!(
            r.read_line(&mut line).unwrap() > 0,
            "the stream closed early"
        );
        if line.starts_with("id: ") {
            ids += 1;
        }
    }
    assert_eq!(with_cookie(p, "POST", "/logout", &cookie).status, 303);
    let step =
        r#"{"event": "pgp_step", "step": "06", "state": "ok", "utcTime": "2026-01-02T10:06:40Z"}"#;
    assert_eq!(post(ip, &path, step, &[]).status, 200);
    let end = Instant::now() + Duration::from_secs(4);
    let mut closed = false;
    while Instant::now() < end {
        line.clear();
        match r.read_line(&mut line) {
            Ok(0) => {
                closed = true;
                break;
            }
            Ok(_) => assert!(
                !line.starts_with("id: "),
                "a message after the session ended: {line}"
            ),
            Err(_) => break,
        }
    }
    assert!(closed, "the stream stayed open after the session ended");
}

#[test]
fn the_plugin_routes_and_metrics_are_behind_the_login() {
    let w = Work::new("plugins-login");
    make_root(&w);
    let mut c = locked_config(&w);
    let tok = w.dir.join("metrics-token");
    write(&tok, "synthetic-metrics-token\n");
    c.metrics_token_file = Some(tok);
    let srv = start(c);
    let p = srv.port();
    for path in [
        format!("/runs/{RUN}/map.json"),
        format!("/runs/{RUN}/results.json"),
        "/family/fam/x.svg".into(),
    ] {
        assert_eq!(get(p, &path).status, 401, "{path}");
    }
    assert_eq!(
        request(
            p,
            "POST",
            &format!("/runs/{RUN}/map.json"),
            &[("Content-Length", "0")]
        )
        .status,
        404
    );
    assert_eq!(get(p, "/metrics").status, 401);
    let bearer = [("Authorization", "Bearer synthetic-metrics-token")];
    assert_eq!(request(p, "GET", "/metrics", &bearer).status, 200);
    let wrong = [("Authorization", "Bearer made-up")];
    assert_eq!(request(p, "GET", "/metrics", &wrong).status, 401);
    let (_, cookie) = login(p, PASSWORD);
    assert_eq!(
        with_cookie(p, "GET", "/metrics", &cookie.unwrap()).status,
        200
    );
    // The assets and the login form need no session.
    assert_eq!(get(p, "/assets/live.css").status, 200);
    assert_eq!(get(p, "/login").status, 200);
}

/// A self-signed certificate for 127.0.0.1 in `dir`, made by openssl;
/// `None` (a skip) without openssl, except in CI.
fn test_certificate(dir: &Path) -> Option<(std::path::PathBuf, std::path::PathBuf)> {
    let (cert, key) = (dir.join("c.pem"), dir.join("k.pem"));
    let made = Command::new("openssl")
        .args([
            "req",
            "-x509",
            "-newkey",
            "rsa:2048",
            "-nodes",
            "-days",
            "1",
            "-subj",
            "/CN=localhost",
        ])
        .arg("-keyout")
        .arg(&key)
        .arg("-out")
        .arg(&cert)
        .output();
    match made {
        Ok(o) if o.status.success() => Some((cert, key)),
        other => {
            assert!(
                std::env::var("CI").as_deref() != Ok("true"),
                "no openssl in CI: {other:?}"
            );
            eprintln!("[SKIP] no openssl to make a test certificate");
            None
        }
    }
}

/// `curl -k` to the site: the status and the headers.
fn curl(args: &[&str], url: &str) -> (String, String) {
    let o = Command::new("curl")
        .args([
            "-sk",
            "--max-time",
            "5",
            "-D",
            "-",
            "-o",
            "/dev/null",
            "-w",
            "\n%{http_code}",
        ])
        .args(args)
        .arg(url)
        .output()
        .expect("curl");
    let text = String::from_utf8_lossy(&o.stdout).into_owned();
    let code = text.lines().last().unwrap_or("").to_owned();
    (code, text)
}

#[test]
fn a_silent_client_does_not_hold_up_https_and_the_cookie_is_secure() {
    let w = Work::new("tls");
    make_root(&w);
    let Some((cert, key)) = test_certificate(&w.dir) else {
        return;
    };
    let mut c = locked_config(&w);
    c.tls = Some((cert, key));
    let srv = start(c);
    let p = srv.port();
    let silent = TcpStream::connect(("127.0.0.1", p)).unwrap();
    thread::sleep(Duration::from_millis(300));
    let url = format!("https://127.0.0.1:{p}/login");
    let (code, _) = curl(&[], &url);
    assert_eq!(code, "200", "a silent client held up an HTTPS request");
    let form = format!("password={}", PASSWORD.replace(' ', "+"));
    let (code, head) = curl(&["--data", &form], &url);
    assert_eq!(code, "303", "{head}");
    let set = head
        .lines()
        .find(|l| l.to_ascii_lowercase().starts_with("set-cookie:"))
        .unwrap_or("");
    assert!(set.contains("; Secure"), "{head}");
    drop(silent);
    // Plain HTTP on the HTTPS port gets no HTTP answer.
    let line = status_line(
        p,
        &format!("GET /login HTTP/1.1\r\nHost: 127.0.0.1:{p}\r\n\r\n"),
        b"",
    );
    assert!(!line.starts_with("HTTP/"), "{line}");
}
