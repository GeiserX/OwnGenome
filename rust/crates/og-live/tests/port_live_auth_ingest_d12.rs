//! The D12 refusals of `og live serve`. Without a password the site stays
//! a loopback site: it refuses to start with the flags of a site reached
//! through a proxy or over the network (`--trust-proxy`, `--allowed-host`,
//! `--tls-cert`, `--tls-key`, `--cookie-secure`), and it refuses a request
//! that says it came through a proxy (`Forwarded`, `X-Forwarded-*`,
//! `X-Real-IP`). The ingest port takes a message only with the token the
//! server made (`STATE_DIR/ingest-token`), and listens on loopback only. A
//! site with a password that serves plain HTTP off loopback warns as it
//! starts. Each check failed on the parity commit, before the refusals.

#[path = "port_live_server_util.rs"]
mod util;

use og_live::auth::{self, Auth};
use og_live::server::{self, Config, Server};
use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::sync::Arc;
use util::*;

const PASSWORD: &str = "synthetic test password 7";

fn args(root: &str, extra: &[&str]) -> Vec<String> {
    let mut a = vec!["--root".to_owned(), root.to_owned()];
    a.extend(extra.iter().map(|x| (*x).to_owned()));
    a
}

fn hash() -> String {
    auth::hash_password(PASSWORD, 1000).unwrap()
}

/// The site of `w` with an ingest listener on a free port.
fn ingest_config(w: &Work) -> Config {
    let mut c = config(w);
    c.ingest = Some(("127.0.0.1".into(), 0));
    c
}

fn post(port: u16, path: &str, body: &str) -> Resp {
    let head = format!(
        "POST {path} HTTP/1.1\r\nHost: 127.0.0.1\r\nContent-Type: application/json\r\n\
         Content-Length: {}\r\nConnection: close\r\n\r\n",
        body.len()
    );
    raw(port, &head, body.as_bytes())
}

fn stored(w: &Work) -> usize {
    fs::read_to_string(events_path(w, RUN))
        .map(|t| t.lines().count())
        .unwrap_or(0)
}

fn token_file(w: &Work) -> std::path::PathBuf {
    w.state().join("ingest-token")
}

const MSG: &str = r#"{"event": "pgp_step", "step": "06", "state": "ok"}"#;

#[test]
fn without_a_password_the_proxy_tls_and_cookie_flags_are_refused() {
    let w = Work::new("d12-flags");
    make_root(&w);
    let root = w.root().to_string_lossy().into_owned();
    for extra in [
        &["--trust-proxy", "127.0.0.1"][..],
        &["--allowed-host", "live.example.org"],
        &["--tls-cert", "/c.pem", "--tls-key", "/k.pem"],
        &["--tls-cert", "/c.pem"],
        &["--cookie-secure"],
    ] {
        match server::configure_with(&args(&root, extra), None) {
            Err((2, m)) => assert!(m.contains("password"), "{extra:?}: {m}"),
            Err((code, m)) => panic!("{extra:?}: exit {code}, not 2: {m}"),
            Ok(_) => panic!("{extra:?} was taken without a password"),
        }
    }
    let conf = w.dir.join("live.conf");
    write(&conf, "cookie_secure = true\n");
    let from_file =
        server::configure_with(&args(&root, &["--config", &conf.to_string_lossy()]), None);
    assert_eq!(
        from_file.err().map(|e| e.0),
        Some(2),
        "cookie_secure from --config"
    );
    // With a password each is taken.
    let h = hash();
    for extra in [
        &["--trust-proxy", "127.0.0.1"][..],
        &["--allowed-host", "live.example.org"],
        &["--cookie-secure"],
    ] {
        assert!(
            server::configure_with(&args(&root, extra), Some(&h)).is_ok(),
            "{extra:?}"
        );
    }
}

#[test]
fn without_a_password_a_forwarded_request_is_refused() {
    let w = Work::new("d12-forwarded");
    make_root(&w);
    let srv = start(config(&w));
    let p = srv.port();
    for (k, v) in [
        ("X-Forwarded-For", "192.0.2.1"),
        ("Forwarded", "for=192.0.2.1"),
        ("X-Forwarded-Host", "live.example.org"),
        ("X-Forwarded-Proto", "https"),
        ("X-Real-IP", "192.0.2.1"),
    ] {
        for path in ["/", &format!("/runs/{RUN}/report.html"), "/runs.json"] {
            let r = request(p, "GET", path, &[(k, v)]);
            assert_eq!(r.status, 403, "{k} on {path}");
            assert!(
                policy_ok(&r) && !r.text().contains("<html"),
                "{k} on {path}"
            );
        }
    }
    assert_eq!(get(p, &format!("/runs/{RUN}/report.html")).status, 200);
    assert!(
        log_rows(&w)
            .iter()
            .any(|r| r[3] == "forwarded-refused" && r[6] == "403")
    );
    // With a password a forwarded request meets the login.
    let w2 = Work::new("d12-forwarded-pw");
    make_root(&w2);
    let mut c = config(&w2);
    c.auth = Arc::new(Auth::new(Some(hash())));
    let srv2 = start(c);
    let r = request(srv2.port(), "GET", "/", &[("X-Forwarded-For", "192.0.2.1")]);
    assert_eq!(
        (r.status, r.header("location")),
        (303, Some("/login?next=/"))
    );
}

#[test]
fn the_ingest_port_takes_only_the_server_made_token() {
    let w = Work::new("d12-token");
    make_root(&w);
    let srv: Server = start(ingest_config(&w));
    let ip = srv.ingest_port().expect("an ingest listener");
    let plain = post(ip, &format!("/ingest/{RUN}"), MSG);
    assert_eq!(
        plain.status,
        403,
        "a message with no token: {}",
        plain.text()
    );
    let made_up = "A".repeat(43);
    assert_eq!(
        post(ip, &format!("/{made_up}/ingest/{RUN}"), MSG).status,
        403
    );
    assert_eq!(post(ip, &format!("/x/ingest/{RUN}"), MSG).status, 403);
    assert_eq!(stored(&w), 0, "fake progress was stored");
    let token = fs::read_to_string(token_file(&w)).expect("STATE_DIR/ingest-token");
    let token = token.trim();
    assert_eq!(token.len(), 43);
    assert!(
        token
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'-' || b == b'_')
    );
    let mode = fs::metadata(token_file(&w)).unwrap().permissions().mode() & 0o777;
    assert_eq!(mode, 0o600, "the token file is the owner's only");
    assert_eq!(post(ip, &format!("/{token}/ingest/{RUN}"), MSG).status, 200);
    // A URL ending in '/' gives a doubled slash before ingest.
    assert_eq!(
        post(ip, &format!("/{token}//ingest/{RUN}"), MSG).status,
        200
    );
    assert_eq!(stored(&w), 2);
    let near: String = token.chars().rev().collect();
    assert_eq!(post(ip, &format!("/{near}/ingest/{RUN}"), MSG).status, 403);
    assert_eq!(stored(&w), 2);
}

#[test]
fn a_restart_keeps_the_token_and_a_token_file_open_to_others_is_replaced() {
    let w = Work::new("d12-restart");
    make_root(&w);
    let first = {
        let _srv = start(ingest_config(&w));
        fs::read_to_string(token_file(&w)).expect("STATE_DIR/ingest-token")
    };
    let second = {
        let _srv = start(ingest_config(&w));
        fs::read_to_string(token_file(&w)).unwrap()
    };
    assert_eq!(
        first, second,
        "a restart must not cut off a run that is posting"
    );
    fs::set_permissions(token_file(&w), fs::Permissions::from_mode(0o644)).unwrap();
    let third = {
        let _srv = start(ingest_config(&w));
        fs::read_to_string(token_file(&w)).unwrap()
    };
    assert_ne!(third, first, "a token others could read was kept");
    let mode = fs::metadata(token_file(&w)).unwrap().permissions().mode() & 0o777;
    assert_eq!(mode, 0o600);
    // A token file someone wrote by hand, not of the server's shape.
    fs::write(token_file(&w), "chosen-by-a-client\n").unwrap();
    fs::set_permissions(token_file(&w), fs::Permissions::from_mode(0o600)).unwrap();
    let srv = start(ingest_config(&w));
    let fourth = fs::read_to_string(token_file(&w)).unwrap();
    assert_ne!(fourth.trim(), "chosen-by-a-client");
    let ip = srv.ingest_port().unwrap();
    assert_eq!(
        post(ip, &format!("/chosen-by-a-client/ingest/{RUN}"), MSG).status,
        403
    );
}

#[test]
fn the_ingest_port_listens_on_loopback_only() {
    let w = Work::new("d12-ingest-bind");
    make_root(&w);
    let root = w.root().to_string_lossy().into_owned();
    let h = hash();
    for pw in [None, Some(h.as_str())] {
        for bind in ["0.0.0.0", "::", "192.0.2.1"] {
            match server::configure_with(&args(&root, &["--ingest-bind", bind]), pw) {
                Err((2, m)) => assert!(m.contains("--ingest-bind"), "{m}"),
                Err((code, m)) => panic!("--ingest-bind {bind}: exit {code}: {m}"),
                Ok(c) => panic!(
                    "--ingest-bind {bind} was taken (notes: {:?})",
                    server::start_notes(&c)
                ),
            }
        }
    }
    for bind in ["127.0.0.1", "::1", "localhost"] {
        assert!(
            server::configure_with(&args(&root, &["--ingest-bind", bind]), None).is_ok(),
            "{bind}"
        );
    }
}

#[test]
fn a_password_site_on_plain_http_off_loopback_warns_as_it_starts() {
    let w = Work::new("d12-cleartext");
    make_root(&w);
    let root = w.root().to_string_lossy().into_owned();
    let h = hash();
    let warned = |extra: &[&str]| {
        let c = server::configure_with(&args(&root, extra), Some(&h)).unwrap();
        server::start_notes(&c)
            .iter()
            .any(|l| l.contains("WARNING") && l.contains("plain HTTP"))
    };
    assert!(
        warned(&["--bind", "0.0.0.0"]),
        "no warning for a cleartext bind off loopback"
    );
    assert!(warned(&["--bind", "0.0.0.0", "--cookie-secure"]));
    assert!(!warned(&[
        "--bind",
        "0.0.0.0",
        "--tls-cert",
        "/c.pem",
        "--tls-key",
        "/k.pem"
    ]));
    assert!(!warned(&["--bind", "127.0.0.1"]));
}
