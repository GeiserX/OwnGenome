//! `og live serve` started the way a person starts it, as a process of its
//! own, so its request threads write to the process's real stderr:
//!
//! 1. `scripts/owngenome-live.sh` with its default logging (no `--quiet`)
//!    answers `GET /` with 200. `og` once held the stderr lock for the
//!    whole command, and the line each request prints blocked the request
//!    for good: the site bound its port and never answered;
//! 2. `og live serve --quiet` answers `family.json` of a family whose
//!    `relationships.tsv` cannot be read with the fixed 500 "family not
//!    readable". That answer prints a note on stderr first, so the held
//!    lock blocked it even with `--quiet`;
//! 3. a post meant for the ingest port that reaches the browser port
//!    prints its request line with the token written `<token>`, never the
//!    token. Control: the same line holds the request's other parts.
//!
//! Every answer is awaited for at most 10 s, so a blocked request fails a
//! test instead of hanging it. `--assets` and `--steps` name the checkout's
//! files, because the test binary's `og` is not inside the checkout.

use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::time::{Duration, Instant};

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");
const WAIT: Duration = Duration::from_secs(10);

fn repo() -> PathBuf {
    Path::new(REPO).canonicalize().unwrap()
}

/// A fresh folder for one case.
fn dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-live-cli-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// A port on 127.0.0.1 that nothing listens on now.
fn free_port() -> u16 {
    TcpListener::bind("127.0.0.1:0")
        .unwrap()
        .local_addr()
        .unwrap()
        .port()
}

/// The checkout's assets, for a binary outside it.
fn asset_flags() -> Vec<String> {
    let r = repo();
    vec![
        "--assets".into(),
        r.join("assets/live").display().to_string(),
        "--steps".into(),
        r.join("assets/steps.tsv").display().to_string(),
    ]
}

/// The site's process, killed when the case ends; its stderr in a file.
struct Site {
    child: Child,
    err: PathBuf,
}

impl Drop for Site {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

impl Site {
    fn start(mut cmd: Command, work: &Path, port: u16) -> Site {
        let err = work.join("site.err");
        let child = cmd
            .stdin(Stdio::null())
            .stdout(Stdio::null())
            .stderr(std::fs::File::create(&err).unwrap())
            .spawn()
            .expect("the site starts");
        let site = Site { child, err };
        let end = Instant::now() + Duration::from_secs(30);
        while TcpStream::connect(("127.0.0.1", port)).is_err() {
            assert!(
                Instant::now() < end,
                "the site never listened: {}",
                site.stderr()
            );
            std::thread::sleep(Duration::from_millis(100));
        }
        site
    }

    fn stderr(&self) -> String {
        std::fs::read_to_string(&self.err).unwrap_or_default()
    }
}

/// The status and body of one request, or why there is none within [`WAIT`].
fn request(port: u16, method: &str, path: &str, body: &str) -> Result<(u16, String), String> {
    let mut s = TcpStream::connect(("127.0.0.1", port)).map_err(|e| e.to_string())?;
    s.set_read_timeout(Some(WAIT)).unwrap();
    write!(
        s,
        "{method} {path} HTTP/1.1\r\nHost: 127.0.0.1:{port}\r\nContent-Type: application/json\r\n\
         Content-Length: {}\r\nConnection: close\r\n\r\n{body}",
        body.len()
    )
    .map_err(|e| e.to_string())?;
    let mut all = Vec::new();
    s.read_to_end(&mut all).map_err(|e| {
        format!(
            "no whole answer within {WAIT:?} ({e}); got {} bytes",
            all.len()
        )
    })?;
    let text = String::from_utf8_lossy(&all).into_owned();
    let status = text
        .split(' ')
        .nth(1)
        .and_then(|c| c.parse().ok())
        .ok_or_else(|| format!("not an HTTP answer: {text:?}"))?;
    let body = text
        .split_once("\r\n\r\n")
        .map(|(_, b)| b.to_owned())
        .unwrap_or_default();
    Ok((status, body))
}

#[test]
fn the_default_launcher_answers_get_root() {
    let w = dir("launcher");
    let genome = w.join("genome");
    std::fs::create_dir_all(&genome).unwrap();
    let (port, ingest) = (free_port(), free_port());
    let mut cmd = Command::new("bash");
    cmd.arg(repo().join("scripts/owngenome-live.sh"))
        .args(asset_flags())
        .env("GENOME_DIR", &genome)
        .env("OG", env!("CARGO_BIN_EXE_og"))
        .env("OWNGENOME_LIVE_PORT", port.to_string())
        .env("OWNGENOME_LIVE_INGEST_PORT", ingest.to_string())
        .env("OWNGENOME_LIVE_VIEWER", "off")
        .env_remove("OWNGENOME_LIVE_STATE_DIR")
        .env_remove("OWNGENOME_LIVE_PASSWORD_FILE");
    let site = Site::start(cmd, &w, port);
    let got = request(port, "GET", "/", "");
    assert!(
        matches!(got, Ok((200, _))),
        "GET / of the default launcher: {:?}; stderr: {}",
        got.map(|(s, _)| s),
        site.stderr()
    );
    // The default logging is on: the request line reached stderr.
    assert!(
        site.stderr().contains("og live serve: GET /"),
        "{}",
        site.stderr()
    );
    drop(site);
    let _ = std::fs::remove_dir_all(&w);
}

/// tests/fixtures/families_live copied into `to`.
fn copy(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for e in std::fs::read_dir(from).unwrap() {
        let e = e.unwrap();
        let (src, dst) = (e.path(), to.join(e.file_name()));
        if e.file_type().unwrap().is_dir() {
            copy(&src, &dst);
        } else {
            std::fs::copy(&src, &dst).unwrap();
        }
    }
}

#[test]
fn an_unreadable_family_is_a_500_through_the_command_line() {
    let w = dir("family");
    let genome = w.join("genome");
    copy(&repo().join("tests/fixtures/families_live"), &genome);
    let rel = genome.join("family_a/relationships.tsv");
    std::fs::set_permissions(&rel, std::fs::Permissions::from_mode(0o000)).unwrap();
    if std::fs::read(&rel).is_ok() {
        eprintln!(
            "note: this account reads a mode 000 file (root), so the unreadable family is not checked"
        );
        let _ = std::fs::remove_dir_all(&w);
        return;
    }
    let (port, ingest) = (free_port(), free_port());
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_og"));
    cmd.args(["live", "serve", "--root"])
        .arg(&genome)
        .arg("--state-dir")
        .arg(w.join("state"))
        .args([
            "--port",
            &port.to_string(),
            "--ingest-port",
            &ingest.to_string(),
            "--quiet",
        ])
        .args(asset_flags());
    let site = Site::start(cmd, &w, port);
    let got = request(port, "GET", "/runs/family_a/family.json", "");
    let err = site.stderr();
    drop(site);
    std::fs::set_permissions(&rel, std::fs::Permissions::from_mode(0o644)).unwrap();
    let _ = std::fs::remove_dir_all(&w);
    assert_eq!(
        got,
        Ok((500, "family not readable\n".to_owned())),
        "family.json of an unreadable family; stderr: {err}"
    );
    assert!(err.contains("family.json of family_a not built"), "{err}");
}

#[test]
fn a_post_to_the_browser_port_never_prints_the_token() {
    let w = dir("wrong-port");
    let genome = w.join("genome");
    std::fs::create_dir_all(&genome).unwrap();
    let (port, ingest) = (free_port(), free_port());
    let state = w.join("state");
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_og"));
    cmd.args(["live", "serve", "--root"])
        .arg(&genome)
        .arg("--state-dir")
        .arg(&state)
        .args([
            "--port",
            &port.to_string(),
            "--ingest-port",
            &ingest.to_string(),
        ])
        .args(asset_flags());
    let site = Site::start(cmd, &w, port);
    let token = std::fs::read_to_string(state.join("ingest-token")).unwrap();
    let token = token.trim();
    assert_eq!(token.len(), 43, "the ingest token has its shape");
    let got = request(
        port,
        "POST",
        &format!("/{token}/ingest/run1"),
        r#"{"event":"started"}"#,
    );
    assert!(got.is_ok(), "the browser port answers the post: {got:?}");
    // Its request line is written before the answer; wait for it.
    let end = Instant::now() + WAIT;
    while !site.stderr().contains("/ingest/run1") && Instant::now() < end {
        std::thread::sleep(Duration::from_millis(50));
    }
    let err = site.stderr();
    drop(site);
    let _ = std::fs::remove_dir_all(&w);
    assert!(
        err.contains("og live serve: POST /<token>/ingest/run1"),
        "control: the request line is there, with the token written <token>: {err}"
    );
    assert!(!err.contains(token), "stderr shows the ingest token: {err}");
}
