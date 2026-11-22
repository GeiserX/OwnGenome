//! The ingest listener: the run's own processes post their progress here
//! (`POST /<token>/ingest/<run>`): a Nextflow weblog message, or a
//! `pgp_step` / `pgp_run_finished` message from `scripts/run-all.sh`. Each
//! message is cut to the fields a progress view needs ([`cut`]) and
//! appended to `STATE_DIR/runs/<run>/events.jsonl`, which the browser side
//! reads. It serves nothing.
//!
//! It listens on loopback only, and it takes a message only with the token
//! the server made ([`token`], `STATE_DIR/ingest-token`, mode 600), so
//! another account on the machine cannot post fake progress. Nextflow's
//! weblog cannot send a header, so the token is the first part of the
//! path: `OWNGENOME_LIVE_URL=http://127.0.0.1:8788/<token>`. The state
//! folder must be the server account's own ([`trusted_dir`]), and no file
//! in it is opened through a symlink.

use crate::auth::{new_token, same_secret};
use crate::runs::valid_run;
use crate::server::{ReadBy, accept_loop, read_head};
use crate::state::{Live, dumps_compact};
use crate::viewer::unquote;
use indexmap::IndexMap;
use pycompat::json::{self, Value};
use rustix::fs::{AtFlags, Mode, OFlags};
use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::net::{Shutdown, TcpListener, TcpStream};
use std::os::unix::fs::{FileExt, MetadataExt};
use std::path::Path;
use std::sync::Arc;
use std::thread;
use std::time::Instant;

/// The weblog's task fields kept. The work directory and the parameters
/// other than [`PARAM_KEYS`] are paths or personal: they never reach a page.
pub const TRACE_KEYS: [&str; 18] = [
    "task_id",
    "hash",
    "name",
    "process",
    "tag",
    "status",
    "exit",
    "attempt",
    "submit",
    "start",
    "complete",
    "realtime",
    "cpus",
    "memory",
    "peak_rss",
    "peak_vmem",
    "rchar",
    "wchar",
];
/// The run parameters kept.
pub const PARAM_KEYS: [&str; 3] = ["outdir", "max_memory", "max_cpus"];
/// The states of a `pgp_step` message.
pub const STEP_STATES: [&str; 3] = ["started", "ok", "failed"];

/// The first `n` characters of `s`.
fn cut_str(s: &str, n: usize) -> String {
    s.chars().take(n).collect()
}

/// `str(v)` of a scalar, as Python writes it.
fn py_str(v: &Value) -> String {
    match v {
        Value::Null => "None".into(),
        Value::Bool(b) => if *b { "True" } else { "False" }.into(),
        Value::Int(i) => i.to_string(),
        Value::Float(f) => pycompat::fmt::repr(*f),
        Value::Str(s) => s.clone(),
        other => dumps_compact(other),
    }
}

/// A scalar as kept: null, a flag and a number as they are, anything else
/// as text cut at `n` characters.
fn scalar(v: &Value, n: usize) -> Value {
    match v {
        Value::Null | Value::Bool(_) | Value::Int(_) | Value::Float(_) => v.clone(),
        other => Value::Str(cut_str(&py_str(other), n)),
    }
}

fn is_nested(v: &Value) -> bool {
    matches!(v, Value::Object(_) | Value::Array(_))
}

/// Whether `v` holds a number that is not finite (`NaN`, `Infinity`), which
/// would make events.jsonl something other than plain JSON.
fn has_non_finite(v: &Value) -> bool {
    match v {
        Value::Float(f) => !f.is_finite(),
        Value::Array(a) => a.iter().any(has_non_finite),
        Value::Object(o) => o.values().any(has_non_finite),
        _ => false,
    }
}

/// The fields of one message the site keeps, or why it is not a message.
pub fn cut(ev: &Value) -> Result<Value, String> {
    let Value::Object(ev) = ev else {
        return Err("not a JSON object".into());
    };
    let kind = match ev.get("event") {
        Some(Value::Str(k))
            if (1..=40).contains(&k.len())
                && k.bytes().all(|b| b.is_ascii_lowercase() || b == b'_') =>
        {
            k.clone()
        }
        _ => return Err("no 'event' name".into()),
    };
    let mut out: IndexMap<String, Value> = IndexMap::new();
    out.insert("event".into(), Value::Str(kind.clone()));
    for k in ["runName", "runId", "utcTime"] {
        if let Some(Value::Str(s)) = ev.get(k) {
            out.insert(k.into(), Value::Str(cut_str(s, 100)));
        }
    }
    if let Some(Value::Object(tr)) = ev.get("trace") {
        let mut t: IndexMap<String, Value> = IndexMap::new();
        for k in TRACE_KEYS {
            if let Some(v) = tr.get(k).filter(|v| !is_nested(v)) {
                t.insert(k.into(), scalar(v, 200));
            }
        }
        // Read as text later: a number there would not be.
        for k in ["name", "process", "tag", "status", "hash"] {
            if let Some(v) = t.get_mut(k).filter(|v| !matches!(v, Value::Null)) {
                *v = Value::Str(cut_str(&py_str(v), 200));
            }
        }
        out.insert("trace".into(), Value::Object(t));
    }
    if let Some(Value::Object(md)) = ev.get("metadata") {
        let mut keep: IndexMap<String, Value> = IndexMap::new();
        if let Some(Value::Object(params)) = md.get("parameters") {
            let mut p: IndexMap<String, Value> = IndexMap::new();
            for k in PARAM_KEYS {
                if let Some(v) = params.get(k).filter(|v| !is_nested(v)) {
                    p.insert(k.into(), scalar(v, 1000));
                }
            }
            if !p.is_empty() {
                keep.insert("parameters".into(), Value::Object(p));
            }
        }
        if let Some(Value::Object(wf)) = md.get("workflow") {
            if let Some(Value::Str(w)) = wf.get("workDir") {
                let mut o = IndexMap::new();
                o.insert("workDir".to_owned(), Value::Str(cut_str(w, 1000)));
                keep.insert("workflow".into(), Value::Object(o));
            }
            if let Some(Value::Str(m)) = wf.get("errorMessage") {
                out.insert("errorMessage".into(), Value::Str(cut_str(m, 500)));
            }
        }
        if !keep.is_empty() {
            out.insert("metadata".into(), Value::Object(keep));
        }
    }
    if let Some(Value::Str(m)) = ev.get("errorMessage") {
        out.insert("errorMessage".into(), Value::Str(cut_str(m, 500)));
    }
    if kind == "pgp_step" {
        let step = match ev.get("step") {
            Some(Value::Str(s)) if step_ok(s) => s.clone(),
            _ => return Err(STEP_NEEDS.into()),
        };
        let state = match ev.get("state") {
            Some(Value::Str(s)) if STEP_STATES.contains(&s.as_str()) => s.clone(),
            _ => return Err(STEP_NEEDS.into()),
        };
        out.insert("step".into(), Value::Str(step));
        out.insert("state".into(), Value::Str(state));
    } else if kind == "pgp_run_finished" {
        match ev.get("exit") {
            Some(Value::Int(n)) => {
                out.insert("exit".into(), Value::Int(*n));
            }
            _ => return Err("pgp_run_finished needs an integer 'exit'".into()),
        }
    }
    Ok(Value::Object(out))
}

const STEP_NEEDS: &str =
    "pgp_step needs 'step' (like 14 or 04b) and 'state' (started, ok or failed)";

/// Two digits and an optional lower-case letter (`14`, `04b`).
fn step_ok(s: &str) -> bool {
    let b = s.as_bytes();
    (b.len() == 2 || (b.len() == 3 && b[2].is_ascii_lowercase()))
        && b[0].is_ascii_digit()
        && b[1].is_ascii_digit()
}

/// The message a body holds, cut ([`cut`]), or why it is not one.
pub fn parse(body: &[u8]) -> Result<Value, String> {
    let text = std::str::from_utf8(body).map_err(|e| e.to_string())?;
    let v = json::loads(text).map_err(|e| format!("{e:?}"))?;
    if has_non_finite(&v) {
        return Err("NaN and Infinity are not JSON".into());
    }
    cut(&v)
}

/// The file of the ingest token, in the state folder.
pub const TOKEN_FILE: &str = "ingest-token";

/// Whether `t` has the shape of a token the server makes: 43 URL-safe
/// characters.
fn token_shape(t: &str) -> bool {
    t.len() == 43
        && t.bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'-' || b == b'_')
}

/// The account the server runs as.
fn me() -> u32 {
    rustix::process::geteuid().as_raw()
}

/// Opens `name` in folder `dir`, never through a symlink: a link there is
/// an error (`ELOOP`), not followed.
fn open_at(dir: &File, name: &str, flags: OFlags, mode: Mode) -> io::Result<File> {
    rustix::fs::openat(dir, name, flags | OFlags::NOFOLLOW | OFlags::CLOEXEC, mode)
        .map(File::from)
        .map_err(io::Error::from)
}

/// Whether the open file `f` is the server account's own: owned by it, and
/// no other account may write to it.
fn own(f: &File) -> io::Result<bool> {
    let m = f.metadata()?;
    Ok(m.uid() == me() && m.mode() & 0o022 == 0)
}

/// The state folder `dir`, opened, when it is the server account's own: a
/// folder owned by the account the server runs as that no other account
/// can write to, so no one else can put a file or a link in it. The checks
/// are on the open folder, so swapping the path afterwards changes nothing.
pub fn trusted_dir(dir: &Path) -> Result<File, String> {
    let f = rustix::fs::open(
        dir,
        OFlags::RDONLY | OFlags::DIRECTORY | OFlags::CLOEXEC,
        Mode::empty(),
    )
    .map(File::from)
    .map_err(|e| format!("{}: {}", dir.display(), io::Error::from(e)))?;
    if own(&f).map_err(|e| format!("{}: {e}", dir.display()))? {
        Ok(f)
    } else {
        Err(format!(
            "the state folder {} must be owned by the account the server runs as, and no other \
             account may write to it, because it holds the ingest token: chown it to this account \
             and run 'chmod go-w {}'",
            dir.display(),
            dir.display()
        ))
    }
}

/// The folder `name` in `parent`, made (mode 755) when missing, never
/// through a symlink.
fn sub_dir(parent: &File, name: &str) -> io::Result<File> {
    match rustix::fs::mkdirat(
        parent,
        name,
        Mode::RWXU | Mode::RGRP | Mode::XGRP | Mode::ROTH | Mode::XOTH,
    ) {
        Ok(()) | Err(rustix::io::Errno::EXIST) => {}
        Err(e) => return Err(e.into()),
    }
    open_at(
        parent,
        name,
        OFlags::RDONLY | OFlags::DIRECTORY,
        Mode::empty(),
    )
}

/// The token the ingest port takes. The one in `STATE_DIR/ingest-token` is
/// kept when it is a plain file the server account owns and only it can
/// read, and has the shape of a token the server makes, so a restart does
/// not cut off a run that is posting. Anything else there (a link, a file
/// others can read or another account owns, a value someone chose) is
/// replaced with a new random token. A state folder that is not the server
/// account's own ([`trusted_dir`]) is an error: the port does not start.
pub fn token(state_dir: &Path) -> Result<String, String> {
    let dir = trusted_dir(state_dir)?;
    let err = |name: &str, e: io::Error| format!("{}: {e}", state_dir.join(name).display());
    // Two servers started at once on one state folder must end with the
    // same token: one reads or writes it at a time. The lock is on the
    // folder itself and ends with its handle. (A lock file would race:
    // macOS can fail an `openat` that creates a name another thread
    // creates at the same moment, with ENOENT.)
    dir.lock().map_err(|e| err("", e))?;
    if let Ok(mut f) = open_at(&dir, TOKEN_FILE, OFlags::RDONLY, Mode::empty())
        && let Ok(meta) = f.metadata()
        && meta.file_type().is_file()
        && meta.uid() == me()
        && meta.mode() & 0o077 == 0
    {
        let mut text = String::new();
        if f.read_to_string(&mut text).is_ok() && token_shape(text.trim()) {
            return Ok(text.trim().to_owned());
        }
    }
    let t = new_token().map_err(|e| err(TOKEN_FILE, e))?;
    let tmp = format!("{TOKEN_FILE}.{}", std::process::id());
    let _ = rustix::fs::unlinkat(&dir, tmp.as_str(), AtFlags::empty());
    let mut f = open_at(
        &dir,
        &tmp,
        OFlags::WRONLY | OFlags::CREATE | OFlags::EXCL,
        Mode::RUSR | Mode::WUSR,
    )
    .map_err(|e| err(&tmp, e))?;
    f.write_all(format!("{t}\n").as_bytes())
        .map_err(|e| err(&tmp, e))?;
    drop(f);
    rustix::fs::renameat(&dir, tmp.as_str(), &dir, TOKEN_FILE)
        .map_err(|e| err(TOKEN_FILE, e.into()))?;
    drop(dir);
    Ok(t)
}

const NEEDS_TOKEN: &str = "an ingest token is needed: post to /<token>/ingest/<run>, with the token of \
                           STATE_DIR/ingest-token";

/// `run`'s events file (`STATE_DIR/runs/<run>/events.jsonl`, the file
/// [`crate::state::Store`] reads), opened to append, and its folder, locked
/// (an exclusive lock that ends with the folder's handle). The state folder
/// must be the server account's own, and neither the folders below it nor
/// the file are followed through a symlink. The lock is taken before the
/// file is opened, so two writers never create it at the same moment
/// (which macOS can fail with ENOENT).
fn open_events(state_dir: &Path, run: &str) -> io::Result<(File, File)> {
    let dir = trusted_dir(state_dir).map_err(io::Error::other)?;
    let runs = sub_dir(&dir, "runs")?;
    let folder = sub_dir(&runs, run)?;
    folder.lock()?;
    let f = open_at(
        &folder,
        "events.jsonl",
        OFlags::RDWR | OFlags::APPEND | OFlags::CREATE,
        Mode::RUSR | Mode::WUSR | Mode::RGRP | Mode::ROTH,
    )?;
    Ok((f, folder))
}

/// Appends message `ev` to `run`'s events file, as a line of its own: a
/// write that fails part way is cut off again, and a file that still ends
/// in the middle of a line (the cut failed too) gets a line break first, so
/// the broken line is skipped on reading and never joins this message.
pub fn append(live: &Live, run: &str, ev: &Value) -> io::Result<()> {
    append_with(&live.cfg.state_dir, run, ev, |f, line| f.write_all(line))
}

/// [`append`], with `write` writing the line. The run's folder stays
/// locked (an exclusive lock every writer takes, in this process or in
/// another `og live serve` on the same state folder) from reading the
/// file's length to the cut after a failed write, so the cut never removes
/// another writer's message.
fn append_with(
    state_dir: &Path,
    run: &str,
    ev: &Value,
    write: impl FnOnce(&mut File, &[u8]) -> io::Result<()>,
) -> io::Result<()> {
    let (mut f, _locked) = open_events(state_dir, run)?;
    let start = f.metadata()?.len();
    let mut line = format!("{}\n", dumps_compact(ev));
    if start > 0 {
        let mut last = [0u8];
        f.read_exact_at(&mut last, start - 1)?;
        if last[0] != b'\n' {
            line.insert(0, '\n');
        }
    }
    write(&mut f, line.as_bytes()).inspect_err(|_| {
        let _ = f.set_len(start);
    })
}

/// The answer: `code` and one line of text.
fn reply(out: &mut dyn Write, code: u16, text: &str) -> io::Result<()> {
    let reason = match code {
        200 => "OK",
        400 => "Bad Request",
        403 => "Forbidden",
        404 => "Not Found",
        411 => "Length Required",
        413 => "Request Entity Too Large",
        _ => "Not Implemented",
    };
    let body = format!("{text}\n");
    let date = jiff::Timestamp::now().strftime("%a, %d %b %Y %H:%M:%S GMT");
    write!(
        out,
        "HTTP/1.1 {code} {reason}\r\nServer: owngenome-live-ingest\r\nDate: {date}\r\n\
         Content-Type: text/plain; charset=utf-8\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{body}",
        body.len()
    )?;
    out.flush()
}

/// The status and text for one request.
fn answer(
    live: &Live,
    method: &str,
    target: &str,
    headers: &[(String, String)],
    body: &mut dyn Read,
) -> (u16, String) {
    let header = |n: &str| {
        headers
            .iter()
            .find(|(k, _)| k == n)
            .map(|(_, v)| v.as_str())
    };
    if method != "POST" {
        return if method == "GET" {
            (404, "not found".into())
        } else {
            (501, "unsupported method".into())
        };
    }
    let path = target.split(['?', '#']).next().unwrap_or("");
    // /<token>/ingest/<run>; a OWNGENOME_LIVE_URL ending in '/' gives
    // '//' before ingest.
    let rest = path.trim_start_matches('/');
    let (given, after) = rest.split_once('/').unwrap_or((rest, ""));
    let run = after.trim_start_matches('/').strip_prefix("ingest/");
    if run.is_none() && given != "ingest" {
        return (404, "not found".into());
    }
    if header("origin").is_some() {
        // Nextflow and curl send no Origin.
        return (403, "a browser cannot post here".into());
    }
    let taken = live
        .cfg
        .ingest_token
        .as_deref()
        .is_some_and(|t| same_secret(t.as_bytes(), given.as_bytes()));
    let Some(run) = run.filter(|_| taken) else {
        return (403, NEEDS_TOKEN.into());
    };
    let run = unquote(run);
    if !valid_run(&run) {
        return (400, "bad run name".into());
    }
    let n = match header("content-length").map(|v| v.trim().parse::<i64>()) {
        None | Some(Err(_)) => return (411, "Content-Length required".into()),
        Some(Ok(n)) => n,
    };
    let Ok(n) = u64::try_from(n) else {
        return (413, "too large".into());
    };
    if n > max_body_bytes(live.cfg.max_body_kb) {
        return (413, "too large".into());
    }
    let mut raw = Vec::new();
    if body.take(n).read_to_end(&mut raw).is_err() || raw.len() as u64 != n {
        return (400, "not a message: the body was cut off".into());
    }
    let ev = match parse(&raw) {
        Ok(ev) => ev,
        Err(e) => return (400, cut_str(&format!("not a message: {e}"), 200)),
    };
    match append(live, &run, &ev) {
        Ok(()) => (200, "ok".into()),
        Err(e) => {
            eprintln!(
                "og live serve: ingest: cannot store a message of {run}: {}",
                e.kind()
            );
            (500, "not stored".into())
        }
    }
}

/// The largest body, in bytes, of a `--max-body-kb` of `kb` (validated by
/// the command line; a value set in code that would overflow means no
/// limit below the largest length).
pub fn max_body_bytes(kb: usize) -> u64 {
    u64::try_from(kb).map_or(u64::MAX, |kb| kb.saturating_mul(1024))
}

fn serve_conn(live: &Live, stream: TcpStream) {
    let timeout = live.cfg.request_timeout;
    let _ = stream.set_write_timeout(Some(timeout));
    let Ok(read_half) = stream.try_clone() else {
        return;
    };
    // The head and the body share one deadline.
    let mut reader = BufReader::new(ReadBy {
        tcp: read_half,
        end: Instant::now() + timeout,
    });
    let mut out = io::BufWriter::new(&stream);
    let (code, text) = match read_head(&mut reader) {
        None => (400, "bad request".to_owned()),
        Some((method, target, headers)) => answer(live, &method, &target, &headers, &mut reader),
    };
    let _ = reply(&mut out, code, &text);
    drop(out);
    let _ = stream.shutdown(Shutdown::Write);
}

/// Serves `listener` on threads until `live` is stopping.
pub fn serve(live: Arc<Live>, listener: TcpListener) -> thread::JoinHandle<()> {
    accept_loop(live, listener, serve_conn)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::os::unix::fs::PermissionsExt;
    use std::process::Command;
    use std::time::Duration;

    /// Set by [`a_failed_append_never_cuts_another_servers_message`] for
    /// the copy of this test binary it starts: the state folder the second
    /// server appends to.
    const OTHER_WRITER: &str = "OG_LIVE_TEST_OTHER_WRITER";

    /// The second server: another process that appends one message to the
    /// run of the state folder in [`OTHER_WRITER`]. Run only by the test
    /// below, as a process of its own.
    #[test]
    #[ignore = "run by a_failed_append_never_cuts_another_servers_message"]
    fn other_writer() {
        let Ok(dir) = std::env::var(OTHER_WRITER) else {
            return;
        };
        let dir = Path::new(&dir);
        fs::write(dir.join("other-is-writing"), "").unwrap();
        append_with(dir, "R1", &obj(r#"{"event":"b"}"#), |f, l| f.write_all(l)).unwrap();
    }

    #[test]
    fn a_failed_append_never_cuts_another_servers_message() {
        let dir = std::env::temp_dir().join(format!("og-live-two-writers-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        fs::set_permissions(&dir, fs::Permissions::from_mode(0o755)).unwrap();
        let mut other = None;
        let got = append_with(&dir, "R1", &obj(r#"{"event":"a"}"#), |f, line| {
            // While this write is under way, another server on the same
            // state folder posts a message of the same run.
            let child = Command::new(std::env::current_exe().unwrap())
                .args(["ingest::tests::other_writer", "--exact", "--ignored"])
                .env(OTHER_WRITER, &dir)
                .stdout(std::process::Stdio::null())
                .spawn()
                .unwrap();
            let began = Instant::now();
            while !dir.join("other-is-writing").exists()
                && began.elapsed() < Duration::from_secs(20)
            {
                thread::sleep(Duration::from_millis(20));
            }
            thread::sleep(Duration::from_millis(500));
            other = Some(child);
            // Then this write fails part way, and its part is cut off.
            f.write_all(&line[..4])?;
            Err(io::Error::other("the disk is full"))
        });
        assert!(got.is_err());
        let status = other.unwrap().wait().unwrap();
        assert!(status.success(), "the second server's append failed");
        let text = fs::read_to_string(dir.join("runs/R1/events.jsonl")).unwrap();
        assert_eq!(
            text, "{\"event\":\"b\"}\n",
            "the cut removed the other server's message"
        );
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn messages_posted_at_once_to_a_new_run_are_all_stored() {
        let dir = std::env::temp_dir().join(format!("og-live-first-posts-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        fs::set_permissions(&dir, fs::Permissions::from_mode(0o755)).unwrap();
        for round in 0..20 {
            let run = format!("R{round}");
            let go = Arc::new(std::sync::Barrier::new(8));
            let done: Vec<io::Result<()>> = (0..8)
                .map(|_| {
                    let (dir, run, go) = (dir.clone(), run.clone(), Arc::clone(&go));
                    thread::spawn(move || {
                        go.wait();
                        append_with(&dir, &run, &obj(r#"{"event":"a"}"#), |f, l| f.write_all(l))
                    })
                })
                .collect::<Vec<_>>()
                .into_iter()
                .map(|h| h.join().unwrap())
                .collect();
            assert!(done.iter().all(Result::is_ok), "{run}: {done:?}");
            let text =
                fs::read_to_string(dir.join("runs").join(&run).join("events.jsonl")).unwrap();
            assert_eq!(text.lines().count(), 8, "{run}: {text:?}");
        }
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn a_state_folder_another_account_can_write_to_is_refused() {
        let dir = std::env::temp_dir().join(format!("og-live-open-state-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        for mode in [0o775, 0o757] {
            fs::set_permissions(&dir, fs::Permissions::from_mode(mode)).unwrap();
            let e = token(&dir).expect_err("a state folder others can write to");
            assert!(e.contains("chmod go-w"), "{e}");
            assert!(!dir.join(TOKEN_FILE).exists());
            let e = append_with(&dir, "R1", &obj(r#"{"event":"a"}"#), |f, l| f.write_all(l));
            assert!(
                e.is_err(),
                "a message was stored in a state folder others can write to"
            );
        }
        fs::set_permissions(&dir, fs::Permissions::from_mode(0o755)).unwrap();
        assert!(token(&dir).is_ok());
        fs::remove_dir_all(&dir).unwrap();
    }

    fn obj(text: &str) -> Value {
        json::loads(text).unwrap()
    }

    #[test]
    fn servers_started_together_on_one_state_folder_get_one_token() {
        let dir = std::env::temp_dir().join(format!("og-live-token-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let go = Arc::new(std::sync::Barrier::new(8));
        let all: Vec<String> = (0..8)
            .map(|_| {
                let (dir, go) = (dir.clone(), Arc::clone(&go));
                thread::spawn(move || {
                    go.wait();
                    token(&dir).unwrap()
                })
            })
            .collect::<Vec<_>>()
            .into_iter()
            .map(|h| h.join().unwrap())
            .collect();
        let on_disk = fs::read_to_string(dir.join(TOKEN_FILE)).unwrap();
        assert!(
            all.iter().all(|t| *t == on_disk.trim()),
            "{all:?} {on_disk}"
        );
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn a_message_keeps_only_what_a_progress_view_needs() {
        let ev = obj(
            r#"{"event":"process_started","runName":"r","trace":{"name":123,"process":"P","tag":null,
            "status":true,"workdir":"/w","script":"bcftools roh","exit":0,"env":{"A":1}},
            "metadata":{"parameters":{"outdir":"/data/genome","hpo":"HP:0001250","vcf":"x.vcf.gz"},
            "workflow":{"workDir":"/w","commandLine":"nextflow run","errorMessage":"boom"}}}"#,
        );
        let got = dumps_compact(&cut(&ev).unwrap());
        assert_eq!(
            got,
            r#"{"event":"process_started","runName":"r","trace":{"name":"123","process":"P","tag":null,"status":"True","exit":0},"errorMessage":"boom","metadata":{"parameters":{"outdir":"/data/genome"},"workflow":{"workDir":"/w"}}}"#
        );
    }

    #[test]
    fn script_messages_need_their_fields() {
        let ok = cut(&obj(r#"{"event":"pgp_step","step":"04b","state":"ok"}"#)).unwrap();
        assert_eq!(
            dumps_compact(&ok),
            r#"{"event":"pgp_step","step":"04b","state":"ok"}"#
        );
        for bad in [
            r#"{"event":"pgp_step","step":"4","state":"ok"}"#,
            r#"{"event":"pgp_step","step":"04","state":"done"}"#,
            r#"{"event":"pgp_run_finished","exit":true}"#,
            r#"{"event":"pgp_run_finished","exit":"0"}"#,
            r#"{"event":"Bad"}"#,
            r#"[1]"#,
        ] {
            assert!(cut(&obj(bad)).is_err(), "{bad}");
        }
        assert!(cut(&obj(r#"{"event":"pgp_run_finished","exit":3}"#)).is_ok());
    }

    #[test]
    fn nan_infinity_and_text_that_is_not_json_are_refused() {
        assert!(parse(br#"{"event":"started","trace":{"attempt":NaN}}"#).is_err());
        assert!(parse(br#"{"event":"started","trace":{"exit":Infinity}}"#).is_err());
        assert!(parse(br#"{"event":"started","x":[1e400]}"#).is_err());
        assert!(parse(b"not json").is_err());
        assert!(parse(b"\xff").is_err());
        assert!(parse(br#"{"event":"started"}"#).is_ok());
    }
}
