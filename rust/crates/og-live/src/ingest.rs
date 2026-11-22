//! The ingest listener: the run's own processes post their progress here
//! (`POST /ingest/<run>`): a Nextflow weblog message, or a `pgp_step` /
//! `pgp_run_finished` message from `scripts/run-all.sh`. Each message is cut
//! to the fields a progress view needs ([`cut`]) and appended to
//! `STATE_DIR/runs/<run>/events.jsonl`, which the browser side reads. It
//! serves nothing.

use crate::runs::valid_run;
use crate::server::{REQUEST_SECONDS, read_head};
use crate::state::{Live, dumps_compact, lock};
use crate::viewer::unquote;
use indexmap::IndexMap;
use pycompat::json::{self, Value};
use std::fs::{self, OpenOptions};
use std::io::{self, BufReader, Read, Write};
use std::net::{Shutdown, TcpListener, TcpStream};
use std::sync::atomic::Ordering;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Duration;

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

/// One writer at a time for every events file.
static APPEND: Mutex<()> = Mutex::new(());

/// Appends message `ev` to `run`'s events file.
pub fn append(live: &Live, run: &str, ev: &Value) -> io::Result<()> {
    let p = live.store.path(run);
    let _g = lock(&APPEND);
    if let Some(d) = p.parent() {
        fs::create_dir_all(d)?;
    }
    let mut f = OpenOptions::new().create(true).append(true).open(p)?;
    f.write_all(format!("{}\n", dumps_compact(ev)).as_bytes())
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
    // A OWNGENOME_LIVE_URL ending in '/' gives '//ingest/...'.
    let path = format!("/{}", path.trim_start_matches('/'));
    let Some(run) = path.strip_prefix("/ingest/") else {
        return (404, "not found".into());
    };
    if header("origin").is_some() {
        // Nextflow and curl send no Origin.
        return (403, "a browser cannot post here".into());
    }
    let run = unquote(run);
    if !valid_run(&run) {
        return (400, "bad run name".into());
    }
    let n = match header("content-length").map(|v| v.trim().parse::<i64>()) {
        None | Some(Err(_)) => return (411, "Content-Length required".into()),
        Some(Ok(n)) => n,
    };
    if n < 0 || n > live.cfg.max_body_kb as i64 * 1024 {
        return (413, "too large".into());
    }
    let mut raw = Vec::with_capacity(n as usize);
    if body.take(n as u64).read_to_end(&mut raw).is_err() {
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

fn serve_conn(live: &Live, stream: TcpStream) {
    let timeout = Some(Duration::from_secs(REQUEST_SECONDS));
    let _ = stream.set_read_timeout(timeout);
    let _ = stream.set_write_timeout(timeout);
    let Ok(read_half) = stream.try_clone() else {
        return;
    };
    let mut reader = BufReader::new(read_half);
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
    thread::spawn(move || {
        for conn in listener.incoming() {
            if live.stopping.load(Ordering::SeqCst) {
                break;
            }
            let Ok(stream) = conn else { continue };
            let l = Arc::clone(&live);
            thread::spawn(move || serve_conn(&l, stream));
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn obj(text: &str) -> Value {
        json::loads(text).unwrap()
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
