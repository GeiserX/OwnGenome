//! One headless Chrome and one page, driven over the DevTools protocol.
//!
//! Chrome runs with `--remote-debugging-pipe`: it reads commands on fd 3
//! and writes answers and events on fd 4, each one JSON message ended by a
//! NUL byte. `/bin/sh` moves our two pipes to those numbers before it
//! starts Chrome, so no port is opened and no unsafe code is needed.

use indexmap::IndexMap;
use pycompat::json::{self, Dump, Value};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::{Child, ChildStdin, Command, Stdio};
use std::sync::mpsc::{self, Receiver, RecvTimeoutError};
use std::thread;
use std::time::{Duration, Instant};

/// How long one call may take before it is an error.
const CALL_TIMEOUT: Duration = Duration::from_secs(30);

/// A headless Chrome with one page attached. Dropping it closes Chrome.
pub struct Browser {
    child: Child,
    to_chrome: ChildStdin,
    from_chrome: Receiver<Value>,
    next_id: i128,
    session: String,
    errors: Vec<String>,
}

impl Browser {
    /// Start `chrome` with its profile in `profile` and attach to a blank
    /// page. A headless Chrome on a busy CI runner sometimes never answers
    /// its first call, so a fresh one (profile `<profile>-2`, `-3`) is
    /// started, up to three times.
    pub fn launch(chrome: &Path, profile: &Path) -> Result<Browser, String> {
        let mut last = String::new();
        for attempt in 1..=3 {
            let dir = if attempt == 1 {
                profile.to_path_buf()
            } else {
                let mut p = profile.as_os_str().to_owned();
                p.push(format!("-{attempt}"));
                p.into()
            };
            match Browser::start(chrome, &dir) {
                Ok(b) => return Ok(b),
                Err(e) => {
                    if attempt < 3 {
                        eprintln!(
                            "Chrome did not start ({e}); trying again ({} of 3)",
                            attempt + 1
                        );
                    }
                    last = e;
                }
            }
        }
        Err(last)
    }

    fn start(chrome: &Path, profile: &Path) -> Result<Browser, String> {
        let mut user_data = std::ffi::OsString::from("--user-data-dir=");
        user_data.push(profile);
        let mut child = Command::new("/bin/sh")
            .arg("-c")
            .arg(r#"exec "$0" "$@" 3<&0 4>&1 </dev/null >/dev/null"#)
            .arg(chrome)
            .args([
                "--headless=new",
                "--disable-gpu",
                "--no-sandbox",
                "--no-first-run",
                "--no-default-browser-check",
            ])
            .arg(user_data)
            .args(["--remote-debugging-pipe", "about:blank"])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .map_err(|e| format!("cannot start {}: {e}", chrome.display()))?;
        let to_chrome = child.stdin.take().expect("stdin is piped");
        let out = child.stdout.take().expect("stdout is piped");
        let (tx, from_chrome) = mpsc::channel();
        thread::spawn(move || {
            let mut r = BufReader::new(out);
            let mut buf = Vec::new();
            loop {
                buf.clear();
                match r.read_until(0, &mut buf) {
                    Ok(0) | Err(_) => break,
                    Ok(_) => {}
                }
                if buf.last() == Some(&0) {
                    buf.pop();
                }
                if let Ok(v) = json::loads(&String::from_utf8_lossy(&buf))
                    && tx.send(v).is_err()
                {
                    break;
                }
            }
        });
        let mut b = Browser {
            child,
            to_chrome,
            from_chrome,
            next_id: 0,
            session: String::new(),
            errors: Vec::new(),
        };
        let attached = (|| {
            let t = b.send(
                "Target.createTarget",
                obj([("url", s("about:blank"))]),
                None,
                CALL_TIMEOUT,
            )?;
            let target = field(&t, "targetId");
            let a = b.send(
                "Target.attachToTarget",
                obj([("targetId", target), ("flatten", Value::Bool(true))]),
                None,
                CALL_TIMEOUT,
            )?;
            b.session = match field(&a, "sessionId") {
                Value::Str(id) => id,
                v => return Err(format!("Target.attachToTarget: no sessionId in {v:?}")),
            };
            b.call("Runtime.enable", obj([]))?;
            b.call("Page.enable", obj([]))
        })();
        match attached {
            Ok(_) => Ok(b),
            Err(e) => {
                b.kill();
                Err(e)
            }
        }
    }

    /// Call `method` on the page with `params` (an object) and return its
    /// result, or the protocol's error.
    pub fn call(&mut self, method: &str, params: Value) -> Result<Value, String> {
        let session = self.session.clone();
        self.send(method, params, Some(&session), CALL_TIMEOUT)
    }

    fn send(
        &mut self,
        method: &str,
        params: Value,
        session: Option<&str>,
        timeout: Duration,
    ) -> Result<Value, String> {
        self.next_id += 1;
        let id = self.next_id;
        let mut msg = IndexMap::new();
        msg.insert("id".to_string(), Value::Int(id));
        msg.insert("method".to_string(), s(method));
        msg.insert("params".to_string(), params);
        if let Some(sid) = session {
            msg.insert("sessionId".to_string(), s(sid));
        }
        let mut line = json::dumps(&Value::Object(msg), &Dump::default())
            .map_err(|e| format!("{method}: {e:?}"))?;
        line.push('\0');
        self.to_chrome
            .write_all(line.as_bytes())
            .and_then(|()| self.to_chrome.flush())
            .map_err(|e| format!("{method}: Chrome closed the pipe ({e})"))?;
        let end = Instant::now() + timeout;
        loop {
            let left = end.saturating_duration_since(Instant::now());
            let m = match self.from_chrome.recv_timeout(left) {
                Ok(m) => m,
                Err(RecvTimeoutError::Timeout) => {
                    return Err(format!("{method}: no answer in {} s", timeout.as_secs()));
                }
                Err(RecvTimeoutError::Disconnected) => {
                    return Err(format!("{method}: Chrome closed the pipe"));
                }
            };
            if field(&m, "id") != Value::Int(id) {
                self.note(&m);
                continue;
            }
            if let Value::Object(e) = field(&m, "error") {
                return Err(format!(
                    "{method}: {}",
                    field(&Value::Object(e), "message").as_text()
                ));
            }
            return Ok(field(&m, "result"));
        }
    }

    /// Keep an event that reports a page error: an uncaught exception or a
    /// `console.error` call.
    fn note(&mut self, m: &Value) {
        let p = field(m, "params");
        match field(m, "method").as_text().as_str() {
            "Runtime.exceptionThrown" => {
                let d = field(&p, "exceptionDetails");
                let desc = field(&field(&d, "exception"), "description").as_text();
                self.errors.push(if desc.is_empty() {
                    field(&d, "text").as_text()
                } else {
                    desc
                });
            }
            "Runtime.consoleAPICalled" if field(&p, "type") == s("error") => {
                let args = match field(&p, "args") {
                    Value::Array(a) => a,
                    _ => Vec::new(),
                };
                let words: Vec<String> = args
                    .iter()
                    .map(|a| match field(a, "value") {
                        Value::Null => field(a, "description").as_text(),
                        v => v.as_text(),
                    })
                    .collect();
                self.errors
                    .push(format!("console.error: {}", words.join(" ")));
            }
            _ => {}
        }
    }

    /// Open `url` in the page. It returns once Chrome has started the
    /// navigation; [`Browser::wait`] for what the page must show.
    pub fn go(&mut self, url: &str) -> Result<(), String> {
        self.call("Page.navigate", obj([("url", s(url))])).map(drop)
    }

    /// The value of the JavaScript expression `expr` in the page (`Null`
    /// for `undefined`), or the exception it threw.
    pub fn eval(&mut self, expr: &str) -> Result<Value, String> {
        let r = self.call(
            "Runtime.evaluate",
            obj([
                ("expression", s(expr)),
                ("returnByValue", Value::Bool(true)),
            ]),
        )?;
        if let Value::Object(d) = field(&r, "exceptionDetails") {
            let d = Value::Object(d);
            let desc = field(&field(&d, "exception"), "description").as_text();
            return Err(if desc.is_empty() {
                field(&d, "text").as_text()
            } else {
                desc
            });
        }
        Ok(field(&field(&r, "result"), "value"))
    }

    /// Evaluate `expr` every 0.2 s until its value is truthy (as in
    /// JavaScript and Python), and return that value; `None` after
    /// `timeout`. An expression that throws counts as not yet.
    pub fn wait(&mut self, expr: &str, timeout: Duration) -> Option<Value> {
        let end = Instant::now() + timeout;
        loop {
            if let Ok(v) = self.eval(expr)
                && pycompat::value::truthy(&v)
            {
                return Some(v);
            }
            if Instant::now() >= end {
                return None;
            }
            thread::sleep(Duration::from_millis(200));
        }
    }

    /// The page's uncaught exceptions and `console.error` calls so far, in
    /// the order they came. Errors of [`Browser::eval`]'s own expressions
    /// are not among them.
    pub fn page_errors(&mut self) -> Vec<String> {
        while let Ok(m) = self.from_chrome.try_recv() {
            self.note(&m);
        }
        self.errors.clone()
    }

    fn kill(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

impl Drop for Browser {
    /// Ask Chrome to close, give it 10 s, then kill it.
    fn drop(&mut self) {
        let _ = self.send("Browser.close", obj([]), None, Duration::from_secs(5));
        let end = Instant::now() + Duration::from_secs(10);
        while Instant::now() < end {
            if matches!(self.child.try_wait(), Ok(Some(_))) {
                return;
            }
            thread::sleep(Duration::from_millis(100));
        }
        self.kill();
    }
}

/// `s` as a JavaScript string literal, to build an expression for
/// [`Browser::eval`] (Python's `json.dumps(s)`).
pub fn js_str(text: &str) -> String {
    let mut out = String::new();
    json::quote(&mut out, text, true);
    out
}

fn s(text: &str) -> Value {
    Value::Str(text.to_string())
}

fn obj<const N: usize>(pairs: [(&str, Value); N]) -> Value {
    Value::Object(pairs.into_iter().map(|(k, v)| (k.to_string(), v)).collect())
}

/// `v[key]`, `Null` when `v` is not an object or has no such key.
fn field(v: &Value, key: &str) -> Value {
    match v {
        Value::Object(m) => m.get(key).cloned().unwrap_or(Value::Null),
        _ => Value::Null,
    }
}

trait Text {
    fn as_text(&self) -> String;
}

impl Text for Value {
    /// A string as it is, anything else as JSON.
    fn as_text(&self) -> String {
        match self {
            Value::Str(t) => t.clone(),
            Value::Null => String::new(),
            v => json::dumps(v, &Dump::default()).unwrap_or_default(),
        }
    }
}
