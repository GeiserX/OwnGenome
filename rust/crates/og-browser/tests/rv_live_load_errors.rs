//! A live page that cannot load its data says so, in a headless Chrome. A small demo site built from
//! assets/live (the home page, DEMO-001, the child DEMO-003 and the family DEMO-FAMILY) is served by a stub
//! that can answer one path with an error, with no answer at all, with a given body or after a delay, with
//! the live site's Content-Security-Policy.
//!
//! - the person page whose results.json answers 500 says the results could not be loaded, with Try again;
//!   once the stub recovers, Try again draws them. An OK answer draws as before.
//! - the home page whose runs.json answers 500, or no answer, says the list of people could not be loaded,
//!   with Try again, and keeps the server's list of runs; Try again draws the page.
//! - the family page whose family.json answers `{}` says the family could not be loaded; Try again draws it.
//! - DEMO-001 with frames.json and results.json at 500: the closed page says "Results not loaded" and never
//!   says child. The control: DEMO-003, a real child with good files, still reads "Results about a child".
//! - frames.json held 15 s: the wire, run health and re-check panels say they wait for the run's data,
//!   never that there is nothing; once the frames arrive the replay shows its steps and the dated re-check.
//!
//! Each case was red on main before the fix (the PR says where). Without Chrome the test says SKIP, and
//! fails in CI.

use og_browser::{Browser, Scratch, Value, chrome_or_skip, js_str, repo_root};
use pycompat::json::{self, Dump};
use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::net::{Shutdown, TcpListener, TcpStream};
use std::path::{Component, Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};

/// The live site's policy (`CSP` of bin/owngenome_live.py; live_replay.rs keeps its copy equal to it).
const CSP: &str = "default-src 'none'; script-src 'self'; connect-src 'self'; style-src 'self'; \
     img-src 'self' data:; form-action 'self'; base-uri 'none'; frame-src 'self'; frame-ancestors 'none'";
const UNREAD: &str = "This page could not read this run's data, so it shows nothing yet.";
const WAITING: &str = "Waiting for the run's data.";
const BANNER: &str = "Database re-check of 2026-01-02T11:05:40Z: ClinVar 2026-01 is newer than the one this run used.";

/// How the stub answers one path.
#[derive(Clone)]
enum Answer {
    /// The file under the site's folder (404 when there is none).
    File,
    /// This status, with a plain-text body.
    Status(u16),
    /// No answer: the connection closes, so fetch rejects (the browser is offline for that file).
    Abort,
    /// 200 with this JSON body.
    Json(&'static str),
    /// The file, after this delay.
    Hold(Duration),
}

/// A site on 127.0.0.1 whose answer per path the test sets with [`Stub::set`]. Threads still holding an
/// answer when it is dropped end on their own.
struct Stub {
    port: u16,
    rules: Arc<Mutex<HashMap<String, Answer>>>,
    stop: Arc<AtomicBool>,
}

impl Stub {
    fn start(dir: &Path) -> Stub {
        let root = Arc::new(dir.canonicalize().unwrap());
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let port = listener.local_addr().unwrap().port();
        let rules: Arc<Mutex<HashMap<String, Answer>>> = Arc::default();
        let stop = Arc::new(AtomicBool::new(false));
        let (r, s) = (Arc::clone(&rules), Arc::clone(&stop));
        thread::spawn(move || {
            for conn in listener.incoming() {
                if s.load(Ordering::SeqCst) {
                    break;
                }
                if let Ok(c) = conn {
                    let (root, rules) = (Arc::clone(&root), Arc::clone(&r));
                    thread::spawn(move || answer(c, &root, &rules));
                }
            }
        });
        Stub { port, rules, stop }
    }

    /// Answer `path` (from the root, like "/runs.json") this way from now on.
    fn set(&self, path: &str, a: Answer) {
        self.rules.lock().unwrap().insert(path.to_string(), a);
    }

    fn url(&self, path: &str) -> String {
        format!("http://127.0.0.1:{}{path}", self.port)
    }
}

impl Drop for Stub {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::SeqCst);
        let _ = TcpStream::connect(("127.0.0.1", self.port));
    }
}

fn answer(mut c: TcpStream, root: &Path, rules: &Mutex<HashMap<String, Answer>>) {
    let _ = c.set_read_timeout(Some(Duration::from_secs(10)));
    let mut r = BufReader::new(&c);
    let mut first = String::new();
    if r.read_line(&mut first).is_err() {
        return;
    }
    let mut line = String::new();
    while r.read_line(&mut line).is_ok_and(|n| n > 2) {
        line.clear();
    }
    let target = first.split_whitespace().nth(1).unwrap_or("/");
    let path = target.split(['?', '#']).next().unwrap_or("/").to_string();
    let rule = rules
        .lock()
        .unwrap()
        .get(&path)
        .cloned()
        .unwrap_or(Answer::File);
    let file = |p: &str| -> Option<(Vec<u8>, &'static str)> {
        let mut rel = p.trim_start_matches('/').to_string();
        if rel.is_empty() || rel.ends_with('/') {
            rel.push_str("index.html");
        }
        let rel = PathBuf::from(rel);
        if !rel.components().all(|c| matches!(c, Component::Normal(_))) {
            return None;
        }
        let ctype = match rel.extension().and_then(|e| e.to_str()) {
            Some("html") => "text/html; charset=utf-8",
            Some("js") => "text/javascript; charset=utf-8",
            Some("css") => "text/css; charset=utf-8",
            Some("json") => "application/json",
            _ => "application/octet-stream",
        };
        fs::read(root.join(rel)).ok().map(|b| (b, ctype))
    };
    let (status, body, ctype) = match rule {
        Answer::Abort => {
            let _ = c.shutdown(Shutdown::Both);
            return;
        }
        Answer::Status(n) => (format!("{n} Stub"), b"stub error\n".to_vec(), "text/plain"),
        Answer::Json(j) => ("200 OK".into(), j.as_bytes().to_vec(), "application/json"),
        Answer::Hold(d) => {
            thread::sleep(d);
            served(file(&path))
        }
        Answer::File => served(file(&path)),
    };
    let head = format!(
        "HTTP/1.1 {status}\r\nContent-Type: {ctype}\r\nContent-Length: {}\r\nCache-Control: no-store\r\n\
         Content-Security-Policy: {CSP}\r\nConnection: close\r\n\r\n",
        body.len()
    );
    let _ = c
        .write_all(head.as_bytes())
        .and_then(|()| c.write_all(&body));
    let _ = c.shutdown(Shutdown::Write);
}

fn served(f: Option<(Vec<u8>, &'static str)>) -> (String, Vec<u8>, &'static str) {
    match f {
        Some((b, t)) => ("200 OK".into(), b, t),
        None => (
            "404 Not Found".into(),
            b"not found\n".to_vec(),
            "text/plain",
        ),
    }
}

/// The frames of tests/fixtures/live/frames.json with `plugins` set on every frame's state.
fn frames_with(plugins: &str) -> String {
    let text = fs::read_to_string(repo_root().join("tests/fixtures/live/frames.json")).unwrap();
    let mut fr = json::loads(&text).unwrap();
    let p = json::loads(plugins).unwrap();
    if let Value::Object(m) = &mut fr
        && let Some(Value::Array(a)) = m.get_mut("frames")
    {
        for f in a.iter_mut() {
            if let Value::Object(f) = f
                && let Some(Value::Object(s)) = f.get_mut("state")
            {
                s.insert("plugins".into(), p.clone());
            }
        }
    }
    json::dumps(&fr, &Dump::default()).unwrap()
}

/// A person page of the demo, as tests/demo/build_live_demo.py run_page fills it, with the ops and
/// re-check plugins.
fn run_page(dir: &Path, name: &str, template: &str, frames: &str, results: Option<&str>) {
    let tags = ["ops.js", "recheck.js"]
        .map(|f| format!("<script src=\"../../assets/plugins/{f}\" defer></script>"))
        .join("\n");
    let page = fs::read_to_string(repo_root().join("assets/live").join(template))
        .unwrap()
        .replace("{{assets}}", "../../assets")
        .replace("{{run}}", name)
        .replace("{{mode}}", "replay")
        .replace(
            "{{plugins}}",
            if template == "index.html" { &tags } else { "" },
        )
        .replace("href=\"/\"", "href=\"../../\"");
    let d = dir.join("runs").join(name);
    fs::create_dir_all(&d).unwrap();
    fs::write(d.join("index.html"), page).unwrap();
    fs::write(d.join("frames.json"), frames).unwrap();
    if let Some(r) = results {
        fs::write(d.join("results.json"), r).unwrap();
    }
}

/// The demo site: assets/, the home page with runs.json, DEMO-001 (an adult), DEMO-003 (a child with a
/// result to act on) and the family page of DEMO-FAMILY.
fn demo_site(dir: &Path) {
    let assets = repo_root().join("assets/live");
    fs::create_dir_all(dir.join("assets/plugins")).unwrap();
    for f in [
        "core.js",
        "live.js",
        "map.js",
        "person.js",
        "home.js",
        "family-page.js",
        "live.css",
        "map.css",
    ] {
        fs::copy(assets.join(f), dir.join("assets").join(f)).unwrap();
    }
    for f in ["ops.js", "recheck.js", "family.js"] {
        fs::copy(
            assets.join("plugins").join(f),
            dir.join("assets/plugins").join(f),
        )
        .unwrap();
    }
    let home = fs::read_to_string(assets.join("home.html"))
        .unwrap()
        .replace("{{assets}}", "assets")
        .replace("{{data}}", "runs.json")
        .replace(
            "{{runs}}",
            "<li><a href=\"runs/DEMO-001/\">DEMO-001</a></li>",
        );
    fs::write(dir.join("index.html"), home).unwrap();
    fs::write(
        dir.join("runs.json"),
        r#"{"runs": [{"name": "DEMO-001", "kind": "person", "url": "runs/DEMO-001/", "state": "finished",
            "started_utc": "2026-01-02T09:00:00Z", "finished_utc": "2026-01-02T10:05:40Z", "planned": 2, "done": 2,
            "running": 0, "failed": 0, "top": []}], "families": []}"#,
    )
    .unwrap();
    let recheck = format!(
        r#"{{"recheck": {{"file": "recheck/state.json", "mtime": 1, "checked_utc": "2026-01-02T11:05:40Z",
            "banner": {}, "steps": {{"03": ["newer ClinVar available"]}}}}}}"#,
        js_str(BANNER)
    );
    let results = r#"{"schema_version": 1, "sample": "DEMO-001", "person": {"sample": "DEMO-001"}, "chapters": [],
        "findings": [], "levels": {}, "note": null, "lead": [],
        "journal": [{"n": 8, "key": "assessed", "title": "Assessed, no finding", "question": "What did the steps look at?",
            "cards": [], "findings": []}]}"#;
    run_page(
        dir,
        "DEMO-001",
        "index.html",
        &frames_with(&recheck),
        Some(results),
    );
    let child = r#"{"family": {"kind": "member", "child": true, "families": [{"name": "DEMO-FAMILY", "url": "../DEMO-FAMILY/"}]}}"#;
    let acts = r#"{"schema_version": 1, "sample": "DEMO-003", "chapters": [], "note": null,
        "findings": [{"id": "act-1", "level": "act", "card": "clinical", "title": "A result to act on"}]}"#;
    run_page(
        dir,
        "DEMO-003",
        "index.html",
        &frames_with(child),
        Some(acts),
    );
    run_page(dir, "DEMO-FAMILY", "family.html", &frames_with("{}"), None);
    fs::write(
        dir.join("runs/DEMO-FAMILY/family.json"),
        r#"{"state": {"kind": "family", "name": "DEMO-FAMILY", "members": [], "children": [], "couples": [],
            "plan": [], "reports": []}, "runs": {}}"#,
    )
    .unwrap();
}

fn text(v: Option<Value>) -> String {
    match v {
        Some(Value::Str(s)) => s,
        Some(v) => format!("{v:?}"),
        None => "None".into(),
    }
}

#[test]
fn a_page_that_cannot_load_its_data_says_so() {
    let Some(chrome) = chrome_or_skip("load the live demo pages") else {
        return;
    };
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "load-errors");
    let site = scratch.path().join("site");
    demo_site(&site);
    let stub = Stub::start(&site);
    let mut b = Browser::launch(&chrome, &scratch.path().join("chrome")).unwrap();
    let mut fails = Vec::new();
    let mut check = |desc: &str, ok: bool, detail: String| {
        println!(
            "[{}] {desc}{}",
            if ok { "PASS" } else { "FAIL" },
            if ok {
                String::new()
            } else {
                format!(" -- {detail}")
            }
        );
        if !ok {
            fails.push(desc.to_string());
        }
    };
    let secs = Duration::from_secs;
    let failed_box = |id: &str| {
        format!(
            "(function () {{ var e = document.getElementById({}); return e ? e.innerText.replace(/\\s+/g, ' ').trim() : ''; }}())",
            js_str(id)
        )
    };

    // The person page: results.json at 500, then an OK answer through Try again.
    stub.set("/runs/DEMO-001/results.json", Answer::Status(500));
    b.go(&stub.url("/runs/DEMO-001/")).unwrap();
    let said = text(b.wait(&failed_box("results-failed"), secs(20)));
    let status = text(
        b.eval("document.documentElement.dataset.pgpStatus || ''")
            .ok(),
    );
    check(
        "person: results.json at 500 says 'The results could not be loaded.' with Try again, and the page is offline",
        said == "The results could not be loaded. Try again" && status == "offline",
        format!("{said:?} {status:?}"),
    );
    stub.set("/runs/DEMO-001/results.json", Answer::File);
    b.eval("var t = document.getElementById('results-retry'); if (t) { t.click(); } !!t")
        .ok();
    let drawn = text(b.wait(
        "!document.getElementById('results-failed') && document.querySelector('#results section.sec h2') && \
         document.querySelector('#results section.sec h2').textContent",
        secs(10),
    ));
    check(
        "person: once results.json answers, Try again draws the results",
        drawn == "Assessed, no finding",
        drawn,
    );
    // An OK answer from the start draws as before.
    b.go(&stub.url("/runs/DEMO-001/?ok")).unwrap();
    let ok = text(b.wait(
        "document.querySelector('#results section.sec h2') && [document.querySelector('#results section.sec h2').textContent, \
         !!document.getElementById('results-failed')].join('|')",
        secs(20),
    ));
    check(
        "person: an OK results.json draws as before, with no failure line",
        ok == "Assessed, no finding|false",
        ok,
    );

    // The home page: runs.json at 500, then no answer at all (offline), then Try again once it answers.
    for (how, a) in [
        ("at 500", Answer::Status(500)),
        ("aborted (offline)", Answer::Abort),
    ] {
        stub.set("/runs.json", a);
        b.go(&stub.url("/")).unwrap();
        let said = text(b.wait(&failed_box("home-failed"), secs(20)));
        let kept = text(
            b.eval("String(document.querySelectorAll('#home ul.runs li a').length)")
                .ok(),
        );
        check(
            &format!(
                "home: runs.json {how} says 'The list of people could not be loaded.' with Try again, above the server's list"
            ),
            said == "The list of people could not be loaded. Try again" && kept == "1",
            format!("{said:?} links {kept}"),
        );
    }
    stub.set("/runs.json", Answer::File);
    b.eval("document.querySelector('#home-failed button').click(); true")
        .ok();
    let people = text(b.wait(
        "!document.getElementById('home-failed') && document.documentElement.dataset.pgpStatus === 'ready' && \
         document.documentElement.dataset.pgpPeople",
        secs(10),
    ));
    check(
        "home: once runs.json answers, Try again draws the people",
        people == "1",
        people,
    );

    // The family page: family.json answers 200 {}, then the family through Try again.
    stub.set("/runs/DEMO-FAMILY/family.json", Answer::Json("{}"));
    b.go(&stub.url("/runs/DEMO-FAMILY/")).unwrap();
    let said = text(b.wait(
        "(function () { var m = document.getElementById('family'); var t = m.innerText.replace(/\\s+/g, ' ').trim(); \
         return m.querySelector('button') ? t : ''; }())",
        secs(20),
    ));
    check(
        "family: family.json answering {} says 'The family could not be loaded.' with Try again",
        said == "The family could not be loaded. Try again",
        said,
    );
    stub.set("/runs/DEMO-FAMILY/family.json", Answer::File);
    b.eval("var t = document.querySelector('#family button'); if (t) { t.click(); } !!t")
        .ok();
    let h1 = text(b.wait(
        "document.querySelector('#family h1') && document.querySelector('#family h1').textContent",
        secs(10),
    ));
    check(
        "family: once family.json answers, Try again draws the family",
        h1 == "Family DEMO-FAMILY",
        h1,
    );

    // A run whose state and results cannot be read: the gate stays closed and names no role.
    stub.set("/runs/DEMO-001/frames.json", Answer::Status(500));
    stub.set("/runs/DEMO-001/results.json", Answer::Status(500));
    b.go(&stub.url("/runs/DEMO-001/?unread")).unwrap();
    let page = text(b.wait(
        "document.getElementById('report-gate') && JSON.stringify([document.querySelector('#report-gate h1').textContent, \
         document.querySelector('#report-gate .lede').textContent, /child/i.test(document.body.innerText), \
         document.querySelectorAll('[data-finding], article.result, #results').length])",
        secs(20),
    ));
    let want = format!("[\"Results not loaded\",{},false,0]", js_str(UNREAD));
    check(
        "gate: DEMO-001 with frames.json and results.json at 500 reads 'Results not loaded', says no 'child' and draws no result",
        page == want,
        page,
    );
    stub.set("/runs/DEMO-001/frames.json", Answer::File);
    stub.set("/runs/DEMO-001/results.json", Answer::File);
    b.go(&stub.url("/runs/DEMO-003/")).unwrap();
    let child = text(b.wait(
        "document.getElementById('report-gate') && document.querySelector('#report-gate h1').textContent",
        secs(20),
    ));
    check(
        "control: DEMO-003, a real child with good files, still reads 'Results about a child'",
        child == "Results about a child",
        child,
    );

    // frames.json held 15 s: the panels wait, then the finished replay draws its steps and the re-check.
    stub.set("/runs/DEMO-001/frames.json", Answer::Hold(secs(15)));
    let t0 = Instant::now();
    b.go(&stub.url("/runs/DEMO-001/?held")).unwrap();
    let panels = "JSON.stringify([document.getElementById('found').innerText.trim(), \
         (document.querySelector('[data-plugin=\"ops\"] p') || {}).textContent || '', \
         (document.querySelector('[data-plugin=\"recheck\"] p') || {}).textContent || ''])";
    let waiting = text(b.wait(
        &format!("document.querySelector('[data-plugin=\"recheck\"] p') && document.querySelector('[data-plugin=\"ops\"] p') && {panels}"),
        secs(10),
    ));
    let absent = text(
        b.eval(
            "['Nothing finished yet.', 'No numbers yet.', 'No re-check yet'].filter(function (s) { \
         return document.body.innerText.indexOf(s) >= 0; }).join('|')",
        )
        .ok(),
    );
    let held = t0.elapsed() < secs(14);
    let want = format!("[{w},{w},{w}]", w = js_str(WAITING));
    check(
        "held: before frames.json answers, the wire, run health and re-check panels read 'Waiting for the run's data.'",
        held && waiting == want && absent.is_empty(),
        format!("{waiting} absent={absent:?} in time={held}"),
    );
    let arrived = text(b.wait(
        "document.documentElement.dataset.pgpFrame && JSON.stringify([document.querySelectorAll('#map g.box').length, \
         (document.querySelector('[data-plugin=\"recheck\"] p') || {}).textContent])",
        secs(30),
    ));
    check(
        "held: once the frames arrive, the finished replay shows its steps and the dated re-check",
        arrived == format!("[2,{}]", js_str(BANNER)),
        arrived,
    );

    drop((b, stub));
    assert!(fails.is_empty(), "failed: {fails:?}");
}
