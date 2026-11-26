//! `og live serve` on the synthetic DEMO-001 run (tests/fixtures/live): the
//! pages, `state.json`, `/runs.json`, the event stream, the run files, the
//! policy on every answer, the access log, the plugins' parts of
//! `state.json` and the family pages. These are the checks of the Python
//! server's tests that need no login; the login, the sessions and the
//! ingest listener are `port_live_auth_ingest*.rs`.

#[path = "port_live_server_util.rs"]
mod util;

use og_live::server::{self, Server};
use og_live::state;
use og_runlog::state::{RunFiles, Sources};
use pycompat::json::Value;
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::net::TcpStream;
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};
use util::*;

fn get_key<'a>(v: &'a Value, k: &str) -> &'a Value {
    match v {
        Value::Object(o) => o.get(k).unwrap_or(&Value::Null),
        _ => &Value::Null,
    }
}

fn s(v: &Value) -> &str {
    match v {
        Value::Str(x) => x,
        _ => "",
    }
}

fn steps_of(v: &Value) -> Vec<(String, String)> {
    match get_key(v, "steps") {
        Value::Array(a) => a
            .iter()
            .map(|x| {
                (
                    s(get_key(x, "step")).to_owned(),
                    s(get_key(x, "state")).to_owned(),
                )
            })
            .collect(),
        _ => Vec::new(),
    }
}

const CLINVAR_DONE: &str = r#"{"event": "process_completed", "utcTime": "2026-01-02T10:06:30Z", "trace": {"task_id": 9, "name": "CLINICAL:CLINVAR_SCREEN (DEMO-001)", "process": "CLINICAL:CLINVAR_SCREEN", "tag": "DEMO-001", "status": "COMPLETED", "exit": 0, "start": 1767348360000, "complete": 1767348390000}}"#;
const STEP06_OK: &str =
    r#"{"event": "pgp_step", "step": "06", "state": "ok", "utcTime": "2026-01-02T10:06:40Z"}"#;

/// The genome folder with DEMO-001's messages stored, and the site on it.
fn site(name: &str) -> (Work, Server) {
    let w = Work::new(name);
    make_root(&w);
    let ev = fixture_events();
    store(&w, RUN, &ev.iter().map(String::as_str).collect::<Vec<_>>());
    store(&w, RUN, &[CLINVAR_DONE, STEP06_OK]);
    let srv = start(config(&w));
    (w, srv)
}

#[test]
fn state_json_is_the_run_state_with_anchors_headlines_and_plugins() {
    let (w, srv) = site("state");
    let p = srv.port();
    let r = get(p, &format!("/runs/{RUN}/state.json"));
    assert_eq!(r.status, 200);
    assert!(policy_ok(&r), "{:?}", r.headers);
    let v = r.json();
    // The same files read by og-runlog directly give the same step states.
    let steps = og_runlog::steps::read_steps(&repo().join("assets/steps.tsv")).unwrap();
    let ev = events_path(&w, RUN);
    let files = RunFiles::read(
        &w.root().join(RUN),
        &Sources {
            traces: &[],
            events: Some(&ev),
            stale: false,
        },
    )
    .unwrap();
    let want = files
        .state(RUN, &steps, state::now(), &jiff::tz::TimeZone::system())
        .to_json();
    assert_eq!(steps_of(&v), steps_of(&want));
    let mine: std::collections::HashMap<_, _> = steps_of(&v).into_iter().collect();
    assert_eq!(
        (
            mine["06"].as_str(),
            mine["07"].as_str(),
            mine["11"].as_str(),
            mine["12"].as_str()
        ),
        ("done", "waiting", "done", "skipped")
    );
    let anchors = get_key(&v, "anchors");
    assert_eq!(s(get_key(anchors, "clinvar")), "clinvar");
    assert_eq!(s(get_key(anchors, "carrier")), "carrier-list");
    assert_eq!(s(get_key(anchors, "roh")), "roh");
    assert_eq!(
        get_key(&v, "events"),
        &Value::Int(fixture_events().len() as i128 + 2)
    );
    let Value::Object(plugins) = get_key(&v, "plugins") else {
        panic!("no plugins: {v:?}")
    };
    let names: Vec<&str> = plugins.keys().map(String::as_str).collect();
    assert_eq!(names, ["charts", "family", "ops", "recheck", "search"]);
    assert_eq!(get_key(get_key(&v, "plugins"), "family"), &Value::Null);
    assert!(matches!(
        get_key(get_key(get_key(&v, "plugins"), "ops"), "run"),
        Value::Object(_)
    ));
    assert!(matches!(get_key(&v, "headlines"), Value::Object(_)));
    assert_eq!(s(get_key(get_key(&v, "run"), "state")), "running");
}

#[test]
fn a_run_is_finished_while_pgp_run_finished_is_its_newest_message() {
    let w = Work::new("finished");
    make_root(&w);
    store(
        &w,
        RUN,
        &[
            r#"{"event": "started"}"#,
            r#"{"event": "pgp_run_finished", "exit": 0}"#,
        ],
    );
    let srv = start(config(&w));
    let state_of = || {
        s(get_key(
            get_key(
                &get(srv.port(), &format!("/runs/{RUN}/state.json")).json(),
                "run",
            ),
            "state",
        ))
        .to_owned()
    };
    assert_eq!(state_of(), "finished");
    store(
        &w,
        RUN,
        &[r#"{"event": "pgp_step", "step": "14", "state": "started"}"#],
    );
    assert_ne!(state_of(), "finished");
}

#[test]
fn a_rerun_is_the_new_run_with_every_message_counted() {
    let w = Work::new("rerun");
    let again = "DEMO-RERUN";
    let d = w.root().join(again);
    let now = state::now() as i64;
    write(
        &d.join("logs/run_status.tsv"),
        &format!("meta\tstarted_epoch\t{now}\nstep\t11\tplanned\n"),
    );
    let roh = r#"{"task_id": 8, "name": "CLINICAL:ROH (DEMO-RERUN)", "process": "CLINICAL:ROH", "tag": "DEMO-RERUN", "status": "COMPLETED", "exit": 0, "submit": 1767348300000, "start": 1767348302000, "complete": 1767348340000}"#;
    let done = format!(r#"{{"event": "process_completed", "trace": {roh}}}"#);
    let old = [
        r#"{"event": "started"}"#,
        done.as_str(),
        r#"{"event": "completed"}"#,
        r#"{"event": "pgp_run_finished", "exit": 0}"#,
        r#"{"event": "started"}"#,
    ];
    store(&w, again, &old);
    let srv = start(config(&w));
    let v = get(srv.port(), &format!("/runs/{again}/state.json")).json();
    let st: std::collections::HashMap<_, _> = steps_of(&v).into_iter().collect();
    assert_eq!(s(get_key(get_key(&v, "run"), "state")), "running");
    assert_ne!(st["11"], "done");
    assert_eq!(get_key(&v, "events"), &Value::Int(5));
    // The control: the whole file read as one run marks step 11 done.
    let all: Vec<Value> = old
        .iter()
        .map(|l| pycompat::json::loads(l).unwrap())
        .collect();
    assert_eq!(state::current_run(&all).len(), 1);
    assert_eq!(state::current_run(&all[..4]).len(), 4);
}

#[test]
fn headlines_read_as_the_report_does_and_say_nothing_is_fine() {
    let sm = |sec: &str| {
        og_render::model::Summary::parse(&format!(
            r#"{{"databases": {{}}, "sections": {{{sec}}}}}"#
        ))
        .unwrap()
    };
    let one = |key: &str, sec: &str| {
        let summary = sm(&format!(r#""{key}": {sec}"#));
        let raw = pycompat::json::loads(sec).unwrap();
        state::headline(key, Some(&raw), &summary)
    };
    let cp = |plp: &str| {
        format!(
            r#"{{"state": "ok", "title": "CPSR", "data": {{"plp": {plp}, "call_floor_dp": 10, "low_confidence": 0}}}}"#
        )
    };
    let tc: Vec<_> = ["[]", r#"[{"gene": "GENEA"}]"#]
        .iter()
        .map(|p| one("cpsr", &cp(p)).unwrap())
        .collect();
    let mito = |n: u8| {
        format!(
            r#"{{"state": "ok", "title": "Mito", "data": {{"annotation": {{"source": "VEP", "clinvar_plp": {n}}}}}}}"#
        )
    };
    let tm: Vec<_> = [0, 2]
        .iter()
        .map(|n| one("mito", &mito(*n)).unwrap())
        .collect();
    let tones: Vec<&str> = tc.iter().chain(&tm).map(|t| t.1).collect();
    assert_eq!(tones, ["info", "attention", "info", "attention"]);
    assert!(
        tc[1].0.contains("Pathogenic/Likely Pathogenic rows: 1"),
        "{}",
        tc[1].0
    );
    assert_eq!(tm[1].0, "ClinVar pathogenic or likely pathogenic: 2");
    let cv: Vec<_> = [0, 3]
        .iter()
        .map(|c| one("clinvar", &format!(r#"{{"state": "ok", "title": "ClinVar screen", "data": {{"count": {c}, "hits": []}}}}"#)).unwrap())
        .collect();
    assert_eq!((cv[0].1, cv[1].1), ("info", "attention"));
    assert!(
        cv.iter()
            .all(|t| t.0.contains(og_render::notes::CLINVAR_FREQ_NOTE))
    );
    assert_eq!(one("roh", r#"{"state": "missing"}"#), None);
    let failed = one(
        "roh",
        r#"{"state": "failed", "title": "ROH", "note": "the step failed"}"#,
    )
    .unwrap();
    assert_eq!(failed, ("ROH: the step failed".to_owned(), "attention"));
    // A helper that cannot read its data: the title and a pointer to the report.
    let odd = one(
        "clinvar",
        r#"{"state": "ok", "title": "ClinVar screen", "data": {"count": "x"}}"#,
    )
    .unwrap();
    assert_eq!(
        odd,
        (
            "ClinVar screen: finished, see the report".to_owned(),
            "info"
        )
    );
}

#[test]
fn the_run_files_the_report_links_and_nothing_else() {
    let (w, srv) = site("files");
    let p = srv.port();
    let mut seen = Vec::new();
    let mut go = |path: &str| {
        let r = get(p, path);
        seen.push((path.to_owned(), r.clone()));
        r
    };
    let r = go(&format!("/runs/{RUN}/report.html"));
    assert_eq!(r.status, 200);
    assert!(r.text().contains("clinvar/hits.tsv"));
    assert_eq!(
        r.header("content-security-policy"),
        Some(og_live::csp::RUN_HTML_CSP)
    );
    let r = go(&format!("/runs/{RUN}/clinvar/hits.tsv?x=1"));
    assert_eq!(
        (r.status, r.header("content-type")),
        (200, Some("text/plain; charset=utf-8"))
    );
    assert!(r.body.starts_with(b"gene\t"));
    for bad in [
        format!("/runs/{RUN}/clinvar/unlinked.tsv"),
        format!("/runs/{RUN}/clinvar/escape.tsv"),
        format!("/runs/{RUN}/logs/run_status.tsv"),
        format!("/runs/{RUN}/..%2F..%2Foutside.txt"),
        format!("/runs/{RUN}/cpsr/{RUN}.cpsr.grch38.classification.tsv.gz"),
        format!("/runs/{RUN}/twopages/"),
        format!("/runs/{RUN}/clinvar/"),
        format!("/runs/{RUN}/twopages/a.html"),
        "/runs/NO-SUCH-RUN/".to_owned(),
        "/runs/../x".to_owned(),
        "/nothing/here".to_owned(),
    ] {
        assert_eq!(go(&bad).status, 404, "{bad}");
    }
    let r = go(&format!("/runs/{RUN}/cpsr/"));
    let page = format!("/runs/{RUN}/cpsr/{RUN}.cpsr.grch38.html");
    assert_eq!((r.status, r.header("location")), (302, Some(page.as_str())));
    let r = go(&page);
    assert_eq!(r.status, 200);
    assert!(r.text().contains("<title>CPSR"));
    assert_eq!(
        r.header("content-security-policy"),
        Some(og_live::csp::RUN_HTML_CSP)
    );
    let r = go(&format!("/runs/{RUN}"));
    assert_eq!(
        (r.status, r.header("location")),
        (301, Some("/runs/DEMO-001/"))
    );
    let r = go("/");
    assert_eq!(r.status, 200);
    assert!(r.text().contains("/runs/DEMO-001/"));
    let r = go(&format!("/runs/{RUN}/"));
    assert_eq!(r.status, 200);
    assert!(r.text().contains("data-mode=\"live\"") || r.text().contains("live"));
    let r = go("/runs.json");
    let v = r.json();
    let Value::Array(people) = get_key(&v, "runs") else {
        panic!("{v:?}")
    };
    assert_eq!(s(get_key(&people[0], "name")), RUN);
    assert_eq!(s(get_key(&people[0], "state")), "running");
    assert_eq!(get_key(&people[0], "attention"), &Value::Null);
    assert_eq!(get_key(&v, "families"), &Value::Array(Vec::new()));
    let r = go("/assets/live.css");
    assert_eq!(
        (r.status, r.header("content-type")),
        (200, Some("text/css; charset=utf-8"))
    );
    assert_eq!(
        go("/assets/index.html").status,
        404,
        "a page template is never an asset"
    );
    assert_eq!(go("/assets/..%2Findex.html").status, 404);
    let r = go("/login?next=%2Fruns%2F");
    assert_eq!(r.status, 200);
    let r = go(&format!("/runs/{RUN}/results.json"));
    assert_eq!(r.status, 200);
    assert_eq!(
        s(get_key(&r.json(), "note")),
        og_live::plugins::NO_RESULTS_NOTE
    );
    assert_eq!(go(&format!("/runs/{RUN}/finding/clinvar-x-1/")).status, 404);
    // Every answer carries the policy, nosniff, no-referrer and no-store.
    let bad: Vec<&str> = seen
        .iter()
        .filter(|(_, r)| {
            !policy_ok(r) && r.header("content-security-policy") != Some(og_live::csp::RUN_HTML_CSP)
        })
        .map(|(p, _)| p.as_str())
        .collect();
    assert!(bad.is_empty(), "no policy on {bad:?}");
    // The access log: a row per page read, no query string.
    let rows = log_rows(&w);
    let has = |a: &str, run: &str, sec: &str, st: &str| {
        rows.iter()
            .any(|r| r.len() == 7 && r[3] == a && r[4] == run && r[5] == sec && r[6] == st)
    };
    assert!(has("report", RUN, "-", "200"), "{rows:?}");
    assert!(has("file", RUN, "clinvar/hits.tsv", "200"), "{rows:?}");
    assert!(has("file", RUN, "clinvar/escape.tsv", "404"), "{rows:?}");
    assert!(has("index", "-", "-", "200") && has("runs", "-", "-", "200"));
    assert!(
        !fs::read_to_string(w.state().join("access.log"))
            .unwrap()
            .contains('?')
    );
}

#[test]
fn a_foreign_host_and_a_post_from_another_site_are_refused() {
    let (w, srv) = site("host");
    let p = srv.port();
    let r = request(p, "GET", "/", &[("Host", "evil.example.net")]);
    assert_eq!(r.status, 403);
    assert!(policy_ok(&r));
    assert_eq!(
        request(p, "GET", "/", &[("Host", "localhost:8443")]).status,
        200
    );
    assert_eq!(request(p, "GET", "/", &[("Host", "[::1]:1")]).status, 200);
    assert_eq!(
        request(p, "GET", "/", &[("Host", "127.0.0.1:x")]).status,
        403
    );
    let evil = request(
        p,
        "POST",
        "/logout",
        &[
            ("Origin", "https://evil.example.net"),
            ("Content-Length", "0"),
        ],
    );
    assert_eq!(evil.status, 403);
    let site = request(
        p,
        "POST",
        "/logout",
        &[("Sec-Fetch-Site", "cross-site"), ("Content-Length", "0")],
    );
    assert_eq!(site.status, 403);
    let own = request(
        p,
        "POST",
        "/logout",
        &[("Origin", "null"), ("Content-Length", "0")],
    );
    assert_eq!((own.status, own.header("location")), (303, Some("/login")));
    assert!(own.header("set-cookie").unwrap_or("").contains("Max-Age=0"));
    let login = request(p, "POST", "/login", &[("Content-Length", "10")]);
    assert_eq!((login.status, login.header("location")), (303, Some("/")));
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
    assert_eq!(request(p, "PUT", "/", &[]).status, 501);
    let head = request(p, "HEAD", "/", &[]);
    assert_eq!(head.status, 200);
    assert!(head.body.is_empty());
    let rows = log_rows(&w);
    assert!(rows.iter().any(|r| r[3] == "host-refused" && r[6] == "403"));
    assert!(
        rows.iter()
            .any(|r| r[3] == "origin-refused" && r[6] == "403")
    );
}

#[test]
fn the_run_page_gets_one_script_tag_per_plugin_script() {
    let w = Work::new("plugin-tags");
    make_root(&w);
    let assets = w.dir.join("assets/live");
    fs::create_dir_all(assets.join("plugins")).unwrap();
    for f in ["index.html", "home.html", "login.html", "live.css"] {
        fs::copy(repo().join("assets/live").join(f), assets.join(f)).unwrap();
    }
    for f in ["a.js", "b.js", "notes.txt", ".hidden.js"] {
        write(&assets.join("plugins").join(f), "// plugin\n");
    }
    let mut c = config(&w);
    c.assets = assets;
    let srv = start(c);
    let page = get(srv.port(), &format!("/runs/{RUN}/")).text();
    let tags: Vec<&str> = page
        .split("<script src=\"/assets/plugins/")
        .skip(1)
        .filter_map(|x| x.split('"').next())
        .collect();
    assert_eq!(tags, ["a.js", "b.js"]);
    assert_eq!(get(srv.port(), "/assets/plugins/a.js").status, 200);
    assert_eq!(get(srv.port(), "/assets/plugins/.hidden.js").status, 404);
}

#[test]
fn search_charts_and_recheck_files_are_served_and_named_in_state_json() {
    let (w, srv) = site("run-globs");
    let d = w.root().join(RUN);
    write(
        &d.join(format!("{RUN}_search_index.json")),
        "{\"genes\": {}}\n",
    );
    write(&d.join(format!("{RUN}_charts.json")), "{}\n");
    write(&d.join("OTHER_charts.json"), "{}\n");
    write(
        &d.join("recheck/state.json"),
        r#"{"checked_utc": "2026-10-06T00:00:00Z", "items": []}"#,
    );
    let p = srv.port();
    let v = get(p, &format!("/runs/{RUN}/state.json")).json();
    let pl = get_key(&v, "plugins");
    assert_eq!(
        s(get_key(get_key(pl, "search"), "file")),
        format!("{RUN}_search_index.json")
    );
    assert_eq!(
        s(get_key(get_key(pl, "charts"), "file")),
        format!("{RUN}_charts.json")
    );
    assert!(matches!(get_key(get_key(pl, "search"), "mtime"), Value::Int(n) if *n > 0));
    assert_eq!(
        s(get_key(get_key(pl, "recheck"), "banner")),
        "Database re-check of 2026-10-06T00:00:00Z: nothing newer than this run"
    );
    for f in [
        format!("{RUN}_search_index.json"),
        format!("{RUN}_charts.json"),
        "OTHER_charts.json".into(),
        "recheck/state.json".into(),
    ] {
        let r = get(p, &format!("/runs/{RUN}/{f}"));
        assert_eq!(
            (r.status, r.header("content-type")),
            (200, Some("application/json")),
            "{f}"
        );
    }
    assert_eq!(get(p, &format!("/runs/{RUN}/recheck/log.tsv")).status, 404);
}

/// An open `/events` stream, read in a thread.
struct Stream {
    items: Arc<Mutex<Vec<(usize, Value)>>>,
    comments: Arc<AtomicUsize>,
    closed: Arc<AtomicBool>,
    status: Arc<Mutex<String>>,
    sock: TcpStream,
}

impl Stream {
    fn open(port: u16, last: Option<usize>) -> Stream {
        let mut sock = TcpStream::connect(("127.0.0.1", port)).unwrap();
        let mut head = format!("GET /runs/{RUN}/events HTTP/1.1\r\nHost: 127.0.0.1:{port}\r\n");
        if let Some(l) = last {
            head.push_str(&format!("Last-Event-ID: {l}\r\n"));
        }
        sock.write_all(format!("{head}\r\n").as_bytes()).unwrap();
        let (items, comments, closed, status) = (
            Arc::new(Mutex::new(Vec::new())),
            Arc::new(AtomicUsize::new(0)),
            Arc::new(AtomicBool::new(false)),
            Arc::new(Mutex::new(String::new())),
        );
        let r = sock.try_clone().unwrap();
        let (i2, c2, d2, s2) = (
            items.clone(),
            comments.clone(),
            closed.clone(),
            status.clone(),
        );
        thread::spawn(move || {
            let mut lines = BufReader::new(r).lines();
            if let Some(Ok(l)) = lines.next() {
                *s2.lock().unwrap() = l;
            }
            for l in lines.by_ref() {
                match l {
                    Ok(l) if l.is_empty() => break,
                    Ok(_) => {}
                    Err(_) => break,
                }
            }
            let mut id = 0;
            let mut data = None;
            for l in lines {
                let Ok(l) = l else { break };
                if let Some(x) = l.strip_prefix("id: ") {
                    id = x.parse().unwrap();
                } else if let Some(x) = l.strip_prefix("data: ") {
                    data = Some(pycompat::json::loads(x).unwrap());
                } else if l.starts_with(": heartbeat") {
                    c2.fetch_add(1, Ordering::SeqCst);
                } else if l.is_empty()
                    && let Some(d) = data.take()
                {
                    i2.lock().unwrap().push((id, d));
                }
            }
            d2.store(true, Ordering::SeqCst);
        });
        Stream {
            items,
            comments,
            closed,
            status,
            sock,
        }
    }

    fn wait(&self, n: usize, t: Duration) -> bool {
        let end = Instant::now() + t;
        while Instant::now() < end && self.items.lock().unwrap().len() < n {
            thread::sleep(Duration::from_millis(20));
        }
        self.items.lock().unwrap().len() >= n
    }

    fn items(&self) -> Vec<(usize, Value)> {
        self.items.lock().unwrap().clone()
    }
}

impl Drop for Stream {
    fn drop(&mut self) {
        let _ = self.sock.shutdown(std::net::Shutdown::Both);
    }
}

#[test]
fn the_event_stream_sends_each_stored_message_and_resumes_after_last_event_id() {
    let w = Work::new("stream");
    make_root(&w);
    let ev = fixture_events();
    store(&w, RUN, &ev.iter().map(String::as_str).collect::<Vec<_>>());
    let mut c = config(&w);
    c.heartbeat = Duration::from_millis(300);
    let srv = start(c);
    let p = srv.port();
    let n0 = ev.len();
    let st = Stream::open(p, None);
    assert!(st.wait(n0, Duration::from_secs(5)), "{:?}", st.items());
    assert!(st.status.lock().unwrap().contains(" 200 "));
    let t0 = Instant::now();
    store(&w, RUN, &[CLINVAR_DONE]);
    assert!(st.wait(n0 + 1, Duration::from_secs(2)), "not within 2 s");
    let last = st.items().last().unwrap().1.clone();
    assert_eq!(
        (s(get_key(&last, "event")), s(get_key(&last, "step"))),
        ("process_completed", "06")
    );
    assert!(t0.elapsed() < Duration::from_secs(2));
    // What the stream sent is messages() of the stored file.
    let steps = og_runlog::steps::read_steps(&repo().join("assets/steps.tsv")).unwrap();
    let stored = og_runlog::files::read_events(Some(&events_path(&w, RUN))).unwrap();
    let want = state::messages(&stored, &steps, 0);
    assert_eq!(st.items(), want);
    let text = state::dumps(&Value::Array(want.iter().map(|(_, v)| v.clone()).collect()));
    for leak in [
        "/data/",
        "workDir",
        "outdir",
        "\"tag\"",
        "HP:0001250",
        "bcftools",
    ] {
        assert!(!text.contains(leak), "{leak} in {text}");
    }
    thread::sleep(Duration::from_millis(800));
    assert!(st.comments.load(Ordering::SeqCst) >= 1, "no heartbeat");
    assert!(
        !st.closed.load(Ordering::SeqCst),
        "a quiet stream was closed"
    );
    let st2 = Stream::open(p, Some(n0 - 1));
    assert!(st2.wait(2, Duration::from_secs(5)));
    let ids: Vec<usize> = st2.items().iter().map(|(i, _)| *i).collect();
    assert_eq!(ids, [n0, n0 + 1]);
    // Every slot taken: 503, and the stream is not opened.
    srv.live().max_streams.store(2, Ordering::SeqCst);
    let r = get(p, &format!("/runs/{RUN}/events"));
    assert_eq!(
        (r.status, r.text().as_str()),
        (503, "too many open event streams\n")
    );
    drop(st);
    drop(st2);
    let end = Instant::now() + Duration::from_secs(5);
    while srv.live().streams.load(Ordering::SeqCst) > 0 && Instant::now() < end {
        thread::sleep(Duration::from_millis(50));
    }
    assert_eq!(
        srv.live().streams.load(Ordering::SeqCst),
        0,
        "a closed stream kept its slot"
    );
    let rows = log_rows(&w);
    assert!(rows.iter().any(|r| r[3] == "events" && r[6] == "503"));
}

#[test]
fn the_command_line_refuses_what_it_cannot_do_yet() {
    let w = Work::new("cli");
    make_root(&w);
    let root = w.root().to_string_lossy().into_owned();
    let run = |a: &[&str]| {
        let args: Vec<String> = a.iter().map(|x| (*x).to_owned()).collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = server::cli(&args, &mut out, &mut err);
        (code, String::from_utf8_lossy(&err).into_owned())
    };
    // Through configure_with with no hash: an exported
    // OWNGENOME_LIVE_PASSWORD_HASH would let cli start and never return.
    let bind_any: Vec<String> = ["--root", &root, "--bind", "0.0.0.0", "--port", "0"]
        .iter()
        .map(|x| (*x).to_owned())
        .collect();
    let (code, err) = server::configure_with(&bind_any, None).err().unwrap();
    assert!(code == 2 && err.contains("password"), "{code} {err}");
    let (code, err) = run(&["--root", &root, "--password-file", "/nonexistent"]);
    assert!(code == 1 && err.contains("password hash"), "{code} {err}");
    let (code, err) = run(&["--port", "0"]);
    assert!(
        code == 2 && err.contains("--root is required"),
        "{code} {err}"
    );
    let (code, err) = run(&["--root", "/nonexistent/genome", "--port", "0"]);
    assert!(code == 1 && err.contains("not a folder"), "{code} {err}");
    let conf = w.dir.join("live.conf");
    write(&conf, &format!("root = {root}\nmax-streams = 3\n"));
    let c = server::configure(&["--config".into(), conf.to_string_lossy().into_owned()]).unwrap();
    assert_eq!((c.root.as_path(), c.max_streams), (Path::new(&root), 3));
    assert_eq!(c.state_dir, w.root().join(".owngenome-live"));
}

#[test]
fn a_family_folder_is_a_run_with_its_own_page_and_files() {
    let w = Work::new("family");
    make_root(&w);
    let fam = "DEMO-FAMILY";
    let fd = w.root().join(fam);
    fs::create_dir_all(&fd).unwrap();
    for e in fs::read_dir(repo().join("tests/fixtures/live/demo").join(fam)).unwrap() {
        let e = e.unwrap();
        fs::copy(e.path(), fd.join(e.file_name())).unwrap();
    }
    let srv = start(config(&w));
    let p = srv.port();
    let v = get(p, "/runs.json").json();
    let Value::Array(fams) = get_key(&v, "families") else {
        panic!("{v:?}")
    };
    assert_eq!(fams.len(), 1);
    assert_eq!(s(get_key(&fams[0], "name")), fam);
    assert_eq!(s(get_key(&fams[0], "url")), "/runs/DEMO-FAMILY/");
    assert!(matches!(get_key(&fams[0], "members"), Value::Array(m) if !m.is_empty()));
    let page = get(p, &format!("/runs/{fam}/")).text();
    let want = fs::read_to_string(repo().join("assets/live/family.html")).unwrap();
    assert!(
        page.starts_with(&want[..want.find("{{").unwrap()]),
        "not the family page"
    );
    let st = get(p, &format!("/runs/{fam}/state.json")).json();
    assert_eq!(
        s(get_key(get_key(get_key(&st, "plugins"), "family"), "kind")),
        "family"
    );
    let r = get(p, &format!("/runs/{fam}/family.json"));
    assert_eq!(r.status, 200, "{}", r.text());
    let page = r.json();
    assert!(
        matches!(get_key(&page, "layout"), Value::Object(_)),
        "{page:?}"
    );
    let runs = get_key(&page, "runs");
    assert_eq!(
        s(get_key(get_key(runs, RUN), "url")),
        "/runs/DEMO-001/",
        "{runs:?}"
    );
    assert_eq!(get_key(get_key(runs, RUN), "child"), &Value::Bool(false));
    let r = get(p, &format!("/family/{fam}/pedigree.svg"));
    assert_eq!(
        (r.status, r.header("content-type")),
        (200, Some("image/svg+xml"))
    );
    assert!(r.text().contains("<svg"));
    let r = get(p, &format!("/family/{fam}/relationships.tsv"));
    assert_eq!(
        (r.status, r.header("content-type")),
        (200, Some("text/plain; charset=utf-8"))
    );
    for bad in [
        format!("/family/{fam}/de_novo.tsv"),
        format!("/family/{fam}/nothing.tsv"),
        format!("/family/{RUN}/pedigree.svg"),
        format!("/runs/{RUN}/family.json"),
        "/family/NO-SUCH/pedigree.svg".into(),
    ] {
        assert_eq!(get(p, &bad).status, 404, "{bad}");
    }
    // A member's page names its family; DEMO-001 has no parents listed.
    let st = get(p, &format!("/runs/{RUN}/state.json")).json();
    let m = get_key(get_key(&st, "plugins"), "family");
    assert_eq!(s(get_key(m, "kind")), "member", "{m:?}");
    assert_eq!(get_key(m, "child"), &Value::Bool(false));
    let Value::Array(listed) = get_key(m, "families") else {
        panic!("{m:?}")
    };
    assert_eq!(s(get_key(&listed[0], "url")), "/runs/DEMO-FAMILY/");
    let rows = log_rows(&w);
    assert!(
        rows.iter()
            .any(|r| r[3] == "family" && r[4] == fam && r[5] == "pedigree.svg" && r[6] == "200")
    );
}

#[test]
fn cli_configuration_serves_the_journal_and_finding_pages() {
    let w = Work::new("cli-results");
    let root = make_root(&w);
    let events = fixture_events();
    store(
        &w,
        RUN,
        &events.iter().map(String::as_str).collect::<Vec<_>>(),
    );
    store(&w, RUN, &[CLINVAR_DONE, STEP06_OK]);
    write(
        &root.join(RUN).join("summary.json"),
        r#"{"sample":"DEMO-001","sections":{"clinvar":{"state":"ok","data":{"hits":[
            {"gene":"GENEA","chrom":"chr1","pos":100,"ref":"A","alt":"G",
             "significance":"Pathogenic","stars":2}
        ]}}}}"#,
    );
    let assets = repo().join("assets");
    let args = [
        "--root",
        root.to_str().unwrap(),
        "--state-dir",
        w.state().to_str().unwrap(),
        "--assets",
        assets.join("live").to_str().unwrap(),
        "--steps",
        assets.join("steps.tsv").to_str().unwrap(),
        "--port",
        "0",
        "--ingest-port",
        "0",
        "--quiet",
    ]
    .map(str::to_owned);
    let srv = start(server::configure_with(&args, None).unwrap());
    let port = srv.port();
    let response = get(port, &format!("/runs/{RUN}/results.json"));
    assert_eq!(response.status, 200);
    let results = response.json();
    assert_eq!(s(get_key(&results, "sample")), RUN);
    assert!(
        matches!(get_key(&results, "chapters"), Value::Array(chapters) if chapters.len() == 7),
        "{results:?}"
    );
    let Value::Array(findings) = get_key(&results, "findings") else {
        panic!("missing journal findings");
    };
    let id = s(get_key(&findings[0], "id"));
    assert_eq!(get(port, &format!("/runs/{RUN}/finding/{id}/")).status, 200);
    assert_eq!(
        get(port, &format!("/runs/{RUN}/finding/unknown/")).status,
        404
    );
    assert_eq!(get(port, &format!("/runs/{RUN}/summary.json")).status, 404);

    let index = get(port, "/runs.json").json();
    let Value::Array(runs) = get_key(&index, "runs") else {
        panic!("missing runs");
    };
    let run = runs
        .iter()
        .find(|run| s(get_key(run, "name")) == RUN)
        .unwrap();
    assert_eq!(get_key(run, "screened"), &Value::Bool(true));
    let rows = log_rows(&w);
    assert_eq!(
        rows.iter()
            .filter(|row| row[3] == "results" && row[4] == RUN)
            .count(),
        1
    );
}
