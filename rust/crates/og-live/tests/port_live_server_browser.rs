//! The live run page served by `og live serve`, in a headless Chrome
//! (og-browser; skipped without Chrome, except in CI): the step map, a box
//! that repaints after a stored message, a stream refused with 503 and
//! opened again, a box's detail with its tasks and its report anchors, and
//! no policy violation.

#[path = "port_live_server_util.rs"]
mod util;

use og_browser::{Browser, Value, chrome_or_skip};
use std::sync::atomic::Ordering;
use std::thread;
use std::time::Duration;
use util::*;

fn text(v: Option<Value>) -> String {
    match v {
        Some(Value::Str(s)) => s,
        Some(Value::Int(n)) => n.to_string(),
        Some(Value::Bool(b)) => b.to_string(),
        _ => String::new(),
    }
}

fn box_state(b: &mut Browser, step: &str) -> String {
    text(b.eval(&format!("(function(){{var g=document.querySelector('g.box[data-step=\"{step}\"]');return g?g.dataset.state:'';}})()")).ok())
}

#[test]
fn the_run_page_draws_the_steps_and_follows_the_stream() {
    let Some(chrome) = chrome_or_skip("load the live run page") else {
        return;
    };
    let w = Work::new("browser");
    make_root(&w);
    let ev = fixture_events();
    store(&w, RUN, &ev.iter().map(String::as_str).collect::<Vec<_>>());
    store(
        &w,
        RUN,
        &[
            r#"{"event": "process_completed", "utcTime": "2026-01-02T10:06:30Z", "trace": {"task_id": 9, "name": "CLINICAL:CLINVAR_SCREEN (DEMO-001)", "process": "CLINICAL:CLINVAR_SCREEN", "tag": "DEMO-001", "status": "COMPLETED", "exit": 0, "start": 1767348360000, "complete": 1767348390000}}"#,
            r#"{"event": "pgp_step", "step": "06", "state": "ok", "utcTime": "2026-01-02T10:06:40Z"}"#,
        ],
    );
    let srv = start(config(&w));
    let base = format!("http://127.0.0.1:{}", srv.port());
    let mut b = Browser::launch(&chrome, &w.dir.join("chrome")).unwrap();
    let wait = |b: &mut Browser, expr: &str, secs: u64| b.wait(expr, Duration::from_secs(secs));

    b.go(&format!("{base}/runs/{RUN}/")).unwrap();
    let ready = wait(
        &mut b,
        "document.documentElement.dataset.pgpStatus === 'ready'",
        30,
    );
    assert!(
        ready.is_some(),
        "the run page did not load: {}",
        text(
            b.eval("location.href + ' ' + document.documentElement.outerHTML.slice(0, 600)")
                .ok()
        )
    );
    let n = text(
        b.eval("document.querySelectorAll('#map g.box').length")
            .ok(),
    );
    assert!(n.parse::<usize>().is_ok_and(|n| n >= 4), "boxes: {n}");
    assert_eq!(box_state(&mut b, "06"), "done");

    // A stored message repaints its box.
    assert_eq!(box_state(&mut b, "07"), "waiting");
    store(
        &w,
        RUN,
        &[
            r#"{"event": "process_started", "utcTime": "2026-01-02T10:07:00Z", "trace": {"task_id": 20, "name": "PGX:PHARMCAT_PREPROCESS (DEMO-001)", "process": "PGX:PHARMCAT_PREPROCESS", "tag": "DEMO-001", "status": "RUNNING", "start": 1767348420000}}"#,
        ],
    );
    let after = wait(
        &mut b,
        "document.querySelector('g.box[data-step=\"07\"]').dataset.state === 'running' && 'running'",
        10,
    );
    assert_eq!(text(after), "running", "07 did not repaint");

    // Every stream slot taken: the stream answers 503, and the page opens it
    // again once a slot is free.
    let live = srv.live();
    let keep = live.max_streams.swap(0, Ordering::SeqCst);
    b.go(&format!("{base}/runs/{RUN}/")).unwrap();
    let refused = wait(
        &mut b,
        "document.documentElement.dataset.pgpStatus === 'offline' || \
         (document.documentElement.dataset.pgpStreams === '1' && \
          document.documentElement.dataset.pgpStatus === 'ready' && window.PGP.state && 'loaded')",
        10,
    );
    thread::sleep(Duration::from_secs(1));
    live.max_streams.store(keep, Ordering::SeqCst);
    let reopened = wait(
        &mut b,
        "Number(document.documentElement.dataset.pgpStreams) >= 2 && 'reopened'",
        15,
    );
    let n0 = text(b.eval("window.PGP.state.events").ok());
    thread::sleep(Duration::from_secs(1));
    store(
        &w,
        RUN,
        &[
            r#"{"event": "pgp_step", "step": "07", "state": "ok", "utcTime": "2026-01-02T10:07:30Z"}"#,
        ],
    );
    let n1 = wait(
        &mut b,
        &format!("window.PGP.state && window.PGP.state.events > {n0} && window.PGP.state.events"),
        10,
    );
    assert!(
        refused.is_some() && text(reopened.clone()) == "reopened" && n1.is_some(),
        "{refused:?} {reopened:?} {n0} {n1:?}"
    );

    // A box's detail: its tasks, and a link to each section's report anchor.
    b.eval("document.querySelector('g.box[data-step=\"11\"]').dispatchEvent(new MouseEvent('click', {bubbles: true})); true").unwrap();
    let detail = text(wait(
        &mut b,
        "document.getElementById('detail').textContent.indexOf('Step 11') >= 0 && document.getElementById('detail').textContent",
        20,
    ));
    assert!(detail.contains("ROH"), "{detail}");
    b.eval("document.querySelector('g.box[data-step=\"06\"]').dispatchEvent(new MouseEvent('click', {bubbles: true})); true").unwrap();
    let links = text(wait(
        &mut b,
        "document.getElementById('detail').textContent.indexOf('Step 06') >= 0 && \
         Array.from(document.querySelectorAll('#detail a')).map(function (a) { return a.getAttribute('href'); }).join(' ')",
        20,
    ));
    assert!(
        links.contains("report.html#clinvar") && links.contains("report.html#carrier-list"),
        "{links}"
    );
    let blocked = text(
        b.eval("document.documentElement.dataset.pgpBlocked || ''")
            .ok(),
    );
    assert_eq!(blocked, "", "the policy blocked what the page needed");
    let errors = b.page_errors();
    assert!(errors.is_empty(), "page errors: {errors:?}");
}
