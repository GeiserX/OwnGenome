//! The live page in replay mode, in a headless Chrome: assets/live's
//! index.html with live.js and live.css, playing tests/fixtures/live/frames.json
//! from static files served with the live site's Content-Security-Policy.
//! Moved here from tests/test_owngenome_live.py, check for check.
//!
//! Red control: the same page without frames.json asks for it, gets a
//! 404 and must fail the first check. Without Chrome the test says SKIP,
//! and fails in CI.

use og_browser::{Browser, Scratch, StaticServer, Value, chrome_or_skip, repo_root};
use pycompat::json::{self, Dump};
use std::fs;
use std::path::Path;
use std::time::Duration;

const RUN: &str = "DEMO-001";
/// The policy the live site sends (`CSP` of og live serve,
/// rust/crates/og-live/src/csp.rs; the test `the_policy_is_the_live_sites`
/// keeps the two equal).
const CSP: &str = "default-src 'none'; script-src 'self'; connect-src 'self'; style-src 'self'; \
     img-src 'self' data:; form-action 'self'; base-uri 'none'; frame-src 'self'; frame-ancestors 'none'";

/// The replay page in `dir`: index.html filled in for replay, live.js,
/// live.css and, when `frames` is given, frames.json.
fn replay_site(dir: &Path, frames: Option<&str>) {
    let assets = repo_root().join("assets/live");
    fs::create_dir_all(dir).unwrap();
    for f in ["live.js", "live.css"] {
        fs::copy(assets.join(f), dir.join(f)).unwrap();
    }
    let page = fs::read_to_string(assets.join("index.html"))
        .unwrap()
        .replace("{{assets}}", ".")
        .replace("{{run}}", RUN)
        .replace("{{mode}}", "replay")
        .replace("{{plugins}}", "");
    fs::write(dir.join("index.html"), page).unwrap();
    if let Some(f) = frames {
        fs::write(dir.join("frames.json"), f).unwrap();
    }
}

fn text(v: Option<Value>) -> String {
    match v {
        Some(Value::Str(s)) => s,
        Some(v) => format!("{v:?}"),
        None => "None".into(),
    }
}

const OPENED: &str = "document.documentElement.dataset.pgpFrame && JSON.stringify([document.documentElement.dataset.pgpFrame, \
     document.querySelector('g.box[data-step=\"11\"]').dataset.state, \
     document.getElementById('replay-time').textContent, document.getElementById('replay-time').title])";

#[test]
fn replay_plays_frames_json() {
    let Some(chrome) = chrome_or_skip("load the replay page") else {
        return;
    };
    // Declared before the server and the browser: dropped after them, also on a panic.
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "replay");
    let work = scratch.path();
    let frames = fs::read_to_string(repo_root().join("tests/fixtures/live/frames.json")).unwrap();
    let rd = work.join("replay");
    replay_site(&rd, Some(&frames));
    // A replay whose last frame is still running (its own folder: no cached frames.json).
    let mut fr = json::loads(&frames).unwrap();
    if let Value::Object(m) = &mut fr
        && let Some(Value::Array(a)) = m.get_mut("frames")
    {
        a.pop();
    }
    replay_site(
        &rd.join("running"),
        Some(&json::dumps(&fr, &Dump::default()).unwrap()),
    );
    // The red control: no frames.json at all.
    replay_site(&rd.join("no-frames"), None);

    let srv = StaticServer::start(&rd, &[("Content-Security-Policy", CSP)]).unwrap();
    let mut b = Browser::launch(&chrome, &work.join("chrome")).unwrap();
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

    b.go(&srv.url("/index.html")).unwrap();
    // A finished run opens at its last frame, its moment short with the full time in the title.
    let opened = text(b.wait(OPENED, secs(20)));
    check(
        "browser: replay of a finished run opens at its end, its time short ('2 Jan 10:05') with the full time in the title",
        opened == r#"["2","done","2 Jan 10:05","2026-01-02T10:05:40Z"]"#,
        opened,
    );
    // Play from the end starts again at the first frame: every frame shown is recorded with step 11's state.
    b.eval(
        "window.seen = []; new MutationObserver(function () { var r = document.documentElement; \
         window.seen.push(r.dataset.pgpFrame + ':' + document.querySelector('g.box[data-step=\"11\"]').dataset.state); })\
         .observe(document.documentElement, {attributes: true, attributeFilter: ['data-pgp-frame']}); \
         document.getElementById('replay-play').click(); true",
    )
    .ok();
    let seen = b.wait(
        "window.seen.length > 1 && window.seen[window.seen.length - 1].indexOf('2:') === 0 && JSON.stringify(window.seen)",
        secs(10),
    );
    let seen: Vec<String> = match seen.as_ref().map(|v| json::loads(&text(Some(v.clone())))) {
        Some(Ok(Value::Array(a))) => a.into_iter().map(|v| text(Some(v))).collect(),
        _ => Vec::new(),
    };
    let found = text(b.eval("document.getElementById('found').textContent").ok());
    check(
        "browser: Play at the end replays frames.json from the first frame (11 waiting at frame 0, done at the last, its headline listed)",
        seen.first().map(String::as_str) == Some("0:waiting")
            && seen.last().map(String::as_str) == Some("2:done")
            && found.contains("Runs of homozygosity"),
        format!("{seen:?} {found}"),
    );
    let aria = text(
        b.eval(
            "var m = document.querySelector('#map svg'), g = document.querySelector('g.box[data-step=\"11\"]'); \
             [m.getAttribute('role'), g.getAttribute('aria-label') || ''].join('|')",
        )
        .ok(),
    );
    check(
        "browser: the map is a group (not an image that hides its buttons) and a box names its step and state",
        aria.starts_with("group|11 ") && aria.ends_with(": done"),
        aria,
    );
    b.eval(
        "document.getElementById('detail').textContent = ''; var li = document.querySelector('#found li[role=button]'); \
         li.focus(); li.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', bubbles: true})); true",
    )
    .ok();
    let keyed = b.wait(
        "document.activeElement && document.activeElement.matches('#found li') && \
         document.getElementById('detail').textContent.indexOf('Step 11') >= 0 && 'opened'",
        secs(5),
    );
    let focused = text(
        b.eval("document.activeElement && document.activeElement.outerHTML.slice(0, 200)")
            .ok(),
    );
    check(
        "browser: a headline takes focus and Enter on it opens its step",
        keyed == Some(Value::Str("opened".into())),
        focused,
    );
    let blocked = text(
        b.eval("document.documentElement.dataset.pgpBlocked || ''")
            .ok(),
    );
    check(
        "browser: the policy blocked nothing in replay",
        blocked.is_empty(),
        blocked,
    );

    // A replay whose last frame is still running opens at its first frame.
    b.go(&srv.url("/running/index.html")).unwrap();
    let at = b.wait(
        "document.documentElement.dataset.pgpFrame && document.getElementById('replay-slider').max === '1' \
         && document.documentElement.dataset.pgpFrame",
        secs(20),
    );
    check(
        "browser: a replay whose last frame is still running opens at its first frame",
        at == Some(Value::Str("0".into())),
        text(at.clone()),
    );

    // The red control: without frames.json the first check must fail. First
    // the proof the page ran and asked: its request for frames.json got a 404.
    b.go(&srv.url("/no-frames/index.html")).unwrap();
    let asked = b.wait(
        "(performance.getEntriesByType('resource').filter(function (e) { \
         return /\\/no-frames\\/frames\\.json$/.test(e.name); })[0] || {}).responseStatus",
        secs(20),
    );
    let none = b.wait(OPENED, secs(5));
    let frame = b.eval("document.documentElement.dataset.pgpFrame || ''");
    check(
        "control: without frames.json the page gets a 404 for it and never opens at a frame",
        asked == Some(Value::Int(404)) && none.is_none() && frame == Ok(Value::Str(String::new())),
        format!("{asked:?} {none:?} {frame:?}"),
    );

    drop((b, srv));
    assert!(fails.is_empty(), "failed: {fails:?}");
}

/// The replay is served with the live site's own policy: `CSP` of og live
/// serve (rust/crates/og-live/src/csp.rs, a Rust string whose lines end in
/// `\` and go on after the next line's indent) must equal [`CSP`].
#[test]
fn the_policy_is_the_live_sites() {
    let src = fs::read_to_string(repo_root().join("rust/crates/og-live/src/csp.rs")).unwrap();
    let open = "pub const CSP: &str = \"";
    let start = src.find(open).expect("csp.rs defines pub const CSP") + open.len();
    let end = start + src[start..].find("\";").expect("pub const CSP ends");
    let mut policy = String::new();
    for (i, piece) in src[start..end].split("\\\n").enumerate() {
        policy.push_str(if i == 0 { piece } else { piece.trim_start() });
    }
    assert_eq!(policy, CSP);
}
