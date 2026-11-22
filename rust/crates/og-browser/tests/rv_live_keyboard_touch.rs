//! The live pages for a keyboard and for a finger, in a headless Chrome, on
//! the demo the docs publish (tests/demo/build_live_demo.py, run with the og
//! this test run built).
//!
//! 1. On the home, person, finding and family pages the first Tab focuses
//!    the skip link, which then shows and points at the page's main; Enter
//!    moves the focus into main.
//! 2. On DEMO-001 at 320x800 and at 390x844, every caveat's More
//!    (`summary.toggle`) has a min-height of at least 44 px and is that
//!    tall, and so are the replay's slider and its Play button. Each More
//!    is named for its own finding, so no two share a name.
//! 3. At 390x844, a tap on the Cyrius row of the step timeline brings the
//!    step's detail into view, clear of the masthead and the tab bar, and
//!    the timeline's box does not scroll inside the page.
//!
//! Each check was red on main before this change (the PR says where).
//! Without Chrome, python3 or a built og the test says SKIP, and fails in CI.

use og_browser::{Browser, Scratch, StaticServer, Value, chrome_or_skip, js_str, repo_root};
use pycompat::json;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;

const RUN: &str = "DEMO-001";
/// A finding of the demo with a page of its own.
const FINDING: &str = "clinvar-demoa-chr2-47403191-c-t";
/// The step of the timeline check 3 taps.
const CYRIUS: &str = "21";
/// The least height a finger can hit, in CSS px.
const TAP: f64 = 44.0;

fn text(v: Option<Value>) -> String {
    match v {
        Some(Value::Str(s)) => s,
        Some(v) => format!("{v:?}"),
        None => "None".into(),
    }
}

fn num(v: &Value) -> f64 {
    match v {
        Value::Int(i) => *i as f64,
        Value::Float(f) => *f,
        _ => f64::NAN,
    }
}

/// A tool this test needs, or `None` after a `[SKIP]` line; in CI it panics.
fn need(found: Option<PathBuf>, what: &str) -> Option<PathBuf> {
    if found.is_some() {
        return found;
    }
    if std::env::var("CI").is_ok_and(|v| v == "true") {
        panic!("[FAIL] no {what}: CI must build the demo for this check");
    }
    println!("[SKIP] no {what} to build the demo");
    None
}

/// The og to build the demo with: OG when set and not empty, else the og
/// beside this test's own folder (target/<profile>/og, which
/// `cargo test --workspace` builds).
fn find_og() -> Option<PathBuf> {
    if let Some(o) = std::env::var_os("OG").filter(|o| !o.is_empty()) {
        let p = PathBuf::from(o);
        return p.is_file().then_some(p);
    }
    let me = std::env::current_exe().ok()?;
    let p = me.parent()?.parent()?.join("og");
    p.is_file().then_some(p)
}

fn find_python() -> Option<PathBuf> {
    let ok = Command::new("python3")
        .arg("-c")
        .arg("import sys; sys.exit(sys.version_info < (3, 10))")
        .status()
        .is_ok_and(|s| s.success());
    ok.then(|| PathBuf::from("python3"))
}

/// The page's size, as a reader's window; on a phone with touch.
fn viewport(b: &mut Browser, width: u32, height: u32, mobile: bool) {
    let p = json::loads(&format!(
        r#"{{"width": {width}, "height": {height}, "deviceScaleFactor": 1, "mobile": {mobile}}}"#
    ))
    .unwrap();
    b.call("Emulation.setDeviceMetricsOverride", p).unwrap();
    let t = json::loads(&format!(r#"{{"enabled": {mobile}, "maxTouchPoints": 1}}"#)).unwrap();
    b.call("Emulation.setTouchEmulationEnabled", t).unwrap();
}

/// One key pressed and let go, as a keyboard sends it.
fn press(b: &mut Browser, key: &str, code: u32, text: &str) {
    let down = if text.is_empty() {
        "rawKeyDown"
    } else {
        "keyDown"
    };
    for (kind, t) in [(down, text), ("keyUp", "")] {
        let p = json::loads(&format!(
            r#"{{"type": "{kind}", "key": {k}, "code": {k}, "windowsVirtualKeyCode": {code}, "nativeVirtualKeyCode": {code}, "text": {t}}}"#,
            k = js_str(key),
            t = js_str(t)
        ))
        .unwrap();
        b.call("Input.dispatchKeyEvent", p).unwrap();
    }
}

/// A finger's tap at (x, y) of the viewport.
fn tap(b: &mut Browser, x: f64, y: f64) {
    for (kind, points) in [
        ("touchStart", format!(r#"[{{"x": {x}, "y": {y}}}]"#)),
        ("touchEnd", "[]".to_string()),
    ] {
        let p = json::loads(&format!(r#"{{"type": "{kind}", "touchPoints": {points}}}"#)).unwrap();
        b.call("Input.dispatchTouchEvent", p).unwrap();
    }
}

/// Open `url` from a blank page and wait until the page has drawn.
fn visit(b: &mut Browser, url: &str) {
    b.go("about:blank").unwrap();
    assert!(
        b.wait("location.href === 'about:blank'", Duration::from_secs(10))
            .is_some(),
        "about:blank never opened"
    );
    b.go(url).unwrap();
    let drawn = "document.documentElement.getAttribute('data-pgp-status') === 'ready' || \
         !!document.querySelector('#family h1') || !!document.querySelector('#finding h1')";
    assert!(
        b.wait(drawn, Duration::from_secs(30)).is_some(),
        "{url} never drew"
    );
}

/// The first Tab stop of the page and where Enter on it moves the focus:
/// [the focused element's tag, its class, its href, whether it shows on
/// screen, main's id, whether the focus is inside main after Enter].
fn skip_link(b: &mut Browser) -> String {
    press(b, "Tab", 9, "");
    let first = text(
        b.eval(
            "(function () { var a = document.activeElement, r = a.getBoundingClientRect(); \
             return JSON.stringify([a.tagName, a.className, a.getAttribute('href') || '', \
             r.top >= 0 && r.bottom <= innerHeight && r.height > 0]); }())",
        )
        .ok(),
    );
    press(b, "Enter", 13, "\r");
    let after = text(b.wait(
        "(function () { var m = document.querySelector('main'); var a = document.activeElement; \
         return m && m.contains(a) && JSON.stringify([m.id, true]); }())",
        Duration::from_secs(3),
    ));
    format!("{first} {after}")
}

/// The sizes a finger meets on the person page: every `summary.toggle`
/// [its height, its min-height, its name, its finding's title], the
/// slider and the Play button.
const SIZES: &str = "(function () { \
     var t = Array.prototype.map.call(document.querySelectorAll('summary.toggle'), function (s) { \
       var row = s.closest('.row'), h = row && row.querySelector('h3'); \
       return [s.getBoundingClientRect().height, parseFloat(getComputedStyle(s).minHeight) || 0, \
         s.getAttribute('aria-label') || s.textContent, h ? h.textContent : '']; }); \
     return JSON.stringify([t, document.getElementById('replay-slider').getBoundingClientRect().height, \
       document.getElementById('replay-play').getBoundingClientRect().height]); }())";

/// Bring the Cyrius row to the middle of the screen (the page's scroll,
/// then the timeline box's own when it has one) and give its place.
const CYRIUS_ROW: &str = "(function () { \
     document.getElementById('progress-fold').open = true; \
     var g = document.querySelector('#map g.box[data-step=\"21\"] rect.hit'); if (!g) { return ''; } \
     var m = document.querySelector('section.map-wrap'), r = g.getBoundingClientRect(); \
     window.scrollBy(0, r.top - innerHeight / 2); r = g.getBoundingClientRect(); \
     if (r.top < 80 || r.bottom > innerHeight - 90) { m.scrollTop += r.top - innerHeight / 2; r = g.getBoundingClientRect(); } \
     return JSON.stringify([r.left, r.top, r.height]); }())";

/// After the tap: the detail's heading, its top and bottom, the masthead's
/// bottom, the tab bar's top, and the timeline box's scroll and client
/// heights.
const AFTER_TAP: &str = "(function () { \
     var d = document.getElementById('detail'), r = d.getBoundingClientRect(), h = d.querySelector('h2'); \
     var m = document.querySelector('section.map-wrap'); \
     return JSON.stringify([h ? h.textContent : '', r.top, r.bottom, \
       document.querySelector('header.mast').getBoundingClientRect().bottom, \
       document.getElementById('tabbar').getBoundingClientRect().top, m.scrollHeight, m.clientHeight]); }())";

#[test]
fn a_keyboard_skips_to_the_results_and_a_finger_hits_the_controls() {
    let Some(chrome) = chrome_or_skip("open the live demo") else {
        return;
    };
    let Some(og) = need(find_og(), "og (OG, or target/<profile>/og)") else {
        return;
    };
    let Some(python) = need(find_python(), "python3 3.10 or newer") else {
        return;
    };
    // Declared before the server and the browser: dropped after them, also on a panic.
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "keyboard-touch");
    let work = scratch.path();
    let site = work.join("site");
    let built = Command::new(&python)
        .arg(repo_root().join("tests/demo/build_live_demo.py"))
        .arg("--out")
        .arg(site.join("demo"))
        .arg("--work")
        .arg(work.join("root"))
        .env("OG", &og)
        .output()
        .unwrap();
    assert!(
        built.status.success(),
        "tests/demo/build_live_demo.py failed:\n{}{}",
        String::from_utf8_lossy(&built.stdout),
        String::from_utf8_lossy(&built.stderr)
    );

    let srv = StaticServer::start(&site, &[]).unwrap();
    let mut b = Browser::launch(&chrome, &work.join("chrome")).unwrap();
    b.call(
        "Emulation.setFocusEmulationEnabled",
        json::loads(r#"{"enabled": true}"#).unwrap(),
    )
    .unwrap();
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
    let run = format!("/demo/runs/{RUN}/");

    // 1. The skip link on each page, on a wide screen.
    viewport(&mut b, 1280, 800, false);
    for (page, path, main) in [
        ("home", "/demo/".to_string(), "home"),
        ("person", run.clone(), "content"),
        ("finding", format!("{run}finding/{FINDING}/"), "finding"),
        ("family", "/demo/runs/DEMO-FAMILY/".to_string(), "family"),
    ] {
        visit(&mut b, &srv.url(&path));
        let got = skip_link(&mut b);
        let want = format!("[\"A\",\"skip\",\"#{main}\",true] [\"{main}\",true]");
        check(
            &format!(
                "{page}: the first Tab focuses the skip link to #{main}, on screen, and Enter moves the focus into main"
            ),
            got == want,
            got,
        );
    }

    // 2. A finger's targets on the person page, at 320 and 390 px wide.
    for (w, h) in [(320, 800), (390, 844)] {
        viewport(&mut b, w, h, true);
        visit(&mut b, &srv.url(&format!("{run}?{w}")));
        let sizes = text(b.eval(SIZES).ok());
        let Ok(Value::Array(a)) = json::loads(&sizes) else {
            check(&format!("{w}x{h}: the sizes can be read"), false, sizes);
            continue;
        };
        let Value::Array(toggles) = &a[0] else {
            panic!("{sizes}")
        };
        let rows: Vec<(f64, f64, String, String)> = toggles
            .iter()
            .map(|t| match t {
                Value::Array(t) => (
                    num(&t[0]),
                    num(&t[1]),
                    text(Some(t[2].clone())),
                    text(Some(t[3].clone())),
                ),
                _ => (0.0, 0.0, String::new(), String::new()),
            })
            .collect();
        check(
            &format!(
                "{w}x{h}: every caveat's More (summary.toggle, at least 2) has min-height 44 px and is 44 px tall"
            ),
            rows.len() >= 2 && rows.iter().all(|r| r.0 >= TAP && r.1 >= TAP),
            sizes.clone(),
        );
        let mut names: Vec<&String> = rows.iter().map(|r| &r.2).collect();
        names.sort();
        names.dedup();
        check(
            &format!("{w}x{h}: each More is named for its own finding, no two alike"),
            names.len() == rows.len()
                && rows
                    .iter()
                    .all(|r| !r.3.is_empty() && r.2 == format!("More of the caveat: {}", r.3)),
            sizes.clone(),
        );
        check(
            &format!("{w}x{h}: the replay's slider and Play button are at least 44 px tall"),
            num(&a[1]) >= TAP && num(&a[2]) >= TAP,
            format!("slider {} play {}", num(&a[1]), num(&a[2])),
        );
    }

    // 3. A tap on the Cyrius row of the step timeline, on a phone.
    viewport(&mut b, 390, 844, true);
    visit(&mut b, &srv.url(&format!("{run}?tap")));
    let row = text(b.wait(CYRIUS_ROW, Duration::from_secs(10)));
    match json::loads(&row) {
        Ok(Value::Array(r)) if r.len() == 3 => {
            tap(&mut b, num(&r[0]) + 60.0, num(&r[1]) + num(&r[2]) / 2.0);
            let heading = format!("Step {CYRIUS}: Cyrius CYP2D6");
            let opened = format!(
                "(function () {{ var h = document.querySelector('#detail h2'); return h && h.textContent === {}; }}())",
                js_str(&heading)
            );
            b.wait(&opened, Duration::from_secs(5));
            let at = text(b.eval(AFTER_TAP).ok());
            let (said, top, bottom, mast, tabs, sh, ch) = match json::loads(&at) {
                Ok(Value::Array(v)) if v.len() == 7 => (
                    text(Some(v[0].clone())),
                    num(&v[1]),
                    num(&v[2]),
                    num(&v[3]),
                    num(&v[4]),
                    num(&v[5]),
                    num(&v[6]),
                ),
                _ => (at.clone(), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            };
            check(
                "390x844: a tap on the Cyrius row shows its detail between the masthead and the tab bar",
                said == heading && top >= mast && bottom <= tabs,
                at.clone(),
            );
            check(
                "390x844: the step timeline does not scroll inside the page (scrollHeight equals clientHeight)",
                sh > 0.0 && sh == ch,
                at,
            );
        }
        _ => check("390x844: the Cyrius row is on the timeline", false, row),
    }

    let errors = b.page_errors();
    check(
        "no page error on any page",
        errors.is_empty(),
        format!("{errors:?}"),
    );
    drop((b, srv));
    assert!(fails.is_empty(), "failed: {fails:?}");
}
