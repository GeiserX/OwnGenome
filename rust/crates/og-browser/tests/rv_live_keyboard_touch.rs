//! The live pages for a keyboard and for a finger, in a headless Chrome, on
//! the demo the docs publish (tests/demo/build_live_demo.py, run with the og
//! this test run built).
//!
//! 1. On the home, person, finding and family pages the first Tab focuses
//!    the skip link, which then shows and points at the page's main; Enter
//!    moves the focus into main, and main's top is not under the masthead.
//! 2. On DEMO-001 at 320x800 and at 390x844, every caveat's More
//!    (`summary.toggle`) has a min-height of at least 44 px, is that tall
//!    and takes a tap over 44 px of its height; so does every card's
//!    caveat summary (`.caveat-full > summary`) once its card is open. No
//!    tap on the control just before or just after one of them (a link, a
//!    fold's summary) lands on it. The replay's slider
//!    and its Play button are 44 px tall. Each More is named for its own
//!    finding, so no two share a name.
//! 3. At 390x844, a tap on the Cyrius row of the step timeline brings the
//!    step's detail into view, clear of the masthead and the tab bar, and
//!    the timeline's box does not scroll inside the page.
//! 4. Two findings with one title (two hits in one gene) get two different
//!    More names: the test adds a second finding with the title of one of
//!    DEMO-002's to the copy of the demo it serves, since the demo has no
//!    such pair. With two detail lines, each More names its own; with no
//!    detail line on either, or the same one on both, each names its place
//!    ("(1 of 2)", "(2 of 2)").
//!
//! Each check was red on main before this change, or on an earlier head of
//! it for the neighbours' check (the PR says where).
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

/// python3 when it is 3.10 or newer, read from `python3 --version`
/// ("Python 3.12.1"; a Python 2 writes it to stderr).
fn find_python() -> Option<PathBuf> {
    let out = Command::new("python3").arg("--version").output().ok()?;
    if !out.status.success() {
        return None;
    }
    let said = format!(
        "{}{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    let mut part = said
        .trim()
        .strip_prefix("Python ")?
        .split('.')
        .map(|p| p.parse::<u32>().ok());
    let (major, minor) = (part.next()??, part.next()??);
    ((major, minor) >= (3, 10)).then(|| PathBuf::from("python3"))
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
/// screen] [main's id, whether the focus is inside main after Enter with
/// main's top at or below the masthead's bottom]. When that never holds,
/// the second part is [main's id, focus inside main, main's top, the
/// masthead's bottom].
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
    let after = match b.wait(
        "(function () { var m = document.querySelector('main'), h = document.querySelector('header.mast'); \
         var a = document.activeElement; \
         return m && h && m.contains(a) && m.getBoundingClientRect().top >= h.getBoundingClientRect().bottom \
           && JSON.stringify([m.id, true]); }())",
        Duration::from_secs(3),
    ) {
        Some(v) => text(Some(v)),
        None => text(
            b.eval(
                "(function () { var m = document.querySelector('main'), h = document.querySelector('header.mast'); \
                 return JSON.stringify([m ? m.id : '', !!m && m.contains(document.activeElement), \
                   m ? m.getBoundingClientRect().top : null, h ? h.getBoundingClientRect().bottom : null]); }())",
            )
            .ok(),
        ),
    };
    format!("{first} {after}")
}

/// A page function `hit(s)`: how many CSS px of `s`'s height take a tap.
/// It brings `s` to the middle of the screen, then asks elementFromPoint at
/// each pixel row of its box, 5 px in from its left edge, whether the tap
/// lands inside `s`.
macro_rules! hit {
    () => {
        "var hit = function (s) { s.scrollIntoView({ block: 'center' }); var r = s.getBoundingClientRect(), n = 0; \
         for (var y = r.top + 0.5; y < r.bottom; y += 1) { var e = document.elementFromPoint(r.left + 5, y); \
           if (e && s.contains(e)) { n += 1; } } return n; }; "
    };
}

/// A page function `beside(s)`: whether `s` (a fold's summary grown to
/// 44 px) takes taps meant for the controls next to it. For the nearest
/// visible control (a link, a summary, a button or a form field) before
/// `s`'s fold and the nearest one after it, it brings the control to the
/// middle of the screen and asks elementFromPoint at each pixel row of the
/// control's first box, from 0.1 px below its top to 0.1 px above its
/// bottom, at its left edge, its middle and its right edge (2 px in).
/// Returns [[what it is, rows, rows where a tap lands on `s`], ...], one
/// entry per neighbour there is.
macro_rules! beside {
    () => {
        "var beside = function (s) { var d = s.parentElement, out = []; \
           var all = Array.prototype.filter.call(document.querySelectorAll('a[href], summary, button, input, select, textarea'), \
             function (c) { return !d.contains(c) && !c.contains(d) && c.checkVisibility() && c.getBoundingClientRect().height > 0; }); \
           var before = all.filter(function (c) { return d.compareDocumentPosition(c) & 2; }).pop(); \
           var after = all.filter(function (c) { return d.compareDocumentPosition(c) & 4; })[0]; \
           [before, after].forEach(function (e) { if (!e) { return; } \
             e.scrollIntoView({ block: 'center' }); var r = e.getClientRects()[0], ys = [], rows = 0, taken = 0; \
             for (var y = r.top + 0.1; y < r.bottom - 0.1; y += 1) { ys.push(y); } ys.push(r.bottom - 0.1); \
             ys.forEach(function (y) { rows += 1; \
               if ([r.left + 2, (r.left + r.right) / 2, r.right - 2].some(function (x) { \
                 var h = document.elementFromPoint(x, y); return h && s.contains(h); })) { taken += 1; } }); \
             out.push([e.tagName + ' ' + (e.textContent || '').trim().slice(0, 40), rows, taken]); }); \
           return out; }; "
    };
}

/// Read `beside`'s answer: (what, rows, rows taken) per neighbour.
fn neighbours(v: &Value) -> Vec<(String, f64, f64)> {
    match v {
        Value::Array(n) => n
            .iter()
            .map(|e| match e {
                Value::Array(e) if e.len() == 3 => {
                    (text(Some(e[0].clone())), num(&e[1]), num(&e[2]))
                }
                _ => (format!("{e:?}"), 0.0, f64::NAN),
            })
            .collect(),
        _ => vec![(format!("{v:?}"), 0.0, f64::NAN)],
    }
}

/// The neighbours that lose part of their tap: (what, rows, rows taken).
fn covered(all: &[Vec<(String, f64, f64)>]) -> Vec<(String, f64, f64)> {
    all.iter()
        .flatten()
        .filter(|n| n.1 <= 0.0 || n.2 != 0.0)
        .cloned()
        .collect()
}

/// The sizes a finger meets on the person page: every `summary.toggle`
/// [its height, its min-height, its name, its finding's title, the px of
/// it that take a tap, its neighbours (`beside`)], the slider and the Play
/// button.
const SIZES: &str = concat!(
    "(function () { ",
    hit!(),
    beside!(),
    "var t = Array.prototype.map.call(document.querySelectorAll('summary.toggle'), function (s) { \
       var row = s.closest('.row'), h = row && row.querySelector('h3'); \
       return [s.getBoundingClientRect().height, parseFloat(getComputedStyle(s).minHeight) || 0, \
         s.getAttribute('aria-label') || s.textContent, h ? h.textContent : '', hit(s), beside(s)]; }); \
     return JSON.stringify([t, document.getElementById('replay-slider').getBoundingClientRect().height, \
       document.getElementById('replay-play').getBoundingClientRect().height]); }())"
);

/// Open the folds around the first card caveat, then every card caveat
/// summary that shows: [its height, its min-height, the px of it that take
/// a tap, its neighbours (`beside`)].
const CAVEATS: &str = concat!(
    "(function () { ",
    hit!(),
    beside!(),
    "var c = document.querySelector('.caveat-full'); if (!c) { return '[]'; } \
     for (var d = c.parentElement; d; d = d.parentElement) { if (d.tagName === 'DETAILS') { d.open = true; } } \
     var s = Array.prototype.filter.call(document.querySelectorAll('.caveat-full > summary'), function (s) { \
       return s.checkVisibility(); }); \
     return JSON.stringify(s.map(function (s) { \
       return [s.getBoundingClientRect().height, parseFloat(getComputedStyle(s).minHeight) || 0, hit(s), beside(s)]; })); }())"
);

/// The More names of a finding and of its twin (the finding the test adds
/// with the same title): [the finding's, the twin's]; '' when one is not
/// on the page.
fn twin_names(id: &str, twin: &str) -> String {
    format!(
        "(function () {{ var name = function (id) {{ var f = document.getElementById('f-' + id); \
           var s = f && f.closest('.row') && f.closest('.row').querySelector('summary.toggle'); \
           return s ? s.getAttribute('aria-label') || '' : ''; }}; \
         return JSON.stringify([name({}), name({})]); }}())",
        js_str(id),
        js_str(twin)
    )
}

/// How the test's twin differs from the finding it copies.
#[derive(Clone, Copy, Debug)]
enum Twin {
    /// Its own detail line: each More names its finding's detail.
    Distinct,
    /// No detail line on either (null on the finding, absent on the twin).
    Missing,
    /// The finding's own detail line on both.
    Repeated,
}

/// Write DEMO-002's results.json (`raw`), in the served copy, with a second
/// finding with the title of one of its findings, as two hits in one gene
/// would: the first finding with a caveat and a detail line in a section the
/// page draws as full articles (section 1, or one of sections 2 to 5 with at
/// most two findings, so at most three with the twin). The twin gets its own
/// id and sits next to it in that section, last in the findings list; its
/// detail line is as `how` says. Returns the finding's and the twin's ids
/// and the More names they must get.
fn add_twin(results: &Path, raw: &str, how: Twin) -> Option<([String; 2], [String; 2])> {
    let mut v = json::loads(raw).ok()?;
    let Value::Object(top) = &mut v else {
        return None;
    };
    let field = |f: &Value, k: &str| match f {
        Value::Object(m) => match m.get(k) {
            Some(Value::Str(s)) => s.clone(),
            _ => String::new(),
        },
        _ => String::new(),
    };
    let ids = |s: &Value| match s {
        Value::Object(m) => match m.get("findings") {
            Some(Value::Array(a)) => a.iter().map(|i| text(Some(i.clone()))).collect(),
            _ => Vec::new(),
        },
        _ => Vec::<String>::new(),
    };
    let n = |s: &Value| match s {
        Value::Object(m) => m.get("n").map(num).unwrap_or(f64::NAN),
        _ => f64::NAN,
    };
    let Some(Value::Array(fs)) = top.get("findings") else {
        return None;
    };
    let usable = |id: &String| {
        fs.iter().find(|f| field(f, "id") == *id).filter(|f| {
            matches!(f, Value::Object(m) if matches!(m.get("caveat"), Some(Value::Array(c)) if !c.is_empty()))
                && !field(f, "title").is_empty()
                && !field(f, "detail").is_empty()
        })
    };
    let Some(Value::Array(journal)) = top.get("journal") else {
        return None;
    };
    let (sec, first) = journal.iter().enumerate().find_map(|(i, s)| {
        let (k, list) = (n(s), ids(s));
        let drawn = k == 1.0 || ((2.0..=5.0).contains(&k) && list.len() <= 2);
        drawn
            .then(|| list.iter().find_map(|id| usable(id).cloned()))
            .flatten()
            .map(|f| (i, f))
    })?;
    let (id, title, detail) = (
        field(&first, "id"),
        field(&first, "title"),
        field(&first, "detail"),
    );
    let twin_id = format!("{id}-twin");
    let mut twin = first;
    if let Value::Object(m) = &mut twin {
        m.insert("id".into(), Value::Str(twin_id.clone()));
        match how {
            Twin::Distinct => {
                m.insert("detail".into(), Value::Str(format!("{detail}, second hit")));
            }
            Twin::Missing => {
                m.shift_remove("detail");
            }
            Twin::Repeated => {}
        }
    }
    if let Some(Value::Array(fs)) = top.get_mut("findings") {
        if let Twin::Missing = how
            && let Some(Value::Object(m)) = fs.iter_mut().find(|f| field(f, "id") == id)
        {
            m.insert("detail".into(), Value::Null);
        }
        fs.push(twin);
    }
    if let Some(Value::Array(journal)) = top.get_mut("journal")
        && let Some(Value::Object(s)) = journal.get_mut(sec)
        && let Some(Value::Array(list)) = s.get_mut("findings")
    {
        let at = list
            .iter()
            .position(|i| *i == Value::Str(id.clone()))
            .map_or(list.len(), |p| p + 1);
        list.insert(at, Value::Str(twin_id.clone()));
    }
    std::fs::write(results, json::dumps(&v, &json::Dump::default()).ok()?).ok()?;
    let more = |tail: &str| format!("More of the caveat: {title} ({tail})");
    let names = match how {
        Twin::Distinct => [more(&detail), more(&format!("{detail}, second hit"))],
        Twin::Missing | Twin::Repeated => [more("1 of 2"), more("2 of 2")],
    };
    Some(([id, twin_id], names))
}

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

    let results = site.join("demo/runs/DEMO-002/results.json");
    let raw = std::fs::read_to_string(&results).unwrap();
    let twin = add_twin(&results, &raw, Twin::Distinct);

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
                "{page}: the first Tab focuses the skip link to #{main}, on screen, and Enter moves the focus into main, below the masthead"
            ),
            got == want,
            got,
        );
    }

    // 2. A finger's targets on the person page, at 320 and 390 px wide.
    for (w, h) in [(320, 800), (390, 844)] {
        viewport(&mut b, w, h, true);
        visit(&mut b, &srv.url(&format!("{run}?{w}")));
        // 'ready' holds once the run is drawn, before results.json is: wait for the findings.
        b.wait(
            "document.querySelectorAll('summary.toggle').length >= 2",
            Duration::from_secs(20),
        );
        let sizes = text(b.eval(SIZES).ok());
        let Ok(Value::Array(a)) = json::loads(&sizes) else {
            check(&format!("{w}x{h}: the sizes can be read"), false, sizes);
            continue;
        };
        let Value::Array(toggles) = &a[0] else {
            panic!("{sizes}")
        };
        type Toggle = (f64, f64, String, String, f64, Vec<(String, f64, f64)>);
        let rows: Vec<Toggle> = toggles
            .iter()
            .map(|t| match t {
                Value::Array(t) if t.len() == 6 => (
                    num(&t[0]),
                    num(&t[1]),
                    text(Some(t[2].clone())),
                    text(Some(t[3].clone())),
                    num(&t[4]),
                    neighbours(&t[5]),
                ),
                _ => (0.0, 0.0, String::new(), String::new(), 0.0, Vec::new()),
            })
            .collect();
        check(
            &format!(
                "{w}x{h}: every caveat's More (summary.toggle, at least 2) has min-height 44 px and is 44 px tall"
            ),
            rows.len() >= 2 && rows.iter().all(|r| r.0 >= TAP && r.1 >= TAP),
            sizes.clone(),
        );
        check(
            &format!("{w}x{h}: a tap lands on each More over at least 44 px of its height"),
            rows.len() >= 2 && rows.iter().all(|r| r.4 >= TAP),
            format!("{:?}", rows.iter().map(|r| r.4).collect::<Vec<_>>()),
        );
        let near: Vec<_> = rows.iter().map(|r| r.5.clone()).collect();
        check(
            &format!(
                "{w}x{h}: no tap on the control just before or just after a More lands on the More"
            ),
            rows.len() >= 2 && near.iter().any(|n| !n.is_empty()) && covered(&near).is_empty(),
            format!("covered {:?} of {near:?}", covered(&near)),
        );
        let mut names: Vec<&String> = rows.iter().map(|r| &r.2).collect();
        names.sort();
        names.dedup();
        check(
            &format!("{w}x{h}: each More is named for its own finding, no two alike"),
            rows.len() >= 2
                && names.len() == rows.len()
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
        let caveats = text(b.eval(CAVEATS).ok());
        type Caveat = (f64, f64, f64, Vec<(String, f64, f64)>);
        let cav: Vec<Caveat> = match json::loads(&caveats) {
            Ok(Value::Array(v)) => v
                .iter()
                .map(|c| match c {
                    Value::Array(c) if c.len() == 4 => {
                        (num(&c[0]), num(&c[1]), num(&c[2]), neighbours(&c[3]))
                    }
                    _ => (0.0, 0.0, 0.0, Vec::new()),
                })
                .collect(),
            _ => Vec::new(),
        };
        check(
            &format!(
                "{w}x{h}: with a card open, every caveat summary (.caveat-full > summary) has min-height 44 px, \
                 is 44 px tall and takes a tap over 44 px"
            ),
            !cav.is_empty() && cav.iter().all(|c| c.0 >= TAP && c.1 >= TAP && c.2 >= TAP),
            caveats.clone(),
        );
        let near: Vec<_> = cav.iter().map(|c| c.3.clone()).collect();
        check(
            &format!(
                "{w}x{h}: no tap on the control just before or just after a caveat summary (such as Open in \
                 the report) lands on the caveat summary"
            ),
            near.iter().any(|n| !n.is_empty()) && covered(&near).is_empty(),
            format!("covered {:?} of {near:?}", covered(&near)),
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

    // 4. Two findings with one title: each More gets a name of its own, with
    //    two detail lines, with none and with the same one on both.
    if twin.is_none() {
        check(
            "DEMO-002 has a finding with a caveat and a detail line to give a twin",
            false,
            String::new(),
        );
    }
    for how in [Twin::Distinct, Twin::Missing, Twin::Repeated] {
        let Some(([id, twin_id], want)) = add_twin(&results, &raw, how) else {
            continue;
        };
        visit(&mut b, &srv.url(&format!("/demo/runs/DEMO-002/?{how:?}")));
        b.wait(
            &format!(
                "!!document.getElementById({})",
                js_str(&format!("f-{twin_id}"))
            ),
            Duration::from_secs(20),
        );
        let got = text(b.eval(&twin_names(&id, &twin_id)).ok());
        let names: Option<Vec<String>> = match json::loads(&got) {
            Ok(Value::Array(v)) => Some(v.into_iter().map(|n| text(Some(n))).collect()),
            _ => None,
        };
        check(
            &format!(
                "two findings with one title ({how:?}): each More gets a name of its own, not empty"
            ),
            names.as_deref() == Some(&want[..]) && want[0] != want[1],
            format!("{got}, want {want:?}"),
        );
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
