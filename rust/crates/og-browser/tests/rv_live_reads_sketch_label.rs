//! The reads offers of the live demo, in a headless Chrome, on the demo the
//! docs publish (og live demo, run with the og this test
//! run built).
//!
//! 1. A finding's reads drawing says it is an illustration before the
//!    picture, and no offer of reads in the demo links to the viewer page,
//!    which only says the demo has no reads: the finding page, a card's
//!    links, a chart mark's panel, a search hit with a place and the map's
//!    card say "Reads: not in the demo" instead. Outside the demo the
//!    finding page keeps its link.
//! 2. After a search for CYP2D6, choosing its finding and selecting the
//!    card's reads button, the reads panel's top is inside the viewport, on
//!    a wide screen and on a phone.
//! 3. Red control: the same demo with a map.js that does not scroll to the
//!    panel leaves it out of view, so check 2 can fail.
//!
//! Without Chrome or a built og the test says SKIP, and fails in CI.

use og_browser::{Browser, Scratch, StaticServer, Value, chrome_or_skip, repo_root};
use pycompat::json;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;

const RUN: &str = "DEMO-001";
/// A finding of the demo with a place: its page draws the reads.
const PLACED: &str = "clinvar-demoa-chr2-47403191-c-t";
/// The line map.js scrolls to the reads panel with: check 3 removes it.
const SCROLL: &str = "      reads.scrollIntoView({ block: 'nearest' });\n";

fn text(v: Option<Value>) -> String {
    match v {
        Some(Value::Str(s)) => s,
        Some(v) => format!("{v:?}"),
        None => "None".into(),
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
    let exe = |p: &Path| p.is_file();
    if let Some(o) = std::env::var_os("OG").filter(|o| !o.is_empty()) {
        let p = PathBuf::from(o);
        return exe(&p).then_some(p);
    }
    let me = std::env::current_exe().ok()?;
    let p = me.parent()?.parent()?.join("og");
    exe(&p).then_some(p)
}

fn copy_dir(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for e in fs::read_dir(from).unwrap() {
        let e = e.unwrap();
        let (src, dst) = (e.path(), to.join(e.file_name()));
        if e.file_type().unwrap().is_dir() {
            copy_dir(&src, &dst);
        } else {
            fs::copy(&src, &dst).unwrap();
        }
    }
}

/// A term of the demo's search index with a mention at a place: its hit
/// offers the reads (the first term in sorted order).
fn located_term(demo: &Path) -> Option<String> {
    let path = demo
        .join("runs")
        .join(RUN)
        .join(format!("{RUN}_search_index.json"));
    let index = json::loads(&fs::read_to_string(path).ok()?).ok()?;
    let Value::Object(top) = index else {
        return None;
    };
    let Some(Value::Object(terms)) = top.get("terms") else {
        return None;
    };
    let placed = |m: &Value| matches!(m, Value::Object(m) if matches!(m.get("locus"), Some(Value::Str(s)) if !s.is_empty()));
    let located = |items: &Value| match items {
        Value::Array(items) => items.iter().any(|it| match it {
            Value::Object(it) => {
                matches!(it.get("mentions"), Some(Value::Array(ms)) if ms.iter().any(placed))
            }
            _ => false,
        }),
        _ => false,
    };
    let mut names: Vec<&String> = terms.keys().collect();
    names.sort();
    names.into_iter().find(|t| located(&terms[*t])).cloned()
}

/// Wide screen or phone: the page's own size, as a reader's window.
fn viewport(b: &mut Browser, width: u32, height: u32, mobile: bool) {
    let p = json::loads(&format!(
        r#"{{"width": {width}, "height": {height}, "deviceScaleFactor": 1, "mobile": {mobile}}}"#
    ))
    .unwrap();
    b.call("Emulation.setDeviceMetricsOverride", p).unwrap();
}

/// Every link on the page whose text names reads and that opens the viewer.
const READS_LINKS: &str = "JSON.stringify(Array.prototype.filter.call(document.querySelectorAll('a[href]'), function (a) { \
     return /reads/i.test(a.textContent) && /(^|\\/)viewer\\//.test(a.getAttribute('href')); })\
     .map(function (a) { return a.textContent + ' -> ' + a.getAttribute('href'); }))";

/// The reader's path to the reads of CYP2D6 on the person page at `url`:
/// the map's search, its finding, the card's reads button. The panel's top
/// against the viewport, the scroll and the button's words.
fn open_cyp2d6_reads(b: &mut Browser, url: &str) -> Result<(f64, f64, String), String> {
    let secs = Duration::from_secs;
    b.go(url)?;
    b.wait(
        "(document.getElementById('genome-map') || {getAttribute: function () {}}).getAttribute('data-pgp-map') === 'drawn' && 'ok'",
        secs(30),
    )
    .ok_or("the map was not drawn")?;
    b.eval(
        "var i = document.getElementById('gm-find'); i.scrollIntoView({block: 'center'}); i.focus(); \
         i.value = 'CYP2D6'; i.dispatchEvent(new Event('input')); true",
    )?;
    let li = text(b.wait(
        "var li = document.querySelector('#gm-find-list li[data-kind=finding]'); li && li.textContent",
        secs(10),
    ));
    if li == "None" {
        return Err("the search offers no CYP2D6 finding".into());
    }
    b.eval(
        "document.querySelector('#gm-find-list li[data-kind=finding]')\
         .dispatchEvent(new MouseEvent('mousedown', {bubbles: true, cancelable: true})); true",
    )?;
    let words = text(b.wait(
        "var c = document.querySelector('.gm-card-reads'); c && !c.closest('[hidden]') && c.textContent",
        secs(10),
    ));
    if words == "None" {
        return Err(format!(
            "choosing {li:?} showed no card with a reads button"
        ));
    }
    // The reader sees the button before selecting it.
    b.eval(
        "var c = document.querySelector('.gm-card-reads'); c.scrollIntoView({block: 'nearest'}); \
         window.scrollBefore = window.scrollY; c.click(); true",
    )?;
    let at = text(b.wait(
        "!document.getElementById('gm-reads').hidden && \
         JSON.stringify([document.getElementById('gm-reads').getBoundingClientRect().top, window.innerHeight])",
        secs(5),
    ));
    match json::loads(&at) {
        Ok(Value::Array(a)) if a.len() == 2 => {
            let n = |v: &Value| match v {
                Value::Int(i) => *i as f64,
                Value::Float(f) => *f,
                _ => f64::NAN,
            };
            Ok((n(&a[0]), n(&a[1]), words))
        }
        _ => Err(format!("the reads panel did not open: {at}")),
    }
}

#[test]
fn the_demo_says_it_has_no_reads_and_the_reads_panel_comes_into_view() {
    let Some(chrome) = chrome_or_skip("open the live demo") else {
        return;
    };
    let Some(og) = need(find_og(), "og (OG, or target/<profile>/og)") else {
        return;
    };
    // Declared before the server and the browser: dropped after them, also on a panic.
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "reads-sketch");
    let work = scratch.path();
    let site = work.join("site");
    let demo = site.join("demo");
    let built = Command::new(&og)
        .args(["live", "demo"])
        .current_dir(repo_root())
        .arg("--out")
        .arg(&demo)
        .arg("--work")
        .arg(work.join("root"))
        .output()
        .unwrap();
    assert!(
        built.status.success(),
        "og live demo failed:\n{}{}",
        String::from_utf8_lossy(&built.stdout),
        String::from_utf8_lossy(&built.stderr)
    );
    // Outside the demo: the same finding page in live mode.
    let fdir = demo.join("runs").join(RUN).join("finding");
    let page = fs::read_to_string(fdir.join(PLACED).join("index.html")).unwrap();
    assert!(
        page.contains("data-mode=\"replay\""),
        "the finding page is not a replay"
    );
    fs::create_dir_all(fdir.join("live")).unwrap();
    fs::write(
        fdir.join("live").join("index.html"),
        page.replace("data-mode=\"replay\"", "data-mode=\"live\""),
    )
    .unwrap();
    // The red control: the demo with a map.js that does not scroll to the panel.
    let planted = site.join("planted");
    copy_dir(&demo, &planted);
    let map_js = planted.join("assets").join("map.js");
    let js = fs::read_to_string(&map_js).unwrap();
    let scrolls = js.matches(SCROLL).count();
    fs::write(&map_js, js.replace(SCROLL, "")).unwrap();

    let srv = StaticServer::start(&site, &[]).unwrap();
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
    check(
        "map.js scrolls to the reads panel with the one line the control removes",
        scrolls == 1,
        format!("{scrolls} copies of {SCROLL:?}"),
    );
    let secs = Duration::from_secs;
    viewport(&mut b, 1280, 800, false);

    // 1. The finding page: the words before the drawing, and no link to the viewer page.
    let run = format!("/demo/runs/{RUN}");
    let ready = "document.documentElement.getAttribute('data-pgp-status') === 'ready' && 'ok'";
    b.go(&srv.url(&format!("{run}/finding/{PLACED}/"))).unwrap();
    b.wait(ready, secs(20));
    let before = text(b.eval(
        "(function () { var f = document.querySelector('figure.reads'); if (!f) { return 'no drawing'; } \
         var t = ''; for (var n = f.firstChild; n && n.nodeName.toLowerCase() !== 'svg'; n = n.nextSibling) { t += n.textContent; } \
         return JSON.stringify([t, f.querySelector('figcaption').textContent]); })()",
    ).ok());
    let (label, caption) = match json::loads(&before) {
        Ok(Value::Array(a)) if a.len() == 2 => (text(Some(a[0].clone())), text(Some(a[1].clone()))),
        _ => (before.clone(), String::new()),
    };
    check(
        "browser: the words before a finding's reads drawing say it is an illustration, not the person's reads",
        label.contains("An illustration, not your reads") && !label.starts_with("The reads at"),
        label.clone(),
    );
    check(
        "browser: in the demo the drawing's caption says the demo has no reads",
        caption == "Reads: not in the demo",
        caption,
    );
    let links = text(b.eval(READS_LINKS).ok());
    check(
        "browser: the demo's finding page links no reads to the viewer page",
        links == "[]",
        links,
    );
    // Every finding page of the demo.
    let mut pages = Vec::new();
    for e in fs::read_dir(&fdir).unwrap() {
        let id = e.unwrap().file_name().to_string_lossy().into_owned();
        if id != "live" {
            pages.push(id);
        }
    }
    pages.sort();
    let mut bad = Vec::new();
    for id in &pages {
        b.go(&srv.url(&format!("{run}/finding/{id}/"))).unwrap();
        let ok = b.wait(ready, secs(20)).is_some();
        let l = text(b.eval(READS_LINKS).ok());
        if !ok || l != "[]" {
            bad.push(format!("{id}: ready {ok}, {l}"));
        }
    }
    check(
        &format!(
            "browser: none of the demo's {} finding pages links reads to the viewer page",
            pages.len()
        ),
        pages.len() > 1 && bad.is_empty(),
        format!("{bad:?}"),
    );
    // Outside the demo the finding page keeps its link to the reads.
    b.go(&srv.url(&format!("{run}/finding/live/"))).unwrap();
    b.wait(ready, secs(20));
    let live = text(b.eval(READS_LINKS).ok());
    check(
        "browser: outside the demo the finding page still links the real reads in the genome viewer",
        live.contains("See the real reads in the genome viewer")
            && live.contains("viewer/#chr2%3A47403191"),
        live,
    );

    // The person page: a card's links and a chart mark's panel.
    b.go(&srv.url(&format!("{run}/"))).unwrap();
    let cards = text(b.wait(
        "Array.prototype.filter.call(document.querySelectorAll('.card-links .viewer'), function (s) { \
         return s.textContent === 'Reads: not in the demo'; }).length > 0 && 'ok'",
        secs(30),
    ));
    check(
        "browser: a card with a place says 'Reads: not in the demo' in its links",
        cards == "ok",
        cards,
    );
    let mark = text(b.wait(
        "(function () { var a = document.querySelector('figure[data-chart=genome] a[data-mark]'); if (!a) { return ''; } \
         a.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true})); \
         var d = document.querySelector('.chart-detail'); return d && !d.hidden && d.textContent; })()",
        secs(20),
    ));
    check(
        "browser: a chart mark's panel says 'Reads: not in the demo'",
        mark.contains("Reads: not in the demo"),
        mark.clone(),
    );
    let links = text(b.eval(READS_LINKS).ok());
    check(
        "browser: the demo's person page links no reads to the viewer page",
        links == "[]",
        links,
    );

    // The search box: a hit with a place says the demo has no reads, as the
    // other offers do (the same search interaction as the demo smoke test).
    let term = located_term(&demo).unwrap_or_default();
    let plain = !term.is_empty()
        && term
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || "._-".contains(c));
    check(
        "the demo's search index has a term with a mention at a place",
        plain,
        format!("{term:?}"),
    );
    b.eval(
        "var s = document.getElementById('replay-slider'); \
         if (s) { s.value = s.max; s.dispatchEvent(new Event('input')); } true",
    )
    .unwrap();
    let field = "document.querySelector('input[type=search]:not(#gm-find)')";
    let open = text(b.wait(
        &format!("var i = {field}; i && !i.disabled && 'enabled'"),
        secs(30),
    ));
    let hit = if open == "enabled" && plain {
        b.eval(&format!(
            "var i = {field}; i.value = '{term}'; i.dispatchEvent(new Event('input')); true"
        ))
        .unwrap();
        text(b.wait(
            "+(document.documentElement.dataset.pgpSearch || 0) > 0 && \
             JSON.stringify(Array.prototype.map.call(document.querySelectorAll('.search-results .viewer'), \
             function (v) { return v.textContent; }))",
            secs(20),
        ))
    } else {
        format!("the search box is {open}")
    };
    check(
        &format!("browser: a search hit for {term} with a place says 'Reads: not in the demo'"),
        hit.contains("Reads: not in the demo") && !hit.contains("reads at "),
        hit,
    );
    let links = text(b.eval(READS_LINKS).ok());
    check(
        "browser: the demo's search results link no reads to the viewer page",
        links == "[]",
        links,
    );

    // 2. The map: search, choose the finding, select the card's reads button.
    for (width, height, mobile, name) in [
        (1280, 800, false, "a wide screen"),
        (390, 844, true, "a phone"),
    ] {
        viewport(&mut b, width, height, mobile);
        match open_cyp2d6_reads(&mut b, &srv.url(&format!("{run}/"))) {
            Ok((top, h, words)) => {
                check(
                    &format!(
                        "browser: on {name} the map card's reads button says 'Reads: not in the demo'"
                    ),
                    words == "Reads: not in the demo",
                    words,
                );
                check(
                    &format!(
                        "browser: on {name} the reads panel of CYP2D6's card opens inside the viewport"
                    ),
                    top >= 0.0 && top < h,
                    format!("panel top {top}, viewport height {h}"),
                );
            }
            Err(e) => check(
                &format!("browser: on {name} the reads of CYP2D6 open from the map"),
                false,
                e,
            ),
        }
    }
    // 3. The red control: without the scroll the same path leaves the panel out of view.
    viewport(&mut b, 1280, 800, false);
    match open_cyp2d6_reads(&mut b, &srv.url(&format!("/planted/runs/{RUN}/"))) {
        Ok((top, h, _)) => check(
            "browser: control: a map.js without the scroll leaves the reads panel out of view (check 2 can fail)",
            !(top >= 0.0 && top < h),
            format!("panel top {top}, viewport height {h}"),
        ),
        Err(e) => check(
            "browser: control: the planted demo opens the reads",
            false,
            e,
        ),
    }
    let errors = b.page_errors();
    check(
        "browser: no page error",
        errors.is_empty(),
        format!("{errors:?}"),
    );
    assert!(fails.is_empty(), "failed: {fails:?}");
}
