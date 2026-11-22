//! The HTML report of the demo people DEMO-001, DEMO-002 and DEMO-003 in a
//! headless Chrome, dark colour scheme, as a reader on a phone or a keyboard
//! meets it. The reports are written by og itself: `og dev demo-sample`, then
//! `og report render` on the sample folder, as step 24 runs it. Each sample
//! folder gets a `run_manifest.tsv` with a 40-character pipeline commit, as
//! a real run writes one, so the header shows the commit.
//!
//! - At 320, 390 and 720 px wide the page is no wider than the screen
//!   (`document.scrollingElement.scrollWidth`), and no card is wider than the
//!   screen; the CPIC (#cpic) and raw-data export (#consumer-export) cards
//!   are there. On main the page was 1195 px at 320 and 390 px, 1171 px at
//!   720 px (1197 and 1173 px on DEMO-003): unbroken paths and wide tables.
//!   The same holds at 320 px for DEMO-001's folder rendered under a sample
//!   name of 57 characters with no place to break: the header, outside every
//!   card, held the commit and the name at their full width.
//! - Every table sits in its own scroll region (`.tscroll`), and no two
//!   regions share a name (a card with several tables names each one).
//! - On a phone the contents are folded: the page opens on the results. At
//!   1440 px they show beside the results, with no summary line to open.
//! - The first Tab reaches "Skip to results"; Enter on it moves the focus to
//!   the results (`main#main`).
//! - The page logs no error (its Content-Security-Policy blocks nothing it needs).
//!
//! The og binary is the `OG` variable, else the `og` beside this test's own
//! target folder (cargo test --workspace builds it). Without og or Chrome the
//! test says SKIP, and fails in CI. Set `OG_BROWSER_SHOTS` to a folder to
//! save a screenshot of each view there. The PR records the failing run on
//! main's report.

use og_browser::{Browser, Scratch, StaticServer, Value, chrome_or_skip, js_str, repo_root};
use pycompat::json;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;

const PEOPLE: [(&str, &[&str]); 3] = [
    ("DEMO-001", &[]),
    ("DEMO-002", &["--sex", "female", "--seed", "2"]),
    ("DEMO-003", &["--seed", "3"]),
];
/// A sample name with no place for a line to break, for DEMO-001's folder.
const LONG_SAMPLE: &str = "DEMO_0001_A_SAMPLE_NAME_WITH_NO_PLACE_TO_BREAK_0123456789";
/// `run_manifest.tsv` rows with a full 40-character pipeline commit.
const MANIFEST: &str =
    "run\twritten_by\trun-all.sh\nrun\tgit_commit\t0123456789abcdef0123456789abcdef01234567\n";
const WIDTHS: [u32; 3] = [320, 390, 720];
const HEIGHT: u32 = 844;
const DESKTOP: u32 = 1440;

/// The og binary, or `None` after a SKIP line (a panic in CI).
fn og_or_skip() -> Option<PathBuf> {
    let from_env = std::env::var_os("OG")
        .filter(|v| !v.is_empty())
        .map(PathBuf::from);
    let beside = std::env::current_exe()
        .ok()
        .and_then(|e| Some(e.parent()?.parent()?.join("og")));
    let found = from_env.or(beside).filter(|p| p.is_file());
    if found.is_none() {
        assert!(
            std::env::var("CI").as_deref() != Ok("true"),
            "[FAIL] no og binary to write the demo reports: CI must run this test (set OG)"
        );
        println!("[SKIP] no og binary to write the demo reports (set OG)");
    }
    found
}

fn og(bin: &Path, args: &[&str]) {
    let out = Command::new(bin)
        .args(args)
        .current_dir(repo_root())
        .output()
        .expect("og runs");
    assert!(
        out.status.success(),
        "og {args:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

fn render(bin: &Path, sample: &str, dir: &Path, out: &Path, assets: &Path) {
    og(
        bin,
        &[
            "report",
            "render",
            "--sample",
            sample,
            "--sample-dir",
            dir.to_str().unwrap(),
            "--assets",
            assets.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ],
    );
}

/// Writes each demo person into `root` with a run manifest and renders
/// `root/<sample>/report.html`; then DEMO-001's folder under [`LONG_SAMPLE`]
/// into `root/long/report.html`.
fn demo_reports(bin: &Path, root: &Path) {
    let repo = repo_root();
    let (r, repo, assets) = (
        root.to_str().unwrap(),
        repo.to_str().unwrap(),
        repo.join("assets"),
    );
    for (sample, extra) in PEOPLE {
        let mut args = vec![
            "dev",
            "demo-sample",
            "--genome-dir",
            r,
            "--records",
            "3000",
            "--repo",
            repo,
            "--sample",
            sample,
        ];
        args.extend(extra);
        og(bin, &args);
        let d = root.join(sample);
        fs::write(d.join("run_manifest.tsv"), MANIFEST).unwrap();
        render(bin, sample, &d, &d.join("report.html"), &assets);
    }
    let long = root.join("long");
    fs::create_dir_all(&long).unwrap();
    render(
        bin,
        LONG_SAMPLE,
        &root.join(PEOPLE[0].0),
        &long.join("report.html"),
        &assets,
    );
    for page in PEOPLE
        .iter()
        .map(|(s, _)| root.join(s))
        .chain([long])
        .map(|d| d.join("report.html"))
    {
        let html = fs::read_to_string(&page).unwrap();
        assert!(
            html.contains("Pipeline commit: 0123456789abcdef0123456789abcdef01234567"),
            "{}: the header has no 40-character commit",
            page.display()
        );
    }
}

/// The page at `width` x `height` CSS pixels, dark; a phone below 900 px.
fn view(b: &mut Browser, width: u32, height: u32) {
    let mobile = width < 900;
    b.call(
        "Emulation.setDeviceMetricsOverride",
        json::loads(&format!(
            r#"{{"width": {width}, "height": {height}, "deviceScaleFactor": 1, "mobile": {mobile}}}"#
        ))
        .unwrap(),
    )
    .unwrap();
    b.call(
        "Emulation.setEmulatedMedia",
        json::loads(r#"{"features":[{"name":"prefers-color-scheme","value":"dark"}]}"#).unwrap(),
    )
    .unwrap();
}

fn open(b: &mut Browser, url: &str) {
    b.go(url).unwrap();
    let ready = format!(
        "document.readyState === 'complete' && location.href === {}",
        js_str(url)
    );
    assert!(
        b.wait(&ready, Duration::from_secs(30)).is_some(),
        "{url} did not load"
    );
}

/// The JSON object the expression returns as a string.
fn probe(b: &mut Browser, expr: &str) -> Value {
    let raw = match b.eval(expr) {
        Ok(Value::Str(s)) => s,
        other => panic!("{expr}: {other:?}"),
    };
    json::loads(&raw).unwrap_or_else(|e| panic!("{expr}: {raw}: {e:?}"))
}

fn get<'a>(v: &'a Value, key: &str) -> &'a Value {
    match v {
        Value::Object(m) => m.get(key).unwrap_or(&Value::Null),
        _ => &Value::Null,
    }
}

fn num(v: &Value) -> f64 {
    match v {
        Value::Int(n) => *n as f64,
        Value::Float(f) => *f,
        other => panic!("not a number: {other:?}"),
    }
}

/// What the page measures: its scroll width, the cards wider than the
/// screen, whether #cpic and #consumer-export are there, the tables outside
/// a scroll region, the region names given to more than one region, whether
/// the contents are open and shown.
const MEASURE: &str = "(function(){var vw=document.documentElement.clientWidth;\
     var seen={},dups=[];[].slice.call(document.querySelectorAll('[role=region]')).forEach(function(r){\
       var n=r.getAttribute('aria-label');if(seen[n]&&dups.indexOf(n)<0)dups.push(n);seen[n]=1;});\
     var wide=[].slice.call(document.querySelectorAll('.card')).filter(function(c){\
       var r=c.getBoundingClientRect();return r.width>vw+0.5||r.right>vw+0.5;})\
       .map(function(c){var s=c.closest('section');return (s&&s.id||'?')+':'+Math.round(c.getBoundingClientRect().width);});\
     var bare=[].slice.call(document.querySelectorAll('table')).filter(function(t){\
       var p=t.parentElement;return !(p&&p.classList.contains('tscroll')&&p.getAttribute('role')==='region'\
       &&p.getAttribute('aria-label'));}).length;\
     var fold=document.querySelector('details.toc-fold');var link=document.querySelector('nav.toc a');\
     var sum=fold&&fold.querySelector('summary');\
     return JSON.stringify({scroll:document.scrollingElement.scrollWidth,vw:vw,wide:wide,\
       cpic:!!document.getElementById('cpic'),consumer:!!document.getElementById('consumer-export'),\
       tables:document.querySelectorAll('table').length,bare:bare,dups:dups,fold:!!fold,open:!!(fold&&fold.open),\
       contents:!!(link&&link.checkVisibility()),summary:!!(sum&&sum.checkVisibility())});})()";

/// Presses one key in the page.
fn key(b: &mut Browser, name: &str, code: u32) {
    for kind in ["keyDown", "keyUp"] {
        b.call(
            "Input.dispatchKeyEvent",
            json::loads(&format!(
                r#"{{"type": "{kind}", "key": "{name}", "code": "{name}", "windowsVirtualKeyCode": {code}}}"#
            ))
            .unwrap(),
        )
        .unwrap();
    }
}

/// Saves a screenshot as `<OG_BROWSER_SHOTS>/report-<name>.png` when that
/// variable names a folder.
fn shot(b: &mut Browser, name: &str) {
    let Some(dir) = std::env::var_os("OG_BROWSER_SHOTS").filter(|d| !d.is_empty()) else {
        return;
    };
    let r = b
        .call(
            "Page.captureScreenshot",
            json::loads(r#"{"format": "png"}"#).unwrap(),
        )
        .unwrap();
    let Value::Str(data) = get(&r, "data") else {
        panic!("no screenshot data")
    };
    fs::write(
        Path::new(&dir).join(format!("report-{name}.png")),
        base64(data),
    )
    .unwrap();
}

fn base64(text: &str) -> Vec<u8> {
    let val = |c: u8| match c {
        b'A'..=b'Z' => c - b'A',
        b'a'..=b'z' => c - b'a' + 26,
        b'0'..=b'9' => c - b'0' + 52,
        b'+' => 62,
        _ => 63,
    };
    let bytes: Vec<u8> = text.bytes().filter(|c| *c != b'=').map(val).collect();
    bytes
        .chunks(4)
        .flat_map(|q| {
            let n = q
                .iter()
                .enumerate()
                .fold(0u32, |n, (i, v)| n | u32::from(*v) << (18 - 6 * i));
            n.to_be_bytes()[1..q.len()].to_vec()
        })
        .collect()
}

#[test]
fn the_demo_reports_fit_a_phone_and_skip_to_the_results() {
    let Some(chrome) = chrome_or_skip("load the demo reports") else {
        return;
    };
    let Some(bin) = og_or_skip() else { return };
    // Declared before the server and the browser: dropped after them, also on a panic.
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "report-frame");
    let runs = scratch.path().join("runs");
    fs::create_dir_all(&runs).unwrap();
    demo_reports(&bin, &runs);
    let srv = StaticServer::start(&runs, &[]).unwrap();
    let mut b = Browser::launch(&chrome, &scratch.path().join("chrome")).unwrap();
    let mut problems = Vec::new();

    for (sample, _) in PEOPLE {
        let url = srv.url(&format!("/{sample}/report.html"));
        for width in WIDTHS {
            view(&mut b, width, HEIGHT);
            open(&mut b, &format!("{url}?w={width}"));
            let m = probe(&mut b, MEASURE);
            let at = format!("{sample} at {width} px");
            let scroll = num(get(&m, "scroll"));
            if scroll > f64::from(width) {
                problems.push(format!("{at}: the page is {scroll} px wide"));
            }
            if !matches!(get(&m, "wide"), Value::Array(a) if a.is_empty()) {
                problems.push(format!(
                    "{at}: cards wider than the screen: {:?}",
                    get(&m, "wide")
                ));
            }
            for (k, sid) in [("cpic", "#cpic"), ("consumer", "#consumer-export")] {
                if !matches!(get(&m, k), Value::Bool(true)) {
                    problems.push(format!("{at}: no {sid} card"));
                }
            }
            if num(get(&m, "tables")) < 10.0 || num(get(&m, "bare")) != 0.0 {
                problems.push(format!(
                    "{at}: {} tables, {} outside a scroll region",
                    num(get(&m, "tables")),
                    num(get(&m, "bare"))
                ));
            }
            if !matches!(get(&m, "dups"), Value::Array(a) if a.is_empty()) {
                problems.push(format!(
                    "{at}: scroll regions that share a name: {:?}",
                    get(&m, "dups")
                ));
            }
            if !matches!(
                (
                    get(&m, "fold"),
                    get(&m, "open"),
                    get(&m, "contents"),
                    get(&m, "summary")
                ),
                (
                    Value::Bool(true),
                    Value::Bool(false),
                    Value::Bool(false),
                    Value::Bool(true)
                )
            ) {
                problems.push(format!(
                    "{at}: the contents are not folded behind their summary: {m:?}"
                ));
            }
            if sample == "DEMO-001" {
                shot(&mut b, &format!("{sample}-{width}"));
            }
            // The keyboard: the first Tab reaches the skip link, which shows;
            // Enter moves the focus to the results.
            if width == 390 {
                key(&mut b, "Tab", 9);
                let f = probe(
                    &mut b,
                    "(function(){var a=document.activeElement;var r=a.getBoundingClientRect();\
                     return JSON.stringify({text:a.textContent,href:a.getAttribute('href'),top:r.top,h:r.height});})()",
                );
                if !matches!(get(&f, "text"), Value::Str(t) if t == "Skip to results")
                    || !matches!(get(&f, "href"), Value::Str(h) if h == "#main")
                    || num(get(&f, "top")) < 0.0
                    || num(get(&f, "h")) <= 0.0
                {
                    problems.push(format!("{at}: the first Tab reaches {f:?}"));
                }
                if sample == "DEMO-001" {
                    shot(&mut b, &format!("{sample}-{width}-skip"));
                }
                key(&mut b, "Enter", 13);
                let moved = b.wait(
                    "document.activeElement && document.activeElement.id === 'main' && location.hash === '#main'",
                    Duration::from_secs(5),
                );
                if moved.is_none() {
                    let now = b.eval(
                        "document.activeElement.outerHTML.slice(0, 80) + ' ' + location.hash",
                    );
                    problems.push(format!(
                        "{at}: Enter on the skip link left the focus on {now:?}"
                    ));
                }
            }
        }
        view(&mut b, DESKTOP, 900);
        open(&mut b, &format!("{url}?w={DESKTOP}"));
        let m = probe(&mut b, MEASURE);
        if !matches!(
            (get(&m, "contents"), get(&m, "summary")),
            (Value::Bool(true), Value::Bool(false))
        ) {
            problems.push(format!(
                "{sample} at {DESKTOP} px: the contents do not show beside the results: {m:?}"
            ));
        }
        if num(get(&m, "scroll")) > f64::from(DESKTOP) {
            problems.push(format!(
                "{sample} at {DESKTOP} px: the page is wider than the screen: {m:?}"
            ));
        }
        if sample == "DEMO-001" {
            shot(&mut b, &format!("{sample}-{DESKTOP}"));
        }
        let errors = b.page_errors();
        if !errors.is_empty() {
            problems.push(format!("{sample}: page errors: {errors:?}"));
        }
    }
    // The header sits outside every card: a long sample name and the
    // 40-character commit must wrap there too.
    let width = WIDTHS[0];
    view(&mut b, width, HEIGHT);
    open(&mut b, &srv.url(&format!("/long/report.html?w={width}")));
    let m = probe(&mut b, MEASURE);
    if num(get(&m, "scroll")) > f64::from(width)
        || !matches!(get(&m, "wide"), Value::Array(a) if a.is_empty())
    {
        problems.push(format!(
            "{LONG_SAMPLE} at {width} px: the page is wider than the screen: {m:?}"
        ));
    }
    drop(b);
    drop(srv);
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}
