//! The live demo's genome map (Figure 1, assets/live/map.js), its In view table, the person page's
//! layout and phone folds, and the live site's genome viewer, in a headless Chrome. Moved here from
//! tests/test_live_demo.py (map_checks, inview_checks, layout_checks, phone_report_checks,
//! fold_checks) and tests/test_live_viewer.py, check for check.
//!
//! Each check runs on the demo as `og live demo` builds it, then once more with its planted mistake
//! (a line of map.js, person.js or live.js changed in the built site, a map.json or results.json
//! changed, or a viewer answer changed by a proxy in front of `og live serve`). The check must fail
//! there: that is its control, printed as `[RED]`. Without Chrome or og every test says SKIP, and
//! fails in CI.

use og_browser::{Browser, Scratch, StaticServer, Value, chrome_or_skip, js_str, repo_root};
use pycompat::json::{Dump, dumps, loads};
use std::collections::{BTreeSet, HashMap};
use std::fs;
use std::io::{BufRead, BufReader, Read, Write};
use std::net::{Shutdown, TcpListener, TcpStream};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};

const S: &str = "DEMO-001";
const PARTNER: &str = "DEMO-002";
const SON: &str = "DEMO-003";

// ---------------------------------------------------------------------------------------------
// Values, files and the report of checks.

static NULL: Value = Value::Null;

fn get<'a>(v: &'a Value, k: &str) -> &'a Value {
    match v {
        Value::Object(o) => o.get(k).unwrap_or(&NULL),
        _ => &NULL,
    }
}

fn idx(v: &Value, i: usize) -> &Value {
    match v {
        Value::Array(a) => a.get(i).unwrap_or(&NULL),
        _ => &NULL,
    }
}

fn items(v: &Value) -> &[Value] {
    match v {
        Value::Array(a) => a,
        _ => &[],
    }
}

fn st(v: &Value) -> String {
    match v {
        Value::Str(s) => s.clone(),
        Value::Null => String::new(),
        Value::Bool(b) => b.to_string(),
        Value::Int(i) => i.to_string(),
        Value::Float(f) => f.to_string(),
        _ => show(v),
    }
}

fn num(v: &Value) -> f64 {
    match v {
        Value::Int(i) => *i as f64,
        Value::Float(f) => *f,
        Value::Str(s) => s.parse().unwrap_or(f64::NAN),
        _ => f64::NAN,
    }
}

fn strs(v: &Value) -> Vec<String> {
    items(v).iter().map(st).collect()
}

fn truthy(v: &Value) -> bool {
    match v {
        Value::Null => false,
        Value::Bool(b) => *b,
        Value::Int(i) => *i != 0,
        Value::Float(f) => *f != 0.0,
        Value::Str(s) => !s.is_empty(),
        Value::Array(a) => !a.is_empty(),
        Value::Object(o) => !o.is_empty(),
    }
}

fn show(v: &Value) -> String {
    let t = dumps(v, &Dump::default()).unwrap_or_default();
    if t.chars().count() > 400 {
        format!("{}…", t.chars().take(400).collect::<String>())
    } else {
        t
    }
}

fn obj(v: &mut Value) -> &mut indexmap::IndexMap<String, Value> {
    match v {
        Value::Object(o) => o,
        _ => panic!("expected an object"),
    }
}

fn arr(v: &mut Value) -> &mut Vec<Value> {
    match v {
        Value::Array(a) => a,
        _ => panic!("expected an array"),
    }
}

fn at<'a>(v: &'a mut Value, path: &[&str]) -> &'a mut Value {
    path.iter().fold(v, |v, k| {
        obj(v)
            .get_mut(*k)
            .unwrap_or_else(|| panic!("no {k} in the document"))
    })
}

fn s(text: &str) -> Value {
    Value::Str(text.to_owned())
}

fn read_map(path: &Path) -> Value {
    let mut text = String::new();
    flate2::read::GzDecoder::new(fs::File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    loads(&text).unwrap()
}

fn write_map(path: &Path, doc: &Value) {
    let text = dumps(
        doc,
        &Dump {
            separators: Some((",", ":")),
            ..Dump::default()
        },
    )
    .unwrap();
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(text.as_bytes()).unwrap();
    fs::write(path, gz.finish().unwrap()).unwrap();
}

/// A file's bytes, written back when dropped (also when a check panics).
struct Keep(PathBuf, Vec<u8>);

impl Keep {
    fn new(path: &Path) -> Keep {
        Keep(path.to_owned(), fs::read(path).unwrap())
    }
    fn text(&self) -> String {
        String::from_utf8(self.1.clone()).unwrap()
    }
}

impl Drop for Keep {
    fn drop(&mut self) {
        let _ = fs::write(&self.0, &self.1);
    }
}

/// `path` with each (from, to, n) replaced, where `from` must be in the file exactly `n` times:
/// a planted mistake, taken out again when the returned value is dropped.
fn plant_n(path: &Path, pairs: &[(&str, &str, usize)]) -> Keep {
    let keep = Keep::new(path);
    let mut text = keep.text();
    for (from, to, n) in pairs {
        assert_eq!(
            keep.text().matches(from).count(),
            *n,
            "the planted line is not {n} times in {}: {from}",
            path.display()
        );
        text = text.replace(from, to);
    }
    fs::write(path, text).unwrap();
    keep
}

fn plant(path: &Path, pairs: &[(&str, &str)]) -> Keep {
    let v: Vec<_> = pairs.iter().map(|(a, z)| (*a, *z, 1)).collect();
    plant_n(path, &v)
}

/// The checks of one test: each prints `[PASS]` or `[FAIL]`, each control `[RED]` (the check failed
/// on its planted mistake, as it must) or `[FAIL]`; the test fails at the end if any did.
#[derive(Default)]
struct Report {
    fails: Vec<String>,
}

impl Report {
    fn check(&mut self, desc: &str, bad: Vec<String>) {
        if bad.is_empty() {
            println!("[PASS] {desc}");
        } else {
            println!("[FAIL] {desc} -- {}", bad.join(" | "));
            self.fails.push(format!("{desc}: {}", bad.join(" | ")));
        }
    }

    fn ok(&mut self, desc: &str, ok: bool, detail: impl std::fmt::Debug) {
        let bad = if ok {
            vec![]
        } else {
            vec![format!("{detail:?}")]
        };
        self.check(desc, bad);
    }

    /// The control `desc`: the check's problems with the mistake planted, which must not be empty.
    fn control(&mut self, desc: &str, bad: Vec<String>) {
        if bad.is_empty() {
            println!("[FAIL] control: {desc}: the check still passes");
            self.fails.push(format!("control passed: {desc}"));
        } else {
            let first = bad[0].chars().take(160).collect::<String>();
            println!("[RED] control: {desc} -- the check fails: {first}");
        }
    }

    fn done(self) {
        assert!(self.fails.is_empty(), "failed:\n{}", self.fails.join("\n"));
    }
}

fn sleep(ms: u64) {
    thread::sleep(Duration::from_millis(ms));
}

fn og() -> Option<PathBuf> {
    let path = std::env::var_os("OG").map(PathBuf::from).or_else(|| {
        let test = std::env::current_exe().ok()?;
        Some(test.parent()?.parent()?.join("og"))
    });
    if path.as_ref().is_some_and(|p| p.is_file()) {
        return path;
    }
    assert_ne!(
        std::env::var("CI").as_deref(),
        Ok("true"),
        "build og or set OG in CI"
    );
    println!("[SKIP] build og or set OG");
    None
}

/// `n` written as Python's `f"{n:,}"`.
fn commas(n: usize) -> String {
    let d = n.to_string();
    let mut out = String::new();
    for (i, c) in d.chars().enumerate() {
        if i > 0 && (d.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

fn count_words(n: usize) -> String {
    if n == 1 {
        "1 feature".into()
    } else {
        format!("{} features", commas(n))
    }
}

/// The text of `pub const NAME: &str = "...";` in a Rust source file of the repository, its
/// continuations joined and its `\"` unescaped: the constant the binary under test uses.
fn source_const(file: &str, name: &str) -> String {
    let src = fs::read_to_string(repo_root().join(file)).unwrap();
    let head = format!("pub const {name}: &str = \"");
    let start = src
        .find(&head)
        .unwrap_or_else(|| panic!("{file} sets no {name}"))
        + head.len();
    let raw = &src[start..start + src[start..].find("\";").unwrap()];
    let joined: String = raw
        .split("\\\n")
        .enumerate()
        .map(|(i, p)| if i == 0 { p } else { p.trim_start() })
        .collect();
    joined.replace("\\\"", "\"")
}

// ---------------------------------------------------------------------------------------------
// The demo site in a Chrome.

struct Demo {
    b: Browser,
    _server: StaticServer,
    base: String,
    site: PathBuf,
    work: PathBuf,
    og: PathBuf,
    r: Report,
    _scratch: Scratch,
}

impl Demo {
    /// `og live demo` built into a scratch folder, served from its root with no cache, and a
    /// Chrome at 1440 x 900 that caches nothing (the controls change files under it).
    fn new(label: &str) -> Option<Demo> {
        let chrome = chrome_or_skip(label)?;
        let og = og()?;
        let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), label);
        let (site, work) = (scratch.path().join("site"), scratch.path().join("work"));
        let built = Command::new(&og)
            .current_dir(repo_root())
            .args(["live", "demo", "--records", "2000", "--out"])
            .arg(&site)
            .arg("--work")
            .arg(&work)
            .output()
            .unwrap();
        assert!(
            built.status.success(),
            "og live demo: {}",
            String::from_utf8_lossy(&built.stderr)
        );
        let server = StaticServer::start(&site, &[("Cache-Control", "no-store")]).unwrap();
        let base = server.url("").trim_end_matches('/').to_owned();
        let b = Browser::launch(&chrome, &scratch.path().join("chrome")).unwrap();
        let mut d = Demo {
            b,
            _server: server,
            base,
            site,
            work,
            og,
            r: Report::default(),
            _scratch: scratch,
        };
        d.call("Page.enable", "{}");
        d.call("Network.enable", "{}");
        d.call("Network.setCacheDisabled", r#"{"cacheDisabled":true}"#);
        d.viewport(1440, 900, 1);
        Some(d)
    }

    fn call(&mut self, method: &str, params: &str) -> Value {
        self.b.call(method, loads(params).unwrap()).unwrap()
    }

    fn viewport(&mut self, w: u32, h: u32, scale: u32) {
        self.call(
            "Emulation.setDeviceMetricsOverride",
            &format!(r#"{{"width":{w},"height":{h},"deviceScaleFactor":{scale},"mobile":false}}"#),
        );
    }

    /// The value of `js`; an exception gives null, as the Python driver did.
    fn ev(&mut self, js: &str) -> Value {
        self.b.eval(js).unwrap_or_else(|e| {
            println!(
                "  (script error: {})",
                e.chars().take(200).collect::<String>()
            );
            Value::Null
        })
    }

    /// The value of `js`, which returns JSON text, parsed.
    fn jv(&mut self, js: &str) -> Value {
        match self.ev(js) {
            Value::Str(t) => loads(&t).unwrap_or(Value::Null),
            v => v,
        }
    }

    /// The value of the promise `js` returns, parsed when it is JSON text.
    fn awaited(&mut self, js: &str) -> Value {
        let params = format!(
            r#"{{"expression":{},"awaitPromise":true,"returnByValue":true}}"#,
            js_str(js)
        );
        let r = self.call("Runtime.evaluate", &params);
        match get(get(&r, "result"), "value") {
            Value::Str(t) => loads(t).unwrap_or(Value::Null),
            v => v.clone(),
        }
    }

    fn wait(&mut self, js: &str, secs: u64) -> bool {
        self.b.wait(js, Duration::from_secs(secs)).is_some()
    }

    fn go(&mut self, path: &str) {
        self.b.go("about:blank").unwrap();
        self.wait("location.href === 'about:blank'", 20);
        self.b.go(&format!("{}{path}", self.base)).unwrap();
    }

    fn map_json(&self) -> PathBuf {
        self.site.join(format!("runs/{S}/map.json"))
    }
}

// ---------------------------------------------------------------------------------------------
// The genome map.

/// What the map shows: its state, the findings it draws (its list of marks), its timings, level and
/// place, its chips, its margin's words and every text and title it holds, and the address.
const MAP_JS: &str = r#"JSON.stringify((function () {
  var f = document.getElementById('genome-map');
  if (!f) { return { state: 'absent', marks: [] }; }
  var t = document.querySelector('.gm-tip');
  return { state: f.dataset.pgpMap || '', level: f.dataset.level || '', chrom: f.dataset.chrom || '', view: f.dataset.view || '',
    parse: Number(f.dataset.pgpMapParseMs), paint: Number(f.dataset.pgpMapPaintMs),
    marks: Array.prototype.map.call(f.querySelectorAll('.gm-marks a'), function (a) { return { id: a.dataset.id,
      level: a.dataset.level, chrom: a.dataset.chrom, s: Number(a.dataset.s), e: Number(a.dataset.e) }; }),
    chips: Array.prototype.map.call(f.querySelectorAll('.gm-chip'), function (c) { return c.textContent; }),
    title: (document.getElementById('gm-title') || {}).textContent || '',
    layers: Array.prototype.map.call(f.querySelectorAll('.gm-layer'), function (l) { return [l.dataset.layer, l.dataset.state,
      l.querySelector('dd').textContent]; }),
    toast: (f.querySelector('.gm-toast') || { hidden: true }).hidden ? '' : f.querySelector('.gm-toast').textContent,
    tip: t && !t.hidden ? t.textContent : '',
    text: [f.textContent].concat(Array.prototype.map.call(f.querySelectorAll('[href],[title],[aria-label]'), function (e) {
      return (e.getAttribute('href') || '') + ' ' + (e.getAttribute('title') || '') + ' ' + (e.getAttribute('aria-label') || ''); })).join(' '),
    hash: location.hash };
}()))"#;

const MAP_READY: &str = "document.documentElement.dataset.pgpResults && (function () { var f = document.getElementById('genome-map'); \
    return f && /^(drawn|failed|nodecompress)$/.test(f.dataset.pgpMap || '') \
    && (f.dataset.pgpMap !== 'drawn' || !!f.dataset.level); }())";

/// The card, the Nearest button, the genes drawn and the search box's lines, read at once.
const NEAR_JS: &str = r#"JSON.stringify((function () {
  var f = document.getElementById('genome-map'), n = f.querySelector('.gm-nearest'), c = f.querySelector('.gm-card');
  return { near: n.hidden ? '' : n.textContent, card: c.hidden ? '' : c.textContent, genes: f.dataset.genes || '',
    crumb: f.querySelector('.gm-crumb').textContent,
    found: Array.prototype.map.call(f.querySelectorAll('.gm-find-list li'), function (li) {
      return { kind: li.dataset.kind || '', text: li.textContent }; }) };
}()))"#;

/// A second figure on DEMO-001's page with a veil that keeps nothing closed: drawn from map.json,
/// then reloaded from another file.
const SECOND_JS: &str = r#"(function () {
  var o = { url: 'map.json', veil: function () { return 'open'; }, now: function () { return true; },
    known: function () { return true; }, stepDone: function () { return true; }, open: function () {},
    close: function () {}, fallback: function () { return []; } };
  var g = window.PGPMap.create(o);
  g.node.id = 'gm-second';
  document.body.appendChild(g.node);
  window.pgpSecond = { g: g, o: o };
  g.refresh('k1');
  return true;
}())"#;

const SECOND_FACTS: &str = r#"JSON.stringify((function () {
  var f = document.getElementById('gm-second');
  var text = function (q) { var n = f.querySelector(q); return n ? n.textContent : ''; }, tip = f.querySelector('.gm-tip');
  return { state: f.dataset.pgpMap, marks: f.querySelectorAll('.gm-marks a').length,
    narr: text('.gm-narration'), layers: text('.gm-layers'), tip: !tip || tip.hidden ? '' : tip.textContent };
}()))"#;

const ROWS_JS: &str = r#"JSON.stringify(Array.prototype.map.call(
  document.querySelectorAll('#gm-second .gm-inview tbody tr[data-layer=findings] .gm-go'),
  function (x) { return x.textContent; }))"#;

/// The closer levels: (address, level, view the address keeps, words the title starts or ends
/// with, genes drawn).
const CLOSER: [(&str, &str, &str, &str, usize); 3] = [
    (
        "#chr11:4.5-6.5M",
        "region",
        "chr11:4.500-6.500M",
        "Chromosome 11, band 11p15",
        10,
    ),
    (
        "#chr11:5.150-5.350M",
        "gene",
        "chr11:5.150-5.350M",
        "and its neighbours",
        10,
    ),
    (
        "#chr11:5225001-5235000",
        "reads",
        "chr11:5225001-5235000",
        "and its neighbours",
        1,
    ),
];

/// (run, id, chrom, mb, level) of each mark the karyogram drew.
fn karyogram_rows() -> Vec<(String, String, String, f64, String)> {
    fs::read_to_string(repo_root().join("tests/fixtures/map/karyogram-marks.tsv"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.starts_with("run\t") && !l.trim().is_empty())
        .map(|l| {
            let p: Vec<_> = l.split('\t').collect();
            (
                p[0].into(),
                p[1].into(),
                format!("chr{}", p[2]),
                p[3].parse().unwrap(),
                p[4].into(),
            )
        })
        .collect()
}

/// The first finding of DEMO-001's map.json that a gate keeps closed: (index, the findings).
fn gated(doc: &Value) -> (usize, Value) {
    let f = get(get(get(doc, "per"), S), "findings").clone();
    let i = items(get(&f, "gate"))
        .iter()
        .position(truthy)
        .expect("DEMO-001 has no gated finding");
    (i, f)
}

impl Demo {
    /// A fresh load of `run`'s page at `hash`: the map's facts once it settled.
    fn map_page(&mut self, run: &str, hash: &str) -> Value {
        self.go(&format!("/runs/{run}/{hash}"));
        if self.wait(MAP_READY, 20) {
            sleep(600); // the margin follows the camera 140 ms after it settles
        }
        self.jv(MAP_JS)
    }

    fn map_now(&mut self) -> Value {
        sleep(1500); // the camera eases for about 0.6 s, then the margin follows 140 ms after it settles
        self.jv(MAP_JS)
    }

    fn near(&mut self) -> Value {
        self.jv(NEAR_JS)
    }

    fn map_key(&mut self, key: &str, canvas: bool) {
        let target = if canvas {
            "document.querySelector('.gm-canvas')"
        } else {
            "document.body"
        };
        self.ev(&format!(
            "{target}.dispatchEvent(new KeyboardEvent('keydown', {{key: {}, bubbles: true}})); true",
            js_str(key)
        ));
    }

    fn slider(&mut self, to: &str) {
        self.ev(&format!(
            "var s = document.getElementById('replay-slider'); s.value = {to}; s.dispatchEvent(new Event('input')); true"
        ));
    }

    fn rect(&mut self, selector: &str) -> (f64, f64, f64, f64) {
        let r = self.jv(&format!(
            "JSON.stringify(document.querySelector({}).getBoundingClientRect())",
            js_str(selector)
        ));
        (
            num(get(&r, "left")),
            num(get(&r, "top")),
            num(get(&r, "width")),
            num(get(&r, "height")),
        )
    }

    fn mouse(&mut self, kind: &str, x: f64, y: f64, more: &str) {
        self.call(
            "Input.dispatchMouseEvent",
            &format!(r#"{{"type":"{kind}","x":{x},"y":{y}{more}}}"#),
        );
    }

    /// A click at (x, y) of the map's canvas, in CSS px from its corner.
    fn click_canvas(&mut self, x: f64, y: f64) {
        self.ev("document.querySelector('.gm-canvas').scrollIntoView({block: 'center'}); true");
        sleep(300);
        let (l, t, _, _) = self.rect(".gm-canvas");
        let (px, py) = (l + x, t + y);
        self.mouse("mouseMoved", px, py, "");
        for kind in ["mousePressed", "mouseReleased"] {
            self.mouse(kind, px, py, r#","button":"left","clickCount":1"#);
        }
        sleep(600);
    }

    fn drag(&mut self, x0: f64, y: f64, x1: f64) {
        self.mouse("mouseMoved", x0, y, "");
        self.mouse(
            "mousePressed",
            x0,
            y,
            r#","button":"left","buttons":1,"clickCount":1"#,
        );
        for i in 1..9 {
            let x = x0 + (x1 - x0) * f64::from(i) / 8.0;
            self.mouse("mouseMoved", x, y, r#","button":"left","buttons":1"#);
            sleep(20);
        }
        self.mouse(
            "mouseReleased",
            x1,
            y,
            r#","button":"left","buttons":0,"clickCount":1"#,
        );
    }

    /// The address shortly after an action: a click or a brush on the strip pushes its place at
    /// once, a pan replaces it 250 ms after the last move.
    fn hash_now(&mut self, ms: u64) -> String {
        sleep(ms);
        st(&self.ev("location.hash"))
    }

    /// The karyogram's marks the map does not place with the same level and chromosome within 0.01 Mb.
    fn parity_problems(&mut self, rows: &[(String, String, String, f64, String)]) -> Vec<String> {
        let mut bad = vec![];
        let runs: BTreeSet<&str> = rows
            .iter()
            .map(|r| r.0.as_str())
            .filter(|r| *r != SON)
            .collect();
        for run in runs {
            let got = self.map_page(run, "");
            if st(get(&got, "state")) != "drawn" {
                bad.push(format!("{run}: map {}", st(get(&got, "state"))));
                continue;
            }
            let by: HashMap<String, &Value> = items(get(&got, "marks"))
                .iter()
                .map(|m| (st(get(m, "id")), m))
                .collect();
            for (_, fid, c, mb, lv) in rows.iter().filter(|r| r.0 == run) {
                match by.get(fid) {
                    None => bad.push(format!("{run}: {fid} not on the map")),
                    Some(g) => {
                        let at = num(get(g, "s")) / 1e6;
                        if st(get(g, "level")) != *lv
                            || st(get(g, "chrom")) != *c
                            || !at.is_finite()
                            || (at - mb).abs() > 0.01
                        {
                            bad.push(format!(
                                "{run}: {fid} at {} {at:.3} Mb {}, want {c} {mb} {lv}",
                                st(get(g, "chrom")),
                                st(get(g, "level"))
                            ));
                        }
                    }
                }
            }
        }
        bad
    }

    /// The genome and chromosome levels, the links that open each, a click on a row, Esc, and a
    /// place that cannot be read.
    fn state_problems(&mut self) -> Vec<String> {
        let mut bad = vec![];
        let g = self.map_page(S, "");
        if !(st(get(&g, "level")) == "genome"
            && st(get(&g, "title")).contains("all 24 chromosomes")
            && st(get(&g, "view")) == "genome")
        {
            bad.push(format!("genome: {}", show(&g)));
        }
        self.click_canvas(200.0, f64::from(9 + 2 * 22 + 11));
        let c = self.map_now();
        if !(st(get(&c, "level")) == "chromosome"
            && st(get(&c, "chrom")) == "chr3"
            && st(get(&c, "hash")) == "#chr3")
        {
            bad.push(format!(
                "a click on chromosome 3's row: {} {} {}",
                st(get(&c, "level")),
                st(get(&c, "chrom")),
                st(get(&c, "hash"))
            ));
        }
        self.map_key("Escape", true);
        let c = self.map_now();
        if st(get(&c, "level")) != "genome" || st(get(&c, "hash")) != "#genome" {
            bad.push(format!(
                "Esc: {} {}",
                st(get(&c, "level")),
                st(get(&c, "hash"))
            ));
        }
        for (link, chrom, view) in [
            ("#chr16", "chr16", "chr16"),
            ("#chr16:10.0-40.0M", "chr16", "chr16:10.0-40.0M"),
            ("#view=chromosome&chr=11", "chr11", "chr11"),
            ("#chr2:20000001-60000000", "chr2", "chr2:20.0-60.0M"),
        ] {
            let c = self.map_page(S, link);
            if !(st(get(&c, "level")) == "chromosome"
                && st(get(&c, "chrom")) == chrom
                && st(get(&c, "view")) == view)
            {
                bad.push(format!(
                    "{link}: {} {} {}",
                    st(get(&c, "level")),
                    st(get(&c, "chrom")),
                    st(get(&c, "view"))
                ));
            } else if items(get(&c, "marks"))
                .iter()
                .any(|m| st(get(m, "chrom")) != chrom)
            {
                bad.push(format!("{link}: marks of another chromosome listed"));
            }
        }
        self.ev("location.hash = '#genome'; true");
        let c = self.map_now();
        if st(get(&c, "level")) != "genome" {
            bad.push(format!("#genome after a place: {}", st(get(&c, "level"))));
        }
        self.ev("location.hash = '#chrQ'; true");
        let c = self.map_now();
        if !st(get(&c, "toast")).contains("Could not read the place")
            || st(get(&c, "level")) != "genome"
        {
            bad.push(format!(
                "#chrQ: toast {:?}, level {}",
                st(get(&c, "toast")),
                st(get(&c, "level"))
            ));
        }
        bad
    }

    /// DEMO-001 at the first frame of its replay, before the traits step: no APOE chip, switch or mark.
    fn replay_apoe_problems(&mut self) -> Vec<String> {
        self.map_page(S, "");
        self.slider("0");
        let c = self.map_now();
        let mut bad = vec![];
        if truthy(&self.ev("!!document.querySelector('details[data-gate=apoe]')")) {
            bad.push(
                "the traits card has its APOE fold at the first frame: the check tests nothing"
                    .into(),
            );
        }
        if strs(get(&c, "chips")).iter().any(|t| t.contains("APOE"))
            || items(get(&c, "marks"))
                .iter()
                .any(|m| st(get(m, "id")) == "apoe")
        {
            bad.push(format!(
                "before the traits step: chips {:?}",
                strs(get(&c, "chips"))
            ));
        }
        bad
    }

    /// A map whose ROH layer is missing says so in words, never a blank.
    fn missing_layer_problems(&mut self) -> Vec<String> {
        let path = self.map_json();
        let _keep = Keep::new(&path);
        let mut doc = read_map(&path);
        let roh = at(&mut doc, &["layers", "roh"]);
        obj(roh).insert("state".into(), s("missing"));
        obj(roh).insert("count".into(), Value::Int(0));
        obj(roh).insert("note".into(), s("planted: no ROH text"));
        obj(at(&mut doc, &["per", S])).shift_remove("roh");
        write_map(&path, &doc);
        let c = self.map_page(S, "#chr1");
        let row: Vec<_> = items(get(&c, "layers"))
            .iter()
            .filter(|x| st(idx(x, 0)) == "roh")
            .collect();
        match row.first() {
            Some(r)
                if st(idx(r, 1)) == "missing"
                    && st(idx(r, 2)).contains("Not on the map: planted: no ROH text") =>
            {
                vec![]
            }
            _ => vec![format!("the missing ROH layer reads {row:?}")],
        }
    }

    /// A replay before the runs-of-homozygosity step ended says nothing about them on an open
    /// chromosome; at its end it does.
    fn runs_sentence_problems(&mut self) -> Vec<String> {
        let path = self.map_json();
        let _keep = Keep::new(&path);
        let mut doc = read_map(&path);
        let text = obj(&mut doc)
            .entry("text".into())
            .or_insert_with(|| Value::Object(Default::default()));
        let chrom = obj(text)
            .entry("chrom".into())
            .or_insert_with(|| Value::Object(Default::default()));
        obj(chrom).insert(
            "chr1".into(),
            loads(r#"[{"t": "PLANTEDRUNS sentence.", "b": false}]"#).unwrap(),
        );
        write_map(&path, &doc);
        self.map_page(S, "#chr1");
        let mut seen = vec![];
        for to in ["0", "s.max"] {
            self.slider(to);
            let c = self.map_now();
            seen.push((
                st(get(&c, "text")).contains("PLANTEDRUNS"),
                format!("{} {}", st(get(&c, "level")), st(get(&c, "chrom"))),
            ));
        }
        if seen[0].1 != "chromosome chr1" {
            return vec![format!(
                "the first frame is at {}: the check tests nothing",
                seen[0].1
            )];
        }
        let mut bad = vec![];
        if seen[0].0 {
            bad.push("the first frame says the runs sentence".into());
        }
        if !seen[1].0 {
            bad.push("the last frame lacks it".into());
        }
        bad
    }

    /// A figure whose reload fails draws nothing more of the document it had.
    fn failed_reload_problems(&mut self) -> Vec<String> {
        self.map_page(S, "");
        self.ev(SECOND_JS);
        if !self.wait(
            "document.getElementById('gm-second').dataset.pgpMap === 'drawn'",
            20,
        ) {
            return vec!["the second figure never drew".into()];
        }
        sleep(600);
        let drawn = self.jv(SECOND_FACTS);
        if !truthy(get(&drawn, "marks")) || !truthy(get(&drawn, "narr")) {
            return vec![format!(
                "the second figure drew no mark or sentence: the check tests nothing ({})",
                show(&drawn)
            )];
        }
        self.ev("pgpSecond.o.url = 'no-such-map.json'; pgpSecond.g.refresh('k2'); true");
        if !self.wait(
            "document.getElementById('gm-second').dataset.pgpMap === 'failed'",
            20,
        ) {
            return vec!["the second figure's reload did not fail".into()];
        }
        sleep(1000);
        let c = self.jv(SECOND_FACTS);
        if st(get(&c, "state")) != "failed" {
            return vec![format!(
                "after the failed reload the figure could not be read: {}",
                show(&c)
            )];
        }
        if ["marks", "narr", "layers", "tip"]
            .iter()
            .any(|k| truthy(get(&c, k)))
        {
            vec![format!("after the failed reload: {}", show(&c))]
        } else {
            vec![]
        }
    }

    /// A browser with no DecompressionStream: the sentence and the findings with a place.
    fn fallback_problems(&mut self) -> Vec<String> {
        let added = self.call(
            "Page.addScriptToEvaluateOnNewDocument",
            r#"{"source":"delete window.DecompressionStream;"}"#,
        );
        let c = self.map_page(S, "");
        let rows = self.ev("document.querySelectorAll('.gm-fallback table tbody tr').length");
        self.call(
            "Page.removeScriptToEvaluateOnNewDocument",
            &format!(
                r#"{{"identifier":{}}}"#,
                js_str(&st(get(&added, "identifier")))
            ),
        );
        let mut bad = vec![];
        if st(get(&c, "state")) != "nodecompress"
            || !st(get(&c, "text")).contains("This browser cannot open the map file")
        {
            bad.push(format!("no DecompressionStream: {}", st(get(&c, "state"))));
        }
        if !truthy(&rows) {
            bad.push("no DecompressionStream: no findings table".into());
        }
        bad
    }

    /// The region, gene and reads levels from their addresses; Esc steps out one level at a time; the
    /// gauge's Region and Gene stops land on their level; + and → zoom and pan.
    fn closer_problems(&mut self) -> Vec<String> {
        let mut bad = vec![];
        for (link, lv, view, words, genes) in CLOSER {
            self.map_page(S, link);
            sleep(800); // the chromosome's gene file arrives after the first frame
            let c = self.map_now();
            let n = self.near();
            if st(get(&c, "level")) != lv || st(get(&c, "view")) != view {
                bad.push(format!(
                    "{link}: level {}, view {}, want {lv} {view}",
                    st(get(&c, "level")),
                    st(get(&c, "view"))
                ));
            }
            let title = st(get(&c, "title"));
            if !(title.starts_with(words) || title.ends_with(words)) {
                bad.push(format!("{link}: title {title:?}"));
            }
            let g = st(get(&n, "genes"));
            if g.parse::<usize>().map_or(true, |g| g < genes) {
                bad.push(format!("{link}: {g:?} genes drawn, want {genes} or more"));
            }
            let crumb = st(get(&n, "crumb"));
            if !crumb.contains("11p15") || !crumb.starts_with("Genome › chr11") {
                bad.push(format!("{link}: path {crumb:?}"));
            }
        }
        self.map_page(S, "#chr11:5.150-5.350M");
        let mut steps = vec![];
        for _ in 0..3 {
            self.map_key("Escape", true);
            let c = self.map_now();
            steps.push(st(get(&c, "level")));
        }
        if steps != ["region", "chromosome", "genome"] {
            bad.push(format!("Esc from the gene level steps through {steps:?}"));
        }
        for stop in ["region", "gene"] {
            self.map_page(S, "");
            self.ev(&format!(
                "document.querySelector('.gm-gauge [data-stop={stop}]').click(); true"
            ));
            let c = self.map_now();
            if st(get(&c, "level")) != stop {
                bad.push(format!(
                    "the {stop} stop from the genome: level {}",
                    st(get(&c, "level"))
                ));
            }
        }
        self.map_page(S, "#chr11:100.0-120.0M");
        self.map_key("+", true);
        let z = self.map_now();
        self.map_key("ArrowRight", true);
        let r = self.map_now();
        if st(get(&z, "view")) != "chr11:105.0-115.0M"
            || st(get(&r, "view")) != "chr11:106.0-116.0M"
        {
            bad.push(format!(
                "+ then → from chr11:100.0-120.0M: {}, then {}",
                st(get(&z, "view")),
                st(get(&r, "view"))
            ));
        }
        bad
    }

    /// `q` typed in the map's search box: its lines once they settled.
    fn map_find(&mut self, q: &str) -> Vec<Value> {
        self.ev(&format!(
            "var i = document.getElementById('gm-find'); i.focus(); i.value = {}; i.dispatchEvent(new Event('input')); true",
            js_str(q)
        ));
        sleep(800);
        self.ev("document.getElementById('gm-find').dispatchEvent(new Event('input')); true");
        sleep(200);
        items(get(&self.near(), "found")).to_vec()
    }

    fn enter_find(&mut self) {
        self.ev("document.getElementById('gm-find').dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', bubbles: true})); true");
    }

    /// The search box: / puts the cursor in it; a gene, a band and a finding are found and Enter
    /// flies there (a finding with its card); APOE, closed, is never a line.
    fn find_problems(&mut self) -> Vec<String> {
        let mut bad = vec![];
        self.map_page(S, "");
        self.map_key("/", false);
        if st(&self.ev("document.activeElement && document.activeElement.id")) != "gm-find" {
            bad.push("/ does not put the cursor in the map's search box".into());
        }
        for (q, kind, words, chrom, lv) in [
            (
                "HBB",
                "gene",
                "HBB · gene · chr11",
                "chr11",
                ["gene", "reads"],
            ),
            (
                "11p15.4",
                "band",
                "11p15.4 · band",
                "chr11",
                ["chromosome", "region"],
            ),
            ("DEMOC", "finding", "DEMOC · ", "chr16", ["gene", "reads"]),
        ] {
            self.map_page(S, "");
            let found = self.map_find(q);
            let first = found.first().cloned().unwrap_or(Value::Null);
            if st(get(&first, "kind")) != kind || !st(get(&first, "text")).starts_with(words) {
                bad.push(format!("{q}: lines {}", show(&Value::Array(found))));
                continue;
            }
            self.enter_find();
            let c = self.map_now();
            let n = self.near();
            if st(get(&c, "chrom")) != chrom || !lv.contains(&st(get(&c, "level")).as_str()) {
                bad.push(format!(
                    "{q}: Enter flies to {} at {}, want {chrom} at {lv:?}",
                    st(get(&c, "chrom")),
                    st(get(&c, "level"))
                ));
            }
            if kind == "finding" && !st(get(&n, "card")).contains("DEMOC") {
                bad.push(format!(
                    "{q}: Enter opens no card for it: {:?}",
                    st(get(&n, "card"))
                ));
            }
        }
        let doc = read_map(&self.map_json());
        let (i, f) = gated(&doc);
        let (fid, gene, title) = (
            st(idx(get(&f, "id"), i)),
            st(idx(get(&f, "gene"), i)),
            st(idx(get(&f, "title"), i)),
        );
        self.map_page(S, "");
        for q in [&gene, &fid, &title] {
            let found = self.map_find(q);
            let named: Vec<_> = found
                .iter()
                .filter(|x| {
                    let t = st(get(x, "text"));
                    st(get(x, "kind")) == "finding"
                        && (t.contains(&gene) || t.contains(&fid) || t.contains(&title))
                })
                .collect();
            if !named.is_empty() {
                bad.push(format!("{gene}, closed, is a line for {q:?}: {named:?}"));
            }
        }
        let found = self.map_find("DEMOC");
        if !found.iter().any(|x| st(get(x, "kind")) == "finding") {
            bad.push("the open finding DEMOC is not found".into());
        }
        bad
    }

    /// Nearest names the closest finding (else the first on the next chromosome) and flies there with
    /// its card; Esc closes the card; it never leads to APOE, closed.
    fn nearest_problems(&mut self) -> Vec<String> {
        let mut bad = vec![];
        for (link, words, chrom) in [
            (
                "#chr11:100.0-101.0M",
                "Nearest: DEMOB at 108.3 Mb ›",
                "chr11",
            ),
            ("#chr3", "Nearest: DEMOB on chr11 ›", "chr11"),
        ] {
            self.map_page(S, link);
            self.wait(
                "!document.querySelector('#genome-map .gm-nearest').hidden",
                5,
            );
            let n = self.near();
            if st(get(&n, "near")) != words {
                bad.push(format!(
                    "{link}: Nearest reads {:?}, want {words:?}",
                    st(get(&n, "near"))
                ));
                continue;
            }
            self.ev("document.querySelector('.gm-nearest').click(); true");
            let demob = |c: &Value| {
                items(get(c, "marks"))
                    .iter()
                    .any(|m| st(get(m, "id")).contains("demob"))
            };
            let mut c = self.map_now();
            for _ in 0..10 {
                if demob(&c) {
                    break;
                }
                c = self.map_now();
            }
            let n = self.near();
            if st(get(&c, "chrom")) != chrom || !demob(&c) || !st(get(&n, "card")).contains("DEMOB")
            {
                bad.push(format!(
                    "{link}: Nearest flew to {} {}, card {:?}",
                    st(get(&c, "chrom")),
                    st(get(&c, "view")),
                    st(get(&n, "card"))
                ));
            }
            self.map_key("Escape", true);
            sleep(200);
            let card = truthy(get(&self.near(), "card"));
            if card || st(get(&self.map_now(), "level")) == "chromosome" {
                bad.push(format!(
                    "{link}: Esc leaves the card open or steps out with a card open"
                ));
            }
        }
        let doc = read_map(&self.map_json());
        let (i, f) = gated(&doc);
        let chrom = st(idx(get(&f, "chrom"), i));
        let si = num(idx(get(&f, "s"), i));
        let gene = st(idx(get(&f, "gene"), i));
        let others: Vec<f64> = (0..items(get(&f, "id")).len())
            .filter(|&j| j != i && st(idx(get(&f, "chrom"), j)) == chrom)
            .map(|j| num(idx(get(&f, "s"), j)))
            .collect();
        let a = (si - 3_000_000.0).max(0.0) as i64;
        if others.iter().any(|x| (x - si).abs() < 6_000_000.0) {
            bad.push(format!(
                "{gene}, closed, has another finding within 6 Mb: no view to check Nearest on"
            ));
        } else {
            self.map_page(S, &format!("#{chrom}:{}-{}", a + 1, a + 2_000_000));
            let n = self.near();
            if st(get(&n, "near")).contains(&gene) || st(get(&n, "card")).contains(&gene) {
                bad.push(format!(
                    "Nearest leads to {gene}, closed: {:?}",
                    st(get(&n, "near"))
                ));
            }
        }
        bad
    }

    /// Paint the copies: Origin is off with its reason, on a run in a family and on one in none.
    fn origin_problems(&mut self) -> Vec<String> {
        let js = "JSON.stringify([document.querySelector('.gm-paint [data-paint=origin]').getAttribute('aria-disabled'), \
                  document.querySelector('.gm-paint-why').textContent])";
        let mut bad = vec![];
        self.map_page(S, "");
        let got = strs(&self.jv(js));
        if got
            != [
                "true",
                "Origin: which parent each stretch came from is painted on the family page.",
            ]
        {
            bad.push(format!("in a family: {got:?}"));
        }
        let path = self.site.join(format!("runs/{S}/frames.json"));
        let keep = Keep::new(&path);
        let mut frames = loads(&keep.text()).unwrap();
        let list = if matches!(frames, Value::Object(_)) {
            at(&mut frames, &["frames"])
        } else {
            &mut frames
        };
        for fr in arr(list) {
            if let Value::Object(o) = fr
                && let Some(Value::Object(state)) = o.get_mut("state")
                && let Some(Value::Object(plugins)) = state.get_mut("plugins")
            {
                plugins.shift_remove("family");
            }
        }
        fs::write(&path, dumps(&frames, &Dump::default()).unwrap()).unwrap();
        self.map_page(S, "");
        let got = strs(&self.jv(js));
        if got
            != [
                "true",
                "Origin needs a family run with a sequenced parent: this run is in none.",
            ]
        {
            bad.push(format!("in no family: {got:?}"));
        }
        drop(keep);
        bad
    }

    fn too_wide(&mut self) -> String {
        st(&self.ev(
            "document.documentElement.scrollWidth > window.innerWidth ? \
                     document.documentElement.scrollWidth + ' > ' + window.innerWidth : ''",
        ))
    }
}

/// Figure 1 on the person pages: drawn from map.json with every mark the karyogram drew; its
/// levels, links and a place it cannot read; APOE before its step; a missing layer; the runs
/// sentence; a failed reload; no DecompressionStream; the closer levels; search, Nearest and
/// Origin; no sideways scroll at 390 px. (First paint and parse times: tests/browser/live_demo_perf.rs.)
#[test]
fn map_levels_search_and_veil() {
    let Some(mut d) = Demo::new("map-levels") else {
        return;
    };
    let rows = karyogram_rows();
    let bad = d.parity_problems(&rows);
    let n = rows.len();
    d.r.check(
        &format!("map: the DEMO-001 to DEMO-003 person pages draw every karyogram mark from map.json ({n} marks: same level and chromosome, within 0.01 Mb)"),
        if n >= 9 { bad } else { vec![format!("only {n} marks")] },
    );
    {
        let path = d.map_json();
        let _keep = Keep::new(&path);
        let mut doc = read_map(&path);
        let want = &rows.iter().find(|r| r.0 == S).unwrap().1;
        let f = at(&mut doc, &["per", S, "findings"]);
        let gone = items(get(f, "id"))
            .iter()
            .position(|v| st(v) == *want)
            .unwrap();
        for (_, col) in obj(f).iter_mut() {
            arr(col).remove(gone);
        }
        write_map(&path, &doc);
        let bad = d.parity_problems(&rows);
        d.r.control("a map.json without one of the karyogram's findings", bad);
    }
    let bad = d.state_problems();
    d.r.check("map: genome and chromosome levels; #chr16, #chr16:10.0-40.0M, chr=11 and a 1-based span open their view; a click on a row opens it, Esc closes it; an unknown place gives a toast", bad);

    let person = d.site.join("assets/person.js");
    let mjs = d.site.join("assets/map.js");
    let bad = d.replay_apoe_problems();
    d.r.check(
        "map: before the traits step a replay has no APOE chip, switch or mark",
        bad,
    );
    {
        let _a = plant(
            &person,
            &[(
                "if (!f) { var t = cards()",
                "if (!f) { return true; var t = cards()",
            )],
        );
        let _b = plant(
            &mjs,
            &[(
                "if (M.gated[id] && opts.now && !opts.now(id)) { return; }",
                "",
            )],
        );
        let bad = d.replay_apoe_problems();
        d.r.control(
            "person.js and map.js that show APOE's chip and row before its card",
            bad,
        );
    }
    let bad = d.missing_layer_problems();
    d.r.check("map: a missing layer says why in words, never a blank", bad);
    let bad = d.runs_sentence_problems();
    d.r.check(
        "map: a replay says an open chromosome's runs sentence only once the runs' step ended",
        bad,
    );
    {
        let _p = plant(&mjs, &[("&& rohLayer() === M.layers.roh", "")]);
        let bad = d.runs_sentence_problems();
        d.r.control("map.js whose runs sentence does not wait for its step", bad);
    }
    let bad = d.failed_reload_problems();
    d.r.check(
        "map: after a failed reload the figure draws nothing more of its old document",
        bad,
    );
    {
        let _p = plant_n(
            &mjs,
            &[
                ("if (!M || failure) { return; }", "if (!M) { return; }", 2),
                ("if (!raf && !failure) {", "if (!raf) {", 1),
                (
                    "window.cancelAnimationFrame(raf);\n      raf = 0;\n      window.clearTimeout(settleT);\n",
                    "",
                    1,
                ),
            ],
        );
        let bad = d.failed_reload_problems();
        d.r.control("map.js that keeps drawing after a failed reload", bad);
    }
    let bad = d.fallback_problems();
    d.r.check("map: a browser with no DecompressionStream gets the sentence and a findings table through the veil", bad);
    let bad = d.closer_problems();
    d.r.check("map: the region, gene and reads levels open from their addresses with their titles, paths and genes; Esc steps out one level at a time; the Region and Gene stops land on their level; + and → zoom and pan", bad);
    {
        let _p = plant(
            &mjs,
            &[(
                "if (bpp > CHROM_BPP) { return 'chromosome'; }",
                "return 'chromosome';",
            )],
        );
        let bad = d.closer_problems();
        d.r.control("map.js whose levels stop at the chromosome", bad);
    }
    let bad = d.find_problems();
    d.r.check("map: / opens the search; a gene, a band and a finding are found and Enter flies there (a finding with its card); APOE, closed, is never a line, by its gene, id or title", bad);
    {
        let _p = plant(
            &mjs,
            &[(
                "var fs = drawn();\n      fs.filter(",
                "var fs = M.findings;\n      fs.filter(",
            )],
        );
        let bad = d.find_problems();
        d.r.control(
            "map.js whose search looks through every finding, closed or not",
            bad,
        );
    }
    let bad = d.nearest_problems();
    d.r.check("map: Nearest names the closest finding (else the first on the next chromosome) and flies there with its card; Esc closes the card first; it never leads to APOE, closed", bad);
    {
        let _p = plant(
            &mjs,
            &[(
                "var ch = chrom(), v = view(), fs = drawn(),",
                "var ch = chrom(), v = view(), fs = M.findings,",
            )],
        );
        let bad = d.nearest_problems();
        d.r.control(
            "map.js whose Nearest leads to every finding, closed or not",
            bad,
        );
    }
    let bad = d.origin_problems();
    d.r.check("map: Paint the copies keeps Origin off with its reason, on a run in a family and on one in none", bad);
    {
        let _p = plant(
            &mjs,
            &[(
                "'aria-pressed': 'false', 'aria-disabled': 'true',",
                "'aria-pressed': 'false',",
            )],
        );
        let bad = d.origin_problems();
        d.r.control("map.js whose Origin button is not off", bad);
    }
    d.viewport(390, 844, 2);
    let mut wide = vec![];
    for link in ["", "#chr16"] {
        d.map_page(S, link);
        let w = d.too_wide();
        if !w.is_empty() {
            wide.push(format!(
                "{}: {w}",
                if link.is_empty() { "genome" } else { link }
            ));
        }
    }
    d.r.check(
        "map: at 390 px the map, open or not, scrolls nothing sideways",
        wide,
    );
    {
        d.ev("var x = document.createElement('div'); x.style.width = '900px'; x.style.height = '4px'; \
              document.getElementById('genome-map').appendChild(x); true");
        let w = d.too_wide();
        d.r.control(
            "a 900 px wide block in the map at 390 px",
            if w.is_empty() { vec![] } else { vec![w] },
        );
    }
    d.viewport(1440, 900, 1);
    d.r.done();
}

// ---------------------------------------------------------------------------------------------
// The In view table, the strip and drawsAt().

const PRED_JS: &str = r#"(async function (runs) {
  var LV = window.PGPMap.LEVELS, on = window.PGPMap.drawsAt, out = {}, genes = [];
  async function gz(u) {
    var r = await fetch(u, { cache: 'no-store' });
    if (!r.ok) { throw new Error(u + ': ' + r.status); }
    return JSON.parse(await new Response(r.body.pipeThrough(new DecompressionStream('gzip'))).text());
  }
  function line(lines, layer, c, s, e, id, row) {
    lines.push([layer, c, s, e, id].concat(LV.map(function (lv) { return on(layer, row, lv) ? 1 : 0; })).join('\t'));
  }
  function each(o, fn) { var at = 0; (o.ds || []).forEach(function (d, i) { at += d; fn(at, at + o.l[i], i); }); }
  var chroms = [];
  for (var i = 1; i <= 22; i++) { chroms.push('chr' + i); }
  chroms.push('chrX', 'chrY');
  for (var c = 0; c < chroms.length; c++) {
    var g = (await gz('/map/genes/' + chroms[c] + '.json')).g || {};
    each(g, function (s, e, j) { line(genes, 'genes', chroms[c], s, e, g.name[j], { name: g.name[j] }); });
  }
  for (var k = 0; k < runs.length; k++) {
    var d = await gz('/runs/' + runs[k] + '/map.json'), per = d.per[runs[k]], f = per.findings, lines = [];
    f.id.forEach(function (id, j) {
      line(lines, 'findings', f.chrom[j], f.s[j], f.e[j], id, { id: id, level: f.level[j], s: f.s[j], e: f.e[j] });
    });
    Object.keys(per.roh || {}).forEach(function (cn) {
      var r = per.roh[cn];
      each(r, function (s, e, j) { line(lines, 'roh', cn, s, e, '', { long: r.long[j] }); });
    });
    out[runs[k]] = lines.concat(genes);
  }
  return JSON.stringify(out);
})(RUNS)"#;

/// The In view table, the canvas and the strip, read at once (Show all first).
const INVIEW_JS: &str = r#"JSON.stringify((function () {
  var f = document.getElementById('genome-map'), cv = f.querySelector('.gm-canvas'), st = f.querySelector('.gm-strip');
  var all = f.querySelector('.gm-showall');
  if (all && all.textContent === 'Show all') { all.click(); }
  var trs = f.querySelectorAll('.gm-table tbody tr[data-key]');
  return { rows: Array.prototype.map.call(trs, function (r) { return r.dataset.key; }), drawn: (cv.pgpDrawn || []).slice(),
    ticks: st.hidden ? null : (st.pgpDrawn || []).slice(), stripH: st.hidden ? 0 : st.clientHeight,
    count: (f.querySelector('.gm-inview-n') || {}).textContent || '', live: f.querySelector('.gm-inview-live').textContent,
    dense: f.querySelector('.gm-inview-dense').textContent, open: f.querySelector('.gm-inview').open,
    summary: f.querySelector('.gm-inview > summary').textContent,
    sort: Array.prototype.map.call(f.querySelectorAll('.gm-table th'), function (t) { return t.getAttribute('aria-sort') || ''; }),
    cells: Array.prototype.map.call(trs, function (r) { return Array.prototype.map.call(r.cells, function (c) { return c.textContent; }); }),
    current: Array.prototype.map.call(f.querySelectorAll('.gm-table tr[aria-current=true]'), function (r) { return r.dataset.key; }),
    card: (f.querySelector('.gm-card') || {}).textContent || '', text: f.querySelector('.gm-inview').textContent };
}()))"#;

const SETTLED_JS: &str = "JSON.stringify((function () { var f = document.getElementById('genome-map'); return [f.dataset.view || '', \
    f.dataset.level || '', (f.querySelector('.gm-inview-n') || {}).textContent || '', \
    (f.querySelector('.gm-canvas').pgpDrawn || []).length, f.dataset.genes || '']; }()))";

const ANNOUNCE_JS: &str = "window.__said = []; new MutationObserver(function () { window.__said.push(document.querySelector(\
    '.gm-inview-live').textContent); }).observe(document.querySelector('.gm-inview-live'), \
    { childList: true, characterData: true, subtree: true }); true";

const SHOWALL_JS: &str = r#"(function () {
  var f = document.getElementById('genome-map'), a = f.querySelector('.gm-showall'), out = [];
  if (!a) { return JSON.stringify(null); }
  for (var i = 0; i < 2; i++) {
    a = f.querySelector('.gm-showall');
    a.focus();
    var was = a.textContent;
    a.click();
    var now = document.activeElement;
    out.push(was + ' > ' + (now && now.classList.contains('gm-showall') ? now.textContent : (now ? now.tagName : '')));
  }
  return JSON.stringify(out);
}())"#;

/// A run of homozygosity of 7 Mb planted in DEMO-001's step 11 text: the demo's own runs are all
/// shorter than 5 Mb, and a long one is the row whose answer differs at the genome level.
const LONG_RG: &str = "RG\tDEMO-001\tchr2\t60000001\t67000000\t7000000\t5000\t80.0\n";

/// Views at every level: a page load for each run's genome level, then the rest by address.
const INVIEW_VIEWS: [(&str, &str, &str); 12] = [
    (S, "", "genome"),
    (S, "#chr11", "chromosome"),
    (S, "#chr22", "chromosome"),
    (S, "#chr22:40.5-43.5M", "region"),
    (S, "#chr11:4.0-7.0M", "region"),
    (S, "#chr19:43.5-46.5M", "region"),
    (S, "#chr11:5.150-5.350M", "gene"),
    (S, "#chr22:42.000-42.200M", "gene"),
    (S, "#chr11:5225001-5250000", "reads"),
    (S, "#chr22:42110001-42140000", "reads"),
    (PARTNER, "", "genome"),
    (PARTNER, "#chr16", "chromosome"),
];

const PLACES: [(&str, &str, &str, &str); 3] = [
    ("", "chr7", "chr7", "chromosome"),
    ("#chr7", "#chr11:4.5-6.5M", "chr11", "region"),
    ("#chr7", "genome", "", "genome"),
];

/// "chrN:a-bM" read as (a, b) in Mb.
fn view_mb(h: &str) -> Option<(f64, f64)> {
    let rest = h.trim_start_matches('#').strip_prefix("chr")?;
    let (_, span) = rest.split_once(':')?;
    let (a, b) = span.strip_suffix('M')?.split_once('-')?;
    Some((a.parse().ok()?, b.parse().ok()?))
}

fn size_bp(words: &str) -> f64 {
    let w = words.replace(',', "");
    let mut p = w.split_whitespace();
    let (Some(n), Some(unit)) = (p.next(), p.next()) else {
        return f64::NAN;
    };
    let mult = match unit {
        "bp" => 1.0,
        "kb" => 1e3,
        "Mb" => 1e6,
        _ => f64::NAN,
    };
    n.parse::<f64>().unwrap_or(f64::NAN) * mult
}

fn key_field(k: &str, i: usize) -> i64 {
    k.split('\t')
        .nth(i)
        .and_then(|x| x.parse().ok())
        .unwrap_or(0)
}

impl Demo {
    fn inview_now(&mut self) -> Value {
        self.jv(INVIEW_JS)
    }

    fn inview_page(&mut self, run: &str, link: &str) -> (Value, Value) {
        let mut c = self.map_page(run, link);
        if ["region", "gene", "reads"].contains(&st(get(&c, "level")).as_str()) {
            sleep(800); // the chromosome's gene file arrives after the first frame
            c = self.map_now();
        }
        let t = self.inview_now();
        (c, t)
    }

    /// Waits until the map has settled on `view` (its address as the map writes it).
    fn settled_at(&mut self, view: &str) -> (Value, Value) {
        let end = Instant::now() + Duration::from_secs(8);
        let mut last = Value::Null;
        while Instant::now() < end {
            sleep(250);
            let now = self.jv(SETTLED_JS);
            if st(idx(&now, 0)) == view && now == last {
                break;
            }
            last = now;
        }
        let c = self.jv(MAP_JS);
        (c, self.inview_now())
    }

    /// The open page moved to `link` without a reload, settled.
    fn go_view(&mut self, link: &str) -> (Value, Value) {
        self.ev(&format!("location.hash = {}; true", js_str(link)));
        self.settled_at(&link[1..])
    }

    fn map_command(&self, run: &str, predicate: bool) -> (bool, Vec<u8>) {
        let results =
            fs::read_to_string(self.site.join(format!("runs/{run}/results.json"))).unwrap();
        let mut cmd = Command::new(&self.og);
        cmd.args(["map", "build", "--root"])
            .arg(&self.work)
            .args(["--run", run, "--assets"])
            .arg(repo_root().join("assets"))
            .arg("--cache-dir")
            .arg(self.work.join(".owngenome-live/map"));
        if predicate {
            cmd.arg("--predicate-table");
        }
        let mut child = cmd
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        write!(child.stdin.take().unwrap(), "{{\"results\": {results}}}").unwrap();
        let out = child.wait_with_output().unwrap();
        (out.status.success(), out.stdout)
    }

    /// drawsAt() of map.js against og map build --predicate-table, row by row: the same rows and the
    /// same answer at each level; the rows cover a finding, a long run, a shorter run and a gene.
    fn predicate_problems(&mut self, runs: &[&str]) -> Vec<String> {
        let mut bad = vec![];
        self.map_page(S, "");
        let list = Value::Array(runs.iter().map(|r| s(r)).collect());
        let js = self.awaited(&PRED_JS.replace("RUNS", &dumps(&list, &Dump::default()).unwrap()));
        let mut seen = BTreeSet::new();
        for run in runs {
            let (ok, out) = self.map_command(run, true);
            let out = String::from_utf8_lossy(&out);
            let mut lines = out.lines();
            if !ok
                || lines.next()
                    != Some("layer\tchrom\ts\te\tid\tgenome\tchromosome\tregion\tgene\treads")
            {
                bad.push(format!(
                    "{run}: og map build --predicate-table failed or has another header"
                ));
                continue;
            }
            let mut want: Vec<String> = lines.map(str::to_owned).collect();
            want.sort();
            let mut got = strs(get(&js, run));
            got.sort();
            if want != got {
                let w: BTreeSet<_> = want.iter().collect();
                let g: BTreeSet<_> = got.iter().collect();
                bad.push(format!(
                    "{run}: og and drawsAt() differ: og only {:?}, drawsAt() only {:?}",
                    w.difference(&g).take(4).collect::<Vec<_>>(),
                    g.difference(&w).take(4).collect::<Vec<_>>()
                ));
            }
            for x in &want {
                let c: Vec<_> = x.split('\t').collect();
                seen.insert((c[0].to_owned(), c[5].to_owned()));
            }
        }
        for want in [
            ("findings", "1"),
            ("roh", "1"),
            ("roh", "0"),
            ("genes", "0"),
        ] {
            if !seen.contains(&(want.0.to_owned(), want.1.to_owned())) {
                bad.push(format!("the rows compared miss {want:?}"));
            }
        }
        bad
    }

    /// On views at every level: the table's rows are exactly what the canvas drew one by one, each
    /// once, its count says how many, and the shorter runs drawn as a band have a line of their own
    /// at the genome level; the strip is hidden at the genome level and over an open chromosome ticks
    /// exactly the findings the genome level draws; T folds the table and opens it again.
    fn inview_problems(&mut self, views: &[(&str, &str, &str)]) -> Vec<String> {
        let mut bad = vec![];
        let mut genome: HashMap<String, Vec<String>> = HashMap::new();
        let mut page = "";
        for &(run, link, lv) in views {
            let (c, t) = if run != page {
                page = run;
                self.inview_page(run, link)
            } else {
                self.go_view(link)
            };
            let wh = format!("{run}{}", if link.is_empty() { " genome" } else { link });
            if st(get(&c, "level")) != lv {
                bad.push(format!("{wh}: level {}, want {lv}", st(get(&c, "level"))));
                continue;
            }
            let rows = strs(get(&t, "rows"));
            let mut a = rows.clone();
            a.sort();
            let mut b = strs(get(&t, "drawn"));
            b.sort();
            let uniq: BTreeSet<_> = rows.iter().collect();
            if a != b || uniq.len() != rows.len() {
                bad.push(format!(
                    "{wh}: the table's rows ({}) are not what the canvas drew ({})",
                    a.len(),
                    b.len()
                ));
            }
            if st(get(&t, "count")) != count_words(rows.len()) {
                bad.push(format!(
                    "{wh}: the count says {:?} for {} rows",
                    st(get(&t, "count")),
                    rows.len()
                ));
            }
            if lv == "genome" && !rows.iter().any(|k| k.starts_with("findings\t")) {
                bad.push(format!("{wh}: no finding in the table"));
            }
            if ["region", "gene", "reads"].contains(&lv)
                && !rows.iter().any(|k| k.starts_with("genes\t"))
            {
                bad.push(format!("{wh}: no gene in the table"));
            }
            let ticks = get(&t, "ticks");
            let strip = num(get(&t, "stripH"));
            if lv == "genome" {
                let mut f: Vec<_> = rows
                    .iter()
                    .filter(|k| k.starts_with("findings\t"))
                    .cloned()
                    .collect();
                f.sort();
                genome.insert(run.to_owned(), f);
                if !matches!(ticks, Value::Null) || strip != 0.0 {
                    bad.push(format!("{wh}: the strip shows at the genome level"));
                }
                if !st(get(&t, "dense"))
                    .contains("shorter runs of homozygosity, drawn as a faint band")
                {
                    bad.push(format!(
                        "{wh}: no line for the shorter runs: {:?}",
                        st(get(&t, "dense"))
                    ));
                }
            } else {
                let mut tk = strs(ticks);
                tk.sort();
                if matches!(ticks, Value::Null)
                    || strip != 34.0
                    || genome.get(run).is_some_and(|g| *g != tk)
                {
                    bad.push(format!(
                        "{wh}: strip {strip} px, ticks {tk:?}, want {:?}",
                        genome.get(run)
                    ));
                }
            }
        }
        if views.len() == INVIEW_VIEWS.len() {
            let mut folds = vec![];
            for _ in 0..2 {
                self.map_key("t", false);
                folds.push(truthy(get(&self.inview_now(), "open")));
            }
            if folds != [false, true] {
                bad.push(format!("T folds and opens the table: {folds:?}"));
            }
        }
        bad
    }

    /// APOE, closed on every page until opened, is no row and no tick on DEMO-001's page.
    fn inview_apoe_problems(&mut self) -> Vec<String> {
        let mut bad = vec![];
        let (_, t) = self.inview_page(S, "");
        if strs(get(&t, "rows")).iter().any(|k| k.contains("\tapoe")) {
            bad.push(format!("{S}: APOE, closed, is a row"));
        }
        let (_, t) = self.go_view("#chr19");
        if strs(get(&t, "ticks"))
            .iter()
            .chain(strs(get(&t, "rows")).iter())
            .any(|k| k.contains("\tapoe"))
        {
            bad.push(format!(
                "{S}: APOE, closed, is a strip tick or a row at chr19"
            ));
        }
        bad
    }

    /// A closed result is never a row or a tick; the son's page, closed whole by the report gate,
    /// has no map, table or strip at all.
    fn inview_veil_problems(&mut self) -> Vec<String> {
        let mut bad = self.inview_apoe_problems();
        self.go(&format!("/runs/{SON}/#chr7"));
        if !self.wait("document.documentElement.dataset.pgpGated === 'true'", 20) {
            bad.push(format!("{SON}: his page is not closed by the report gate"));
        } else if truthy(
            &self.ev("!!document.querySelector('#genome-map, .gm-inview, .gm-strip, canvas')"),
        ) {
            bad.push(format!(
                "{SON}: his closed page has a map, a table or a strip"
            ));
        }
        bad
    }

    /// By place with the findings first (aria-sort ascending on Where); Where again reverses it;
    /// Size sorts by size, smallest first.
    fn sort_problems(&mut self, default_only: bool) -> Vec<String> {
        let mut bad = vec![];
        let (_, t) = self.inview_page(S, "#chr22:40.5-43.5M");
        let keys = strs(get(&t, "rows"));
        let mut place = keys.clone();
        place.sort_by(|a, b| {
            (key_field(a, 2), key_field(a, 3), a).cmp(&(key_field(b, 2), key_field(b, 3), b))
        });
        let first: Vec<_> = place
            .iter()
            .filter(|k| k.starts_with("findings\t"))
            .cloned()
            .collect();
        let want: Vec<_> = first
            .iter()
            .cloned()
            .chain(place.iter().filter(|k| !first.contains(k)).cloned())
            .collect();
        if keys != want || place.first().is_some_and(|p| first.contains(p)) {
            bad.push(format!("default order {:?}: want the findings first, then by place (and a gene before the finding)", &keys[..keys.len().min(4)]));
        }
        let sort = strs(get(&t, "sort"));
        if sort != ["", "", "ascending", "", "", "", ""] {
            bad.push(format!("default aria-sort {sort:?}"));
        }
        if default_only {
            return bad;
        }
        self.ev("document.querySelector('.gm-table th button[data-col=where]').click(); true");
        let r = self.inview_now();
        let mut rev = keys.clone();
        rev.reverse();
        if strs(get(&r, "rows")) != rev || st(idx(get(&r, "sort"), 2)) != "descending" {
            bad.push(format!(
                "Where again: aria-sort {:?}",
                strs(get(&r, "sort"))
            ));
        }
        self.ev("document.querySelector('.gm-table th button[data-col=size]').click(); true");
        let r = self.inview_now();
        let sizes: Vec<f64> = items(get(&r, "cells"))
            .iter()
            .map(|c| size_bp(&st(idx(c, 3))))
            .collect();
        let sorted = sizes.windows(2).all(|w| w[0] <= w[1]);
        if !sorted || strs(get(&r, "sort")) != ["", "", "", "ascending", "", "", ""] {
            bad.push(format!(
                "Size: {:?}, aria-sort {:?}",
                &sizes[..sizes.len().min(5)],
                strs(get(&r, "sort"))
            ));
        }
        bad
    }

    /// The table announces its count once when it changes and never when it does not.
    fn announce_problems(&mut self, pans_only: bool) -> Vec<String> {
        let mut bad = vec![];
        let (_, t) = self.inview_page(S, "#chr3:100.0-101.0M");
        self.ev(ANNOUNCE_JS);
        for _ in 0..3 {
            self.map_key("ArrowRight", true);
            sleep(400);
        }
        self.settled_at("chr3:100.300-101.300M");
        let said = st(&self.ev("JSON.stringify(window.__said)"));
        if st(get(&t, "count")) != "0 features" || said != "[]" {
            bad.push(format!(
                "panning an empty stretch ({}) announced {said}",
                st(get(&t, "count"))
            ));
        }
        if pans_only {
            return bad;
        }
        let (_, t) = self.go_view("#chr11");
        let said = strs(&self.jv("JSON.stringify(window.__said)"));
        let count = st(get(&t, "count"));
        if said != [format!("{count} in view")] || count != "4 features" {
            bad.push(format!("#chr11 ({count}) announced {said:?}"));
        }
        bad
    }

    /// The keyboard stays on Show all when the table is built again, and the count's live line sits
    /// outside the foldable table.
    fn showall_problems(&mut self) -> Vec<String> {
        let mut bad = vec![];
        self.inview_page(S, "#chr22:40.5-43.5M");
        if truthy(&self.ev(
            "(function () { var l = document.querySelector('#genome-map .gm-inview-live'); \
                            return !l || !!l.closest('details'); }())",
        )) {
            bad.push("the count's live line is missing or inside the foldable table".into());
        }
        let said = self.jv(SHOWALL_JS);
        if matches!(said, Value::Null) {
            bad.push("no Show all at chr22:40.5-43.5M: the check tests nothing".into());
            return bad;
        }
        let mut said = strs(&said);
        said.sort();
        if said != ["Show all > Show fewer", "Show fewer > Show all"] {
            bad.push(format!(
                "after Show all and Show fewer the keyboard is on {said:?}"
            ));
        }
        bad
    }

    fn click_row(&mut self, key: &str) {
        self.ev(&format!(
            "document.querySelector({}).click(); true",
            js_str(&format!(".gm-table tr[data-key=\"{key}\"] .gm-go"))
        ));
    }

    /// A row's button flies the map to its feature: a finding's row with its card, the row then
    /// current; a run's row to a view that holds the run, with no card.
    fn fly_problems(&mut self, finding_only: bool) -> Vec<String> {
        let mut bad = vec![];
        let (_, t) = self.inview_page(S, "");
        let Some(key) = strs(get(&t, "rows"))
            .into_iter()
            .find(|k| k.starts_with("findings\tchr16\t"))
        else {
            return vec!["no chr16 finding row at the genome level".into()];
        };
        self.click_row(&key);
        let mut c = Value::Null;
        for _ in 0..6 {
            c = self.map_now();
            if st(get(&c, "chrom")) == "chr16" {
                break;
            }
        }
        let t = self.inview_now();
        if st(get(&c, "chrom")) != "chr16"
            || !strs(get(&t, "drawn")).contains(&key)
            || !st(get(&t, "card")).contains("DEMOC")
            || strs(get(&t, "current")) != [key.clone()]
        {
            bad.push(format!(
                "the DEMOC row flew to {} {}, current {:?}",
                st(get(&c, "chrom")),
                st(get(&c, "view")),
                strs(get(&t, "current"))
            ));
        }
        if finding_only {
            return bad;
        }
        let (_, t) = self.go_view("#chr11");
        let Some(key) = strs(get(&t, "rows"))
            .into_iter()
            .find(|k| k.starts_with("roh\t"))
        else {
            bad.push("no run of homozygosity row at chr11".into());
            return bad;
        };
        self.click_row(&key);
        let c = self.map_now();
        let t = self.inview_now();
        let (s0, e0) = (key_field(&key, 2) as f64, key_field(&key, 3) as f64);
        let view = st(get(&c, "view"));
        let ok = view.starts_with("chr11:")
            && view_mb(&view).is_some_and(|(a, b)| a * 1e6 <= s0 + 1e5 && b * 1e6 >= e0 - 1e5);
        if !ok || truthy(get(&t, "card")) || strs(get(&t, "current")) != [key.clone()] {
            bad.push(format!(
                "the run's row {key:?} flew to {view}, current {:?}",
                strs(get(&t, "current"))
            ));
        }
        bad
    }

    /// Where `frac` of `chrom` sits on the strip, in page px (map.js stripGeo() and sx(); `x0` and
    /// `pad` are a phone's 64 and 10).
    fn strip_x(
        &mut self,
        doc: &Value,
        chrom: &str,
        frac: f64,
        x0: f64,
        pad: f64,
    ) -> (f64, f64, f64) {
        let (l, t, w, h) = self.rect(".gm-strip");
        let no_y = matches!(get(idx(get(doc, "people"), 0), "has_y"), Value::Bool(false));
        let rows: Vec<_> = items(get(get(doc, "ref"), "chroms"))
            .iter()
            .filter(|c| !(st(get(c, "n")) == "chrY" && no_y))
            .collect();
        let total: f64 =
            rows.iter().map(|c| num(get(c, "len"))).sum::<f64>() + 5e6 * (rows.len() - 1) as f64;
        let mut off = 0.0;
        for c in rows {
            let len = num(get(c, "len"));
            if st(get(c, "n")) == chrom {
                let y = t + h / 2.0 - if h > 20.0 { 2.0 } else { 0.0 };
                return (l + x0 + (off + frac * len) * (w - pad - x0) / total, y, len);
            }
            off += len + 5e6;
        }
        panic!("{chrom} is not in the map's reference")
    }

    /// The strip: a click opens a chromosome whole, a drag brushes a new view, a drag on the brush
    /// pans; at 390 px it shows 12 px high and a click there does nothing, and the table starts
    /// folded with its count in its summary.
    fn strip_problems(&mut self, phone_only: bool, click_only: bool) -> Vec<String> {
        let mut bad = vec![];
        let doc = read_map(&self.map_json());
        if !phone_only {
            self.map_page(S, "#chr11");
            self.ev("document.querySelector('.gm-strip').scrollIntoView({block: 'center'}); true");
            sleep(300);
            let (x, y, _) = self.strip_x(&doc, "chr16", 0.5, 128.0, 22.0);
            self.drag(x, y, x);
            let h = self.hash_now(350);
            if h != "#chr16" {
                bad.push(format!("a click on chr16 of the strip opens {h}"));
            }
            if click_only {
                return bad;
            }
            sleep(1200);
            let (x0, y, n) = self.strip_x(&doc, "chr2", 0.2, 128.0, 22.0);
            let x1 = self.strip_x(&doc, "chr2", 0.4, 128.0, 22.0).0;
            self.drag(x0, y, x1);
            let h = self.hash_now(350);
            let v = view_mb(&h);
            let near = v.is_some_and(|(a, b)| {
                (a - 0.2 * n / 1e6).abs() <= 4.0 && (b - 0.4 * n / 1e6).abs() <= 4.0
            });
            if !h.starts_with("#chr2:") || !near {
                bad.push(format!("a brush over chr2 20-40% opens {h}"));
            }
            if let Some((a, b)) = v {
                sleep(1200);
                let mid = self
                    .strip_x(&doc, "chr2", (a + b) / 2.0 * 1e6 / n, 128.0, 22.0)
                    .0;
                self.drag(mid, y, mid + 12.0);
                let h = self.hash_now(600);
                let p = view_mb(&h);
                if !p.is_some_and(|p| p.0 - a > 5.0) {
                    bad.push(format!(
                        "a drag on the brush from {a}-{b} moved the view to {p:?}"
                    ));
                }
            }
        }
        self.viewport(390, 844, 2);
        self.map_page(S, "#chr11");
        let t = self.inview_now();
        let summary = st(get(&t, "summary"));
        if num(get(&t, "stripH")) != 12.0
            || truthy(get(&t, "open"))
            || !summary.starts_with("In view · ")
            || !summary.ends_with("features")
        {
            bad.push(format!(
                "at 390 px: strip {} px, table open {}, summary {summary:?}",
                num(get(&t, "stripH")),
                truthy(get(&t, "open"))
            ));
        }
        self.ev("document.querySelector('.gm-strip').scrollIntoView({block: 'center'}); true");
        sleep(300);
        let (x, y, _) = self.strip_x(&doc, "chr16", 0.5, 64.0, 10.0);
        self.drag(x, y, x);
        let h = self.hash_now(350);
        if h != "#chr11" {
            bad.push(format!(
                "at 390 px a click on chr16 of the strip moved the map to {h}"
            ));
        }
        self.viewport(1440, 900, 1);
        bad
    }

    /// The search box takes a place as the address writes it, and Enter flies there.
    fn place_find_problems(&mut self, places: &[(&str, &str, &str, &str)]) -> Vec<String> {
        let mut bad = vec![];
        for &(link, q, chrom, lv) in places {
            self.map_page(S, link);
            let found = self.map_find(q);
            if found.first().is_none_or(|f| st(get(f, "kind")) != "place") {
                bad.push(format!("{q}: lines {}", show(&Value::Array(found))));
                continue;
            }
            self.enter_find();
            let mut c = Value::Null;
            for _ in 0..5 {
                c = self.map_now();
                if st(get(&c, "level")) == lv && (chrom.is_empty() || st(get(&c, "chrom")) == chrom)
                {
                    break;
                }
            }
            if st(get(&c, "level")) != lv || (!chrom.is_empty() && st(get(&c, "chrom")) != chrom) {
                bad.push(format!(
                    "{q}: Enter flies to {} at {}, want {chrom} at {lv}",
                    st(get(&c, "chrom")),
                    st(get(&c, "level"))
                ));
            }
        }
        bad
    }

    /// Tab from a finding's mark to its card's link: the view settles again and the keyboard stays
    /// on the link.
    fn card_focus_problems(&mut self) -> Vec<String> {
        self.map_page(S, "#chr11");
        if !truthy(&self.ev(
            "(function () { var a = document.querySelector('#genome-map .gm-marks a[data-id]'); \
                             if (!a) { return false; } a.focus(); return true; }())",
        )) {
            return vec!["no finding's mark on chr11: the check tests nothing".into()];
        }
        sleep(600);
        if !truthy(&self.ev("(function () { var l = document.querySelector('#genome-map .gm-card:not([hidden]) a'); \
                             if (!l) { return false; } l.focus(); return document.activeElement === l; }())"))
        {
            return vec!["the finding's card has no link to focus: the check tests nothing".into()];
        }
        sleep(1000);
        if truthy(&self.ev(
            "!!(document.activeElement && document.activeElement.closest && \
                            document.activeElement.closest('#genome-map .gm-card'))",
        )) {
            vec![]
        } else {
            vec![format!(
                "the keyboard left the card's link once the view settled: {}",
                st(&self.ev("document.activeElement ? document.activeElement.tagName : ''"))
            )]
        }
    }

    /// A reload whose document gives a finding new words under the same place and id: the table
    /// shows the new words.
    fn reload_table_problems(&mut self) -> Vec<String> {
        let other = self.site.join(format!("runs/{S}/map-reload.json"));
        let bad = (|| {
            self.map_page(S, "");
            self.ev(SECOND_JS);
            if !self.wait(
                "document.getElementById('gm-second').dataset.pgpMap === 'drawn'",
                20,
            ) {
                return vec!["the second figure never drew".into()];
            }
            sleep(600);
            let before = strs(&self.jv(ROWS_JS));
            let mut doc = read_map(&self.map_json());
            let titles = arr(at(&mut doc, &["per", S, "findings", "title"]));
            let first = before.first().cloned().unwrap_or_default();
            if before.is_empty() || titles.iter().filter(|t| st(t) == first).count() != 1 {
                return vec![format!(
                    "no finding row to rename: the check tests nothing ({before:?})"
                )];
            }
            let i = titles.iter().position(|t| st(t) == first).unwrap();
            titles[i] = s("Planted new words for this finding");
            write_map(&other, &doc);
            self.ev("pgpSecond.o.url = 'map-reload.json'; pgpSecond.g.refresh('k2'); true");
            let probe = format!(
                "{}.indexOf('Planted new words for this finding') >= 0",
                ROWS_JS.replace("JSON.stringify", "")
            );
            if self.wait(&probe, 5) {
                return vec![];
            }
            vec![format!(
                "after the reload the table still reads {:?}",
                strs(&self.jv(ROWS_JS))
            )]
        })();
        let _ = fs::remove_file(&other);
        bad
    }
}

/// The whole-genome strip, the In view table and drawsAt(), each check with its planted mistake.
/// (Show all's window and the frame budget: tests/browser/live_demo_perf.rs.)
#[test]
fn map_in_view_table_and_strip() {
    let Some(mut d) = Demo::new("map-in-view") else {
        return;
    };
    let mjs = d.site.join("assets/map.js");
    {
        let roh = d.work.join(format!("{S}/vcf/{S}_roh.txt"));
        let (_k1, _k2) = (Keep::new(&roh), Keep::new(&d.map_json()));
        fs::OpenOptions::new()
            .append(true)
            .open(&roh)
            .unwrap()
            .write_all(LONG_RG.as_bytes())
            .unwrap();
        let (ok, map) = d.map_command(S, false);
        assert!(ok, "og map build failed on the planted long run");
        fs::write(d.map_json(), map).unwrap();
        let bad = d.predicate_problems(&[S, PARTNER, SON]);
        d.r.check("map: drawsAt() in map.js gives og map build --predicate-table's answer for every row of the three demo runs and every gene, at every level", bad);
        let _p = plant(
            &mjs,
            &[("return i > 0 || !!(row && row.long);", "return i > 0;")],
        );
        let bad = d.predicate_problems(&[S]);
        d.r.control(
            "map.js whose drawsAt() drops the long runs at the genome level",
            bad,
        );
    }
    let bad = d.inview_problems(&INVIEW_VIEWS);
    d.r.check("map: on views at every level the In view table lists exactly what the canvas draws one by one, with its count and a line for the runs drawn as a band; the strip shows only over an open chromosome and ticks the findings the genome level draws; T folds and opens the table", bad);
    {
        let _p = plant(
            &mjs,
            &[(
                "if (drawsAt('roh', { long: R.o.long[i] }, lv)) { rows.push(",
                "if (true) { rows.push(",
            )],
        );
        let bad = d.inview_problems(&INVIEW_VIEWS[..1]);
        d.r.control(
            "map.js whose table lists every run of homozygosity, drawn one by one or not",
            bad,
        );
    }
    {
        let _p = plant(
            &mjs,
            &[(
                "        if (!drawsAt('findings', r, 'genome')) { return; }",
                "        if (!drawsAt('findings', r, 'genome') || r.chrom !== chrom().n) { return; }",
            )],
        );
        let bad = d.inview_problems(&INVIEW_VIEWS[..2]);
        d.r.control(
            "map.js whose strip ticks the findings of the open chromosome only",
            bad,
        );
    }
    let bad = d.inview_veil_problems();
    d.r.check("map: a closed result is no row of the In view table and no tick of the strip (APOE on an adult's page); the son's page, closed by the report gate, has no map", bad);
    {
        let _p = plant(
            &mjs,
            &[(
                "drawn().forEach(function (r) { if (meets(",
                "M.findings.forEach(function (r) { if (meets(",
            )],
        );
        let bad = d.inview_apoe_problems();
        d.r.control("map.js whose table lists every finding, closed or not", bad);
    }
    {
        let _p = plant(
            &mjs,
            &[(
                "      drawn().forEach(function (r) {\n        if (!drawsAt('findings', r, 'genome'))",
                "      M.findings.forEach(function (r) {\n        if (!drawsAt('findings', r, 'genome'))",
            )],
        );
        let bad = d.inview_apoe_problems();
        d.r.control("map.js whose strip ticks every finding, closed or not", bad);
    }
    let bad = d.sort_problems(false);
    d.r.check("map: the In view table sorts by place with the findings first (aria-sort ascending on Where), Where again reverses it, and Size sorts by size", bad);
    {
        let _p = plant(
            &mjs,
            &[(
                "return (a.layer === 'findings' ? 0 : 1) - (b.layer === 'findings' ? 0 : 1) || place(a, b);",
                "return place(a, b);",
            )],
        );
        let bad = d.sort_problems(true);
        d.r.control("map.js whose table sorts by place alone", bad);
    }
    let bad = d.announce_problems(false);
    d.r.check("map: the In view table announces its count once when it changes, never while it stays the same", bad);
    {
        let _p = plant(
            &mjs,
            &[(
                "if (n !== announced) { announced = n; inLive",
                "if (true) { announced = n; inLive",
            )],
        );
        let bad = d.announce_problems(true);
        d.r.control(
            "map.js whose table announces its count on every settle",
            bad,
        );
    }
    let bad = d.showall_problems();
    d.r.check("map: the keyboard stays on Show all and Show fewer, and the count is said outside the foldable table", bad);
    {
        let _p = plant(
            &mjs,
            &[("if (more) { more.focus(", "if (false) { more.focus(")],
        );
        let bad = d.showall_problems();
        d.r.control("map.js whose Show all drops the keyboard", bad);
    }
    {
        let _p = plant(
            &mjs,
            &[(
                "main.insertBefore(inLive, inview.nextSibling);",
                "inview.insertBefore(inLive, inWrap);",
            )],
        );
        let bad = d.showall_problems();
        d.r.control("map.js whose count's live line folds with the table", bad);
    }
    let bad = d.fly_problems(false);
    d.r.check("map: a row of the In view table flies the map to its feature (a finding with its card) and becomes the current row", bad);
    {
        let _p = plant(
            &mjs,
            &[(
                "if (x.layer === 'findings') { pick(x.r.id); } else { pick(null); }",
                "pick(null);",
            )],
        );
        let bad = d.fly_problems(true);
        d.r.control("map.js whose rows fly without the finding's card", bad);
    }
    let bad = d.strip_problems(false, false);
    d.r.check("map: the strip opens a chromosome on a click, brushes a new view and pans on a drag; at 390 px it is a 12 px picture and the table starts folded", bad);
    {
        let _p = plant(
            &mjs,
            &[(
                "if (i >= 0) { openChrom(i, undefined, undefined, true); }\n      } else if (b0.mode",
                "void i;\n      } else if (b0.mode",
            )],
        );
        let bad = d.strip_problems(false, true);
        d.r.control("map.js whose strip ignores a click", bad);
    }
    {
        let _p = plant(
            &mjs,
            &[
                (
                    "if (!M || phone || !b0) { return; }",
                    "if (!M || !b0) { return; }",
                ),
                (
                    "if (!M || phone || e.button !== 0 || level() === 'genome') { return; }",
                    "if (!M || e.button !== 0 || level() === 'genome') { return; }",
                ),
            ],
        );
        let bad = d.strip_problems(true, false);
        d.r.control("map.js whose strip takes clicks on a phone", bad);
    }
    let bad = d.place_find_problems(&PLACES);
    d.r.check("map: the search box takes a place as the address writes it (chr7, #chr11:4.5-6.5M, genome) and Enter flies there", bad);
    {
        let _p = plant(
            &mjs,
            &[(
                "var p = /^(genome|chr[0-9XYMT]+(:.*)?)$/i.test(pq) ? parsePlace(pq) : null;",
                "var p = /^chr/i.test(pq) && /[:]/.test(pq) ? parsePlace(pq) : null;",
            )],
        );
        let bad = d.place_find_problems(&PLACES[..1]);
        d.r.control("map.js whose search takes a place only with a colon", bad);
    }
    let bad = d.card_focus_problems();
    d.r.check(
        "map: the keyboard stays on the card's link when the view settles again",
        bad,
    );
    {
        let _p = plant(
            &mjs,
            &[("if (sig === cardSig && !card.hidden) { return; }", "")],
        );
        let bad = d.card_focus_problems();
        d.r.control("map.js whose card is built again on every settle", bad);
    }
    let bad = d.reload_table_problems();
    d.r.check("map: a reloaded document with new words for a finding under the same key shows them in the In view table", bad);
    {
        let _p = plant(
            &mjs,
            &[(
                "        tSig = '';  // a new document may",
                "        void 0;  // a new document may",
            )],
        );
        let bad = d.reload_table_problems();
        d.r.control("map.js whose table keeps its rows across a reload", bad);
    }
    d.r.done();
}

// ---------------------------------------------------------------------------------------------
// The person page: layout, phone report, folds.

const LAYOUT_JS: &str = r#"JSON.stringify((function () {
  var main = document.getElementById('content'), fold = document.getElementById('progress-fold');
  var band = document.getElementById('runhero'), lead = document.getElementById('lead-slot');
  var wire = Array.prototype.map.call(document.querySelectorAll('#found > li[role=button]'), function (li) {
    return (li.childNodes[1] || {}).textContent || ''; });
  var hla = document.getElementById('card-hla');
  var shownLines = hla ? hla.querySelectorAll('.card-text .lines-box > ul.lines > li').length : -1;
  return {
    last: main.lastElementChild && main.lastElementChild.id, folded: !!fold && !fold.open,
    band: !band.hidden, bandFirst: !!(band.compareDocumentPosition(lead) & Node.DOCUMENT_POSITION_FOLLOWING),
    wire: wire, hla: shownLines, hlaFold: !!(hla && hla.querySelector('details.lines-more')),
    margin: !!(hla && hla.querySelector('.card-margin .caveat .caveat-lead'))
  };
}()))"#;

const PHONE_JS: &str = r#"JSON.stringify((function () {
  // The content of a closed details keeps a layout box (content-visibility: hidden): ask whether it shows.
  var shown = function (e) { return !!e && e.checkVisibility() && e.getBoundingClientRect().height > 0; };
  var loud = Array.prototype.filter.call(document.querySelectorAll('#results [data-finding]'), function (e) {
    return /^(act|medicine|carrier|review)$/.test(e.dataset.level); });
  return { height: document.documentElement.scrollHeight, loud: loud.length,
    hidden: loud.filter(function (e) { return !shown(e); }).map(function (e) { return e.dataset.finding; }),
    deep: shown(document.getElementById('card-' + 'roh')),
    sec: shown(document.querySelector('#s7 details.secfold details.pfold > summary')),
    folds: document.querySelectorAll('#results details.pfold').length };
}()))"#;

/// CSS px: the Journal design's phone page is about 4,600; every card open was about 37,000.
const PHONE_MAX: f64 = 8000.0;
/// A card of section 5 with no finding of its own, flagged for the fold checks.
const FLAGGED: &str = "callability";

impl Demo {
    /// The layout facts of DEMO-001's page at a replay frame (mid, or the end).
    fn person_layout(&mut self, mid: bool) -> Value {
        self.go(&format!("/runs/{S}/"));
        self.wait("document.documentElement.dataset.pgpStatus === 'ready' && document.documentElement.dataset.pgpFrame", 20);
        self.slider(if mid {
            "Math.floor(s.max / 2)"
        } else {
            "s.max"
        });
        let frame = if mid {
            "String(Math.floor(document.getElementById('replay-slider').max / 2))"
        } else {
            "document.getElementById('replay-slider').max"
        };
        self.wait(&format!("document.documentElement.dataset.pgpResults && document.documentElement.dataset.pgpFrame === {frame}"), 20);
        sleep(500);
        self.jv(LAYOUT_JS)
    }

    fn layout_problems(&mut self) -> Vec<String> {
        let (mid, end) = (self.person_layout(true), self.person_layout(false));
        let mut out = vec![];
        for (name, got) in [("mid", &mid), ("end", &end)] {
            if st(get(got, "last")) != "progress" || !truthy(get(got, "folded")) {
                out.push(format!(
                    "{name}: the step timeline is not the last block, folded: {}, {}",
                    st(get(got, "last")),
                    truthy(get(got, "folded"))
                ));
            }
            let wire = strs(get(got, "wire"));
            if wire.is_empty() || wire.len() > 8 {
                out.push(format!(
                    "{name}: the wire shows {} lines, want 1 to 8",
                    wire.len()
                ));
            }
            let plain = wire
                .iter()
                .filter(|w| w.contains("finished, see the report"))
                .count();
            if plain > 1 {
                out.push(format!(
                    "{name}: {plain} plain 'finished' lines on the wire, want them folded into one"
                ));
            }
            let long: Vec<_> = wire.iter().filter(|w| w.chars().count() > 141).collect();
            if !long.is_empty() {
                out.push(format!(
                    "{name}: wire lines longer than one sentence: {:?}",
                    &long[..long.len().min(2)]
                ));
            }
        }
        if !truthy(get(&mid, "band")) || !truthy(get(&mid, "bandFirst")) {
            out.push("mid: the run band is not shown above the lead".into());
        }
        if truthy(get(&end, "band")) {
            out.push("end: the run band still shows on a finished run".into());
        }
        let hla = num(get(&end, "hla"));
        if !(0.0..=6.0).contains(&hla) || !truthy(get(&end, "hlaFold")) {
            out.push(format!(
                "end: the HLA card shows {hla} lines, the absent ones not folded ({})",
                truthy(get(&end, "hlaFold"))
            ));
        }
        if !truthy(get(&end, "margin")) {
            out.push("end: the HLA card's caveat is not in its margin".into());
        }
        out
    }

    /// DEMO-002's page at 390 px: its height, its findings to act on, a deep link to a card two
    /// folds deep and one to a folded section, each before and after the link.
    fn phone_report_problems(&mut self) -> Vec<String> {
        self.viewport(390, 844, 2);
        self.go(&format!("/runs/{PARTNER}/"));
        self.wait("document.documentElement.dataset.pgpResults", 20);
        sleep(1000);
        let pre = self.jv(PHONE_JS);
        self.ev("location.hash = '#s/roh'; true");
        sleep(500);
        let card = self.jv(PHONE_JS);
        self.ev("location.hash = '#s7'; true");
        sleep(500);
        let section = self.jv(PHONE_JS);
        self.viewport(1440, 900, 1);
        let mut out = vec![];
        let h = num(get(&pre, "height"));
        if !(h > 0.0 && h < PHONE_MAX) {
            out.push(format!(
                "the page is {h} CSS px tall at 390 px, want under {PHONE_MAX}"
            ));
        }
        if !truthy(get(&pre, "loud")) || truthy(get(&pre, "hidden")) {
            out.push(format!("findings to act on, for a prescriber, for the family or worth a look closed behind a tap: {:?}", strs(get(&pre, "hidden"))));
        }
        if truthy(get(&pre, "deep")) || truthy(get(&pre, "sec")) {
            out.push("card roh or section 7 starts open on a phone".into());
        }
        if !truthy(get(&card, "deep")) {
            out.push("#s/roh leaves its card closed".into());
        }
        if !truthy(get(&section, "sec")) {
            out.push("#s7 leaves section 7 closed".into());
        }
        out
    }

    /// DEMO-002's page at 390 px with its section 5 card FLAGGED flagged.
    fn fold_problems(&mut self) -> Vec<String> {
        self.viewport(390, 844, 2);
        self.go(&format!("/runs/{PARTNER}/"));
        self.wait("document.documentElement.dataset.pgpResults && document.querySelector('#results details.pfold')", 20);
        sleep(500);
        let flagged = truthy(&self.ev(&format!(
            "var e = document.getElementById('card-{FLAGGED}'); !!e && e.checkVisibility() && e.getBoundingClientRect().height > 0"
        )));
        let opened = st(&self.ev(r#"(function () {
          var s6 = document.querySelector('#s6 details.secfold'), roh = document.querySelector('details[data-card="fold-roh"]');
          var more = document.querySelector('#results article.result details.lines-more');
          if (!s6 || !roh || !more) { return ''; }
          s6.open = true; roh.open = true;
          for (var d = more; d; d = d.parentElement && d.parentElement.closest('details')) { d.open = true; }
          return more.closest('article').id; }())"#));
        self.viewport(1280, 844, 2);
        self.wait(
            "!document.querySelector('#results details.pfold') && 'wide'",
            10,
        );
        sleep(500);
        let wide = truthy(&self.ev(&format!(
            "var a = document.getElementById({}); !!a && !!a.querySelector('details.lines-more[open]')",
            js_str(&opened)
        )));
        self.viewport(390, 844, 2);
        self.wait(
            "document.querySelector('#results details.pfold') && 'phone'",
            10,
        );
        sleep(500);
        let back = self.jv(&format!(
            r#"JSON.stringify([!!document.querySelector('#s6 details.secfold[open]'),
          !!document.querySelector('details[data-card="fold-roh"][open]'),
          !!document.querySelector('#' + {} + ' details.lines-more[open]')])"#,
            js_str(if opened.is_empty() { "x" } else { &opened })
        ));
        self.ev("location.hash = '#s8'; true");
        self.wait(
            "document.querySelector('#s8 details.secfold[open]') && 'open'",
            5,
        );
        self.ev(r##"window.addEventListener('click', function (e) { e.preventDefault(); });
          document.querySelector('#s8 details.secfold').open = false;
          document.querySelector('#tabbar a[href="#s8"]').dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true, ctrlKey: true}));
          true"##);
        sleep(300);
        let ctrl = self.ev("document.querySelector('#s8 details.secfold').open");
        self.ev(r##"document.querySelector('#tabbar a[href="#s8"]').dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true})); true"##);
        sleep(300);
        let plain = self.ev("document.querySelector('#s8 details.secfold').open");
        self.viewport(1440, 900, 1);
        let mut out = vec![];
        if !flagged {
            out.push(format!("the flagged section 5 card {FLAGGED} with no finding of its own is folded on a phone"));
        }
        if opened.is_empty() {
            out.push("no section 6 fold, ROH fold or card with other lines to open".into());
        }
        if !wide {
            out.push(format!(
                "the other lines opened in {opened} closed when the page turned wide"
            ));
        }
        if strs(&back) != ["true", "true", "true"] {
            out.push(format!("section 6, the ROH card and the other lines opened on the phone, open after wide and back: {}", show(&back)));
        }
        if ctrl != Value::Bool(false) {
            out.push(format!(
                "a Ctrl+click on a link to the current place opened its fold: {}",
                show(&ctrl)
            ));
        }
        if plain != Value::Bool(true) {
            out.push(format!(
                "a plain click on a link to the current place left its fold closed: {}",
                show(&plain)
            ));
        }
        out
    }
}

/// The person page answers "is there anything serious?" first, and on a phone it is short: the
/// findings that matter open, the rest one tap away, and folds the reader opened stay open.
#[test]
fn person_page_layout_and_phone_folds() {
    let Some(mut d) = Demo::new("person-layout") else {
        return;
    };
    let live = d.site.join("assets/live.js");
    let person = d.site.join("assets/person.js");
    let bad = d.layout_problems();
    d.r.check("browser: the person page leads with the report, a short run band while it runs, the step timeline last and folded, a short wire, and each card's caveat in its margin with its quiet lines folded", bad);
    for (file, a, z, desc) in [
        (
            &live,
            "said.slice(0, WIRE_SHOWN)",
            "said",
            "a wire that shows every line",
        ),
        (
            &live,
            "if (m && hl[k].tone !== 'attention')",
            "if (false)",
            "a wire that lists every plain finish",
        ),
        (
            &person,
            "var LINES_SHOWN = 6;",
            "var LINES_SHOWN = 99;",
            "a card that shows every line",
        ),
        (
            &person,
            "folded = true; fold.open = ",
            "folded = true; fold.open = true || ",
            "a step timeline left open",
        ),
    ] {
        let _p = plant(file, &[(a, z)]);
        let bad = d.layout_problems();
        d.r.control(desc, bad);
    }

    let bad = d.phone_report_problems();
    d.r.check(&format!("browser: at 390 px {PARTNER}'s page is under {PHONE_MAX} CSS px with every finding to act on, for a prescriber, for the family or worth a look open, and a link to a folded card or section opens it"), bad);
    d.go(&format!("/runs/{PARTNER}/"));
    d.wait("document.documentElement.dataset.pgpResults", 20);
    sleep(500);
    let folds = d.ev("document.querySelectorAll('#results details.pfold').length");
    d.r.ok(
        "browser: wider than a phone the person page folds no card and no section",
        folds == Value::Int(0),
        &folds,
    );
    for (a, z, desc) in [
        (
            "var phoneQuery = window.matchMedia ?",
            "var phoneQuery = false ?",
            "a phone page with every card open",
        ),
        (
            "d = d.parentElement && d.parentElement.closest('details')",
            "d = null",
            "a link that opens only the nearest fold",
        ),
    ] {
        let _p = plant(&person, &[(a, z)]);
        let bad = d.phone_report_problems();
        d.r.control(desc, bad);
    }

    let rj = d.site.join(format!("runs/{PARTNER}/results.json"));
    let keep = Keep::new(&rj);
    let mut res = loads(&keep.text()).unwrap();
    let mut found = 0;
    for ch in arr(at(&mut res, &["chapters"])) {
        for c in arr(at(ch, &["sections"])) {
            if st(get(c, "id")) == FLAGGED {
                if !truthy(get(c, "findings")) {
                    found += 1;
                }
                obj(c).insert("attention".into(), Value::Bool(true));
            }
        }
    }
    let sec5 = items(get(&res, "journal"))
        .iter()
        .find(|x| num(get(x, "n")) == 5.0)
        .cloned();
    let in5 = sec5
        .as_ref()
        .is_some_and(|x| strs(get(x, "cards")).iter().any(|c| c == FLAGGED));
    d.r.ok(&format!("browser: {FLAGGED} is a section 5 card with no finding of its own, so it can be flagged alone"), found == 1 && in5, (found, in5));
    fs::write(&rj, dumps(&res, &Dump::default()).unwrap()).unwrap();
    let bad = d.fold_problems();
    d.r.check("browser: at 390 px a flagged section 5 card stays open, folds opened on a phone stay open after the page turns wide and back (other lines too), and a Ctrl+click leaves a fold as it is", bad);
    for (a, z, desc) in [
        (
            "function flaggedAlone(c) { return ",
            "function flaggedAlone(c) { return false && ",
            "a flagged section 5 card folded with its section",
        ),
        (
            "function restoreFolds() { foldKeys()",
            "function restoreFolds() { return; foldKeys()",
            "a redraw that forgets the folds",
        ),
        (
            "if (e.button !== 0 || e.ctrlKey",
            "if (false && e.ctrlKey",
            "a Ctrl+click that opens the fold",
        ),
    ] {
        let _p = plant(&person, &[(a, z)]);
        let bad = d.fold_problems();
        d.r.control(desc, bad);
    }
    drop(keep);
    d.r.done();
}

// ---------------------------------------------------------------------------------------------
// The live site's genome viewer (`og live serve --viewer-fasta --viewer-igv-js`).

const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03, 0, 0,
    0, 0, 0, 0, 0, 0, 0,
];

fn bgzf_block(data: &[u8]) -> Vec<u8> {
    let mut e = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::new(6));
    e.write_all(data).unwrap();
    let body = e.finish().unwrap();
    let mut crc = flate2::Crc::new();
    crc.update(data);
    let mut out = vec![0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff];
    out.extend(6u16.to_le_bytes());
    out.extend(b"BC");
    out.extend(2u16.to_le_bytes());
    out.extend(u16::try_from(body.len() + 25).unwrap().to_le_bytes());
    out.extend(&body);
    out.extend(crc.sum().to_le_bytes());
    out.extend(u32::try_from(data.len()).unwrap().to_le_bytes());
    out
}

fn reg2bin(beg: i64, end: i64) -> u32 {
    let end = end - 1;
    for (shift, off) in [(14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)] {
        if beg >> shift == end >> shift {
            return (off + (beg >> shift)) as u32;
        }
    }
    0
}

/// The per-reference part of a BAI or TBI: every record (ref, beg, end) sits in one BGZF block
/// that starts at `rec_voff`.
fn binning_index(
    n_ref: usize,
    recs: &[(usize, i64, i64)],
    rec_voff: i64,
    end_voff: i64,
) -> Vec<u8> {
    let mut out = Vec::new();
    for r in 0..n_ref {
        let mine: Vec<_> = recs.iter().filter(|x| x.0 == r).collect();
        if mine.is_empty() {
            out.extend(0i32.to_le_bytes());
            out.extend(0i32.to_le_bytes());
            continue;
        }
        let mut bins: Vec<u32> = mine.iter().map(|&&(_, b, e)| reg2bin(b, e)).collect();
        bins.sort_unstable();
        bins.dedup();
        out.extend((bins.len() as i32).to_le_bytes());
        for b in bins {
            out.extend(b.to_le_bytes());
            out.extend(1i32.to_le_bytes());
            out.extend(rec_voff.to_le_bytes());
            out.extend(end_voff.to_le_bytes());
        }
        let last = mine.iter().map(|x| x.2).max().unwrap();
        let n_intv = ((last - 1) >> 14) + 1;
        out.extend((n_intv as i32).to_le_bytes());
        for _ in 0..n_intv {
            out.extend(rec_voff.to_le_bytes());
        }
    }
    out
}

/// A BAM of `reads` (ref index, 0-based pos, sequence), sorted, 50M each, with its .bai.
fn write_bam(path: &Path, contigs: &[(String, String)], reads: &[(usize, i64, String)]) {
    let mut text = "@HD\tVN:1.6\tSO:coordinate\n".to_owned();
    for (n, s) in contigs {
        text.push_str(&format!("@SQ\tSN:{n}\tLN:{}\n", s.len()));
    }
    text.push_str("@RG\tID:S\tSM:S\n");
    let mut head = b"BAM\x01".to_vec();
    head.extend((text.len() as i32).to_le_bytes());
    head.extend(text.as_bytes());
    head.extend((contigs.len() as i32).to_le_bytes());
    for (n, s) in contigs {
        head.extend((n.len() as i32 + 1).to_le_bytes());
        head.extend(n.as_bytes());
        head.push(0);
        head.extend((s.len() as i32).to_le_bytes());
    }
    let code = |c: u8| "=ACMGRSVTWYHKDBN".bytes().position(|x| x == c).unwrap() as u8;
    let (mut body, mut recs) = (Vec::new(), Vec::new());
    for (i, (r, pos, seq)) in reads.iter().enumerate() {
        let name = format!("r{i}\0");
        let sb = seq.as_bytes();
        let packed: Vec<u8> = sb
            .chunks(2)
            .map(|c| (code(c[0]) << 4) | c.get(1).map_or(0, |&x| code(x)))
            .collect();
        let len = sb.len() as i64;
        let mut rec = Vec::new();
        rec.extend((*r as i32).to_le_bytes());
        rec.extend((*pos as i32).to_le_bytes());
        rec.push(name.len() as u8);
        rec.push(60);
        rec.extend((reg2bin(*pos, pos + len) as u16).to_le_bytes());
        rec.extend(1u16.to_le_bytes());
        rec.extend(0u16.to_le_bytes());
        rec.extend((len as i32).to_le_bytes());
        for v in [-1i32, -1, 0] {
            rec.extend(v.to_le_bytes());
        }
        rec.extend(name.as_bytes());
        rec.extend(((len as u32) << 4).to_le_bytes());
        rec.extend(&packed);
        rec.extend(std::iter::repeat_n(30u8, sb.len()));
        rec.extend(b"RGZS\0");
        body.extend((rec.len() as i32).to_le_bytes());
        body.extend(rec);
        recs.push((*r, *pos, pos + len));
    }
    let (b0, b1) = (bgzf_block(&head), bgzf_block(&body));
    let mut all = b0.clone();
    all.extend(&b1);
    all.extend(BGZF_EOF);
    fs::write(path, all).unwrap();
    let mut bai = b"BAI\x01".to_vec();
    bai.extend((contigs.len() as i32).to_le_bytes());
    bai.extend(binning_index(
        contigs.len(),
        &recs,
        (b0.len() as i64) << 16,
        ((b0.len() + b1.len()) as i64) << 16,
    ));
    fs::write(format!("{}.bai", path.display()), bai).unwrap();
}

/// A bgzipped VCF of `records` (chrom, 1-based pos, ref, alt, info), with a .tbi when `index`.
fn write_vcf(
    path: &Path,
    contigs: &[(String, String)],
    records: &[(&str, i64, String, String, &str)],
    index: bool,
) {
    let mut head = "##fileformat=VCFv4.2\n".to_owned();
    for (n, s) in contigs {
        head.push_str(&format!("##contig=<ID={n},length={}>\n", s.len()));
    }
    head.push_str("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">\n");
    head.push_str("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n");
    head.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    head.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n");
    let body: String = records
        .iter()
        .map(|(c, p, r, a, i)| format!("{c}\t{p}\t.\t{r}\t{a}\t50\tPASS\t{i}\tGT\t0/1\n"))
        .collect();
    let (b0, b1) = (bgzf_block(head.as_bytes()), bgzf_block(body.as_bytes()));
    let mut all = b0.clone();
    all.extend(&b1);
    all.extend(BGZF_EOF);
    fs::write(path, all).unwrap();
    if index {
        let order: Vec<&str> = contigs.iter().map(|(n, _)| n.as_str()).collect();
        let mut names: Vec<&str> = records.iter().map(|r| r.0).collect();
        names.sort_by_key(|n| order.iter().position(|o| o == n));
        names.dedup();
        let nm: Vec<u8> = names.iter().flat_map(|n| n.bytes().chain([0])).collect();
        let recs: Vec<_> = records
            .iter()
            .map(|(c, p, r, _, _)| {
                (
                    names.iter().position(|n| n == c).unwrap(),
                    p - 1,
                    p - 1 + r.len() as i64,
                )
            })
            .collect();
        let mut tbi = b"TBI\x01".to_vec();
        for v in [
            names.len() as i32,
            2,
            1,
            2,
            0,
            i32::from(b'#'),
            0,
            nm.len() as i32,
        ] {
            tbi.extend(v.to_le_bytes());
        }
        tbi.extend(&nm);
        tbi.extend(binning_index(
            names.len(),
            &recs,
            (b0.len() as i64) << 16,
            ((b0.len() + b1.len()) as i64) << 16,
        ));
        let mut out = bgzf_block(&tbi);
        out.extend(BGZF_EOF);
        fs::write(format!("{}.tbi", path.display()), out).unwrap();
    }
}

/// A random sequence, the same on every run.
fn sequence(seed: &mut u64, n: usize) -> String {
    (0..n)
        .map(|_| {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(*seed >> 62) as usize] as char
        })
        .collect()
}

fn write_fasta(fa: &Path, contigs: &[(String, String)]) {
    let (mut fasta, mut fai, mut off) = (String::new(), String::new(), 0);
    for (n, s) in contigs {
        off += n.len() + 2;
        fai.push_str(&format!("{n}\t{}\t{off}\t60\t61\n", s.len()));
        let lines: Vec<&str> = s
            .as_bytes()
            .chunks(60)
            .map(|c| std::str::from_utf8(c).unwrap())
            .collect();
        fasta.push_str(&format!(">{n}\n{}\n", lines.join("\n")));
        off += s.len() + lines.len();
    }
    fs::write(fa, fasta).unwrap();
    fs::write(format!("{}.fai", fa.display()), fai).unwrap();
}

/// The genome folder `gd`: reference/synthetic.fasta (chr1 60 kb, so a 49 kb view holds reads),
/// and run S with its BAM, VCF, SV VCF and a report whose Viewer links point at the local viewer.
fn make_sample(gd: &Path, viewer_url: &str, viewer_note: &str) {
    let mut seed = 7;
    let mut chr1 = sequence(&mut seed, 3000);
    let chr2 = sequence(&mut seed, 1000);
    chr1.push_str(&"ACGT".repeat(14250));
    let contigs = vec![("chr1".to_owned(), chr1), ("chr2".to_owned(), chr2)];
    fs::create_dir_all(gd.join("reference")).unwrap();
    write_fasta(&gd.join("reference/synthetic.fasta"), &contigs);
    let s = gd.join("S");
    for d in ["aligned", "vcf", "sv_merged", "nextflow"] {
        fs::create_dir_all(s.join(d)).unwrap();
    }
    let seq = &contigs[0].1;
    let reads: Vec<_> = (400..700)
        .step_by(10)
        .map(|p| (0, p as i64, seq[p..p + 50].to_owned()))
        .collect();
    write_bam(&s.join("aligned/S_sorted.bam"), &contigs, &reads);
    let r = &seq[450..451];
    let alt = "ACGT".replace(r, "")[..1].to_owned();
    write_vcf(
        &s.join("vcf/S.vcf.gz"),
        &contigs,
        &[("chr1", 451, r.to_owned(), alt, ".")],
        true,
    );
    write_vcf(
        &s.join("sv_merged/S_sv_consensus.vcf.gz"),
        &contigs,
        &[(
            "chr1",
            1000,
            seq[999..1000].to_owned(),
            "<DEL>".into(),
            "SVTYPE=DEL;END=1500",
        )],
        false,
    );
    fs::write(
        s.join("S_report.html"),
        format!(
            "<!doctype html><html><body><table><tr><th>Where</th><th>Viewer</th></tr>\n\
             <tr><td>chr1:451</td><td><a href=\"{viewer_url}#chr1:451\">reads</a></td></tr></table>\n{viewer_note}\n</body></html>\n"
        ),
    )
    .unwrap();
}

/// The igv.js release versions.env pins: (version, sha256).
fn igv_pin() -> (String, String) {
    let env = fs::read_to_string(repo_root().join("versions.env")).unwrap();
    let pin = |k: &str| {
        env.lines()
            .find_map(|l| l.strip_prefix(&format!("{k}=\"")))
            .and_then(|v| v.strip_suffix('"'))
            .unwrap_or_else(|| panic!("versions.env sets no {k}"))
            .to_owned()
    };
    (pin("IGV_JS_VERSION"), pin("IGV_JS_SHA256"))
}

fn sha256(p: &Path) -> String {
    let out = Command::new("shasum")
        .args(["-a", "256"])
        .arg(p)
        .output()
        .or_else(|_| Command::new("sha256sum").arg(p).output())
        .expect("shasum or sha256sum");
    String::from_utf8_lossy(&out.stdout)
        .split_whitespace()
        .next()
        .unwrap_or_default()
        .to_owned()
}

struct Resp {
    code: u16,
    head: Vec<(String, String)>,
    body: Vec<u8>,
}

impl Resp {
    fn h(&self, k: &str) -> String {
        self.head
            .iter()
            .find(|(n, _)| n.eq_ignore_ascii_case(k))
            .map(|(_, v)| v.clone())
            .unwrap_or_default()
    }
    fn text(&self) -> String {
        String::from_utf8_lossy(&self.body).into_owned()
    }
}

fn parse_resp(raw: &[u8]) -> Resp {
    let split = raw
        .windows(4)
        .position(|w| w == b"\r\n\r\n")
        .map_or(raw.len(), |i| i + 4);
    let head = String::from_utf8_lossy(&raw[..split]);
    let mut lines = head.split("\r\n");
    let code = lines
        .next()
        .and_then(|l| l.split_whitespace().nth(1))
        .and_then(|c| c.parse().ok())
        .unwrap_or(0);
    let head = lines
        .filter_map(|l| l.split_once(':'))
        .map(|(k, v)| (k.trim().to_owned(), v.trim().to_owned()))
        .collect();
    Resp {
        code,
        head,
        body: raw[split..].to_vec(),
    }
}

/// One HTTP/1.1 request to 127.0.0.1:`port`, the connection closed after the answer.
fn request(port: u16, method: &str, path: &str, headers: &[(&str, &str)], body: &[u8]) -> Resp {
    let mut c = TcpStream::connect(("127.0.0.1", port)).unwrap();
    c.set_read_timeout(Some(Duration::from_secs(20))).unwrap();
    let mut req =
        format!("{method} {path} HTTP/1.1\r\nHost: 127.0.0.1:{port}\r\nConnection: close\r\n");
    for (k, v) in headers {
        req.push_str(&format!("{k}: {v}\r\n"));
    }
    if method == "POST" {
        req.push_str(&format!("Content-Length: {}\r\n", body.len()));
    }
    req.push_str("\r\n");
    c.write_all(req.as_bytes()).unwrap();
    c.write_all(body).unwrap();
    let mut raw = vec![];
    let _ = c.read_to_end(&mut raw);
    parse_resp(&raw)
}

fn free_port() -> u16 {
    TcpListener::bind("127.0.0.1:0")
        .unwrap()
        .local_addr()
        .unwrap()
        .port()
}

/// `og live serve` on a free port, with the repository's assets and steps.
struct Site {
    child: Child,
    port: u16,
}

impl Site {
    fn start(og: &Path, args: &[std::ffi::OsString]) -> Site {
        let root = repo_root();
        let mut site = Site {
            child: Command::new(og)
                .args([
                    "live",
                    "serve",
                    "--port",
                    "0",
                    "--ingest-port",
                    "0",
                    "--quiet",
                    "--assets",
                ])
                .arg(root.join("assets/live"))
                .arg("--steps")
                .arg(root.join("assets/steps.tsv"))
                .args(args)
                .env_remove("OWNGENOME_LIVE_PASSWORD_HASH")
                .stdin(Stdio::null())
                .stdout(Stdio::null())
                .stderr(Stdio::piped())
                .spawn()
                .unwrap(),
            port: 0,
        };
        let mut lines = BufReader::new(site.child.stderr.take().unwrap()).lines();
        while let Some(line) = lines.next() {
            let line = line.unwrap();
            if let Some(url) = line.strip_prefix("og live serve: browser ") {
                site.port = url
                    .split_whitespace()
                    .next()
                    .and_then(|u| u.trim_end_matches('/').rsplit(':').next())
                    .and_then(|p| p.parse().ok())
                    .expect("og live serve printed no port");
                // Keep reading what the server says, so it never writes to a
                // closed pipe.
                std::thread::spawn(move || {
                    for line in lines.map_while(Result::ok) {
                        println!("og live serve: {line}");
                    }
                });
                return site;
            }
            println!("og live serve: {line}");
        }
        panic!("og live serve exited before listening");
    }
}

impl Drop for Site {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

/// A mistake a proxy in front of `og live serve` plants in the viewer's answers: the planted
/// mistakes of tests/test_live_viewer.py, made in the answer since the server is compiled.
#[derive(Clone, Copy, PartialEq, Debug)]
enum Mistake {
    Nothing,
    /// A viewer that serves any file of the run folder.
    AnyFile,
    /// A report filter that does nothing.
    NoFilter,
    /// A viewer that compares the still-quoted file name.
    Quoted,
    /// A viewer that serves a track index symlinked out of GENOME_DIR.
    IndexOut,
}

struct ProxyState {
    upstream: u16,
    mistake: Mistake,
    block_frames: bool,
    csp: String,
    gd: PathBuf,
    map_doc: Vec<u8>,
}

/// A proxy on its own port: `/map-1g/` (the map's reads fixture page) and `/map-1g.json` (its gzip
/// document); every other request goes to `og live serve` (`upstream`), with `mistake` planted.
struct Proxy {
    port: u16,
    state: Arc<Mutex<ProxyState>>,
}

const MAP_1G_PAGE: &str = r#"<!doctype html><html data-theme="dark"><head><meta name="viewport" content="width=device-width, initial-scale=1"><link rel="stylesheet" href="/assets/live.css"><link rel="stylesheet" href="/assets/map.css"><script src="/assets/core.js" defer></script><script src="/assets/map.js" defer></script></head><body></body></html>"#;

const MAP_1G_DOC: &str = r##"{"v": 1, "people": [{"id": "S"}], "dict": {"veil": ["closed"]},
 "ref": {"chroms": [{"n": "chr1", "len": 60000, "masked": []}], "bands": {}},
 "per": {"S": {"findings": {
   "id": ["open", "closed"], "level": ["review", "act"], "gate": ["", ""],
   "chrom": ["chr1", "chr1"], "pos": [451, 551], "s": [450, 550], "e": [451, 551],
   "gene": ["", ""], "title": ["Synthetic finding", "Closed finding"],
   "href": ["#open", "#closed"], "phrase": ["Synthetic row", "Closed row"], "veil": [-1, 0]}}},
 "layers": {"findings": {"state": "landed"}, "roh": {"state": "not_run"}}, "text": {}}"##;

fn reply(c: &mut TcpStream, code: u16, ctype: &str, body: &[u8], extra: &str) {
    let _ = write!(
        c,
        "HTTP/1.1 {code} X\r\nContent-Type: {ctype}\r\nContent-Length: {}\r\n{extra}Connection: close\r\n\r\n",
        body.len()
    );
    let _ = c.write_all(body);
}

fn find_named(dir: &Path, name: &str) -> Option<PathBuf> {
    for e in fs::read_dir(dir).ok()?.flatten() {
        let p = e.path();
        if p.is_dir() {
            if let Some(f) = find_named(&p, name) {
                return Some(f);
            }
        } else if e.file_name().to_string_lossy() == name && p.is_file() {
            return Some(p);
        }
    }
    None
}

impl Proxy {
    fn start(upstream: u16, gd: &Path) -> Proxy {
        let csp = request(upstream, "GET", "/", &[], b"").h("content-security-policy");
        assert!(
            csp.contains("frame-src 'self'; "),
            "the live policy has no frame-src 'self': {csp}"
        );
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gz.write_all(MAP_1G_DOC.as_bytes()).unwrap();
        let state = Arc::new(Mutex::new(ProxyState {
            upstream,
            mistake: Mistake::Nothing,
            block_frames: false,
            csp,
            gd: gd.to_owned(),
            map_doc: gz.finish().unwrap(),
        }));
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let port = listener.local_addr().unwrap().port();
        let shared = Arc::clone(&state);
        thread::spawn(move || {
            for c in listener.incoming().flatten() {
                let shared = Arc::clone(&shared);
                thread::spawn(move || Proxy::answer(c, &shared));
            }
        });
        Proxy { port, state }
    }

    fn set(&self, f: impl FnOnce(&mut ProxyState)) {
        f(&mut self.state.lock().unwrap());
    }

    /// The request (`head` without its blank line, then `body`) sent on to `og live serve` as
    /// 127.0.0.1:`upstream`, the connection closed after the answer.
    fn send(head: &str, body: &[u8], upstream: u16) -> Option<TcpStream> {
        let mut up = TcpStream::connect(("127.0.0.1", upstream)).ok()?;
        let mut out = String::new();
        for (i, l) in head.split("\r\n").enumerate() {
            let lower = l.to_ascii_lowercase();
            if i > 0 && (lower.starts_with("host:") || lower.starts_with("connection:")) {
                continue;
            }
            out.push_str(l);
            out.push_str("\r\n");
        }
        out.push_str(&format!(
            "Host: 127.0.0.1:{upstream}\r\nConnection: close\r\n\r\n"
        ));
        up.write_all(out.as_bytes()).ok()?;
        up.write_all(body).ok()?;
        Some(up)
    }

    fn forward(head: &str, body: &[u8], upstream: u16) -> Vec<u8> {
        let mut raw = vec![];
        if let Some(mut up) = Proxy::send(head, body, upstream) {
            let _ = up.read_to_end(&mut raw);
        }
        raw
    }

    fn answer(mut c: TcpStream, shared: &Mutex<ProxyState>) {
        let _ = c.set_read_timeout(Some(Duration::from_secs(30)));
        let mut raw = vec![];
        let mut buf = [0u8; 8192];
        let end = loop {
            if let Some(i) = raw.windows(4).position(|w| w == b"\r\n\r\n") {
                break i + 4;
            }
            match c.read(&mut buf) {
                Ok(0) | Err(_) => return,
                Ok(n) => raw.extend(&buf[..n]),
            }
        };
        let head = String::from_utf8_lossy(&raw[..end]).into_owned();
        let len: usize = head
            .lines()
            .find_map(|l| {
                l.to_ascii_lowercase()
                    .strip_prefix("content-length:")
                    .map(|v| v.trim().parse().unwrap_or(0))
            })
            .unwrap_or(0);
        while raw.len() < end + len {
            match c.read(&mut buf) {
                Ok(0) | Err(_) => break,
                Ok(n) => raw.extend(&buf[..n]),
            }
        }
        let target = head.split_whitespace().nth(1).unwrap_or("/").to_owned();
        let path = target.split(['?', '#']).next().unwrap_or("/").to_owned();
        let (upstream, mistake, csp, gd, doc) = {
            let g = shared.lock().unwrap();
            let csp = if g.block_frames {
                g.csp.replace("frame-src 'self'; ", "")
            } else {
                g.csp.clone()
            };
            (g.upstream, g.mistake, csp, g.gd.clone(), g.map_doc.clone())
        };
        if path == "/map-1g/" {
            return reply(
                &mut c,
                200,
                "text/html; charset=utf-8",
                MAP_1G_PAGE.as_bytes(),
                &format!("Content-Security-Policy: {csp}\r\nCache-Control: no-store\r\n"),
            );
        }
        if path == "/map-1g.json" {
            return reply(
                &mut c,
                200,
                "application/gzip",
                &doc,
                "Cache-Control: no-store\r\n",
            );
        }
        let parts: Vec<&str> = path.trim_start_matches('/').splitn(4, '/').collect();
        let viewer_rest = (parts.len() == 4 && parts[0] == "runs" && parts[2] == "viewer")
            .then(|| (parts[1], parts[3]));
        match (mistake, viewer_rest) {
            // Any file inside the genome folder: a path that leaves it, by
            // `..` or by a symlink, is not served.
            (Mistake::AnyFile, Some((run, rest)))
                if !rest.is_empty()
                    && !format!("{run}/{rest}")
                        .split('/')
                        .any(|p| p.is_empty() || p == "..")
                    && gd.join(run).join(rest).is_file()
                    && fs::canonicalize(gd.join(run).join(rest))
                        .is_ok_and(|p| p.starts_with(fs::canonicalize(&gd).unwrap())) =>
            {
                let body = fs::read(gd.join(run).join(rest)).unwrap();
                return reply(&mut c, 200, "application/octet-stream", &body, "");
            }
            (Mistake::Quoted, Some((_, rest)))
                if rest.starts_with("data/") && rest.contains('%') =>
            {
                return reply(&mut c, 404, "text/plain", b"not found\n", "");
            }
            (Mistake::IndexOut, Some((run, rest))) if rest.starts_with("data/") => {
                let answer = Proxy::forward(&head[..end - 4], &raw[end..], upstream);
                let name = rest.rsplit('/').next().unwrap_or_default();
                // Only a file that lies outside GENOME_DIR (a symlink out of it), as the planted
                // viewer kept such a track index in its list of files.
                let root = gd.canonicalize().unwrap_or_default();
                if parse_resp(&answer).code == 404
                    && let Some(f) = find_named(&gd.join(run), name)
                    && !f.canonicalize().unwrap_or_default().starts_with(&root)
                {
                    return reply(
                        &mut c,
                        200,
                        "application/octet-stream",
                        &fs::read(f).unwrap(),
                        "",
                    );
                }
                let _ = c.write_all(&answer);
                return;
            }
            _ => {}
        }
        if mistake == Mistake::NoFilter
            && parts.len() == 3
            && parts[0] == "runs"
            && parts[2] == "report.html"
        {
            let body = fs::read(gd.join(parts[1]).join(format!("{}_report.html", parts[1])))
                .unwrap_or_default();
            return reply(&mut c, 200, "text/html; charset=utf-8", &body, "");
        }
        let Some(mut up) = Proxy::send(&head[..end - 4], &raw[end..], upstream) else {
            return reply(&mut c, 502, "text/plain", b"no upstream\n", "");
        };
        let _ = std::io::copy(&mut up, &mut c);
        let _ = c.shutdown(Shutdown::Both);
    }
}

/// The paths the viewer, the report and the outside checks read, run with `log` collecting
/// each check's outcome: (description, problems).
struct Viewer<'a> {
    port: u16,
    auth: Vec<(&'a str, String)>,
    gd: &'a Path,
    viewer_csp: &'a str,
    viewer_url: &'a str,
    viewer_note: &'a str,
}

impl Viewer<'_> {
    fn get(&self, path: &str, extra: &[(&str, &str)]) -> Resp {
        let mut h: Vec<(&str, &str)> = self.auth.iter().map(|(k, v)| (*k, v.as_str())).collect();
        h.extend_from_slice(extra);
        request(self.port, "GET", path, &h, b"")
    }

    /// Check 2: what the viewer serves to a signed-in browser.
    fn viewer_checks(&self, out: &mut Vec<(String, Vec<String>)>) {
        let mut ck = |d: &str, ok: bool, detail: String| {
            out.push((d.to_owned(), if ok { vec![] } else { vec![detail] }))
        };
        let r = self.get("/runs/S/viewer/", &[]);
        let csp = r.h("content-security-policy");
        let page = r.text();
        ck(
            "page: 200 with the viewer's policy, served by the site",
            r.code == 200
                && csp == self.viewer_csp
                && page.contains("served by the run site")
                && !page.contains("served from this computer only"),
            format!("{} {csp:?}", r.code),
        );
        let igv = self.get("/runs/S/viewer/igv.min.js", &[]);
        let js = self.get("/runs/S/viewer/viewer.js", &[]);
        ck(
            "page: igv.min.js and viewer.js are served beside it",
            igv.code == 200
                && igv.body.len() > 100_000
                && js.code == 200
                && js.text().contains("config.json"),
            format!("{} {}", igv.code, js.code),
        );
        let r = self.get("/runs/S/viewer/config.json", &[]);
        let cfg = if r.code == 200 {
            loads(&r.text()).unwrap_or(Value::Null)
        } else {
            Value::Null
        };
        let tracks = items(get(&cfg, "tracks"));
        let names: Vec<String> = tracks
            .iter()
            .map(|t| {
                st(get(t, "url"))
                    .rsplit('/')
                    .next()
                    .unwrap_or_default()
                    .to_owned()
            })
            .collect();
        ck(
            "config: the run's BAM, VCF and SV VCF, by relative URL",
            names == ["S_sorted.bam", "S.vcf.gz", "S_sv_consensus.vcf.gz"]
                && tracks
                    .iter()
                    .all(|t| st(get(t, "url")).starts_with("data/")),
            format!("{} {names:?}", r.code),
        );
        let raw = fs::read(self.gd.join("S/aligned/S_sorted.bam")).unwrap();
        let bam = tracks
            .first()
            .map_or("/runs/S/viewer/data/x".to_owned(), |t| {
                format!("/runs/S/viewer/{}", st(get(t, "url")))
            });
        let r = self.get(&bam, &[("Range", "bytes=10-109")]);
        ck(
            "reads: a BAM range is 206 with exactly those bytes",
            r.code == 206
                && r.body == raw[10..110]
                && r.h("content-range") == format!("bytes 10-109/{}", raw.len()),
            format!("{} {:?}", r.code, r.h("content-range")),
        );
        let past = format!("bytes={}-", raw.len() + 5);
        let r = self.get(&bam, &[("Range", &past)]);
        ck(
            "reads: a range past the end is 416",
            r.code == 416,
            format!("{}", r.code),
        );
        let bad: Vec<u16> = [
            "/runs/S/viewer/data/2/S_report.html",
            "/runs/S/viewer/data/99/S_sorted.bam",
            "/runs/S/viewer/S_report.html",
            "/runs/S/viewer/data/2/..%2F..%2FS_report.html",
            "/runs/S/viewer/serve.py",
            "/runs/NOPE/viewer/config.json",
        ]
        .iter()
        .map(|p| self.get(p, &[]).code)
        .collect();
        ck(
            "a file the config does not name, and a run that does not exist, are 404",
            bad.iter().all(|c| *c == 404),
            format!("{bad:?}"),
        );
    }

    /// Check 3: the served report's links point at the site's viewer.
    fn report_checks(&self, out: &mut Vec<(String, Vec<String>)>) {
        let r = self.get("/runs/S/report.html", &[]);
        let page = r.text();
        let ok = r.code == 200
            && page.contains("href=\"/runs/S/viewer/#chr1:451\"")
            && !page.contains(self.viewer_url)
            && page.contains("this site's genome viewer")
            && !page.contains(self.viewer_note);
        out.push(("report: positions link to /runs/S/viewer/ with the site's note, none to 127.0.0.1:8765".into(),
            if ok { vec![] } else { vec![format!("{} {:?}", r.code, page.lines().filter(|l| l.contains("viewer")).take(4).collect::<Vec<_>>())] }));
    }

    /// Check 5: files out of GENOME_DIR are left out; an odd reference name is served.
    fn outside_checks(&self, out: &mut Vec<(String, Vec<String>)>) {
        let mut ck = |d: &str, ok: bool, detail: String| {
            out.push((d.to_owned(), if ok { vec![] } else { vec![detail] }))
        };
        let r = self.get("/runs/T/viewer/config.json", &[]);
        let cfg = loads(&r.text()).unwrap_or(Value::Null);
        ck(
            "a BAM symlinked out of GENOME_DIR is left out of the config",
            r.code == 200 && get(&cfg, "tracks") == &Value::Array(vec![]),
            format!("{} {}", r.code, show(&cfg)),
        );
        let r = self.get("/runs/U/viewer/config.json", &[]);
        let cfg = if r.code == 200 {
            loads(&r.text()).unwrap_or(Value::Null)
        } else {
            Value::Null
        };
        let c2 = self
            .get("/runs/U/viewer/data/2/U.vcf.gz", &[("Range", "bytes=0-9")])
            .code;
        let c3 = self
            .get(
                "/runs/U/viewer/data/3/U.vcf.gz.tbi",
                &[("Range", "bytes=0-9")],
            )
            .code;
        ck(
            "a VCF whose index is symlinked out of GENOME_DIR: no track, neither file served",
            r.code == 200 && get(&cfg, "tracks") == &Value::Array(vec![]) && c2 == 404 && c3 == 404,
            format!("{} {} {c2} {c3}", r.code, show(get(&cfg, "tracks"))),
        );
        let url = st(get(get(&cfg, "reference"), "fastaURL"));
        let r = self.get(&format!("/runs/U/viewer/{url}"), &[("Range", "bytes=0-4")]);
        ck(
            "a reference named with a space, '#' or '?' is served from its encoded URL",
            url == "data/0/my%20ref%231%3F.fasta" && r.code == 206 && r.body == b">chr1",
            format!("{url} {}", r.code),
        );
    }

    fn all(&self) -> Vec<(String, Vec<String>)> {
        let mut out = vec![];
        self.viewer_checks(&mut out);
        self.report_checks(&mut out);
        self.outside_checks(&mut out);
        out
    }
}

const PASSWORD: &str = "synthetic test password 7";
/// `og live serve hash-password` of PASSWORD, with a fixed salt.
const PASSWORD_HASH: &str = "pbkdf2_sha256$1000$00112233445566778899aabbccddeeff$2d0f8cdad5ceae90d47e704316f775d61bd586f3c81c2a18614c665f0cbc80b2";

fn login(port: u16) -> Option<String> {
    let r = request(
        port,
        "POST",
        "/login",
        &[("Content-Type", "application/x-www-form-urlencoded")],
        format!("password={}&next=%2F", PASSWORD.replace(' ', "+")).as_bytes(),
    );
    let cookie = r.h("set-cookie");
    let value = cookie
        .strip_prefix("owngenome_live=")?
        .split(';')
        .next()?
        .to_owned();
    (r.code == 303).then(|| format!("owngenome_live={value}"))
}

/// A process stopped and reaped when dropped, also when a check panics.
struct Reap(Child);

impl Drop for Reap {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}

fn wait_port(port: u16, child: &mut Child, secs: u64) -> bool {
    let end = Instant::now() + Duration::from_secs(secs);
    while Instant::now() < end {
        if child.try_wait().ok().flatten().is_some() {
            return false;
        }
        if TcpStream::connect(("127.0.0.1", port)).is_ok() {
            return true;
        }
        sleep(500);
    }
    false
}

fn access_rows(path: &Path) -> Vec<String> {
    fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .filter(|r| r.contains("\tviewer\t"))
        .map(str::to_owned)
        .collect()
}

/// The genome viewer on the live site, on a synthetic sample: scripts/owngenome-live.sh downloads
/// the pinned igv.js and starts the site with the viewer behind its login (`og live serve`,
/// through bin/og); the viewer's page, policy, config and byte ranges; the report's links; the
/// access log; files out of GENOME_DIR; a site without the viewer; the page in Chrome; and the
/// map's reads panel. Checks 2, 3 and 5 also run against planted mistakes (the ones the removed
/// Python tests planted), made by a proxy, and must fail there.
#[test]
fn live_site_genome_viewer() {
    let Some(chrome) = chrome_or_skip("open the live site's genome viewer") else {
        return;
    };
    let Some(og) = og() else { return };
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "live-viewer");
    let work = scratch.path();
    let gd = work.join("genome");
    let viewer_url = source_const("rust/crates/og-render/src/html.rs", "VIEWER_URL");
    let viewer_note = source_const("rust/crates/og-render/src/notes.rs", "VIEWER_NOTE");
    let viewer_csp = source_const("rust/crates/og-live/src/viewer.rs", "LIVE_CSP");
    make_sample(&gd, &viewer_url, &viewer_note);
    let state = work.join("state");
    fs::create_dir_all(&state).unwrap();
    fs::write(state.join("password"), format!("{PASSWORD_HASH}\n")).unwrap();
    let mut r = Report::default();

    // 1. owngenome-live.sh, with a password.
    let (port, iport) = (free_port(), free_port());
    let log = work.join("site.log");
    let mut script = Reap(
        Command::new("bash")
            .arg(repo_root().join("scripts/owngenome-live.sh"))
            .arg("--quiet")
            .env("GENOME_DIR", &gd)
            .env("REF_FASTA", "reference/synthetic.fasta")
            .env("OWNGENOME_LIVE_PORT", port.to_string())
            .env("OWNGENOME_LIVE_INGEST_PORT", iport.to_string())
            .env("OWNGENOME_LIVE_STATE_DIR", &state)
            .env_remove("OWNGENOME_LIVE_PASSWORD_HASH")
            .stdin(Stdio::null())
            .stdout(fs::File::create(&log).unwrap())
            .stderr(fs::File::create(work.join("site.err")).unwrap())
            .spawn()
            .unwrap(),
    );
    let up = wait_port(port, &mut script.0, 120);
    let said = || {
        format!(
            "{}{}",
            fs::read_to_string(work.join("site.log")).unwrap_or_default(),
            fs::read_to_string(work.join("site.err")).unwrap_or_default()
        )
    };
    r.ok("owngenome-live.sh: the site starts", up, said());
    assert!(up, "the site did not start: {}", said());
    let (version, sha) = igv_pin();
    let igv = gd.join(format!("tools/igv.js-{version}/igv.min.js"));
    r.ok(
        "owngenome-live.sh: the igv.js pinned in versions.env is downloaded, with its sha256",
        igv.is_file() && sha256(&igv) == sha,
        said(),
    );
    let p1 = request(port, "GET", "/runs/S/viewer/", &[], b"");
    let c2 = request(port, "GET", "/runs/S/viewer/config.json", &[], b"").code;
    let c3 = request(
        port,
        "GET",
        "/runs/S/viewer/data/2/S_sorted.bam",
        &[("Range", "bytes=0-9")],
        b"",
    )
    .code;
    r.ok(
        "no login: the page sends the browser to /login, the config and the reads are 401",
        p1.code == 303
            && p1.h("location") == "/login?next=/runs/S/viewer/"
            && c2 == 401
            && c3 == 401,
        (p1.code, p1.h("location"), c2, c3),
    );
    let auth = login(port);
    r.ok(
        "login: the password opens a session",
        auth.is_some(),
        "no 303 with a cookie",
    );
    let auth = auth.unwrap_or_default();

    // 2 and 3: the viewer and the served report.
    let acc = state.join("access.log");
    let before = access_rows(&acc).len();
    let v = Viewer {
        port,
        auth: vec![("Cookie", auth.clone())],
        gd: &gd,
        viewer_csp: &viewer_csp,
        viewer_url: &viewer_url,
        viewer_note: &viewer_note,
    };
    let mut out = vec![];
    v.viewer_checks(&mut out);
    v.report_checks(&mut out);
    for (desc, bad) in out {
        r.check(&desc, bad);
    }
    let red = request(port, "GET", "/runs/S/viewer", &[("Cookie", &auth)], b"");
    r.ok(
        "/runs/S/viewer redirects to /runs/S/viewer/",
        red.code == 301 && red.h("location") == "/runs/S/viewer/",
        (red.code, red.h("location")),
    );

    // 4: the access log.
    sleep(500);
    let mine: Vec<String> = access_rows(&acc).into_iter().skip(before).collect();
    r.ok("access log: one row for the page, none for config.json, igv.js or the byte ranges (and one for the run that does not exist)",
        mine.len() == 2
            && mine.iter().filter(|x| x.ends_with("\tviewer\tS\t-\t200")).count() == 1
            && mine.iter().filter(|x| x.ends_with("\tviewer\t-\t-\t404")).count() == 1,
        &mine);
    drop(script);
    // The script's exec replaced bash with og live serve, so stopping it stopped the site.

    // 5: a BAM, and a VCF index, symlinked out of GENOME_DIR; an odd reference name.
    let (sd, td, ud) = (gd.join("S"), gd.join("T"), gd.join("U"));
    for d in [
        td.join("aligned"),
        td.join("nextflow"),
        ud.join("vcf"),
        ud.join("nextflow"),
    ] {
        fs::create_dir_all(d).unwrap();
    }
    fs::copy(sd.join("aligned/S_sorted.bam"), work.join("outside.bam")).unwrap();
    fs::copy(
        sd.join("aligned/S_sorted.bam.bai"),
        work.join("outside.bam.bai"),
    )
    .unwrap();
    std::os::unix::fs::symlink(work.join("outside.bam"), td.join("aligned/T_sorted.bam")).unwrap();
    std::os::unix::fs::symlink(
        work.join("outside.bam.bai"),
        td.join("aligned/T_sorted.bam.bai"),
    )
    .unwrap();
    fs::copy(sd.join("vcf/S.vcf.gz"), ud.join("vcf/U.vcf.gz")).unwrap();
    fs::copy(sd.join("vcf/S.vcf.gz.tbi"), work.join("outside.vcf.gz.tbi")).unwrap();
    std::os::unix::fs::symlink(work.join("outside.vcf.gz.tbi"), ud.join("vcf/U.vcf.gz.tbi"))
        .unwrap();
    let odd = gd.join("reference/my ref#1?.fasta");
    fs::copy(gd.join("reference/synthetic.fasta"), &odd).unwrap();
    fs::copy(
        gd.join("reference/synthetic.fasta.fai"),
        gd.join("reference/my ref#1?.fasta.fai"),
    )
    .unwrap();
    let site = Site::start(
        &og,
        &[
            "--root".into(),
            gd.clone().into(),
            "--state-dir".into(),
            work.join("state2").into(),
            "--viewer-fasta".into(),
            odd.clone().into(),
            "--viewer-igv-js".into(),
            igv.clone().into(),
        ],
    );
    let v = Viewer {
        port: site.port,
        auth: vec![],
        gd: &gd,
        viewer_csp: &viewer_csp,
        viewer_url: &viewer_url,
        viewer_note: &viewer_note,
    };
    let mut out = vec![];
    v.outside_checks(&mut out);
    for (desc, bad) in out {
        r.check(&desc, bad);
    }
    for path in [
        "/runs/S/viewer/",
        "/runs/S/viewer/viewer.js",
        "/runs/S/viewer/config.json",
        "/runs/S/viewer/data/2/S_sorted.bam",
    ] {
        let a = request(site.port, "GET", path, &[], b"");
        r.ok(
            &format!("same-origin viewer framing: {path}"),
            a.code == 200
                && a.h("content-security-policy")
                    .contains("frame-ancestors 'self'")
                && a.h("x-frame-options") == "SAMEORIGIN",
            (
                a.code,
                a.h("content-security-policy"),
                a.h("x-frame-options"),
            ),
        );
    }
    for path in [
        "/",
        "/login",
        "/assets/map.js",
        "/runs/S/",
        "/runs/S/state.json",
        "/runs/S/report.html",
        "/not-found",
    ] {
        let a = request(site.port, "GET", path, &[], b"");
        r.ok(
            &format!("non-viewer framing stays closed: {path}"),
            a.h("content-security-policy")
                .contains("frame-ancestors 'none'")
                && a.h("x-frame-options") == "DENY",
            (
                a.code,
                a.h("content-security-policy"),
                a.h("x-frame-options"),
            ),
        );
    }

    // 6: a site without the viewer flags.
    let bare = Site::start(
        &og,
        &[
            "--root".into(),
            gd.clone().into(),
            "--state-dir".into(),
            work.join("state3").into(),
        ],
    );
    let a = request(bare.port, "GET", "/runs/S/viewer/", &[], b"");
    let rep = request(bare.port, "GET", "/runs/S/report.html", &[], b"").text();
    r.ok(
        "without the flags: the viewer is 404 with the reason, the report keeps its local links",
        a.code == 404 && a.text().contains("not set up") && rep.contains(&viewer_url),
        (a.code, a.text()),
    );

    // The planted mistakes: checks 2, 3 and 5 pass through a proxy that changes nothing, and each
    // mistake the proxy plants makes one of them fail.
    let proxy = Proxy::start(site.port, &gd);
    let pv = Viewer {
        port: proxy.port,
        auth: vec![],
        gd: &gd,
        viewer_csp: &viewer_csp,
        viewer_url: &viewer_url,
        viewer_note: &viewer_note,
    };
    let through: Vec<String> = pv
        .all()
        .into_iter()
        .flat_map(|(d, b)| b.into_iter().map(move |x| format!("{d}: {x}")))
        .collect();
    r.check(
        "checks 2, 3 and 5 pass through the proxy that plants nothing",
        through,
    );
    for (mistake, desc) in [
        (
            Mistake::AnyFile,
            "a viewer that serves any file of the run folder",
        ),
        (Mistake::NoFilter, "a report filter that does nothing"),
        (
            Mistake::Quoted,
            "a viewer that compares the still-quoted name",
        ),
        (
            Mistake::IndexOut,
            "a viewer that serves an index out of GENOME_DIR",
        ),
    ] {
        proxy.set(|s| s.mistake = mistake);
        let bad: Vec<String> = pv
            .all()
            .into_iter()
            .flat_map(|(d, b)| b.into_iter().map(move |x| format!("{d}: {x}")))
            .collect();
        r.control(
            &format!("checks 2, 3 and 5 on the planted mistake '{desc}'"),
            bad,
        );
    }
    proxy.set(|s| s.mistake = Mistake::Nothing);

    // 7: the page from the site, in Chrome.
    let mut b = Browser::launch(&chrome, &work.join("chrome")).unwrap();
    b.go(&format!(
        "http://127.0.0.1:{}/runs/S/viewer/#chr1:451",
        site.port
    ))
    .unwrap();
    let ready = b.wait(
        "document.documentElement.dataset.pgpStatus==='ready' && document.documentElement.dataset.pgpTracks==='3'",
        Duration::from_secs(30),
    );
    r.ok(
        "browser: the page from the site starts igv.js with the run's 3 tracks",
        ready.is_some(),
        b.eval("document.documentElement.dataset.pgpStatus").ok(),
    );
    let blocked = b
        .eval("document.documentElement.dataset.pgpBlocked || ''")
        .unwrap_or(Value::Null);
    r.ok(
        "browser: the policy blocked nothing the page needed",
        blocked == s(""),
        &blocked,
    );
    drop(b);

    // 8: the map's reads panel.
    map_reads_checks(&mut r, &chrome, work, &proxy, site.port, bare.port);
    drop(proxy);
    drop(site);
    drop(bare);
    r.done();
}

/// The real map and viewer, on invented map rows and the synthetic reads: the reads pill, the
/// same-origin iframe that follows the camera, closed findings, Escape, R, the Reads gauge, the
/// card's reads, a site without the viewer, the planted page policy without frame-src, the demo.
fn map_reads_checks(
    r: &mut Report,
    chrome: &Path,
    work: &Path,
    proxy: &Proxy,
    viewer_port: u16,
    bare_port: u16,
) {
    let mut b = Browser::launch(chrome, &work.join("chrome-map-1g")).unwrap();
    let call = |b: &mut Browser, m: &str, p: &str| b.call(m, loads(p).unwrap()).unwrap();
    call(
        &mut b,
        "Emulation.setEmulatedMedia",
        r#"{"features":[{"name":"prefers-color-scheme","value":"dark"},{"name":"prefers-reduced-motion","value":"reduce"}]}"#,
    );
    call(
        &mut b,
        "Emulation.setDeviceMetricsOverride",
        r#"{"width":1440,"height":900,"deviceScaleFactor":1,"mobile":false}"#,
    );
    let w = |b: &mut Browser, js: &str, secs: u64| b.wait(js, Duration::from_secs(secs)).is_some();
    let e = |b: &mut Browser, js: &str| b.eval(js).unwrap_or(Value::Null);
    let mut visit = 0;
    let mut mount = |b: &mut Browser, demo: bool| -> bool {
        visit += 1;
        b.go(&format!(
            "http://127.0.0.1:{}/map-1g/?visit={visit}#chr1:1-49000",
            proxy.port
        ))
        .unwrap();
        if !w(b, "!!(window.PGPMap && window.PGPUI)", 20) {
            return false;
        }
        e(
            b,
            &format!(
                "window.closedResult = true;document.documentElement.dataset.run='S'; \
             document.documentElement.dataset.mode = {}; \
             window.testMap = PGPMap.create({{url:'/map-1g.json',\
             veil:function(id){{return id==='closed' && window.closedResult ? 'closed':'open';}},\
             known:function(){{return true;}},now:function(){{return true;}}}});\
             document.body.appendChild(testMap.node);testMap.refresh('one'); true",
                js_str(if demo { "replay" } else { "live" })
            ),
        );
        w(b, "!!document.querySelector('.gmap[data-level=reads]')", 20)
    };
    let loaded = |b: &mut Browser| {
        w(
            b,
            "(function(){var f=document.querySelector('.gm-reads iframe');\
        return f && f.contentDocument && f.contentDocument.documentElement.dataset.pgpStatus==='ready';})()",
            30,
        )
    };
    let locus = |b: &mut Browser| {
        st(&e(
            b,
            "savedFrame.contentDocument.getElementById('igv').shadowRoot.querySelector('input.igv-search-input').value",
        ))
    };
    let resized = |b: &mut Browser| {
        call(
            b,
            "Runtime.evaluate",
            r#"{"expression":"new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))","awaitPromise":true}"#,
        );
    };

    let m = mount(&mut b, false);
    r.ok("map: opens at a 49 kb view", m, "");
    let pill = truthy(&e(
        &mut b,
        "!!document.querySelector('.gm-reads-pill:not([hidden])')",
    ));
    r.ok("map: under 50 kb offers the reads pill", pill, "");
    assert!(
        truthy(&e(&mut b, "!!document.querySelector('.gm-reads-pill')")),
        "no reads pill to go on with"
    );
    e(
        &mut b,
        "document.querySelector('.gm-reads-pill').click(); true",
    );
    let ok = loaded(&mut b);
    r.ok(
        "map: the reads panel loads the viewer in an iframe at the nearest open finding",
        ok,
        e(&mut b, "document.querySelector('.gm-reads').textContent"),
    );
    let hash = e(
        &mut b,
        "document.querySelector('.gm-reads iframe').contentWindow.location.hash",
    );
    r.ok(
        "map: the iframe locus is 1-based, finding plus or minus 50",
        hash == s("#chr1:401-501"),
        &hash,
    );
    let closed = e(
        &mut b,
        "!document.querySelector('.gm-marks [data-id=closed]') && !document.querySelector('.gm-card-reads[data-id=closed]')",
    );
    r.ok(
        "map: a closed result has no card or reads link",
        truthy(&closed),
        "",
    );
    e(
        &mut b,
        "window.savedFrame=document.querySelector('.gm-reads iframe');window.savedViewerDocument=savedFrame.contentDocument;\
              window.savedHistoryLength=history.length;location.hash='chr1:1-49001'; true",
    );
    let ok = w(
        &mut b,
        "savedFrame===document.querySelector('.gm-reads iframe') && savedFrame.contentWindow.location.hash==='#chr1:1-49001'",
        20,
    );
    r.ok(
        "map: moving updates the same iframe to the map's locus",
        ok,
        "",
    );
    let ok = truthy(&e(
        &mut b,
        "savedViewerDocument===savedFrame.contentDocument && history.length===savedHistoryLength+1",
    ));
    r.ok(
        "map: a camera move keeps the viewer document and adds no iframe history entry",
        ok,
        "",
    );
    let ok = w(
        &mut b,
        "(function(){var h=savedFrame.contentDocument.getElementById('igv');\
        var v=h && h.shadowRoot && h.shadowRoot.querySelector('.igv-viewport[data-track-type=alignment]');\
        var n=v && v.querySelector('.igv-zoom-in-notice-container');\
        return !!(v && v.querySelector('canvas') && n && n.style.display!=='flex');})()",
        20,
    );
    r.ok("viewer: alignments remain visible at 49 kb", ok, "");
    let expected = locus(&mut b);
    e(&mut b, "savedFrame.style.width='50%'; true");
    resized(&mut b);
    sleep(500); // past viewer.js's 100 ms resize timer
    let same = format!(
        "savedFrame.contentDocument.getElementById('igv').shadowRoot.querySelector('input.igv-search-input').value==={}",
        js_str(&expected)
    );
    let ok = w(&mut b, &same, 20);
    r.ok(
        "viewer: resizing keeps the full requested locus",
        ok,
        locus(&mut b),
    );
    e(&mut b, "savedFrame.style.width=''; true");
    resized(&mut b);
    sleep(500);
    let ok = w(&mut b, &same, 20);
    r.ok(
        "viewer: restoring the width keeps the full requested locus",
        ok,
        locus(&mut b),
    );
    // A locus the reader types in igv stays through a resize: only a new map locus moves igv again.
    e(
        &mut b,
        "var i=savedFrame.contentDocument.getElementById('igv').shadowRoot.querySelector('input.igv-search-input'); i.focus();i.select(); true",
    );
    call(&mut b, "Input.insertText", r#"{"text":"chr1:2001-3000"}"#);
    for kind in ["keyDown", "keyUp"] {
        call(
            &mut b,
            "Input.dispatchKeyEvent",
            &format!(
                r#"{{"type":"{kind}","key":"Enter","code":"Enter","windowsVirtualKeyCode":13,"text":"\r"}}"#
            ),
        );
    }
    let typed = w(
        &mut b,
        "/^chr1:2,0/.test(savedFrame.contentDocument.getElementById('igv').shadowRoot.querySelector('input.igv-search-input').value)",
        20,
    );
    e(&mut b, "savedFrame.style.width='50%'; true");
    resized(&mut b);
    sleep(500); // past viewer.js's 100 ms resize timer
    let now = locus(&mut b);
    let span = now
        .strip_prefix("chr1:")
        .and_then(|x| x.split_once('-'))
        .map(|(a, z)| {
            (
                a.replace(',', "").parse::<u64>().unwrap_or(0),
                z.replace(',', "").parse::<u64>().unwrap_or(u64::MAX),
            )
        });
    r.ok(
        "viewer: a locus typed in igv stays through a resize",
        typed && span.is_some_and(|(a, z)| (1500..=3000).contains(&a) && z <= 3500),
        &now,
    );
    e(&mut b, "savedFrame.style.width=''; true");
    resized(&mut b);
    e(
        &mut b,
        "savedFrame.contentDocument.dispatchEvent(new KeyboardEvent('keydown',{key:'Escape',bubbles:true})); true",
    );
    let ok = truthy(&e(
        &mut b,
        "document.querySelector('.gm-reads').hidden && document.activeElement===document.querySelector('.gm-canvas') && \
        document.querySelector('.gmap').dataset.view==='chr1:1-49001'",
    ));
    r.ok(
        "map: Escape closes the reads and returns focus without moving the map",
        ok,
        "",
    );
    e(
        &mut b,
        "document.querySelector('.gm-canvas').dispatchEvent(new KeyboardEvent('keydown',{key:'r',bubbles:true})); true",
    );
    let ok = loaded(&mut b);
    r.ok("map: R reopens the panel", ok, "");
    e(
        &mut b,
        "document.querySelector('.gm-reads-close').click();window.gaugeHistoryLength=history.length;\
              document.querySelector('[data-stop=reads]').click(); true",
    );
    let ok = w(
        &mut b,
        "!document.querySelector('.gm-reads').hidden && document.querySelector('.gmap').dataset.level==='reads' && \
        (function(){var m=/:(\\d+)-(\\d+)$/.exec(document.querySelector('.gmap').dataset.view||'');return !!m&&m[2]-m[1]+1===2000;})()",
        20,
    );
    let view = e(&mut b, "document.querySelector('.gmap').dataset.view");
    r.ok("map: Reads gauge opens a 2 kb view and the panel", ok, view);
    let ok = truthy(&e(&mut b, "history.length===gaugeHistoryLength+1"));
    r.ok("map: the Reads gauge adds one camera history entry", ok, "");
    e(
        &mut b,
        "document.querySelector('.gm-reads-close').click();location.hash='chr1:1-2000'; true",
    );
    w(
        &mut b,
        "document.querySelector('.gmap').dataset.view==='chr1:1-2000'",
        20,
    );
    e(
        &mut b,
        "document.querySelector('.gm-marks [data-id=open]').focus();document.querySelector('.gm-card-reads').click(); true",
    );
    let ok = w(
        &mut b,
        "document.querySelector('.gm-reads iframe').contentWindow.location.hash==='#chr1:401-501'",
        20,
    );
    r.ok("map: the card opens reads at its finding", ok, "");
    e(
        &mut b,
        "window.closedResult=false;testMap.refresh('one');document.querySelector('.gm-marks [data-id=closed]').focus();\
              document.querySelector('.gm-card-reads').click(); true",
    );
    let ok = w(
        &mut b,
        "document.querySelector('.gm-reads iframe').contentWindow.location.hash==='#chr1:501-601'",
        20,
    );
    r.ok("map: an opened finding offers reads", ok, "");
    e(
        &mut b,
        "window.closedResult=true;testMap.refresh('one'); true",
    );
    let ok = truthy(&e(
        &mut b,
        "document.querySelector('.gm-reads').hidden && !document.querySelector('.gm-card-reads') && \
        !document.querySelector('.gm-reads a[href]') && !document.querySelector('.gm-reads').textContent.includes('Closed finding')",
    ));
    r.ok(
        "map: closing a result removes its open reads panel and link",
        ok,
        "",
    );

    // No finding in view: the view itself; a site with no viewer flags: the server's reason.
    e(&mut b, "location.hash='chr1:10000-59000'; true");
    w(
        &mut b,
        "document.querySelector('.gmap').dataset.view==='chr1:10000-59000'",
        20,
    );
    proxy.set(|s| s.upstream = bare_port);
    e(
        &mut b,
        "document.querySelector('.gm-reads-pill').click(); true",
    );
    let ok = w(
        &mut b,
        "document.querySelector('.gm-reads').textContent.includes('not set up')",
        20,
    ) && truthy(&e(&mut b, "!document.querySelector('.gm-reads iframe')"));
    r.ok(
        "map: a missing viewer shows the server's reason, without a broken frame",
        ok,
        e(&mut b, "document.querySelector('.gm-reads').textContent"),
    );
    proxy.set(|s| s.upstream = viewer_port);

    // The control: the page policy without frame-src blocks the reads iframe, and the check sees it.
    proxy.set(|s| s.block_frames = true);
    mount(&mut b, false);
    e(
        &mut b,
        "window.blockedFrames=0;document.addEventListener('securitypolicyviolation',function(e){\
              if(e.violatedDirective==='frame-src')window.blockedFrames++;});document.querySelector('.gm-reads-pill').click(); true",
    );
    let blocked = w(&mut b, "window.blockedFrames>0", 20)
        && truthy(&e(
            &mut b,
            "(function(){var f=document.querySelector('.gm-reads iframe');return !f || !f.contentDocument || \
            f.contentDocument.documentElement.dataset.pgpStatus!=='ready';})()",
        ));
    r.control(
        "a page policy without frame-src: the reads iframe must not load",
        if blocked {
            vec!["the policy blocked the reads iframe".into()]
        } else {
            vec![]
        },
    );
    proxy.set(|s| s.block_frames = false);
    mount(&mut b, true);
    e(
        &mut b,
        "document.querySelector('.gm-reads-pill').click(); true",
    );
    let ok = w(
        &mut b,
        "document.querySelector('.gm-reads').textContent.includes('The genome viewer is not in the demo')",
        20,
    ) && truthy(&e(&mut b, "!document.querySelector('.gm-reads iframe')"));
    r.ok("map: the static demo explains why it has no viewer", ok, "");
}
