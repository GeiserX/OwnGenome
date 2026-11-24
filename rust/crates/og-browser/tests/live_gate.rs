//! Chrome coverage of report gate isolation, transparency, and failure paths.
//!
//! `demo_report_gate` runs on the static demo: the opaque closed handles, the
//! order of the closed controls, each guard of the gate with its planted
//! mistake, the transparency of the gate on adult and no-act pages, the slow
//! release, a closed child's family data that leaves the family page as it
//! was and two children with the same act IDs. `live_serve_report_gate`
//! runs on `og live serve`: the finding route of a child and of an adult,
//! the live-mode requests, a parent's page whatever card holds a child's
//! results, a closed child's side by side, the family page's closed and
//! opened lines, the cross-family partner, an unreadable family file and a
//! failed family fetch. `families_live_page` runs `og live serve` on
//! tests/fixtures/families_live: the family card, the pedigree, the couple
//! screen's risk drawing, the trio and compound-het lines and the side by
//! side. Moved here from tests/test_live_demo.py and
//! tests/test_owngenome_live_family.py; every check whose Python had a
//! planted mistake runs again on it here and must fail ("[RED] control").

use og_browser::{Browser, Scratch, StaticServer, Value, chrome_or_skip, js_str, repo_root};
use pycompat::json::{dumps, loads};
use std::fs;
use std::io::{BufRead, BufReader, Read, Write};
use std::net::TcpListener;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::thread;
use std::time::Duration;

const ADULT: &str = "DEMO-001";
// On og live serve DEMO-001 is the demo's replayed run, with no step done
// yet; DEMO-002, the other parent, has its report.
const LIVE_ADULT: &str = "DEMO-002";
const CHILD: &str = "DEMO-003";
const OTHER: &str = "DEMO-004";
const FAMILY: &str = "DEMO-FAMILY";
const PARTNERS: &str = "DEMO-A-PARTNERS";
const PARTNERS_LAST: &str = "DEMO-Z-PARTNERS";
const CHILD_URL: &str = "/runs/DEMO-003/";
const READY: &str = "document.documentElement.dataset.pgpStatus === 'ready'";
const CLOSED: &str = "document.documentElement.dataset.pgpGated === 'true'";
const OPEN: &str = "!document.documentElement.dataset.pgpGated && !!document.querySelector('#results .finding[data-level=act]')";
const WAIT: Duration = Duration::from_secs(16);

struct Fixture {
    browser: Browser,
    scratch: Scratch,
    og: PathBuf,
    site: PathBuf,
    work: PathBuf,
}

impl Fixture {
    fn build(label: &str) -> Option<Self> {
        let chrome = chrome_or_skip(label)?;
        let og = og_binary()?;
        let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), label);
        let site = scratch.path().join("site");
        let work = scratch.path().join("work");
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
        let mut browser = Browser::launch(&chrome, &scratch.path().join("chrome")).unwrap();
        call(&mut browser, "Network.enable", "{}");
        call(
            &mut browser,
            "Network.setCacheDisabled",
            r#"{"cacheDisabled":true}"#,
        );
        Some(Self {
            browser,
            scratch,
            og,
            site,
            work,
        })
    }
}

/// The og binary built beside this test, or `$OG`; None (a skip) when it is
/// missing outside CI.
fn og_binary() -> Option<PathBuf> {
    let og = std::env::var_os("OG").map(PathBuf::from).or_else(|| {
        let exe = std::env::current_exe().ok()?;
        Some(exe.parent()?.parent()?.join("og"))
    })?;
    if !og.is_file() {
        assert_ne!(std::env::var("CI").as_deref(), Ok("true"), "build og in CI");
        eprintln!("[SKIP] build og or set OG");
        return None;
    }
    Some(og)
}

fn call(b: &mut Browser, method: &str, params: &str) -> Value {
    b.call(method, loads(params).unwrap()).unwrap()
}

fn eval(b: &mut Browser, js: &str) -> Value {
    b.eval(js).unwrap_or_else(|e| panic!("{js}: {e}"))
}

fn wait(b: &mut Browser, js: &str) -> Value {
    b.wait(js, WAIT).unwrap_or_else(|| {
        let page = b.eval("[location.href, Array.from(document.documentElement.attributes, a => a.name + '=' + a.value.slice(0, 60)), (window.ogRequests || []).map(r => r.join(' ')), document.body.textContent.slice(0,350)]");
        panic!("{js} did not settle: {page:?}; {:?}", b.page_errors());
    })
}

/// Open `url` from about:blank. A release's own location.replace can land
/// after about:blank, so about:blank is asked for again until it stays.
fn go(b: &mut Browser, url: &str) {
    for _ in 0..3 {
        b.go("about:blank").unwrap();
        if b.wait("location.href === 'about:blank'", Duration::from_secs(5))
            .is_some()
        {
            thread::sleep(Duration::from_millis(100));
            if b.eval("location.href").ok() == Some(Value::Str("about:blank".into())) {
                b.go(url).unwrap();
                return;
            }
        }
    }
    panic!("about:blank never stayed before {url}");
}

fn run_url(base: &str, run: &str) -> String {
    format!("{base}/runs/{run}/")
}

fn install(b: &mut Browser, source: &str) -> String {
    let installed = call(
        b,
        "Page.addScriptToEvaluateOnNewDocument",
        &format!(r#"{{"source":{}}}"#, js_str(source)),
    );
    text(field(&installed, "identifier")).to_owned()
}

fn uninstall(b: &mut Browser, id: &str) {
    call(
        b,
        "Page.removeScriptToEvaluateOnNewDocument",
        &format!(r#"{{"identifier":{}}}"#, js_str(id)),
    );
}

fn doc(path: &Path) -> Value {
    loads(&fs::read_to_string(path).unwrap()).unwrap()
}

fn save(path: &Path, value: &Value) {
    fs::write(path, dumps(value, &Default::default()).unwrap()).unwrap();
}

fn obj(value: &mut Value) -> &mut indexmap::IndexMap<String, Value> {
    let Value::Object(value) = value else {
        panic!("expected object")
    };
    value
}

fn field<'a>(value: &'a Value, key: &str) -> &'a Value {
    let Value::Object(value) = value else {
        panic!("expected object with {key}")
    };
    &value[key]
}

fn text(value: &Value) -> &str {
    let Value::Str(text) = value else {
        panic!("expected text, got {value:?}")
    };
    text
}

fn array(value: &mut Value) -> &mut Vec<Value> {
    let Value::Array(value) = value else {
        panic!("expected array")
    };
    value
}

fn items(value: &Value) -> &[Value] {
    let Value::Array(value) = value else {
        panic!("expected array, got {value:?}")
    };
    value
}

fn int(value: &Value) -> i64 {
    match value {
        Value::Int(n) => i64::try_from(*n).unwrap(),
        other => panic!("expected a number, got {other:?}"),
    }
}

fn s(value: &str) -> Value {
    Value::Str(value.into())
}

fn is_act(finding: &Value) -> bool {
    text(field(finding, "level")) == "act"
}

fn replace_once(source: &str, from: &str, to: &str) -> String {
    assert_eq!(source.matches(from).count(), 1, "control anchor: {from}");
    source.replacen(from, to, 1)
}

fn pass(what: &str) {
    println!("[PASS] {what}");
}

fn red(what: &str, seen: &dyn std::fmt::Debug) {
    let seen = format!("{seen:?}");
    let cut: String = seen.chars().take(300).collect();
    println!("[RED] control: {what} -- the check fails: {cut}");
}

/// `og live results -` with `summary` on stdin and `args` after it.
fn og_results(og: &Path, summary: &str, args: &[&str]) -> Value {
    let mut child = Command::new(og)
        .args(["live", "results", "-"])
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(summary.as_bytes())
        .unwrap();
    let out = child.wait_with_output().unwrap();
    assert!(
        out.status.success(),
        "og live results: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    loads(&String::from_utf8(out.stdout).unwrap()).unwrap()
}

// Draws the same numbers in every load (DRAW for each), so which act a closed
// control opens is fixed for the load.
const SEEDED: &str = "crypto.getRandomValues = function (a) { for (var i = 0; i < a.length; i++) { a[i] = DRAW; } return a; };";
const SEARCH_CLOSED: &str =
    "Results this page keeps closed are left out of the search until they are opened.";

// Each network read a page script starts, in the order it starts it, below
// every wrapper of the page's own: the path, and for fetch its cache mode and
// credentials. Installed before any script of the page runs.
const REQUEST_LOG: &str = r#"(function () {
  var log = window.ogRequests = [], f = window.fetch, E = window.EventSource;
  window.fetch = function (u, o) {
    log.push(['fetch', new URL(u && u.url || u, location.href).pathname, (o && o.cache) || 'default', (o && o.credentials) || 'same-origin']);
    return f.apply(this, arguments);
  };
  if (E) {
    var W = function (u, o) { log.push(['events', new URL(u, location.href).pathname]); return new E(u, o); };
    W.prototype = E.prototype; W.CONNECTING = E.CONNECTING; W.OPEN = E.OPEN; W.CLOSED = E.CLOSED;
    window.EventSource = W;
  }
}());"#;

/// core.js with the report gate cut out: the transparent reference a page
/// without the gate is compared to.
fn bare_core(source: &str) -> String {
    let marker = "  // --- the report gate";
    assert_eq!(source.matches(marker).count(), 1, "gate marker moved");
    let (prefix, _) = source.split_once(marker).unwrap();
    format!(
        "{prefix}
  PGPUI.ready = Promise.resolve();
  PGPUI.gateOwnsPage = function () {{ return false; }};
  PGPUI.openedActs = function () {{ return {{}}; }};
  PGPUI.closeAct = function () {{ return false; }};
}}());
"
    )
}

fn resources(b: &mut Browser, name: &str) -> i64 {
    int(&eval(
        b,
        &format!(
            "performance.getEntriesByType('resource').filter(function (r) {{ return r.name.endsWith({}); }}).length",
            js_str(name)
        ),
    ))
}

/// Type `term` in the header search once its index has loaded; the text of
/// what it shows.
fn type_search(b: &mut Browser, term: &str) -> String {
    wait(
        b,
        "document.querySelector('input[type=search]') && !document.querySelector('input[type=search]').disabled",
    );
    eval(
        b,
        &format!(
            "var s=document.querySelector('input[type=search]'); s.value={}; s.dispatchEvent(new Event('input')); true",
            js_str(term)
        ),
    );
    wait(
        b,
        "document.documentElement.dataset.pgpSearch !== undefined",
    );
    text(&eval(
        b,
        "document.querySelector('.search-results').textContent",
    ))
    .to_owned()
}

/// Open every closed act of the page. The last opening loads the page again.
fn open_gate(b: &mut Browser) {
    wait(b, READY);
    let closed = "document.querySelectorAll('[data-opening^=closed]').length";
    for _ in 0..12 {
        let n = match b.eval(closed) {
            Ok(Value::Int(n)) if n > 0 => n,
            _ => return,
        };
        let _ = b.eval("document.querySelector('[data-opening^=closed]').click(); true");
        if n == 1 {
            return;
        }
        wait(b, &format!("{closed} < {n}"));
    }
    panic!("the report never opened");
}

/// Open the act whose own record has `title`. A closed control opens an act
/// the reader cannot choose, so each is tried in turn and a wrong one is
/// closed again. The number of clicks.
fn open_act(b: &mut Browser, title: &str, scope: &str) -> i64 {
    let closed = js_str(&format!("{scope}[data-opening^=closed]"));
    let opened = js_str(&format!("{scope}[data-opening^=opened]"));
    let mut clicks = 0;
    for k in 0..10 {
        let n = int(&eval(
            b,
            &format!("document.querySelectorAll({closed}).length"),
        ));
        assert!(n > k, "no closed control opens {title}");
        let m = int(&eval(
            b,
            &format!("document.querySelectorAll({opened}).length"),
        )) + 1;
        eval(
            b,
            &format!("document.querySelectorAll({closed})[{k}].click(); true"),
        );
        clicks += 1;
        if n == 1 {
            return clicks;
        }
        let got = js_str(&format!("{scope}[data-finding=opened-{m}] .opened-act h2"));
        wait(b, &format!("!!document.querySelector({got})"));
        if text(&eval(
            b,
            &format!("document.querySelector({got}).textContent"),
        )) == title
        {
            return clicks;
        }
        let again = js_str(&format!("{scope}[data-opening=opened-{m}]"));
        eval(b, &format!("document.querySelector({again}).click(); true"));
        clicks += 1;
        wait(
            b,
            &format!("document.querySelectorAll({closed}).length === {n}"),
        );
    }
    panic!("no closed control opens {title}");
}

/// The child's results.json rewritten with `acts` as its only acts (before
/// the other findings when `first`, else after them), listed in its veil.
fn write_acts(path: &Path, data: &Value, acts: &[Value], first: bool) {
    let mut d = data.clone();
    let o = obj(&mut d);
    let findings = array(o.get_mut("findings").unwrap());
    findings.retain(|f| !is_act(f));
    let at = if first { 0 } else { findings.len() };
    findings.splice(at..at, acts.iter().cloned());
    if !matches!(o.get("veil"), Some(Value::Object(_))) {
        o.insert("veil".into(), Value::Object(Default::default()));
    }
    let ids = acts.iter().map(|f| field(f, "id").clone()).collect();
    obj(o.get_mut("veil").unwrap()).insert("acts".into(), Value::Array(ids));
    save(path, &d);
}

fn acts_of(data: &Value) -> Vec<Value> {
    items(field(data, "findings"))
        .iter()
        .filter(|f| is_act(f))
        .cloned()
        .collect()
}

/// The child's results with three acts of their own titles: the demo's
/// first two and an unplaced two-gene one. The saved file and the acts.
fn three_acts(site: &Path) -> (PathBuf, Vec<u8>, Vec<Value>) {
    let path = site.join(format!("runs/{CHILD}/results.json"));
    let saved = fs::read(&path).unwrap();
    let data = doc(&path);
    let mut acts = acts_of(&data);
    let mut third = acts[0].clone();
    let t = obj(&mut third);
    t.insert("id".into(), s("gate-third"));
    t.insert("gene".into(), s("HBA1"));
    t.insert("genes".into(), Value::Array(vec![s("HBA1"), s("HBA2")]));
    t.insert("viewer".into(), Value::Null);
    t.insert("title".into(), s("Gate fixture result in HBA1"));
    t.insert("detail".into(), s("HBA1 and HBA2"));
    acts.truncate(2);
    acts.push(third);
    let titles: std::collections::BTreeSet<_> =
        acts.iter().map(|f| text(field(f, "title"))).collect();
    assert_eq!(titles.len(), 3, "the acts need their own titles");
    write_acts(&path, &data, &acts, false);
    (path, saved, acts)
}

/// test_live_demo.py's gate handles and the end of its gate addresses: a
/// closed page's address names an opaque handle whatever was typed, the
/// handle follows its act when another act opens, and the opened report
/// keeps #f/<id> as the full report does.
fn gate_handles(f: &mut Fixture, base: &str, core: &Path, source: &str) {
    let (rj, saved, acts) = three_acts(&f.site);
    let b = &mut f.browser;
    let child = run_url(base, CHILD);
    let id0 = text(field(&acts[0], "id")).to_owned();
    let guesses = |b: &mut Browser| -> Vec<Value> {
        let typed = [
            format!("#f/{id0}"),
            "#f/closed-99".into(),
            "#f/opened-1".into(),
            "#f/closed-1".into(),
            "#s/exomiser".into(),
            "?chr=17&start=43050000#f/closed-99".into(),
        ];
        typed
            .iter()
            .map(|h| {
                go(b, &format!("{child}{h}"));
                wait(
                    b,
                    &format!("document.readyState === 'complete' && ({READY}) && ({CLOSED})"),
                );
                eval(b, "[location.hash, location.search]")
            })
            .collect()
    };
    let first = Value::Array(vec![s("#f/closed-1"), s("")]);
    let got = guesses(b);
    assert!(
        got.iter().all(|g| *g == first),
        "gate handles: typed and guessed ids must share the first closed handle: {got:?}"
    );
    pass(
        "gate handles: typed and guessed ids, a query and a section address share the first opaque closed handle",
    );
    fs::write(
        core,
        replace_once(
            source,
            "var to = found ? location.hash : all.length ? '#f/' + all[0].handle : '';",
            "var to = location.hash;",
        ),
    )
    .unwrap();
    let got = guesses(b);
    assert!(
        !got.iter().all(|g| *g == first),
        "the kept-hash control went unnoticed"
    );
    red("a closed page that keeps the typed address", &got);
    fs::write(core, source).unwrap();

    let follow = |b: &mut Browser| -> Value {
        go(b, &child);
        wait(b, CLOSED);
        for h in ["closed-2", "closed-3"] {
            eval(b, &format!("location.hash = 'f/{h}'; true"));
            wait(b, &format!("location.hash === '#f/{h}'"));
        }
        eval(b, "history.back(); true");
        wait(b, "location.hash === '#f/closed-2'");
        eval(b, "history.forward(); true");
        wait(b, "location.hash === '#f/closed-3'");
        eval(
            b,
            "document.querySelector('[data-opening=closed-1]').click(); true",
        );
        wait(b, "!!document.querySelector('[data-opening=opened-1]')");
        eval(b, "location.hash")
    };
    let got = follow(b);
    assert_eq!(text(&got), "#f/closed-2", "gate handles follow their act");
    pass("gate handles follow their act after another act opens");
    fs::write(core, replace_once(source, "if (at) {", "if (false) {")).unwrap();
    let got = follow(b);
    assert_ne!(
        text(&got),
        "#f/closed-2",
        "the handle control went unnoticed"
    );
    red(
        "an address that keeps its handle when another act opens",
        &got,
    );
    fs::write(core, source).unwrap();

    go(b, &child);
    open_gate(b);
    wait(b, "!!document.querySelector('#results .finding')");
    let links = eval(
        b,
        &format!(
            r##"[!!document.querySelector('a[href="#f/{id0}"]'), !!document.querySelector('a[href^="#f/opened-"],a[href^="#f/closed-"]')]"##
        ),
    );
    eval(b, &format!("location.hash = 'f/{id0}'; true"));
    let kept = b.wait(
        &format!("location.hash === '#f/{id0}' && !!document.getElementById('f-{id0}')"),
        Duration::from_secs(5),
    );
    assert!(
        links == Value::Array(vec![Value::Bool(true), Value::Bool(false)]) && kept.is_some(),
        "the opened report must link and keep #f/<id>: {links:?} {kept:?}"
    );
    pass("gate addresses: the opened report links and keeps #f/<id>, as the full report does");
    fs::write(rj, saved).unwrap();
}

/// test_live_demo.py's gate order: only an unopened act changes (its rule
/// and gene, so its place in the producer's order); the opened act is the
/// same. The pages must not differ, and which act a closed control opens
/// comes from a draw of this page load.
fn gate_order(f: &mut Fixture, base: &str, core: &Path, source: &str) {
    let rj = f.site.join(format!("runs/{CHILD}/results.json"));
    let raw = fs::read(&rj).unwrap();
    let data = doc(&rj);
    let acts = acts_of(&data);
    let (opened, closed) = (acts[0].clone(), acts[1].clone());
    let mut tp53 = closed.clone();
    let t = obj(&mut tp53);
    t.insert("id".into(), s("cpsr-tp53-c-743g-a"));
    t.insert("gene".into(), s("TP53"));
    t.insert("genes".into(), Value::Array(vec![s("TP53")]));
    t.insert("title".into(), s("Cancer predisposition (TP53)"));
    t.insert("detail".into(), s("TP53 c.743G>A"));
    t.insert("viewer".into(), s("chr17:7674220"));
    let title = text(field(&opened, "title")).to_owned();
    let b = &mut f.browser;
    let routes = [
        (run_url(base, CHILD), String::new()),
        (run_url(base, FAMILY), format!("[data-gate-child={CHILD}] ")),
    ];
    let mut pages = Vec::new();
    for order in [vec![opened.clone(), closed], vec![tp53, opened]] {
        write_acts(&rj, &data, &order, true);
        for (url, scope) in &routes {
            go(b, url);
            wait(
                b,
                &format!(
                    "!!document.querySelector({})",
                    js_str(&format!("{scope}[data-opening^=closed]"))
                ),
            );
            let completed = resources(b, &format!("/{CHILD}/results.json"));
            let clicks = open_act(b, &title, scope);
            wait(
                b,
                &format!(
                    "!!document.querySelector({})",
                    js_str(&format!("{scope}.opened-act"))
                ),
            );
            if !scope.is_empty() {
                wait(
                    b,
                    &format!(
                        "performance.getEntriesByType('resource').filter(function(r){{return r.name.endsWith('/{CHILD}/results.json');}}).length >= {}",
                        completed + clicks
                    ),
                );
            }
            thread::sleep(Duration::from_millis(300));
            pages.push(eval(
                b,
                "[document.documentElement.outerHTML, document.body.innerText, location.href, document.title]",
            ));
        }
    }
    for (i, (url, _)) in routes.iter().enumerate() {
        assert!(
            pages[i] == pages[i + 2],
            "gate order: an unopened act earlier in the report changed {url}"
        );
    }
    pass(
        "gate order: an unopened act earlier or later in the report leaves the page with an opened act unchanged",
    );
    let first_opens = |b: &mut Browser, draw: &str| -> String {
        let seed = install(b, &SEEDED.replace("DRAW", draw));
        go(b, &run_url(base, CHILD));
        wait(b, "!!document.querySelector('[data-opening=closed-1]')");
        eval(
            b,
            "document.querySelector('[data-opening=closed-1]').click(); true",
        );
        let got = wait(
            b,
            "document.querySelector('[data-finding=opened-1] .opened-act h2')?.textContent",
        );
        uninstall(b, &seed);
        text(&got).to_owned()
    };
    let draws = |b: &mut Browser| {
        std::collections::BTreeSet::from([first_opens(b, "0"), first_opens(b, "4294967295")])
    };
    let got = draws(b);
    assert_eq!(got.len(), 2, "gate order: the draw chose nothing: {got:?}");
    pass(
        "gate order: which act the first closed control opens comes from the page load, not the report",
    );
    fs::write(
        core,
        replace_once(source, ": shuffled(ids);", ": ids.slice();"),
    )
    .unwrap();
    let got = draws(b);
    assert_eq!(got.len(), 1, "the producer-order control went unnoticed");
    red("closed controls in the producer order", &got);
    fs::write(core, source).unwrap();
    fs::write(rj, raw).unwrap();
}

/// What `<html>` may keep once a child's open report is cleared: the gate's
/// own two attributes.
const LEFT: &str = "Array.prototype.map.call(document.documentElement.attributes, function (a) { return a.name; }).filter(function (n) { return (n.indexOf('data-pgp-') === 0 && n !== 'data-pgp-status' && n !== 'data-pgp-gated') || n === 'data-finding'; })";

/// The edits that plant one mistake in core.js: (anchor, replacement).
type Edits<'a> = Vec<(&'a str, &'a str)>;

#[derive(Clone, Copy, Debug)]
enum Guard {
    Stale,
    RefreshCleanup,
    Address,
    Partial,
    RawView,
    Handoff,
    FamilyRows,
    Cleanup,
    Ticket,
    CloseAttrs,
    Pagehide,
    Back,
}

struct Guards<'a> {
    base: &'a str,
    core: &'a Path,
    rj: &'a Path,
    data: &'a Value,
    fid: String,
}

impl Guards<'_> {
    fn visit_child(&self, b: &mut Browser) {
        go(b, &run_url(self.base, CHILD));
        wait(b, CLOSED);
    }

    fn opened(&self, b: &mut Browser) {
        self.visit_child(b);
        open_gate(b);
        wait(
            b,
            "!!document.querySelector('#results .finding') && !!document.documentElement.dataset.pgpResults",
        );
    }

    /// The child's results with every act's text changed by `mark`.
    fn changed(&self, mark: &str, first_only: bool) {
        let mut d = self.data.clone();
        let mut done = false;
        for f in array(obj(&mut d).get_mut("findings").unwrap()) {
            if is_act(f) && !(first_only && done) {
                let o = obj(f);
                let t = format!("{}{mark}", text(&o["text"]));
                o.insert("text".into(), s(&t));
                done = true;
            }
        }
        save(self.rj, &d);
    }

    /// Whether the guard holds, and what was seen.
    fn probe(&self, b: &mut Browser, guard: Guard) -> (bool, Value) {
        let fid = js_str(&self.fid);
        match guard {
            Guard::Partial => {
                self.visit_child(b);
                eval(
                    b,
                    "document.querySelector('[data-opening^=closed]').click(); true",
                );
                let ok = b
                    .wait(
                        &format!("({CLOSED}) && document.querySelectorAll('[data-opening^=closed]').length === 1"),
                        Duration::from_secs(2),
                    )
                    .is_some();
                let released = !ok
                    && b.wait(
                        "!!document.documentElement.dataset.pgpResults && !document.documentElement.dataset.pgpGated",
                        Duration::from_secs(10),
                    )
                    .is_some();
                (ok, Value::Bool(released))
            }
            Guard::RawView => {
                self.visit_child(b);
                b.wait("document.readyState === 'complete'", Duration::from_secs(2));
                if fs::read_to_string(self.core)
                    .unwrap()
                    .contains("data-gate-raw-read")
                {
                    assert!(
                        b.wait(
                            "document.documentElement.dataset.gateRawRead === 'done'",
                            Duration::from_secs(2)
                        )
                        .is_some(),
                        "raw probe did not settle"
                    );
                }
                let leaked = eval(
                    b,
                    "document.documentElement.outerHTML.includes('FUTURE_GATE_SENTINEL')",
                );
                (leaked == Value::Bool(false), leaked)
            }
            Guard::Handoff => {
                // A view that asks a closed run for its files gets no answer
                // at all, whatever it reads the body with. An answer is
                // stored as [its length], so an empty one counts too.
                self.visit_child(b);
                eval(
                    b,
                    "window.gateRaw = null; fetch('results.json').then(function (r) { return r.text(); }).then(function (t) { window.gateRaw = [t.length]; }); true",
                );
                let got = b.wait("window.gateRaw", Duration::from_secs(2));
                (got.is_none(), got.unwrap_or(Value::Null))
            }
            Guard::FamilyRows => {
                go(b, &run_url(self.base, FAMILY));
                wait(
                    b,
                    &format!("!!document.querySelector('[data-gate-child={CHILD}]')"),
                );
                // What the family view gets of the closed child: its trio,
                // relationship and shared-segment rows, and its entries of
                // the per-child tables (trio status, inherited stretches).
                eval(
                    b,
                    &format!(
                        "window.gateRows = null; PGPUI.getJSON('family.json').then(function (d) {{ var c = {}; function mine(r) {{ return [r.child, r.sample_a, r.sample_b, r.a, r.b].indexOf(c) >= 0; }} var rows = ['de_novo', 'compound_het', 'relationships', 'ibd_pairs'].reduce(function (n, k) {{ return n + (d[k] || []).filter(mine).length; }}, 0); var keys = ['trio_status', 'segments'].filter(function (k) {{ return !!d[k] && Object.prototype.hasOwnProperty.call(d[k], c); }}).length; window.gateRows = rows + keys; }}); true",
                        js_str(CHILD)
                    ),
                );
                wait(b, "window.gateRows !== null");
                let n = eval(b, "window.gateRows");
                (n == Value::Int(0), n)
            }
            Guard::Cleanup => {
                self.opened(b);
                type_search(b, "KCNH2");
                // Read in the same turn as the close: the document
                // replacement that follows would hide a kept query.
                let got = eval(
                    b,
                    &format!(
                        "PGPUI.closeAct({fid}); var s = document.querySelector('input[type=search]'); [PGPUI.gateOwnsPage() && !document.querySelector('input[type=search],.search-results,.search-status,canvas,svg,img'), s ? s.value : '']"
                    ),
                );
                (items(&got)[0] == Value::Bool(true), got)
            }
            Guard::Address => {
                go(
                    b,
                    &format!(
                        "{}?chr=17&start=43050000#f/closed-99",
                        run_url(self.base, CHILD)
                    ),
                );
                wait(b, CLOSED);
                let got = eval(b, "[location.search, location.hash]");
                (got == Value::Array(vec![s(""), s("#f/closed-1")]), got)
            }
            Guard::Stale => {
                go(b, &run_url(self.base, FAMILY));
                let control = format!("[data-gate-child={CHILD}] [data-opening^=closed]");
                wait(
                    b,
                    &format!("!!document.querySelector({})", js_str(&control)),
                );
                eval(
                    b,
                    &format!(
                        "window.gateOldControl = document.querySelector({}); true",
                        js_str(&control)
                    ),
                );
                self.changed(" REPLACEMENT_ACT_RECORD", false);
                eval(
                    b,
                    "window.gateRefreshDone = false; PGPUI.getJSON('family.json').then(function () { window.gateRefreshDone = true; }); true",
                );
                wait(b, "window.gateRefreshDone");
                eval(b, "gateOldControl.click(); true");
                // Read the outcome while the replacement is still what the
                // server answers: the redraw the click starts reads
                // results.json again.
                let act = js_str(&format!("[data-gate-child={CHILD}] .opened-act"));
                let opened = b.wait(
                    &format!("!!document.querySelector({act})"),
                    Duration::from_millis(1500),
                );
                let replacement = opened.is_some()
                    && eval(
                        b,
                        &format!(
                            "document.querySelector({act}).textContent.includes('REPLACEMENT_ACT_RECORD')"
                        ),
                    ) == Value::Bool(true);
                save(self.rj, self.data);
                (opened.is_none(), Value::Bool(replacement))
            }
            Guard::RefreshCleanup => {
                self.opened(b);
                // Hold document replacement so the guard must remove the old
                // report synchronously.
                eval(
                    b,
                    "window.gateTransferStarted = false; crypto.subtle.digest = function () { window.gateTransferStarted = true; return new Promise(function () {}); }; true",
                );
                self.changed(" REFRESHED_ACT_RECORD", true);
                eval(b, "PGPUI.getJSON('results.json'); true");
                wait(b, "window.gateTransferStarted");
                let ok = eval(
                    b,
                    "PGPUI.gateOwnsPage() && !document.querySelector('#results,canvas')",
                );
                let report = eval(b, "!!document.querySelector('#results .finding')");
                save(self.rj, self.data);
                (ok == Value::Bool(true), report)
            }
            Guard::Ticket => {
                // The opening ticket crosses one load only: after closing an
                // opened report again, the closed page keeps no ticket.
                self.opened(b);
                eval(b, &format!("PGPUI.closeAct({fid}); true"));
                wait(
                    b,
                    &format!("({CLOSED}) && document.readyState === 'complete'"),
                );
                let left = eval(
                    b,
                    "[sessionStorage.getItem('og-report-opening'), window.name]",
                );
                let l = items(&left);
                let ok = l[0] == Value::Null && !text(&l[1]).starts_with("og-report-opening:");
                (ok, left)
            }
            Guard::CloseAttrs => {
                // Read in the same turn as the close, before the document is
                // replaced.
                self.opened(b);
                let attrs = eval(b, &format!("PGPUI.closeAct({fid}); {LEFT}"));
                (items(&attrs).is_empty(), attrs)
            }
            Guard::Pagehide => {
                self.opened(b);
                let attrs = eval(
                    b,
                    &format!(
                        "window.dispatchEvent(new PageTransitionEvent('pagehide', {{ persisted: false }})); {LEFT}"
                    ),
                );
                let gone = eval(b, "!document.querySelector('#results')");
                (items(&attrs).is_empty() && gone == Value::Bool(true), attrs)
            }
            Guard::Back => {
                // Leave the open report for another page, then Back: from the
                // back-forward cache the page loads afresh, closed. The
                // static server sends no no-store, so Chrome may keep the
                // page in that cache.
                self.opened(b);
                // A page restored from that cache keeps its script state and
                // gets pageshow with persisted true; a reloaded one has neither.
                eval(
                    b,
                    "window.gateRestored = false; window.addEventListener('pageshow', function (e) { window.gateRestored = e.persisted; }); true",
                );
                b.go(&format!("{}/", self.base)).unwrap();
                wait(
                    b,
                    &format!("({READY}) && !!document.documentElement.dataset.pgpPeople"),
                );
                eval(b, "history.back(); true");
                b.wait(
                    &format!(
                        "location.pathname === '{CHILD_URL}' && document.readyState === 'complete'"
                    ),
                    WAIT,
                );
                let got = b
                    .wait(
                        &format!("({CLOSED}) && performance.getEntriesByType('navigation')[0].type"),
                        Duration::from_secs(10),
                    )
                    .unwrap_or_else(|| {
                        b.eval("[performance.getEntriesByType('navigation')[0].type, location.pathname, document.documentElement.dataset.pgpGated || '', window.gateRestored === true]")
                            .unwrap_or(Value::Null)
                    });
                (got == s("reload"), got)
            }
        }
    }
}

/// Whether a guard's probe saw the outcome its planted mistake causes, so a
/// control cannot go red for another reason (a timeout on a busy runner).
fn broke(guard: Guard, seen: &Value) -> bool {
    match guard {
        // The replacement act shows; the old report stays; one act opens.
        Guard::Stale | Guard::RefreshCleanup | Guard::Partial | Guard::RawView => {
            *seen == Value::Bool(true)
        }
        // The closed address keeps the map query.
        Guard::Address => text(&items(seen)[0]).contains("chr=17"),
        // The closed run's results.json reaches the view.
        Guard::Handoff => {
            matches!(seen, Value::Array(n) if matches!(n[..], [Value::Int(len)] if len > 0))
        }
        // Its rows or its per-child keys reach the family view.
        Guard::FamilyRows => matches!(seen, Value::Int(n) if *n > 0),
        // The search box keeps the query after the close.
        Guard::Cleanup => items(seen)[1] == s("KCNH2"),
        // The opening ticket is still stored.
        Guard::Ticket => items(seen)[0] != Value::Null,
        // Report attributes are left on the page.
        Guard::CloseAttrs | Guard::Pagehide => !items(seen).is_empty(),
        // Back reaches the child's page, which comes back from the
        // back-forward cache (pageshow persisted) open, not reloaded and
        // closed.
        Guard::Back => {
            matches!(&seen, Value::Array(a) if a.len() == 4 && a[0] != s("reload") && a[1] == s(CHILD_URL) && a[2] != s("true") && a[3] == Value::Bool(true))
        }
    }
}

/// test_live_demo.py's gate controls: each guard of the gate holds, and
/// fails when its implementation in core.js is deliberately broken, with
/// the outcome that mistake causes.
fn gate_controls(f: &mut Fixture, base: &str, core: &Path, source: &str) {
    let rj = f.site.join(format!("runs/{CHILD}/results.json"));
    let raw = fs::read(&rj).unwrap();
    let mut data = doc(&rj);
    obj(&mut data).insert("future_payload_field".into(), s("FUTURE_GATE_SENTINEL"));
    let acts = acts_of(&data);
    write_acts(&rj, &data, &acts[..2], false);
    let data = doc(&rj);
    let guards = Guards {
        base,
        core,
        rj: &rj,
        data: &data,
        fid: text(field(&acts[0], "id")).to_owned(),
    };
    let gate_open =
        "if (ownsPage || closedRun(requestURL)) { return new Promise(function () {}); }";
    let gate_late = "return got.then(function (response) {\n        if (ownsPage) { return new Promise(function () {}); }";
    let mutations: [(&str, Edits, Guard); 14] = [
        (
            "a stale control opens a replacement act",
            vec![(
                "if (!same(byId(active, id), byId(rec, id))) { changed(); return; }",
                "",
            )],
            Guard::Stale,
        ),
        (
            "a refresh leaves the old report visible",
            vec![(
                "if (shut(current)) { ownsPage = true; eventSources.forEach(function (es) { es.close(); }); cleanPage(); transfer();",
                "if (shut(current)) { transfer();",
            )],
            Guard::RefreshCleanup,
        ),
        (
            "the closed address keeps a map query",
            vec![(
                "location.pathname + to",
                "location.pathname + location.search + to",
            )],
            Guard::Address,
        ),
        (
            "one act is released early",
            vec![(
                "rec.acts.some(function (f) { return !rec.opened[f.id]; })",
                "rec.acts.filter(function (f) { return !rec.opened[f.id]; }).length > 1",
            )],
            Guard::Partial,
        ),
        (
            "one view reads the raw payload",
            vec![(
                "add(document.body, wrap);",
                "add(document.body, wrap); jsonAt(current.url + 'results.json').then(function (r) { add(wrap, el('p', {}, JSON.stringify(r))); document.documentElement.setAttribute('data-gate-raw-read', 'done'); });",
            )],
            Guard::RawView,
        ),
        (
            "a view receives a closed run's file",
            vec![
                (gate_open, ""),
                (gate_late, "return got.then(function (response) {"),
            ],
            Guard::Handoff,
        ),
        (
            "the family view receives a gated child's trio rows",
            vec![(
                "copy((raw[key] || []).filter(function (r) { return !gated(r.child); }))",
                "copy(raw[key] || [])",
            )],
            Guard::FamilyRows,
        ),
        (
            "the family view receives a gated child's relationship and shared-segment rows",
            vec![(FREE_CLOSED, FREE_LEAKED)],
            Guard::FamilyRows,
        ),
        (
            "the family view receives a gated child's trio status",
            vec![(
                "key === 'trio_status' ? gated(k) || parentGated(k) : gated(k)",
                "key === 'trio_status' ? false : gated(k)",
            )],
            Guard::FamilyRows,
        ),
        (
            "closing skips the cleanup",
            vec![(
                "ownsPage = true; eventSources.forEach(function (es) { es.close(); }); cleanPage(); transfer(); return true;",
                "return true;",
            )],
            Guard::Cleanup,
        ),
        (
            "the opening ticket is kept after it is read",
            vec![("    sessionStorage.removeItem(TRANSFER);\n", "")],
            Guard::Ticket,
        ),
        (
            "closing keeps the page attributes",
            vec![(
                "if (a.name.indexOf('data-pgp-') === 0 || a.name === 'data-finding') { root.removeAttribute(a.name); }",
                "",
            )],
            Guard::CloseAttrs,
        ),
        (
            "pagehide keeps the open report",
            vec![(
                "ownsPage = true; cleaned = true; eventSources.forEach(function (es) { es.close(); }); cleanPage(); }",
                "ownsPage = true; cleaned = true; eventSources.forEach(function (es) { es.close(); }); }",
            )],
            Guard::Pagehide,
        ),
        (
            "the back-forward cache shows the cleared report again",
            vec![("if (e.persisted && cleaned) { location.reload(); }", "")],
            Guard::Back,
        ),
    ];
    let b = &mut f.browser;
    for (what, edits, guard) in mutations {
        let (ok, seen) = guards.probe(b, guard);
        assert!(ok, "gate guard {guard:?} ({what}) failed: {seen:?}");
        pass(&format!("gate guard: {what} never happens"));
        let mut planted = source.to_owned();
        for (from, to) in edits {
            planted = replace_once(&planted, from, to);
        }
        fs::write(core, planted).unwrap();
        let (ok, seen) = guards.probe(b, guard);
        assert!(!ok, "the gate control went unnoticed: {what}");
        assert!(
            broke(guard, &seen),
            "the gate control failed without its broken outcome: {what}: {seen:?}"
        );
        red(&format!("core.js where {what}"), &seen);
        fs::write(core, source).unwrap();
    }
    fs::write(rj, raw).unwrap();
}

fn add_second_child(site: &Path) -> usize {
    let family_file = site.join(format!("runs/{FAMILY}/family.json"));
    let mut family = doc(&family_file);
    let f = obj(&mut family);
    let state = obj(f.get_mut("state").unwrap());
    for (key, name) in [("members", "sample"), ("children", "child")] {
        let entries = array(state.get_mut(key).unwrap());
        let mut cloned = entries
            .iter()
            .find(|x| text(field(x, name)) == CHILD)
            .unwrap()
            .clone();
        obj(&mut cloned).insert(name.into(), Value::Str(OTHER.into()));
        entries.push(cloned);
    }
    let runs = obj(f.get_mut("runs").unwrap());
    let mut run = runs[CHILD].clone();
    obj(&mut run).insert("url".into(), Value::Str(format!("../../runs/{OTHER}/")));
    runs.insert(OTHER.into(), run);
    let side = obj(f.get_mut("side_by_side").unwrap());
    let members = array(side.get_mut("members").unwrap());
    let index = members
        .iter()
        .position(|sample| text(sample) == CHILD)
        .unwrap();
    members.push(Value::Str(OTHER.into()));
    for row in array(side.get_mut("rows").unwrap()) {
        for key in ["values", "stale"] {
            let values = array(obj(row).get_mut(key).unwrap());
            values.push(values[index].clone());
        }
    }
    save(&family_file, &family);

    // Start from the same document so both children have identical finding IDs.
    let mut results = doc(&site.join(format!("runs/{CHILD}/results.json")));
    let results_obj = obj(&mut results);
    results_obj.insert("sample".into(), Value::Str(OTHER.into()));
    obj(results_obj.get_mut("person").unwrap()).insert("sample".into(), Value::Str(OTHER.into()));
    let acts = array(results_obj.get_mut("findings").unwrap())
        .iter()
        .filter(|f| text(field(f, "level")) == "act")
        .count();
    assert!(acts >= 2);
    let second = site.join(format!("runs/{OTHER}"));
    fs::create_dir(&second).unwrap();
    save(&second.join("results.json"), &results);
    acts
}

// The child's own values on the family page's side by side, in the core.js
// that leaves them in.
const SIDE_CLOSED: &str =
    "values: row.values.map(function (v, i) { return gated(sb.members[i]) ? 'closed' : v; }),";
const SIDE_LEAKED: &str = "values: row.values.map(function (v, i) { return v; }),";
// The relationship and shared-segment rows that name a gated child, in the
// core.js that leaves them in.
const FREE_CLOSED: &str = "function free(r) { return !gated(r.sample_a) && !gated(r.sample_b) && !gated(r.a) && !gated(r.b); }";
const FREE_LEAKED: &str = "function free(r) { return true; }";
const CHANGED: &str = "changed value";

/// Every value inside `v` changed: text gets a suffix, numbers grow by 37.
fn change(v: &mut Value) {
    match v {
        Value::Str(x) => x.push_str(" changed"),
        Value::Int(n) => *n += 37,
        Value::Float(x) => *x += 37.0,
        Value::Array(xs) => xs.iter_mut().for_each(change),
        Value::Object(o) => o.values_mut().for_each(change),
        Value::Null | Value::Bool(_) => {}
    }
}

/// Every field of `row` changed but those named in `keep`.
fn change_but(row: &mut Value, keep: &[&str]) {
    for (k, v) in obj(row).iter_mut() {
        if !keep.contains(&k.as_str()) {
            change(v);
        }
    }
}

/// The family page while DEMO-003 is closed: the child's cells of the side
/// by side, and the whole page (every chapter, and the lead that sums them
/// up). The closed controls are drawn in the same order in every load.
fn closed_family(b: &mut Browser, base: &str) -> (Value, Value) {
    let seed = install(b, &SEEDED.replace("DRAW", "0"));
    go(b, &run_url(base, FAMILY));
    wait(
        b,
        &format!(
            "({READY}) && document.querySelectorAll('[data-gate-child={CHILD}] [data-opening^=closed]').length > 0 && !!document.querySelector('#ch-side-by-side td')"
        ),
    );
    // Let the last redraw settle before the page is read.
    thread::sleep(Duration::from_millis(500));
    let got = eval(
        b,
        &format!(
            "[[...document.querySelectorAll('#ch-side-by-side td[data-label=\"{CHILD}\"]')].map(td => td.textContent), document.getElementById('family').outerHTML]"
        ),
    );
    uninstall(b, &seed);
    let got = items(&got);
    (got[0].clone(), got[1].clone())
}

/// Every cell of the child's column reads "closed", and there is one.
fn all_closed(cells: &Value) -> bool {
    !items(cells).is_empty() && items(cells).iter().all(|c| *c == s("closed"))
}

/// The family page while every field of a gated child in family.json
/// changes, and fields a producer may add tomorrow appear: the child's
/// side-by-side cells read "closed", and no chapter of the page changes
/// (test_live_demo.py gate_checks, its family half, on main eba8400b).
fn gate_family_payload(f: &mut Fixture, base: &str, core: &Path, source: &str) {
    let family_file = f.site.join(format!("runs/{FAMILY}/family.json"));
    let saved = fs::read(&family_file).unwrap();
    let b = &mut f.browser;
    let (cells, page) = closed_family(b, base);
    assert!(
        all_closed(&cells),
        "family page: the closed child's side-by-side cells: {cells:?}"
    );

    let mut family = doc(&family_file);
    let fam = obj(&mut family);
    let side = obj(fam.get_mut("side_by_side").unwrap());
    let at = items(&side["members"])
        .iter()
        .position(|m| *m == s(CHILD))
        .unwrap();
    let width = items(&side["members"]).len();
    let rows = array(side.get_mut("rows").unwrap());
    for row in rows.iter_mut() {
        let row = obj(row);
        array(row.get_mut("values").unwrap())[at] = s(CHANGED);
        array(row.get_mut("stale").unwrap())[at] = Value::Bool(true);
    }
    // A row only the child has a value in.
    let mut only = vec![s("."); width];
    only[at] = s(CHANGED);
    rows.push(
        loads(&format!(
            r#"{{"group": "Blood groups", "label": "Changed system", "values": {}, "stale": {}}}"#,
            dumps(&Value::Array(only), &Default::default()).unwrap(),
            dumps(
                &Value::Array(vec![Value::Bool(false); width]),
                &Default::default()
            )
            .unwrap()
        ))
        .unwrap(),
    );
    obj(fam.get_mut("trio_status").unwrap()).insert(
        CHILD.into(),
        loads(r#"{"de_novo": {"assessed": false, "why": "changed value"}, "compound_het": {"assessed": null, "why": "changed value"}}"#).unwrap(),
    );
    for key in ["de_novo", "compound_het"] {
        let trio = array(fam.get_mut(key).unwrap());
        let mut extra = None;
        for row in trio.iter_mut() {
            if *field(row, "child") != s(CHILD) {
                continue;
            }
            for (k, v) in obj(row).iter_mut() {
                if k.as_str() != "child"
                    && let Value::Str(x) = v
                {
                    x.push_str(" changed");
                }
            }
            extra = Some(row.clone());
        }
        trio.push(extra.unwrap_or_else(|| panic!("the demo has no {key} row of {CHILD}")));
    }
    array(
        obj(fam.get_mut("segments").unwrap())
            .get_mut(CHILD)
            .unwrap(),
    )
    .truncate(1);
    // The child's relationship and shared-segment rows, its trio files, its
    // entry among the children, the plan rows and reports that name it, its
    // run past what any run's page shows, and fields added tomorrow.
    let mut named = 0;
    for key in ["relationships", "ibd_pairs"] {
        for row in array(fam.get_mut(key).unwrap()) {
            if ["sample_a", "sample_b", "a", "b"]
                .iter()
                .any(|k| obj(row).get(*k) == Some(&s(CHILD)))
            {
                change_but(row, &["sample_a", "sample_b", "a", "b"]);
                named += 1;
            }
        }
    }
    assert!(named > 0, "the demo has no relationship row of {CHILD}");
    // Which trio files the folder holds, by file name: the one child's.
    for key in ["trio_files", "trio_file_state"] {
        change(fam.get_mut(key).unwrap());
    }
    for key in ["de_novo", "compound_het"] {
        array(fam.get_mut(key).unwrap()).push(
            loads(&format!(
                r#"{{"child": "{CHILD}", "future_child_field": "new child record"}}"#
            ))
            .unwrap(),
        );
    }
    let state = obj(fam.get_mut("state").unwrap());
    for child in array(state.get_mut("children").unwrap()) {
        if *field(child, "child") == s(CHILD) {
            change_but(child, &["child"]);
        }
    }
    let plan = array(state.get_mut("plan").unwrap());
    for row in plan.iter_mut() {
        if matches!(obj(row).get("members"), Some(Value::Str(m))
            if m.split([' ', ',', ';', '\t']).any(|m| m == CHILD))
        {
            change_but(row, &["members"]);
        }
    }
    plan.push(
        loads(&format!(
            r#"{{"analysis": "future_child_analysis", "members": "{CHILD}", "decision": "runs", "reason": "new"}}"#
        ))
        .unwrap(),
    );
    array(state.get_mut("reports").unwrap()).push(
        loads(&format!(
            r#"{{"url": "{CHILD}_future_report.html", "title": "{CHILD} future report"}}"#
        ))
        .unwrap(),
    );
    let run = obj(fam.get_mut("runs").unwrap()).get_mut(CHILD).unwrap();
    let progress = [
        "name",
        "kind",
        "url",
        "state",
        "started_utc",
        "last_heard_utc",
        "finished_utc",
        "planned",
        "done",
        "running",
        "failed",
        "child",
    ];
    change_but(run, &progress);
    obj(run).insert("future_child_field".into(), s("new child value"));
    fam.insert("future_child_payload".into(), s("new child value"));
    save(&family_file, &family);

    let (changed_cells, changed_page) = closed_family(b, base);
    assert!(
        all_closed(&changed_cells) && changed_page == page,
        "family page: a closed child's family.json data changed the page: {changed_cells:?}"
    );
    pass(
        "family page: a closed child's side-by-side cells read closed, and every field of it in family.json (side by side, trio status and rows, segments, relationships, trio files, plan, reports, run, new fields) leaves every chapter as it was",
    );

    fs::write(core, replace_once(source, FREE_CLOSED, FREE_LEAKED)).unwrap();
    let (_, leaked_page) = closed_family(b, base);
    fs::write(core, source).unwrap();
    assert!(
        leaked_page != page,
        "the relationship-row control went unnoticed"
    );
    red(
        "core.js that keeps a closed child's relationship and shared-segment rows on the family page",
        &"the page changed",
    );

    fs::write(core, replace_once(source, SIDE_CLOSED, SIDE_LEAKED)).unwrap();
    let (leaked, leaked_page) = closed_family(b, base);
    fs::write(core, source).unwrap();
    assert!(
        items(&leaked).contains(&s(CHANGED)) && leaked_page != page,
        "the side-by-side control went unnoticed: {leaked:?}"
    );
    red(
        "core.js that shows a closed child's side-by-side values on the family page",
        &leaked,
    );
    fs::write(&family_file, saved).unwrap();
}

fn gate_family_independence(f: &mut Fixture, base: &str, core: &Path, source: &str) {
    let acts = add_second_child(&f.site);
    let b = &mut f.browser;
    let other = "document.querySelector('[data-gate-child=DEMO-004]')?.outerHTML || ''";
    let buttons = "document.querySelectorAll('[data-gate-child=DEMO-003] [data-opening^=closed]')";
    let mut exercise = || {
        go(b, &run_url(base, FAMILY));
        wait(b, &format!("{buttons}.length === {acts}"));
        let before = wait(b, other);
        for remaining in (0..acts).rev() {
            eval(b, &format!("{buttons}[0].click(); true"));
            wait(b, &format!("{buttons}.length === {remaining}"));
        }
        wait(
            b,
            "document.querySelector('[data-gate-child=DEMO-003]') === null",
        );
        (before, eval(b, other))
    };
    let (before, after) = exercise();
    assert_eq!(before, after, "opening DEMO-003 changed DEMO-004");
    pass("gate: two children with the same act IDs; opening one leaves the other closed");
    // Plant a shared record for the two children, modeling an opening keyed by finding ID alone.
    let planted = replace_once(
        source,
        "familyGates = recs; return projectFamily(data, recs);",
        "recs['DEMO-004'] = recs['DEMO-003']; familyGates = recs; return projectFamily(data, recs);",
    );
    fs::write(core, planted).unwrap();
    let (before, after) = exercise();
    assert_ne!(
        before, after,
        "the shared finding ID control must change DEMO-004"
    );
    red(
        "a family opening keyed by finding ID alone",
        &"DEMO-004 changed",
    );
    fs::write(core, source).unwrap();
}

// A summary that gives an adult a result to act on (BRCA1) and two symptom
// rankings, one of the same gene.
const SYNTH: &str = r#"{"sections": {"cpsr": {"state": "ok", "data": {"plp": [{"gene": "BRCA1", "variant": "c.68_69del", "final": "Pathogenic", "clinvar": "Pathogenic"}]}}, "exomiser": {"state": "ok", "data": {"genes": [{"gene": "BRCA1", "rank": 1, "exomiser_score": 0.9, "moi": "AD"}, {"gene": "TTN", "rank": 2, "exomiser_score": 0.8, "moi": "AR"}]}}, "clinvar": {"state": "ok", "data": {"hits": []}}, "coverage": {"state": "ok", "data": {"mean_depth": 31.2}}}}"#;
const NO_FINDING: &str = "No finding in what these steps assessed needs a doctor.";

// The gate transparency capture: the DOM with each element's attributes
// sorted by name (the page's scripts set them as their files arrive, in any
// order), the address, the canvases, the resources, every read with its
// options, and the reads of the files the gate reads first, in their order.
const CAPTURE: &str = r#"(() => {
  const clone = document.documentElement.cloneNode(true);
  clone.querySelectorAll('[data-pgp-map-paint-ms],[data-pgp-map-parse-ms]').forEach(node => {
    node.removeAttribute('data-pgp-map-paint-ms');
    node.removeAttribute('data-pgp-map-parse-ms');
  });
  for (const node of [clone, ...clone.querySelectorAll('*')]) {
    const attrs = Array.from(node.attributes).sort((a, b) => a.name < b.name ? -1 : a.name > b.name ? 1 : 0);
    for (const attr of attrs) node.removeAttributeNode(attr);
    for (const attr of attrs) node.setAttributeNodeNS(attr);
  }
  const reads = window.ogRequests || [];
  const first = reads.filter(r => /\/(state|frames|family|results)\.json$|\/events$/.test(r[1]));
  return [clone.outerHTML, location.href, Array.from(document.querySelectorAll('canvas'), c => c.toDataURL()),
    performance.getEntriesByType('resource').map(r => new URL(r.name).pathname)
      .filter(path => path !== '/favicon.ico').sort(),
    reads.map(r => JSON.stringify(r)).sort(), first];
})()"#;

#[derive(Default)]
struct How {
    apoe: bool,
    search: Option<String>,
    no_family: bool,
}

fn capture(
    b: &mut Browser,
    base: &str,
    route: &str,
    how: &How,
    timeout: Duration,
) -> Option<Vec<Value>> {
    go(b, &format!("{base}{route}"));
    if route == "/" {
        b.wait(
            &format!("({READY}) && !!document.documentElement.dataset.pgpPeople"),
            timeout,
        )?;
    } else if route.contains("/finding/") {
        b.wait(
            &format!("({READY}) && !!document.querySelector('h1')"),
            timeout,
        )?;
    } else {
        b.wait(
            "document.documentElement.dataset.pgpResults && !!document.querySelector('[data-pgp-map=drawn]')",
            timeout,
        )?;
        b.wait(
            &format!(
                "({} || !!document.querySelector('#rail a[href*={FAMILY}]')) && !!document.documentElement.dataset.pgpCharts",
                how.no_family
            ),
            timeout,
        )?;
    }
    if how.apoe {
        // A redraw while the address names a place: the contents rail must
        // come back as it did.
        eval(
            b,
            "document.querySelector('details[data-gate=apoe] summary').click(); true",
        );
    }
    if let Some(term) = &how.search {
        type_search(b, term);
    }
    // Let the last redraw settle before the page is read.
    thread::sleep(Duration::from_millis(500));
    let Value::Array(got) = eval(b, CAPTURE) else {
        panic!("missing capture")
    };
    Some(got)
}

/// A term of the demo's search index with a mention at a locus.
fn located_term(site: &Path) -> String {
    let index = doc(&site.join(format!("runs/{ADULT}/{ADULT}_search_index.json")));
    let Value::Object(terms) = field(&index, "terms") else {
        panic!("no terms")
    };
    let mut names: Vec<&String> = terms.keys().collect();
    names.sort();
    let located = |item: &Value| match item {
        Value::Object(o) => match o.get("mentions") {
            Some(Value::Array(ms)) => ms.iter().any(|m| match m {
                Value::Object(m) => m.get("locus").is_some_and(pycompat::value::truthy),
                _ => false,
            }),
            _ => false,
        },
        _ => false,
    };
    names
        .into_iter()
        .find(|t| items(&terms[*t]).iter().any(located))
        .expect("a located term")
        .clone()
}

#[derive(Clone, Copy, PartialEq)]
enum Prep {
    Nothing,
    NoAct,
    Synthetic,
    NoFamily,
    NoActNoFamily,
}

struct Files {
    adult: PathBuf,
    child: PathBuf,
    runs: PathBuf,
    family: PathBuf,
    saved: Vec<(PathBuf, Vec<u8>)>,
}

impl Files {
    fn new(site: &Path) -> Self {
        let adult = site.join(format!("runs/{ADULT}/results.json"));
        let child = site.join(format!("runs/{CHILD}/results.json"));
        let runs = site.join("runs.json");
        let family = site.join(format!("runs/{FAMILY}/family.json"));
        let saved = [&adult, &child, &runs, &family]
            .iter()
            .map(|p| ((*p).clone(), fs::read(p).unwrap()))
            .collect();
        Files {
            adult,
            child,
            runs,
            family,
            saved,
        }
    }

    fn restore(&self) {
        for (p, raw) in &self.saved {
            fs::write(p, raw).unwrap();
        }
    }

    fn run_entry<'a>(runs: &'a mut Value, name: &str) -> &'a mut indexmap::IndexMap<String, Value> {
        let all = array(obj(runs).get_mut("runs").unwrap());
        obj(all
            .iter_mut()
            .find(|r| text(field(r, "name")) == name)
            .unwrap())
    }

    /// The child's results and its runs.json entry with no act.
    fn no_act(&self) {
        let mut d = doc(&self.child);
        let o = obj(&mut d);
        array(o.get_mut("findings").unwrap()).retain(|f| !is_act(f));
        obj(o.get_mut("levels").unwrap()).insert("act".into(), Value::Int(0));
        save(&self.child, &d);
        let mut idx = doc(&self.runs);
        let r = Self::run_entry(&mut idx, CHILD);
        obj(r.get_mut("levels").unwrap()).insert("act".into(), Value::Int(0));
        array(r.get_mut("top").unwrap()).retain(|t| !is_act(t));
        r.insert("lead".into(), s(NO_FINDING));
        r.insert("lead_closed".into(), s(NO_FINDING));
        save(&self.runs, &idx);
    }

    /// A fresh producer run for the adult: results and home entry from SYNTH.
    fn synthetic(&self, og: &Path) {
        self.no_act();
        save(&self.adult, &og_results(og, SYNTH, &["--sample", ADULT]));
        let o = og_results(og, SYNTH, &["--overview"]);
        let mut idx = doc(&self.runs);
        let r = Self::run_entry(&mut idx, ADULT);
        for key in [
            "attention",
            "levels",
            "top",
            "lead",
            "lead_closed",
            "screened",
        ] {
            r.insert(key.into(), field(&o, key).clone());
        }
        save(&self.runs, &idx);
    }

    fn prepare(&self, prep: Prep, og: &Path) {
        self.restore();
        match prep {
            Prep::Nothing => {}
            Prep::NoAct => self.no_act(),
            Prep::Synthetic => self.synthetic(og),
            // A family that cannot be read (its family.json answers 404).
            Prep::NoFamily => fs::remove_file(&self.family).unwrap(),
            Prep::NoActNoFamily => {
                self.no_act();
                fs::remove_file(&self.family).unwrap();
            }
        }
    }
}

/// Where two captures differ, a little of each side.
fn capture_difference(x: &[Value], y: &[Value]) -> String {
    let mut out = Vec::new();
    for (i, label) in [
        "DOM",
        "address",
        "canvas",
        "resources",
        "reads",
        "first reads",
    ]
    .iter()
    .enumerate()
    {
        if x[i] == y[i] {
            continue;
        }
        if i == 0 {
            let at = text(&x[0])
                .split("><")
                .zip(text(&y[0]).split("><"))
                .find(|(a, b)| a != b);
            out.push(format!("{label}: {at:?}"));
        } else if i >= 3 {
            let (a, b) = (items(&x[i]), items(&y[i]));
            let only = |p: &[Value], q: &[Value]| -> Vec<Value> {
                p.iter().filter(|v| !q.contains(v)).cloned().collect()
            };
            out.push(format!(
                "{label}: without the gate only {:?}, with it only {:?}, {a:?} / {b:?}",
                only(a, b),
                only(b, a)
            ));
        } else {
            out.push(format!("{label}: {:?} / {:?}", x[i], y[i]));
        }
    }
    out.join("; ")
}

/// test_live_demo.py's gate transparency: the gate leaves adult and no-act
/// pages equal to the renderer without it: DOM, attributes, address, canvas
/// and the requests, counted, each with its fetch options, and the files the
/// gate reads first in their order.
fn gate_transparency(f: &mut Fixture, base: &str, core: &Path, source: &str) {
    let files = Files::new(&f.site);
    let bare = bare_core(source);
    let term = located_term(&f.site);
    let first = items(field(&doc(&files.adult), "findings"))
        .iter()
        .find(|x| !is_act(x))
        .map(|x| text(field(x, "id")).to_owned())
        .unwrap();
    let adult = format!("/runs/{ADULT}/");
    let child = format!("/runs/{CHILD}/");
    let finding = format!("/runs/{ADULT}/finding/{first}/");
    let cases = [
        ("an adult", Prep::Nothing, adult.clone(), How::default()),
        (
            "an adult, a redraw at a section address",
            Prep::Nothing,
            format!("{adult}#s3"),
            How {
                apoe: true,
                ..How::default()
            },
        ),
        (
            "an adult's finding page",
            Prep::Nothing,
            finding.clone(),
            How::default(),
        ),
        (
            "a child with no act, a search typed",
            Prep::NoAct,
            child.clone(),
            How {
                search: Some(term.clone()),
                ..How::default()
            },
        ),
        (
            "the home page with no gated child",
            Prep::NoAct,
            "/".into(),
            How::default(),
        ),
        (
            "a fresh adult with an act and a same-gene ranking",
            Prep::Synthetic,
            adult.clone(),
            How::default(),
        ),
        (
            "an adult whose family cannot be read",
            Prep::NoFamily,
            adult.clone(),
            How {
                no_family: true,
                ..How::default()
            },
        ),
        (
            "an adult's finding page whose family cannot be read",
            Prep::NoFamily,
            finding,
            How::default(),
        ),
        (
            "a child with no act whose family cannot be read",
            Prep::NoActNoFamily,
            child.clone(),
            How {
                no_family: true,
                ..How::default()
            },
        ),
    ];
    let b = &mut f.browser;
    let log = install(b, REQUEST_LOG);
    for (label, prep, route, how) in &cases {
        files.prepare(*prep, &f.og);
        let mut pair = Vec::new();
        for core_text in [bare.as_str(), source] {
            fs::write(core, core_text).unwrap();
            pair.push(
                capture(b, base, route, how, WAIT)
                    .unwrap_or_else(|| panic!("{label}: {route} never drew")),
            );
        }
        assert!(
            pair[0] == pair[1],
            "gate transparency: {label} differs: {}",
            capture_difference(&pair[0], &pair[1])
        );
        pass(&format!(
            "gate transparency: {label} has identical DOM, attributes, address, canvas and requests"
        ));
        if how.search.is_some() {
            // A child's page keeps results closed until opened, with or
            // without results to act on, and says so under the hits.
            assert!(
                text(&pair[1][0]).contains(SEARCH_CLOSED),
                "the search on a child with no act lost its closed note ({term})"
            );
            pass(
                "gate transparency: a search on a child with no act says that closed results are left out",
            );
        }
    }
    files.prepare(Prep::NoAct, &f.og);
    let reference = capture(b, base, &child, &How::default(), WAIT);
    assert!(reference.is_some(), "the no-act child never drew");
    for (what, from, to, draws) in [
        (
            "closing an adult or no-act report",
            "ownsPage = shut(current);",
            "ownsPage = true;",
            false,
        ),
        (
            "reading a file twice",
            "var kept = shared[requestURL];",
            "var kept = null;",
            true,
        ),
    ] {
        fs::write(core, replace_once(source, from, to)).unwrap();
        let planted = capture(b, base, &child, &How::default(), Duration::from_secs(4));
        assert!(
            planted != reference,
            "transparency control went unnoticed: {to}"
        );
        // A page that reads a file twice still draws: what differs is its
        // reads, more of them than the reference made.
        if draws {
            let want = items(&reference.as_ref().unwrap()[4]).len();
            let got = planted.as_ref().map(|p| items(&p[4]).len());
            assert!(
                got.is_some_and(|n| n > want),
                "transparency control {what}: {got:?} reads, the reference {want}"
            );
        }
        let seen = match &planted {
            None => "the page never drew".to_owned(),
            Some(p) => capture_difference(reference.as_ref().unwrap(), p),
        };
        red(&format!("gate transparency: core.js {what}"), &seen);
    }
    uninstall(b, &log);
    files.restore();
    fs::write(core, source).unwrap();
}

fn slow_open(b: &mut Browser, base: &str) -> bool {
    go(b, &run_url(base, CHILD));
    wait(b, CLOSED);
    let closed = "document.querySelectorAll('[data-opening^=closed]')";
    for _ in 0..10 {
        if eval(b, &format!("{closed}.length")) == Value::Int(1) {
            break;
        }
        let Value::Int(before) = eval(b, &format!("{closed}.length")) else {
            panic!("missing closed act count")
        };
        eval(b, &format!("{closed}[0].click(); true"));
        wait(b, &format!("{closed}.length === {before} - 1"));
    }
    wait(b, &format!("{closed}.length === 1"));
    call(
        b,
        "Network.emulateNetworkConditionsByRule",
        &format!(
            r#"{{"offline":false,"matchedNetworkConditions":[
            {{"urlPattern":{},"latency":4000,"downloadThroughput":-1,"uploadThroughput":-1}},
            {{"urlPattern":{},"latency":6000,"downloadThroughput":-1,"uploadThroughput":-1}}
        ]}}"#,
            js_str(&format!("{base}{CHILD_URL}")),
            js_str(&format!("{base}{CHILD_URL}results.json")),
        ),
    );
    eval(
        b,
        &format!("setTimeout(() => {closed}[0].click(), 0); true"),
    );
    let outcome = b.wait(
        &format!("location.pathname === '{CHILD_URL}' && (({OPEN}) ? 'open' :           (performance.getEntriesByType('navigation')[0]?.responseStart >= 3500 && ({CLOSED}) ? 'closed' : ''))"),
        Duration::from_secs(25),
    );
    call(
        b,
        "Network.emulateNetworkConditionsByRule",
        r#"{"offline":false,"matchedNetworkConditions":[]}"#,
    );
    let outcome = outcome.expect("the slow reload did not reach a gate decision");
    assert_eq!(
        eval(
            b,
            r#"performance.getEntriesByType('navigation')[0].responseStart >= 3500 &&
        performance.getEntriesByType('resource').some(r => r.name.endsWith('/results.json') && r.responseStart - r.startTime >= 5500)"#
        ),
        Value::Bool(true),
        "Chrome must observe both delays"
    );
    text(&outcome) == "open"
}

fn gate_slow_release(b: &mut Browser, base: &str, core: &Path, source: &str) {
    assert!(
        slow_open(b, base),
        "the slow reload did not release the report"
    );
    pass(
        "gate release: a reload whose page takes 4 s and whose results.json takes 6 s still releases the report",
    );
    let planted = replace_once(
        source,
        "sessionStorage.setItem(TRANSFER, JSON.stringify({ path: target,",
        "sessionStorage.setItem(TRANSFER, JSON.stringify({ at: Date.now(), path: target,",
    );
    let planted = replace_once(
        &planted,
        "t.path === location.pathname) { ticket = t; }",
        "t.path === location.pathname && Date.now() - t.at < 1500) { ticket = t; }",
    );
    fs::write(core, planted).unwrap();
    assert!(
        !slow_open(b, base),
        "the 1.5 s ticket control still released the report"
    );
    red(
        "an opening accepted only within 1.5 s",
        &"the report stayed closed",
    );
    fs::write(core, source).unwrap();
}

#[test]
fn demo_report_gate() {
    let Some(mut f) = Fixture::build("demo-gate") else {
        return;
    };
    let server = StaticServer::start(&f.site, &[]).unwrap();
    let base = server.url("").trim_end_matches('/').to_owned();
    let core = f.site.join("assets/core.js");
    let source = fs::read_to_string(&core).unwrap();
    gate_handles(&mut f, &base, &core, &source);
    gate_order(&mut f, &base, &core, &source);
    gate_controls(&mut f, &base, &core, &source);
    gate_transparency(&mut f, &base, &core, &source);
    gate_slow_release(&mut f.browser, &base, &core, &source);
    gate_family_payload(&mut f, &base, &core, &source);
    gate_family_independence(&mut f, &base, &core, &source);
}

struct LiveServe {
    child: Child,
    base: String,
}

impl LiveServe {
    fn start(og: &Path, root: &Path, state: &Path, assets: &Path) -> Self {
        let mut server = Self {
            child: Command::new(og)
                .args(["live", "serve", "--port", "0", "--quiet", "--root"])
                .arg(root)
                .arg("--state-dir")
                .arg(state)
                .arg("--assets")
                .arg(assets)
                .arg("--steps")
                .arg(repo_root().join("assets/steps.tsv"))
                .env_remove("OWNGENOME_LIVE_PASSWORD_HASH")
                .stdin(Stdio::null())
                .stdout(Stdio::null())
                .stderr(Stdio::piped())
                .spawn()
                .unwrap(),
            base: String::new(),
        };
        let stderr = server.child.stderr.take().unwrap();
        let mut lines = BufReader::new(stderr).lines();
        while let Some(line) = lines.next() {
            let line = line.unwrap();
            if let Some(url) = line.strip_prefix("og live serve: browser ") {
                server.base = url
                    .split_whitespace()
                    .next()
                    .unwrap()
                    .trim_end_matches('/')
                    .to_owned();
                // Keep reading what the server says, so it never writes to a
                // closed pipe.
                thread::spawn(move || {
                    for line in lines.map_while(Result::ok) {
                        eprintln!("{line}");
                    }
                });
                return server;
            }
            eprintln!("{line}");
        }
        panic!("og live serve exited before listening");
    }
}

impl Drop for LiveServe {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

/// One GET to `og live serve`: the status and the body.
fn http_get(base: &str, path: &str) -> (u16, Vec<u8>) {
    let host = base.strip_prefix("http://").unwrap();
    let mut stream = std::net::TcpStream::connect(host).unwrap();
    stream.set_read_timeout(Some(WAIT)).unwrap();
    write!(
        stream,
        "GET {path} HTTP/1.1\r\nHost: {host}\r\nConnection: close\r\n\r\n"
    )
    .unwrap();
    let mut raw = Vec::new();
    stream.read_to_end(&mut raw).unwrap();
    let end = raw
        .windows(4)
        .position(|w| w == b"\r\n\r\n")
        .unwrap_or_else(|| panic!("GET {path}: no end of headers"));
    let head = String::from_utf8_lossy(&raw[..end]);
    let status = head.split(' ').nth(1).unwrap().parse().unwrap();
    (status, raw[end + 4..].to_vec())
}

/// The ids of a run's findings, from its results.json on `og live serve`,
/// results to act on first.
fn finding_ids(base: &str, run: &str) -> Vec<String> {
    let (status, body) = http_get(base, &format!("/runs/{run}/results.json"));
    assert_eq!(status, 200, "{run}'s results.json");
    let results = loads(&String::from_utf8(body).unwrap()).unwrap();
    let mut found = items(field(&results, "findings")).to_vec();
    found.sort_by_key(|f| !is_act(f));
    found
        .iter()
        .map(|f| text(field(f, "id")).to_owned())
        .collect()
}

// A well-formed finding id that names no finding of the demo, and nonsense.
const GUESSED: &str = "two-hits-cftr";
const NONSENSE: &str = "zz-not-a-finding";

/// A run's finding pages for `ids`: the status, and the body with the id it
/// was asked for replaced (the page names it, as its address does).
fn finding_answers(base: &str, run: &str, ids: &[&str]) -> Vec<(u16, String)> {
    ids.iter()
        .map(|id| {
            let (status, body) = http_get(base, &format!("/runs/{run}/finding/{id}/"));
            (status, String::from_utf8_lossy(&body).replace(id, "ID"))
        })
        .collect()
}

/// The same page, 200, for every id.
fn same_page(answers: &[(u16, String)]) -> bool {
    answers.iter().all(|a| a.0 == 200 && a.1 == answers[0].1)
}

fn statuses(answers: &[(u16, String)]) -> Vec<u16> {
    answers.iter().map(|a| a.0).collect()
}

/// The finding route of `og live serve` (test_owngenome_live_family.py
/// finding_route_checks, on main 79b5382c): a child with a result to act on
/// gets the same page for its own act, a guessed act and nonsense, so a
/// guessed address never tells whether the child has that result; an adult
/// keeps the 404 for an id that is not its finding.
fn finding_route_checks(base: &str) -> String {
    let child = finding_ids(base, CHILD);
    let adult = finding_ids(base, LIVE_ADULT);
    assert!(
        !child.is_empty() && !adult.is_empty(),
        "the demo's child and adult have findings"
    );
    for id in [GUESSED, NONSENSE] {
        assert!(
            !child.iter().chain(&adult).any(|x| x == id),
            "{id} must be no finding of the demo"
        );
    }
    let act = child[0].clone();
    let answers = finding_answers(base, CHILD, &[act.as_str(), GUESSED, NONSENSE]);
    assert!(
        same_page(&answers),
        "a child's finding route tells its act, a guessed act and nonsense apart: {:?}",
        statuses(&answers)
    );
    pass("finding route: a child's act, a guessed act and nonsense answer the same page");
    let answers = finding_answers(base, LIVE_ADULT, &[adult[0].as_str(), GUESSED, NONSENSE]);
    assert_eq!(
        statuses(&answers),
        [200, 404, 404],
        "an adult's finding route"
    );
    pass("finding route: an adult's own finding answers its page, any other id 404");
    act
}

/// A family folder's table with every row that names DEMO-003 changed: each
/// cell that names no person grows by 37 when it is a number and gets a
/// suffix when it is text. How many rows changed.
fn change_tsv(path: &Path) -> usize {
    let mut rows = 0;
    let mut out = String::new();
    for (i, line) in fs::read_to_string(path).unwrap().lines().enumerate() {
        let cells: Vec<&str> = line.split('\t').collect();
        if i == 0 || !cells.iter().any(|c| c.split(',').any(|m| m == CHILD)) {
            out.push_str(line);
        } else {
            rows += 1;
            let changed: Vec<String> = cells
                .iter()
                .map(|c| {
                    if c.split(',').all(|m| m.starts_with("DEMO-")) {
                        (*c).to_owned()
                    } else if let Ok(n) = c.parse::<i64>() {
                        (n + 37).to_string()
                    } else if let Ok(x) = c.parse::<f64>() {
                        (x + 37.0).to_string()
                    } else {
                        format!("{c} changed")
                    }
                })
                .collect();
            out.push_str(&changed.join("\t"));
        }
        out.push('\n');
    }
    fs::write(path, out).unwrap();
    rows
}

/// The family page on `og live serve` while DEMO-003 is closed: the child's
/// side-by-side cells read "closed", and its summary.json, which the side by
/// side reads, and every row of the family folder that names it (the
/// relationship, shared-segment, trio and plan rows, its inherited
/// stretches) leave every chapter of the page as it was.
fn live_family_payload(f: &mut Fixture, base: &str, core: &Path, source: &str) {
    let child_summary = f.work.join(CHILD).join("summary.json");
    let saved = fs::read(&child_summary).unwrap();
    let tables: Vec<PathBuf> = [
        "relationships.tsv",
        "ibd_pairs.tsv",
        "family_plan.tsv",
        "de_novo.tsv",
        "compound_het.tsv",
        format!("{CHILD}_segments.tsv").as_str(),
    ]
    .iter()
    .map(|name| f.work.join(FAMILY).join(name))
    .collect();
    let saved_tables: Vec<Vec<u8>> = tables.iter().map(|t| fs::read(t).unwrap()).collect();
    let b = &mut f.browser;
    let (cells, page) = closed_family(b, base);
    assert!(
        all_closed(&cells),
        "live family page: the closed child's side-by-side cells: {cells:?}"
    );
    let mut changed = doc(&child_summary);
    let sections = obj(obj(&mut changed).get_mut("sections").unwrap());
    let blood = obj(obj(sections.get_mut("bloodtype").unwrap())
        .get_mut("data")
        .unwrap());
    for system in array(blood.get_mut("systems").unwrap()) {
        obj(system).insert("phenotype".into(), s(CHANGED));
    }
    save(&child_summary, &changed);
    for table in &tables {
        assert!(change_tsv(table) > 0, "{table:?} has no row of {CHILD}");
    }
    let (changed_cells, changed_page) = closed_family(b, base);
    assert!(
        all_closed(&changed_cells) && changed_page == page,
        "live family page: a closed child's data changed the page: {changed_cells:?}"
    );
    pass(
        "live family page: a closed child's side-by-side cells read closed, and its blood groups and every family row naming it leave every chapter as it was",
    );
    fs::write(core, replace_once(source, FREE_CLOSED, FREE_LEAKED)).unwrap();
    let (_, leaked_page) = closed_family(b, base);
    fs::write(core, source).unwrap();
    assert!(
        leaked_page != page,
        "the live relationship-row control went unnoticed"
    );
    red(
        "live: core.js that keeps a closed child's relationship and shared-segment rows on the family page",
        &"the page changed",
    );
    for (table, raw) in tables.iter().zip(saved_tables) {
        fs::write(table, raw).unwrap();
        touch(table, None);
    }
    fs::write(core, replace_once(source, SIDE_CLOSED, SIDE_LEAKED)).unwrap();
    let (leaked, _) = closed_family(b, base);
    fs::write(core, source).unwrap();
    assert!(
        items(&leaked).contains(&s(CHANGED)),
        "the live side-by-side control went unnoticed: {leaked:?}"
    );
    red(
        "live: core.js that shows a closed child's side-by-side values on the family page",
        &leaked,
    );
    fs::write(&child_summary, saved).unwrap();
    touch(&child_summary, None);
}

fn copy_tree(from: &Path, to: &Path, rename: impl Fn(&str) -> String + Copy) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let dest = to.join(rename(entry.file_name().to_str().unwrap()));
        if entry.file_type().unwrap().is_dir() {
            copy_tree(&entry.path(), &dest, rename);
        } else {
            fs::copy(entry.path(), &dest).unwrap();
        }
    }
}

// A person's page as text, without the timings that change on every load.
const PAGE_JS: &str = r#"(function () {
  var clone = document.documentElement.cloneNode(true);
  clone.querySelectorAll('[data-pgp-map-paint-ms],[data-pgp-map-parse-ms]').forEach(function (n) {
    n.removeAttribute('data-pgp-map-paint-ms'); n.removeAttribute('data-pgp-map-parse-ms'); });
  clone.querySelectorAll('.ops-disk, #timing .v, #timing .r').forEach(function (n) { n.textContent = ''; });
  return clone.outerHTML;
}())"#;
const MAP_SETTLED: &str = "/^(drawn|failed|nodecompress|unavailable)$/.test((document.getElementById('genome-map') || {dataset: {}}).dataset.pgpMap || '')";

// A copied answer read as JSON (Response.clone().json()) resolves 1 s late.
const SLOW_CLONE: &str = r#"(function () {
  var clone = Response.prototype.clone;
  Response.prototype.clone = function () {
    var copy = clone.call(this), json = copy.json.bind(copy);
    copy.json = function () { return json().then(function (v) { return new Promise(function (done) { setTimeout(function () { done(v); }, 1000); }); }); };
    return copy;
  };
}())"#;

/// A page script that holds back the answer of every fetch whose path ends
/// with a key of `held` that many milliseconds: fixes which of two answers
/// the page reads first.
fn hold(held: &[(&str, u32)]) -> String {
    let list: Vec<String> = held
        .iter()
        .map(|(k, ms)| format!("[{}, {ms}]", js_str(k)))
        .collect();
    format!(
        r#"(function () {{
  var f = window.fetch, held = [{}];
  window.fetch = function (u) {{
    var p = new URL(u && u.url || u, location.href).pathname, ms = 0;
    held.forEach(function (h) {{ if (p.slice(-h[0].length) === h[0]) {{ ms += h[1]; }} }});
    var got = f.apply(this, arguments);
    return ms ? got.then(function (r) {{ return new Promise(function (done) {{ setTimeout(function () {{ done(r); }}, ms); }}); }}) : got;
  }};
}}());"#,
        list.join(", ")
    )
}

// The fixture summaries, as og live results reads them.
const PLP_BRCA1: &str =
    r#"{"gene": "BRCA1", "variant": "c.68_69del", "final": "Pathogenic", "clinvar": "Pathogenic"}"#;
const SV_DEL: &str = r#"{"state": "ok", "data": {"count": 1, "consensus": 1, "single_caller": 0, "exon_edge_checked": false, "exon_edge": 0, "common": 0, "rows": [{"sv": "chr17:43044295-43125483", "type": "DEL", "gt": "0/1", "length": -81188, "acmg_class": 5, "genes": ["BRCA1"]}]}}"#;
const STALE_SECTIONS: &str = r#""bloodtype": {"state": "stale", "note": "older than its input", "data": {"abo": "A1", "rhd": "D+", "systems": [{"system": "ABO", "phenotype": "A1"}, {"system": "RhD", "phenotype": "D+"}, {"system": "Kell", "phenotype": "K-k+"}, {"system": "Duffy", "phenotype": "Fy(a+b-)"}]}}, "haplogroup": {"state": "ok", "data": {"haplogroup": "H1"}}"#;

fn two_in_one_card() -> String {
    format!(
        r#"{{"sections": {{"cpsr": {{"state": "ok", "data": {{"plp": [{PLP_BRCA1}, {{"gene": "TP53", "variant": "c.743G>A", "final": "Pathogenic", "clinvar": "Pathogenic"}}]}}}}}}}}"#
    )
}

fn two_cards_sections() -> String {
    format!(r#""cpsr": {{"state": "ok", "data": {{"plp": [{PLP_BRCA1}]}}}}, "annotsv": {SV_DEL}"#)
}

fn two_cards() -> String {
    format!(r#"{{"sections": {{{}}}}}"#, two_cards_sections())
}

/// Write a run's summary.json. The family view reads it as it is; the run's
/// own results.json reads its step files, as the demo's steps have no end
/// time (`og live serve` reads summary.json only when it is newer than the
/// end of every done step).
fn summary(work: &Path, run: &str, body: &str) {
    fs::write(work.join(run).join("summary.json"), body).unwrap();
}

/// Take the child's two results to act on (KCNH2, from its ClinVar hits and
/// from AnnotSV) out of its step files, and touch its summary.json so the
/// server reads the run again. The files as they were.
fn without_acts(work: &Path) -> Vec<(PathBuf, Vec<u8>)> {
    let d = work.join(CHILD);
    let files = [
        d.join(format!("clinvar/{CHILD}_clinvar_hits.vcf")),
        d.join(format!("annotsv/{CHILD}_sv_class4_5_checked.tsv")),
    ];
    let saved = files
        .iter()
        .map(|p| (p.clone(), fs::read(p).unwrap()))
        .collect();
    for p in &files {
        let all = fs::read_to_string(p).unwrap();
        let kept: String = all
            .lines()
            .filter(|l| l.starts_with('#') || !l.contains("KCNH2"))
            .map(|l| format!("{l}\n"))
            .collect();
        assert!(kept.len() < all.len(), "no KCNH2 row in {p:?}");
        fs::write(p, kept).unwrap();
    }
    touch(&d.join("summary.json"), None);
    saved
}

fn put_back(work: &Path, saved: Vec<(PathBuf, Vec<u8>)>) {
    for (p, raw) in saved {
        fs::write(p, raw).unwrap();
    }
    touch(&work.join(CHILD).join("summary.json"), None);
}

/// Set a file's time to now, or to `stamp` (touch -t).
fn touch(path: &Path, stamp: Option<&str>) {
    let mut cmd = Command::new("touch");
    if let Some(stamp) = stamp {
        cmd.args(["-t", stamp]);
    }
    assert!(cmd.arg(path).status().unwrap().success(), "touch {path:?}");
}

/// A live page from about:blank once the family and the map have settled:
/// the page, its address, its reads in order and the map's state.
fn capture_live(
    b: &mut Browser,
    url: &str,
    core: &Path,
    core_text: &str,
    held: &[(&str, u32)],
    slow_clone: bool,
) -> Vec<Value> {
    fs::write(core, core_text).unwrap();
    let mut ids = Vec::new();
    if slow_clone {
        ids.push(install(b, SLOW_CLONE));
    }
    ids.push(install(b, &hold(held)));
    ids.push(install(b, REQUEST_LOG));
    go(b, url);
    wait(
        b,
        &format!("({READY}) && !!document.querySelector('#rail a[href*={FAMILY}]')"),
    );
    assert!(
        b.wait(MAP_SETTLED, Duration::from_secs(15)).is_some(),
        "the map never settled: {url}"
    );
    thread::sleep(Duration::from_secs(1));
    let got = eval(
        b,
        &format!(
            "[{PAGE_JS}, location.href, (window.ogRequests || []).map(function (r) {{ return r.join(' '); }}), (document.getElementById('genome-map') || {{dataset: {{}}}}).dataset.pgpMap || '']"
        ),
    );
    for id in ids {
        uninstall(b, &id);
    }
    items(&got).to_vec()
}

/// Up to four places where two pages differ, a little of each side.
fn differences(x: &str, y: &str) -> Vec<(String, String)> {
    x.split("><")
        .zip(y.split("><"))
        .filter(|(a, b)| a != b)
        .take(4)
        .map(|(a, b)| (a.chars().take(160).collect(), b.chars().take(160).collect()))
        .collect()
}

/// test_owngenome_live_family.py's live browser checks: in live mode
/// (state.json, then the event stream) an adult's and a no-act child's page
/// make the requests of the page without the gate, in its order; the marks a
/// section address sets follow main's outcome in each answer order; a
/// parent's page is the same whichever card holds a child's results.
fn live_mode_checks(f: &mut Fixture, base: &str, core: &Path, source: &str) {
    let child_summary = f.work.join(CHILD).join("summary.json");
    let saved = fs::read(&child_summary).unwrap();
    let steps = without_acts(&f.work);
    let b = &mut f.browser;
    let bare = bare_core(source);
    let map_last = [("/map.json", 500)];
    let mut adult_bare = Vec::new();
    for (label, run) in [("an adult", LIVE_ADULT), ("a no-act child", CHILD)] {
        let url = run_url(base, run);
        let pair = [
            capture_live(b, &url, core, &bare, &map_last, false),
            capture_live(b, &url, core, source, &map_last, false),
        ];
        assert!(
            pair[0][..3] == pair[1][..3],
            "live mode: {label}'s page differs with the gate: {:?} / {:?}; {:?}",
            pair[0][2],
            pair[1][2],
            differences(text(&pair[0][0]), text(&pair[1][0]))
        );
        pass(&format!(
            "live mode: the page of {label} makes the requests of the page without the gate, in its order, and is the same page"
        ));
        if run == LIVE_ADULT {
            assert_eq!(text(&pair[1][3]), "drawn", "the adult's genome map");
            pass("live mode: the adult's page compared has its genome map drawn");
            adult_bare = items(&pair[0][2]).to_vec();
        }
    }
    let first: Vec<&str> = adult_bare
        .iter()
        .take(3)
        .map(|r| {
            text(r)
                .split(' ')
                .nth(1)
                .unwrap()
                .rsplit('/')
                .next()
                .unwrap()
        })
        .collect();
    assert_eq!(
        first,
        ["state.json", "events", "family.json"],
        "live mode: the first reads without the gate: {adult_bare:?}"
    );
    pass(
        "live mode: without the gate the page asks for state.json, then the event stream, then family.json",
    );

    // A section address on an adult's and a no-act child's page: the marks it
    // sets depend on whether the family or the results answer first. Each
    // order is forced by holding the other answer back, twice each, and the
    // map's document comes last. The outcomes are the ones main's own view
    // files give: the marks stay when the family answers first and the
    // family's redraw clears them when the results answer first, and the
    // map's tooltip stays hidden="hidden". A copied answer read as JSON
    // resolves late, so a gate that checked a child's results on such a copy
    // would give the family-first outcome in the results-first order.
    let tip = regex::Regex::new(r#"<div class="gm-tip"[^>]*>"#).unwrap();
    let target = regex::Regex::new(r#"<article class="[^"]*\btarget\b"#).unwrap();
    let mut marks = Vec::new();
    for (first, held) in [
        ("results.json", "/family.json"),
        ("family.json", "/results.json"),
    ] {
        for who in [LIVE_ADULT, CHILD] {
            for place in ["#s3", "#s/clinvar"] {
                for _ in 0..2 {
                    let url = format!("{}{place}", run_url(base, who));
                    let page = capture_live(
                        b,
                        &url,
                        core,
                        source,
                        &[(held, 500), ("/map.json", 1000)],
                        true,
                    );
                    let page = text(&page[0]).to_owned();
                    let got = (
                        page.matches(r#"aria-current="true""#).count(),
                        target.is_match(&page),
                        tip.find(&page)
                            .is_some_and(|t| t.as_str().contains(r#"hidden="hidden""#)),
                    );
                    let want = if first == "results.json" {
                        (0, false, true)
                    } else {
                        (1, place == "#s/clinvar", true)
                    };
                    marks.push((first, who, place, got, got == want));
                }
            }
        }
    }
    let bad: Vec<_> = marks.iter().filter(|m| !m.4).collect();
    assert!(
        bad.is_empty(),
        "live mode: section marks in each order: {bad:?}"
    );
    pass(
        "live mode: an adult's and a no-act child's page at a section address give main's outcome in each order",
    );

    put_back(&f.work, steps);

    // A parent's page, whichever card holds the child's two results to act on.
    // The two summaries must give different flagged-card counts, or the
    // comparison could not see a page that tells them apart.
    let attention = |body: &str| {
        int(field(
            &og_results(&f.og, body, &["--overview"]),
            "attention",
        ))
    };
    let counts = (attention(&two_in_one_card()), attention(&two_cards()));
    assert_ne!(
        counts.0, counts.1,
        "live mode: the two summaries give the same flagged-card count"
    );
    let mut seen = Vec::new();
    for body in [two_in_one_card(), two_cards()] {
        summary(&f.work, CHILD, &body);
        seen.push(
            capture_live(
                b,
                &run_url(base, LIVE_ADULT),
                core,
                source,
                &map_last,
                false,
            )[0]
            .clone(),
        );
    }
    assert!(
        seen[0] == seen[1],
        "live mode: a parent's page differs by the card that holds the child's results: {:?}",
        differences(text(&seen[0]), text(&seen[1]))
    );
    pass("live mode: a parent's page is the same whichever card holds a child's results to act on");
    fs::write(&child_summary, saved).unwrap();
    touch(&child_summary, None);
    fs::write(core, source).unwrap();
}

const MEMBERS_HEADER: &str = "sample\tsex\tfather\tmother\tpartner\tin_pedigree\tsequenced\n";
const FAMILY_CARD: &str = "[...document.querySelectorAll('.person-card')].find(c => c.querySelector('.person-name').textContent === 'DEMO-003')";
const PARTNER_GATE: &str = "document.querySelector('[data-gate-child=DEMO-003]')";

fn partner_view(b: &mut Browser, url: &str) -> Value {
    go(b, url);
    wait(
        b,
        &format!("{READY} && !!document.querySelector('.person-card')"),
    );
    eval(
        b,
        &format!(
            r#"[
        {FAMILY_CARD}?.querySelector('.person-look')?.textContent,
        {PARTNER_GATE}?.parentElement.closest('section')?.id,
        document.querySelector('#ch-children')?.textContent.includes('The pedigree names no child.'),
        document.querySelector('#ch-couples')?.textContent.includes('DEMO-003 and DEMO-005')
    ]"#
        ),
    )
}

/// The member line of `who` on a family page, and its pedigree flag.
fn member_line(b: &mut Browser, who: &str) -> (Value, Value) {
    let line = eval(
        b,
        &format!(
            "[...document.querySelectorAll('.person-card')].filter(c => c.querySelector('.person-name').textContent === '{who}').map(c => {{ const l = c.querySelector('.person-look'); return l ? l.textContent : ''; }})[0]"
        ),
    );
    let flag = eval(
        b,
        &format!(
            "[...document.querySelectorAll('.ped-person')].filter(g => g.querySelector('.ped-name').textContent === '{who}').map(g => {{ const f = g.querySelector('.ped-flag'); return f ? f.textContent : ''; }})[0]"
        ),
    );
    (line, flag)
}

fn family_gate_checks(f: &mut Fixture, base: &str, core: &Path, source: &str) {
    // The child's summary.json is older than the results the page reads: its
    // count is the one family.json gives without the gate.
    let child_summary = f.work.join(CHILD).join("summary.json");
    let adult_summary = f.work.join(LIVE_ADULT).join("summary.json");
    let (child_saved, adult_saved) = (
        fs::read_to_string(&child_summary).unwrap(),
        fs::read(&adult_summary).unwrap(),
    );
    let stale = format!(
        r#"{{"sections": {{{}, {STALE_SECTIONS}}}}}"#,
        two_cards_sections()
    );
    fs::write(&child_summary, &stale).unwrap();
    touch(&child_summary, Some("200001010000"));
    let attention = |body: &str| {
        int(field(
            &og_results(&f.og, body, &["--overview"]),
            "attention",
        ))
    };
    let (n, now) = (attention(&stale), attention(&child_saved));
    assert!(n > 0 && n != now, "the stale count must differ: {n} {now}");
    fs::write(&adult_summary, "{not json").unwrap();

    let b = &mut f.browser;
    go(b, &run_url(base, FAMILY));
    let closed_controls =
        "document.querySelectorAll('[data-gate-child=DEMO-003] [data-opening^=closed]')";
    wait(b, &format!("{closed_controls}.length"));
    let closed = member_line(b, CHILD);
    let adult = member_line(b, LIVE_ADULT).0;
    for _ in 0..12 {
        let left = int(&eval(b, &format!("{closed_controls}.length")));
        if left == 0 {
            break;
        }
        eval(b, &format!("{closed_controls}[0].click(); true"));
        wait(
            b,
            &format!("!{PARTNER_GATE} || {closed_controls}.length < {left}"),
        );
    }
    wait(
        b,
        &format!(
            "!{PARTNER_GATE} && /worth a look|Nothing flagged/.test({FAMILY_CARD}?.querySelector('.person-look')?.textContent || '')"
        ),
    );
    let opened = member_line(b, CHILD);
    assert!(
        closed == (s("Report closed"), s(""))
            && opened == (s(&format!("{n} worth a look")), s(&format!("▲ {n}"))),
        "family page lines: closed {closed:?}, opened {opened:?}, summary.json {n}, results {now}"
    );
    pass(
        "family page: a gated child's member line is the fixed closed line with no flag; once opened, summary.json's count (not the newer results') shows on the line and the pedigree",
    );
    assert_eq!(
        text(&adult),
        "Nothing flagged; see the caveats",
        "an adult whose summary.json cannot be read"
    );
    pass("family page: an adult whose summary.json cannot be read reads as without the gate");
    fs::write(&child_summary, child_saved).unwrap();
    fs::write(&adult_summary, adult_saved).unwrap();

    // A person who is a child in another family is closed on this family's
    // page too, keeps the partner role and the couple's line, whichever
    // family sorts first (a child in any family is one).
    let expected = loads(r#"["Report closed","ch-members",true,true]"#).unwrap();
    let planted = replace_once(
        source,
        "return childRun(f, m.sample);",
        "return childRun(f, m.sample) && (f.state.children || []).some(c => c.child === m.sample);",
    );
    for partners in [PARTNERS, PARTNERS_LAST] {
        let dir = f.work.join(partners);
        fs::create_dir(&dir).unwrap();
        fs::write(dir.join("members.tsv"), format!(
            "{MEMBERS_HEADER}DEMO-003\tmale\t\t\tDEMO-005\tyes\tyes\nDEMO-005\tfemale\t\t\tDEMO-003\tyes\tno\n"
        )).unwrap();
        let url = run_url(base, partners);
        assert_eq!(partner_view(&mut f.browser, &url), expected, "{partners}");
        pass(&format!(
            "family page {partners}: a child in another family is closed here too, keeps the partner role and the couple's line"
        ));
        fs::write(core, &planted).unwrap();
        let got = partner_view(&mut f.browser, &url);
        assert_ne!(got, expected);
        red(
            &format!("{partners}: deciding from this pedigree alone"),
            &got,
        );
        fs::write(core, source).unwrap();
        if partners == PARTNERS_LAST {
            // Who is a child is the server's decision (family.json marks each
            // member run), over every family. A server where the family that
            // sorts last decides would mark DEMO-003 no child here, as this
            // family names it a partner; on DEMO-A-PARTNERS it would still
            // answer right. Its answer, planted in the page:
            let last = install(
                &mut f.browser,
                &format!(
                    r#"(function () {{
  var f = window.fetch;
  window.fetch = function (u) {{
    var p = new URL(u && u.url || u, location.href).pathname;
    var got = f.apply(this, arguments);
    if (p !== '/runs/{PARTNERS_LAST}/family.json') {{ return got; }}
    return got.then(function (r) {{ return r.json().then(function (d) {{
      if (d.runs && d.runs['{CHILD}']) {{ d.runs['{CHILD}'].child = false; }}
      return new Response(JSON.stringify(d), {{ status: r.status, headers: r.headers }});
    }}); }});
  }};
}}());"#
                ),
            );
            let got = partner_view(&mut f.browser, &url);
            uninstall(&mut f.browser, &last);
            assert_ne!(
                got, expected,
                "the last-family-decides control went unnoticed"
            );
            red(
                &format!(
                    "{partners}: a server where the family that sorts last decides who is a child"
                ),
                &got,
            );
        }
        fs::remove_dir_all(dir).unwrap();
    }
}

fn run_state(b: &mut Browser, base: &str, name: &str) -> Value {
    go(b, &run_url(base, name));
    wait(
        b,
        &format!("({CLOSED}) ? 'closed' : ({OPEN}) ? 'open' : ''"),
    )
}

fn home_states(b: &mut Browser, base: &str) -> Value {
    go(b, &format!("{base}/"));
    wait(
        b,
        r#"(() => {
        const tiles = ['DEMO-003', 'DEMO-004'].map(n => document.querySelector('.ptile[data-sample=' + n + ']'));
        return tiles.every(Boolean) && tiles.map(t => t.textContent.includes(PGPUI.closedText) ? 'closed' : 'open');
    })()"#,
    )
}

fn failure_browser_checks(f: &mut Fixture, base: &str, act: &str) {
    let family = f.work.join(FAMILY);
    // The result remains actionable, but the pedigree now lists DEMO-003 as an adult.
    fs::write(family.join("members.tsv"), format!(
        "{MEMBERS_HEADER}DEMO-003\tmale\t\t\tDEMO-002\tyes\tyes\nDEMO-002\tfemale\t\t\tDEMO-003\tyes\tyes\n"
    )).unwrap();
    copy_tree(&f.work.join(CHILD), &f.work.join(OTHER), |name| {
        name.replace(CHILD, OTHER)
    });
    let relationships = family.join("relationships.tsv");
    let mode = fs::metadata(&relationships).unwrap().permissions();
    fs::set_permissions(&relationships, fs::Permissions::from_mode(0o000)).unwrap();
    assert!(
        fs::File::open(&relationships).is_err(),
        "chmod must make the fixture unreadable"
    );

    let closed_home = loads(r#"["closed","open"]"#).unwrap();
    assert_eq!(text(&run_state(&mut f.browser, base, CHILD)), "closed");
    assert_eq!(text(&run_state(&mut f.browser, base, OTHER)), "open");
    assert_eq!(home_states(&mut f.browser, base), closed_home);
    pass(
        "a family file that cannot be read keeps the runs it lists closed, page and tile, and no other run",
    );
    let ids = [act, GUESSED, NONSENSE];
    let listed = finding_answers(base, CHILD, &ids);
    let other = finding_answers(base, OTHER, &ids);
    assert!(
        same_page(&listed) && statuses(&other) == [200, 404, 404],
        "finding route while a family file cannot be read: {CHILD} {:?}, {OTHER} {:?}",
        statuses(&listed),
        statuses(&other)
    );
    pass(
        "finding route while a family file cannot be read: a run it lists answers the same page for any id, another run 404 for a guessed one",
    );

    fs::set_permissions(relationships, mode).unwrap();
    assert_eq!(text(&run_state(&mut f.browser, base, CHILD)), "open");
    assert_eq!(
        home_states(&mut f.browser, base),
        loads(r#"["open","open"]"#).unwrap()
    );
    pass("the child's page and tile reopened once the family file reads");
    // members.tsv now names DEMO-003 an adult: the same-page check must fail
    // for it, with a guessed id and nonsense answering 404.
    let adult = finding_answers(base, CHILD, &ids);
    assert!(
        !same_page(&adult) && statuses(&adult) == [200, 404, 404],
        "the finding route control went unnoticed: {:?}",
        statuses(&adult)
    );
    red(
        "the finding route of DEMO-003 made an adult in members.tsv",
        &statuses(&adult),
    );

    // While the family's members cannot be read, every run is closed and its
    // family is unknown, never "in no family" (test_owngenome_live_family.py,
    // on main 79b5382c). Before, every run is an adult: open, and known.
    let runs = person_runs(base);
    assert!(
        [CHILD, OTHER, LIVE_ADULT]
            .iter()
            .all(|r| runs.iter().any(|x| x == r)),
        "the live runs: {runs:?}"
    );
    // The runs with a result to act on, fixed while the family reads, so
    // none can drop out of the finding-route check later.
    let acts = act_runs(base, &runs);
    assert!(
        acts.iter().any(|r| r == CHILD),
        "the runs with an act: {acts:?}"
    );
    let open = member_answers(base, &runs, &acts);
    assert!(
        !all_unknown(&open)
            && open
                .iter()
                .all(|a| a.1 == Value::Null && a.2.is_none_or(|n| n == 404)),
        "every run is an adult while the family reads: {open:?}"
    );
    let members = family.join("members.tsv");
    let mode = fs::metadata(&members).unwrap().permissions();
    fs::set_permissions(&members, fs::Permissions::from_mode(0o000)).unwrap();
    assert!(
        fs::File::open(&members).is_err(),
        "chmod must make the fixture unreadable"
    );
    let unread = member_answers(base, &runs, &acts);
    let pages = [
        run_state(&mut f.browser, base, CHILD),
        run_state(&mut f.browser, base, OTHER),
    ];
    let home = home_states(&mut f.browser, base);
    fs::set_permissions(&members, mode).unwrap();
    assert!(
        all_unknown(&unread),
        "while members.tsv cannot be read, a run is open or known: {unread:?}"
    );
    assert_eq!(pages, [s("closed"), s("closed")], "closed pages");
    assert_eq!(
        home,
        loads(r#"["closed","closed"]"#).unwrap(),
        "closed tiles"
    );
    pass(
        "while a family's members cannot be read, every run is closed (state.json, the finding route, page and tile) and its family is unknown",
    );
    let back = member_answers(base, &runs, &acts);
    assert_eq!(back, open, "the runs once members.tsv reads again");
    assert_eq!(text(&run_state(&mut f.browser, base, OTHER)), "open");
    pass("every run reopened once members.tsv reads");
    red(
        "the same check on the readable members.tsv, where every run is an adult",
        &open,
    );
}

/// The person runs `og live serve` lists in /runs.json.
fn person_runs(base: &str) -> Vec<String> {
    let (status, body) = http_get(base, "/runs.json");
    assert_eq!(status, 200, "/runs.json");
    let index = loads(&String::from_utf8(body).unwrap()).unwrap();
    items(field(&index, "runs"))
        .iter()
        .filter(|r| *field(r, "kind") == s("sample"))
        .map(|r| text(field(r, "name")).to_owned())
        .collect()
}

/// The runs whose results.json holds a result to act on.
fn act_runs(base: &str, runs: &[String]) -> Vec<String> {
    runs.iter()
        .filter(|run| {
            let (status, body) = http_get(base, &format!("/runs/{run}/results.json"));
            // A run whose results are a note alone has no findings.
            status == 200
                && matches!(
                    loads(&String::from_utf8(body).unwrap()).unwrap(),
                    Value::Object(o) if matches!(o.get("findings"), Some(Value::Array(f)) if f.iter().any(is_act))
                )
        })
        .cloned()
        .collect()
}

/// Each run's family plugin in its state.json (whether its family is
/// unknown) and, for a run in `acts`, what a guessed finding id answers:
/// 200 while the run is closed, 404 while it is an adult's. A run with no
/// result to act on answers 404 either way: it has nothing to hide.
fn member_answers(
    base: &str,
    runs: &[String],
    acts: &[String],
) -> Vec<(String, Value, Option<u16>)> {
    runs.iter()
        .map(|run| {
            let (status, body) = http_get(base, &format!("/runs/{run}/state.json"));
            assert_eq!(status, 200, "{run}'s state.json");
            let st = loads(&String::from_utf8(body).unwrap()).unwrap();
            let family = field(field(&st, "plugins"), "family");
            let unknown = match family {
                Value::Object(o) => o.get("unknown").cloned().unwrap_or(Value::Null),
                _ => Value::Null,
            };
            let guessed = acts
                .contains(run)
                .then(|| http_get(base, &format!("/runs/{run}/finding/{GUESSED}/")).0);
            (run.clone(), unknown, guessed)
        })
        .collect()
}

/// Every run's family unknown, and every run with a result to act on
/// closed on the finding route.
fn all_unknown(answers: &[(String, Value, Option<u16>)]) -> bool {
    answers.iter().any(|a| a.2.is_some())
        && answers
            .iter()
            .all(|a| a.1 == Value::Bool(true) && a.2.is_none_or(|n| n == 200))
}

const LIVE_RESULTS: &str = r#"[
    document.querySelector('#results').innerHTML,
    Boolean(PGP.searchHidden(null)),
    [...document.querySelectorAll('.gm-chip')].map(node => node.textContent),
    document.querySelectorAll('[data-opening]').length
]"#;

fn adult_results(b: &mut Browser, base: &str, expected_family: bool) -> Value {
    go(b, &run_url(base, CHILD));
    wait(
        b,
        &format!(
            r#"({OPEN}) && document.querySelector('[data-pgp-map=drawn]')
        && Boolean(document.querySelector('#rail a[href*=DEMO-FAMILY]')) === {expected_family}
        && (!window.familyFault || document.documentElement.dataset.familyFault === window.familyFault)"#
        ),
    );
    eval(b, LIVE_RESULTS)
}

// Answer one request with a 502, or read it and close the socket without a reply.
fn family_fault(listener: TcpListener, status: &str) {
    listener.set_nonblocking(true).unwrap();
    let until = std::time::Instant::now() + WAIT;
    let mut connection = loop {
        match listener.accept() {
            Ok((connection, _)) => break connection,
            Err(e) if e.kind() == std::io::ErrorKind::WouldBlock => {
                assert!(
                    std::time::Instant::now() < until,
                    "family fault was never requested"
                );
                thread::sleep(Duration::from_millis(10));
            }
            Err(e) => panic!("cannot accept the family fault request: {e}"),
        }
    };
    connection.set_read_timeout(Some(WAIT)).unwrap();
    let mut reader = BufReader::new(&mut connection);
    let mut line = String::new();
    reader.read_line(&mut line).unwrap();
    assert!(line.starts_with("GET /family.json "), "{line}");
    while line != "\r\n" {
        line.clear();
        assert_ne!(reader.read_line(&mut line).unwrap(), 0);
    }
    if status == "502" {
        write!(connection, "HTTP/1.1 502 Bad Gateway\r\nAccess-Control-Allow-Origin: *\r\nContent-Length: 0\r\nConnection: close\r\n\r\n").unwrap();
    }
}

fn family_fetch_failure_checks(f: &mut Fixture, base: &str, person_js: &Path) {
    let b = &mut f.browser;
    let healthy = adult_results(b, base, true);
    let Value::Array(metrics) = &healthy else {
        panic!("missing adult report metrics")
    };
    assert_eq!(metrics[1], Value::Bool(false));
    assert_eq!(metrics[3], Value::Int(0));
    let Value::Array(chips) = &metrics[2] else {
        panic!("missing map chips")
    };
    assert!(chips.iter().all(|chip| !text(chip).contains("Opened")));
    let source = fs::read_to_string(person_js).unwrap();
    // Let the page reach the second loopback port. Only the family request is redirected.
    call(b, "Page.setBypassCSP", r#"{"enabled":true}"#);
    for (fault, planted) in [("502", false), ("drop", false), ("502", true)] {
        if planted {
            fs::write(
                person_js,
                replace_once(
                    &source,
                    "function closedPage() {",
                    "function closedPage() { if (familyStatus === 'failed') { return true; }",
                ),
            )
            .unwrap();
        }
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let address = listener.local_addr().unwrap();
        let redirect = install(
            b,
            &format!(
                r#"(() => {{
            window.familyFault = '{fault}';
            const original = window.fetch;
            window.fetch = function(url, opts) {{
                const path = new URL(url?.url || url, location.href).pathname;
                if (path !== '/runs/DEMO-FAMILY/family.json') return original.apply(this, arguments);
                return original('http://{address}/family.json', {{credentials:'omit'}}).then(response => {{
                    document.documentElement.dataset.familyFault = String(response.status);
                    return response;
                }}, error => {{
                    document.documentElement.dataset.familyFault = 'drop';
                    throw error;
                }});
            }};
        }})()"#
            ),
        );
        thread::scope(|scope| {
            let served = scope.spawn(|| family_fault(listener, fault));
            let failed = adult_results(b, base, false);
            served.join().unwrap();
            if planted {
                assert!(
                    failed != healthy,
                    "the faulty 502 handling changed no adult controls"
                );
                red(
                    "person.js that treats a failed family fetch as closed",
                    &"the adult report changed",
                );
            } else {
                assert!(failed == healthy, "{fault} changed the adult report");
                pass(&format!(
                    "family.json {fault}: the adult report stayed open"
                ));
            }
        });
        uninstall(b, &redirect);
    }
    call(b, "Page.setBypassCSP", r#"{"enabled":false}"#);

    fs::write(person_js, source).unwrap();
}

#[test]
fn live_serve_report_gate() {
    let Some(mut f) = Fixture::build("serve-gate") else {
        return;
    };
    let assets = f.scratch.path().join("assets/live");
    copy_tree(&repo_root().join("assets/live"), &assets, str::to_owned);
    std::os::unix::fs::symlink(
        repo_root().join("assets/reference"),
        assets.with_file_name("reference"),
    )
    .unwrap();
    let core = assets.join("core.js");
    let source = fs::read_to_string(&core).unwrap();
    let server = LiveServe::start(
        &f.og,
        &f.work,
        &f.scratch.path().join("server-state"),
        &assets,
    );
    let act = finding_route_checks(&server.base);
    live_mode_checks(&mut f, &server.base, &core, &source);
    live_family_payload(&mut f, &server.base, &core, &source);
    family_gate_checks(&mut f, &server.base, &core, &source);
    failure_browser_checks(&mut f, &server.base, &act);
    family_fetch_failure_checks(&mut f, &server.base, &assets.join("person.js"));
}

// child_a's summary.json for the side by side: a stale blood group with
// four systems, and a haplogroup that is not stale.
const STALE_SUMMARY: &str = r#"{"sections": {"bloodtype": {"state": "stale", "note": "older than its input", "data": {"abo": "A1", "rhd": "D+", "systems": [{"system": "ABO", "phenotype": "A1"}, {"system": "RhD", "phenotype": "D+"}, {"system": "Kell", "phenotype": "K-k+"}, {"system": "Duffy", "phenotype": "Fy(a+b-)"}]}}, "haplogroup": {"state": "ok", "data": {"haplogroup": "H1"}}}}"#;
const FAMILY_READY: &str = "document.documentElement.dataset.pgpStatus === 'ready'";
const COMPOUND_LINE: &str = "var p = document.querySelector('details[data-child=child_a] p[data-trio=compound_het]'); p && p.textContent";

/// child_a's compound-het line on a fresh load of family_a's page.
fn compound_line(b: &mut Browser, base: &str) -> String {
    go(b, &run_url(base, "family_a"));
    wait(
        b,
        &format!(
            "{FAMILY_READY} && !!document.querySelector('details[data-child=child_a] p[data-trio=compound_het]')"
        ),
    );
    words(&eval(b, COMPOUND_LINE))
}

/// The family page of tests/fixtures/families_live: the card, the pedigree,
/// a child whose results cannot be read, then the couple screen's risk
/// drawing, a child's trio files (one missing, one linking out of the
/// folder, then its compound-het line from no candidates, newer candidates
/// and older ones) and the side by side (a stale value, every blood group
/// system). Moved here from browser_checks and depth_checks of
/// tests/test_owngenome_live_family.py, which ran them on the Python site;
/// those had no planted mistake, so every check here names the exact words
/// or count it needs.
#[test]
fn families_live_page() {
    let Some(chrome) = chrome_or_skip("open the family page") else {
        return;
    };
    let Some(og) = og_binary() else { return };
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "families-live");
    let root = scratch.path().join("genome");
    copy_tree(
        &repo_root().join("tests/fixtures/families_live"),
        &root,
        str::to_owned,
    );
    let server = LiveServe::start(
        &og,
        &root,
        &scratch.path().join("state"),
        &repo_root().join("assets/live"),
    );
    let base = server.base.clone();
    let b = &mut Browser::launch(&chrome, &scratch.path().join("chrome")).unwrap();

    go(b, &run_url(&base, "family_a"));
    let card = wait(
        b,
        "var c = document.querySelector('section[data-plugin=family]'); c && !c.hidden && c.querySelector('img.pedigree') && c.textContent",
    );
    assert!(words(&card).contains("Couple carrier screen"), "{card:?}");
    pass("browser: the family card shows on family_a's page");
    wait(
        b,
        "var i = document.querySelector('img.pedigree'); i && i.complete && i.naturalWidth > 0",
    );
    pass("browser: the pedigree image loads");
    // child_a has a run whose results cannot be read: closed whole, its
    // chromosome image too. Its parents have no run here, so nothing of
    // theirs is closed: their couple screen shows.
    let got = wait(
        b,
        "!!document.querySelector('[data-gate-child=child_a]') && [!!document.querySelector('img.segments'), !!document.querySelector('#ch-couples table')].join()",
    );
    assert_eq!(words(&got), "false,true");
    pass(
        "browser: child_a's unreadable report stays closed whole; its parents' couple screen shows",
    );
    let blocked = eval(b, "document.documentElement.dataset.pgpBlocked || ''");
    assert_eq!(words(&blocked), "", "the policy blocked {blocked:?}");
    pass("browser: the policy blocked nothing on the family page");

    let fam = root.join("family_a");
    let tsv = fam.join("couple_father_a_mother_a.tsv");
    let mut text = fs::read_to_string(&tsv)
        .unwrap()
        .replace("\tshared\t25%\t", "\tshared\t1 in 4 per pregnancy\t");
    text.push_str("GENE4\tAR\tACMG\tinvented condition four\tshared\t1 in 2 per pregnancy\tfather_a\tinvented finding\tmother_a\tinvented finding\tinvented\n");
    fs::write(&tsv, text).unwrap();
    let secret = scratch.path().join("outside-folder-name.tsv");
    fs::write(
        &secret,
        "child\tgene\tvariant_a\tvariant_b\ta_from\tb_from\tphase\n",
    )
    .unwrap();
    std::os::unix::fs::symlink(&secret, fam.join("compound_het.tsv")).unwrap();
    fs::write(root.join("child_a/summary.json"), STALE_SUMMARY).unwrap();
    // Readable no-act results for child_a and its parents: results that
    // cannot be read never open.
    let log = root.join("child_a/logs/run_status.tsv");
    let old_log = fs::read_to_string(&log).unwrap();
    fs::write(&log, format!("{old_log}step\t06\tok\n")).unwrap();
    for parent in ["father_a", "mother_a"] {
        fs::create_dir_all(root.join(parent).join("logs")).unwrap();
        fs::write(
            root.join(parent).join("logs/run_status.tsv"),
            "step\t06\tok\n",
        )
        .unwrap();
    }

    go(b, &run_url(&base, "family_a"));
    wait(
        b,
        &format!("{FAMILY_READY} && !!document.querySelector('div[data-shared=GENE4]')"),
    );
    let four = eval(
        b,
        "var k = document.querySelector('div[data-shared=GENE1] .kids4 .grid4'); k ? k.children.length : -1",
    );
    assert_eq!(four, Value::Int(4));
    pass("browser: a 1-in-4 shared row draws the four children");
    let two = words(&eval(
        b,
        "var d = document.querySelector('div[data-shared=GENE4]'); d ? [!!d.querySelector('.kids4'), d.textContent].join('|') : ''",
    ));
    assert!(
        two.starts_with("false|") && two.contains("1 in 2 per pregnancy"),
        "{two}"
    );
    pass("browser: a 1-in-2 shared row gives the chance in words, without the drawing");
    let de_novo = words(&eval(
        b,
        "var p = document.querySelector('details[data-child=child_a] p[data-trio=de_novo]'); p && p.textContent",
    ));
    let compound = words(&eval(b, COMPOUND_LINE));
    let page = words(&eval(b, "document.body.textContent"));
    assert!(
        de_novo.contains("Not assessed")
            && de_novo.contains("de_novo.tsv is not here")
            && compound.contains("link out of the family folder")
            && !page.contains("outside-folder-name"),
        "{de_novo:?} {compound:?}"
    );
    pass(
        "browser: a missing trio file and one linking out of the folder are told apart, the link's target unshown",
    );
    let rows = words(&eval(
        b,
        "Array.from(document.querySelectorAll('#ch-side-by-side tbody tr')).map(function (r) { return Array.from(r.cells).map(function (c) { return c.textContent; }).join('|'); }).join('\\n')",
    ));
    let cell = |system: &str| {
        rows.lines()
            .find(|r| r.contains(&format!("|{system}|")))
            .and_then(|r| r.split('|').nth(2))
            .unwrap_or("")
            .to_owned()
    };
    let legend = eval(
        b,
        "var p = document.querySelector('#ch-side-by-side p[data-stale]'); p && p.dataset.stale",
    );
    assert!(
        rows.contains("|Duffy|")
            && cell("Kell") == "K-k+ (older than its input)"
            && cell("Mitochondrial") == "H1"
            && words(&legend) == "4",
        "{rows}\n{legend:?}"
    );
    pass(
        "browser: side by side has a row per blood group system, a stale value marked in its cell and no other",
    );

    fs::remove_file(fam.join("compound_het.tsv")).unwrap();
    fs::copy(&secret, fam.join("compound_het.tsv")).unwrap();
    let line = compound_line(b, &base);
    assert!(
        line.starts_with("Not assessed") && line.contains("wrote no candidates for child_a"),
        "{line}"
    );
    pass(
        "browser: a child whose run has no step 23 candidates is not assessed for compound-het pairs",
    );
    let clinical = root.join("child_a/clinical");
    fs::create_dir_all(&clinical).unwrap();
    let vcf = clinical.join("child_a_clinical.vcf.gz");
    fs::write(&vcf, "").unwrap();
    touch(&fam.join("compound_het.tsv"), Some("202601010000"));
    let line = compound_line(b, &base);
    assert!(
        line.starts_with("Not known whether assessed")
            && line.contains("newer than compound_het.tsv"),
        "{line}"
    );
    pass("browser: candidates newer than compound_het.tsv leave the compound-het line not known");
    touch(&vcf, Some("202501010000"));
    let line = compound_line(b, &base);
    assert!(
        line.starts_with("Assessed (") && line.contains("step 23"),
        "{line}"
    );
    pass(
        "browser: an assessed child's empty compound-het line says pairs come only from step 23's candidates",
    );

    // child_a's results unreadable again: its page shows the fixed gate.
    fs::write(&log, old_log).unwrap();
    fs::remove_file(root.join("child_a/summary.json")).unwrap();
    go(b, &run_url(&base, "child_a"));
    wait(b, CLOSED);
    pass("browser: an unreadable child's page shows the fixed gate");
    drop(server);
}

/// A page value as text: a string as it is, anything else (null, false) as
/// nothing.
fn words(value: &Value) -> String {
    match value {
        Value::Str(text) => text.clone(),
        _ => String::new(),
    }
}
