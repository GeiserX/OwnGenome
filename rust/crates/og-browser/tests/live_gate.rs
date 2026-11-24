//! Chrome coverage of report gate isolation, transparency, and failure paths.

use og_browser::{Browser, Scratch, StaticServer, Value, chrome_or_skip, js_str, repo_root};
use pycompat::json::{dumps, loads};
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::net::TcpListener;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::thread;
use std::time::Duration;

const ADULT: &str = "DEMO-001";
const CHILD: &str = "DEMO-003";
const OTHER: &str = "DEMO-004";
const FAMILY: &str = "DEMO-FAMILY";
const PARTNERS: &str = "DEMO-A-PARTNERS";
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
        let og = std::env::var_os("OG").map(PathBuf::from).or_else(|| {
            let exe = std::env::current_exe().ok()?;
            Some(exe.parent()?.parent()?.join("og"))
        })?;
        if !og.is_file() {
            assert_ne!(std::env::var("CI").as_deref(), Ok("true"), "build og in CI");
            eprintln!("[SKIP] build og or set OG");
            return None;
        }
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

fn call(b: &mut Browser, method: &str, params: &str) -> Value {
    b.call(method, loads(params).unwrap()).unwrap()
}

fn eval(b: &mut Browser, js: &str) -> Value {
    b.eval(js).unwrap_or_else(|e| panic!("{js}: {e}"))
}

fn wait(b: &mut Browser, js: &str) -> Value {
    b.wait(js, WAIT).unwrap_or_else(|| {
        let page = b.eval("[location.href, document.documentElement.dataset, document.body.textContent.slice(0,350)]");
        panic!("{js} did not settle: {page:?}; {:?}", b.page_errors());
    })
}

fn go(b: &mut Browser, url: &str) {
    b.go("about:blank").unwrap();
    wait(b, "location.href === 'about:blank'");
    b.go(url).unwrap();
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

fn replace_once(source: &str, from: &str, to: &str) -> String {
    assert_eq!(source.matches(from).count(), 1, "control anchor: {from}");
    source.replacen(from, to, 1)
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
    println!("red control: shared family opening changed DEMO-004");
    fs::write(core, source).unwrap();
}

const CAPTURE: &str = r#"(() => {
  const clone = document.documentElement.cloneNode(true);
  clone.querySelectorAll('[data-pgp-map-paint-ms],[data-pgp-map-parse-ms]').forEach(node => {
    node.removeAttribute('data-pgp-map-paint-ms');
    node.removeAttribute('data-pgp-map-parse-ms');
  });
  for (const node of [clone, ...clone.querySelectorAll('*')]) {
    const attrs = Array.from(node.attributes).sort((a, b) => a.name.localeCompare(b.name));
    for (const attr of attrs) node.removeAttributeNode(attr);
    for (const attr of attrs) node.setAttributeNodeNS(attr);
  }
  const resources = performance.getEntriesByType('resource').map(r => new URL(r.name).pathname)
    .filter(path => path !== '/favicon.ico').sort();
  const reads = ogRequests.map(r => JSON.stringify(r)).sort();
  const order = ogRequests.filter(r => /\/(frames|family|results)\.json$/.test(r[0]));
  return [clone.outerHTML, location.href, Array.from(document.querySelectorAll('canvas'), c => c.toDataURL()),
    resources, reads, order];
})()"#;

const REQUESTS: &str = r#"(() => {
  const original = window.fetch;
  window.ogRequests = [];
  window.fetch = function(url, options) {
    ogRequests.push([new URL(url?.url || url, location.href).pathname,
      options?.cache || 'default', options?.credentials || 'same-origin']);
    return original.apply(this, arguments);
  };
})()"#;

fn demo_snapshot(b: &mut Browser, url: &str, timeout: Duration) -> Option<Vec<Value>> {
    go(b, url);
    b.wait(
        "document.readyState === 'complete' && document.documentElement.dataset.pgpResults          && document.querySelector('[data-pgp-map=drawn] canvas')          && document.querySelector('#rail a[href*=DEMO-FAMILY]')          && document.documentElement.dataset.pgpCharts",
        timeout,
    )?;
    eval(b, "location.hash = '#s3'; true");
    b.wait(
        "document.querySelector('#rail [aria-current=true]')?.getAttribute('data-toc') === '3'",
        timeout,
    )?;
    let Value::Array(snapshot) = eval(b, CAPTURE) else {
        panic!("missing demo snapshot")
    };
    Some(snapshot)
}

fn gate_transparency(f: &mut Fixture, base: &str, core: &Path, source: &str) {
    let results = f.site.join(format!("runs/{CHILD}/results.json"));
    let saved = fs::read(&results).unwrap();
    let mut no_acts = doc(&results);
    let obj = obj(&mut no_acts);
    array(obj.get_mut("findings").unwrap())
        .retain(|finding| text(field(finding, "level")) != "act");
    if let Some(Value::Object(levels)) = obj.get_mut("levels") {
        levels.insert("act".into(), Value::Int(0));
    }
    save(&results, &no_acts);
    let (prefix, _) = source.split_once("  // --- the report gate").unwrap();
    let bare = format!(
        "{prefix}
  PGPUI.ready = Promise.resolve();
  PGPUI.gateOwnsPage = function () {{ return false; }};
  PGPUI.openedActs = function () {{ return {{}}; }};
  PGPUI.closeAct = function () {{ return false; }};
}}());"
    );
    let b = &mut f.browser;
    let log = install(b, REQUESTS);
    let mut reference = None;
    for run in [ADULT, CHILD] {
        let url = run_url(base, run);
        fs::write(core, &bare).unwrap();
        let bare = demo_snapshot(b, &url, WAIT)
            .unwrap_or_else(|| panic!("{run}: bare core did not render"));
        fs::write(core, source).unwrap();
        let gated =
            demo_snapshot(b, &url, WAIT).unwrap_or_else(|| panic!("{run}: core did not render"));
        for (i, label) in [
            "DOM and attributes",
            "address",
            "canvas",
            "resources",
            "fetches",
            "fetch order",
        ]
        .iter()
        .enumerate()
        {
            if i == 0 {
                // Rendered pages are large. Print just the first differing element.
                let diff = text(&bare[i])
                    .split("><")
                    .zip(text(&gated[i]).split("><"))
                    .find(|(a, b)| a != b)
                    .map(|(a, b)| (a.to_owned(), b.to_owned()));
                assert!(bare[i] == gated[i], "{run}: {label} differs: {diff:?}");
            } else {
                assert!(bare[i] == gated[i], "{run}: {label} differs");
            }
        }
        if run == CHILD {
            reference = Some(gated);
        }
    }
    let url = run_url(base, CHILD);
    for (from, to) in [
        ("ownsPage = shut(current);", "ownsPage = true;"),
        ("var kept = shared[requestURL];", "var kept = null;"),
    ] {
        fs::write(core, replace_once(source, from, to)).unwrap();
        let planted = demo_snapshot(b, &url, Duration::from_secs(5));
        assert!(
            planted != reference,
            "transparency control went unnoticed: {to}"
        );
        println!("red control: {to} changed the no-act child's render or requests");
    }
    uninstall(b, &log);
    fs::write(results, saved).unwrap();
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
    println!("red control: a 1.5 s ticket left the report closed");
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
    gate_transparency(&mut f, &base, &core, &source);
    gate_slow_release(&mut f.browser, &base, &core, &source);
    gate_family_independence(&mut f, &base, &core, &source);
}

struct LiveServe {
    child: Child,
    base: String,
}

impl LiveServe {
    fn start(f: &Fixture, assets: &Path) -> Self {
        let mut server = Self {
            child: Command::new(&f.og)
                .args(["live", "serve", "--port", "0", "--quiet", "--root"])
                .arg(&f.work)
                .arg("--state-dir")
                .arg(f.scratch.path().join("server-state"))
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
        for line in BufReader::new(stderr).lines() {
            let line = line.unwrap();
            if let Some(url) = line.strip_prefix("og live serve: browser ") {
                server.base = url
                    .split_whitespace()
                    .next()
                    .unwrap()
                    .trim_end_matches('/')
                    .to_owned();
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

fn family_gate_checks(f: &mut Fixture, base: &str, core: &Path, source: &str) {
    let partners = f.work.join(PARTNERS);
    fs::create_dir(&partners).unwrap();
    fs::write(partners.join("members.tsv"), format!(
        "{MEMBERS_HEADER}DEMO-003\tmale\t\t\tDEMO-005\tyes\tyes\nDEMO-005\tfemale\t\t\tDEMO-003\tyes\tno\n"
    )).unwrap();

    go(&mut f.browser, &run_url(base, FAMILY));
    wait(&mut f.browser, &format!("!!{PARTNER_GATE}"));
    assert_eq!(
        text(&eval(
            &mut f.browser,
            &format!("{FAMILY_CARD}.querySelector('.person-look').textContent")
        )),
        "Report closed"
    );
    // The child remains closed on the partner's family page even though it sorts first.
    let url = run_url(base, PARTNERS);
    let expected = loads(r#"["Report closed","ch-members",true,true]"#).unwrap();
    assert_eq!(partner_view(&mut f.browser, &url), expected);
    let planted = replace_once(
        source,
        "return childRun(f, m.sample);",
        "return childRun(f, m.sample) && (f.state.children || []).some(c => c.child === m.sample);",
    );
    fs::write(core, planted).unwrap();
    assert_ne!(partner_view(&mut f.browser, &url), expected);
    println!("red control: deciding from this pedigree alone exposed the cross-family partner");
    fs::write(core, source).unwrap();
    fs::remove_dir_all(partners).unwrap();
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

fn failure_browser_checks(f: &mut Fixture, base: &str) {
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

    fs::set_permissions(relationships, mode).unwrap();
    assert_eq!(text(&run_state(&mut f.browser, base, CHILD)), "open");
    assert_eq!(
        home_states(&mut f.browser, base),
        loads(r#"["open","open"]"#).unwrap()
    );
    println!(
        "red control: letting the bad family file read makes the closed-page and tile assertions fail"
    );
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
                println!(
                    "red control: treating a failed family fetch as closed changed the adult report"
                );
            } else {
                assert!(failed == healthy, "{fault} changed the adult report");
                println!("family.json {fault}: the adult report stayed open");
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
    let server = LiveServe::start(&f, &assets);
    family_gate_checks(&mut f, &server.base, &core, &source);
    failure_browser_checks(&mut f, &server.base);
    family_fetch_failure_checks(&mut f, &server.base, &assets.join("person.js"));
}
