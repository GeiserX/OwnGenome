//! Exercise the child report gate with the result and index files the demo
//! builds. Check what Chrome presents on the person, finding, home and
//! family pages, including when an unopened result's data changes.

use og_browser::{Browser, Scratch, StaticServer, Value, chrome_or_skip, js_str, repo_root};
use pycompat::json::{dumps, loads};
use std::{fs, path::Path, path::PathBuf, process::Command, time::Duration};

const CHILD: &str = "/runs/DEMO-003/";
const FAMILY: &str = "/runs/DEMO-FAMILY/";
const CLOSED: &str = "[data-opening^=closed]";
const GATED: &str = "document.documentElement.dataset.pgpGated === 'true'";
const OPEN: &str =
    "!document.documentElement.dataset.pgpGated && !!document.querySelector('#results .finding')";
const WAIT: Duration = Duration::from_secs(15);

fn command(mut cmd: Command) {
    let out = cmd.output().expect("start fixture command");
    assert!(
        out.status.success(),
        "{cmd:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

fn value(path: &Path) -> Value {
    loads(&fs::read_to_string(path).unwrap()).unwrap()
}

fn object(v: &mut Value) -> &mut indexmap::IndexMap<String, Value> {
    let Value::Object(o) = v else {
        panic!("expected object")
    };
    o
}

fn write(path: &Path, doc: &Value) {
    fs::write(path, dumps(doc, &Default::default()).unwrap()).unwrap();
}

fn field<'a>(v: &'a Value, key: &str) -> &'a Value {
    let Value::Object(o) = v else {
        panic!("expected object with {key}");
    };
    o.get(key).unwrap_or_else(|| panic!("missing {key}"))
}

fn text<'a>(v: &'a Value, key: &str) -> &'a str {
    let Value::Str(s) = field(v, key) else {
        panic!("expected text at {key}")
    };
    s
}

fn wait(b: &mut Browser, expr: &str) -> Value {
    b.wait(expr, WAIT)
        .unwrap_or_else(|| panic!("waiting for {expr}: {:?}", b.page_errors()))
}

fn go(b: &mut Browser, srv: &StaticServer, route: &str) {
    b.go("about:blank").unwrap();
    wait(b, "location.href === 'about:blank'");
    b.go(&srv.url(route)).unwrap();
    wait(
        b,
        "document.readyState === 'complete' && document.documentElement.dataset.pgpStatus === 'ready'",
    );
}

fn count(b: &mut Browser, selector: &str, expected: usize) {
    wait(
        b,
        &format!(
            "document.querySelectorAll({}).length === {expected}",
            js_str(selector)
        ),
    );
}

fn click(b: &mut Browser, selector: &str, index: usize) {
    b.eval(&format!(
        "document.querySelectorAll({})[{index}].click(); true",
        js_str(selector)
    ))
    .unwrap();
}

fn closed_view(b: &mut Browser, count_expected: usize) -> Value {
    wait(b, GATED);
    count(b, CLOSED, count_expected);
    count(
        b,
        "input[type=search], #results, #genome-map, canvas, svg, img, #report-link",
        0,
    );
    b.eval("JSON.stringify([location.pathname, location.search, location.hash, document.body.innerHTML, sessionStorage.length, window.name])").unwrap()
}

fn closed_snapshot(b: &mut Browser, expected: usize) -> (Value, Vec<Value>) {
    closed_view(b, expected);
    b.eval(include_str!("../demo/gate-property.js")).unwrap();
    let frame = b
        .eval(
            r#"(() => {
      const visible = GateProperty.snapshot();
      visible.requests = visible.requests.filter(path => path !== '/favicon.ico');
      return visible;
    })()"#,
        )
        .unwrap();
    let tree = b
        .call("Accessibility.getFullAXTree", loads("{}").unwrap())
        .unwrap();
    let Value::Array(nodes) = field(&tree, "nodes") else {
        panic!("expected accessibility tree");
    };
    let names = nodes
        .iter()
        .filter_map(|node| {
            if field(node, "ignored") == &Value::Bool(true) {
                return None;
            }
            let Value::Object(node) = node else {
                return None;
            };
            node.get("name").map(|name| field(name, "value").clone())
        })
        .collect();
    (frame, names)
}

fn require_index_tags(index: &Value, acts: &[Value]) {
    let Value::Array(items) = field(field(index, "terms"), "KCNH2") else {
        panic!("the planted gene has no search entries")
    };
    assert!(!items.is_empty());
    for item in items {
        let id = text(item, "veil");
        assert!(
            acts.iter().any(|act| text(act, "id") == id),
            "search tag {id} must name a finding to act on"
        );
    }
}

#[test]
fn child_results_stay_closed_until_every_finding_is_opened() {
    let Some(chrome) = chrome_or_skip("check the live results veil") else {
        return;
    };
    let Some(og) = std::env::var_os("OG").map(PathBuf::from).or_else(|| {
        let test = std::env::current_exe().ok()?;
        let og = test.parent()?.parent()?.join("og");
        og.is_file().then_some(og)
    }) else {
        assert_ne!(std::env::var("CI").as_deref(), Ok("true"), "build og in CI");
        println!("[SKIP] no og to build the live demo");
        return;
    };
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "results-veil");
    let site = scratch.path().join("site");
    let work = scratch.path().join("work");
    fs::create_dir(&work).unwrap();
    let mut build = Command::new("python3");
    build
        .arg(repo_root().join("tests/demo/build_live_demo.py"))
        .args(["--records", "2000", "--out"])
        .arg(&site)
        .arg("--work")
        .arg(&work)
        .arg("--og")
        .arg(&og);
    command(build);

    let result_path = site.join("runs/DEMO-003/results.json");
    let original = value(&result_path);
    let Value::Array(findings) = field(&original, "findings") else {
        panic!("the demo has no findings")
    };
    let (mut acts, others): (Vec<_>, Vec<_>) = findings
        .iter()
        .cloned()
        .partition(|f| text(f, "level") == "act");
    let act_count = acts.len();
    assert!(act_count >= 2, "the demo needs two independent acts");
    require_index_tags(
        &value(&site.join("runs/DEMO-003/DEMO-003_search_index.json")),
        &acts,
    );
    // Exercise both native index entry points as well as the legacy demo's.
    let run = work.join("DEMO-003");
    let index = scratch.path().join("index.json");
    let mut index_build = Command::new(&og);
    index_build
        .args(["report", "index", "build", "--summary"])
        .arg(run.join("summary.json"))
        .arg("--sample-dir")
        .arg(&run)
        .arg("-o")
        .arg(&index);
    command(index_build);
    require_index_tags(&value(&index), &acts);
    let mut report = Command::new(&og);
    report
        .args(["report", "render", "--summary"])
        .arg(run.join("summary.json"))
        .arg("--sample-dir")
        .arg(&run)
        .arg("-o")
        .arg(scratch.path().join("report.txt"))
        .arg("--search-index")
        .arg(&index);
    command(report);
    require_index_tags(&value(&index), &acts);

    let mut third = acts[0].clone();
    object(&mut third).extend([
        ("id".into(), Value::Str("gate-fixture-third".into())),
        (
            "title".into(),
            Value::Str("Separate synthetic result".into()),
        ),
        ("viewer".into(), Value::Null),
    ]);
    acts.truncate(2);
    acts.push(third);
    let fid = text(&acts[0], "id");
    let finding = format!("{CHILD}finding/{fid}/");
    let guessed = format!("{CHILD}finding/gate-fixture-unknown/");
    let guessed_page = site.join(guessed.trim_start_matches('/'));
    fs::create_dir_all(&guessed_page).unwrap();
    let finding_page = site
        .join(finding.trim_start_matches('/'))
        .join("index.html");
    let html = fs::read_to_string(finding_page).unwrap();
    fs::write(
        guessed_page.join("index.html"),
        html.replace(fid, "gate-fixture-unknown"),
    )
    .unwrap();
    let srv = StaticServer::start(&site, &[]).unwrap();
    let mut b = Browser::launch(&chrome, &scratch.path().join("chrome")).unwrap();
    b.call("Network.enable", loads("{}").unwrap()).unwrap();
    b.call(
        "Network.setCacheDisabled",
        loads(r#"{"cacheDisabled":true}"#).unwrap(),
    )
    .unwrap();

    for (width, media) in [(1440, "screen"), (390, "screen"), (1440, "print")] {
        b.call(
            "Emulation.setDeviceMetricsOverride",
            loads(&format!(
                r#"{{"width":{width},"height":900,"deviceScaleFactor":1,"mobile":false}}"#
            ))
            .unwrap(),
        )
        .unwrap();
        b.call(
            "Emulation.setEmulatedMedia",
            loads(&format!(r#"{{"media":"{media}"}}"#)).unwrap(),
        )
        .unwrap();
        for route in [
            CHILD.to_owned(),
            finding.clone(),
            format!("{CHILD}?chr=7#f/{fid}"),
        ] {
            let mut altered = original.clone();
            let secret = "unopened fixture payload";
            let changed = |mut f: Value| {
                object(&mut f).insert("title".into(), Value::Str(secret.into()));
                object(&mut f).insert("future_field".into(), Value::Array(vec![Value::Bool(true)]));
                f
            };
            let Value::Array(altered_findings) = &mut object(&mut altered)["findings"] else {
                panic!("the demo has no findings")
            };
            for f in altered_findings {
                *f = changed(f.clone());
            }
            object(&mut altered).insert("future_field".into(), Value::Str(secret.into()));
            write(&result_path, &original);
            go(&mut b, &srv, &route);
            let before = closed_snapshot(&mut b, act_count);
            write(&result_path, &altered);
            go(&mut b, &srv, &route);
            assert_eq!(
                closed_snapshot(&mut b, act_count),
                before,
                "{width} {media} {route}"
            );
        }
    }
    b.call(
        "Emulation.setEmulatedMedia",
        loads(r#"{"media":"screen"}"#).unwrap(),
    )
    .unwrap();

    // The indices give every order in which the remaining controls can be selected.
    let orders: &[&[usize]] = &[
        &[0],
        &[0, 0],
        &[1, 0],
        &[0, 0, 0],
        &[0, 1, 0],
        &[1, 0, 0],
        &[1, 1, 0],
        &[2, 0, 0],
        &[2, 1, 0],
    ];
    for order in orders {
        let total = order.len();
        let mut doc = original.clone();
        object(&mut doc).insert(
            "findings".into(),
            Value::Array(others.iter().chain(&acts[..total]).cloned().collect()),
        );
        write(&result_path, &doc);
        go(&mut b, &srv, CHILD);
        for (stage, &which) in order.iter().enumerate() {
            let left = total - stage;
            let before = closed_view(&mut b, left);
            click(&mut b, CLOSED, which);
            if left > 1 {
                count(&mut b, ".opened-act", stage + 1);
                closed_view(&mut b, left - 1);
                click(&mut b, "[data-opening^=opened]", stage);
                assert_eq!(closed_view(&mut b, left), before, "close again");
                click(&mut b, CLOSED, which);
                closed_view(&mut b, left - 1);
            }
        }
        wait(&mut b, OPEN);
        wait(
            &mut b,
            "!!document.querySelector('input[type=search]:not(:disabled)')",
        );
        b.eval("var search = document.querySelector('input[type=search]'); search.value = 'KCNH2'; search.dispatchEvent(new Event('input')); true").unwrap();
        wait(
            &mut b,
            "document.querySelector('.search-results').textContent.includes('Results this page keeps closed are left out of the search')",
        );
        b.eval(&format!("PGPUI.closeAct({}); true", js_str(fid)))
            .unwrap();
        closed_view(&mut b, 1);
        count(&mut b, ".opened-act", total - 1);
        b.call("Page.reload", loads("{}").unwrap()).unwrap();
        closed_view(&mut b, total);
    }

    write(&result_path, &original);
    go(&mut b, &srv, &finding);
    let known = closed_view(&mut b, act_count);
    go(&mut b, &srv, &guessed);
    assert_eq!(closed_view(&mut b, act_count), known);
    for left in (1..=act_count).rev() {
        count(&mut b, CLOSED, left);
        click(&mut b, CLOSED, 0);
    }
    wait(&mut b, OPEN);
    assert_eq!(
        b.eval("location.pathname").unwrap(),
        Value::Str(CHILD.into())
    );
    go(&mut b, &srv, FAMILY);
    wait(
        &mut b,
        "!!document.querySelector('[data-gate-child=DEMO-003]')",
    );
    assert_eq!(
        b.eval("document.body.innerText.includes('KCNH2')").unwrap(),
        Value::Bool(false)
    );
    go(&mut b, &srv, "/");
    assert_eq!(
        b.eval("document.body.innerText.includes('KCNH2')").unwrap(),
        Value::Bool(false)
    );
    go(&mut b, &srv, "/runs/DEMO-002/");
    wait(&mut b, OPEN);
    for widget in ["str", "ci", "vaf", "copies"] {
        wait(
            &mut b,
            &format!("!!document.querySelector('[data-widget={widget}]')"),
        );
    }
}
