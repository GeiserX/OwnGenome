//! Static-site contract and red controls for the public demo.

use flate2::read::GzDecoder;
use og_browser::{Scratch, asset_problems, repo_root};
use og_live::{csp, static_files};
use pycompat::json::{Value, loads};
use pycompat::value::truthy;
use regex::Regex;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const PEOPLE: [&str; 3] = ["DEMO-001", "DEMO-002", "DEMO-003"];
const FAMILY: &str = "DEMO-FAMILY";

fn og() -> Option<PathBuf> {
    let path = std::env::var_os("OG").map(PathBuf::from).or_else(|| {
        let exe = std::env::current_exe().ok()?;
        Some(exe.parent()?.parent()?.join("og"))
    });
    if path.as_ref().is_some_and(|p| p.is_file()) {
        return path;
    }
    assert_ne!(
        std::env::var("CI").as_deref(),
        Ok("true"),
        "build og or set OG"
    );
    eprintln!("[SKIP] build og or set OG to check the public demo");
    None
}

fn command(og: &Path, argv: &[&str], paths: &[&Path]) -> Output {
    Command::new(og)
        .args(argv)
        .args(paths)
        .current_dir(repo_root())
        .output()
        .unwrap()
}

fn pass(result: Output) {
    assert!(
        result.status.success(),
        "{}{}",
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
}

fn json(path: &Path) -> Value {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    loads(&text).unwrap()
}

fn field<'a>(v: &'a Value, key: &str) -> &'a Value {
    if let Value::Object(o) = v {
        o.get(key).unwrap_or(&Value::Null)
    } else {
        &Value::Null
    }
}

fn list(v: &Value) -> &[Value] {
    if let Value::Array(a) = v { a } else { &[] }
}

fn text(v: &Value) -> Option<&str> {
    if let Value::Str(s) = v { Some(s) } else { None }
}

fn step_state(fr: &Value, name: &str) -> Option<&'static str> {
    let step = list(field(fr, "steps"))
        .iter()
        .find(|step| text(field(step, "step")) == Some(name))?;
    match text(field(step, "state")) {
        Some("done") => Some("done"),
        Some("waiting") => Some("waiting"),
        _ => None,
    }
}

fn frames_of(path: &Path) -> Vec<Value> {
    let doc = json(path);
    assert_eq!(field(&doc, "interval_ms"), &Value::Int(600));
    list(field(&doc, "frames")).to_vec()
}

fn plugins(frame: &Value) -> &Value {
    field(field(frame, "state"), "plugins")
}

fn missing_hooks(last: &Value, family: &Value) -> Vec<&'static str> {
    let last = plugins(last);
    ["charts", "ops", "recheck", "search", "family"]
        .into_iter()
        .filter(|name| {
            let value = field(last, name);
            !truthy(value)
                || truthy(field(value, "error"))
                || (*name == "family"
                    && (text(field(value, "kind")) != Some("member")
                        || text(field(field(plugins(family), "family"), "kind")) != Some("family")))
        })
        .collect()
}

fn follow_problems(frames: &[Value], planned: &BTreeSet<String>) -> Vec<String> {
    let mut bad = Vec::new();
    let first = field(&frames[0], "state");
    let last = field(frames.last().unwrap(), "state");
    for (state, expected) in [(first, "running"), (last, "finished")] {
        if text(field(field(state, "run"), "state")) != Some(expected) {
            bad.push(format!("the replay should be {expected}"));
        }
    }
    for name in planned {
        if step_state(last, name) != Some("done") {
            bad.push(format!("step {name} does not end done"));
        }
    }
    for (i, frame) in frames.iter().enumerate() {
        let report = step_state(field(frame, "state"), "24") == Some("done");
        for name in ["charts", "search"] {
            if truthy(field(field(plugins(frame), name), "file")) != report {
                bad.push(format!("frame {i} shows {name} at the wrong step"));
            }
        }
        if truthy(field(plugins(frame), "recheck")) != (i + 1 == frames.len()) {
            bad.push(format!("frame {i} shows the re-check at the wrong time"));
        }
    }
    if !truthy(field(
        field(plugins(frames.last().unwrap()), "recheck"),
        "steps",
    )) {
        bad.push("the last frame has no re-check badge".into());
    }
    bad
}

fn inside_site(root: &Path, base: &Path, link: &str) -> bool {
    if link.starts_with("https://") {
        return true;
    }
    if link.starts_with('/') {
        return false;
    }
    let path = base.join(link);
    let file = if link.ends_with('/') {
        path.join("index.html")
    } else {
        path
    };
    file.is_file() && file.canonicalize().is_ok_and(|p| p.starts_with(root))
}

fn links<'a>(value: &'a Value, output: &mut Vec<&'a str>) {
    match value {
        Value::Object(o) => {
            for (key, v) in o {
                if matches!(
                    key.as_str(),
                    "url" | "pedigree" | "segments" | "report" | "table" | "file"
                ) && let Some(s) =
                    text(v).filter(|s| !s.is_empty() && !s.bytes().all(|b| b.is_ascii_digit()))
                {
                    output.push(s);
                } else {
                    links(v, output);
                }
            }
        }
        Value::Array(a) => a.iter().for_each(|v| links(v, output)),
        _ => {}
    }
}

fn problems(root: &Path, base: &Path, doc: &Value) -> Vec<String> {
    let mut urls = Vec::new();
    links(doc, &mut urls);
    for chapter in list(field(doc, "chapters")) {
        for section in list(field(chapter, "sections")) {
            urls.extend(list(field(section, "links")).iter().filter_map(text));
        }
    }
    urls.into_iter()
        .filter(|u| !inside_site(root, base, u))
        .map(str::to_owned)
        .collect()
}

fn no_host_path(frames: &str, paths: &[&Path]) -> bool {
    paths.iter().all(|p| !frames.contains(p.to_str().unwrap()))
}

fn collect_files(dir: &Path, files: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(dir).unwrap() {
        let entry = entry.unwrap();
        if entry.file_type().unwrap().is_dir() {
            collect_files(&entry.path(), files);
        } else {
            files.push(entry.path());
        }
    }
}

fn report_problems(root: &Path, run: &Path, html: &str) -> Vec<String> {
    let href = Regex::new(r#"(?i)\bhref\s*=\s*["']([^"']+)["']"#).unwrap();
    let scheme = Regex::new(r"^[A-Za-z][A-Za-z0-9+.-]*:").unwrap();
    href.captures_iter(html)
        .map(|m| static_files::unescape(m[1].split(['#', '?']).next().unwrap()))
        .map(|link| og_live::viewer::unquote(&link))
        .filter(|link| !link.is_empty() && !scheme.is_match(link))
        .filter(|link| !inside_site(root, run, link))
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

fn object(value: &mut Value) -> &mut indexmap::IndexMap<String, Value> {
    let Value::Object(fields) = value else {
        panic!("expected a JSON object");
    };
    fields
}

fn chart_fixture(og: &Path, run: &Path, summary: &Path, html: &Path, charts: &Path) {
    pass(
        Command::new(og)
            .args(["report", "render", "--summary"])
            .arg(summary)
            .arg("--sample-dir")
            .arg(run)
            .arg("-o")
            .arg(html)
            .arg("--charts-json")
            .arg(charts)
            .current_dir(repo_root())
            .output()
            .unwrap(),
    );
}

fn fragment_links(html: &str) -> BTreeSet<String> {
    let id = Regex::new(r#"\sid="([^"]+)""#).unwrap();
    let ids: BTreeSet<_> = id.captures_iter(html).map(|m| m[1].to_owned()).collect();
    let href = Regex::new(r##"href="#([^"]+)""##).unwrap();
    href.captures_iter(html)
        .map(|m| static_files::unescape(&m[1]))
        .filter(|href| !ids.contains(href))
        .collect()
}

fn charts_contract(og: &Path, work: &Path, scratch: &Path) {
    let run = work.join(PEOPLE[0]);
    let summary = run.join("summary.json");
    let mut rendered = Vec::new();
    for i in 1..=2 {
        let html = scratch.join(format!("charts-{i}.html"));
        let charts = scratch.join(format!("charts-{i}.json"));
        chart_fixture(og, &run, &summary, &html, &charts);
        let mut doc = json(&charts);
        object(&mut doc).shift_remove("generated_utc");
        rendered.push(doc);
    }
    assert_eq!(
        rendered[0], rendered[1],
        "the same summary generated different charts"
    );
    assert_eq!(field(&rendered[0], "schema_version"), &Value::Int(1));
    let ids: Vec<_> = list(field(&rendered[0], "charts"))
        .iter()
        .filter_map(|chart| text(field(chart, "id")))
        .collect();
    assert_eq!(
        ids,
        ["clinvar-stars", "genome", "mito", "ancestry-fine", "prs"]
    );
    let mut planted = rendered[1].clone();
    let Some(Value::Array(charts)) = object(&mut planted).get_mut("charts") else {
        panic!("no charts")
    };
    let Some(Value::Array(marks)) = object(&mut charts[0]).get_mut("marks") else {
        panic!("no ClinVar marks")
    };
    object(&mut marks[0]).insert("value".into(), Value::Int(-1));
    assert_ne!(
        rendered[0], planted,
        "a changed chart value must be visible"
    );

    // The larger fixture includes a hit after the report table's row cap.
    let mut merged = json(&summary);
    let mut fixture = json(&repo_root().join("tests/fixtures/report_charts/sections.json"));
    let fixture = object(&mut fixture).shift_remove("sections").unwrap();
    let Value::Object(sections) = fixture else {
        panic!("no chart sections")
    };
    object(object(&mut merged).get_mut("sections").unwrap()).extend(sections);
    let overlaid = scratch.join("charts-summary.json");
    fs::write(&overlaid, og_live::state::dumps_compact(&merged)).unwrap();
    let html = scratch.join("charts-overlay.html");
    chart_fixture(
        og,
        &run,
        &overlaid,
        &html,
        &scratch.join("charts-overlay.json"),
    );
    let html = fs::read_to_string(html).unwrap();
    assert!(!html.to_ascii_lowercase().contains("<script"));
    assert_eq!(fragment_links(&html), BTreeSet::new());
    assert_eq!(
        fragment_links(&html.replace(r#"id="clinvar-hits""#, r#"id="gone""#)),
        BTreeSet::from(["clinvar-hits".into()])
    );
    for (id, caveat) in [
        ("clinvar", og_render::notes::CLINVAR_FREQ_NOTE),
        ("roh", og_render::notes::ROH_NOTE),
        ("mito", og_render::notes::MITO_NOTE),
        ("ancestry-fine", og_render::notes::ANCESTRY_FINE_NOTE),
        ("prs", og_render::notes::PRS_ANCESTRY_NOTE),
    ] {
        let section = Regex::new(&format!(
            r#"(?s)<section class="sec[^"]*" id="{id}">(.*?)\n  </section>"#
        ))
        .unwrap();
        let body = &section
            .captures(&html)
            .unwrap_or_else(|| panic!("missing chart section {id}"))[1];
        let figure = r#"<figure class="chart">"#;
        assert_eq!(
            body.matches(figure).count(),
            1,
            "{id} should contain one figure"
        );
        assert!(
            body.find(&static_files::esc(caveat)).unwrap() < body.find(figure).unwrap(),
            "{id} shows its chart before the caveat"
        );
    }
    assert_eq!(html.matches(r#"<a id="clinvar-hits"></a>"#).count(), 1);
    assert!(html.matches("<tr><td>chr").count() >= 50);
}

#[test]
fn live_charts_render_after_sign_in() {
    use og_browser::{Browser, chrome_or_skip, js_str};
    use og_live::auth::{Auth, hash_password};
    use og_live::server::{Config, Server};
    use std::sync::Arc;
    use std::time::{Duration, UNIX_EPOCH};

    let Some(chrome) = chrome_or_skip("read charts on the live page") else {
        return;
    };
    let Some(og) = og() else { return };
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "demo-live-charts");
    let site = scratch.path().join("site");
    let work = scratch.path().join("work");
    pass(
        Command::new(&og)
            .args(["live", "demo", "--records", "2000", "--out"])
            .arg(&site)
            .arg("--work")
            .arg(&work)
            .current_dir(repo_root())
            .output()
            .unwrap(),
    );

    let run = work.join(PEOPLE[0]);
    let charts_path = run.join(format!("{}_charts.json", PEOPLE[0]));
    let test_run = scratch.path().join(PEOPLE[0]);
    fs::create_dir(&test_run).unwrap();
    assert!(!truthy(&og_live::charts::state(&test_run)));
    let copy = test_run.join(charts_path.file_name().unwrap());
    fs::copy(charts_path, &copy).unwrap();
    let plugin = og_live::charts::state(&test_run);
    assert_eq!(text(field(&plugin, "file")), Some("DEMO-001_charts.json"));
    let Value::Int(stamp) = field(&plugin, "mtime") else {
        panic!("plugin has no chart mtime")
    };
    fs::File::open(copy)
        .unwrap()
        .set_modified(UNIX_EPOCH + Duration::from_nanos(u64::try_from(*stamp + 1000).unwrap()))
        .unwrap();
    assert_eq!(
        field(&og_live::charts::state(&test_run), "mtime"),
        &Value::Int(*stamp + 1000)
    );

    let repo = repo_root();
    let assets = site.join("assets");
    // Use the same unfilled pages the live server loads, beside the site's
    // existing script and style files.
    for page in static_files::PAGES {
        fs::copy(repo.join("assets/live").join(page), assets.join(page)).unwrap();
    }
    let route = format!("/runs/{}/", PEOPLE[0]);
    let password = "synthetic chart test password";
    for (tag, show) in [("charts", true), ("no-charts", false)] {
        if !show {
            fs::remove_file(assets.join("plugins/charts.js")).unwrap();
        }
        let mut cfg = Config::new(&work);
        cfg.repo_assets = repo.join("assets");
        cfg.assets = assets.clone();
        cfg.steps = cfg.repo_assets.join("steps.tsv");
        cfg.port = 0;
        cfg.quiet = true;
        cfg.og = og.clone();
        cfg.auth = Arc::new(Auth::new(Some(hash_password(password, 1000).unwrap())));
        cfg.results = Some(og_live::results::FNS);
        let server = Server::start(cfg).unwrap();
        let base = format!("http://127.0.0.1:{}", server.port());
        let mut b = Browser::launch(&chrome, &scratch.path().join(tag)).unwrap();
        b.call(
            "Emulation.setEmulatedMedia",
            loads(r#"{"features":[{"name":"prefers-color-scheme","value":"dark"}]}"#).unwrap(),
        )
        .unwrap();
        b.go(&format!("{base}{route}")).unwrap();
        let wait = Duration::from_secs(15);
        assert!(
            b.wait(
                "location.pathname === '/login' && !!document.getElementById('password')",
                wait
            )
            .is_some(),
            "no login form"
        );
        b.eval(&format!("document.getElementById('password').value = {}; document.querySelector('form').submit(); true",
            js_str(password))).unwrap();
        assert!(b.wait(&format!("location.pathname === {} && document.documentElement.dataset.pgpStatus === 'ready'",
            js_str(&route)), wait).is_some(), "the signed-in run page did not load");
        if show {
            assert!(
                b.wait("document.documentElement.dataset.pgpCharts === '5'", wait)
                    .is_some(),
                "the live page did not draw five charts"
            );
            b.eval(r#"document.querySelector('figure[data-chart="genome"] a[data-mark^="hit-"]')
                    .dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true})); true"#).unwrap();
            let detail = format!(
                r#"(() => {{
                const panel = document.querySelector('.chart-detail');
                return panel && !panel.hidden && panel.textContent.includes({})
                    && /^viewer\/#chr[0-9XY]+:[0-9]+$/.test(panel.querySelector('a.viewer')?.getAttribute('href'))
                    && location.pathname === {};
            }})()"#,
                js_str(og_render::notes::CLINVAR_FREQ_NOTE),
                js_str(&route)
            );
            assert!(
                b.wait(&detail, wait).is_some(),
                "the live chart click did not open its caveat and viewer link"
            );
        } else {
            assert!(
                b.wait(
                    "Boolean(document.documentElement.dataset.pgpCharts)",
                    Duration::from_millis(700)
                )
                .is_none(),
                "the planted demo without charts.js drew charts"
            );
        }
        assert_eq!(b.eval(r#"performance.getEntriesByType('resource').every(entry => new URL(entry.name).origin === location.origin)"#).unwrap(), Value::Bool(true));
        assert_eq!(
            b.eval("document.documentElement.dataset.pgpBlocked || ''")
                .unwrap(),
            Value::Str(String::new())
        );
    }
}

#[test]
fn the_weblog_matches_its_catalogue() {
    let Some(og) = og() else { return };
    pass(command(&og, &["dev", "demo-events", "--check"], &[]));
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "demo-events");
    let path = scratch.path().join("truncated.jsonl");
    let fixture =
        fs::read_to_string(repo_root().join("tests/fixtures/live/demo-001.weblog.jsonl")).unwrap();
    let (truncated, _) = fixture.trim_end().rsplit_once('\n').unwrap();
    fs::write(&path, format!("{truncated}\n")).unwrap();
    let bad = command(&og, &["dev", "demo-events", "--check", "--out"], &[&path]);
    assert_eq!(bad.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&bad.stderr).contains("differs from the generated weblog"));
    pass(command(&og, &["dev", "demo-events", "--out"], &[&path]));
    assert_eq!(fs::read_to_string(path).unwrap(), fixture);
}

#[test]
fn the_built_site_keeps_the_demo_contract() {
    let Some(og) = og() else { return };
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "live-demo");
    let site = scratch.path().join("site");
    let work = scratch.path().join("run");
    let result = Command::new(&og)
        .args(["live", "demo", "--records", "2000", "--out"])
        .arg(&site)
        .arg("--work")
        .arg(&work)
        .current_dir(repo_root())
        .output()
        .unwrap();
    pass(result);
    let site = site.canonicalize().unwrap();
    let primary = site.join("runs").join(PEOPLE[0]);
    let family = site.join("runs").join(FAMILY);
    let frames_file = primary.join("frames.json");
    let frames = frames_of(&frames_file);
    let first = field(&frames[0], "state");
    let last = frames.last().unwrap();
    let family_frames = frames_of(&family.join("frames.json"));

    let raw = fs::read_to_string(&frames_file).unwrap();
    let repository = repo_root().canonicalize().unwrap();
    let work_real = work.canonicalize().unwrap();
    let home = PathBuf::from(std::env::var_os("HOME").expect("test runner has HOME"));
    let paths = [scratch.path(), &work, &work_real, &repository, &home];
    assert!(
        no_host_path(&raw, &paths),
        "a frame contains a build-machine path"
    );
    assert!(!no_host_path(&format!("{raw} {}", work.display()), &paths));
    let ops = field(plugins(last), "ops");
    assert_eq!(
        field(field(ops, "memory"), "host_bytes"),
        &Value::Int(64 << 30)
    );
    assert_eq!(
        field(&list(field(ops, "disk"))[0], "size_bytes"),
        &Value::Int(2 << 40)
    );

    assert_eq!(step_state(first, "27"), Some("waiting"));
    let fixture =
        fs::read_to_string(repository.join("tests/fixtures/live/demo-001.weblog.jsonl")).unwrap();
    let events: Vec<_> = fixture.lines().map(|line| loads(line).unwrap()).collect();
    let start = field(&events[0], "utcTime");
    for frame in &frames {
        assert_eq!(
            field(field(field(frame, "state"), "run"), "started_utc"),
            start
        );
    }
    let processes: BTreeSet<_> = events
        .iter()
        .filter_map(|ev| {
            let name = text(field(field(ev, "trace"), "name"))?;
            og_live::state::task_name(name).0
        })
        .collect();
    let planned: BTreeSet<_> = og_runlog::steps::read_steps(&repository.join("assets/steps.tsv"))
        .unwrap()
        .into_iter()
        .filter(|st| st.processes.iter().any(|p| processes.contains(p)))
        .map(|st| st.step)
        .collect();
    assert!(planned.len() > 30);
    assert_eq!(follow_problems(&frames, &planned), Vec::<String>::new());
    let mut premature = frames.clone();
    let Value::Object(head) = &mut premature[0] else {
        unreachable!()
    };
    let Value::Object(st) = head.get_mut("state").unwrap() else {
        unreachable!()
    };
    let Value::Object(plugin) = st.get_mut("plugins").unwrap() else {
        unreachable!()
    };
    plugin.insert("charts".into(), field(plugins(last), "charts").clone());
    assert!(
        follow_problems(&premature, &planned)
            .contains(&"frame 0 shows charts at the wrong step".to_owned())
    );
    assert_eq!(missing_hooks(last, &family_frames[0]), Vec::<&str>::new());
    let mut dropped = last.clone();
    let Value::Object(last_frame) = &mut dropped else {
        unreachable!()
    };
    let Value::Object(st) = last_frame.get_mut("state").unwrap() else {
        unreachable!()
    };
    let Value::Object(plugin) = st.get_mut("plugins").unwrap() else {
        unreachable!()
    };
    plugin.shift_remove("charts");
    assert_eq!(missing_hooks(&dropped, &family_frames[0]), ["charts"]);

    let mut broken = Vec::new();
    for name in PEOPLE {
        let path = site.join("runs").join(name);
        let report = fs::read_to_string(path.join("report.html")).unwrap();
        let warning =
            format!("Invented data: {name} is not a real person; nothing here is a result");
        assert_eq!(
            report.lines().next().unwrap(),
            format!("<!-- {warning} -->")
        );
        let heading = report
            .lines()
            .find(|line| line.contains("Pipeline commit:"))
            .unwrap();
        assert!(!heading.contains("unknown"), "{heading}");
        assert!(
            !report.contains("No run_manifest.tsv"),
            "{name} has no manifest"
        );
        let manifest = fs::read_to_string(path.join("run_manifest.tsv")).unwrap();
        assert!(
            manifest.contains("run\tdemo\tinvented data\n"),
            "{name} has no demo marker"
        );
        assert!(
            manifest.contains("version\tPCGR_DATA_BUNDLE\t"),
            "{name} has no version pins"
        );
        let mut files = Vec::new();
        collect_files(&path, &mut files);
        let texts: Vec<_> = files
            .iter()
            .filter(|p| p.extension().is_some_and(|e| e == "txt"))
            .collect();
        assert!(
            !texts.is_empty(),
            "{name} is missing report text attachments"
        );
        for file in texts {
            assert_eq!(
                fs::read_to_string(file).unwrap().lines().next(),
                Some(warning.as_str()),
                "{}",
                file.display()
            );
        }
        for link in report_problems(&site, &path, &report) {
            broken.push(format!("{name}: {link}"));
        }
        for doc in ["results.json", "frames.json"] {
            let content = json(&path.join(doc));
            let bad = if doc == "frames.json" {
                list(field(&content, "frames"))
                    .iter()
                    .flat_map(|frame| problems(&site, &path, plugins(frame)))
                    .collect()
            } else {
                problems(&site, &path, &content)
            };
            for link in bad {
                broken.push(format!("{name}/{doc}: {link}"));
            }
        }
    }
    assert!(broken.is_empty(), "unresolved links: {}", broken.join(", "));
    for name in PEOPLE {
        let path = site.join("runs").join(name);
        let cpsr = path.join("cpsr/index.html");
        assert!(
            fs::read_to_string(&cpsr)
                .unwrap()
                .contains(&format!("CPSR {name}"))
        );
        let report = fs::read_to_string(path.join("report.html")).unwrap();
        fs::rename(&cpsr, cpsr.with_extension("held")).unwrap();
        assert!(report_problems(&site, &path, &report).contains(&"cpsr/".into()));
        fs::rename(cpsr.with_extension("held"), &cpsr).unwrap();
    }
    let mut index = json(&site.join("runs.json"));
    assert!(problems(&site, &site, &index).is_empty());
    let Value::Object(o) = &mut index else {
        unreachable!()
    };
    o.insert("url".into(), Value::Str("/runs/DEMO-001/".into()));
    assert_eq!(problems(&site, &site, &index), ["/runs/DEMO-001/"]);
    let family_json = json(&family.join("family.json"));
    assert!(problems(&site, &family, &family_json).is_empty());
    let lead = fs::read_to_string(family.join("frames.json")).unwrap();
    let planted = loads(&lead.replace("../../family/", "/family/")).unwrap();
    assert!(!problems(&site, &family, plugins(&list(field(&planted, "frames"))[0])).is_empty());

    let mut genes = String::new();
    GzDecoder::new(fs::File::open(site.join("map/genes/chr16.json")).unwrap())
        .read_to_string(&mut genes)
        .unwrap();
    let map = loads(&genes).unwrap();
    assert!(
        list(field(field(&map, "g"), "name"))
            .iter()
            .any(|v| text(v) == Some("DEMOC")),
        "chr16 slice has no DEMOC"
    );

    let mut files = Vec::new();
    collect_files(&site, &mut files);
    let mut texts = BTreeMap::new();
    let mut pages = 0;
    let policy = csp::CSP.rsplit_once("; ").unwrap().0;
    let meta = format!(r#"<meta http-equiv="Content-Security-Policy" content="{policy}">"#);
    for path in files {
        let ext = path.extension().and_then(|e| e.to_str());
        if !matches!(ext, Some("html" | "js" | "css")) {
            continue;
        }
        let raw = fs::read_to_string(&path).unwrap();
        let name = path.strip_prefix(&site).unwrap().to_str().unwrap();
        if ext == Some("html") {
            if raw.starts_with("<!-- Invented data:") {
                let limited = csp::RUN_HTML_CSP
                    .split("; ")
                    .skip(1)
                    .take_while(|p| !p.starts_with("frame-ancestors"))
                    .collect::<Vec<_>>()
                    .join("; ");
                assert!(
                    static_files::unescape(&raw).contains(&format!("content=\"{limited}\"")),
                    "{name} has no report policy"
                );
                continue;
            }
            pages += 1;
            assert!(
                static_files::unescape(&raw).contains(&meta),
                "{name} has no live policy"
            );
        }
        texts.insert(name.to_owned(), raw);
    }
    assert!(pages >= 4);
    assert_eq!(format!("{policy}; frame-ancestors 'none'"), csp::CSP);
    let external = asset_problems(&texts);
    assert!(external.is_empty(), "{external:?}");
    texts.insert(
        "planted.html".into(),
        r#"<script src="https://example.invalid/plant.js"></script>"#.into(),
    );
    assert!(!asset_problems(&texts).is_empty());
    charts_contract(&og, &work, scratch.path());
}
