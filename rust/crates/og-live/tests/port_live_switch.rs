//! The family checks of the Python site that only its server tests held,
//! now on `og live serve`, the site in production (moved from
//! tests/test_owngenome_live_family.py, which went with the Python server):
//!
//! 1. a child with a result to act on: its finding route answers the same
//!    page for its own finding, a guessed finding id and nonsense, so an
//!    address cannot tell whether the child has a result; an adult keeps
//!    the 404 for an id that is not its finding; the access log names no
//!    id. Control: the same requests with the child made an adult in the
//!    members file tell the ids apart, so the check can fail;
//! 2. family.json of the family: no flagged-card count for a child with a
//!    result to act on, its count for a child without one and for an
//!    adult, and each member run's `child` flag by the one rule;
//! 3. a family file that cannot be read: family.json answers the fixed 500
//!    "family not readable", and while it cannot be read the adults the
//!    family lists are closed like children on the finding route (fail
//!    closed); they reopen once it reads again;
//! 4. family.js: nothing from outside, no innerHTML, no eval or new
//!    Function, and `node --check` reads it (node is needed on CI). A
//!    control plants each fault in a copy.

#[path = "port_live_server_util.rs"]
mod util;

use og_live::plugins::ResultsFns;
use pycompat::json::{Value, loads};
use std::collections::BTreeSet;
use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::Path;
use std::process::Command;
use util::*;

const STEP06_OK: &str =
    r#"{"event": "pgp_step", "step": "06", "state": "ok", "utcTime": "2026-01-02T10:06:40Z"}"#;

/// The one finding to act on of a summary that has `"planted_act": true`.
const ACT: &str = "sf-planted-act";
/// A well-formed id that is no finding of any run here.
const GUESSED: &str = "two-hits-cftr";
const NONSENSE: &str = "zz-not-a-finding";

fn field<'a>(v: &'a Value, k: &str) -> Option<&'a Value> {
    match v {
        Value::Object(o) => o.get(k),
        _ => None,
    }
}

fn planted_act(summary: &Value) -> bool {
    matches!(field(summary, "planted_act"), Some(Value::Bool(true)))
}

/// The stub results: a finding to act on when the summary plants one,
/// otherwise one finding for information.
fn build(name: &str, summary: &Value, _: Option<&BTreeSet<String>>) -> Result<Value, String> {
    let (id, level) = if planted_act(summary) {
        (ACT, "act")
    } else {
        ("info-planted", "info")
    };
    loads(&format!(
        r#"{{"schema_version": 1, "sample": {}, "chapters": [], "note": null,
            "findings": [{{"id": "{id}", "gene": "PLANTGENE1", "title": "planted finding", "level": "{level}"}}]}}"#,
        og_live::state::dumps(&Value::Str(name.into()))
    ))
    .map_err(|e| e.to_string())
}

/// The stub overview: the summary's planted attention count, and one result
/// to act on when it plants one.
fn overview(summary: &Value) -> Result<Value, String> {
    let n = match field(summary, "planted_attention") {
        Some(Value::Int(n)) => *n,
        _ => 0,
    };
    let act = i32::from(planted_act(summary));
    loads(&format!(
        r#"{{"attention": {n}, "levels": {{"act": {act}}}}}"#
    ))
    .map_err(|e| e.to_string())
}

/// tests/fixtures/families_live as the genome folder, with a run for
/// child_a, father_a (a child of gf_a and gm_a), gf_a (an adult) and loner
/// (in no family), each with step 06 done and a summary.json.
struct Fam {
    w: Work,
}

impl Fam {
    fn new(name: &str) -> Fam {
        let w = Work::new(name);
        copy(&repo().join("tests/fixtures/families_live"), &w.root());
        let f = Fam { w };
        f.run(
            "child_a",
            r#"{"sections": {}, "planted_act": true, "planted_attention": 1}"#,
        );
        f.run("father_a", r#"{"sections": {}, "planted_attention": 3}"#);
        f.run(
            "gf_a",
            r#"{"sections": {}, "planted_act": true, "planted_attention": 2}"#,
        );
        f.run(
            "loner",
            r#"{"sections": {}, "planted_act": true, "planted_attention": 4}"#,
        );
        f
    }

    fn run(&self, name: &str, summary: &str) {
        let d = self.w.root().join(name);
        write(
            &d.join("logs/run_status.tsv"),
            "meta\tstarted_epoch\t1767347940\nstep\t06\tok\n",
        );
        // Step 06 ends with its ClinVar task, which gives the step its end.
        let task = format!(
            r#"{{"runName": "r", "event": "process_completed", "utcTime": "2026-01-02T10:06:40Z", "trace": {{"task_id": 9,
                "hash": "cd/000009", "name": "CLINICAL:CLINVAR_SCREEN ({name})", "process": "CLINICAL:CLINVAR_SCREEN",
                "tag": "{name}", "status": "COMPLETED", "exit": 0, "submit": 1767348359000, "start": 1767348360000,
                "complete": 1767348400000}}}}"#
        );
        store(
            &self.w,
            name,
            &[task.replace('\n', " ").as_str(), STEP06_OK],
        );
        // Newer than the end of step 06, so the site reads it rather than
        // collecting a summary of the folder.
        let p = d.join("summary.json");
        write(&p, summary);
        let later = std::time::SystemTime::now() + std::time::Duration::from_secs(3600);
        fs::File::options()
            .write(true)
            .open(&p)
            .unwrap()
            .set_modified(later)
            .unwrap();
    }

    fn config(&self) -> og_live::server::Config {
        let mut c = config(&self.w);
        c.results = Some(ResultsFns { build, overview });
        c
    }

    fn family_file(&self, f: &str) -> std::path::PathBuf {
        self.w.root().join("family_a").join(f)
    }
}

fn copy(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for e in fs::read_dir(from).unwrap() {
        let e = e.unwrap();
        let (src, dst) = (e.path(), to.join(e.file_name()));
        if e.file_type().unwrap().is_dir() {
            copy(&src, &dst);
        } else {
            fs::copy(&src, &dst).unwrap();
        }
    }
}

/// The finding route of `run` for its finding, a guessed id and nonsense:
/// (status, body with the asked id written as ID).
fn answers(port: u16, run: &str) -> Vec<(u16, String)> {
    [ACT, GUESSED, NONSENSE]
        .iter()
        .map(|id| {
            let r = get(port, &format!("/runs/{run}/finding/{id}/"));
            (r.status, r.text().replace(id, "ID"))
        })
        .collect()
}

fn statuses(a: &[(u16, String)]) -> Vec<u16> {
    a.iter().map(|(s, _)| *s).collect()
}

/// Whether the answers are one page for every id: what a child's route must give.
fn alike(a: &[(u16, String)]) -> bool {
    statuses(a) == [200; 3] && a.iter().all(|(_, b)| *b == a[0].1)
}

#[test]
fn a_childs_finding_route_cannot_tell_a_guessed_id_apart() {
    let f = Fam::new("switch-finding");
    let srv = start(f.config());
    let p = srv.port();
    let child = answers(p, "child_a");
    assert!(alike(&child), "{:?}", statuses(&child));
    assert_eq!(statuses(&answers(p, "loner")), [200, 404, 404]);
    assert_eq!(statuses(&answers(p, "gf_a")), [200, 404, 404]);
    let log = fs::read_to_string(f.w.state().join("access.log")).unwrap();
    for id in [ACT, GUESSED, NONSENSE] {
        assert!(!log.contains(id), "{id} in the access log");
    }
    drop(srv);

    // Control: child_a without its parents in the members file is an adult,
    // and the same requests tell its finding from the guessed ids.
    let members = f.family_file("members.tsv");
    let text = fs::read_to_string(&members).unwrap();
    let adult = text.replace("child_a\tfemale\tfather_a\tmother_a", "child_a\tfemale\t\t");
    assert_ne!(adult, text);
    fs::write(&members, adult).unwrap();
    let srv = start(f.config());
    let control = answers(srv.port(), "child_a");
    assert!(
        !alike(&control),
        "control: an adult's route answered every id alike"
    );
    assert_eq!(statuses(&control), [200, 404, 404]);
}

#[test]
fn family_json_gives_no_count_for_a_gated_child_and_the_child_flags() {
    let f = Fam::new("switch-family-json");
    let srv = start(f.config());
    let r = get(srv.port(), "/runs/family_a/family.json");
    assert_eq!(r.status, 200, "{}", r.text());
    let v = r.json();
    let runs = field(&v, "runs").expect("family.json has runs");
    let run = |n: &str| field(runs, n).unwrap_or_else(|| panic!("no {n} in family.json"));
    let attention = |n: &str| field(run(n), "attention").cloned();
    let child = |n: &str| field(run(n), "child").cloned();
    // A child with a result to act on: no flagged-card count, as its page is closed.
    assert_eq!(attention("child_a"), Some(Value::Null));
    // A child with none keeps its count; so does an adult with one.
    assert_eq!(attention("father_a"), Some(Value::Int(3)));
    assert_eq!(attention("gf_a"), Some(Value::Int(2)));
    assert_eq!(child("child_a"), Some(Value::Bool(true)));
    assert_eq!(child("father_a"), Some(Value::Bool(true)));
    assert_eq!(child("gf_a"), Some(Value::Bool(false)));
    // loner is in no family, so not in its page.
    assert!(field(runs, "loner").is_none());
}

#[test]
fn an_unreadable_family_file_is_a_fixed_500_and_closes_its_adults() {
    let f = Fam::new("switch-unreadable");
    let rel = f.family_file("relationships.tsv");
    let srv = start(f.config());
    let p = srv.port();
    let gf = || statuses(&answers(p, "gf_a"));
    assert_eq!(
        gf(),
        [200, 404, 404],
        "gf_a is an adult while every file reads"
    );
    let mode = fs::metadata(&rel).unwrap().permissions().mode();
    fs::set_permissions(&rel, fs::Permissions::from_mode(0o000)).unwrap();
    if fs::read(&rel).is_ok() {
        fs::set_permissions(&rel, fs::Permissions::from_mode(mode)).unwrap();
        eprintln!(
            "note: this account reads a mode 000 file (root), so the unreadable family is not checked"
        );
        return;
    }
    // Before any family.json request: closed, as a child's route is.
    let closed = answers(p, "gf_a");
    let loner = statuses(&answers(p, "loner"));
    let r = get(p, "/runs/family_a/family.json");
    fs::set_permissions(&rel, fs::Permissions::from_mode(mode)).unwrap();
    assert!(
        alike(&closed),
        "gf_a while family_a cannot be read: {:?}",
        statuses(&closed)
    );
    assert_eq!(loner, [200, 404, 404], "a run in no family stays as it was");
    assert_eq!(r.status, 500);
    assert_eq!(r.text(), "family not readable\n");
    assert!(
        r.header("content-type")
            .unwrap_or("")
            .starts_with("text/plain")
    );
    assert!(
        log_rows(&f.w).iter().any(|row| row[3..].starts_with(&[
            "family".into(),
            "family_a".into(),
            "family.json".into(),
            "500".into()
        ])),
        "{:?}",
        log_rows(&f.w)
    );
    assert_eq!(
        gf(),
        [200, 404, 404],
        "gf_a reopens once the file reads again"
    );
}

/// What family.js must not do: each fault found in `js`.
fn script_faults(js: &str) -> Vec<String> {
    let mut bad: Vec<String> = regex::Regex::new(r"(?:https?:)?//[A-Za-z0-9.-]+\.[a-z]{2,}")
        .unwrap()
        .find_iter(js)
        .map(|m| format!("loads from outside: {}", m.as_str()))
        .collect();
    if js.contains("innerHTML") {
        bad.push("sets innerHTML".into());
    }
    if regex::Regex::new(r"\beval\s*\(").unwrap().is_match(js) {
        bad.push("calls eval".into());
    }
    if regex::Regex::new(r"new\s+Function\s*\(")
        .unwrap()
        .is_match(js)
    {
        bad.push("makes a new Function".into());
    }
    bad
}

#[test]
fn family_js_loads_nothing_from_outside_and_sets_no_html() {
    let path = repo().join("assets/live/plugins/family.js");
    let js = fs::read_to_string(&path).unwrap();
    assert_eq!(script_faults(&js), Vec::<String>::new());
    // Control: each fault planted in a copy is found.
    let planted = format!(
        "{js}\nx.innerHTML = y;\neval (s);\nnew Function(\"a\");\nfetch(\"https://cdn.example.org/a.js\");\n"
    );
    assert_eq!(
        script_faults(&planted).len(),
        4,
        "{:?}",
        script_faults(&planted)
    );
    match Command::new("node").arg("--check").arg(&path).output() {
        Ok(out) => {
            assert!(
                out.status.success(),
                "node --check family.js: {}",
                String::from_utf8_lossy(&out.stderr)
            );
            // Control: node --check refuses a copy cut in the middle of a function.
            let w = Work::new("switch-node");
            let cut = w.dir.join("family.js");
            fs::write(&cut, format!("{js}\nfunction planted() {{\n")).unwrap();
            let out = Command::new("node")
                .arg("--check")
                .arg(&cut)
                .output()
                .unwrap();
            assert!(
                !out.status.success(),
                "control: node --check read a cut copy"
            );
        }
        Err(_) => {
            assert!(
                std::env::var_os("CI").is_none(),
                "node is not on the PATH, and CI needs it for this check"
            );
            eprintln!("note: node is not on the PATH, so family.js is not syntax-checked");
        }
    }
}
