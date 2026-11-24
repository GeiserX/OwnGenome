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
        || matches!(field(summary, "planted_acts"), Some(Value::Int(n)) if *n > 0)
}

/// The sample of a summary `og live serve` collected itself from the run
/// folder (`og summary collect` names it; the planted summary.json files of
/// these tests never do). The stubs give child_a's collected summary two
/// results to act on and every other one none, so a test can tell which
/// summary a count came from.
fn collected(summary: &Value) -> Option<&str> {
    match field(summary, "sample") {
        Some(Value::Str(s)) => Some(s.as_str()),
        _ => None,
    }
}

/// The stub results: a finding to act on when the summary plants one (or
/// is child_a's collected summary), otherwise one finding for information.
fn build(name: &str, summary: &Value, _: Option<&BTreeSet<String>>) -> Result<Value, String> {
    let act = match collected(summary) {
        Some(s) => s == "child_a",
        None => planted_act(summary),
    };
    let (id, level) = if act {
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

/// The stub overview: the summary's planted attention count, top and lead,
/// and its results to act on (`planted_acts`, or one for `planted_act`).
/// A collected summary: child_a's has two results to act on and 9 cards,
/// any other none, 7 cards, top "fresh top" and lead "fresh lead".
fn overview(summary: &Value) -> Result<Value, String> {
    if let Some(sample) = collected(summary) {
        return loads(if sample == "child_a" {
            r#"{"attention": 9, "levels": {"act": 2}, "top": ["fresh act"], "lead": "fresh act lead"}"#
        } else {
            r#"{"attention": 7, "levels": {"act": 0}, "top": ["fresh top"], "lead": "fresh lead"}"#
        })
        .map_err(|e| e.to_string());
    }
    let n = match field(summary, "planted_attention") {
        Some(Value::Int(n)) => *n,
        _ => 0,
    };
    let act = match field(summary, "planted_acts") {
        Some(Value::Int(n)) => *n,
        _ => i128::from(planted_act(summary)),
    };
    let mut o = format!(r#"{{"attention": {n}, "levels": {{"act": {act}}}"#);
    for k in ["top", "lead"] {
        if let Some(v) = field(summary, &format!("planted_{k}")) {
            o.push_str(&format!(r#", "{k}": {}"#, og_live::state::dumps(v)));
        }
    }
    o.push('}');
    loads(&o).map_err(|e| e.to_string())
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

    /// `name`'s summary.json made older than the end of its step 06, so
    /// results.json reads the summary the site collects from the folder.
    fn stale(&self, name: &str) {
        let p = self.w.root().join(name).join("summary.json");
        fs::File::options()
            .write(true)
            .open(&p)
            .unwrap()
            .set_modified(std::time::UNIX_EPOCH + std::time::Duration::from_secs(946_684_800))
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

// --- the home page and family.json over the summary results.json reads ----------------------
// (moved from tests/test_owngenome_live_family.py's gate_server_checks)

/// `/runs.json`'s entry of `name`.
fn index_entry(port: u16, name: &str) -> Value {
    let r = get(port, "/runs.json");
    assert_eq!(r.status, 200, "{}", r.text());
    let v = r.json();
    let Some(Value::Array(runs)) = field(&v, "runs") else {
        panic!("/runs.json has no runs: {}", r.text())
    };
    runs.iter()
        .find(|e| matches!(field(e, "name"), Some(Value::Str(n)) if n == name))
        .cloned()
        .unwrap_or_else(|| panic!("no {name} in /runs.json"))
}

fn acts(entry: &Value) -> Option<Value> {
    field(entry, "levels")
        .and_then(|l| field(l, "act"))
        .cloned()
}

/// family.json's entry of member run `name`.
fn family_entry(port: u16, name: &str) -> Value {
    let r = get(port, "/runs/family_a/family.json");
    assert_eq!(r.status, 200, "{}", r.text());
    let v = r.json();
    field(&v, "runs")
        .and_then(|runs| field(runs, name))
        .cloned()
        .unwrap_or_else(|| panic!("no {name} in family.json: {}", r.text()))
}

#[test]
fn home_counts_a_childs_results_to_act_on_from_the_summary_results_json_reads() {
    // child_a's summary.json has no result to act on and is older than its
    // step 06; the summary the site collects (results.json's) has two.
    let f = Fam::new("switch-stale-child");
    f.run("child_a", r#"{"sections": {}, "planted_attention": 5}"#);
    f.stale("child_a");
    let srv = start(f.config());
    let e = index_entry(srv.port(), "child_a");
    assert_eq!(
        acts(&e),
        Some(Value::Int(2)),
        "home counts the stale summary.json: {e:?}"
    );
    assert_eq!(
        field(&family_entry(srv.port(), "child_a"), "attention"),
        Some(&Value::Null),
        "family.json gives a count for a child whose results hold one to act on"
    );
    drop(srv);
    // Control: the same summary.json, current: the home has no result to act on and the count.
    let f = Fam::new("switch-fresh-child");
    f.run("child_a", r#"{"sections": {}, "planted_attention": 5}"#);
    let srv = start(f.config());
    let e = index_entry(srv.port(), "child_a");
    assert_eq!(acts(&e), Some(Value::Int(0)), "control: {e:?}");
    assert_eq!(
        field(&family_entry(srv.port(), "child_a"), "attention"),
        Some(&Value::Int(5))
    );
}

#[test]
fn a_child_with_no_result_to_act_on_keeps_the_numbers_of_its_summary_json() {
    // father_a is a child (of gf_a and gm_a). Its summary.json is older than
    // step 06 and has no result to act on, nor has the collected summary,
    // whose numbers differ: the home keeps summary.json's, as no gate applies.
    let f = Fam::new("switch-stale-noact");
    f.run(
        "father_a",
        r#"{"sections": {}, "planted_attention": 3, "planted_top": ["stale top"], "planted_lead": "stale lead"}"#,
    );
    f.stale("father_a");
    // A finished run: while a run goes, the home counts from the summary its
    // page last read (a rule of its own, not the gate's).
    let status = f.w.root().join("father_a/logs/run_status.tsv");
    let mut text = fs::read_to_string(&status).unwrap();
    text.push_str("meta\tfinished_utc\t2026-01-02 11:00:00 UTC\n");
    fs::write(&status, text).unwrap();
    let srv = start(f.config());
    let want = |e: &Value| ["attention", "levels", "top", "lead"].map(|k| field(e, k).cloned());
    let stale = [
        Some(Value::Int(3)),
        loads(r#"{"act": 0}"#).ok(),
        loads(r#"["stale top"]"#).ok(),
        Some(Value::Str("stale lead".into())),
    ];
    for round in ["first", "second"] {
        // The second request finds the collected summary cached by the first.
        let e = index_entry(srv.port(), "father_a");
        assert_eq!(want(&e), stale, "{round} /runs.json: {e:?}");
    }
    assert_eq!(
        field(&family_entry(srv.port(), "father_a"), "attention"),
        Some(&Value::Int(3))
    );
}

#[test]
fn a_childs_family_entry_is_the_same_whichever_card_holds_its_results_to_act_on() {
    // Two results to act on in one card (0 more cards worth a look) or in two
    // (one more): the counts differ, father_a's family.json entry must not.
    let mut seen = Vec::new();
    for (k, attention) in [("one-card", 0), ("two-cards", 1)] {
        let f = Fam::new(&format!("switch-cards-{k}"));
        f.run(
            "father_a",
            &format!(
                r#"{{"sections": {{}}, "planted_acts": 2, "planted_attention": {attention}}}"#
            ),
        );
        let srv = start(f.config());
        seen.push(family_entry(srv.port(), "father_a"));
    }
    assert_eq!(
        seen[0], seen[1],
        "father_a's entry tells the two reports apart"
    );
    assert_eq!(field(&seen[0], "attention"), Some(&Value::Null));
    // Control: an adult's entries (gf_a's place) do differ, so the comparison can fail.
    let mut adult = Vec::new();
    for (k, attention) in [("one-card", 0), ("two-cards", 1)] {
        let f = Fam::new(&format!("switch-cards-adult-{k}"));
        f.run(
            "gf_a",
            &format!(
                r#"{{"sections": {{}}, "planted_acts": 2, "planted_attention": {attention}}}"#
            ),
        );
        let srv = start(f.config());
        adult.push(family_entry(srv.port(), "gf_a"));
    }
    assert_ne!(
        adult[0], adult[1],
        "control: an adult's two entries are alike"
    );
}

// --- the run index: a run that went quiet, and one that comes back ---------------------------
// (moved from tests/test_run_state_stopped.py's site checks)

const DAY: u64 = 86_400;

/// A run folder `name` under `root` with two planned steps and one trace
/// row in the shared `root/pipeline_info`, every file written `age` seconds ago.
fn aged_run(root: &Path, name: &str, age: u64) -> std::time::SystemTime {
    let t = std::time::SystemTime::now() - std::time::Duration::from_secs(age);
    let secs = t.duration_since(std::time::UNIX_EPOCH).unwrap().as_secs();
    let status = root.join(name).join("logs/run_status.tsv");
    write(
        &status,
        &format!(
            "meta\tstarted_epoch\t{}\nstep\t16\tplanned\nstep\t03\tplanned\n",
            secs - 3600
        ),
    );
    let when = jiff::Timestamp::from_second(secs as i64 - 600)
        .unwrap()
        .strftime("%Y-%m-%d %H:%M:%S.000")
        .to_string();
    let trace = root.join("pipeline_info").join(format!("trace_{name}.txt"));
    write(
        &trace,
        &format!(
            "task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\trealtime\n\
             1\tab/cdef12\t1\tFASTP ({name})\tCOMPLETED\t0\t{when}\t5m\t5m\n"
        ),
    );
    for p in [status, trace] {
        fs::File::options()
            .write(true)
            .open(&p)
            .unwrap()
            .set_modified(t)
            .unwrap();
    }
    t
}

fn run_state_of(e: &Value) -> (Option<Value>, Option<Value>) {
    (
        field(e, "state").cloned(),
        field(e, "last_heard_utc").cloned(),
    )
}

#[test]
fn the_index_reads_a_quiet_run_stopped_and_a_message_brings_it_back() {
    let w = Work::new("switch-stopped");
    let old = aged_run(&w.root(), "old1", 3 * DAY);
    aged_run(&w.root(), "new1", 60);
    let srv = start(config(&w));
    let p = srv.port();
    let (state, heard) = run_state_of(&index_entry(p, "old1"));
    // The day of the instant the fixture used, so a UTC midnight in between changes nothing.
    let three_days = jiff::Timestamp::try_from(old)
        .unwrap()
        .strftime("%Y-%m-%d")
        .to_string();
    assert_eq!(state, Some(Value::Str("stopped".into())));
    assert!(
        matches!(&heard, Some(Value::Str(h)) if h.starts_with(&three_days)),
        "old1 was last heard from three days ago: {heard:?}"
    );
    // new1 writes its trace into the same pipeline_info: that keeps new1
    // running, never old1.
    assert_eq!(
        run_state_of(&index_entry(p, "new1")).0,
        Some(Value::Str("running".into()))
    );
    // A new message of old1: running again.
    let now = jiff::Timestamp::now()
        .strftime("%Y-%m-%dT%H:%M:%SZ")
        .to_string();
    store(
        &w,
        "old1",
        &[&format!(
            r#"{{"runName": "r", "event": "process_started", "utcTime": "{now}", "trace": {{"task_id": 2,
               "name": "BWA_MEM2 (old1)", "process": "BWA_MEM2", "tag": "old1", "status": "RUNNING"}}}}"#
        )
        .replace('\n', " ")],
    );
    assert_eq!(
        run_state_of(&index_entry(p, "old1")).0,
        Some(Value::Str("running".into())),
        "a stopped run that gets a message"
    );
}

// --- the run state opens no result file -------------------------------------------------------
// (moved from tests/test_run_state.py check 9)

/// A FIFO at `p`: opening it to read blocks until a writer comes, so a
/// reader that opens it never returns.
fn fifo(p: &Path) {
    fs::create_dir_all(p.parent().unwrap()).unwrap();
    let st = Command::new("mkfifo").arg(p).status().expect("mkfifo runs");
    assert!(st.success(), "mkfifo {}", p.display());
}

/// Whether `f` returns within 10 s.
fn returns<T: Send + 'static>(f: impl FnOnce() -> T + Send + 'static) -> Option<T> {
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = tx.send(f());
    });
    rx.recv_timeout(std::time::Duration::from_secs(10)).ok()
}

#[test]
fn the_run_state_reads_its_files_and_opens_no_result_file() {
    let w = Work::new("switch-no-result-open");
    let d = w.root().join("demo1");
    write(
        &d.join("logs/run_status.tsv"),
        "meta\tstarted_epoch\t1767347940\nstep\t27\tok\nstep\t03\tplanned\n",
    );
    write(&d.join("summary.json"), r#"{"sections": {}}"#);
    let (vcf, cpic) = (d.join("vcf/demo1.vcf.gz"), d.join("cpic/demo1_cpic.tsv"));
    fifo(&vcf);
    fifo(&cpic);
    // Control: a reader that opens a planted result file does not return.
    let probe = vcf.clone();
    assert!(
        returns(move || fs::File::open(probe).is_ok()).is_none(),
        "control: opening the FIFO returned, so the check cannot tell a reader that opens it"
    );
    // Let the control's open return, so its thread ends.
    let _ = fs::OpenOptions::new().write(true).open(&vcf);
    let steps = og_runlog::steps::read_steps(&repo().join("assets/steps.tsv")).unwrap();
    let dd = d.clone();
    let got = returns(move || {
        let src = og_runlog::state::Sources {
            traces: &[],
            events: None,
            stale: true,
        };
        let files = og_runlog::state::RunFiles::read(&dd, &src).map_err(|e| e.to_string())?;
        let st = files.state("demo1", &steps, 1_767_400_000.0, &jiff::tz::TimeZone::UTC);
        Ok::<String, String>(og_live::state::dumps(&st.to_json()))
    });
    let state = got.expect("the run state opened a planted result file (it never returned)");
    let state = state.expect("the run state reads its files");
    assert!(
        state.contains(r#""step": "27""#) || state.contains(r#""step":"27""#),
        "{state}"
    );
    // The site's state.json, the same way.
    let srv = start(config(&w));
    let p = srv.port();
    let r = returns(move || get(p, "/runs/demo1/state.json").status);
    assert_eq!(
        r,
        Some(200),
        "state.json of a run with planted result FIFOs"
    );
}

// --- the ingest port takes the token as HTTP Basic -------------------------------------------

#[test]
fn the_ingest_port_takes_the_token_as_the_http_basic_password_of_og() {
    let w = Work::new("switch-basic");
    let mut c = config(&w);
    c.ingest = Some(("127.0.0.1".into(), 0));
    let srv = start(c);
    let ip = srv.ingest_port().expect("an ingest listener");
    let token = fs::read_to_string(w.state().join("ingest-token")).unwrap();
    let token = token.trim().to_owned();
    // Standard base64 with padding, as HTTP Basic sends the credentials.
    let b64 = |t: &str| {
        const ABC: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        let mut out = String::new();
        for c in t.as_bytes().chunks(3) {
            let n = c
                .iter()
                .enumerate()
                .fold(0u32, |n, (i, b)| n | (u32::from(*b) << (16 - 8 * i)));
            for i in 0..4 {
                out.push(if i <= c.len() {
                    char::from(ABC[(n >> (18 - 6 * i)) as usize & 63])
                } else {
                    '='
                });
            }
        }
        out
    };
    assert_eq!(b64("og:abc"), "b2c6YWJj");
    assert_eq!(b64("og:ab"), "b2c6YWI=");
    let body = r#"{"event":"pgp_step","step":"14","state":"started"}"#;
    let post = |path: &str, auth: Option<String>| {
        let mut head = format!(
            "POST {path} HTTP/1.1\r\nHost: 127.0.0.1\r\nContent-Type: application/json\r\nContent-Length: {}\r\n",
            body.len()
        );
        if let Some(a) = auth {
            head.push_str(&format!("Authorization: {a}\r\n"));
        }
        head.push_str("Connection: close\r\n\r\n");
        raw(ip, &head, body.as_bytes()).status
    };
    let stored = || {
        fs::read_to_string(events_path(&w, "basic1"))
            .map(|t| t.lines().count())
            .unwrap_or(0)
    };
    assert_eq!(post("/ingest/basic1", None), 403, "no token");
    assert_eq!(
        post(
            "/ingest/basic1",
            Some(format!("Basic {}", b64(&format!("og:{}x", &token[1..]))))
        ),
        403,
        "another token"
    );
    assert_eq!(
        post(
            "/ingest/basic1",
            Some(format!("Basic {}", b64(&format!("me:{token}"))))
        ),
        403,
        "another user"
    );
    assert_eq!(
        post("/ingest/basic1", Some(format!("Bearer {token}"))),
        403,
        "a bearer token"
    );
    assert_eq!(stored(), 0, "a refused post was stored");
    assert_eq!(
        post(
            "/ingest/basic1",
            Some(format!("Basic {}", b64(&format!("og:{token}"))))
        ),
        200
    );
    assert_eq!(stored(), 1);
    // The token as the first part of the path still works.
    assert_eq!(post(&format!("/{token}/ingest/basic1"), None), 200);
    assert_eq!(stored(), 2);
}

// --- trace files: read on each call, with the right rows ------------------------------------
// (moved from tests/test_run_state_trace_cache.py; the cache itself was not
// ported: in a release build one /runs.json over eight runs, each reading
// 16,000 trace rows, took 53 ms on a test machine and one state.json 21 ms,
// under the 100 ms that would call for it; see
// `measure_the_run_index_over_big_traces`)

const TRACE_HEAD: &str = "task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\trealtime\t%cpu\tpeak_rss\tpeak_vmem\trchar\twchar\n";

fn trace_row(id: usize, sample: &str, status: &str, duration: &str) -> String {
    format!(
        "{id}\tab/{id:06}\t{}\tPGX:CLINVAR_SCREEN ({sample})\t{status}\t0\t2026-01-02 10:{:02}:00.000\t{duration}\t4s\t50.0%\t1 MB\t1 MB\t1 MB\t1 MB\n",
        1000 + id,
        id % 60
    )
}

fn trace_ids(p: &Path, sample: &str) -> Vec<(String, String, bool, bool)> {
    og_runlog::tasks::trace_tasks(&[p], sample, &jiff::tz::TimeZone::UTC)
        .into_iter()
        .map(|t| (t.id, t.status, t.start.is_some(), t.end.is_some()))
        .collect()
}

#[test]
fn a_trace_gives_each_sample_its_own_rows_and_a_bad_duration_costs_only_its_end() {
    let w = Work::new("switch-trace-rows");
    let p = w.dir.join("trace.txt");
    let text = format!(
        "{TRACE_HEAD}{}{}{}{}{}",
        trace_row(1, "demo1", "COMPLETED", "5s"),
        trace_row(2, "demo2", "COMPLETED", "5s"),
        trace_row(3, "demo1", "COMPLETED", "1.2.3s"),
        trace_row(4, "demo1", "FAILED", ".s"),
        trace_row(5, "demo1", "COMPLETED", "1m 5s"),
    );
    write(&p, &text);
    let one = trace_ids(&p, "demo1");
    let two = trace_ids(&p, "demo2");
    let row = |id: &str, status: &str, end: bool| (id.to_owned(), status.to_owned(), true, end);
    assert_eq!(
        one,
        [
            row("1", "COMPLETED", true),
            row("3", "COMPLETED", false),
            row("4", "FAILED", false),
            row("5", "COMPLETED", true),
        ],
        "demo1's rows: a duration that is no number costs that row its end only"
    );
    assert_eq!(
        two,
        [row("2", "COMPLETED", true)],
        "demo2 gets only its own row"
    );
}

#[test]
fn a_trace_that_grew_or_was_replaced_is_read_again() {
    let w = Work::new("switch-trace-change");
    let p = w.dir.join("trace.txt");
    write(
        &p,
        &format!("{TRACE_HEAD}{}", trace_row(1, "demo1", "COMPLETED", "5s")),
    );
    assert_eq!(trace_ids(&p, "demo1").len(), 1);
    // A task row appended: the next read has it.
    let mut f = fs::OpenOptions::new().append(true).open(&p).unwrap();
    std::io::Write::write_all(&mut f, trace_row(99, "demo1", "COMPLETED", "5s").as_bytes())
        .unwrap();
    drop(f);
    let grown = trace_ids(&p, "demo1");
    assert_eq!(grown.last().map(|r| r.0.as_str()), Some("99"), "{grown:?}");
    // Replaced by a file of the same size and modification time: the next read has its rows.
    let before = fs::metadata(&p).unwrap();
    let text = fs::read_to_string(&p).unwrap();
    let swapped = text.replacen("\tCOMPLETED\t0\t", "\tFAILED...\t1\t", 1);
    assert_eq!(swapped.len(), text.len());
    let new = w.dir.join("trace.new");
    write(&new, &swapped);
    fs::File::options()
        .write(true)
        .open(&new)
        .unwrap()
        .set_modified(before.modified().unwrap())
        .unwrap();
    fs::rename(&new, &p).unwrap();
    let after = fs::metadata(&p).unwrap();
    assert_eq!(
        (after.len(), after.modified().unwrap()),
        (before.len(), before.modified().unwrap())
    );
    let statuses: Vec<String> = trace_ids(&p, "demo1").into_iter().map(|r| r.1).collect();
    assert_eq!(
        statuses,
        ["FAILED...", "COMPLETED"],
        "the replaced file's rows"
    );
}

/// The cost of the run index over big traces: eight runs whose traces sit in
/// one pipeline_info, 2,000 task rows each, so each run's state reads all
/// 16,000 rows. Prints the time per /runs.json; run it with --ignored (and
/// --release for the server's build).
#[test]
#[ignore = "a measurement: run with --ignored"]
fn measure_the_run_index_over_big_traces() {
    let w = Work::new("switch-trace-measure");
    let samples: Vec<String> = (1..=8).map(|i| format!("s{i}")).collect();
    for s in &samples {
        write(
            &w.root().join(s).join("logs/run_status.tsv"),
            "meta\tstarted_epoch\t1767347940\nstep\t06\tplanned\nstep\t03\tplanned\n",
        );
        let mut t = String::from(TRACE_HEAD);
        for i in 0..2000 {
            t.push_str(&trace_row(
                i,
                s,
                if i % 7 == 0 { "FAILED" } else { "COMPLETED" },
                "1m 5s",
            ));
        }
        write(
            &w.root()
                .join("pipeline_info")
                .join(format!("trace_{s}.txt")),
            &t,
        );
    }
    let rows: usize = og_runlog::tasks::trace_tasks(
        &samples
            .iter()
            .map(|s| {
                w.root()
                    .join("pipeline_info")
                    .join(format!("trace_{s}.txt"))
            })
            .collect::<Vec<_>>(),
        "s1",
        &jiff::tz::TimeZone::UTC,
    )
    .len();
    assert_eq!(rows, 2000);
    let srv = start(config(&w));
    let p = srv.port();
    assert_eq!(get(p, "/runs.json").status, 200);
    let n = 10;
    let t0 = std::time::Instant::now();
    for _ in 0..n {
        assert_eq!(get(p, "/runs.json").status, 200);
    }
    let per = t0.elapsed() / n;
    let t1 = std::time::Instant::now();
    for _ in 0..n {
        assert_eq!(get(p, "/runs/s1/state.json").status, 200);
    }
    let one = t1.elapsed() / n;
    println!(
        "MEASURE /runs.json over 8 runs x 16,000 trace rows: {per:?} per call; state.json of one run: {one:?}"
    );
}
