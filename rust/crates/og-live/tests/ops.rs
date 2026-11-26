//! og-live's ops card and `/metrics` text on the made-up run of
//! `tests/fixtures/live/ops`, with a fake disk and a fake host; then the same
//! inputs through the Python plugin it replaces (`bin/owngenome_live_ops.py`,
//! while the Python site still runs it), which must give the same card by
//! meaning and the same metrics text byte for byte. A second parity case
//! runs on the traces and messages of `tests/fixtures/run_eta` and
//! `tests/fixtures/run_headroom`.
//!
//! With `--features fault-inject` each `planted_*` test turns on one planted
//! mistake and runs the behaviour checks: every one of them must fail
//! (Guard's rust job requires it). A check that cannot fail is no check.

use jiff::tz::TimeZone;
use og_live::ops::{
    self, Asker, Ctx, METRICS, OpsConfig, OpsState, RunOps, current_run, metrics_allowed,
    metrics_text, run_ops,
};
use og_runlog::host::{FsStat, Host};
use og_runlog::state::{Inputs, RunState, State, build};
use og_runlog::steps::{Step, read_steps};
use pycompat::json::{Value, loads};
use std::collections::{BTreeMap, BTreeSet};
use std::io;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Duration, SystemTime};

const SAMPLE: &str = "OPS-001";
/// 2026-01-02T10:00:00Z, the fixture's start.
const T0: f64 = 1_767_348_000.0;
const G: u128 = 1 << 30;
const UTC: TimeZone = TimeZone::UTC;
/// The last commit whose `bin/` still held the Python plugin, for the parity
/// cases once the Python site is gone.
const ORACLE: &str = "1ec630690e63796b4410d23ec65768200adb8dc0";

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn fix(name: &str) -> PathBuf {
    repo().join("tests/fixtures/live/ops").join(name)
}

fn steps() -> &'static [Step] {
    static S: OnceLock<Vec<Step>> = OnceLock::new();
    S.get_or_init(|| read_steps(&repo().join("assets/steps.tsv")).expect("assets/steps.tsv reads"))
}

/// A new empty folder of this test run.
fn tmp() -> PathBuf {
    static N: AtomicUsize = AtomicUsize::new(0);
    let d = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!(
        "ops-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn set_mtime(p: &Path, t: f64) {
    std::fs::File::options()
        .write(true)
        .open(p)
        .unwrap()
        .set_modified(SystemTime::UNIX_EPOCH + Duration::from_secs_f64(t))
        .unwrap();
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

/// A host whose every file system is `ratio` used; device ids and
/// existence are the real ones.
struct Fake {
    ratio: f64,
    phys: u128,
}

impl Host for Fake {
    fn statvfs(&self, _: &Path) -> io::Result<FsStat> {
        let blocks = 1_000_000u64;
        let free = (blocks as f64 * (1.0 - self.ratio)) as u64;
        Ok(FsStat {
            frsize: 4096,
            bsize: 4096,
            blocks,
            bfree: free,
            bavail: free,
            files: 1000,
            ffree: 900,
        })
    }
    fn dev(&self, p: &Path) -> io::Result<u64> {
        use std::os::unix::fs::MetadataExt;
        Ok(std::fs::metadata(p)?.dev())
    }
    fn exists(&self, p: &Path) -> bool {
        p.exists()
    }
    fn physical_memory(&self) -> Option<u128> {
        Some(self.phys)
    }
}

/// The genome folder of one case.
struct Tree {
    root: PathBuf,
    run_dir: PathBuf,
    events: Vec<Value>,
    events_file: PathBuf,
    work: PathBuf,
    out: PathBuf,
    sample: String,
}

/// GENOME_DIR in a new folder: the history of box `opsbox` and, for a
/// run-all.sh run, the run folder; `extra` messages after the fixture's.
fn make_tree(run_all: bool, extra: &[Value]) -> Tree {
    let tmp = tmp();
    let root = tmp.join("genome");
    let (work, out) = (tmp.join("nfwork"), tmp.join("nfout"));
    let hist = root.join("pipeline_info/step_history");
    for d in [&root, &work, &out, &hist] {
        std::fs::create_dir_all(d).unwrap();
    }
    std::fs::copy(fix("step_history.tsv"), hist.join("opsbox.tsv")).unwrap();
    let text = std::fs::read_to_string(fix("events.jsonl"))
        .unwrap()
        .replace("@WORK@", &s(&work))
        .replace("@OUT@", &s(&out));
    let mut events: Vec<Value> = text
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| loads(l).unwrap())
        .collect();
    events.extend(extra.iter().cloned());
    let run_dir = root.join(SAMPLE);
    if run_all {
        std::fs::create_dir_all(run_dir.join("logs")).unwrap();
        std::fs::create_dir_all(run_dir.join("nextflow/work")).unwrap();
        let status = run_dir.join("logs/run_status.tsv");
        std::fs::copy(fix("run_status.tsv"), &status).unwrap();
        set_mtime(&status, T0);
    }
    let events_file = tmp.join("events.jsonl");
    write_events(&events_file, &events);
    Tree {
        root,
        run_dir,
        events,
        events_file,
        work,
        out,
        sample: SAMPLE.into(),
    }
}

fn write_events(path: &Path, events: &[Value]) {
    let opts = pycompat::json::Dump::default();
    let text: String = events
        .iter()
        .map(|e| pycompat::json::dumps(e, &opts).unwrap() + "\n")
        .collect();
    std::fs::write(path, text).unwrap();
}

/// How one case runs.
#[derive(Clone)]
struct Case {
    now: f64,
    physical: u128,
    disk: f64,
    stopped: bool,
    stall_minutes: f64,
    disk_warn_percent: Option<&'static str>,
}

fn case(now: f64) -> Case {
    Case {
        now,
        physical: 64 * G,
        disk: 0.50,
        stopped: false,
        stall_minutes: 60.0,
        disk_warn_percent: None,
    }
}

/// The warn percents are given (90 unless the case says otherwise), so an
/// `OWNGENOME_*_WARN_PERCENT` of the shell running the tests changes nothing.
fn cfg_for(t: &Tree, c: &Case) -> OpsConfig {
    OpsConfig {
        root: Some(t.root.clone()),
        eta_box: Some("opsbox".into()),
        eta_cores: Some(8),
        disk_warn_percent: Some(c.disk_warn_percent.unwrap_or("90").to_owned()),
        mem_warn_percent: Some("90".into()),
        stall_minutes: Some(c.stall_minutes),
    }
}

/// The state the site builds for the tree at `c.now`.
fn state_of(t: &Tree, c: &Case) -> State<'static> {
    let status = og_runlog::files::read_run_status(&t.run_dir.join("logs/run_status.tsv")).unwrap();
    let cur = current_run(&t.events);
    let mut st = build(
        &t.sample,
        steps(),
        &Inputs {
            status: status.as_ref(),
            event_tasks: og_runlog::tasks::event_tasks(cur, &t.sample),
            events: cur,
            ..Default::default()
        },
        c.now,
        &UTC,
    );
    if c.stopped {
        st.run = RunState::Stopped;
    }
    st
}

fn ops_of(t: &Tree, c: &Case) -> RunOps {
    let host = Fake {
        ratio: c.disk,
        phys: c.physical,
    };
    run_ops(
        Some(&t.run_dir),
        &state_of(t, c),
        &t.events,
        &cfg_for(t, c),
        &Ctx {
            now: c.now,
            steps: steps(),
            host: &host,
            tz: &UTC,
        },
    )
    .expect("the warn percents are good")
}

fn ops_at(c: &Case, run_all: bool, extra: &[Value]) -> RunOps {
    ops_of(&make_tree(run_all, extra), c)
}

/// `{kind: level}` of the alerts.
fn kinds(o: &RunOps) -> BTreeMap<&'static str, &'static str> {
    o.alerts
        .iter()
        .map(|a| (a.kind.as_str(), a.level.as_str()))
        .collect()
}

/// The series (name and labels, no value) of a Prometheus text.
fn series(text: &str) -> BTreeSet<String> {
    text.lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| l.rsplit_once(' ').unwrap().0.to_owned())
        .collect()
}

fn ev(text: &str) -> Value {
    loads(text).unwrap()
}

fn finished(code: i64) -> Value {
    ev(&format!(
        r#"{{"event":"pgp_run_finished","exit":{code},"utcTime":"2026-01-02T10:06:00Z"}}"#
    ))
}

fn failed_task() -> Value {
    ev(&format!(
        r#"{{"runName":"ops_demo","event":"process_completed","utcTime":"2026-01-02T10:06:00Z","trace":{{"task_id":4,"hash":"ab/000004","name":"CLINICAL:ROH ({SAMPLE})","process":"CLINICAL:ROH","tag":"{SAMPLE}","status":"FAILED","exit":1,"attempt":1,"complete":{}}}}}"#,
        (T0 as i64 + 360) * 1000
    ))
}

// --- the behaviour checks ------------------------------------------------------------------

#[test]
fn half_a_minute_into_step_11_not_late_not_stalled_and_both_ends() {
    let o = ops_at(&case(T0 + 330.0), true, &[]);
    let k = kinds(&o);
    let h = o.hours_left.as_ref().expect("the run is going");
    assert!(
        !k.contains_key("late")
            && !k.contains_key("stalled")
            && !o.run.stalled
            && matches!(h.high_s, Some(hi) if h.low_s <= hi),
        "{o:?}"
    );
}

#[test]
fn two_hours_later_and_quiet_step_11_is_late_and_the_run_stalled() {
    let o = ops_at(&case(T0 + 302.0 + 7200.0), true, &[]);
    let h = o.hours_left.as_ref().expect("the run is going");
    let k = kinds(&o);
    assert!(
        h.late_steps == ["11"]
            && h.high_s.is_none()
            && k.get("late") == Some(&"warning")
            && k.get("stalled") == Some(&"warning"),
        "{:?}",
        o.alerts
    );
}

#[test]
fn the_same_quiet_run_read_stopped_is_still_stalled() {
    let mut c = case(T0 + 302.0 + 7200.0);
    c.stopped = true;
    let o = ops_at(&c, true, &[]);
    assert!(
        o.run.stalled && kinds(&o).get("stalled") == Some(&"warning"),
        "{:?} {:?}",
        o.run,
        o.alerts
    );
}

#[test]
fn quiet_past_stall_minutes_but_step_11_inside_its_usual_time_is_not_stalled() {
    let mut c = case(T0 + 302.0 + 90.0);
    c.stall_minutes = 1.0;
    let o = ops_at(&c, true, &[]);
    let k = kinds(&o);
    assert!(
        o.run.seconds_since_activity.is_some_and(|q| q >= 60.0)
            && !k.contains_key("stalled")
            && !k.contains_key("late"),
        "{:?} {:?}",
        o.run,
        o.alerts
    );
}

#[test]
fn stall_minutes_300_keeps_a_2_hour_quiet_run_from_being_stalled() {
    let mut c = case(T0 + 302.0 + 7200.0);
    c.stall_minutes = 300.0;
    let o = ops_at(&c, true, &[]);
    let k = kinds(&o);
    assert!(
        !k.contains_key("stalled") && k.contains_key("late"),
        "{:?}",
        o.alerts
    );
}

#[test]
fn a_failed_task_fails_the_run_exit_1_fails_it_exit_0_does_not() {
    let c = case(T0 + 400.0);
    let f1 = ops_at(&c, true, &[failed_task()]);
    let f2 = ops_at(&c, true, &[finished(1)]);
    let f0 = ops_at(&c, true, &[finished(0)]);
    assert!(
        f1.run.failed
            && kinds(&f1).get("failed") == Some(&"critical")
            && f2.run.failed
            && f2.run.exit == Value::Int(1)
            && !f0.run.failed
            && !kinds(&f0).contains_key("failed"),
        "{:?} {:?} {:?}",
        f1.run,
        f2.run,
        f0.run
    );
}

#[test]
fn one_retry_after_exit_137_and_7_6_of_8_gb_asked_is_a_hint() {
    let o = ops_at(&case(T0 + 330.0), true, &[]);
    let m = &o.memory;
    assert!(
        m.retries == 1 && m.near_request == 1 && kinds(&o).get("memory") == Some(&"warning"),
        "{m:?} {:?}",
        o.alerts
    );
}

#[test]
fn max_memory_30_gb_on_a_16_gb_host_warns_on_64_it_does_not() {
    let big = ops_at(&case(T0 + 330.0), true, &[]);
    let mut c = case(T0 + 330.0);
    c.physical = 16 * G;
    let small = ops_at(&c, true, &[]);
    let text = small
        .alerts
        .iter()
        .find(|a| a.kind.as_str() == "memory")
        .map(|a| a.text.clone())
        .unwrap_or_default();
    assert!(
        small.memory.max_memory_above_ram
            && !big.memory.max_memory_above_ram
            && text.contains("max_memory"),
        "{:?}",
        small.memory
    );
}

#[test]
fn a_disk_95_used_warns_at_90_not_at_disk_warn_percent_96() {
    let early = ops_at(&case(T0 + 330.0), true, &[]);
    let mut c = case(T0 + 330.0);
    c.disk = 0.95;
    let full = ops_at(&c, true, &[]);
    c.disk_warn_percent = Some("96");
    let full96 = ops_at(&c, true, &[]);
    assert!(
        kinds(&full).get("disk") == Some(&"warning")
            && !kinds(&full96).contains_key("disk")
            && !kinds(&early).contains_key("disk")
            && full
                .disk
                .first()
                .and_then(|d| d.used_ratio)
                .is_some_and(|r| (r - 0.95).abs() < 0.001),
        "{:?} {:?}",
        full.disk,
        full96.alerts
    );
}

#[test]
fn a_bare_nextflow_run_gives_the_same_series_as_the_run_all_run() {
    let early = ops_at(&case(T0 + 330.0), true, &[]);
    let bare = ops_at(&case(T0 + 330.0), false, &[]);
    let a = series(&metrics_text(&[(SAMPLE, &early)], Some(0.9), None, 0));
    let b = series(&metrics_text(&[(SAMPLE, &bare)], Some(0.9), None, 0));
    assert!(
        a == b && a.len() > 15,
        "{:?}",
        a.symmetric_difference(&b).collect::<Vec<_>>()
    );
}

/// Checks 1 to 9: what each planted mistake must break.
#[cfg_attr(not(feature = "fault-inject"), allow(dead_code))]
fn behaviour() {
    half_a_minute_into_step_11_not_late_not_stalled_and_both_ends();
    two_hours_later_and_quiet_step_11_is_late_and_the_run_stalled();
    the_same_quiet_run_read_stopped_is_still_stalled();
    quiet_past_stall_minutes_but_step_11_inside_its_usual_time_is_not_stalled();
    stall_minutes_300_keeps_a_2_hour_quiet_run_from_being_stalled();
    a_failed_task_fails_the_run_exit_1_fails_it_exit_0_does_not();
    one_retry_after_exit_137_and_7_6_of_8_gb_asked_is_a_hint();
    max_memory_30_gb_on_a_16_gb_host_warns_on_64_it_does_not();
    a_disk_95_used_warns_at_90_not_at_disk_warn_percent_96();
    a_bare_nextflow_run_gives_the_same_series_as_the_run_all_run();
}

// --- no path, two labels, a valid text ------------------------------------------------------

fn strings(v: &Value, out: &mut Vec<String>) {
    match v {
        Value::Str(x) => out.push(x.clone()),
        Value::Array(a) => a.iter().for_each(|x| strings(x, out)),
        Value::Object(o) => o.values().for_each(|x| strings(x, out)),
        _ => {}
    }
}

/// A sample line: `pgp_name{labels} value` or `pgp_name value`.
fn valid_sample(line: &str) -> bool {
    let Some((head, value)) = line.rsplit_once(' ') else {
        return false;
    };
    let name = head.split_once('{').map_or(head, |(n, _)| n);
    let labels_ok = head
        .split_once('{')
        .is_none_or(|(_, rest)| rest.ends_with('}') && !rest[..rest.len() - 1].contains('}'));
    name.starts_with("pgp_")
        && name.chars().all(|c| c.is_ascii_lowercase() || c == '_')
        && labels_ok
        && !value.is_empty()
        && value
            .chars()
            .all(|c| c.is_ascii_digit() || ".e+-".contains(c))
}

fn label_names(text: &str) -> BTreeSet<String> {
    let mut out = BTreeSet::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        if let Some((_, rest)) = line.split_once('{')
            && let Some((labels, _)) = rest.split_once('}')
        {
            out.extend(
                labels
                    .split(',')
                    .filter_map(|kv| kv.split_once('='))
                    .map(|(k, _)| k.to_owned()),
            );
        }
    }
    out
}

#[test]
fn no_path_in_the_card_and_the_metrics_carry_only_run_and_fs_in_a_valid_text() {
    let t = make_tree(true, &[]);
    let o = ops_of(&t, &case(T0 + 330.0));
    let mut all = Vec::new();
    strings(&o.to_json(), &mut all);
    let leaked: Vec<&String> = all
        .iter()
        .filter(|x| {
            x.starts_with('/') || [&t.root, &t.work, &t.out].iter().any(|p| x.contains(&s(p)))
        })
        .collect();
    let text = metrics_text(&[(SAMPLE, &o)], Some(0.9), Some(64 * G), 0);
    let bad: Vec<&str> = text
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#') && !valid_sample(l))
        .collect();
    let labels = label_names(&text);
    assert!(
        leaked.is_empty()
            && labels == BTreeSet::from(["fs".to_owned(), "run".to_owned()])
            && bad.is_empty(),
        "{leaked:?} {labels:?} {bad:?}"
    );
    // Control: a path label and a broken line must be seen.
    let planted = format!("{text}pgp_runs{{path=\"/x\"}} 1\npgp_runs one\n");
    assert!(label_names(&planted).contains("path"));
    assert!(
        planted
            .lines()
            .any(|l| !l.starts_with('#') && !valid_sample(l))
    );
}

// --- who may read /metrics ------------------------------------------------------------------

fn who(login: bool, logged_in: bool, auth: Option<&str>) -> Asker<'_> {
    Asker {
        login_required: login,
        logged_in,
        authorization: auth,
    }
}

fn tok(t: Option<&str>) -> impl FnOnce() -> io::Result<Option<String>> + '_ {
    move || Ok(t.map(str::to_owned))
}

fn unreadable() -> io::Result<Option<String>> {
    Err(io::Error::new(io::ErrorKind::NotFound, "no token file"))
}

#[test]
fn metrics_open_without_a_password_and_with_one_a_login_or_the_bearer_token_only() {
    let yes = |w: &Asker<'_>, t: Option<&str>| metrics_allowed(w, tok(t)).unwrap();
    assert!(yes(&who(false, false, None), None));
    assert!(!yes(&who(true, false, None), None));
    assert!(!yes(&who(true, false, Some("Bearer s3cret")), None));
    assert!(yes(
        &who(true, false, Some("Bearer s3cret")),
        Some("s3cret")
    ));
    assert!(yes(
        &who(true, false, Some("bearer  s3cret ")),
        Some("s3cret")
    ));
    assert!(!yes(
        &who(true, false, Some("Bearer wrong")),
        Some("s3cret")
    ));
    assert!(!yes(
        &who(true, false, Some("Bearer s3cre")),
        Some("s3cret")
    ));
    assert!(!yes(
        &who(true, false, Some("Basic s3cret")),
        Some("s3cret")
    ));
    assert!(yes(&who(true, true, None), None));
    // The token file is read only when a token is needed; unreadable is an error.
    assert!(metrics_allowed(&who(true, false, Some("Bearer s3cret")), unreadable).is_err());
    assert!(metrics_allowed(&who(true, true, None), unreadable).unwrap());
    assert!(metrics_allowed(&who(false, false, None), unreadable).unwrap());
}

#[test]
fn the_token_file_is_read_stripped_and_an_empty_one_is_no_token() {
    let d = tmp();
    std::fs::write(d.join("t"), " s3cret\n").unwrap();
    std::fs::write(d.join("e"), "\n").unwrap();
    assert_eq!(
        ops::read_token(Some(&d.join("t"))).unwrap().as_deref(),
        Some("s3cret")
    );
    assert_eq!(ops::read_token(Some(&d.join("e"))).unwrap(), None);
    assert_eq!(ops::read_token(None).unwrap(), None);
    assert!(ops::read_token(Some(&d.join("missing"))).is_err());
}

// --- the docs and the example rules ---------------------------------------------------------

/// The `pgp_*` names the rules' expressions use.
fn rule_metrics(rules: &str) -> BTreeSet<String> {
    let mut out = BTreeSet::new();
    for line in rules.lines() {
        if !(line.trim_start().starts_with("expr:") || line.starts_with(&" ".repeat(10))) {
            continue;
        }
        let code = line.split('#').next().unwrap_or("");
        out.extend(
            code.split(|c: char| !(c.is_ascii_alphanumeric() || c == '_'))
                .filter(|w| w.starts_with("pgp_"))
                .map(str::to_owned),
        );
    }
    out
}

fn missing_from_doc(doc: &str) -> Vec<&'static str> {
    METRICS
        .iter()
        .map(|(n, _, _)| *n)
        .filter(|n| !doc.contains(&format!("`{n}`")))
        .collect()
}

fn unknown(used: &BTreeSet<String>) -> Vec<&String> {
    used.iter()
        .filter(|n| !METRICS.iter().any(|(m, _, _)| m == n))
        .collect()
}

#[test]
fn every_metric_is_in_the_docs_and_every_metric_the_rules_use_exists() {
    let doc = std::fs::read_to_string(repo().join("docs/alerts.md")).unwrap();
    let rules = std::fs::read_to_string(repo().join("assets/monitoring/pgp-alerts.yml")).unwrap();
    let used = rule_metrics(&rules);
    assert!(
        missing_from_doc(&doc).is_empty() && unknown(&used).is_empty() && !used.is_empty(),
        "{:?} {:?}",
        missing_from_doc(&doc),
        unknown(&used)
    );
    // Control: a doc without one name, and a rule with a name /metrics never
    // writes, must each be seen.
    assert_eq!(
        missing_from_doc(&doc.replace("`pgp_run_stalled`", "")),
        ["pgp_run_stalled"]
    );
    let planted = format!("{rules}\n        expr: pgp_run_nonexistent > 0\n");
    assert_eq!(
        unknown(&rule_metrics(&planted)),
        [&"pgp_run_nonexistent".to_owned()]
    );
}

// --- read only, and the state hook ----------------------------------------------------------

fn mtimes(root: &Path) -> BTreeMap<PathBuf, f64> {
    let mut out = BTreeMap::new();
    let mut todo = vec![root.to_path_buf()];
    while let Some(d) = todo.pop() {
        for e in std::fs::read_dir(&d).unwrap() {
            let p = e.unwrap().path();
            if p.is_dir() {
                todo.push(p);
            } else {
                out.insert(p.clone(), pycompat::mtime::mtime(&p).unwrap());
            }
        }
    }
    out
}

#[test]
fn the_step_history_is_read_never_written() {
    let t = make_tree(true, &[]);
    let hist = t.root.join("pipeline_info/step_history/opsbox.tsv");
    std::fs::write(
        t.root.join("pipeline_info/trace_2026-01-02_10-00-00.txt"),
        "task_id\thash\tname\tstatus\texit\n1\tzz/1\tROH (OPS-001)\tCOMPLETED\t0\n",
    )
    .unwrap();
    let before = mtimes(&t.root);
    let text = std::fs::read_to_string(&hist).unwrap();
    ops_of(&t, &case(T0 + 330.0));
    assert_eq!(mtimes(&t.root), before);
    assert_eq!(std::fs::read_to_string(&hist).unwrap(), text);
}

fn hook_of(t: &Tree, cfg: &OpsConfig) -> OpsState {
    let c = case(T0 + 330.0);
    let host = Fake {
        ratio: 0.5,
        phys: 64 * G,
    };
    ops::hook(
        Some(&t.run_dir),
        &state_of(t, &c),
        &t.events,
        cfg,
        &Ctx {
            now: c.now,
            steps: steps(),
            host: &host,
            tz: &UTC,
        },
    )
}

/// The Python plugin turned any exception into an error card; the one left
/// here is a warn percent outside 1 to 100. A history with no `trace_file`
/// column (the Python card's KeyError) now gives the card from what the
/// history holds.
#[test]
fn the_state_hook_gives_the_card_or_one_error_line() {
    let t = make_tree(true, &[]);
    std::fs::write(
        t.root.join("pipeline_info/step_history/opsbox.tsv"),
        "step\tseconds\nROH\t60\n",
    )
    .unwrap();
    let c = case(T0 + 330.0);
    let ok = hook_of(&t, &cfg_for(&t, &c));
    assert!(
        matches!(&ok, OpsState::Ok(o) if o.hours_left.is_some()),
        "{ok:?}"
    );
    let mut bad = cfg_for(&t, &c);
    bad.disk_warn_percent = Some("0".into());
    let err = hook_of(&t, &bad);
    assert_eq!(
        err.to_json(),
        ev(r#"{"error":"warn percent must be a number from 1 to 100, got '0'"}"#)
    );
}

// --- the Python live site, while it still runs the plugin -----------------------------------

/// A free local port.
fn free_port() -> u16 {
    std::net::TcpListener::bind("127.0.0.1:0")
        .and_then(|l| l.local_addr())
        .expect("a local port")
        .port()
}

/// One GET on the local site: status, content type and body; `None` while
/// it does not answer.
fn get(port: u16, path: &str) -> Option<(u16, String, String)> {
    use std::io::{Read, Write};
    let mut c = std::net::TcpStream::connect(("127.0.0.1", port)).ok()?;
    c.set_read_timeout(Some(Duration::from_secs(10))).ok()?;
    write!(c, "GET {path} HTTP/1.0\r\nHost: 127.0.0.1\r\n\r\n").ok()?;
    let mut raw = Vec::new();
    c.read_to_end(&mut raw).ok()?;
    let raw = String::from_utf8(raw).ok()?;
    let (head, body) = raw.split_once("\r\n\r\n")?;
    let status = head.split(' ').nth(1)?.parse().ok()?;
    let ctype = head
        .lines()
        .find_map(|l| {
            let (k, v) = l.split_once(':')?;
            k.eq_ignore_ascii_case("content-type")
                .then(|| v.trim().to_owned())
        })
        .unwrap_or_default();
    Some((status, ctype, body.to_owned()))
}

/// The site's process, stopped when the test ends, pass or fail.
struct Site(std::process::Child);

impl Drop for Site {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}

/// What `bin/owngenome_live.py` serves on the ops tree.
struct Served {
    /// `/runs/OPS-001/state.json`.
    state: Value,
    /// `/metrics`: content type and text.
    metrics: (String, String),
    /// `state.json` again, once the history has lost its `trace_file` column.
    broken: Option<Value>,
}

fn key<'a>(v: &'a Value, k: &str) -> Option<&'a Value> {
    match v {
        Value::Object(m) => m.get(k),
        _ => None,
    }
}

/// Starts the Python site on a new ops tree (`plugin_dir` in place of
/// `bin/` when given) and reads it.
fn python_site(plugin_dir: Option<&Path>) -> Result<Served, String> {
    let t = make_tree(true, &[]);
    let base = t.root.parent().unwrap().to_path_buf();
    let state_dir = base.join("state");
    let runs = state_dir.join("runs").join(SAMPLE);
    std::fs::create_dir_all(&runs).unwrap();
    write_events(&runs.join("events.jsonl"), &t.events);
    let (port, iport) = (free_port(), free_port());
    let exe = std::env::var("OG_PYTHON").unwrap_or_else(|_| "python3".into());
    let err_log = base.join("site.stderr");
    let mut cmd = Command::new(&exe);
    cmd.arg(repo().join("bin/owngenome_live.py"))
        .arg("--root")
        .arg(&t.root)
        .arg("--state-dir")
        .arg(&state_dir)
        .args([
            "--port",
            &port.to_string(),
            "--ingest-port",
            &iport.to_string(),
        ])
        .args(["--quiet", "--stall-minutes", "30"])
        .args(["--disk-warn-percent", "90", "--mem-warn-percent", "90"])
        .env("PGP_BOX", "opsbox")
        .env_remove("OWNGENOME_LIVE_PASSWORD_HASH")
        .stdout(std::process::Stdio::null())
        .stderr(std::fs::File::create(&err_log).unwrap());
    if let Some(d) = plugin_dir {
        cmd.arg("--plugin-dir").arg(d);
    }
    let mut site = Site(
        cmd.spawn()
            .map_err(|e| format!("cannot run {exe}: {e}; set OG_PYTHON"))?,
    );
    let why = || std::fs::read_to_string(&err_log).unwrap_or_default();
    let state_path = format!("/runs/{SAMPLE}/state.json");
    let read_state = || match get(port, &state_path) {
        Some((200, _, body)) => loads(&body).map_err(|e| format!("state.json: {e:?}")),
        Some((code, _, body)) => Err(format!("state.json answered {code}: {body}")),
        None => Err("state.json: no answer".into()),
    };
    let mut state = Err(String::new());
    for _ in 0..200 {
        if let Some(code) = site.0.try_wait().unwrap() {
            return Err(format!("the site exited ({code}): {}", why()));
        }
        state = read_state();
        if state.is_ok() {
            break;
        }
        std::thread::sleep(Duration::from_millis(100));
    }
    let state = state.map_err(|e| format!("{e}; {}", why()))?;
    let (_, ctype, text) = get(port, "/metrics").ok_or("/metrics: no answer")?;
    std::fs::write(
        t.root.join("pipeline_info/step_history/opsbox.tsv"),
        "step\tseconds\nROH\t60\n",
    )
    .unwrap();
    let broken = read_state().ok();
    Ok(Served {
        state,
        metrics: (ctype, text),
        broken,
    })
}

/// Check 14 and 15 of the Python suite this file replaces: the site loads
/// the plugin (state.json has `plugins.ops`, `/metrics` answers in the
/// Prometheus text format), and an error in the plugin's hook becomes the
/// card's error field while the rest of state.json stays.
fn site_serves_the_plugin(s: &Served) -> Result<(), String> {
    let ops = key(&s.state, "plugins").and_then(|p| key(p, "ops"));
    if ops.and_then(|o| key(o, "alerts")).is_none() {
        return Err(format!("state.json has no plugins.ops card: {ops:?}"));
    }
    let (ctype, text) = &s.metrics;
    let running = format!("pgp_run_running{{run=\"{SAMPLE}\"}} 1");
    if !ctype.contains("version=0.0.4")
        || !text.lines().any(|l| l == "pgp_runs 1")
        || !text.lines().any(|l| l == running)
    {
        return Err(format!("/metrics: {ctype} {text}"));
    }
    let broken = s
        .broken
        .as_ref()
        .ok_or("state.json fails on a bad history")?;
    let err = key(broken, "plugins")
        .and_then(|p| key(p, "ops"))
        .and_then(|o| key(o, "error"));
    match err {
        Some(Value::Str(e))
            if e.contains("could not be read") && key(broken, "steps").is_some() =>
        {
            Ok(())
        }
        _ => Err(format!("a hook error gave {err:?}")),
    }
}

/// While `bin/owngenome_live.py` lives, it must load the plugin and serve
/// `/metrics`; the check goes with the Python server. Two controls run the
/// site on a copy of the plugin with its hook unregistered, and with its
/// catch-all cut: both must fail the check.
#[test]
fn the_python_site_loads_the_plugin_serves_metrics_and_keeps_a_hook_error_on_the_card() {
    let bin = repo().join("bin");
    let plugin = bin.join("owngenome_live_ops.py");
    if !bin.join("owngenome_live.py").is_file() || !plugin.is_file() {
        eprintln!("bin/owngenome_live.py is gone: the Python site check retired with it");
        return;
    }
    let served = python_site(None).unwrap_or_else(|e| panic!("{e}"));
    site_serves_the_plugin(&served).unwrap_or_else(|e| panic!("{e}"));
    let src = std::fs::read_to_string(&plugin).unwrap();
    let catch_all = "        except Exception as e:  # anything else: the card says so, the rest of the run's state stays\n            return {\"error\": f\"run health could not be read ({type(e).__name__})\"}\n";
    for (what, from, to) in [
        ("unregistered", "    app.state_hook(\"ops\", hook)\n", ""),
        ("no catch-all", catch_all, ""),
    ] {
        assert!(
            src.contains(from),
            "the control found no line to cut ({what})"
        );
        let dir = tmp().join("plugins");
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("owngenome_live_ops.py"), src.replacen(from, to, 1)).unwrap();
        let red = python_site(Some(&dir)).and_then(|s| site_serves_the_plugin(&s));
        assert!(
            red.is_err(),
            "the site check passed a plugin with its hook {what}"
        );
    }
}

// --- parity with the Python plugin ----------------------------------------------------------

fn git(args: &[&str]) -> std::process::Output {
    Command::new("git")
        .arg("-C")
        .arg(repo())
        .args(args)
        .output()
        .expect("git runs")
}

/// The folder holding `owngenome_live_ops.py` and the modules it imports:
/// `bin/` while the Python site lives, else those files at [`ORACLE`].
fn oracle_bin() -> &'static Path {
    static DIR: OnceLock<PathBuf> = OnceLock::new();
    DIR.get_or_init(|| {
        let bin = repo().join("bin");
        if bin.join("owngenome_live_ops.py").is_file() {
            return bin;
        }
        if !git(&["cat-file", "-e", &format!("{ORACLE}^{{commit}}")])
            .status
            .success()
        {
            let f = git(&["fetch", "-q", "--depth=1", "origin", ORACLE]);
            assert!(f.status.success(), "commit {ORACLE} could not be fetched");
        }
        let root = tmp().join("oracle");
        std::fs::create_dir_all(root.join("bin")).unwrap();
        std::fs::create_dir_all(root.join("assets")).unwrap();
        for f in ["owngenome_live_ops", "run_eta", "run_headroom", "run_state"] {
            let o = git(&["show", &format!("{ORACLE}:bin/{f}.py")]);
            assert!(o.status.success(), "git show {ORACLE}:bin/{f}.py");
            std::fs::write(root.join("bin").join(format!("{f}.py")), o.stdout).unwrap();
        }
        std::fs::copy(
            repo().join("assets/steps.tsv"),
            root.join("assets/steps.tsv"),
        )
        .unwrap();
        root.join("bin")
    })
}

const PY: &str = r#"
import json, os, sys, types
a = json.loads(sys.argv[1])
sys.path.insert(0, a["bin"])
import run_headroom, run_state
import owngenome_live_ops as ops

class FakeFs:
    def __init__(self, used_ratio):
        self.f_frsize = self.f_bsize = 4096
        self.f_blocks = 1000000
        self.f_bfree = self.f_bavail = int(self.f_blocks * (1 - used_ratio))
        self.f_files, self.f_ffree = 1000, 900

run_headroom._statvfs = lambda p: FakeFs(a["disk"])
events = [json.loads(l) for l in open(a["events"]) if l.strip()]
steps = run_state.read_steps()
status_path = os.path.join(a["run_dir"], "logs", "run_status.tsv")
status = run_state.read_run_status(status_path) if os.path.isfile(status_path) else None
cur = ops.current_run(events)
state = run_state.build(a["sample"], steps, status, None, [], run_state.event_tasks(cur, a["sample"]), cur, None, a["now"])
if a["stopped"]:
    state["run"]["state"] = "stopped"
cfg = types.SimpleNamespace(root=a["root"], steps=run_state.STEPS_TSV, eta_box="opsbox", eta_cores=8,
                            disk_warn_percent=a["dwp"], mem_warn_percent="90", stall_minutes=a["stall"])
o = ops.run_ops(a["run_dir"], state, events, cfg, now=a["now"], steps=steps, physical=a["physical"])
sys.stdout.write(json.dumps(o) + "\n---\n" + ops.metrics_text([(a["sample"], o)], 0.9, a["physical"], 0))
"#;

/// `text` with the tree's folder written as `{TMP}` and the repository as
/// `{REPO}`, so an answer does not hold this run's paths.
fn norm(t: &Tree, text: &str) -> String {
    text.replace(&t.root.parent().unwrap().display().to_string(), "{TMP}")
        .replace(&repo().display().to_string(), "{REPO}")
}

/// RECORD (this commit only): with `OG_RECORD_VECTORS=<dir>`, appends the
/// plugin's answer to `<dir>/live_ops.jsonl` under `label`.
fn record(label: &str, t: &Tree, out: &str) {
    use std::io::Write;
    static LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
    let Some(dir) = std::env::var_os("OG_RECORD_VECTORS") else {
        return;
    };
    let line = pycompat::json::dumps(
        &Value::Array(vec![Value::Str(label.to_owned()), Value::Str(norm(t, out))]),
        &pycompat::json::Dump::default(),
    )
    .unwrap();
    let _g = LOCK.lock().unwrap();
    std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(PathBuf::from(dir).join("live_ops.jsonl"))
        .unwrap()
        .write_all(format!("{line}\n").as_bytes())
        .unwrap();
}

/// The Python card and metrics text of the tree at `c`.
fn py_ops(label: &str, t: &Tree, c: &Case) -> (Value, String) {
    let mut args = indexmap::IndexMap::new();
    let mut put = |k: &str, v: Value| {
        args.insert(k.to_owned(), v);
    };
    put("bin", Value::Str(s(oracle_bin())));
    put("events", Value::Str(s(&t.events_file)));
    put("run_dir", Value::Str(s(&t.run_dir)));
    put("root", Value::Str(s(&t.root)));
    put("sample", Value::Str(t.sample.clone()));
    put("now", Value::Float(c.now));
    put("disk", Value::Float(c.disk));
    put("physical", Value::Int(c.physical as i128));
    put("stopped", Value::Bool(c.stopped));
    put("stall", Value::Float(c.stall_minutes));
    put(
        "dwp",
        Value::Str(c.disk_warn_percent.unwrap_or("90").into()),
    );
    let arg =
        pycompat::json::dumps(&Value::Object(args), &pycompat::json::Dump::default()).unwrap();
    let exe = std::env::var("OG_PYTHON").unwrap_or_else(|_| "python3".into());
    let o = Command::new(&exe)
        .args(["-c", PY, &arg])
        .env("TZ", "UTC")
        .env_remove("PGP_BOX")
        .env_remove("OWNGENOME_DISK_WARN_PERCENT")
        .env_remove("OWNGENOME_MEM_WARN_PERCENT")
        .output()
        .unwrap_or_else(|e| panic!("cannot run the Python oracle {exe}: {e}; set OG_PYTHON"));
    assert!(
        o.status.success(),
        "the Python plugin failed: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    let out = String::from_utf8(o.stdout).unwrap();
    record(label, t, &out);
    let (card, text) = out.split_once("\n---\n").expect("card and text");
    (loads(card).unwrap(), text.to_owned())
}

/// The same by meaning: numbers by value, the rest equal, keys in order.
fn same(a: &Value, b: &Value, at: &str) -> Result<(), String> {
    let num = |v: &Value| match v {
        Value::Int(i) => Some(*i as f64),
        Value::Float(x) => Some(*x),
        _ => None,
    };
    match (a, b) {
        (Value::Object(x), Value::Object(y)) => {
            let (kx, ky): (Vec<_>, Vec<_>) = (x.keys().collect(), y.keys().collect());
            if kx != ky {
                return Err(format!("{at}: keys {kx:?} != {ky:?}"));
            }
            x.iter()
                .try_for_each(|(k, v)| same(v, &y[k], &format!("{at}.{k}")))
        }
        (Value::Array(x), Value::Array(y)) if x.len() == y.len() => x
            .iter()
            .zip(y)
            .enumerate()
            .try_for_each(|(i, (v, w))| same(v, w, &format!("{at}[{i}]"))),
        _ => match (num(a), num(b)) {
            (Some(x), Some(y)) if (x - y).abs() <= 1e-9 * x.abs().max(1.0) => Ok(()),
            (Some(_), Some(_)) => Err(format!("{at}: {a:?} != {b:?}")),
            _ if a == b => Ok(()),
            _ => Err(format!("{at}: {a:?} != {b:?}")),
        },
    }
}

fn parity(t: &Tree, c: &Case, what: &str) {
    let o = ops_of(t, c);
    let (card, text) = py_ops(what, t, c);
    if let Err(e) = same(&o.to_json(), &card, what) {
        panic!("og and the Python plugin differ: {e}");
    }
    let mine = metrics_text(&[(&t.sample, &o)], Some(0.9), Some(c.physical), 0);
    assert_eq!(mine, text, "{what}: the metrics text differs");
}

#[test]
fn the_card_and_metrics_match_the_python_plugin_on_the_ops_fixture() {
    let quiet = T0 + 302.0 + 7200.0;
    let mut cases: Vec<(&str, bool, Vec<Value>, Case)> = vec![
        ("early", true, vec![], case(T0 + 330.0)),
        ("late and stalled", true, vec![], case(quiet)),
        ("failed task", true, vec![failed_task()], case(T0 + 400.0)),
        ("exit 1", true, vec![finished(1)], case(T0 + 400.0)),
        ("exit 0", true, vec![finished(0)], case(T0 + 400.0)),
        (
            "error then exit 0",
            true,
            vec![
                ev(r#"{"event":"error","utcTime":"2026-01-02T10:05:50Z"}"#),
                finished(0),
            ],
            case(T0 + 400.0),
        ),
        ("bare nextflow run", false, vec![], case(T0 + 330.0)),
    ];
    let mut c = case(quiet);
    c.stopped = true;
    cases.push(("stopped", true, vec![], c));
    let mut c = case(T0 + 302.0 + 90.0);
    c.stall_minutes = 1.0;
    cases.push(("inside its usual time", true, vec![], c));
    let mut c = case(quiet);
    c.stall_minutes = 300.0;
    cases.push(("stall 300", true, vec![], c));
    let mut c = case(T0 + 330.0);
    c.physical = 16 * G;
    cases.push(("16 GB host", true, vec![], c.clone()));
    c.physical = 64 * G;
    c.disk = 0.95;
    cases.push(("disk 95%", true, vec![], c.clone()));
    c.disk = 0.985;
    cases.push(("disk 98.5%", true, vec![], c.clone()));
    c.disk = 0.95;
    c.disk_warn_percent = Some("96");
    cases.push(("disk warn 96", true, vec![], c));
    for (what, run_all, extra, c) in &cases {
        parity(&make_tree(*run_all, extra), c, what);
    }
}

#[test]
fn the_card_and_metrics_match_the_python_plugin_on_the_eta_and_headroom_fixtures() {
    let d = tmp();
    let root = d.join("genome");
    let (work, out) = (d.join("nfwork"), d.join("nfout"));
    let info = root.join("pipeline_info");
    let run_dir = root.join("demo1");
    for p in [
        &info,
        &work,
        &out,
        &run_dir.join("logs"),
        &run_dir.join("nextflow/work"),
    ] {
        std::fs::create_dir_all(p).unwrap();
    }
    // An earlier run's trace makes this box's history.
    let old = info.join("trace_2026-01-01_09-59-00.txt");
    std::fs::copy(repo().join("tests/fixtures/run_eta/trace_new.txt"), &old).unwrap();
    set_mtime(&old, T0 - 86400.0);
    og_runlog::eta::ingest(&root, "opsbox", steps(), 16, &UTC).unwrap();
    // This run: started at T0, its trace and its messages.
    let new = info.join("trace_2026-01-02_10-00-00.txt");
    std::fs::copy(repo().join("tests/fixtures/run_headroom/trace.txt"), &new).unwrap();
    set_mtime(&new, T0 + 60.0);
    let status = run_dir.join("logs/run_status.tsv");
    std::fs::write(
        &status,
        format!(
            "meta\tstarted_epoch\t{}\nmeta\tstarted_utc\t2026-01-02 10:00:00 UTC\nstep\t01b\tplanned\nstep\t02\tplanned\nstep\t03\tplanned\nstep\t13\tplanned\n",
            T0 as i64
        ),
    )
    .unwrap();
    set_mtime(&status, T0);
    let mut events = vec![ev(&format!(
        r#"{{"event":"started","utcTime":"2026-01-02T10:00:00Z","metadata":{{"parameters":{{"outdir":"{}","max_memory":"64.GB"}},"workflow":{{"workDir":"{}"}}}}}}"#,
        s(&out),
        s(&work)
    ))];
    events.extend(
        og_runlog::files::read_events(Some(
            &repo().join("tests/fixtures/run_headroom/events.jsonl"),
        ))
        .unwrap(),
    );
    let events_file = d.join("events.jsonl");
    write_events(&events_file, &events);
    let t = Tree {
        root,
        run_dir,
        events,
        events_file,
        work,
        out,
        sample: "demo1".into(),
    };
    let mut c = case(T0 + 3600.0);
    c.physical = 32 * G;
    parity(&t, &c, "eta and headroom fixtures");
    let o = ops_of(&t, &c);
    // The case is worth comparing: hours left from the box, a retry, a host warning.
    assert!(
        o.hours_left.is_some() && o.memory.retries == 1 && o.memory.max_memory_above_ram,
        "{o:?}"
    );
    // Control: a card that differs in one number must be seen.
    let (card, _) = py_ops("eta and headroom fixtures", &t, &c);
    let mut changed = o.clone();
    changed.memory.retries = 2;
    assert!(same(&changed.to_json(), &card, "control").is_err());
}

// --- the planted mistakes: each test below must FAIL --------------------------------------

#[cfg(feature = "fault-inject")]
mod planted {
    use og_runlog::eta::fault::plant;

    macro_rules! planted {
        ($($name:ident: $fault:literal => $check:path;)*) => {$(
            #[test]
            fn $name() {
                plant($fault);
                $check();
            }
        )*};
    }

    planted! {
        planted_ops_stall_ignores_quiet: "ops-stall-ignores-quiet" => super::behaviour;
        planted_ops_stopped_not_stalled: "ops-stopped-not-stalled" => super::behaviour;
        planted_ops_stall_ignores_steps: "ops-stall-ignores-steps" => super::behaviour;
        planted_ops_late_dropped: "ops-late-dropped" => super::behaviour;
        planted_ops_exit_ignored: "ops-exit-ignored" => super::behaviour;
        planted_ops_disk_warn_ignored: "ops-disk-warn-ignored" => super::behaviour;
        planted_ops_retries_uncounted: "ops-retries-uncounted" => super::behaviour;
        planted_ops_bare_run_no_work: "ops-bare-run-no-work" => super::behaviour;
        // The parity cases must see a port that drifts from the Python plugin.
        planted_ops_parity_retries: "ops-retries-uncounted" => super::the_card_and_metrics_match_the_python_plugin_on_the_ops_fixture;
        planted_ops_parity_late: "ops-late-dropped" => super::the_card_and_metrics_match_the_python_plugin_on_the_ops_fixture;
        planted_ops_parity_headroom: "headroom-exits-143-only" => super::the_card_and_metrics_match_the_python_plugin_on_the_eta_and_headroom_fixtures;
    }
}
