//! og-runlog's headroom, as the live site's ops card reads it: synthetic
//! disks on a fake host, the tasks of `tests/fixtures/run_headroom`, the
//! weblog metadata, and the retried exits against `conf/base.config`.
//!
//! With `--features fault-inject` each `planted_*` test turns on one planted
//! mistake and runs the check that must catch it: every one of them must
//! fail (Guard's rust job requires it). A check that cannot fail is no check.

use jiff::tz::TimeZone;
use og_runlog::eta::{TraceTask, read_trace};
use og_runlog::headroom::{
    self, Disk, MEMORY_EXITS, Request, Snapshot, TaskMemory, TaskRow, event_rows, metadata,
    task_memory, warn_percent,
};
use og_runlog::host::{FsStat, Host, System};
use og_runlog::steps::{Step, read_steps};
use pycompat::json::{Value, loads};
use std::cell::RefCell;
use std::io;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const G: i128 = 1 << 30;

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn fixture(name: &str) -> PathBuf {
    repo().join("tests/fixtures/run_headroom").join(name)
}

fn steps() -> Vec<Step> {
    read_steps(&repo().join("assets/steps.tsv")).expect("assets/steps.tsv reads")
}

/// A new empty folder of this test run.
fn tmp() -> PathBuf {
    static N: AtomicUsize = AtomicUsize::new(0);
    let d = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!(
        "ops-headroom-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// The warn percents the checks expect come from no environment variable.
fn no_warn_env() {
    for k in ["OWNGENOME_DISK_WARN_PERCENT", "OWNGENOME_MEM_WARN_PERCENT"] {
        assert!(
            std::env::var_os(k).is_none(),
            "unset {k}: the checks test the default warn percents"
        );
    }
}

/// A host whose every file system is `used_pct`% used; its device ids are
/// `devs`, else 1; existence is the real one.
struct Fake {
    used_pct: u64,
    files: u64,
    ffree: u64,
    devs: Vec<(PathBuf, u64)>,
    phys: Option<u128>,
    seen: RefCell<Vec<PathBuf>>,
}

fn fake(used_pct: u64) -> Fake {
    Fake {
        used_pct,
        files: 1000,
        ffree: 900,
        devs: Vec::new(),
        phys: Some(64 << 30),
        seen: RefCell::new(Vec::new()),
    }
}

impl Host for Fake {
    fn statvfs(&self, path: &Path) -> io::Result<FsStat> {
        self.seen.borrow_mut().push(path.to_path_buf());
        let free = 1000 - self.used_pct * 10;
        Ok(FsStat {
            frsize: 4096,
            bsize: 4096,
            blocks: 1000,
            bfree: free,
            bavail: free,
            files: self.files,
            ffree: self.ffree,
        })
    }
    fn dev(&self, path: &Path) -> io::Result<u64> {
        Ok(self
            .devs
            .iter()
            .find(|(p, _)| p == path)
            .map_or(1, |(_, d)| *d))
    }
    fn exists(&self, path: &Path) -> bool {
        path.exists()
    }
    fn physical_memory(&self) -> Option<u128> {
        self.phys
    }
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

fn default_warn() -> f64 {
    warn_percent(None, None).unwrap()
}

fn snapshot_of(host: &dyn Host, r: &Request<'_>) -> Snapshot {
    headroom::snapshot(host, &steps(), r).expect("the warn percents are good")
}

fn request<'a>(work: Option<&'a str>, outdir: Option<&'a str>) -> Request<'a> {
    Request {
        events: &[],
        trace_rows: &[],
        work_dir: work,
        outdir,
        max_memory: None,
        disk_warn_percent: None,
        mem_warn_percent: None,
    }
}

fn by_process(ts: Vec<TaskMemory>) -> std::collections::HashMap<String, TaskMemory> {
    ts.into_iter().map(|t| (t.process.clone(), t)).collect()
}

fn close(a: Option<f64>, b: f64) -> bool {
    a.is_some_and(|a| (a - b).abs() < 1e-6)
}

fn trace(name: &str) -> Vec<TraceTask> {
    read_trace(&fixture(name), &TimeZone::UTC)
}

fn events() -> Vec<Value> {
    og_runlog::files::read_events(Some(&fixture("events.jsonl"))).unwrap()
}

fn trace_memory(rows: &[TraceTask]) -> Vec<TaskMemory> {
    let rows: Vec<TaskRow<'_>> = rows.iter().map(TaskRow::Trace).collect();
    task_memory(&rows, default_warn(), &steps())
}

fn weblog_memory(rows: &[&Value]) -> Vec<TaskMemory> {
    let rows: Vec<TaskRow<'_>> = rows.iter().map(|v| TaskRow::Weblog(v)).collect();
    task_memory(&rows, default_warn(), &steps())
}

// --- disks --------------------------------------------------------------------------------

#[test]
fn disk_91_warns_at_the_default_89_does_not() {
    let d = tmp();
    let hi = headroom::disk(&fake(91), &[("work", Some(s(&d)))], default_warn());
    let lo = headroom::disk(&fake(89), &[("work", Some(s(&d)))], default_warn());
    assert!(
        hi.len() == 1
            && hi[0].warn
            && close(hi[0].used_ratio, 0.91)
            && lo.len() == 1
            && !lo[0].warn,
        "{hi:?} {lo:?}"
    );
}

/// The environment half runs in a child process (the test binary itself):
/// a test may not set a variable of its own process.
#[test]
fn disk_warn_percent_80_and_its_environment_variable_warn_at_85() {
    no_warn_env();
    let d = tmp();
    let host = fake(85);
    let plain = headroom::lines(&snapshot_of(&host, &request(Some(s(&d)), None)));
    let mut with = request(Some(s(&d)), None);
    with.disk_warn_percent = Some("80");
    let arg = headroom::lines(&snapshot_of(&host, &with)).join("\n");
    let mut child = Command::new(std::env::current_exe().unwrap());
    child
        .args([
            "child_disk_warn_from_the_environment",
            "--exact",
            "--ignored",
            "--quiet",
        ])
        .env("OWNGENOME_DISK_WARN_PERCENT", "80")
        .env("OPS_CHILD_DIR", &d);
    #[cfg(feature = "fault-inject")]
    if let Some(f) = og_runlog::eta::fault::planted() {
        child.env("OG_PLANTED", f);
    }
    let out = child.output().unwrap();
    assert!(
        !plain[0].contains("WARN") && arg.contains("WARN: at or above 80%") && out.status.success(),
        "{plain:?} {arg} {}",
        String::from_utf8_lossy(&out.stdout)
    );
}

#[test]
#[ignore = "run by disk_warn_percent_80_and_its_environment_variable_warn_at_85"]
fn child_disk_warn_from_the_environment() {
    #[cfg(feature = "fault-inject")]
    if let Ok(f) = std::env::var("OG_PLANTED") {
        og_runlog::eta::fault::plant(Box::leak(f.into_boxed_str()));
    }
    let d = std::env::var("OPS_CHILD_DIR").expect("run by its parent test");
    let text = headroom::lines(&snapshot_of(&fake(85), &request(Some(&d), None))).join("\n");
    assert!(text.contains("WARN: at or above 80%"), "{text}");
}

#[test]
fn a_missing_work_folder_is_measured_on_its_parent() {
    let d = tmp();
    let host = fake(50);
    let miss = headroom::disk(
        &host,
        &[("work", Some(s(&d.join("not/yet"))))],
        default_warn(),
    );
    let seen = host.seen.borrow().clone();
    assert!(
        miss.len() == 1 && miss[0].missing && seen == [d.clone()],
        "{miss:?} {seen:?}"
    );
}

#[test]
fn one_device_reports_once_two_report_twice() {
    let d = tmp();
    let (a, b) = (d.join("a"), d.join("b"));
    std::fs::create_dir_all(&a).unwrap();
    std::fs::create_dir_all(&b).unwrap();
    let paths = [("work", Some(s(&a))), ("results", Some(s(&b)))];
    let mut host = fake(50);
    host.devs = vec![(a.clone(), 7), (b.clone(), 7)];
    let one = headroom::disk(&host, &paths, default_warn());
    host.devs = vec![(a.clone(), 7), (b.clone(), 8)];
    let two = headroom::disk(&host, &paths, default_warn());
    let fs = |v: &[Disk]| v.iter().map(|d| d.fs.clone()).collect::<Vec<_>>();
    assert!(
        fs(&one) == ["work and results"] && fs(&two) == ["work", "results"],
        "{one:?} {two:?}"
    );
}

#[test]
fn f_files_0_gives_no_inode_ratio() {
    let d = tmp();
    let mut none = fake(50);
    none.files = 0;
    none.ffree = 0;
    let noi = headroom::disk(&none, &[("work", Some(s(&d)))], default_warn());
    let withi = headroom::disk(&fake(50), &[("work", Some(s(&d)))], default_warn());
    assert!(
        noi[0].inodes_used_ratio.is_none() && close(withi[0].inodes_used_ratio, 0.1),
        "{noi:?} {withi:?}"
    );
}

#[test]
fn a_real_statvfs_gives_a_ratio_in_0_to_1_and_a_disk_line() {
    no_warn_env();
    let d = tmp();
    let real = headroom::disk(&System, &[("work", Some(s(&d)))], default_warn());
    let lines = headroom::lines(&snapshot_of(&System, &request(Some(s(&d)), Some(s(&d)))));
    assert!(
        real.len() == 1
            && real[0].used_ratio.is_some_and(|r| (0.0..=1.0).contains(&r))
            && real[0].size_bytes > 0
            && lines[0].starts_with("disk work and results: "),
        "{real:?} {lines:?}"
    );
}

// --- memory -------------------------------------------------------------------------------

#[test]
fn seven_point_six_of_8_gb_warns_6_of_8_does_not() {
    let t = by_process(trace_memory(&trace("trace.txt")));
    assert!(
        t["FASTP"].warn
            && close(t["FASTP"].ratio, 0.95)
            && !t["VEP"].warn
            && close(t["VEP"].ratio, 0.75),
        "{t:?}"
    );
}

#[test]
fn a_trace_size_and_weblog_bytes_give_the_same_ratio() {
    let t = by_process(trace_memory(&trace("trace.txt")));
    let evs = events();
    let e = by_process(weblog_memory(&event_rows(&evs)));
    assert!(
        t["FASTP"].ratio.is_some()
            && close(e["FASTP"].ratio, t["FASTP"].ratio.unwrap())
            && e["FASTP"].step.as_deref() == Some("01b"),
        "{:?} {:?}",
        t["FASTP"],
        e["FASTP"]
    );
}

#[test]
fn no_memory_column_gives_no_ratio() {
    let t = trace_memory(&trace("trace_no_memory.txt"));
    assert!(t.len() == 1 && t[0].ratio.is_none() && !t[0].warn, "{t:?}");
}

#[test]
fn exit_137_then_attempt_2_is_one_retry_and_a_running_task_shows_its_request() {
    let evs = events();
    let e = by_process(weblog_memory(&event_rows(&evs)));
    let (v, f, dv) = (&e["VEP"], &e["FASTP"], &e["DEEPVARIANT"]);
    assert!(
        v.oom_retry == 1
            && v.attempt == 2
            && close(v.ratio, 9.0 / 16.0)
            && f.oom_retry == 0
            && dv.state == "running"
            && dv.ratio.is_none()
            && dv.request_bytes == Some(16 * G)
            && dv.peak_rss_bytes.is_none(),
        "{e:?}"
    );
}

fn task(hash: &str, status: &str, exit: &str, peak: Option<i128>) -> TraceTask {
    TraceTask {
        trace_file: String::new(),
        hash: hash.into(),
        attempt: "1".into(),
        process: Some("FASTP".into()),
        tag: Some("demo".into()),
        status: status.into(),
        exit: exit.into(),
        cpus: None,
        memory_bytes: Some(8 * G),
        peak_rss_bytes: peak,
        realtime_s: None,
        start: None,
        complete: None,
    }
}

#[test]
fn a_running_row_after_a_finished_one_of_the_same_attempt_keeps_the_finished_ratio() {
    let rows = [
        task("ab/1", "COMPLETED", "0", Some(6 * G)),
        task("ab/2", "RUNNING", "", None),
    ];
    let m = trace_memory(&rows);
    assert!(
        m.len() == 1 && m[0].state == "completed" && close(m[0].ratio, 0.75),
        "{m:?}"
    );
}

#[test]
fn a_running_message_after_a_finished_one_of_the_same_task_id_keeps_the_finished_ratio() {
    let fin = loads(&format!(
        r#"{{"task_id":7,"process":"VEP","name":"VEP (demo)","attempt":1,"status":"COMPLETED","exit":0,"memory":{},"peak_rss":{}}}"#,
        8 * G,
        6 * G
    ))
    .unwrap();
    let run = loads(&format!(
        r#"{{"task_id":7,"process":"VEP","name":"VEP (demo)","attempt":1,"status":"RUNNING","memory":{}}}"#,
        8 * G
    ))
    .unwrap();
    let m = weblog_memory(&[&fin, &run]);
    assert!(
        m.len() == 1 && m[0].state == "completed" && close(m[0].ratio, 0.75),
        "{m:?}"
    );
}

#[test]
fn max_memory_64_gb_on_a_32_gb_host_warns_on_128_it_does_not() {
    no_warn_env();
    let mut r = request(None, None);
    r.max_memory = Some("64.GB");
    let mut host = fake(50);
    host.phys = Some(32 << 30);
    let small = snapshot_of(&host, &r).host;
    host.phys = Some(128 << 30);
    let big = snapshot_of(&host, &r).host;
    assert!(
        small.warn && small.max_memory_bytes == Some(64 * G) && !big.warn,
        "{small:?} {big:?}"
    );
}

// --- the snapshot and the weblog metadata ---------------------------------------------------

fn strings(snap: &Snapshot) -> Vec<String> {
    let mut out: Vec<String> = snap.disk.iter().map(|d| d.fs.clone()).collect();
    for t in &snap.tasks {
        out.extend([t.process.clone(), t.state.clone()]);
        out.extend(t.step.clone());
    }
    out
}

#[test]
fn the_snapshot_holds_no_path_and_takes_its_folders_from_the_metadata() {
    no_warn_env();
    let d = tmp();
    let (work, res) = (d.join("work"), d.join("results"));
    std::fs::create_dir_all(&work).unwrap();
    std::fs::create_dir_all(&res).unwrap();
    let meta = loads(&format!(
        r#"{{"event":"started","metadata":{{"workflow":{{"workDir":"{}"}},"parameters":{{"outdir":"{}","max_memory":"64.GB"}}}}}}"#,
        s(&work),
        s(&res)
    ))
    .unwrap();
    let mut evs = vec![meta];
    evs.extend(events());
    let rows = trace("trace.txt");
    let mut host = fake(50);
    host.devs = vec![(work.clone(), 1), (res.clone(), 2)];
    host.phys = Some(32 << 30);
    let mut r = request(None, None);
    r.events = &evs;
    r.trace_rows = &rows;
    let snap = snapshot_of(&host, &r);
    let bad: Vec<String> = strings(&snap)
        .into_iter()
        .filter(|x| x.contains('/'))
        .collect();
    let fs: Vec<&str> = snap.disk.iter().map(|d| d.fs.as_str()).collect();
    let mut procs: Vec<&str> = snap.tasks.iter().map(|t| t.process.as_str()).collect();
    procs.sort_unstable();
    procs.dedup();
    assert!(
        bad.is_empty()
            && fs == ["work", "results"]
            && snap.host.warn
            && snap.warn_used_ratio == 0.9
            && snap.mem_warn_ratio == 0.9
            && procs == ["DEEPVARIANT", "FASTP", "VEP"],
        "{bad:?} {snap:?}"
    );
}

fn md(text: &str) -> (Option<Value>, Option<Value>, Option<Value>) {
    metadata(&[loads(text).unwrap()])
}

fn st(v: &str) -> Option<Value> {
    Some(Value::Str(v.into()))
}

#[test]
fn a_relative_outdir_is_under_launch_dir_an_absolute_one_is_kept() {
    let rel = md(
        r#"{"metadata":{"workflow":{"workDir":"/x/work","launchDir":"/x"},"parameters":{"outdir":"./results"}}}"#,
    );
    let ab = md(
        r#"{"metadata":{"workflow":{"workDir":"/x/work","launchDir":"/x"},"parameters":{"outdir":"/y/results"}}}"#,
    );
    assert_eq!(rel, (st("/x/work"), st("/x/./results"), None));
    assert_eq!(ab, (st("/x/work"), st("/y/results"), None));
}

#[test]
fn metadata_under_params_gives_the_outdir_and_max_memory() {
    let alt = md(
        r#"{"metadata":{"workflow":{"workDir":"/x/work"},"params":{"outdir":"/y/r","max_memory":"8.GB"}}}"#,
    );
    assert_eq!(alt, (st("/x/work"), st("/y/r"), st("8.GB")));
}

/// The exits of `errorStrategy = { task.exitStatus in [...] ? 'retry' ...`.
fn retried_exits(config: &str) -> Option<Vec<String>> {
    let at = config.find("errorStrategy")?;
    let rest = &config[at..];
    let open = rest.find("task.exitStatus in [")? + "task.exitStatus in [".len();
    let close = rest[open..].find(']')? + open;
    Some(
        rest[open..close]
            .split(',')
            .map(|x| x.trim().to_owned())
            .collect(),
    )
}

#[test]
fn memory_exits_are_the_exits_base_config_retries() {
    let base = std::fs::read_to_string(repo().join("conf/base.config")).unwrap();
    let mut want: Vec<String> = MEMORY_EXITS.iter().map(|x| (*x).to_owned()).collect();
    want.sort();
    let sorted = |v: Option<Vec<String>>| {
        v.map(|mut v| {
            v.sort();
            v
        })
    };
    assert_eq!(sorted(retried_exits(&base)), Some(want.clone()));
    // Control: the same reading of a config without 137 must differ.
    let broken = base.replace("[143,137,", "[143,");
    assert_ne!(broken, base, "the control found no exit list to break");
    assert_ne!(sorted(retried_exits(&broken)), Some(want));
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
        planted_headroom_default_95: "headroom-default-95" => super::disk_91_warns_at_the_default_89_does_not;
        planted_headroom_percent_ignored: "headroom-percent-ignored" => super::disk_warn_percent_80_and_its_environment_variable_warn_at_85;
        planted_headroom_missing_not_parent: "headroom-missing-not-parent" => super::a_missing_work_folder_is_measured_on_its_parent;
        planted_headroom_no_merge: "headroom-no-merge" => super::one_device_reports_once_two_report_twice;
        planted_headroom_inodes_always: "headroom-inodes-always" => super::f_files_0_gives_no_inode_ratio;
        planted_headroom_used_above_size: "headroom-used-above-size" => super::a_real_statvfs_gives_a_ratio_in_0_to_1_and_a_disk_line;
        planted_headroom_warn_strict: "headroom-warn-strict" => super::seven_point_six_of_8_gb_warns_6_of_8_does_not;
        planted_headroom_weblog_decimal_peak: "headroom-weblog-decimal-peak" => super::a_trace_size_and_weblog_bytes_give_the_same_ratio;
        planted_headroom_request_default_8g: "headroom-request-default-8g" => super::no_memory_column_gives_no_ratio;
        planted_headroom_exits_143_only: "headroom-exits-143-only" => super::exit_137_then_attempt_2_is_one_retry_and_a_running_task_shows_its_request;
        planted_headroom_running_wins: "headroom-running-wins" => super::a_running_row_after_a_finished_one_of_the_same_attempt_keeps_the_finished_ratio;
        planted_headroom_replaces_always: "headroom-replaces-always" => super::a_running_message_after_a_finished_one_of_the_same_task_id_keeps_the_finished_ratio;
        planted_headroom_host_never_warns: "headroom-host-never-warns" => super::max_memory_64_gb_on_a_32_gb_host_warns_on_128_it_does_not;
        planted_headroom_label_is_path: "headroom-label-is-path" => super::the_snapshot_holds_no_path_and_takes_its_folders_from_the_metadata;
        planted_headroom_no_launch_dir: "headroom-no-launch-dir" => super::a_relative_outdir_is_under_launch_dir_an_absolute_one_is_kept;
        planted_headroom_parameters_only: "headroom-parameters-only" => super::metadata_under_params_gives_the_outdir_and_max_memory;
    }
}
