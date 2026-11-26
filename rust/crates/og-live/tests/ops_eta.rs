//! og-runlog's step history and hours left, as the live site's ops card
//! reads them: the synthetic traces of `tests/fixtures/run_eta`, synthetic
//! histories and synthetic run states; `nextflow.config`'s trace columns;
//! and the old Python command lines that must refuse to run.
//!
//! With `--features fault-inject` each `planted_*` test turns on one planted
//! mistake and runs the check that must catch it: every one of them must
//! fail (Guard's rust job requires it). A check that cannot fail is no check.

use indexmap::{IndexMap, IndexSet};
use jiff::tz::TimeZone;
use og_runlog::eta::{
    self, EtaArgs, HISTORY_COLUMNS, HistoryRow, Range, Source, TraceTask, box_name, history_path,
    history_rows, ingest, read_trace, remaining, step_ranges,
};
use og_runlog::state::{RunState, State, StepOut, StepState};
use og_runlog::steps::{Step, process_steps, read_steps};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::sync::atomic::{AtomicUsize, Ordering};

const H: f64 = 3600.0;
/// 2026-01-02T10:00:00Z.
const NOW: f64 = 1_767_348_000.0;
const UTC: TimeZone = TimeZone::UTC;

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn steps() -> &'static [Step] {
    static S: OnceLock<Vec<Step>> = OnceLock::new();
    S.get_or_init(|| read_steps(&repo().join("assets/steps.tsv")).expect("assets/steps.tsv reads"))
}

/// A new empty folder of this test run.
fn tmp() -> PathBuf {
    static N: AtomicUsize = AtomicUsize::new(0);
    let d = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!(
        "ops-eta-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// A copy written now (a new modification time, as Python's shutil.copy).
fn copy(from: &Path, to: &Path) {
    std::fs::write(to, std::fs::read(from).unwrap()).unwrap();
}

/// A fresh outdir whose `pipeline_info` holds the named fixture traces.
fn outdir(traces: &[(&str, &str)]) -> PathBuf {
    let d = tmp();
    std::fs::create_dir_all(d.join("pipeline_info")).unwrap();
    for (name, fx) in traces {
        copy(
            &repo().join("tests/fixtures/run_eta").join(fx),
            &d.join("pipeline_info").join(name),
        );
    }
    d
}

fn ingest_as(d: &Path, box_: &str) -> Vec<HistoryRow> {
    ingest(d, &box_name(Some(box_)), steps(), 16, &UTC).unwrap()
}

/// `{step: sorted spans}` of history rows.
fn spans(rows: &[HistoryRow]) -> BTreeMap<String, Vec<f64>> {
    let mut out: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for r in rows {
        if let Ok(v) = r[10].parse::<f64>() {
            let e = out.entry(r[3].clone()).or_default();
            if !e.contains(&v) {
                e.push(v);
            }
        }
    }
    for v in out.values_mut() {
        v.sort_by(f64::total_cmp);
    }
    out
}

const OLD: &str = "trace_2026-01-02_09-59-00.txt";

// --- the history ---------------------------------------------------------------------------

#[test]
fn an_old_and_a_new_trace_give_the_same_spans() {
    let old = ingest_as(&outdir(&[(OLD, "trace_old.txt")]), "t");
    let new = ingest_as(&outdir(&[(OLD, "trace_new.txt")]), "t");
    let want: BTreeMap<String, Vec<f64>> = [
        ("01b", vec![1200.0]),
        ("02", vec![120.0, 4.0 * H]),
        ("03", vec![11.0 * H + 10.0]),
        ("04", vec![1800.0]),
    ]
    .into_iter()
    .map(|(k, v)| (k.to_owned(), v))
    .collect();
    assert_eq!(spans(&old), want);
    assert_eq!(spans(&new), want);
}

#[test]
fn human_and_raw_units_read_the_same() {
    use og_nf::units::{duration, memory, time};
    let raw_18 = ((18.3 * 1024f64.powi(3)) as i128).to_string();
    assert_eq!(duration(Some("1h 2m 3s")), Some(3723.0));
    assert_eq!(duration(Some("3723000")), Some(3723.0));
    assert_eq!(duration(Some("850ms")), Some(0.85));
    assert_eq!(duration(Some("850")), Some(0.85));
    assert_eq!(memory(Some("18.3 GB")), memory(Some(&raw_18)));
    assert_eq!(memory(Some("40 MB")), memory(Some(&(40 << 20).to_string())));
    assert_eq!(time(Some("2026-01-02 10:00:00.000"), &UTC), Some(NOW));
    assert_eq!(time(Some("1767348000000"), &UTC), Some(NOW));
    // And as the history reads a trace: a size and its bytes are one value.
    let d = tmp();
    let head = "task_id\thash\tname\tstatus\texit\tmemory\n";
    std::fs::write(
        d.join("a.txt"),
        format!("{head}1\tab/1\tFASTP (s)\tCOMPLETED\t0\t18.3 GB\n"),
    )
    .unwrap();
    std::fs::write(
        d.join("b.txt"),
        format!("{head}1\tab/1\tFASTP (s)\tCOMPLETED\t0\t{raw_18}\n"),
    )
    .unwrap();
    let (a, b) = (
        read_trace(&d.join("a.txt"), &UTC),
        read_trace(&d.join("b.txt"), &UTC),
    );
    assert!(a[0].memory_bytes.is_some());
    assert_eq!(a[0].memory_bytes, b[0].memory_bytes);
}

#[test]
fn cached_failed_and_unknown_process_tasks_never_enter_the_history() {
    let rows = ingest_as(&outdir(&[(OLD, "trace_new.txt")]), "t");
    let mut hashes: Vec<&str> = rows.iter().map(|r| r[1].as_str()).collect();
    hashes.sort_unstable();
    assert_eq!(
        hashes,
        ["aa000001", "aa000002", "aa000003", "aa000006", "aa000007"]
    );
}

#[test]
fn ingesting_twice_or_a_copied_trace_adds_nothing() {
    let d = outdir(&[(OLD, "trace_new.txt")]);
    ingest_as(&d, "t");
    let path = history_path(&d, "t");
    let first = std::fs::read_to_string(&path).unwrap();
    ingest_as(&d, "t");
    copy(
        &repo().join("tests/fixtures/run_eta/trace_new.txt"),
        &d.join("pipeline_info/trace_2026-01-03_08-00-00.txt"),
    );
    ingest_as(&d, "t");
    assert_eq!(std::fs::read_to_string(&path).unwrap(), first);
}

#[test]
fn two_boxes_write_two_files() {
    let d = outdir(&[(OLD, "trace_new.txt")]);
    ingest_as(&d, "t");
    ingest_as(&d, "box-b");
    let mut files: Vec<String> = std::fs::read_dir(d.join("pipeline_info/step_history"))
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|f| f.ends_with(".tsv"))
        .collect();
    files.sort();
    assert_eq!(files, ["box-b.tsv", "t.tsv"]);
}

#[test]
fn the_history_holds_no_sample_id_and_no_slash() {
    let d = outdir(&[(OLD, "trace_new.txt")]);
    ingest_as(&d, "t");
    let body = std::fs::read_to_string(history_path(&d, "t")).unwrap();
    let head: Vec<&str> = body.lines().next().unwrap().split('\t').collect();
    assert!(
        !body.contains("demo1") && !body.contains('/') && head == HISTORY_COLUMNS,
        "{body}"
    );
}

#[test]
fn the_box_name_keeps_only_letters_digits_underscore_and_dash() {
    assert_eq!(box_name(Some("my box.local")), "my_box_local");
}

// --- the ranges ---------------------------------------------------------------------------

/// History rows of one step: one run per span in hours.
fn hist(step: &str, runs: &[f64]) -> Vec<HistoryRow> {
    runs.iter()
        .enumerate()
        .map(|(i, hours)| {
            [
                format!("trace_{i}.txt"),
                format!("h{i}"),
                "1".into(),
                step.into(),
                "X".into(),
                "no".into(),
                "8".into(),
                String::new(),
                String::new(),
                String::new(),
                format!("{:?}", hours * H),
                format!("2026-01-{:02}T00:00:00Z", i + 1),
                "16".into(),
            ]
        })
        .collect()
}

fn ranges(history: &[HistoryRow], cores: usize) -> IndexMap<String, Range> {
    step_ranges(history, steps(), cores, &IndexSet::new())
}

fn step(id: &str) -> &'static Step {
    steps()
        .iter()
        .find(|s| s.step == id)
        .expect("a catalogue step")
}

#[test]
fn three_runs_give_min_max_with_n_3() {
    let r = &ranges(&hist("02", &[4.0, 3.0, 5.0]), 16)["02"];
    assert_eq!(
        (r.low_s, r.high_s, r.n, r.source),
        (Some(3.0 * H), Some(5.0 * H), 3, Source::Box)
    );
}

#[test]
fn only_the_last_5_runs_count() {
    let r = &ranges(&hist("02", &[9.0, 3.0, 4.0, 5.0, 6.0, 7.0]), 16)["02"];
    assert_eq!((r.low_s, r.high_s, r.n), (Some(3.0 * H), Some(7.0 * H), 5));
}

fn docs_hours(id: &str) -> (f64, f64) {
    let s = step(id);
    (
        s.hours.unwrap() * H,
        s.hours_low.parse::<f64>().unwrap() * H,
    )
}

#[test]
fn an_8_core_box_takes_the_measured_low_end_and_doubles_the_high_end() {
    let (dv, low) = docs_hours("03");
    let r = &ranges(&[], 8)["03"];
    assert!(
        r.source == Source::Docs
            && r.low_s == Some(low)
            && r.high_s == Some(2.0 * dv)
            && r.label.contains("16 threads"),
        "{r:?}"
    );
}

#[test]
fn a_32_core_box_is_not_scaled_below_the_docs_hours() {
    let (dv, low) = docs_hours("03");
    let r = &ranges(&[], 32)["03"];
    assert!(r.low_s == Some(low) && r.high_s == Some(dv), "{r:?}");
}

#[test]
fn on_16_and_32_cores_every_measured_step_has_low_below_high() {
    let measured: Vec<&str> = steps()
        .iter()
        .filter(|s| s.hours_from == "measured")
        .map(|s| s.step.as_str())
        .collect();
    assert!(measured.len() >= 25, "{measured:?}");
    for cores in [16, 32] {
        let r = ranges(&[], cores);
        let one: Vec<&&str> = measured
            .iter()
            .filter(|k| {
                let x = &r[**k];
                !(x.source == Source::Docs
                    && matches!((x.low_s, x.high_s), (Some(l), Some(h)) if l < h))
            })
            .collect();
        assert!(one.is_empty(), "{cores} cores: {one:?}");
    }
}

/// A run state with these steps: (step, state, started_utc); upstream from
/// the catalogue.
fn state(rows: &[(&str, StepState, Option<&str>)], run: RunState) -> State<'static> {
    State {
        sample: "s1".into(),
        generated: NOW,
        run,
        started: None,
        finished_utc: None,
        declared_sex: None,
        status_file: true,
        last_heard: None,
        steps: rows
            .iter()
            .map(|(k, st, t)| StepOut {
                step: step(k),
                state: *st,
                reason: None,
                note: None,
                earlier: None,
                started_utc: t.map(str::to_owned),
                ended_utc: None,
                processes: Vec::new(),
            })
            .collect(),
    }
}

use StepState::{Done, Running, Waiting};

#[test]
fn a_step_the_table_only_estimates_gets_no_number() {
    let est = ranges(&[], 16);
    let unk = remaining(
        &state(
            &[("03", Done, None), ("14", Waiting, None)],
            RunState::Running,
        ),
        &est,
        NOW,
        16,
    );
    assert!(
        est["14"].source == Source::None
            && est["04b"].source == Source::None
            && est["14"].low_s.is_none()
            && unk.unknown_steps == ["14"]
            && unk.low_s == 0.0,
        "{:?} {:?} {unk:?}",
        est["14"],
        est["04b"]
    );
}

fn task(
    hash: &str,
    process: &str,
    status: &str,
    exit: &str,
    start: f64,
    end: f64,
    attempt: &str,
) -> TraceTask {
    TraceTask {
        trace_file: "trace_x.txt".into(),
        hash: hash.into(),
        attempt: attempt.into(),
        process: Some(process.into()),
        tag: Some("s1".into()),
        status: status.into(),
        exit: exit.into(),
        cpus: Some(1),
        memory_bytes: None,
        peak_rss_bytes: None,
        realtime_s: Some(end - start),
        start: Some(start),
        complete: Some(end),
    }
}

#[test]
fn a_step_with_a_task_that_failed_or_was_stopped_gives_no_span() {
    // Step 07: the preprocessor ended, then PharmCAT was stopped. Step 03: a
    // failed first attempt that a second one replaced still counts.
    let rows = history_rows(
        &[
            task(
                "p1",
                "PHARMCAT_PREPROCESS",
                "COMPLETED",
                "0",
                0.0,
                40.0,
                "1",
            ),
            task("p2", "PHARMCAT", "ABORTED", "-", 41.0, 50.0, "1"),
            task("d1", "DEEPVARIANT", "FAILED", "137", 0.0, 100.0, "1"),
            task("d2", "DEEPVARIANT", "COMPLETED", "0", 110.0, 300.0, "2"),
        ],
        &process_steps(steps()),
        16,
    );
    let got: BTreeMap<(String, String), Option<f64>> = rows
        .iter()
        .map(|r| ((r[3].clone(), r[4].clone()), r[10].parse::<f64>().ok()))
        .collect();
    let want: BTreeMap<(String, String), Option<f64>> = [
        (("07".to_owned(), "PHARMCAT_PREPROCESS".to_owned()), None),
        (("03".to_owned(), "DEEPVARIANT".to_owned()), Some(300.0)),
    ]
    .into_iter()
    .collect();
    assert_eq!(got, want);
}

// --- the hours left -----------------------------------------------------------------------

fn range(low: Option<f64>, high: Option<f64>, cpus: Option<i128>) -> Range {
    Range {
        low_s: low,
        high_s: high,
        n: 1,
        source: Source::Box,
        label: String::new(),
        cpus,
    }
}

fn fixed_ranges() -> IndexMap<String, Range> {
    [
        ("02", range(Some(4.0 * H), Some(5.0 * H), Some(8))),
        ("03", range(Some(10.0 * H), Some(12.0 * H), Some(8))),
        ("13", range(Some(4.0 * H), Some(5.0 * H), Some(8))),
        ("23", range(Some(0.1 * H), Some(0.2 * H), Some(1))),
        ("04", range(Some(0.5 * H), Some(H), Some(8))),
        ("33", range(None, None, None)),
    ]
    .into_iter()
    .map(|(k, v)| (k.to_owned(), v))
    .collect()
}

fn close(a: Option<f64>, b: f64) -> bool {
    a.is_some_and(|a| (a - b).abs() < 1e-6)
}

fn chain_state() -> State<'static> {
    state(
        &[
            ("02", Done, None),
            ("03", Waiting, None),
            ("13", Waiting, None),
            ("23", Waiting, None),
            ("04", Waiting, None),
        ],
        RunState::Running,
    )
}

#[test]
fn the_hours_left_are_the_longest_chain() {
    // 03 -> 13 -> 23 is the longest chain: low 10 + 4 + 0.1, high 12 + 5 + 0.2.
    let big = remaining(&chain_state(), &fixed_ranges(), NOW, 64);
    assert!(
        close(Some(big.low_s), 14.1 * H) && close(big.high_s, 17.2 * H),
        "{big:?}"
    );
}

#[test]
fn the_cpu_work_bound_when_it_is_larger() {
    // CPU work: 12*8 + 5*8 + 0.2*1 + 1*8 = 144.2 CPU-hours; over 4 cores 36.05 h.
    let small = remaining(&chain_state(), &fixed_ranges(), NOW, 4);
    assert!(
        close(Some(small.low_s), 14.1 * H) && close(small.high_s, 36.05 * H),
        "{small:?}"
    );
}

#[test]
fn a_running_step_past_its_high_end_makes_the_high_end_unknown_and_says_so() {
    let late = remaining(
        &state(
            &[
                ("02", Done, None),
                ("03", Running, Some("2026-01-01T21:00:00Z")),
                ("13", Waiting, None),
            ],
            RunState::Running,
        ),
        &fixed_ranges(),
        NOW,
        64,
    );
    assert!(
        late.high_s.is_none()
            && late.late_steps == ["03"]
            && late.note.as_deref() == Some("step 03 is taking longer than usual")
            && close(Some(late.low_s), 4.0 * H),
        "{late:?}"
    );
}

#[test]
fn a_running_step_counts_only_what_is_left_of_it() {
    let going = remaining(
        &state(
            &[
                ("02", Done, None),
                ("03", Running, Some("2026-01-02T06:00:00Z")),
                ("13", Waiting, None),
            ],
            RunState::Running,
        ),
        &fixed_ranges(),
        NOW,
        64,
    );
    assert!(
        close(Some(going.low_s), 10.0 * H) && close(going.high_s, 13.0 * H),
        "{going:?}"
    );
}

#[test]
fn a_step_with_no_time_lands_in_unknown_steps_never_0() {
    let r = ranges(&[], 16);
    let unk = remaining(
        &state(
            &[("02", Done, None), ("33", Waiting, None)],
            RunState::Running,
        ),
        &r,
        NOW,
        16,
    );
    assert!(
        unk.unknown_steps == ["33"] && r["33"].low_s.is_none() && r["33"].source == Source::None,
        "{unk:?} {:?}",
        r["33"]
    );
}

#[test]
fn a_finished_run_has_nothing_left() {
    let done = remaining(
        &state(&[("03", Waiting, None)], RunState::Finished),
        &fixed_ranges(),
        NOW,
        16,
    );
    assert!(done.low_s == 0.0 && done.high_s == Some(0.0), "{done:?}");
}

// --- the estimate of an outdir ------------------------------------------------------------

fn now() -> f64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_secs_f64()
}

fn estimate(d: &Path, cores: usize) -> eta::Estimate<'static> {
    eta::estimate(
        d,
        steps(),
        &EtaArgs {
            sample: None,
            cores: Some(cores),
            box_: Some("ci"),
            now: now(),
            tz: &UTC,
        },
    )
    .unwrap()
}

#[test]
fn the_estimate_gives_the_ranges_and_writes_the_history() {
    let d = outdir(&[(OLD, "trace_new.txt")]);
    let est = estimate(&d, 8);
    assert!(
        est.box_ == "ci"
            && est.ranges["02"].source == Source::Box
            && est.ranges["02"].high_s == Some(4.0 * H)
            && est.ranges["13"].source == Source::Docs
            && est.remaining.is_none()
            && history_path(&d, "ci").is_file(),
        "{} {:?}",
        est.box_,
        est.ranges["02"]
    );
}

/// An outdir with a run that run-all.sh started a minute ago.
fn going(d: &Path) -> PathBuf {
    let logs = d.join("demo1/logs");
    std::fs::create_dir_all(&logs).unwrap();
    let status = logs.join("run_status.tsv");
    std::fs::write(
        &status,
        format!(
            "meta\tstarted_epoch\t{}\nstep\t02\tplanned\n",
            now() as i64 - 60
        ),
    )
    .unwrap();
    status
}

#[test]
fn without_a_sample_the_trace_of_a_run_not_finished_is_left_out() {
    let d = outdir(&[(OLD, "trace_new.txt")]);
    let status = going(&d);
    let before = estimate(&d, 16).ranges["02"].source;
    let mut text = std::fs::read_to_string(&status).unwrap();
    text.push_str("meta\tfinished_utc\t2026-01-02 12:00:00 UTC\n");
    std::fs::write(&status, text).unwrap();
    let after = estimate(&d, 16).ranges["02"].source;
    assert_eq!((before, after), (Source::Docs, Source::Box));
}

#[test]
fn an_unreadable_run_status_leaves_out_the_newest_trace() {
    use std::os::unix::fs::PermissionsExt;
    let d = outdir(&[(OLD, "trace_new.txt")]);
    let status = going(&d);
    let mut text = std::fs::read_to_string(&status).unwrap();
    text.push_str("meta\tfinished_utc\t2026-01-02 12:00:00 UTC\n");
    std::fs::write(&status, text).unwrap();
    std::fs::set_permissions(&status, std::fs::Permissions::from_mode(0o000)).unwrap();
    assert!(
        std::fs::read(&status).is_err(),
        "run this test as a user that file modes bind (not root)"
    );
    let got = estimate(&d, 16).ranges["02"].source;
    std::fs::set_permissions(&status, std::fs::Permissions::from_mode(0o644)).unwrap();
    assert_eq!(got, Source::Docs);
}

// --- what the Python modules and nextflow.config must still say ---------------------------

/// `trace.fields` of nextflow.config.
fn trace_fields(config: &str) -> Option<Vec<String>> {
    let start = config.find("\ntrace {\n")?;
    let block = &config[start..];
    let block = &block[..block.find("\n}")?];
    let line = block
        .lines()
        .find(|l| l.trim_start().starts_with("fields"))?;
    let value = line.split_once('\'')?.1.split_once('\'')?.0;
    Some(value.split(',').map(str::to_owned).collect())
}

/// Nextflow's 14 default trace columns, then the 7 the history reads.
const TRACE_FIELDS: [&str; 21] = [
    "task_id",
    "hash",
    "native_id",
    "name",
    "status",
    "exit",
    "submit",
    "duration",
    "realtime",
    "%cpu",
    "peak_rss",
    "peak_vmem",
    "rchar",
    "wchar",
    "process",
    "tag",
    "attempt",
    "cpus",
    "memory",
    "start",
    "complete",
];

#[test]
fn nextflow_config_traces_the_14_default_columns_then_the_7_added() {
    let cfg = std::fs::read_to_string(repo().join("nextflow.config")).unwrap();
    assert_eq!(
        trace_fields(&cfg),
        Some(TRACE_FIELDS.map(str::to_owned).to_vec())
    );
    // Control: the same reading of a config without memory must differ.
    let broken = cfg.replace(",memory,start", ",start");
    assert_ne!(broken, cfg, "the control found no trace.fields to break");
    assert_ne!(
        trace_fields(&broken),
        Some(TRACE_FIELDS.map(str::to_owned).to_vec())
    );
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
        planted_eta_old_start_is_submit: "eta-old-start-is-submit" => super::an_old_and_a_new_trace_give_the_same_spans;
        planted_eta_memory_decimal_gb: "eta-memory-decimal-gb" => super::human_and_raw_units_read_the_same;
        planted_eta_forced_ok: "eta-forced-ok" => super::cached_failed_and_unknown_process_tasks_never_enter_the_history;
        planted_eta_ingest_adds_again: "eta-ingest-adds-again" => super::ingesting_twice_or_a_copied_trace_adds_nothing;
        planted_eta_box_same: "eta-box-same" => super::two_boxes_write_two_files;
        planted_eta_history_path_in_name: "eta-history-path-in-name" => super::the_history_holds_no_sample_id_and_no_slash;
        planted_eta_history_sample_in_process: "eta-history-sample-in-process" => super::the_history_holds_no_sample_id_and_no_slash;
        planted_eta_box_raw: "eta-box-raw" => super::the_box_name_keeps_only_letters_digits_underscore_and_dash;
        planted_eta_last_runs_2: "eta-last-runs-2" => super::three_runs_give_min_max_with_n_3;
        planted_eta_last_runs_6: "eta-last-runs-6" => super::only_the_last_5_runs_count;
        planted_eta_docs_threads_8: "eta-docs-threads-8" => super::an_8_core_box_takes_the_measured_low_end_and_doubles_the_high_end;
        planted_eta_docs_scale_down: "eta-docs-scale-down" => super::a_32_core_box_is_not_scaled_below_the_docs_hours;
        planted_eta_no_hours_low: "eta-no-hours-low" => super::on_16_and_32_cores_every_measured_step_has_low_below_high;
        planted_eta_estimates_count: "eta-estimates-count" => super::a_step_the_table_only_estimates_gets_no_number;
        planted_eta_finished_always: "eta-finished-always" => super::a_step_with_a_task_that_failed_or_was_stopped_gives_no_span;
        planted_eta_chain_ignores_upstream: "eta-chain-ignores-upstream" => super::the_hours_left_are_the_longest_chain;
        planted_eta_no_cpu_work: "eta-no-cpu-work" => super::the_cpu_work_bound_when_it_is_larger;
        planted_eta_late_not_flagged: "eta-late-not-flagged" => super::a_running_step_past_its_high_end_makes_the_high_end_unknown_and_says_so;
        planted_eta_gone_ignored: "eta-gone-ignored" => super::a_running_step_counts_only_what_is_left_of_it;
        planted_eta_none_is_zero: "eta-none-is-zero" => super::a_step_with_no_time_lands_in_unknown_steps_never_0;
        planted_eta_finished_ignored: "eta-finished-ignored" => super::a_finished_run_has_nothing_left;
        planted_eta_box_same_estimate: "eta-box-same" => super::the_estimate_gives_the_ranges_and_writes_the_history;
        planted_eta_live_traces_none: "eta-live-traces-none" => super::without_a_sample_the_trace_of_a_run_not_finished_is_left_out;
        planted_eta_unreadable_status_ignored: "eta-unreadable-status-ignored" => super::an_unreadable_run_status_leaves_out_the_newest_trace;
    }
}
