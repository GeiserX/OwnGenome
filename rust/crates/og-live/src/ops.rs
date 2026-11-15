//! How a run is doing, on the page and in Prometheus: the live site's ops
//! card and `GET /metrics`.
//!
//! Per run, [`run_ops`] gives the hours left (a range, from
//! [`og_runlog::eta`]), the memory and disk headroom (from
//! [`og_runlog::headroom`]), when the run last showed activity, and five
//! alerts: failed, stalled, late, memory and disk. Two outputs, one source:
//! [`RunOps::to_json`] is `plugins.ops` of `state.json`, drawn by
//! `assets/live/plugins/ops.js`, and [`metrics_text`] writes the same numbers
//! in the Prometheus text format for the example rules in
//! `assets/monitoring/pgp-alerts.yml`.
//!
//! Nothing here writes: the step history is read, never updated (`og run
//! eta` adds a run to it), so the genome folder can stay read-only. The paths
//! it measures are the run folder, the run's work folder (`nextflow/work`, or
//! an absolute `workDir` from the weblog) and an absolute outdir from the
//! weblog, as this machine sees them. No path, task name or parameter reaches
//! the card or a metric; the labels are `run` (the run's folder name) and
//! `fs` (`work`, `results` or `work and results`).

use jiff::tz::TimeZone;
use og_nf::time::{ISO_Z, iso, parse_utc};
use og_runlog::eta::{self, fault};
use og_runlog::headroom::{self, BadPercent, Disk, TaskMemory};
use og_runlog::host::{Host, cpu_count};
use og_runlog::state::{RunState, State, StepState, trace_files};
use og_runlog::steps::Step;
use pycompat::json::Value;
use pycompat::value::py_str;
use std::fmt::Write as _;
use std::io;
use std::path::{Path, PathBuf};

/// A run quiet this long, with no running step inside its usual time, is
/// stalled (`--stall-minutes`).
pub const DEFAULT_STALL_MINUTES: f64 = 60.0;
/// The content type of `GET /metrics`.
pub const CONTENT_TYPE: &str = "text/plain; version=0.0.4; charset=utf-8";

/// Every metric `/metrics` can write, in the order it writes them: name,
/// type, help. `docs/alerts.md` lists the same names.
pub const METRICS: [(&str, &str, &str); 21] = [
    ("pgp_runs", "gauge", "Runs the site knows."),
    (
        "pgp_disk_warn_ratio",
        "gauge",
        "The disk warn share (--disk-warn-percent / 100).",
    ),
    (
        "pgp_host_memory_bytes",
        "gauge",
        "RAM of the machine the site runs on.",
    ),
    (
        "pgp_metrics_errors",
        "gauge",
        "Runs whose numbers could not be read in this scrape.",
    ),
    ("pgp_run_running", "gauge", "1 while the run is going."),
    (
        "pgp_run_failed",
        "gauge",
        "1 when the run ended badly: a failed step, a non-zero exit or a Nextflow error.",
    ),
    (
        "pgp_run_failed_steps",
        "gauge",
        "Steps of the run in state failed.",
    ),
    (
        "pgp_run_stalled",
        "gauge",
        "1 when the run is going but has shown no activity for --stall-minutes.",
    ),
    (
        "pgp_run_seconds_since_activity",
        "gauge",
        "Seconds since the run's last message or file change.",
    ),
    (
        "pgp_run_late_steps",
        "gauge",
        "Running steps already past the high end of their usual time.",
    ),
    (
        "pgp_run_unknown_steps",
        "gauge",
        "Steps left with no usual time.",
    ),
    ("pgp_run_seconds_left_low", "gauge", "Low end of the time left."),
    (
        "pgp_run_seconds_left_high",
        "gauge",
        "High end of the time left (absent when a step is late).",
    ),
    (
        "pgp_run_tasks_near_memory_request",
        "gauge",
        "Finished tasks whose peak_rss reached the memory warn share of their request.",
    ),
    (
        "pgp_run_memory_retries",
        "gauge",
        "Task attempts retried after a memory or time kill.",
    ),
    ("pgp_run_max_memory_bytes", "gauge", "The run's max_memory."),
    (
        "pgp_run_max_memory_above_ram",
        "gauge",
        "1 when the run's max_memory is above this machine's RAM.",
    ),
    (
        "pgp_disk_used_ratio",
        "gauge",
        "Used share of the filesystem.",
    ),
    (
        "pgp_disk_inodes_used_ratio",
        "gauge",
        "Used share of the filesystem's inodes.",
    ),
    ("pgp_disk_free_bytes", "gauge", "Bytes free on the filesystem."),
    ("pgp_disk_size_bytes", "gauge", "Size of the filesystem."),
];

/// The site's ops options (also `KEY=VALUE` lines in `--config`).
#[derive(Debug, Clone, Default)]
pub struct OpsConfig {
    /// The genome folder; else the run folder's parent.
    pub root: Option<PathBuf>,
    /// The step history to read (`--eta-box`); else `PGP_BOX` or the short
    /// host name.
    pub eta_box: Option<String>,
    /// The cores of the box that runs the pipeline (`--eta-cores`); else
    /// this machine's.
    pub eta_cores: Option<usize>,
    /// `--disk-warn-percent`, as given; else `OWNGENOME_DISK_WARN_PERCENT`
    /// or 90.
    pub disk_warn_percent: Option<String>,
    /// `--mem-warn-percent`, as given; else `OWNGENOME_MEM_WARN_PERCENT` or
    /// 90.
    pub mem_warn_percent: Option<String>,
    /// `--stall-minutes`; none or 0 is [`DEFAULT_STALL_MINUTES`].
    pub stall_minutes: Option<f64>,
}

/// What [`run_ops`] reads besides the run.
pub struct Ctx<'a> {
    pub now: f64,
    pub steps: &'a [Step],
    pub host: &'a dyn Host,
    /// The zone the trace files' times are in.
    pub tz: &'a TimeZone,
}

/// How bad an alert is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Level {
    Info,
    Warning,
    Critical,
}

impl Level {
    pub fn as_str(self) -> &'static str {
        match self {
            Level::Info => "info",
            Level::Warning => "warning",
            Level::Critical => "critical",
        }
    }
}

/// What an alert is about.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlertKind {
    Failed,
    Stalled,
    Late,
    Memory,
    Disk,
}

impl AlertKind {
    pub fn as_str(self) -> &'static str {
        match self {
            AlertKind::Failed => "failed",
            AlertKind::Stalled => "stalled",
            AlertKind::Late => "late",
            AlertKind::Memory => "memory",
            AlertKind::Disk => "disk",
        }
    }
}

/// One alert of the card.
#[derive(Debug, Clone, PartialEq)]
pub struct Alert {
    pub kind: AlertKind,
    pub level: Level,
    pub text: String,
}

/// Whether the run is going, and how it went.
#[derive(Debug, Clone, PartialEq)]
pub struct RunHealth {
    pub state: RunState,
    pub failed: bool,
    pub failed_steps: Vec<String>,
    /// The exit of the last `pgp_run_finished` message, as it came (null
    /// when there is none).
    pub exit: Value,
    pub stalled: bool,
    pub seconds_since_activity: Option<f64>,
    pub stall_seconds: f64,
}

/// The hours left of a run going.
#[derive(Debug, Clone, PartialEq)]
pub struct HoursLeft {
    pub low_s: f64,
    /// `None` when a running step is past its usual high end.
    pub high_s: Option<f64>,
    pub late_steps: Vec<String>,
    pub unknown_steps: Vec<String>,
    pub basis: String,
    pub cores: usize,
}

/// The run's memory: the machine against `max_memory`, and the tasks near
/// their request.
#[derive(Debug, Clone, PartialEq)]
pub struct Memory {
    pub host_bytes: Option<u128>,
    pub max_memory_bytes: Option<i128>,
    pub max_memory_above_ram: bool,
    pub warn_ratio: f64,
    pub near_request: usize,
    pub retries: usize,
    /// The tasks at or above the warn share of their request.
    pub tasks: Vec<TaskMemory>,
}

/// The ops card of one run.
#[derive(Debug, Clone, PartialEq)]
pub struct RunOps {
    pub run: RunHealth,
    /// Only while the run is going.
    pub hours_left: Option<HoursLeft>,
    pub memory: Memory,
    pub disk: Vec<Disk>,
    pub disk_warn_ratio: f64,
    pub alerts: Vec<Alert>,
}

/// The card, or the one line it shows instead.
#[derive(Debug, Clone, PartialEq)]
pub enum OpsState {
    Ok(Box<RunOps>),
    Error(String),
}

/// The messages of the newest run: after the last `pgp_run_finished` that
/// is not the newest message, and from the last weblog `started`.
pub fn current_run(events: &[Value]) -> &[Value] {
    let mut start = 0;
    for (i, ev) in events.iter().enumerate().take(events.len().saturating_sub(1)) {
        if is_event(ev, "pgp_run_finished") {
            start = i + 1;
        }
    }
    if let Some(i) = (start..events.len())
        .rev()
        .find(|&i| is_event(&events[i], "started"))
    {
        start = i;
    }
    &events[start..]
}

fn field<'v>(ev: &'v Value, k: &str) -> Option<&'v Value> {
    match ev {
        Value::Object(o) => o.get(k),
        _ => None,
    }
}

fn is_event(ev: &Value, kind: &str) -> bool {
    matches!(field(ev, "event"), Some(Value::Str(e)) if e == kind)
}

/// An ISO time such as `2026-01-02T10:00:00Z`, in seconds.
fn epoch_of(v: Option<&Value>) -> Option<f64> {
    match v {
        Some(Value::Str(s)) => parse_utc(s, &[ISO_Z]),
        _ => None,
    }
}

fn is_dir(p: Option<&Path>) -> Option<&Path> {
    p.filter(|p| p.is_dir())
}

/// An absolute path from the weblog, else nothing.
fn abs_dir(v: Option<Value>) -> Option<String> {
    match v {
        Some(Value::Str(s)) if s.starts_with('/') && !fault::on("ops-bare-run-no-work") => Some(s),
        _ => None,
    }
}

/// The newest sign of life: a message's `utcTime`, the run's
/// `logs/run_status.tsv`, or a trace file of the run.
pub fn last_activity(run_dir: Option<&Path>, events: &[Value], started: Option<f64>) -> Option<f64> {
    let mut times: Vec<f64> = events
        .iter()
        .filter_map(|ev| epoch_of(field(ev, "utcTime")))
        .collect();
    if let Some(dir) = is_dir(run_dir) {
        times.extend(pycompat::mtime::mtime(dir.join("logs").join("run_status.tsv")).ok());
        times.extend(
            trace_files(dir, started)
                .iter()
                .filter_map(|p| pycompat::mtime::mtime(p).ok()),
        );
    }
    times.into_iter().reduce(f64::max)
}

/// `(work, results)` to measure: the run folder's `nextflow/work`, else an
/// absolute `workDir` from the weblog; the run folder, else an absolute
/// outdir from the weblog (a relative one is under the weblog's
/// `launchDir`).
pub fn folders(run_dir: Option<&Path>, events: &[Value]) -> (Option<String>, Option<String>) {
    let (m_work, m_out, _) = headroom::metadata(events);
    let have = is_dir(run_dir);
    let text = |p: &Path| p.to_string_lossy().into_owned();
    let work = match have.map(|d| d.join("nextflow").join("work")) {
        Some(w) if w.is_dir() => Some(text(&w)),
        _ => abs_dir(m_work),
    };
    let results = match have {
        Some(d) => Some(text(d)),
        None => match m_out {
            Some(Value::Str(s)) if s.starts_with('/') => Some(s),
            _ => None,
        },
    };
    (work, results)
}

/// The run's started time as the history reads it (whole seconds of its ISO
/// form).
fn started_of(state: &State<'_>) -> Option<f64> {
    state
        .started
        .and_then(iso)
        .and_then(|v| parse_utc(&v, &[ISO_Z]))
}

/// This run's own tasks from its trace files.
fn trace_rows(run_dir: Option<&Path>, sample: &str, started: Option<f64>, tz: &TimeZone) -> Vec<eta::TraceTask> {
    let Some(dir) = is_dir(run_dir) else {
        return Vec::new();
    };
    trace_files(dir, started)
        .iter()
        .flat_map(|p| eta::read_trace(p, tz))
        .filter(|r| og_runlog::tasks::is_sample(r.tag.as_deref(), sample))
        .collect()
}

/// A weblog value as text (`max_memory` may come as a number).
fn value_text(v: &Value) -> Option<String> {
    match v {
        Value::Null => None,
        Value::Str(s) => Some(s.clone()),
        v => py_str(v).ok(),
    }
}

/// An exit that means the run ended well, or none: null, 0 (or 0.0, false).
fn exit_ok(v: &Value) -> bool {
    match v {
        Value::Null => true,
        Value::Int(i) => *i == 0,
        Value::Float(x) => *x == 0.0,
        Value::Bool(b) => !b,
        _ => false,
    }
}

fn plural(n: usize, one: &str, many: &str) -> String {
    if n > 1 { many.to_owned() } else { one.to_owned() }
}

/// The ops card of the run in `run_dir` (none or missing for a run known
/// from its messages only), whose state is `state` and whose stored
/// messages are `events`. A warn percent outside 1 to 100 is the error.
pub fn run_ops(
    run_dir: Option<&Path>,
    state: &State<'_>,
    events: &[Value],
    cfg: &OpsConfig,
    ctx: &Ctx<'_>,
) -> Result<RunOps, BadPercent> {
    let now = ctx.now;
    let root = cfg
        .root
        .clone()
        .or_else(|| run_dir.and_then(Path::parent).map(Path::to_path_buf));
    let cur = current_run(events);
    let sample = if state.sample.is_empty() {
        run_dir
            .and_then(Path::file_name)
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default()
    } else {
        state.sample.clone()
    };
    // run-all.sh posts pgp_run_finished as its last message; the site marks
    // the run finished on it too.
    let ended = cur.last().is_some_and(|ev| is_event(ev, "pgp_run_finished"));
    let run_st = if ended && state.run == RunState::Running {
        RunState::Finished
    } else {
        state.run
    };
    let running = run_st == RunState::Running;
    let started = started_of(state);

    // Hours left: the ranges from the history as it is, without adding to it.
    let cores = cfg.eta_cores.filter(|c| *c != 0).unwrap_or_else(cpu_count);
    let history = root
        .as_deref()
        .map(|r| eta::read_history(&eta::history_path(r, &eta::box_name(cfg.eta_box.as_deref()))))
        .unwrap_or_default();
    let mut exclude = root.as_deref().map(eta::live_traces).unwrap_or_default();
    if let Some(dir) = is_dir(run_dir)
        && state.run == RunState::Running
    {
        exclude.extend(trace_files(dir, started).iter().filter_map(|p| {
            p.file_name().map(|n| n.to_string_lossy().into_owned())
        }));
    }
    let ranges = eta::step_ranges(&history, ctx.steps, cores, &exclude);
    let rem = running.then(|| eta::remaining(state, &ranges, now, cores));

    // Headroom on the folders found here. The weblog metadata is taken out
    // of the messages first, so the snapshot cannot fall back to a weblog
    // path that folders() left out.
    let (work, results) = folders(run_dir, cur);
    let (_, _, cap) = headroom::metadata(cur);
    let plain: Vec<Value> = cur
        .iter()
        .map(|ev| match ev {
            Value::Object(o) => Value::Object(
                o.iter()
                    .filter(|(k, _)| *k != "metadata")
                    .map(|(k, v)| (k.clone(), v.clone()))
                    .collect(),
            ),
            v => v.clone(),
        })
        .collect();
    let rows = trace_rows(run_dir, &sample, started, ctx.tz);
    let cap = cap.as_ref().and_then(value_text);
    let snap = headroom::snapshot(
        ctx.host,
        ctx.steps,
        &headroom::Request {
            events: &plain,
            trace_rows: &rows,
            work_dir: work.as_deref(),
            outdir: results.as_deref(),
            max_memory: cap.as_deref(),
            disk_warn_percent: cfg.disk_warn_percent.as_deref(),
            mem_warn_percent: cfg.mem_warn_percent.as_deref(),
        },
    )?;

    // What went wrong or is going wrong.
    let failed_steps: Vec<String> = state
        .steps
        .iter()
        .filter(|s| s.state == StepState::Failed)
        .map(|s| s.step.step.clone())
        .collect();
    let exit = cur
        .iter()
        .rfind(|ev| is_event(ev, "pgp_run_finished"))
        .map_or(Value::Null, |ev| field(ev, "exit").cloned().unwrap_or(Value::Null));
    let error = cur.iter().any(|ev| is_event(ev, "error"));
    let finished = run_st == RunState::Finished;
    let bad_exit = !exit_ok(&exit) && !fault::on("ops-exit-ignored");
    let failed = !failed_steps.is_empty() || (finished && (bad_exit || error));
    let quiet = last_activity(run_dir, cur, started).map(|a| (now - a).max(0.0));
    let stall_s = cfg.stall_minutes.filter(|m| *m != 0.0).unwrap_or(DEFAULT_STALL_MINUTES) * 60.0;
    // Running steps still inside their usual high end: a long task posts
    // nothing meanwhile.
    let covered = state.steps.iter().any(|s| {
        let high = ranges.get(&s.step.step).and_then(|r| r.high_s);
        let st = s.started_utc.as_deref().and_then(|v| parse_utc(v, &[ISO_Z]));
        s.state == StepState::Running
            && matches!((high, st), (Some(h), Some(st)) if now - st < h)
            && !fault::on("ops-stall-ignores-steps")
    });
    // A run with no end record that the state calls stopped has been quiet
    // even longer: still stalled, so the alert holds for a run that died
    // without an end instead of resolving.
    let alive = run_st == RunState::Running
        || (run_st == RunState::Stopped && !fault::on("ops-stopped-not-stalled"));
    let stalled = alive
        && quiet.is_some_and(|q| q >= stall_s || fault::on("ops-stall-ignores-quiet"))
        && !covered;
    let late: Vec<String> = match &rem {
        Some(r) if !fault::on("ops-late-dropped") => r.late_steps.clone(),
        _ => Vec::new(),
    };
    let mem_hints: Vec<&TaskMemory> = snap.tasks.iter().filter(|t| t.warn).collect();
    let retries: usize = if fault::on("ops-retries-uncounted") {
        0
    } else {
        snap.tasks.iter().map(|t| t.oom_retry).sum()
    };
    let disk_warn: Vec<&Disk> = snap
        .disk
        .iter()
        .filter(|d| d.warn && !fault::on("ops-disk-warn-ignored"))
        .collect();

    let mut alerts = Vec::new();
    if failed {
        let text = if !failed_steps.is_empty() {
            format!(
                "{} {} failed",
                plural(failed_steps.len(), "step", "steps"),
                failed_steps.join(", ")
            )
        } else if !exit_ok(&exit) {
            format!("the run ended with exit {}", py_str(&exit).unwrap_or_default())
        } else {
            "Nextflow reported an error".to_owned()
        };
        alerts.push(Alert {
            kind: AlertKind::Failed,
            level: Level::Critical,
            text,
        });
    }
    if stalled {
        let minutes = pycompat::fmt::fixed(quiet.unwrap_or(0.0) / 60.0, 0);
        alerts.push(Alert {
            kind: AlertKind::Stalled,
            level: Level::Warning,
            text: format!(
                "no activity for {minutes} minutes, and no running step is inside its usual time"
            ),
        });
    }
    if !late.is_empty() {
        let text = late
            .iter()
            .map(|k| format!("step {k} is taking longer than usual"))
            .collect::<Vec<_>>()
            .join("; ");
        alerts.push(Alert {
            kind: AlertKind::Late,
            level: Level::Warning,
            text,
        });
    }
    if retries > 0 || snap.host.warn {
        let mut parts = Vec::new();
        if retries > 0 {
            parts.push(format!(
                "{retries} task {} retried after a memory or time kill",
                plural(retries, "attempt", "attempts")
            ));
        }
        if snap.host.warn {
            parts.push("max_memory is above this machine's memory".to_owned());
        }
        alerts.push(Alert {
            kind: AlertKind::Memory,
            level: Level::Warning,
            text: parts.join("; "),
        });
    } else if !mem_hints.is_empty() {
        alerts.push(Alert {
            kind: AlertKind::Memory,
            level: Level::Info,
            text: format!(
                "{} {} reached {}% of the memory asked for (peak_rss counts shared pages, so this may be fine)",
                mem_hints.len(),
                plural(mem_hints.len(), "task", "tasks"),
                pycompat::fmt::fixed(snap.mem_warn_ratio * 100.0, 0)
            ),
        });
    }
    for d in &disk_warn {
        alerts.push(Alert {
            kind: AlertKind::Disk,
            level: if d.used_ratio.unwrap_or(0.0) >= 0.98 {
                Level::Critical
            } else {
                Level::Warning
            },
            text: format!(
                "disk {} is at or above {}% used",
                d.fs,
                pycompat::fmt::fixed(snap.warn_used_ratio * 100.0, 0)
            ),
        });
    }

    Ok(RunOps {
        run: RunHealth {
            state: run_st,
            failed,
            failed_steps,
            exit,
            stalled,
            seconds_since_activity: quiet,
            stall_seconds: stall_s,
        },
        hours_left: rem.map(|r| HoursLeft {
            low_s: r.low_s,
            high_s: r.high_s,
            late_steps: late,
            unknown_steps: r.unknown_steps,
            basis: r.basis,
            cores,
        }),
        memory: Memory {
            host_bytes: snap.host.memory_bytes,
            max_memory_bytes: snap.host.max_memory_bytes,
            max_memory_above_ram: snap.host.warn,
            warn_ratio: snap.mem_warn_ratio,
            near_request: mem_hints.len(),
            retries,
            tasks: mem_hints.into_iter().cloned().collect(),
        },
        disk: snap.disk,
        disk_warn_ratio: snap.warn_used_ratio,
        alerts,
    })
}

/// The state hook: the card, or the line it shows instead. Only a warn
/// percent outside 1 to 100 fails it (the readers take a missing or
/// half-written file as no value), and the rest of the run's state stays.
pub fn hook(run_dir: Option<&Path>, state: &State<'_>, events: &[Value], cfg: &OpsConfig, ctx: &Ctx<'_>) -> OpsState {
    match run_ops(run_dir, state, events, cfg, ctx) {
        Ok(o) => OpsState::Ok(Box::new(o)),
        Err(BadPercent(m)) => OpsState::Error(m.chars().take(200).collect()),
    }
}

fn int(v: impl TryInto<i128>) -> Value {
    v.try_into().map_or(Value::Null, Value::Int)
}

fn opt_int(v: Option<impl TryInto<i128>>) -> Value {
    v.map_or(Value::Null, int)
}

fn opt_float(v: Option<f64>) -> Value {
    v.map_or(Value::Null, Value::Float)
}

fn strs(v: &[String]) -> Value {
    Value::Array(v.iter().map(|s| Value::Str(s.clone())).collect())
}

fn obj<const N: usize>(pairs: [(&str, Value); N]) -> Value {
    Value::Object(pairs.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}

impl RunOps {
    /// `plugins.ops` of `state.json`, key for key.
    pub fn to_json(&self) -> Value {
        let r = &self.run;
        let m = &self.memory;
        obj([
            (
                "run",
                obj([
                    ("state", Value::Str(r.state.as_str().to_owned())),
                    ("failed", Value::Bool(r.failed)),
                    ("failed_steps", strs(&r.failed_steps)),
                    ("exit", r.exit.clone()),
                    ("stalled", Value::Bool(r.stalled)),
                    ("seconds_since_activity", opt_float(r.seconds_since_activity)),
                    ("stall_seconds", Value::Float(r.stall_seconds)),
                ]),
            ),
            (
                "hours_left",
                self.hours_left.as_ref().map_or(Value::Null, |h| {
                    obj([
                        ("low_s", Value::Float(h.low_s)),
                        ("high_s", opt_float(h.high_s)),
                        ("late_steps", strs(&h.late_steps)),
                        ("unknown_steps", strs(&h.unknown_steps)),
                        ("basis", Value::Str(h.basis.clone())),
                        ("cores", int(h.cores)),
                    ])
                }),
            ),
            (
                "memory",
                obj([
                    ("host_bytes", opt_int(m.host_bytes)),
                    ("max_memory_bytes", opt_int(m.max_memory_bytes)),
                    ("max_memory_above_ram", Value::Bool(m.max_memory_above_ram)),
                    ("warn_ratio", Value::Float(m.warn_ratio)),
                    ("near_request", int(m.near_request)),
                    ("retries", int(m.retries)),
                    (
                        "tasks",
                        Value::Array(
                            m.tasks
                                .iter()
                                .map(|t| {
                                    obj([
                                        ("step", t.step.clone().map_or(Value::Null, Value::Str)),
                                        ("process", Value::Str(t.process.clone())),
                                        ("state", Value::Str(t.state.clone())),
                                        ("request_bytes", opt_int(t.request_bytes)),
                                        ("peak_rss_bytes", opt_int(t.peak_rss_bytes)),
                                        ("ratio", opt_float(t.ratio)),
                                        ("warn", Value::Bool(t.warn)),
                                        ("oom_retry", int(t.oom_retry)),
                                    ])
                                })
                                .collect(),
                        ),
                    ),
                ]),
            ),
            (
                "disk",
                Value::Array(
                    self.disk
                        .iter()
                        .map(|d| {
                            obj([
                                ("fs", Value::Str(d.fs.clone())),
                                ("used_ratio", opt_float(d.used_ratio)),
                                ("free_bytes", int(d.free_bytes)),
                                ("size_bytes", int(d.size_bytes)),
                                ("inodes_used_ratio", opt_float(d.inodes_used_ratio)),
                                ("missing", Value::Bool(d.missing)),
                                ("warn", Value::Bool(d.warn)),
                            ])
                        })
                        .collect(),
                ),
            ),
            ("disk_warn_ratio", Value::Float(self.disk_warn_ratio)),
            (
                "alerts",
                Value::Array(
                    self.alerts
                        .iter()
                        .map(|a| {
                            obj([
                                ("kind", Value::Str(a.kind.as_str().to_owned())),
                                ("level", Value::Str(a.level.as_str().to_owned())),
                                ("text", Value::Str(a.text.clone())),
                            ])
                        })
                        .collect(),
                ),
            ),
        ])
    }
}

impl OpsState {
    /// `plugins.ops` of `state.json`: the card, or `{"error": ...}`.
    pub fn to_json(&self) -> Value {
        match self {
            OpsState::Ok(o) => o.to_json(),
            OpsState::Error(e) => obj([("error", Value::Str(e.clone()))]),
        }
    }
}

/// A number as `/metrics` writes it: a flag as 1 or 0, a whole number as
/// itself, any other as Python's `repr` (`3600.0`, `1e-05`).
#[derive(Debug, Clone, Copy)]
enum Num {
    Flag(bool),
    Whole(i128),
    Real(f64),
}

impl std::fmt::Display for Num {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Num::Flag(b) => f.write_str(if *b { "1" } else { "0" }),
            Num::Whole(i) => write!(f, "{i}"),
            Num::Real(x) => f.write_str(&pycompat::fmt::repr(*x)),
        }
    }
}

fn whole(v: impl TryInto<i128>) -> Option<Num> {
    v.try_into().ok().map(Num::Whole)
}

/// A label value, escaped as the text format asks.
fn label(v: &str) -> String {
    v.replace('\\', "\\\\")
        .replace('\n', "\\n")
        .replace('"', "\\\"")
}

/// The disk warn share `/metrics` reports: `--disk-warn-percent` (else
/// `OWNGENOME_DISK_WARN_PERCENT`, else 90) over 100; none when it is not a
/// number from 1 to 100.
pub fn disk_warn_ratio(cfg: &OpsConfig) -> Option<f64> {
    headroom::warn_percent(cfg.disk_warn_percent.as_deref(), Some("OWNGENOME_DISK_WARN_PERCENT"))
        .ok()
        .map(|p| p / 100.0)
}

/// The Prometheus text of `runs` (run name, card). A value that is not
/// known is left out, never written as 0.
pub fn metrics_text(
    runs: &[(&str, &RunOps)],
    disk_warn_ratio: Option<f64>,
    host_bytes: Option<u128>,
    errors: usize,
) -> String {
    let mut samples: Vec<Vec<String>> = vec![Vec::new(); METRICS.len()];
    let mut put = |name: &str, value: Option<Num>, labels: &[(&str, &str)]| {
        let Some(value) = value else { return };
        let i = METRICS
            .iter()
            .position(|(n, _, _)| *n == name)
            .expect("every metric written is in METRICS");
        let lab = labels
            .iter()
            .map(|(k, v)| format!("{k}=\"{}\"", label(v)))
            .collect::<Vec<_>>()
            .join(",");
        samples[i].push(if lab.is_empty() {
            format!("{name} {value}")
        } else {
            format!("{name}{{{lab}}} {value}")
        });
    };
    put("pgp_runs", whole(runs.len()), &[]);
    put("pgp_disk_warn_ratio", disk_warn_ratio.map(Num::Real), &[]);
    put("pgp_host_memory_bytes", host_bytes.and_then(whole), &[]);
    put("pgp_metrics_errors", whole(errors), &[]);
    for (name, o) in runs {
        let run = [("run", *name)];
        let (r, m) = (&o.run, &o.memory);
        put("pgp_run_running", Some(Num::Flag(r.state == RunState::Running)), &run);
        put("pgp_run_failed", Some(Num::Flag(r.failed)), &run);
        put("pgp_run_failed_steps", whole(r.failed_steps.len()), &run);
        put("pgp_run_stalled", Some(Num::Flag(r.stalled)), &run);
        put("pgp_run_seconds_since_activity", r.seconds_since_activity.map(Num::Real), &run);
        if let Some(h) = &o.hours_left {
            put("pgp_run_late_steps", whole(h.late_steps.len()), &run);
            put("pgp_run_unknown_steps", whole(h.unknown_steps.len()), &run);
            put("pgp_run_seconds_left_low", Some(Num::Real(h.low_s)), &run);
            put("pgp_run_seconds_left_high", h.high_s.map(Num::Real), &run);
        }
        put("pgp_run_tasks_near_memory_request", whole(m.near_request), &run);
        put("pgp_run_memory_retries", whole(m.retries), &run);
        put("pgp_run_max_memory_bytes", m.max_memory_bytes.and_then(whole), &run);
        if m.max_memory_bytes.is_some() {
            put("pgp_run_max_memory_above_ram", Some(Num::Flag(m.max_memory_above_ram)), &run);
        }
        for d in &o.disk {
            let at = [("run", *name), ("fs", d.fs.as_str())];
            put("pgp_disk_used_ratio", d.used_ratio.map(Num::Real), &at);
            put("pgp_disk_inodes_used_ratio", d.inodes_used_ratio.map(Num::Real), &at);
            put("pgp_disk_free_bytes", whole(d.free_bytes), &at);
            put("pgp_disk_size_bytes", whole(d.size_bytes), &at);
        }
    }
    let mut out = String::new();
    for ((name, typ, help), lines) in METRICS.iter().zip(&samples) {
        if lines.is_empty() {
            continue;
        }
        let _ = writeln!(out, "# HELP {name} {help}");
        let _ = writeln!(out, "# TYPE {name} {typ}");
        for l in lines {
            out.push_str(l);
            out.push('\n');
        }
    }
    if out.is_empty() {
        out.push('\n');
    }
    out
}

/// The token in `path` (`--metrics-token-file`), surrounding white space
/// cut; none for no file given or an empty one.
pub fn read_token(path: Option<&Path>) -> io::Result<Option<String>> {
    let Some(path) = path else { return Ok(None) };
    let text = std::fs::read_to_string(path)?;
    let tok = pycompat::text::strip(&text, None);
    Ok((!tok.is_empty()).then(|| tok.to_owned()))
}

/// Who asks for `/metrics`.
#[derive(Debug, Clone, Copy, Default)]
pub struct Asker<'a> {
    /// The site has a password.
    pub password_set: bool,
    /// The request carries a valid session.
    pub logged_in: bool,
    /// Its `Authorization` header.
    pub authorization: Option<&'a str>,
}

/// Equal bytes, in time that does not depend on where they differ.
fn same_secret(a: &[u8], b: &[u8]) -> bool {
    a.len() == b.len() && a.iter().zip(b).fold(0u8, |acc, (x, y)| acc | (x ^ y)) == 0
}

/// Whether `/metrics` answers `who`. Without a password the site is
/// loopback only and `/metrics` is open, like every page. With one: a
/// logged-in browser, or `Authorization: Bearer <token>`. `token` reads the
/// token, and is called only when one is needed (an error: the token file
/// cannot be read).
pub fn metrics_allowed(
    who: &Asker<'_>,
    token: impl FnOnce() -> io::Result<Option<String>>,
) -> io::Result<bool> {
    if !who.password_set || who.logged_in {
        return Ok(true);
    }
    let auth = who.authorization.unwrap_or("");
    let (scheme, given) = auth.split_once(' ').unwrap_or((auth, ""));
    let Some(token) = token()? else {
        return Ok(false);
    };
    Ok(pycompat::text::lower(scheme) == "bearer"
        && same_secret(
            pycompat::text::strip(given, None).as_bytes(),
            token.as_bytes(),
        ))
}
