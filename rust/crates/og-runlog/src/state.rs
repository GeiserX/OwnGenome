//! Where one sample's run stands, step by step: `run_state.json`.
//!
//! Per step, run-all.sh's last row decides (`ok`, `failed`, `skipped
//! (reason)`, `not finished`). A `planned` step, or one with no row, takes
//! its state from main.nf's process table when the run has ended, else from
//! its tasks so far; a planned step with neither is waiting (or not
//! finished once the run has ended), a step with no row and no task is not
//! planned. A waiting or running step whose summary section was written
//! before this run started carries `earlier_result`: that file is from an
//! earlier run and is being redone.

use crate::files::{self, ProcessCounts, RunStatus};
use crate::steps::{Step, stale_after};
use crate::tasks::{self, Task, done_bad};
use crate::util::{Obj, list, opt_f, opt_s, py_max, py_min, s};
use indexmap::IndexMap;
use jiff::tz::TimeZone;
use og_nf::time::{RUN_ALL_UTC, iso, parse_iso8601, parse_utc};
use pycompat::json::Value;
use std::io;
use std::path::{Path, PathBuf};

pub const SCHEMA_VERSION: i128 = 1;

/// A step's state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepState {
    Waiting,
    Running,
    Done,
    Failed,
    Skipped,
    NotFinished,
    NotPlanned,
}

impl StepState {
    pub fn as_str(self) -> &'static str {
        match self {
            StepState::Waiting => "waiting",
            StepState::Running => "running",
            StepState::Done => "done",
            StepState::Failed => "failed",
            StepState::Skipped => "skipped",
            StepState::NotFinished => "not finished",
            StepState::NotPlanned => "not planned",
        }
    }
}

/// The run's state. `Stopped`: no end record, and nothing heard from it for
/// longer than [`stale_after`]; only when the caller asks for it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunState {
    NotStarted,
    Running,
    Stopped,
    Finished,
}

impl RunState {
    pub fn as_str(self) -> &'static str {
        match self {
            RunState::NotStarted => "not started",
            RunState::Running => "running",
            RunState::Stopped => "stopped",
            RunState::Finished => "finished",
        }
    }
}

/// How many of a process's tasks are in each state.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Counts {
    pub waiting: i128,
    pub running: i128,
    pub done: i128,
    pub cached: i128,
    pub failed: i128,
}

/// A section file an earlier run left, which this run is redoing.
#[derive(Debug, Clone, PartialEq)]
pub struct Earlier {
    pub section: String,
    pub source: String,
    pub file_date: Option<String>,
}

/// One step of the state.
#[derive(Debug, Clone, PartialEq)]
pub struct StepOut<'a> {
    pub step: &'a Step,
    pub state: StepState,
    pub reason: Option<String>,
    pub note: Option<String>,
    pub earlier: Option<Earlier>,
    pub started_utc: Option<String>,
    pub ended_utc: Option<String>,
    /// The step's processes that have tasks, in catalogue order.
    pub processes: Vec<(&'a str, Counts)>,
}

/// The whole state.
#[derive(Debug, Clone, PartialEq)]
pub struct State<'a> {
    pub sample: String,
    pub generated: f64,
    pub run: RunState,
    pub started: Option<f64>,
    /// run-all.sh's `meta finished_utc`, as ISO 8601 when it has its usual
    /// form, else as it is.
    pub finished_utc: Option<String>,
    pub declared_sex: Option<String>,
    pub status_file: bool,
    /// When the run was last heard from; only when the caller asked.
    pub last_heard: Option<f64>,
    pub steps: Vec<StepOut<'a>>,
}

/// What [`build`] reads.
#[derive(Default)]
pub struct Inputs<'a> {
    pub status: Option<&'a RunStatus>,
    pub process_status: Option<&'a IndexMap<String, ProcessCounts>>,
    pub trace_tasks: Vec<Task>,
    pub event_tasks: Vec<Task>,
    /// The weblog messages, for the run's own start and end.
    pub events: &'a [Value],
    /// `{section: (source, file time)}` of the last summary.
    pub sources: IndexMap<String, (String, f64)>,
    /// When the run's own files last changed. Without it a run is never
    /// called stopped.
    pub heard: Option<f64>,
}

pub fn counts_of<'t>(ts: impl IntoIterator<Item = &'t Task>) -> Counts {
    let mut c = Counts::default();
    for t in ts {
        match t.status.as_str() {
            "CACHED" => c.cached += 1,
            "COMPLETED" => c.done += 1,
            s if done_bad(s) => c.failed += 1,
            "RUNNING" => c.running += 1,
            _ => c.waiting += 1,
        }
    }
    c
}

fn state_from_tasks(ts: &[&Task]) -> Option<StepState> {
    if ts.is_empty() {
        return None;
    }
    let c = counts_of(ts.iter().copied());
    Some(if c.failed > 0 {
        StepState::Failed
    } else if c.running > 0 || (c.waiting > 0 && (c.done > 0 || c.cached > 0)) {
        StepState::Running
    } else if c.waiting > 0 {
        StepState::Waiting
    } else {
        StepState::Done
    })
}

/// HAPLOCHECK gets no task when step 20 did not run or failed, as run-all.sh
/// reads it.
fn no_haplocheck_task_expected(ps: &IndexMap<String, ProcessCounts>) -> bool {
    ps.get("MITO_VARIANTS")
        .is_none_or(|mv| mv.failed + mv.ignored > 0)
}

/// The step from main.nf's end-of-run table, `rows` being its processes
/// that have a row. A process with no task had not started, so the step is
/// not finished, except HAPLOCHECK when step 20 gave it nothing to check.
fn state_from_process_status(
    rows: &[(&str, &ProcessCounts)],
    all: &IndexMap<String, ProcessCounts>,
) -> StepState {
    if rows.iter().any(|(_, r)| r.failed + r.ignored != 0) {
        return StepState::Failed;
    }
    if rows.iter().any(|(_, r)| r.unfinished != 0) {
        return StepState::NotFinished;
    }
    for (p, r) in rows {
        let excused = *p == "HAPLOCHECK" && no_haplocheck_task_expected(all);
        if !(r.any() || excused) {
            return StepState::NotFinished;
        }
    }
    if rows.iter().any(|(_, r)| r.succeeded + r.cached != 0) {
        StepState::Done
    } else {
        StepState::NotFinished
    }
}

/// What run-all.sh's last row of a step says.
enum Row {
    Planned,
    Is(StepState, Option<String>),
    Unknown,
}

fn from_row(row: Option<&str>) -> Row {
    match row {
        None => Row::Unknown,
        Some("planned") => Row::Planned,
        Some("ok") => Row::Is(StepState::Done, None),
        Some(r) if r.starts_with("failed") => Row::Is(StepState::Failed, None),
        Some("not finished") => Row::Is(StepState::NotFinished, None),
        Some("skipped") => Row::Is(StepState::Skipped, None),
        Some(r) => match r
            .strip_prefix("skipped (")
            .and_then(|x| x.strip_suffix(')'))
            .filter(|x| !x.contains('\n'))
        {
            Some(reason) => Row::Is(StepState::Skipped, Some(reason.to_owned())),
            None => Row::Unknown,
        },
    }
}

/// The run state at `now`, from what the run files hold.
pub fn build<'a>(
    sample: &str,
    steps: &'a [Step],
    inp: &Inputs<'_>,
    now: f64,
    tz: &TimeZone,
) -> State<'a> {
    let empty = RunStatus::default();
    let status = inp.status.unwrap_or(&empty);
    let started = status.started();
    let by_events = inp.events.iter().any(|ev| {
        matches!(ev, Value::Object(o) if matches!(o.get("event"), Some(Value::Str(e)) if e == "completed" || e == "error"))
    });
    let finished = status.meta.contains_key("finished_utc")
        || (inp.process_status.is_some() && inp.status.is_none())
        || by_events;
    let tasks = tasks::latest_tasks(inp.trace_tasks.clone(), inp.event_tasks.clone());
    let mut out = Vec::with_capacity(steps.len());
    for s in steps {
        let mine: Vec<&Task> = tasks
            .iter()
            .filter(|t| s.processes.contains(&t.process))
            .collect();
        let (mut state, reason) = match from_row(status.steps.get(&s.step).map(String::as_str)) {
            Row::Is(st, r) => (Some(st), r),
            Row::Planned | Row::Unknown => (None, None),
        };
        if state.is_none() {
            let planned = matches!(
                from_row(status.steps.get(&s.step).map(String::as_str)),
                Row::Planned
            );
            let ps: Vec<(&str, &ProcessCounts)> = inp.process_status.map_or_else(Vec::new, |all| {
                s.processes
                    .iter()
                    .filter_map(|p| all.get(p).map(|c| (p.as_str(), c)))
                    .collect()
            });
            let live = match inp.process_status {
                Some(all) if !ps.is_empty() => Some(state_from_process_status(&ps, all)),
                _ => state_from_tasks(&mine),
            };
            state = Some(match live {
                Some(l) => l,
                None if planned && finished => StepState::NotFinished,
                None if planned => StepState::Waiting,
                None => StepState::NotPlanned,
            });
        }
        let state = state.unwrap_or(StepState::NotPlanned);
        let (mut earlier, mut note) = (None, None);
        if matches!(state, StepState::Waiting | StepState::Running)
            && let Some(start) = started
        {
            for k in &s.sections {
                if let Some((src, t)) = inp.sources.get(k)
                    && *t < start
                {
                    let date = iso(*t);
                    note = Some(format!(
                        "an earlier result ({}) is being redone",
                        date.as_deref().unwrap_or("None")
                    ));
                    earlier = Some(Earlier {
                        section: k.clone(),
                        source: src.clone(),
                        file_date: date,
                    });
                    break;
                }
            }
        }
        let starts = py_min(mine.iter().filter_map(|t| t.start));
        let ends = py_max(mine.iter().filter_map(|t| t.end));
        let processes = s
            .processes
            .iter()
            .filter(|p| mine.iter().any(|t| &t.process == *p))
            .map(|p| {
                (
                    p.as_str(),
                    counts_of(mine.iter().copied().filter(|t| &t.process == p)),
                )
            })
            .collect();
        out.push(StepOut {
            step: s,
            state,
            reason,
            note,
            earlier,
            started_utc: starts.and_then(iso),
            ended_utc: ends
                .filter(|_| matches!(state, StepState::Done | StepState::Failed))
                .and_then(iso),
            processes,
        });
    }
    let mut run = if finished {
        RunState::Finished
    } else if inp.status.is_some() || !tasks.is_empty() {
        RunState::Running
    } else {
        RunState::NotStarted
    };
    let last_heard = inp.heard.map(|h| {
        let msgs = inp.events.iter().map(|ev| match ev {
            Value::Object(o) => match o.get("utcTime") {
                Some(Value::Str(t)) => parse_iso8601(t, tz).unwrap_or(0.0),
                _ => 0.0,
            },
            _ => 0.0,
        });
        let task_times = tasks.iter().flat_map(|t| [t.start, t.end]).flatten();
        py_max(
            [h, started.filter(|x| *x != 0.0).unwrap_or(0.0)]
                .into_iter()
                .chain(msgs)
                .chain(task_times),
        )
        .unwrap_or(h)
    });
    if let Some(h) = last_heard
        && run == RunState::Running
        && now - h > stale_after(steps)
    {
        run = RunState::Stopped;
    }
    let finished_utc = status.meta.get("finished_utc").map(|v| {
        parse_utc(v, &[RUN_ALL_UTC])
            .and_then(iso)
            .unwrap_or_else(|| v.clone())
    });
    State {
        sample: sample.to_owned(),
        generated: now,
        run,
        started,
        finished_utc,
        declared_sex: status.meta.get("declared_sex").cloned(),
        status_file: inp.status.is_some(),
        last_heard,
        steps: out,
    }
}

impl State<'_> {
    /// The JSON of `run_state.json`, key for key.
    pub fn to_json(&self) -> Value {
        let mut run = Obj::new();
        run.insert("state".into(), s(self.run.as_str()));
        run.insert("started_utc".into(), opt_s(self.started.and_then(iso)));
        run.insert("finished_utc".into(), opt_s(self.finished_utc.clone()));
        run.insert("declared_sex".into(), opt_s(self.declared_sex.clone()));
        run.insert("status_file".into(), Value::Bool(self.status_file));
        if let Some(h) = self.last_heard {
            run.insert("last_heard_utc".into(), opt_s(iso(h)));
        }
        let steps = self
            .steps
            .iter()
            .map(|o| {
                let mut m = Obj::new();
                m.insert("step".into(), s(o.step.step.as_str()));
                m.insert("label".into(), s(o.step.label.as_str()));
                m.insert(
                    "tool".into(),
                    opt_s(Some(o.step.tool.as_str()).filter(|t| !t.is_empty())),
                );
                m.insert("state".into(), s(o.state.as_str()));
                m.insert("reason".into(), opt_s(o.reason.clone()));
                m.insert("note".into(), opt_s(o.note.clone()));
                m.insert(
                    "earlier_result".into(),
                    o.earlier.as_ref().map_or(Value::Null, |e| {
                        let mut x = Obj::new();
                        x.insert("section".into(), s(e.section.as_str()));
                        x.insert("source".into(), s(e.source.as_str()));
                        x.insert("file_date".into(), opt_s(e.file_date.clone()));
                        Value::Object(x)
                    }),
                );
                m.insert("started_utc".into(), opt_s(o.started_utc.clone()));
                m.insert("ended_utc".into(), opt_s(o.ended_utc.clone()));
                let procs = o
                    .processes
                    .iter()
                    .map(|(p, c)| {
                        let mut x = Obj::new();
                        for (k, v) in [
                            ("waiting", c.waiting),
                            ("running", c.running),
                            ("done", c.done),
                            ("cached", c.cached),
                            ("failed", c.failed),
                        ] {
                            x.insert(k.into(), Value::Int(v));
                        }
                        ((*p).to_owned(), Value::Object(x))
                    })
                    .collect();
                m.insert("processes".into(), Value::Object(procs));
                m.insert("sections".into(), list(&o.step.sections));
                m.insert("upstream".into(), list(&o.step.upstream));
                m.insert("hours".into(), opt_f(o.step.hours));
                m.insert(
                    "hours_from".into(),
                    opt_s(Some(o.step.hours_from.as_str()).filter(|t| !t.is_empty())),
                );
                Value::Object(m)
            })
            .collect();
        let mut top = Obj::new();
        top.insert("schema_version".into(), Value::Int(SCHEMA_VERSION));
        top.insert("sample".into(), s(self.sample.as_str()));
        top.insert("generated_utc".into(), opt_s(iso(self.generated)));
        top.insert("run".into(), Value::Object(run));
        top.insert("steps".into(), Value::Array(steps));
        Value::Object(top)
    }
}

/// The trace files of a run: those under `SAMPLE_DIR/../pipeline_info`
/// written since it started, oldest first, or only the newest when the
/// start is not known.
pub fn trace_files(sample_dir: &Path, started: Option<f64>) -> Vec<PathBuf> {
    let parent = crate::util::abspath(sample_dir);
    let parent = parent.parent().unwrap_or(Path::new("/"));
    let pattern = parent.join("pipeline_info").join("trace_*.txt");
    let mut found: Vec<(f64, PathBuf)> = pycompat::mtime::glob(pattern.as_os_str(), false)
        .into_iter()
        .filter_map(|p| pycompat::mtime::mtime(&p).ok().map(|t| (t, p)))
        .collect();
    found.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| pycompat::mtime::py_cmp(&a.1, &b.1))
    });
    match started {
        None => found.pop().map(|(_, p)| p).into_iter().collect(),
        Some(st) => found
            .into_iter()
            .map(|(_, p)| p)
            .filter(|p| pycompat::mtime::mtime(p).unwrap_or(0.0) >= st)
            .collect(),
    }
}

/// The newest time any of the run's own files changed: `logs/run_status.tsv`,
/// `nextflow/process_status.tsv` and `extra` (its events file). Not the
/// trace files: every sample's run writes them, so another run would keep
/// this one alive.
pub fn last_heard(sample_dir: &Path, extra: &[&Path]) -> Option<f64> {
    let own = [
        sample_dir.join("logs").join("run_status.tsv"),
        sample_dir.join("nextflow").join("process_status.tsv"),
    ];
    py_max(
        own.iter()
            .map(PathBuf::as_path)
            .chain(extra.iter().copied())
            .filter_map(|p| pycompat::mtime::mtime(p).ok()),
    )
}

/// What [`RunFiles::read`] reads beyond the sample folder.
pub struct Sources<'a> {
    /// Trace files; empty means the run's own under `../pipeline_info`.
    pub traces: &'a [PathBuf],
    pub events: Option<&'a Path>,
    /// Also call a quiet run with no end record `stopped` (the live site's
    /// view). The guards that refuse to rerun under a run going on leave it
    /// off: they must keep reading `running` for a run alive but quiet.
    pub stale: bool,
}

/// The run files of `sample_dir`, read for [`build`]. A process table older
/// than the run start is an earlier run's and is left out.
pub struct RunFiles {
    pub status: Option<RunStatus>,
    pub process_status: Option<IndexMap<String, ProcessCounts>>,
    pub events: Vec<Value>,
    pub traces: Vec<PathBuf>,
    pub sources: IndexMap<String, (String, f64)>,
    pub heard: Option<f64>,
}

impl RunFiles {
    pub fn read(sample_dir: &Path, src: &Sources<'_>) -> io::Result<RunFiles> {
        let status = files::read_run_status(&sample_dir.join("logs").join("run_status.tsv"))?;
        let started = status.as_ref().and_then(RunStatus::started);
        let pt_path = sample_dir.join("nextflow").join("process_status.tsv");
        let mut process_status = files::read_process_status(&pt_path)?;
        if process_status.is_some()
            && let Some(st) = started
            && pycompat::mtime::mtime(&pt_path)? < st
        {
            process_status = None;
        }
        Ok(RunFiles {
            status,
            process_status,
            events: files::read_events(src.events)?,
            traces: if src.traces.is_empty() {
                trace_files(sample_dir, started)
            } else {
                src.traces.to_vec()
            },
            sources: files::summary_sources(sample_dir),
            heard: if src.stale {
                last_heard(sample_dir, &src.events.into_iter().collect::<Vec<_>>())
            } else {
                None
            },
        })
    }

    /// The state of these files at `now`.
    pub fn state<'s>(&self, sample: &str, steps: &'s [Step], now: f64, tz: &TimeZone) -> State<'s> {
        let inp = Inputs {
            status: self.status.as_ref(),
            process_status: self.process_status.as_ref(),
            trace_tasks: tasks::trace_tasks(&self.traces, sample, tz),
            event_tasks: tasks::event_tasks(&self.events, sample),
            events: &self.events,
            sources: self.sources.clone(),
            heard: self.heard,
        };
        build(sample, steps, &inp, now, tz)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn step(id: &str, procs: &[&str], hours: Option<f64>) -> Step {
        Step {
            step: id.into(),
            label: format!("step {id}"),
            tool: String::new(),
            processes: procs.iter().map(|p| (*p).to_owned()).collect(),
            sections: vec![format!("sec{id}")],
            upstream: vec![],
            hours,
            hours_from: String::new(),
            hours_low: String::new(),
            reads: vec![],
        }
    }

    fn task(id: &str, proc: &str, status: &str, start: Option<f64>) -> Task {
        Task {
            id: id.into(),
            name: format!("W:{proc} (S1)"),
            process: proc.into(),
            status: status.into(),
            start,
            end: start.map(|s| s + 10.0),
        }
    }

    fn states(st: &State<'_>) -> Vec<&'static str> {
        st.steps.iter().map(|s| s.state.as_str()).collect()
    }

    #[test]
    fn rows_then_tables_then_tasks() {
        let steps = [
            step("01", &["A"], Some(1.0)),
            step("02", &["B"], None),
            step("03", &["C"], None),
            step("04", &["D"], None),
            step("05", &["E"], None),
        ];
        let mut status = RunStatus::default();
        status.meta.insert("started_epoch".into(), "100".into());
        for (k, v) in [
            ("01", "ok"),
            ("02", "planned"),
            ("03", "skipped (no BAM)"),
            ("04", "planned"),
        ] {
            status.steps.insert(k.into(), v.into());
        }
        let inp = Inputs {
            status: Some(&status),
            event_tasks: vec![
                task("1", "B", "RUNNING", Some(150.0)),
                task("2", "E", "COMPLETED", Some(160.0)),
            ],
            sources: [("sec02".to_owned(), ("02/x.tsv".to_owned(), 50.0))]
                .into_iter()
                .collect(),
            ..Inputs::default()
        };
        let st = build("S1", &steps, &inp, 1000.0, &TimeZone::UTC);
        assert_eq!(st.run, RunState::Running);
        assert_eq!(
            states(&st),
            ["done", "running", "skipped", "waiting", "done"]
        );
        assert_eq!(st.steps[2].reason.as_deref(), Some("no BAM"));
        assert_eq!(st.steps[1].earlier.as_ref().unwrap().source, "02/x.tsv");
        assert_eq!(
            st.steps[1].note.as_deref(),
            Some("an earlier result (1970-01-01T00:00:50Z) is being redone")
        );
        assert_eq!(
            st.steps[4].ended_utc.as_deref(),
            Some("1970-01-01T00:02:50Z")
        );
        assert_eq!(st.steps[1].ended_utc, None);
    }

    #[test]
    fn a_finished_run_leaves_planned_steps_not_finished() {
        let steps = [step("01", &["A"], None), step("02", &["HAPLOCHECK"], None)];
        let mut status = RunStatus::default();
        status
            .meta
            .insert("finished_utc".into(), "2026-01-03 05:12:00 UTC".into());
        status.steps.insert("01".into(), "planned".into());
        let mut ps = IndexMap::new();
        ps.insert("HAPLOCHECK".to_owned(), ProcessCounts::default());
        let inp = Inputs {
            status: Some(&status),
            process_status: Some(&ps),
            ..Inputs::default()
        };
        let st = build("S1", &steps, &inp, 0.0, &TimeZone::UTC);
        assert_eq!(st.run, RunState::Finished);
        // No MITO_VARIANTS row: step 20 did not run, so no HAPLOCHECK task
        // was due and an empty row is not "not finished".
        assert_eq!(states(&st), ["not finished", "not finished"]);
        assert_eq!(st.finished_utc.as_deref(), Some("2026-01-03T05:12:00Z"));
    }

    #[test]
    fn a_quiet_run_reads_stopped_only_when_asked() {
        let steps = [step("01", &["A"], Some(1.0))];
        let tasks = vec![task("1", "A", "RUNNING", Some(10.0))];
        let mut inp = Inputs {
            event_tasks: tasks,
            ..Inputs::default()
        };
        let late = 10.0 + 6.0 * 3600.0 + 11.0;
        assert_eq!(
            build("S1", &steps, &inp, late, &TimeZone::UTC).run,
            RunState::Running
        );
        inp.heard = Some(5.0);
        let st = build("S1", &steps, &inp, late, &TimeZone::UTC);
        assert_eq!(st.run, RunState::Stopped);
        assert_eq!(st.last_heard, Some(20.0));
        assert_eq!(
            build("S1", &steps, &inp, late - 2.0, &TimeZone::UTC).run,
            RunState::Running
        );
    }
}
