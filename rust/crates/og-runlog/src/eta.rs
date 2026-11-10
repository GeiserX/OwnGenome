//! How long each step usually takes on this box, and how many hours a run
//! has left.
//!
//! Each call first reads the trace files `OUTDIR/pipeline_info/trace_*.txt`
//! into this box's history, `OUTDIR/pipeline_info/step_history/<box>.tsv`.
//! The history keeps one row per task that ended COMPLETED with exit 0,
//! with no sample id and no path. A step's range is the shortest and the
//! longest time it took in the last 5 runs on this box that hold it; with
//! no history it takes the catalogue's measured hours (16 threads), the
//! high end scaled up on a box with fewer cores. The hours left: the low
//! end is the longest chain of waiting and running steps at their low ends;
//! the high end is the longer of that chain at the high ends and the
//! steps' CPU work over the cores; a running step past its high end makes
//! the high end unknown.

use crate::state::{RunFiles, Sources, State, trace_files};
use crate::steps::{Step, process_steps};
use crate::util::{Obj, basename, opt_f, opt_s, py_max, py_min, s};
use indexmap::{IndexMap, IndexSet};
use jiff::tz::TimeZone;
use og_nf::time::{ISO_Z, iso, parse_utc};
use og_nf::units::{blank, duration, memory, time};
use pycompat::json::Value;
use std::fs::{self, OpenOptions};
use std::io::{self, Write};
use std::path::{Path, PathBuf};

pub const HISTORY_COLUMNS: [&str; 13] = [
    "trace_file",
    "hash",
    "attempt",
    "step",
    "process",
    "run_level",
    "cpus",
    "memory_bytes",
    "peak_rss_bytes",
    "realtime_s",
    "step_span_s",
    "finished_utc",
    "box_cpus",
];
/// The runs per step the range looks back on.
pub const LAST_RUNS: usize = 5;
/// The threads the catalogue's hours were measured on.
pub const DOCS_THREADS: f64 = 16.0;

/// One trace row, as the history and the headroom read it.
#[derive(Debug, Clone, PartialEq)]
pub struct TraceTask {
    pub trace_file: String,
    /// The task hash without its `/`.
    pub hash: String,
    /// As the trace spells it; `1` when empty.
    pub attempt: String,
    pub process: Option<String>,
    pub tag: Option<String>,
    pub status: String,
    pub exit: String,
    pub cpus: Option<i128>,
    pub memory_bytes: Option<i128>,
    pub peak_rss_bytes: Option<i128>,
    pub realtime_s: Option<f64>,
    pub start: Option<f64>,
    pub complete: Option<f64>,
}

/// Every row of a trace file, old (the default columns) or new (with the
/// columns nextflow.config adds), human or raw units. An old trace has no
/// start or complete: the task ended at submit + duration and ran for its
/// realtime before that. A trace that cannot be read gives nothing.
pub fn read_trace(path: &Path, tz: &TimeZone) -> Vec<TraceTask> {
    let mut out = Vec::new();
    let file = basename(path);
    let _ = og_nf::for_each_row(path, |r| {
        if r.is_blank() {
            return;
        }
        let (mut proc, mut tag) = r
            .get("name")
            .and_then(og_nf::task_name)
            .map_or((None, None), |(p, t)| {
                (Some(p.to_owned()), t.map(str::to_owned))
            });
        if !blank(r.get("process")) {
            let p = r.get("process").unwrap_or("");
            proc = Some(p.rsplit(':').next().unwrap_or(p).to_owned());
        }
        if let Some(t) = r.get("tag") {
            tag = (!blank(Some(t))).then(|| t.to_owned());
        }
        let submit = time(r.get("submit"), tz);
        let dur = duration(r.get("duration"));
        let real = duration(r.get("realtime"));
        let complete = time(r.get("complete"), tz).or_else(|| submit.zip(dur).map(|(s, d)| s + d));
        let start = time(r.get("start"), tz).or_else(|| {
            complete.and_then(|c| {
                if let Some(x) = real {
                    Some(c - x)
                } else {
                    submit
                }
            })
        });
        let cpus = if blank(r.get("cpus")) {
            None
        } else {
            r.get("cpus").and_then(|v| pycompat::num::parse_int(v).ok())
        };
        out.push(TraceTask {
            trace_file: file.clone(),
            hash: r.get("hash").unwrap_or("").replace('/', ""),
            attempt: if blank(r.get("attempt")) {
                "1"
            } else {
                r.get("attempt").unwrap_or("1")
            }
            .to_owned(),
            process: proc,
            tag,
            status: pycompat::text::upper(r.get("status").unwrap_or("")),
            exit: pycompat::text::strip(r.get("exit").unwrap_or(""), None).to_owned(),
            cpus,
            memory_bytes: memory(r.get("memory")),
            peak_rss_bytes: memory(r.get("peak_rss")),
            realtime_s: real,
            start,
            complete,
        });
    });
    out
}

fn shared(t: &TraceTask) -> bool {
    t.process
        .as_deref()
        .is_some_and(|p| crate::tasks::SHARED_PROCESSES.contains(&p))
}

fn run_level(t: &TraceTask) -> bool {
    shared(t) || t.tag.is_none()
}

/// True when every task of a step's group ended well: per process and tag,
/// the attempt that ended last is COMPLETED with exit 0.
fn finished(ts: &[&TraceTask]) -> bool {
    type Last<'t> = (f64, usize, &'t TraceTask);
    let mut last: IndexMap<(Option<&str>, Option<&str>), Last<'_>> = IndexMap::new();
    for (i, t) in ts.iter().enumerate() {
        let k = (t.process.as_deref(), t.tag.as_deref());
        let c = t.complete.unwrap_or(f64::NEG_INFINITY);
        if last.get(&k).is_none_or(|(oc, oi, _)| (c, i) > (*oc, *oi)) {
            last.insert(k, (c, i, t));
        }
    }
    last.values()
        .all(|(_, _, t)| t.status == "COMPLETED" && t.exit == "0")
}

/// `fmt()` of the history: a float to 3 decimals, trailing zeros and point
/// dropped (correctly rounded, as Python's `.3f`).
pub fn fmt_float(v: f64) -> String {
    let s = pycompat::fmt::fixed(v, 3);
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_owned()
    } else {
        s
    }
}

fn fmt_opt_f(v: Option<f64>) -> String {
    v.map_or_else(String::new, fmt_float)
}

fn fmt_opt_i(v: Option<i128>) -> String {
    v.map_or_else(String::new, |v| v.to_string())
}

/// One history row, its cells already written as the file holds them.
pub type HistoryRow = [String; 13];

/// The history rows of one trace. A step's span is its first start to its
/// last end among one sample's tasks (or the run-level tasks); a span with
/// a cached task in it, or a task whose last attempt failed, is not the
/// step's full time and is left empty.
pub fn history_rows(
    tasks: &[TraceTask],
    proc_step: &IndexMap<&str, &str>,
    box_cpus: usize,
) -> Vec<HistoryRow> {
    let mut groups: IndexMap<(&str, String), Vec<&TraceTask>> = IndexMap::new();
    for t in tasks {
        let Some(step) = t.process.as_deref().and_then(|p| proc_step.get(p)) else {
            continue;
        };
        let g = if run_level(t) {
            String::new()
        } else {
            let tag = t.tag.as_deref().unwrap_or("");
            tag.split_once(':').map_or(tag, |(a, _)| a).to_owned()
        };
        groups.entry((step, g)).or_default().push(t);
    }
    let mut span_of: IndexMap<(&str, String), Option<f64>> = IndexMap::new();
    for (key, ts) in &groups {
        let span = if ts.iter().any(|t| t.status == "CACHED") || !finished(ts) {
            None
        } else {
            let lo = py_min(ts.iter().filter_map(|t| t.start));
            let hi = py_max(ts.iter().filter_map(|t| t.complete));
            lo.zip(hi).map(|(a, b)| b - a)
        };
        span_of.insert(key.clone(), span);
    }
    let mut rows = Vec::new();
    let mut seen: IndexSet<(&str, &str)> = IndexSet::new();
    for (key, ts) in &groups {
        for t in ts {
            if t.status != "COMPLETED" || t.exit != "0" || t.hash.is_empty() {
                continue;
            }
            // The trace can list a finished task twice.
            if !seen.insert((t.hash.as_str(), t.attempt.as_str())) {
                continue;
            }
            rows.push([
                t.trace_file.clone(),
                t.hash.clone(),
                t.attempt.clone(),
                key.0.to_owned(),
                t.process.clone().unwrap_or_default(),
                if key.1.is_empty() { "yes" } else { "no" }.to_owned(),
                fmt_opt_i(t.cpus),
                fmt_opt_i(t.memory_bytes),
                fmt_opt_i(t.peak_rss_bytes),
                fmt_opt_f(t.realtime_s),
                fmt_opt_f(span_of[key]),
                t.complete.and_then(iso).unwrap_or_default(),
                box_cpus.to_string(),
            ]);
        }
    }
    rows
}

/// The history's name for this box: `name`, else `PGP_BOX`, else the short
/// host name, with every character other than letters, digits, `_` and
/// `-` turned into `_`.
pub fn box_name(name: Option<&str>) -> String {
    let env = std::env::var("PGP_BOX").ok();
    let raw = name
        .filter(|n| !n.is_empty())
        .map(str::to_owned)
        .or(env.filter(|e| !e.is_empty()))
        .unwrap_or_else(crate::host::short_hostname);
    let clean: String = raw
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' || c == '-' {
                c
            } else {
                '_'
            }
        })
        .collect();
    if clean.is_empty() {
        "box".to_owned()
    } else {
        clean
    }
}

pub fn history_path(outdir: &Path, box_: &str) -> PathBuf {
    outdir
        .join("pipeline_info")
        .join("step_history")
        .join(format!("{box_}.tsv"))
}

/// The history file's rows, every column of [`HISTORY_COLUMNS`] filled
/// (empty when the file lacks it).
pub fn read_history(path: &Path) -> Vec<HistoryRow> {
    let Ok(f) = fs::File::open(path) else {
        return Vec::new();
    };
    let mut lines = og_io::Lines::new(io::BufReader::new(f), og_io::Utf8::Replace);
    let Some(Ok(head)) = lines.next() else {
        return Vec::new();
    };
    let header: Vec<&str> = head.split('\t').collect();
    let mut out = Vec::new();
    for line in lines {
        let Ok(line) = line else { break };
        if line.chars().all(pycompat::text::is_space) {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        let get = |k: &str| {
            (0..header.len().min(cells.len()))
                .rev()
                .find(|&i| header[i] == k)
                .map_or("", |i| cells[i])
        };
        out.push(HISTORY_COLUMNS.map(|k| get(k).to_owned()));
    }
    out
}

/// Reads `OUTDIR/pipeline_info/trace_*.txt` into this box's history and
/// returns the history. A trace already read is read again and its rows
/// replaced, so a trace still being written ends with its full rows; a
/// task the history holds from another trace is not added twice; rows of a
/// deleted trace stay. The file is rewritten only when its text changes,
/// under an exclusive `flock(2)` on `<history>.lock`, the lock `flock(1)`
/// takes too.
pub fn ingest(
    outdir: &Path,
    box_: &str,
    steps: &[Step],
    box_cpus: usize,
    tz: &TimeZone,
) -> io::Result<Vec<HistoryRow>> {
    let path = history_path(outdir, box_);
    fs::create_dir_all(path.parent().expect("history has a folder"))?;
    let proc_step = process_steps(steps);
    let mut lock_name = path.clone().into_os_string();
    lock_name.push(".lock");
    let lock = OpenOptions::new()
        .append(true)
        .create(true)
        .open(PathBuf::from(lock_name))?;
    let _ = &lock;
    let mut by_file: IndexMap<String, Vec<HistoryRow>> = IndexMap::new();
    for r in read_history(&path) {
        by_file.entry(r[0].clone()).or_default().push(r);
    }
    let pattern = outdir.join("pipeline_info").join("trace_*.txt");
    let mut traces = pycompat::mtime::glob(pattern.as_os_str(), false);
    pycompat::mtime::sort(&mut traces);
    for tp in traces {
        let name = basename(&tp);
        let others: IndexSet<(String, String)> = by_file
            .iter()
            .filter(|(f, _)| **f != name)
            .flat_map(|(_, rs)| rs.iter().map(|r| (r[1].clone(), r[2].clone())))
            .collect();
        let new: Vec<HistoryRow> = history_rows(&read_trace(&tp, tz), &proc_step, box_cpus)
            .into_iter()
            .filter(|r| !others.contains(&(r[1].clone(), r[2].clone())))
            .collect();
        if !new.is_empty() || by_file.contains_key(&name) {
            by_file.insert(name, new);
        }
    }
    let mut files: Vec<&String> = by_file.keys().collect();
    files.sort();
    let rows: Vec<HistoryRow> = files
        .into_iter()
        .flat_map(|f| by_file[f].iter().cloned())
        .collect();
    let mut text = HISTORY_COLUMNS.join("\t");
    text.push('\n');
    for r in &rows {
        text.push_str(&r.join("\t"));
        text.push('\n');
    }
    if crate::util::read_text(&path).ok().as_deref() != Some(text.as_str()) {
        let mut f = og_io::AtomicFile::create(&path)?;
        f.write_all(text.as_bytes())?;
        f.commit()?;
    }
    Ok(rows)
}

/// `float(v)` of a history cell; `None` for an empty or bad one.
fn num(v: &str) -> Option<f64> {
    if v.is_empty() {
        None
    } else {
        pycompat::num::parse_float(v).ok()
    }
}

/// The catalogue's hours were measured on 16 threads: fewer cores scale
/// them up, more never scale them below that.
pub fn docs_scale(cores: usize) -> f64 {
    let x = DOCS_THREADS / cores.max(1) as f64;
    if x > 1.0 { x } else { 1.0 }
}

/// Where a step's range comes from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Source {
    /// The last runs on this box.
    Box,
    /// The catalogue's measured hours.
    Docs,
    /// No number.
    None,
}

/// One step's usual time.
#[derive(Debug, Clone, PartialEq)]
pub struct Range {
    pub low_s: Option<f64>,
    pub high_s: Option<f64>,
    pub n: usize,
    pub source: Source,
    pub label: String,
    /// The most CPUs its tasks asked for (box), 1 (docs), none (no number).
    pub cpus: Option<i128>,
}

/// `{step: range}` in catalogue order, leaving out the history of the
/// traces in `exclude` (runs still going).
pub fn step_ranges(
    history: &[HistoryRow],
    steps: &[Step],
    cores: usize,
    exclude: &IndexSet<String>,
) -> IndexMap<String, Range> {
    struct Run {
        spans: Vec<f64>,
        end: String,
        cpus: i128,
    }
    let mut runs: IndexMap<&str, IndexMap<&str, Run>> = IndexMap::new();
    for r in history {
        if exclude.contains(&r[0]) {
            continue;
        }
        let run = runs
            .entry(r[3].as_str())
            .or_default()
            .entry(r[0].as_str())
            .or_insert(Run {
                spans: Vec::new(),
                end: String::new(),
                cpus: 0,
            });
        if let Some(span) = num(&r[10]) {
            run.spans.push(span);
        }
        if r[11] > run.end {
            run.end.clone_from(&r[11]);
        }
        let cpus = num(&r[6])
            .filter(|c| *c != 0.0 && c.is_finite())
            .map_or(0, |c| c.trunc() as i128);
        run.cpus = run.cpus.max(cpus);
    }
    let scale = docs_scale(cores);
    let mut out = IndexMap::new();
    for s in steps {
        let mut mine: Vec<&Run> = runs
            .get(s.step.as_str())
            .map(|m| m.values().filter(|v| !v.spans.is_empty()).collect())
            .unwrap_or_default();
        mine.sort_by(|a, b| a.end.cmp(&b.end));
        let mine = &mine[mine.len().saturating_sub(LAST_RUNS)..];
        let range = if !mine.is_empty() {
            let spans: Vec<f64> = mine
                .iter()
                .filter_map(|v| py_max(v.spans.iter().copied()))
                .collect();
            let n = mine.len();
            Range {
                low_s: py_min(spans.iter().copied()),
                high_s: py_max(spans.iter().copied()),
                n,
                source: Source::Box,
                label: format!("on this box, {n} run{}", if n > 1 { "s" } else { "" }),
                cpus: Some(
                    mine.iter()
                        .map(|v| v.cpus)
                        .max()
                        .filter(|c| *c != 0)
                        .unwrap_or(1),
                ),
            }
        } else if let (Some(h), "measured") = (s.hours, s.hours_from.as_str()) {
            let h = h * 3600.0;
            let lo = num(&s.hours_low);
            Range {
                low_s: Some(lo.map_or(h, |l| l * 3600.0)),
                high_s: Some(h * scale),
                n: 0,
                source: Source::Docs,
                label: format!(
                    "from the docs, measured on 16 threads{}",
                    if scale > 1.0 {
                        format!(", scaled to {cores} cores")
                    } else {
                        String::new()
                    }
                ),
                cpus: Some(1),
            }
        } else {
            Range {
                low_s: None,
                high_s: None,
                n: 0,
                source: Source::None,
                label: "not measured".to_owned(),
                cpus: None,
            }
        };
        out.insert(s.step.clone(), range);
    }
    out
}

/// The hours left of a run.
#[derive(Debug, Clone, PartialEq)]
pub struct Remaining {
    pub low_s: f64,
    /// `None` when a running step is past its high end.
    pub high_s: Option<f64>,
    pub unknown_steps: Vec<String>,
    pub late_steps: Vec<String>,
    pub note: Option<String>,
    pub basis: String,
}

/// The hours left of the run in `state` at `now`.
pub fn remaining(
    state: &State<'_>,
    ranges: &IndexMap<String, Range>,
    now: f64,
    cores: usize,
) -> Remaining {
    use crate::state::{RunState, StepState};
    let mut steps: IndexMap<&str, &crate::state::StepOut<'_>> = IndexMap::new();
    for s in &state.steps {
        steps.insert(s.step.step.as_str(), s);
    }
    let left: Vec<&str> = if state.run == RunState::Finished {
        Vec::new()
    } else {
        steps
            .iter()
            .filter(|(_, s)| matches!(s.state, StepState::Waiting | StepState::Running))
            .map(|(k, _)| *k)
            .collect()
    };
    let mut dur: IndexMap<&str, (f64, f64)> = IndexMap::new();
    let (mut unknown, mut late) = (Vec::new(), Vec::new());
    for &k in &left {
        let r = ranges.get(k);
        let (Some(lo), Some(hi)) = (r.and_then(|r| r.low_s), r.and_then(|r| r.high_s)) else {
            unknown.push(k.to_owned());
            dur.insert(k, (0.0, 0.0));
            continue;
        };
        let (mut lo, mut hi) = (lo, hi);
        let started = steps[k]
            .started_utc
            .as_deref()
            .and_then(|v| parse_utc(v, &[ISO_Z]));
        if steps[k].state == StepState::Running
            && let Some(st) = started
        {
            let gone = if now - st > 0.0 { now - st } else { 0.0 };
            let l = lo - gone;
            lo = if l > 0.0 { l } else { 0.0 };
            hi -= gone;
            if hi < 0.0 {
                hi = 0.0;
                late.push(k.to_owned());
            }
        }
        dur.insert(k, (lo, hi));
    }
    fn chain<'k>(
        k: &'k str,
        i: usize,
        steps: &IndexMap<&'k str, &crate::state::StepOut<'k>>,
        dur: &IndexMap<&str, (f64, f64)>,
        memo: &mut IndexMap<&'k str, f64>,
        path: &mut Vec<&'k str>,
    ) -> f64 {
        if let Some(v) = memo.get(k) {
            return *v;
        }
        if path.contains(&k) {
            return 0.0;
        }
        path.push(k);
        let mut up = Vec::new();
        for u in &steps[k].step.upstream {
            if let Some((u, _)) = dur.get_key_value(u.as_str())
                && let Some((&u, _)) = steps.get_key_value(*u)
            {
                up.push(chain(u, i, steps, dur, memo, path));
            }
        }
        path.pop();
        let d = dur[k];
        let v = if i == 0 { d.0 } else { d.1 } + py_max(up).unwrap_or(0.0);
        memo.insert(k, v);
        v
    }
    let longest = |i: usize| {
        py_max(
            left.iter()
                .map(|&k| chain(k, i, &steps, &dur, &mut IndexMap::new(), &mut Vec::new())),
        )
        .unwrap_or(0.0)
    };
    let low = longest(0);
    let high_path = longest(1);
    let work = left
        .iter()
        .filter(|k| !unknown.iter().any(|u| u == *k))
        .fold(0.0, |acc, k| {
            let cpus = ranges
                .get(*k)
                .and_then(|r| r.cpus)
                .filter(|c| *c != 0)
                .unwrap_or(1);
            acc + dur[k].1 * cpus as f64
        });
    let cores_f = cores.max(1) as f64;
    let bound = work / cores_f;
    let high = if late.is_empty() {
        py_max([high_path, bound, low])
    } else {
        None
    };
    let note = late
        .iter()
        .map(|k| format!("step {k} is taking longer than usual"))
        .collect::<Vec<_>>()
        .join("; ");
    Remaining {
        low_s: low,
        high_s: high,
        note: (!note.is_empty()).then_some(note),
        basis: format!(
            "low: the longest chain of waiting and running steps at their low ends; high: the longer of that chain at the high ends and their CPU work ({} CPU-hours) over {cores} cores{}",
            pycompat::fmt::fixed(work / 3600.0, 1),
            if unknown.is_empty() {
                ""
            } else {
                "; steps with no number add nothing"
            }
        ),
        unknown_steps: unknown,
        late_steps: late,
    }
}

/// The trace files of every run that run-all.sh started under `outdir` and
/// has not ended: its `run_status.tsv` has `meta started_epoch` and no `meta
/// finished_utc`. Their partial spans are not history yet.
pub fn live_traces(outdir: &Path) -> IndexSet<String> {
    let mut out = IndexSet::new();
    let pattern = outdir.join("*").join("logs").join("run_status.tsv");
    for p in pycompat::mtime::glob(pattern.as_os_str(), false) {
        let sample_dir = p.parent().and_then(Path::parent).unwrap_or(Path::new("/"));
        let meta = match crate::files::read_run_status(&p) {
            Ok(st) => st.unwrap_or_default().meta,
            Err(_) => {
                out.extend(trace_files(sample_dir, None).iter().map(|t| basename(t)));
                continue;
            }
        };
        if meta.contains_key("finished_utc") {
            continue;
        }
        let Some(started) = meta
            .get("started_epoch")
            .and_then(|v| pycompat::num::parse_float(v).ok())
        else {
            continue;
        };
        out.extend(
            trace_files(sample_dir, Some(started))
                .iter()
                .map(|t| basename(t)),
        );
    }
    out
}

/// The estimate: this box's ranges and, with a sample, its hours left.
pub struct Estimate<'a> {
    pub box_: String,
    pub cores: usize,
    pub generated: f64,
    pub ranges: IndexMap<String, Range>,
    pub remaining: Option<(State<'a>, Remaining)>,
}

/// What [`estimate`] needs beyond the outdir.
pub struct EtaArgs<'a> {
    pub sample: Option<&'a str>,
    pub cores: Option<usize>,
    pub box_: Option<&'a str>,
    pub now: f64,
    pub tz: &'a TimeZone,
}

pub fn estimate<'s>(outdir: &Path, steps: &'s [Step], a: &EtaArgs<'_>) -> io::Result<Estimate<'s>> {
    let cpus = crate::host::cpu_count();
    let cores = a.cores.filter(|c| *c != 0).unwrap_or(cpus);
    let box_ = box_name(a.box_);
    let history = ingest(outdir, &box_, steps, cpus, a.tz)?;
    let mut exclude = live_traces(outdir);
    let state = match a.sample {
        Some(sample) => {
            let dir = outdir.join(sample);
            let files = RunFiles::read(
                &dir,
                &Sources {
                    traces: &[],
                    events: None,
                    stale: false,
                },
            )?;
            let st = files.state(sample, steps, a.now, a.tz);
            if st.run == crate::state::RunState::Running {
                let started = st
                    .started
                    .and_then(iso)
                    .and_then(|v| parse_utc(&v, &[ISO_Z]));
                exclude.extend(trace_files(&dir, started).iter().map(|t| basename(t)));
            }
            Some(st)
        }
        None => None,
    };
    let ranges = step_ranges(&history, steps, cores, &exclude);
    let remaining = state.map(|st| {
        let r = remaining(&st, &ranges, a.now, cores);
        (st, r)
    });
    Ok(Estimate {
        box_,
        cores,
        generated: a.now,
        ranges,
        remaining,
    })
}

impl Estimate<'_> {
    pub fn to_json(&self) -> Value {
        let steps = self
            .ranges
            .iter()
            .map(|(k, r)| {
                let mut m = Obj::new();
                m.insert("low_s".into(), opt_f(r.low_s));
                m.insert("high_s".into(), opt_f(r.high_s));
                m.insert("n".into(), Value::Int(r.n as i128));
                m.insert(
                    "source".into(),
                    s(match r.source {
                        Source::Box => "box",
                        Source::Docs => "docs",
                        Source::None => "none",
                    }),
                );
                m.insert("label".into(), s(r.label.as_str()));
                m.insert("cpus".into(), r.cpus.map_or(Value::Null, Value::Int));
                (k.clone(), Value::Object(m))
            })
            .collect();
        let mut top = Obj::new();
        top.insert("box".into(), s(self.box_.as_str()));
        top.insert("cores".into(), Value::Int(self.cores as i128));
        top.insert("generated_utc".into(), opt_s(iso(self.generated)));
        top.insert("steps".into(), Value::Object(steps));
        top.insert(
            "remaining".into(),
            self.remaining.as_ref().map_or(Value::Null, |(_, r)| {
                let strs = |v: &[String]| Value::Array(v.iter().map(|x| s(x.as_str())).collect());
                let mut m = Obj::new();
                m.insert("low_s".into(), Value::Float(r.low_s));
                m.insert("high_s".into(), opt_f(r.high_s));
                m.insert("unknown_steps".into(), strs(&r.unknown_steps));
                m.insert("late_steps".into(), strs(&r.late_steps));
                m.insert("note".into(), opt_s(r.note.clone()));
                m.insert("basis".into(), s(r.basis.as_str()));
                Value::Object(m)
            }),
        );
        Value::Object(top)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn history_numbers() {
        assert_eq!(fmt_float(1210.0), "1210");
        assert_eq!(fmt_float(0.0), "0");
        assert_eq!(fmt_float(2.5), "2.5");
        assert_eq!(fmt_float(0.0005), "0.001");
        assert_eq!(fmt_float(1.0005), "1");
        assert_eq!(fmt_float(100.25), "100.25");
        assert_eq!(fmt_float(-0.0001), "-0");
    }

    #[test]
    fn scale() {
        assert_eq!(docs_scale(8), 2.0);
        assert_eq!(docs_scale(32), 1.0);
        assert_eq!(docs_scale(0), 16.0);
    }

    #[test]
    fn box_names() {
        assert_eq!(box_name(Some("my box.1")), "my_box_1");
        assert_eq!(box_name(Some("ñ")), "_");
    }
}
