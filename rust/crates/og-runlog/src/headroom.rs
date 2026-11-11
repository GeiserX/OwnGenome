//! How close each task came to the memory it asked for, and how full the
//! work and results disks are.
//!
//! Per task, the ratio is its peak_rss over its request. At or above the
//! warn percent it is a hint, not a fault: peak_rss counts shared and
//! file-mapped pages. An attempt that ended with an exit conf/base.config
//! retries (killed for memory or time) and that a later attempt replaced
//! counts as a retry. The disks are the work folder and the results folder,
//! given or from the weblog metadata; a folder not created yet is measured
//! on its nearest existing parent, and two folders on one file system are
//! reported once. Labels are `work` and `results`, never a path.

use crate::eta::TraceTask;
use crate::host::Host;
use crate::steps::{Step, process_steps};
use crate::util::abspath;
use indexmap::IndexMap;
use pycompat::json::Value;
use pycompat::value::{py_repr, py_str, truthy};
use std::path::{Path, PathBuf};

pub const DEFAULT_WARN_PERCENT: f64 = 90.0;
/// The exits conf/base.config's errorStrategy retries: usually killed for
/// memory or time.
pub const MEMORY_EXITS: [&str; 6] = ["143", "137", "104", "134", "139", "247"];

fn finished(status: &str) -> bool {
    matches!(status, "COMPLETED" | "FAILED" | "ABORTED" | "CACHED")
}

/// A warn percent that is not a number from 1 to 100; its text is the
/// message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BadPercent(pub String);

impl std::fmt::Display for BadPercent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

/// The warn percent: `value`, else the environment variable `env`, else 90.
pub fn warn_percent(value: Option<&str>, env: Option<&str>) -> Result<f64, BadPercent> {
    let from_env = env.and_then(|e| std::env::var(e).ok());
    let raw = value
        .filter(|v| !v.is_empty())
        .map(str::to_owned)
        .or(from_env);
    let Some(raw) = raw.filter(|r| !r.is_empty()) else {
        return Ok(DEFAULT_WARN_PERCENT);
    };
    let bad = || {
        let shown = py_repr(&Value::Str(raw.clone())).unwrap_or_else(|_| format!("'{raw}'"));
        BadPercent(format!(
            "warn percent must be a number from 1 to 100, got {shown}"
        ))
    };
    let p = pycompat::num::parse_float(&raw).map_err(|_| bad())?;
    if (1.0..=100.0).contains(&p) {
        Ok(p)
    } else {
        Err(bad())
    }
}

/// A task row: a trace row or a weblog message's trace object.
pub enum TaskRow<'a> {
    Trace(&'a TraceTask),
    Weblog(&'a Value),
}

/// Identifies one attempt's rows across messages: the weblog's task id, or
/// the hash and attempt, compared as Python compares dict keys.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum Key {
    One(String),
    Pair(String, String),
}

/// A value as a dict key: numbers by value (`1`, `1.0` and `True` are one
/// key), strings apart from numbers.
fn key_part(v: Option<&Value>) -> String {
    match v {
        None | Some(Value::Null) => "N".into(),
        Some(Value::Bool(b)) => format!("n{}", u8::from(*b)),
        Some(Value::Int(i)) => format!("n{i}"),
        Some(Value::Float(x)) if x.fract() == 0.0 && x.abs() < 1e30 => format!("n{}", *x as i128),
        Some(Value::Str(s)) => format!("s{s}"),
        Some(other) => format!("x{}", py_repr(other).unwrap_or_default()),
    }
}

struct Row {
    key: Key,
    order: usize,
    process: Option<String>,
    tag: Option<String>,
    attempt: i128,
    status: String,
    exit: String,
    memory: Option<i128>,
    peak: Option<i128>,
}

/// `str(v)` of a weblog field, `None` for a missing or null one.
fn text(v: Option<&Value>) -> Option<String> {
    match v {
        None | Some(Value::Null) => None,
        Some(v) => py_str(v).ok(),
    }
}

/// `int(v or 1)`; 1 for a value that is not a whole number.
fn attempt_of(v: Option<&Value>) -> i128 {
    match v.filter(|v| truthy(v)) {
        None => 1,
        Some(Value::Int(i)) => *i,
        Some(Value::Bool(_)) => 1,
        Some(Value::Float(x)) if x.is_finite() => x.trunc() as i128,
        Some(Value::Str(s)) => pycompat::num::parse_int(s).unwrap_or(1),
        Some(_) => 1,
    }
}

fn row(r: &TaskRow<'_>, order: usize) -> Row {
    match r {
        TaskRow::Trace(t) => Row {
            key: Key::Pair(format!("s{}", t.hash), format!("s{}", t.attempt)),
            order,
            process: t.process.clone(),
            tag: t.tag.clone(),
            attempt: pycompat::num::parse_int(&t.attempt).unwrap_or(1),
            status: t.status.clone(),
            exit: t.exit.clone(),
            memory: t.memory_bytes,
            peak: t.peak_rss_bytes,
        },
        TaskRow::Weblog(v) => {
            let f = |k: &str| match v {
                Value::Object(o) => o.get(k),
                _ => None,
            };
            let (mut proc, mut tag) = match f("name") {
                Some(Value::Str(n)) => og_nf::task_name(n).map_or((None, None), |(p, t)| {
                    (Some(p.to_owned()), t.map(str::to_owned))
                }),
                _ => (None, None),
            };
            let blank = |k: &str| og_nf::units::blank(text(f(k)).as_deref());
            if !blank("process") {
                let p = text(f("process")).unwrap_or_default();
                proc = Some(p.rsplit(':').next().unwrap_or(&p).to_owned());
            }
            if !blank("tag") {
                tag = text(f("tag"));
            }
            let key = match f("task_id") {
                Some(id) if *id != Value::Null => Key::One(key_part(Some(id))),
                _ => Key::Pair(key_part(f("hash")), key_part(f("attempt"))),
            };
            Row {
                key,
                order,
                process: proc,
                tag,
                attempt: attempt_of(f("attempt")),
                status: pycompat::text::upper(
                    &f("status")
                        .filter(|v| truthy(v))
                        .and_then(|v| py_str(v).ok())
                        .unwrap_or_default(),
                ),
                exit: if blank("exit") {
                    String::new()
                } else {
                    pycompat::text::strip(&text(f("exit")).unwrap_or_default(), None).to_owned()
                },
                memory: og_nf::units::memory(text(f("memory")).as_deref()),
                peak: og_nf::units::memory(text(f("peak_rss")).as_deref()),
            }
        }
    }
}

/// One task's memory.
#[derive(Debug, Clone, PartialEq)]
pub struct TaskMemory {
    pub step: Option<String>,
    pub process: String,
    pub attempt: i128,
    /// `running` while it runs, else its status in lower case.
    pub state: String,
    pub request_bytes: Option<i128>,
    pub peak_rss_bytes: Option<i128>,
    pub ratio: Option<f64>,
    pub warn: bool,
    /// Earlier attempts killed for memory or time.
    pub oom_retry: usize,
}

/// One entry per task (process and tag), from its last attempt: the latest
/// row of an attempt wins, unless it says running after a finished one.
pub fn task_memory(rows: &[TaskRow<'_>], warn: f64, steps: &[Step]) -> Vec<TaskMemory> {
    let proc_step = process_steps(steps);
    let mut latest: IndexMap<Key, Row> = IndexMap::new();
    for (i, r) in rows.iter().enumerate() {
        let t = row(r, i);
        if t.process.as_deref().is_none_or(str::is_empty) {
            continue;
        }
        let replaces = latest.get(&t.key).is_none_or(|old| {
            (t.attempt, finished(&t.status)) >= (old.attempt, finished(&old.status))
        });
        if replaces {
            latest.insert(t.key.clone(), t);
        }
    }
    let mut groups: IndexMap<(String, Option<String>), Vec<Row>> = IndexMap::new();
    for t in latest.into_values() {
        groups
            .entry((t.process.clone().unwrap_or_default(), t.tag.clone()))
            .or_default()
            .push(t);
    }
    let mut groups: Vec<_> = groups.into_iter().collect();
    groups.sort_by_key(|(_, ts)| ts.iter().map(|t| t.order).min());
    let mut out = Vec::new();
    for ((proc, _), mut ts) in groups {
        ts.sort_by_key(|t| (t.attempt, finished(&t.status), t.order));
        let last = ts.last().expect("a group has a row");
        let mut retried: Vec<i128> = ts
            .iter()
            .filter(|t| t.attempt < last.attempt && MEMORY_EXITS.contains(&t.exit.as_str()))
            .map(|t| t.attempt)
            .collect();
        retried.sort_unstable();
        retried.dedup();
        let running = !finished(&last.status);
        let ratio = match (last.memory, last.peak) {
            (Some(m), Some(p)) if !running && m != 0 => Some(p as f64 / m as f64),
            _ => None,
        };
        out.push(TaskMemory {
            step: proc_step.get(proc.as_str()).map(|s| (*s).to_owned()),
            attempt: last.attempt,
            state: if running {
                "running".to_owned()
            } else {
                pycompat::text::lower(&last.status)
            },
            request_bytes: last.memory,
            peak_rss_bytes: if running { None } else { last.peak },
            warn: ratio.is_some_and(|r| r * 100.0 >= warn),
            ratio,
            oom_retry: retried.len(),
            process: proc,
        });
    }
    out
}

/// The machine's memory against the pipeline's `max_memory`.
#[derive(Debug, Clone, PartialEq)]
pub struct HostMemory {
    pub memory_bytes: Option<u128>,
    pub max_memory_bytes: Option<i128>,
    /// `max_memory` is above the machine's memory.
    pub warn: bool,
}

/// One file system's use.
#[derive(Debug, Clone, PartialEq)]
pub struct Disk {
    /// `work`, `results`, or both joined with ` and `.
    pub fs: String,
    pub used_ratio: Option<f64>,
    pub free_bytes: u128,
    pub size_bytes: u128,
    pub inodes_used_ratio: Option<f64>,
    /// A folder is not created yet: measured on its nearest parent.
    pub missing: bool,
    pub warn: bool,
}

/// `path`, or its nearest parent that exists.
fn nearest_existing(host: &dyn Host, path: &Path) -> PathBuf {
    let mut p = abspath(path);
    while !host.exists(&p) {
        match p.parent() {
            Some(up) if up != p => p = up.to_path_buf(),
            _ => break,
        }
    }
    p
}

/// The disks of `paths` (label, folder); a pair with no folder is left out,
/// and so is a folder whose file system cannot be read.
pub fn disk(host: &dyn Host, paths: &[(&str, Option<&str>)], warn: f64) -> Vec<Disk> {
    let mut out: Vec<Disk> = Vec::new();
    let mut by_dev: IndexMap<u64, usize> = IndexMap::new();
    for (label, path) in paths {
        let Some(path) = path.filter(|p| !p.is_empty()) else {
            continue;
        };
        let at = nearest_existing(host, Path::new(path));
        let Ok(dev) = host.dev(&at) else { continue };
        let missing = at != abspath(Path::new(path));
        if let Some(&i) = by_dev.get(&dev) {
            let e = &mut out[i];
            e.fs = format!("{} and {label}", e.fs);
            e.missing |= missing;
            continue;
        }
        let Ok(s) = host.statvfs(&at) else { continue };
        let frsize = u128::from(if s.frsize != 0 { s.frsize } else { s.bsize });
        let used = u128::from(s.blocks.saturating_sub(s.bfree)) * frsize;
        let free = u128::from(s.bavail) * frsize;
        let size = u128::from(s.blocks) * frsize;
        let inodes =
            (s.files > 0).then(|| (s.files - s.ffree.min(s.files)) as f64 / s.files as f64);
        let used_ratio = (used + free > 0).then(|| used as f64 / (used + free) as f64);
        by_dev.insert(dev, out.len());
        out.push(Disk {
            fs: (*label).to_owned(),
            used_ratio,
            free_bytes: free,
            size_bytes: size,
            inodes_used_ratio: inodes,
            missing,
            warn: [used_ratio, inodes]
                .into_iter()
                .flatten()
                .any(|r| r * 100.0 >= warn),
        });
    }
    out
}

/// `(workDir, outdir, max_memory)` from the last weblog messages that carry
/// each. A relative outdir (the default `./results`) is relative to
/// Nextflow's launch folder, so it is joined to the message's `launchDir`.
pub fn metadata(events: &[Value]) -> (Option<Value>, Option<Value>, Option<Value>) {
    let (mut work, mut outdir, mut cap, mut launch): (
        Option<Value>,
        Option<Value>,
        Option<Value>,
        Option<Value>,
    ) = (None, None, None, None);
    let get = |d: Option<&Value>, k: &str| -> Option<Value> {
        match d {
            Some(Value::Object(o)) => o.get(k).filter(|v| truthy(v)).cloned(),
            _ => None,
        }
    };
    for ev in events {
        let Value::Object(e) = ev else { continue };
        let Some(md @ Value::Object(_)) = e.get("metadata") else {
            continue;
        };
        let wf = get(Some(md), "workflow");
        let params = get(Some(md), "parameters").or_else(|| get(Some(md), "params"));
        work = get(wf.as_ref(), "workDir").or(work);
        launch = get(wf.as_ref(), "launchDir").or(launch);
        outdir = get(params.as_ref(), "outdir").or(outdir);
        cap = get(params.as_ref(), "max_memory").or(cap);
    }
    if let (Some(o), Some(l)) = (&outdir, &launch) {
        let o = py_str(o).unwrap_or_default();
        if !o.starts_with('/') {
            let l = py_str(l).unwrap_or_default();
            outdir = Some(Value::Str(
                Path::new(&l).join(&o).to_string_lossy().into_owned(),
            ));
        }
    }
    (work, outdir, cap)
}

/// The weblog `process_*` messages' trace objects.
pub fn event_rows(events: &[Value]) -> Vec<&Value> {
    events
        .iter()
        .filter_map(|ev| match ev {
            Value::Object(o) => {
                let kind = o
                    .get("event")
                    .map(|k| py_str(k).unwrap_or_default())
                    .unwrap_or_default();
                match o.get("trace") {
                    Some(tr @ Value::Object(_)) if kind.starts_with("process_") => Some(tr),
                    _ => None,
                }
            }
            _ => None,
        })
        .collect()
}

/// One snapshot of the headroom.
#[derive(Debug, Clone, PartialEq)]
pub struct Snapshot {
    pub disk: Vec<Disk>,
    pub tasks: Vec<TaskMemory>,
    pub host: HostMemory,
    pub warn_used_ratio: f64,
    pub mem_warn_ratio: f64,
}

/// What [`snapshot`] reads.
pub struct Request<'a> {
    pub events: &'a [Value],
    pub trace_rows: &'a [TraceTask],
    pub work_dir: Option<&'a str>,
    pub outdir: Option<&'a str>,
    pub max_memory: Option<&'a str>,
    pub disk_warn_percent: Option<&'a str>,
    pub mem_warn_percent: Option<&'a str>,
}

/// Folders and max_memory not given come from the weblog metadata; tasks
/// come from the trace rows, then the messages.
pub fn snapshot(host: &dyn Host, steps: &[Step], r: &Request<'_>) -> Result<Snapshot, BadPercent> {
    let dw = warn_percent(r.disk_warn_percent, Some("OWNGENOME_DISK_WARN_PERCENT"))?;
    let mw = warn_percent(r.mem_warn_percent, Some("OWNGENOME_MEM_WARN_PERCENT"))?;
    let (m_work, m_out, m_cap) = metadata(r.events);
    let as_path = |v: Option<Value>| match v {
        Some(Value::Str(s)) => Some(s),
        _ => None,
    };
    let (m_work, m_out) = (as_path(m_work), as_path(m_out));
    let work = r.work_dir.filter(|w| !w.is_empty()).or(m_work.as_deref());
    let out = r.outdir.filter(|o| !o.is_empty()).or(m_out.as_deref());
    let mut rows: Vec<TaskRow<'_>> = r.trace_rows.iter().map(TaskRow::Trace).collect();
    rows.extend(event_rows(r.events).into_iter().map(TaskRow::Weblog));
    let cap_text = match r.max_memory.filter(|m| !m.is_empty()) {
        Some(m) => Some(m.to_owned()),
        None => text(m_cap.as_ref()),
    };
    let phys = host.physical_memory();
    let cap = cap_text.and_then(|c| og_nf::units::memory(Some(&c)));
    Ok(Snapshot {
        disk: disk(host, &[("work", work), ("results", out)], dw),
        tasks: task_memory(&rows, mw, steps),
        host: HostMemory {
            memory_bytes: phys,
            max_memory_bytes: cap,
            warn: matches!((phys, cap), (Some(p), Some(c)) if c > 0 && c as u128 > p),
        },
        warn_used_ratio: dw / 100.0,
        mem_warn_ratio: mw / 100.0,
    })
}

/// A size as `512 B`, `1.5 GB`, ... up to TB; `unknown` for none.
pub fn size(b: Option<f64>) -> String {
    let Some(mut b) = b else {
        return "unknown".to_owned();
    };
    for unit in ["B", "KB", "MB", "GB", "TB"] {
        if b < 1024.0 || unit == "TB" {
            return if unit == "B" {
                format!("{} {unit}", pycompat::fmt::fixed(b, 0))
            } else {
                format!("{} {unit}", pycompat::fmt::fixed(b, 1))
            };
        }
        b /= 1024.0;
    }
    unreachable!("TB ends the loop")
}

fn pct(r: f64) -> String {
    format!("{}%", pycompat::fmt::fixed(r * 100.0, 0))
}

/// The plain lines `og run headroom` prints.
pub fn lines(snap: &Snapshot) -> Vec<String> {
    let mut out = Vec::new();
    for d in &snap.disk {
        let mut s = format!(
            "disk {}: {} used, {} free of {}",
            d.fs,
            d.used_ratio.map_or_else(|| "unknown".to_owned(), pct),
            size(Some(d.free_bytes as f64)),
            size(Some(d.size_bytes as f64))
        );
        if let Some(i) = d.inodes_used_ratio {
            s += &format!(", inodes {} used", pct(i));
        }
        if d.missing {
            s += " (folder not created yet: measured on the disk that will hold it)";
        }
        if d.warn {
            s += &format!(
                "; WARN: at or above {}%",
                pycompat::fmt::fixed(snap.warn_used_ratio * 100.0, 0)
            );
        }
        out.push(s);
    }
    if snap.disk.is_empty() {
        out.push(
            "disk: no work or results folder given (--work-dir, --outdir, or weblog metadata)"
                .to_owned(),
        );
    }
    let h = &snap.host;
    let mut s = format!("memory host: {}", size(h.memory_bytes.map(|m| m as f64)));
    if let Some(m) = h.max_memory_bytes {
        s += &format!(", max_memory {}", size(Some(m as f64)));
    }
    if h.warn {
        s += "; WARN: max_memory is above this machine's memory, so a task can ask for more than there is (lower --max_memory; under WSL2 this machine is the VM, and Docker Desktop's containers get only the memory set in its settings)";
    }
    out.push(s);
    for t in &snap.tasks {
        let name = format!(
            "task {}{}",
            t.process,
            t.step
                .as_deref()
                .filter(|s| !s.is_empty())
                .map_or_else(String::new, |s| format!(" (step {s})"))
        );
        let mut s = if t.state == "running" {
            format!(
                "{name}: running, asked for {}",
                size(t.request_bytes.map(|m| m as f64))
            )
        } else if let Some(r) = t.ratio {
            let mut s = format!(
                "{name}: peak_rss {} of {} asked ({})",
                size(t.peak_rss_bytes.map(|m| m as f64)),
                size(t.request_bytes.map(|m| m as f64)),
                pct(r)
            );
            if t.warn {
                s += "; hint: close to its request (peak_rss counts shared pages, so this may be fine)";
            }
            s
        } else {
            format!(
                "{name}: {}, no ratio (request or peak_rss unknown)",
                t.state
            )
        };
        if t.oom_retry > 0 {
            s += &format!(
                "; retried {} time{} after a memory or time kill",
                t.oom_retry,
                if t.oom_retry > 1 { "s" } else { "" }
            );
        }
        out.push(s);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::host::FsStat;
    use pycompat::json::loads;
    use std::io;

    struct Fake {
        dirs: Vec<(&'static str, u64)>,
        ratio: f64,
    }

    impl Host for Fake {
        fn statvfs(&self, _: &Path) -> io::Result<FsStat> {
            let blocks = 1000;
            let used = (blocks as f64 * self.ratio) as u64;
            Ok(FsStat {
                frsize: 4096,
                bsize: 4096,
                blocks,
                bfree: blocks - used,
                bavail: blocks - used,
                files: 100,
                ffree: 50,
            })
        }
        fn dev(&self, p: &Path) -> io::Result<u64> {
            self.dirs
                .iter()
                .find(|(d, _)| Path::new(d) == p)
                .map(|(_, dev)| *dev)
                .ok_or_else(|| io::Error::other("no dev"))
        }
        fn exists(&self, p: &Path) -> bool {
            self.dirs.iter().any(|(d, _)| Path::new(d) == p)
        }
        fn physical_memory(&self) -> Option<u128> {
            Some(16 << 30)
        }
    }

    #[test]
    fn percents() {
        assert_eq!(warn_percent(Some("95"), None), Ok(95.0));
        assert_eq!(warn_percent(None, None), Ok(90.0));
        assert_eq!(
            warn_percent(Some("0"), None).unwrap_err().0,
            "warn percent must be a number from 1 to 100, got '0'"
        );
        assert!(warn_percent(Some("nan"), None).is_err());
        assert!(warn_percent(Some("x"), None).is_err());
    }

    #[test]
    fn sizes() {
        assert_eq!(size(Some(512.0)), "512 B");
        assert_eq!(size(Some(1536.0)), "1.5 KB");
        assert_eq!(size(Some(8.0 * 1024f64.powi(4) * 2048.0)), "16384.0 TB");
        assert_eq!(size(None), "unknown");
    }

    #[test]
    fn disks_share_a_file_system_and_measure_a_missing_folder_on_its_parent() {
        let host = Fake {
            dirs: vec![("/", 1), ("/w", 2), ("/r", 2)],
            ratio: 0.95,
        };
        let d = disk(
            &host,
            &[
                ("work", Some("/w/new/sub")),
                ("results", Some("/r")),
                ("x", None),
            ],
            90.0,
        );
        assert_eq!(d.len(), 1);
        assert_eq!(d[0].fs, "work and results");
        assert!(d[0].missing && d[0].warn);
        assert_eq!(d[0].inodes_used_ratio, Some(0.5));
        let ok = Fake {
            dirs: vec![("/", 1), ("/w", 2)],
            ratio: 0.5,
        };
        assert!(!disk(&ok, &[("work", Some("/w"))], 90.0)[0].warn);
    }

    #[test]
    fn the_last_attempt_counts_and_a_memory_kill_is_a_retry() {
        let evs: Vec<Value> = [
            r#"{"event":"process_completed","trace":{"task_id":1,"name":"W:A (s)","attempt":1,"status":"FAILED","exit":137,"memory":100,"peak_rss":100}}"#,
            r#"{"event":"process_completed","trace":{"task_id":2,"name":"W:A (s)","attempt":2,"status":"COMPLETED","exit":0,"memory":200,"peak_rss":190}}"#,
            r#"{"event":"process_started","trace":{"task_id":2,"name":"W:A (s)","attempt":2,"status":"RUNNING"}}"#,
            r#"{"event":"process_started","trace":{"task_id":3,"name":"W:B (s)","status":"RUNNING","memory":"2 GB"}}"#,
        ]
        .iter()
        .map(|s| loads(s).unwrap())
        .collect();
        let rows: Vec<TaskRow<'_>> = event_rows(&evs).into_iter().map(TaskRow::Weblog).collect();
        let t = task_memory(&rows, 90.0, &[]);
        assert_eq!(t.len(), 2);
        assert_eq!(
            (t[0].attempt, t[0].oom_retry, t[0].ratio, t[0].warn),
            (2, 1, Some(0.95), true)
        );
        assert_eq!(
            (t[1].state.as_str(), t[1].request_bytes),
            ("running", Some(2 << 30))
        );
    }

    #[test]
    fn a_relative_outdir_is_under_the_launch_folder() {
        let evs = vec![
            loads(r#"{"metadata":{"parameters":{"outdir":"results","max_memory":"30.GB"},"workflow":{"workDir":"/w","launchDir":"/l"}}}"#).unwrap(),
        ];
        let (w, o, c) = metadata(&evs);
        assert_eq!(w, Some(Value::Str("/w".into())));
        assert_eq!(o, Some(Value::Str("/l/results".into())));
        assert_eq!(c, Some(Value::Str("30.GB".into())));
    }
}
