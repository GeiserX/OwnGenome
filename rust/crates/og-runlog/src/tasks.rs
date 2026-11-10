//! One sample's tasks, from the trace and from weblog messages.

use indexmap::IndexMap;
use jiff::tz::TimeZone;
use og_nf::task_name;
use pycompat::json::Value;
use pycompat::value::{py_str, truthy};
use std::path::Path;

/// Processes that run once per run with a tag that is not a sample id (the
/// reference basename, the score set): never a sample's task, even when
/// the tag equals the id.
pub const SHARED_PROCESSES: [&str; 2] = ["MINIMAP2_INDEX", "PRS_PREPARE"];

/// A task status that means it ended well.
pub fn done_ok(s: &str) -> bool {
    matches!(s, "COMPLETED" | "CACHED")
}

/// A task status that means it ended badly.
pub fn done_bad(s: &str) -> bool {
    matches!(s, "FAILED" | "ABORTED")
}

/// One task as the run state sees it.
#[derive(Debug, Clone, PartialEq)]
pub struct Task {
    /// The task id (or its hash, or its name when neither is known).
    pub id: String,
    /// The full name, `WORKFLOW:PROCESS (tag)`: attempts of one task share it.
    pub name: String,
    pub process: String,
    /// Nextflow's status, upper case: SUBMITTED, RUNNING, COMPLETED, CACHED,
    /// FAILED, ABORTED.
    pub status: String,
    /// Epoch seconds.
    pub start: Option<f64>,
    pub end: Option<f64>,
}

/// A task is the sample's when its tag is the sample id, or the id and
/// `:<part>`.
pub fn is_sample(tag: Option<&str>, sample: &str) -> bool {
    tag.is_some_and(|t| t == sample || t.strip_prefix(sample).is_some_and(|r| r.starts_with(':')))
}

/// The rows of the trace files tagged with the sample. A trace start is the
/// `submit` time in the zone `tz`; its end is that plus the `duration`. A
/// trace that cannot be read gives nothing.
pub fn trace_tasks(paths: &[impl AsRef<Path>], sample: &str, tz: &TimeZone) -> Vec<Task> {
    let mut tasks = Vec::new();
    for p in paths {
        let _ = og_nf::for_each_row(p.as_ref(), |r| {
            let Some(name) = r.get("name") else { return };
            let Some((proc, tag)) = task_name(name) else {
                return;
            };
            if SHARED_PROCESSES.contains(&proc) || !is_sample(tag, sample) {
                return;
            }
            let start = r
                .get("submit")
                .and_then(|v| og_nf::time::parse_local(v, tz));
            let dur = r.get("duration").and_then(og_nf::units::duration_text);
            let id = [r.get("task_id"), r.get("hash"), Some(name)]
                .into_iter()
                .flatten()
                .find(|v| !v.is_empty())
                .unwrap_or(name);
            tasks.push(Task {
                id: id.to_owned(),
                name: name.to_owned(),
                process: proc.to_owned(),
                status: pycompat::text::upper(r.get("status").unwrap_or("")),
                start,
                end: start.zip(dur).map(|(s, d)| s + d),
            });
        });
    }
    tasks
}

/// `d.get(key)` of an object, `None` for anything else.
fn field<'a>(d: &'a Value, key: &str) -> Option<&'a Value> {
    match d {
        Value::Object(o) => o.get(key),
        _ => None,
    }
}

/// A positive weblog time in milliseconds as epoch seconds.
fn ms(v: Option<&Value>) -> Option<f64> {
    let x = match v? {
        Value::Int(i) => *i as f64,
        Value::Float(x) => *x,
        Value::Bool(b) => f64::from(u8::from(*b)),
        _ => return None,
    };
    (x > 0.0).then(|| x / 1000.0)
}

/// `str(v)` of a field, `"None"` for a missing one.
fn text(v: Option<&Value>) -> String {
    v.map_or_else(|| "None".to_owned(), |v| py_str(v).unwrap_or_default())
}

/// The weblog `process_*` messages of the sample's tasks, in file order.
/// A message's status is the event (`process_submitted`, `process_started`)
/// or its trace's status; a `process_completed` with no final status is
/// COMPLETED on exit 0, else FAILED.
pub fn event_tasks(events: &[Value], sample: &str) -> Vec<Task> {
    let mut tasks = Vec::new();
    for ev in events {
        let kind = field(ev, "event");
        let Some(tr @ Value::Object(_)) = field(ev, "trace") else {
            continue;
        };
        if !text(kind).starts_with("process_") {
            continue;
        }
        let name = match field(tr, "name").filter(|v| truthy(v)) {
            Some(Value::Str(n)) => n.clone(),
            Some(_) => continue,
            None => match field(tr, "tag").filter(|v| truthy(v)) {
                Some(tag) => format!("{} ({})", text(field(tr, "process")), text(Some(tag))),
                None => match field(tr, "process") {
                    Some(Value::Str(p)) => p.clone(),
                    _ => continue,
                },
            },
        };
        let Some((proc, tag)) = task_name(&name) else {
            continue;
        };
        if SHARED_PROCESSES.contains(&proc) || !is_sample(tag, sample) {
            continue;
        }
        let kind = match kind {
            Some(Value::Str(k)) => k.as_str(),
            _ => "",
        };
        let mut st = match kind {
            "process_submitted" => "SUBMITTED".to_owned(),
            "process_started" => "RUNNING".to_owned(),
            _ => pycompat::text::upper(
                &field(tr, "status")
                    .filter(|v| truthy(v))
                    .map(|v| py_str(v).unwrap_or_default())
                    .unwrap_or_default(),
            ),
        };
        if kind == "process_completed" && !done_ok(&st) && !done_bad(&st) {
            st = if text(field(tr, "exit")) == "0" {
                "COMPLETED"
            } else {
                "FAILED"
            }
            .to_owned();
        }
        let id = [field(tr, "task_id"), field(tr, "hash")]
            .into_iter()
            .flatten()
            .find(|v| truthy(v))
            .map_or_else(|| name.clone(), |v| py_str(v).unwrap_or_default());
        tasks.push(Task {
            id,
            process: proc.to_owned(),
            status: st,
            start: ms(field(tr, "start")).or_else(|| ms(field(tr, "submit"))),
            end: ms(field(tr, "complete")),
            name,
        });
    }
    tasks
}

/// One task per task id (the trace first, then the messages in order: the
/// later wins, except that a late `started` does not undo a finished
/// task), then one per task name: the attempt with the highest id is the
/// task's outcome (a retry replaces a failed attempt).
pub fn latest_tasks(trace: Vec<Task>, events: Vec<Task>) -> Vec<Task> {
    let mut by_id: IndexMap<String, Task> = IndexMap::new();
    for t in trace.into_iter().chain(events) {
        let merged = match by_id.get(&t.id) {
            Some(old) => {
                let finished = |s: &str| done_ok(s) || done_bad(s);
                if finished(&old.status) && !finished(&t.status) {
                    continue;
                }
                Task {
                    start: t.start.or(old.start),
                    end: t.end.or(old.end),
                    ..t
                }
            }
            None => t,
        };
        by_id.insert(merged.id.clone(), merged);
    }
    let mut all: Vec<Task> = by_id.into_values().collect();
    // Numeric ids first, in number order, then the rest in text order.
    all.sort_by(|a, b| {
        let key = |t: &Task| pycompat::num::parse_int(&t.id).ok();
        match (key(a), key(b)) {
            (Some(x), Some(y)) => x.cmp(&y),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => a.id.cmp(&b.id),
        }
    });
    let mut by_name: IndexMap<String, Task> = IndexMap::new();
    for t in all {
        by_name.insert(t.name.clone(), t);
    }
    by_name.into_values().collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use pycompat::json::loads;

    fn task(id: &str, name: &str, status: &str) -> Task {
        Task {
            id: id.into(),
            name: name.into(),
            process: "P".into(),
            status: status.into(),
            start: None,
            end: None,
        }
    }

    #[test]
    fn sample_tags() {
        assert!(is_sample(Some("S1"), "S1"));
        assert!(is_sample(Some("S1:chr1"), "S1"));
        assert!(!is_sample(Some("S10"), "S1"));
        assert!(!is_sample(None, "S1"));
    }

    #[test]
    fn the_retry_with_the_higher_id_wins() {
        let t = latest_tasks(
            vec![
                task("10", "P (s)", "COMPLETED"),
                task("9", "P (s)", "FAILED"),
            ],
            vec![],
        );
        assert_eq!(t.len(), 1);
        assert_eq!(t[0].status, "COMPLETED");
        // 10 > 9 as numbers, though "10" < "9" as text.
    }

    #[test]
    fn a_late_start_does_not_undo_an_end() {
        let mut done = task("1", "P (s)", "COMPLETED");
        done.end = Some(5.0);
        let t = latest_tasks(vec![done], vec![task("1", "P (s)", "RUNNING")]);
        assert_eq!(t[0].status, "COMPLETED");
        let mut run = task("1", "P (s)", "RUNNING");
        run.start = Some(1.0);
        let t = latest_tasks(vec![run], vec![task("1", "P (s)", "COMPLETED")]);
        assert_eq!((t[0].status.as_str(), t[0].start), ("COMPLETED", Some(1.0)));
    }

    #[test]
    fn weblog_messages() {
        let evs: Vec<Value> = [
            r#"{"event":"process_submitted","trace":{"task_id":1,"name":"W:A (S1)","submit":1000}}"#,
            r#"{"event":"process_completed","trace":{"task_id":1,"process":"W:A","tag":"S1","exit":0,"complete":5000}}"#,
            r#"{"event":"process_completed","trace":{"task_id":2,"name":"W:B (S2)","exit":0}}"#,
            r#"{"event":"process_completed","trace":{"task_id":3,"name":"W:MINIMAP2_INDEX (S1)"}}"#,
            r#"{"event":"started","trace":{"task_id":4,"name":"W:C (S1)"}}"#,
        ]
        .iter()
        .map(|s| loads(s).unwrap())
        .collect();
        let t = event_tasks(&evs, "S1");
        assert_eq!(t.len(), 2);
        assert_eq!(
            (t[0].id.as_str(), t[0].status.as_str(), t[0].start),
            ("1", "SUBMITTED", Some(1.0))
        );
        assert_eq!(
            (t[1].name.as_str(), t[1].status.as_str(), t[1].end),
            ("W:A (S1)", "COMPLETED", Some(5.0))
        );
    }
}
