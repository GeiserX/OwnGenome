//! Regenerate the public demo's invented Nextflow weblog from assets/steps.tsv.

use crate::demo::{RUN_STARTED, utc};
use og_runlog::steps::{Step, read_steps};
use pycompat::json::{Dump, Value, dumps, loads};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

const FIXTURE: &str = "tests/fixtures/live/demo-001.weblog.jsonl";
const USAGE: &str = "usage: og dev demo-events [--check] [--out FILE] [--repo DIR]";
const GIB: i64 = 1 << 30;

fn duration(step: &Step) -> i64 {
    if let Some(high) = step.hours
        && step.hours_from == "measured"
    {
        let low = step.hours_low.parse::<f64>().unwrap_or(high);
        return (((low + high) / 2.0) * 3600.0).round_ties_even() as i64;
    }
    600
}

fn finish<'a>(
    step: &'a Step,
    steps: &BTreeMap<&'a str, &'a Step>,
    ends: &mut BTreeMap<&'a str, i64>,
    active: &mut BTreeSet<&'a str>,
) -> Result<i64, String> {
    let name = step.step.as_str();
    if let Some(end) = ends.get(name) {
        return Ok(*end);
    }
    if !active.insert(name) {
        return Err(format!("the demo step graph has a cycle at step {name}"));
    }
    let mut start = RUN_STARTED;
    for (&other, upstream) in steps {
        // The report reads all the final sections. Step 38 is a separate family run.
        if (name == "24" && other != name) || step.upstream.iter().any(|s| s == other) {
            start = start.max(finish(upstream, steps, ends, active)?);
        }
    }
    active.remove(name);
    let end = start + duration(step);
    ends.insert(name, end);
    Ok(end)
}

fn iso(t: i64) -> String {
    utc(t).replace(" UTC", "Z").replace(' ', "T")
}

fn object<const N: usize>(fields: [(&str, Value); N]) -> Value {
    Value::Object(fields.into_iter().map(|(k, v)| (k.into(), v)).collect())
}

fn text(s: &str) -> Value {
    Value::Str(s.to_owned())
}

fn number(n: i64) -> Value {
    Value::Int(i128::from(n))
}

fn message(event: &str, t: i64, detail: (&str, Value)) -> Value {
    object([
        ("runName", text("demo_run")),
        ("runId", text("00000000-0000-0000-0000-00000000d001")),
        ("event", text(event)),
        ("utcTime", Value::Str(iso(t))),
        detail,
    ])
}

struct Tasks {
    next_id: i64,
    events: Vec<(i64, Value)>,
}

impl Tasks {
    fn task(
        &mut self,
        process: &str,
        submit: i64,
        start: i64,
        end: i64,
        attempt: i64,
        failed: bool,
    ) {
        self.next_id += 1;
        let id = self.next_id;
        let long = end - start > 3600;
        let memory = if long { 8 * GIB } else { 2 * GIB };
        let mut trace = indexmap::IndexMap::from([
            ("task_id".into(), number(id)),
            ("hash".into(), Value::Str(format!("de/{id:06}"))),
            (
                "name".into(),
                Value::Str(format!("OWNGENOME:{process} (DEMO-001)")),
            ),
            ("process".into(), Value::Str(format!("OWNGENOME:{process}"))),
            ("tag".into(), text("DEMO-001")),
            ("attempt".into(), number(attempt)),
            ("cpus".into(), number(if long { 4 } else { 2 })),
            ("memory".into(), number(memory)),
            ("submit".into(), number(submit * 1000)),
        ]);
        trace.insert("status".into(), text("SUBMITTED"));
        self.events.push((
            submit,
            message(
                "process_submitted",
                submit,
                ("trace", Value::Object(trace.clone())),
            ),
        ));
        let started = start + 1;
        trace.insert("status".into(), text("RUNNING"));
        trace.insert("start".into(), number(started * 1000));
        self.events.push((
            started,
            message(
                "process_started",
                started,
                ("trace", Value::Object(trace.clone())),
            ),
        ));
        let peak_percent = match (failed, process) {
            (true, _) => 99,
            (false, "HLA_TYPING") => 92,
            _ => 40,
        };
        trace.insert(
            "status".into(),
            text(if failed { "FAILED" } else { "COMPLETED" }),
        );
        trace.insert("exit".into(), number(if failed { 137 } else { 0 }));
        trace.insert("complete".into(), number(end * 1000));
        trace.insert("realtime".into(), number((end - started) * 1000));
        trace.insert("peak_rss".into(), number(memory * peak_percent / 100));
        self.events.push((
            end,
            message("process_completed", end, ("trace", Value::Object(trace))),
        ));
    }
}

fn render(catalogue: &Path) -> Result<String, String> {
    let steps = read_steps(catalogue).map_err(|e| e.to_string())?;
    let steps: BTreeMap<_, _> = steps
        .iter()
        .filter(|s| !s.processes.is_empty() && s.step != "38")
        .map(|s| (s.step.as_str(), s))
        .collect();
    let mut ends = BTreeMap::new();
    let mut schedule = Vec::new();
    let mut active = BTreeSet::new();
    for step in steps.values() {
        let end = finish(step, &steps, &mut ends, &mut active)?;
        schedule.push((end - duration(step), step.step.as_str(), end));
    }
    schedule.sort_unstable();
    let mut done = BTreeSet::new();
    let mut tasks = Tasks {
        next_id: 0,
        events: Vec::new(),
    };
    for (start, step, end) in schedule {
        // PRS is shared by steps 25 and 26 and runs only once.
        let processes: Vec<_> = steps[step]
            .processes
            .iter()
            .filter(|p| done.insert(*p))
            .collect();
        let count = processes.len() as i64;
        if count == 0 {
            continue;
        }
        let span = (end - start) / count;
        for (i, process) in processes.into_iter().enumerate() {
            let a = start + i as i64 * span;
            let b = if i as i64 + 1 == count { end } else { a + span };
            if process == "CLINVAR_SCREEN" {
                let retry = a + (b - a) / 3;
                tasks.task(process, start, a, retry, 1, true);
                tasks.task(process, retry, retry + 5, b, 2, false);
            } else {
                tasks.task(process, start, a, b, 1, false);
            }
        }
    }
    // A failed attempt and its resubmission may share a second. Keep their order.
    tasks.events.sort_by_key(|(time, _)| *time);
    let end = tasks
        .events
        .last()
        .ok_or("the step catalogue schedules no demo tasks")?
        .0;
    let started = loads(r#"{"parameters":{"outdir":"/genome","max_cpus":16,"max_memory":"12.GB"},"workflow":{"commandLine":"nextflow run main.nf -resume"}}"#)
        .expect("valid invented metadata");
    let completed = object([("workflow", object([("success", Value::Bool(true))]))]);
    let events = std::iter::once(message("started", RUN_STARTED, ("metadata", started)))
        .chain(tasks.events.into_iter().map(|(_, event)| event))
        .chain([
            message("completed", end + 30, ("metadata", completed)),
            object([
                ("event", text("pgp_run_finished")),
                ("exit", number(0)),
                ("utcTime", Value::Str(iso(end + 60))),
            ]),
        ]);
    let options = Dump {
        sort_keys: true,
        separators: Some((",", ":")),
        ..Dump::default()
    };
    let mut output = String::new();
    for event in events {
        output.push_str(&dumps(&event, &options).map_err(|e| e.to_string())?);
        output.push('\n');
    }
    Ok(output)
}

struct Args {
    repo: PathBuf,
    output: PathBuf,
    check: bool,
}

impl Args {
    fn parse(args: &[String]) -> Result<Self, String> {
        let mut repo = None;
        let mut output = None;
        let mut check = false;
        let mut args = args.iter();
        while let Some(flag) = args.next() {
            match flag.as_str() {
                "--check" => check = true,
                "--out" | "--repo" => {
                    let path = args
                        .next()
                        .filter(|v| !v.starts_with("--"))
                        .ok_or_else(|| format!("{flag} takes a path"))?;
                    if flag == "--out" {
                        output = Some(PathBuf::from(path));
                    } else {
                        repo = Some(PathBuf::from(path));
                    }
                }
                _ => return Err(format!("unknown argument {flag}")),
            }
        }
        let repo = repo
            .or_else(super::demo::find_repo)
            .ok_or("cannot find the checkout; give --repo")?;
        let output = output.unwrap_or_else(|| repo.join(FIXTURE));
        Ok(Self {
            repo,
            output,
            check,
        })
    }

    fn run(&self, out: &mut dyn Write) -> Result<(), String> {
        let generated = render(&self.repo.join("assets/steps.tsv"))?;
        if self.check {
            if fs::read(&self.output).ok().as_deref() != Some(generated.as_bytes()) {
                return Err(format!(
                    "{} differs from the generated weblog; run og dev demo-events",
                    self.output.display()
                ));
            }
            writeln!(out, "{}: up to date", self.output.display())
        } else {
            fs::write(&self.output, &generated)
                .map_err(|e| format!("{}: {e}", self.output.display()))?;
            writeln!(
                out,
                "Written: {} ({} messages)",
                self.output.display(),
                generated.lines().count()
            )
        }
        .map_err(|e| e.to_string())
    }
}

/// Write the deterministic public fixture, or check it against the catalogue.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if matches!(args.first().map(String::as_str), Some("-h" | "--help")) {
        let _ = writeln!(out, "{USAGE}");
        return 0;
    }
    let args = match Args::parse(args) {
        Ok(args) => args,
        Err(e) => {
            let _ = writeln!(err, "{USAGE}\nog dev demo-events: {e}");
            return 2;
        }
    };
    match args.run(out) {
        Ok(()) => 0,
        Err(e) => {
            let _ = writeln!(err, "og dev demo-events: {e}");
            1
        }
    }
}
