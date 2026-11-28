//! `og dev steps-check`: `assets/steps.tsv` against the code that runs the
//! steps. Each row of the catalogue must agree with what runs, and each step
//! that runs must have its row:
//!
//! - every step `scripts/run-all.sh` plans has its row, with the label and
//!   the `--tools` name it plans it under;
//! - the rows marked `recorded=yes` are the steps run-all.sh writes into
//!   `logs/run_status.tsv`, and the summary's list of those steps
//!   (`RUN_ALL_STEPS` of og-summary) is the same set;
//! - every summary.json section belongs to the row of the step that writes
//!   it, and no row names a section the summary does not have;
//! - every process of `modules/local` (and every name a workflow includes a
//!   process under) belongs to exactly one row, and every process a row
//!   names is one of them;
//! - step 17's cDNA read check, which `workflows/clinical.nf` includes under
//!   its own name, counts under step 17, and step 06's under step 06;
//! - the rows' tools are `main.nf`'s known tools;
//! - every measured row has an `hours_low` below its `hours`, and no other
//!   row has one;
//! - one row per step, and every upstream step has a row;
//! - every row is a step something runs: run-all.sh plans or records it, or
//!   it names a process. A row for a step that nothing runs would show a
//!   step on the live page that can never start.
//!
//! These were the catalogue checks of the Python test of the run state
//! reader; the last one is new.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use regex::Regex;

const USAGE: &str = "usage: og dev steps-check [--root DIR]\n";

/// One row of the catalogue, as the checks read it.
#[derive(Debug, Clone)]
pub struct Row {
    pub step: String,
    pub label: String,
    pub tool: String,
    pub recorded: bool,
    pub processes: Vec<String>,
    pub sections: Vec<String>,
    pub upstream: Vec<String>,
    pub hours: Option<f64>,
    pub hours_from: String,
    pub hours_low: String,
}

/// What the checks compare the catalogue with: the code's side.
#[derive(Debug, Clone)]
pub struct Code {
    /// `scripts/run-all.sh`.
    pub run_all: String,
    /// `main.nf`.
    pub main_nf: String,
    /// `workflows/clinical.nf`.
    pub clinical: String,
    /// The process names of `modules/local`, and the names `main.nf` and
    /// `workflows/*.nf` include a process under.
    pub processes: BTreeSet<String>,
    /// Each summary.json section and the step it belongs to.
    pub sections: Vec<(String, String)>,
    /// The steps the summary treats as recorded by run-all.sh.
    pub run_all_steps: Vec<String>,
}

fn re(p: &str) -> Regex {
    Regex::new(p).expect("a valid regex")
}

fn read(root: &Path, rel: &str) -> Result<String, String> {
    fs::read_to_string(root.join(rel)).map_err(|e| format!("{rel}: {e}"))
}

/// The catalogue at `path`: og run state's reading of it, plus the
/// `recorded` column.
pub fn read_rows(path: &Path) -> Result<Vec<Row>, String> {
    let typed = og_runlog::steps::read_steps(path).map_err(|e| e.to_string())?;
    let text = fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut lines = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty());
    let head: Vec<&str> = lines.next().unwrap_or_default().split('\t').collect();
    let rec = head
        .iter()
        .position(|h| *h == "recorded")
        .ok_or_else(|| format!("{}: no column 'recorded'", path.display()))?;
    let recorded: Vec<bool> = lines
        .map(|l| l.split('\t').nth(rec) == Some("yes"))
        .collect();
    if recorded.len() != typed.len() {
        return Err(format!(
            "{}: {} rows read, {} recorded cells",
            path.display(),
            typed.len(),
            recorded.len()
        ));
    }
    Ok(typed
        .into_iter()
        .zip(recorded)
        .map(|(s, recorded)| Row {
            step: s.step,
            label: s.label,
            tool: s.tool,
            recorded,
            processes: s.processes,
            sections: s.sections,
            upstream: s.upstream,
            hours: s.hours,
            hours_from: s.hours_from,
            hours_low: s.hours_low,
        })
        .collect())
}

fn nf_files(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(rd) = fs::read_dir(dir) else { return };
    let mut entries: Vec<PathBuf> = rd.filter_map(|e| e.ok().map(|e| e.path())).collect();
    entries.sort();
    for p in entries {
        if p.is_dir() {
            nf_files(&p, out);
        } else if p.extension().is_some_and(|x| x == "nf") {
            out.push(p);
        }
    }
}

/// The processes `modules/local` defines, and the names a workflow
/// includes one under (`include { A as B }` runs, and is counted, as B).
fn module_processes(root: &Path) -> Result<BTreeSet<String>, String> {
    let process = re(r"(?m)^process\s+([A-Z0-9_]+)");
    let include = re(r"(?m)^\s*include\s*\{([^}]*)\}");
    let alias = re(r"\w+\s+as\s+([A-Z0-9_]+)");
    let mut names = BTreeSet::new();
    let mut modules = Vec::new();
    nf_files(&root.join("modules/local"), &mut modules);
    for f in modules {
        let t = fs::read_to_string(&f).map_err(|e| format!("{}: {e}", f.display()))?;
        names.extend(process.captures_iter(&t).map(|c| c[1].to_owned()));
    }
    let mut wf = vec![root.join("main.nf")];
    let mut more = Vec::new();
    nf_files(&root.join("workflows"), &mut more);
    wf.extend(
        more.into_iter()
            .filter(|p| p.parent() == Some(&root.join("workflows"))),
    );
    for f in wf {
        let t = fs::read_to_string(&f).map_err(|e| format!("{}: {e}", f.display()))?;
        for inc in include.captures_iter(&t) {
            names.extend(alias.captures_iter(&inc[1]).map(|c| c[1].to_owned()));
        }
    }
    Ok(names)
}

/// Each summary.json section and its step: the sections of the summary
/// og summary collect writes for an empty sample folder.
fn summary_sections(root: &Path) -> Result<Vec<(String, String)>, String> {
    use pycompat::json::Value;
    let dir = std::env::temp_dir().join(format!("og-steps-check-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    let summary = og_summary::collect::post::collect(&og_summary::collect::Input {
        sample: "S".into(),
        sample_dir: dir.to_string_lossy().into_owned(),
        declared_sex: None,
        assets: root.join("assets").to_string_lossy().into_owned(),
        group: og_summary::collect::Group::All,
        tables: root.join("bin").to_string_lossy().into_owned(),
    });
    let _ = fs::remove_dir_all(&dir);
    let summary = summary?;
    let field = |v: &Value, k: &str| match v {
        Value::Object(m) => m.get(k).cloned(),
        _ => None,
    };
    let Some(Value::Object(sections)) = field(&summary, "sections") else {
        return Err("og summary collect wrote no sections".into());
    };
    sections
        .iter()
        .map(|(k, v)| match field(v, "step") {
            Some(Value::Str(s)) => Ok((k.clone(), s)),
            _ => Err(format!("summary section {k} has no step")),
        })
        .collect()
}

/// The code's side of the checks, read from the checkout at `root`.
pub fn read_code(root: &Path) -> Result<Code, String> {
    Ok(Code {
        run_all: read(root, "scripts/run-all.sh")?,
        main_nf: read(root, "main.nf")?,
        clinical: read(root, "workflows/clinical.nf")?,
        processes: module_processes(root)?,
        sections: summary_sections(root)?,
        run_all_steps: og_summary::collect::model::RUN_ALL_STEPS
            .iter()
            .map(|s| (*s).to_owned())
            .collect(),
    })
}

/// (step, label, tool) of every step run-all.sh plans, and the other steps
/// it records: those it marks skipped by hand, the scripts it runs by name
/// (`TIDDIT:04a-tiddit`) and the four variant callers it globs. A line
/// commented out is no step, so full-line comments are dropped first.
fn run_all_steps(src: &str) -> (BTreeSet<(String, String, String)>, BTreeSet<String>) {
    let code: Vec<&str> = src
        .lines()
        .filter(|l| !l.trim_start().starts_with('#'))
        .collect();
    let src = code.join("\n");
    let src = src.as_str();
    let mut planned = BTreeSet::new();
    for c in re(r#""(\d+b?) ([^":]+):(\w+)""#).captures_iter(src) {
        planned.insert((c[1].to_owned(), c[2].to_owned(), c[3].to_owned()));
    }
    for c in re(r#"plan "(\d+b?) ([^"]+)" (\w+)"#).captures_iter(src) {
        planned.insert((c[1].to_owned(), c[2].to_owned(), c[3].to_owned()));
    }
    let steps: BTreeSet<&str> = planned.iter().map(|(s, _, _)| s.as_str()).collect();
    let mut other: BTreeSet<String> = re(r"printf 'step\\t(\d+b?)\\tskipped")
        .captures_iter(src)
        .map(|c| c[1].to_owned())
        .filter(|s| !steps.contains(s.as_str()))
        .collect();
    other.extend(
        re(r"\b\w+:(\d+[ab]?)-[\w-]+\b")
            .captures_iter(src)
            .map(|c| c[1].to_owned()),
    );
    if re(r#"compgen -G "03\[([a-d-]+)\]"#)
        .captures(src)
        .is_some_and(|c| &c[1] == "a-d")
    {
        other.extend(["03a", "03b", "03c", "03d"].map(str::to_owned));
    }
    (planned, other)
}

/// The quoted words of `def known_tools = [ ... ]` in main.nf.
fn known_tools(main_nf: &str) -> BTreeSet<String> {
    let Some(block) = re(r"(?s)def known_tools = \[(.*?)\]").captures(main_nf) else {
        return BTreeSet::new();
    };
    re(r"'(\w+)'")
        .captures_iter(&block[1])
        .map(|c| c[1].to_owned())
        .collect()
}

/// A check's name and, when it fails, why.
pub type Outcome = (&'static str, Result<(), String>);

fn outcome(name: &'static str, ok: bool, why: impl FnOnce() -> String) -> Outcome {
    (name, if ok { Ok(()) } else { Err(why()) })
}

/// Every check of the catalogue `rows` against `code`, in order.
pub fn checks(rows: &[Row], code: &Code) -> Vec<Outcome> {
    let mut out = Vec::new();
    let row: BTreeMap<&str, &Row> = rows.iter().map(|r| (r.step.as_str(), r)).collect();
    let (planned, other) = run_all_steps(&code.run_all);

    let bad: Vec<String> = planned
        .iter()
        .filter(|(n, l, t)| {
            row.get(n.as_str())
                .is_none_or(|r| r.label != *l || r.tool != *t)
        })
        .map(|(n, l, t)| {
            let r = row.get(n.as_str());
            format!(
                "run-all.sh plans {n} '{l}' {t}, the row says {:?} {:?}",
                r.map(|r| &r.label),
                r.map(|r| &r.tool)
            )
        })
        .collect();
    out.push(outcome(
        "every run-all.sh plan step has its row, label and tool",
        planned.len() >= 40 && bad.is_empty(),
        || {
            if bad.is_empty() {
                format!("only {} plan steps read in run-all.sh", planned.len())
            } else {
                bad.join("; ")
            }
        },
    ));

    let recorded: BTreeSet<&str> = rows
        .iter()
        .filter(|r| r.recorded)
        .map(|r| r.step.as_str())
        .collect();
    let want: BTreeSet<&str> = planned
        .iter()
        .map(|(s, _, _)| s.as_str())
        .chain(other.iter().map(String::as_str))
        .collect();
    out.push(outcome(
        "the recorded steps are what run-all.sh records",
        recorded == want,
        || {
            format!(
                "differ in {:?}",
                recorded.symmetric_difference(&want).collect::<Vec<_>>()
            )
        },
    ));
    let summary_steps: BTreeSet<&str> = code.run_all_steps.iter().map(String::as_str).collect();
    out.push(outcome(
        "the summary's recorded steps are the recorded rows",
        summary_steps == recorded,
        || {
            format!(
                "differ in {:?}",
                summary_steps
                    .symmetric_difference(&recorded)
                    .collect::<Vec<_>>()
            )
        },
    ));

    // Every row that names a section, so a wrong row is seen even when the
    // right row names it too.
    let mut owner: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for r in rows {
        for k in &r.sections {
            owner.entry(k.as_str()).or_default().push(r.step.as_str());
        }
    }
    let keys: BTreeMap<&str, &str> = code
        .sections
        .iter()
        .map(|(k, s)| (k.as_str(), s.as_str()))
        .collect();
    let wrong: Vec<String> = keys
        .iter()
        .filter(|(k, s)| owner.get(*k).map(Vec::as_slice) != Some(&[**s][..]))
        .map(|(k, s)| format!("{k} is step {s}'s, the catalogue says {:?}", owner.get(k)))
        .collect();
    let unknown: Vec<&str> = owner
        .keys()
        .filter(|k| !keys.contains_key(*k))
        .copied()
        .collect();
    out.push(outcome(
        "every summary section belongs to the row of its step",
        wrong.is_empty() && unknown.is_empty() && !keys.is_empty(),
        || format!("{wrong:?}; sections the summary does not have: {unknown:?}"),
    ));

    let mut procs: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for r in rows {
        for p in &r.processes {
            procs.entry(p.as_str()).or_default().push(r.step.as_str());
        }
    }
    let mods: BTreeSet<&str> = code.processes.iter().map(String::as_str).collect();
    let lost: Vec<&str> = mods
        .iter()
        .filter(|p| !procs.contains_key(*p))
        .copied()
        .collect();
    // Step 26 runs inside PRS, so rows 25 and 26 both name it, and no
    // other row may.
    let twice: Vec<String> = procs
        .iter()
        .filter(|(p, v)| v.len() > 1 && !(**p == "PRS" && v[..] == ["25", "26"]))
        .map(|(p, v)| format!("{p} in {v:?}"))
        .collect();
    let extra: Vec<&str> = procs
        .keys()
        .filter(|p| !mods.contains(*p))
        .copied()
        .collect();
    out.push(outcome(
        "every module process belongs to one step",
        mods.len() > 100 && lost.is_empty() && twice.is_empty() && extra.is_empty(),
        || {
            format!(
                "{} processes; in no row: {lost:?}; in two rows: {twice:?}; no such process: {extra:?}",
                mods.len()
            )
        },
    ));

    // Step 17 runs step 06's read check on CPSR's rows: workflows/clinical.nf
    // includes it under its own name, so that task counts under step 17, and
    // step 06's own call under step 06.
    let step_of: BTreeMap<&str, &str> = rows
        .iter()
        .flat_map(|r| {
            r.processes
                .iter()
                .map(move |p| (p.as_str(), r.step.as_str()))
        })
        .collect();
    let of = |task: &str| og_nf::task_name(task).and_then(|(p, _)| step_of.get(p).copied());
    let got = (
        of("CLINICAL:CPSR_CDNA_CHECK (S1)"),
        of("PGX:CLINVAR_CDNA_CHECK (S1)"),
    );
    let included = re(r"(?m)^include \{ CLINVAR_CDNA_CHECK as CPSR_CDNA_CHECK \}")
        .is_match(&code.clinical)
        && code.clinical.contains("CPSR_CDNA_CHECK(");
    out.push(outcome(
        "CPSR's read check counts under step 17, ClinVar's under step 06",
        got == (Some("17"), Some("06")) && included,
        || format!("steps {got:?}, clinical.nf includes and calls it: {included}"),
    ));

    let tools: BTreeSet<String> = rows
        .iter()
        .filter(|r| !r.tool.is_empty())
        .map(|r| r.tool.clone())
        .collect();
    let known = known_tools(&code.main_nf);
    out.push(outcome(
        "the tools are main.nf's known tools",
        tools == known && !known.is_empty(),
        || {
            format!(
                "differ in {:?}",
                tools.symmetric_difference(&known).collect::<Vec<_>>()
            )
        },
    ));

    let low_ok = |r: &Row| {
        if r.hours_from != "measured" {
            return r.hours_low.is_empty();
        }
        match (pycompat::num::parse_float(&r.hours_low), r.hours) {
            (Ok(lo), Some(h)) => 0.0 < lo && lo < h,
            _ => false,
        }
    };
    let bad_low: Vec<String> = rows
        .iter()
        .filter(|r| !low_ok(r))
        .map(|r| {
            format!(
                "{} low {:?} hours {:?} ({})",
                r.step, r.hours_low, r.hours, r.hours_from
            )
        })
        .collect();
    let measured = rows.iter().filter(|r| r.hours_from == "measured").count();
    out.push(outcome(
        "every measured row has hours_low below hours, no other row has one",
        measured >= 25 && bad_low.is_empty(),
        || format!("{measured} measured rows; {bad_low:?}"),
    ));

    let nums: Vec<&str> = rows.iter().map(|r| r.step.as_str()).collect();
    let dup: BTreeSet<&str> = nums
        .iter()
        .filter(|n| nums.iter().filter(|m| m == n).count() > 1)
        .copied()
        .collect();
    let orphan: Vec<String> = rows
        .iter()
        .flat_map(|r| r.upstream.iter().map(move |u| (r, u)))
        .filter(|(_, u)| !row.contains_key(u.as_str()))
        .map(|(r, u)| format!("{} reads {u}", r.step))
        .collect();
    out.push(outcome(
        "one row per step, every upstream step has a row",
        dup.is_empty() && orphan.is_empty(),
        || format!("twice: {dup:?}; no row: {orphan:?}"),
    ));

    let idle: Vec<&str> = rows
        .iter()
        .filter(|r| !want.contains(r.step.as_str()) && r.processes.is_empty())
        .map(|r| r.step.as_str())
        .collect();
    out.push(outcome(
        "every row is a step run-all.sh runs or a process of main.nf",
        idle.is_empty(),
        || format!("nothing runs step {idle:?}"),
    ));
    out
}

/// Every check on the checkout at `root`.
pub fn check(root: &Path) -> Result<Vec<Outcome>, String> {
    let rows = read_rows(&root.join("assets/steps.tsv"))?;
    let code = read_code(root)?;
    Ok(checks(&rows, &code))
}

/// `og dev steps-check ARGS`: 0 when the catalogue agrees with the code, 1
/// when a check fails, 2 on a usage or read error.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let mut root = PathBuf::from(".");
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match (a.as_str(), it.clone().next()) {
            ("--root", Some(v)) => {
                root = PathBuf::from(v);
                it.next();
            }
            ("-h" | "--help", _) => {
                let _ = write!(out, "{USAGE}");
                return 0;
            }
            _ => {
                let _ = write!(err, "{USAGE}");
                let _ = writeln!(err, "og dev steps-check: error: unexpected argument '{a}'");
                return 2;
            }
        }
    }
    let results = match check(&root) {
        Ok(r) => r,
        Err(e) => {
            let _ = writeln!(err, "og dev steps-check: {e}");
            return 2;
        }
    };
    let mut failed = false;
    for (name, r) in results {
        let _ = match r {
            Ok(()) => writeln!(out, "ok      {name}"),
            Err(why) => {
                failed = true;
                writeln!(out, "FAIL    {name}: {why}")
            }
        };
    }
    u8::from(failed)
}
