//! `--rerun`, `--refresh-clinvar` and `--launch-dir`.
//!
//! The locks are `flock(2)` locks, the kind `flock(1)` takes, so they
//! interlock with `scripts/run-all.sh` and step 06.

use crate::files::append_log;
use crate::{CLINVAR_LOCK, Check, now_epoch, now_utc};
use indexmap::IndexMap;
use og_runlog::state::{RunFiles, RunState, Sources, State};
use og_runlog::steps::Step;
use pycompat::json::Value;
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::process::Command;

fn open_lock(path: &Path) -> std::io::Result<File> {
    if let Some(d) = path.parent() {
        std::fs::create_dir_all(d)?;
    }
    OpenOptions::new().append(true).create(true).open(path)
}

fn files(sample_dir: &Path) -> std::io::Result<RunFiles> {
    RunFiles::read(
        sample_dir,
        &Sources {
            traces: &[],
            events: None,
            stale: false,
        },
    )
}

fn state<'s>(files: &RunFiles, sample: &str, steps: &'s [Step]) -> State<'s> {
    files.state(sample, steps, now_epoch(), &jiff::tz::TimeZone::system())
}

/// Whether a run of the sample is going on (no end record).
pub fn running(sample: &str, sample_dir: &Path, steps: &[Step]) -> bool {
    files(sample_dir).is_ok_and(|f| state(&f, sample, steps).run == RunState::Running)
}

/// The steps of which a task ran in this run (done, not taken from the
/// cache).
pub fn redone_steps(st: &State<'_>) -> Vec<String> {
    st.steps
        .iter()
        .filter(|s| s.processes.iter().any(|(_, c)| c.done > 0))
        .map(|s| s.step.step.clone())
        .collect()
}

/// What the rerun redid against what the check named.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Comparison {
    pub predicted: Vec<String>,
    pub redone: Vec<String>,
    pub predicted_not_redone: Vec<String>,
    pub redone_not_predicted: Vec<String>,
    /// Named script steps (no Nextflow process): they skip when their output
    /// exists, so a rerun does not redo them.
    pub not_redone_script_steps: Vec<String>,
}

impl Comparison {
    pub fn to_json(&self) -> Value {
        let a = |v: &[String]| Value::Array(v.iter().map(|s| Value::Str(s.clone())).collect());
        Value::Object(IndexMap::from([
            ("predicted".to_owned(), a(&self.predicted)),
            ("redone".to_owned(), a(&self.redone)),
            (
                "predicted_not_redone".to_owned(),
                a(&self.predicted_not_redone),
            ),
            (
                "redone_not_predicted".to_owned(),
                a(&self.redone_not_predicted),
            ),
            (
                "not_redone_script_steps".to_owned(),
                a(&self.not_redone_script_steps),
            ),
        ]))
    }
}

pub fn compare(predicted: &[String], redone: Vec<String>, steps: &[Step]) -> Comparison {
    let script: Vec<String> = predicted
        .iter()
        .filter(|n| {
            steps
                .iter()
                .find(|s| s.step == **n)
                .is_none_or(|s| s.processes.is_empty())
        })
        .cloned()
        .collect();
    Comparison {
        predicted: predicted.to_vec(),
        predicted_not_redone: predicted
            .iter()
            .filter(|n| !redone.contains(n) && !script.contains(n))
            .cloned()
            .collect(),
        redone_not_predicted: redone
            .iter()
            .filter(|n| !predicted.contains(n))
            .cloned()
            .collect(),
        redone,
        not_redone_script_steps: script,
    }
}

/// How a replay ended.
pub enum Rerun {
    /// It ran; its exit code.
    Ran(i32),
    /// Refused before anything ran: the exit code and why.
    Refused(u8, String),
}

/// Replays the last `run-all.sh` command of the sample (`meta command` of
/// `logs/run_status.tsv`) with the same switches: the switch variables the
/// run read (`meta switches`) are dropped from the environment, so the
/// command's own assignments decide. One `--rerun` per sample at a time:
/// `recheck/rerun.lock` is held from the running check through the end of
/// the replay. Records the rerun in `state` and the log.
pub fn rerun(
    sample: &str,
    sample_dir: &Path,
    steps: &[Step],
    state_: &mut Check,
    log: &Path,
) -> std::io::Result<Rerun> {
    let lock = open_lock(&sample_dir.join("recheck").join("rerun.lock"))?;
    if lock.try_lock().is_err() {
        return Ok(Rerun::Refused(
            3,
            format!("another --rerun of {sample} is going on: not rerunning."),
        ));
    }
    let before = files(sample_dir)?;
    if state(&before, sample, steps).run == RunState::Running {
        return Ok(Rerun::Refused(
            3,
            format!(
                "a run of {sample} is going on (logs/run_status.tsv has no finished_utc): not rerunning. \
                 When it has ended, run this again."
            ),
        ));
    }
    let meta = before.status.map(|s| s.meta).unwrap_or_default();
    let Some(cmd) = meta.get("command").filter(|c| !c.is_empty()).cloned() else {
        return Ok(Rerun::Refused(
            2,
            "logs/run_status.tsv names no run-all.sh command (a run of an older version, or a bare \
             Nextflow run): run the same command again by hand."
                .into(),
        ));
    };
    let switches = meta.get("switches").cloned().unwrap_or_default();
    let started = now_utc();
    // After --refresh-clinvar the new file is installed, so it is rerun-ready.
    let predicted = state_.rerun_ready_steps.clone();
    println!("Rerunning: {cmd}");
    let _ = std::io::stdout().flush();
    append_log(log, &[[&started, "rerun-start", "", "", &cmd]])?;
    let mut c = Command::new("bash");
    c.arg("-c").arg(&cmd);
    for name in switches.split_whitespace() {
        c.env_remove(name);
    }
    let rc = c.status().map(|s| s.code().unwrap_or(1)).unwrap_or(127);
    let after = files(sample_dir)?;
    let cmp = compare(
        &predicted,
        redone_steps(&state(&after, sample, steps)),
        steps,
    );
    let finished = now_utc();
    state_.rerun = Some(Value::Object(IndexMap::from([
        ("started_utc".to_owned(), Value::Str(started)),
        ("finished_utc".to_owned(), Value::Str(finished.clone())),
        ("rc".to_owned(), Value::Int(i128::from(rc))),
        ("command".to_owned(), Value::Str(cmd)),
        ("predicted_vs_redone".to_owned(), cmp.to_json()),
    ])));
    append_log(
        log,
        &[[
            &finished,
            "rerun-end",
            "",
            &format!("rc {rc}"),
            &format!(
                "predicted {}; redone {}",
                cmp.predicted.join(","),
                cmp.redone.join(",")
            ),
        ]],
    )?;
    let none = |v: &[String]| {
        if v.is_empty() {
            "none".to_owned()
        } else {
            v.join(", ")
        }
    };
    println!(
        "Predicted steps: {}; redone: {}",
        none(&cmp.predicted),
        none(&cmp.redone)
    );
    for (list, what) in [
        (&cmp.predicted_not_redone, "named but not redone"),
        (&cmp.redone_not_predicted, "redone but not named"),
        (
            &cmp.not_redone_script_steps,
            "script steps not redone (they skip when their output exists; delete it and run the script)",
        ),
    ] {
        if !list.is_empty() {
            println!("  {what}: {}", list.join(", "));
        }
    }
    drop(lock);
    Ok(Rerun::Ran(rc))
}

/// The last command in `DIR/.nextflow/history`, with `-resume`.
pub fn launch_dir_command(dir: &Path) -> Option<String> {
    let raw = std::fs::read(dir.join(".nextflow").join("history")).ok()?;
    let text = String::from_utf8_lossy(&raw);
    let last = text.lines().rfind(|l| !l.trim().is_empty())?;
    let cmd = last.rsplit('\t').next().unwrap_or(last).trim();
    let words =
        shlex::split(cmd).unwrap_or_else(|| cmd.split_whitespace().map(str::to_owned).collect());
    Some(if words.iter().any(|w| w == "-resume") {
        cmd.to_owned()
    } else {
        format!("{cmd} -resume")
    })
}

/// The samples of `genome_dir` (a folder with `logs/run_status.tsv`)
/// whose run is going on.
pub fn running_samples(genome_dir: &Path, steps: &[Step]) -> Vec<String> {
    let Ok(rd) = std::fs::read_dir(genome_dir) else {
        return vec![];
    };
    let mut names: Vec<String> = rd
        .flatten()
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
        .into_iter()
        .filter(|n| {
            let d = genome_dir.join(n);
            d.join("logs").join("run_status.tsv").is_file() && running(n, &d, steps)
        })
        .collect()
}

/// `setup.sh --refresh clinvar GENOME_DIR` under the ClinVar lock of
/// GENOME_DIR. The ClinVar files are shared by every sample, so a run of any
/// of them blocks the refresh: `run-all.sh` and step 06 hold the lock shared
/// while they read the files, and wait while this holds it. A run with no
/// lock (a bare Nextflow run, an older checkout) is caught by its
/// `logs/run_status.tsv`. The lock is let go before `--rerun` starts
/// `run-all.sh`, which takes it shared. Returns the exit code and, on a
/// refusal, why.
pub fn refresh_clinvar(
    genome_dir: &Path,
    setup: &Path,
    steps: &[Step],
) -> std::io::Result<(i32, Option<String>)> {
    let path = genome_dir.join(CLINVAR_LOCK);
    let lock = open_lock(&path)?;
    if lock.try_lock().is_err() {
        return Ok((
            3,
            Some(format!(
                "a run or step 06 holds {}: not refreshing the ClinVar files it reads.",
                path.display()
            )),
        ));
    }
    let busy = running_samples(genome_dir, steps);
    if !busy.is_empty() {
        return Ok((
            3,
            Some(format!(
                "a run of {} is going on: not refreshing the ClinVar files it reads.",
                busy.join(", ")
            )),
        ));
    }
    let rc = Command::new(setup)
        .args(["--refresh", "clinvar"])
        .arg(genome_dir)
        .status()
        .map(|s| s.code().unwrap_or(1))
        .unwrap_or(127);
    drop(lock);
    Ok((
        rc,
        (rc != 0).then(|| "setup.sh --refresh clinvar failed: see above.".to_owned()),
    ))
}
