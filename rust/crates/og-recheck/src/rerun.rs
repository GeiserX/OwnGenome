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

/// The switches `scripts/run-all.sh` records (its `SWITCHES` array; a unit
/// test keeps the two equal). A replay sets or drops only these, whatever
/// `meta switches` says, so an edited row can never hand `BASH_ENV`, `PATH`,
/// `LD_PRELOAD` or `DYLD_*` to the `bash` it starts.
pub const SWITCHES: &[&str] = &[
    "GENOME_DIR",
    "SKIP_VALIDATION",
    "THREADS",
    "SKIP_TRIM",
    "INTERVALS",
    "DEEPVARIANT_DISABLE_SMALL_MODEL",
    "ALIGN_DIR",
    "TOOLS",
    "KIR",
    "HLA_QC",
    "HPO",
    "REF_FASTA",
    "EH_CATALOG",
    "MANTA_CALL_REGIONS",
    "PARASCOPY_POPULATION",
    "PARAPHASE_REGIONS",
    "ANCESTRY_PANEL",
    "ANCESTRY_FINE_DIR",
    "ARCHAIC_DIR",
    "ARCHAIC_POP",
    "PGSC_CALC_DIR",
    "ANCIENT_PANEL_DIR",
    "VIRAL_REF",
    "MEGANE_KMERS",
    "MOCHA_DATA",
    "GRIDSS",
    "IMPUTATION",
    "SOMATIC",
    "EXTRA_CALLERS",
    "BENCHMARK",
    "FAMILIES",
    "OWNGENOME_LIVE_URL",
    "LOCITYPER_DATA",
    "LOCITYPER_BG_REGION",
    "NXF_VER",
    "NXF_WORK",
    "NXF_HOME",
    "NXF_OFFLINE",
];

/// The names of `meta switches` that are also in [`SWITCHES`].
fn known(switches: &str) -> Vec<&str> {
    switches
        .split_whitespace()
        .filter(|n| SWITCHES.contains(n))
        .collect()
}

/// One shell word of a recorded command.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Token {
    /// A word: its text, and how many of its first bytes stood unquoted.
    Word(String, usize),
    /// The `&&` between the `cd` and the call.
    And,
}

/// Bytes `printf %q` leaves unquoted anywhere in a word, and that mean
/// nothing to the shell there. `#` and `~` do too, except first in a word
/// (and `~` after an unquoted `=` or `:`, see [`words`]).
fn plain(b: u8, first: bool) -> bool {
    b.is_ascii_alphanumeric()
        || b"%+-./:=@_,".contains(&b)
        || (!first && (b == b'#' || b == b'~'))
        || b >= 0x80
}

/// The text of a `$'...'` quote from just after its `$'`; `None` for an
/// escape `printf %q` never writes. Moves `i` past the closing quote.
fn ansi_c(s: &[u8], i: &mut usize, out: &mut Vec<u8>) -> Option<()> {
    loop {
        let b = *s.get(*i)?;
        *i += 1;
        match b {
            b'\'' => return Some(()),
            b'\\' => {
                let e = *s.get(*i)?;
                *i += 1;
                out.push(match e {
                    b'a' => 7,
                    b'b' => 8,
                    b'e' | b'E' => 27,
                    b'f' => 12,
                    b'n' => b'\n',
                    b'r' => b'\r',
                    b't' => b'\t',
                    b'v' => 11,
                    b'\\' | b'\'' | b'"' | b'?' => e,
                    b'0'..=b'7' => {
                        let mut v = u32::from(e - b'0');
                        for _ in 0..2 {
                            match s.get(*i) {
                                Some(d @ b'0'..=b'7') => {
                                    v = v * 8 + u32::from(d - b'0');
                                    *i += 1;
                                }
                                _ => break,
                            }
                        }
                        u8::try_from(v).ok()?
                    }
                    b'x' => {
                        let mut v = 0u32;
                        let mut n = 0;
                        while n < 2
                            && let Some(d) = s.get(*i).and_then(|d| char::from(*d).to_digit(16))
                        {
                            v = v * 16 + d;
                            *i += 1;
                            n += 1;
                        }
                        if n == 0 {
                            return None;
                        }
                        u8::try_from(v).ok()?
                    }
                    _ => return None,
                });
            }
            _ => out.push(b),
        }
    }
}

/// The words of a command as `run-all.sh` writes it with `printf %q`:
/// plain bytes, `\x`, `'...'` and `$'...'`, words split on spaces and tabs,
/// and the `&&` operator. `None` for anything else a shell would act on (a
/// `;`, a pipe, a redirect, an expansion, a glob, a double quote, a line
/// break, a `~` first in a word or after an unquoted `=` or `:`, where bash
/// expands it) and for text that is not UTF-8.
fn words(cmd: &str) -> Option<Vec<Token>> {
    let s = cmd.as_bytes();
    let (mut out, mut i) = (Vec::new(), 0);
    while i < s.len() {
        if matches!(s[i], b' ' | b'\t') {
            i += 1;
            continue;
        }
        if s[i..].starts_with(b"&&") && s.get(i + 2).is_none_or(|b| matches!(*b, b' ' | b'\t')) {
            out.push(Token::And);
            i += 2;
            continue;
        }
        let (mut w, mut lead, mut quoted) = (Vec::new(), 0, false);
        // The byte before, when it stood unquoted.
        let mut prev: Option<u8> = None;
        while i < s.len() && !matches!(s[i], b' ' | b'\t') {
            let b = s[i];
            if b == b'\\' {
                let ch = cmd.get(i + 1..)?.chars().next()?;
                if ch == '\n' {
                    return None;
                }
                let mut buf = [0u8; 4];
                w.extend_from_slice(ch.encode_utf8(&mut buf).as_bytes());
                i += 1 + ch.len_utf8();
                quoted = true;
                prev = None;
            } else if b == b'\'' {
                let end = s[i + 1..].iter().position(|&c| c == b'\'')? + i + 1;
                w.extend_from_slice(&s[i + 1..end]);
                i = end + 1;
                quoted = true;
                prev = None;
            } else if b == b'$' && s.get(i + 1) == Some(&b'\'') {
                i += 2;
                ansi_c(s, &mut i, &mut w)?;
                quoted = true;
                prev = None;
            } else if plain(b, w.is_empty() && !quoted)
                && !(b == b'~' && matches!(prev, Some(b'=' | b':')))
            {
                w.push(b);
                prev = Some(b);
                i += 1;
                if !quoted {
                    lead = w.len();
                }
            } else {
                return None;
            }
        }
        out.push(Token::Word(String::from_utf8(w).ok()?, lead));
    }
    Some(out)
}

/// A recorded `run-all.sh` call as typed parts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Call {
    /// The directory it was run from.
    pub dir: String,
    /// The switches it set, as `NAME`, value.
    pub env: Vec<(String, String)>,
    /// The `run-all.sh` it ran.
    pub script: String,
    /// Its arguments.
    pub args: Vec<String>,
}

/// `meta command` as a typed call, or `None` when it is not the shape
/// `run-all.sh` writes: `cd DIR && [NAME=VALUE ...] .../scripts/run-all.sh
/// ARG ...`, every `NAME` one of `switches` and of [`SWITCHES`]. Nothing of
/// it reaches a shell.
pub fn parse_call(cmd: &str, switches: &str) -> Option<Call> {
    let names = known(switches);
    let mut t = words(cmd)?.into_iter();
    let mut word = || match t.next() {
        Some(Token::Word(w, lead)) => Some((w, lead)),
        _ => None,
    };
    let (cd, lead) = word()?;
    if cd != "cd" || lead != 2 {
        return None;
    }
    let (dir, _) = word()?;
    if !matches!(t.next(), Some(Token::And)) {
        return None;
    }
    let mut env = Vec::new();
    let script = loop {
        let (w, lead) = match t.next() {
            Some(Token::Word(w, lead)) => (w, lead),
            _ => return None,
        };
        match w.find('=') {
            Some(eq) if eq < lead && names.contains(&&w[..eq]) => {
                env.push((w[..eq].to_owned(), w[eq + 1..].to_owned()));
            }
            _ => break w,
        }
    };
    if !script.ends_with("/scripts/run-all.sh") {
        return None;
    }
    let args = t
        .map(|x| match x {
            Token::Word(w, _) => Some(w),
            Token::And => None,
        })
        .collect::<Option<Vec<String>>>()?;
    Some(Call {
        dir,
        env,
        script,
        args,
    })
}

/// Replays the last `run-all.sh` command of the sample (`meta command` of
/// `logs/run_status.tsv`) with the same switches: the switch variables the
/// run read (`meta switches`, only names in [`SWITCHES`]) are dropped from
/// the environment, so the
/// command's own assignments decide. The command is read as a typed call
/// ([`parse_call`]) and `run-all.sh` is started by `bash` with its
/// arguments one by one; a command of any other shape is refused and
/// nothing runs. One `--rerun` per sample at a time: `recheck/rerun.lock`
/// is held from the running check through the end of the replay. Records
/// the rerun in `state` and the log.
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
    let Some(call) = parse_call(&cmd, &switches) else {
        return Ok(Rerun::Refused(
            2,
            "the command in logs/run_status.tsv is not a run-all.sh call as run-all.sh records it \
             (cd DIR && [SWITCH=VALUE ...] .../scripts/run-all.sh ARGS): not running it. Check it, \
             then run it by hand."
                .into(),
        ));
    };
    if !Path::new(&call.dir).is_dir() {
        return Ok(Rerun::Refused(
            2,
            format!(
                "the run's directory {} is gone: run the command by hand from where it should run.",
                call.dir
            ),
        ));
    }
    let started = now_utc();
    // After --refresh-clinvar the new file is installed, so it is rerun-ready.
    let predicted = state_.rerun_ready_steps.clone();
    println!("Rerunning: {cmd}");
    let _ = std::io::stdout().flush();
    append_log(log, &[[&started, "rerun-start", "", "", &cmd]])?;
    let mut c = Command::new("bash");
    c.arg(&call.script).args(&call.args).current_dir(&call.dir);
    for name in known(&switches) {
        c.env_remove(name);
    }
    for (name, value) in &call.env {
        c.env(name, value);
    }
    c.env("PWD", &call.dir);
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
