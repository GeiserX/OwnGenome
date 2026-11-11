//! `og dev report-only [--repo DIR]`: every process whose output only the
//! report reads has the report-only error strategy. It replaces
//! `scripts/ci/check-report-only.py`.
//!
//! `nextflow.config` gives report-only processes a `withName` block whose
//! `errorStrategy` closure retries a task killed for memory or time once
//! and otherwise returns `'ignore'`, so one odd input cannot stop the run
//! before DeepVariant and the report. A process left out keeps
//! `conf/base.config`'s stop.
//!
//! A process counts as report-only when, in its `workflows/<name>.nf`,
//! every reference to `<PROCESS>.out.<channel>` is either `ch_versions =
//! ch_versions.mix(<PROCESS>.out.versions)` or a plain `ch_x =
//! <PROCESS>.out.<y>`, that `ch_x` is used nowhere else in the file but in
//! its own assignments and in the emit block, and `main.nf` uses every such
//! emit only in a line `.join(<WORKFLOW>.out.<emit>... remainder: true)`.
//! The check errs toward silence: it never asks for a process that feeds
//! another step. [`KEEP_STOP`] names the report-only processes that keep the
//! stop by choice; a name there the scan no longer finds fails, so the list
//! cannot go stale.
//!
//! Exit 0 when every report-only process is covered, 1 when one is not, 2
//! when a file cannot be read or an argument is wrong. `--repo` defaults to
//! the current folder.
//!
//! Where it differs from the Python on purpose: comments are blanked with
//! strings respected (a `//` inside a string is not a comment), and a
//! `withName` block is read to its own closing brace, so a nested brace in
//! the `errorStrategy` closure, or another setting before it, no longer
//! hides the block.

use crate::py_os_error;
use regex::Regex;
use std::collections::BTreeSet;
use std::fs;
use std::io::Write;
use std::path::Path;
use std::sync::LazyLock;

/// Report-only by their channels, kept fatal until a decision says
/// otherwise: ROH (step 17's runs of homozygosity) and CPIC_LOOKUP (the
/// CPIC recommendations) have not been through the report-only review yet.
pub const KEEP_STOP: [&str; 2] = ["CPIC_LOOKUP", "ROH"];

const USAGE: &str = "usage: og dev report-only [--repo DIR]\n";

fn read(path: &Path) -> Result<String, String> {
    let s = fs::read_to_string(path).map_err(|e| py_os_error(&e, &path.to_string_lossy()))?;
    Ok(s.replace("\r\n", "\n").replace('\r', "\n"))
}

fn re(s: &str) -> Regex {
    Regex::new(s).expect("a valid built-in pattern")
}

static PROC_REF: LazyLock<Regex> = LazyLock::new(|| re(r"\b([A-Z][A-Z0-9_]+)\.out\."));
static VERSIONS: LazyLock<Regex> = LazyLock::new(|| {
    re(r"^\s*ch_versions\s*=\s*ch_versions\.mix\(\s*(\w+)\.out\.versions\s*\)\s*$")
});
static ALIAS: LazyLock<Regex> = LazyLock::new(|| re(r"^\s*(ch_\w+)\s*=\s*(\w+)\.out\.\w+\s*$"));
static ASSIGNED: LazyLock<Regex> = LazyLock::new(|| re(r"^\s*(\w+)\s*="));
static PASSED_ON: LazyLock<Regex> = LazyLock::new(|| re(r"^\s*\w+\s*=\s*(\w+)\s*$"));
static REMAINDER: LazyLock<Regex> = LazyLock::new(|| re(r"remainder:\s*true"));

/// Python's `\w` for one character.
fn is_word(c: char) -> bool {
    c.is_alphanumeric() || c == '_'
}

/// The byte offsets where `token` stands in `line` with no word character
/// right before it (`\btoken`), and with none right after it when `whole`
/// (`\btoken\b`).
fn token_at<'a>(line: &'a str, token: &'a str, whole: bool) -> impl Iterator<Item = usize> + 'a {
    line.match_indices(token).map(|(i, _)| i).filter(move |&i| {
        let before = line[..i].chars().next_back().is_none_or(|c| !is_word(c));
        let after = line[i + token.len()..].chars().next();
        before && (!whole || after.is_none_or(|c| !is_word(c)))
    })
}

/// True when `line` holds `\b<proc>.out.\w`.
fn refers_to(line: &str, proc: &str) -> bool {
    let token = format!("{proc}.out.");
    token_at(line, &token, false)
        .any(|i| line[i + token.len()..].chars().next().is_some_and(is_word))
}

/// True when `line` holds `.join( <target>` with `remainder: true` after it.
fn report_join(line: &str, target: &str) -> bool {
    line.match_indices(".join(").any(|(i, _)| {
        let rest = line[i + ".join(".len()..].trim_start_matches(pycompat::text::is_space);
        rest.strip_prefix(target).is_some_and(|after| {
            after.chars().next().is_none_or(|c| !is_word(c)) && REMAINDER.is_match(after)
        })
    })
}

/// A report-only process: its workflow file and the emits it reaches
/// `main.nf` through.
#[derive(Debug, PartialEq, Eq)]
pub struct Found {
    pub file: String,
    pub process: String,
    pub emits: Vec<String>,
}

/// Every report-only process of `repo`, by workflow file name, then process
/// name.
pub fn report_only_processes(repo: &Path) -> Result<Vec<Found>, String> {
    let main = og_nf::groovy::strip_comments(&read(&repo.join("main.nf"))?);
    let main_lines: Vec<&str> = main.split('\n').collect();
    let wdir = repo.join("workflows");
    let mut names: Vec<String> = fs::read_dir(&wdir)
        .map_err(|e| py_os_error(&e, &wdir.to_string_lossy()))?
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with(".nf"))
        .collect();
    names.sort();
    let mut found = Vec::new();
    for fname in names {
        let Some(wf) = og_nf::workflow::read(&read(&wdir.join(&fname))?) else {
            continue;
        };
        let lines: Vec<&str> = wf.text.split('\n').collect();
        let procs: BTreeSet<&str> = PROC_REF
            .captures_iter(&wf.text)
            .map(|c| c.get(1).expect("group").as_str())
            .collect();
        for proc in procs {
            let mut used = Vec::new();
            let mut ok = true;
            for line in lines.iter().filter(|l| refers_to(l, proc)) {
                if VERSIONS.captures(line).is_some_and(|c| &c[1] == proc) {
                    continue;
                }
                let Some(ch) = ALIAS
                    .captures(line)
                    .filter(|c| &c[2] == proc)
                    .map(|c| c[1].to_owned())
                else {
                    ok = false;
                    break;
                };
                let own = |l: &str| {
                    ASSIGNED.captures(l).is_some_and(|c| c[1] == ch)
                        || PASSED_ON.captures(l).is_some_and(|c| c[1] == ch)
                };
                let others = lines
                    .iter()
                    .any(|l| token_at(l, &ch, true).next().is_some() && !own(l));
                let emits: Vec<&String> = wf
                    .emits
                    .iter()
                    .filter(|(_, c)| *c == ch)
                    .map(|(e, _)| e)
                    .collect();
                if others || emits.is_empty() {
                    ok = false;
                    break;
                }
                for emit in emits {
                    let target = format!("{}.out.{emit}", wf.name);
                    let in_main: Vec<&&str> = main_lines
                        .iter()
                        .filter(|l| token_at(l, &target, true).next().is_some())
                        .collect();
                    if in_main.is_empty() || !in_main.iter().all(|l| report_join(l, &target)) {
                        ok = false;
                    }
                    used.push(emit.clone());
                }
            }
            if ok && !used.is_empty() {
                found.push(Found {
                    file: fname.clone(),
                    process: proc.to_owned(),
                    emits: used,
                });
            }
        }
    }
    Ok(found)
}

/// Runs `og dev report-only ARGS`; returns the exit code.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let repo = match args {
        [] => ".".to_owned(),
        [flag, dir] if flag == "--repo" => dir.clone(),
        [flag] if flag.starts_with("--repo=") => flag["--repo=".len()..].to_owned(),
        [flag] if flag == "-h" || flag == "--help" => {
            let _ = write!(out, "{USAGE}");
            return 0;
        }
        _ => {
            let _ = write!(err, "{USAGE}");
            let _ = writeln!(
                err,
                "og dev report-only: error: unrecognized arguments: {}",
                args.join(" ")
            );
            return 2;
        }
    };
    let repo = Path::new(&repo);
    let config = match read(&repo.join("nextflow.config")) {
        Ok(c) => c,
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            return 2;
        }
    };
    // PLANTED: the Python's comment cut at any // and its one-level closure regex.
    let config: String = config
        .lines()
        .map(|l| l.split("//").next().unwrap_or(""))
        .collect::<Vec<_>>()
        .join("\n");
    let ignore =
        Regex::new(r"withName:\s*'([^']*)'\s*\{\s*errorStrategy\s*=\s*\{[^}]*'ignore'[^}]*\}")
            .unwrap();
    let planted: Vec<String> = ignore
        .captures_iter(&config)
        .map(|c| c[1].to_owned())
        .collect();
    let mut patterns = Vec::new();
    for w in og_nf::config::with_names(&config) {
        if !planted.contains(&w.selector) {
            continue;
        }
        if !w.error_strategy.is_some_and(|s| s.contains("'ignore'")) {
            continue;
        }
        match Regex::new(&format!("^(?:{})$", w.selector)) {
            Ok(r) => patterns.push(r),
            Err(e) => {
                let _ = writeln!(
                    err,
                    "ERROR: withName '{}' is not a pattern: {e}",
                    w.selector
                );
                return 2;
            }
        }
    }
    if patterns.is_empty() {
        let _ = writeln!(
            err,
            "ERROR: no report-only withName block in nextflow.config"
        );
        return 1;
    }
    let found = match report_only_processes(repo) {
        Ok(f) => f,
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            return 2;
        }
    };
    let mut bad = 0;
    for f in &found {
        let covered = patterns.iter().any(|p| p.is_match(&f.process));
        if KEEP_STOP.contains(&f.process.as_str()) {
            let _ = writeln!(
                out,
                "keep-stop: {} (workflows/{}), by choice",
                f.process, f.file
            );
        } else if covered {
            let _ = writeln!(out, "ok: {} (workflows/{})", f.process, f.file);
        } else {
            let _ = writeln!(
                err,
                "ERROR: {} (workflows/{}) feeds only the report ({}, joined with remainder: true in main.nf) but no report-only withName block in nextflow.config names it, so one failed task stops the whole run",
                f.process,
                f.file,
                f.emits.join(", ")
            );
            bad += 1;
        }
    }
    for k in KEEP_STOP {
        if !found.iter().any(|f| f.process == k) {
            let _ = writeln!(
                err,
                "ERROR: KEEP_STOP names {k}, which the scan no longer finds as report-only; drop it from the list"
            );
            bad += 1;
        }
    }
    if found.is_empty() {
        let _ = writeln!(
            err,
            "ERROR: the scan found no report-only process at all; the patterns no longer match the code"
        );
        bad += 1;
    }
    u8::from(bad > 0)
}
