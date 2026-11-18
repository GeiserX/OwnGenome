//! The help rules every `og` command keeps (docs/rust.md, "Command
//! lines"):
//!
//! - every command on a line of `og --help` or `og dev --help` answers
//!   `--help` and `-h` with its usage on stdout and exit 0, and writes no
//!   file while doing so;
//! - `og --help` starts with how to run a genome, and the `og dev` lines
//!   come only in their own block at the end, which `og dev --help` prints
//!   alone;
//! - `og clinvar --help` names the defaults of its options.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

struct Out {
    code: Option<i32>,
    stdout: String,
    stderr: String,
}

/// Runs `og ARGS` in `cwd` with no stdin; killed and failed after 30 s.
fn og_in(cwd: &Path, args: &[&str]) -> Out {
    let mut child = Command::new(env!("CARGO_BIN_EXE_og"))
        .args(args)
        .current_dir(cwd)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("og starts");
    let mut so = child.stdout.take().expect("piped stdout");
    let mut se = child.stderr.take().expect("piped stderr");
    let ro = std::thread::spawn(move || {
        let mut b = Vec::new();
        so.read_to_end(&mut b).map(|_| b)
    });
    let re = std::thread::spawn(move || {
        let mut b = Vec::new();
        se.read_to_end(&mut b).map(|_| b)
    });
    let start = Instant::now();
    let status = loop {
        if let Some(status) = child.try_wait().expect("og waits") {
            break status;
        }
        if start.elapsed() > Duration::from_secs(30) {
            let _ = child.kill();
            let _ = child.wait();
            panic!("og {} still ran after 30 s", args.join(" "));
        }
        std::thread::sleep(Duration::from_millis(10));
    };
    let text = |b: Vec<u8>| String::from_utf8_lossy(&b).into_owned();
    Out {
        code: status.code(),
        stdout: text(ro.join().expect("reader").expect("stdout")),
        stderr: text(re.join().expect("reader").expect("stderr")),
    }
}

/// A new empty directory of this test.
fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-cli-conventions-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).expect("scratch dir");
    d
}

/// The `og ...` lines of a help text, without `usage:` and indentation.
fn og_lines(text: &str) -> Vec<String> {
    text.lines()
        .map(|l| l.trim_start_matches("usage:").trim())
        .filter(|l| l.starts_with("og "))
        .map(str::to_owned)
        .collect()
}

fn is_verb(w: &str) -> bool {
    w.starts_with(|c: char| c.is_ascii_lowercase())
        && w.chars()
            .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '-')
}

/// The commands a usage line names: the tool alone, and the tool with each
/// of the verbs of its second word.
fn commands(line: &str) -> Vec<Vec<String>> {
    let w: Vec<&str> = line.split_whitespace().collect();
    // A root option with its own line (`og --tree-hash`) is a command too.
    let Some(tool) = w
        .get(1)
        .filter(|t| is_verb(t) || t.strip_prefix("--").is_some_and(is_verb))
    else {
        return vec![];
    };
    let mut out = vec![vec![(*tool).to_owned()]];
    if let Some(verbs) = w.get(2).filter(|v| v.split('|').all(is_verb)) {
        for v in verbs.split('|') {
            out.push(vec![(*tool).to_owned(), v.to_owned()]);
        }
    }
    out
}

#[test]
fn every_command_answers_help_on_stdout_and_writes_nothing() {
    let here = scratch("help");
    let root = og_in(&here, &["--help"]);
    let dev = og_in(&here, &["dev", "--help"]);
    assert_eq!(root.code, Some(0), "og --help: {}", root.stderr);
    let mut lines = og_lines(&root.stdout);
    lines.extend(og_lines(&dev.stdout));
    let mut seen = std::collections::BTreeSet::new();
    let mut failures = Vec::new();
    for line in &lines {
        for cmd in commands(line) {
            if !seen.insert(cmd.clone()) {
                continue;
            }
            for flag in ["--help", "-h"] {
                let mut args: Vec<&str> = cmd.iter().map(String::as_str).collect();
                args.push(flag);
                let o = og_in(&here, &args);
                if o.code != Some(0) || o.stdout.trim().is_empty() {
                    failures.push(format!(
                        "og {}: exit {:?}, {} bytes on stdout, stderr {:?}",
                        args.join(" "),
                        o.code,
                        o.stdout.len(),
                        o.stderr.lines().last().unwrap_or("")
                    ));
                }
                let left: Vec<_> = std::fs::read_dir(&here)
                    .expect("scratch dir")
                    .map(|e| e.expect("entry").file_name())
                    .collect();
                if !left.is_empty() {
                    failures.push(format!("og {} wrote {left:?}", args.join(" ")));
                    let _ = std::fs::remove_dir_all(&here);
                    std::fs::create_dir_all(&here).expect("scratch dir");
                }
            }
        }
    }
    assert!(
        seen.len() > 100,
        "only {} commands read from the usage",
        seen.len()
    );
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    let _ = std::fs::remove_dir_all(&here);
}

#[test]
fn the_help_fallback_names_the_verb_and_an_unknown_tool_is_still_refused() {
    let here = scratch("fallback");
    // A verb answered from the root usage gets its own line, not the tool's.
    let o = og_in(&here, &["sv", "consensus-check", "--help"]);
    assert_eq!(o.code, Some(0));
    assert_eq!(o.stdout, "usage: og sv consensus-check [args]\n");
    // A help word after `--` is an argument, not a request.
    let o = og_in(&here, &["sv", "consensus-check", "--", "--help"]);
    assert_ne!(o.code, Some(0));
    // The control: no line, no help.
    let o = og_in(&here, &["no-such-tool", "--help"]);
    assert_eq!(o.code, Some(2));
    assert!(
        o.stderr.contains("unknown tool 'no-such-tool'"),
        "{}",
        o.stderr
    );
    let _ = std::fs::remove_dir_all(&here);
}

#[test]
fn og_help_starts_with_how_to_run_a_genome() {
    let here = scratch("root");
    let root = og_in(&here, &["--help"]).stdout;
    let head: Vec<&str> = root.lines().take(15).collect();
    assert!(
        head.iter()
            .any(|l| l.contains("./scripts/run-all.sh SAMPLE")),
        "the first 15 lines of og --help do not say how to start a run:\n{}",
        head.join("\n")
    );
    assert!(
        head.iter()
            .any(|l| l.contains("./scripts/validate-setup.sh"))
    );
    let run = root
        .lines()
        .position(|l| l.trim().starts_with("og run state"))
        .expect("an og run line");
    let what = root.lines().nth(run + 1).unwrap_or("");
    assert!(
        what.contains("progress") && what.contains("never starts a run"),
        "og run is not described as inspecting a run: {what:?}"
    );
    // Every og dev line is in the last block, and og dev --help prints that
    // block alone.
    let block = root
        .lines()
        .position(|l| l.starts_with("Developer and CI tools"))
        .expect("a developer block");
    let before: Vec<&str> = root
        .lines()
        .take(block)
        .filter(|l| l.trim().starts_with("og dev "))
        .collect();
    assert!(
        before.is_empty(),
        "og dev lines before their block: {before:?}"
    );
    let in_root: Vec<String> = og_lines(&root)
        .into_iter()
        .filter(|l| l.starts_with("og dev "))
        .collect();
    let dev = og_in(&here, &["dev", "--help"]);
    assert_eq!(dev.code, Some(0));
    let in_dev = og_lines(&dev.stdout);
    assert!(in_root.len() > 20, "{in_root:?}");
    let missing: Vec<&String> = in_root.iter().filter(|l| !in_dev.contains(l)).collect();
    assert!(missing.is_empty(), "not in og dev --help: {missing:?}");
    let other: Vec<&String> = in_dev
        .iter()
        .filter(|l| !l.starts_with("og dev "))
        .collect();
    assert!(
        other.is_empty(),
        "og dev --help shows other tools: {other:?}"
    );
    let _ = std::fs::remove_dir_all(&here);
}

#[test]
fn og_clinvar_help_names_its_defaults() {
    let here = scratch("clinvar");
    for args in [
        &["clinvar", "--help"][..],
        &["clinvar", "hits", "--help"],
        &["clinvar", "deep-intronic", "-h"],
    ] {
        let o = og_in(&here, args);
        assert_eq!(o.code, Some(0), "og {}: {}", args.join(" "), o.stderr);
        for want in [
            "og clinvar join-adjacent",
            "og clinvar hits",
            "og clinvar deep-intronic",
        ] {
            assert!(o.stdout.contains(want), "{want} not in:\n{}", o.stdout);
        }
        for (opt, default) in [
            ("--min-gq N", "30"),
            ("--min-dp N", "10"),
            ("--window N", "10"),
        ] {
            let line = o
                .stdout
                .lines()
                .find(|l| l.trim_start().starts_with(opt))
                .unwrap_or_else(|| panic!("no {opt} line in:\n{}", o.stdout));
            assert!(line.ends_with(&format!("(default {default})")), "{line}");
        }
    }
    let _ = std::fs::remove_dir_all(&here);
}
