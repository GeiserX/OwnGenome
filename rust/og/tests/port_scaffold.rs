//! Every verb the 1.0.0 plan adds is registered in `og`: listed on its own
//! line of `og --help` and routed to its module, never to the "unknown"
//! error of a router. A port replaces a stub's body, so this test keeps
//! holding after it; a verb that a later merge drops turns it red.
//!
//! While a verb is still a stub it must say exactly `og <verb>: not ported
//! yet` on stderr and exit 2, so no caller mistakes it for a result.

use std::io::Read;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

/// Every verb the 1.0.0 plan adds, as its words after `og`.
const PLANNED: &[&str] = &[
    "report render",
    "report index",
    "report charts",
    "live serve",
    "live results",
    "live family",
    "live demo",
    "view serve",
    "couple screen",
    "couple genes",
    "vcf index",
    "vcf shift-back",
    "setup chrm-shift",
    "setup revel",
    "setup roh-af",
    "sv survivor-prep",
    "sv consensus-check",
    "archaic compare",
    "archaic genomes",
    "ancestry-fine freqs",
    "ancestry-fine project",
    "ancestry-fine paint",
    "ancestry-fine simulate",
    "longread clair3-model",
    "map svg",
    "dev demo-sample",
    "dev demo-events",
    "dev params-doc",
    "dev python-gate",
    "dev everything-check",
    "dev licences",
    "dev docs-check",
    "dev stub-http",
    "dev steps-check",
    "dev build roh-af",
    "dev build acmg-lof",
];

/// An option no verb takes, so a ported verb stops at its argument check
/// instead of doing its work.
const PROBE: &str = "--og-scaffold-probe";

/// Runs `og ARGS` with no stdin; the exit code and stderr. A verb still
/// running after 30 s is killed and fails the test.
fn og(args: &[&str]) -> (Option<i32>, String) {
    let mut child = Command::new(env!("CARGO_BIN_EXE_og"))
        .args(args)
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .expect("og starts");
    let mut stderr = child.stderr.take().expect("piped stderr");
    let reader = std::thread::spawn(move || {
        let mut text = String::new();
        stderr.read_to_string(&mut text).map(|_| text)
    });
    let start = Instant::now();
    let status = loop {
        if let Some(status) = child.try_wait().expect("og waits") {
            break status;
        }
        if start.elapsed() > Duration::from_secs(30) {
            let _ = child.kill();
            let _ = child.wait();
            panic!("og {} {PROBE} still ran after 30 s", args.join(" "));
        }
        std::thread::sleep(Duration::from_millis(20));
    };
    let text = reader.join().expect("reader").expect("utf-8 stderr");
    (status.code(), text)
}

#[test]
fn every_planned_verb_has_its_own_usage_line() {
    let out = Command::new(env!("CARGO_BIN_EXE_og"))
        .arg("--help")
        .output()
        .expect("og --help runs");
    assert!(out.status.success());
    let usage = String::from_utf8(out.stdout).expect("utf-8 usage");
    let lines: Vec<&str> = usage
        .lines()
        .map(|l| l.trim_start_matches("usage:").trim())
        .collect();
    let missing: Vec<&&str> = PLANNED
        .iter()
        .filter(|verb| {
            let want = format!("og {verb}");
            !lines
                .iter()
                .any(|l| l == &want || l.starts_with(&format!("{want} ")))
        })
        .collect();
    assert!(missing.is_empty(), "no usage line for: {missing:?}");
}

#[test]
fn every_planned_verb_is_routed_and_a_stub_says_so() {
    let mut failures = Vec::new();
    for verb in PLANNED {
        let words: Vec<&str> = verb.split(' ').collect();
        let mut args = words.clone();
        args.push(PROBE);
        let (code, stderr) = og(&args);
        let last = words[words.len() - 1];
        let tool = words[words.len() - 2];
        let refusals = [
            format!("unknown tool '{tool}'"),
            format!("unknown {tool} tool '{last}'"),
            format!("unknown verb '{last}'"),
            format!("unknown builder '{last}'"),
            format!("invalid choice: '{last}'"),
        ];
        if let Some(r) = refusals.iter().find(|r| stderr.contains(r.as_str())) {
            failures.push(format!("og {verb}: not routed ({r})"));
        } else if stderr.contains("not ported yet") {
            let want = format!("og {verb}: not ported yet\n");
            if stderr != want || code != Some(2) {
                failures.push(format!(
                    "og {verb}: a stub must print {want:?} and exit 2, got {stderr:?} and {code:?}"
                ));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn an_unplanned_verb_is_still_refused() {
    // The control for the test above: the router's refusal is what it looks for.
    let (code, stderr) = og(&["sv", "no-such-verb", PROBE]);
    assert_eq!(code, Some(2));
    assert!(
        stderr.contains("unknown sv tool 'no-such-verb'"),
        "{stderr}"
    );
    let (code, stderr) = og(&["dev", "build", "no-such-builder", PROBE]);
    assert_eq!(code, Some(2));
    assert!(
        stderr.contains("unknown builder 'no-such-builder'"),
        "{stderr}"
    );
}
