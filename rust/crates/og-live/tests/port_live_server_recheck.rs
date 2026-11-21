//! The database re-check on the live site, and the ClinVar lock that
//! `og recheck --refresh-clinvar` interlocks with. The check itself (its
//! classes, steps, exit codes and files, `--affected`) is tested by `cargo
//! test -p og-recheck`.
//!
//! 1. `plugins.recheck` of `state.json` gives each step its badges ("newer
//!    ClinVar available", "redone <date>") and the banner, `{}` without a
//!    check, nothing from a rerun block that is not an object, and nothing
//!    from a check older than the run's `run_manifest.tsv`;
//! 2. the ClinVar lock block of `scripts/run-all.sh` and of step 06 waits
//!    while a refresh holds the lock and goes on once it ends (needs
//!    `flock`; CI fails without it). Its control: the block with its `flock`
//!    call taken out must fail the check.

#[path = "port_live_server_util.rs"]
mod util;

use og_live::recheck;
use pycompat::json::{Value, loads};
use std::fs;
use std::path::Path;
use std::process::{Command, Stdio};
use std::thread;
use std::time::{Duration, Instant};
use util::*;

fn steps(v: &Value) -> Value {
    match v {
        Value::Object(o) => o.get("steps").cloned().unwrap_or(Value::Null),
        _ => Value::Null,
    }
}

fn banner(v: &Value) -> String {
    match v {
        Value::Object(o) => match o.get("banner") {
            Some(Value::Str(b)) => b.clone(),
            _ => String::new(),
        },
        _ => String::new(),
    }
}

fn empty() -> Value {
    Value::Object(indexmap::IndexMap::new())
}

#[test]
fn each_step_gets_its_badges_and_the_banner() {
    let w = Work::new("recheck-badges");
    let run = w.dir.join("S1");
    fs::create_dir_all(&run).unwrap();
    assert_eq!(recheck::state(&run), empty());
    write(
        &run.join("recheck/state.json"),
        r#"{"checked_utc": "2026-10-06T00:00:00Z", "items": [
            {"item": "clinvar_file_date", "label": "ClinVar (step 06's file)", "class": "refresh-available", "steps": ["06", "24"]},
            {"item": "VEP_CACHE_RELEASE", "label": "Ensembl / VEP cache", "class": "no-newer-release", "steps": ["13"]}],
            "rerun": {"finished_utc": "2026-10-05T10:00:00Z", "predicted_vs_redone": {"redone": ["06"]}}}"#,
    );
    let got = recheck::state(&run);
    assert_eq!(
        steps(&got),
        loads(r#"{"06": ["newer ClinVar available", "redone 2026-10-05"], "24": ["newer ClinVar available"]}"#).unwrap()
    );
    assert_eq!(
        banner(&got),
        "Database re-check of 2026-10-06T00:00:00Z: newer ClinVar available"
    );
    let Value::Object(o) = &got else { panic!() };
    assert_eq!(o["file"], Value::Str("recheck/state.json".into()));
    assert_eq!(o["checked_utc"], Value::Str("2026-10-06T00:00:00Z".into()));
    assert!(matches!(o["mtime"], Value::Int(n) if n > 0));
}

#[test]
fn a_rerun_block_that_is_not_an_object_gives_no_badge() {
    let w = Work::new("recheck-odd");
    let run = w.dir.join("S1");
    for rr in [
        "[1]",
        r#"{"finished_utc": "2026-10-05T10:00:00Z", "predicted_vs_redone": ["06"]}"#,
    ] {
        write(
            &run.join("recheck/state.json"),
            &format!(r#"{{"checked_utc": "2026-10-06T00:00:00Z", "items": [], "rerun": {rr}}}"#),
        );
        let got = recheck::state(&run);
        assert_eq!(steps(&got), empty(), "{rr}");
        assert!(banner(&got).contains("nothing newer"), "{rr}");
    }
    // Other classes are named in the banner, with no badge.
    write(
        &run.join("recheck/state.json"),
        r#"{"checked_utc": "2026-10-06", "items": [{"item": "PGS", "class": "error", "steps": ["30"]}, 5]}"#,
    );
    let got = recheck::state(&run);
    assert_eq!(
        (steps(&got), banner(&got)),
        (
            empty(),
            "Database re-check of 2026-10-06: PGS: error".to_owned()
        )
    );
    write(&run.join("recheck/state.json"), "not json");
    assert_eq!(recheck::state(&run), empty());
}

#[test]
fn a_check_older_than_the_runs_manifest_is_left_out() {
    let w = Work::new("recheck-stale");
    let run = w.dir.join("S1");
    write(
        &run.join("recheck/state.json"),
        r#"{"checked_utc": "2026-10-06T00:00:00Z", "items": [
            {"item": "clinvar_file_date", "label": "ClinVar", "class": "refresh-available", "steps": ["06"]}]}"#,
    );
    write(
        &run.join("run_manifest.tsv"),
        "run\twritten_utc\t2026-10-07 09:00:00 UTC\r\n",
    );
    assert_eq!(recheck::state(&run), empty());
    write(
        &run.join("run_manifest.tsv"),
        "run\twritten_utc\t2026-10-05 09:00:00 UTC\n",
    );
    assert_eq!(
        steps(&recheck::state(&run)),
        loads(r#"{"06": ["newer ClinVar available"]}"#).unwrap()
    );
}

/// The ClinVar lock block of a script: from `CVLOCK=` through its closing
/// `fi`.
fn lock_block(path: &Path) -> String {
    let text = fs::read_to_string(path).unwrap();
    let re = regex::Regex::new(r"(?ms)^CVLOCK=.*?\n  fi\nfi\n").unwrap();
    re.find(&text)
        .map(|m| m.as_str().to_owned())
        .unwrap_or_default()
}

/// `block`, run in bash while another process holds the lock exclusive: it
/// must wait, then go on once the lock is let go.
fn lock_waits(block: &str, name: &str) -> (bool, String) {
    let w = Work::new(name);
    let g = w.dir.join("genome");
    fs::create_dir_all(g.join("clinvar")).unwrap();
    let lock = g.join("clinvar/.refresh.lock");
    fs::write(&lock, "").unwrap();
    let mut holder = Command::new("flock")
        .arg("-x")
        .arg(&lock)
        .args(["sleep", "3"])
        .spawn()
        .unwrap();
    // Until the holder has it: a shared try fails while it is held.
    let end = Instant::now() + Duration::from_secs(5);
    while Instant::now() < end
        && Command::new("flock")
            .args(["-s", "-n"])
            .arg(&lock)
            .arg("true")
            .status()
            .is_ok_and(|s| s.success())
    {
        thread::sleep(Duration::from_millis(50));
    }
    let mut p = Command::new("bash")
        .args(["-c", &format!("set -euo pipefail\n{block}echo through\n")])
        .env("G", &g)
        .env("GENOME_DIR", &g)
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    thread::sleep(Duration::from_millis(1500));
    let waited = p.try_wait().unwrap().is_none();
    let _ = holder.wait();
    let out = p.wait_with_output().unwrap();
    let text = String::from_utf8_lossy(&out.stdout).into_owned();
    let ok = !block.is_empty()
        && waited
        && out.status.success()
        && text.contains("Waiting")
        && text.trim_end().ends_with("through");
    (ok, format!("waited {waited}, {:?}, {text}", out.status))
}

#[test]
fn the_clinvar_lock_block_waits_for_a_refresh_then_goes_on() {
    let have = Command::new("flock")
        .arg("-h")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok();
    if !have {
        let ci = std::env::var("CI").is_ok_and(|v| v == "true");
        assert!(
            !ci,
            "[FAIL] flock is not installed: CI must run the ClinVar lock checks"
        );
        println!("[SKIP] flock is not installed: the ClinVar lock block is not run here");
        return;
    }
    let flock_line = regex::Regex::new(r"(?m)^    flock -s -n 7 .*$").unwrap();
    for (i, script) in ["scripts/run-all.sh", "scripts/06-clinvar-screen.sh"]
        .iter()
        .enumerate()
    {
        let block = lock_block(&repo().join(script));
        assert!(!block.is_empty(), "no ClinVar lock block in {script}");
        let (ok, detail) = lock_waits(&block, &format!("cvlock-{i}"));
        assert!(ok, "{script}: {detail}");
        let planted = flock_line.replace(&block, "    :").into_owned();
        assert_ne!(planted, block, "the control's line is not in {script}");
        let (still, _) = lock_waits(&planted, &format!("cvlock-planted-{i}"));
        assert!(
            !still,
            "control: {script}'s block with no flock call passed the check"
        );
    }
}
