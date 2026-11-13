//! `og recheck`.
//!
//! ```text
//! og recheck --sample S [--genome-dir DIR] [--offline FILE] [--owngenome-repo OWNER/NAME] [--rerun [--refresh-clinvar]]
//! og recheck --affected KEY[,KEY...]      the steps a change of these manifest keys redoes
//! og recheck --launch-dir DIR             the last bare `nextflow run` there, with -resume
//! ```
//!
//! Exit codes: 0, 1 when a lookup failed, 2 on a usage error or a sample
//! with no manifest, 3 when a lock is held or a run is going on, else the
//! replay's own exit code.

use crate::answers::{self, Answers, Memo, Offline};
use crate::files::{last_rerun, log_check, write_check};
use crate::rerun::{Rerun, launch_dir_command, refresh_clinvar, rerun};
use crate::{DEFAULT_REPO, Here, RAW, affected_steps, check, print_check};
use og_fresh::pins::{Pins, hla_db_release};
use og_runlog::steps::read_steps;
use std::io::Write;
use std::path::{Path, PathBuf};

pub const USAGE: &str =
    "usage: og recheck --sample S [--genome-dir DIR] [--offline FILE] [--owngenome-repo OWNER/NAME]
                  [--rerun [--refresh-clinvar]]
       og recheck --affected KEY[,KEY...]
       og recheck --launch-dir DIR
";

const HELP: &str =
    "Say which of a finished sample's results a newer ClinVar or catalog would change.

  --sample S              the sample folder under GENOME_DIR
  --genome-dir DIR        default: $GENOME_DIR
  --offline FILE          read the upstream answers from this JSON file
  --owngenome-repo O/N    the pipeline's repository (default GeiserX/OwnGenome)
  --rerun                 replay the last run-all.sh command
  --refresh-clinvar       with --rerun: setup.sh --refresh clinvar first
  --affected KEYS         print the steps a change of these manifest keys redoes
  --launch-dir DIR        print the last nextflow run of this launch directory with -resume

docs/recheck.md has the classes, the files it writes and how to schedule it.
";

fn usage_error(err: &mut dyn Write, msg: &str) -> u8 {
    let _ = write!(err, "{USAGE}");
    let _ = writeln!(err, "og recheck: error: {msg}");
    2
}

fn code(rc: i32) -> u8 {
    u8::try_from(rc.clamp(0, 255)).unwrap_or(1)
}

/// `og recheck` on this machine.
pub fn main(args: &[String]) -> u8 {
    run(
        args,
        &mut std::io::stdout().lock(),
        &mut std::io::stderr().lock(),
    )
}

pub fn run(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let p = match og_fresh::cli::parse(
        args,
        &[
            "sample",
            "genome-dir",
            "offline",
            "owngenome-repo",
            "pgs-api",
            "raw",
            "affected",
            "launch-dir",
            "steps",
            "setup-script",
        ],
        &["rerun", "refresh-clinvar", "help"],
    ) {
        Ok(p) => p,
        Err(e) => return usage_error(err, &e),
    };
    let one = |k: &str| og_fresh::cli::one(&p, k);
    let flag = |k: &str| og_fresh::cli::flag(&p, k);
    if flag("help") {
        let _ = write!(out, "{USAGE}\n{HELP}");
        return 0;
    }
    let root = og_fresh::checkout_root();
    let steps_path = match (one("steps"), &root) {
        (Some(s), _) => PathBuf::from(s),
        (None, Some(r)) => r.join("assets/steps.tsv"),
        (None, None) => {
            let _ = writeln!(
                err,
                "ERROR: cannot find assets/steps.tsv above og or the working directory; pass --steps"
            );
            return 2;
        }
    };
    let steps = match read_steps(&steps_path) {
        Ok(s) => s,
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            return 1;
        }
    };
    if let Some(keys) = one("affected") {
        let keys: Vec<&str> = keys.split(',').filter(|k| !k.is_empty()).collect();
        let _ = writeln!(out, "{}", affected_steps(&steps, &keys).join(" "));
        return 0;
    }
    if let Some(d) = one("launch-dir") {
        let Some(cmd) = launch_dir_command(Path::new(d)) else {
            let _ = writeln!(err, "ERROR: no .nextflow/history in {d}");
            return 2;
        };
        let _ = writeln!(
            out,
            "The last run there, with -resume (not run; check it, then run it from that directory):"
        );
        let _ = writeln!(out, "{cmd}");
        return 0;
    }
    let genome = one("genome-dir")
        .map(str::to_owned)
        .or_else(|| std::env::var("GENOME_DIR").ok().filter(|g| !g.is_empty()));
    let (Some(sample), Some(genome)) = (one("sample").filter(|s| !s.is_empty()), genome) else {
        return usage_error(err, "give --sample and --genome-dir (or GENOME_DIR)");
    };
    if flag("refresh-clinvar") && !flag("rerun") {
        return usage_error(err, "--refresh-clinvar goes with --rerun");
    }
    let Some(root) = root else {
        let _ = writeln!(
            err,
            "ERROR: cannot find the OwnGenome checkout above og or the working directory"
        );
        return 2;
    };
    let genome = PathBuf::from(genome);
    let sample_dir = genome.join(sample);
    if !sample_dir.join("run_manifest.tsv").is_file() {
        let _ = writeln!(
            err,
            "ERROR: no {}/run_manifest.tsv: run the pipeline on {sample} first.",
            sample_dir.display()
        );
        return 2;
    }
    if flag("refresh-clinvar") {
        let setup =
            one("setup-script").map_or_else(|| root.join("scripts/setup.sh"), PathBuf::from);
        let _ = out.flush();
        match refresh_clinvar(&genome, &setup, &steps) {
            Ok((0, _)) => {}
            Ok((rc, why)) => {
                if let Some(w) = why {
                    let _ = writeln!(err, "ERROR: {w}");
                }
                return code(rc);
            }
            Err(e) => {
                let _ = writeln!(err, "ERROR: {e}");
                return 1;
            }
        }
    }
    let this_pins = match Pins::read(&root.join("versions.env")) {
        Ok(p) => p.values,
        Err(e) => {
            let _ = writeln!(err, "ERROR: cannot read versions.env: {e}");
            return 1;
        }
    };
    let this_hla = hla_db_release(&root.join("scripts/lib/common.sh")).ok();
    let here = Here {
        pins: &this_pins,
        hla: this_hla.as_deref(),
    };
    let inner: Box<dyn Answers> = match one("offline") {
        Some(f) => match Offline::read(Path::new(f)) {
            Ok(o) => Box::new(o),
            Err(e) => {
                let _ = writeln!(err, "ERROR: {e}");
                return 1;
            }
        },
        None => {
            let man = og_runlog::files::read_manifest(Some(&sample_dir.join("run_manifest.tsv")))
                .unwrap_or_default();
            let pins = Pins {
                values: man
                    .iter()
                    .filter(|((sec, _), _)| sec == "version")
                    .map(|((_, k), v)| (k.clone(), v.clone()))
                    .collect(),
                notes: Default::default(),
            };
            let raw = one("raw")
                .map(str::to_owned)
                .or_else(|| std::env::var("OWNGENOME_RECHECK_RAW").ok())
                .unwrap_or_else(|| RAW.into());
            Box::new(answers::online(
                pins,
                one("owngenome-repo").unwrap_or(DEFAULT_REPO),
                one("pgs-api").unwrap_or(answers::PGS_API),
                &raw,
                &root,
            ))
        }
    };
    let answers = Memo::new(inner);
    let d = sample_dir.join("recheck");
    let (state_path, log) = (d.join("state.json"), d.join("log.tsv"));
    let mut state = check(sample, &genome, &answers, &steps, &here);
    // The last rerun stays until another one replaces it.
    state.rerun = last_rerun(&state_path);
    let io = |e: std::io::Error, err: &mut dyn Write| {
        let _ = writeln!(err, "ERROR: cannot write {}: {e}", d.display());
        1
    };
    if let Err(e) = log_check(&log, &state).and_then(|()| write_check(&state_path, &state)) {
        return io(e, err);
    }
    let _ = print_check(&state, out);
    let _ = writeln!(out, "Written: {}", state_path.display());
    let mut errors = state.errors();
    if flag("rerun") {
        let _ = out.flush();
        let rc = match rerun(sample, &sample_dir, &steps, &mut state, &log) {
            Err(e) => return io(e, err),
            Ok(Rerun::Refused(rc, why)) => {
                let _ = writeln!(err, "ERROR: {why}");
                return rc;
            }
            Ok(Rerun::Ran(rc)) => rc,
        };
        // The replay wrote a new run_manifest.tsv: check again against it,
        // with the same answers, so state.json no longer calls ready what
        // the rerun has just used.
        let mut fresh = check(sample, &genome, &answers, &steps, &here);
        fresh.rerun = state.rerun.take();
        if let Err(e) = log_check(&log, &fresh).and_then(|()| write_check(&state_path, &fresh)) {
            return io(e, err);
        }
        if rc != 0 {
            return code(rc);
        }
        // What state.json now says: the replay's manifest may name a score
        // the first check did not ask about.
        errors = fresh.errors();
    }
    u8::from(errors > 0)
}
