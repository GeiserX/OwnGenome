//! `og run` against the Python it replaces, case by case, on the same
//! files: `bin/run_state.py`, `bin/run_eta.py`, `bin/run_headroom.py` and
//! `bin/run_delta.py` still sit in the checkout (the live server and the
//! summary import them). Outputs must agree byte for byte, except
//! `generated_utc` (the clock), the disk numbers of `headroom` (free space
//! moves between two calls) and the host memory line where Python can read
//! the RAM and og cannot (not Linux).
//!
//! The Python is `OG_PYTHON`, else `python3`; it must be 3.9 or newer. A
//! missing Python fails the test: a parity test that skips proves nothing.

use jiff::tz::TimeZone;
use og_runlog::cli::{Env, run};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const TZ: &str = "Europe/Madrid";
/// 2026-01-03 12:00:00 UTC: after every fixture time.
const NOW: f64 = 1_767_441_600.0;

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .canonicalize()
        .unwrap()
}

fn steps() -> String {
    repo()
        .join("assets/steps.tsv")
        .to_string_lossy()
        .into_owned()
}

fn tmp(name: &str) -> PathBuf {
    static N: AtomicUsize = AtomicUsize::new(0);
    let d = std::env::temp_dir().join(format!(
        "og-runlog-parity-{}-{name}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d.canonicalize().unwrap()
}

fn copy_dir(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for e in fs::read_dir(from).unwrap() {
        let e = e.unwrap();
        let t = to.join(e.file_name());
        if e.file_type().unwrap().is_dir() {
            copy_dir(&e.path(), &t);
        } else {
            fs::copy(e.path(), t).unwrap();
        }
    }
}

struct Out {
    code: i32,
    stdout: String,
    stderr: String,
}

/// The Python tool `bin/<tool>.py` with `args`, in the zone [`TZ`].
fn py(tool: &str, args: &[&str]) -> Out {
    let exe = std::env::var("OG_PYTHON").unwrap_or_else(|_| "python3".into());
    let o = Command::new(&exe)
        .arg(repo().join("bin").join(format!("{tool}.py")))
        .args(args)
        .env("TZ", TZ)
        .env_remove("PGP_BOX")
        .env_remove("OWNGENOME_DISK_WARN_PERCENT")
        .env_remove("OWNGENOME_MEM_WARN_PERCENT")
        .output()
        .unwrap_or_else(|e| panic!("cannot run the Python oracle {exe}: {e}; set OG_PYTHON"));
    Out {
        code: o.status.code().unwrap_or(-1),
        stdout: String::from_utf8(o.stdout).unwrap(),
        stderr: String::from_utf8(o.stderr).unwrap(),
    }
}

/// Python code run with `bin/` importable, in the zone [`TZ`].
fn py_code(code: &str, args: &[&str]) -> Out {
    let exe = std::env::var("OG_PYTHON").unwrap_or_else(|_| "python3".into());
    let o = Command::new(&exe)
        .arg("-c")
        .arg(format!(
            "import sys; sys.path.insert(0, {:?})\n{code}",
            repo().join("bin").to_string_lossy()
        ))
        .args(args)
        .env("TZ", TZ)
        .env_remove("PGP_BOX")
        .output()
        .unwrap();
    Out {
        code: o.status.code().unwrap_or(-1),
        stdout: String::from_utf8(o.stdout).unwrap(),
        stderr: String::from_utf8(o.stderr).unwrap(),
    }
}

/// `og run <args>` in the zone [`TZ`] at [`NOW`].
fn og(args: &[&str]) -> Out {
    let env = Env {
        tz: TimeZone::get(TZ).expect("the system has zone data"),
        now: NOW,
    };
    let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
    let (mut o, mut e) = (Vec::new(), Vec::new());
    let code = run(&args, &env, &mut o, &mut e);
    Out {
        code: i32::from(code),
        stdout: String::from_utf8(o).unwrap(),
        stderr: String::from_utf8(e).unwrap(),
    }
}

fn mask_generated(text: &str) -> String {
    text.lines()
        .map(|l| {
            if l.trim_start().starts_with("\"generated_utc\"") {
                " \"generated_utc\": MASKED,"
            } else {
                l
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
}

fn same(what: &str, oracle: &str, new: &str) {
    if oracle != new {
        let at = oracle.lines().zip(new.lines()).position(|(a, b)| a != b);
        panic!(
            "{what}: og differs from the Python{}\n--- Python\n{oracle}\n--- og\n{new}",
            at.map_or_else(String::new, |i| format!(" at line {}", i + 1))
        );
    }
}

// ---- og run state ------------------------------------------------------

/// A sample folder with run-all.sh rows of each kind, a trace under
/// ../pipeline_info with this sample's and another's tasks, weblog
/// messages, and a summary section an earlier run wrote.
fn state_case(dir: &Path) -> (PathBuf, PathBuf) {
    let s = dir.join("S1");
    fs::create_dir_all(s.join("logs")).unwrap();
    fs::create_dir_all(dir.join("pipeline_info")).unwrap();
    fs::create_dir_all(s.join("clinvar")).unwrap();
    fs::write(s.join("clinvar/S1_hits.tsv"), "x\n").unwrap();
    // The section file is from before the run started (2026-01-02 09:00 UTC).
    let old = std::time::UNIX_EPOCH + std::time::Duration::from_secs(1_767_000_000);
    fs::File::options()
        .write(true)
        .open(s.join("clinvar/S1_hits.tsv"))
        .unwrap()
        .set_modified(old)
        .unwrap();
    fs::write(
        s.join("summary.json"),
        r#"{"sections": {"clinvar": {"state": "ok", "source": "clinvar/S1_hits.tsv"}, "cpsr": {"source": "gone.tsv"}, "x": 3}}"#,
    )
    .unwrap();
    fs::write(
        s.join("logs/run_status.tsv"),
        "# run-all.sh\nmeta\tstarted_epoch\t1767344400\nmeta\tdeclared_sex\tfemale\n\
         step\t01b\tplanned\nstep\t02\tplanned\nstep\t03\tplanned\nstep\t06\tplanned\nstep\t07\tok\n\
         step\t14\tskipped (no BAM)\nstep\t16\tfailed (rc=3)\nstep\t17\tnot finished\nstep\t18\tskipped\nstep\t19\tweird\n",
    )
    .unwrap();
    let head = "task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\trealtime\n";
    fs::write(
        dir.join("pipeline_info/trace_2026-01-02_10-00-00.txt"),
        format!(
            "{head}1\taa/01\t1\tW:FASTP (S1)\tCOMPLETED\t0\t2026-01-02 10:00:00.000\t20m 10s\t20m\n\
             2\taa/02\t2\tW:MINIMAP2_INDEX (S1)\tCOMPLETED\t0\t2026-01-02 10:00:00.000\t2m\t2m\n\
             3\taa/03\t3\tW:ALIGN_MINIMAP2 (S1)\tFAILED\t137\t2026-01-02 10:21:00.000\t1h\t1h\n\
             4\taa/04\t4\tW:FASTP (S2)\tCOMPLETED\t0\t2026-01-02 10:00:00.000\t1m\t1m\n\
             \n5\taa/05\t5\tbroken name here\tCOMPLETED\t0\t-\t-\t-\n"
        ),
    )
    .unwrap();
    let events = dir.join("events.jsonl");
    fs::write(
        &events,
        concat!(
            r#"{"event":"started","utcTime":"2026-01-02T09:00:00Z"}"#, "\n",
            r#"{"event":"process_submitted","trace":{"task_id":10,"name":"W:ALIGN_MINIMAP2 (S1)","status":"SUBMITTED","submit":1767349000000}}"#, "\n",
            r#"{"event":"process_started","trace":{"task_id":10,"name":"W:ALIGN_MINIMAP2 (S1)","status":"RUNNING","start":1767349005000}}"#, "\n",
            r#"{"event":"process_completed","trace":{"task_id":11,"process":"W:DEEPVARIANT","tag":"S1","exit":0,"start":1767349100000,"complete":1767352000000}}"#, "\n",
            r#"{"event":"process_completed","trace":{"task_id":12,"name":"W:VCF_PRECHECK (S1:chr1)","status":"FAILED","exit":1}}"#, "\n",
            "not json\n",
            r#"{"event":"process_started","trace":{"task_id":1,"name":"W:FASTP (S1)","status":"RUNNING"}}"#, "\n",
        ),
    )
    .unwrap();
    (s, events)
}

fn state_parity(sample: &str, dir: &Path, extra: &[&str]) {
    let st = steps();
    let py_out = dir.join("py.json");
    let og_out = dir.join("og.json");
    let d = dir.to_string_lossy().into_owned();
    let mut base = vec!["--sample", sample, "--sample-dir", &d, "--steps", &st];
    base.extend_from_slice(extra);
    let po = py_out.to_string_lossy().into_owned();
    let oo = og_out.to_string_lossy().into_owned();
    let mut pa = base.clone();
    pa.extend(["--out", &po]);
    let mut oa = vec!["state"];
    oa.extend(base);
    oa.extend(["--out", &oo]);
    let (p, o) = (py("run_state", &pa), og(&oa));
    assert_eq!((p.code, o.code), (0, 0), "{}{}", p.stderr, o.stderr);
    same(
        "stdout",
        &p.stdout.replace(&po, "OUT"),
        &o.stdout.replace(&oo, "OUT"),
    );
    same(
        "run_state.json",
        &mask_generated(&fs::read_to_string(&py_out).unwrap()),
        &mask_generated(&fs::read_to_string(&og_out).unwrap()),
    );
}

#[test]
fn state_a_run_going() {
    let dir = tmp("state");
    let (s, events) = state_case(&dir);
    state_parity("S1", &s, &["--events", &events.to_string_lossy()]);
    // Without the messages, from the trace alone.
    state_parity("S1", &s, &[]);
}

#[test]
fn state_a_run_that_ended() {
    let dir = tmp("state-end");
    let (s, _) = state_case(&dir);
    fs::create_dir_all(s.join("nextflow")).unwrap();
    fs::write(
        s.join("nextflow/process_status.tsv"),
        "process\tsucceeded\tcached\tfailed\tignored\tunfinished\n\
         W:FASTP\t1\t0\t0\t0\t0\nW:ALIGN_MINIMAP2\t0\t0\t0\t0\t1\nW:CRAM_TO_BAM\t0\t0\t0\t0\t0\nHAPLOCHECK\t0\t0\t0\t0\t0\nbad\trow\n",
    )
    .unwrap();
    let mut rs = fs::read_to_string(s.join("logs/run_status.tsv")).unwrap();
    rs.push_str("meta\tfinished_utc\t2026-01-03 05:12:00 UTC\n");
    fs::write(s.join("logs/run_status.tsv"), rs).unwrap();
    state_parity("S1", &s, &[]);
    // A run with no run-all.sh file: the process table alone says it ended.
    fs::remove_file(s.join("logs/run_status.tsv")).unwrap();
    state_parity("S1", &s, &[]);
}

#[test]
fn state_not_started_and_the_live_fixture() {
    let dir = tmp("state-empty");
    fs::create_dir_all(dir.join("S9")).unwrap();
    state_parity("S9", &dir.join("S9"), &[]);
    let live = tmp("state-live");
    copy_dir(
        &repo().join("tests/fixtures/live/DEMO-001"),
        &live.join("DEMO-001"),
    );
    let ev = repo().join("tests/fixtures/live/demo-001.weblog.jsonl");
    state_parity(
        "DEMO-001",
        &live.join("DEMO-001"),
        &["--events", &ev.to_string_lossy()],
    );
    let ops = tmp("state-ops");
    fs::create_dir_all(ops.join("OPS-001/logs")).unwrap();
    fs::copy(
        repo().join("tests/fixtures/live/ops/run_status.tsv"),
        ops.join("OPS-001/logs/run_status.tsv"),
    )
    .unwrap();
    let ev = repo().join("tests/fixtures/live/ops/events.jsonl");
    state_parity(
        "OPS-001",
        &ops.join("OPS-001"),
        &["--events", &ev.to_string_lossy()],
    );
}

// ---- og run eta --------------------------------------------------------

/// The Python estimate at [`NOW`], printed as `run_eta.py` prints it.
fn py_eta(outdir: &Path, sample: Option<&str>, cores: usize, box_: &str) -> Out {
    py_code(
        "import json, run_eta\n\
         a = sys.argv[1:]\n\
         est = run_eta.estimate(a[0], a[1] or None, now=float(a[2]), cores=int(a[3]), box=a[4])\n\
         sys.stdout.write(json.dumps(est, indent=1) + '\\n')",
        &[
            &outdir.to_string_lossy(),
            sample.unwrap_or(""),
            &NOW.to_string(),
            &cores.to_string(),
            box_,
        ],
    )
}

fn og_eta(outdir: &Path, sample: Option<&str>, cores: usize, box_: &str) -> Out {
    let (o, c, st) = (
        outdir.to_string_lossy().into_owned(),
        cores.to_string(),
        steps(),
    );
    let mut a = vec!["eta", &o, "--cores", &c, "--box", box_, "--steps", &st];
    if let Some(s) = sample {
        a.extend(["--sample", s]);
    }
    og(&a)
}

fn eta_outdir(name: &str) -> PathBuf {
    let d = tmp(name);
    fs::create_dir_all(d.join("pipeline_info")).unwrap();
    for t in ["trace_new.txt", "trace_old.txt"] {
        fs::copy(
            repo().join("tests/fixtures/run_eta").join(t),
            d.join("pipeline_info").join(format!(
                "trace_2026-01-0{}.txt",
                if t == "trace_new.txt" { 2 } else { 1 }
            )),
        )
        .unwrap();
    }
    d
}

fn history(d: &Path, box_: &str) -> PathBuf {
    d.join("pipeline_info/step_history")
        .join(format!("{box_}.tsv"))
}

#[test]
fn eta_history_and_ranges() {
    for cores in [8, 32] {
        let (a, b) = (eta_outdir("eta-py"), eta_outdir("eta-og"));
        let (p, o) = (
            py_eta(&a, None, cores, "parity"),
            og_eta(&b, None, cores, "parity"),
        );
        assert_eq!((p.code, o.code), (0, 0), "{}{}", p.stderr, o.stderr);
        same("eta JSON", &p.stdout, &o.stdout);
        same(
            "step history",
            &fs::read_to_string(history(&a, "parity")).unwrap(),
            &fs::read_to_string(history(&b, "parity")).unwrap(),
        );
        // A second call on the same traces leaves the history as it is, so
        // og on a history Python wrote does not rewrite it.
        let h = history(&a, "parity");
        let before = fs::metadata(&h).unwrap();
        std::thread::sleep(std::time::Duration::from_millis(20));
        let again = og_eta(&a, None, cores, "parity");
        assert_eq!(again.code, 0);
        let after = fs::metadata(&h).unwrap();
        use std::os::unix::fs::MetadataExt;
        assert_eq!(
            (before.ino(), before.mtime_nsec(), before.mtime()),
            (after.ino(), after.mtime_nsec(), after.mtime()),
            "og rewrote a history whose text had not changed"
        );
    }
}

#[test]
fn eta_hours_left_of_a_run_going() {
    let mk = |name: &str| {
        let d = eta_outdir(name);
        // A box history of earlier runs, then a run of demo1 still going.
        fs::create_dir_all(d.join("pipeline_info/step_history")).unwrap();
        fs::copy(
            repo().join("tests/fixtures/live/ops/step_history.tsv"),
            history(&d, "opsbox"),
        )
        .unwrap();
        fs::create_dir_all(d.join("demo1/logs")).unwrap();
        fs::write(
            d.join("demo1/logs/run_status.tsv"),
            "meta\tstarted_epoch\t1767340000\nstep\t01b\tok\nstep\t02\tplanned\nstep\t03\tplanned\nstep\t06\tplanned\nstep\t11\tplanned\nstep\t25\tplanned\n",
        )
        .unwrap();
        // Another sample's run, ended: its traces count.
        fs::create_dir_all(d.join("other/logs")).unwrap();
        fs::write(
            d.join("other/logs/run_status.tsv"),
            "meta\tstarted_epoch\t1\nmeta\tfinished_utc\tx\n",
        )
        .unwrap();
        d
    };
    for cores in [4, 16] {
        let (a, b) = (mk("eta-run-py"), mk("eta-run-og"));
        let (p, o) = (
            py_eta(&a, Some("demo1"), cores, "opsbox"),
            og_eta(&b, Some("demo1"), cores, "opsbox"),
        );
        assert_eq!((p.code, o.code), (0, 0), "{}{}", p.stderr, o.stderr);
        same("eta JSON with hours left", &p.stdout, &o.stdout);
        same(
            "step history",
            &fs::read_to_string(history(&a, "opsbox")).unwrap(),
            &fs::read_to_string(history(&b, "opsbox")).unwrap(),
        );
    }
}

// ---- og run headroom ---------------------------------------------------

/// Lines with the numbers of each disk line and, off Linux, the host
/// memory line masked.
fn mask_host(text: &str) -> String {
    text.lines()
        .map(|l| {
            if l.starts_with("disk ") {
                l.chars()
                    .map(|c| if c.is_ascii_digit() { '#' } else { c })
                    .collect()
            } else if l.starts_with("memory host:") && !cfg!(target_os = "linux") {
                "memory host: MASKED".to_owned()
            } else {
                l.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
}

fn headroom_parity(args: &[&str]) {
    let st = steps();
    let mut oa = vec!["headroom", "--steps", &st];
    oa.extend_from_slice(args);
    let (p, o) = (py("run_headroom", args), og(&oa));
    assert_eq!(p.code, o.code, "{}{}", p.stderr, o.stderr);
    if p.code == 0 {
        same(
            "headroom lines",
            &mask_host(&p.stdout),
            &mask_host(&o.stdout),
        );
    } else {
        let last = |e: &str| {
            e.lines()
                .last()
                .unwrap_or("")
                .split_once("error: ")
                .map(|x| x.1.to_owned())
        };
        assert_eq!(last(&p.stderr), last(&o.stderr), "{}{}", p.stderr, o.stderr);
    }
}

#[test]
fn headroom_tasks_hosts_and_disks() {
    let fx = repo().join("tests/fixtures/run_headroom");
    let (tr, tr0, ev) = (
        fx.join("trace.txt").to_string_lossy().into_owned(),
        fx.join("trace_no_memory.txt")
            .to_string_lossy()
            .into_owned(),
        fx.join("events.jsonl").to_string_lossy().into_owned(),
    );
    headroom_parity(&["--trace", &tr, "--trace", &tr0]);
    headroom_parity(&["--events", &ev, "--max-memory", "64.GB"]);
    headroom_parity(&[
        "--events",
        &ev,
        "--trace",
        &tr,
        "--mem-warn-percent",
        "95",
        "--max-memory",
        "1 TB",
    ]);
    let d = tmp("headroom");
    let (w, r) = (
        d.join("work").to_string_lossy().into_owned(),
        d.join("not/yet").to_string_lossy().into_owned(),
    );
    fs::create_dir_all(&w).unwrap();
    headroom_parity(&["--work-dir", &w, "--outdir", &r, "--disk-warn-percent", "1"]);
    headroom_parity(&["--work-dir", &w, "--disk-warn-percent", "100"]);
    // The weblog metadata names the folders: a relative outdir is under
    // the launch folder.
    let ops = d.join("ops.jsonl");
    fs::write(
        &ops,
        fs::read_to_string(repo().join("tests/fixtures/live/ops/events.jsonl"))
            .unwrap()
            .replace("@WORK@", &w)
            .replace("\"@OUT@\"", "\"results\"")
            .replacen(
                "\"workflow\":{",
                &format!("\"workflow\":{{\"launchDir\":\"{w}\","),
                1,
            ),
    )
    .unwrap();
    headroom_parity(&["--events", &ops.to_string_lossy()]);
    headroom_parity(&[]);
    for bad in ["0", "101", "x", "nan"] {
        headroom_parity(&["--disk-warn-percent", bad]);
        headroom_parity(&["--mem-warn-percent", bad]);
    }
}

// ---- og run delta ------------------------------------------------------

fn summary(clinvar: &str, prs: &str, extra: &str) -> String {
    format!(
        r#"{{"generated_utc": "2026-01-0{extra}T10:00:00Z", "sections": {{
  "clinvar": {{"state": "ok", "title": "ClinVar", "data": {{"hits": [{clinvar}]}}}},
  "cpsr": {{"state": "ok", "data": {{"classification": {{"Pathogenic": 1, "VUS": {extra}}}}}}},
  "cpic": {{"state": "stale", "data": {{"genes": [{{"gene": "CYP2C19", "diplotype": "*1/*{extra}", "phenotype": "x"}}, {{"gene": null, "diplotype": "?"}}]}}}},
  "pypgx": {{"state": "ok", "data": {{"cyp2d6": {{"call": "*1/*{extra}", "n": [1, 2.5, null, true]}}}}}},
  "prs": {{"state": "ok", "data": {{"scores": [{prs}]}}}},
  "roh": {{"state": "{roh}", "title": "Runs <of> homozygosity"}}
}}, "cyp2d6": {{"consensus": {{"result": {cons}}}}}}}"#,
        roh = if extra == "1" { "ok" } else { "missing" },
        cons = if extra == "1" { "null" } else { "\"*1/*2\"" },
    )
}

#[test]
fn delta_every_kind_of_change() {
    let d = tmp("delta");
    let s = d.join("S1");
    let prev = s.join("history/2026-01-01T100000Z");
    fs::create_dir_all(&prev).unwrap();
    fs::create_dir_all(s.join("recheck")).unwrap();
    fs::write(
        prev.join("summary.json"),
        summary(
            r#"{"chrom": "chr1", "pos": 1000, "ref": "A", "alt": "G", "gene": "GENEA", "significance": "Pathogenic", "stars": 2, "review_status": "a", "genotype": "0/1"},
               {"chrom": "chr2", "pos": 2000, "ref": "C", "alt": "T", "gene": "GENEB", "significance": "Likely pathogenic", "stars": 1},
               {"chrom": "chr3", "pos": 3, "ref": "C", "alt": "T"}"#,
            r#"{"pgs_id": "PGS1", "condition": "One", "percentile": 0}, {"pgs_id": "PGS2", "condition": "Two", "percentile": "50.0"},
               {"pgs_id": "PGS3", "condition": "Three", "percentile": 10}, {"pgs_id": "PGS4", "condition": "Four"}"#,
            "1",
        ),
    )
    .unwrap();
    fs::write(
        prev.join("run_manifest.tsv"),
        "# manifest\nrun\twritten_utc\t2026-01-01 10:00:00 UTC\nrun\tgit_commit\tabc\nversion\tPYTHON_IMAGE\tpython:3.15\n\
         image\tpython:3.15\tpython@sha256:0000\nversion\tCLINVAR_DATE\t20251201\ndata\tclinvar_file_date\t20251201\tx\n",
    )
    .unwrap();
    let cur = s.join("summary.json");
    fs::write(
        &cur,
        summary(
            r#"{"chrom": "chr1", "pos": 1000, "ref": "A", "alt": "G", "gene": "GENEA", "significance": "Pathogenic", "stars": 3, "review_status": "b", "genotype": "1/1"},
               {"chrom": "chr2", "pos": 2000, "ref": "C", "alt": "T", "gene": "GENEB", "significance": "Uncertain significance", "stars": 1},
               {"chrom": "chr4", "pos": 4, "ref": "G", "alt": "A", "gene": "GENÉ"}"#,
            r#"{"pgs_id": "PGS1", "condition": "One", "percentile": 7.5}, {"pgs_id": "PGS2", "condition": "Two", "percentile": "53.0"},
               {"pgs_id": "PGS3", "condition": "Three"}, {"pgs_id": "PGS5", "condition": "Five", "percentile": 1}"#,
            "2",
        ),
    )
    .unwrap();
    let man = s.join("run_manifest.tsv");
    fs::write(
        &man,
        "run\twritten_utc\t2026-01-02 10:00:00 UTC\nrun\tgit_commit\tdef\nversion\tPYTHON_IMAGE\tpython:3.15\n\
         image\tpython:3.15\tpython@sha256:1111\nversion\tCLINVAR_DATE\t20260101\ndata\tclinvar_file_date\t20260101\n",
    )
    .unwrap();
    fs::write(
        s.join("recheck/state.json"),
        r#"{"checked_utc": "2026-01-03T00:00:00Z", "items": [{"item": "clinvar", "label": "ClinVar", "class": "rerun-ready", "action": "rerun"}, {"item": "vep", "class": "no-newer-release"}, 3], "rerun": {"rc": 0}}"#,
    )
    .unwrap();
    let (c, m, p) = (
        cur.to_string_lossy().into_owned(),
        man.to_string_lossy().into_owned(),
        prev.to_string_lossy().into_owned(),
    );
    delta_parity(
        &d,
        &[
            "--sample",
            "S1",
            "--summary",
            &c,
            "--manifest",
            &m,
            "--previous",
            &p,
        ],
    );
    // The same run again: no change. No previous run. An old re-check.
    delta_parity(
        &d,
        &[
            "--sample",
            "S1",
            "--summary",
            &c,
            "--manifest",
            &m,
            "--previous",
            &format!("{p}/"),
        ],
    );
    fs::copy(&cur, prev.join("summary.json")).unwrap();
    fs::copy(&man, prev.join("run_manifest.tsv")).unwrap();
    delta_parity(
        &d,
        &[
            "--sample",
            "S1",
            "--summary",
            &c,
            "--manifest",
            &m,
            "--previous",
            &p,
        ],
    );
    delta_parity(&d, &["--sample", "S1", "--summary", &c, "--manifest", &m]);
    fs::write(&man, "run\twritten_utc\t2026-01-04 10:00:00 UTC\n").unwrap();
    delta_parity(&d, &["--sample", "S1", "--summary", &c]);
    // A summary that cannot be read exits 2.
    let bad = d.join("bad.json");
    fs::write(&bad, "{").unwrap();
    let out = d.join("x").to_string_lossy().into_owned();
    let (pb, ob) = (
        py(
            "run_delta",
            &[
                "--sample",
                "S1",
                "--summary",
                &bad.to_string_lossy(),
                "--out",
                &out,
            ],
        ),
        og(&[
            "delta",
            "--sample",
            "S1",
            "--summary",
            &bad.to_string_lossy(),
            "--out",
            &out,
        ]),
    );
    assert_eq!((pb.code, ob.code), (2, 2), "{}{}", pb.stderr, ob.stderr);
    assert!(
        ob.stderr.starts_with("ERROR: cannot read a summary: "),
        "{}",
        ob.stderr
    );
}

fn delta_parity(d: &Path, args: &[&str]) {
    let (po, oo) = (
        d.join("py_delta").to_string_lossy().into_owned(),
        d.join("og_delta").to_string_lossy().into_owned(),
    );
    let mut pa = args.to_vec();
    pa.extend(["--out", &po]);
    let mut oa = vec!["delta"];
    oa.extend_from_slice(args);
    oa.extend(["--out", &oo]);
    let (p, o) = (py("run_delta", &pa), og(&oa));
    assert_eq!((p.code, o.code), (0, 0), "{}{}", p.stderr, o.stderr);
    same("delta stdout", &p.stdout, &o.stdout);
    for ext in ["txt", "html", "json"] {
        same(
            &format!("delta .{ext}"),
            &fs::read_to_string(format!("{po}.{ext}")).unwrap(),
            &fs::read_to_string(format!("{oo}.{ext}")).unwrap(),
        );
    }
}

/// The history lock is the one `run_eta.py` and `flock(1)` take: while
/// another process holds `<history>.lock`, og waits for it.
#[test]
fn eta_waits_for_the_history_lock() {
    use std::io::{BufRead, BufReader};
    let d = eta_outdir("eta-lock");
    let lock = d.join("pipeline_info/step_history/held.tsv.lock");
    fs::create_dir_all(lock.parent().unwrap()).unwrap();
    let exe = std::env::var("OG_PYTHON").unwrap_or_else(|_| "python3".into());
    let mut holder = Command::new(&exe)
        .args([
            "-c",
            "import fcntl, sys, time\nf = open(sys.argv[1], 'a')\nfcntl.flock(f, fcntl.LOCK_EX)\nprint('held', flush=True)\ntime.sleep(1.5)",
        ])
        .arg(&lock)
        .stdout(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    let mut line = String::new();
    BufReader::new(holder.stdout.take().unwrap())
        .read_line(&mut line)
        .unwrap();
    assert_eq!(line, "held\n");
    let t = std::time::Instant::now();
    let o = og_eta(&d, None, 8, "held");
    let waited = t.elapsed();
    holder.wait().unwrap();
    assert_eq!(o.code, 0, "{}", o.stderr);
    assert!(
        waited >= std::time::Duration::from_millis(1000),
        "og did not wait for the lock: {waited:?}"
    );
    assert!(history(&d, "held").is_file());
}
