//! `og run state|tasks|eta|headroom|delta`. Options keep the spelling of the
//! Python tools they replace; an abbreviated option is refused. A usage
//! error prints `usage:` and the error on stderr and exits 2.

use crate::eta::{EtaArgs, estimate, read_trace};
use crate::state::{RunFiles, Sources};
use crate::steps::{Step, default_path, read_steps};
use jiff::tz::TimeZone;
use pycompat::json::{Dump, Value, dumps};
use std::io::Write;
use std::path::{Path, PathBuf};

pub const USAGE: &str = "usage: og run state --sample S --sample-dir DIR [--trace FILE]... [--events FILE] [--steps FILE] [--out FILE]
       og run tasks --trace FILE [--trace FILE]... [--steps FILE]
       og run eta OUTDIR [--sample S] [--cores N] [--box NAME] [--out FILE] [--steps FILE]
       og run headroom [--events FILE] [--trace FILE]... [--work-dir DIR] [--outdir DIR] [--max-memory SIZE]
                       [--disk-warn-percent N] [--mem-warn-percent N] [--steps FILE]
       og run delta --sample S --summary FILE [--manifest FILE] [--previous DIR | --newest-history SAMPLE_DIR]
                    --out PREFIX [--recheck FILE]
";

/// The clock and zone a call runs with.
pub struct Env {
    /// The zone naive trace times are in.
    pub tz: TimeZone,
    /// Now, as epoch seconds.
    pub now: f64,
}

impl Env {
    /// This machine's zone (`TZ`, else `/etc/localtime`, else UTC) and clock.
    pub fn system() -> Env {
        let now = jiff::Timestamp::now();
        Env {
            tz: TimeZone::system(),
            now: now.as_second() as f64 + f64::from(now.subsec_microsecond()) / 1e6,
        }
    }
}

/// `og run ...` on this machine.
pub fn main(args: &[String]) -> u8 {
    run(
        args,
        &Env::system(),
        &mut std::io::stdout().lock(),
        &mut std::io::stderr().lock(),
    )
}

/// The options one verb takes: (name, repeatable), and its positionals.
/// A single option given twice keeps the last value, as argparse does.
struct Spec {
    verb: &'static str,
    opts: &'static [(&'static str, bool)],
    positionals: &'static [&'static str],
    required: &'static [&'static str],
}

/// A parsed command line.
#[derive(Default)]
struct Parsed {
    opts: Vec<(&'static str, String)>,
    positionals: Vec<String>,
}

impl Parsed {
    fn one(&self, k: &str) -> Option<&str> {
        self.opts
            .iter()
            .rev()
            .find(|(n, _)| *n == k)
            .map(|(_, v)| v.as_str())
    }
    fn all(&self, k: &str) -> Vec<&str> {
        self.opts
            .iter()
            .filter(|(n, _)| *n == k)
            .map(|(_, v)| v.as_str())
            .collect()
    }
}

fn parse(spec: &Spec, args: &[String]) -> Result<Parsed, String> {
    let mut p = Parsed::default();
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        if let Some(body) = a.strip_prefix("--").filter(|b| !b.is_empty()) {
            let (name, inline) = match body.split_once('=') {
                Some((n, v)) => (n, Some(v.to_owned())),
                None => (body, None),
            };
            let Some(&(opt, _)) = spec.opts.iter().find(|(o, _)| *o == name) else {
                return Err(format!("unrecognized arguments: {a}"));
            };
            let value = match inline {
                Some(v) => v,
                None => {
                    i += 1;
                    args.get(i)
                        .cloned()
                        .ok_or_else(|| format!("argument --{opt}: expected one argument"))?
                }
            };
            p.opts.push((opt, value));
        } else if p.positionals.len() < spec.positionals.len() {
            p.positionals.push(a.clone());
        } else {
            return Err(format!("unrecognized arguments: {a}"));
        }
        i += 1;
    }
    let missing: Vec<String> = spec
        .required
        .iter()
        .filter(|o| p.one(o).is_none())
        .map(|o| format!("--{o}"))
        .chain(
            spec.positionals[p.positionals.len()..]
                .iter()
                .map(|x| (*x).to_owned()),
        )
        .collect();
    if !missing.is_empty() {
        return Err(format!(
            "the following arguments are required: {}",
            missing.join(", ")
        ));
    }
    Ok(p)
}

const STATE: Spec = Spec {
    verb: "state",
    opts: &[
        ("sample", false),
        ("sample-dir", false),
        ("trace", true),
        ("events", false),
        ("steps", false),
        ("out", false),
    ],
    positionals: &[],
    required: &["sample", "sample-dir"],
};
const TASKS: Spec = Spec {
    verb: "tasks",
    opts: &[("trace", true), ("steps", false)],
    positionals: &[],
    required: &["trace"],
};
const ETA: Spec = Spec {
    verb: "eta",
    opts: &[
        ("sample", false),
        ("cores", false),
        ("box", false),
        ("out", false),
        ("steps", false),
    ],
    positionals: &["outdir"],
    required: &[],
};
const HEADROOM: Spec = Spec {
    verb: "headroom",
    opts: &[
        ("events", false),
        ("trace", true),
        ("work-dir", false),
        ("outdir", false),
        ("max-memory", false),
        ("disk-warn-percent", false),
        ("mem-warn-percent", false),
        ("steps", false),
    ],
    positionals: &[],
    required: &[],
};
const DELTA: Spec = Spec {
    verb: "delta",
    opts: &[
        ("sample", false),
        ("summary", false),
        ("manifest", false),
        ("previous", false),
        ("newest-history", false),
        ("out", false),
        ("recheck", false),
    ],
    positionals: &[],
    required: &["sample", "summary", "out"],
};

fn usage_error(err: &mut dyn Write, verb: &str, msg: &str) -> u8 {
    let _ = write!(err, "{USAGE}");
    let _ = writeln!(err, "og run {verb}: error: {msg}");
    2
}

fn fail(err: &mut dyn Write, msg: impl std::fmt::Display) -> u8 {
    let _ = writeln!(err, "ERROR: {msg}");
    1
}

fn steps_of(p: &Parsed) -> Result<Vec<Step>, String> {
    let path = match p.one("steps") {
        Some(s) => PathBuf::from(s),
        None => default_path()
            .ok_or("cannot find assets/steps.tsv above the og binary; pass --steps")?,
    };
    read_steps(&path).map_err(|e| e.to_string())
}

fn json_text(v: &Value) -> String {
    dumps(
        v,
        &Dump {
            indent: Some(1),
            ..Dump::default()
        },
    )
    .expect("allow_nan is on")
        + "\n"
}

fn write_atomic(path: &Path, text: &str) -> std::io::Result<()> {
    let mut f = og_io::AtomicFile::create(path)?;
    f.write_all(text.as_bytes())?;
    f.commit()
}

/// Runs one `og run` verb with `env`'s clock and zone.
pub fn run(args: &[String], env: &Env, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let Some(verb) = args.first() else {
        return usage_error(err, "", "the following arguments are required: verb");
    };
    let rest = &args[1..];
    if rest.iter().any(|a| a == "-h" || a == "--help") {
        let _ = write!(out, "{USAGE}");
        return 0;
    }
    let spec = match verb.as_str() {
        "state" => &STATE,
        "tasks" => &TASKS,
        "eta" => &ETA,
        "headroom" => &HEADROOM,
        "delta" => &DELTA,
        "-h" | "--help" => {
            let _ = write!(out, "{USAGE}");
            return 0;
        }
        other => return usage_error(err, "", &format!("unknown verb '{other}'")),
    };
    let p = match parse(spec, rest) {
        Ok(p) => p,
        Err(e) => return usage_error(err, spec.verb, &e),
    };
    match spec.verb {
        "state" => state(&p, env, out, err),
        "tasks" => tasks(&p, out, err),
        "eta" => eta(&p, env, out, err),
        "headroom" => headroom(&p, env, out, err),
        _ => delta(&p, out, err),
    }
}

fn state(p: &Parsed, env: &Env, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let steps = match steps_of(p) {
        Ok(s) => s,
        Err(e) => return fail(err, e),
    };
    let sample = p.one("sample").unwrap_or_default();
    let dir = PathBuf::from(p.one("sample-dir").unwrap_or_default());
    let traces: Vec<PathBuf> = p.all("trace").into_iter().map(PathBuf::from).collect();
    let events = p.one("events").map(Path::new);
    let files = match RunFiles::read(
        &dir,
        &Sources {
            traces: &traces,
            events,
            stale: false,
        },
    ) {
        Ok(f) => f,
        Err(e) => return fail(err, e),
    };
    let st = files.state(sample, &steps, env.now, &env.tz);
    let target = p
        .one("out")
        .map_or_else(|| dir.join("run_state.json"), PathBuf::from);
    if let Err(e) = write_atomic(&target, &json_text(&st.to_json())) {
        return fail(err, format!("cannot write {}: {e}", target.display()));
    }
    let _ = writeln!(out, "Run state written: {}", target.display());
    0
}

fn tasks(p: &Parsed, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let steps = match steps_of(p) {
        Ok(s) => s,
        Err(e) => return fail(err, e),
    };
    let proc_step = crate::steps::process_steps(&steps);
    let mut text = String::new();
    for t in p.all("trace") {
        let r = og_nf::for_each_row(Path::new(t), |r| {
            if r.is_blank() {
                return;
            }
            let proc = r
                .get("name")
                .and_then(og_nf::task_name)
                .map_or("", |(p, _)| p);
            let step = proc_step.get(proc).copied().unwrap_or("");
            text.push_str(&format!(
                "{step}\t{proc}\t{}\n",
                r.get("status").unwrap_or("")
            ));
        });
        if let Err(e) = r {
            return fail(err, format!("cannot read {t}: {e}"));
        }
    }
    let _ = out.write_all(text.as_bytes());
    0
}

fn eta(p: &Parsed, env: &Env, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let cores = match p.one("cores").map(pycompat::num::parse_int) {
        None => None,
        Some(Ok(c)) if c >= 0 => Some(usize::try_from(c).unwrap_or(usize::MAX)),
        Some(_) => {
            return usage_error(
                err,
                "eta",
                &format!(
                    "argument --cores: invalid int value: '{}'",
                    p.one("cores").unwrap_or_default()
                ),
            );
        }
    };
    let steps = match steps_of(p) {
        Ok(s) => s,
        Err(e) => return fail(err, e),
    };
    let outdir = PathBuf::from(&p.positionals[0]);
    let est = match estimate(
        &outdir,
        &steps,
        &EtaArgs {
            sample: p.one("sample"),
            cores,
            box_: p.one("box"),
            now: env.now,
            tz: &env.tz,
        },
    ) {
        Ok(e) => e,
        Err(e) => return fail(err, e),
    };
    let text = json_text(&est.to_json());
    match p.one("out") {
        Some(o) => {
            if let Err(e) = write_atomic(Path::new(o), &text) {
                return fail(err, format!("cannot write {o}: {e}"));
            }
        }
        None => {
            let _ = out.write_all(text.as_bytes());
        }
    }
    0
}

fn headroom(p: &Parsed, env: &Env, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let steps = match steps_of(p) {
        Ok(s) => s,
        Err(e) => return fail(err, e),
    };
    let events = match crate::files::read_events(p.one("events").map(Path::new)) {
        Ok(e) => e,
        Err(e) => return fail(err, e),
    };
    let rows: Vec<_> = p
        .all("trace")
        .into_iter()
        .flat_map(|t| read_trace(Path::new(t), &env.tz))
        .collect();
    let req = crate::headroom::Request {
        events: &events,
        trace_rows: &rows,
        work_dir: p.one("work-dir"),
        outdir: p.one("outdir"),
        max_memory: p.one("max-memory"),
        disk_warn_percent: p.one("disk-warn-percent"),
        mem_warn_percent: p.one("mem-warn-percent"),
    };
    match crate::headroom::snapshot(&crate::host::System, &steps, &req) {
        Ok(snap) => {
            let text: String = crate::headroom::lines(&snap)
                .iter()
                .map(|l| format!("{l}\n"))
                .collect();
            let _ = out.write_all(text.as_bytes());
            0
        }
        Err(e) => usage_error(err, "headroom", &e.0),
    }
}

/// `--previous DIR` names the previous run's folder; `--newest-history
/// SAMPLE_DIR` lets og pick it by [`crate::delta::newest_history`], so the
/// report step and the summary share one rule. Neither: no previous run.
fn delta(p: &Parsed, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let previous: Option<PathBuf> = match (p.one("previous"), p.one("newest-history")) {
        (Some(_), Some(_)) => {
            return usage_error(
                err,
                "delta",
                "argument --newest-history: not allowed with argument --previous",
            );
        }
        (Some(d), None) => Some(PathBuf::from(d)),
        (None, Some(sample_dir)) => crate::delta::newest_history(Path::new(sample_dir)),
        (None, None) => None,
    };
    let a = crate::delta::DeltaArgs {
        sample: p.one("sample").unwrap_or_default(),
        summary: Path::new(p.one("summary").unwrap_or_default()),
        manifest: p.one("manifest").map(Path::new),
        previous: previous.as_deref(),
        out: p.one("out").unwrap_or_default(),
        recheck: p.one("recheck").map(Path::new),
    };
    match crate::delta::run(&a) {
        Ok(text) => {
            let _ = out.write_all(text.as_bytes());
            0
        }
        Err(crate::delta::DeltaError::Summary(e)) => {
            let _ = writeln!(err, "ERROR: cannot read a summary: {e}");
            2
        }
        Err(crate::delta::DeltaError::Write(path, e)) => {
            fail(err, format!("cannot write {path}: {e}"))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn call(args: &[&str]) -> (u8, String, String) {
        let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
        let env = Env {
            tz: TimeZone::UTC,
            now: 0.0,
        };
        let (mut o, mut e) = (Vec::new(), Vec::new());
        let rc = run(&args, &env, &mut o, &mut e);
        (
            rc,
            String::from_utf8(o).unwrap(),
            String::from_utf8(e).unwrap(),
        )
    }

    #[test]
    fn usage_errors_exit_2() {
        let (rc, _, e) = call(&["state", "--sample", "S"]);
        assert_eq!(rc, 2);
        assert!(
            e.starts_with("usage:") && e.contains("required: --sample-dir"),
            "{e}"
        );
        assert_eq!(call(&["state", "--samp", "S", "--sample-dir", "d"]).0, 2);
        assert_eq!(call(&["nope"]).0, 2);
        assert_eq!(call(&["eta"]).0, 2);
        assert_eq!(call(&["eta", "o", "--cores", "x"]).0, 2);
        assert_eq!(call(&["delta", "--help"]).0, 0);
    }

    #[test]
    fn newest_history_picks_the_newest_folder_with_a_summary() {
        let dir = std::env::temp_dir().join(format!("og-run-newest-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let h = dir.join("S/history");
        // <stamp>-02 sorts after <stamp>; the newest folder lost its
        // summary.json half way through a copy, so it is passed over.
        for (name, summary) in [
            ("2026-01-01T000000Z", true),
            ("2026-01-01T000000Z-02", true),
            ("2026-01-02T000000Z", false),
        ] {
            std::fs::create_dir_all(h.join(name)).unwrap();
            if summary {
                std::fs::write(
                    h.join(name).join("summary.json"),
                    r#"{"generated_utc": "before", "sections": {}}"#,
                )
                .unwrap();
            }
        }
        std::fs::write(
            dir.join("S/summary.json"),
            r#"{"generated_utc": "now", "sections": {}}"#,
        )
        .unwrap();
        let (sd, sum, outp) = (
            dir.join("S").to_string_lossy().into_owned(),
            dir.join("S/summary.json").to_string_lossy().into_owned(),
            dir.join("S/S_delta").to_string_lossy().into_owned(),
        );
        let base = ["delta", "--sample", "S", "--summary", &sum, "--out", &outp];
        let previous = || {
            let j = std::fs::read_to_string(dir.join("S/S_delta.json")).unwrap();
            match pycompat::json::loads(&j).unwrap() {
                Value::Object(o) => o.get("previous").cloned(),
                _ => None,
            }
        };
        let mut a = base.to_vec();
        a.extend(["--newest-history", &sd]);
        let (rc, _, e) = call(&a);
        assert_eq!(rc, 0, "{e}");
        assert_eq!(previous(), Some(Value::Str("2026-01-01T000000Z-02".into())));
        // No history folder at all: no previous run, not an error.
        std::fs::remove_dir_all(&h).unwrap();
        let (rc, o, e) = call(&a);
        assert_eq!(rc, 0, "{e}");
        assert_eq!(previous(), Some(Value::Null));
        assert!(
            o.contains("No previous run of this sample on record"),
            "{o}"
        );
        // --previous and --newest-history together is a usage error.
        a.extend(["--previous", &sd]);
        let (rc, _, e) = call(&a);
        assert_eq!(rc, 2);
        assert!(e.contains("not allowed with argument --previous"), "{e}");
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn tasks_name_the_step_that_owns_each_process() {
        let dir = std::env::temp_dir().join(format!("og-run-tasks-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let t = dir.join("trace_1.txt");
        std::fs::write(
            &t,
            "task_id\tname\tstatus\n1\tW:PRS (s1)\tCOMPLETED\n\n2\tW:SUB:ALIGN_MINIMAP2 (s1)\tCACHED\n3\tW:NOT_A_STEP\tCOMPLETED\n",
        )
        .unwrap();
        let steps = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../assets/steps.tsv");
        let (rc, out, e) = call(&[
            "tasks",
            "--trace",
            &t.to_string_lossy(),
            "--steps",
            &steps.to_string_lossy(),
        ]);
        assert_eq!(rc, 0, "{e}");
        // PRS is listed by steps 25 and 26: the first owns it.
        assert_eq!(
            out,
            "25\tPRS\tCOMPLETED\n02\tALIGN_MINIMAP2\tCACHED\n\tNOT_A_STEP\tCOMPLETED\n"
        );
        let (rc, _, e) = call(&[
            "tasks",
            "--trace",
            &dir.join("none.txt").to_string_lossy(),
            "--steps",
            &steps.to_string_lossy(),
        ]);
        assert_eq!(rc, 1, "{e}");
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
