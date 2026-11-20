//! `og report render`: the text and HTML reports of one run, and the files
//! written beside them (summary.json, the search index, the charts).
//!
//! ```text
//! og report render (--summary IN.json | --sample S --sample-dir DIR [--json OUT.json]
//!                   [--declared-sex SEX] [--assets DIR])
//!                  [-o REPORT.txt] [-o REPORT.html]
//!                  [--search-index INDEX.json] [--charts-json CHARTS.json]
//! ```
//!
//! With `--summary` it renders that summary. Without it, it collects the
//! summary of the sample folder first (`og summary collect`) and keeps it in
//! `--json` (default `DIR/summary.json`). `-o` writes the text report for a
//! `.txt` name and the HTML report for `.html`; it repeats. The HTML report
//! closes with every gene and rsID the results name, from the search index
//! `og report index` builds; `--sample-dir` adds the gene columns of the
//! step 23 and 31 tables to it. `--search-index` also writes that index, and
//! `--charts-json` the charts the HTML report embeds (`og report charts`).
//! Each file is written under a temporary name and renamed into place.
//!
//! A section that cannot be written (its data has a shape the report does
//! not read, or its writer stops on an internal error) says so in its place
//! in both reports, and the rest of the report is still written: one
//! warning on stderr counts those parts and names them, and the exit code
//! stays 0. The last line on stdout counts the summary's sections by state.

use std::io::Write;

use pycompat::json::{Dump, Value, dumps, loads};

use crate::chapters::{Failed, WriterFn, guarded, writer};
use crate::model::{Summary, py};
use crate::{charts, html, index, text};

const USAGE: &str = "usage: og report render (--summary IN.json | --sample SAMPLE --sample-dir SAMPLE_DIR \
[--json OUT.json] [--declared-sex SEX] [--assets ASSETS]) [-o REPORT.txt] [-o REPORT.html] \
[--search-index INDEX.json] [--charts-json CHARTS.json]";

#[derive(Debug, Default)]
struct Args {
    summary: Option<String>,
    sample: Option<String>,
    sample_dir: Option<String>,
    json: Option<String>,
    declared_sex: Option<String>,
    assets: Option<String>,
    out: Vec<String>,
    search_index: Option<String>,
    charts_json: Option<String>,
}

/// Reads the arguments, or says why they are wrong.
fn parse(args: &[String]) -> Result<Args, String> {
    let mut a = Args::default();
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        let (name, inline) = match arg.split_once('=') {
            Some((n, v)) if n.starts_with("--") => (n, Some(v.to_owned())),
            _ if arg.starts_with("-o") && arg.len() > 2 && !arg.starts_with("--") => {
                ("-o", Some(arg[2..].to_owned()))
            }
            _ => (arg.as_str(), None),
        };
        let mut value = || {
            inline
                .clone()
                .or_else(|| it.next().cloned())
                .ok_or_else(|| format!("argument {name}: expected one argument"))
        };
        match name {
            "--summary" => a.summary = Some(value()?),
            "--sample" => a.sample = Some(value()?),
            "--sample-dir" => a.sample_dir = Some(value()?),
            "--json" => a.json = Some(value()?),
            "--declared-sex" => a.declared_sex = Some(value()?),
            "--assets" => a.assets = Some(value()?),
            "-o" | "--out" => a.out.push(value()?),
            "--search-index" => a.search_index = Some(value()?),
            "--charts-json" => a.charts_json = Some(value()?),
            n => return Err(format!("unrecognized arguments: {n}")),
        }
    }
    for o in &a.out {
        if !(o.ends_with(".txt") || o.ends_with(".html")) {
            return Err(format!("-o {o}: the name must end in .txt or .html"));
        }
    }
    if a.summary.is_none() && (a.sample.is_none() || a.sample_dir.is_none()) {
        return Err("give --summary, or --sample and --sample-dir".to_owned());
    }
    Ok(a)
}

/// `Sections: missing 28, ok 20`: the sections counted by state, by name.
fn state_counts(s: &Summary) -> Result<String, String> {
    let mut counts: Vec<(String, usize)> = Vec::new();
    let Some(Value::Object(secs)) = s.get("sections") else {
        return Err("summary.sections is not an object".to_owned());
    };
    for (key, sec) in secs {
        let state = crate::model::need(sec, "state", &format!("sections.{key}"))
            .and_then(py)
            .map_err(|e| e.to_string())?;
        match counts.iter_mut().find(|(k, _)| *k == state) {
            Some((_, n)) => *n += 1,
            None => counts.push((state, 1)),
        }
    }
    counts.sort();
    let parts: Vec<String> = counts.iter().map(|(k, n)| format!("{k} {n}")).collect();
    Ok(format!("Sections: {}", parts.join(", ")))
}

fn write(path: &str, body: &str) -> std::io::Result<()> {
    let mut f = og_io::AtomicFile::create(path)?;
    f.write_all(body.as_bytes())?;
    f.commit()
}

/// A JSON value as the reports' side files hold it: one-space indent and a
/// closing newline.
fn json_text(v: &Value) -> Result<String, String> {
    let opts = Dump {
        indent: Some(1),
        ..Dump::default()
    };
    dumps(v, &opts)
        .map(|b| b + "\n")
        .map_err(|e| format!("{e:?}"))
}

/// Collects the summary of the sample folder with `og summary collect`,
/// keeps it in `dest` and returns its text.
fn collect(a: &Args, dest: &str, err: &mut dyn Write) -> Result<String, String> {
    let sample = a.sample.clone().unwrap_or_default();
    let mut argv = vec![
        "--sample".to_owned(),
        sample.clone(),
        "--sample-dir".to_owned(),
        a.sample_dir.clone().unwrap_or_default(),
        "--out".to_owned(),
        dest.to_owned(),
    ];
    if let Some(assets) = &a.assets {
        argv.extend(["--assets".to_owned(), assets.clone()]);
    }
    // An empty sex (a samplesheet row without one) is no sex.
    if let Some(sex) = a.declared_sex.as_ref().filter(|s| !s.is_empty()) {
        argv.extend(["--declared-sex".to_owned(), sex.clone()]);
    }
    let rc = og_summary::collect::cli(&argv, &mut std::io::sink(), &mut &mut *err);
    if rc != 0 {
        return Err(format!(
            "ERROR: og summary collect wrote no summary of {sample} (exit {rc})"
        ));
    }
    std::fs::read_to_string(dest).map_err(|e| format!("{dest}: {e}"))
}

/// `og report render ARGS`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    cli_using(args, out, err, writer)
}

/// `og report render ARGS` with the section writers of `writers` (a test
/// plants a fault through them).
pub fn cli_using(
    args: &[String],
    out: &mut dyn Write,
    err: &mut dyn Write,
    writers: WriterFn,
) -> u8 {
    if args.iter().any(|a| a == "-h" || a == "--help") {
        let _ = writeln!(out, "{USAGE}");
        return 0;
    }
    let a = match parse(args) {
        Ok(a) => a,
        Err(why) => {
            let _ = writeln!(err, "{USAGE}\nog report render: error: {why}");
            return 2;
        }
    };
    let fail = |err: &mut dyn Write, why: String| {
        let _ = writeln!(err, "og report render: {why}");
        1
    };
    let (path, body) = match &a.summary {
        Some(p) => (
            p.clone(),
            std::fs::read_to_string(p).map_err(|e| format!("{p}: {e}")),
        ),
        None => {
            let dest = a.json.clone().unwrap_or_else(|| {
                std::path::Path::new(a.sample_dir.as_deref().unwrap_or_default())
                    .join("summary.json")
                    .display()
                    .to_string()
            });
            let body = collect(&a, &dest, err);
            (dest, body)
        }
    };
    let body = match body {
        Ok(b) => b,
        Err(why) if why.starts_with("ERROR: ") => {
            let _ = writeln!(err, "{why}");
            return 1;
        }
        Err(why) => return fail(err, why),
    };
    let s = match Summary::parse(&body) {
        Ok(s) => s,
        Err(why) => return fail(err, format!("{path}: {why}")),
    };
    let mut failed: Vec<Failed> = Vec::new();
    // The search index, for the HTML report's appendix and --search-index;
    // a text report alone does not read the step 23 and 31 tables.
    let need_index = a.search_index.is_some() || a.out.iter().any(|o| o.ends_with(".html"));
    let idx = if need_index {
        let built = guarded(|| {
            let root = loads(&body).map_err(|e| crate::model::Malformed(e.to_string()))?;
            index::build(&root, a.sample_dir.as_deref(), &[], &mut &mut *err)
                .map_err(crate::model::Malformed)
        });
        match built {
            Ok(i) => Some(i),
            Err(why) => {
                failed.push(Failed {
                    part: "terms".to_owned(),
                    why: why.to_string(),
                });
                None
            }
        }
    } else {
        None
    };
    for o in &a.out {
        let made = if o.ends_with(".txt") {
            text::report_using(&s, writers)
        } else {
            // Without an index the appendix is left out, which `failed` names.
            html::report_using(&s, idx.as_ref(), writers)
        };
        let done = made.map_err(|e| format!("{path}: {e}")).and_then(|(b, f)| {
            write(o, &b).map_err(|e| format!("{o}: {e}"))?;
            Ok(f)
        });
        match done {
            Ok(f) => {
                for x in f {
                    if !failed.contains(&x) {
                        failed.push(x);
                    }
                }
            }
            Err(why) => return fail(err, why),
        }
    }
    if let Some(dest) = &a.search_index {
        let Some(i) = &idx else {
            return fail(err, format!("{dest}: the search index could not be built"));
        };
        if let Err(why) = json_text(i).and_then(|b| write(dest, &b).map_err(|e| e.to_string())) {
            return fail(err, format!("{dest}: {why}"));
        }
    }
    if let Some(dest) = &a.charts_json {
        let made = guarded(|| charts::json(&s))
            .map_err(|e| format!("{path}: {e}"))
            .and_then(|v| json_text(&v))
            .and_then(|b| write(dest, &b).map_err(|e| format!("{dest}: {e}")));
        if let Err(why) = made {
            return fail(err, why);
        }
    }
    if !failed.is_empty() {
        let n = failed.len();
        let named: Vec<String> = failed
            .iter()
            .map(|f| format!("{} ({})", f.part, f.why))
            .collect();
        let _ = writeln!(
            err,
            "og report render: WARNING: {n} part{} of the report could not be written and say{} so in \
             place: {}",
            if n == 1 { "" } else { "s" },
            if n == 1 { "s" } else { "" },
            named.join("; ")
        );
    }
    match state_counts(&s) {
        Ok(line) => {
            let _ = writeln!(out, "{line}");
            0
        }
        Err(why) => fail(err, format!("{path}: {why}")),
    }
}
