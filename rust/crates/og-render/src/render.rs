//! `og report render`: the text and HTML reports from summary.json, which
//! `bin/render_report.py` writes today.
//!
//! ```text
//! og report render --summary IN.json -o REPORT.txt -o REPORT.html
//! ```
//!
//! `-o` writes the text report for a `.txt` name and the HTML report for
//! `.html`; it repeats. Each file is written under a temporary name and
//! renamed into place. The last line on stdout counts the summary's sections
//! by state. No caller uses the verb yet: it renders the frame and the
//! sections ported so far, and refuses the options of the parts still in
//! Python (collecting the summary, the search index, the charts).

use std::io::Write;

use crate::model::{Summary, py};
use crate::{html, text};

const USAGE: &str = "usage: og report render --summary IN.json [-o REPORT.txt] [-o REPORT.html]";

/// The options `bin/render_report.py` takes that this verb does not yet.
const NOT_PORTED: [&str; 7] = [
    "--sample",
    "--sample-dir",
    "--json",
    "--declared-sex",
    "--assets",
    "--search-index",
    "--charts-json",
];

#[derive(Debug, Default)]
struct Args {
    summary: Option<String>,
    out: Vec<String>,
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
            "-o" | "--out" => a.out.push(value()?),
            n if NOT_PORTED.contains(&n) => {
                return Err(format!(
                    "{n} is not ported yet: bin/render_report.py still collects the summary and writes \
                     the search index and the charts"
                ));
            }
            n => return Err(format!("unrecognized arguments: {n}")),
        }
    }
    for o in &a.out {
        if !(o.ends_with(".txt") || o.ends_with(".html")) {
            return Err(format!("-o {o}: the name must end in .txt or .html"));
        }
    }
    if a.summary.is_none() {
        return Err("give --summary".to_owned());
    }
    Ok(a)
}

/// `Sections: missing 28, ok 20`: the sections counted by state, by name.
fn state_counts(s: &Summary) -> Result<String, String> {
    let mut counts: Vec<(String, usize)> = Vec::new();
    let Some(pycompat::json::Value::Object(secs)) = s.get("sections") else {
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

/// `og report render ARGS`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
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
    let path = a.summary.as_deref().unwrap_or_default();
    let fail = |err: &mut dyn Write, why: String| {
        let _ = writeln!(err, "og report render: {why}");
        1
    };
    let s = match std::fs::read_to_string(path)
        .map_err(|e| e.to_string())
        .and_then(|t| Summary::parse(&t).map_err(|e| e.to_string()))
    {
        Ok(s) => s,
        Err(why) => return fail(err, format!("{path}: {why}")),
    };
    for o in &a.out {
        let body = if o.ends_with(".txt") {
            text::report(&s)
        } else {
            html::report(&s)
        };
        let done = body
            .map_err(|e| format!("{path}: {e}"))
            .and_then(|b| write(o, &b).map_err(|e| format!("{o}: {e}")));
        if let Err(why) = done {
            return fail(err, why);
        }
    }
    match state_counts(&s) {
        Ok(line) => {
            let _ = writeln!(out, "{line}");
            0
        }
        Err(why) => fail(err, format!("{path}: {why}")),
    }
}
