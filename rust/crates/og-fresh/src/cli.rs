//! `og dev freshness`: the monthly report of pins that fell behind.

use crate::http::Client;
use crate::lookups::Ctx;
use crate::pins::Pins;
use crate::registry::Online;
use crate::report::{SECTIONS, build, render};
use std::io::Write;
use std::path::PathBuf;

pub const USAGE: &str =
    "usage: og dev freshness [--out FILE] [--sections a,b] [--versions F] [--setup F]
       og dev freshness --update-issue BODY_FILE [--repo OWNER/NAME]
       og dev freshness --self-test
  sections: databases, clinvar (the ClinVar row alone), images, digests
";

fn usage_error(err: &mut dyn Write, msg: &str) -> u8 {
    let _ = write!(err, "{USAGE}");
    let _ = writeln!(err, "og dev freshness: error: {msg}");
    2
}

/// Options as `--name value` or `--name=value`; flags take no value. A
/// repeated option keeps its last value; an unknown one is refused.
pub fn parse(
    args: &[String],
    opts: &[&'static str],
    flags: &[&'static str],
) -> Result<Vec<(&'static str, Option<String>)>, String> {
    let mut out = vec![];
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        let Some(body) = a.strip_prefix("--").filter(|b| !b.is_empty()) else {
            return Err(format!("unrecognized arguments: {a}"));
        };
        let (name, inline) = match body.split_once('=') {
            Some((n, v)) => (n, Some(v.to_owned())),
            None => (body, None),
        };
        if let Some(f) = flags.iter().find(|f| **f == name) {
            if inline.is_some() {
                return Err(format!("argument --{f}: ignored explicit argument"));
            }
            out.push((*f, None));
        } else if let Some(o) = opts.iter().find(|o| **o == name) {
            let v = match inline {
                Some(v) => v,
                None => {
                    i += 1;
                    args.get(i)
                        .cloned()
                        .ok_or_else(|| format!("argument --{o}: expected one argument"))?
                }
            };
            out.push((*o, Some(v)));
        } else {
            return Err(format!("unrecognized arguments: {a}"));
        }
        i += 1;
    }
    Ok(out)
}

pub fn one<'a>(p: &'a [(&'static str, Option<String>)], k: &str) -> Option<&'a str> {
    p.iter()
        .rev()
        .find(|(n, _)| *n == k)
        .and_then(|(_, v)| v.as_deref())
}

pub fn flag(p: &[(&'static str, Option<String>)], k: &str) -> bool {
    p.iter().any(|(n, _)| *n == k)
}

/// `og dev freshness` on this machine.
pub fn main(args: &[String]) -> u8 {
    run(
        args,
        &mut std::io::stdout().lock(),
        &mut std::io::stderr().lock(),
    )
}

pub fn run(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let p = match parse(
        args,
        &[
            "versions",
            "setup",
            "sections",
            "out",
            "update-issue",
            "repo",
        ],
        &["self-test", "help"],
    ) {
        Ok(p) => p,
        Err(e) => return usage_error(err, &e),
    };
    if flag(&p, "help") {
        let _ = write!(out, "{USAGE}");
        return 0;
    }
    if flag(&p, "self-test") {
        return crate::selftest::run(out);
    }
    let client = Client::default();
    if let Some(body) = one(&p, "update-issue") {
        let repo = one(&p, "repo").map(str::to_owned).or_else(|| {
            std::env::var("GITHUB_REPOSITORY")
                .ok()
                .filter(|r| !r.is_empty())
        });
        let Some(repo) = repo else {
            return usage_error(err, "--update-issue needs --repo or GITHUB_REPOSITORY");
        };
        let text = match std::fs::read_to_string(body) {
            Ok(t) => t,
            Err(e) => {
                let _ = writeln!(err, "ERROR: cannot read {body}: {e}");
                return 1;
            }
        };
        return match crate::issue::update(&client, &repo, &text) {
            Ok(_) => 0,
            Err(e) => {
                let _ = writeln!(err, "ERROR: {e}");
                1
            }
        };
    }
    let sections: Vec<&str> = one(&p, "sections").map_or_else(
        || SECTIONS.to_vec(),
        |s| s.split(',').filter(|x| !x.is_empty()).collect(),
    );
    let bad: Vec<&str> = sections
        .iter()
        .copied()
        .filter(|s| !SECTIONS.contains(s))
        .collect();
    if !bad.is_empty() || sections.is_empty() {
        return usage_error(
            err,
            &format!(
                "unknown section(s) {bad:?}; choose from {}",
                SECTIONS.join(", ")
            ),
        );
    }
    let root = crate::checkout_root();
    let path = |opt: &str, rel: &str| -> Option<PathBuf> {
        one(&p, opt)
            .map(PathBuf::from)
            .or_else(|| root.as_ref().map(|r| r.join(rel)))
    };
    let versions = path("versions", "versions.env");
    // common.sh sits in the checkout of the versions.env read.
    let common = root.as_ref().map_or_else(
        || {
            versions
                .as_ref()
                .and_then(|v| v.parent())
                .map(|d| d.join("scripts/lib/common.sh"))
        },
        |r| Some(r.join("scripts/lib/common.sh")),
    );
    let (Some(versions), Some(setup), Some(common)) =
        (versions, path("setup", "scripts/setup.sh"), common)
    else {
        let _ = writeln!(
            err,
            "ERROR: cannot find the OwnGenome checkout above og or the working directory"
        );
        return 2;
    };
    let pins = match Pins::read(&versions) {
        Ok(p) => p,
        Err(e) => {
            let _ = writeln!(err, "ERROR: cannot read {}: {e}", versions.display());
            return 1;
        }
    };
    let ctx = Ctx {
        fetch: &client,
        registry: &Online(&client),
        pins: &pins,
        versions: &versions,
        setup: &setup,
        common: &common,
        now: jiff::Timestamp::now(),
    };
    let report = build(&ctx, &sections);
    let text = render(&report, &sections);
    if let Some(o) = one(&p, "out")
        && let Err(e) = std::fs::write(o, format!("{text}\n"))
    {
        let _ = writeln!(err, "ERROR: cannot write {o}: {e}");
        return 1;
    }
    let _ = writeln!(out, "{text}");
    if !report.errors.is_empty() {
        let _ = writeln!(err, "ERROR: {} lookup(s) failed.", report.errors.len());
        return 1;
    }
    0
}
