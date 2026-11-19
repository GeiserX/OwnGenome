//! `og dev freshness`: the monthly report of pins that fell behind.

use crate::http::Client;
use crate::lookups::Ctx;
use crate::pins::Pins;
use crate::registry::Online;
use crate::report::{Report, SECTIONS, build, render_with};
use std::io::Write;
use std::path::PathBuf;

pub const USAGE: &str =
    "usage: og dev freshness [--out FILE] [--sections a,b] [--versions F] [--setup F]
                        [--unreachable-warns]
       og dev freshness --update-issue BODY_FILE [--repo OWNER/NAME]
       og dev freshness --self-test
  sections: databases, clinvar (the ClinVar row alone), images, digests (default: all)
  --out FILE        also write the report to FILE; it is always printed on stdout
  --versions F      the pins to check (default: versions.env of the checkout above og)
  --setup F         the setup script with the download URLs (default: scripts/setup.sh)
  --update-issue B  put the report file B in the freshness issue of --repo
                    (default: $GITHUB_REPOSITORY)
  --unreachable-warns
                    an upstream that does not answer (no connection, a timeout,
                    HTTP 5xx or 429 after the retries, or HTTP 403) is a
                    warning, not a failure; every other failed lookup still
                    exits 1
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
    // A help word after `--` is an argument, as in og's own help rule.
    if args
        .iter()
        .take_while(|a| *a != "--")
        .any(|a| a == "-h" || a == "--help")
    {
        let _ = write!(out, "{USAGE}");
        return 0;
    }
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
        &["self-test", "unreachable-warns"],
    ) {
        Ok(p) => p,
        Err(e) => return usage_error(err, &e),
    };
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
    finish(
        &report,
        &sections,
        flag(&p, "unreachable-warns"),
        one(&p, "out"),
        out,
        err,
    )
}

/// Writes the report to `out_file` and `out`, and gives the exit status:
/// 1 when a lookup failed, not counting, with `warns`, the upstreams that
/// did not answer.
pub fn finish(
    report: &Report,
    sections: &[&str],
    warns: bool,
    out_file: Option<&str>,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> u8 {
    let text = render_with(report, sections, warns);
    if let Some(o) = out_file
        && let Err(e) = std::fs::write(o, format!("{text}\n"))
    {
        let _ = writeln!(err, "ERROR: cannot write {o}: {e}");
        return 1;
    }
    let _ = writeln!(out, "{text}");
    let hard = report.hard_errors(warns);
    let soft = report.errors.len() - hard;
    if soft > 0 {
        let _ = writeln!(
            err,
            "WARNING: {soft} upstream lookup(s) got no answer; not counted (--unreachable-warns)."
        );
    }
    if hard > 0 {
        let _ = writeln!(err, "ERROR: {hard} lookup(s) failed.");
        return 1;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::LookupError;

    fn report(errors: &[LookupError]) -> Report {
        let mut r = Report::default();
        r.db_error("Ensembl / VEP cache", &errors[0]);
        for (i, e) in errors.iter().enumerate().skip(1) {
            r.error(&format!("IMAGE_{i}"), e);
        }
        r
    }

    fn exit(r: &Report, warns: bool) -> (u8, String, String) {
        let (mut o, mut e) = (Vec::new(), Vec::new());
        let rc = finish(r, &["databases", "images"], warns, None, &mut o, &mut e);
        let s = |b: Vec<u8>| String::from_utf8(b).unwrap();
        (rc, s(o), s(e))
    }

    fn down() -> LookupError {
        LookupError::Unreachable("https://rest.ensembl.org/x: no answer (timed out)".into())
    }

    #[test]
    fn an_outage_alone_warns_only_when_asked() {
        let r = report(&[down()]);
        let (rc, text, err) = exit(&r, true);
        assert_eq!(rc, 0, "{err}");
        assert!(!text.contains("### Errors"), "{text}");
        let warn = text
            .split_once("### Upstreams that did not answer")
            .expect("the warning heading")
            .1;
        assert!(warn.contains("| Ensembl / VEP cache | https://rest.ensembl.org/x: no answer"));
        assert!(text.contains("| NO ANSWER |"), "{text}");
        assert!(err.contains("WARNING: 1 upstream"), "{err}");
        // Without the flag the run is what it always was.
        let (rc, text, err) = exit(&r, false);
        assert_eq!(rc, 1);
        assert!(
            text.contains("### Errors") && text.contains("| ERROR | lookup failed, see Errors |")
        );
        assert!(!text.contains("did not answer"), "{text}");
        assert!(err.contains("ERROR: 1 lookup(s) failed."), "{err}");
        assert_eq!(
            text.trim_end(),
            crate::report::render(&r, &["databases", "images"]).trim_end()
        );
    }

    #[test]
    fn a_wrong_answer_still_fails_with_the_flag() {
        for wrong in [
            LookupError::Failed("https://x/y: HTTP 404".into()),
            LookupError::Shape("no key \"data\"".into()),
        ] {
            let r = report(&[down(), wrong]);
            for warns in [true, false] {
                let (rc, text, err) = exit(&r, warns);
                assert_eq!(rc, 1, "warns={warns}: {text}");
                let errors = text.split_once("### Errors").expect("an Errors table").1;
                let errors = errors.split("### Upstreams").next().unwrap_or("");
                assert!(errors.contains("| IMAGE_1 |"), "{text}");
                assert_eq!(errors.contains("Ensembl"), !warns, "{text}");
                assert!(
                    err.contains("ERROR: 1 lookup(s) failed.") || !warns,
                    "{err}"
                );
            }
        }
    }

    #[test]
    fn the_flag_is_an_option() {
        let p = parse(
            &["--unreachable-warns".into()],
            &["out"],
            &["unreachable-warns"],
        )
        .unwrap();
        assert!(flag(&p, "unreachable-warns"));
        assert!(
            parse(
                &["--unreachable-warns=1".into()],
                &[],
                &["unreachable-warns"]
            )
            .is_err()
        );
    }
}
