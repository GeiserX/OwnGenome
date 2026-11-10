//! `og map <verb>`: the genome map (`rust/crates/og-map`).
//!
//! `og map build --root DIR --run RUN --assets DIR --cache-dir DIR < facts.json > map.json.gz`
//!
//! Exit codes: 0 the document is on stdout; 2 usage; 3 no map for this run
//! (the reason on stderr); 1 an internal error. A layer that cannot be read
//! is not an exit code: it is that layer's `unreadable` state.

use og_map::{BuildError, Request};
use std::io::{Read, Write};
use std::path::Path;
use std::process::ExitCode;

const USAGE: &str = "usage: og map build --root DIR --run RUN --assets DIR --cache-dir DIR < facts.json > map.json.gz\n";

fn usage(msg: &str) -> ExitCode {
    eprint!("{USAGE}");
    eprintln!("og map: error: {msg}");
    ExitCode::from(2)
}

/// The values of `build`'s four options, in order.
#[derive(Debug, Default, PartialEq, Eq)]
struct Build {
    root: String,
    run: String,
    assets: String,
    cache_dir: String,
}

fn parse_build(args: &[String]) -> Result<Build, String> {
    let mut b = Build::default();
    let mut seen = [false; 4];
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let (key, inline) = match a.split_once('=') {
            Some((k, v)) if k.starts_with("--") => (k, Some(v.to_owned())),
            _ => (a.as_str(), None),
        };
        let k = match key {
            "--root" => 0,
            "--run" => 1,
            "--assets" => 2,
            "--cache-dir" => 3,
            _ => return Err(format!("unrecognized arguments: {a}")),
        };
        let v = match inline {
            Some(v) => v,
            None => it
                .next()
                .cloned()
                .ok_or_else(|| format!("argument {key}: expected one argument"))?,
        };
        seen[k] = true;
        *[&mut b.root, &mut b.run, &mut b.assets, &mut b.cache_dir][k] = v;
    }
    let names = ["--root", "--run", "--assets", "--cache-dir"];
    let missing: Vec<&str> = names
        .iter()
        .zip(seen)
        .filter(|(_, s)| !s)
        .map(|(n, _)| *n)
        .collect();
    if !missing.is_empty() {
        return Err(format!(
            "the following arguments are required: {}",
            missing.join(", ")
        ));
    }
    Ok(b)
}

pub fn run(args: &[String], tree_hash: &str) -> ExitCode {
    match args.first().map(String::as_str) {
        Some("-h" | "--help") => {
            print!("{USAGE}");
            ExitCode::SUCCESS
        }
        Some("build") => build(&args[1..], tree_hash),
        Some(other) => usage(&format!("unknown verb '{other}'")),
        None => usage("the following arguments are required: verb"),
    }
}

fn build(args: &[String], tree_hash: &str) -> ExitCode {
    let b = match parse_build(args) {
        Ok(b) => b,
        Err(e) => return usage(&e),
    };
    let mut stdin = Vec::new();
    if let Err(e) = std::io::stdin().lock().read_to_end(&mut stdin) {
        eprintln!("og map build: stdin: {e}");
        return ExitCode::FAILURE;
    }
    let req = Request {
        root: Path::new(&b.root),
        run: &b.run,
        assets: Path::new(&b.assets),
        made_by: tree_hash,
        stdin: &stdin,
    };
    match og_map::build_cached(&req, Path::new(&b.cache_dir)) {
        Ok(doc) => {
            let mut out = std::io::stdout().lock();
            match out.write_all(&doc).and_then(|()| out.flush()) {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    eprintln!("og map build: stdout: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        Err(BuildError::Usage(e)) => usage(&e),
        Err(BuildError::NoMap(e)) => {
            eprintln!("og map build: {e}");
            ExitCode::from(3)
        }
        Err(BuildError::Internal(e)) => {
            eprintln!("og map build: {e}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(s: &str) -> Vec<String> {
        s.split_whitespace().map(str::to_owned).collect()
    }

    #[test]
    fn build_arguments() {
        let b = parse_build(&args("--root /g --run DEMO-001 --assets=/a --cache-dir /c")).unwrap();
        assert_eq!(
            b,
            Build {
                root: "/g".into(),
                run: "DEMO-001".into(),
                assets: "/a".into(),
                cache_dir: "/c".into()
            }
        );
        assert_eq!(
            parse_build(&args("--root /g --run X")).unwrap_err(),
            "the following arguments are required: --assets, --cache-dir"
        );
        assert_eq!(
            parse_build(&args("--root")).unwrap_err(),
            "argument --root: expected one argument"
        );
        assert_eq!(
            parse_build(&args("--family F")).unwrap_err(),
            "unrecognized arguments: --family"
        );
    }
}
