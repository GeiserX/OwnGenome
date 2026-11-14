//! `og map <verb>`: the genome map (`rust/crates/og-map`).
//!
//! `og map build --root DIR --run RUN --assets DIR --cache-dir DIR < facts.json > map.json.gz`
//! `og map build ... --predicate-table < facts.json > table.tsv`: whether each row of the
//! document (and each installed gene) draws one by one at each level (`og_map::predicate`),
//! for the browser check that holds the page's `drawsAt()` to the same answers
//! `og map genes --root DIR --cache-dir DIR (--chrom CHROM | --names) > genes.json.gz`
//!
//! Exit codes: 0 the document is on stdout; 2 usage; 3 no map for this run,
//! or no gene files (the reason on stderr); 1 an internal error. A layer that cannot be read
//! is not an exit code: it is that layer's `unreadable` state.

use og_map::{BuildError, Request};
use std::io::{Read, Write};
use std::path::Path;
use std::process::ExitCode;

const USAGE: &str = "usage: og map build --root DIR --run RUN --assets DIR --cache-dir DIR [--predicate-table] < facts.json > map.json.gz
       og map genes --root DIR --cache-dir DIR (--chrom CHROM | --names) > genes.json.gz
";

fn usage(msg: &str) -> ExitCode {
    eprint!("{USAGE}");
    eprintln!("og map: error: {msg}");
    ExitCode::from(2)
}

/// Parses `--name VALUE`, `--name=VALUE` and bare flags, argparse's way.
/// Returns each option's value in the order of `valued`, and whether each
/// of `flags` was given.
fn options(
    args: &[String],
    valued: &[&str],
    flags: &[&str],
) -> Result<(Vec<Option<String>>, Vec<bool>), String> {
    let mut values = vec![None; valued.len()];
    let mut set = vec![false; flags.len()];
    let mut it = args.iter();
    while let Some(a) = it.next() {
        if let Some(k) = flags.iter().position(|f| f == a) {
            set[k] = true;
            continue;
        }
        let (key, inline) = match a.split_once('=') {
            Some((k, v)) if k.starts_with("--") => (k, Some(v.to_owned())),
            _ => (a.as_str(), None),
        };
        let Some(k) = valued.iter().position(|n| *n == key) else {
            return Err(format!("unrecognized arguments: {a}"));
        };
        let v = match inline {
            Some(v) => v,
            None => it
                .next()
                .cloned()
                .ok_or_else(|| format!("argument {key}: expected one argument"))?,
        };
        values[k] = Some(v);
    }
    let missing: Vec<&str> = valued
        .iter()
        .zip(&values)
        .filter(|(_, v)| v.is_none())
        .map(|(n, _)| *n)
        .collect();
    if !missing.is_empty() {
        return Err(format!(
            "the following arguments are required: {}",
            missing.join(", ")
        ));
    }
    Ok((values, set))
}

/// The values of `build`'s four options, in order, and its one flag.
#[derive(Debug, Default, PartialEq, Eq)]
struct Build {
    root: String,
    run: String,
    assets: String,
    cache_dir: String,
    /// `--predicate-table`: the predicate's answers instead of the document.
    predicate_table: bool,
}

fn parse_build(args: &[String]) -> Result<Build, String> {
    let (v, set) = options(
        args,
        &["--root", "--run", "--assets", "--cache-dir"],
        &["--predicate-table"],
    )?;
    let [root, run, assets, cache_dir] = <[Option<String>; 4]>::try_from(v)
        .expect("one value per option")
        .map(Option::unwrap_or_default);
    Ok(Build {
        root,
        run,
        assets,
        cache_dir,
        predicate_table: set[0],
    })
}

/// What `genes` was asked for.
#[derive(Debug, PartialEq, Eq)]
struct Genes {
    root: String,
    cache_dir: String,
    /// The chromosome, or `None` for the names file.
    chrom: Option<String>,
}

/// Reads the options of `og map genes`.
fn parse_genes(args: &[String]) -> Result<Genes, String> {
    // --chrom is optional, so it is parsed apart from the two required options.
    let mut chrom = None;
    let mut rest = Vec::with_capacity(args.len());
    let mut it = args.iter();
    while let Some(a) = it.next() {
        if a == "--chrom" {
            let v = it.next().ok_or("argument --chrom: expected one argument")?;
            chrom = Some(v.clone());
        } else if let Some(v) = a.strip_prefix("--chrom=") {
            chrom = Some(v.to_owned());
        } else {
            rest.push(a.clone());
        }
    }
    let (v, set) = options(&rest, &["--root", "--cache-dir"], &["--names"])?;
    let [root, cache_dir] = <[Option<String>; 2]>::try_from(v)
        .expect("one value per option")
        .map(Option::unwrap_or_default);
    match (chrom.is_some(), set[0]) {
        (true, true) => Err("argument --names: not allowed with argument --chrom".to_owned()),
        (false, false) => Err("one of the arguments --chrom --names is required".to_owned()),
        _ => Ok(Genes {
            root,
            cache_dir,
            chrom,
        }),
    }
}

pub fn run(args: &[String], tree_hash: &str) -> ExitCode {
    match args.first().map(String::as_str) {
        Some("-h" | "--help") => {
            print!("{USAGE}");
            ExitCode::SUCCESS
        }
        Some("build") => build(&args[1..], tree_hash),
        Some("genes") => genes(&args[1..], tree_hash),
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
    if b.predicate_table {
        return finish(
            "build",
            og_map::predicate::build_table(&req).map(String::into_bytes),
        );
    }
    finish("build", og_map::build_cached(&req, Path::new(&b.cache_dir)))
}

/// Runs `og map genes`: writes the requested gene file to stdout.
fn genes(args: &[String], tree_hash: &str) -> ExitCode {
    let g = match parse_genes(args) {
        Ok(g) => g,
        Err(e) => return usage(&e),
    };
    let want = match g.chrom.as_deref() {
        Some(c) => og_map::genes::Want::chrom(c),
        None => Ok(og_map::genes::Want::Names),
    };
    let bytes = want.and_then(|w| {
        og_map::genes::genes_cached(Path::new(&g.root), Path::new(&g.cache_dir), tree_hash, w)
    });
    finish("genes", bytes)
}

/// Writes the bytes to stdout, or reports why there are none, with the
/// verb's exit code.
fn finish(verb: &str, result: Result<Vec<u8>, BuildError>) -> ExitCode {
    match result {
        Ok(doc) => {
            let mut out = std::io::stdout().lock();
            match out.write_all(&doc).and_then(|()| out.flush()) {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    eprintln!("og map {verb}: stdout: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        Err(BuildError::Usage(e)) => usage(&e),
        Err(BuildError::NoMap(e)) => {
            eprintln!("og map {verb}: {e}");
            ExitCode::from(3)
        }
        Err(BuildError::Internal(e)) => {
            eprintln!("og map {verb}: {e}");
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
                cache_dir: "/c".into(),
                predicate_table: false,
            }
        );
        let t = parse_build(&args(
            "--predicate-table --root /g --run R --assets /a --cache-dir /c",
        ))
        .unwrap();
        assert!(t.predicate_table);
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

    #[test]
    fn genes_arguments() {
        assert_eq!(
            parse_genes(&args("--root /g --cache-dir /c --chrom chr7")).unwrap(),
            Genes {
                root: "/g".into(),
                cache_dir: "/c".into(),
                chrom: Some("chr7".into())
            }
        );
        assert_eq!(
            parse_genes(&args("--names --root=/g --cache-dir=/c")).unwrap(),
            Genes {
                root: "/g".into(),
                cache_dir: "/c".into(),
                chrom: None
            }
        );
        assert_eq!(
            parse_genes(&args("--root /g --cache-dir /c")).unwrap_err(),
            "one of the arguments --chrom --names is required"
        );
        assert_eq!(
            parse_genes(&args("--root /g --cache-dir /c --names --chrom chr1")).unwrap_err(),
            "argument --names: not allowed with argument --chrom"
        );
        assert_eq!(
            parse_genes(&args("--names --cache-dir /c")).unwrap_err(),
            "the following arguments are required: --root"
        );
        assert_eq!(
            parse_genes(&args("--root /g --cache-dir /c --chrom")).unwrap_err(),
            "argument --chrom: expected one argument"
        );
    }
}
