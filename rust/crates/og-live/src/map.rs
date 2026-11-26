//! The genome map's files on the live site.
//!
//! ```text
//! GET /runs/<run>/map.json        `og map build` of a person's run: stdin is {"results": the
//!                                 body /runs/<run>/results.json sends}, the answer is og's
//!                                 stdout (gzip)
//! GET /map/genes/<chrom>.json     `og map genes --chrom <chrom>`: one chromosome's genes (gzip)
//! GET /map/genes/names.json       `og map genes --names`: every gene's name and place (gzip)
//! ```
//!
//! The site decides nothing about the map: og builds each file and caches it
//! under `STATE_DIR/map`, and the site sends og's bytes unchanged, with its
//! own headers. og writes a file only once it is complete and its exit code
//! sets the status, so the site waits for og to end: 0 is 200; 3 (no map for
//! this run, no gene files) is 404 with og's reason; 127 (og cannot run on
//! this host) is 503; anything else is 500. On 127 and 500, og's first stderr
//! line goes to the server's stderr. Each request writes one access-log row
//! `map`, with the run (or `genes` as the section for a gene file) and never
//! a finding, a gene or a chromosome.
//!
//! The map uses the run's live results, including the veil for a child's
//! closed results. With no configured results producer, it returns 503
//! before calling `og map build`.

use crate::server::Req;
use crate::state::dumps;
use pycompat::json::Value;
use std::io::{self, Read, Write};
use std::path::Path;
use std::process::{Command, Stdio};
use std::thread;
use std::time::{Duration, Instant};

/// How long og may take.
pub const OG_SECONDS: u64 = 300;
pub const GZIP: &str = "application/gzip";
const TEXT: &str = "text/plain; charset=utf-8";
/// The answer when og cannot run on this host.
pub const NO_OG: &str = "The map needs og on this host: run scripts/setup.sh --og\n";
/// The answer when og failed.
pub const FAILED: &str = "The map could not be built: the server's log says why.\n";
/// The answer while the site has no results to build a map from.
pub const NO_RESULTS: &str =
    "The map needs the results, which og live serve does not serve yet: the report still works.\n";

/// How og ended: its exit code, or why it has none.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Ended {
    Code(i32),
    /// Not found, or not runnable here.
    Missing,
    Timeout,
    /// Killed by a signal, or another error starting it.
    Other(String),
}

/// Runs `og ARGS` with `stdin`, up to [`OG_SECONDS`]: how it ended, its
/// stdout and its first stderr line.
pub fn run_og(
    og: &Path,
    args: &[String],
    stdin: &[u8],
    limit: Duration,
) -> (Ended, Vec<u8>, String) {
    let child = Command::new(og)
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn();
    let mut child = match child {
        Ok(c) => c,
        Err(e) if e.kind() == io::ErrorKind::NotFound => {
            return (Ended::Missing, Vec::new(), String::new());
        }
        Err(e) => {
            return (
                Ended::Other(e.kind().to_string()),
                Vec::new(),
                e.to_string(),
            );
        }
    };
    let mut pipe_in = child.stdin.take();
    let input = stdin.to_vec();
    let writer = thread::spawn(move || {
        if let Some(p) = pipe_in.as_mut() {
            let _ = p.write_all(&input);
        }
    });
    let mut so = child.stdout.take();
    let mut se = child.stderr.take();
    let out_t = thread::spawn(move || {
        let mut b = Vec::new();
        if let Some(p) = so.as_mut() {
            let _ = p.read_to_end(&mut b);
        }
        b
    });
    let err_t = thread::spawn(move || {
        let mut b = Vec::new();
        if let Some(p) = se.as_mut() {
            let _ = p.read_to_end(&mut b);
        }
        b
    });
    let end = Instant::now() + limit;
    let status = loop {
        match child.try_wait() {
            Ok(Some(s)) => break Some(s),
            Ok(None) if Instant::now() >= end => {
                let _ = child.kill();
                let _ = child.wait();
                break None;
            }
            Ok(None) => thread::sleep(Duration::from_millis(20)),
            Err(_) => break None,
        }
    };
    let _ = writer.join();
    let out = out_t.join().unwrap_or_default();
    let err = err_t.join().unwrap_or_default();
    let first = String::from_utf8_lossy(&err)
        .trim()
        .lines()
        .next()
        .unwrap_or("")
        .to_owned();
    let Some(status) = status else {
        return (
            Ended::Timeout,
            Vec::new(),
            format!("og did not end within {} s", limit.as_secs()),
        );
    };
    let ended = match status.code() {
        Some(c) => Ended::Code(c),
        None => {
            use std::os::unix::process::ExitStatusExt;
            Ended::Other(format!("-{}", status.signal().unwrap_or(0)))
        }
    };
    (ended, out, first)
}

/// Runs og with `args` and `stdin` and sends its answer; one access-log row.
pub fn relay(
    req: &mut Req<'_>,
    args: &[String],
    stdin: &[u8],
    run: &str,
    section: &str,
) -> io::Result<()> {
    let og = req.live.cfg.og.clone();
    let (ended, out, err) = run_og(&og, args, stdin, Duration::from_secs(OG_SECONDS));
    let verb = args.get(1).map_or("", String::as_str);
    match ended {
        Ended::Code(0) => {
            req.log("map", run, section, 200);
            req.send_body(200, GZIP, out)
        }
        Ended::Code(3) => {
            req.log("map", run, section, 404);
            req.send_body(404, TEXT, format!("{err}\n"))
        }
        Ended::Code(127) | Ended::Missing => {
            if !err.is_empty() {
                eprintln!("og live serve: og map {verb} unavailable: {err}");
            }
            req.log("map", run, section, 503);
            req.send_body(503, TEXT, NO_OG)
        }
        other => {
            let code = match other {
                Ended::Code(c) => c.to_string(),
                Ended::Timeout => "timeout".into(),
                Ended::Other(w) => w,
                Ended::Missing => unreachable!("answered above"),
            };
            eprintln!("og live serve: og map {verb} failed ({code}): {err}");
            req.log("map", run, section, 500);
            req.send_body(500, TEXT, FAILED)
        }
    }
}

/// Where og keeps the map files it built: `STATE_DIR/map`.
fn cache_dir(req: &Req<'_>) -> String {
    req.live
        .cfg
        .state_dir
        .join("map")
        .to_string_lossy()
        .into_owned()
}

/// `GET /runs/<run>/map.json`: og map build of a person's run.
pub fn build(req: &mut Req<'_>, name: &str) -> io::Result<()> {
    let live = req.live;
    let Some(run) = live.run(name).filter(|r| !r.is_family()) else {
        req.log("map", name, "-", 404);
        return req.send_body(404, TEXT, "not found\n");
    };
    let Some(fns) = live.cfg.results else {
        req.log("map", &run.name, "-", 503);
        return req.send_body(503, TEXT, NO_RESULTS);
    };
    let body = match crate::plugins::results_body(live, &run, &fns) {
        Ok(b) => b,
        Err(e) => {
            eprintln!(
                "og live serve: results of {} not built for the map: {e}",
                run.name
            );
            req.log("map", &run.name, "-", 500);
            return req.send_body(500, TEXT, FAILED);
        }
    };
    let mut o = indexmap::IndexMap::new();
    o.insert("results".to_owned(), body);
    let stdin = dumps(&Value::Object(o));
    let args: Vec<String> = [
        "map",
        "build",
        "--root",
        &live.root.to_string_lossy(),
        "--run",
        &run.name,
        "--assets",
        &live.cfg.repo_assets.to_string_lossy(),
        "--cache-dir",
        &cache_dir(req),
    ]
    .iter()
    .map(|s| (*s).to_owned())
    .collect();
    relay(req, &args, stdin.as_bytes(), &run.name, "-")
}

/// Whether `which` is a gene file the map has: `chr` and one or two letters
/// or digits, or `names`.
pub fn valid_genes(which: &str) -> bool {
    which == "names"
        || which.strip_prefix("chr").is_some_and(|c| {
            (1..=2).contains(&c.len()) && c.bytes().all(|b| b.is_ascii_alphanumeric())
        })
}

/// `GET /map/genes/<chrom>.json` and `/map/genes/names.json`: og map genes.
pub fn genes(req: &mut Req<'_>, which: &str) -> io::Result<()> {
    let mut args: Vec<String> = vec![
        "map".into(),
        "genes".into(),
        "--root".into(),
        req.live.root.to_string_lossy().into_owned(),
        "--cache-dir".into(),
        cache_dir(req),
    ];
    if which == "names" {
        args.push("--names".into());
    } else {
        args.push("--chrom".into());
        args.push(which.into());
    }
    relay(req, &args, b"", "-", "genes")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gene_files() {
        assert!(
            valid_genes("names")
                && valid_genes("chr1")
                && valid_genes("chrX")
                && valid_genes("chr22")
        );
        assert!(
            !valid_genes("chr")
                && !valid_genes("chr123")
                && !valid_genes("chr1.x")
                && !valid_genes("x1")
        );
    }

    #[test]
    fn og_missing_or_failing() {
        let (e, _, _) = run_og(
            Path::new("/nonexistent/og"),
            &[],
            b"",
            Duration::from_secs(5),
        );
        assert_eq!(e, Ended::Missing);
        let (e, out, err) = run_og(
            Path::new("/bin/sh"),
            &[
                "-c".into(),
                "cat; echo why >&2; echo more >&2; exit 3".into(),
            ],
            b"in",
            Duration::from_secs(5),
        );
        assert_eq!(
            (e, out.as_slice(), err.as_str()),
            (Ended::Code(3), &b"in"[..], "why")
        );
        let (e, _, err) = run_og(
            Path::new("/bin/sh"),
            &["-c".into(), "sleep 5".into()],
            b"",
            Duration::from_millis(200),
        );
        assert_eq!(e, Ended::Timeout);
        assert!(err.contains("did not end"));
    }
}
