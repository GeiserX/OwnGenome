//! `og dev build`: the builders of the asset tables OwnGenome commits.
//!
//! Each verb rebuilds one committed file from pinned public downloads that
//! a shell wrapper in `scripts/ci/` has already fetched and checked against
//! its sha256 (`scripts/ci/lib/fetch-pinned.sh`): `og` itself never
//! downloads. The committed file is the oracle, byte for byte.
//!
//! - `og dev build cn-intervals` ([`cn_intervals`]): the interval rows of
//!   `bin/cn_loci.tsv` (`scripts/ci/build-cn-intervals.sh`).
//! - `og dev build callability` ([`callability`]):
//!   `assets/callability_genes.tsv` (`scripts/ci/build-callability.sh`).
//! - `og dev build sv-rule-regions` ([`sv_rules`]):
//!   `assets/sv_rule_regions.tsv` (`scripts/ci/build-sv-rule-regions.sh`).
//!
//! Exit 0 when the file was written (or `--check` found it right); 1 when
//! `--check` found it different, or when the callability or SV rule
//! downloads hold something their rules refuse (a CMRG gene with no GENCODE
//! gene, a symbol listed twice, not one MHC region), as the Python builders
//! exited; 2 for a usage error or an input that cannot be read, and for
//! every cn-intervals input the rules refuse, as its Python builder exited.

pub mod callability;
pub mod cn_intervals;
mod gtf;
pub mod sv_rules;

use og_io::{Lines, Utf8};
use std::fmt;
use std::io::{self, BufRead, Write};
use std::path::Path;

const USAGE: &str = "usage: og dev build cn-intervals --clingen BED --segdups BED --gaps BED --satellites BED --parascopy TXT (--catalog FILE [--check] | --stdout)\n       og dev build callability --work DIR --repo DIR --out FILE\n       og dev build sv-rule-regions --work DIR --out FILE\n";

/// GRCh38's primary chromosomes in their sort order: chr1-22, X, Y.
pub(crate) const CHROMS: [&str; 24] = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
    "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
    "chr22", "chrX", "chrY",
];

/// The index of a primary chromosome in [`CHROMS`], from its name's bytes.
pub(crate) fn chrom_index(name: impl AsRef<[u8]>) -> Option<usize> {
    let rest = name.as_ref().strip_prefix(b"chr")?;
    match rest {
        b"X" => Some(22),
        b"Y" => Some(23),
        [d @ b'1'..=b'9'] => Some(usize::from(d - b'1')),
        [a @ b'1'..=b'2', b @ b'0'..=b'9'] => {
            let n = usize::from(a - b'0') * 10 + usize::from(b - b'0');
            (n <= 22).then(|| n - 1)
        }
        _ => None,
    }
}

/// Why a build stopped: exit 1 for downloads the rules refuse, exit 2 for
/// an input that cannot be read.
#[derive(Debug)]
pub enum BuildError {
    /// The inputs were read, but they hold what the rules refuse.
    Data(String),
    /// An input could not be read or parsed.
    Input(String),
}

impl BuildError {
    fn code(&self) -> u8 {
        match self {
            BuildError::Data(_) => 1,
            BuildError::Input(_) => 2,
        }
    }
}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BuildError::Data(m) | BuildError::Input(m) => f.write_str(m),
        }
    }
}

/// An input that cannot be read, named by its path.
pub(crate) fn input_error(path: &Path, e: &io::Error) -> BuildError {
    BuildError::Input(crate::py_os_error(e, &path.to_string_lossy()))
}

/// The lines of a plain or gzip file, without their endings, as Python's
/// text mode reads them.
pub(crate) fn lines(path: &Path, utf8: Utf8) -> Result<Lines<Box<dyn BufRead>>, BuildError> {
    let r = og_io::open(path).map_err(|e| input_error(path, &e))?;
    Ok(Lines::new(r, utf8))
}

/// Reads each line of `path` into `f`, stopping at the first error.
pub(crate) fn each_line(
    path: &Path,
    utf8: Utf8,
    mut f: impl FnMut(&str) -> Result<(), BuildError>,
) -> Result<(), BuildError> {
    for line in lines(path, utf8)? {
        f(&line.map_err(|e| input_error(path, &e))?)?;
    }
    Ok(())
}

/// Reads each line of `path` (plain or gzip) into `f` as bytes, without
/// its `\n` or `\r\n`, stopping at the first error. For the big files,
/// where a text line per line would cost more than the reading.
pub(crate) fn each_raw_line(
    path: &Path,
    mut f: impl FnMut(&[u8]) -> Result<(), BuildError>,
) -> Result<(), BuildError> {
    let mut r = og_io::open(path).map_err(|e| input_error(path, &e))?;
    let mut buf = Vec::with_capacity(1 << 12);
    loop {
        buf.clear();
        if r.read_until(b'\n', &mut buf)
            .map_err(|e| input_error(path, &e))?
            == 0
        {
            return Ok(());
        }
        let mut line = buf.as_slice();
        if let Some(l) = line.strip_suffix(b"\n") {
            line = l.strip_suffix(b"\r").unwrap_or(l);
        }
        f(line)?;
    }
}

/// A whole number of a text column, as Python's `int()` reads one.
pub(crate) fn int(cell: &str, what: &str) -> Result<i64, BuildError> {
    pycompat::text::strip(cell, None)
        .parse()
        .map_err(|_| BuildError::Input(format!("{what}: not a whole number: {cell:?}")))
}

/// Takes `--flag VALUE` pairs and bare `--switch`es from `args`.
pub(crate) struct Args {
    values: Vec<(String, String)>,
    switches: Vec<String>,
}

impl Args {
    pub(crate) fn parse(
        args: &[String],
        flags: &[&str],
        switches: &[&str],
    ) -> Result<Args, String> {
        let mut out = Args {
            values: Vec::new(),
            switches: Vec::new(),
        };
        let mut it = args.iter();
        while let Some(a) = it.next() {
            if flags.contains(&a.as_str()) {
                let v = it
                    .next()
                    .ok_or_else(|| format!("argument {a}: expected one argument"))?;
                if out.values.iter().any(|(k, _)| k == a) {
                    return Err(format!("argument {a}: given twice"));
                }
                out.values.push((a.clone(), v.clone()));
            } else if switches.contains(&a.as_str()) {
                out.switches.push(a.clone());
            } else {
                return Err(format!("unrecognized arguments: {a}"));
            }
        }
        Ok(out)
    }

    pub(crate) fn get(&self, flag: &str) -> Option<&str> {
        self.values
            .iter()
            .find(|(k, _)| k == flag)
            .map(|(_, v)| v.as_str())
    }

    pub(crate) fn need(&self, flag: &str) -> Result<&str, String> {
        self.get(flag)
            .ok_or_else(|| format!("the following arguments are required: {flag}"))
    }

    pub(crate) fn has(&self, switch: &str) -> bool {
        self.switches.iter().any(|s| s == switch)
    }
}

/// `og dev build VERB ARGS`; returns the exit code.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let result = match args.first().map(String::as_str) {
        Some("cn-intervals") => cn_intervals::cli(&args[1..], out, err),
        Some("callability") => callability::cli(&args[1..], err),
        Some("sv-rule-regions") => sv_rules::cli(&args[1..], err),
        Some("-h" | "--help") => {
            let _ = out.write_all(USAGE.as_bytes());
            return 0;
        }
        Some(other) => Err(Failure::Usage(format!("unknown builder '{other}'"))),
        None => Err(Failure::Usage(
            "the following arguments are required: builder".to_owned(),
        )),
    };
    match result {
        Ok(code) => code,
        Err(Failure::Usage(m)) => {
            let _ = writeln!(err, "{USAGE}og dev build: error: {m}");
            2
        }
        Err(Failure::Build(e)) => {
            let _ = writeln!(err, "ERROR: {e}");
            e.code()
        }
    }
}

/// How a builder's command line ended badly.
pub(crate) enum Failure {
    Usage(String),
    Build(BuildError),
}

impl From<BuildError> for Failure {
    fn from(e: BuildError) -> Self {
        Failure::Build(e)
    }
}

impl From<String> for Failure {
    fn from(m: String) -> Self {
        Failure::Usage(m)
    }
}

/// Writes `bytes` to `path` through a temporary file renamed over it.
pub(crate) fn write_atomic(path: &Path, bytes: &[u8]) -> Result<(), BuildError> {
    let mut f = og_io::AtomicFile::create(path).map_err(|e| input_error(path, &e))?;
    f.write_all(bytes).map_err(|e| input_error(path, &e))?;
    f.commit().map_err(|e| input_error(path, &e))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(v: &[&str]) -> Vec<String> {
        v.iter().map(|x| (*x).to_owned()).collect()
    }

    #[test]
    fn unknown_builder_is_a_usage_error() {
        let (mut o, mut e) = (Vec::new(), Vec::new());
        assert_eq!(cli(&s(&["nope"]), &mut o, &mut e), 2);
        assert!(
            String::from_utf8(e)
                .unwrap()
                .contains("unknown builder 'nope'")
        );
        assert_eq!(cli(&[], &mut o, &mut Vec::new()), 2);
    }

    #[test]
    fn args_refuse_a_repeated_or_unknown_flag() {
        assert!(Args::parse(&s(&["--out", "a", "--out", "b"]), &["--out"], &[]).is_err());
        assert!(Args::parse(&s(&["--what"]), &["--out"], &[]).is_err());
        assert!(Args::parse(&s(&["--out"]), &["--out"], &[]).is_err());
        let a = Args::parse(&s(&["--out", "a", "--check"]), &["--out"], &["--check"]).unwrap();
        assert_eq!(a.get("--out"), Some("a"));
        assert!(a.has("--check"));
    }

    #[test]
    fn chromosome_names() {
        for (i, c) in CHROMS.iter().enumerate() {
            assert_eq!(chrom_index(c), Some(i));
        }
        for c in [
            "chr0",
            "chr23",
            "chr01",
            "chrM",
            "chr1_KI270706v1_random",
            "1",
            "chrx",
            "",
        ] {
            assert_eq!(chrom_index(c), None, "{c}");
        }
    }

    #[test]
    fn int_reads_as_python_does() {
        assert_eq!(int(" 12 ", "x").unwrap(), 12);
        assert!(int("1.5", "x").is_err());
    }
}

#[cfg(all(test, unix))]
mod fetch_tests {
    use std::io::{BufRead, BufReader, Write};
    use std::net::TcpListener;
    use std::process::Command;
    use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
    use std::time::Duration;

    static NEXT: AtomicU64 = AtomicU64::new(0);
    const BODY: &str = "hello\n";
    const SHA: &str = "5891b5b522d5df086d0ff0b110fbd9d21bb4fc7163af34d08286a2e846f6be03";

    // Real HTTP responses exercise curl's own retries too. Shorten only the
    // test delays; the production script still supplies its stall limits.
    fn fetch(statuses: &[u16], sha: &str, cached: bool, succeeds: bool, expected: &[u16]) {
        let dir = std::env::temp_dir().join(format!(
            "og-fetch-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).unwrap();
        if cached {
            std::fs::write(dir.join("asset"), BODY).unwrap();
        }
        let server = TcpListener::bind("127.0.0.1:0").unwrap();
        server.set_nonblocking(true).unwrap();
        let url = format!("http://{}/asset", server.local_addr().unwrap());
        let stop = AtomicBool::new(false);
        let (output, seen) = std::thread::scope(|scope| {
            let worker = scope.spawn(|| {
                let mut seen = Vec::new();
                while !stop.load(Ordering::Relaxed) {
                    match server.accept() {
                        Ok((mut stream, _)) => {
                            stream.set_read_timeout(Some(Duration::from_secs(5))).unwrap();
                            let mut reader = BufReader::new(&stream);
                            let mut line = String::new();
                            loop {
                                line.clear();
                                if reader.read_line(&mut line).unwrap() == 0 || line == "\r\n" {
                                    break;
                                }
                            }
                            let code = statuses.get(seen.len()).copied().unwrap_or(200);
                            seen.push(code);
                            write!(stream,
                                "HTTP/1.1 {code} Test\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{BODY}",
                                BODY.len()).unwrap();
                        }
                        Err(e) if e.kind() == std::io::ErrorKind::WouldBlock => {
                            std::thread::sleep(Duration::from_millis(5));
                        }
                        Err(e) => panic!("HTTP test server: {e}"),
                    }
                }
                seen
            });
            let output = Command::new("bash")
                .args([
                    "-c",
                    r#"
                    source "$1"
                    sleep() { :; }
                    curl() { command curl "$@" --noproxy '*' --retry-delay 1 --max-time 5; }
                    fetch_pinned "$2" asset "$3" "$4"
                "#,
                    "_",
                ])
                .arg(concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/../../../scripts/ci/lib/fetch-pinned.sh"
                ))
                .arg(&dir)
                .arg(url)
                .arg(sha)
                .output()
                .unwrap();
            stop.store(true, Ordering::Relaxed);
            (output, worker.join().unwrap())
        });
        let file = std::fs::read_to_string(dir.join("asset")).ok();
        let files: Vec<_> = std::fs::read_dir(&dir)
            .unwrap()
            .map(|e| e.unwrap().file_name())
            .collect();
        std::fs::remove_dir_all(&dir).unwrap();
        assert_eq!(seen, expected, "HTTP sequence {statuses:?}");
        assert_eq!(
            output.status.success(),
            succeeds,
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(file.as_deref(), succeeds.then_some(BODY));
        assert_eq!(
            files.len(),
            usize::from(succeeds),
            "temporary download leaked: {files:?}"
        );
        if let Some(code @ (403 | 404)) = expected.last() {
            assert!(String::from_utf8_lossy(&output.stderr).contains(&format!("HTTP {code}")));
        }
    }

    #[test]
    fn http_403_stops_after_any_attempt() {
        fetch(&[403], SHA, false, false, &[403]);
        fetch(&[503, 403], SHA, false, false, &[503, 403]);
    }

    #[test]
    fn http_404_stops_after_any_attempt() {
        fetch(&[404], SHA, false, false, &[404]);
        fetch(&[503, 404], SHA, false, false, &[503, 404]);
    }

    #[test]
    fn other_failures_retry_but_stop_at_the_bound() {
        fetch(&[503, 503, 200], SHA, false, true, &[503, 503, 200]);
        fetch(&[503; 7], SHA, false, false, &[503; 7]);
    }

    #[test]
    fn downloads_are_pinned_and_cached() {
        fetch(&[200], SHA, false, true, &[200]);
        fetch(&[200], "wrong-sha256", false, false, &[200]);
        fetch(&[], SHA, true, true, &[]);
    }
}
