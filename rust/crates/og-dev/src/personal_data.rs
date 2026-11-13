//! `og dev scan personal-data`: fail when a tracked text file holds
//! something that looks like personal data. It replaces the Python inside
//! `scripts/ci/check-personal-data.sh`, which keeps the lists of documented
//! examples (`--allow`) and skipped files (`--skip`) and passes them in.
//!
//! ```text
//! og dev scan personal-data --root DIR --allow TEXT --skip TEXT
//!                           [--patterns FILE] [--cached] [--] [FILE...]
//! ```
//!
//! With no FILE, every file `git ls-files` lists under DIR. `--cached`
//! reads each file's staged blob (one `git cat-file --batch` process)
//! instead of the working copy. `--patterns` adds one regular expression
//! per line (kind `own`).
//!
//! Each line of a text file (no NUL byte in its first 8 KiB) is matched
//! against the rules below; a match inside an allowed example of that file
//! is skipped, the rest of the line still counts. Output: one
//! `FAIL <kind> <escaped-file>:<line>` per rule and line, then
//! `scanned N text files`. Exit 0 clean, 1 on a finding or when no text
//! file was scanned, 2 on a usage error or a file or pattern that cannot be
//! read.
//!
//! Where it differs from the Python on purpose: a scan that read no text
//! file fails (`FAIL: no text file was scanned`) instead of passing, and a
//! file is placed under DIR after resolving symbolic links in both, so
//! DIR given as `/tmp/x` still finds `/private/tmp/x/CLAUDE.md`'s examples.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::LazyLock;

use pycompat::mtime::fnmatch;
use pycompat::num::py_strip;

use crate::lookaround::PythonRule as Rule;

const USAGE: &str = "usage: og dev scan personal-data --root DIR --allow TEXT --skip TEXT [--patterns FILE] [--cached] [--] [FILE...]\n";

/// This file's own path: it holds every allowed text, so all of them apply.
const SELF: &str = "scripts/ci/check-personal-data.sh";

/// How many leading bytes are searched for a NUL to call a file binary.
const BINARY_PROBE: usize = 8192;

const STR_GENES: &str = "HTT|FMR1|C9orf72|C9ORF72|ATXN[0-9]+|DMPK|FXN|CNBP|RFC1|PABPN1|CACNA1A|TCF4|AR|PPP2R2B|TBP|JPH3|NOP56|ATN1|CSTB|AFF2|FMR2|NOTCH2NLC|DAB1|BEAN1|STARD7|YEATS2";
const OCTET: &str = "(?:25[0-5]|2[0-4][0-9]|1[0-9][0-9]|[1-9]?[0-9])";

/// Python has three distinct line operations here. Pattern files use text
/// mode universal newlines; exception tables use str.splitlines(); scanned
/// blobs are decoded from bytes and split on LF, preserving the last empty
/// line. Keep the distinction in this one helper, including CRLF as one row.
#[derive(Clone, Copy)]
enum Lines {
    Blob,
    TextFile,
    ExceptionTable,
}

fn lines(text: &str, mode: Lines) -> impl Iterator<Item = &str> {
    let (mut rest, mut finished) = (text, false);
    std::iter::from_fn(move || {
        if finished || (rest.is_empty() && !matches!(mode, Lines::Blob)) {
            return None;
        }
        let at = rest.char_indices().find(|&(_, c)| match mode {
            Lines::Blob => c == '\n',
            Lines::TextFile => matches!(c, '\n' | '\r'),
            Lines::ExceptionTable => matches!(
                c,
                '\n' | '\r'
                    | '\u{b}'
                    | '\u{c}'
                    | '\u{1c}'
                    | '\u{1d}'
                    | '\u{1e}'
                    | '\u{85}'
                    | '\u{2028}'
                    | '\u{2029}'
            ),
        });
        match at {
            Some((i, c)) => {
                let line = &rest[..i];
                let mut end = i + c.len_utf8();
                if c == '\r' && rest[end..].starts_with('\n') {
                    end += 1;
                }
                rest = &rest[end..];
                Some(line)
            }
            None => {
                finished = true;
                Some(rest)
            }
        }
    })
}

/// A rule: its kind, its pattern, and a pattern the line must also match.
#[derive(Clone)]
struct Check {
    kind: &'static str,
    rule: Rule,
    cond: Option<&'static [Rule]>,
}

impl Check {
    /// Necessary literal bytes only. This can admit extra lines, but cannot
    /// decide a match or discard a candidate span. Custom rules always run.
    fn candidate(&self, line: &str) -> bool {
        match self.kind {
            "home" => {
                line.contains("/home/") || line.contains("/mnt/use") || line.contains("Users")
            }
            "address" => line.contains('.'),
            "diplotype" => line.contains('*') && line.contains('/'),
            "rsid" => line.contains("rs"),
            "repeat" => line.contains('/'),
            _ => true,
        }
    }
}

fn rule(pattern: &str) -> Rule {
    Rule::compile(pattern, true).expect("a valid built-in pattern")
}

static LABEL: LazyLock<[Rule; 1]> = LazyLock::new(|| {
    [rule(
        r"(?:(?i:(?<![A-Za-z])(?:sample|patient|subject|proband|participant|individual|donor)s?(?:[_-]?[0-9]+)?(?![A-Za-z]))|\b(?:NA|HG)[0-9]{3,5}\b)",
    )]
});
static REPCN: LazyLock<[Rule; 1]> = LazyLock::new(|| [rule(r"\bREPCN\b")]);

/// Built-ins use the same compiler and iterator as custom patterns.
fn builtin() -> Vec<Check> {
    static RULES: LazyLock<Vec<Check>> = LazyLock::new(|| {
        let gt = r"\(?[ACGTDI][/|;]?[ACGTDI]\)?(?![A-Za-z0-9])";
        vec![
            Check {
                kind: "home",
                rule: rule(
                    r"(?<![A-Za-z0-9_])(?:/home/[A-Za-z]|/Users/[A-Za-z]|/mnt/use[r]\b|[A-Za-z]:\\\\?Users\\\\?[A-Za-z])",
                ),
                cond: None,
            },
            Check {
                kind: "address",
                rule: rule(&format!(
                    r"(?<![0-9.])(?:10\.{OCTET}|172\.(?:1[6-9]|2[0-9]|3[01])|192\.168|100\.(?:6[4-9]|[7-9][0-9]|1[01][0-9]|12[0-7]))\.{OCTET}\.{OCTET}(?![0-9]|\.[0-9])"
                )),
                cond: None,
            },
            Check {
                kind: "address",
                rule: rule(r"[A-Za-z0-9-]\.ts\.net\b"),
                cond: None,
            },
            Check {
                kind: "diplotype",
                rule: rule(r"\*[0-9]+[A-Za-z0-9.]*/\*[0-9]+"),
                cond: Some(&*LABEL),
            },
            Check {
                kind: "rsid",
                rule: rule(&format!(r#"\brs[0-9]+\b[\s:=,;|`'"\[\]-]{{1,6}}{gt}"#)),
                cond: None,
            },
            Check {
                kind: "rsid",
                rule: rule(r"\brs[0-9]+\t[0-9XYMT]{1,2}\t[0-9]+\t[ACGTDI-]{1,2}\b"),
                cond: None,
            },
            Check {
                kind: "repeat",
                rule: rule(r"(?<![0-9./])[0-9]{1,4}/[0-9]{1,4}(?![0-9/])"),
                cond: Some(&*REPCN),
            },
            Check {
                kind: "repeat",
                rule: rule(&format!(
                    r"\b(?:{STR_GENES})\b[^/\n]{{0,40}}?(?<![0-9.])[0-9]{{1,4}}/[0-9]{{1,4}}(?![0-9/])"
                )),
                cond: None,
            },
        ]
    });
    RULES.clone()
}

/// The arguments.
struct Args {
    root: PathBuf,
    allow: Vec<(String, String)>,
    skip: Vec<String>,
    patterns: Option<String>,
    cached: bool,
    files: Vec<PathBuf>,
}

fn parse(args: &[String]) -> Result<Args, String> {
    let (mut root, mut allow, mut skip, mut patterns, mut cached) = (None, None, None, None, false);
    let mut files = Vec::new();
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let mut value = |name: &str| {
            it.next()
                .cloned()
                .ok_or_else(|| format!("{name} needs a value"))
        };
        match a.as_str() {
            "--root" => root = Some(value("--root")?),
            "--allow" => allow = Some(value("--allow")?),
            "--skip" => skip = Some(value("--skip")?),
            "--patterns" => patterns = Some(value("--patterns")?),
            "--cached" => cached = true,
            "--" => {
                files.extend(it.by_ref().map(PathBuf::from));
                break;
            }
            s if s.starts_with('-') => return Err(format!("unknown option {s}")),
            s => files.push(PathBuf::from(s)),
        }
    }
    let root = root.ok_or("--root is required")?;
    let allow_text = allow.ok_or("--allow is required")?;
    let allow = lines(&allow_text, Lines::ExceptionTable)
        .filter(|l| !py_strip(l).is_empty())
        .map(|l| {
            let l = py_strip(l);
            let (path, text) = l.split_once(' ').unwrap_or((l, ""));
            (path.to_owned(), py_strip(text).to_owned())
        })
        .collect();
    let skip_text = skip.ok_or("--skip is required")?;
    let skip = lines(&skip_text, Lines::ExceptionTable)
        .filter_map(|l| {
            pycompat::text::split(l, None, None)
                .first()
                .map(|s| (*s).to_owned())
        })
        .collect();
    Ok(Args {
        root: PathBuf::from(root),
        allow,
        skip,
        patterns: patterns.filter(|p| !p.is_empty()),
        cached,
        files,
    })
}

/// The user's own patterns, one per line; `#` lines and blank lines skipped.
fn own_rules(path: &str) -> Result<Vec<Check>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    let mut out = Vec::new();
    for (line, l) in lines(&text, Lines::TextFile).enumerate() {
        let p = py_strip(l);
        if p.is_empty() || l.starts_with('#') {
            continue;
        }
        // Rule reports capture 1 for the built-in look-around patterns.
        // A user's captures must not narrow the span checked against ALLOW.
        let rule = Rule::compile(p, false).map_err(|_| {
            format!(
                "{}:{}: invalid or unsupported custom pattern",
                path.escape_debug(),
                line + 1
            )
        })?;
        out.push(Check {
            kind: "own",
            rule,
            cond: None,
        });
    }
    Ok(out)
}

/// `path` relative to `root` when it lies below it (both with symbolic
/// links resolved), else `path` as given.
fn relative(path: &Path, root: &Path) -> Result<String, String> {
    utf8_path(path)?;
    let abs = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir().map_or_else(|_| path.to_path_buf(), |d| d.join(path))
    };
    // The folder is resolved, not the file: a tracked symbolic link keeps
    // its own name.
    let placed = match (abs.parent(), abs.file_name()) {
        (Some(dir), Some(name)) => dir
            .canonicalize()
            .map_or_else(|_| abs.clone(), |d| d.join(name)),
        _ => abs.clone(),
    };
    match placed.strip_prefix(root) {
        Ok(rel) if !rel.as_os_str().is_empty() => utf8_path(rel).map(str::to_owned),
        _ => utf8_path(path).map(str::to_owned),
    }
}

fn utf8_path(path: &Path) -> Result<&str, String> {
    path.to_str()
        .ok_or_else(|| "path is not valid UTF-8".into())
}

fn git_files(root: &Path) -> Result<Vec<PathBuf>, String> {
    let out = Command::new("git")
        .arg("-C")
        .arg(root)
        .args(["ls-files", "-z"])
        .stderr(Stdio::piped())
        .output()
        .map_err(|e| format!("cannot run git: {e}"))?;
    if !out.status.success() {
        return Err(format!("git ls-files failed in {}", root.display()));
    }
    out.stdout
        .split(|&b| b == 0)
        .filter(|f| !f.is_empty())
        .map(|f| {
            std::str::from_utf8(f)
                .map(|s| root.join(s))
                .map_err(|_| "path is not valid UTF-8".into())
        })
        .collect::<Result<Vec<_>, String>>()
}

/// The staged blob of each path (relative to `root`), from one
/// `git cat-file --batch`. A path with no blob in the index is absent.
fn staged(root: &Path, rels: &[String]) -> Result<HashMap<String, Vec<u8>>, String> {
    let wanted: Vec<&String> = rels.iter().collect();
    let mut blobs = HashMap::new();
    if wanted.is_empty() {
        return Ok(blobs);
    }
    let mut child = Command::new("git")
        .arg("-C")
        .arg(root)
        .args(["cat-file", "--batch", "-Z"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .map_err(|e| format!("cannot run git: {e}"))?;
    let mut stdin = child.stdin.take().ok_or("no stdin for git")?;
    let input: Vec<u8> = wanted
        .iter()
        .flat_map(|r| format!(":{r}\0").into_bytes())
        .collect();
    let writer = std::thread::spawn(move || stdin.write_all(&input));
    let mut out = BufReader::new(child.stdout.take().ok_or("no stdout for git")?);
    for rel in &wanted {
        let mut header = Vec::new();
        if out.read_until(0, &mut header).map_err(|e| e.to_string())? == 0 {
            return Err("git cat-file stopped early".into());
        }
        if header.pop() != Some(0) {
            return Err("git cat-file: incomplete header".into());
        }
        let head =
            std::str::from_utf8(&header).map_err(|_| "git cat-file: invalid header encoding")?;
        let fields: Vec<&str> = head.split(' ').collect();
        let [_, kind, size] = fields[..] else {
            continue; // "<name> missing" or "ambiguous"
        };
        let size: usize = size
            .parse()
            .map_err(|_| "git cat-file: invalid object size".to_owned())?;
        let mut data = vec![0; size + 1];
        out.read_exact(&mut data).map_err(|e| e.to_string())?;
        data.pop();
        if kind == "blob" {
            blobs.insert((*rel).clone(), data);
        }
    }
    writer
        .join()
        .map_err(|_| "git input writer stopped")?
        .map_err(|e| e.to_string())?;
    if !child.wait().map_err(|e| e.to_string())?.success() {
        return Err("git cat-file failed".into());
    }
    Ok(blobs)
}

/// Reports metadata only. Escape control characters in paths so each
/// finding occupies one line, including filenames containing newlines.
fn report(out: &mut dyn Write, kind: &str, rel: &str, n: usize) {
    let _ = writeln!(out, "FAIL {kind:<9} {}:{n}", rel.escape_debug());
}

/// Scans one file's text; true when something was reported.
fn scan_text(
    out: &mut dyn Write,
    rules: &[Check],
    texts: &[&str],
    rel: &str,
    data: &[u8],
) -> Result<bool, ()> {
    let mut bad = false;
    // Intentional content decoding, matching Python errors="replace". Paths
    // and Git protocol headers are validated strictly before any conversion.
    let text = String::from_utf8_lossy(data);
    for (i, line) in lines(&text, Lines::Blob).enumerate() {
        let spans: Vec<(usize, usize)> = texts
            .iter()
            .flat_map(|t| line.match_indices(*t).map(|(a, s)| (a, a + s.len())))
            .collect();
        for c in rules {
            if !c.candidate(line) {
                continue;
            }
            if let Some(conditions) = c.cond {
                let mut matched = false;
                for condition in conditions {
                    matched |= condition.is_match(line)?;
                }
                if !matched {
                    continue;
                }
            }
            for found in c.rule.find_iter(line) {
                let (s, e) = found?;
                if !spans.iter().any(|&(a, b)| a <= s && e <= b) {
                    report(out, c.kind, rel, i + 1);
                    bad = true;
                    break;
                }
            }
        }
    }
    Ok(bad)
}

fn run(a: &Args, out: &mut dyn Write) -> Result<bool, String> {
    utf8_path(&a.root)?;
    let root = a
        .root
        .canonicalize()
        .map_err(|e| format!("--root {}: {e}", a.root.display()))?;
    utf8_path(&root)?;
    let mut rules = builtin();
    if let Some(p) = &a.patterns {
        rules.extend(own_rules(p)?);
    }
    let files = if a.files.is_empty() {
        git_files(&root)?
    } else {
        a.files.clone()
    };
    if files.is_empty() {
        let _ = writeln!(out, "FAIL: no file to scan");
        return Ok(false);
    }
    let todo: Vec<(&PathBuf, String)> = files
        .iter()
        .map(|f| relative(f, &root).map(|rel| (f, rel)))
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .filter(|(_, rel)| !a.skip.iter().any(|g| fnmatch(rel, g)))
        .collect();
    let blobs = if a.cached {
        staged(
            &root,
            &todo.iter().map(|(_, r)| r.clone()).collect::<Vec<_>>(),
        )?
    } else {
        HashMap::new()
    };
    let (mut bad, mut scanned) = (false, 0usize);
    for (path, rel) in &todo {
        let data = if a.cached {
            match blobs.get(rel) {
                Some(d) => std::borrow::Cow::Borrowed(d.as_slice()),
                None => continue,
            }
        } else {
            match std::fs::read(path) {
                Ok(d) => std::borrow::Cow::Owned(d),
                Err(e)
                    if matches!(
                        e.kind(),
                        std::io::ErrorKind::NotFound | std::io::ErrorKind::IsADirectory
                    ) =>
                {
                    continue;
                }
                Err(e) => return Err(format!("cannot read {}: {e}", path.display())),
            }
        };
        if data[..data.len().min(BINARY_PROBE)].contains(&0) {
            continue;
        }
        scanned += 1;
        let texts: Vec<&str> = a
            .allow
            .iter()
            .filter(|(p, _)| rel == SELF || fnmatch(rel, p))
            .map(|(_, t)| t.as_str())
            .collect();
        bad |= scan_text(out, &rules, &texts, rel, &data).map_err(|()| {
            format!(
                "{}: pattern matching could not complete",
                rel.escape_debug()
            )
        })?;
    }
    let _ = writeln!(out, "scanned {scanned} text files");
    if scanned == 0 {
        let _ = writeln!(out, "FAIL: no text file was scanned");
        return Ok(false);
    }
    Ok(!bad)
}

/// `og dev scan personal-data ARGS`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let a = match parse(args) {
        Ok(a) => a,
        Err(e) => {
            let _ = write!(err, "{USAGE}");
            let _ = writeln!(err, "og dev scan personal-data: error: {e}");
            return 2;
        }
    };
    match run(&a, out) {
        Ok(true) => 0,
        Ok(false) => 1,
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            2
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::OsStr;
    use std::os::unix::ffi::OsStrExt;

    fn kinds(line: &str) -> Vec<&'static str> {
        let rules = builtin();
        let mut out = Vec::new();
        for c in &rules {
            if c.cond
                .is_some_and(|cs| !cs.iter().any(|r| r.is_match(line).unwrap()))
            {
                continue;
            }
            if c.rule.is_match(line).unwrap() {
                out.push(c.kind);
            }
        }
        out
    }

    // Every planted value is split by "~", which is dropped, so this file
    // holds none of them whole.
    fn p(s: &str) -> String {
        s.replace('~', "")
    }

    #[test]
    fn each_kind_is_caught() {
        for (line, kind) in [
            ("cp results.html /home/~alice/genome/", "home"),
            ("C:\\~Users\\bob", "home"),
            ("scp out.vcf.gz root@192.~168.7.20:/data/", "address"),
            ("ssh nas.tail1234.ts.~net", "address"),
            ("| sample7 | CYP2D6 | *1/~*4 |", "diplotype"),
            ("sample_01 CYP2D6 *1/~*4", "diplotype"),
            ("HG0~02 *1/~*2", "diplotype"),
            ("rs4244285~ AG", "rsid"),
            ("rs123: (C;~T)", "rsid"),
            ("rs4477212\t1\t82154\t~AA", "rsid"),
            ("HTT CAG repeats: 17/~19", "repeat"),
            ("GT:REPCN 0~/1:17/~19", "repeat"),
        ] {
            let line = p(line);
            assert!(kinds(&line).contains(&kind), "{line}: {:?}", kinds(&line));
        }
    }

    #[test]
    fn edges_as_python_reads_them() {
        for line in [
            "- `*1/*4` -- Intermediate metabolizer",
            "samtools 1.21, GATK 4.6.2.0",
            "https://example.org/home/page",
            "CYP2C19*2 is rs4244285 (G>A)",
            // An address with a fifth part, or a version that looks like one.
            "v10.1.2.3.4",
            "1.10.1.2.3",
            "10.1.2.256",
            // The word must stand alone: "samples" is a label, "sampler" is not.
            "sampler *1/*4",
            // A repeat count with another digit or slash on either side.
            "REPCN 1.17/19",
            "REPCN 17/19/3",
            "/home/1x",
            "my/home/user",
        ] {
            assert!(kinds(line).is_empty(), "{line}: {:?}", kinds(line));
        }
        // Python takes the label "sample" out of "sample7a" by backing off
        // the number: the line still counts as labelled.
        assert_eq!(kinds(&p("sample7a *1/~*4")), ["diplotype"]);
        assert_eq!(kinds(&p("samples: *1/~*4")), ["diplotype"]);
        assert_eq!(kinds(&p("10.1.2.~3.")), ["address"]);
        assert_eq!(kinds(&p("FMR1 has 29/~31")), ["repeat"]);
    }

    #[test]
    fn an_allowed_example_covers_only_its_own_text() {
        let rules = builtin();
        let ex = p("(`/mnt/u~ser/`, `/ho~me/username/`, etc.)");
        let mut out = Vec::new();
        assert!(!scan_text(&mut out, &rules, &[ex.as_str()], "x.md", ex.as_bytes()).unwrap());
        let line = p(&format!("{ex} see /home/~bob/x"));
        assert!(scan_text(&mut out, &rules, &[ex.as_str()], "x.md", line.as_bytes()).unwrap());
        let got = String::from_utf8(out).unwrap();
        assert!(got.starts_with("FAIL home      x.md:1"), "{got}");
    }

    #[test]
    fn own_patterns_use_the_complete_match() {
        let path =
            std::env::temp_dir().join(format!("og-personal-patterns-{}", std::process::id()));
        std::fs::write(
            &path,
            "PREFIX (PUBLIC) SUFFIX\nPREFIX (?P<label>PUBLIC) SUFFIX\n(?i)prefix (public) suffix\nPREFIX (?:PUBLIC) SUFFIX\n(?x)PREFIX\\s+(PUBLIC)\\s+SUFFIX # trailing comment\n",
        )
        .unwrap();
        let rules = own_rules(path.to_str().unwrap()).unwrap();
        std::fs::remove_file(path).unwrap();

        let line = b"PREFIX PUBLIC SUFFIX";
        let mut out = Vec::new();
        assert!(scan_text(&mut out, &rules, &["PUBLIC"], "x.md", line).unwrap());
        let got = String::from_utf8(out).unwrap();
        assert_eq!(got.lines().count(), 5, "{got}");
        assert!(
            !scan_text(
                &mut Vec::new(),
                &rules,
                &["PREFIX PUBLIC SUFFIX"],
                "x.md",
                line,
            )
            .unwrap()
        );
    }

    struct Fixture(PathBuf);
    impl Fixture {
        fn new(name: &str) -> Self {
            let root = std::env::temp_dir().join(format!("og-scan-{name}-{}", std::process::id()));
            std::fs::create_dir_all(&root).unwrap();
            Self(root.canonicalize().unwrap())
        }
        fn args(&self, files: &[&str]) -> Vec<String> {
            let mut args = vec![
                "--root".into(),
                self.0.display().to_string(),
                "--allow".into(),
                "".into(),
                "--skip".into(),
                "".into(),
            ];
            args.extend(files.iter().map(|f| self.0.join(f).display().to_string()));
            args
        }
    }
    impl Drop for Fixture {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    // Execute the removed scanner itself. Only its reporting is adapted to
    // structured metadata, because printing findings is the bug being fixed.
    fn oracle(f: &Fixture, files: &[&str], cached: bool) -> String {
        oracle_options(f, files, cached, "", "", "")
    }

    fn oracle_options(
        f: &Fixture,
        files: &[&str],
        cached: bool,
        allow: &str,
        skip: &str,
        patterns: &str,
    ) -> String {
        static PYTHON: LazyLock<String> = LazyLock::new(|| {
            const SHA: &str = "4104a93f2121b457413336574b8563fe85291fd7";
            let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..");
            let read = || {
                Command::new("git")
                    .arg("-C")
                    .arg(&repo)
                    .args(["show", &format!("{SHA}:scripts/ci/check-personal-data.sh")])
                    .output()
                    .unwrap()
            };
            let mut source = read();
            if !source.status.success() {
                assert!(
                    Command::new("git")
                        .arg("-C")
                        .arg(&repo)
                        .args(["fetch", "--quiet", "--depth=1", "origin", SHA])
                        .status()
                        .unwrap()
                        .success()
                );
                source = read();
            }
            assert!(source.status.success());
            let source = String::from_utf8(source.stdout).unwrap();
            let python = source
                .split("<<'PY'\n")
                .nth(1)
                .unwrap()
                .split("\nPY\n")
                .next()
                .unwrap();
            python
                .replace(
                    "import fnmatch, os, re, subprocess, sys",
                    "import fnmatch, os, re, subprocess, sys, json",
                )
                .replace(
                    "print(\"FAIL %-9s %s:%d: %s\" % (kind, rel, n, line.strip()[:160]))",
                    "print(json.dumps([kind, rel, n]))",
                )
        });
        let mut cmd = Command::new("python3");
        cmd.args(["-I", "-c", &PYTHON]).arg(&f.0).args([
            allow,
            skip,
            patterns,
            if cached { "1" } else { "" },
        ]);
        for file in files {
            cmd.arg(f.0.join(file));
        }
        let output = cmd.output().unwrap();
        assert!(
            output.status.code().is_some_and(|c| c == 0 || c == 1),
            "oracle failed"
        );
        assert!(output.stderr.is_empty(), "oracle diagnostic");
        String::from_utf8(output.stdout).unwrap()
    }

    fn assert_oracle_findings(f: &Fixture, files: &[&str], cached: bool) {
        assert_options(f, files, cached, "", "", "", true);
    }

    fn assert_options(
        f: &Fixture,
        files: &[&str],
        cached: bool,
        allow: &str,
        skip: &str,
        patterns: &str,
        require_finding: bool,
    ) {
        let expected = oracle_options(f, files, cached, allow, skip, patterns);
        let mut args = f.args(files);
        args[3] = allow.into();
        args[5] = skip.into();
        if !patterns.is_empty() {
            args.extend(["--patterns".into(), patterns.into()]);
        }
        if cached {
            args.insert(0, "--cached".into());
        }
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = cli(&args, &mut out, &mut err);
        let out = String::from_utf8(out).unwrap();
        assert!(err.is_empty(), "unexpected scanner diagnostic");
        let reported: std::collections::HashSet<_> = out.lines().collect();
        let mut count = 0;
        let mut required_reports = Vec::new();
        for line in expected.lines().filter(|l| l.starts_with('[')) {
            let pycompat::json::Value::Array(row) = pycompat::json::loads(line).unwrap() else {
                panic!("metadata")
            };
            let [
                pycompat::json::Value::Str(kind),
                pycompat::json::Value::Str(rel),
                pycompat::json::Value::Int(n),
            ] = row.as_slice()
            else {
                panic!("metadata fields")
            };
            // Keep the required report spelling independent of report().
            let finding = format!("FAIL {kind:<9} {}:{n}", rel.escape_debug());
            assert!(
                reported.contains(finding.as_str()),
                "missing oracle finding #{count}"
            );
            required_reports.push(finding);
            count += 1;
        }
        let mut actual_reports: Vec<_> = out
            .lines()
            .filter(|l| l.starts_with("FAIL "))
            .map(str::to_owned)
            .collect();
        actual_reports.sort();
        required_reports.sort();
        assert_eq!(
            actual_reports, required_reports,
            "scanner metadata multiset differs"
        );
        if require_finding {
            assert!(count > 0, "oracle control was clean");
        }
        assert_eq!(code, u8::from(count > 0));
        assert_eq!(
            out.lines().filter(|l| l.starts_with("FAIL ")).count(),
            count
        );
    }

    #[test]
    fn privacy_report_never_contains_values_or_patterns() {
        let f = Fixture::new("redaction");
        let secret = p("/home/~synthetic_private_marker/path");
        std::fs::write(f.0.join("probe.md"), &secret).unwrap();
        let expected = oracle(&f, &["probe.md"], false);
        assert!(expected.contains("home"));
        let (mut out, mut err) = (Vec::new(), Vec::new());
        assert_eq!(cli(&f.args(&["probe.md"]), &mut out, &mut err), 1);
        for stream in [&out, &err] {
            assert!(
                !String::from_utf8_lossy(stream).contains(&secret),
                "reported private value"
            );
        }
        assert_oracle_findings(&f, &["probe.md"], false);
        let pattern = "CUSTOM_PRIVATE_PATTERN_MARKER[";
        std::fs::write(f.0.join("patterns"), pattern).unwrap();
        let mut args = f.args(&["probe.md"]);
        args.extend([
            "--patterns".into(),
            f.0.join("patterns").display().to_string(),
        ]);
        out.clear();
        err.clear();
        assert_eq!(cli(&args, &mut out, &mut err), 2);
        for stream in [&out, &err] {
            assert!(
                !String::from_utf8_lossy(stream).contains("CUSTOM_PRIVATE_PATTERN_MARKER"),
                "reported custom pattern"
            );
        }
    }

    #[test]
    fn privacy_python_boundaries_differential() {
        let f = Fixture::new("boundaries");
        let text = p(
            "\u{301}rs123~ AG\nHG0~02\u{301} *1/~*4\nHTT\u{301} 17/~19\nREPCN\u{301} 17/~19\n/mnt/u~ser\u{301}\nx.ts.~net\u{301}\nrs123\t1\t9\t~AA\u{301}\n",
        );
        std::fs::write(f.0.join("probe.md"), text).unwrap();
        assert_oracle_findings(&f, &["probe.md"], false);
    }

    #[test]
    fn privacy_python_casefold_differential() {
        let f = Fixture::new("casefold");
        let text = p(
            "patİent *1/~*4\npatıent *1/~*4\nſample *1/~*4\nİpatient *1/~*4\nıpatient *1/~*4\nKpatient *1/~*4\n",
        );
        std::fs::write(f.0.join("probe.md"), text).unwrap();
        assert_oracle_findings(&f, &["probe.md"], false);
    }

    #[test]
    fn privacy_staged_newline_differential() {
        let f = Fixture::new("staged");
        let git = |args: &[&str]| {
            assert!(
                Command::new("git")
                    .arg("-C")
                    .arg(&f.0)
                    .args(args)
                    .status()
                    .unwrap()
                    .success()
            )
        };
        git(&["init", "-q"]);
        std::fs::write(f.0.join("clean.md"), "clean\n").unwrap();
        std::fs::write(
            f.0.join("odd\nname.md"),
            p("/home/~synthetic_private_marker/path\n"),
        )
        .unwrap();
        git(&["add", "--", "clean.md", "odd\nname.md"]);
        std::fs::write(f.0.join("odd\nname.md"), "clean\n").unwrap();
        assert_oracle_findings(&f, &["clean.md", "odd\nname.md"], true);
    }

    #[test]
    fn privacy_fixed_seed_differential_fuzz() {
        let f = Fixture::new("fuzz");
        let vocabulary = [
            "rs123~ AG",
            "HG0~02 *1/~*4",
            "patient *1/~*4",
            "patİent *1/~*4",
            "patıent *1/~*4",
            "ſample *1/~*4",
            "HTT 17/~19",
            "REPCN 17/~19",
            "/home/~synthetic/path",
            "/mnt/u~ser",
            "node.ts.~net",
            "10.2.~3.4",
            "rs123\t1\t9\t~AA",
            "clean",
        ];
        let noise = [
            "", "\u{301}", "\u{345}", "\u{200d}", "\u{200b}", "İ", "ı", "ſ", "K", "é", "²", "_",
            "\u{1c}", "\u{a0}", "\r", "\t",
        ];
        let mut seed = 0x5ca1ab1eu64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed as usize
        };
        let mut data = Vec::new();
        for _ in 0..6144 {
            let value = p(vocabulary[next() % vocabulary.len()]);
            let before = noise[next() % noise.len()];
            let after = noise[next() % noise.len()];
            // Noise at every internal boundary, not only outside the line.
            let value = value.replace(' ', &format!("{after} {before}"));
            data.extend_from_slice(format!("{before}{value}{after}\r\n").as_bytes());
        }
        data.extend_from_slice(&[0xff, 0xfe]);
        data.extend_from_slice(p("rs123~ AG\n").as_bytes());
        data.extend(std::iter::repeat_n(b'x', 100_000));
        data.extend_from_slice(p(" HTT 17/~19\n\0REPCN 17/~19\n").as_bytes());
        let names = [
            "noise.md",
            "odd\nname.md",
            "tab\tname.md",
            "-dash.md",
            "unicode-\u{301}.md",
            "binary.md",
        ];
        for name in &names[..5] {
            std::fs::write(f.0.join(name), &data).unwrap();
        }
        std::fs::write(f.0.join(names[5]), b"\0binary").unwrap();
        assert_oracle_findings(&f, &names, false);
    }

    #[test]
    fn round2_boundary_backtracking() {
        let f = Fixture::new("r2-boundary");
        std::fs::write(f.0.join("chip.md"), "rs987654\t1\t9\tA-\n").unwrap();
        assert_oracle_findings(&f, &["chip.md"], false);
        assert!(
            Command::new("git")
                .arg("-C")
                .arg(&f.0)
                .args(["init", "-q"])
                .status()
                .unwrap()
                .success()
        );
        assert!(
            Command::new("git")
                .arg("-C")
                .arg(&f.0)
                .args(["add", "chip.md"])
                .status()
                .unwrap()
                .success()
        );
        assert_oracle_findings(&f, &["chip.md"], true);
        assert_oracle_findings(&f, &[], false);
        assert_oracle_findings(&f, &[], true);
    }

    #[test]
    fn round2_custom_unicode() {
        let f = Fixture::new("r2-custom");
        for (pattern, text) in [
            (r"\bZetaSecret\b", "\u{301}ZetaSecret"),
            (r"(?i)patient", "patİent"),
            (r"(?i)[A-Z]+", "İı"),
            (r"\s+", "\u{1c}"),
            (r"\W+", "\u{301}"),
            (r"[^\w]+", "\u{301}"),
        ] {
            std::fs::write(f.0.join("patterns"), pattern).unwrap();
            std::fs::write(f.0.join("own.md"), text).unwrap();
            assert_options(
                &f,
                &["own.md"],
                false,
                "",
                "",
                f.0.join("patterns").to_str().unwrap(),
                true,
            );
        }
    }

    #[test]
    fn round2_cr_patterns() {
        let f = Fixture::new("r2-pattern-lines");
        for prefix in ["# a comment", "Unrelated"] {
            std::fs::write(
                f.0.join("patterns"),
                format!("{prefix}\rZETA_PRIVATE_MARKER\r"),
            )
            .unwrap();
            std::fs::write(f.0.join("own.md"), "ZETA_PRIVATE_MARKER\n").unwrap();
            assert_options(
                &f,
                &["own.md"],
                false,
                "",
                "",
                f.0.join("patterns").to_str().unwrap(),
                true,
            );
        }
    }

    #[test]
    fn round2_empty_alternative() {
        let f = Fixture::new("r2-empty");
        std::fs::write(f.0.join("patterns"), "^|PUBLIC_PRIVATE\n").unwrap();
        std::fs::write(f.0.join("own.md"), "PUBLIC_PRIVATE").unwrap();
        assert_options(
            &f,
            &["own.md"],
            false,
            "own.md PUBLIC",
            "",
            f.0.join("patterns").to_str().unwrap(),
            true,
        );
    }

    #[test]
    fn round2_exception_lines() {
        let f = Fixture::new("r2-allow");
        let text = p("safe\rOTHER.md /home/~synthetic_private_marker/path\n");
        std::fs::write(f.0.join("builtin-allow.md"), &text).unwrap();
        assert_options(
            &f,
            &["builtin-allow.md"],
            false,
            &format!("builtin-allow.md {text}"),
            "",
            "",
            true,
        );
    }

    #[test]
    fn round2_generated_patterns_and_tables() {
        let f = Fixture::new("r2-grammar");
        std::fs::write(f.0.join("clean.md"), "clean").unwrap();
        let atoms = [
            "PUBLIC",
            "patient",
            r"\w",
            r"\W",
            r"\s",
            r"\S",
            r"\d",
            r"\D",
            "[A-Z]",
            r"[^\w]",
            r"\bZetaSecret\b",
            r"\BZetaSecret",
        ];
        let quantifiers = ["", "?", "+", "*", "{1,2}", "{0,2}?"];
        let noise = [
            "", "\u{301}", "\u{200c}", "\u{200d}", "\u{1c}", "\u{85}", "\u{2028}", "\u{2029}",
            "İı", "ſK", "²", "_", "\t",
        ];
        let separators = [
            "\n", "\r\n", "\r", "\u{85}", "\u{2028}", "\u{2029}", "\u{b}", "\u{1c}",
        ];
        let mut seed = 0x2d5ca1ab1eu64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed as usize
        };
        for case in 0..96 {
            let atom = if case < atoms.len() {
                atoms[case]
            } else {
                atoms[next() % atoms.len()]
            };
            let quantified = format!("(?:{atom}){}", quantifiers[next() % quantifiers.len()]);
            let expression = if case < atoms.len() {
                atom.to_owned()
            } else {
                match next() % 4 {
                    0 => quantified,
                    1 => format!("^(?:{quantified})"),
                    2 => format!("(?:{quantified})$"),
                    _ => format!("^|{quantified}|PUBLIC_PRIVATE"),
                }
            };
            let pattern = format!("{}{expression}", if next() % 2 == 0 { "(?i)" } else { "" });
            // Unicode separators are pattern characters, not text-mode
            // newlines. Only CR, CRLF and LF separate pattern-file rows.
            let sep = if case < atoms.len() {
                "\n"
            } else {
                separators[case % separators.len()]
            };
            std::fs::write(
                f.0.join("patterns"),
                format!("# comment{sep}{pattern}{sep}Unrelated{sep}"),
            )
            .unwrap();
            let mut text = String::from("PUBLIC_PRIVATE\n");
            for value in noise {
                text.push_str(value);
                text.push('\n');
            }
            for _ in 0..24 {
                text.push_str(noise[next() % noise.len()]);
                text.push_str(
                    [
                        "ZetaSecret",
                        "PUBLIC_PRIVATE",
                        "patİent",
                        "İı",
                        "a9",
                        "\u{1c}",
                    ][next() % 6],
                );
                text.push_str(noise[next() % noise.len()]);
                text.push_str(separators[next() % separators.len()]);
            }
            for ending in ["A-", "-A", "--", "AA", "A"] {
                text.push_str(&format!("rs987654\t1\t9\t{ending}\n"));
            }
            std::fs::write(f.0.join("own.md"), text).unwrap();
            let allow = format!("own.md PUBLIC{sep}OTHER.md PRIVATE{sep}own.md Zeta");
            let skip = format!("missing.md reason{sep}OTHER.md reason");
            assert_options(
                &f,
                &["own.md", "clean.md"],
                false,
                &allow,
                &skip,
                f.0.join("patterns").to_str().unwrap(),
                false,
            );
        }
        // Cross-file scope and SKIP parsing, with a scanned clean companion.
        for sep in separators {
            let secret = p("safe /home/~synthetic/path");
            std::fs::write(f.0.join("own.md"), &secret).unwrap();
            let allow = format!("own.md safe{sep}OTHER.md {}", p("/home/~synthetic/path"));
            assert_options(&f, &["own.md", "clean.md"], false, &allow, "", "", true);
            let skip = format!("missing.md reason{sep}own.md reason");
            assert_options(&f, &["own.md", "clean.md"], false, "", &skip, "", false);
        }
    }

    #[test]
    fn round2_unsupported_and_exhausted_patterns_fail_closed() {
        let f = Fixture::new("r2-unsupported");
        std::fs::write(f.0.join("own.md"), format!("{}!", "a".repeat(64))).unwrap();
        for (index, pattern) in [
            r"\b*",
            r"(?i)*",
            r"(?P<a>x)(?P<a>y)",
            r"(?=PRIVATE)",
            r"(PRIVATE)\1",
            r"\p{L}",
            r"[[:alpha:]]",
            r"(?a)PRIVATE",
            r"\N{LATIN CAPITAL LETTER A}",
            r"(a|\Ba)+$",
        ]
        .into_iter()
        .enumerate()
        {
            std::fs::write(f.0.join("patterns"), pattern).unwrap();
            let mut args = f.args(&["own.md"]);
            args.extend([
                "--patterns".into(),
                f.0.join("patterns").display().to_string(),
            ]);
            let (mut out, mut err) = (Vec::new(), Vec::new());
            assert_eq!(cli(&args, &mut out, &mut err), 2, "control {index}");
            assert!(!String::from_utf8_lossy(&err).contains("PRIVATE"));
            assert!(!String::from_utf8_lossy(&out).contains("PRIVATE"));
        }
    }

    #[test]
    fn round3_nullable_repetitions_fail_closed() {
        let f = Fixture::new("r3-nullable");
        assert!(
            Command::new("git")
                .arg("-C")
                .arg(&f.0)
                .args(["init", "-q"])
                .status()
                .unwrap()
                .success()
        );
        for (pattern, text, allow) in [
            ("(a*|b*?){1,3}", "ba", "own.md b\nown.md a"),
            (
                "(?:a*|[ab]{0,2}){0,3}",
                "baa",
                "own.md\nown.md a\nown.md ba",
            ),
        ] {
            // Reconfirm that a partial ALLOW does not cover Python's match.
            std::fs::write(f.0.join("patterns"), pattern).unwrap();
            std::fs::write(f.0.join("own.md"), text).unwrap();
            assert!(
                Command::new("git")
                    .arg("-C")
                    .arg(&f.0)
                    .args(["add", "own.md"])
                    .status()
                    .unwrap()
                    .success()
            );
            for cached in [false, true] {
                for files in [&["own.md"][..], &[][..]] {
                    let expected = oracle_options(
                        &f,
                        files,
                        cached,
                        allow,
                        "",
                        f.0.join("patterns").to_str().unwrap(),
                    );
                    assert!(
                        expected.contains("[\"own\""),
                        "partial ALLOW must leave a Python finding"
                    );
                    let mut args = f.args(files);
                    args[3] = allow.into();
                    args.extend([
                        "--patterns".into(),
                        f.0.join("patterns").to_str().unwrap().into(),
                    ]);
                    if cached {
                        args.push("--cached".into());
                    }
                    let (mut out, mut err) = (Vec::new(), Vec::new());
                    assert_eq!(
                        cli(&args, &mut out, &mut err),
                        2,
                        "nullable repetition must be refused"
                    );
                    assert!(!String::from_utf8(err).unwrap().contains(pattern));
                    assert!(out.is_empty());
                }
            }
        }
    }

    #[test]
    fn round3_raw_byte_paths_fail_closed() {
        let f = Fixture::new("r3-raw-path");
        let git = |args: &[&str], input: &[u8]| {
            let mut child = Command::new("git")
                .arg("-C")
                .arg(&f.0)
                .args(args)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .spawn()
                .unwrap();
            child.stdin.take().unwrap().write_all(input).unwrap();
            let output = child.wait_with_output().unwrap();
            assert!(output.status.success());
            output.stdout
        };
        git(&["init", "-q"], b"");
        std::fs::write(f.0.join("clean.md"), "clean").unwrap();
        git(&["add", "clean.md"], b"");
        let hash = git(
            &["hash-object", "-w", "--stdin"],
            p("/home/~reviewer_three_synthetic/path").as_bytes(),
        );
        let mut entry = b"100644 ".to_vec();
        entry.extend_from_slice(hash.trim_ascii());
        entry.extend_from_slice(b"\tbad\xff.md\0");
        git(&["update-index", "-z", "--index-info"], &entry);
        for cached in [false, true] {
            let mut args = f.args(&[]);
            if cached {
                args.push("--cached".into());
            }
            let (mut out, mut err) = (Vec::new(), Vec::new());
            assert_eq!(
                cli(&args, &mut out, &mut err),
                2,
                "invalid tracked path must stop the whole scan"
            );
            assert!(out.is_empty());
            assert_eq!(
                String::from_utf8(err).unwrap(),
                "ERROR: path is not valid UTF-8\n"
            );
            // Explicit paths reach the same validation before SKIP or I/O.
            let mut a = parse(&f.args(&["clean.md"])).unwrap();
            a.cached = cached;
            a.files.push(f.0.join(OsStr::from_bytes(b"bad\xff.md")));
            a.skip.push("*".into());
            assert_eq!(
                run(&a, &mut Vec::new()).unwrap_err(),
                "path is not valid UTF-8"
            );
        }
    }
    #[test]
    fn round4_unicode_versions_fail_closed() {
        let f = Fixture::new("r4-unicode");
        assert!(
            Command::new("git")
                .arg("-C")
                .arg(&f.0)
                .args(["init", "-q"])
                .status()
                .unwrap()
                .success()
        );
        let cases = [
            (p("\u{105c0}rs555~ AG"), None),
            (p("HG004\u{105c0} *2/~*8"), None),
            (p("HTT\u{105c0} 18/~32"), None),
            ("\u{105c0}".into(), Some(r"\W+")),
            ("\u{1e6c0}".into(), Some(r"\w+")),
            ("\u{16ebb}".into(), Some(r"(?i)\U00016ea0")),
        ];
        let mut failures = Vec::new();
        for (case, (text, pattern)) in cases.iter().enumerate() {
            std::fs::write(f.0.join("own.md"), text).unwrap();
            std::fs::write(f.0.join("patterns"), pattern.unwrap_or("# none")).unwrap();
            assert!(
                Command::new("git")
                    .arg("-C")
                    .arg(&f.0)
                    .args(["add", "own.md"])
                    .status()
                    .unwrap()
                    .success()
            );
            for cached in [false, true] {
                for explicit in [false, true] {
                    let mut args = f.args(if explicit { &["own.md"] } else { &[] });
                    if cached {
                        args.push("--cached".into());
                    }
                    if pattern.is_some() {
                        args.extend([
                            "--patterns".into(),
                            f.0.join("patterns").display().to_string(),
                        ]);
                    }
                    let (mut out, mut err) = (Vec::new(), Vec::new());
                    let code = cli(&args, &mut out, &mut err);
                    if code != 2 {
                        failures.push((case, cached, explicit, code));
                    }
                    assert!(!String::from_utf8_lossy(&out).contains(text));
                    assert!(!String::from_utf8_lossy(&err).contains(text));
                }
            }
        }
        assert!(
            failures.is_empty(),
            "version-dependent text must fail closed: {failures:?}"
        );
    }

    #[test]
    fn round4_symlink_targets_fail_closed() {
        let f = Fixture::new("r4-symlinks");
        assert!(
            Command::new("git")
                .arg("-C")
                .arg(&f.0)
                .args(["init", "-q"])
                .status()
                .unwrap()
                .success()
        );
        std::fs::write(f.0.join("clean.md"), "clean").unwrap();
        std::os::unix::fs::symlink(
            OsStr::from_bytes(b"private_target_\xff.md"),
            f.0.join("alias.md"),
        )
        .unwrap();
        assert!(
            Command::new("git")
                .arg("-C")
                .arg(&f.0)
                .args(["add", "alias.md", "clean.md"])
                .status()
                .unwrap()
                .success()
        );
        let mut failures = Vec::new();
        for cached in [false, true] {
            for explicit in [false, true] {
                let mut args = f.args(if explicit {
                    &["alias.md", "clean.md"]
                } else {
                    &[]
                });
                if cached {
                    args.push("--cached".into());
                }
                let (mut out, mut err) = (Vec::new(), Vec::new());
                let code = cli(&args, &mut out, &mut err);
                if code != 2 {
                    failures.push((cached, explicit, code));
                }
                if code == 2 {
                    assert!(out.is_empty());
                    assert_eq!(err, b"ERROR: path is not valid UTF-8\n");
                }
            }
        }
        assert!(
            failures.is_empty(),
            "raw symlink target must stop scan: {failures:?}"
        );
    }
}
