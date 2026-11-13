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
use std::ffi::OsStr;
use std::io::{BufRead, BufReader, Read, Write};
use std::os::unix::ffi::OsStrExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::LazyLock;

use pycompat::mtime::fnmatch;
use pycompat::num::py_strip;

use crate::lookaround::Rule;

const USAGE: &str = "usage: og dev scan personal-data --root DIR --allow TEXT --skip TEXT [--patterns FILE] [--cached] [--] [FILE...]\n";

/// This file's own path: it holds every allowed text, so all of them apply.
const SELF: &str = "scripts/ci/check-personal-data.sh";

/// How many leading bytes are searched for a NUL to call a file binary.
const BINARY_PROBE: usize = 8192;

fn not_word(c: char) -> bool {
    c.is_ascii_alphanumeric() || c == '_'
}

fn digit_or_dot(c: char) -> bool {
    c.is_ascii_digit() || c == '.'
}

fn digit_dot_slash(c: char) -> bool {
    c.is_ascii_digit() || c == '.' || c == '/'
}

/// `[A-Za-z]` under Python's `(?i)`, which also folds to the long s and
/// the Kelvin sign, plus dotted and dotless I.
fn letter_ci(c: char) -> bool {
    c.is_ascii_alphabetic() || c == '\u{130}' || c == '\u{131}' || c == '\u{17F}' || c == '\u{212A}'
}

const STR_GENES: &str = "HTT|FMR1|C9orf72|C9ORF72|ATXN[0-9]+|DMPK|FXN|CNBP|RFC1|PABPN1|CACNA1A|TCF4|AR|PPP2R2B|TBP|JPH3|NOP56|ATN1|CSTB|AFF2|FMR2|NOTCH2NLC|DAB1|BEAN1|STARD7|YEATS2";
const OCTET: &str = "(?:25[0-5]|2[0-4][0-9]|1[0-9][0-9]|[1-9]?[0-9])";
/// Python's `\s` also holds the four separators U+001C to U+001F.
const SPACE: &str = r"\s\x1c-\x1f";

/// A rule: its kind, its pattern, and a pattern the line must also match.
struct Check {
    kind: &'static str,
    rule: Rule,
    cond: Option<&'static [Rule]>,
}

fn rule(pattern: &str, not_after: Option<fn(char) -> bool>) -> Rule {
    Rule::python(pattern, not_after).expect("a valid built-in pattern")
}

/// A sample label, which may carry a number after a separator: sample7,
/// sample_01, donor-3; or a public sample id such as NA12878 or HG002.
static LABEL: LazyLock<[Rule; 2]> = LazyLock::new(|| {
    [
        rule(
            r"(?i)((?:sample|pat[iİı]ent|subject|proband|part[iİı]c[iİı]pant|[iİı]nd[iİı]v[iİı]dual|donor)s?(?:[_-]?[0-9]+)?)(?:[^A-Za-zİı]|$)",
            Some(letter_ci),
        ),
        rule(r"\b(?:NA|HG)[0-9]{3,5}\b", None),
    ]
});

static REPCN: LazyLock<[Rule; 1]> = LazyLock::new(|| [rule(r"\bREPCN\b", None)]);

/// The built-in rules, in the order they are reported. Each pattern marks
/// only the leaked value, so an allowed example covers just its span.
fn builtin() -> Vec<Check> {
    let gt = r"\(?[ACGTDI][/|;]?[ACGTDI]\)?";
    vec![
        Check {
            kind: "home",
            rule: rule(
                r"/home/[A-Za-z]|/Users/[A-Za-z]|/mnt/use[r]\b|[A-Za-z]:\\\\?Users\\\\?[A-Za-z]",
                Some(not_word),
            ),
            cond: None,
        },
        Check {
            kind: "address",
            rule: rule(
                &format!(
                    r"((?:10\.{OCTET}|172\.(?:1[6-9]|2[0-9]|3[01])|192\.168|100\.(?:6[4-9]|[7-9][0-9]|1[01][0-9]|12[0-7]))\.{OCTET}\.{OCTET})(?:$|\.$|\.[^0-9]|[^0-9.])"
                ),
                Some(digit_or_dot),
            ),
            cond: None,
        },
        Check {
            kind: "address",
            rule: rule(r"[A-Za-z0-9-]\.ts\.net\b", None),
            cond: None,
        },
        Check {
            kind: "diplotype",
            rule: rule(r"\*[0-9]+[A-Za-z0-9.]*/\*[0-9]+", None),
            cond: Some(&*LABEL),
        },
        Check {
            kind: "rsid",
            rule: rule(
                &format!(r#"(\brs[0-9]+\b[{SPACE}:=,;|`'"\[\]-]{{1,6}}{gt})(?:[^A-Za-z0-9]|$)"#),
                None,
            ),
            cond: None,
        },
        Check {
            kind: "rsid",
            rule: rule(
                r"\brs[0-9]+\t[0-9XYMT]{1,2}\t[0-9]+\t[ACGTDI-]{1,2}\b",
                None,
            ),
            cond: None,
        },
        Check {
            kind: "repeat",
            rule: rule(
                r"([0-9]{1,4}/[0-9]{1,4})(?:[^0-9/]|$)",
                Some(digit_dot_slash),
            ),
            cond: Some(&*REPCN),
        },
        // Python: \b(?:GENES)\b[^/\n]{0,40}?(?<![0-9.])N/N(?![0-9/]). The
        // gap after the gene is never empty (a digit right after a gene
        // name is no word boundary), so its last character carries the
        // look-behind.
        Check {
            kind: "repeat",
            rule: rule(
                &format!(
                    r"(\b(?:{STR_GENES})\b[^/\n]{{0,39}}?[^/\n0-9.][0-9]{{1,4}}/[0-9]{{1,4}})(?:[^0-9/]|$)"
                ),
                None,
            ),
            cond: None,
        },
    ]
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
    let allow = allow
        .ok_or("--allow is required")?
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let l = py_strip(l);
            let (path, text) = l.split_once(' ').unwrap_or((l, ""));
            (path.to_owned(), py_strip(text).to_owned())
        })
        .collect();
    let skip = skip
        .ok_or("--skip is required")?
        .lines()
        .filter_map(|l| l.split_whitespace().next().map(str::to_owned))
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
    for (line, l) in text.lines().enumerate() {
        let p = py_strip(l);
        if p.is_empty() || l.starts_with('#') {
            continue;
        }
        // Rule reports capture 1 for the built-in look-around patterns.
        // A user's captures must not narrow the span checked against ALLOW.
        let rule = Rule::whole_match(p).map_err(|_| {
            format!(
                "{}:{}: invalid custom pattern",
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
fn relative(path: &Path, root: &Path) -> String {
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
        Ok(rel) if !rel.as_os_str().is_empty() => rel.to_string_lossy().into_owned(),
        _ => path.to_string_lossy().into_owned(),
    }
}

fn git_files(root: &Path) -> Result<Vec<PathBuf>, String> {
    let out = Command::new("git")
        .arg("-C")
        .arg(root)
        .args(["ls-files", "-z"])
        .stderr(Stdio::inherit())
        .output()
        .map_err(|e| format!("cannot run git: {e}"))?;
    if !out.status.success() {
        return Err(format!("git ls-files failed in {}", root.display()));
    }
    Ok(out
        .stdout
        .split(|&b| b == 0)
        .filter(|f| !f.is_empty())
        .map(|f| root.join(OsStr::from_bytes(f)))
        .collect())
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
        let head = String::from_utf8_lossy(&header);
        let fields: Vec<&str> = head.split(' ').collect();
        let [_, kind, size] = fields[..] else {
            continue; // "<name> missing" or "ambiguous"
        };
        let size: usize = size.parse().map_err(|_| format!("git cat-file: {head}"))?;
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
fn scan_text(out: &mut dyn Write, rules: &[Check], texts: &[&str], rel: &str, data: &[u8]) -> bool {
    let mut bad = false;
    let text = String::from_utf8_lossy(data);
    for (i, line) in text.split('\n').enumerate() {
        let spans: Vec<(usize, usize)> = texts
            .iter()
            .flat_map(|t| line.match_indices(*t).map(|(a, s)| (a, a + s.len())))
            .collect();
        for c in rules {
            if c.cond
                .is_some_and(|cs| !cs.iter().any(|r| r.is_match(line)))
            {
                continue;
            }
            let leaked = c
                .rule
                .find_iter(line)
                .find(|&(s, e)| !spans.iter().any(|&(a, b)| a <= s && e <= b));
            if leaked.is_some() {
                report(out, c.kind, rel, i + 1);
                bad = true;
            }
        }
    }
    bad
}

fn run(a: &Args, out: &mut dyn Write) -> Result<bool, String> {
    let root = a
        .root
        .canonicalize()
        .map_err(|e| format!("--root {}: {e}", a.root.display()))?;
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
        .map(|f| (f, relative(f, &root)))
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
        bad |= scan_text(out, &rules, &texts, rel, &data);
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

    fn kinds(line: &str) -> Vec<&'static str> {
        let rules = builtin();
        let mut out = Vec::new();
        for c in &rules {
            if c.cond
                .is_some_and(|cs| !cs.iter().any(|r| r.is_match(line)))
            {
                continue;
            }
            if c.rule.is_match(line) {
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
        assert!(!scan_text(
            &mut out,
            &rules,
            &[ex.as_str()],
            "x.md",
            ex.as_bytes()
        ));
        let line = p(&format!("{ex} see /home/~bob/x"));
        assert!(scan_text(
            &mut out,
            &rules,
            &[ex.as_str()],
            "x.md",
            line.as_bytes()
        ));
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
        assert!(scan_text(&mut out, &rules, &["PUBLIC"], "x.md", line));
        let got = String::from_utf8(out).unwrap();
        assert_eq!(got.lines().count(), 5, "{got}");
        assert!(!scan_text(
            &mut Vec::new(),
            &rules,
            &["PREFIX PUBLIC SUFFIX"],
            "x.md",
            line,
        ));
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
        let python = python
            .replace(
                "import fnmatch, os, re, subprocess, sys",
                "import fnmatch, os, re, subprocess, sys, json",
            )
            .replace(
                "print(\"FAIL %-9s %s:%d: %s\" % (kind, rel, n, line.strip()[:160]))",
                "print(json.dumps([kind, rel, n]))",
            );
        let mut cmd = Command::new("python3");
        cmd.args(["-I", "-c", &python])
            .arg(&f.0)
            .args(["", "", "", if cached { "1" } else { "" }]);
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
        let expected = oracle(f, files, cached);
        let mut args = f.args(files);
        if cached {
            args.insert(0, "--cached".into());
        }
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = cli(&args, &mut out, &mut err);
        let out = String::from_utf8(out).unwrap();
        assert!(err.is_empty(), "unexpected scanner diagnostic");
        let reported: std::collections::HashSet<_> = out.lines().collect();
        let mut count = 0;
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
                reported.iter().any(|l| *l == finding
                    || l.strip_prefix(&finding)
                        .is_some_and(|tail| tail.starts_with(": "))),
                "missing oracle finding #{count}"
            );
            count += 1;
        }
        assert!(count > 0, "oracle control was clean");
        assert_eq!(code, 1);
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
}
