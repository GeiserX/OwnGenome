//! `og dev python-gate [--root DIR] [--allow FILE] [--list]`: no owned
//! Python is left in the repository.
//!
//! Every file `git ls-files` lists (tracked, plus new files that are not
//! ignored) is read, and these are the findings:
//!
//! - `file`: a `*.py`, `*.pyw` or `*.pyi` file;
//! - `shebang`: a `#!` line whose program is python, python3 or
//!   pythonX.Y, directly or through `env`: the first line of any file, or
//!   any line of a shell, Nextflow or workflow file (a script a heredoc or a
//!   Nextflow `script:` block writes);
//! - `call`: the word python, python3 or pythonX.Y outside a comment in a
//!   shell script, a Nextflow file or config, a GitHub workflow or action,
//!   or a `command` column of a TSV command table (`tests/smoke/commands.tsv`):
//!   a `python3 -c`, a `python3 - <<EOF`, a `python3 tool` run in an image;
//! - `heredoc`: a heredoc whose delimiter says it holds Python (`<<'PY'`,
//!   `<<PYTHON`, `<<EOF_PYTHON`);
//! - `rust`: a `Command::new("python3")` (any python program, any path)
//!   or the variable that named a Python oracle, in Rust code outside
//!   comments.
//!
//! A comment, a doc sentence and a Markdown, text, JSON or image file are
//! prose: they never count. `python:3.15` (an image), `python=...` (a
//! stamp), `PYTHON_IMAGE`, `setup-python` and a `site-packages/python3.12/`
//! path are not calls.
//!
//! `tests/python-allow.tsv` names, one row per file, the files whose
//! findings are third-party Python: a tool run inside its own pinned image
//! (PharmCAT's preprocessor, eKLIPse, Cyrius in `PYTHON_IMAGE`), or the
//! docs builder. Each row is `file<TAB>reason`, and the gate fails when
//!
//! - a finding is in a file no row names;
//! - a row names a missing file, a directory or a glob, or has no reason;
//! - a reason names an `*_IMAGE` variable that `versions.env` does not set;
//! - a row covers nothing (the Python left that file: delete the row);
//! - a row names a Python file (a `file` or `shebang` finding): every Python
//!   file in the repository is owned code, so none can be allowed, and an
//!   allowed file can never import a module from `bin/` or the demo
//!   builders;
//! - a file is named twice.
//!
//! `--list` prints every finding, allowed or not. Exit 0 when all is
//! allowed, 1 when a check fails, 2 on a usage or read error.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::versions_env;

const USAGE: &str = "usage: og dev python-gate [--root DIR] [--allow FILE] [--list]\n";

/// The allow list, relative to the root.
pub const ALLOW: &str = "tests/python-allow.tsv";

/// The variable a Python oracle was named by. Written in two halves so this
/// file does not hold it as one word and the gate can scan itself.
fn oracle_var() -> String {
    ["OG_", "PYTHON"].concat()
}

/// One thing the gate found.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Finding {
    /// The path, relative to the root, with `/`.
    pub file: String,
    /// The 1-based line, 0 for the file itself.
    pub line: usize,
    /// `file`, `shebang`, `call`, `heredoc` or `rust`.
    pub kind: &'static str,
    /// The line, trimmed and cut to 120 characters.
    pub text: String,
}

/// How the gate reads a file's lines.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Lang {
    /// Shell, workflow YAML: `#` comments.
    Shell,
    /// Nextflow and Groovy: `//` and `/* */` comments, and `#` comments in
    /// the shell of `script:` blocks.
    Nextflow,
    /// A TSV command table: `#` lines are comments, the rest is shell.
    CommandTable,
    /// Rust.
    Rust,
    /// Anything else: only the file name and the first line count.
    Other,
}

fn lang(rel: &str, text: &str) -> Lang {
    let first_line = text.lines().next().unwrap_or("");
    let name = rel.rsplit('/').next().unwrap_or(rel);
    let ext = name.rsplit_once('.').map(|(_, e)| e).unwrap_or("");
    match ext {
        "sh" | "bash" | "ksh" | "zsh" => return Lang::Shell,
        "nf" | "config" | "groovy" => return Lang::Nextflow,
        "rs" => return Lang::Rust,
        "yml" | "yaml" if rel.starts_with(".github/") => return Lang::Shell,
        "tsv" if is_command_table(text) => return Lang::CommandTable,
        _ => {}
    }
    if name.starts_with("Dockerfile") || name.starts_with("Containerfile") || name == "Makefile" {
        return Lang::Shell;
    }
    let shell_script = shebang_program(first_line)
        .is_some_and(|p| matches!(p.as_str(), "sh" | "bash" | "dash" | "ksh" | "zsh"));
    if ext.is_empty() && shell_script {
        return Lang::Shell;
    }
    Lang::Other
}

/// Does this TSV have a `command` column? Its header is one of the leading
/// `#` lines (`#var<TAB>opts<TAB>needs<TAB>command<TAB>expect`) or the
/// first line after them.
fn is_command_table(text: &str) -> bool {
    let mut header = Vec::new();
    for line in text.lines().take(20) {
        header.push(line);
        if !line.starts_with('#') {
            break;
        }
    }
    header.iter().any(|l| {
        l.trim_start_matches('#')
            .split('\t')
            .any(|c| c.trim() == "command")
    })
}

/// The program a `#!` line runs: the base name of its first word, or of the
/// first word after `env` and its options and assignments.
fn shebang_program(line: &str) -> Option<String> {
    let rest = line.trim_start().strip_prefix("#!")?;
    let mut words = rest.split_whitespace();
    let base = |w: &str| w.rsplit('/').next().unwrap_or(w).to_owned();
    let first = base(words.next()?);
    if first != "env" {
        return Some(first);
    }
    let mut skip_value = false;
    for w in words {
        if skip_value {
            skip_value = false;
            continue;
        }
        if w == "-u" || w == "-C" || w == "--unset" || w == "--chdir" {
            skip_value = true;
            continue;
        }
        if w.starts_with('-') || w.contains('=') {
            continue;
        }
        return Some(base(w));
    }
    None
}

/// python, python3, python3.12: not python3.12/ (a path part) and not
/// python. (a sentence).
fn is_python_program(word: &str) -> bool {
    let Some(rest) = word.strip_prefix("python") else {
        return false;
    };
    rest.bytes().all(|b| b.is_ascii_digit() || b == b'.') && !rest.ends_with('.')
}

fn shebang_is_python(line: &str) -> bool {
    shebang_program(line).is_some_and(|p| is_python_program(&p))
}

/// The line without its comment: `#` at the start or after a blank, outside
/// quotes opened on this line. `nf` also strips `//` after a blank or at the
/// start. A `#!` line is kept whole (the caller reads it first).
fn strip_line_comment(line: &str, nf: bool) -> &str {
    let b = line.as_bytes();
    let (mut single, mut double) = (false, false);
    let mut i = 0;
    while i < b.len() {
        let c = b[i];
        let after_blank = i == 0 || b[i - 1] == b' ' || b[i - 1] == b'\t';
        match c {
            b'\\' if !single => {
                i += 2;
                continue;
            }
            b'\'' if !double => single = !single,
            b'"' if !single => double = !double,
            b'#' if !single && !double && after_blank => return &line[..i],
            b'/' if nf && !single && !double && after_blank && b.get(i + 1) == Some(&b'/') => {
                return &line[..i];
            }
            _ => {}
        }
        i += 1;
    }
    line
}

/// The python program words of a code line, as calls.
fn python_words(code: &str) -> bool {
    let b = code.as_bytes();
    let mut from = 0;
    while let Some(k) = code[from..].find("python") {
        let s = from + k;
        from = s + 6;
        let before_ok = s == 0 || {
            let p = b[s - 1];
            p.is_ascii_whitespace() || b"\"'`(;|&{/".contains(&p)
        };
        if !before_ok {
            continue;
        }
        let mut e = s + 6;
        while e < b.len() && (b[e].is_ascii_digit() || b[e] == b'.') {
            e += 1;
        }
        let after_ok = e == b.len() || {
            let n = b[e];
            n.is_ascii_whitespace() || b"\"'`;|&)}".contains(&n)
        };
        if after_ok && is_python_program(&code[s..e]) {
            return true;
        }
    }
    false
}

/// Does this code line open a heredoc whose delimiter says Python?
fn python_heredoc(code: &str) -> bool {
    let mut rest = code;
    while let Some(k) = rest.find("<<") {
        rest = &rest[k + 2..];
        let t = rest.trim_start_matches(['-', '~']).trim_start();
        let t = t.trim_start_matches(['\'', '"']);
        let word: String = t
            .chars()
            .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
            .collect();
        let up = word.to_ascii_uppercase();
        if up.starts_with("PY") || up.contains("PYTHON") || up == "EOPY" {
            return true;
        }
    }
    false
}

fn excerpt(line: &str) -> String {
    let t = line.trim();
    if t.chars().count() <= 120 {
        t.to_owned()
    } else {
        let cut: String = t.chars().take(117).collect();
        format!("{cut}...")
    }
}

/// The findings in one file's text.
pub fn scan_text(rel: &str, text: &str) -> Vec<Finding> {
    let mut out = Vec::new();
    let name = rel.rsplit('/').next().unwrap_or(rel);
    if [".py", ".pyw", ".pyi"].iter().any(|e| name.ends_with(e)) {
        out.push(Finding {
            file: rel.to_owned(),
            line: 0,
            kind: "file",
            text: "a Python file".into(),
        });
    }
    let first = text.lines().next().unwrap_or("");
    if shebang_is_python(first) {
        out.push(Finding {
            file: rel.to_owned(),
            line: 1,
            kind: "shebang",
            text: excerpt(first),
        });
    }
    let push = |out: &mut Vec<Finding>, n: usize, kind, line: &str| {
        out.push(Finding {
            file: rel.to_owned(),
            line: n,
            kind,
            text: excerpt(line),
        });
    };
    match lang(rel, text) {
        Lang::Other => {}
        Lang::Rust => scan_rust(rel, text, &mut out),
        l @ (Lang::Shell | Lang::Nextflow | Lang::CommandTable) => {
            let nf = l == Lang::Nextflow;
            let mut block = false;
            for (i, line) in text.lines().enumerate() {
                let n = i + 1;
                if n > 1 && shebang_is_python(line) {
                    push(&mut out, n, "shebang", line);
                    continue;
                }
                if l == Lang::CommandTable && line.starts_with('#') {
                    continue;
                }
                let mut code = line.to_owned();
                if nf {
                    code = strip_block_comments(&code, &mut block);
                }
                let code = strip_line_comment(&code, nf);
                if python_heredoc(code) {
                    push(&mut out, n, "heredoc", line);
                } else if python_words(code) {
                    push(&mut out, n, "call", line);
                }
            }
        }
    }
    out
}

/// A Nextflow line without its `/* */` parts; `open` carries a block that
/// spans lines.
fn strip_block_comments(line: &str, open: &mut bool) -> String {
    let mut out = String::new();
    let mut rest = line;
    loop {
        if *open {
            match rest.find("*/") {
                Some(k) => {
                    rest = &rest[k + 2..];
                    *open = false;
                }
                None => return out,
            }
        } else {
            // Only a `/*` outside quotes on this line opens a block.
            match block_start(rest) {
                Some(k) => {
                    out.push_str(&rest[..k]);
                    rest = &rest[k + 2..];
                    *open = true;
                }
                None => {
                    out.push_str(rest);
                    return out;
                }
            }
        }
    }
}

fn block_start(s: &str) -> Option<usize> {
    let b = s.as_bytes();
    let (mut single, mut double) = (false, false);
    let mut i = 0;
    while i + 1 < b.len() {
        match b[i] {
            b'\\' => i += 1,
            b'\'' if !double => single = !single,
            b'"' if !single => double = !double,
            b'/' if !single && !double && b[i + 1] == b'*' => return Some(i),
            _ => {}
        }
        i += 1;
    }
    None
}

/// What each byte of Rust source is.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Part {
    Code,
    Str,
    Comment,
}

/// Marks each byte of Rust source as code, string or comment: line and
/// block comments (nested), strings with escapes, raw strings, byte
/// strings, and char literals told from lifetimes.
fn rust_parts(src: &str) -> Vec<Part> {
    let b = src.as_bytes();
    let mut part = vec![Part::Code; b.len()];
    let mut i = 0;
    while i < b.len() {
        let c = b[i];
        if c == b'/' && b.get(i + 1) == Some(&b'/') {
            while i < b.len() && b[i] != b'\n' {
                part[i] = Part::Comment;
                i += 1;
            }
        } else if c == b'/' && b.get(i + 1) == Some(&b'*') {
            let mut depth = 0usize;
            while i < b.len() {
                if b[i] == b'/' && b.get(i + 1) == Some(&b'*') {
                    depth += 1;
                    part[i] = Part::Comment;
                    part[i + 1] = Part::Comment;
                    i += 2;
                } else if b[i] == b'*' && b.get(i + 1) == Some(&b'/') {
                    depth -= 1;
                    part[i] = Part::Comment;
                    part[i + 1] = Part::Comment;
                    i += 2;
                    if depth == 0 {
                        break;
                    }
                } else {
                    part[i] = Part::Comment;
                    i += 1;
                }
            }
        } else if (c == b'r' || (c == b'b' && b.get(i + 1) == Some(&b'r')))
            && (i == 0 || !(b[i - 1].is_ascii_alphanumeric() || b[i - 1] == b'_'))
            && {
                let j = if c == b'b' { i + 2 } else { i + 1 };
                let mut k = j;
                while b.get(k) == Some(&b'#') {
                    k += 1;
                }
                b.get(k) == Some(&b'"')
            }
        {
            let j = if c == b'b' { i + 2 } else { i + 1 };
            let mut k = j;
            while b[k] == b'#' {
                k += 1;
            }
            let hashes = k - j;
            let close: Vec<u8> = std::iter::once(b'"')
                .chain(std::iter::repeat_n(b'#', hashes))
                .collect();
            // The prefix and the opening quote are code, the body a string.
            let body = k + 1;
            let end = b[body..]
                .windows(close.len())
                .position(|w| w == close.as_slice())
                .map_or(b.len(), |p| body + p);
            for p in &mut part[body..end] {
                *p = Part::Str;
            }
            i = (end + close.len()).min(b.len());
        } else if c == b'"' {
            i += 1;
            while i < b.len() && b[i] != b'"' {
                if b[i] == b'\\' {
                    part[i] = Part::Str;
                    i += 1;
                }
                if i < b.len() {
                    part[i] = Part::Str;
                    i += 1;
                }
            }
            i += 1;
        } else if c == b'\'' {
            // 'x', '\n', '\u{..}' are chars; 'a in &'a str is a lifetime.
            if b.get(i + 1) == Some(&b'\\') {
                let mut k = i + 2;
                while k < b.len() && b[k] != b'\'' {
                    k += 1;
                }
                i = k + 1;
            } else if b.get(i + 2) == Some(&b'\'') {
                i += 3;
            } else if let Some(ch) = src[i + 1..].chars().next() {
                let w = ch.len_utf8();
                i += if b.get(i + 1 + w) == Some(&b'\'') {
                    2 + w
                } else {
                    1
                };
            } else {
                i += 1;
            }
        } else {
            i += 1;
        }
    }
    part
}

fn scan_rust(rel: &str, src: &str, out: &mut Vec<Finding>) {
    let part = rust_parts(src);
    let line_of = |pos: usize| src[..pos].bytes().filter(|&c| c == b'\n').count() + 1;
    let line_text = |pos: usize| {
        let s = src[..pos].rfind('\n').map_or(0, |k| k + 1);
        let e = src[pos..].find('\n').map_or(src.len(), |k| pos + k);
        excerpt(&src[s..e])
    };
    let mut seen = BTreeSet::new();
    let mut hit = |pos: usize, out: &mut Vec<Finding>| {
        let n = line_of(pos);
        if seen.insert(n) {
            out.push(Finding {
                file: rel.to_owned(),
                line: n,
                kind: "rust",
                text: line_text(pos),
            });
        }
    };
    // Command::new("…python3") with the call in code and the program a
    // string literal.
    let needle = "Command::new(";
    let mut from = 0;
    while let Some(k) = src[from..].find(needle) {
        let s = from + k;
        from = s + needle.len();
        if part[s] != Part::Code {
            continue;
        }
        let rest = &src[from..];
        let skip = rest.len() - rest.trim_start().len();
        let q = from + skip;
        if src.as_bytes().get(q) != Some(&b'"') || part[q] != Part::Code {
            continue;
        }
        let Some(len) = src[q + 1..].find('"') else {
            continue;
        };
        let prog = &src[q + 1..q + 1 + len];
        let base = prog.rsplit('/').next().unwrap_or(prog);
        if is_python_program(base) {
            hit(s, out);
        }
    }
    // The oracle variable, in code or in a string (env::var("…")).
    let var = oracle_var();
    let mut from = 0;
    while let Some(k) = src[from..].find(&var) {
        let s = from + k;
        from = s + var.len();
        let word = |c: Option<&u8>| c.is_some_and(|c| c.is_ascii_alphanumeric() || *c == b'_');
        if part[s] != Part::Comment
            && !word(s.checked_sub(1).and_then(|p| src.as_bytes().get(p)))
            && !word(src.as_bytes().get(from))
        {
            hit(s, out);
        }
    }
}

/// One row of the allow list.
#[derive(Debug)]
struct Row {
    file: String,
    reason: String,
    line: usize,
}

fn read_allow(path: &Path) -> Result<(Vec<Row>, Vec<String>), String> {
    let text =
        fs::read_to_string(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let mut rows = Vec::new();
    let mut bad = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() != 2 {
            bad.push(format!(
                "{}:{}: a row is file<TAB>reason, this one has {} columns",
                path.display(),
                i + 1,
                cols.len()
            ));
            continue;
        }
        rows.push(Row {
            file: cols[0].trim().to_owned(),
            reason: cols[1].trim().to_owned(),
            line: i + 1,
        });
    }
    Ok((rows, bad))
}

/// The files to scan: tracked, plus new files that are not ignored.
fn repo_files(root: &Path) -> Result<Vec<String>, String> {
    let out = Command::new("git")
        .arg("-C")
        .arg(root)
        .args([
            "ls-files",
            "-z",
            "--cached",
            "--others",
            "--exclude-standard",
        ])
        .output()
        .map_err(|e| format!("cannot run git ls-files in {}: {e}", root.display()))?;
    if !out.status.success() {
        return Err(format!(
            "git ls-files failed in {}: {}",
            root.display(),
            String::from_utf8_lossy(&out.stderr).trim()
        ));
    }
    let mut files: BTreeSet<String> = BTreeSet::new();
    for p in out.stdout.split(|&c| c == 0) {
        if p.is_empty() {
            continue;
        }
        files.insert(String::from_utf8_lossy(p).into_owned());
    }
    Ok(files.into_iter().collect())
}

/// What a run of the gate found.
pub struct Report {
    pub findings: Vec<Finding>,
    pub failures: Vec<String>,
    pub scanned: usize,
    pub code_files: usize,
    pub rust_files: usize,
    pub rows: usize,
}

/// Runs the gate on ROOT with the allow list ALLOW.
pub fn check(root: &Path, allow: &Path) -> Result<Report, String> {
    let files = repo_files(root)?;
    let mut findings = Vec::new();
    let (mut scanned, mut code_files, mut rust_files) = (0, 0, 0);
    for rel in &files {
        let p = root.join(rel);
        // A deleted file still in the index, a symlink, a submodule: skip.
        let Ok(meta) = fs::symlink_metadata(&p) else {
            continue;
        };
        if !meta.is_file() {
            continue;
        }
        let bytes = fs::read(&p).map_err(|e| format!("cannot read {rel}: {e}"))?;
        let text = String::from_utf8_lossy(&bytes);
        scanned += 1;
        match lang(rel, &text) {
            Lang::Rust => rust_files += 1,
            Lang::Shell | Lang::Nextflow | Lang::CommandTable => code_files += 1,
            Lang::Other => {}
        }
        findings.extend(scan_text(rel, &text));
    }
    if code_files == 0 || rust_files == 0 {
        return Err(format!(
            "{} holds {code_files} shell, Nextflow or workflow files and {rust_files} Rust files: \
             the gate read nothing it could check",
            root.display()
        ));
    }
    let (rows, mut failures) = read_allow(allow)?;
    let pins: BTreeSet<String> = versions_env::read(&root.join("versions.env"))?
        .into_iter()
        .map(|p| p.name)
        .collect();
    let mut by_file: BTreeMap<&str, Vec<&Finding>> = BTreeMap::new();
    for f in &findings {
        by_file.entry(f.file.as_str()).or_default().push(f);
    }
    let mut named: BTreeMap<String, usize> = BTreeMap::new();
    let al = allow.display();
    for r in &rows {
        let at = format!("{al}:{}", r.line);
        if let Some(first) = named.insert(r.file.clone(), r.line) {
            failures.push(format!(
                "{at}: {} is named twice (also line {first})",
                r.file
            ));
            continue;
        }
        if r.reason.is_empty() {
            failures.push(format!("{at}: {} has no reason", r.file));
        }
        if r.file.contains(['*', '?', '[']) || r.file.ends_with('/') {
            failures.push(format!("{at}: {} is a pattern; name one file", r.file));
            continue;
        }
        let p = root.join(&r.file);
        if p.is_dir() {
            failures.push(format!("{at}: {} is a directory; name one file", r.file));
            continue;
        }
        if !p.is_file() {
            failures.push(format!("{at}: {} does not exist", r.file));
            continue;
        }
        for var in image_vars(&r.reason) {
            if !pins.contains(&var) {
                failures.push(format!(
                    "{at}: the reason names {var}, which versions.env does not set"
                ));
            }
        }
        let fs_here = by_file.get(r.file.as_str()).cloned().unwrap_or_default();
        if fs_here.is_empty() {
            failures.push(format!(
                "{at}: {} holds no Python any more: delete the row",
                r.file
            ));
        }
        if fs_here
            .iter()
            .any(|f| f.kind == "file" || f.kind == "shebang")
        {
            failures.push(format!(
                "{at}: {} is a Python file, which is owned code: it cannot be allowed",
                r.file
            ));
        }
    }
    for f in &findings {
        if !named.contains_key(&f.file) {
            failures.push(format!(
                "{}:{}: {}: {} (not in {})",
                f.file, f.line, f.kind, f.text, ALLOW
            ));
        }
    }
    Ok(Report {
        findings,
        failures,
        scanned,
        code_files,
        rust_files,
        rows: rows.len(),
    })
}

/// The `*_IMAGE` variable names a reason mentions.
fn image_vars(reason: &str) -> Vec<String> {
    reason
        .split(|c: char| !(c.is_ascii_alphanumeric() || c == '_'))
        .filter(|w| {
            w.ends_with("_IMAGE")
                && w.len() > "_IMAGE".len()
                && w.bytes()
                    .all(|b| b.is_ascii_uppercase() || b.is_ascii_digit() || b == b'_')
        })
        .map(str::to_owned)
        .collect()
}

/// `og dev python-gate ARGS`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let mut root = PathBuf::from(".");
    let mut allow: Option<PathBuf> = None;
    let mut list = false;
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--root" | "--allow" => {
                let Some(v) = it.next() else {
                    let _ = write!(err, "{USAGE}");
                    let _ = writeln!(err, "og dev python-gate: error: {a} needs a value");
                    return 2;
                };
                if a == "--root" {
                    root = PathBuf::from(v);
                } else {
                    allow = Some(PathBuf::from(v));
                }
            }
            "--list" => list = true,
            "-h" | "--help" => {
                let _ = write!(out, "{USAGE}");
                return 0;
            }
            _ => {
                let _ = write!(err, "{USAGE}");
                let _ = writeln!(err, "og dev python-gate: error: unexpected argument '{a}'");
                return 2;
            }
        }
    }
    let allow = allow.unwrap_or_else(|| root.join(ALLOW));
    let r = match check(&root, &allow) {
        Ok(r) => r,
        Err(e) => {
            let _ = writeln!(err, "og dev python-gate: {e}");
            return 2;
        }
    };
    if list {
        for f in &r.findings {
            let _ = writeln!(out, "{}\t{}\t{}\t{}", f.file, f.line, f.kind, f.text);
        }
    }
    for f in &r.failures {
        let _ = writeln!(out, "FAIL {f}");
    }
    let files: BTreeSet<&str> = r.findings.iter().map(|f| f.file.as_str()).collect();
    if r.failures.is_empty() {
        let _ = writeln!(
            out,
            "OK: {} files read ({} shell, Nextflow or workflow, {} Rust); {} Python uses in {} files, \
             each third-party and named in {} ({} rows)",
            r.scanned,
            r.code_files,
            r.rust_files,
            r.findings.len(),
            files.len(),
            ALLOW,
            r.rows
        );
        0
    } else {
        let _ = writeln!(
            out,
            "{} problems: owned Python, or a row of {} that is wrong",
            r.failures.len(),
            ALLOW
        );
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn kinds(rel: &str, text: &str) -> Vec<(usize, &'static str)> {
        scan_text(rel, text)
            .into_iter()
            .map(|f| (f.line, f.kind))
            .collect()
    }

    #[test]
    fn python_files_and_shebangs() {
        assert_eq!(kinds("bin/x.py", "print(1)\n"), [(0, "file")]);
        assert_eq!(kinds("tests/y.pyi", ""), [(0, "file")]);
        assert_eq!(
            kinds("bin/tool", "#!/usr/bin/env python3\nprint(1)\n"),
            [(1, "shebang")]
        );
        assert_eq!(
            kinds("bin/tool", "#!/usr/bin/python3.12 -u\n"),
            [(1, "shebang")]
        );
        assert_eq!(
            kinds("bin/tool", "#!/usr/bin/env -S python3 -I\n"),
            [(1, "shebang")]
        );
        assert_eq!(
            kinds("bin/tool", "#!/usr/bin/env PYTHONPATH=x python\n"),
            [(1, "shebang")]
        );
        // The word counts anywhere in shell code, not only first on a line:
        // that is how `docker run IMG python3 x` is caught.
        assert_eq!(
            kinds("bin/tool", "#!/bin/sh\necho python3\n"),
            [(2, "call")]
        );
        assert!(kinds("bin/og", "#!/bin/sh\nexec og \"$@\"\n").is_empty());
        // A Nextflow script block that is a Python script.
        let nf = "process P {\n    script:\n    \"\"\"\n    #!/usr/bin/env python3\n    print(1)\n    \"\"\"\n}\n";
        assert_eq!(kinds("modules/local/p/main.nf", nf), [(4, "shebang")]);
    }

    #[test]
    fn shell_calls_count_and_prose_does_not() {
        let sh = "#!/usr/bin/env bash\n\
                  # No python3 is needed.\n\
                  echo hi  # python3 here is a comment\n\
                  python3 -c 'print(1)'\n\
                  out=$(python3 - <<'EOF'\nprint(2)\nEOF\n)\n\
                  PYTHONPATH=\"$1\" python3 -m cyrius\n\
                  docker run --rm \"$IMG\" /usr/bin/python3 x\n\
                  stamp=\"python=${PYTHON_IMAGE}\"\n\
                  image=python:3.15\n\
                  ls /usr/local/lib/python3.12/site-packages\n\
                  echo \"${#arr[@]}\" | python\n";
        assert_eq!(
            kinds("scripts/x.sh", sh),
            [
                (4, "call"),
                (5, "call"),
                (9, "call"),
                (10, "call"),
                (14, "call")
            ]
        );
    }

    #[test]
    fn heredocs_with_a_python_delimiter() {
        assert_eq!(
            kinds("s.sh", "cat > x <<'PY'\nprint(1)\nPY\n"),
            [(1, "heredoc")]
        );
        assert_eq!(kinds("s.sh", "run <<-EOF_PYTHON\n"), [(1, "heredoc")]);
        assert!(kinds("s.sh", "cat <<'EOF'\nhello\nEOF\n").is_empty());
        assert!(kinds("s.sh", "cat <<COPY\nx\nCOPY\n").is_empty());
    }

    #[test]
    fn nextflow_comments_are_prose() {
        let nf = "/*\n  CNVPYTOR (python3 successor)\n*/\nprocess P {\n    // python3 in a comment\n    \
                  container 'python:3.15'\n    script:\n    \"\"\"\n    # python3 in a shell comment\n    \
                  python3 /pharmcat/pre \\\\\n    \"\"\"\n}\n";
        assert_eq!(kinds("modules/local/p/main.nf", nf), [(10, "call")]);
    }

    #[test]
    fn workflows_and_command_tables() {
        let yml = "jobs:\n  a:\n    steps:\n      - uses: actions/setup-python@abc # v7\n        with:\n          \
                   python-version: '3.15'\n      # python -m x in a comment\n      - run: python -m zensical build\n";
        assert_eq!(kinds(".github/workflows/docs.yml", yml), [(8, "call")]);
        // A YAML file outside .github is data.
        assert!(kinds("zensical.yml", "run: python -m x\n").is_empty());
        let tsv = "# Rows, tab separated:\n#var\topts\tneeds\tcommand\texpect\n# python3 in a comment\nX_IMAGE\t-\t-\tpython -c 'import x'\tok\n";
        assert_eq!(kinds("tests/smoke/commands.tsv", tsv), [(4, "call")]);
        // A TSV without a command column is data.
        assert!(
            kinds(
                "scripts/ci/renovate-pins.tsv",
                "#dep\tfile\npython\tdocs.yml\n"
            )
            .is_empty()
        );
    }

    #[test]
    fn markdown_and_text_are_prose() {
        assert!(kinds("docs/x.md", "Run `python3 -c 'x'` first.\n").is_empty());
        assert!(kinds("notes.txt", "python3 x.py\n").is_empty());
    }

    #[test]
    fn rust_commands_and_the_oracle_variable() {
        let var = oracle_var();
        let rs = format!(
            "// Command::new(\"python3\") in a comment\n\
             /* {var} in a block */\n\
             fn a() {{ let _ = std::process::Command::new(\"python3\"); }}\n\
             fn b() {{ let _ = Command::new( \"/usr/bin/python3.12\" ); }}\n\
             fn c() {{ let _ = std::env::var(\"{var}\"); }}\n\
             fn d() {{ let _ = Command::new(\"og\"); let s = \"Command::new(\\\"python3\\\")\"; }}\n\
             fn e<'a>(x: &'a str) -> char {{ let _ = x; 'p' }}\n\
             fn f() {{ let _ = \"{var}S\"; let _ = r#\"Command::new(\"python3\")\"#; }}\n"
        );
        assert_eq!(
            kinds("rust/crates/x/src/a.rs", &rs),
            [(3, "rust"), (4, "rust"), (5, "rust")]
        );
    }

    #[test]
    fn image_vars_in_a_reason() {
        assert_eq!(
            image_vars("PHARMCAT_IMAGE: its preprocessor"),
            ["PHARMCAT_IMAGE"]
        );
        assert!(image_vars("the docs builder, no _IMAGE").is_empty());
    }
}
