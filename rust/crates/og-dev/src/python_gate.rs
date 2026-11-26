//! `og dev python-gate [--root DIR] [--allow FILE] [--list]`: no owned
//! Python is left in the repository.
//!
//! Every file `git ls-files` lists (tracked, plus new files that are not
//! ignored) is read, and these are the findings:
//!
//! - `file`: a `*.py`, `*.pyw`, `*.pyi` or `*.ipynb` file;
//! - `shebang`: a `#!` line whose program is python, python3 or
//!   pythonX.Y, directly or through `env`: the first line of any file, or
//!   any line of a shell, Nextflow or workflow file (a script a heredoc or a
//!   Nextflow `script:` block writes);
//! - `call`: a Python program word outside a comment in a shell script (by
//!   its suffix, by a shell `#!` line whatever the suffix, a Makefile, a
//!   `*.mk`, a justfile or a Dockerfile), a Nextflow file or config, a
//!   GitHub workflow or action (not its `name:` and `description:`
//!   labels), JavaScript or TypeScript, or the `command` field of a TSV
//!   command table (`tests/smoke/commands.tsv`): a `python3 -c`, a
//!   `python3 - <<EOF`, a `python3 tool` run in an image, a `PY=python3`
//!   or `${PYTHON:-python3}` that picks the interpreter. The program words
//!   are python, python3, pythonX.Y and the package runners pip, pip3,
//!   pipx, uv and uvx. A workflow step that `uses:` a `setup-python` or
//!   `setup-uv` action is a call too;
//! - `heredoc`: a heredoc whose delimiter says it holds Python (`<<'PY'`,
//!   `<<PYTHON`, `<<EOF_PYTHON`), not a `<<PY` inside quotes;
//! - `rust`: in Rust code outside comments, a `Command::new("python3")`
//!   (any python program, any path, an ordinary or raw string), the
//!   variable that named a Python oracle, or a string literal that the
//!   shell rules above find a call, heredoc or shebang in
//!   (`.arg("python3")`, `"sh -c 'python3 x'"`, `const PY = "python3"`).
//!
//! A comment, a doc sentence and a Markdown, text, JSON or image file are
//! prose: they never count. `python:3.15` (an image), `python=...` (a
//! stamp), `PYTHON_IMAGE`, `setup-python` outside a `uses:` line and a
//! `site-packages/python3.12/` path are not calls.
//!
//! Known limits: a Python program name built at run time
//! (`format!("python{}", 3)`, `["py", "thon"].concat()`) is not seen, and
//! YAML outside `.github/workflows/` and `.github/actions/` is data (a
//! `.pre-commit-config.yaml` `entry: python3`, a compose `command:`).
//! Review catches those.
//!
//! `tests/python-allow.tsv` names, one row per file, the files whose
//! findings are third-party Python: a tool run inside its own pinned image
//! (PharmCAT's preprocessor, eKLIPse, Cyrius and the docs builder Zensical
//! in `PYTHON_IMAGE`). Each row is `file<TAB>uses<TAB>reason[<TAB>marks]`,
//! where `uses` is the number of lines with a finding and `marks`, if
//! given, is a `|`-separated list of texts one of which every such line
//! must hold (`/pharmcat/pharmcat_vcf_preprocessor`), so a use swapped for
//! another at the same count fails too. The gate fails when
//!
//! - a finding is in a file no row names;
//! - a row names a missing file, a directory or a glob, or has no reason;
//! - a line with a finding holds none of its row's marks, or a mark is
//!   empty;
//! - a reason names an `*_IMAGE` variable that `versions.env` does not set,
//!   or sets to "";
//! - a row covers nothing (the Python left that file: delete the row);
//! - the file holds more or fewer uses than the row says (a new use is
//!   checked against the reason before the count changes);
//! - a row names a Python file (a `file` or `shebang` finding): every Python
//!   file in the repository is owned code, so none can be allowed, and an
//!   allowed file can never import a module from `bin/` or the demo
//!   builders;
//! - a file is named twice.
//!
//! A symlink is not read, but one named like a Python file is a `file`
//! finding. A path git lists that is not UTF-8 is an error.
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

/// The program words, each in two halves for the same reason.
const PYTHON: &str = concat!("pyt", "hon");
const PIP: &str = concat!("pi", "p");
const RUNNERS: [&str; 3] = [concat!("pi", "px"), concat!("u", "v"), concat!("u", "vx")];

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
    /// The whole line, trimmed: what a row's marks are looked for in.
    pub code: String,
}

/// How the gate reads a file's lines.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Lang {
    /// Shell: `#` comments.
    Shell,
    /// A GitHub workflow or action: shell, but a `name:` or `description:`
    /// line is a label, not code.
    Workflow,
    /// Nextflow and Groovy: `//` and `/* */` comments, and `#` comments in
    /// the shell of `script:` blocks.
    Nextflow,
    /// A TSV command table: `#` lines are comments, the rest is shell.
    CommandTable,
    /// JavaScript and TypeScript: `//` and `/* */` comments.
    Script,
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
        "sh" | "bash" | "ksh" | "zsh" | "mk" => return Lang::Shell,
        "nf" | "config" | "groovy" => return Lang::Nextflow,
        "rs" => return Lang::Rust,
        "js" | "mjs" | "cjs" | "jsx" | "ts" | "mts" | "cts" | "tsx" => return Lang::Script,
        // Workflows and actions run code; the rest of .github (Dependabot,
        // issue forms) is data.
        "yml" | "yaml"
            if rel.starts_with(".github/workflows/") || rel.starts_with(".github/actions/") =>
        {
            return Lang::Workflow;
        }
        "tsv" if is_command_table(text) => return Lang::CommandTable,
        _ => {}
    }
    let make_or_just = matches!(
        name,
        "Makefile" | "makefile" | "GNUmakefile" | "justfile" | "Justfile" | ".justfile"
    );
    if name.starts_with("Dockerfile") || name.starts_with("Containerfile") || make_or_just {
        return Lang::Shell;
    }
    // A shell shebang makes a shell script, whatever the name's suffix
    // (`bin/report.local`).
    let shell_script = shebang_program(first_line)
        .is_some_and(|p| matches!(p.as_str(), "sh" | "bash" | "dash" | "ksh" | "zsh"));
    if shell_script {
        return Lang::Shell;
    }
    Lang::Other
}

/// Does this TSV have a `command` column?
fn is_command_table(text: &str) -> bool {
    command_column(text).is_some()
}

/// The index of a TSV's `command` column. Its header is one of the leading
/// `#` lines (`#var<TAB>opts<TAB>needs<TAB>command<TAB>expect`) or the
/// first line after them.
fn command_column(text: &str) -> Option<usize> {
    for line in text.lines().take(20) {
        let col = line
            .trim_start_matches('#')
            .split('\t')
            .position(|c| c.trim() == "command");
        if col.is_some() || !line.starts_with('#') {
            return col;
        }
    }
    None
}

/// A workflow line that is a label: `name:` or `description:`, also as
/// the first key of a list item.
fn yaml_label(line: &str) -> bool {
    let t = line.trim_start();
    let t = t.strip_prefix("- ").unwrap_or(t).trim_start();
    t.starts_with("name:") || t.starts_with("description:")
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
    word.strip_prefix(PYTHON).is_some_and(version_suffix)
}

/// "", "3" or "3.12": not "3." (a sentence).
fn version_suffix(rest: &str) -> bool {
    rest.bytes().all(|b| b.is_ascii_digit() || b == b'.') && !rest.ends_with('.')
}

/// A package runner that runs Python: pip, pip3, pip3.12, pipx, uv, uvx.
fn is_runner_program(word: &str) -> bool {
    RUNNERS.contains(&word) || word.strip_prefix(PIP).is_some_and(version_suffix)
}

fn shebang_is_python(line: &str) -> bool {
    shebang_program(line).is_some_and(|p| is_python_program(&p))
}

/// Does a `#` at `i` start a shell comment? It does at the start of a word:
/// at the start of the line, after a blank, or after `;`, `&` or `|`.
fn comment_start(b: &[u8], i: usize) -> bool {
    i == 0 || matches!(b[i - 1], b' ' | b'\t' | b';' | b'&' | b'|')
}

/// A shell line without its comment, outside quotes opened on this line.
fn strip_line_comment(line: &str) -> &str {
    let b = line.as_bytes();
    let (mut single, mut double) = (false, false);
    let mut i = 0;
    while i < b.len() {
        match b[i] {
            b'\\' if !single => {
                i += 2;
                continue;
            }
            b'\'' if !double => single = !single,
            b'"' if !single => double = !double,
            b'#' if !single && !double && comment_start(b, i) => return &line[..i],
            _ => {}
        }
        i += 1;
    }
    line
}

/// A Nextflow line without its comments, read in order: `/* */` (`open`
/// carries a block that spans lines), `//` at the start or after a blank,
/// and, with `hash`, a `#` shell comment in a `script:` block (JavaScript
/// has none). Only a marker outside quotes counts, and a `/*` after `//` is part of the line comment.
fn strip_nf_comments(line: &str, open: &mut bool, hash: bool) -> String {
    let b = line.as_bytes();
    let mut out = String::new();
    let (mut single, mut double) = (false, false);
    let (mut i, mut seg) = (0, 0);
    while i < b.len() {
        if *open {
            if b[i] == b'*' && b.get(i + 1) == Some(&b'/') {
                *open = false;
                seg = i + 2;
                i += 2;
            } else {
                i += 1;
            }
            continue;
        }
        let quoted = single || double;
        let after_blank = i == 0 || b[i - 1] == b' ' || b[i - 1] == b'\t';
        match b[i] {
            b'\\' if !single => {
                i += 2;
                continue;
            }
            b'\'' if !double => single = !single,
            b'"' if !single => double = !double,
            b'/' if !quoted && b.get(i + 1) == Some(&b'*') => {
                out.push_str(&line[seg..i]);
                *open = true;
                i += 2;
                continue;
            }
            b'/' if !quoted && after_blank && b.get(i + 1) == Some(&b'/') => {
                out.push_str(&line[seg..i]);
                return out;
            }
            b'#' if hash && !quoted && comment_start(b, i) => {
                out.push_str(&line[seg..i]);
                return out;
            }
            _ => {}
        }
        i += 1;
    }
    if !*open {
        out.push_str(&line[seg.min(line.len())..]);
    }
    out
}

/// Can a program word start at `s`? At the line start, after a blank, a
/// quote, a bracket, a shell operator, a path `/`, a list `,` or an
/// assignment `=` (`PY=python3`), and after the `-` or `+` of
/// `${NAME-…}`, `${NAME:-…}`, `${NAME+…}` and `${NAME:+…}`.
fn word_start(b: &[u8], s: usize) -> bool {
    if s == 0 {
        return true;
    }
    let p = b[s - 1];
    if p.is_ascii_whitespace() || b"\"'`(;|&{/[,=".contains(&p) {
        return true;
    }
    (p == b'-' || p == b'+') && {
        let head = &b[..s - 1];
        let head = head.strip_suffix(b":").unwrap_or(head);
        let name = head
            .iter()
            .rev()
            .take_while(|c| c.is_ascii_alphanumeric() || **c == b'_')
            .count();
        name > 0 && head[..head.len() - name].ends_with(b"${")
    }
}

/// The Python program words of a code line (python, python3, pythonX.Y,
/// pip, pip3, pipx, uv, uvx), as calls.
fn python_words(code: &str) -> bool {
    let b = code.as_bytes();
    for s in 0..b.len() {
        if !matches!(b[s], b'p' | b'u') || !word_start(b, s) {
            continue;
        }
        let mut e = s;
        while e < b.len() && (b[e].is_ascii_lowercase() || b[e].is_ascii_digit() || b[e] == b'.') {
            e += 1;
        }
        // A redirection may follow at once (`python3</dev/null`); a `/`
        // may not (`site-packages/python3.12/`). A `\` ends a word only
        // before an escaped quote (`"sh -c \"python3\""`), not in `\n` or `\.`.
        let after_ok = e == b.len() || {
            let n = b[e];
            n.is_ascii_whitespace()
                || b"\"'`;|&)}<>,".contains(&n)
                || (n == b'\\' && matches!(b.get(e + 1), Some(b'"' | b'\'')))
        };
        // `{python}`, `${python}`: a variable or a format placeholder.
        let braced = s > 0 && b[s - 1] == b'{' && b.get(e) == Some(&b'}');
        let word = &code[s..e];
        if after_ok && !braced && (is_python_program(word) || is_runner_program(word)) {
            return true;
        }
    }
    false
}

/// A workflow step that installs a Python: `uses:` a `setup-python` or
/// `setup-uv` action, at any owner and ref.
fn setup_python(code: &str) -> bool {
    let t = code.trim_start();
    let t = t.strip_prefix("- ").unwrap_or(t).trim_start();
    t.strip_prefix("uses:").is_some_and(|r| {
        let action = r.trim().trim_matches(['"', '\'']);
        let action = action.split('@').next().unwrap_or("");
        let name = action.rsplit('/').next().unwrap_or("");
        name == ["setup-", PYTHON].concat() || name == ["setup-", RUNNERS[1]].concat()
    })
}

/// Does this code line open a heredoc whose delimiter says Python? Only a
/// `<<` outside quotes is an operator (`printf 'x <<PY'` prints text); the
/// delimiter after it may be quoted (`<<'PY'`). `<<<` is a here-string. A
/// Groovy `"""` or `'''` opens a string whose body is shell, so it does not
/// count as a quote.
fn python_heredoc(code: &str) -> bool {
    let b = code.as_bytes();
    let (mut single, mut double) = (false, false);
    let mut i = 0;
    while i < b.len() {
        if b[i..].starts_with(b"\"\"\"") || b[i..].starts_with(b"'''") {
            i += 3;
            continue;
        }
        match b[i] {
            b'\\' if !single => {
                i += 2;
                continue;
            }
            b'\'' if !double => single = !single,
            b'"' if !single => double = !double,
            b'<' if !single && !double && b.get(i + 1) == Some(&b'<') => {
                if b.get(i + 2) == Some(&b'<') {
                    i += 3;
                    continue;
                }
                let t = code[i + 2..].trim_start_matches(['-', '~']).trim_start();
                let t = t.trim_start_matches(['\'', '"']);
                let word: String = t
                    .chars()
                    .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
                    .collect();
                let up = word.to_ascii_uppercase();
                if up.starts_with("PY") || up.contains("PYTHON") || up == "EOPY" {
                    return true;
                }
                i += 2;
                continue;
            }
            _ => {}
        }
        i += 1;
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

/// A finding on line `line` (1-based, 0 for the file) of `rel`.
fn finding(rel: &str, line: usize, kind: &'static str, text: &str) -> Finding {
    Finding {
        file: rel.to_owned(),
        line,
        kind,
        text: excerpt(text),
        code: text.trim().to_owned(),
    }
}

/// The `file` finding of a path whose name is a Python file's.
fn python_name(rel: &str) -> Option<Finding> {
    let name = rel.rsplit('/').next().unwrap_or(rel);
    [".py", ".pyw", ".pyi", ".ipynb"]
        .iter()
        .any(|e| name.ends_with(e))
        .then(|| finding(rel, 0, "file", "a Python file"))
}

/// The findings in one file's text.
pub fn scan_text(rel: &str, text: &str) -> Vec<Finding> {
    let mut out: Vec<Finding> = python_name(rel).into_iter().collect();
    let first = text.lines().next().unwrap_or("");
    if shebang_is_python(first) {
        out.push(finding(rel, 1, "shebang", first));
    }
    let push = |out: &mut Vec<Finding>, n: usize, kind, line: &str| {
        out.push(finding(rel, n, kind, line));
    };
    match lang(rel, text) {
        Lang::Other => {}
        Lang::Rust => scan_rust(rel, text, &mut out),
        l @ (Lang::Shell | Lang::Workflow | Lang::Nextflow | Lang::CommandTable | Lang::Script) => {
            // `//` and `/* */` comments; Nextflow's script blocks also have
            // shell `#` comments, JavaScript has none.
            let slashes = matches!(l, Lang::Nextflow | Lang::Script);
            // A command table: only its `command` field is shell.
            let column = if l == Lang::CommandTable {
                command_column(text)
            } else {
                None
            };
            let mut block = false;
            for (i, line) in text.lines().enumerate() {
                let n = i + 1;
                if l == Lang::CommandTable && line.starts_with('#') {
                    continue;
                }
                if l == Lang::Workflow && yaml_label(line) {
                    continue;
                }
                let field = match column {
                    Some(c) => line.split('\t').nth(c).unwrap_or(""),
                    None => line,
                };
                // A `#!` line inside a Nextflow block comment is prose.
                if n > 1 && !block && shebang_is_python(field) {
                    push(&mut out, n, "shebang", line);
                    continue;
                }
                let stripped;
                let code = if slashes {
                    stripped = strip_nf_comments(field, &mut block, l == Lang::Nextflow);
                    stripped.as_str()
                } else {
                    strip_line_comment(field)
                };
                if python_heredoc(code) {
                    push(&mut out, n, "heredoc", line);
                } else if python_words(code) || (l == Lang::Workflow && setup_python(code)) {
                    push(&mut out, n, "call", line);
                }
            }
        }
    }
    out
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
        &src[s..e]
    };
    let mut seen = BTreeSet::new();
    let mut hit = |pos: usize, out: &mut Vec<Finding>| {
        let n = line_of(pos);
        if seen.insert(n) {
            out.push(finding(rel, n, "rust", line_text(pos)));
        }
    };
    // Command::new("…python3") with the call in code and the program a
    // string literal, ordinary or raw (`r"…"`, `r#"…"#`), after any blanks
    // and comments.
    let b = src.as_bytes();
    let needle = "Command::new(";
    let mut from = 0;
    while let Some(k) = src[from..].find(needle) {
        let s = from + k;
        from = s + needle.len();
        if part[s] != Part::Code {
            continue;
        }
        let mut q = from;
        while q < b.len() && (b[q].is_ascii_whitespace() || part[q] == Part::Comment) {
            q += 1;
        }
        let mut open = q;
        if b.get(open) == Some(&b'r') {
            open += 1;
            while b.get(open) == Some(&b'#') {
                open += 1;
            }
        }
        if b.get(open) != Some(&b'"') || part[q] != Part::Code {
            continue;
        }
        let body = open + 1;
        let mut end = body;
        while end < b.len() && part[end] == Part::Str {
            end += 1;
        }
        let prog = &src[body..end];
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
    // Every string literal, read with the shell rules on each of its lines
    // (a real line break or a `\n` escape): a program word (`.arg("python3")`,
    // `"sh -c 'python3 x'"`), a Python heredoc or a python `#!` line.
    let mut i = 0;
    while i < b.len() {
        if part[i] != Part::Str {
            i += 1;
            continue;
        }
        let start = i;
        while i < b.len() && part[i] == Part::Str {
            i += 1;
        }
        let mut at = start;
        for real in src[start..i].split('\n') {
            let mut piece_at = at;
            for piece in real.split("\\n") {
                let code = strip_line_comment(piece);
                if shebang_is_python(piece) || python_heredoc(code) || python_words(code) {
                    hit(piece_at, out);
                }
                piece_at += piece.len() + 2;
            }
            at += real.len() + 1;
        }
    }
    out.sort_by_key(|f| f.line);
}

/// One row of the allow list.
#[derive(Debug)]
struct Row {
    file: String,
    /// How many Python uses (lines with a finding) the file holds.
    uses: usize,
    reason: String,
    /// Texts one of which every line with a finding must hold; empty: any.
    marks: Vec<String>,
    line: usize,
}

fn read_allow(path: &Path, label: &str) -> Result<(Vec<Row>, Vec<String>), String> {
    let text =
        fs::read_to_string(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let mut rows = Vec::new();
    let mut bad = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if !(3..=4).contains(&cols.len()) {
            bad.push(format!(
                "{label}:{}: a row is file<TAB>uses<TAB>reason[<TAB>marks], this one has {} columns",
                i + 1,
                cols.len()
            ));
            continue;
        }
        let Some(uses) = cols[1].trim().parse::<usize>().ok().filter(|&n| n > 0) else {
            bad.push(format!(
                "{label}:{}: uses is the number of Python uses in {}, a whole number above 0, not '{}'",
                i + 1,
                cols[0].trim(),
                cols[1].trim()
            ));
            continue;
        };
        rows.push(Row {
            file: cols[0].trim().to_owned(),
            uses,
            reason: cols[2].trim().to_owned(),
            marks: cols
                .get(3)
                .map(|m| m.split('|').map(|s| s.trim().to_owned()).collect())
                .unwrap_or_default(),
            line: i + 1,
        });
    }
    Ok((rows, bad))
}

/// A path from `git ls-files -z`, which must be UTF-8: a path the gate
/// cannot name exactly is an error, never a replaced byte.
fn git_path(p: &[u8]) -> Result<String, String> {
    String::from_utf8(p.to_vec()).map_err(|_| {
        format!(
            "git lists a path that is not UTF-8 ({}): rename it",
            String::from_utf8_lossy(p)
        )
    })
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
        files.insert(git_path(p)?);
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
    /// The allow list, as the messages name it.
    pub allow: String,
}

/// Runs the gate on ROOT with the allow list ALLOW.
pub fn check(root: &Path, allow: &Path) -> Result<Report, String> {
    let al = if allow == root.join(ALLOW) {
        ALLOW.to_owned()
    } else {
        allow.display().to_string()
    };
    let files = repo_files(root)?;
    let mut findings = Vec::new();
    let (mut scanned, mut code_files, mut rust_files) = (0, 0, 0);
    for rel in &files {
        let p = root.join(rel);
        // A deleted file still in the index, a submodule: skip. A symlink's
        // target is not read, but a symlink named like a Python file is one.
        let Ok(meta) = fs::symlink_metadata(&p) else {
            continue;
        };
        if !meta.is_file() {
            if meta.file_type().is_symlink() {
                findings.extend(python_name(rel));
            }
            continue;
        }
        let bytes = fs::read(&p).map_err(|e| format!("cannot read {rel}: {e}"))?;
        let text = String::from_utf8_lossy(&bytes);
        scanned += 1;
        match lang(rel, &text) {
            Lang::Rust => rust_files += 1,
            Lang::Shell | Lang::Workflow | Lang::Nextflow | Lang::CommandTable | Lang::Script => {
                code_files += 1
            }
            Lang::Other => {}
        }
        findings.extend(scan_text(rel, &text));
    }
    if code_files == 0 || rust_files == 0 {
        return Err(format!(
            "{} holds {code_files} shell, Nextflow, workflow or JavaScript files and {rust_files} Rust files: \
             the gate read nothing it could check",
            root.display()
        ));
    }
    let (rows, mut failures) = read_allow(allow, &al)?;
    // A variable set to "" names no image.
    let pins: BTreeSet<String> = versions_env::read(&root.join("versions.env"))?
        .into_iter()
        .filter(|p| !p.value.trim().is_empty())
        .map(|p| p.name)
        .collect();
    let mut by_file: BTreeMap<&str, Vec<&Finding>> = BTreeMap::new();
    for f in &findings {
        by_file.entry(f.file.as_str()).or_default().push(f);
    }
    let mut named: BTreeMap<String, usize> = BTreeMap::new();
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
        } else if fs_here.len() != r.uses {
            // A use added to (or gone from) an allowed file is reviewed
            // against the row's reason before the count is changed.
            failures.push(format!(
                "{at}: {} holds {} Python uses, the row says {}: check each use against the reason, then set the count",
                r.file,
                fs_here.len(),
                r.uses
            ));
        }
        if r.marks.iter().any(String::is_empty) {
            failures.push(format!(
                "{at}: {} has an empty mark: marks are texts separated by |",
                r.file
            ));
        } else if !r.marks.is_empty() {
            // A use swapped for another at the same count still fails.
            for f in &fs_here {
                if !r.marks.iter().any(|m| f.code.contains(m.as_str())) {
                    failures.push(format!(
                        "{at}: {}:{} holds none of the row's marks ({}): {}",
                        f.file,
                        f.line,
                        r.marks.join("|"),
                        f.text
                    ));
                }
            }
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
                "{}:{}: {}: {} (not in {al})",
                f.file, f.line, f.kind, f.text
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
        allow: al,
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
            "OK: {} files read ({} shell, Nextflow, workflow or JavaScript, {} Rust); {} Python uses in {} files, \
             each third-party and named in {} ({} rows)",
            r.scanned,
            r.code_files,
            r.rust_files,
            r.findings.len(),
            files.len(),
            r.allow,
            r.rows
        );
        0
    } else {
        let _ = writeln!(
            out,
            "{} problems: owned Python, or a row of {} that is wrong",
            r.failures.len(),
            r.allow
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
        assert_eq!(k("bin/x.py", "print(1)\n"), [(0, "file")]);
        assert_eq!(k("tests/y.pyi", ""), [(0, "file")]);
        assert_eq!(
            k("bin/tool", "#!/usr/bin/env @p3\nprint(1)\n"),
            [(1, "shebang")]
        );
        assert_eq!(k("bin/tool", "#!/usr/bin/@p3.12 -u\n"), [(1, "shebang")]);
        assert_eq!(
            k("bin/tool", "#!/usr/bin/env -S @p3 -I\n"),
            [(1, "shebang")]
        );
        assert_eq!(
            k("bin/tool", "#!/usr/bin/env @PTHONPATH=x @p\n"),
            [(1, "shebang")]
        );
        // The word counts anywhere in shell code, not only first on a line:
        // that is how `docker run IMG python3 x` is caught.
        assert_eq!(k("bin/tool", "#!/bin/sh\necho @p3\n"), [(2, "call")]);
        assert!(k("bin/og", "#!/bin/sh\nexec og \"$@\"\n").is_empty());
        // A Nextflow script block that is a Python script.
        let nf = "process P {\n    script:\n    \"\"\"\n    #!/usr/bin/env @p3\n    print(1)\n    \"\"\"\n}\n";
        assert_eq!(k("modules/local/p/main.nf", nf), [(4, "shebang")]);
    }

    #[test]
    fn shell_calls_count_and_prose_does_not() {
        let sh = "#!/usr/bin/env bash\n\
                  # No @p3 is needed.\n\
                  echo hi  # @p3 here is a comment\n\
                  @p3 -c 'print(1)'\n\
                  out=$(@p3 - <<'EOF'\nprint(2)\nEOF\n)\n\
                  @PTHONPATH=\"$1\" @p3 -m cyrius\n\
                  docker run --rm \"$IMG\" /usr/bin/@p3 x\n\
                  stamp=\"@p=${@PTHON_IMAGE}\"\n\
                  image=@p:3.15\n\
                  ls /usr/local/lib/@p3.12/site-packages\n\
                  echo \"${#arr[@]}\" | @p\n";
        assert_eq!(
            k("scripts/x.sh", sh),
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
            k("s.sh", "cat > x <<'@P'\nprint(1)\n@P\n"),
            [(1, "heredoc")]
        );
        assert_eq!(k("s.sh", "run <<-EOF_@PTHON\n"), [(1, "heredoc")]);
        assert!(k("s.sh", "cat <<'EOF'\nhello\nEOF\n").is_empty());
        assert!(k("s.sh", "cat <<CO@P\nx\nCO@P\n").is_empty());
    }

    #[test]
    fn nextflow_comments_are_prose() {
        let nf = "/*\n  CNV@PTOR (@p3 successor)\n*/\nprocess P {\n    // @p3 in a comment\n    \
                  container '@p:3.15'\n    script:\n    \"\"\"\n    # @p3 in a shell comment\n    \
                  @p3 /pharmcat/pre \\\\\n    \"\"\"\n}\n";
        assert_eq!(k("modules/local/p/main.nf", nf), [(10, "call")]);
    }

    #[test]
    fn workflows_and_command_tables() {
        let yml = "jobs:\n  a:\n    steps:\n      - uses: actions/setup-@p@abc # v7\n        with:\n          \
                   @p-version: '3.15'\n      # @p -m x in a comment\n      - run: @p -m zensical build\n";
        assert_eq!(
            k(".github/workflows/docs.yml", yml),
            [(4, "call"), (8, "call")]
        );
        // A YAML file outside .github is data.
        assert!(k("zensical.yml", "run: @p -m x\n").is_empty());
        let tsv = "# Rows, tab separated:\n#var\topts\tneeds\tcommand\texpect\n# @p3 in a comment\nX_IMAGE\t-\t-\t@p -c 'import x'\tok\n";
        assert_eq!(k("tests/smoke/commands.tsv", tsv), [(4, "call")]);
        // A TSV without a command column is data.
        assert!(k("scripts/ci/renovate-pins.tsv", "#dep\tfile\n@p\tdocs.yml\n").is_empty());
    }

    #[test]
    fn markdown_and_text_are_prose() {
        assert!(k("docs/x.md", "Run `@p3 -c 'x'` first.\n").is_empty());
        assert!(k("notes.txt", "@p3 x.py\n").is_empty());
    }

    #[test]
    fn rust_commands_and_the_oracle_variable() {
        let var = oracle_var();
        let rs = format!(
            "// Command::new(\"@p3\") in a comment\n\
             /* {var} in a block */\n\
             fn a() {{ let _ = std::process::Command::new(\"@p3\"); }}\n\
             fn b() {{ let _ = Command::new( \"/usr/bin/@p3.12\" ); }}\n\
             fn c() {{ let _ = std::env::var(\"{var}\"); }}\n\
             fn d() {{ let _ = Command::new(\"og\"); let s = \"Command::new(\\\"@p3\\\")\"; }}\n\
             fn e<'a>(x: &'a str) -> char {{ let _ = x; 'p' }}\n\
             fn f() {{ let _ = \"{var}S\"; let _ = r#\"Command::new(\"@p3\")\"#; }}\n"
        );
        // A string that holds a python word counts too (lines 6 and 8),
        // even when it is Rust source in a string.
        assert_eq!(
            k("rust/crates/x/src/a.rs", &rs),
            [
                (3, "rust"),
                (4, "rust"),
                (5, "rust"),
                (6, "rust"),
                (8, "rust")
            ]
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

    #[test]
    fn shell_shebangs_comments_heredocs_and_redirections() {
        // A shell #! line makes a shell script, whatever the suffix.
        assert_eq!(k("bin/report.local", "#!/bin/sh\n@p3 x\n"), [(2, "call")]);
        // `#` after `;`, `&` or `|` starts a comment.
        assert!(k("s.sh", ":;# @p3 is not needed\n").is_empty());
        assert!(k("s.sh", "true &# @p3 later\n").is_empty());
        assert!(k("s.sh", "echo \"${#arr[@]}\"\n").is_empty());
        // A <<PY inside quotes is text; a quoted delimiter after a real
        // operator is a heredoc, also in a Groovy """ line.
        assert!(k("s.sh", "printf '%s\\n' 'expected <<@P'\n").is_empty());
        assert!(k("s.sh", "echo \"a <<@PTHON b\"\n").is_empty());
        assert_eq!(k("s.sh", "cat <<\"@P\"\n"), [(1, "heredoc")]);
        assert_eq!(k("m.nf", "    \"\"\"cat > x <<'@P'\n"), [(1, "heredoc")]);
        assert!(k("s.sh", "cat <<<@P\n").is_empty());
        // A redirection right after the program; a path part stays out.
        assert_eq!(k("s.sh", "@p3</dev/null\n"), [(1, "call")]);
        assert_eq!(k("s.sh", "x=$(@p3>/dev/null)\n"), [(1, "call")]);
        assert!(k("s.sh", "ls lib/@p3.12/x\n").is_empty());
    }

    #[test]
    fn nextflow_comments_in_source_order() {
        // `// … /*` is a line comment: it opens no block, so the call
        // below is still read.
        let nf = "// Example: /*\nprocess P {\n    script:\n    \"\"\"\n    @p3 x\n    \"\"\"\n}\n";
        assert_eq!(k("m.nf", nf), [(5, "call")]);
        // A #! line inside a block comment is prose; after it, code again.
        let nf = "/*\n#!/usr/bin/env @p3\n*/\nprocess P {\n    script:\n    \"\"\"\n    #!/usr/bin/env @p3\n    \"\"\"\n}\n";
        assert_eq!(k("m.nf", nf), [(7, "shebang")]);
        // A block that closes and a call after it on the same line.
        assert_eq!(k("m.nf", "/* x */ @p3 y\n"), [(1, "call")]);
        assert!(k("m.nf", "x = 1 /* @p3 */\n").is_empty());
    }

    #[test]
    fn command_fields_and_workflow_labels() {
        let tsv = "#var\topts\tneeds\tcommand\texpect\n\
                   X_IMAGE\t-\t-\tprintf ok\trequires @p3 on the old runner\n\
                   Y_IMAGE\t-\t-\t@p3 -V\tPython 3\n";
        assert_eq!(k("tests/smoke/commands.tsv", tsv), [(3, "call")]);
        let yml = "name: @p3 checks\non: push\njobs:\n  a:\n    steps:\n      - name: No @p3 is needed\n        run: echo ok\n      - name: x\n        description: @p3 only in a label\n        run: @p3 -V\n";
        assert_eq!(k(".github/workflows/x.yml", yml), [(10, "call")]);
    }

    #[test]
    fn rust_raw_strings_and_comments_before_the_program() {
        let rs = "fn a() { let _ = Command::new(r\"@p3\"); }\n\
                  fn b() { let _ = Command::new(/* the program */ \"@p3\"); }\n\
                  fn c() { let _ = Command::new(r#\"/usr/bin/@p3\"#); }\n\
                  fn d() { let _ = Command::new(// next line\n\"@p3\"); }\n\
                          // Line 4 opens the call; its program string is on line 5, which
        // the string rule finds as well.
fn e() { let _ = Command::new(r\"og\"); }\n";
        assert_eq!(
            k("x.rs", rs),
            [
                (1, "rust"),
                (2, "rust"),
                (3, "rust"),
                (4, "rust"),
                (5, "rust")
            ]
        );
    }

    #[test]
    fn a_path_that_is_not_utf8_is_an_error() {
        assert_eq!(git_path(b"bin/x.py").as_deref(), Ok("bin/x.py"));
        let e = git_path(b"bin/x\xff.py").unwrap_err();
        assert!(e.contains("not UTF-8"), "{e}");
    }

    /// A git repository in a temporary folder with these files, and the
    /// gate's failures on it with the allow list `allow`.
    fn gate(name: &str, files: &[(&str, &str)], allow: &str) -> Vec<String> {
        let d = std::env::temp_dir().join(format!("og-@p-gate-{}-{name}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        let ok = Command::new("git")
            .args(["init", "-q"])
            .current_dir(&d)
            .status()
            .unwrap()
            .success();
        assert!(ok, "git init");
        let base = [
            (
                "versions.env",
                "TOOL_IMAGE=\"t/tool:1\"\nEMPTY_IMAGE=\"\"\n",
            ),
            ("a.rs", "fn main() {}\n"),
            ("b.sh", "echo ok\n"),
        ];
        for (rel, text) in base.iter().chain(files) {
            let p = d.join(rel);
            fs::create_dir_all(p.parent().unwrap()).unwrap();
            fs::write(p, t(text)).unwrap();
        }
        fs::write(d.join("allow.tsv"), allow).unwrap();
        #[cfg(unix)]
        if name == "symlink" {
            std::os::unix::fs::symlink("tool", d.join("bin/tool.py")).unwrap();
        }
        let r = check(&d, &d.join("allow.tsv")).unwrap();
        let _ = fs::remove_dir_all(&d);
        r.failures
    }

    #[test]
    fn allow_rows_count_uses_and_name_set_images() {
        let one = [("s.sh", "docker run x @p3 a\n")];
        let two = [("s.sh", "docker run x @p3 a\n@p3 -c 'print(1)'\n")];
        let row = "s.sh\t1\tTOOL_IMAGE: the tool's own Python\n";
        assert_eq!(gate("ok", &one, row), Vec::<String>::new());
        // A new use in an allowed file fails until the count changes.
        let f = gate("added", &two, row);
        assert!(
            f.len() == 1 && f[0].contains("s.sh holds 2 Python uses, the row says 1"),
            "{f:?}"
        );
        assert!(gate("counted", &two, "s.sh\t2\tTOOL_IMAGE: x\n").is_empty());
        // The old two-column row, a count that is not a number.
        let f = gate("cols", &one, "s.sh\tTOOL_IMAGE: x\n");
        assert!(
            f.iter().any(|f| f.contains("file<TAB>uses<TAB>reason")),
            "{f:?}"
        );
        let f = gate("nan", &one, "s.sh\tsome\tTOOL_IMAGE: x\n");
        assert!(f.iter().any(|f| f.contains("not 'some'")), "{f:?}");
        // An image variable set to "" is not set.
        let f = gate("empty", &one, "s.sh\t1\tEMPTY_IMAGE: x\n");
        assert!(
            f.iter()
                .any(|f| f.contains("names EMPTY_IMAGE, which versions.env does not set")),
            "{f:?}"
        );
        // A finding names the allow list the gate read.
        let f = gate("other", &one, "");
        assert!(
            f.len() == 1 && f[0].ends_with("allow.tsv)") && !f[0].contains(ALLOW),
            "{f:?}"
        );
    }

    /// Test text with `@p` for the word python, `@P` for PY, `@i` for pip
    /// and `@u` for uv, so no string in this file holds one of them and the
    /// gate scans this file like any other Rust file.
    fn t(s: &str) -> String {
        s.replace("@p", concat!("pyt", "hon"))
            .replace("@P", concat!("P", "Y"))
            .replace("@i", concat!("pi", "p"))
            .replace("@u", concat!("u", "v"))
    }

    fn k(rel: &str, text: &str) -> Vec<(usize, &'static str)> {
        kinds(rel, &t(text))
    }

    #[test]
    fn a_python_picked_into_a_variable() {
        assert_eq!(k("s.sh", "PY=@p3; $PY -c \"print(1)\"\n"), [(1, "call")]);
        assert_eq!(
            k("s.sh", "PY=${PYTHON:-@p3}; \"$PY\" -c 1\n"),
            [(1, "call")]
        );
        assert_eq!(k("s.sh", "PY=${PYTHON-@p3}\n"), [(1, "call")]);
        assert_eq!(k("s.sh", ": \"${PY:=@p3}\"\n"), [(1, "call")]);
        assert_eq!(k("s.sh", ": \"${PY=@p3.12}\"\n"), [(1, "call")]);
        assert_eq!(k("s.sh", "exec env PY=@p \"$@\"\n"), [(1, "call")]);
        // Not a word start: an image, a stamp, a dash inside a word.
        assert!(k("s.sh", "image=@p:3.15\n").is_empty());
        assert!(k("s.sh", "stamp=\"@p=${PYTHON_IMAGE}\"\n").is_empty());
        assert!(k("s.sh", "x=a-@p3\n--with-@p=3\n").is_empty());
    }

    #[test]
    fn brackets_escapes_and_braces_that_are_not_calls() {
        // An awk character class, a python shebang printed into a stub file,
        // an escaped dot in a path pattern, a variable and a placeholder.
        assert!(k("s.sh", "awk '/-O[uv]/ { next }'\n").is_empty());
        assert!(k("s.sh", "printf '#!/usr/bin/env @p\\n' > \"$d/x\"\n").is_empty());
        assert!(k("s.sh", "MNT=/usr/lib/@p3\\\\.12/x\n").is_empty());
        assert!(k("s.sh", "echo \"${@p}\"\n").is_empty());
        assert!(k("x.rs", "fn a() { println!(\"Python={@p}\"); }\n").is_empty());
        // An escaped quote ends the word.
        assert_eq!(
            k("x.rs", "const C: &str = \"sh -c \\\"@p3\\\"\";\n"),
            [(1, "rust")]
        );
    }

    #[test]
    fn python_run_from_a_rust_string() {
        let rs = t(
            r##"fn a() { Command::new("sh").args(["-c", "@p3 -c 'print(1)'"]); }
fn b() { Command::new("env").arg("@p3"); }
const PY: &str = "@p3";
fn c() { let _ = "#!/usr/bin/env @p3\nprint(1)"; }
fn d() { let _ = r#"cat > x <<'@P'"#; }
fn e() { let _ = "sh -c \"@i install x\""; }
fn f() { let _ = ["lib/@p3.12/x", "@p:3.15", "--@p-free", "PYTHON_IMAGE"]; }
"##,
        );
        assert_eq!(
            kinds("rust/crates/x/src/a.rs", &rs),
            [
                (1, "rust"),
                (2, "rust"),
                (3, "rust"),
                (4, "rust"),
                (5, "rust"),
                (6, "rust")
            ]
        );
    }

    #[test]
    fn package_runners_and_setup_python() {
        let sh = "@i install x\n@i3 install -r r.txt\n@ix run black\n@u run t.py\n\
                  @ux ruff\n/opt/v/bin/@i3.12 install y\necho @i.\n@iline\n@uvar=1\n";
        assert_eq!(
            k("s.sh", sh),
            [
                (1, "call"),
                (2, "call"),
                (3, "call"),
                (4, "call"),
                (5, "call"),
                (6, "call")
            ]
        );
        let yml = "jobs:\n  a:\n    steps:\n      - uses: actions/setup-@p@abc # v7\n      \
                   - uses: 'astral-sh/setup-@u@def'\n      - run: echo ok\n      \
                   # - uses: actions/setup-@p@abc\n";
        assert_eq!(
            k(".github/workflows/x.yml", yml),
            [(4, "call"), (5, "call")]
        );
        // Dependabot's config is data: it runs nothing.
        assert!(k(".github/dependabot.yml", "  - package-ecosystem: @i\n").is_empty());
    }

    #[test]
    fn makefiles_notebooks_brackets_and_javascript() {
        assert_eq!(
            k("justfile", "serve:\n    @p3 -m http.server\n"),
            [(2, "call")]
        );
        assert_eq!(
            k("assets/stub/db/megane.mk", "all:\n\t@p3 -c 1\n"),
            [(2, "call")]
        );
        assert_eq!(k("GNUmakefile", "x:\n\t@p3 -V\n"), [(2, "call")]);
        assert_eq!(k("docs/x.ipynb", "{}\n"), [(0, "file")]);
        assert_eq!(
            k(
                ".github/workflows/x.yml",
                "      - run: docker run --entrypoint [@p3,-c] img\n"
            ),
            [(1, "call")]
        );
        let js = "// @p3 in a comment\nexecSync('@p3 -c 1');\n/* @p3\n */ const a = '#x';\n";
        assert_eq!(k("tests/demo/gate-property.js", js), [(2, "call")]);
        assert_eq!(
            k("x.ts", "spawn(\"@p3\", [\"-c\", \"1\"]);\n"),
            [(1, "call")]
        );
    }

    #[test]
    fn a_mark_column_pins_what_each_use_runs() {
        let one = [("s.sh", t("docker run x @p3 /tool/pre a\n"))];
        let swapped = [("s.sh", t("@p3 -c 'print(1)'\n"))];
        let g = |name: &str, f: &[(&str, String)], allow: &str| {
            let refs: Vec<(&str, &str)> = f.iter().map(|(a, b)| (*a, b.as_str())).collect();
            gate(name, &refs, allow)
        };
        let row = "s.sh\t1\tTOOL_IMAGE: x\t/tool/pre\n";
        assert_eq!(g("mark-ok", &one, row), Vec::<String>::new());
        // A use swapped for a host one at the same count fails.
        let f = g("mark-swapped", &swapped, row);
        assert!(
            f.len() == 1 && f[0].contains("s.sh:1 holds none of the row's marks (/tool/pre)"),
            "{f:?}"
        );
        // Any one of several marks is enough.
        assert!(
            g(
                "mark-any",
                &one,
                "s.sh\t1\tTOOL_IMAGE: x\t/other|/tool/pre\n"
            )
            .is_empty()
        );
        // No mark column: the count alone, as before.
        assert!(g("mark-none", &swapped, "s.sh\t1\tTOOL_IMAGE: x\n").is_empty());
        // Five columns, or an empty mark, is a bad row.
        let f = g(
            "mark-five",
            &one,
            "s.sh\t1\tTOOL_IMAGE: x\t/tool/pre\textra\n",
        );
        assert!(f.iter().any(|f| f.contains("has 5 columns")), "{f:?}");
        let f = g("mark-empty", &one, "s.sh\t1\tTOOL_IMAGE: x\t/tool/pre||\n");
        assert!(f.iter().any(|f| f.contains("an empty mark")), "{f:?}");
    }

    #[cfg(unix)]
    #[test]
    fn a_symlink_named_like_a_python_file_is_one() {
        let f = gate("symlink", &[("bin/tool", "x\n")], "");
        assert!(
            f.len() == 1 && f[0].starts_with("bin/tool.py:0: file:"),
            "{f:?}"
        );
    }
}
