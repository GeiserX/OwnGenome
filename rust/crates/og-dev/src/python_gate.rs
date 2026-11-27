//! `og dev python-gate [--root DIR] [--allow FILE] [--list]`: no owned
//! Python is left in the repository.
//!
//! Every file `git ls-files` lists (tracked, plus new files that are not
//! ignored) is read, and these are the findings:
//!
//! - `file`: a `*.py`, `*.pyw`, `*.pyi`, `*.ipynb`, Cython `*.pyx` or
//!   `*.pxd`, a `*.py.in` template, or a Snakemake `Snakefile` or `*.smk`;
//! - `shebang`: a `#!` line whose program is python, python3, pythonX.Y,
//!   pypy, pypy3 or pypyX.Y, directly or through `env`: the first line of
//!   any file, or any line of a shell, Nextflow or workflow file (a script
//!   a heredoc or a Nextflow `script:` block writes);
//! - `call`: a Python program word outside a comment in a shell script (by
//!   its suffix, `*.bats` among them, by a shell or bats `#!` line whatever
//!   the suffix, a Makefile, a `*.mk`, a justfile, or a Dockerfile or
//!   Containerfile by any name: `Dockerfile.x`, `x.dockerfile`,
//!   `x.Dockerfile`), a Nextflow file or config, a workflow (GitHub, Gitea
//!   or Forgejo) or an action (`.github/actions/`, or any `action.yml`;
//!   not its `name:` and `description:` labels), JavaScript or TypeScript,
//!   or the `command` field of a TSV command table
//!   (`tests/smoke/commands.tsv`): a `python3 -c`, a `python3 - <<EOF`, a
//!   `python3 tool` run in an image, a `\python3` past an alias, a
//!   `PY=python3` or `${PYTHON:-python3}` that picks the interpreter. The
//!   program words are python, python3, pythonX.Y, pypy, pypy3, pypyX.Y and
//!   the package runners pip, pip3, pipx, uv and uvx. A run of the python
//!   image (`PYTHON_IMAGE`, `python:3.15`, `library/python`) by `docker
//!   run`, `podman run` or `run_in` with no program after the image is a
//!   call too: the image's own python3 runs what comes on stdin
//!   (`docker run -i "$PYTHON_IMAGE" <<'END'`, `... < tool.txt`), wherever
//!   the redirection stands; a command continued with `\` is read whole,
//!   in a Rust string too. A workflow step that `uses:` a
//!   `setup-python` or `setup-uv` action (in a block or a flow mapping) is
//!   a call too. A shell
//!   quote may span lines (a `#` inside it is text) but never opens inside
//!   a heredoc body, where an apostrophe is data. A `//` comment may follow
//!   a token with no blank (`"ok";// note`). A Nextflow `"""` or
//!   `'''` body is shell (`rm cache/*` is a glob, not a comment), and a
//!   Make or just recipe drops its `@`, `-` and `+` prefixes first;
//! - `heredoc`: a heredoc whose delimiter says it holds Python (`<<'PY'`,
//!   `<<PYTHON`, `<<EOF_PYTHON`), not a `<<PY` inside quotes;
//! - `rust`: in Rust code outside comments, a `Command::new("python3")`
//!   (any python program, any path, an ordinary or raw string), the
//!   variable that named a Python oracle, or a string literal that the
//!   shell rules above find a call, heredoc or shebang in
//!   (`.arg("python3")`, `"sh -c 'python3 x'"`, `const PY = "python3"`).
//!   An ordinary literal is read with its escapes decoded, as the program
//!   gets it: `"python\x33"` and `"python\u{33}"` are python3.
//!
//! A comment, a doc sentence and a Markdown, text, JSON or image file are
//! prose: they never count. `python:3.15` (an image), `python=...` (a
//! stamp), `PYTHON_IMAGE`, `setup-python` outside a `uses:` line and a
//! `site-packages/python3.12/` path are not calls.
//!
//! Known limits: a Python program name built at run time
//! (`format!("python{}", 3)`, `["py", "thon"].concat()`) or by shell
//! quoting (`py""thon3`) is not seen; nor is a program word with a suffix
//! of its own (`python3-dbg`), a shell file with no suffix and no `#!`
//! line that another script sources, a python image named through a
//! variable of its own (`img=$PYTHON_IMAGE; docker run -i "$img" < x`) or
//! through separate Rust strings (`.args(["run", "-i", &image])`), and
//! YAML other than workflows and actions is data (a
//! `.pre-commit-config.yaml` `entry: python3`, a compose `command:`).
//! Review catches those.
//!
//! `tests/python-allow.tsv` names, one row per file, the files whose
//! findings are not owned Python: a third-party tool run inside its own
//! pinned image (PharmCAT's preprocessor, eKLIPse, Cyrius and the docs
//! builder Zensical in `PYTHON_IMAGE`), or a line that runs no Python (a
//! check of a log, a stub, an image's name). Each row is `file<TAB>uses<TAB>reason[<TAB>marks]`,
//! where `uses` is the number of findings (one per use: a line with two
//! Python commands holds two) and `marks`, if given, is a `|`-separated
//! list of texts one of which the code of every such line, without its
//! comments, must hold (`/pharmcat/pharmcat_vcf_preprocessor`), so a use swapped for
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
const PYPY: &str = concat!("py", "py");
const PIP: &str = concat!("pi", "p");
/// The python image's variable in `versions.env`.
const PYTHON_IMAGE: &str = concat!("PYT", "HON_IMAGE");
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
        "sh" | "bash" | "ksh" | "zsh" | "mk" | "bats" => return Lang::Shell,
        // `og.dockerfile`, `og.Dockerfile`.
        "dockerfile" | "Dockerfile" => return Lang::Shell,
        "nf" | "config" | "groovy" => return Lang::Nextflow,
        "rs" => return Lang::Rust,
        "js" | "mjs" | "cjs" | "jsx" | "ts" | "mts" | "cts" | "tsx" => return Lang::Script,
        // Workflows and actions run code; the rest of .github (Dependabot,
        // issue forms) is data. A local action may live anywhere
        // (`uses: ./tools/act`), and Gitea and Forgejo read workflows of
        // their own.
        "yml" | "yaml"
            if [
                ".github/workflows/",
                ".github/actions/",
                ".gitea/workflows/",
                ".forgejo/workflows/",
            ]
            .iter()
            .any(|d| rel.starts_with(d))
                || matches!(name, "action.yml" | "action.yaml") =>
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
        .is_some_and(|p| matches!(p.as_str(), "sh" | "bash" | "dash" | "ksh" | "zsh" | "bats"));
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

/// python, python3, python3.12, and PyPy's pypy, pypy3, pypy3.10: not
/// python3.12/ (a path part) and not python. (a sentence).
fn is_python_program(word: &str) -> bool {
    [PYTHON, PYPY]
        .iter()
        .any(|p| word.strip_prefix(p).is_some_and(version_suffix))
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

/// Shell quotes a line leaves open, carried to the next line.
#[derive(Clone, Copy, Default)]
struct Quotes {
    single: bool,
    double: bool,
}

/// A shell line without its comment. `q` carries the quotes a line leaves
/// open, so a `#` inside a quoted string that spans lines is text, and the
/// command after the quote closes is still read.
fn strip_line_comment_in<'a>(line: &'a str, q: &mut Quotes) -> &'a str {
    let b = line.as_bytes();
    let mut i = 0;
    while i < b.len() {
        match b[i] {
            b'\\' if !q.single => {
                i += 2;
                continue;
            }
            b'\'' if !q.double => q.single = !q.single,
            b'"' if !q.single => q.double = !q.double,
            b'#' if !q.single && !q.double && comment_start(b, i) => return &line[..i],
            _ => {}
        }
        i += 1;
    }
    line
}

/// A shell line without its comment, outside quotes opened on this line.
fn strip_line_comment(line: &str) -> &str {
    strip_line_comment_in(line, &mut Quotes::default())
}

/// What a Nextflow or JavaScript file leaves open at the end of a line.
#[derive(Clone, Copy, Default)]
struct Open {
    /// A `/* */` block comment.
    block: bool,
    /// A Nextflow `"""` or `'''` string: the shell body of a `script:`
    /// block, where `/*` and `//` are a glob and a path, not comments.
    body: bool,
    /// A JavaScript backtick template literal, which may span lines.
    tick: bool,
}

/// A Nextflow line without its comments, read in order: `/* */` (`st`
/// carries a block that spans lines), `//` anywhere outside quotes (a JavaScript backtick quotes too),
/// and, with `nf`, a `#` shell comment. With `nf`, a `"""` or `'''` opens
/// and closes a shell body, where only `#` comments count (`rm cache/*` is
/// a glob). Only a marker outside quotes counts, and a `/*` after `//` is
/// part of the line comment. JavaScript (`nf` false) has neither.
fn strip_nf_comments(line: &str, st: &mut Open, nf: bool) -> String {
    let b = line.as_bytes();
    let mut out = String::new();
    // JavaScript also quotes with backticks: a template literal
    // (`http://x`), carried in `st` because it may span lines.
    let (mut single, mut double) = (false, false);
    let (mut i, mut seg) = (0, 0);
    while i < b.len() {
        if st.block {
            if b[i] == b'*' && b.get(i + 1) == Some(&b'/') {
                st.block = false;
                seg = i + 2;
                i += 2;
            } else {
                i += 1;
            }
            continue;
        }
        let quoted = single || double || st.tick;
        if nf && !quoted && (b[i..].starts_with(b"\"\"\"") || b[i..].starts_with(b"'''")) {
            st.body = !st.body;
            i += 3;
            continue;
        }
        match b[i] {
            b'\\' if !single => {
                i += 2;
                continue;
            }
            b'\'' if !double && !st.tick => single = !single,
            b'"' if !single && !st.tick => double = !double,
            b'`' if !nf && !single && !double => st.tick = !st.tick,
            b'/' if !quoted && !st.body && b.get(i + 1) == Some(&b'*') => {
                out.push_str(&line[seg..i]);
                st.block = true;
                i += 2;
                continue;
            }
            // `//` right after a token is a comment too (`"ok";// note`).
            b'/' if !quoted && !st.body && b.get(i + 1) == Some(&b'/') => {
                out.push_str(&line[seg..i]);
                return out;
            }
            b'#' if nf && !quoted && comment_start(b, i) => {
                out.push_str(&line[seg..i]);
                return out;
            }
            _ => {}
        }
        i += 1;
    }
    if !st.block {
        out.push_str(&line[seg.min(line.len())..]);
    }
    out
}

/// Can a program word start at `s`? At the line start, after a blank, a
/// quote, a bracket, a shell operator, a path `/`, a list `,`, an
/// assignment `=` (`PY=python3`) or a `\` (`\python3` skips an alias and
/// runs python3), and after the `-` or `+` of `${NAME-…}`, `${NAME:-…}`,
/// `${NAME+…}` and `${NAME:+…}`.
fn word_start(b: &[u8], s: usize) -> bool {
    if s == 0 {
        return true;
    }
    let p = b[s - 1];
    if p.is_ascii_whitespace() || b"\"'`(;|&{/[,=\\".contains(&p) {
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

/// How many Python program words a code line holds (python, python3,
/// pythonX.Y, pip, pip3, pipx, uv, uvx): each is a use.
fn python_words(code: &str) -> usize {
    let b = code.as_bytes();
    let mut uses = 0;
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
            uses += 1;
        }
    }
    uses
}

/// A workflow step that installs a Python: its `uses:` key, in a block
/// mapping (`- uses: x`) or a flow mapping (`- {name: n, uses: x}`), names
/// a `setup-python` or `setup-uv` action, at any owner and ref.
fn setup_python(code: &str) -> bool {
    let mut from = 0;
    while let Some(k) = code[from..].find("uses:") {
        let s = from + k;
        from = s + 5;
        // The key starts the line (after `- `) or follows `{` or `,`: not a
        // `- uses:` inside a printed string.
        let before = code[..s].trim();
        if !(before.is_empty() || before == "-" || before.ends_with('{') || before.ends_with(',')) {
            continue;
        }
        let value = code[from..].trim_start();
        let value = value.trim_start_matches(['"', '\'']);
        let action = value
            .split(['@', ',', '}', '"', '\'', ' ', '\t'])
            .next()
            .unwrap_or("");
        let name = action.rsplit('/').next().unwrap_or("");
        if name == ["setup-", PYTHON].concat() || name == ["setup-", RUNNERS[1]].concat() {
            return true;
        }
    }
    false
}

/// The heredocs a code line opens, in order: each delimiter with its quotes
/// and backslashes removed (`<<'EOF'` ends at `EOF`), and whether `<<-`
/// lets the closing line start with tabs. Only a `<<` outside quotes is an
/// operator (`printf 'x <<PY'` prints text), and `<<<` is a here-string. A
/// Groovy `"""` or `'''` opens a string whose body is shell, so it does not
/// count as a quote.
fn heredoc_delims(code: &str) -> Vec<(String, bool)> {
    let b = code.as_bytes();
    let mut out = Vec::new();
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
                i += 2;
                let dash = b.get(i) == Some(&b'-');
                while i < b.len() && (b[i] == b'-' || b[i] == b'~' || b[i] == b' ' || b[i] == b'\t')
                {
                    i += 1;
                }
                let mut word = String::new();
                let mut q = 0u8;
                while i < b.len() {
                    let c = b[i];
                    if q != 0 {
                        if c == q {
                            q = 0;
                        } else {
                            word.push(c as char);
                        }
                    } else if c == b'\'' || c == b'"' {
                        q = c;
                    } else if c == b'\\' {
                        if let Some(&n) = b.get(i + 1) {
                            word.push(n as char);
                        }
                        i += 1;
                    } else if c.is_ascii_whitespace() || b";|&<>()".contains(&c) {
                        break;
                    } else {
                        word.push(c as char);
                    }
                    i += 1;
                }
                if !word.is_empty() {
                    out.push((word, dash));
                }
                continue;
            }
            _ => {}
        }
        i += 1;
    }
    out
}

/// Does this code line open a heredoc whose delimiter says Python
/// (`<<'PY'`, `<<-EOF_PYTHON`, `<<EOPY`)?
fn python_heredoc(code: &str) -> bool {
    heredoc_delims(code).iter().any(|(d, _)| {
        let word: String = d
            .chars()
            .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
            .collect();
        let up = word.to_ascii_uppercase();
        up.starts_with("PY") || up.contains("PYTHON") || up == "EOPY"
    })
}

/// The words and operators of a shell command, roughly: quotes are dropped
/// (`"$PYTHON_IMAGE"` is the word `$PYTHON_IMAGE`), a `${…}` stays one word
/// with its blanks, a `\` keeps the character after it, and a run of
/// `; | & ( ) < >` is one operator.
fn shell_tokens(code: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut it = code.chars().peekable();
    while let Some(c) = it.next() {
        match c {
            '\'' | '"' | '`' => {}
            '\\' => cur.extend(it.next()),
            '$' if it.peek() == Some(&'{') => {
                cur.push(c);
                let mut depth = 0usize;
                for d in it.by_ref() {
                    cur.push(d);
                    if d == '{' {
                        depth += 1;
                    } else if d == '}' {
                        depth -= 1;
                        if depth == 0 {
                            break;
                        }
                    }
                }
            }
            ';' | '|' | '&' | '(' | ')' | '<' | '>' => {
                if !cur.is_empty() {
                    out.push(std::mem::take(&mut cur));
                }
                let mut op = String::from(c);
                while let Some(n) = it.next_if(|n| ";|&()<>".contains(*n)) {
                    op.push(n);
                }
                out.push(op);
            }
            c if c.is_whitespace() => {
                if !cur.is_empty() {
                    out.push(std::mem::take(&mut cur));
                }
            }
            _ => cur.push(c),
        }
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}

/// Is this word the python image: `$PYTHON_IMAGE`, `${PYTHON_IMAGE}` (with
/// any `:-`, `:?` tail), or a `python` or `…/library/python` reference with
/// any tag or digest (`python:3.15`, `docker.io/library/python@sha256:…`)?
fn python_image_ref(word: &str) -> bool {
    if let Some(r) = word.strip_prefix("${") {
        return r
            .strip_prefix(PYTHON_IMAGE)
            .is_some_and(|rest| rest.starts_with(['}', ':', '-', '+', '=', '?']));
    }
    if let Some(r) = word.strip_prefix('$') {
        return r == PYTHON_IMAGE;
    }
    let name = word.split('@').next().unwrap_or(word);
    let slash = name.rfind('/').map_or(0, |k| k + 1);
    let repo = &name[..name[slash..].find(':').map_or(name.len(), |k| slash + k)];
    let library = ["library/", PYTHON].concat();
    repo == PYTHON || repo == library || repo.ends_with(&format!("/{library}"))
}

/// Does a command run the python image with no program after the image?
/// Then the image's own python3, its default command, runs what comes on
/// stdin (`docker run -i "$PYTHON_IMAGE" <<'END'`, `… < tool.txt`). The run
/// is `docker run`, `podman run` or `nerdctl run` (also `container run`),
/// or the scripts' `run_in`. A redirection may stand anywhere in the
/// command (`docker run -i < tool.txt "$PYTHON_IMAGE"`): it and its file
/// are skipped, and a word after it is still the program (`... < x sh`).
/// A bare `python` image is a program word already, so it is not counted
/// twice.
fn python_image_default(code: &str) -> bool {
    let toks = shell_tokens(code);
    // `;`, `|`, `&`, `&&`, `||`, `(` and `)` end the command; `<`, `>`,
    // `<<`, `>&`, `&>` and the like are redirections, each with one word
    // after it.
    let control = |t: &str| {
        t.starts_with([';', '|', '(', ')']) || (t.starts_with('&') && !t.starts_with("&>"))
    };
    let redirect = |t: &str| t.starts_with(['<', '>']) || t.starts_with("&>");
    // A file descriptor before a redirection (`2>/dev/null`) is no word of
    // the command.
    let fd = |p: usize| {
        toks[p].bytes().all(|c| c.is_ascii_digit()) && toks.get(p + 1).is_some_and(|n| redirect(n))
    };
    let mut i = 0;
    while i < toks.len() {
        let base = toks[i].rsplit('/').next().unwrap_or("");
        let mut from = None;
        if matches!(base, "docker" | "podman" | "nerdctl") {
            let mut j = i + 1;
            if toks.get(j).is_some_and(|t| t == "container") {
                j += 1;
            }
            if toks.get(j).is_some_and(|t| t == "run") {
                from = Some(j + 1);
            }
        } else if toks[i] == "run_in" {
            from = Some(i + 1);
        }
        if let Some(mut j) = from {
            while j < toks.len() && !control(&toks[j]) {
                if redirect(&toks[j]) {
                    j += 2;
                    continue;
                }
                if python_image_ref(&toks[j]) {
                    if is_python_program(&toks[j]) {
                        break;
                    }
                    // The program: the first word after the image that is no
                    // redirection, no file of one and no descriptor.
                    let mut p = j + 1;
                    loop {
                        match toks.get(p) {
                            None => return true,
                            Some(t) if control(t) => return true,
                            Some(t) if redirect(t) => p += 2,
                            Some(_) if fd(p) => p += 1,
                            Some(_) => break,
                        }
                    }
                    break;
                }
                j += 1;
            }
            i = j;
        }
        i += 1;
    }
    false
}

/// A shell command continued with `\` from the lines before, so a run of
/// the python image is read whole.
#[derive(Default)]
struct Continued(String);

impl Continued {
    /// Reads the next line's code. On the command's last line (one that
    /// does not end in `\`): does the whole command run the python image's
    /// own python3 (`python_image_default`)?
    fn ends_python_image_run(&mut self, code: &str) -> bool {
        let whole = if self.0.is_empty() {
            code.to_owned()
        } else {
            format!("{} {code}", self.0)
        };
        match whole.trim_end().strip_suffix('\\') {
            Some(head) => {
                self.0 = head.to_owned();
                false
            }
            None => {
                self.0.clear();
                python_image_default(&whole)
            }
        }
    }
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

/// A finding on line `line` (1-based, 0 for the file) of `rel`: `text` is
/// the whole line, `code` the part a row's marks are looked for in (the
/// line without its comments).
fn finding(rel: &str, line: usize, kind: &'static str, text: &str, code: &str) -> Finding {
    Finding {
        file: rel.to_owned(),
        line,
        kind,
        text: excerpt(text),
        code: code.trim().to_owned(),
    }
}

/// The `file` finding of a path whose name is a Python file's: Python,
/// a notebook, Cython (`.pyx`, `.pxd`), a template of a Python file
/// (`.py.in`) or Snakemake, which is Python (`Snakefile`, `.smk`).
fn python_name(rel: &str) -> Option<Finding> {
    let name = rel.rsplit('/').next().unwrap_or(rel);
    ([
        ".py", ".pyw", ".pyi", ".ipynb", ".pyx", ".pxd", ".py.in", ".smk",
    ]
    .iter()
    .any(|e| name.ends_with(e))
        || matches!(name, "Snakefile" | "snakefile"))
    .then(|| finding(rel, 0, "file", "a Python file", ""))
}

/// A Makefile or a justfile: an indented recipe line may start with `@`,
/// `-` or `+`, which make and just take off before the shell runs it.
fn recipe_file(rel: &str) -> bool {
    let name = rel.rsplit('/').next().unwrap_or(rel);
    name.ends_with(".mk")
        || matches!(
            name,
            "Makefile" | "makefile" | "GNUmakefile" | "justfile" | "Justfile" | ".justfile"
        )
}

/// The findings in one file's text: one per use, so a line with two
/// Python commands has two.
pub fn scan_text(rel: &str, text: &str) -> Vec<Finding> {
    let mut out: Vec<Finding> = python_name(rel).into_iter().collect();
    let first = text.lines().next().unwrap_or("");
    if shebang_is_python(first) {
        out.push(finding(rel, 1, "shebang", first, first));
    }
    match lang(rel, text) {
        Lang::Other => {}
        Lang::Rust => scan_rust(rel, text, &mut out),
        l @ (Lang::Shell | Lang::Workflow | Lang::Nextflow | Lang::CommandTable | Lang::Script) => {
            // `//` and `/* */` comments; Nextflow's script blocks also have
            // shell `#` comments, JavaScript has none.
            let slashes = matches!(l, Lang::Nextflow | Lang::Script);
            let recipes = recipe_file(rel);
            // A command table: only its `command` field is shell.
            let column = if l == Lang::CommandTable {
                command_column(text)
            } else {
                None
            };
            let mut open = Open::default();
            // A shell quote can span lines; the other kinds start each line
            // outside quotes.
            let mut quotes = Quotes::default();
            // Heredoc bodies still to read, in order. A body is data: its quotes
            // stay on their own line and never open a quote in the shell
            // around it (`It's ready.` inside `cat <<EOF`).
            let mut bodies: std::collections::VecDeque<(String, bool)> =
                std::collections::VecDeque::new();
            let mut command = Continued::default();
            for (i, line) in text.lines().enumerate() {
                let n = i + 1;
                if l == Lang::CommandTable && line.starts_with('#') {
                    continue;
                }
                if l == Lang::Workflow && yaml_label(line) {
                    continue;
                }
                let mut field = match column {
                    Some(c) => line.split('\t').nth(c).unwrap_or(""),
                    None => line,
                };
                if recipes && field.starts_with([' ', '\t']) {
                    field = field.trim_start().trim_start_matches(['@', '-', '+']);
                }
                // A `#!` line inside a Nextflow block comment is prose.
                if n > 1
                    && !open.block
                    && !quotes.single
                    && !quotes.double
                    && shebang_is_python(field)
                {
                    out.push(finding(rel, n, "shebang", line, field));
                    continue;
                }
                let stripped;
                let code = if slashes {
                    stripped = strip_nf_comments(field, &mut open, l == Lang::Nextflow);
                    stripped.as_str()
                } else if l == Lang::Shell && !bodies.is_empty() {
                    strip_line_comment(field)
                } else if l == Lang::Shell {
                    strip_line_comment_in(field, &mut quotes)
                } else {
                    strip_line_comment(field)
                };
                if l == Lang::Shell {
                    match bodies.front() {
                        Some((end, dash)) => {
                            let close = if *dash {
                                field.trim_start_matches('\t')
                            } else {
                                field
                            };
                            if close == end {
                                bodies.pop_front();
                                continue;
                            }
                        }
                        None => bodies.extend(heredoc_delims(code)),
                    }
                }
                // The marks of a command-table row may name its image column.
                let marked = match column {
                    Some(_) => format!("{}\t{code}", line.split('\t').next().unwrap_or("")),
                    None => code.to_owned(),
                };
                // A heredoc line holds at least the heredoc; a second command
                // beside it is a second use.
                let default_python = command.ends_python_image_run(code);
                // A Python heredoc fed to the image is that one use.
                let uses = if python_heredoc(code) {
                    (python_words(code).max(1), "heredoc")
                } else {
                    let words = python_words(code);
                    let setup = usize::from(l == Lang::Workflow && setup_python(code));
                    (words + setup + usize::from(default_python), "call")
                };
                for _ in 0..uses.0 {
                    out.push(finding(rel, n, uses.1, line, &marked));
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

/// The lines of a Rust string literal's body as the program gets them,
/// each with the source position it starts at (`start` is the body's). In
/// an ordinary literal the escapes are decoded (`\x33` and `\u{33}` are
/// `3`, `\"` is `"`, a `\` at the end of a line joins the next one and the
/// blanks it starts with); a raw literal keeps its text. A real line break
/// or a `\n` ends a line in both, as it does for the shell that may run it.
fn literal_lines(body: &str, start: usize, raw: bool) -> Vec<(String, usize)> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut at = start;
    let mut it = body.char_indices().peekable();
    while let Some((j, c)) = it.next() {
        if c == '\n' {
            out.push((std::mem::take(&mut cur), at));
            at = start + j + 1;
            continue;
        }
        if c != '\\' {
            cur.push(c);
            continue;
        }
        if raw {
            if it.next_if(|&(_, n)| n == 'n').is_some() {
                out.push((std::mem::take(&mut cur), at));
                at = start + j + 2;
            } else {
                cur.push('\\');
            }
            continue;
        }
        let Some((k, e)) = it.next() else { break };
        match e {
            'n' => {
                out.push((std::mem::take(&mut cur), at));
                at = start + k + 1;
            }
            't' => cur.push('\t'),
            'r' => cur.push('\r'),
            '0' => cur.push('\0'),
            'x' => {
                let mut hex = String::new();
                while hex.len() < 2 {
                    match it.next_if(|&(_, h)| h.is_ascii_hexdigit()) {
                        Some((_, h)) => hex.push(h),
                        None => break,
                    }
                }
                cur.extend(u8::from_str_radix(&hex, 16).ok().map(char::from));
            }
            'u' if it.next_if(|&(_, h)| h == '{').is_some() => {
                let mut hex = String::new();
                for (_, h) in it.by_ref() {
                    if h == '}' {
                        break;
                    }
                    if h != '_' {
                        hex.push(h);
                    }
                }
                cur.extend(u32::from_str_radix(&hex, 16).ok().and_then(char::from_u32));
            }
            '\n' => while it.next_if(|&(_, w)| w.is_whitespace()).is_some() {},
            other => cur.push(other),
        }
    }
    out.push((cur, at));
    out
}

fn scan_rust(rel: &str, src: &str, out: &mut Vec<Finding>) {
    let part = rust_parts(src);
    let line_of = |pos: usize| src[..pos].bytes().filter(|&c| c == b'\n').count() + 1;
    let line_text = |pos: usize| {
        let s = src[..pos].rfind('\n').map_or(0, |k| k + 1);
        let e = src[pos..].find('\n').map_or(src.len(), |k| pos + k);
        &src[s..e]
    };
    // The line without its comments, where a row's marks are looked for.
    let line_code = |pos: usize| {
        let s = src[..pos].rfind('\n').map_or(0, |k| k + 1);
        let e = src[pos..].find('\n').map_or(src.len(), |k| pos + k);
        let kept: Vec<u8> = (s..e)
            .filter(|&p| part[p] != Part::Comment)
            .map(|p| src.as_bytes()[p])
            .collect();
        String::from_utf8_lossy(&kept).into_owned()
    };
    // Uses per line (at a position on it): each python word of a string, a
    // heredoc or #! line, the oracle variable. A Command::new of a python
    // program is one use at least on its line; its program string is that
    // same use, not a second one.
    let mut uses: BTreeMap<usize, (usize, usize)> = BTreeMap::new();
    let mut calls: BTreeMap<usize, usize> = BTreeMap::new();
    let mut add = |pos: usize, k: usize| {
        uses.entry(line_of(pos)).or_insert((0, pos)).0 += k;
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
        // `"python\x33"` is python3; a raw string keeps its backslashes.
        let prog: Vec<String> = literal_lines(&src[body..end], body, open > q)
            .into_iter()
            .map(|(s, _)| s)
            .collect();
        let prog = prog.join("\n");
        let base = prog.rsplit('/').next().unwrap_or(&prog);
        if is_python_program(base) {
            calls.entry(line_of(s)).or_insert(s);
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
            add(s, 1);
        }
    }
    // Every string literal, read with the shell rules on each of its lines
    // (a real line break or a `\n` escape) as the program gets them (an
    // ordinary literal's escapes decoded): a program word (`.arg("python3")`,
    // `"sh -c 'python3 x'"`), a Python heredoc, a python `#!` line, or a run
    // of the python image with no program after it.
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
        // The body of `r"…"` or `r#"…"#` follows `r"` or `#"`; that of an
        // ordinary or byte string follows `"` or `b"`.
        let raw = start >= 2 && b[start - 1] == b'"' && matches!(b[start - 2], b'r' | b'#');
        // Each literal starts a command of its own.
        let mut command = Continued::default();
        for (piece, at) in literal_lines(&src[start..i], start, raw) {
            let code = strip_line_comment(&piece);
            let default_python = command.ends_python_image_run(code);
            let k = if shebang_is_python(&piece) {
                1
            } else if python_heredoc(code) {
                python_words(code).max(1)
            } else {
                python_words(code) + usize::from(default_python)
            };
            if k > 0 {
                add(at, k);
            }
        }
    }
    for (n, pos) in calls {
        let e = uses.entry(n).or_insert((0, pos));
        e.0 = e.0.max(1);
    }
    for (n, (k, pos)) in uses {
        for _ in 0..k {
            out.push(finding(rel, n, "rust", line_text(pos), &line_code(pos)));
        }
    }
}

/// One row of the allow list.
#[derive(Debug)]
struct Row {
    file: String,
    /// How many Python uses (findings) the file holds.
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

    /// Test text with `@p` for the word python, `@P` for PY, `@i` for pip,
    /// `@u` for uv, `@y` for pypy and `@V` for the python image's variable,
    /// so no string in this file holds one of them and the gate scans this
    /// file like any other Rust file.
    fn t(s: &str) -> String {
        s.replace("@p", concat!("pyt", "hon"))
            .replace("@P", concat!("P", "Y"))
            .replace("@i", concat!("pi", "p"))
            .replace("@u", concat!("u", "v"))
            .replace("@y", concat!("py", "py"))
            .replace("@V", concat!("PYT", "HON_IMAGE"))
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
        // A second Python on an allowed line is a second use.
        let more = [("s.sh", t("docker run x @p3 /tool/pre a; @p3 -V\n"))];
        let f = g("mark-more", &more, row);
        assert!(
            f.len() == 1 && f[0].contains("s.sh holds 2 Python uses, the row says 1"),
            "{f:?}"
        );
        // A mark in a comment does not count: it must be in the code.
        let commented = [("s.sh", t("@p3 -V # /tool/pre\n"))];
        let f = g("mark-comment", &commented, row);
        assert!(
            f.len() == 1 && f[0].contains("s.sh:1 holds none of the row's marks"),
            "{f:?}"
        );
    }

    #[test]
    fn each_program_word_on_a_line_is_a_use() {
        assert_eq!(
            k("s.sh", "@p3 a; @p3 b | @i install c\n"),
            [(1, "call"), (1, "call"), (1, "call")]
        );
        let rs = t("fn a() { let _ = [\"@p3 a\", \"sh -c '@p3 b; @ux c'\"]; }\n");
        assert_eq!(kinds("x.rs", &rs), [(1, "rust"), (1, "rust"), (1, "rust")]);
        // Command::new("python3") is one use, not one for the call and one
        // for its string.
        assert_eq!(
            k("x.rs", "fn a() { Command::new(\"@p3\"); }\n"),
            [(1, "rust")]
        );
        // A Python heredoc line counts its own program words too: a second
        // command beside the heredoc is a second use.
        assert_eq!(
            k("s.sh", "docker run tool @p3 /tool/pre a; @p3 -V <<@P\n@P\n"),
            [(1, "heredoc"), (1, "heredoc")]
        );
        assert_eq!(k("s.sh", "@p3 - <<@P\n@P\n"), [(1, "heredoc")]);
        assert_eq!(k("s.sh", "cat <<@P\n@P\n"), [(1, "heredoc")]);
        let rs = t("fn a() { let _ = \"@p3 x; @p3 - <<@P\"; }\n");
        assert_eq!(kinds("x.rs", &rs), [(1, "rust"), (1, "rust")]);
    }

    #[test]
    fn a_quote_that_spans_lines_hides_no_command() {
        // The # on line 2 is inside the quote opened on line 1; python3
        // runs after the quote closes.
        let sh = "printf '%s\\n' 'first\n# second'; @p3 -V\n";
        assert_eq!(k("s.sh", sh), [(2, "call")]);
        // A real comment after a closed quote is still a comment.
        assert!(k("s.sh", "echo 'a\nb' # @p3 later\n").is_empty());
    }

    #[test]
    fn a_heredoc_body_opens_no_shell_quote() {
        // The apostrophe in the body is data; python3 after the body runs.
        let sh = "cat <<'EOF'\nIt's ready.\nEOF\nprintf '%s' 'a # b'; @p3 -V\n";
        assert_eq!(k("s.sh", sh), [(4, "call")]);
        // <<- lets the closing line start with tabs; a substitution in the
        // body still runs.
        let sh = "cat <<-EOF\n\tdon't\n\tv=$(@p3 -V)\n\tEOF\necho 'a # b'; @p3 -V\n";
        assert_eq!(k("s.sh", sh), [(3, "call"), (5, "call")]);
        // Two heredocs on one line are read in order.
        let sh = "cat <<A <<\"B\"\nit's\nA\nit's\nB\necho 'x # y'; @p3 -V\n";
        assert_eq!(k("s.sh", sh), [(6, "call")]);
    }

    #[test]
    fn a_javascript_line_comment_right_after_a_token() {
        let js = "const label = \"ok\";// Example: /*\nspawn(\"@p3\", [\"-V\"]);\n";
        assert_eq!(k("a.js", js), [(2, "call")]);
        // After a colon too: `key:// note` is a comment in an object literal.
        let js = "const c = { key:// Example: /*\n  spawn(\"@p3\", [\"-V\"]) };\n";
        assert_eq!(k("a.js", js), [(2, "call")]);
        // A // inside a template literal is text.
        assert_eq!(
            k("a.js", "const u = `http://x`; spawn(\"@p3\");\n"),
            [(1, "call")]
        );
        // A template literal may span lines.
        let js = "const t = `first\nhttp://x`; spawn(\"@p3\", [\"-V\"]);\n";
        assert_eq!(k("a.js", js), [(2, "call")]);
    }

    #[test]
    fn a_glob_in_a_nextflow_shell_body_opens_no_comment() {
        let nf =
            "process P {\n    script:\n    \"\"\"\n    rm -f cache/*\n    @p3 -V\n    \"\"\"\n}\n";
        assert_eq!(k("m.nf", nf), [(5, "call")]);
        // Outside the body, /* still opens a Groovy block comment.
        let nf = "/* cache/\n@p3 -V\n*/\nprocess P {\n    script:\n    '''\n    ls a//b\n    @p3 x\n    '''\n}\n";
        assert_eq!(k("m.nf", nf), [(8, "call")]);
    }

    #[test]
    fn setup_python_in_a_flow_mapping() {
        let yml = "jobs:\n  a:\n    steps:\n      - {uses: actions/setup-@p@v7}\n      \
                   - {name: x, uses: 'astral-sh/setup-@u@v1'}\n      - {name: setup-@p, run: echo}\n";
        assert_eq!(
            k(".github/workflows/x.yml", yml),
            [(4, "call"), (5, "call")]
        );
        // A step printed into a file by a control runs nothing.
        let yml = "        run: |\n          printf '      - uses: actions/setup-@p@v7\\n' >> f.yml\n          check '- uses: actions/setup-@p@v7'\n";
        assert!(k(".github/workflows/x.yml", yml).is_empty());
    }

    #[test]
    fn make_and_just_recipe_prefixes() {
        // Make drops @, - and + before it hands the line to the shell; so
        // does just with @ and -.
        assert_eq!(
            k("Makefile", "x:\n\t@@p3 -V\n\t-@p3 a\n\t+@@i install x\n"),
            [(2, "call"), (3, "call"), (4, "call")]
        );
        assert_eq!(
            k("justfile", "serve:\n    @@p3 -m http.server\n"),
            [(2, "call")]
        );
        // In a shell script, -python3 stays an option name.
        assert!(k("s.sh", "\t-@p3 a\n").is_empty());
    }

    #[test]
    fn python_fed_to_the_python_image_on_stdin() {
        // The python image's default command is python3: with no program
        // after the image, it runs what comes on stdin.
        let sh = ". ./versions.env\ndocker run -i --rm \"$@V\" <<'END'\nimport os\nEND\n";
        assert_eq!(k("s.sh", sh), [(2, "call")]);
        assert_eq!(
            k("s.sh", "docker run -i --rm \"${@V}\" < tool.txt\n"),
            [(1, "call")]
        );
        assert_eq!(
            k("s.sh", "podman run -i @p:3.15 2>/dev/null <x\n"),
            [(1, "call")]
        );
        assert_eq!(
            k(
                "s.sh",
                "out=$(docker container run -i docker.io/library/@p:3@sha256:ab < x)\n"
            ),
            [(1, "call")]
        );
        assert_eq!(k("s.sh", "run_in \"$@V\" < x\n"), [(1, "call")]);
        // A continued command: the image ends one line, the redirection
        // starts the next.
        assert_eq!(
            k("s.sh", "docker run -i \\\n  \"$@V\" \\\n  < x\n"),
            [(3, "call")]
        );
        // A program after the image: python3 is counted by its word, bash
        // and sh are no Python, also on the next line of the command.
        assert_eq!(
            k("s.sh", "docker run --rm \"$@V\" @p3 -m cyrius\n"),
            [(1, "call")]
        );
        assert!(
            k(
                "s.sh",
                "run_in -w /tmp \\\n  \"${@V}\" \\\n  bash -c 'og x'\n"
            )
            .is_empty()
        );
        assert!(
            k(
                "s.sh",
                "docker run \"${@V:?is not set in versions.env}\" sh -euc 'ls'\n"
            )
            .is_empty()
        );
        // Not a run of the image: a pull, an echo, a stamp, another image.
        assert!(
            k(
                "s.sh",
                "docker pull \"$@V\"\necho \"$@V\"\nimg=@p:3.15\ndocker run -i \"$OG_IMAGE\" < x\n"
            )
            .is_empty()
        );
        // A Python heredoc to the image is one use, not two.
        assert_eq!(
            k("s.sh", "docker run -i \"$@V\" <<'@P'\n@P\n"),
            [(1, "heredoc")]
        );
        // From a Rust string and a workflow step too.
        let rs = t("fn a() { let _ = \"docker run -i @p:3.15 < x\"; }\n");
        assert_eq!(kinds("x.rs", &rs), [(1, "rust")]);
        assert_eq!(
            k(
                ".github/workflows/x.yml",
                "      - run: docker run -i \"$@V\" < x\n"
            ),
            [(1, "call")]
        );
    }

    #[test]
    fn a_redirection_anywhere_in_a_run_of_the_python_image() {
        // A redirection may stand before the image too.
        assert_eq!(
            k("s.sh", "docker run --rm -i < tool.txt \"$@V\"\n"),
            [(1, "call")]
        );
        assert_eq!(
            k(
                "s.sh",
                "docker run -i 2>/dev/null \"$@V\" >out.txt <in.txt\n"
            ),
            [(1, "call")]
        );
        // A program after a redirection replaces the default command.
        assert!(k("s.sh", "docker run -i \"$@V\" < tool.txt sh\n").is_empty());
        assert!(k("s.sh", "docker run -i \"$@V\" 2>&1 bash -c 'og x'\n").is_empty());
    }

    #[test]
    fn a_continued_command_in_a_rust_string() {
        let rs = t("fn a() { let s = r#\"docker run -i \\\n  \"$@V\" < tool.txt\"#; }\n");
        assert_eq!(kinds("x.rs", &rs), [(2, "rust")]);
        // Each literal starts a command of its own.
        let rs = t("fn a() { let _ = [\"docker run -i \\\\\", \"\\\"$@V\\\" < x\"]; }\n");
        assert!(kinds("x.rs", &rs).is_empty());
    }

    #[test]
    fn a_backslash_before_the_program_still_runs_it() {
        // `\python3` skips an alias and runs python3.
        assert_eq!(k("s.sh", "\\@p3 -c 1\n"), [(1, "call")]);
        assert_eq!(k("s.sh", "x; \\@i install y\n"), [(1, "call")]);
    }

    #[test]
    fn pypy_is_a_python() {
        assert_eq!(
            k("s.sh", "@y3 -c 1\n@y -V\nX=@y3.10\n"),
            [(1, "call"), (2, "call"), (3, "call")]
        );
        assert_eq!(k("bin/tool", "#!/usr/bin/env @y3\n"), [(1, "shebang")]);
        assert!(k("s.sh", "ls lib/@y3.10/site-packages\n").is_empty());
    }

    #[test]
    fn bats_files_are_shell() {
        assert_eq!(
            k("tests/zz.bats", "@test \"x\" {\n  run @p3 -c 1\n}\n"),
            [(2, "call")]
        );
        assert_eq!(
            k("tests/zz", "#!/usr/bin/env bats\nrun @p3 -c 1\n"),
            [(2, "call")]
        );
    }

    #[test]
    fn dockerfiles_by_any_name_are_shell() {
        assert_eq!(k("docker/og.dockerfile", "RUN @p3 -c 1\n"), [(1, "call")]);
        assert_eq!(k("docker/og.Dockerfile", "RUN @p3 -c 1\n"), [(1, "call")]);
        assert_eq!(k("Containerfile.dev", "RUN @p3 -c 1\n"), [(1, "call")]);
    }

    #[test]
    fn actions_anywhere_and_other_forges_are_workflows() {
        let steps = "runs:\n  steps:\n    - name: @p3 in a label\n      run: @p3 -c 1\n";
        for rel in [
            "tools/act/action.yml",
            "action.yaml",
            ".gitea/workflows/a.yml",
            ".forgejo/workflows/a.yaml",
        ] {
            assert_eq!(k(rel, steps), [(4, "call")], "{rel}");
        }
        // Other YAML stays data.
        assert!(k("tools/act/config.yml", steps).is_empty());
    }

    #[test]
    fn snakemake_cython_and_templates_are_python_files() {
        for rel in [
            "Snakefile",
            "workflow/snakefile",
            "rules/a.smk",
            "src/a.pyx",
            "src/a.pxd",
            "bin/tool.py.in",
        ] {
            assert_eq!(k(rel, "x\n"), [(0, "file")], "{rel}");
        }
        assert!(k("bin/tool.in", "x\n").is_empty());
    }

    #[test]
    fn rust_escapes_are_read_as_the_program_gets_them() {
        // `\x33` and `\u{33}` are 3: the program is python3.
        assert_eq!(
            k("x.rs", "fn a() { Command::new(\"@p\\x33\"); }\n"),
            [(1, "rust")]
        );
        assert_eq!(
            k("x.rs", "fn a() { let _ = \"sh -c '@p\\u{33} x'\"; }\n"),
            [(1, "rust")]
        );
        // A `\` at the end of a line joins the next one.
        assert_eq!(
            k("x.rs", "fn a() { let _ = \"sh -c '@p\\\n    3 x'\"; }\n"),
            [(1, "rust")]
        );
        // A raw string keeps its backslashes: no escape runs there.
        assert!(k("x.rs", "fn a() { let _ = r\"@p\\x33\"; }\n").is_empty());
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
