//! `og dev helpers`: the helper binaries (`bcftools`, `bgzip`, `tabix`,
//! `samtools`) each script and module calls, and the image each call runs
//! in. `scripts/ci/check-container-helpers.sh` runs it; that wrapper keeps
//! the docker probe (rule 2), the orad download and the self-test.
//!
//! The table is derived, never kept by hand:
//!
//! - `scripts/*.sh`: every `docker run … IMAGE cmd` command and every call
//!   of the wrapper `run_in [--net] [--root] [--rw DIR] … IMAGE cmd`. The
//!   in-image commands are `cmd`, or for `bash -c '…'` the command word of
//!   every command in the quoted body.
//! - `modules/local/*/main.nf`: every process, its `container` (or its
//!   `withName` selector in `conf/containers.config`) and the command word
//!   of every command in the triple-quoted strings of its `script:`,
//!   `shell:` and `exec:` sections. The `stub:` section never runs a tool.
//!
//! A call by path (`"$d/bgzip"`) runs that file, not the helper on `PATH`,
//! and is left out. An image word is named by its `versions.env` variable
//! when one matches; `versions.env` is read, never executed.
//!
//! Modes (each takes `--root DIR`, default the current folder):
//!
//! - `--list`: the table, one `VAR bin file:line` row per call.
//! - `--static`: rule 1, `bgzip` and `tabix` are never called in
//!   `BCFTOOLS_IMAGE` or `SAMTOOLS_IMAGE`; and every call's image is known.
//! - `--images`: the table grouped by image for the docker probe, as tab
//!   separated `image NAME VALUE BINS` and `call NAME VALUE BIN LOCATION`
//!   rows.
//! - `--orad-calls`: is there an `"$ORAD" …` call in `scripts/*.sh`?
//! - `--orad --help-file F [--orad-version V]`: rule 3, every option a
//!   script passes to `"$ORAD"` is listed in orad's `--help` text `F`.
//!
//! Exit 0 clean, 1 a finding, 2 a usage error or a tree the readers no
//! longer understand (no helper call derived from the scripts, no orad
//! call, a `versions.env` line that is not `NAME="value"`): a check that
//! checked nothing never passes.
//!
//! Where it differs from the Python it replaces, on purpose: helper calls
//! in `shell:` and `exec:` sections are read (the Python read only
//! `script:`); a here-document ends only on its exact delimiter line, as in
//! bash (the Python also ended it on a line with blanks around the
//! delimiter); two `*_IMAGE` variables with the same value are both named,
//! joined by `|` (the Python kept the alphabetically last).

use crate::shlex::{self, Command};
use regex::Regex;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::LazyLock;

const USAGE: &str = "usage: og dev helpers --list|--static|--images|--orad-calls [--root DIR]\n       og dev helpers --orad --help-file FILE [--orad-version V] [--root DIR]\n";

/// The helper binaries this check follows.
const HELPERS: [&str; 4] = ["bcftools", "bgzip", "tabix", "samtools"];
/// Images that must never be asked for `bgzip` or `tabix` (rule 1).
const NO_BGZIP_TABIX: [&str; 2] = ["BCFTOOLS_IMAGE", "SAMTOOLS_IMAGE"];

/// `docker run` options that take a value as the next word.
const VALUED: &[&str] = &[
    "-a",
    "--attach",
    "--add-host",
    "-c",
    "--cpu-shares",
    "--cap-add",
    "--cap-drop",
    "--cidfile",
    "--cpus",
    "--cpuset-cpus",
    "--device",
    "--dns",
    "-e",
    "--env",
    "--env-file",
    "--entrypoint",
    "--expose",
    "--gpus",
    "--group-add",
    "-h",
    "--hostname",
    "--ipc",
    "-l",
    "--label",
    "--label-file",
    "--link",
    "--log-driver",
    "--log-opt",
    "-m",
    "--memory",
    "--memory-reservation",
    "--memory-swap",
    "--mount",
    "--name",
    "--net",
    "--network",
    "-p",
    "--publish",
    "--pid",
    "--platform",
    "--pull",
    "--restart",
    "--runtime",
    "--security-opt",
    "--shm-size",
    "--stop-signal",
    "--tmpfs",
    "-u",
    "--user",
    "--ulimit",
    "--userns",
    "--uts",
    "-v",
    "--volume",
    "--volumes-from",
    "-w",
    "--workdir",
];

fn re(s: &str) -> Regex {
    Regex::new(s).expect("a valid built-in pattern")
}

static VERSIONS_LINE: LazyLock<Regex> =
    LazyLock::new(|| re(r#"^([A-Z][A-Z0-9_]*)="([^"$`]*)"\s*(?:#.*)?$"#));
static SELECTOR: LazyLock<Regex> = LazyLock::new(|| {
    re(r#"withName:\s*['"]?(\w+)['"]?\s*\{[^}]*?container\s*=\s*['"]([^'"]+)['"]"#)
});
static PROCESS: LazyLock<Regex> = LazyLock::new(|| re(r"(?m)^process\s+(\w+)\s*\{"));
static CONTAINER: LazyLock<Regex> = LazyLock::new(|| re(r#"(?m)^\s*container\s+['"]([^'"]+)['"]"#));
static SECTION: LazyLock<Regex> = LazyLock::new(|| re(r"(?m)^\s*(script|shell|exec|stub):\s*$"));
static PLAIN: LazyLock<Regex> = LazyLock::new(|| re(r"^[A-Za-z0-9._/:@+-]+$"));
static ORAD_HELP: LazyLock<Regex> = LazyLock::new(|| re(r"(?m)^\s*-h, --help\b"));

/// The `*_IMAGE` variables of `versions.env`, and each value's names.
struct Versions {
    by_name: BTreeMap<String, String>,
    by_value: BTreeMap<String, String>,
}

impl Versions {
    /// Reads `versions.env` strictly: every line that is neither blank nor a
    /// comment must be `NAME="value"` with no expansion in the value, or the
    /// file is refused (nothing is executed, so nothing may need to be).
    fn read(path: &Path) -> Result<Self, String> {
        let text =
            fs::read_to_string(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
        let mut by_name = BTreeMap::new();
        for (n, line) in text.lines().enumerate() {
            let t = line.trim();
            if t.is_empty() || t.starts_with('#') {
                continue;
            }
            let c = VERSIONS_LINE.captures(line).ok_or_else(|| {
                format!(
                    "{}:{}: not NAME=\"value\" (no $, ` or quote in the value): {line}",
                    path.display(),
                    n + 1
                )
            })?;
            if c[1].ends_with("_IMAGE") {
                by_name.insert(c[1].to_owned(), c[2].to_owned());
            }
        }
        let mut names: BTreeMap<String, Vec<&str>> = BTreeMap::new();
        for (k, v) in &by_name {
            names.entry(v.clone()).or_default().push(k);
        }
        let by_value = names.into_iter().map(|(v, ks)| (v, ks.join("|"))).collect();
        Ok(Versions { by_name, by_value })
    }

    /// The variable name(s) of an image value, or the value itself.
    fn name_of(&self, value: &str) -> String {
        self.by_value
            .get(value)
            .cloned()
            .unwrap_or_else(|| value.to_owned())
    }
}

/// One helper call: the image's name and value, the helper, where.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct Call {
    name: String,
    value: String,
    bin: String,
    loc: String,
}

/// What the scripts and modules of a tree call.
#[derive(Default)]
struct Table {
    calls: Vec<Call>,
    /// Locations whose image could not be determined.
    unresolved: Vec<String>,
    has_scripts: bool,
}

impl Table {
    /// No call from `scripts/` was derived although there are scripts: the
    /// parser no longer sees how the scripts start a container.
    fn blind(&self) -> bool {
        self.has_scripts
            && !self.calls.iter().any(|c| c.loc.starts_with("scripts/"))
            && !self.unresolved.iter().any(|l| l.starts_with("scripts/"))
    }
}

/// Python's `\w` for one character.
fn is_word(c: char) -> bool {
    c.is_alphanumeric() || c == '_'
}

/// The first offset at or after `from` where `token` stands with none of
/// `before` (and no word character) right before it, and no word character
/// or `-` right after it.
fn find_token(text: &str, token: &str, from: usize, before: &[char]) -> Option<usize> {
    let hay = text.get(from..)?;
    hay.match_indices(token).map(|(i, _)| from + i).find(|&i| {
        let prev = text[..i].chars().next_back();
        let next = text[i + token.len()..].chars().next();
        prev.is_none_or(|c| !is_word(c) && !before.contains(&c))
            && next.is_none_or(|c| !is_word(c) && c != '-')
    })
}

fn line_of(text: &str, offset: usize) -> usize {
    text.as_bytes()[..offset.min(text.len())]
        .iter()
        .filter(|&&c| c == b'\n')
        .count()
        + 1
}

/// The files `dir/*SUFFIX`, sorted by name, hidden ones left out, as glob
/// lists them.
fn files(dir: &Path, suffix: &str) -> Vec<PathBuf> {
    let Ok(rd) = fs::read_dir(dir) else {
        return Vec::new();
    };
    let mut out: Vec<PathBuf> = rd
        .filter_map(Result::ok)
        .filter(|e| {
            let n = e.file_name();
            let n = n.to_string_lossy();
            !n.starts_with('.') && n.ends_with(suffix) && e.path().is_file()
        })
        .map(|e| e.path())
        .collect();
    out.sort();
    out
}

/// `modules/local/*/main.nf`, sorted by module folder.
fn module_files(root: &Path) -> Vec<PathBuf> {
    let Ok(rd) = fs::read_dir(root.join("modules").join("local")) else {
        return Vec::new();
    };
    let mut out: Vec<PathBuf> = rd
        .filter_map(Result::ok)
        .filter(|e| !e.file_name().to_string_lossy().starts_with('.'))
        .map(|e| e.path().join("main.nf"))
        .filter(|p| p.is_file())
        .collect();
    out.sort();
    out
}

fn read_lossy(path: &Path) -> Result<String, String> {
    fs::read(path)
        .map(|b| String::from_utf8_lossy(&b).into_owned())
        .map_err(|e| format!("cannot read {}: {e}", path.display()))
}

/// The words of a command after its keywords and assignments.
fn command_words(cmd: &Command) -> &[String] {
    let k = cmd
        .words
        .iter()
        .position(|w| !shlex::is_keyword(w) && !shlex::is_assignment(w))
        .unwrap_or(cmd.words.len());
    &cmd.words[k..]
}

/// `(image word, in-image argv)` for `docker run …` or `run_in …`.
fn docker_run(words: &[String]) -> Option<(&str, Vec<String>)> {
    let mut k = match words {
        [d, r, ..] if d == "docker" && r == "run" => 2,
        [w, ..] if w == "run_in" => {
            let mut k = 1;
            while let Some(w) = words.get(k) {
                match w.as_str() {
                    "--rw" => k += 2,
                    "--net" | "--root" => k += 1,
                    _ => break,
                }
            }
            k
        }
        _ => return None,
    };
    let mut entry = Vec::new();
    while let Some(w) = words.get(k).filter(|w| w.starts_with('-')) {
        if w == "--entrypoint"
            && let Some(e) = words.get(k + 1)
        {
            entry = vec![e.clone()];
        }
        k += if VALUED.contains(&w.as_str()) && !w.contains('=') {
            2
        } else {
            1
        };
    }
    let image = words.get(k)?;
    entry.extend(words[k + 1..].iter().cloned());
    Some((image.as_str(), entry))
}

/// The binaries a container command runs: `argv[0]`, or every command word
/// of a `bash -c` / `sh -c` body. Option groups such as `-euo pipefail` are
/// understood: each `o` or `O` takes the next word, and a `c` anywhere in a
/// group makes the first word after the options the command body.
fn in_image_heads(argv: &[String]) -> Vec<String> {
    let Some(first) = argv.first() else {
        return Vec::new();
    };
    let base = first.rsplit('/').next().unwrap_or(first);
    if base != "bash" && base != "sh" {
        return if first.contains('/') {
            Vec::new()
        } else {
            vec![first.clone()]
        };
    }
    let (mut k, mut want_c, mut takes) = (1, false, 0usize);
    while let Some(w) = argv.get(k) {
        if takes > 0 {
            takes -= 1;
        } else if w == "--" {
            k += 1;
            break;
        } else if w.starts_with("--") {
        } else if w.len() > 1 && (w.starts_with('-') || w.starts_with('+')) {
            let opts = &w[1..];
            want_c = want_c || opts.contains('c');
            takes = opts.matches(['o', 'O']).count();
        } else {
            break;
        }
        k += 1;
    }
    match argv.get(k) {
        Some(body) if want_c => shlex::heads(&shlex::lex(body)),
        _ => Vec::new(),
    }
}

/// The first `$NAME` or `${NAME` in a word.
fn first_var(word: &str) -> Option<&str> {
    let b = word.as_bytes();
    for (i, _) in word.match_indices('$') {
        let mut s = i + 1;
        if b.get(s) == Some(&b'{') {
            s += 1;
        }
        if !b
            .get(s)
            .is_some_and(|c| c.is_ascii_alphabetic() || *c == b'_')
        {
            continue;
        }
        let mut e = s + 1;
        while e < b.len() && (b[e].is_ascii_alphanumeric() || b[e] == b'_') {
            e += 1;
        }
        return Some(&word[s..e]);
    }
    None
}

/// The value of `NAME=VALUE` (after optional `export`, `local` or
/// `readonly`) at the start of a line: the word up to the first blank.
fn assigned<'a>(line: &'a str, name: &str) -> Option<&'a str> {
    let at = |s: &'a str| {
        s.strip_prefix(name)?
            .strip_prefix('=')
            .map(|v| v.split(char::is_whitespace).next().unwrap_or(""))
    };
    let s = line.trim_start();
    for kw in ["export", "local", "readonly"] {
        if let Some(r) = s.strip_prefix(kw)
            && r.starts_with(char::is_whitespace)
            && let Some(v) = at(r.trim_start())
        {
            return Some(v);
        }
    }
    at(s)
}

/// `(name, value)` for an image word of a script. A script-local
/// assignment counts only when its value is a plain literal or one
/// `versions.env` variable; anything else leaves the value unknown.
fn resolve(word: &str, script: &str, v: &Versions) -> (String, Option<String>) {
    let Some(name) = first_var(word) else {
        return (v.name_of(word), Some(word.to_owned()));
    };
    if let Some(val) = v.by_name.get(name) {
        return (name.to_owned(), Some(val.clone()));
    }
    for line in script.lines() {
        let Some(raw) = assigned(line, name) else {
            continue;
        };
        let mut val = raw.trim_end_matches(';');
        let b = val.as_bytes();
        if b.len() >= 2 && b[0] == b[b.len() - 1] && (b[0] == b'"' || b[0] == b'\'') {
            val = &val[1..val.len() - 1];
        }
        let inner = val.strip_prefix('$').map(|r| {
            let r = r.strip_prefix('{').unwrap_or(r);
            r.strip_suffix('}').unwrap_or(r)
        });
        if let Some(r) = inner
            && !r.is_empty()
            && r.bytes()
                .next()
                .is_some_and(|c| c.is_ascii_alphabetic() || c == b'_')
            && r.bytes().all(|c| c.is_ascii_alphanumeric() || c == b'_')
            && let Some(val) = v.by_name.get(r)
        {
            return (r.to_owned(), Some(val.clone()));
        }
        if PLAIN.is_match(val) {
            return (name.to_owned(), Some(val.to_owned()));
        }
        return (name.to_owned(), None);
    }
    (name.to_owned(), None)
}

/// A Groovy string's text as the shell sees it: escapes resolved, and with
/// `interpolate` every `${…}` or `$name` made the word `GROOVY`.
fn groovy_to_shell(s: &str, interpolate: bool) -> String {
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        let c = b[i];
        if c == b'\\' && i + 1 < b.len() {
            let nxt = b[i + 1];
            if b"$\\\"'".contains(&nxt) {
                out.push(nxt);
            } else {
                out.extend_from_slice(&[c, nxt]);
            }
            i += 2;
            continue;
        }
        if c == b'$' && interpolate {
            if b.get(i + 1) == Some(&b'{') {
                i = shlex::match_brace(b, i + 1) + 1;
                out.extend_from_slice(b"GROOVY");
                continue;
            }
            if b.get(i + 1)
                .is_some_and(|c| c.is_ascii_alphabetic() || *c == b'_')
            {
                i += 2;
                while i < b.len() && (b[i].is_ascii_alphanumeric() || b"_.".contains(&b[i])) {
                    i += 1;
                }
                out.extend_from_slice(b"GROOVY");
                continue;
            }
        }
        out.push(c);
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// The triple-quoted strings of `s`: `(interpolating, content start,
/// content end)`, leftmost first.
fn triple_quoted(s: &str) -> Vec<(bool, usize, usize)> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < s.len() {
        let rest = &s[i..];
        let delim = ["\"\"\"", "'''"].into_iter().find(|d| rest.starts_with(d));
        if let Some(d) = delim
            && let Some(close) = s[i + 3..].find(d)
        {
            out.push((d == "\"\"\"", i + 3, i + 3 + close));
            i += 3 + close + 3;
            continue;
        }
        i += rest.chars().next().map_or(1, char::len_utf8);
    }
    out
}

fn collect(root: &Path, v: &Versions) -> Result<Table, String> {
    let mut table = Table::default();
    let scripts = files(&root.join("scripts"), ".sh");
    table.has_scripts = !scripts.is_empty();
    for path in &scripts {
        let rel = format!(
            "scripts/{}",
            path.file_name().unwrap_or_default().to_string_lossy()
        );
        let text = read_lossy(path)?;
        for cmd in shlex::lex(&text) {
            let Some((image, argv)) = docker_run(command_words(&cmd)) else {
                continue;
            };
            let (name, value) = resolve(image, &text, v);
            for bin in in_image_heads(&argv) {
                if !HELPERS.contains(&bin.as_str()) {
                    continue;
                }
                let at = find_token(&text, &bin, cmd.offset, &['/', '-']).unwrap_or(cmd.offset);
                let loc = format!("{rel}:{}", line_of(&text, at));
                match &value {
                    None => table.unresolved.push(format!("{loc} ({image})")),
                    Some(value) => table.calls.push(Call {
                        name: name.clone(),
                        value: value.clone(),
                        bin,
                        loc,
                    }),
                }
            }
        }
    }

    let mut selectors = BTreeMap::new();
    let cfg = root.join("conf").join("containers.config");
    if cfg.exists() {
        let text = read_lossy(&cfg)?;
        for c in SELECTOR.captures_iter(&text) {
            selectors.insert(c[1].to_owned(), c[2].to_owned());
        }
    }

    for path in module_files(root) {
        let folder = path
            .parent()
            .and_then(Path::file_name)
            .unwrap_or_default()
            .to_string_lossy()
            .into_owned();
        let rel = format!("modules/local/{folder}/main.nf");
        let text = read_lossy(&path)?;
        let procs: Vec<(usize, usize, String)> = PROCESS
            .captures_iter(&text)
            .map(|c| {
                let m = c.get(0).expect("match");
                (m.start(), m.end(), c[1].to_owned())
            })
            .collect();
        for (k, (_, start, pname)) in procs.iter().enumerate() {
            let end = procs.get(k + 1).map_or(text.len(), |p| p.0);
            let body = &text[*start..end];
            // A withName selector in the config overrides the directive.
            let image = selectors
                .get(pname)
                .cloned()
                .or_else(|| CONTAINER.captures(body).map(|c| c[1].to_owned()));
            let heads: Vec<(&str, usize, usize)> = SECTION
                .captures_iter(body)
                .map(|c| {
                    let m = c.get(0).expect("match");
                    (c.get(1).expect("kind").as_str(), m.start(), m.end())
                })
                .collect();
            for (j, (kind, _, sec)) in heads.iter().enumerate() {
                if *kind == "stub" {
                    continue;
                }
                let sec_start = start + sec;
                let sec_end = start + heads.get(j + 1).map_or(body.len(), |h| h.1);
                let section = &text[sec_start..sec_end];
                for (interp, qs, qe) in triple_quoted(section) {
                    let shell = groovy_to_shell(&section[qs..qe], interp);
                    let qstart = sec_start + qs;
                    for bin in shlex::heads(&shlex::lex(&shell)) {
                        if !HELPERS.contains(&bin.as_str()) {
                            continue;
                        }
                        let at = find_token(&text[..sec_end], &bin, qstart, &['/', '-'])
                            .unwrap_or(qstart);
                        let loc = format!("{rel}:{} ({pname})", line_of(&text, at));
                        match &image {
                            None => table.unresolved.push(loc),
                            Some(img) => table.calls.push(Call {
                                name: v.name_of(img),
                                value: img.clone(),
                                bin,
                                loc,
                            }),
                        }
                    }
                }
            }
        }
    }

    let mut seen = BTreeSet::new();
    table.calls.retain(|c| seen.insert(c.clone()));
    Ok(table)
}

const BLIND: &str = "ERROR: derived no helper call in scripts/*.sh; the parser does not see how the scripts start containers.";

fn list(t: &Table, out: &mut dyn Write) -> std::io::Result<u8> {
    for c in &t.calls {
        writeln!(out, "{:<16} {:<9} {}", c.name, c.bin, c.loc)?;
    }
    if t.blind() {
        writeln!(out, "{BLIND}")?;
        return Ok(2);
    }
    Ok(0)
}

/// The derived table is usable: some call was derived, and from the
/// scripts too. Prints why not and returns `false` otherwise.
fn usable(t: &Table, out: &mut dyn Write) -> std::io::Result<bool> {
    if t.blind() {
        writeln!(out, "{BLIND}")?;
        return Ok(false);
    }
    if t.calls.is_empty() {
        writeln!(
            out,
            "ERROR: derived no helper calls at all; the parser is broken."
        )?;
        return Ok(false);
    }
    Ok(true)
}

fn rule1(t: &Table, out: &mut dyn Write) -> std::io::Result<u8> {
    if !usable(t, out)? {
        return Ok(2);
    }
    let mut failed = false;
    if !t.unresolved.is_empty() {
        failed = true;
        writeln!(
            out,
            "FAIL: helper calls whose image could not be determined:"
        )?;
        for loc in &t.unresolved {
            writeln!(out, "  {loc}")?;
        }
    }
    let bad: Vec<&Call> = t
        .calls
        .iter()
        .filter(|c| {
            (c.bin == "bgzip" || c.bin == "tabix")
                && c.name.split('|').any(|n| NO_BGZIP_TABIX.contains(&n))
        })
        .collect();
    if bad.is_empty() {
        writeln!(
            out,
            "OK: no bgzip/tabix call in {}.",
            NO_BGZIP_TABIX.join(" or ")
        )?;
    } else {
        failed = true;
        writeln!(
            out,
            "FAIL: bgzip/tabix called in an image that is not meant to provide them"
        )?;
        writeln!(
            out,
            "      (use `bcftools view -Oz -o X.vcf.gz` and `bcftools index -t X.vcf.gz`):"
        )?;
        for c in bad {
            writeln!(out, "  {:<9} in {:<15} {}", c.bin, c.name, c.loc)?;
        }
    }
    Ok(u8::from(failed))
}

fn images(t: &Table, out: &mut dyn Write) -> std::io::Result<u8> {
    if !usable(t, out)? {
        return Ok(2);
    }
    let mut by_image: BTreeMap<(&str, &str), BTreeMap<&str, Vec<&str>>> = BTreeMap::new();
    for c in &t.calls {
        by_image
            .entry((&c.name, &c.value))
            .or_default()
            .entry(&c.bin)
            .or_default()
            .push(&c.loc);
    }
    for ((name, value), bins) in &by_image {
        let names: Vec<&str> = bins.keys().copied().collect();
        writeln!(out, "image\t{name}\t{value}\t{}", names.join(" "))?;
        for (bin, locs) in bins {
            for loc in locs {
                writeln!(out, "call\t{name}\t{value}\t{bin}\t{loc}")?;
            }
        }
    }
    Ok(0)
}

/// `(location, options)` for every `"$ORAD" …` command of `scripts/*.sh`.
fn orad_calls(root: &Path) -> Result<Vec<(String, Vec<String>)>, String> {
    let mut found = Vec::new();
    for path in files(&root.join("scripts"), ".sh") {
        let rel = format!(
            "scripts/{}",
            path.file_name().unwrap_or_default().to_string_lossy()
        );
        let text = read_lossy(&path)?;
        for cmd in shlex::lex(&text) {
            let words = command_words(&cmd);
            let Some(first) = words.first() else { continue };
            let bare = first.strip_prefix('$').map(|r| {
                let r = r.strip_prefix('{').unwrap_or(r);
                r.strip_suffix('}').unwrap_or(r)
            });
            if bare != Some("ORAD") {
                continue;
            }
            let opts = words[1..]
                .iter()
                .filter(|w| w.len() > 1 && w.starts_with('-'))
                .map(|w| w.split('=').next().unwrap_or(w).to_owned())
                .collect();
            found.push((format!("{rel}:{}", line_of(&text, cmd.offset)), opts));
        }
    }
    Ok(found)
}

const NO_ORAD: &str = "ERROR: found no \"$ORAD\" call in scripts/*.sh; the parser does not see how step 01 runs orad.";

fn orad(root: &Path, help_file: &Path, version: &str, out: &mut dyn Write) -> Result<u8, String> {
    let calls = orad_calls(root)?;
    let w = |e: std::io::Error| e.to_string();
    if calls.is_empty() {
        writeln!(out, "{NO_ORAD}").map_err(w)?;
        return Ok(2);
    }
    let text = read_lossy(help_file)?;
    // An empty or foreign help text must not pass as "every option listed".
    if !ORAD_HELP.is_match(&text) {
        writeln!(out, "FAIL: orad --help printed no option list:").map_err(w)?;
        writeln!(out, "{text}").map_err(w)?;
        return Ok(1);
    }
    let mut failed = false;
    for (loc, opts) in &calls {
        for o in opts {
            if find_token(&text, o, 0, &['-']).is_some() {
                writeln!(out, "OK: {o:<20} {loc}").map_err(w)?;
            } else {
                failed = true;
                writeln!(
                    out,
                    "FAIL: {loc} passes {o}, which orad {version} does not list in --help"
                )
                .map_err(w)?;
            }
        }
    }
    if failed {
        writeln!(out, "orad {version} --help:").map_err(w)?;
        writeln!(out, "{text}").map_err(w)?;
    }
    Ok(u8::from(failed))
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    List,
    Static,
    Images,
    OradCalls,
    Orad,
}

/// `og dev helpers ARGS`; returns the exit code.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let usage = |err: &mut dyn Write, msg: &str| {
        let _ = write!(err, "{USAGE}");
        let _ = writeln!(err, "og dev helpers: error: {msg}");
        2
    };
    let mut mode = None;
    let mut root = PathBuf::from(".");
    let mut help_file = None;
    let mut version = String::from("?");
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let m = match a.as_str() {
            "--list" => Mode::List,
            "--static" => Mode::Static,
            "--images" => Mode::Images,
            "--orad-calls" => Mode::OradCalls,
            "--orad" => Mode::Orad,
            "--root" | "--help-file" | "--orad-version" => {
                let Some(v) = it.next() else {
                    return usage(err, &format!("{a} needs a value"));
                };
                match a.as_str() {
                    "--root" => root = PathBuf::from(v),
                    "--help-file" => help_file = Some(PathBuf::from(v)),
                    _ => version.clone_from(v),
                }
                continue;
            }
            "-h" | "--help" => {
                let _ = write!(out, "{USAGE}");
                return 0;
            }
            other => return usage(err, &format!("unknown argument '{other}'")),
        };
        if mode.replace(m).is_some_and(|old| old != m) {
            return usage(err, "give one mode");
        }
    }
    let Some(mode) = mode else {
        return usage(err, "give one mode");
    };
    let fail = |err: &mut dyn Write, msg: String| {
        let _ = writeln!(err, "og dev helpers: error: {msg}");
        2
    };
    let io = |r: std::io::Result<u8>| r.unwrap_or(2);
    match mode {
        Mode::OradCalls => match orad_calls(&root) {
            Ok(c) if c.is_empty() => {
                let _ = writeln!(out, "{NO_ORAD}");
                2
            }
            Ok(_) => 0,
            Err(e) => fail(err, e),
        },
        Mode::Orad => {
            let Some(help) = help_file else {
                return usage(err, "--orad needs --help-file FILE");
            };
            orad(&root, &help, &version, out).unwrap_or_else(|e| fail(err, e))
        }
        Mode::List | Mode::Static | Mode::Images => {
            let table = Versions::read(&root.join("versions.env")).and_then(|v| collect(&root, &v));
            match table {
                Err(e) => fail(err, e),
                Ok(t) => io(match mode {
                    Mode::List => list(&t, out),
                    Mode::Static => rule1(&t, out),
                    _ => images(&t, out),
                }),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tree(files: &[(&str, &str)]) -> tempdir::Dir {
        let d = tempdir::Dir::new();
        for (p, text) in files {
            let path = d.0.join(p);
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            fs::write(path, text).unwrap();
        }
        d
    }

    mod tempdir {
        use std::path::PathBuf;
        use std::sync::atomic::{AtomicUsize, Ordering};
        pub struct Dir(pub PathBuf);
        static N: AtomicUsize = AtomicUsize::new(0);
        impl Dir {
            pub fn new() -> Self {
                let p = std::env::temp_dir().join(format!(
                    "og-helpers-{}-{}",
                    std::process::id(),
                    N.fetch_add(1, Ordering::SeqCst)
                ));
                std::fs::create_dir_all(&p).unwrap();
                Dir(p)
            }
        }
        impl Drop for Dir {
            fn drop(&mut self) {
                let _ = std::fs::remove_dir_all(&self.0);
            }
        }
    }

    const VERSIONS: &str = "# images\nBCFTOOLS_IMAGE=\"example/bcftools:1.0\"\nSAMTOOLS_IMAGE=\"example/samtools:1.0\"  # pinned\nTOOL_IMAGE=\"example/tool:1.0\"\nOTHER=\"x\"\n";

    fn run(args: &[&str], root: &Path) -> (u8, String) {
        let mut a: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
        a.push("--root".into());
        a.push(root.to_string_lossy().into_owned());
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = cli(&a, &mut out, &mut err);
        let mut s = String::from_utf8(out).unwrap();
        s.push_str(&String::from_utf8(err).unwrap());
        (code, s)
    }

    fn module(name: &str, image: &str, section: &str, body: &str) -> String {
        format!(
            "process {name} {{\n    container '{image}'\n    {section}:\n    \"\"\"\n    {body}\n    \"\"\"\n}}\n"
        )
    }

    #[test]
    fn script_calls_behind_every_prefix() {
        let script = concat!(
            "#!/usr/bin/env bash\n",
            "if ! docker run --rm -v \"$G:/genome\" \"${BCFTOOLS_IMAGE}\" bash -c \"bcftools view x | bgzip -c > y.gz\"; then\n",
            "OUT=$(time docker run --rm \"${BCFTOOLS_IMAGE}\" tabix z.gz)\n",
            "VAR=x docker run --rm \"${SAMTOOLS_IMAGE}\" bgzip w\n",
            "docker run --rm \"${BCFTOOLS_IMAGE}\" bash -euo pipefail -c \"bgzip c.vcf\"\n",
            "IMG=\"${BCFTOOLS_IMAGE}\"\n",
            "docker run --rm \"$IMG\" tabix i.vcf.gz\n",
            "run_in --net --root --rw \"$G/ref\" --cpus 2 \"${SAMTOOLS_IMAGE}\" tabix w.vcf.gz\n",
            "run_in \"${BCFTOOLS_IMAGE}\" bash -c \"\\\"\\$d/bgzip\\\" -c by_path.vcf\"\n",
            "docker run --rm \"${TOOL_IMAGE}\" bgzip allowed.vcf\n",
        );
        let d = tree(&[("versions.env", VERSIONS), ("scripts/a.sh", script)]);
        let (code, out) = run(&["--list"], &d.0);
        assert_eq!(code, 0, "{out}");
        let rows: Vec<&str> = out.lines().collect();
        let want = [
            "BCFTOOLS_IMAGE   bcftools  scripts/a.sh:2",
            "BCFTOOLS_IMAGE   bgzip     scripts/a.sh:2",
            "BCFTOOLS_IMAGE   tabix     scripts/a.sh:3",
            "SAMTOOLS_IMAGE   bgzip     scripts/a.sh:4",
            "BCFTOOLS_IMAGE   bgzip     scripts/a.sh:5",
            "BCFTOOLS_IMAGE   tabix     scripts/a.sh:7",
            "SAMTOOLS_IMAGE   tabix     scripts/a.sh:8",
            "TOOL_IMAGE       bgzip     scripts/a.sh:10",
        ];
        assert_eq!(rows, want);
        let (code, out) = run(&["--static"], &d.0);
        assert_eq!(code, 1, "{out}");
        assert!(
            out.contains("  bgzip     in SAMTOOLS_IMAGE  scripts/a.sh:4\n"),
            "{out}"
        );
        assert!(!out.contains("TOOL_IMAGE"), "{out}");
    }

    #[test]
    fn shell_and_exec_sections_are_read_and_stub_is_not() {
        let d = tree(&[
            ("versions.env", VERSIONS),
            (
                "scripts/a.sh",
                "docker run --rm \"${TOOL_IMAGE}\" bcftools view x\n",
            ),
            (
                "modules/local/s/main.nf",
                &module(
                    "S",
                    "example/bcftools:1.0",
                    "shell",
                    "bgzip -c !{vcf} > o.gz",
                ),
            ),
            (
                "modules/local/e/main.nf",
                &module("E", "example/samtools:1.0", "exec", "tabix w.vcf.gz"),
            ),
            (
                "modules/local/x/main.nf",
                "process X {\n    container 'example/bcftools:1.0'\n    script:\n    \"\"\"\n    bcftools view ${vcf} -Oz -o o.gz\n    \"\"\"\n    stub:\n    \"\"\"\n    tabix stub_only.vcf.gz\n    \"\"\"\n}\n",
            ),
        ]);
        let (code, out) = run(&["--static"], &d.0);
        assert_eq!(code, 1, "{out}");
        assert!(
            out.contains("bgzip     in BCFTOOLS_IMAGE  modules/local/s/main.nf:5 (S)"),
            "{out}"
        );
        assert!(
            out.contains("tabix     in SAMTOOLS_IMAGE  modules/local/e/main.nf:5 (E)"),
            "{out}"
        );
        assert!(
            !out.contains("stub_only") && !out.contains(":9 (X)"),
            "{out}"
        );
        let (_, list) = run(&["--list"], &d.0);
        assert!(
            list.contains("BCFTOOLS_IMAGE   bcftools  modules/local/x/main.nf:5 (X)"),
            "{list}"
        );
    }

    #[test]
    fn a_selector_overrides_the_directive() {
        let d = tree(&[
            ("versions.env", VERSIONS),
            ("scripts/a.sh", "run_in \"${TOOL_IMAGE}\" samtools view x\n"),
            (
                "conf/containers.config",
                "process { withName: 'Y' { container = 'example/bcftools:1.0' } }\n",
            ),
            (
                "modules/local/y/main.nf",
                &module("Y", "example/tool:1.0", "script", "tabix y.vcf.gz"),
            ),
        ]);
        let (code, out) = run(&["--static"], &d.0);
        assert_eq!(code, 1);
        assert!(
            out.contains("tabix     in BCFTOOLS_IMAGE  modules/local/y/main.nf:5 (Y)"),
            "{out}"
        );
    }

    #[test]
    fn unknown_images_and_blind_trees() {
        let d = tree(&[
            ("versions.env", VERSIONS),
            (
                "scripts/a.sh",
                "ODD=$(pick)\ndocker run --rm \"$ODD\" samtools view x\nrun_in \"${TOOL_IMAGE}\" tabix y\n",
            ),
        ]);
        let (code, out) = run(&["--static"], &d.0);
        assert_eq!(code, 1);
        assert!(out.contains("  scripts/a.sh:2 ($ODD)\n"), "{out}");

        let blind = tree(&[
            ("versions.env", VERSIONS),
            (
                "scripts/a.sh",
                "start_container \"${BCFTOOLS_IMAGE}\" bcftools view x\n",
            ),
            (
                "modules/local/x/main.nf",
                &module("X", "example/bcftools:1.0", "script", "bcftools view x"),
            ),
        ]);
        for mode in ["--static", "--list", "--images"] {
            let (code, out) = run(&[mode], &blind.0);
            assert_eq!(code, 2, "{mode}: {out}");
            assert!(out.contains("no helper call in scripts/"), "{out}");
        }
    }

    #[test]
    fn versions_env_is_read_strictly() {
        let d = tree(&[
            ("versions.env", "A_IMAGE=\"x/a:1\"\nB_IMAGE=\"$A_IMAGE\"\n"),
            ("scripts/a.sh", "run_in \"${A_IMAGE}\" bcftools view\n"),
        ]);
        let (code, out) = run(&["--list"], &d.0);
        assert_eq!(code, 2);
        assert!(out.contains("versions.env:2: not NAME=\"value\""), "{out}");
    }

    #[test]
    fn a_shared_image_value_names_both_variables() {
        let d = tree(&[
            (
                "versions.env",
                "BCFTOOLS_IMAGE=\"x/b:1\"\nZ_IMAGE=\"x/b:1\"\n",
            ),
            ("scripts/a.sh", "docker run --rm x/b:1 tabix y\n"),
        ]);
        let (code, out) = run(&["--static"], &d.0);
        assert_eq!(code, 1, "{out}");
        assert!(
            out.contains("tabix     in BCFTOOLS_IMAGE|Z_IMAGE scripts/a.sh:1"),
            "{out}"
        );
    }

    #[test]
    fn images_groups_the_calls() {
        let d = tree(&[
            ("versions.env", VERSIONS),
            (
                "scripts/a.sh",
                "run_in \"${TOOL_IMAGE}\" tabix a\nrun_in \"${TOOL_IMAGE}\" bgzip b\nrun_in \"${TOOL_IMAGE}\" tabix c\n",
            ),
        ]);
        let (code, out) = run(&["--images"], &d.0);
        assert_eq!(code, 0);
        assert_eq!(
            out,
            "image\tTOOL_IMAGE\texample/tool:1.0\tbgzip tabix\n\
             call\tTOOL_IMAGE\texample/tool:1.0\tbgzip\tscripts/a.sh:2\n\
             call\tTOOL_IMAGE\texample/tool:1.0\ttabix\tscripts/a.sh:1\n\
             call\tTOOL_IMAGE\texample/tool:1.0\ttabix\tscripts/a.sh:3\n"
        );
    }

    #[test]
    fn orad_options_against_help() {
        let help = "-P  path, --path  path   : output\n-h, --help               : help\n--ora-reference          : ref\n";
        let bad = tree(&[
            ("help", help),
            (
                "scripts/b.sh",
                "\"$ORAD\" \\\n  --ora-reference \"$R\" \\\n  --output-directory \"$O\" \"$F\"\n",
            ),
        ]);
        let hf = bad.0.join("help");
        let hf = hf.to_str().unwrap();
        let (code, out) = run(
            &["--orad", "--help-file", hf, "--orad-version", "9"],
            &bad.0,
        );
        assert_eq!(code, 1, "{out}");
        assert!(
            out.contains("OK: --ora-reference      scripts/b.sh:1"),
            "{out}"
        );
        assert!(
            out.contains("FAIL: scripts/b.sh:1 passes --output-directory, which orad 9"),
            "{out}"
        );

        let good = tree(&[
            ("help", help),
            (
                "scripts/b.sh",
                "time \"${ORAD}\" --ora-reference=\"$R\" -P \"$O\" --path \"$O\" \"$F\"\n",
            ),
        ]);
        let hf = good.0.join("help");
        let (code, out) = run(&["--orad", "--help-file", hf.to_str().unwrap()], &good.0);
        assert_eq!(code, 0, "{out}");
        assert_eq!(out.lines().filter(|l| l.starts_with("OK: ")).count(), 3);
        assert_eq!(run(&["--orad-calls"], &good.0).0, 0);

        fs::write(&hf, "").unwrap();
        assert_eq!(
            run(&["--orad", "--help-file", hf.to_str().unwrap()], &good.0).0,
            1
        );

        let none = tree(&[("scripts/a.sh", "echo hi\n")]);
        let (code, out) = run(&["--orad-calls"], &none.0);
        assert_eq!(code, 2);
        assert!(out.contains("no \"$ORAD\" call"), "{out}");
    }

    #[test]
    fn usage_errors() {
        let d = tree(&[]);
        assert_eq!(run(&[], &d.0).0, 2);
        assert_eq!(run(&["--list", "--static"], &d.0).0, 2);
        assert_eq!(run(&["--orad"], &d.0).0, 2);
        assert_eq!(run(&["--bogus"], &d.0).0, 2);
    }

    #[test]
    fn groovy_strings_become_shell() {
        assert_eq!(
            groovy_to_shell("bgzip ${x} $y.z \\$HOME \\n", true),
            "bgzip GROOVY GROOVY $HOME \\n"
        );
        assert_eq!(groovy_to_shell("a ${x}", false), "a ${x}");
        assert_eq!(
            triple_quoted("x \"\"\"a\"\"\" '''b''' \"\"\"open"),
            [(true, 5, 6), (false, 13, 14)]
        );
    }
}
