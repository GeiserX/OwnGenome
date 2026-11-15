//! `og dev licences`: the notices of the third-party code in the og binary.
//!
//! ```text
//! og dev licences notices [--manifest-path FILE] [--package NAME]
//!     (--target TRIPLE | --target-no-default TRIPLE)... [--extra NAME=FILE]...
//!     [--out FILE]
//! ```
//!
//! The og binaries link about seventy crates statically, and the Linux pair
//! also links htslib and htscodecs (C, built by `hts-sys`). Their licences
//! (MIT, BSD-3-Clause, Unicode-3.0 and others) ask that a binary carry their
//! notices, so every release ships the file this verb writes
//! (`rust/og.sh build`, `.github/workflows/release.yml`).
//!
//! For each target it asks `cargo metadata --filter-platform TRIPLE` for the
//! resolved graph (`--target-no-default` adds `--no-default-features`, as the
//! macOS build is made) and walks it from the package (default `og`) along
//! normal edges only: build and dev dependencies, and procedural macros with
//! everything below them, run at compile time and leave no code in the
//! binary. Each crate reached outside the workspace must have a licence
//! expression the allow list ([`ALLOWED`]) satisfies, and a licence text:
//! the LICENSE, LICENCE, COPYING, NOTICE, UNLICENSE or COPYRIGHT files at its
//! root and its `license-file`. A crate that builds C (it sets `links`) also
//! contributes every such file below its root, so the vendored htslib and
//! htscodecs texts come from the crate itself and cannot go stale on a bump.
//! A crate that ships no text, under a licence that MIT satisfies, gets the
//! MIT text with its authors; any other crate without a text is an error.
//!
//! The binary also holds code that is no cargo package: the Rust standard
//! library and, in a musl build, the musl C library. `--extra NAME=FILE`
//! adds such a part with its licence file (an `.html` file, as the Rust
//! toolchain ships its notices, is turned into plain text); `rust/og.sh`
//! passes both from the builder image. A licence file that cannot be read,
//! a crate's or an extra one, is an error.
//!
//! The file lists the extra parts, then every crate with its version,
//! licence and, when it is not in every binary, the targets it is in; then
//! each distinct licence text
//! once, headed by the crates and files it covers. Errors name every crate
//! at fault and exit 1; usage errors exit 2.

use pycompat::json::{self, Value};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

/// The licences a linked crate may carry. Each is permissive and
/// compatible with GPL-3.0-or-later; an expression passes when these alone
/// satisfy it (an `OR` needs one side, an `AND` both).
pub const ALLOWED: &[&str] = &[
    "0BSD",
    "Apache-2.0",
    "Apache-2.0 WITH LLVM-exception",
    "BSD-3-Clause",
    "ISC",
    "MIT",
    "Unicode-3.0",
    "Unlicense",
    "Zlib",
];

const USAGE: &str = "usage: og dev licences notices [--manifest-path FILE] [--package NAME] (--target TRIPLE | --target-no-default TRIPLE)... [--extra NAME=FILE]... [--out FILE]";

/// The MIT text, for a crate under MIT that ships no licence file.
const MIT_TEXT: &str =
    "Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the \"Software\"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
";

/// `og dev licences ARGS`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    match args.first().map(String::as_str) {
        Some("notices") => notices_cli(&args[1..], out, err),
        _ => {
            let _ = writeln!(err, "{USAGE}");
            2
        }
    }
}

/// One cargo build of the binary: a target triple and its feature set.
struct Build {
    triple: String,
    default_features: bool,
}

struct Opts {
    manifest: Option<String>,
    package: String,
    builds: Vec<Build>,
    /// Parts that are no cargo package: (name, licence file).
    extras: Vec<(String, String)>,
    out: Option<String>,
}

fn parse(args: &[String]) -> Result<Opts, String> {
    let mut o = Opts {
        manifest: None,
        package: "og".into(),
        builds: Vec::new(),
        extras: Vec::new(),
        out: None,
    };
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let mut value = || {
            it.next()
                .cloned()
                .ok_or_else(|| format!("{a} needs a value"))
        };
        match a.as_str() {
            "--manifest-path" => o.manifest = Some(value()?),
            "--package" => o.package = value()?,
            "--out" => o.out = Some(value()?),
            "--extra" => {
                let v = value()?;
                match v.split_once('=') {
                    Some((name, file)) if !name.is_empty() && !file.is_empty() => {
                        o.extras.push((name.into(), file.into()));
                    }
                    _ => return Err(format!("--extra wants NAME=FILE, got {v:?}")),
                }
            }
            "--target" | "--target-no-default" => o.builds.push(Build {
                default_features: a == "--target",
                triple: value()?,
            }),
            _ => return Err(format!("unknown argument {a:?}")),
        }
    }
    if o.builds.is_empty() {
        return Err("name at least one --target or --target-no-default".into());
    }
    Ok(o)
}

fn notices_cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let opts = match parse(args) {
        Ok(o) => o,
        Err(e) => {
            let _ = writeln!(err, "og dev licences notices: {e}\n{USAGE}");
            return 2;
        }
    };
    match notices(&opts) {
        Ok(text) => {
            let written = match &opts.out {
                Some(path) => std::fs::write(path, &text).map_err(|e| format!("{path}: {e}")),
                None => out.write_all(text.as_bytes()).map_err(|e| e.to_string()),
            };
            match written {
                Ok(()) => 0,
                Err(e) => {
                    let _ = writeln!(err, "og dev licences notices: {e}");
                    1
                }
            }
        }
        Err(problems) => {
            for p in problems {
                let _ = writeln!(err, "og dev licences notices: {p}");
            }
            1
        }
    }
}

/// A third-party package the binary links, gathered over every build.
struct Crate {
    name: String,
    version: String,
    license: Option<String>,
    license_file: Option<String>,
    authors: Vec<String>,
    dir: PathBuf,
    builds_c: bool,
    triples: BTreeSet<String>,
}

fn notices(opts: &Opts) -> Result<String, Vec<String>> {
    let mut crates: BTreeMap<(String, String), Crate> = BTreeMap::new();
    for b in &opts.builds {
        let meta = cargo_metadata(opts, b).map_err(|e| vec![e])?;
        for c in linked(&meta, &opts.package).map_err(|e| vec![e])? {
            crates
                .entry((c.name.clone(), c.version.clone()))
                .or_insert(c)
                .triples
                .insert(b.triple.clone());
        }
    }
    let mut problems = Vec::new();
    // Each distinct text once, in the order its first part comes.
    let mut texts: Vec<(String, Vec<String>)> = Vec::new();
    let mut index: HashMap<String, usize> = HashMap::new();
    let mut add = |who: String, text: String| match index.get(&text) {
        Some(&i) => texts[i].1.push(who),
        None => {
            index.insert(text.clone(), texts.len());
            texts.push((text, vec![who]));
        }
    };
    for (name, file) in &opts.extras {
        match std::fs::read(file) {
            Ok(bytes) => {
                let raw = String::from_utf8_lossy(&bytes);
                let text = if file.ends_with(".html") {
                    html_text(&raw)
                } else {
                    raw.into_owned()
                };
                let base = Path::new(file)
                    .file_name()
                    .map_or(file.as_str(), |n| n.to_str().unwrap_or(file.as_str()));
                add(format!("{name}: {base}"), text);
            }
            Err(e) => problems.push(format!("{name}: cannot read {file}: {e}")),
        }
    }
    for c in crates.values() {
        let id = format!("{} {}", c.name, c.version);
        let Some(expr) = &c.license else {
            problems.push(format!("{id} declares no licence expression"));
            continue;
        };
        match satisfied(expr, &|l| ALLOWED.contains(&l)) {
            Ok(true) => {}
            Ok(false) => {
                problems.push(format!(
                    "{id} is licensed {expr:?}, which the allow list ({}) does not satisfy",
                    ALLOWED.join(", ")
                ));
                continue;
            }
            Err(e) => {
                problems.push(format!("{id}: cannot read its licence {expr:?}: {e}"));
                continue;
            }
        }
        let mut found = match licence_files(c) {
            Ok(f) => f,
            Err(e) => {
                problems.push(format!("{id}: {e}"));
                continue;
            }
        };
        if found.is_empty() {
            if satisfied(expr, &|l| l == "MIT") == Ok(true) {
                let who = if c.authors.is_empty() {
                    format!("the {} authors", c.name)
                } else {
                    c.authors.join(", ")
                };
                found.push((
                    "MIT, no licence file in the crate".into(),
                    format!("Copyright (c) {who}\n\n{MIT_TEXT}"),
                ));
            } else {
                problems.push(format!(
                    "{id} ships no licence file, and {expr:?} cannot be met with the MIT text"
                ));
                continue;
            }
        }
        for (label, text) in found {
            add(format!("{id}: {label}"), text);
        }
    }
    if !problems.is_empty() {
        return Err(problems);
    }
    Ok(render(opts, &crates, &texts))
}

fn render(
    opts: &Opts,
    crates: &BTreeMap<(String, String), Crate>,
    texts: &[(String, Vec<String>)],
) -> String {
    let all: BTreeSet<String> = opts.builds.iter().map(|b| b.triple.clone()).collect();
    let mut s = String::new();
    s.push_str("Third-party licences of the og binaries\n");
    s.push_str("=======================================\n\n");
    s.push_str("OwnGenome is GPL-3.0-or-later (LICENSE in its source). The og binaries\n");
    s.push_str("also hold the code of the parts and crates below, linked statically. A\n");
    s.push_str("crate that builds C carries that library's licence in its own folder:\n");
    s.push_str("hts-sys gives htslib/LICENSE and htslib/htscodecs/LICENSE.md, the licences\n");
    s.push_str("of htslib and htscodecs. Each licence text comes once below, after the\n");
    s.push_str("parts, crates and files it covers. `og dev licences notices` writes this\n");
    s.push_str("file from cargo metadata and the files each crate ships.\n\n");
    s.push_str(&format!(
        "Targets: {}\n\n",
        all.iter().cloned().collect::<Vec<_>>().join(", ")
    ));
    if !opts.extras.is_empty() {
        s.push_str("Parts that are no crate\n-----------------------\n");
        for (name, _) in &opts.extras {
            s.push_str(&format!("{name}\n"));
        }
        s.push('\n');
    }
    let title = format!("Crates ({})", crates.len());
    s.push_str(&format!("{title}\n{}\n", "-".repeat(title.len())));
    for c in crates.values() {
        s.push_str(&format!(
            "{} {}  {}",
            c.name,
            c.version,
            c.license.as_deref().unwrap_or("")
        ));
        if c.triples != all {
            let only: Vec<&str> = c.triples.iter().map(String::as_str).collect();
            s.push_str(&format!("  (only {})", only.join(", ")));
        }
        s.push('\n');
    }
    let title = format!("Licence texts ({})", texts.len());
    s.push_str(&format!("\n{title}\n{}\n", "-".repeat(title.len())));
    for (i, (text, who)) in texts.iter().enumerate() {
        s.push_str(&format!("\n=== {} ===\n", i + 1));
        for w in who {
            s.push_str(&format!("For {w}\n"));
        }
        s.push('\n');
        s.push_str(text.trim_end());
        s.push('\n');
    }
    s
}

fn cargo_metadata(opts: &Opts, b: &Build) -> Result<Value, String> {
    let cargo = std::env::var_os("CARGO").unwrap_or_else(|| "cargo".into());
    let mut cmd = Command::new(cargo);
    cmd.args([
        "metadata",
        "--format-version",
        "1",
        "--locked",
        "--filter-platform",
    ])
    .arg(&b.triple);
    if !b.default_features {
        cmd.arg("--no-default-features");
    }
    if let Some(m) = &opts.manifest {
        cmd.arg("--manifest-path").arg(m);
    }
    let run = cmd
        .output()
        .map_err(|e| format!("cargo metadata did not start: {e}"))?;
    if !run.status.success() {
        return Err(format!(
            "cargo metadata for {} failed: {}",
            b.triple,
            String::from_utf8_lossy(&run.stderr).trim()
        ));
    }
    let text =
        String::from_utf8(run.stdout).map_err(|_| "cargo metadata wrote non-UTF-8".to_string())?;
    json::loads(&text)
        .map_err(|e| format!("cargo metadata wrote unreadable JSON: {}", e.exception()))
}

fn field<'a>(v: &'a Value, key: &str) -> Option<&'a Value> {
    match v {
        Value::Object(m) => m.get(key),
        _ => None,
    }
}

fn text<'a>(v: &'a Value, key: &str) -> Option<&'a str> {
    match field(v, key) {
        Some(Value::Str(s)) => Some(s),
        _ => None,
    }
}

fn list<'a>(v: &'a Value, key: &str) -> &'a [Value] {
    match field(v, key) {
        Some(Value::Array(a)) => a,
        _ => &[],
    }
}

/// The third-party packages `package` links on this build's target.
fn linked(meta: &Value, package: &str) -> Result<Vec<Crate>, String> {
    let packages: HashMap<&str, &Value> = list(meta, "packages")
        .iter()
        .filter_map(|p| Some((text(p, "id")?, p)))
        .collect();
    let members: HashSet<&str> = list(meta, "workspace_members")
        .iter()
        .filter_map(|v| match v {
            Value::Str(s) => Some(s.as_str()),
            _ => None,
        })
        .collect();
    let resolve = field(meta, "resolve").ok_or("cargo metadata gave no resolve graph")?;
    let nodes: HashMap<&str, &Value> = list(resolve, "nodes")
        .iter()
        .filter_map(|n| Some((text(n, "id")?, n)))
        .collect();
    let root = members
        .iter()
        .copied()
        .find(|id| packages.get(id).and_then(|p| text(p, "name")) == Some(package))
        .ok_or_else(|| format!("no workspace package is named {package}"))?;
    let mut seen = HashSet::from([root]);
    let mut stack = vec![root];
    while let Some(id) = stack.pop() {
        let Some(node) = nodes.get(id) else { continue };
        for dep in list(node, "deps") {
            let Some(pkg) = text(dep, "pkg") else {
                continue;
            };
            // kind null is a normal dependency; "build" and "dev" are not linked.
            let normal = list(dep, "dep_kinds")
                .iter()
                .any(|k| matches!(field(k, "kind"), None | Some(Value::Null)));
            let proc_macro = packages.get(pkg).is_some_and(|p| {
                list(p, "targets").iter().any(|t| {
                    list(t, "kind")
                        .iter()
                        .any(|k| matches!(k, Value::Str(s) if s == "proc-macro"))
                })
            });
            if normal && !proc_macro && seen.insert(pkg) {
                stack.push(pkg);
            }
        }
    }
    let mut out = Vec::new();
    for id in seen {
        if members.contains(id) {
            continue;
        }
        let p = packages.get(id).ok_or_else(|| {
            format!("cargo metadata lists {id} in its graph but not as a package")
        })?;
        let manifest =
            text(p, "manifest_path").ok_or_else(|| format!("{id} has no manifest_path"))?;
        out.push(Crate {
            name: text(p, "name").unwrap_or(id).into(),
            version: text(p, "version").unwrap_or("").into(),
            license: text(p, "license").map(str::to_owned),
            license_file: text(p, "license_file").map(str::to_owned),
            authors: list(p, "authors")
                .iter()
                .filter_map(|a| match a {
                    Value::Str(s) => Some(s.clone()),
                    _ => None,
                })
                .collect(),
            dir: Path::new(manifest)
                .parent()
                .unwrap_or(Path::new("."))
                .to_path_buf(),
            builds_c: text(p, "links").is_some(),
            triples: BTreeSet::new(),
        });
    }
    Ok(out)
}

fn is_licence_name(name: &str) -> bool {
    let n = name.to_ascii_lowercase();
    [
        "license",
        "licence",
        "copying",
        "notice",
        "unlicense",
        "copyright",
    ]
    .iter()
    .any(|p| n.starts_with(p))
}

/// (label, text) of every licence file of a crate, label = path in the
/// crate. A file found but not readable is an error, never a silent gap.
fn licence_files(c: &Crate) -> Result<Vec<(String, String)>, String> {
    let mut paths: BTreeSet<PathBuf> = BTreeSet::new();
    collect(&c.dir, c.builds_c, &mut paths);
    if let Some(f) = &c.license_file {
        let p = c.dir.join(f);
        if p.is_file() {
            paths.insert(p);
        }
    }
    paths
        .into_iter()
        .map(|p| {
            let label = p.strip_prefix(&c.dir).unwrap_or(&p).display().to_string();
            let bytes = std::fs::read(&p).map_err(|e| format!("cannot read {label}: {e}"))?;
            Ok((label, String::from_utf8_lossy(&bytes).into_owned()))
        })
        .collect()
}

/// The text of an HTML notice file: tags dropped, the five escapes HTML
/// writes undone, runs of blank lines folded to one.
fn html_text(html: &str) -> String {
    let mut bare = String::with_capacity(html.len());
    let mut in_tag = false;
    for ch in html.chars() {
        match ch {
            '<' => in_tag = true,
            '>' if in_tag => in_tag = false,
            c if !in_tag => bare.push(c),
            _ => {}
        }
    }
    let bare = bare
        .replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", "\"")
        .replace("&#39;", "'")
        .replace("&#x27;", "'")
        .replace("&amp;", "&");
    let mut out = String::with_capacity(bare.len());
    let mut blank = true;
    for line in bare.lines().map(str::trim_end) {
        if line.is_empty() {
            if !blank {
                out.push('\n');
            }
            blank = true;
        } else {
            out.push_str(line);
            out.push('\n');
            blank = false;
        }
    }
    out.truncate(out.trim_end().len());
    out.push('\n');
    out
}

/// Licence files in `dir`, and below it when `deep` (a crate that builds C).
fn collect(dir: &Path, deep: bool, out: &mut BTreeSet<PathBuf>) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    for e in entries.flatten() {
        let Ok(kind) = e.file_type() else { continue };
        let path = e.path();
        if kind.is_file() && is_licence_name(&e.file_name().to_string_lossy()) {
            out.insert(path);
        } else if kind.is_dir() && deep {
            collect(&path, deep, out);
        }
    }
}

/// Whether an SPDX expression is met by the licences `ok` accepts: `OR`
/// (or the old `/`) needs one side, `AND` both, `X WITH Y` is one licence.
pub fn satisfied(expr: &str, ok: &dyn Fn(&str) -> bool) -> Result<bool, String> {
    let spaced = expr
        .replace('(', " ( ")
        .replace(')', " ) ")
        .replace('/', " OR ");
    let tokens: Vec<&str> = spaced.split_whitespace().collect();
    let mut p = Parser {
        tokens: &tokens,
        at: 0,
        ok,
    };
    let v = p.or()?;
    if p.at != tokens.len() {
        return Err(format!("unexpected {:?}", tokens[p.at]));
    }
    Ok(v)
}

struct Parser<'a> {
    tokens: &'a [&'a str],
    at: usize,
    ok: &'a dyn Fn(&str) -> bool,
}

impl<'a> Parser<'a> {
    fn peek(&self) -> Option<&'a str> {
        self.tokens.get(self.at).copied()
    }
    fn or(&mut self) -> Result<bool, String> {
        let mut v = self.and()?;
        while self.peek().is_some_and(|t| t.eq_ignore_ascii_case("OR")) {
            self.at += 1;
            v |= self.and()?;
        }
        Ok(v)
    }
    fn and(&mut self) -> Result<bool, String> {
        let mut v = self.one()?;
        while self.peek().is_some_and(|t| t.eq_ignore_ascii_case("AND")) {
            self.at += 1;
            v &= self.one()?;
        }
        Ok(v)
    }
    fn one(&mut self) -> Result<bool, String> {
        let t = self.peek().ok_or("the expression ends early")?;
        self.at += 1;
        if t == "(" {
            let v = self.or()?;
            if self.peek() != Some(")") {
                return Err("a ( is not closed".into());
            }
            self.at += 1;
            return Ok(v);
        }
        if [")", "AND", "OR", "WITH"]
            .iter()
            .any(|k| t.eq_ignore_ascii_case(k))
        {
            return Err(format!("unexpected {t:?}"));
        }
        if self.peek().is_some_and(|w| w.eq_ignore_ascii_case("WITH")) {
            let exception = self
                .tokens
                .get(self.at + 1)
                .ok_or("WITH needs an exception")?;
            self.at += 2;
            return Ok((self.ok)(&format!("{t} WITH {exception}")));
        }
        Ok((self.ok)(t))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn allowed(expr: &str) -> Result<bool, String> {
        satisfied(expr, &|l| ALLOWED.contains(&l))
    }

    #[test]
    fn expressions_the_locked_crates_use_pass() {
        for e in [
            "MIT",
            "MIT OR Apache-2.0",
            "MIT/Apache-2.0",
            "(MIT OR Apache-2.0) AND Unicode-3.0",
            "Apache-2.0 WITH LLVM-exception OR Apache-2.0 OR MIT",
            "MIT OR Apache-2.0 OR LGPL-2.1-or-later",
            "0BSD OR MIT OR Apache-2.0",
            "Unlicense OR MIT",
        ] {
            assert_eq!(allowed(e), Ok(true), "{e}");
        }
    }

    #[test]
    fn a_licence_outside_the_list_fails() {
        for e in [
            "SSPL-1.0",
            "MIT AND SSPL-1.0",
            "LGPL-2.1-or-later",
            "GPL-2.0 WITH Classpath-exception-2.0",
            "(MIT OR Apache-2.0) AND MPL-2.0",
        ] {
            assert_eq!(allowed(e), Ok(false), "{e}");
        }
    }

    #[test]
    fn html_notices_become_text() {
        let html =
            "<h1>Notices</h1>\n<p>A &amp; B &lt;x&gt; &quot;q&quot;</p>\n\n\n<pre>line\n</pre>\n";
        assert_eq!(html_text(html), "Notices\nA & B <x> \"q\"\n\nline\n");
    }

    #[test]
    fn an_unreadable_expression_is_an_error() {
        for e in ["", "MIT OR", "(MIT", "MIT)", "AND MIT", "MIT WITH"] {
            assert!(allowed(e).is_err(), "{e:?} should not parse");
        }
    }
}
