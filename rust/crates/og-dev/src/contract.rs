//! `og dev contract`: the 1.0 output contract, `assets/output_contract.tsv`
//! (`docs/output-contract.md`), against the code and against a run. It
//! replaces `scripts/ci/check-output-contract.py`.
//!
//! - `--static [--repo DIR] [--contract FILE]`: with no run, every param of
//!   `nextflow_schema.json`, every `--tools` name of `main.nf`
//!   (`known_tools`) and every `run-all.sh` switch (`SWITCHES`) has a row
//!   and every row still names one; every file row of a process names a
//!   file that process still writes and publishes, in the folder of a
//!   `publishDir` that keeps it; every file a module publishes has a row; a
//!   `nextflow` row (not internal) whose step script writes the same path is
//!   refused (its legs are both); every default of `nextflow_schema.json`
//!   equals `nextflow.config`'s.
//! - `--outputs DIR --leg bash|nextflow --sample S [--observe]`: every
//!   stable row whose `e2e` column names the leg is there under DIR with its
//!   required columns (tsv, csv), INFO/FORMAT/FILTER fields (vcf) or
//!   top-level keys (json). `--observe` also lists every file under
//!   `DIR/S` with what it has.
//!
//! Exit 0 when nothing is broken, 1 when something is (or the contract is
//! malformed: `ERROR: …`), 2 on a usage error. `--repo` defaults to the
//! current folder.
//!
//! Where it differs from the Python on purpose:
//! - a schema default and a `nextflow.config` value of different JSON types
//!   differ: `true` is not `1` and `false` is not `0` (Python's `True == 1`
//!   let that drift pass); `1` and `1.0` stay equal.
//! - a file the run check cannot read (a truncated gzip) is a `FAIL` naming
//!   it, not a traceback.
//! - Groovy is read with comments blanked and strings respected, so a
//!   commented-out `path(...)` is not an output and a `}` inside a string
//!   does not end a `saveAs` closure.

use crate::py_os_error;
use og_io::{Lines, Utf8};
use og_nf::module::Process;
use pycompat::json::{self, Value};
use pycompat::mtime::{fnmatch, glob};
use pycompat::value::py_repr;
use regex::Regex;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs;
use std::io::{self, BufRead, Read, Write};
use std::path::Path;
use std::sync::LazyLock;

const USAGE: &str = "usage: og dev contract [--contract CONTRACT] [--repo REPO] [--static] [--outputs OUTPUTS]\n                       [--leg {bash,nextflow}] [--sample SAMPLE] [--observe]\n";

const COLUMNS: [&str; 9] = [
    "kind",
    "name",
    "step",
    "process",
    "legs",
    "format",
    "stability",
    "e2e",
    "required",
];
const LEGS: [&str; 3] = ["bash", "nextflow", "both"];
const E2E: [&str; 4] = ["-", "bash", "nextflow", "both"];
const FORMATS: [&str; 18] = [
    "tsv", "csv", "vcf", "json", "text", "html", "bam", "cram", "index", "bed", "dir", "zip",
    "svg", "png", "md", "yml", "fasta", "bin",
];

/// What a contract row promises.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Kind {
    File,
    Param,
    Tool,
    Switch,
}

impl Kind {
    fn parse(s: &str) -> Option<Kind> {
        Some(match s {
            "file" => Kind::File,
            "param" => Kind::Param,
            "tool" => Kind::Tool,
            "switch" => Kind::Switch,
            _ => return None,
        })
    }

    fn word(self) -> &'static str {
        match self {
            Kind::File => "file",
            Kind::Param => "param",
            Kind::Tool => "tool",
            Kind::Switch => "switch",
        }
    }
}

/// How long a row's promise holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stability {
    Stable,
    Experimental,
    Internal,
}

impl Stability {
    fn parse(s: &str) -> Option<Stability> {
        Some(match s {
            "stable" => Stability::Stable,
            "experimental" => Stability::Experimental,
            "internal" => Stability::Internal,
            _ => return None,
        })
    }
}

/// One row of the contract. `legs`, `format` and `e2e` are checked against
/// their word lists on file rows; `step` and `process` are `-` when unset.
#[derive(Debug, Clone)]
pub struct Row {
    pub kind: Kind,
    pub name: String,
    pub step: String,
    pub process: String,
    pub legs: String,
    pub format: String,
    pub stability: Stability,
    pub e2e: String,
    pub required: String,
    pub line: usize,
}

/// The rows of the contract at `path`, or why it is malformed.
pub fn read_contract(path: &str) -> Result<Vec<Row>, String> {
    let file = fs::File::open(path).map_err(|e| py_os_error(&e, path))?;
    let mut rows = Vec::new();
    let mut header = false;
    for (i, line) in Lines::new(io::BufReader::new(file), Utf8::Strict).enumerate() {
        let n = i + 1;
        let line = line.map_err(|e| format!("{path}:{n}: {e}"))?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        if !header {
            if cells != COLUMNS {
                return Err(format!(
                    "{path}:{n}: header must be {}, got {}",
                    COLUMNS.join(" "),
                    cells.join(" ")
                ));
            }
            header = true;
            continue;
        }
        if cells.len() != COLUMNS.len() {
            return Err(format!(
                "{path}:{n}: {} cells, want {}",
                cells.len(),
                COLUMNS.len()
            ));
        }
        let mut bad = Vec::new();
        let kind = Kind::parse(cells[0]);
        let stability = Stability::parse(cells[6]);
        if kind.is_none() {
            bad.push(format!("kind {}", cells[0]));
        }
        if stability.is_none() {
            bad.push(format!("stability {}", cells[6]));
        }
        if kind == Some(Kind::File) {
            if !LEGS.contains(&cells[4]) {
                bad.push(format!("legs {}", cells[4]));
            }
            if !FORMATS.contains(&cells[5]) {
                bad.push(format!("format {}", cells[5]));
            }
            if !E2E.contains(&cells[7]) {
                bad.push(format!("e2e {}", cells[7]));
            }
            if cells[4] == "bash" && cells[3] != "-" {
                bad.push("a bash-only row names a process".to_owned());
            }
        }
        let (Some(kind), Some(stability)) = (kind, stability) else {
            return Err(format!("{path}:{n}: bad {}", bad.join(", ")));
        };
        if !bad.is_empty() {
            return Err(format!("{path}:{n}: bad {}", bad.join(", ")));
        }
        let s = |k: usize| cells[k].to_owned();
        rows.push(Row {
            kind,
            name: s(1),
            step: s(2),
            process: s(3),
            legs: s(4),
            format: s(5),
            stability,
            e2e: s(7),
            required: s(8),
            line: n,
        });
    }
    if !header {
        return Err(format!("{path}: no header line"));
    }
    let mut seen: HashMap<(Kind, &str, &str), usize> = HashMap::new();
    for r in &rows {
        let legs = if r.kind == Kind::File {
            r.legs.as_str()
        } else {
            ""
        };
        if let Some(first) = seen.insert((r.kind, &r.name, legs), r.line) {
            return Err(format!(
                "{path}:{}: {} {} is listed twice (line {first})",
                r.line,
                r.kind.word(),
                r.name
            ));
        }
    }
    Ok(rows)
}

// ---------------------------------------------------------------- the code side

fn re(s: &str) -> Regex {
    Regex::new(s).expect("a valid built-in pattern")
}

static BRACES: LazyLock<Regex> = LazyLock::new(|| re(r"\{([^{}]*)\}"));
static FAMILY: LazyLock<Regex> = LazyLock::new(|| re(r"\$\{params\.family_name[^}]*\}"));
static SAVE_AS: LazyLock<Regex> = LazyLock::new(|| {
    re(
        r"(?s)\A\s*(\w+)\s*->\s*(\w+)\s*(?:\.(endsWith|startsWith)\(\s*'([^']*)'\s*\)|==\s*'([^']*)')\s*\?\s*([^?:]+?)\s*:\s*([^?:]+?)\s*\z",
    )
});
static KNOWN_TOOLS: LazyLock<Regex> = LazyLock::new(|| re(r"(?s)def known_tools = \[(.*?)\]"));
static TOOL: LazyLock<Regex> = LazyLock::new(|| re(r"'([a-z0-9_]+)'"));
static SWITCHES: LazyLock<Regex> = LazyLock::new(|| re(r"(?ms)^SWITCHES=\((.*?)\)"));
static VAR: LazyLock<Regex> = LazyLock::new(|| re(r"\$\{(\w+)\}|\$(\w+)"));
static ASSIGN: LazyLock<Regex> = LazyLock::new(|| {
    re(
        r#"(?m)^\s*(?:local\s+|export\s+|readonly\s+)?([A-Za-z_]\w*)=("[^"\n]*"|'[^'\n]*'|[^\s;#()]+)"#,
    )
});
static VCF_FIELD: LazyLock<Regex> = LazyLock::new(|| re(r"^##(INFO|FORMAT|FILTER)=<ID=([^,>]+)"));

/// One level of `{a,b}` alternatives expanded, as Nextflow globs do.
pub fn braces(pattern: &str) -> Vec<String> {
    let Some(c) = BRACES.captures(pattern) else {
        return vec![pattern.to_owned()];
    };
    let m = c.get(0).expect("match");
    c[1].split(',')
        .flat_map(|alt| {
            braces(&format!(
                "{}{alt}{}",
                &pattern[..m.start()],
                &pattern[m.end()..]
            ))
        })
        .collect()
}

/// A Nextflow path with `${meta.id}` and the family name filled in.
pub fn concrete(pattern: &str) -> String {
    let p = pattern
        .replace("${meta.id}", "SAMPLE")
        .replace("{sample}", "SAMPLE");
    let p = FAMILY
        .replace_all(&p, "FAMILY")
        .replace("{family}", "FAMILY");
    let p = p
        .replace("${params.outdir}/", "")
        .replace("${params.outdir}", "");
    p.trim_matches('/').to_owned()
}

/// What a `saveAs` closure does to a task file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SaveAs {
    /// Published under its own name (also: no closure).
    Keep,
    /// The closure returns `null`: not published.
    Drop,
    /// Published under another name, or by a closure this reader does not
    /// follow; matched by basename only.
    Rename,
}

/// `x -> x.endsWith('…') ? A : B` (or `startsWith`, or `x == '…'`)
/// evaluated for `name`; any other closure is [`SaveAs::Rename`].
pub fn save_as(closure: Option<&str>, name: &str) -> SaveAs {
    let Some(closure) = closure else {
        return SaveAs::Keep;
    };
    let Some(c) = SAVE_AS.captures(closure) else {
        return SaveAs::Rename;
    };
    let var = &c[1];
    if &c[2] != var {
        return SaveAs::Rename;
    }
    let cond = match c.get(3).map(|m| m.as_str()) {
        Some("endsWith") => name.ends_with(&c[4]),
        Some(_) => name.starts_with(&c[4]),
        None => name == &c[5],
    };
    let expr = if cond { &c[6] } else { &c[7] }.trim();
    if expr == "null" {
        SaveAs::Drop
    } else if expr == var {
        SaveAs::Keep
    } else {
        SaveAs::Rename
    }
}

/// A module process with its file, as the static check uses it.
struct Proc {
    file: String,
    /// `(concrete dir, expanded patterns, saveAs body)` of each publishDir.
    publish: Vec<(String, Option<Vec<String>>, Option<String>)>,
    outputs: Vec<String>,
}

fn read_text(repo: &str, rel: &str) -> Result<String, String> {
    let path = Path::new(repo).join(rel);
    let mut s = String::new();
    fs::File::open(&path)
        .and_then(|mut f| f.read_to_string(&mut s))
        .map_err(|e| py_os_error(&e, &path.to_string_lossy()))?;
    // Python's text mode: \r\n and a lone \r read as \n.
    Ok(if s.contains('\r') {
        s.replace("\r\n", "\n").replace('\r', "\n")
    } else {
        s
    })
}

fn module_processes(repo: &str) -> Result<BTreeMap<String, Proc>, String> {
    let pattern = Path::new(repo).join("modules/local/*/main.nf");
    let mut files = glob(pattern.as_os_str(), false);
    files.sort();
    let mut procs = BTreeMap::new();
    for f in files {
        let rel = f
            .strip_prefix(repo)
            .unwrap_or(&f)
            .to_string_lossy()
            .trim_start_matches('/')
            .to_owned();
        let src = read_text(repo, &rel)?;
        for Process {
            name,
            publish,
            outputs,
        } in og_nf::module::scan(&src)
        {
            let publish = publish
                .into_iter()
                .map(|p| {
                    (
                        concrete(&p.dir),
                        p.pattern.as_deref().map(braces),
                        p.save_as,
                    )
                })
                .collect();
            procs.insert(
                name,
                Proc {
                    file: rel.clone(),
                    publish,
                    outputs,
                },
            );
        }
    }
    Ok(procs)
}

/// `(contracted-path glob, basename only)` of each file a process
/// publishes.
fn published(p: &Proc) -> Vec<(String, bool)> {
    let mut items = Vec::new();
    for (d, pats, closure) in &p.publish {
        for g in &p.outputs {
            let c = concrete(g);
            if let Some(pats) = pats {
                let x = c.replace('*', "X");
                if !pats.iter().any(|pat| fnmatch(&x, pat) || fnmatch(&c, pat)) {
                    continue;
                }
            }
            let action = save_as(closure.as_deref(), &c);
            if action == SaveAs::Drop {
                continue;
            }
            let path = if d.is_empty() { c } else { format!("{d}/{c}") };
            items.push((path, action == SaveAs::Rename));
        }
    }
    items
}

fn basename(p: &str) -> &str {
    p.rsplit('/').next().unwrap_or(p)
}

/// Equal as the contract means it: by JSON type (a bool never equals a
/// number), an int equal to a float of the same value.
fn same(a: &Value, b: &Value) -> bool {
    match (a, b) {
        (Value::Int(x), Value::Float(y)) | (Value::Float(y), Value::Int(x)) => *x as f64 == *y,
        _ => a == b,
    }
}

fn repr(v: &Value) -> String {
    py_repr(v).unwrap_or_else(|_| json::dumps(v, &json::Dump::default()).unwrap_or_default())
}

/// The lines of a step script with its shell variables filled in: SAMPLE
/// as `{sample}`, GENOME_DIR as `G`, every other variable as each value it
/// is assigned (at most 4 values, 64 expansions per line, depth 24).
fn bash_paths(src: &str) -> HashSet<String> {
    let mut assign: HashMap<&str, Vec<String>> = HashMap::new();
    for c in ASSIGN.captures_iter(src) {
        let v = c.get(1).expect("group").as_str();
        assign
            .entry(v)
            .or_default()
            .push(c[2].trim_matches(['"', '\'']).to_owned());
    }
    let nul = vec!["\0".to_owned()];
    let fixed = |name: &str| -> Option<&'static str> {
        match name {
            "SAMPLE" | "1" => Some("{sample}"),
            "GENOME_DIR" => Some("G"),
            _ => None,
        }
    };
    fn expand(
        s: &str,
        depth: usize,
        values: &dyn Fn(&str) -> Vec<String>,
        out: &mut Vec<String>,
    ) -> usize {
        let start = out.len();
        let Some(c) = VAR.captures(s) else {
            out.push(s.to_owned());
            return 1;
        };
        if depth > 24 {
            out.push(s.to_owned());
            return 1;
        }
        let m = c.get(0).expect("match");
        let name = c.get(1).or(c.get(2)).expect("a branch").as_str();
        for val in values(name).iter().take(4) {
            let next = format!("{}{val}{}", &s[..m.start()], &s[m.end()..]);
            expand(&next, depth + 1, values, out);
            if out.len() - start >= 64 {
                break;
            }
        }
        out.truncate(start + 64);
        out.len() - start
    }
    let values = |name: &str| -> Vec<String> {
        match fixed(name) {
            Some(v) => vec![v.to_owned()],
            None => assign.get(name).unwrap_or(&nul).clone(),
        }
    };
    let mut set = HashSet::new();
    let mut buf = Vec::new();
    for line in src.split('\n') {
        if pycompat::text::lstrip(line, None).starts_with('#') {
            continue;
        }
        buf.clear();
        expand(line, 0, &values, &mut buf);
        set.extend(buf.drain(..));
    }
    set
}

/// Python's `\w` for one character.
fn is_word(c: char) -> bool {
    c.is_alphanumeric() || c == '_'
}

type ScriptCache = HashMap<String, (Option<String>, HashSet<String>)>;

/// The step script (`scripts/<step>-*.sh`) that writes this per-sample path
/// under the same name, the whole path (`x.tsv.gz` is not `x.tsv`).
pub fn bash_writes(
    repo: &str,
    step: &str,
    name: &str,
    cache: &mut HashMap<String, (Option<String>, HashSet<String>)>,
) -> Result<Option<String>, String> {
    if !name.starts_with("{sample}/") || name.contains('*') || step == "-" {
        return Ok(None);
    }
    if !cache.contains_key(step) {
        let mut files = glob(
            Path::new(repo)
                .join("scripts")
                .join(format!("{step}-*.sh"))
                .as_os_str(),
            false,
        );
        files.sort();
        let mut src = String::new();
        let mut first = None;
        for f in &files {
            let rel = f
                .strip_prefix(repo)
                .unwrap_or(f)
                .to_string_lossy()
                .trim_start_matches('/')
                .to_owned();
            src.push_str(&read_text(repo, &rel)?);
            first.get_or_insert(rel);
        }
        let lines = if first.is_some() {
            bash_paths(&src)
        } else {
            HashSet::new()
        };
        cache.insert(step.to_owned(), (first, lines));
    }
    let (first, lines) = &cache[step];
    let needle = format!("G/{name}");
    let hit = lines.iter().any(|l| {
        l.match_indices(&needle).any(|(i, _)| {
            l[i + needle.len()..]
                .chars()
                .next()
                .is_none_or(|c| !(is_word(c) || c == '.' || c == '-'))
        })
    });
    Ok(if hit { first.clone() } else { None })
}

fn check_static(repo: &str, rows: &[Row]) -> Result<Vec<String>, String> {
    let schema_text = read_text(repo, "nextflow_schema.json")?;
    let schema = json::loads(&schema_text)
        .map_err(|e| format!("{repo}/nextflow_schema.json: not JSON: {e:?}"))?;
    let schema_params = og_nf::schema_params::read(&schema);
    let main_nf = read_text(repo, "main.nf")?;
    let run_all = read_text(repo, "scripts/run-all.sh")?;
    let config = read_text(repo, "nextflow.config")?;

    let tools: BTreeSet<String> = KNOWN_TOOLS
        .captures(&main_nf)
        .map(|c| TOOL.captures_iter(&c[1]).map(|t| t[1].to_owned()).collect())
        .unwrap_or_default();
    let switches: BTreeSet<String> = SWITCHES
        .captures(&run_all)
        .map(|c| c[1].split_whitespace().map(str::to_owned).collect())
        .unwrap_or_default();
    let params: BTreeSet<String> = schema_params.iter().map(|(n, _)| n.clone()).collect();

    let mut problems = Vec::new();
    for (kind, have, where_) in [
        (Kind::Param, &params, "nextflow_schema.json"),
        (Kind::Tool, &tools, "main.nf known_tools"),
        (Kind::Switch, &switches, "scripts/run-all.sh SWITCHES"),
    ] {
        let k = kind.word();
        if have.is_empty() {
            problems.push(format!(
                "found no {k}s in {where_}: the reader no longer matches the code"
            ));
            continue;
        }
        let listed: BTreeMap<&str, &Row> = rows
            .iter()
            .filter(|r| r.kind == kind)
            .map(|r| (r.name.as_str(), r))
            .collect();
        for (name, r) in &listed {
            if !have.contains(*name) {
                let what = if r.stability == Stability::Stable {
                    "BREAKING: the stable"
                } else {
                    "the"
                };
                problems.push(format!(
                    "{what} {k} {name} (contract line {}) is gone from {where_}",
                    r.line
                ));
            }
        }
        for name in have {
            if !listed.contains_key(name.as_str()) {
                problems.push(format!(
                    "{k} {name} of {where_} has no row in the contract: add one"
                ));
            }
        }
    }

    // The last default and the last config value of a name win, as a dict.
    let defaults: BTreeMap<String, Value> = schema_params
        .into_iter()
        .filter_map(|(n, d)| d.map(|d| (n, d)))
        .collect();
    let config: HashMap<String, Value> = og_nf::config::params(&config).into_iter().collect();
    if config.is_empty() {
        problems.push(
            "found no params in nextflow.config: the reader no longer matches the code".to_owned(),
        );
    }
    for (name, val) in &defaults {
        if let Some(c) = config.get(name)
            && !same(c, val)
        {
            problems.push(format!(
                "param {name}: nextflow_schema.json's default {} is not nextflow.config's {}",
                repr(val),
                repr(c)
            ));
        }
    }

    let procs = module_processes(repo)?;
    if procs.is_empty() {
        problems.push("found no processes under modules/local".to_owned());
    }
    let file_rows: Vec<&Row> = rows
        .iter()
        .filter(|r| r.kind == Kind::File && r.process != "-")
        .collect();
    for r in &file_rows {
        let Some(p) = procs.get(&r.process) else {
            let what = if r.stability == Stability::Stable {
                "BREAKING: "
            } else {
                ""
            };
            problems.push(format!(
                "{what}{} (contract line {}): process {} is not in modules/local",
                r.name, r.line, r.process
            ));
            continue;
        };
        let want = concrete(&r.name);
        let ok = if r.stability == Stability::Internal {
            p.outputs
                .iter()
                .any(|g| fnmatch(basename(&want), basename(&concrete(g))))
        } else {
            published(p).iter().any(|(path, renamed)| {
                if *renamed {
                    fnmatch(basename(&want), basename(path))
                } else {
                    fnmatch(&want, path)
                }
            })
        };
        if !ok {
            let what = if r.stability == Stability::Stable {
                "BREAKING: the stable file"
            } else {
                "the file"
            };
            problems.push(format!(
                "{what} {} (contract line {}) is not what {} publishes ({}): renamed, moved or no longer written",
                r.name, r.line, r.process, p.file
            ));
        }
    }
    let mut cache: ScriptCache = HashMap::new();
    for r in rows {
        if r.kind == Kind::File
            && r.legs == "nextflow"
            && r.stability != Stability::Internal
            && let Some(script) = bash_writes(repo, &r.step, &r.name, &mut cache)?
        {
            problems.push(format!(
                "{} (contract line {}) has legs nextflow, but {script} writes it under the same name too: its legs are both",
                r.name, r.line
            ));
        }
    }
    for (name, p) in &procs {
        if p.publish.is_empty() {
            continue;
        }
        let rows_p: Vec<String> = file_rows
            .iter()
            .filter(|r| &r.process == name)
            .map(|r| concrete(&r.name))
            .collect();
        for g in &p.outputs {
            let c = concrete(g);
            let c = basename(&c);
            if !rows_p.iter().any(|w| fnmatch(basename(w), c)) {
                problems.push(format!(
                    "{name} ({}) writes {g} and has no contract row for it: add one (internal when it is not published)",
                    p.file
                ));
            }
        }
    }
    Ok(problems)
}

// ---------------------------------------------------------------- the run side

/// The lines of a plain or gzip file, invalid UTF-8 replaced.
fn lines(path: &Path) -> io::Result<Lines<Box<dyn BufRead>>> {
    Ok(Lines::new(og_io::open(path)?, Utf8::Replace))
}

/// What a file has: its columns, its INFO/FORMAT/FILTER ids or its
/// top-level keys; `None` for a JSON file that does not parse.
fn fields(path: &Path, fmt: &str) -> io::Result<Option<Vec<String>>> {
    match fmt {
        "tsv" | "csv" => {
            let sep = if fmt == "tsv" { '\t' } else { ',' };
            for (i, line) in lines(path)?.enumerate() {
                if i > 200 {
                    break;
                }
                let line = line?;
                if line.starts_with("##") || !line.contains(sep) {
                    continue;
                }
                return Ok(Some(
                    line.trim_start_matches('#')
                        .split(sep)
                        .map(|c| pycompat::text::strip(c, None).to_owned())
                        .collect(),
                ));
            }
            Ok(Some(Vec::new()))
        }
        "vcf" => {
            let mut out = Vec::new();
            for line in lines(path)? {
                let line = line?;
                if !line.starts_with("##") {
                    break;
                }
                if let Some(c) = VCF_FIELD.captures(&line) {
                    out.push(format!("{}/{}", &c[1], &c[2]));
                }
            }
            Ok(Some(out))
        }
        "json" => {
            let mut raw = Vec::new();
            og_io::open(path)?.read_to_end(&mut raw)?;
            let text = String::from_utf8_lossy(&raw);
            let text = text.replace("\r\n", "\n").replace('\r', "\n");
            Ok(match json::loads(&text) {
                Ok(Value::Object(o)) => Some(o.keys().cloned().collect()),
                Ok(_) => Some(Vec::new()),
                Err(_) => None,
            })
        }
        _ => Ok(Some(Vec::new())),
    }
}

fn guess_format(name: &str) -> &'static str {
    for (ext, fmt) in [
        (".vcf.gz", "vcf"),
        (".vcf", "vcf"),
        (".tsv.gz", "tsv"),
        (".tsv", "tsv"),
        (".csv", "csv"),
        (".json", "json"),
        (".txt", "tsv"),
    ] {
        if name.ends_with(ext) {
            return fmt;
        }
    }
    "-"
}

/// `os.path.join(a, b)` for a relative `b`.
fn join(a: &str, b: &str) -> String {
    if a.is_empty() || a.ends_with('/') {
        format!("{a}{b}")
    } else {
        format!("{a}/{b}")
    }
}

/// `os.path.relpath(path, base)` for a path made by `join(base, …)`.
fn rel_to(path: &str, base: &str) -> String {
    let rest = path
        .strip_prefix(base)
        .map(|r| r.trim_start_matches('/'))
        .unwrap_or(path);
    rest.strip_prefix("./").unwrap_or(rest).to_owned()
}

fn check_outputs(
    rows: &[Row],
    outdir: &str,
    leg: &str,
    sample: &str,
    observe: bool,
    out: &mut dyn Write,
) -> (usize, Vec<String>) {
    let mut problems = Vec::new();
    let mut checked = 0;
    for r in rows {
        if r.kind != Kind::File
            || r.stability != Stability::Stable
            || !(r.e2e == leg || r.e2e == "both")
        {
            continue;
        }
        if !(r.legs == leg || r.legs == "both") {
            problems.push(format!(
                "contract line {}: e2e {} on a {} row",
                r.line, r.e2e, r.legs
            ));
            continue;
        }
        checked += 1;
        let name = r.name.replace("{sample}", sample);
        let mut hits: Vec<String> = glob(join(outdir, &name), false)
            .into_iter()
            .map(|p| p.to_string_lossy().into_owned())
            .collect();
        hits.sort();
        if hits.is_empty() {
            problems.push(format!(
                "BREAKING: the stable file {name} (contract line {}, step {}) is missing from {outdir}: renamed, moved or not written",
                r.line, r.step
            ));
            continue;
        }
        let req: Vec<&str> = if r.required == "-" {
            Vec::new()
        } else {
            r.required.split(',').collect()
        };
        if r.format == "dir" {
            if !hits.iter().all(|h| Path::new(h).is_dir()) {
                problems.push(format!("{name} (contract line {}) is not a folder", r.line));
            }
            continue;
        }
        for h in &hits {
            let rel = rel_to(h, outdir);
            if Path::new(h).is_dir() {
                problems.push(format!(
                    "{rel} (contract line {}) is a folder, want a file",
                    r.line
                ));
                continue;
            }
            if req.is_empty() {
                continue;
            }
            let have = match fields(Path::new(h), &r.format) {
                Ok(Some(h)) => h,
                Ok(None) => {
                    problems.push(format!(
                        "{rel} (contract line {}) is not valid {}",
                        r.line, r.format
                    ));
                    continue;
                }
                Err(e) => {
                    problems.push(format!(
                        "{rel} (contract line {}) is not readable {}: {e}",
                        r.line, r.format
                    ));
                    continue;
                }
            };
            for col in &req {
                if !have.iter().any(|h| h == col) {
                    let list = if have.is_empty() {
                        "nothing".to_owned()
                    } else {
                        have.join(",")
                    };
                    problems.push(format!(
                        "BREAKING: {rel} (contract line {}) has no {col} (it has: {list})",
                        r.line
                    ));
                }
            }
        }
    }
    if observe {
        let base = join(outdir, sample);
        let _ = writeln!(
            out,
            "# observed under {base} (leg {leg}): path, format, fields"
        );
        walk(&base, 0, &mut |p: &str| {
            let rel = rel_to(p, outdir).replace(sample, "{sample}");
            let fmt = guess_format(basename(p));
            let have = if fmt == "-" {
                Vec::new()
            } else {
                match fs::metadata(p) {
                    Ok(m) if m.len() < 2 << 30 => {
                        fields(Path::new(p), fmt).ok().flatten().unwrap_or_default()
                    }
                    _ => Vec::new(),
                }
            };
            let line = format!("OBSERVED\t{leg}\t{rel}\t{fmt}\t{}", have.join(","));
            let line: String = line.chars().take(4000).collect();
            let _ = writeln!(out, "{line}");
        });
    }
    (checked, problems)
}

/// `os.walk(root)` as the observation uses it: the files of each folder in
/// name order, then its subfolders in name order (not `work`, `nf-work`,
/// `.nextflow`, nor a link), four levels below the root at most. The Python
/// stops below a folder whose `relpath(folder, root)` holds three
/// separators; the root (`.`) and its children both hold none, so that is
/// the folder four levels down, and its files are still listed.
fn walk(root: &str, depth: usize, f: &mut dyn FnMut(&str)) {
    let Ok(entries) = fs::read_dir(root) else {
        return;
    };
    let (mut dirs, mut files) = (Vec::new(), Vec::new());
    for e in entries.flatten() {
        let name = e.file_name().to_string_lossy().into_owned();
        if fs::metadata(e.path()).is_ok_and(|m| m.is_dir()) {
            let link = e.file_type().is_ok_and(|t| t.is_symlink());
            dirs.push((name, link));
        } else {
            files.push(name);
        }
    }
    files.sort();
    for name in &files {
        f(&join(root, name));
    }
    if depth >= 4 {
        return;
    }
    dirs.retain(|(d, _)| !matches!(d.as_str(), "work" | "nf-work" | ".nextflow"));
    dirs.sort();
    for (d, link) in dirs {
        if !link {
            walk(&join(root, &d), depth + 1, f);
        }
    }
}

// ---------------------------------------------------------------- the command

#[derive(Default)]
struct Args {
    contract: Option<String>,
    repo: Option<String>,
    static_: bool,
    outputs: Option<String>,
    leg: Option<String>,
    sample: Option<String>,
    observe: bool,
}

fn parse(args: &[String]) -> Result<Option<Args>, String> {
    let mut a = Args::default();
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        let (flag, inline) = match arg.split_once('=') {
            Some((f, v)) if f.starts_with("--") => (f, Some(v.to_owned())),
            _ => (arg.as_str(), None),
        };
        let mut value = |name: &str| -> Result<String, String> {
            match inline.clone() {
                Some(v) => Ok(v),
                None => it
                    .next()
                    .filter(|v| !v.starts_with("--"))
                    .cloned()
                    .ok_or_else(|| format!("argument {name}: expected one argument")),
            }
        };
        match flag {
            "-h" | "--help" => return Ok(None),
            "--contract" => a.contract = Some(value(flag)?),
            "--repo" => a.repo = Some(value(flag)?),
            "--outputs" => a.outputs = Some(value(flag)?),
            "--sample" => a.sample = Some(value(flag)?),
            "--leg" => {
                let v = value(flag)?;
                if v != "bash" && v != "nextflow" {
                    return Err(format!(
                        "argument --leg: invalid choice: {} (choose from bash, nextflow)",
                        repr(&Value::Str(v))
                    ));
                }
                a.leg = Some(v);
            }
            "--static" if inline.is_none() => a.static_ = true,
            "--observe" if inline.is_none() => a.observe = true,
            _ => return Err(format!("unrecognized arguments: {arg}")),
        }
    }
    if !a.static_ && a.outputs.is_none() {
        return Err("give --static or --outputs".to_owned());
    }
    if a.outputs.is_some() && (a.leg.is_none() || a.sample.is_none()) {
        return Err("--outputs needs --leg and --sample".to_owned());
    }
    Ok(Some(a))
}

/// Runs `og dev contract ARGS`; returns the exit code.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let a = match parse(args) {
        Ok(Some(a)) => a,
        Ok(None) => {
            let _ = write!(out, "{USAGE}");
            return 0;
        }
        Err(e) => {
            let _ = write!(err, "{USAGE}");
            let _ = writeln!(err, "og dev contract: error: {e}");
            return 2;
        }
    };
    let repo = a.repo.as_deref().unwrap_or(".");
    let contract = a
        .contract
        .clone()
        .unwrap_or_else(|| join(repo, "assets/output_contract.tsv"));
    let rows = match read_contract(&contract) {
        Ok(r) => r,
        Err(e) => {
            let _ = writeln!(out, "ERROR: {e}");
            return 1;
        }
    };
    let mut problems = Vec::new();
    if a.static_ {
        match check_static(repo, &rows) {
            Ok(p) => {
                let _ = writeln!(
                    out,
                    "static: {} contract rows, {} problems",
                    rows.len(),
                    p.len()
                );
                problems.extend(p);
            }
            Err(e) => {
                let _ = writeln!(out, "ERROR: {e}");
                return 1;
            }
        }
    }
    if let (Some(dir), Some(leg), Some(sample)) = (&a.outputs, &a.leg, &a.sample) {
        let (n, mut p) = check_outputs(&rows, dir, leg, sample, a.observe, out);
        let _ = writeln!(
            out,
            "outputs ({leg}): {n} stable files checked under {dir}, {} problems",
            p.len()
        );
        if n == 0 {
            p.push(format!(
                "no stable row has e2e {leg}: the check checked nothing"
            ));
        }
        problems.extend(p);
    }
    for p in &problems {
        let _ = writeln!(out, "FAIL {p}");
    }
    u8::from(!problems.is_empty())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn braces_and_concrete() {
        assert_eq!(braces("a.{tsv,txt}"), ["a.tsv", "a.txt"]);
        assert_eq!(braces("{a,b}{c,d}"), ["ac", "ad", "bc", "bd"]);
        assert_eq!(braces("none"), ["none"]);
        assert_eq!(
            concrete("${params.outdir}/${meta.id}/fam/${params.family_name ?: 'x'}/{family}.tsv/"),
            "SAMPLE/fam/FAMILY/FAMILY.tsv"
        );
    }

    #[test]
    fn save_as_closures() {
        let c = Some(" fn -> fn.endsWith('.log') ? null : fn ");
        assert_eq!(save_as(c, "a.log"), SaveAs::Drop);
        assert_eq!(save_as(c, "a.tsv"), SaveAs::Keep);
        assert_eq!(
            save_as(Some("x -> x == 'a' ? x : \"b/$x\""), "a"),
            SaveAs::Keep
        );
        assert_eq!(
            save_as(Some("x -> x == 'a' ? x : \"b/$x\""), "c"),
            SaveAs::Rename
        );
        assert_eq!(
            save_as(Some("x -> y.endsWith('a') ? x : null"), "a"),
            SaveAs::Rename
        );
        assert_eq!(save_as(Some("fn ->\n  if (fn) fn\n"), "a"), SaveAs::Rename);
        assert_eq!(save_as(None, "a"), SaveAs::Keep);
    }

    #[test]
    fn a_closure_after_saveas_options_is_read_to_its_own_brace() {
        let stmt = "publishDir \"x\", saveAs: { fn -> fn.endsWith('.log') ? null : fn }, enabled: { true }";
        let body = og_nf::module::closure_body(stmt).unwrap();
        assert!(!body.contains("enabled"));
        assert_eq!(save_as(Some(&body), "a.log"), SaveAs::Drop);
    }

    #[test]
    fn defaults_compare_by_json_type() {
        assert!(!same(&Value::Bool(true), &Value::Int(1)));
        assert!(!same(&Value::Int(0), &Value::Bool(false)));
        assert!(same(&Value::Int(1), &Value::Float(1.0)));
        assert!(same(&Value::Str("a".into()), &Value::Str("a".into())));
    }

    #[test]
    fn a_longer_path_is_not_the_path() {
        let dir = std::env::temp_dir().join(format!("og-contract-bw-{}", std::process::id()));
        fs::create_dir_all(dir.join("scripts")).unwrap();
        fs::write(
            dir.join("scripts/11-roh.sh"),
            "S=\"${GENOME_DIR}/${SAMPLE}\"\nbgzip -c x > \"${S}/roh/${SAMPLE}_roh.txt.gz\"\n",
        )
        .unwrap();
        let repo = dir.to_str().unwrap();
        let mut cache = HashMap::new();
        let short = bash_writes(repo, "11", "{sample}/roh/{sample}_roh.txt", &mut cache).unwrap();
        let full = bash_writes(repo, "11", "{sample}/roh/{sample}_roh.txt.gz", &mut cache).unwrap();
        fs::remove_dir_all(&dir).unwrap();
        assert_eq!(short, None);
        assert_eq!(full.as_deref(), Some("scripts/11-roh.sh"));
    }

    #[test]
    fn the_walk_lists_files_four_folders_down_and_no_deeper() {
        let dir = std::env::temp_dir().join(format!("og-contract-walk-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(dir.join("S1/a/b/c/d/e")).unwrap();
        for p in [
            "S1/top.tsv",
            "S1/a/1.tsv",
            "S1/a/b/c/3.tsv",
            "S1/a/b/c/d/4.tsv",
            "S1/a/b/c/d/e/5.tsv",
        ] {
            fs::write(dir.join(p), "x\n").unwrap();
        }
        let mut seen = Vec::new();
        walk(dir.join("S1").to_str().unwrap(), 0, &mut |p: &str| {
            seen.push(p.rsplit('/').next().unwrap_or(p).to_string());
        });
        fs::remove_dir_all(&dir).unwrap();
        // What os.walk with the oracle's depth cut lists, in its order.
        assert_eq!(seen, ["top.tsv", "1.tsv", "3.tsv", "4.tsv"]);
    }
}
