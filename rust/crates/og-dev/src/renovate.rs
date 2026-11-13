//! `og dev renovate-check LOG [--pins FILE] [--repo DIR]`: checks what a
//! Renovate dry run (`renovate --platform=local --dry-run=lookup`, JSON
//! debug log) found against every pin in the repository. It replaces the
//! Python the `renovate-dry-run.yml` workflow wrote into a file and ran.
//!
//! The checks, each with a planted fault in the self-test that runs after
//! them:
//!
//! - `count`: every `*_IMAGE` line of `versions.env` is a detected docker
//!   dependency;
//! - `names`: no docker package name carries a `:` or `@`;
//! - `twins`: `docs/versions.md` holds exactly the images of `versions.env`,
//!   and `conf/containers.config` holds only those;
//! - `one-pr`: a dependency gets the same proposals in every file it is in;
//! - `pins`: every non-image pin of the pins table
//!   (`scripts/ci/renovate-pins.tsv`) is detected as often as its file pins
//!   it;
//! - `lookups`: no lookup failed, no tag was unreadable, and the only
//!   skipped dependencies are the ones a rule in `renovate.json` disables.
//!
//! The pins table has three tab-separated columns: the dependency, a file
//! or a glob of files (relative to the repository; `*` and `?` in the last
//! part), and a regular expression for one pin of it in that text, read
//! with `^` and `$` at line ends. A file that is missing, or a glob that
//! matches nothing, is an error (exit 2): a pin counted in no file could
//! never fail.
//!
//! Output: the Markdown table of every dependency and its proposed update
//! (also appended to `$GITHUB_STEP_SUMMARY`), one `FAIL <check>: ...` line
//! per problem, the self-test lines, then `OK: N dependencies; every check
//! passed`. Exit 0 when clean, 1 on a finding, a self-test miss, or a log
//! without its `packageFiles with updates` record, 2 on a usage error or a
//! pins table it cannot use.

use indexmap::IndexMap;
use pycompat::json::{self, Value};
use pycompat::value::{py_repr, py_str, truthy};
use regex::Regex;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs;
use std::io::Write;
use std::path::Path;

const USAGE: &str = "usage: og dev renovate-check LOG [--pins FILE] [--repo DIR]\n";
const ENV: &str = "versions.env";
const DOCS: &str = "docs/versions.md";
const CONFIG: &str = "conf/containers.config";
const TWINS: [&str; 2] = [DOCS, CONFIG];
const RULE_SKIPS: [&str; 2] = ["disabled", "package-rules"];
const PINS: &str = "scripts/ci/renovate-pins.tsv";

/// One dependency Renovate found: the file it is in and its record.
type Dep = (String, IndexMap<String, Value>);

/// The pins table, expanded: per (dependency, file) how many pins the text
/// holds, and per dependency the files listed for it.
struct Pins {
    want: BTreeMap<(String, String), usize>,
    files: HashMap<String, BTreeSet<String>>,
    setup_nextflow: Option<Regex>,
}

/// Runs `og dev renovate-check ARGS`; returns the exit code.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let summary = std::env::var_os("GITHUB_STEP_SUMMARY").filter(|p| !p.is_empty());
    cli_with_summary(args, out, err, summary.as_deref().map(Path::new))
}

fn cli_with_summary(
    args: &[String],
    out: &mut dyn Write,
    err: &mut dyn Write,
    summary_path: Option<&Path>,
) -> u8 {
    let (mut log, mut pins, mut repo) = (None, PINS.to_owned(), ".".to_owned());
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "-h" | "--help" => {
                let _ = out.write_all(USAGE.as_bytes());
                return 0;
            }
            "--pins" | "--repo" => {
                let Some(v) = it.next() else {
                    let _ = writeln!(
                        err,
                        "{USAGE}og dev renovate-check: error: {a} needs a value"
                    );
                    return 2;
                };
                if a == "--pins" {
                    pins = v.clone();
                } else {
                    repo = v.clone();
                }
            }
            _ if log.is_none() => log = Some(a.clone()),
            _ => {
                let _ = writeln!(
                    err,
                    "{USAGE}og dev renovate-check: error: unexpected argument '{a}'"
                );
                return 2;
            }
        }
    }
    let Some(log) = log else {
        let _ = writeln!(
            err,
            "{USAGE}og dev renovate-check: error: the following arguments are required: LOG"
        );
        return 2;
    };
    let repo = Path::new(&repo);
    let pins = match read_pins(repo, &repo.join(&pins)) {
        Ok(p) => p,
        Err(e) => {
            let _ = writeln!(err, "og dev renovate-check: {e}");
            return 2;
        }
    };
    let env_text = match crate::evidence::read_text(&repo.join(ENV)) {
        Ok(t) => t,
        Err(e) => {
            let _ = writeln!(err, "og dev renovate-check: {ENV}: {e}");
            return 2;
        }
    };
    let deps = match load(&log) {
        Ok(d) => d,
        Err(e) => {
            let _ = writeln!(err, "{e}");
            return 1;
        }
    };
    let table = summary(&deps);
    let _ = writeln!(out, "{table}");
    if let Some(path) = summary_path {
        let appended = (|| -> Result<(), String> {
            let (pending, mut file) = crate::evidence::PendingOutput::new(path)?;
            match fs::File::open(path) {
                Ok(mut old) => {
                    std::io::copy(&mut old, &mut file).map_err(|e| e.to_string())?;
                }
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
                Err(e) => return Err(e.to_string()),
            }
            writeln!(file, "## Renovate dry run\n\n{table}").map_err(|e| e.to_string())?;
            file.flush().map_err(|e| e.to_string())?;
            drop(file);
            pending.commit()
        })();
        if let Err(e) = appended {
            let _ = writeln!(err, "og dev renovate-check: GITHUB_STEP_SUMMARY: {e}");
            return 1;
        }
    }
    let mut fails = check(&deps, &env_text, &pins);
    for x in &fails {
        let _ = writeln!(out, "FAIL {x}");
    }
    if !self_test(&deps, &env_text, &pins, out) {
        fails.push("the self-test did not catch a planted fault".to_owned());
        let _ = writeln!(out, "FAIL the self-test did not catch a planted fault");
    }
    if !fails.is_empty() {
        return 1;
    }
    let _ = writeln!(out, "OK: {} dependencies; every check passed", deps.len());
    0
}

/// `fnmatch` for `*` and `?`, case-sensitive.
fn fnmatch(pat: &[u8], name: &[u8]) -> bool {
    match pat.split_first() {
        None => name.is_empty(),
        Some((b'*', rest)) => (0..=name.len()).any(|i| fnmatch(rest, &name[i..])),
        Some((b'?', rest)) => !name.is_empty() && fnmatch(rest, &name[1..]),
        Some((c, rest)) => name.first() == Some(c) && fnmatch(rest, &name[1..]),
    }
}

/// The files of `pattern` under `repo`, relative and sorted, as Python's
/// `sorted(glob.glob(pattern))` lists them: a wildcard never matches a name
/// that starts with a dot unless the pattern does.
fn glob(repo: &Path, pattern: &str) -> Result<Vec<String>, String> {
    let (dir, last) = pattern.rsplit_once('/').unwrap_or(("", pattern));
    if dir.contains(['*', '?']) {
        return Err(format!(
            "{pattern}: a wildcard is allowed only in the last part"
        ));
    }
    if !last.contains(['*', '?']) {
        return if repo.join(pattern).is_file() {
            Ok(vec![pattern.to_owned()])
        } else {
            Err(format!("{pattern} is not a file"))
        };
    }
    let mut out = Vec::new();
    let entries = fs::read_dir(repo.join(if dir.is_empty() { "." } else { dir }))
        .map_err(|e| format!("{pattern}: {e}"))?;
    for e in entries {
        let name = e.map_err(|e| format!("{pattern}: {e}"))?.file_name();
        let name = name.to_string_lossy();
        if (name.starts_with('.') && !last.starts_with('.'))
            || !fnmatch(last.as_bytes(), name.as_bytes())
        {
            continue;
        }
        out.push(if dir.is_empty() {
            name.into_owned()
        } else {
            format!("{dir}/{name}")
        });
    }
    out.sort();
    if out.is_empty() {
        return Err(format!("{pattern} matches no file"));
    }
    Ok(out)
}

fn read_pins(repo: &Path, path: &Path) -> Result<Pins, String> {
    let text = crate::evidence::read_text(path)?;
    let mut pins = Pins {
        want: BTreeMap::new(),
        files: HashMap::new(),
        setup_nextflow: None,
    };
    let mut texts: HashMap<String, String> = HashMap::new();
    for (i, line) in text.lines().enumerate() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let at = |msg: &str| format!("{}:{}: {msg}", path.display(), i + 1);
        let cols: Vec<&str> = line.split('\t').collect();
        let [dep, files, pattern] = cols[..] else {
            return Err(at(
                "want three tab-separated columns: dependency, file or glob, pattern",
            ));
        };
        let re = Regex::new(&format!("(?m){pattern}")).map_err(|e| at(&e.to_string()))?;
        if pattern.contains("setup-nextflow") {
            pins.setup_nextflow = Some(re.clone());
        }
        for f in glob(repo, files).map_err(|e| at(&e))? {
            if !texts.contains_key(&f) {
                let t = crate::evidence::read_text(&repo.join(&f)).map_err(|e| at(&e))?;
                texts.insert(f.clone(), t);
            }
            let n = re.find_iter(&texts[&f]).count();
            if n > 0 {
                *pins.want.entry((dep.to_owned(), f.clone())).or_default() += n;
            }
            pins.files.entry(dep.to_owned()).or_default().insert(f);
        }
    }
    Ok(pins)
}

/// Every dependency of the log's last `packageFiles with updates` record.
fn load(log: &str) -> Result<Vec<Dep>, String> {
    let mut files = None;
    crate::evidence::lines(log, |line| {
        // The old load() skipped only JSON parse errors, not parsed scalars.
        let Ok(value) = json::loads(std::str::from_utf8(line).expect("strict reader")) else {
            return Ok(());
        };
        let Value::Object(mut rec) = value else {
            return Err(format!("FAIL: {log}: a log record is not an object"));
        };
        if rec.get("msg") == Some(&Value::Str("packageFiles with updates".to_owned())) {
            files = Some(
                rec.swap_remove("config")
                    .ok_or_else(|| format!("FAIL: {log}: missing config"))?,
            );
        }
        Ok(())
    })?;
    let Some(files) = files else {
        return Err(format!(
            "FAIL: no 'packageFiles with updates' record in {log}"
        ));
    };
    let shape = || {
        format!(
            "FAIL: the 'packageFiles with updates' record of {log} is not a map of package file lists"
        )
    };
    let Value::Object(files) = files else {
        return Err(shape());
    };
    let mut deps = Vec::new();
    for entries in files.into_values() {
        let Value::Array(entries) = entries else {
            return Err(shape());
        };
        for pf in entries {
            let Value::Object(mut pf) = pf else {
                return Err(shape());
            };
            let ds = match pf.swap_remove("deps") {
                None => Vec::new(),
                Some(Value::Array(ds)) => ds,
                _ => return Err(shape()),
            };
            let name = match pf.swap_remove("packageFile") {
                None => None,
                Some(Value::Str(name)) => Some(name),
                // An unused field may be absent, but a present wrong type is
                // still a deliberate input error.
                _ => return Err(shape()),
            };
            if ds.is_empty() {
                continue;
            }
            let name = name.ok_or_else(shape)?;
            for d in ds {
                let Value::Object(d) = d else {
                    return Err(shape());
                };
                // Optional means absent. Null or a scalar is not an empty list.
                if let Some(u) = d.get("updates") {
                    let Value::Array(us) = u else {
                        return Err(shape());
                    };
                    for u in us {
                        let Value::Object(u) = u else {
                            return Err(shape());
                        };
                        for k in ["newValue", "newDigest", "updateType"] {
                            if u.get(k)
                                .is_some_and(|v| !matches!(v, Value::Null | Value::Str(_)))
                            {
                                return Err(shape());
                            }
                        }
                    }
                }
                for k in [
                    "datasource",
                    "depName",
                    "currentValue",
                    "currentDigest",
                    "skipReason",
                ] {
                    if d.get(k)
                        .is_some_and(|v| !matches!(v, Value::Null | Value::Str(_)))
                    {
                        return Err(shape());
                    }
                }
                if let Some(ws) = d.get("warnings").filter(|v| truthy(v)) {
                    let Value::Array(ws) = ws else {
                        return Err(shape());
                    };
                    if ws.iter().any(|v| !matches!(v, Value::Object(_))) {
                        return Err(shape());
                    }
                }
                deps.push((name.clone(), d));
            }
        }
    }
    Ok(deps)
}

fn get<'a>(d: &'a IndexMap<String, Value>, k: &str) -> &'a Value {
    d.get(k).unwrap_or(&Value::Null)
}

/// Python's `"%s" % v`.
fn s(v: &Value) -> String {
    py_str(v).unwrap_or_else(|_| format!("{v:?}"))
}

/// Python's `repr(v)`.
fn r(v: &Value) -> String {
    py_repr(v).unwrap_or_else(|_| format!("{v:?}"))
}

/// `repr` of a tuple of values.
fn tuple(vs: &[&Value]) -> String {
    let items: Vec<String> = vs.iter().map(|v| r(v)).collect();
    format!("({})", items.join(", "))
}

fn is_docker(d: &IndexMap<String, Value>) -> bool {
    get(d, "datasource") == &Value::Str("docker".to_owned())
}

/// The image names of `versions.env`: `name:tag@sha256:...`,
/// `name@sha256:...` and `name:tag` all give `name`.
fn env_images(text: &str) -> Vec<String> {
    let line = Regex::new(r#"(?m)^[A-Z0-9_]+_IMAGE="([^"]*)""#).expect("a valid pattern");
    let tag = Regex::new(r":[^:/]*$").expect("a valid pattern");
    line.captures_iter(text)
        .map(|c| {
            let v = c[1].split('@').next().unwrap_or("");
            tag.replace(v, "").into_owned()
        })
        .collect()
}

/// `(depName, currentValue, currentDigest)` as Python prints the tuple.
fn key(d: &IndexMap<String, Value>) -> String {
    tuple(&[
        get(d, "depName"),
        get(d, "currentValue"),
        get(d, "currentDigest"),
    ])
}

/// `v or ""`.
fn or_empty(v: &Value) -> Value {
    if truthy(v) {
        v.clone()
    } else {
        Value::Str(String::new())
    }
}

/// The updates as sorted `(newValue or "", newDigest or "")` pairs, printed
/// as Python prints the list.
fn proposals(d: &IndexMap<String, Value>) -> String {
    let mut ps: Vec<(Value, Value)> = updates(d)
        .iter()
        .map(|u| match u {
            Value::Object(u) => (or_empty(get(u, "newValue")), or_empty(get(u, "newDigest"))),
            _ => (Value::Null, Value::Null),
        })
        .collect();
    ps.sort_by(|a, b| (s(&a.0), s(&a.1)).cmp(&(s(&b.0), s(&b.1))));
    let items: Vec<String> = ps.iter().map(|(a, b)| tuple(&[a, b])).collect();
    format!("[{}]", items.join(", "))
}

fn updates(d: &IndexMap<String, Value>) -> &[Value] {
    match d.get("updates") {
        Some(Value::Array(u)) => u,
        _ => &[],
    }
}

fn counter<'a>(it: impl Iterator<Item = &'a IndexMap<String, Value>>) -> BTreeMap<String, usize> {
    let mut c = BTreeMap::new();
    for d in it {
        *c.entry(key(d)).or_default() += 1;
    }
    c
}

/// The docker dependencies of one file.
fn in_file<'a>(
    docker: &'a [&'a Dep],
    f: &'a str,
) -> impl Iterator<Item = &'a IndexMap<String, Value>> {
    docker.iter().filter(move |(g, _)| g == f).map(|(_, d)| d)
}

/// Every problem, as `check: ...` lines.
fn check(deps: &[Dep], env_text: &str, pins: &Pins) -> Vec<String> {
    let mut fails = Vec::new();
    let docker: Vec<&Dep> = deps.iter().filter(|(_, d)| is_docker(d)).collect();
    let want = env_images(env_text);
    let got: Vec<&Value> = in_file(&docker, ENV).map(|d| get(d, "depName")).collect();
    if got.len() != want.len() {
        fails.push(format!(
            "count: versions.env has {} *_IMAGE lines, Renovate found {} docker dependencies there",
            want.len(),
            got.len()
        ));
    }
    let want_set: BTreeSet<&String> = want.iter().collect();
    for name in want_set {
        if !got.contains(&&Value::Str(name.clone())) {
            fails.push(format!(
                "count: {name} (versions.env) is not a detected dependency"
            ));
        }
    }
    for (f, d) in &docker {
        let name = get(d, "depName");
        let text = if truthy(name) {
            s(name)
        } else {
            ":".to_owned()
        };
        if text.contains([':', '@']) {
            fails.push(format!("names: {f}: package name {}", r(name)));
        }
    }
    let env_keys = counter(in_file(&docker, ENV));
    for f in TWINS {
        let keys = counter(in_file(&docker, f));
        if f == DOCS && keys != env_keys {
            let all: BTreeSet<&String> = keys.keys().chain(env_keys.keys()).collect();
            for k in all {
                if keys.get(k) != env_keys.get(k) {
                    fails.push(format!(
                        "twins: docs/versions.md and versions.env differ on {k}"
                    ));
                }
            }
        } else {
            for k in keys.keys().filter(|k| !env_keys.contains_key(*k)) {
                fails.push(format!("twins: {f} has {k}, which versions.env does not"));
            }
        }
    }
    let mut seen: HashMap<String, (&str, String)> = HashMap::new();
    for (f, d) in deps {
        let name = tuple(&[
            get(d, "datasource"),
            get(d, "depName"),
            get(d, "currentValue"),
            get(d, "currentDigest"),
        ]);
        let p = proposals(d);
        match seen.get(&name) {
            Some((first, q)) if *q != p => fails.push(format!(
                "one-pr: {} gets {q} in {first} but {p} in {f}",
                s(get(d, "depName"))
            )),
            Some(_) => {}
            None => {
                seen.insert(name, (f, p));
            }
        }
    }
    for (f, d) in deps {
        let warnings = get(d, "warnings");
        if truthy(warnings) {
            let msgs: Vec<String> = match warnings {
                Value::Array(ws) => ws
                    .iter()
                    .map(|w| match w {
                        Value::Object(o) if o.contains_key("message") => s(get(o, "message")),
                        w => s(w),
                    })
                    .collect(),
                w => vec![s(w)],
            };
            fails.push(format!(
                "lookups: {} ({f}): {}",
                s(get(d, "depName")),
                msgs.join("; ")
            ));
        }
        let skip = get(d, "skipReason");
        if truthy(skip) && !RULE_SKIPS.contains(&s(skip).as_str()) {
            fails.push(format!(
                "lookups: {} ({f}) skipped: {}",
                s(get(d, "depName")),
                s(skip)
            ));
        }
    }
    let mut got: BTreeMap<(String, String), usize> = BTreeMap::new();
    for (f, d) in deps {
        if let Value::Str(name) = get(d, "depName")
            && pins.files.get(name).is_some_and(|fs| fs.contains(f))
        {
            *got.entry((name.clone(), f.clone())).or_default() += 1;
        }
    }
    let all: BTreeSet<&(String, String)> = pins.want.keys().chain(got.keys()).collect();
    for k in all {
        let (w, g) = (
            pins.want.get(k).copied().unwrap_or(0),
            got.get(k).copied().unwrap_or(0),
        );
        if w != g {
            fails.push(format!(
                "pins: {} in {}: the text pins it {w} times, Renovate found {g}",
                k.0, k.1
            ));
        }
    }
    fails
}

/// Plants one fault of each kind in a copy of the dependencies and shows
/// that the check it belongs to catches it; true when every one is caught.
fn self_test(deps: &[Dep], env_text: &str, pins: &Pins, out: &mut dyn Write) -> bool {
    let mut ok = true;
    let base: BTreeSet<String> = check(deps, env_text, pins).into_iter().collect();
    let first = |d: &[Dep], f: &str| d.iter().position(|(g, x)| g == f && is_docker(x));
    type Plant = fn(&mut Vec<Dep>, usize);
    let plants: [(&str, &str, Option<&str>, Plant); 6] = [
        ("count", "count:", Some(ENV), |d, i| {
            d.remove(i);
        }),
        ("names", "names:", Some(ENV), |d, i| {
            let name = format!("{}@sha256", s(get(&d[i].1, "depName")));
            d[i].1.insert("depName".to_owned(), Value::Str(name));
        }),
        ("twins", "twins:", Some(DOCS), |d, i| {
            d[i].1.insert(
                "currentValue".to_owned(),
                Value::Str("0.0.0-planted".to_owned()),
            );
        }),
        ("one-pr", "one-pr:", Some(CONFIG), |d, i| {
            let u = IndexMap::from([(
                "newValue".to_owned(),
                Value::Str("9.9.9-planted".to_owned()),
            )]);
            d[i].1
                .insert("updates".to_owned(), Value::Array(vec![Value::Object(u)]));
        }),
        ("pins", "pins:", None, |d, i| {
            d.remove(i);
        }),
        ("lookups", "lookups:", Some(ENV), |d, i| {
            d[i].1.insert(
                "skipReason".to_owned(),
                Value::Str("invalid-value".to_owned()),
            );
        }),
    ];
    for (label, pattern, file, plant) in plants {
        let at = match file {
            Some(f) => first(deps, f),
            None => deps.iter().position(|(g, x)| {
                g == ENV && get(x, "depName") == &Value::Str("nextflow-io/nextflow".to_owned())
            }),
        };
        let Some(i) = at else {
            let _ = writeln!(
                out,
                "self-test {label:<8} NOT PLANTED: no dependency to plant it on"
            );
            ok = false;
            continue;
        };
        let mut d = deps.to_vec();
        plant(&mut d, i);
        // Only a failure the plant added counts as caught.
        let hit = check(&d, env_text, pins)
            .into_iter()
            .find(|x| x.starts_with(pattern) && !base.contains(x));
        match hit {
            Some(h) => {
                let _ = writeln!(out, "self-test {label:<8} caught: {h}");
            }
            None => {
                let _ = writeln!(out, "self-test {label:<8} NOT CAUGHT");
                ok = false;
            }
        }
    }
    // The setup-nextflow pattern itself, on planted workflow text.
    let Some(setup) = &pins.setup_nextflow else {
        let _ = writeln!(
            out,
            "self-test pin-text  NOT PLANTED: the pins table has no setup-nextflow pattern"
        );
        return false;
    };
    let computed = "${{ steps.nf.outputs.version }}";
    let cases: [(&str, &str, &str, &str, usize); 7] = [
        ("computed", "", "", computed, 0),
        ("computed, uses quoted", "\"", "", computed, 0),
        ("computed after a blank line", "", "\n", computed, 0),
        ("literal", "", "", "\"25.10.8\"", 1),
        ("literal after a blank line", "", "\n", "\"25.10.8\"", 1),
        ("single-quoted literal", "'", "", "'25.10.8'", 1),
        ("bare literal, uses quoted", "\"", "", "25.10.8", 1),
    ];
    for (label, quote, gap, version, want) in cases {
        let text = format!(
            "      - uses: {quote}nf-core/setup-nextflow@{a}{quote} # v3.0.1\n{gap}        if: steps.scope.outputs.run == 'true'\n        with:\n          version: {version}\n\n      - uses: actions/setup-java@{b} # v6\n        with:\n          java-version: \"17\"\n",
            a = "a".repeat(40),
            b = "b".repeat(40),
        );
        let n = setup.find_iter(&text).count();
        let _ = writeln!(
            out,
            "self-test pin-text {label:<27} counted {n}, want {want}"
        );
        ok = ok && n == want;
    }
    let clean = check(deps, env_text, pins);
    let verdict = clean
        .first()
        .map_or_else(|| "passes".to_owned(), |c| format!("fails: {c}"));
    let _ = writeln!(out, "self-test clean    {verdict}");
    ok
}

/// The Markdown table of every dependency, once per datasource and name.
fn summary(deps: &[Dep]) -> String {
    let id = |d: &IndexMap<String, Value>| tuple(&[get(d, "datasource"), get(d, "depName")]);
    let mut files: HashMap<String, BTreeSet<&str>> = HashMap::new();
    for (f, d) in deps {
        files.entry(id(d)).or_default().insert(f);
    }
    let mut rows = vec![
        "| Dependency | Current | Proposed update | Skipped | Files |".to_owned(),
        "|---|---|---|---|---|".to_owned(),
    ];
    let mut seen = BTreeSet::new();
    for (_, d) in deps {
        let k = id(d);
        if !seen.insert(k.clone()) {
            continue;
        }
        let cur = [get(d, "currentValue"), get(d, "currentDigest")]
            .into_iter()
            .find(|v| truthy(v))
            .map_or_else(String::new, s);
        let ups: Vec<String> = updates(d)
            .iter()
            .map(|u| {
                let u = match u {
                    Value::Object(u) => u,
                    _ => return s(u),
                };
                let new = if truthy(get(u, "newValue")) {
                    get(u, "newValue")
                } else {
                    get(u, "newDigest")
                };
                format!("{} ({})", s(new), s(get(u, "updateType")))
            })
            .collect();
        let ups = if ups.is_empty() {
            "none".to_owned()
        } else {
            ups.join(", ")
        };
        let skip = get(d, "skipReason");
        let note = if truthy(skip) { s(skip) } else { String::new() };
        rows.push(format!(
            "| `{}` | `{cur}` | {ups} | {note} | {} |",
            s(get(d, "depName")),
            files[&k].len()
        ));
    }
    rows.join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fnmatch_and_env_images() {
        assert!(
            fnmatch(b"*.yml", b"e2e.yml")
                && !fnmatch(b"*.yml", b"e2e.yaml")
                && fnmatch(b"a?c", b"abc")
        );
        let env = "A_IMAGE=\"x/y:1.2@sha256:ab\"\nB_IMAGE=\"z@sha256:cd\"\nC_IMAGE=\"host:5000/w:3\"\nD=\"1\"\n";
        assert_eq!(env_images(env), ["x/y", "z", "host:5000/w"]);
    }

    #[test]
    fn keys_and_proposals_print_as_python() {
        let d: IndexMap<String, Value> = match json::loads(
            r#"{"depName":"it's","currentValue":"1","updates":[{"newValue":"3"},{"newValue":"2","newDigest":"d"}]}"#,
        ) {
            Ok(Value::Object(o)) => o,
            other => panic!("{other:?}"),
        };
        assert_eq!(key(&d), r#"("it's", '1', None)"#);
        assert_eq!(proposals(&d), "[('2', 'd'), ('3', '')]");
    }

    /// A repository of the files the check reads, and a log that matches it.
    fn planted(name: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("og-renovate-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        for sub in ["docs", "conf", "scripts/ci", ".github/workflows"] {
            fs::create_dir_all(d.join(sub)).expect("mkdir");
        }
        let w = |p: &str, t: &str| fs::write(d.join(p), t).expect("write");
        w(
            ENV,
            "NEXTFLOW_VERSION=\"25.10.8\"\nA_IMAGE=\"x/a:1@sha256:aa\"\nB_IMAGE=\"x/b:2\"\n",
        );
        w(DOCS, "| a | x/a:1 |\n");
        w(CONFIG, "container = 'x/a:1'\n");
        w(
            ".github/workflows/ci.yml",
            "      - uses: nf-core/setup-nextflow@abc # v3\n        with:\n          version: \"25.10.8\"\n",
        );
        w(
            ".github/workflows/.hidden.yml",
            "      - uses: nf-core/setup-nextflow@abc # v3\n        with:\n          version: \"25.10.8\"\n",
        );
        w(
            PINS,
            "# comment\nnextflow-io/nextflow\tversions.env\t^NEXTFLOW_VERSION=\"\nnextflow-io/nextflow\t.github/workflows/*.yml\tuses: [\"\\x27]?nf-core/setup-nextflow@[^\\n]*\\n(?:[ \\t]*(?:[^\\s-][^\\n]*)?\\n)*?[ \\t]+version: [\"\\x27]?\\d\n",
        );
        let dep = |n: &str, v: &str, ds: &str| {
            format!(r#"{{"depName":"{n}","currentValue":"{v}","datasource":"{ds}"}}"#)
        };
        let log = format!(
            "not json\n{{\"msg\":\"other\"}}\n{{\"msg\":\"packageFiles with updates\",\"config\":{{\"regex\":[{{\"packageFile\":\"versions.env\",\"deps\":[{},{},{}]}},{{\"packageFile\":\"docs/versions.md\",\"deps\":[{},{}]}},{{\"packageFile\":\"conf/containers.config\",\"deps\":[{}]}},{{\"packageFile\":\".github/workflows/ci.yml\",\"deps\":[{}]}}]}}}}\n",
            dep("x/a", "1", "docker"),
            dep("x/b", "2", "docker"),
            dep("nextflow-io/nextflow", "25.10.8", "github-releases"),
            dep("x/a", "1", "docker"),
            dep("x/b", "2", "docker"),
            dep("x/a", "1", "docker"),
            dep("nextflow-io/nextflow", "25.10.8", "github-releases"),
        );
        w("renovate.log", &log);
        d
    }

    fn run(d: &Path) -> (u8, String) {
        let args: Vec<String> = [
            d.join("renovate.log").to_string_lossy().into_owned(),
            "--repo".to_owned(),
            d.to_string_lossy().into_owned(),
        ]
        .into();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = cli(&args, &mut out, &mut err);
        (
            code,
            String::from_utf8_lossy(&out).into_owned() + &String::from_utf8_lossy(&err),
        )
    }

    #[test]
    fn a_clean_tree_passes_and_every_plant_is_caught() {
        let d = planted("clean");
        let (code, out) = run(&d);
        assert_eq!(code, 0, "{out}");
        assert!(
            out.contains("OK: 7 dependencies; every check passed"),
            "{out}"
        );
        assert!(out.contains("| `x/a` | `1` | none |  | 3 |"), "{out}");
        assert_eq!(out.matches("caught: ").count(), 6, "{out}");
        let _ = fs::remove_dir_all(&d);
    }

    #[test]
    fn a_pin_renovate_misses_fails() {
        let d = planted("pin");
        fs::write(
            d.join(".github/workflows/more.yml"),
            "      - uses: 'nf-core/setup-nextflow@abc'\n\n        with:\n          version: '26.04.0'\n",
        )
        .expect("write");
        let (code, out) = run(&d);
        assert_eq!(code, 1, "{out}");
        assert!(
            out.contains("FAIL pins: nextflow-io/nextflow in .github/workflows/more.yml: the text pins it 1 times, Renovate found 0"),
            "{out}"
        );
        let _ = fs::remove_dir_all(&d);
    }

    #[test]
    fn an_image_renovate_misses_and_a_stale_twin_fail() {
        let d = planted("image");
        let log = fs::read_to_string(d.join("renovate.log")).expect("read");
        let log = log.replacen(
            r#"{"depName":"x/b","currentValue":"2","datasource":"docker"},"#,
            "",
            1,
        );
        let log = log.replacen(
            r#""docs/versions.md","deps":[{"depName":"x/a","currentValue":"1""#,
            r#""docs/versions.md","deps":[{"depName":"x/a","currentValue":"0""#,
            1,
        );
        fs::write(d.join("renovate.log"), log).expect("write");
        let (code, out) = run(&d);
        assert_eq!(code, 1, "{out}");
        assert!(
            out.contains("FAIL count: x/b (versions.env) is not a detected dependency"),
            "{out}"
        );
        assert!(
            out.contains(
                "FAIL twins: docs/versions.md and versions.env differ on ('x/a', '0', None)"
            ),
            "{out}"
        );
        let _ = fs::remove_dir_all(&d);
    }

    #[test]
    fn a_glob_that_matches_nothing_is_refused() {
        let d = planted("glob");
        fs::write(d.join(PINS), "aldy\tdocs/*.txt\taldy==\\d\n").expect("write");
        let (code, out) = run(&d);
        assert_eq!(code, 2, "{out}");
        assert!(out.contains("docs/*.txt matches no file"), "{out}");
        let _ = fs::remove_dir_all(&d);
    }

    #[test]
    fn configured_summary_write_failure_fails_the_check() {
        let d = planted("summary");
        let args = vec![
            d.join("renovate.log").to_str().unwrap().to_owned(),
            "--repo".to_owned(),
            d.to_str().unwrap().to_owned(),
        ];
        let (mut out, mut err) = (Vec::new(), Vec::new());
        assert_eq!(cli_with_summary(&args, &mut out, &mut err, Some(&d)), 1);
        assert!(!String::from_utf8_lossy(&out).contains("every check passed"));
        fs::remove_dir_all(d).unwrap();
    }

    #[test]
    fn malformed_log_shapes_and_text_fail() {
        for (label, bad) in [
            ("null-deps", r#"{"packageFile":"unused","deps":null}"#),
            (
                "null-updates",
                r#"{"packageFile":"unused","deps":[{"updates":null}]}"#,
            ),
            (
                "scalar-update",
                r#"{"packageFile":"unused","deps":[{"updates":[42]}]}"#,
            ),
            ("scalar-dep", r#"{"packageFile":"unused","deps":[42]}"#),
            ("scalar-file", r#"{"packageFile":[],"deps":[]}"#),
            (
                "map-updates",
                r#"{"packageFile":"unused","deps":[{"updates":{}}]}"#,
            ),
            (
                "number-update",
                r#"{"packageFile":"unused","deps":[{"updates":[{"newValue":0}]}]}"#,
            ),
            (
                "number-name",
                r#"{"packageFile":"unused","deps":[{"depName":0}]}"#,
            ),
            (
                "number-current",
                r#"{"packageFile":"unused","deps":[{"currentValue":0}]}"#,
            ),
        ] {
            let d = planted(label);
            let log = fs::read_to_string(d.join("renovate.log")).unwrap();
            let log = log.replacen(r#""regex":["#, &format!(r#""regex":[{bad},"#), 1);
            fs::write(d.join("renovate.log"), log).unwrap();
            let (code, out) = run(&d);
            assert_eq!(code, 1, "{label}: {out}");
            assert!(
                out.contains("not a map of package file lists"),
                "{label}: {out}"
            );
            fs::remove_dir_all(d).unwrap();
        }
        for (label, prefix) in [
            ("utf8", b"\xff\n".as_slice()),
            ("array", b"[]\n".as_slice()),
            ("null", b"null\n".as_slice()),
        ] {
            let d = planted(label);
            let log = [prefix, &fs::read(d.join("renovate.log")).unwrap()].concat();
            fs::write(d.join("renovate.log"), log).unwrap();
            let (code, out) = run(&d);
            assert_eq!(code, 1, "{label}: {out}");
            fs::remove_dir_all(d).unwrap();
        }
    }

    #[test]
    fn repository_sources_use_universal_newlines() {
        for name in [ENV, PINS, ".github/workflows/ci.yml"] {
            for ending in ["\r", "\r\n"] {
                let d = planted(&format!(
                    "newlines-{}-{}",
                    name.replace('/', "-"),
                    ending.len()
                ));
                let baseline = run(&d);
                let text = fs::read_to_string(d.join(name))
                    .unwrap()
                    .replace('\n', ending);
                fs::write(d.join(name), text).unwrap();
                assert_eq!(run(&d), baseline, "{name}/{ending:?}");
                fs::remove_dir_all(d).unwrap();
            }
        }
    }

    #[test]
    fn real_log_falsy_warnings_and_absent_empty_packages_are_accepted() {
        for entry in [
            r#"{}"#,
            r#"{"deps":[]}"#,
            r#"{"packageFile":"unused","deps":[{"warnings":null}]}"#,
            r#"{"packageFile":"unused","deps":[{"warnings":false}]}"#,
            r#"{"packageFile":"unused","deps":[{"warnings":""}]}"#,
            r#"{"packageFile":"unused","deps":[{"warnings":{}}]}"#,
        ] {
            let d = planted("falsy");
            let log = fs::read_to_string(d.join("renovate.log")).unwrap();
            fs::write(
                d.join("renovate.log"),
                log.replacen(r#""regex":["#, &format!(r#""regex":[{entry},"#), 1),
            )
            .unwrap();
            let (code, out) = run(&d);
            assert_eq!(code, 0, "{entry}: {out}");
            fs::remove_dir_all(d).unwrap();
        }
    }

    #[test]
    fn repository_sources_reject_gzip_even_under_plain_names() {
        for name in [ENV, PINS, ".github/workflows/ci.yml"] {
            let d = planted(&format!("compressed-{}", name.replace('/', "-")));
            let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
            gz.write_all(&fs::read(d.join(name)).unwrap()).unwrap();
            fs::write(d.join(name), gz.finish().unwrap()).unwrap();
            let (code, out) = run(&d);
            assert_eq!(code, 2, "{out}");
            assert!(out.contains("wrong compression"), "{out}");
            fs::remove_dir_all(d).unwrap();
        }
    }

    #[test]
    fn renovate_requires_plain_compression_and_pins_record_limit() {
        let d = planted("compression-limit");
        let text = fs::read(d.join("renovate.log")).unwrap();
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gz.write_all(&text).unwrap();
        fs::write(d.join("renovate.log"), gz.finish().unwrap()).unwrap();
        let (code, out) = run(&d);
        assert_eq!(code, 1);
        assert!(out.contains("compression"), "{out}");
        fs::write(d.join("renovate.log"), vec![b'x'; (16 << 20) + 1]).unwrap();
        let (code, out) = run(&d);
        assert_eq!(code, 1);
        assert!(out.contains("record exceeds 16777216 bytes"), "{out}");
        fs::remove_dir_all(d).unwrap();
    }

    #[test]
    fn a_log_without_the_record_fails() {
        let d = planted("nolog");
        fs::write(d.join("renovate.log"), "{\"msg\":\"other\"}\n").expect("write");
        let (code, out) = run(&d);
        assert_eq!(code, 1, "{out}");
        assert!(
            out.contains("FAIL: no 'packageFiles with updates' record in"),
            "{out}"
        );
        let _ = fs::remove_dir_all(&d);
    }
}
