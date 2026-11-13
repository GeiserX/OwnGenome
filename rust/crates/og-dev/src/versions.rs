//! `og dev versions check-run` and `og dev versions expand`: the version
//! checks of `scripts/ci/gen-containers-config.sh --check-versions` and
//! `scripts/ci/check-links.sh` that were Python.
//!
//! ```text
//! og dev versions check-run --config conf/containers.config FILE [PROCESS...]
//! og dev versions expand --env versions.env < URL<TAB>LOC lines
//! ```
//!
//! `check-run` proves a run used the images `conf/containers.config` gives
//! each process. Each module writes, as the first line of its block in
//! `pipeline_info/software_versions.yml` (FILE), the tag of the image it ran
//! in; a process that runs on the host (an `ext.pipeline_version` selector)
//! writes the release of the pipeline it starts. One `OK   PROC key: value`
//! or `FAIL PROC ...` line per process of PROCESS... (default: every
//! selector) and of FILE. Exit 0 when all agree, 1 otherwise, 2 on a usage
//! error or a file that cannot be read.
//!
//! `expand` fills `${NAME}` and `$NAME` in each URL from `versions.env` (a
//! missing file fills nothing), drops trailing `.,;:]` and a `}` that closes
//! nothing, and prints the lines back. Exit 2 on a line with no tab or a
//! `versions.env` the strict reader refuses: on purpose, a line it cannot
//! read is no longer skipped ([`crate::versions_env`]).

use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, Write};
use std::path::Path;
use std::sync::LazyLock;

use regex::{Captures, Regex};

use crate::versions_env;

const USAGE: &str = "usage: og dev versions check-run --config FILE VERSIONS.yml [PROCESS...]\n       og dev versions expand --env versions.env < LINES\n       og dev versions coupled --env versions.env\n";

fn re(s: &str) -> Regex {
    Regex::new(s).expect("a valid built-in pattern")
}

static CONTAINER: LazyLock<Regex> =
    LazyLock::new(|| re(r"withName: '(\w+)' \{ container = '([^']+)' \}"));
static PIPELINE: LazyLock<Regex> =
    LazyLock::new(|| re(r"withName: '(\w+)' \{ ext\.pipeline_version = '([^']+)';"));
static BLOCK: LazyLock<Regex> = LazyLock::new(|| re(r#"^"([^"]+)":\s*$"#));
static ENTRY: LazyLock<Regex> = LazyLock::new(|| re(r"^\s+([\w.-]+):\s*(.*?)\s*$"));
static TAG_PREFIX: LazyLock<Regex> = LazyLock::new(|| re(r"^[^:@]+[:@]"));

/// The image (or pipeline release) of every process in the config text.
pub fn selectors(config: &str) -> BTreeMap<String, String> {
    let mut images = BTreeMap::new();
    for re in [&*CONTAINER, &*PIPELINE] {
        for c in re.captures_iter(config) {
            images.insert(c[1].to_owned(), c[2].to_owned());
        }
    }
    images
}

/// The first `key: value` under each `"PATH:PROC":` block of a
/// `software_versions.yml`.
pub fn first_entries(yml: &str) -> BTreeMap<String, (String, String)> {
    let mut first = BTreeMap::new();
    let mut proc: Option<String> = None;
    for line in yml.split_inclusive('\n') {
        let line = line.strip_suffix('\n').unwrap_or(line);
        if let Some(c) = BLOCK.captures(line) {
            proc = c[1].rsplit(':').next().map(str::to_owned);
            continue;
        }
        if let (Some(c), Some(p)) = (ENTRY.captures(line), &proc)
            && !p.is_empty()
            && !first.contains_key(p)
        {
            first.insert(p.clone(), (c[1].to_owned(), c[2].to_owned()));
        }
    }
    first
}

/// The check on parsed inputs. `expect` empty means every selector.
pub fn check_run(
    images: &BTreeMap<String, String>,
    first: &BTreeMap<String, (String, String)>,
    yml_name: &str,
    expect: &[String],
    out: &mut dyn Write,
) -> bool {
    let mut ok = true;
    if first.is_empty() {
        let _ = writeln!(out, "FAIL: {yml_name} has no process block");
        ok = false;
    }
    let wanted: BTreeSet<&str> = if expect.is_empty() {
        images.keys().map(String::as_str).collect()
    } else {
        expect.iter().map(String::as_str).collect()
    };
    let all: BTreeSet<&str> = wanted
        .iter()
        .copied()
        .chain(first.keys().map(String::as_str))
        .collect();
    for p in all {
        match (images.get(p), first.get(p)) {
            (None, _) => {
                let _ = writeln!(
                    out,
                    "FAIL {p:<20} has no selector in conf/containers.config"
                );
                ok = false;
            }
            (Some(_), None) => {
                if wanted.contains(p) {
                    let _ = writeln!(out, "FAIL {p:<20} missing from {yml_name}");
                    ok = false;
                }
            }
            (Some(image), Some((key, value))) => {
                let tag = TAG_PREFIX.replace(image, "");
                if *value == tag {
                    let _ = writeln!(out, "OK   {p:<20} {key}: {value}");
                } else {
                    let _ = writeln!(out, "FAIL {p:<20} {key}: {value}, want {tag} ({image})");
                    ok = false;
                }
            }
        }
    }
    ok
}

fn read(path: &str) -> Result<String, String> {
    std::fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))
}

/// `og dev versions check-run ARGS`.
pub fn cli_check_run(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let [flag, config, yml, expect @ ..] = args else {
        let _ = write!(err, "{USAGE}");
        return 2;
    };
    if flag != "--config" {
        let _ = write!(err, "{USAGE}");
        return 2;
    }
    let texts = read(config).and_then(|c| Ok((c, read(yml)?)));
    match texts {
        Ok((c, y)) => u8::from(!check_run(
            &selectors(&c),
            &first_entries(&y),
            yml,
            expect,
            out,
        )),
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            2
        }
    }
}

static VAR: LazyLock<Regex> = LazyLock::new(|| re(r"\$\{([A-Z][A-Z0-9_]*)\}|\$([A-Z][A-Z0-9_]*)"));

/// One URL with its variables filled and its trailing punctuation dropped.
pub fn expand(url: &str, vals: &BTreeMap<String, String>) -> String {
    let filled = VAR.replace_all(url, |c: &Captures| {
        let name = c.get(1).or_else(|| c.get(2)).map_or("", |m| m.as_str());
        vals.get(name).cloned().unwrap_or_else(|| c[0].to_owned())
    });
    let mut u = filled.into_owned();
    // Trailing punctuation, and the brace that closes ${VAR:-https://...}.
    while let Some(last) = u.chars().next_back() {
        let unbalanced = last == '}' && u.matches('}').count() > u.matches('{').count();
        if ".,;:]".contains(last) || unbalanced {
            u.pop();
        } else {
            break;
        }
    }
    u
}

/// `og dev versions expand --env FILE`.
pub fn cli_expand(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let [flag, env] = args else {
        let _ = write!(err, "{USAGE}");
        return 2;
    };
    if flag != "--env" {
        let _ = write!(err, "{USAGE}");
        return 2;
    }
    let pins = match std::fs::read_to_string(Path::new(env)) {
        Ok(text) => versions_env::parse(&text),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(Vec::new()),
        Err(e) => Err(format!("cannot read {env}: {e}")),
    };
    let vals: BTreeMap<String, String> = match pins {
        Ok(p) => p.into_iter().map(|p| (p.name, p.value)).collect(),
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            return 2;
        }
    };
    for (i, line) in std::io::stdin().lock().lines().enumerate() {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                let _ = writeln!(err, "ERROR: cannot read line {}: {e}", i + 1);
                return 2;
            }
        };
        let Some((url, loc)) = line.split_once('\t') else {
            let _ = writeln!(err, "ERROR: line {} has no tab: {line}", i + 1);
            return 2;
        };
        let _ = writeln!(out, "{}\t{loc}", expand(url, &vals));
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    const CONFIG: &str = "process {\n    withName: 'ROH' { container = 'quay.io/x/bcftools:1.22--h0_0' }\n    withName: 'VEP' { container = 'ensemblorg/ensembl-vep@sha256:abc' }\n    withName: 'PRS' { ext.pipeline_version = 'v2.1.0'; ext.x = 'y' }\n}\n";

    fn run(yml: &str, expect: &[&str]) -> (bool, String) {
        let expect: Vec<String> = expect.iter().map(|s| (*s).to_owned()).collect();
        let mut out = Vec::new();
        let ok = check_run(
            &selectors(CONFIG),
            &first_entries(yml),
            "v.yml",
            &expect,
            &mut out,
        );
        (ok, String::from_utf8(out).unwrap())
    }

    #[test]
    fn right_tags_pass() {
        let (ok, out) = run(
            "\"PGX:ROH\":\n    bcftools: 1.22--h0_0\n    other: 9\n\"VEP\":\n    ensemblvep: sha256:abc\n\"CLINICAL:PRS\":\n    pgsc_calc: v2.1.0\n",
            &[],
        );
        assert!(ok, "{out}");
        assert!(
            out.contains("OK   ROH                  bcftools: 1.22--h0_0\n"),
            "{out}"
        );
    }

    #[test]
    fn a_wrong_missing_or_unknown_process_fails() {
        let (ok, out) = run(
            "\"PGX:ROH\":\n    bcftools: 0.0.0-planted\n\"NEW_TOOL\":\n    tool: 1.0\n",
            &["ROH", "VEP"],
        );
        assert!(!ok);
        assert!(out.contains("FAIL ROH                  bcftools: 0.0.0-planted, want 1.22--h0_0 (quay.io/x/bcftools:1.22--h0_0)"), "{out}");
        assert!(
            out.contains("FAIL NEW_TOOL             has no selector"),
            "{out}"
        );
        assert!(
            out.contains("FAIL VEP                  missing from v.yml"),
            "{out}"
        );
        assert!(!out.contains("PRS"), "{out}");
        let (ok, out) = run("nothing: here\n", &["ROH"]);
        assert!(!ok);
        assert!(out.contains("FAIL: v.yml has no process block"), "{out}");
    }

    #[test]
    fn urls_expand_and_lose_trailing_punctuation() {
        let vals = BTreeMap::from([("REL".to_owned(), "1.2".to_owned())]);
        assert_eq!(
            expand("https://x.org/${REL}/a.gz", &vals),
            "https://x.org/1.2/a.gz"
        );
        assert_eq!(
            expand("https://x.org/$REL/a.gz).", &vals),
            "https://x.org/1.2/a.gz)"
        );
        assert_eq!(
            expand("https://x.org/${NOPE}/a.", &vals),
            "https://x.org/${NOPE}/a"
        );
        assert_eq!(
            expand("${U:-https://x.org/a}", &vals),
            "${U:-https://x.org/a}"
        );
        assert_eq!(expand("https://x.org/a}", &vals), "https://x.org/a");
        assert_eq!(expand("https://x.org/a];:", &vals), "https://x.org/a");
    }
}
