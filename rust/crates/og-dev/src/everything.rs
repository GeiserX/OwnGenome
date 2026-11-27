//! `og dev everything-check [--repo DIR]`: `assets/everything.tsv` against
//! every tool, switch, installer and optional image.
//!
//! `EVERYTHING=true` (scripts/run-all.sh) and `-profile everything`
//! (conf/everything.config) run every optional step from that one list, so
//! an opt-in that has no row there would be missed with no word. The check
//! reads, as text:
//!
//! - the names of `known_tools` in `main.nf`;
//! - the parameters of `nextflow.config`'s `params` block that default to
//!   `false`, and those whose comment says `opt-in`;
//! - the words of `SWITCHES=(...)` in `scripts/run-all.sh`;
//! - the flags of og's setup table (`og setup flags`,
//!   rust/crates/og-setup/src/components.rs), which must equal the first
//!   `case "${1:-}" in` of `scripts/setup.sh` while setup.sh still has one;
//! - the images of `versions.env` whose comment starts with `optional`;
//!
//! and fails when one of them has no row, when a row names an item that is
//! gone, when `conf/everything.config` sets other tools or parameters than
//! the rows say, or when `conf/test_all.config`'s tools miss one of them (the
//! stub run of every process would then miss it too).
//!
//! Exit 0 when everything agrees, 1 when a check fails (one `FAIL` line per
//! problem, naming the item), 2 when a file cannot be read or an argument is
//! wrong. `--repo` defaults to the current folder.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::Path;

const USAGE: &str = "usage: og dev everything-check [--repo DIR]\n";
const LIST: &str = "assets/everything.tsv";
const KINDS: [&str; 6] = ["tool", "switch", "setup", "input", "image", "excluded"];

/// One row of `assets/everything.tsv`.
struct Row {
    kind: String,
    name: String,
    value: String,
    nextflow: String,
    line: usize,
}

/// What the repository declares, read from its files.
#[derive(Default)]
struct Items {
    tools: BTreeSet<String>,
    nf_params: BTreeSet<String>,
    switches: BTreeSet<String>,
    setup_flags: BTreeSet<String>,
    /// The flags of setup.sh's case, while it has one.
    setup_case: Option<BTreeSet<String>>,
    images: BTreeSet<String>,
}

fn read(repo: &Path, rel: &str) -> Result<String, String> {
    fs::read_to_string(repo.join(rel)).map_err(|e| format!("cannot read {rel}: {e}"))
}

fn read_rows(text: &str) -> Result<Vec<Row>, String> {
    let mut rows = Vec::new();
    let mut header = false;
    for (i, line) in text.lines().enumerate() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if !header {
            if c != ["kind", "name", "value", "nextflow", "reason"] {
                return Err(format!(
                    "{LIST} line {}: the header must be kind, name, value, nextflow, reason",
                    i + 1
                ));
            }
            header = true;
            continue;
        }
        let [kind, name, value, nextflow, reason] = c.as_slice() else {
            return Err(format!(
                "{LIST} line {}: {} fields, want 5: kind, name, value, nextflow, reason",
                i + 1,
                c.len()
            ));
        };
        if [kind, name, value, nextflow, reason]
            .iter()
            .any(|f| f.trim().is_empty())
        {
            return Err(format!("{LIST} line {}: an empty field", i + 1));
        }
        rows.push(Row {
            kind: (*kind).to_owned(),
            name: (*name).to_owned(),
            value: (*value).to_owned(),
            nextflow: (*nextflow).to_owned(),
            line: i + 1,
        });
    }
    if rows.is_empty() {
        return Err(format!("{LIST}: no rows"));
    }
    Ok(rows)
}

/// The quoted words of `def known_tools = [ ... ]` in `main.nf`.
fn known_tools(main_nf: &str) -> Result<BTreeSet<String>, String> {
    let start = main_nf
        .find("def known_tools = [")
        .ok_or("main.nf: no `def known_tools = [`")?;
    let rest = &main_nf[start + "def known_tools = [".len()..];
    let end = rest.find(']').ok_or("main.nf: known_tools has no `]`")?;
    let names: BTreeSet<String> = rest[..end]
        .split('\'')
        .skip(1)
        .step_by(2)
        .map(str::to_owned)
        .collect();
    if names.len() < 10 {
        return Err(format!(
            "main.nf: known_tools read as {} names: that is the list read wrong",
            names.len()
        ));
    }
    Ok(names)
}

/// A line split into its code and its `//` comment, outside quotes.
fn split_comment(line: &str) -> (&str, &str) {
    let mut quote = None;
    let bytes = line.as_bytes();
    for (i, &b) in bytes.iter().enumerate() {
        match (quote, b) {
            (None, b'\'' | b'"') => quote = Some(b),
            (Some(q), _) if b == q => quote = None,
            (None, b'/') if bytes.get(i + 1) == Some(&b'/') => {
                return (&line[..i], &line[i + 2..]);
            }
            _ => {}
        }
    }
    (line, "")
}

/// `name = value` of one config line, with its comment.
fn assignment(line: &str) -> Option<(&str, &str, &str)> {
    let (code, comment) = split_comment(line);
    let (name, value) = code.split_once('=')?;
    let name = name.trim();
    let ok = !name.is_empty()
        && name
            .chars()
            .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_');
    ok.then(|| (name, value.trim(), comment))
}

/// The lines of the first top-level `params {` block of a config.
fn params_block(config: &str) -> Option<Vec<&str>> {
    let mut lines = config.lines();
    lines.find(|l| l.trim_end() == "params {")?;
    Some(lines.take_while(|l| !l.starts_with('}')).collect())
}

/// The parameters of `nextflow.config` that default to `false` or whose
/// comment says `opt-in`.
fn optional_params(config: &str) -> Result<BTreeSet<String>, String> {
    let block = params_block(config).ok_or("nextflow.config: no `params {` block")?;
    let mut out = BTreeSet::new();
    let mut seen = 0;
    for line in block {
        if let Some((name, value, comment)) = assignment(line) {
            seen += 1;
            if value == "false" || comment.contains("opt-in") {
                out.insert(name.to_owned());
            }
        }
    }
    if seen < 20 {
        return Err(format!(
            "nextflow.config: the params block read as {seen} parameters: that is the block read wrong"
        ));
    }
    Ok(out)
}

/// The words of `SWITCHES=( ... )` in `run-all.sh`.
fn switches(run_all: &str) -> Result<BTreeSet<String>, String> {
    let start = run_all
        .find("SWITCHES=(")
        .ok_or("scripts/run-all.sh: no `SWITCHES=(`")?;
    let rest = &run_all[start + "SWITCHES=(".len()..];
    let end = rest
        .find(')')
        .ok_or("scripts/run-all.sh: SWITCHES has no `)`")?;
    let words: BTreeSet<String> = rest[..end].split_whitespace().map(str::to_owned).collect();
    if words.len() < 10 {
        return Err(format!(
            "scripts/run-all.sh: SWITCHES read as {} names: that is the list read wrong",
            words.len()
        ));
    }
    Ok(words)
}

/// The flags of the first `case "${1:-}" in` of `setup.sh`: every `--flag`
/// of each pattern line, up to its `esac`; `None` once setup.sh has no such
/// case (a wrapper that hands its arguments to og).
fn setup_flags(setup: &str) -> Result<Option<BTreeSet<String>>, String> {
    let mut lines = setup.lines();
    if !lines.any(|l| l.trim() == r#"case "${1:-}" in"#) {
        return Ok(None);
    }
    let mut out = BTreeSet::new();
    for line in lines.take_while(|l| l.trim() != "esac") {
        let t = line.trim();
        let Some(pattern) = t.strip_suffix(')') else {
            continue;
        };
        if pattern.starts_with("--") {
            out.extend(pattern.split('|').map(str::to_owned));
        }
    }
    if out.len() < 5 {
        return Err(format!(
            "scripts/setup.sh: its flags read as {}: that is the case read wrong",
            out.len()
        ));
    }
    Ok(Some(out))
}

/// The images of `versions.env` whose comment starts with `optional`.
fn optional_images(env: &str) -> BTreeSet<String> {
    env.lines()
        .filter_map(|l| {
            let (name, rest) = l.split_once('=')?;
            let ok = name.ends_with("_IMAGE")
                && name
                    .chars()
                    .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit() || c == '_');
            let (_, comment) = rest.split_once("  # ")?;
            (ok && comment.trim_start().starts_with("optional")).then(|| name.to_owned())
        })
        .collect()
}

/// The comma-separated names of a config's `tools = '...'`.
fn config_tools(config: &str, rel: &str) -> Result<BTreeSet<String>, String> {
    let block = params_block(config).ok_or(format!("{rel}: no `params {{` block"))?;
    let value = block
        .iter()
        .find_map(|l| assignment(l).filter(|(n, _, _)| *n == "tools"))
        .map(|(_, v, _)| v)
        .ok_or(format!("{rel}: no `tools =` in its params block"))?;
    Ok(value
        .trim_matches(|c| c == '\'' || c == '"')
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_owned)
        .collect())
}

fn read_items(repo: &Path) -> Result<Items, String> {
    let env = read(repo, "versions.env")?;
    let images = optional_images(&env);
    if images.len() < 5 {
        return Err(format!(
            "versions.env: {} optional images: that is the file read wrong",
            images.len()
        ));
    }
    Ok(Items {
        tools: known_tools(&read(repo, "main.nf")?)?,
        nf_params: optional_params(&read(repo, "nextflow.config")?)?,
        switches: switches(&read(repo, "scripts/run-all.sh")?)?,
        setup_flags: og_setup::components::flags().map(str::to_owned).collect(),
        setup_case: setup_flags(&read(repo, "scripts/setup.sh")?)?,
        images,
    })
}

/// Does the repository still declare the item a row names?
fn row_item_exists(row: &Row, items: &Items) -> bool {
    let n = row.name.as_str();
    let param = n.strip_prefix("--");
    match row.kind.as_str() {
        "tool" => items.tools.contains(n),
        "switch" => items.switches.contains(n),
        "setup" => items.setup_flags.contains(n),
        "image" => items.images.contains(n),
        _ if n.ends_with("_IMAGE") => items.images.contains(n),
        _ if param.is_some() => {
            items.setup_flags.contains(n) || param.is_some_and(|p| items.nf_params.contains(p))
        }
        _ if n.chars().any(|c| c.is_ascii_uppercase()) => items.switches.contains(n),
        "input" => true,
        _ => items.tools.contains(n),
    }
}

/// Every problem, one line each; empty when the list agrees with the repo.
fn problems(repo: &Path) -> Result<(Vec<String>, usize), String> {
    let rows = read_rows(&read(repo, LIST)?)?;
    let items = read_items(repo)?;
    let mut bad = Vec::new();
    let mut by_name: BTreeMap<&str, &Row> = BTreeMap::new();
    for r in &rows {
        if !KINDS.contains(&r.kind.as_str()) {
            bad.push(format!(
                "{LIST} line {}: kind '{}' is none of {}",
                r.line,
                r.kind,
                KINDS.join(", ")
            ));
        }
        if by_name.insert(&r.name, r).is_some() {
            bad.push(format!("{LIST} line {}: '{}' has two rows", r.line, r.name));
        }
    }
    for r in &rows {
        if KINDS.contains(&r.kind.as_str()) && !row_item_exists(r, &items) {
            bad.push(format!(
                "{LIST} line {}: {} row '{}' names nothing the repository declares (renamed or removed?)",
                r.line, r.kind, r.name
            ));
        }
        if matches!(r.kind.as_str(), "setup" | "image") {
            for target in r.value.split(',') {
                if !by_name.contains_key(target) {
                    bad.push(format!(
                        "{LIST} line {}: {} row '{}' serves '{target}', which has no row",
                        r.line, r.kind, r.name
                    ));
                }
            }
        }
        let nf_ok = if r.kind == "switch" {
            r.nextflow == "-"
                || r.nextflow
                    .strip_prefix("--")
                    .and_then(|s| s.split_once(' '))
                    .is_some_and(|(p, v)| !p.is_empty() && !v.is_empty() && !v.contains(' '))
        } else {
            r.nextflow == "-"
        };
        if !nf_ok {
            bad.push(format!(
                "{LIST} line {}: nextflow '{}' of '{}': a switch row takes '--param value' or '-', any other row '-'",
                r.line, r.nextflow, r.name
            ));
        }
    }
    // Every item the repository declares has a row.
    let has = |name: &str| by_name.contains_key(name);
    let nf_named: BTreeSet<&str> = rows
        .iter()
        .filter_map(|r| r.nextflow.strip_prefix("--")?.split(' ').next())
        .collect();
    for t in &items.tools {
        if !has(t) {
            bad.push(format!(
                "main.nf: tool '{t}' of known_tools has no row in {LIST}"
            ));
        }
    }
    for p in &items.nf_params {
        if !has(&format!("--{p}")) && !nf_named.contains(p.as_str()) {
            bad.push(format!(
                "nextflow.config: params.{p} (false by default or opt-in) has no row in {LIST}"
            ));
        }
    }
    for s in &items.switches {
        if !has(s) {
            bad.push(format!(
                "scripts/run-all.sh: switch {s} of SWITCHES has no row in {LIST}"
            ));
        }
    }
    // og's table is the list of setup's flags; while setup.sh parses its
    // own flags, the two must agree.
    if let Some(case) = &items.setup_case {
        for f in case.difference(&items.setup_flags) {
            bad.push(format!(
                "scripts/setup.sh: flag {f} is not in og's setup table (rust/crates/og-setup/src/components.rs)"
            ));
        }
        for f in items.setup_flags.difference(case) {
            bad.push(format!(
                "rust/crates/og-setup/src/components.rs: flag {f} is not a flag of scripts/setup.sh"
            ));
        }
    }
    for f in &items.setup_flags {
        if !has(f) {
            bad.push(format!("og setup flags: flag {f} has no row in {LIST}"));
        }
    }
    for i in &items.images {
        if !has(i) {
            bad.push(format!(
                "versions.env: optional image {i} has no row in {LIST}"
            ));
        }
    }
    // conf/everything.config sets exactly the tool rows and the switch rows'
    // Nextflow parameters.
    let want_tools: BTreeSet<String> = rows
        .iter()
        .filter(|r| r.kind == "tool")
        .map(|r| r.name.clone())
        .collect();
    let ev_rel = "conf/everything.config";
    let ev = read(repo, ev_rel)?;
    let ev_tools = config_tools(&ev, ev_rel)?;
    for t in want_tools.difference(&ev_tools) {
        bad.push(format!("{ev_rel}: tools lacks '{t}', a tool row of {LIST}"));
    }
    for t in ev_tools.difference(&want_tools) {
        bad.push(format!(
            "{ev_rel}: tools names '{t}', which is no tool row of {LIST}"
        ));
    }
    let want_params: BTreeMap<&str, &str> = rows
        .iter()
        .filter(|r| r.kind == "switch")
        .filter_map(|r| r.nextflow.strip_prefix("--")?.split_once(' '))
        .collect();
    let mut set: BTreeMap<String, String> = BTreeMap::new();
    for line in params_block(&ev).unwrap_or_default() {
        if let Some((name, value, _)) = assignment(line)
            && name != "tools"
        {
            set.insert(name.to_owned(), value.trim_matches('\'').to_owned());
        }
    }
    for (p, v) in &want_params {
        match set.get(*p) {
            Some(got) if got == v => {}
            Some(got) => bad.push(format!(
                "{ev_rel}: {p} = {got}, but its switch row says {v}"
            )),
            None => bad.push(format!(
                "{ev_rel}: does not set {p} = {v} of its switch row"
            )),
        }
    }
    for p in set.keys() {
        if !want_params.contains_key(p.as_str()) {
            bad.push(format!(
                "{ev_rel}: sets {p}, which no switch row of {LIST} names"
            ));
        }
    }
    // The stub run of every process covers every tool of the everything profile.
    let ta_rel = "conf/test_all.config";
    let ta_tools = config_tools(&read(repo, ta_rel)?, ta_rel)?;
    for t in want_tools.difference(&ta_tools) {
        bad.push(format!(
            "{ta_rel}: tools lacks '{t}', so the stub run of every tool misses it"
        ));
    }
    Ok((bad, rows.len()))
}

/// `og dev everything-check ARGS`: 0 when the list covers every optional
/// item, 1 when it does not, 2 on a usage or read error.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let mut repo = std::path::PathBuf::from(".");
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match (a.as_str(), it.clone().next()) {
            ("--repo", Some(v)) => {
                repo = v.into();
                it.next();
            }
            ("-h" | "--help", _) => {
                let _ = write!(out, "{USAGE}");
                return 0;
            }
            _ => {
                let _ = write!(err, "{USAGE}");
                let _ = writeln!(
                    err,
                    "og dev everything-check: error: unexpected argument '{a}'"
                );
                return 2;
            }
        }
    }
    match problems(&repo) {
        Ok((bad, n)) if bad.is_empty() => {
            let _ = writeln!(
                out,
                "ok      {LIST}: {n} rows cover every tool, switch, setup flag and optional image; conf/everything.config and conf/test_all.config agree"
            );
            0
        }
        Ok((bad, _)) => {
            for b in &bad {
                let _ = writeln!(out, "FAIL    {b}");
            }
            let _ = writeln!(
                out,
                "{} problem(s): give each new opt-in a row in {LIST} (kind tool, switch, setup, input, image or excluded with the reason)",
                bad.len()
            );
            1
        }
        Err(e) => {
            let _ = writeln!(err, "og dev everything-check: {e}");
            2
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    const FILES: [&str; 8] = [
        LIST,
        "main.nf",
        "nextflow.config",
        "scripts/run-all.sh",
        "scripts/setup.sh",
        "versions.env",
        "conf/everything.config",
        "conf/test_all.config",
    ];

    fn repo_root() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
    }

    /// A copy of the files the check reads, in a folder of its own.
    fn copy(name: &str) -> PathBuf {
        let d =
            std::env::temp_dir().join(format!("og-dev-everything-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        for f in FILES {
            let to = d.join(f);
            fs::create_dir_all(to.parent().unwrap()).unwrap();
            fs::copy(repo_root().join(f), &to).unwrap();
        }
        d
    }

    fn edit(d: &Path, rel: &str, from: &str, to: &str) {
        let p = d.join(rel);
        let text = fs::read_to_string(&p).unwrap();
        assert!(text.contains(from), "{rel} has no {from:?}");
        fs::write(&p, text.replacen(from, to, 1)).unwrap();
    }

    fn run(d: &Path) -> (u8, String) {
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let args = ["--repo".to_owned(), d.display().to_string()];
        let code = cli(&args, &mut out, &mut err);
        let mut text = String::from_utf8(out).unwrap();
        text.push_str(&String::from_utf8(err).unwrap());
        (code, text)
    }

    fn assert_fails(d: &Path, want: &str) {
        let (code, text) = run(d);
        let _ = fs::remove_dir_all(d);
        assert_eq!(code, 1, "{text}");
        assert!(text.contains(want), "no {want:?} in:\n{text}");
    }

    #[test]
    fn the_repository_agrees() {
        let d = copy("ok");
        let (code, text) = run(&d);
        let _ = fs::remove_dir_all(&d);
        assert_eq!(code, 0, "{text}");
        assert!(text.starts_with("ok "), "{text}");
    }

    #[test]
    fn a_new_tool_in_known_tools_fails() {
        let d = copy("tool");
        edit(
            &d,
            "main.nf",
            "'locityper',\n",
            "'locityper',\n        'newtool',\n",
        );
        assert_fails(&d, "main.nf: tool 'newtool' of known_tools has no row");
    }

    #[test]
    fn a_new_false_param_fails() {
        let d = copy("param");
        edit(
            &d,
            "nextflow.config",
            "    kir                = false",
            "    foo                = false      // a new switch\n    kir                = false",
        );
        assert_fails(
            &d,
            "nextflow.config: params.foo (false by default or opt-in) has no row",
        );
    }

    #[test]
    fn a_new_opt_in_param_fails() {
        let d = copy("optin");
        edit(
            &d,
            "nextflow.config",
            "    pedigree           = null",
            "    bar                = null       // opt-in: a new input\n    pedigree           = null",
        );
        assert_fails(
            &d,
            "nextflow.config: params.bar (false by default or opt-in) has no row",
        );
    }

    #[test]
    fn a_new_setup_flag_fails() {
        let d = copy("setup");
        edit(
            &d,
            "scripts/setup.sh",
            "|--archaic)",
            "|--archaic|--new-data)",
        );
        assert_fails(
            &d,
            "scripts/setup.sh: flag --new-data is not in og's setup table",
        );
    }

    #[test]
    fn a_flag_gone_from_setup_sh_fails() {
        let d = copy("setup-gone");
        edit(&d, "scripts/setup.sh", "|--archaic)", ")");
        assert_fails(
            &d,
            "rust/crates/og-setup/src/components.rs: flag --archaic is not a flag of scripts/setup.sh",
        );
    }

    #[test]
    fn a_table_flag_with_no_row_fails() {
        let d = copy("setup-row");
        let p = d.join(LIST);
        let text = fs::read_to_string(&p).unwrap();
        let kept: Vec<&str> = text
            .lines()
            .filter(|l| l.split('\t').nth(1) != Some("--kir-data"))
            .collect();
        assert_ne!(
            kept.len(),
            text.lines().count(),
            "{LIST} has no --kir-data row"
        );
        fs::write(&p, kept.join("\n") + "\n").unwrap();
        assert_fails(&d, "og setup flags: flag --kir-data has no row in");
    }

    #[test]
    fn a_new_optional_image_fails() {
        let d = copy("image");
        edit(
            &d,
            "versions.env",
            "TIDDIT_IMAGE=",
            "NEWTOOL_IMAGE=\"quay.io/biocontainers/newtool:1.0--0\"  # optional: step 99\nTIDDIT_IMAGE=",
        );
        assert_fails(&d, "versions.env: optional image NEWTOOL_IMAGE has no row");
    }

    #[test]
    fn a_new_run_all_switch_fails() {
        let d = copy("switch");
        edit(
            &d,
            "scripts/run-all.sh",
            "SWITCHES=(",
            "SWITCHES=(NEWSWITCH ",
        );
        assert_fails(
            &d,
            "scripts/run-all.sh: switch NEWSWITCH of SWITCHES has no row",
        );
    }

    #[test]
    fn a_stale_row_fails() {
        let d = copy("stale");
        edit(&d, "main.nf", "        'mocha',\n", "");
        assert_fails(&d, "tool row 'mocha' names nothing the repository declares");
    }

    #[test]
    fn the_profile_and_the_stub_run_must_follow_the_rows() {
        let d = copy("profile");
        edit(&d, "conf/everything.config", "hla_qc", "hla_qc_off");
        assert_fails(&d, "conf/everything.config: does not set hla_qc = true");
        let d = copy("profile-tool");
        edit(&d, "conf/everything.config", ",locityper", "");
        assert_fails(&d, "conf/everything.config: tools lacks 'locityper'");
        let d = copy("test-all");
        edit(&d, "conf/test_all.config", ",locityper", "");
        assert_fails(&d, "conf/test_all.config: tools lacks 'locityper'");
    }

    #[test]
    fn a_bad_argument_is_refused() {
        let (mut out, mut err) = (Vec::new(), Vec::new());
        assert_eq!(
            cli(&["--og-scaffold-probe".to_owned()], &mut out, &mut err),
            2
        );
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let missing = ["--repo".to_owned(), "/nonexistent-og-repo".to_owned()];
        assert_eq!(cli(&missing, &mut out, &mut err), 2);
    }
}
