//! `og dev params-doc [--repo DIR]`: `docs/parameters.md` lists every
//! parameter of `nextflow_schema.json`, and nothing else.
//!
//! It replaces `tests/test_parameters_doc.py`. Three checks:
//!
//! 1. every parameter of the schema (the top-level `properties` and those
//!    of each group under `definitions`, or `$defs` when there is no
//!    `definitions`) has a row that starts ``| `--<name>` |`` in the page;
//! 2. every such row names a schema parameter, or one of [`DOC_ONLY`] that
//!    `main.nf` still reads (`params.<name>` or `'<name>'`), so a removed or
//!    renamed parameter cannot keep its row;
//! 3. no parameter has two rows.
//!
//! A schema with 20 parameters or fewer is refused too: that is a schema
//! read wrong, not a page that agrees with it.
//!
//! Exit 0 when the page agrees, 1 when a check fails (each failure names
//! the parameters), 2 when a file cannot be read or an argument is wrong.
//! `--repo` defaults to the current folder. An unknown option is refused
//! before any file is read.

use crate::py_os_error;
use pycompat::json::{Value, loads};
use std::collections::BTreeSet;
use std::fs;
use std::io::Write;
use std::path::Path;

/// Read by `main.nf` but not in the schema yet; the page writes their rows
/// by hand. Each keeps its row only while `main.nf` reads it.
pub const DOC_ONLY: [&str; 1] = ["pedigree"];

const USAGE: &str = "usage: og dev params-doc [--repo DIR]\n";

/// The parameter names of a parsed `nextflow_schema.json`.
pub fn schema_params(schema: &Value) -> BTreeSet<String> {
    let mut names = BTreeSet::new();
    let Value::Object(top) = schema else {
        return names;
    };
    let mut add = |v: Option<&Value>| {
        if let Some(Value::Object(props)) = v {
            names.extend(props.keys().cloned());
        }
    };
    add(top.get("properties"));
    let groups = top.get("definitions").or_else(|| top.get("$defs"));
    if let Some(Value::Object(groups)) = groups {
        for g in groups.values() {
            if let Value::Object(g) = g {
                add(g.get("properties"));
            }
        }
    }
    names
}

/// The parameter each row of the page names, in page order: a line that
/// starts ``| `--NAME` |`` with NAME of ASCII letters, digits and `_`.
pub fn rows(page: &str) -> Vec<&str> {
    page.lines()
        .filter_map(|l| {
            let rest = l.strip_prefix("| `--")?;
            let end = rest.find('`')?;
            let name = &rest[..end];
            let ok = !name.is_empty()
                && name.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_')
                && rest[end..].starts_with("` |");
            ok.then_some(name)
        })
        .collect()
}

/// `main.nf` reads `name`: `params.name` as a whole word, or `'name'`.
fn main_reads(main_nf: &str, name: &str) -> bool {
    let word = |c: char| c.is_alphanumeric() || c == '_';
    let dotted = format!("params.{name}");
    let whole = main_nf.match_indices(&dotted).any(|(i, _)| {
        let before = main_nf[..i].chars().next_back().is_none_or(|c| !word(c));
        let after = main_nf[i + dotted.len()..].chars().next();
        before && after.is_none_or(|c| !word(c))
    });
    whole || main_nf.contains(&format!("'{name}'"))
}

/// What the page gets wrong: (schema parameters with no row, rows that
/// name no parameter, parameters with more than one row), each sorted.
pub fn problems(
    params: &BTreeSet<String>,
    page: &str,
    main_nf: &str,
) -> (Vec<String>, Vec<String>, Vec<String>) {
    let all = rows(page);
    let named: BTreeSet<&str> = all.iter().copied().collect();
    let missing = params
        .iter()
        .filter(|p| !named.contains(p.as_str()))
        .cloned()
        .collect();
    let extra = named
        .iter()
        .filter(|r| !params.contains(**r))
        .filter(|r| !(DOC_ONLY.contains(r) && main_reads(main_nf, r)))
        .map(|r| (*r).to_owned())
        .collect();
    let twice = named
        .iter()
        .filter(|r| all.iter().filter(|x| x == r).count() > 1)
        .map(|r| (*r).to_owned())
        .collect();
    (missing, extra, twice)
}

fn read(path: &Path) -> Result<String, String> {
    fs::read_to_string(path).map_err(|e| py_os_error(&e, &path.to_string_lossy()))
}

/// `og dev params-doc ARGS`; returns the exit code.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let repo = match args {
        [] => ".".to_owned(),
        [flag, dir] if flag == "--repo" => dir.clone(),
        [flag] if flag.starts_with("--repo=") => flag["--repo=".len()..].to_owned(),
        [flag] if flag == "-h" || flag == "--help" => {
            let _ = write!(out, "{USAGE}");
            return 0;
        }
        _ => {
            let _ = write!(err, "{USAGE}");
            let _ = writeln!(
                err,
                "og dev params-doc: error: unrecognized arguments: {}",
                args.join(" ")
            );
            return 2;
        }
    };
    let repo = Path::new(&repo);
    let texts = ["nextflow_schema.json", "docs/parameters.md", "main.nf"]
        .map(|f| read(&repo.join(f)))
        .into_iter()
        .collect::<Result<Vec<_>, _>>();
    let [schema, page, main_nf] = match texts {
        Ok(t) => <[String; 3]>::try_from(t).expect("three files"),
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            return 2;
        }
    };
    let schema = match loads(&schema) {
        Ok(v) => v,
        Err(e) => {
            let _ = writeln!(err, "ERROR: nextflow_schema.json: {e}");
            return 2;
        }
    };
    let params = schema_params(&schema);
    let (missing, extra, twice) = problems(&params, &page, &main_nf);
    let mut failures = 0;
    let mut check = |desc: &str, ok: bool, detail: String| {
        let _ = if ok {
            writeln!(out, "[PASS] {desc}")
        } else {
            failures += 1;
            writeln!(out, "[FAIL] {desc} -- {detail}")
        };
    };
    check(
        "the schema has parameters to compare",
        params.len() > 20,
        format!("{} parameters", params.len()),
    );
    check(
        "every schema parameter has a row in docs/parameters.md",
        missing.is_empty(),
        format!("no row for: {}", missing.join(", ")),
    );
    check(
        "every row of docs/parameters.md is a schema parameter or one main.nf reads",
        extra.is_empty(),
        format!("not in the schema: {}", extra.join(", ")),
    );
    check(
        "each parameter has one row",
        twice.is_empty(),
        format!("more than one row: {}", twice.join(", ")),
    );
    let _ = writeln!(
        out,
        "RESULT: {}",
        if failures == 0 {
            "ALL PASS".to_owned()
        } else {
            format!("{failures} FAILED")
        }
    );
    u8::from(failures > 0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn repo() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
    }

    fn run(args: &[&str]) -> (u8, String, String) {
        let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = cli(&args, &mut out, &mut err);
        (
            code,
            String::from_utf8(out).unwrap(),
            String::from_utf8(err).unwrap(),
        )
    }

    /// A scratch copy of the three files the check reads, `page` edited.
    fn scratch(name: &str, edit: impl Fn(&str) -> String) -> PathBuf {
        let d = std::env::temp_dir().join(format!("og-params-doc-{}-{name}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(d.join("docs")).unwrap();
        for f in ["nextflow_schema.json", "main.nf"] {
            fs::copy(repo().join(f), d.join(f)).unwrap();
        }
        let page = fs::read_to_string(repo().join("docs/parameters.md")).unwrap();
        fs::write(d.join("docs/parameters.md"), edit(&page)).unwrap();
        d
    }

    fn real() -> (BTreeSet<String>, String, String) {
        let schema = fs::read_to_string(repo().join("nextflow_schema.json")).unwrap();
        (
            schema_params(&loads(&schema).unwrap()),
            fs::read_to_string(repo().join("docs/parameters.md")).unwrap(),
            fs::read_to_string(repo().join("main.nf")).unwrap(),
        )
    }

    /// The repository's own page agrees with its schema: this is the check
    /// Guard runs through `cargo test`.
    #[test]
    fn the_repository_page_lists_every_schema_parameter_once() {
        let (code, out, err) = run(&["--repo", repo().to_str().unwrap()]);
        assert_eq!(code, 0, "{out}{err}");
        assert!(out.ends_with("RESULT: ALL PASS\n"), "{out}");
        assert_eq!(out.matches("[PASS]").count(), 4, "{out}");
    }

    #[test]
    fn a_page_without_one_row_fails_and_names_that_parameter() {
        let (params, page, main_nf) = real();
        let victim = params.iter().next().unwrap().clone();
        let row = format!("| `--{victim}` |");
        let dropped: String = page
            .lines()
            .filter(|l| !l.starts_with(&row))
            .map(|l| format!("{l}\n"))
            .collect();
        assert_ne!(dropped, page, "the row was there to drop");
        let (missing, extra, twice) = problems(&params, &dropped, &main_nf);
        assert_eq!(
            (missing, extra, twice),
            (vec![victim.clone()], vec![], vec![])
        );
        let d = scratch("dropped", |p| {
            p.lines()
                .filter(|l| !l.starts_with(&row))
                .map(|l| format!("{l}\n"))
                .collect()
        });
        let (code, out, _) = run(&["--repo", d.to_str().unwrap()]);
        assert_eq!(code, 1, "{out}");
        assert!(out.contains(&format!("no row for: {victim}\n")), "{out}");
        fs::remove_dir_all(d).unwrap();
    }

    #[test]
    fn a_planted_row_for_a_parameter_the_schema_lacks_fails() {
        let d = scratch("planted", |p| {
            format!("{p}| `--no_such_parameter` | none | planted |\n")
        });
        let (code, out, _) = run(&["--repo", d.to_str().unwrap()]);
        assert_eq!(code, 1, "{out}");
        assert!(
            out.contains("not in the schema: no_such_parameter\n"),
            "{out}"
        );
        fs::remove_dir_all(d).unwrap();
    }

    #[test]
    fn a_row_given_twice_fails() {
        let (params, page, main_nf) = real();
        let first = params.iter().next().unwrap();
        let twice = format!("{page}| `--{first}` | none | again |\n");
        assert_eq!(problems(&params, &twice, &main_nf).2, vec![first.clone()]);
    }

    #[test]
    fn a_doc_only_row_stays_only_while_main_nf_reads_it() {
        let params: BTreeSet<String> = ["input".to_owned()].into();
        let page = "| `--input` | none | x |\n| `--pedigree` | none | y |\n";
        assert!(
            problems(&params, page, "def p = params.pedigree\n")
                .1
                .is_empty()
        );
        assert!(problems(&params, page, "x = ['pedigree']\n").1.is_empty());
        assert_eq!(
            problems(&params, page, "def p = params.pedigree_file\n").1,
            vec!["pedigree".to_owned()]
        );
    }

    #[test]
    fn groups_under_defs_count_when_there_is_no_definitions() {
        let s = loads(r#"{"properties": {"a": {}}, "$defs": {"g": {"properties": {"b": {}}}}}"#)
            .unwrap();
        assert_eq!(schema_params(&s), ["a".to_owned(), "b".to_owned()].into());
    }

    #[test]
    fn only_a_real_row_counts() {
        let page = "| `--ok_1` | x |\n|`--no` | x |\n| `--bad-name` | x |\n| `--x`|\nsee `--y` |\n";
        assert_eq!(rows(page), vec!["ok_1"]);
    }

    #[test]
    fn an_unknown_option_is_refused_before_any_file_is_read() {
        let (code, out, err) = run(&["--og-scaffold-probe"]);
        assert_eq!(code, 2);
        assert!(out.is_empty());
        assert!(
            err.ends_with(
                "og dev params-doc: error: unrecognized arguments: --og-scaffold-probe\n"
            ),
            "{err}"
        );
    }

    #[test]
    fn a_missing_file_exits_2() {
        let d = std::env::temp_dir().join(format!("og-params-doc-{}-empty", std::process::id()));
        fs::create_dir_all(&d).unwrap();
        let (code, _, err) = run(&["--repo", d.to_str().unwrap()]);
        assert_eq!(code, 2);
        assert!(err.starts_with("ERROR: [Errno 2]"), "{err}");
        fs::remove_dir_all(d).unwrap();
    }
}
