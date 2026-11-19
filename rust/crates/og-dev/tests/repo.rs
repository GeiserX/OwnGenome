//! `og dev contract` and `og dev report-only` on this checkout, and on
//! copies of it with one planted break each. They replace
//! `tests/test_output_contract.py` and `tests/test_contract_report_only.py`.
//!
//! 1. The static contract check passes on this checkout; each planted break
//!    of the code (a renamed published file, a removed param, a removed
//!    switch, a new `--tools` name with no row, a new published file with
//!    no row, a schema default that is not `nextflow.config`'s, also one of
//!    another JSON type) or of the contract (a file row in the folder a
//!    `saveAs` drops it from, a `nextflow` row a step script writes too) is
//!    refused, the stable ones as BREAKING.
//! 2. The run check passes on a results folder planted from the contract's
//!    own e2e rows, for both legs; a dropped column, VCF field or JSON key,
//!    a renamed file and an unreadable gzip are each refused and named; an
//!    experimental row is never enforced; an empty folder fails.
//! 3. A contract with a row listed twice or a bad stability is refused.
//! 4. Against a release's contract (`--baseline`, this checkout's standing
//!    in for it): a stable output renamed in the code and in its row passes
//!    the static check and fails here; a removed stable param or tool, a
//!    stable row made experimental or internal, another format, a dropped
//!    required column and a dropped leg each fail; a new row, a new leg, a
//!    new required column, an `e2e` or `process` change and any change to
//!    an experimental row pass; a baseline with no stable row fails.
//! 5. The report-only check passes on this checkout and names both cDNA
//!    checks; dropping a process from the report-only list fails naming
//!    it; a second consumer makes a process a feeder; a stale KEEP_STOP
//!    name fails; a `//` inside a string and a nested brace in the strategy
//!    closure do not hide a block.

use std::fs;
use std::path::{Path, PathBuf};

const SAMPLE: &str = "S1";

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn s(p: &Path) -> String {
    p.to_str().expect("a UTF-8 path").to_owned()
}

/// A fresh folder under the system temp folder, removed on drop.
struct Tmp(PathBuf);

impl Tmp {
    fn new(tag: &str) -> Tmp {
        use std::sync::atomic::{AtomicUsize, Ordering};
        static N: AtomicUsize = AtomicUsize::new(0);
        let p = std::env::temp_dir().join(format!(
            "og-dev-{tag}-{}-{}",
            std::process::id(),
            N.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = fs::remove_dir_all(&p);
        fs::create_dir_all(&p).unwrap();
        Tmp(p)
    }
}

impl Drop for Tmp {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn copy_dir(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for e in fs::read_dir(from).unwrap() {
        let e = e.unwrap();
        let dst = to.join(e.file_name());
        if e.file_type().unwrap().is_dir() {
            copy_dir(&e.path(), &dst);
        } else {
            fs::copy(e.path(), dst).unwrap();
        }
    }
}

fn edit(path: &Path, old: &str, new: &str) {
    let src = fs::read_to_string(path).unwrap();
    assert!(src.contains(old), "{old:?} not in {}", path.display());
    fs::write(path, src.replacen(old, new, 1)).unwrap();
}

/// Runs a verb; returns its exit code and stdout plus stderr.
fn run(
    verb: fn(&[String], &mut dyn std::io::Write, &mut dyn std::io::Write) -> u8,
    args: &[&str],
) -> (u8, String) {
    let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let rc = verb(&args, &mut out, &mut err);
    out.extend(err);
    (rc, String::from_utf8(out).unwrap())
}

fn contract(args: &[&str]) -> (u8, String) {
    run(og_dev::contract::cli, args)
}

fn report_only(dir: &Path) -> (u8, String) {
    run(og_dev::report_only::cli, &["--repo", &s(dir)])
}

fn contract_path() -> PathBuf {
    repo().join("assets/output_contract.tsv")
}

/// The files the static check reads, copied, so a control can change them.
fn code_copy(tmp: &Tmp) -> PathBuf {
    let r = repo();
    let dst = tmp.0.join("repo");
    copy_dir(&r.join("modules"), &dst.join("modules"));
    fs::create_dir_all(dst.join("scripts")).unwrap();
    for f in ["main.nf", "nextflow_schema.json", "nextflow.config"] {
        fs::copy(r.join(f), dst.join(f)).unwrap();
    }
    for e in fs::read_dir(r.join("scripts")).unwrap() {
        let e = e.unwrap();
        if e.file_name().to_string_lossy().ends_with(".sh") {
            fs::copy(e.path(), dst.join("scripts").join(e.file_name())).unwrap();
        }
    }
    dst
}

/// The contract with the cells of the file row called `name` changed.
fn contract_copy(tmp: &Tmp, name: &str, change: impl Fn(&mut Vec<String>)) -> PathBuf {
    let path = tmp.0.join("contract.tsv");
    let mut hit = false;
    let mut out = String::new();
    for line in fs::read_to_string(contract_path()).unwrap().lines() {
        let mut c: Vec<String> = line.split('\t').map(str::to_owned).collect();
        if c.len() > 1 && c[0] == "file" && c[1] == name {
            change(&mut c);
            hit = true;
        }
        out.push_str(&c.join("\t"));
        out.push('\n');
    }
    assert!(hit, "no file row {name}");
    fs::write(&path, out).unwrap();
    path
}

fn has(out: &str, re: &str) -> bool {
    regex::Regex::new(re).unwrap().is_match(out)
}

#[test]
fn static_check_passes_on_this_checkout() {
    let (rc, out) = contract(&["--static", "--repo", &s(&repo())]);
    assert!(rc == 0 && out.contains(" 0 problems"), "{out}");
}

#[test]
fn static_controls_are_refused() {
    type Plant = fn(&Path);
    let cases: [(&str, Plant, &[&str]); 7] = [
        (
            "a renamed published file (ROH's summary)",
            |d| {
                edit(
                    &d.join("modules/local/roh/main.nf"),
                    "path(\"*_roh_summary.txt\")",
                    "path(\"*_roh_segments.txt\")",
                )
            },
            &[
                r"BREAKING: the stable file \{sample\}/roh/\{sample\}_roh_summary\.txt",
                r"writes \*_roh_segments\.txt and has no contract row",
            ],
        ),
        (
            "a removed param (--intervals)",
            |d| drop_param(d, "intervals"),
            &[r"BREAKING: the stable param intervals"],
        ),
        (
            "a removed run-all.sh switch (SKIP_TRIM)",
            |d| edit(&d.join("scripts/run-all.sh"), " SKIP_TRIM ", " "),
            &[r"BREAKING: the stable switch SKIP_TRIM"],
        ),
        (
            "a new --tools name with no row",
            |d| {
                edit(
                    &d.join("main.nf"),
                    "'pharmcat', 'cpic',",
                    "'pharmcat', 'newtool', 'cpic',",
                )
            },
            &[r"tool newtool of main.nf known_tools has no row"],
        ),
        (
            "a new published file with no row",
            |d| {
                edit(
                    &d.join("modules/local/roh/main.nf"),
                    "path(\"*_roh.txt\"),",
                    "path(\"*_roh.txt\"),\n    tuple val(meta), path(\"*_roh_extra.tsv\"), emit: extra\n",
                )
            },
            &[r"ROH .* writes \*_roh_extra\.tsv and has no contract row"],
        ),
        (
            "a schema default that is not nextflow.config's (traits left out of --tools)",
            |d| {
                edit(
                    &d.join("nextflow_schema.json"),
                    "mito_variants,traits,html_report",
                    "mito_variants,html_report",
                )
            },
            &[r"param tools: nextflow_schema\.json's default .* is not nextflow\.config's"],
        ),
        (
            // Red on the Python, which compared with False == 0.
            "a schema default of another JSON type (0 for nextflow.config's false)",
            |d| set_default(d, "allow_unfiltered", "0"),
            &[
                r"param allow_unfiltered: nextflow_schema\.json's default 0 is not nextflow\.config's False",
            ],
        ),
    ];
    for (desc, plant, want) in cases {
        let tmp = Tmp::new("static");
        let d = code_copy(&tmp);
        plant(&d);
        let (rc, out) = contract(&[
            "--static",
            "--repo",
            &s(&d),
            "--contract",
            &s(&contract_path()),
        ]);
        let missing: Vec<&&str> = want.iter().filter(|w| !has(&out, w)).collect();
        assert!(
            rc == 1 && missing.is_empty(),
            "{desc}: rc={rc}, missing {missing:?}: {out}"
        );
    }
}

/// Removes param `name` from every group of the copy's schema.
fn drop_param(d: &Path, name: &str) {
    let path = d.join("nextflow_schema.json");
    let mut v = pycompat::json::loads(&fs::read_to_string(&path).unwrap()).unwrap();
    if let pycompat::json::Value::Object(top) = &mut v
        && let Some(pycompat::json::Value::Object(groups)) = top.get_mut("definitions")
    {
        for g in groups.values_mut() {
            if let pycompat::json::Value::Object(g) = g
                && let Some(pycompat::json::Value::Object(props)) = g.get_mut("properties")
            {
                props.shift_remove(name);
            }
        }
    }
    let text = pycompat::json::dumps(&v, &pycompat::json::Dump::default()).unwrap();
    fs::write(path, text).unwrap();
}

/// Sets the schema default of param `name` to the JSON text `value`.
fn set_default(d: &Path, name: &str, value: &str) {
    let path = d.join("nextflow_schema.json");
    let mut v = pycompat::json::loads(&fs::read_to_string(&path).unwrap()).unwrap();
    let new = pycompat::json::loads(value).unwrap();
    let mut hit = false;
    if let pycompat::json::Value::Object(top) = &mut v
        && let Some(pycompat::json::Value::Object(groups)) = top.get_mut("definitions")
    {
        for g in groups.values_mut() {
            if let pycompat::json::Value::Object(g) = g
                && let Some(pycompat::json::Value::Object(props)) = g.get_mut("properties")
                && let Some(pycompat::json::Value::Object(p)) = props.get_mut(name)
            {
                p.insert("default".to_owned(), new.clone());
                hit = true;
            }
        }
    }
    assert!(hit, "no param {name}");
    let text = pycompat::json::dumps(&v, &pycompat::json::Dump::default()).unwrap();
    fs::write(path, text).unwrap();
}

#[test]
fn contract_controls_are_refused() {
    type Change = fn(&mut Vec<String>);
    let cases: [(&str, &str, Change, &str); 3] = [
        (
            "the read count in aligned/, whose publishDir drops it (it is published to qc/)",
            "{sample}/qc/{sample}_read_count.tsv",
            |c| c[1] = "{sample}/aligned/{sample}_read_count.tsv".to_owned(),
            r"BREAKING: the stable file \{sample\}/aligned/\{sample\}_read_count\.tsv .* is not what ALIGN_MINIMAP2 publishes",
        ),
        (
            "legs nextflow on a file step 49 writes under the same name",
            "{sample}/chip/{sample}_chip.tsv",
            |c| c[4] = "nextflow".to_owned(),
            r"\{sample\}/chip/\{sample\}_chip\.tsv .* has legs nextflow, but scripts/49-chip\.sh writes it",
        ),
        (
            "legs nextflow on a file step 50 writes on a line of seven shell variables",
            "{sample}/mei/{sample}_megane.log",
            |c| c[4] = "nextflow".to_owned(),
            r"\{sample\}/mei/\{sample\}_megane\.log .* has legs nextflow, but scripts/50-mei\.sh writes it",
        ),
    ];
    for (desc, name, change, want) in cases {
        let tmp = Tmp::new("contract");
        let c = contract_copy(&tmp, name, change);
        let (rc, out) = contract(&["--static", "--repo", &s(&repo()), "--contract", &s(&c)]);
        assert!(rc == 1 && has(&out, want), "{desc}: rc={rc}: {out}");
    }
}

// ---------------------------------------------------------------- the release baseline

/// The contract with `change` applied to its lines, the cells of each split
/// on tabs; `tag` names the copy.
fn contract_edit(tmp: &Tmp, tag: &str, change: impl Fn(&mut Vec<Vec<String>>)) -> PathBuf {
    let path = tmp.0.join(format!("{tag}.tsv"));
    let mut lines: Vec<Vec<String>> = fs::read_to_string(contract_path())
        .unwrap()
        .lines()
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect();
    let before = lines.clone();
    change(&mut lines);
    assert!(lines != before, "{tag}: the change changed nothing");
    let text: Vec<String> = lines.iter().map(|c| c.join("\t")).collect();
    fs::write(&path, text.join("\n") + "\n").unwrap();
    path
}

/// The line of `kind` row `name`.
fn row_mut<'a>(lines: &'a mut [Vec<String>], kind: &str, name: &str) -> &'a mut Vec<String> {
    lines
        .iter_mut()
        .find(|c| c.len() > 1 && c[0] == kind && c[1] == name)
        .unwrap_or_else(|| panic!("no {kind} row {name}"))
}

fn drop_row(lines: &mut Vec<Vec<String>>, kind: &str, name: &str) {
    let n = lines.len();
    lines.retain(|c| !(c.len() > 1 && c[0] == kind && c[1] == name));
    assert_eq!(lines.len(), n - 1, "no {kind} row {name}");
}

const ROH: &str = "{sample}/roh/{sample}_roh_summary.txt";
const HITS: &str = "{sample}/clinvar/{sample}_clinvar_hits.tsv";

#[test]
fn the_contract_keeps_its_own_baseline() {
    let (rc, out) = contract(&["--repo", &s(&repo()), "--baseline", &s(&contract_path())]);
    assert!(
        rc == 0
            && has(
                &out,
                r"baseline: \d{2,} stable rows of .* checked, 0 problems"
            ),
        "{out}"
    );
}

/// The red-first control: a stable output renamed in the code and in its row
/// in one change passes the static check, which reads only this checkout,
/// and fails against the release's contract.
#[test]
fn a_stable_output_renamed_with_its_row_fails_against_the_baseline() {
    let tmp = Tmp::new("baseline-roh");
    let d = code_copy(&tmp);
    edit(
        &d.join("modules/local/roh/main.nf"),
        "path(\"*_roh_summary.txt\")",
        "path(\"*_roh_segments.txt\")",
    );
    let c = contract_edit(&tmp, "renamed", |l| {
        row_mut(l, "file", ROH)[1] = "{sample}/roh/{sample}_roh_segments.txt".to_owned();
    });
    let (rc, out) = contract(&["--static", "--repo", &s(&d), "--contract", &s(&c)]);
    assert!(rc == 0 && out.contains(" 0 problems"), "{out}");
    let base = s(&contract_path());
    let (rc, out) = contract(&[
        "--static",
        "--repo",
        &s(&d),
        "--contract",
        &s(&c),
        "--baseline",
        &base,
    ]);
    assert!(
        rc == 1
            && has(
                &out,
                r"FAIL BREAKING: the stable file \{sample\}/roh/\{sample\}_roh_summary\.txt \(baseline line \d+\) has no row in the contract"
            )
            && out.matches("FAIL ").count() == 1,
        "{out}"
    );
}

#[test]
fn baseline_controls_are_refused() {
    type Change = fn(&mut Vec<Vec<String>>);
    let cases: [(&str, Change, &str); 7] = [
        (
            "a stable param removed with its row",
            |l| drop_row(l, "param", "intervals"),
            r"BREAKING: the stable param intervals \(baseline line \d+\) has no row in the contract",
        ),
        (
            "a stable file made experimental",
            |l| row_mut(l, "file", ROH)[6] = "experimental".to_owned(),
            r"BREAKING: the stable file \{sample\}/roh/\{sample\}_roh_summary\.txt \(baseline line \d+\) is experimental now \(contract line \d+\)",
        ),
        (
            "a stable switch made internal",
            |l| row_mut(l, "switch", "SKIP_TRIM")[6] = "internal".to_owned(),
            r"BREAKING: the stable switch SKIP_TRIM \(baseline line \d+\) is internal now",
        ),
        (
            "a stable file of another format",
            |l| row_mut(l, "file", ROH)[5] = "csv".to_owned(),
            r"BREAKING: the stable file \{sample\}/roh/\{sample\}_roh_summary\.txt \(contract line \d+\) is csv, it was tsv at the release",
        ),
        (
            "a required column dropped",
            |l| {
                let r = row_mut(l, "file", HITS);
                r[8] = r[8].replace(",clnsig,", ",");
            },
            r"BREAKING: the stable file \{sample\}/clinvar/\{sample\}_clinvar_hits\.tsv \(contract line \d+\) no longer requires clnsig",
        ),
        (
            "a leg dropped (both to nextflow)",
            |l| row_mut(l, "file", HITS)[4] = "nextflow".to_owned(),
            r"BREAKING: the stable file \{sample\}/clinvar/\{sample\}_clinvar_hits\.tsv \(baseline line \d+, legs both\) is no longer written by the bash leg",
        ),
        (
            "a stable tool removed with its row",
            |l| drop_row(l, "tool", "pharmcat"),
            r"BREAKING: the stable tool pharmcat \(baseline line \d+\) has no row in the contract",
        ),
    ];
    let base = s(&contract_path());
    for (desc, change, want) in cases {
        let tmp = Tmp::new("baseline");
        let c = contract_edit(&tmp, "changed", change);
        let (rc, out) = contract(&["--contract", &s(&c), "--baseline", &base]);
        assert!(
            rc == 1 && has(&out, want) && out.matches("FAIL ").count() == 1,
            "{desc}: rc={rc}: {out}"
        );
    }
}

#[test]
fn what_the_baseline_allows() {
    type Change = fn(&mut Vec<Vec<String>>);
    let cases: [(&str, Change); 6] = [
        ("a new stable row", |l| {
            l.push(
                "file\t{sample}/new/{sample}_new.tsv\t-\t-\tbash\ttsv\tstable\t-\t-"
                    .split('\t')
                    .map(str::to_owned)
                    .collect(),
            )
        }),
        ("a leg added (nextflow to both)", |l| {
            row_mut(l, "file", ROH)[4] = "both".to_owned()
        }),
        ("a new required column", |l| {
            row_mut(l, "file", HITS)[8].push_str(",new_column")
        }),
        ("an e2e cell and a process changed", |l| {
            let r = row_mut(l, "file", HITS);
            r[7] = "-".to_owned();
            r[3] = "OTHER".to_owned();
        }),
        ("an experimental row removed", |l| {
            drop_row(l, "file", "{sample}/bloodtype/{sample}_bloodtype.tsv")
        }),
        ("an experimental row made stable", |l| {
            row_mut(l, "file", "{sample}/bloodtype/{sample}_bloodtype.tsv")[6] = "stable".to_owned()
        }),
    ];
    let base = s(&contract_path());
    for (desc, change) in cases {
        let tmp = Tmp::new("baseline-ok");
        let c = contract_edit(&tmp, "changed", change);
        let (rc, out) = contract(&["--contract", &s(&c), "--baseline", &base]);
        assert!(
            rc == 0 && out.contains(" 0 problems"),
            "{desc}: rc={rc}: {out}"
        );
    }
}

#[test]
fn a_baseline_that_checks_nothing_or_does_not_parse_fails() {
    let tmp = Tmp::new("baseline-bad");
    let empty = tmp.0.join("empty.tsv");
    fs::write(
        &empty,
        "kind\tname\tstep\tprocess\tlegs\tformat\tstability\te2e\trequired\n",
    )
    .unwrap();
    let (rc, out) = contract(&["--repo", &s(&repo()), "--baseline", &s(&empty)]);
    assert!(
        rc == 1 && out.contains("has no stable row: the check checked nothing"),
        "{out}"
    );
    let bad = tmp.0.join("bad.tsv");
    fs::write(&bad, "not a contract\n").unwrap();
    let (rc, out) = contract(&["--repo", &s(&repo()), "--baseline", &s(&bad)]);
    assert!(rc == 1 && out.contains("ERROR: baseline: "), "{out}");
}

#[test]
fn pedigree_is_a_declared_param() {
    let config = fs::read_to_string(repo().join("nextflow.config")).unwrap();
    assert!(has(&config, r"(?m)^\s*pedigree\s*=\s*null\b"));
    let schema =
        pycompat::json::loads(&fs::read_to_string(repo().join("nextflow_schema.json")).unwrap())
            .unwrap();
    assert!(
        og_nf::schema_params::read(&schema)
            .iter()
            .any(|(n, _)| n == "pedigree")
    );
}

struct ContractRow {
    name: String,
    format: String,
    stability: String,
    e2e: String,
    required: Vec<String>,
}

fn rows() -> Vec<ContractRow> {
    let text = fs::read_to_string(contract_path()).unwrap();
    let mut out = Vec::new();
    let mut header = false;
    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if !header {
            header = true;
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c[0] != "file" {
            continue;
        }
        out.push(ContractRow {
            name: c[1].to_owned(),
            format: c[5].to_owned(),
            stability: c[6].to_owned(),
            e2e: c[7].to_owned(),
            required: if c[8] == "-" {
                Vec::new()
            } else {
                c[8].split(',').map(str::to_owned).collect()
            },
        });
    }
    out
}

fn gzip(text: &str) -> Vec<u8> {
    use std::io::Write;
    let mut g = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    g.write_all(text.as_bytes()).unwrap();
    g.finish().unwrap()
}

/// A results folder with every stable e2e file of a leg, each with its
/// required fields.
fn plant(root: &Path, leg: &str) {
    for r in rows() {
        if r.stability != "stable" || !(r.e2e == leg || r.e2e == "both") {
            continue;
        }
        let path = root.join(r.name.replace("{sample}", SAMPLE).replace('*', "x"));
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        if r.format == "dir" {
            fs::create_dir_all(&path).unwrap();
            continue;
        }
        let text = match r.format.as_str() {
            "vcf" => {
                let mut t = "##fileformat=VCFv4.2\n".to_owned();
                for f in &r.required {
                    let (kind, id) = f.split_once('/').unwrap();
                    t.push_str(&format!(
                        "##{kind}=<ID={id},Number=1,Type=String,Description=\"x\">\n"
                    ));
                }
                t + &format!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}\n")
            }
            "json" => {
                let keys: Vec<String> = r
                    .required
                    .iter()
                    .map(|k| format!("\"{k}\": null"))
                    .collect();
                format!("{{{}}}", keys.join(", "))
            }
            "tsv" | "csv" if !r.required.is_empty() => {
                let sep = if r.format == "tsv" { "\t" } else { "," };
                let ones: Vec<&str> = r.required.iter().map(|_| "1").collect();
                format!(
                    "# a comment line\n{}\n{}\n",
                    r.required.join(sep),
                    ones.join(sep)
                )
            }
            _ => "x\n".to_owned(),
        };
        if path.to_string_lossy().ends_with(".gz") {
            fs::write(&path, gzip(&text)).unwrap();
        } else {
            fs::write(&path, text).unwrap();
        }
    }
}

fn outputs(root: &Path, leg: &str, extra: &[&str]) -> (u8, String) {
    let root = s(root);
    let r = s(&repo());
    let mut args = vec![
        "--repo",
        &r,
        "--outputs",
        &root,
        "--leg",
        leg,
        "--sample",
        SAMPLE,
    ];
    args.extend_from_slice(extra);
    contract(&args)
}

#[test]
fn run_controls() {
    let tmp = Tmp::new("run");
    let checked = regex::Regex::new(r"(\d+) stable files checked").unwrap();
    for leg in ["bash", "nextflow"] {
        let root = tmp.0.join(leg);
        plant(&root, leg);
        let (rc, out) = outputs(&root, leg, &[]);
        let n: usize = checked
            .captures(&out)
            .map(|c| c[1].parse().unwrap())
            .unwrap_or(0);
        assert!(rc == 0 && n >= 5, "a planted {leg} folder: {out}");
    }

    let nf = tmp.0.join("nextflow");
    let hits = nf.join(format!("{SAMPLE}/clinvar/{SAMPLE}_clinvar_hits.tsv"));
    let text = fs::read_to_string(&hits).unwrap();
    let lines: Vec<&str> = text.split('\n').collect();
    let i = lines[1].split('\t').position(|c| c == "clnsig").unwrap();
    let dropped: Vec<String> = lines
        .iter()
        .map(|l| {
            l.split('\t')
                .enumerate()
                .filter(|(j, _)| *j != i)
                .map(|(_, c)| c)
                .collect::<Vec<_>>()
                .join("\t")
        })
        .collect();
    fs::write(&hits, dropped.join("\n")).unwrap();
    let (rc, out) = outputs(&nf, "nextflow", &[]);
    assert!(
        rc == 1
            && has(
                &out,
                r"BREAKING: S1/clinvar/S1_clinvar_hits\.tsv .* has no clnsig"
            ),
        "{out}"
    );

    let roh = nf.join(format!("{SAMPLE}/roh/{SAMPLE}_roh_summary.txt"));
    fs::rename(
        &roh,
        nf.join(format!("{SAMPLE}/roh/{SAMPLE}_roh_segments.txt")),
    )
    .unwrap();
    let (rc, out) = outputs(&nf, "nextflow", &[]);
    assert!(
        rc == 1
            && has(
                &out,
                r"BREAKING: the stable file S1/roh/S1_roh_summary\.txt .* is missing"
            ),
        "{out}"
    );

    let ann = nf.join(format!("{SAMPLE}/vep/{SAMPLE}_annotated.vcf.gz"));
    let mut text = String::new();
    std::io::Read::read_to_string(&mut og_io::open(&ann).unwrap(), &mut text).unwrap();
    let text = regex::Regex::new(r"##INFO=<ID=REVEL,[^\n]*\n")
        .unwrap()
        .replace_all(&text, "");
    fs::write(&ann, gzip(&text)).unwrap();
    let (rc, out) = outputs(&nf, "nextflow", &[]);
    assert!(
        rc == 1 && has(&out, r"S1_annotated\.vcf\.gz .* has no INFO/REVEL"),
        "{out}"
    );

    // A truncated gzip is a named failure, not a crash.
    let bytes = fs::read(&ann).unwrap();
    fs::write(&ann, &bytes[..bytes.len() / 2]).unwrap();
    let (rc, out) = outputs(&nf, "nextflow", &[]);
    assert!(
        rc == 1
            && has(
                &out,
                r"FAIL S1/vep/S1_annotated\.vcf\.gz \(contract line \d+\) is not readable vcf"
            ),
        "{out}"
    );

    let bash = tmp.0.join("bash");
    let summ = bash.join(format!("{SAMPLE}/summary.json"));
    edit(&summ, "\"not_assessed\": null", "\"other\": null");
    let (rc, out) = outputs(&bash, "bash", &[]);
    assert!(
        rc == 1 && has(&out, r"summary\.json .* has no not_assessed"),
        "{out}"
    );

    // An experimental row with e2e set is still not enforced.
    let exp = contract_copy(&tmp, "{sample}/bloodtype/{sample}_bloodtype.tsv", |c| {
        c[7] = "both".to_owned()
    });
    let fresh = tmp.0.join("fresh");
    plant(&fresh, "bash");
    let (rc, out) = outputs(&fresh, "bash", &["--contract", &s(&exp)]);
    assert_eq!(rc, 0, "{out}");
    let (rc, out) = outputs(&tmp.0.join("empty"), "bash", &[]);
    assert!(rc == 1 && out.contains("is missing"), "{out}");
}

#[test]
fn malformed_contracts_are_refused() {
    let tmp = Tmp::new("malformed");
    let text = fs::read_to_string(contract_path()).unwrap();
    let lines: Vec<&str> = text.split('\n').collect();
    let first = lines.iter().position(|l| l.starts_with("file\t")).unwrap();
    let dup = tmp.0.join("dup.tsv");
    let mut d = lines.clone();
    d.insert(first + 1, lines[first]);
    fs::write(&dup, d.join("\n")).unwrap();
    let (rc, out) = contract(&["--static", "--repo", &s(&repo()), "--contract", &s(&dup)]);
    assert!(rc == 1 && out.contains("listed twice"), "{out}");
    let bad = tmp.0.join("bad.tsv");
    let mut b = lines.clone();
    let solid = lines[first].replace("\tstable\t", "\tsolid\t");
    b[first] = &solid;
    fs::write(&bad, b.join("\n")).unwrap();
    let (rc, out) = contract(&["--static", "--repo", &s(&repo()), "--contract", &s(&bad)]);
    assert!(rc == 1 && out.contains("bad stability solid"), "{out}");
    let (rc, out) = contract(&["--outputs", "x"]);
    assert!(
        rc == 2 && out.contains("--outputs needs --leg and --sample"),
        "{out}"
    );
}

// ---------------------------------------------------------------- report-only

fn nf_copy(tmp: &Tmp) -> PathBuf {
    let d = tmp.0.join("repo");
    fs::create_dir_all(&d).unwrap();
    for f in ["main.nf", "nextflow.config"] {
        fs::copy(repo().join(f), d.join(f)).unwrap();
    }
    copy_dir(&repo().join("workflows"), &d.join("workflows"));
    d
}

#[test]
fn report_only_passes_on_this_checkout() {
    let (rc, out) = report_only(&repo());
    assert_eq!(rc, 0, "{out}");
    assert!(
        out.contains("ok: CLINVAR_CDNA_CHECK") && out.contains("ok: CPSR_CDNA_CHECK"),
        "{out}"
    );
}

#[test]
fn report_only_controls() {
    for (name, old, new) in [
        ("CLINVAR_CDNA_CHECK", "|CLINVAR_CDNA_CHECK|", "|"),
        ("CPSR_CDNA_CHECK", "|CPSR_CDNA_CHECK|", "|"),
        (
            "MITO_DELETIONS",
            "withName: 'MITO_DELETIONS'",
            "withName: 'NO_SUCH_PROCESS'",
        ),
    ] {
        let tmp = Tmp::new("ro");
        let d = nf_copy(&tmp);
        edit(&d.join("nextflow.config"), old, new);
        let (rc, out) = report_only(&d);
        assert!(
            rc == 1 && out.contains(&format!("ERROR: {name} (workflows/")),
            "{name}: {out}"
        );
    }

    let tmp = Tmp::new("ro");
    let d = nf_copy(&tmp);
    edit(&d.join("nextflow.config"), "|CLINVAR_CDNA_CHECK|", "|");
    edit(
        &d.join("workflows/pgx.nf"),
        "    clinvar_cdna         = ch_clinvar_cdna",
        "    clinvar_cdna_more    = ch_clinvar_cdna.map { it }\n    clinvar_cdna         = ch_clinvar_cdna",
    );
    let (rc, out) = report_only(&d);
    assert!(
        rc == 0 && !out.contains("CLINVAR_CDNA_CHECK"),
        "a second consumer: {out}"
    );

    let tmp = Tmp::new("ro");
    let d = nf_copy(&tmp);
    let clin = d.join("workflows/clinical.nf");
    let text = fs::read_to_string(&clin).unwrap();
    let m = regex::Regex::new(r"(?m)^\s*roh_regions\s*=\s*(ch_\w+)").unwrap();
    let c = m
        .captures(&text)
        .expect("ROH has its emit in workflows/clinical.nf");
    let whole = c[0].to_owned();
    edit(
        &clin,
        &whole,
        &format!("    roh_more = {}.map {{ it }}\n{whole}", &c[1]),
    );
    let (rc, out) = report_only(&d);
    assert!(rc == 1 && out.contains("KEEP_STOP names ROH"), "{out}");
}

/// Red on the Python: it cut every line at `//`, strings included, and read
/// a strategy closure only up to its first `}`.
#[test]
fn report_only_reads_strings_and_nested_braces() {
    // A // inside a string before the selector: the Python cut the line
    // there and lost the block, so MITO_DELETIONS read as uncovered.
    let tmp = Tmp::new("ro");
    let d = nf_copy(&tmp);
    edit(
        &d.join("nextflow.config"),
        "    withName: 'MITO_DELETIONS' {",
        "    ext.doc = 'https://example.org'; withName: 'MITO_DELETIONS' {",
    );
    let (rc, out) = report_only(&d);
    assert!(rc == 0 && out.contains("ok: MITO_DELETIONS"), "{out}");

    // A nested brace in the closure before 'ignore'.
    let tmp = Tmp::new("ro");
    let d = nf_copy(&tmp);
    let config = d.join("nextflow.config");
    let text = fs::read_to_string(&config).unwrap();
    let block = "withName: 'MITO_DELETIONS' {\n        errorStrategy = { task.exitStatus in [143,137,104,134,139,247] && task.attempt == 1 ? 'retry' : 'ignore' }";
    assert!(text.contains(block));
    fs::write(
        &config,
        text.replacen(
            block,
            "withName: 'MITO_DELETIONS' {\n        errorStrategy = { def why = { it }; task.exitStatus in [143,137,104,134,139,247] && task.attempt == 1 ? 'retry' : 'ignore' }",
            1,
        ),
    )
    .unwrap();
    let (rc, out) = report_only(&d);
    assert!(rc == 0 && out.contains("ok: MITO_DELETIONS"), "{out}");
}
