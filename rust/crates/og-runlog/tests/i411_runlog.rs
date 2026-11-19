//! `og run delta` words the reason it gives for a changed result as a
//! guess. The delta only knows which inputs changed and which inputs feed
//! a section; it cannot know that one of them caused the change. So:
//!
//! 1. a changed result whose section has changed inputs (each kind of
//!    result: ClinVar, CPSR, CPIC, pypgx, the CYP2D6 consensus, PRS) gets a
//!    reason that starts with "likely", in the JSON, the text and the HTML,
//!    and still names each changed input with its old and new value;
//! 2. a changed result with none of its inputs changed gets the line that
//!    names no cause, and that line claims no input;
//! 3. an input change itself carries no reason.

use og_runlog::cli::{Env, run};
use pycompat::json::{Value, loads};
use std::fs;
use std::path::{Path, PathBuf};

/// One summary with a result in each area the delta compares; `v` moves
/// every result between the two runs.
fn summary(v: &str) -> String {
    let (sig, vus, dip, cyp, pct, cons) = if v == "1" {
        ("Pathogenic", 3, "*1/*1", "*1/*2", "40.0", "*1/*2")
    } else {
        ("Benign", 4, "*1/*2", "*1/*4", "60.0", "*1/*4")
    };
    format!(
        r#"{{"generated_utc": "2026-0{v}-01T00:00:00Z", "sections": {{
  "clinvar": {{"state": "ok", "title": "ClinVar screen", "data": {{"hits": [
    {{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G", "gene": "GENEA", "significance": "{sig}", "stars": 2, "genotype": "het"}}]}}}},
  "cpsr": {{"state": "ok", "data": {{"classification": {{"Benign": 10, "VUS": {vus}}}}}}},
  "cpic": {{"state": "ok", "data": {{"genes": [{{"gene": "CYP2C19", "diplotype": "{dip}", "phenotype": "x"}}]}}}},
  "pypgx": {{"state": "ok", "data": {{"cyp2d6": "{cyp}"}}}},
  "prs": {{"state": "ok", "data": {{"scores": [{{"pgs_id": "PGS000001", "condition": "Trait", "percentile": "{pct}"}}]}}}}
}}, "cyp2d6": {{"consensus": {{"result": "{cons}"}}}}}}"#
    )
}

fn manifest(commit: &str, vep: &str) -> String {
    format!(
        "run\twritten_utc\t2026-01-01 00:00:00 UTC\nrun\twritten_by\trun-all.sh\nrun\tgit_commit\t{commit}\n\
         version\tVEP_CACHE_RELEASE\t{vep}\n"
    )
}

fn write_run(d: &Path, summary: &str, manifest: &str) {
    fs::create_dir_all(d).unwrap();
    fs::write(d.join("summary.json"), summary).unwrap();
    fs::write(d.join("run_manifest.tsv"), manifest).unwrap();
}

fn tmp(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-i411-runlog-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

/// `og run delta` of `cur` against `prev`: its JSON, text and HTML.
fn delta(root: &Path, prev: &Path, cur: &Path) -> (Value, String, String) {
    let out = root.join("delta");
    let args: Vec<String> = [
        "delta",
        "--sample",
        "S",
        "--summary",
        cur.join("summary.json").to_str().unwrap(),
        "--manifest",
        cur.join("run_manifest.tsv").to_str().unwrap(),
        "--previous",
        prev.to_str().unwrap(),
        "--out",
        out.to_str().unwrap(),
    ]
    .map(str::to_owned)
    .to_vec();
    let (mut o, mut e) = (Vec::new(), Vec::new());
    let code = run(&args, &Env::system(), &mut o, &mut e);
    assert_eq!(code, 0, "og run delta exits 0: {}", String::from_utf8_lossy(&e));
    let read = |ext: &str| fs::read_to_string(out.with_extension(ext)).unwrap();
    (loads(&read("json")).unwrap(), read("txt"), read("html"))
}

fn get<'a>(v: &'a Value, k: &str) -> Option<&'a Value> {
    match v {
        Value::Object(o) => o.get(k),
        _ => None,
    }
}

fn text(v: Option<&Value>) -> &str {
    match v {
        Some(Value::Str(x)) => x,
        _ => "",
    }
}

fn changes(d: &Value) -> Vec<Value> {
    match get(d, "changes") {
        Some(Value::Array(a)) => a.clone(),
        other => panic!("changes is not a list: {other:?}"),
    }
}

/// The result changes (every area but the inputs and the sections).
fn results(d: &Value) -> Vec<Value> {
    changes(d)
        .into_iter()
        .filter(|c| !matches!(text(get(c, "kind")), "input" | "section"))
        .collect()
}

const AREAS: [&str; 6] = [
    "ClinVar",
    "CPSR",
    "PGx (CPIC)",
    "PGx (pypgx)",
    "PGx (CYP2D6 consensus)",
    "PRS",
];

#[test]
fn every_reason_that_names_inputs_starts_with_likely() {
    let root = tmp("likely");
    let (prev, cur) = (root.join("prev"), root.join("cur"));
    write_run(&prev, &summary("1"), &manifest("abc", "115"));
    write_run(&cur, &summary("2"), &manifest("def", "116"));
    let (d, txt, page) = delta(&root, &prev, &cur);

    let res = results(&d);
    let mut areas: Vec<&str> = res.iter().map(|c| text(get(c, "area"))).collect();
    areas.sort_unstable();
    areas.dedup();
    let mut want = AREAS.to_vec();
    want.sort_unstable();
    assert_eq!(areas, want, "a changed result in each area: {d:?}");

    for c in &res {
        let why = text(get(c, "why"));
        assert!(
            why.starts_with("likely"),
            "{} {}: the reason is worded as a guess: {why:?}",
            text(get(c, "area")),
            text(get(c, "item"))
        );
        assert!(
            why.contains("git_commit abc -> def"),
            "{}: the reason still names the changed commit: {why:?}",
            text(get(c, "area"))
        );
        assert!(
            txt.contains(&format!("    why: {why}\n")),
            "the text prints the reason as written: {why:?}\n{txt}"
        );
        assert!(
            page.contains(&format!("<td>{}</td>", why.replace('>', "&gt;"))),
            "the HTML prints the reason as written: {why:?}"
        );
    }
    let cpsr = res
        .iter()
        .find(|c| text(get(c, "area")) == "CPSR")
        .unwrap();
    let why = text(get(cpsr, "why"));
    assert!(
        why.contains("VEP_CACHE_RELEASE 115 -> 116"),
        "CPSR: every changed input it reads is named: {why:?}"
    );

    for c in changes(&d)
        .iter()
        .filter(|c| text(get(c, "kind")) == "input")
    {
        assert!(
            get(c, "why").is_none(),
            "an input change carries no reason: {c:?}"
        );
    }
    let _ = fs::remove_dir_all(&root);
}

#[test]
fn a_result_with_no_changed_input_names_no_cause() {
    let root = tmp("same");
    let (prev, cur) = (root.join("prev"), root.join("cur"));
    write_run(&prev, &summary("1"), &manifest("abc", "115"));
    write_run(&cur, &summary("2"), &manifest("abc", "115"));
    let (d, _, _) = delta(&root, &prev, &cur);
    let res = results(&d);
    assert_eq!(res.len(), AREAS.len(), "one change per area: {d:?}");
    for c in &res {
        let why = text(get(c, "why"));
        assert_eq!(
            why,
            "the inputs this delta tracks for it are the same; the call itself differs",
            "{}: no changed input, so no cause is named",
            text(get(c, "area"))
        );
        assert!(!why.contains("likely") && !why.contains("->"));
    }
    let _ = fs::remove_dir_all(&root);
}
