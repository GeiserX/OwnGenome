//! The PGx, HLA and mitochondrial sections of the reports (og report
//! render), and the summary data they show (og summary collect), on
//! synthetic inputs. These are the checks of the old
//! tests/test_hla_report.py and tests/test_mito_extended.py, and the rules
//! of the CYP2D6 verdict line and the KIR block.
//!
//! Each detector is shown red: the table checks on a planted bad table, the
//! derivations on a swapped allele, the report lines on the summary without
//! the input that makes them.

use std::fs;
use std::path::{Path, PathBuf};

use og_render::model::Summary;
use og_render::sections::{hla, pgx};
use og_render::{html, text};
use og_summary::collect::{Group, Input, post};
use pycompat::json::{Dump, Value, dumps, loads};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

/// A fresh scratch folder for one test.
fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-render-pgx-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).expect("scratch folder");
    d
}

/// og summary collect on one sample folder, every section, with ASSETS.
fn collect_with(sample: &str, dir: &Path, assets: &Path) -> String {
    let input = Input {
        sample: sample.to_owned(),
        sample_dir: dir.display().to_string(),
        declared_sex: None,
        assets: assets.display().to_string(),
        group: Group::All,
        tables: root().join("bin").display().to_string(),
    };
    let v = post::collect(&input).expect("a summary");
    dumps(
        &v,
        &Dump {
            indent: Some(1),
            ..Dump::default()
        },
    )
    .expect("json")
}

fn collect(sample: &str, dir: &Path) -> String {
    collect_with(sample, dir, &root().join("assets"))
}

fn reports(json: &str) -> (String, String) {
    let s = Summary::parse(json).expect("a summary");
    (
        text::report(&s).expect("the text report"),
        html::report(&s).expect("the HTML report"),
    )
}

/// `doc[a][b]...`, or `Null` where a key is missing.
fn at<'a>(doc: &'a Value, path: &[&str]) -> &'a Value {
    let mut v = doc;
    for k in path {
        v = match v {
            Value::Object(o) => o.get(*k).unwrap_or(&Value::Null),
            _ => &Value::Null,
        };
    }
    v
}

fn s(v: &Value) -> &str {
    match v {
        Value::Str(s) => s,
        _ => "",
    }
}

/// The errors og dev schema finds in a summary.
fn schema_errors(json: &str) -> Vec<String> {
    let schema = loads(
        &fs::read_to_string(root().join("tests/schema/summary.schema.json")).expect("schema"),
    )
    .expect("schema json");
    og_dev::schema::validate(&schema, &loads(json).expect("summary json")).expect("a schema")
}

// --- HLA ------------------------------------------------------------------------------

const COLUMNS: [&str; 5] = ["allele", "kind", "association", "pmids", "citation"];
const REQUIRED: [&str; 14] = [
    "A*31:01",
    "B*15:02",
    "B*57:01",
    "B*58:01",
    "B*27",
    "C*06:02",
    "DQB1*06:02",
    "DRB1*15:01",
    "DRB1*01:01",
    "DRB1*04:01",
    "DRB1*04:04",
    "DRB1*04:05",
    "DRB1*04:08",
    "DRB1*14:02",
];

/// `[A-Z]+[0-9]?\*[0-9]{2}(:[0-9]{2,3})?`: an allele group at a field
/// boundary.
fn allele_group(a: &str) -> bool {
    let Some((gene, fields)) = a.split_once('*') else {
        return false;
    };
    let letters = gene.trim_end_matches(|c: char| c.is_ascii_digit());
    let digits = gene.len() - letters.len();
    let f: Vec<&str> = fields.split(':').collect();
    !letters.is_empty()
        && letters.chars().all(|c| c.is_ascii_uppercase())
        && digits <= 1
        && f[0].len() == 2
        && f.iter()
            .all(|x| !x.is_empty() && x.chars().all(|c| c.is_ascii_digit()))
        && match f.len() {
            1 => true,
            2 => (2..=3).contains(&f[1].len()),
            _ => false,
        }
}

/// Every schema problem of an allele table.
fn table_errors(text: &str) -> Vec<String> {
    let (mut errs, mut head, mut seen) = (Vec::new(), None::<Vec<String>>, Vec::<String>::new());
    for (n, line) in text.lines().enumerate().map(|(i, l)| (i + 1, l)) {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if head.is_none() {
            if c != COLUMNS {
                errs.push(format!("header {c:?}"));
            }
            head = Some(c.iter().map(|x| (*x).to_owned()).collect());
            continue;
        }
        if c.len() != COLUMNS.len() {
            errs.push(format!("line {n}: {} columns", c.len()));
            continue;
        }
        let (allele, kind, assoc, pmids, cite) = (c[0], c[1], c[2], c[3], c[4]);
        if !allele_group(allele) {
            errs.push(format!("line {n}: allele {allele:?}"));
        }
        if !matches!(kind, "drug" | "disease") {
            errs.push(format!("line {n}: kind {kind:?}"));
        }
        if assoc.trim().is_empty() {
            errs.push(format!("line {n}: no association"));
        }
        if pmids.is_empty()
            || !pmids
                .split(',')
                .all(|p| !p.is_empty() && p.chars().all(|c| c.is_ascii_digit()))
        {
            errs.push(format!("line {n}: pmids {pmids:?}"));
        }
        if cite.trim().is_empty() {
            errs.push(format!("line {n}: no citation"));
        }
        if seen.iter().any(|x| x == allele) {
            errs.push(format!("line {n}: duplicate {allele}"));
        }
        seen.push(allele.to_owned());
    }
    errs.extend(
        REQUIRED
            .iter()
            .filter(|a| !seen.iter().any(|x| x == *a))
            .map(|a| format!("missing {a}")),
    );
    errs
}

/// One allele of a T1K row: (allele, quality), or `None` for T1K's
/// one-allele row.
type Pair<'a> = Option<(&'a str, i32)>;

/// T1K's genotype TSV.
fn t1k(rows: &[(&str, &str, i32, Pair<'_>)]) -> String {
    rows.iter()
        .map(|(gene, a1, q1, two)| match two {
            None => format!("{gene}\t1\t{a1}\t30\t{q1}\t.\t0\t-1\n"),
            Some((a2, q2)) => format!("{gene}\t2\t{a1}\t30\t{q1}\t{a2}\t20\t{q2}\n"),
        })
        .collect()
}

/// The summary of a sample folder that holds only T1K's genotype TSV.
fn hla_summary(work: &Path, name: &str, rows: &[(&str, &str, i32, Pair<'_>)]) -> String {
    let d = work.join(name).join("S");
    fs::create_dir_all(d.join("hla_t1k")).expect("folder");
    fs::write(d.join("hla_t1k/S_hla_genotype.tsv"), t1k(rows)).expect("t1k");
    collect("S", &d)
}

fn coeliac(json: &str) -> Value {
    at(
        &loads(json).expect("json"),
        &["sections", "hla", "data", "derived", "coeliac"],
    )
    .clone()
}

fn state_of(json: &str, allele: &str) -> String {
    let v = loads(json).expect("json");
    let Value::Array(rows) = at(&v, &["sections", "hla", "data", "derived", "alleles"]) else {
        panic!("no allele table");
    };
    rows.iter()
        .find(|r| s(at(r, &["allele"])) == allele)
        .map(|r| s(at(r, &["state"])).to_owned())
        .unwrap_or_default()
}

#[test]
fn the_hla_allele_table_passes_its_schema_and_each_planted_problem_is_caught() {
    let good = fs::read_to_string(root().join("bin/hla_alleles.tsv")).expect("table");
    assert_eq!(table_errors(&good), Vec::<String>::new());
    for (what, bad) in [
        ("a wrong kind", good.replacen("\tdrug\t", "\tdrugs\t", 1)),
        (
            "a PubMed id that is not a number",
            good.replace("\t21428769\t", "\tPMID21428769\t"),
        ),
        (
            "a duplicate allele",
            format!("{good}B*27\tdisease\tx\t1\ty\n"),
        ),
        (
            "a missing required allele",
            good.replace("C*06:02\t", "C*07:02\t"),
        ),
        (
            "a three-field allele",
            good.replace("B*58:01\t", "B*58:01:01\t"),
        ),
    ] {
        assert!(!table_errors(&bad).is_empty(), "not caught: {what}");
    }
    // og summary collect reads the same rows, in order.
    let work = scratch("table");
    let json = hla_summary(
        &work,
        "table",
        &[("HLA-B", "B*57:01:01", 40, Some(("B*08:01:01", 40)))],
    );
    let v = loads(&json).expect("json");
    let Value::Array(rows) = at(&v, &["sections", "hla", "data", "derived", "alleles"]) else {
        panic!("no allele table");
    };
    let got: Vec<&str> = rows.iter().map(|r| s(at(r, &["allele"]))).collect();
    let want: Vec<&str> = good
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .skip(1)
        .map(|l| l.split('\t').next().unwrap_or(""))
        .collect();
    assert_eq!(got, want);
    let _ = fs::remove_dir_all(work);
}

#[test]
fn the_coeliac_dq_types_and_the_dqb1_02_dose_follow_the_alleles() {
    let work = scratch("coeliac");
    let cis = hla_summary(
        &work,
        "cis",
        &[
            ("HLA-DQA1", "DQA1*05:01:01", 40, Some(("DQA1*01:01:01", 40))),
            ("HLA-DQB1", "DQB1*02:01:01", 40, Some(("DQB1*05:01:01", 40))),
            ("HLA-DRB1", "DRB1*03:01:01", 40, Some(("DRB1*01:01:01", 40))),
            ("HLA-B", "B*57:01:01", 40, Some(("B*08:01:01", 40))),
        ],
    );
    let c = coeliac(&cis);
    assert_eq!(s(at(&c, &["status"])), "derived");
    assert_eq!(s(at(&c, &["dq2_5"])), "present");
    assert!(
        s(at(&c, &["dq2_5_phase"])).starts_with("likely cis"),
        "{c:?}"
    );
    assert_eq!(
        (s(at(&c, &["dq2_2"])), s(at(&c, &["dq8"]))),
        ("absent", "absent")
    );
    assert_eq!(at(&c, &["dqb1_02_dose"]), &Value::Int(1));
    assert_eq!(schema_errors(&cis), Vec::<String>::new());
    let trans = coeliac(&hla_summary(
        &work,
        "trans",
        &[
            ("HLA-DQA1", "DQA1*05:05:01", 40, Some(("DQA1*02:01:01", 40))),
            ("HLA-DQB1", "DQB1*03:01:01", 40, Some(("DQB1*02:02:01", 40))),
            ("HLA-DRB1", "DRB1*11:01:01", 40, Some(("DRB1*07:01:01", 40))),
        ],
    ));
    assert_eq!(s(at(&trans, &["dq2_5"])), "present");
    assert!(s(at(&trans, &["dq2_5_phase"])).starts_with("possible trans"));
    assert_eq!(s(at(&trans, &["dq2_2"])), "present");
    let dq22 = coeliac(&hla_summary(
        &work,
        "dq22",
        &[
            ("HLA-DQA1", "DQA1*02:01:01", 40, Some(("DQA1*01:02:01", 40))),
            ("HLA-DQB1", "DQB1*02:02:01", 40, Some(("DQB1*06:02:01", 40))),
        ],
    ));
    assert_eq!(
        (
            s(at(&dq22, &["dq2_5"])),
            s(at(&dq22, &["dq2_2"])),
            s(at(&dq22, &["dq8"]))
        ),
        ("absent", "present", "absent")
    );
    assert_eq!(at(&dq22, &["dq2_5_phase"]), &Value::Null);
    let dq8 = coeliac(&hla_summary(
        &work,
        "dq8",
        &[
            ("HLA-DQA1", "DQA1*03:01:01", 40, Some(("DQA1*01:01:01", 40))),
            ("HLA-DQB1", "DQB1*03:02:01", 40, Some(("DQB1*05:01:01", 40))),
        ],
    ));
    assert_eq!(
        (s(at(&dq8, &["dq8"])), s(at(&dq8, &["dq2_5"]))),
        ("present", "absent")
    );
    assert_eq!(at(&dq8, &["dqb1_02_dose"]), &Value::Int(0));
    // Red control: DQA1*05:01 swapped for DQA1*01:01, so no DQ2.5.
    let red = coeliac(&hla_summary(
        &work,
        "red",
        &[
            ("HLA-DQA1", "DQA1*01:01:01", 40, Some(("DQA1*01:02:01", 40))),
            ("HLA-DQB1", "DQB1*02:01:01", 40, Some(("DQB1*05:01:01", 40))),
            ("HLA-DRB1", "DRB1*03:01:01", 40, Some(("DRB1*01:01:01", 40))),
        ],
    ));
    assert_eq!(s(at(&red, &["dq2_5"])), "absent");
    assert_eq!(at(&red, &["dq2_5_phase"]), &Value::Null);
    let low = coeliac(&hla_summary(
        &work,
        "low",
        &[
            ("HLA-DQA1", "DQA1*05:01:01", 0, Some(("DQA1*01:01:01", 40))),
            ("HLA-DQB1", "DQB1*02:01:01", 40, Some(("DQB1*05:01:01", 40))),
        ],
    ));
    assert_eq!(s(at(&low, &["status"])), "not_determined");
    assert_eq!(at(&low, &["dq2_5"]), &Value::Null);
    assert!(s(at(&low, &["reason"])).contains("HLA-DQA1 low confidence"));
    assert_eq!(at(&low, &["dqb1_02_dose"]), &Value::Int(1));
    let two = coeliac(&hla_summary(
        &work,
        "two",
        &[("HLA-DQB1", "DQB1*02:01:01", 40, Some(("DQB1*02:02:01", 40)))],
    ));
    assert_eq!(at(&two, &["dqb1_02_dose"]), &Value::Int(2));
    let hom = coeliac(&hla_summary(
        &work,
        "hom",
        &[("HLA-DQB1", "DQB1*02:01:01", 40, None)],
    ));
    assert_eq!(
        at(&hom, &["dqb1_02_dose"]),
        &Value::Int(2),
        "T1K's one homozygous allele"
    );
    assert_eq!(s(at(&hom, &["reason"])), "HLA-DQA1 not typed");
    let _ = fs::remove_dir_all(work);
}

#[test]
fn the_allele_table_reads_present_absent_uncertain_and_not_determined() {
    let work = scratch("alleles");
    let cis = hla_summary(
        &work,
        "cis",
        &[("HLA-B", "B*57:01:01", 40, Some(("B*08:01:01", 40)))],
    );
    assert_eq!(state_of(&cis, "B*57:01"), "present");
    let b7 = hla_summary(
        &work,
        "b7",
        &[("HLA-B", "B*07:02:01", 40, Some(("B*08:01:01", 40)))],
    );
    assert_eq!(state_of(&b7, "B*57:01"), "absent");
    assert_eq!(state_of(&b7, "B*58:01"), "absent");
    let grp = hla_summary(
        &work,
        "b27",
        &[
            (
                "HLA-B",
                "B*27:05:02",
                40,
                Some(("B*57:03:01,B*57:01:01", 40)),
            ),
            ("HLA-C", "C*06:02:01", 0, Some(("C*07:01:01", 40))),
        ],
    );
    assert_eq!(state_of(&grp, "B*27"), "present");
    assert_eq!(state_of(&grp, "B*57:01"), "uncertain");
    assert_eq!(state_of(&grp, "C*06:02"), "not_determined");
    assert_eq!(state_of(&grp, "DQB1*06:02"), "not_determined");
    let _ = fs::remove_dir_all(work);
}

#[test]
fn both_reports_show_the_hla_rows_the_derived_lines_and_the_note() {
    let work = scratch("hla-report");
    let many: Vec<String> = (1..=12).map(|i| format!("DRB5*01:{i:02}")).collect();
    let many = many.join(",");
    let cap = hla_summary(
        &work,
        "cap",
        &[
            ("HLA-DRB5", &many, 0, None),
            ("HLA-HFE", "HFE*001:01", 40, Some(("HFE*002:01", 40))),
            ("HLA-DQA1", "DQA1*05:01:01", 40, Some(("DQA1*01:01:01", 40))),
            ("HLA-DQB1", "DQB1*02:01:01", 40, Some(("DQB1*05:01:01", 40))),
            ("HLA-B", "B*57:01:01", 40, Some(("B*08:01:01", 40))),
        ],
    );
    let v = loads(&cap).expect("json");
    let Value::Array(loci) = at(&v, &["sections", "hla", "data", "loci"]) else {
        panic!("no loci");
    };
    let drb5 = loci
        .iter()
        .find(|l| s(at(l, &["gene"])) == "HLA-DRB5")
        .expect("DRB5");
    assert_eq!(
        at(drb5, &["alleles"]),
        &Value::Array(vec![Value::Str("DRB5*01:01,DRB5*01:02,DRB5*01:03".into())])
    );
    assert_eq!(
        at(drb5, &["candidates"]),
        &Value::Array(vec![Value::Int(12)])
    );
    let hfe = loci
        .iter()
        .find(|l| s(at(l, &["gene"])) == "HLA-HFE")
        .expect("HFE");
    assert!(
        s(at(hfe, &["label"])).starts_with("T1K sequence allele, not the HFE C282Y/H63D genotype")
    );
    assert_eq!(schema_errors(&cap), Vec::<String>::new());
    let (txt, page) = reports(&cap);
    assert!(
        txt.contains("HLA-DRB5 no confident call (12 candidates, quality 0); low confidence"),
        "{txt}"
    );
    assert!(!txt.contains("DRB5*01:04") && !page.contains("DRB5*01:04"));
    assert!(page.contains("no confident call (12 candidates, quality 0)"));
    assert_eq!(txt.matches("not the HFE C282Y/H63D genotype").count(), 1);
    assert_eq!(page.matches("not the HFE C282Y/H63D genotype").count(), 1);
    let want = [
        "Coeliac DQ: DQ2.5 present, phase unknown (HLA-DRB1 not typed); DQ2.2 absent; DQ8 absent",
        "DQB1*02 copies: 1",
        "B*57:01 present: abacavir hypersensitivity; flucloxacillin liver injury",
        "C*06:02 not determined (HLA-C not typed)",
    ];
    for w in want {
        assert!(txt.contains(w), "text: {w}");
        assert!(
            page.contains(&og_render::escape::e(w).into_owned()),
            "HTML: {w}"
        );
    }
    assert!(txt.contains("not a diagnosis"));
    let note = og_render::notes::HLA_DERIVED_NOTE;
    assert!(page.contains(&format!(
        "<p class=\"caveat\">{}</p>",
        og_render::escape::e(note)
    )));
    // No clinical advice in the derived lines or the note.
    let advice = |t: &str| {
        let t = t.to_lowercase();
        ["you should", "you must", "consult", "treat"]
            .iter()
            .any(|w| t.contains(w))
    };
    let derived = hla::hla_derived_lines(Some(at(&v, &["sections", "hla", "data", "derived"])))
        .expect("lines");
    assert!(!advice(&format!("{} {note}", derived.join(" "))));
    assert!(
        advice("you should consult"),
        "red control of the advice check"
    );
    // Red control: without T1K's table the section is not run, and neither
    // report prints a derived line.
    let none = collect("S", &work.join("none"));
    let (txt, page) = reports(&none);
    assert!(!txt.contains("Coeliac DQ") && !page.contains("Coeliac DQ"));
    let _ = fs::remove_dir_all(work);
}

/// The KIR block of the text report, as one string.
fn kir_block(txt: &str) -> String {
    let start = txt
        .find("  KIR genes (T1K, KIR=true):")
        .expect("a KIR block");
    let rest = &txt[start..];
    let end = rest.find("\n\n").unwrap_or(rest.len());
    rest[..end].to_owned()
}

#[test]
fn the_kir_block_lists_each_typed_gene_then_the_genes_with_no_allele() {
    let kir = loads(
        r#"{"reason": "", "typed": [{"gene": "KIR2DL1", "alleles": ["KIR2DL1*001", "KIR2DL1*002"],
            "quality": ["30", ""]}, {"gene": "KIR3DL2", "alleles": ["KIR3DL2*007"], "quality": [12, "x"]}],
            "not_called": ["KIR2DS1", "KIR3DS1"]}"#,
    )
    .expect("json");
    assert_eq!(
        hla::kir_lines(&kir).expect("lines"),
        [
            (
                "KIR2DL1".to_owned(),
                "KIR2DL1*001 (quality 30) / KIR2DL1*002".to_owned()
            ),
            ("KIR3DL2".to_owned(), "KIR3DL2*007 (quality 12)".to_owned()),
            (
                String::new(),
                "No allele called (absent, or too few reads): KIR2DS1, KIR3DS1".to_owned()
            ),
        ]
    );
    let why = loads(r#"{"reason": "KIR=false", "typed": [], "not_called": ["X"]}"#).expect("json");
    assert_eq!(
        hla::kir_lines(&why).expect("lines"),
        [(String::new(), "KIR=false".to_owned())]
    );
    let summary = synthetic(&format!(
        r#""hla": {{"title": "HLA", "state": "ok", "source": "x", "file_date": null, "note": null,
            "data": {{"loci": [], "kir": {}}}}}"#,
        dumps(&kir, &Dump::default()).expect("json")
    ));
    let (txt, page) = reports(&summary);
    assert_eq!(
        kir_block(&txt),
        "  KIR genes (T1K, KIR=true):\n    KIR2DL1  KIR2DL1*001 (quality 30) / KIR2DL1*002\n    \
         KIR3DL2  KIR3DL2*007 (quality 12)\n    No allele called (absent, or too few reads): KIR2DS1, \
         KIR3DS1\n  KIR genotypes as T1K typed them. No step interprets them: they matter for transplant \
         matching and some immune conditions, not for drug guidance."
    );
    assert!(page.contains(
        "<span class=\"label\">KIR2DL1</span><span class=\"value\">KIR2DL1*001 (quality 30) / KIR2DL1*002</span>"
    ));
    assert!(page.contains(
        "<p class=\"note\">No allele called (absent, or too few reads): KIR2DS1, KIR3DS1</p>"
    ));
}

// --- CYP2D6 ---------------------------------------------------------------------------

/// A summary with the frame's keys and the sections and top-level keys of
/// `sections` (a JSON object body).
fn synthetic(sections: &str) -> String {
    synthetic_with(sections, r#""cyp2d6": null"#)
}

fn synthetic_with(sections: &str, top: &str) -> String {
    format!(
        r#"{{"schema_version": 2, "sample": "S1", "generated_utc": "2026-11-06 12:00 UTC",
 "run": {{"started_utc": null, "declared_sex": null}},
 "databases": {{"clinvar_release": null, "cpsr_clinvar_release": null, "hla_database": "IPD-IMGT/HLA 3.60"}},
 "manifest": [], "sections": {{{sections}}}, {top}, "not_run": [], "not_assessed": []}}"#
    )
}

#[test]
fn the_cyp2d6_verdict_line_says_where_the_call_went_and_why() {
    let c =
        |cons: &str| loads(&format!(r#"{{"calls": {{}}, "consensus": {cons}}}"#)).expect("json");
    assert_eq!(
        pgx::cyp2d6_verdict(&c(
            r#"{"result": "*1/*2", "passed_to_pharmcat": true, "reason": "pypgx and Cyrius agree"}"#
        ))
        .expect("line"),
        "*1/*2, passed to PharmCAT (pypgx and Cyrius agree)"
    );
    assert_eq!(
        pgx::cyp2d6_verdict(&c(r#"{"result": null, "passed_to_pharmcat": false}"#)).expect("line"),
        "?, not passed to PharmCAT (no reason given)"
    );
    assert_eq!(pgx::cyp2d6_verdict(&c("null")).expect("line"), "");
    let top = r#""cyp2d6": {"calls": {"PharmCAT": "*1/*4", "pypgx": "*1/*2", "Cyrius": null},
        "agree": false, "consensus": {"result": "*1/*4", "passed_to_pharmcat": false,
        "reason": "the callers disagree"}, "probable": null}"#;
    let (txt, page) = reports(&synthetic_with("", top));
    assert!(
        txt.contains(
            "  PharmCAT (step 36's call)  *1/*4\n  pypgx                      *1/*2\n  Cyrius                     not run\n  \
             pypgx and Cyrius agree: NO, review before acting on CYP2D6\n  \
             Step 36: *1/*4, not passed to PharmCAT (the callers disagree)\n"
        ),
        "{txt}"
    );
    assert!(page.contains(
        "<span class=\"label\">Step 36</span><span class=\"value\"><span class=\"badge badge-yellow\">*1/*4</span>"
    ));
    assert!(page.contains("<p class=\"note\">Not passed to PharmCAT: the callers disagree</p>"));
    // Red control: with no call and no verdict the text report has no block.
    let empty = r#""cyp2d6": {"calls": {"PharmCAT": null}, "agree": null, "consensus": null, "probable": null}"#;
    let (txt, _) = reports(&synthetic_with("", empty));
    assert!(!txt.contains("CYP2D6 across callers"));
    // A summary without the block (collected by part) shows none.
    let (txt, page) = reports(&synthetic(""));
    assert!(!txt.contains("CYP2D6 across callers") && !page.contains("CYP2D6 Across Callers"));
}

// --- mitochondrial variants --------------------------------------------------------

const EK_HEAD: &str = "\"Title\";\"5' breakpoint\";\"3' breakpoint\";\"Freq\";\"Freq For\";\"Freq Rev\";\
\"5' Blast\";\"3' Blast\";\"5' Depth\";\"3' Depth\";\"Repetition\"\n";

fn run_mito(verb: &str, args: &[&str]) -> u8 {
    let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
    og_summary::mito::cli(Some(verb), &args, &mut Vec::new(), &mut Vec::new())
}

/// A sample folder with a chrM VCF of (pos, ref, alt, AF), and eKLIPse's
/// deletions and the idxstats when given.
fn mito_dir(
    work: &Path,
    name: &str,
    records: &[(u32, &str, &str, &str)],
    deletions: Option<&[(u32, u32, f64)]>,
    idxstats: Option<String>,
) -> PathBuf {
    let d = work.join(name);
    fs::create_dir_all(d.join("mito")).expect("folder");
    let mut vcf = "##fileformat=VCFv4.2\n##contig=<ID=chrM,length=16569>\n\
                   #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        .to_owned();
    for (pos, r, a, af) in records {
        vcf += &format!("chrM\t{pos}\t.\t{r}\t{a}\t.\tPASS\t.\tGT:AF\t0/1:{af}\n");
    }
    fs::write(d.join(format!("mito/{name}_chrM_filtered.vcf.gz")), vcf).expect("vcf");
    if let Some(rows) = deletions {
        let mut csv = EK_HEAD.to_owned();
        for (b5, b3, lvl) in rows {
            let v = lvl.to_string().replace('.', ",");
            csv += &format!(
                "\"{name}\";\"{b5}\";\"{b3}\";\"{v}\";\"{v}\";\"{v}\";\"12\";\"14\";\"390.0\";\"402.0\";\"none\"\n"
            );
        }
        let p = work.join(format!("{name}_eklipse.csv"));
        fs::write(&p, csv).expect("csv");
        let out = d.join(format!("mito/{name}_chrM_deletions.tsv"));
        assert_eq!(
            run_mito(
                "deletions",
                &[
                    "--eklipse",
                    &p.display().to_string(),
                    "--out",
                    &out.display().to_string()
                ]
            ),
            0
        );
    }
    if let Some(text) = idxstats {
        let p = work.join(format!("{name}.idxstats"));
        fs::write(&p, text).expect("idxstats");
        let out = d.join(format!("mito/{name}_mtdna_copies.tsv"));
        assert_eq!(
            run_mito(
                "copies",
                &[
                    "--idxstats",
                    &p.display().to_string(),
                    "--out",
                    &out.display().to_string()
                ]
            ),
            0
        );
    }
    d
}

fn mito_data(json: &str) -> Value {
    at(&loads(json).expect("json"), &["sections", "mito", "data"]).clone()
}

fn idxstats(chrm: u64, per_chrom: u64) -> String {
    let mut s: String = (1..=22)
        .map(|i| format!("chr{i}\t50000000\t{per_chrom}\t0\n"))
        .collect();
    s += &format!("chrX\t156040895\t{per_chrom}\t0\nchrM\t16569\t{chrm}\t0\n*\t0\t0\t500\n");
    s
}

#[test]
fn the_mito_list_has_its_columns_an_id_on_every_row_and_its_first_rows() {
    let text = fs::read_to_string(root().join("assets/mito_confirmed.tsv")).expect("list");
    let errors = |t: &str| -> Vec<String> {
        let mut rows = t.lines().filter(|l| !l.starts_with('#'));
        let head = rows.next().unwrap_or("");
        let cols = [
            "kind",
            "pos",
            "end",
            "ref",
            "alt",
            "hgvs",
            "gene",
            "disease",
            "clinvar_id",
            "pmids",
            "heteroplasmy_note",
        ];
        if head.split('\t').collect::<Vec<_>>() != cols {
            return vec!["header".into()];
        }
        let digits = |x: &str| !x.is_empty() && x.chars().all(|c| c.is_ascii_digit());
        let mut out = Vec::new();
        for r in rows {
            let c: Vec<&str> = r.split('\t').collect();
            if c.len() != cols.len() {
                out.push(format!("{r}: columns"));
                continue;
            }
            let (cv, pm) = (!matches!(c[8], "" | "."), !matches!(c[9], "" | "."));
            if !cv && !pm {
                out.push(format!("{}: neither a ClinVar id nor a PubMed id", c[5]));
            }
            if cv && !digits(c[8]) {
                out.push(format!("{}: ClinVar id", c[5]));
            }
            if pm && !c[9].split(';').all(digits) {
                out.push(format!("{}: PubMed ids", c[5]));
            }
            if !matches!(c[0], "snv" | "deletion") {
                out.push(format!("{}: kind", c[5]));
            } else {
                let (p, e) = (c[1].parse::<u32>().ok(), c[2].parse::<u32>().ok());
                if !matches!((p, e), (Some(p), Some(e)) if 1 <= p && p <= e && e <= 16569) {
                    out.push(format!("{}: positions", c[5]));
                }
            }
        }
        out
    };
    assert_eq!(errors(&text), Vec::<String>::new());
    for want in [
        "m.3243A>G",
        "m.8344A>G",
        "m.11778G>A",
        "m.3460G>A",
        "m.14484T>C",
    ] {
        assert!(text.contains(&format!("\t{want}\t")), "{want}");
    }
    assert!(
        text.lines()
            .any(|l| l.starts_with("deletion\t8470\t13446\t"))
    );
    let bad =
        format!("{text}snv\t3271\t3271\tT\tC\tm.3271T>C\tMT-TL1\tMELAS\t.\t.\tno id at all\n");
    assert!(
        errors(&bad)
            .iter()
            .any(|e| e.contains("neither a ClinVar id nor a PubMed id")),
        "red control: a planted row with neither id"
    );
}

#[test]
fn a_confirmed_variant_shows_with_its_level_in_the_summary_and_both_reports() {
    let work = scratch("confirmed");
    let d = mito_dir(
        &work,
        "S3243",
        &[
            (73, "A", "G", "1.00"),
            (3243, "A", "G", "0.12"),
            (16519, "T", "C", "0.998"),
        ],
        None,
        None,
    );
    let json = collect("S3243", &d);
    let data = mito_data(&json);
    let Value::Array(conf) = at(&data, &["confirmed"]) else {
        panic!("no confirmed list: {data:?}");
    };
    assert_eq!(conf.len(), 1);
    assert_eq!(s(at(&conf[0], &["hgvs"])), "m.3243A>G");
    assert_eq!(s(at(&conf[0], &["af"])), "0.12");
    let (txt, page) = reports(&json);
    assert!(
        txt.contains("m.3243A>G (MT-TL1)  AF 0.12  filter PASS  MELAS"),
        "{txt}"
    );
    assert!(page.contains(
        "Confirmed disease variants</span><span class=\"value\"><span class=\"badge badge-yellow\">1<"
    ));
    assert!(page.contains("m.3243A&gt;G (MT-TL1)  AF 0.12"));
    // Red control: the same VCF without m.3243A>G lists none.
    let d = mito_dir(
        &work,
        "Snone",
        &[(73, "A", "G", "1.00"), (16519, "T", "C", "0.998")],
        None,
        None,
    );
    let json = collect("Snone", &d);
    assert_eq!(
        at(&mito_data(&json), &["confirmed"]),
        &Value::Array(Vec::new())
    );
    let (txt, _) = reports(&json);
    assert!(txt.contains("Confirmed disease variants (assets/mito_confirmed.tsv): 0"));
    // Without the list there is no count, only why.
    let no_assets = work.join("no-assets");
    fs::create_dir_all(&no_assets).expect("folder");
    let json = collect_with("Snone", &d, &no_assets);
    let data = mito_data(&json);
    assert_eq!(at(&data, &["confirmed"]), &Value::Null);
    assert!(s(at(&data, &["confirmed_note"])).contains("mito_confirmed.tsv is not in"));
    let (txt, page) = reports(&json);
    assert!(!txt.contains("Confirmed disease variants (assets"));
    assert!(!page.contains("Confirmed disease variants</span>"));
    assert!(txt.contains("Confirmed disease variants: not checked"));
    let _ = fs::remove_dir_all(work);
}

#[test]
fn the_copy_number_shows_in_both_reports_or_why_it_is_missing() {
    let work = scratch("copies");
    let per = 200 * 50_000;
    let chrm = (100.0 * 16569.0 * per as f64 / 50_000_000.0).round() as u64;
    let d = mito_dir(&work, "Snocn", &[(73, "A", "G", "1.00")], None, None);
    let json = collect("Snocn", &d);
    let data = mito_data(&json);
    assert_eq!(at(&data, &["copies"]), &Value::Null);
    assert!(s(at(&data, &["copies_note"])).contains("mtdna_copies.tsv"));
    let (txt, page) = reports(&json);
    assert!(txt.contains("mtDNA copy number: not counted ("));
    assert!(page.contains("mtDNA copy number: not counted ("));
    let d = mito_dir(
        &work,
        "Scn",
        &[(73, "A", "G", "1.00")],
        None,
        Some(idxstats(chrm, per)),
    );
    let (txt, page) = reports(&collect("Scn", &d));
    let n = |t: &str, lead: &str| -> f64 {
        let at = t.find(lead).expect(lead) + lead.len();
        t[at..]
            .split(|c: char| !(c.is_ascii_digit() || c == '.'))
            .next()
            .and_then(|x| x.parse().ok())
            .unwrap_or(0.0)
    };
    let got = n(&txt, "mtDNA copy number: ");
    assert!((200.0..210.0).contains(&got), "{got}");
    let got = n(&page, "mtDNA copies per cell</span><span class=\"value\">");
    assert!((200.0..210.0).contains(&got), "{got}");
    let _ = fs::remove_dir_all(work);
}

#[test]
fn the_large_deletions_are_counted_listed_and_named() {
    let work = scratch("deletions");
    let d = mito_dir(
        &work,
        "Sdel",
        &[(73, "A", "G", "1.00")],
        Some(&[
            (8469, 13447, 10.2),
            (8482, 13460, 2.5),
            (9000, 14000, 3.1),
            (5000, 7000, 0.4),
        ]),
        None,
    );
    let json = collect("Sdel", &d);
    let dl = at(&mito_data(&json), &["deletions"]).clone();
    let Value::Array(listed) = at(&dl, &["listed"]) else {
        panic!("no deletions: {dl:?}");
    };
    let num = |v: &Value| match v {
        Value::Int(i) => *i,
        _ => -1,
    };
    let named: Vec<(i128, i128)> = listed
        .iter()
        .filter(|r| !matches!(at(r, &["known"]), Value::Null))
        .map(|r| (num(at(r, &["bkp5"])), num(at(r, &["bkp3"]))))
        .collect();
    assert_eq!(named, [(8469, 13447), (8482, 13460)]);
    assert!(
        listed
            .iter()
            .any(|r| num(at(r, &["bkp5"])) == 9000 && matches!(at(r, &["known"]), Value::Null))
    );
    assert_eq!((at(&dl, &["total"]), listed.len()), (&Value::Int(4), 3));
    let (txt, page) = reports(&json);
    assert!(page.contains(">Large deletions (eKLIPse)</span><span class=\"value\">4<"));
    assert!(page.contains(
        "Named in the confirmed list</span><span class=\"value\"><span class=\"badge badge-yellow\">2<"
    ));
    let line = txt
        .lines()
        .find(|l| l.contains("m.8470_13446del (4977 bp)  level 10.2%"))
        .expect("the common deletion");
    assert!(
        line.contains("The common 4977 bp deletion") && line.contains("ClinVar 430680"),
        "{line}"
    );
    assert!(txt.contains("often found in muscle and not in blood"));
    // Red control: a table with no deletion lists none.
    let d = mito_dir(&work, "Snodel", &[(73, "A", "G", "1.00")], Some(&[]), None);
    let dl = at(&mito_data(&collect("Snodel", &d)), &["deletions"]).clone();
    assert_eq!(at(&dl, &["total"]), &Value::Int(0));
    assert_eq!(at(&dl, &["listed"]), &Value::Array(Vec::new()));
    let _ = fs::remove_dir_all(work);
}

#[test]
fn the_eklipse_image_is_pinned_by_digest_and_names_its_upstream_commit() {
    let env = fs::read_to_string(root().join("versions.env")).expect("versions.env");
    let line = env
        .lines()
        .find(|l| l.starts_with("EKLIPSE_IMAGE="))
        .expect("EKLIPSE_IMAGE");
    let image = line
        .strip_prefix("EKLIPSE_IMAGE=\"")
        .and_then(|r| r.split('"').next())
        .unwrap_or("");
    let (name, digest) = image.split_once("@sha256:").unwrap_or(("", ""));
    assert!(name.contains(':') && digest.len() == 64, "{line}");
    assert!(
        digest
            .chars()
            .all(|c| c.is_ascii_hexdigit() && !c.is_ascii_uppercase())
    );
    assert!(
        line.contains("3606cb2edac983d2623ddc667b49206c3d01373c"),
        "{line}"
    );
}
