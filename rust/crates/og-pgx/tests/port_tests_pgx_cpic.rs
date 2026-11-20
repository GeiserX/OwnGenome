//! `og pgx cpic-report` and `og pgx step36-outside`, the one PharmCAT
//! reader (step 27 and the CPIC_LOOKUP module run the same verbs), run
//! through [`og_pgx::main`] as the `og` binary runs them. The
//! function-level rules (the phenotype classes, the layouts, the activity
//! score bands) are this crate's unit tests; these run the verbs.
//!
//! 1. The real PharmCAT `report.json` of the HG002 fixture, from each
//!    PharmCAT release in [`REAL_REPORTS`]
//!    (`tests/fixtures/pharmcat/report-<version>.json`), parses to genes; a
//!    gene PharmCAT lists several possible diplotypes for, with different
//!    phenotypes (CYP2C19, 528 of them), is 'ambiguous' and never read as
//!    its first diplotype; a report that parses to zero genes exits 1 with a
//!    PARSING FAILED report, never an all-clear.
//! 2. On PharmCAT's own example report exactly 3 genes count as non-normal
//!    (CYP2D6, HLA-B, NAT2); the 7 whose phenotype is n/a, 'no phenotype
//!    assigned', a genotype or 'Uncertain Susceptibility' are
//!    'unclassified', listed on their own and given no drug list; every
//!    non-normal gene gets a drug list, with the recommendation for its own
//!    diplotype; synthetic genes check the fallbacks: the report's `drugs`
//!    section first, then `relatedDrugs`, then the static table, and
//!    otherwise a line saying it is not in the drug table.
//! 3. The flat (3.x), nested (2.x) and list layouts all parse.
//! 4. A gene PharmCAT could not call (or could not resolve) while pypgx did
//!    prints a warning naming every drug PharmCAT links to the gene, or a
//!    note when pypgx's call is normal, and the comparison marks it 'pypgx
//!    only'; an unreadable pypgx summary is reported, not fatal, and makes
//!    no comparison.
//! 5. A gene PharmCAT leaves 'Indeterminate' at gene level while its drug
//!    annotations give it one phenotype shows that phenotype, with a note,
//!    and takes its status; with no annotation, two that disagree, one that
//!    is Indeterminate itself, or possible diplotypes whose guidance
//!    differs, it stays Indeterminate and unclassified.
//! 6. A single-position gene PharmCAT labels 'rsX reference (C)/rsX variant
//!    (T)' matches pypgx's 'Reference/rsX', and a real difference still does
//!    not; the comparison carries both tools' phenotypes; a CPIC
//!    recommendation over 200 characters prints whole; CYP4F2 (and VKORC1)
//!    with no PharmCAT result takes the genotype PharmCAT read at the *3
//!    (-1639) position as a labelled single-position fallback, and a
//!    position with no genotype is 'not callable'.
//! 7. An MT-RNR1 outside call is not checked against step 36's consensus
//!    table, which never lists it, while an unconfirmed CYP2D6 outside call
//!    is still held back; `step36-outside` lists the step 36 genes the
//!    report calls from an outside call and never MT-RNR1.
//! 8. Phenoconversion (opt-in `--medications` with
//!    `assets/pgx_inhibitors.tsv`): the table's shape, a strong and a
//!    moderate CYP2D6 inhibitor, CYP2C19, CYP2C9 and CYP3A notes, a name not
//!    in the table, an empty list, and the control: the inhibitor's factor
//!    set to 1 must not give a poor metabolizer.
//! 9. Disease-variant contraindications (`--contraindications` with
//!    `assets/pgx_contraindications.tsv` and step 06's `--clinvar-hits`):
//!    the table's shape, a pathogenic POLG hit and the hits that add
//!    nothing, no or an unreadable hits table, and the control: the table
//!    without its POLG row gives no valproic acid line for the same hit.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use pycompat::json::{Dump, Value, dumps, loads};

/// One real HG002 report per PharmCAT release the reader is tested on: the
/// version, and the genes it lists more than one possible diplotype for,
/// with how many.
const REAL_REPORTS: [(&str, [(&str, usize); 2]); 2] = [
    ("3.2.0", [("CYP2B6", 4), ("CYP2C19", 528)]),
    ("3.4.0", [("CYP2B6", 4), ("CYP2C19", 528)]),
];
const PHENOCONVERSION_GENES: [&str; 5] = ["CYP2D6", "CYP2C19", "CYP2C9", "CYP3A4", "CYP3A5"];
const CONTRA_COLUMNS: [&str; 8] = [
    "Gene",
    "Drug",
    "Names",
    "CPIC_level",
    "Note",
    "Pair_source",
    "Label_source",
    "Report_if",
];
const FDA_TABLE: &str = "https://www.fda.gov/drugs/drug-interactions-labeling/\
                         healthcare-professionals-fdas-examples-drugs-interact-cyp-enzymes-and-transporter-systems";

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn repo_path(rel: &str) -> String {
    repo().join(rel).to_string_lossy().into_owned()
}

fn inhibitors() -> String {
    repo_path("assets/pgx_inhibitors.tsv")
}

fn contra() -> String {
    repo_path("assets/pgx_contraindications.tsv")
}

fn fixture(name: &str) -> Value {
    loads(&fs::read_to_string(repo().join("tests/fixtures/pharmcat").join(name)).unwrap()).unwrap()
}

/// A scratch directory of one test: emptied when made, removed when
/// dropped.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Self {
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("pgx-cpic-{name}"));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).to_string_lossy().into_owned()
    }

    fn file(&self, name: &str, text: &str) -> String {
        let p = self.path(name);
        fs::write(&p, text).unwrap();
        p
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

// --- JSON built as Python's dicts are -----------------------------------------

fn s(x: &str) -> Value {
    Value::Str(x.to_owned())
}

fn o(pairs: Vec<(&str, Value)>) -> Value {
    Value::Object(pairs.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}

fn a(items: Vec<Value>) -> Value {
    Value::Array(items)
}

/// `dict(v, key=val)`.
fn with(mut v: Value, key: &str, val: Value) -> Value {
    if let Value::Object(m) = &mut v {
        m.insert(key.to_owned(), val);
    }
    v
}

fn at<'a>(v: &'a mut Value, key: &str) -> &'a mut Value {
    match v {
        Value::Object(m) => m.get_mut(key).unwrap(),
        other => panic!("not an object: {other:?}"),
    }
}

fn dip(a1: &str, a2: &str, phen: &str) -> Value {
    o(vec![(
        "sourceDiplotypes",
        a(vec![o(vec![
            ("allele1", o(vec![("name", s(a1))])),
            ("allele2", o(vec![("name", s(a2))])),
            ("label", s(&format!("{a1}/{a2}"))),
            ("phenotypes", a(vec![s(phen)])),
        ])]),
    )])
}

fn genes(pairs: Vec<(&str, Value)>) -> Value {
    o(vec![("genes", o(pairs))])
}

/// `{"CPIC Guideline Annotation": {drug: {"name": drug, "guidelines":
/// [{"annotations": anns}]}}}`.
fn cpic_drugs(drug: &str, anns: Vec<Value>) -> Value {
    o(vec![(
        "CPIC Guideline Annotation",
        o(vec![(
            drug,
            o(vec![
                ("name", s(drug)),
                ("guidelines", a(vec![o(vec![("annotations", a(anns))])])),
            ]),
        )]),
    )])
}

// --- one run of the verb --------------------------------------------------------

#[derive(Default)]
struct Opts<'a> {
    pypgx: Option<Vec<(&'a str, &'a str, &'a str)>>,
    consensus: Option<Vec<(&'a str, &'a str, &'a str)>>,
    medications: Option<&'a str>,
    inhibitors: Option<String>,
    contraindications: Option<String>,
    hits: Option<String>,
}

/// What one run wrote: exit code, recommendations text, phenotypes rows,
/// comparison rows.
struct Report {
    rc: u8,
    rec: String,
    rows: Vec<Vec<String>>,
    comp: Vec<Vec<String>>,
}

fn og(args: &[&str]) -> u8 {
    let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
    og_pgx::main(&args)
}

fn rows_of(path: &str) -> Vec<Vec<String>> {
    fs::read_to_string(path)
        .map(|t| {
            t.lines()
                .skip(1)
                .map(|l| l.split('\t').map(str::to_owned).collect())
                .collect()
        })
        .unwrap_or_default()
}

fn run_report(work: &Scratch, name: &str, data: &Value, o: &Opts) -> Report {
    let d = work.0.join(name);
    fs::create_dir_all(&d).unwrap();
    let p = |n: &str| d.join(n).to_string_lossy().into_owned();
    let rep = p("report.json");
    fs::write(&rep, dumps(data, &Dump::default()).unwrap()).unwrap();
    let mut args: Vec<String> = [
        "cpic-report",
        "--sample",
        name,
        "--report",
        &rep,
        "--outdir",
    ]
    .map(str::to_owned)
    .to_vec();
    args.push(d.to_string_lossy().into_owned());
    let mut opt = |k: &str, v: String| {
        args.push(k.to_owned());
        args.push(v);
    };
    if let Some(rows) = &o.consensus {
        let mut text = "Gene\tResult\tOutside_call\tReason\tEvidence\n".to_owned();
        for (g, res, out) in rows {
            text.push_str(&format!("{g}\t{res}\t{out}\ttest\ttest\n"));
        }
        fs::write(p("consensus.tsv"), text).unwrap();
        opt("--consensus", p("consensus.tsv"));
    }
    if let Some(m) = o.medications {
        fs::write(p("medications.txt"), m).unwrap();
        opt("--medications", p("medications.txt"));
        opt(
            "--inhibitors",
            o.inhibitors.clone().unwrap_or_else(inhibitors),
        );
    }
    if let Some(c) = &o.contraindications {
        opt("--contraindications", c.clone());
    }
    if let Some(h) = &o.hits {
        fs::write(p("clinvar_hits.tsv"), h).unwrap();
        opt("--clinvar-hits", p("clinvar_hits.tsv"));
    }
    let comp = p("comparison.tsv");
    if let Some(rows) = &o.pypgx {
        let mut text = "Gene\tDiplotype\tPhenotype\tCNV_call\tSource\n".to_owned();
        for (g, dl, ph) in rows {
            text.push_str(&format!("{g}\t{dl}\t{ph}\t.\tbam\n"));
        }
        fs::write(p("pypgx.tsv"), text).unwrap();
        opt("--pypgx", p("pypgx.tsv"));
        opt("--comparison", comp.clone());
    }
    let rc = og_pgx::main(&args);
    Report {
        rc,
        rec: fs::read_to_string(p(&format!("{name}_cpic_recommendations.txt"))).unwrap(),
        rows: rows_of(&p(&format!("{name}_phenotypes.tsv"))),
        comp: rows_of(&comp),
    }
}

fn plain(work: &Scratch, name: &str, data: &Value) -> Report {
    run_report(work, name, data, &Opts::default())
}

fn row<'a>(rows: &'a [Vec<String>], gene: &str) -> Option<&'a Vec<String>> {
    rows.iter().find(|r| r[0] == gene)
}

fn cell(rows: &[Vec<String>], gene: &str, i: usize) -> String {
    row(rows, gene)
        .and_then(|r| r.get(i).cloned())
        .unwrap_or_default()
}

fn strs(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| (*s).to_owned()).collect()
}

/// The medications block of one gene in the recommendations text.
fn block(rec: &str, gene: &str) -> String {
    let mut lines = Vec::new();
    let mut on = false;
    for line in rec.lines() {
        if line.starts_with(&format!("  {gene} -- ")) {
            on = true;
        }
        if on {
            lines.push(line);
            if line.contains("Action:") {
                break;
            }
        }
    }
    lines.join("\n")
}

/// The text after the first `start` and before the first `end` after it.
fn between<'a>(text: &'a str, start: &str, end: &str) -> &'a str {
    match text.split_once(start) {
        Some((_, rest)) => rest.split(end).next().unwrap_or_default(),
        None => "",
    }
}

/// The Affected Medications section of a recommendations text (the whole
/// text when it has none).
fn medications(text: &str) -> &str {
    let rest = text
        .split_once("Affected Medications:")
        .map_or(text, |(_, r)| r);
    rest.split("Calls From Other Tools")
        .next()
        .unwrap_or_default()
}

fn genes_named(rows: &[Vec<String>], status: &str) -> Vec<String> {
    let mut g: Vec<String> = rows
        .iter()
        .filter(|r| r[3] == status)
        .map(|r| r[0].clone())
        .collect();
    g.sort();
    g
}

// --- 1. the real reports, one per PharmCAT release -------------------------------

#[test]
fn real_hg002_reports_parse_and_ambiguous_genes_stay_ambiguous() {
    let work = Scratch::new("real");
    for (version, ambiguous) in REAL_REPORTS {
        let real = fixture(&format!("report-{version}.json"));
        let tag = format!("{version}:");
        let r = plain(&work, &format!("HG002_{version}"), &real);
        let Value::Object(m) = &real else {
            panic!("not an object")
        };
        assert_eq!(
            m.get("pharmcatVersion"),
            Some(&s(version)),
            "{tag} the report is PharmCAT {version}"
        );
        assert_eq!(r.rc, 0, "{tag} cpic-report on the HG002 report exits 0");
        assert!(
            r.rows.len() >= 10,
            "{tag} the HG002 report yields genes: {}",
            r.rows.len()
        );
        assert!(
            r.rows.iter().any(|r| r[3] != "not called"),
            "{tag} the HG002 report calls at least one gene"
        );
        let mut names: Vec<&str> = r.rows.iter().map(|r| r[0].as_str()).collect();
        names.sort_unstable();
        names.dedup();
        assert!(
            names.len() == r.rows.len() && r.rows.iter().all(|r| r.len() == 4),
            "{tag} one phenotypes row per gene, with a Status"
        );
        let mut several: Vec<(String, usize)> = r
            .rows
            .iter()
            .filter(|r| r[1].contains(" possible)"))
            .map(|r| {
                let (_, n) = r[1].rsplit_once("(1 of ").unwrap();
                (
                    r[0].clone(),
                    n.split_whitespace().next().unwrap().parse().unwrap(),
                )
            })
            .collect();
        several.sort();
        let want: Vec<(String, usize)> = ambiguous
            .iter()
            .map(|(g, n)| ((*g).to_owned(), *n))
            .collect();
        assert_eq!(
            several, want,
            "{tag} the genes with more than one possible diplotype are CYP2B6, CYP2C19"
        );
        for (gene, n) in ambiguous {
            let g = row(&r.rows, gene);
            assert_eq!(
                g.map(|g| g[3].as_str()),
                Some("ambiguous"),
                "{tag} HG002 {gene}: {n} possible diplotypes, status ambiguous"
            );
            assert!(
                g.is_some_and(|g| g[2].starts_with("ambiguous: ")
                    && g[1].contains(&format!("(1 of {n} possible)"))),
                "{tag} HG002 {gene}: not reported as its first diplotype's phenotype: {g:?}"
            );
            assert!(
                r.rec.contains("Genes With More Than One Possible Result:")
                    && r.rec
                        .contains(&format!("  {gene} -- {n} possible diplotypes")),
                "{tag} the recommendations list {gene} with the ambiguous genes"
            );
        }
    }
}

#[test]
fn a_report_with_no_gene_fails_loud() {
    let work = Scratch::new("empty");
    let real = fixture("report-3.4.0.json");
    let names: Vec<String> = match &real {
        Value::Object(m) => match m.get("genes") {
            Some(Value::Object(g)) => g.keys().cloned().collect(),
            _ => panic!("no genes"),
        },
        _ => panic!("not an object"),
    };
    let stripped = with(
        real,
        "genes",
        Value::Object(
            names
                .iter()
                .map(|g| (g.clone(), o(vec![("geneSymbol", s(g))])))
                .collect(),
        ),
    );
    let r = plain(&work, "EMPTY", &stripped);
    assert_eq!(r.rc, 1, "a report whose genes carry no diplotype exits 1");
    assert!(
        r.rec.contains("PARSING FAILED")
            && !r.rec.contains("all genes were successfully called")
            && r.rows.is_empty(),
        "and writes PARSING FAILED, never an all-clear"
    );
    let r = plain(&work, "NOKEYS", &o(vec![("metadata", o(vec![]))]));
    assert!(
        r.rc == 1 && r.rec.contains("PARSING FAILED"),
        "a report with no gene section at all exits 1"
    );
}

// --- 2. PharmCAT's example report and the drug-list fallbacks -----------------------

const UNCLASSIFIED: [&str; 7] = [
    "CACNA1S", "CFTR", "CYP4F2", "IFNL3", "MT-RNR1", "RYR1", "VKORC1",
];

#[test]
fn the_example_report_has_three_non_normal_and_seven_unclassified_genes() {
    let work = Scratch::new("example");
    let example = fixture("pharmcat-docs-example.json");
    let r = plain(&work, "EXAMPLE", &example);
    assert_eq!(r.rc, 0, "cpic-report on the example report exits 0");
    assert_eq!(
        genes_named(&r.rows, "non-normal"),
        ["CYP2D6", "HLA-B", "NAT2"],
        "example report: 3 genes with a non-normal phenotype (CYP2D6, HLA-B, NAT2)"
    );
    assert_eq!(
        genes_named(&r.rows, "unclassified"),
        UNCLASSIFIED,
        "example report: n/a, 'no phenotype assigned', a genotype and 'Uncertain Susceptibility' are unclassified"
    );
    let after = r.rec.split("Affected Medications:").nth(1).unwrap();
    let meds = after
        .split_once("\n\n")
        .unwrap()
        .1
        .split("Called Genes Without")
        .next()
        .unwrap();
    let mut listed: Vec<&str> = meds
        .lines()
        .filter(|l| l.starts_with("  ") && l.contains(" -- ") && !l.starts_with("   "))
        .map(|l| l.split(" -- ").next().unwrap().trim())
        .collect();
    listed.sort_unstable();
    assert_eq!(
        listed,
        ["CYP2D6", "HLA-B", "NAT2"],
        "example: the medications section lists only the 3 non-normal genes"
    );
    let own = r
        .rec
        .rsplit("Called Genes Without a Function Phenotype:")
        .next()
        .unwrap()
        .split("Uncallable Genes:")
        .next()
        .unwrap();
    assert!(
        UNCLASSIFIED
            .iter()
            .all(|g| own.contains(&format!("  {g} -- ")))
            && own.contains("  VKORC1 -- -1639 GG ("),
        "example: the unclassified genes are listed in their own section, VKORC1 with its genotype"
    );
    assert_eq!(
        genes_named(&r.rows, "unclassified"),
        UNCLASSIFIED,
        "example: the phenotypes table says unclassified for them"
    );
    for gene in ["CYP2D6", "HLA-B", "NAT2"] {
        let b = block(&r.rec, gene);
        assert!(
            b.contains("    Drugs") && !b.contains("not in the drug table"),
            "example {gene}: a drug list: {b}"
        );
    }
    assert!(
        block(&r.rec, "CYP2D6").contains(
            "- codeine [CPIC, Moderate]: Use codeine label recommended age- or weight-specific dosing. \
             If no response and opioid use is warranted, consider a non-tramadol opioid."
        ),
        "example CYP2D6 *1/*3 (Intermediate): codeine's recommendation for that diplotype"
    );

    let mut ga = example;
    let first = match at(at(at(&mut ga, "genes"), "VKORC1"), "sourceDiplotypes") {
        Value::Array(v) => &mut v[0],
        other => panic!("{other:?}"),
    };
    *at(first, "label") = s("rs9923231 reference (C)/rs9923231 variant (T)");
    *at(first, "phenotypes") = a(vec![s("-1639 GA")]);
    let r = plain(&work, "EXAMPLE_GA", &ga);
    assert_eq!(
        cell(&r.rows, "VKORC1", 3),
        "non-normal",
        "example with VKORC1 -1639 GA: VKORC1 counts as non-normal"
    );
    let b = block(&r.rec, "VKORC1");
    assert!(
        r.rc == 0 && b.contains("    Drugs") && b.contains("warfarin"),
        "example with VKORC1 -1639 GA: its drug list names warfarin: {b}"
    );
}

#[test]
fn drug_list_fallbacks() {
    let work = Scratch::new("fallback");
    let data = with(
        genes(vec![
            ("CYP2C19", dip("*2", "*2", "Poor Metabolizer")),
            (
                "ABCG2",
                with(
                    dip(
                        "rs2231142 variant (T)",
                        "rs2231142 variant (T)",
                        "Poor Function",
                    ),
                    "relatedDrugs",
                    a(vec![o(vec![("name", s("rosuvastatin"))])]),
                ),
            ),
            ("NEWGENE", dip("*1", "*7", "Decreased Function")),
            ("CYP4F2", dip("*1", "*3", "Decreased Function")),
            ("SLCO1B1", dip("*1", "*1", "Normal Function")),
        ]),
        "drugs",
        cpic_drugs(
            "warfarin",
            vec![o(vec![
                ("phenotypes", o(vec![("CYP4F2", s("Decreased Function"))])),
                ("classification", s("Optional")),
                ("drugRecommendation", s("Increase the dose by 5-10%.")),
            ])],
        ),
    );
    let r = plain(&work, "FALLBACK", &data);
    assert_eq!(r.rc, 0, "fallbacks: exits 0");
    assert!(
        block(&r.rec, "CYP4F2")
            .contains("- warfarin [CPIC, Optional]: Increase the dose by 5-10%."),
        "CYP4F2: the report's own guidance (warfarin, CPIC Optional)"
    );
    let b = block(&r.rec, "ABCG2");
    assert!(
        b.contains("Drugs PharmCAT links to ABCG2") && b.contains("rosuvastatin"),
        "ABCG2: relatedDrugs when the report has no guidance for it"
    );
    let b = block(&r.rec, "CYP2C19");
    assert!(
        b.contains("pipeline fallback table") && b.contains("clopidogrel"),
        "CYP2C19: the static table when the report names no drug"
    );
    assert!(
        block(&r.rec, "NEWGENE")
            .contains("NEWGENE has a non-normal phenotype but is not in the drug table"),
        "NEWGENE: a line saying it is not in the drug table, not silence"
    );
    assert_eq!(
        block(&r.rec, "SLCO1B1"),
        "",
        "SLCO1B1 normal: no medications block"
    );
}

/// PharmCAT lists an annotation for every diplotype the sample may have:
/// the recommendation shown must be the one for the called diplotype.
#[test]
fn the_recommendation_is_the_called_diplotypes() {
    let work = Scratch::new("match");
    let ann = |label: &str, phen: &str, text: &str| {
        o(vec![
            ("phenotypes", o(vec![("CYP2C19", s(phen))])),
            ("classification", s("Strong")),
            ("drugRecommendation", s(text)),
            (
                "genotypes",
                a(vec![o(vec![(
                    "diplotypes",
                    a(vec![o(vec![("gene", s("CYP2C19")), ("label", s(label))])]),
                )])]),
            ),
        ])
    };
    let data = with(
        genes(vec![(
            "CYP2C19",
            dip("*1", "*2", "Intermediate Metabolizer"),
        )]),
        "drugs",
        cpic_drugs(
            "clopidogrel",
            vec![
                ann("*1/*1", "Normal Metabolizer", "Standard dose."),
                ann(
                    "*1/*2",
                    "Intermediate Metabolizer",
                    "Use an alternative antiplatelet.",
                ),
                ann("*2/*2", "Poor Metabolizer", "Avoid clopidogrel."),
            ],
        ),
    );
    let r = plain(&work, "MATCH", &data);
    let b = block(&r.rec, "CYP2C19");
    assert!(
        b.contains("- clopidogrel [CPIC, Strong]: Use an alternative antiplatelet.")
            && !b.contains("Standard dose.")
            && !b.contains("Avoid clopidogrel."),
        "the recommendation is the called diplotype's (*1/*2), not another possible one's: {b}"
    );

    // The phenotype fallback (annotations naming no diplotype) matches the
    // called phenotype exactly.
    let pann = |phen: &str, text: &str| {
        o(vec![
            ("phenotypes", o(vec![("CYP2C19", s(phen))])),
            ("classification", s("Strong")),
            ("drugRecommendation", s(text)),
        ])
    };
    let data = with(
        genes(vec![(
            "CYP2C19",
            dip("*1", "*9", "Likely Intermediate Metabolizer"),
        )]),
        "drugs",
        cpic_drugs(
            "clopidogrel",
            vec![
                pann("Intermediate Metabolizer", "Guidance for IM."),
                pann("Likely Intermediate Metabolizer", "Guidance for likely IM."),
            ],
        ),
    );
    let r = plain(&work, "PHENMATCH", &data);
    let b = block(&r.rec, "CYP2C19");
    assert!(
        b.contains("Guidance for likely IM.") && !b.contains("Guidance for IM."),
        "phenotype fallback: the exact phenotype's guidance only, not a substring match: {b}"
    );
}

// --- 3. layouts ----------------------------------------------------------------

#[test]
fn flat_nested_list_and_recommendation_layouts_parse() {
    let work = Scratch::new("layouts");
    let g2 = || {
        o(vec![
            ("CYP2C19", dip("*1", "*2", "Intermediate Metabolizer")),
            ("CYP2D6", dip("*1", "*1", "Normal Metabolizer")),
        ])
    };
    let cases = [
        ("flat 3.x", o(vec![("genes", g2())]), 2),
        (
            "nested 2.x",
            o(vec![("genes", o(vec![("CPIC", g2()), ("DPWG", g2())]))]),
            2,
        ),
        (
            "list",
            o(vec![(
                "genes",
                a(vec![with(
                    dip("*1", "*2", "Intermediate Metabolizer"),
                    "gene",
                    s("CYP2C19"),
                )]),
            )]),
            1,
        ),
        (
            "recommendationDiplotypes",
            genes(vec![(
                "DPYD",
                o(vec![(
                    "recommendationDiplotypes",
                    a(vec![o(vec![
                        ("label", s("c.1905+1G>A/Reference")),
                        ("phenotypes", a(vec![s("Intermediate Metabolizer")])),
                    ])]),
                )]),
            )]),
            1,
        ),
    ];
    for (label, d, want) in cases {
        let r = plain(&work, &format!("LAYOUT{}", label.len()), &d);
        assert!(
            r.rc == 0 && r.rows.len() == want,
            "layout {label}: {want} gene(s): {} {:?}",
            r.rc,
            r.rows
        );
    }
}

// --- 4. PharmCAT no result, a pypgx call -------------------------------------------

#[test]
fn pypgx_call_where_pharmcat_has_no_result() {
    let work = Scratch::new("pypgx");
    let data = with(
        genes(vec![
            (
                "CYP2D6",
                with(
                    dip("Unknown", "Unknown", "No Result"),
                    "relatedDrugs",
                    a(vec![
                        o(vec![("name", s("codeine"))]),
                        o(vec![("name", s("tramadol"))]),
                    ]),
                ),
            ),
            ("CYP2C9", dip("Unknown", "Unknown", "No Result")),
            ("CYP2C19", dip("*1", "*1", "Normal Metabolizer")),
        ]),
        // a two-gene annotation that names CYP2D6 for one drug only
        "drugs",
        cpic_drugs(
            "amitriptyline",
            vec![o(vec![
                (
                    "phenotypes",
                    o(vec![
                        ("CYP2C19", s("Normal Metabolizer")),
                        ("CYP2D6", s("No Result")),
                    ]),
                ),
                ("drugRecommendation", s("x")),
                ("classification", s("Optional")),
            ])],
        ),
    );
    let r = run_report(
        &work,
        "PYPGX",
        &data,
        &Opts {
            pypgx: Some(vec![
                ("CYP2D6", "*1/*4", "Intermediate Metabolizer"),
                ("CYP2C9", "*1/*1", "Normal Metabolizer"),
                ("CYP2C19", "*1/*1", "Normal Metabolizer"),
            ]),
            ..Opts::default()
        },
    );
    assert_eq!(r.rc, 0, "pypgx: exits 0");
    assert!(
        r.rec.contains(
            "WARNING: PharmCAT has no result for CYP2D6, but pypgx (step 32) called *1/*4"
        ),
        "pypgx: the warning names CYP2D6 and pypgx's call"
    );
    assert!(
        r.rec
            .contains("Drugs affected by CYP2D6: codeine, tramadol"),
        "pypgx: the warning lists every drug PharmCAT links to CYP2D6, not only the one another gene matched"
    );
    assert!(
        r.rec.contains(
            "NOTE: PharmCAT has no result for CYP2C9; pypgx (step 32) called *1/*1 (Normal Metabolizer)."
        ) && !r.rec.contains("WARNING: PharmCAT has no result for CYP2C9"),
        "pypgx: a normal pypgx call is a NOTE, not a WARNING"
    );
    assert!(
        r.rec.contains("CYP2D6 -- No Result (not callable"),
        "pypgx: CYP2D6 is in the uncallable list"
    );
    assert!(
        r.comp.contains(&strs(&[
            "CYP2D6",
            "Unknown/Unknown",
            "*1/*4",
            "pypgx only",
            "pypgx only",
            "No Result",
            "Intermediate Metabolizer"
        ])) && r.comp.contains(&strs(&[
            "CYP2C19",
            "*1/*1",
            "*1/*1",
            "Yes",
            "both",
            "Normal Metabolizer",
            "Normal Metabolizer"
        ])),
        "pypgx: comparison marks CYP2D6 'pypgx only' and CYP2C19 concordant: {:?}",
        r.comp
    );
}

/// An ambiguous PharmCAT gene is not compared as if it were its first
/// diplotype.
#[test]
fn an_ambiguous_gene_is_compared_as_ambiguous() {
    let work = Scratch::new("ambig");
    let one = |a2: &str, phen: &str| {
        o(vec![
            ("allele1", o(vec![("name", s("*1"))])),
            ("allele2", o(vec![("name", s(a2))])),
            ("label", s(&format!("*1/{a2}"))),
            ("phenotypes", a(vec![s(phen)])),
        ])
    };
    let data = genes(vec![(
        "CYP2C19",
        o(vec![(
            "sourceDiplotypes",
            a(vec![
                one("*1", "Normal Metabolizer"),
                one("*2", "Intermediate Metabolizer"),
            ]),
        )]),
    )]);
    let r = run_report(
        &work,
        "AMBIG",
        &data,
        &Opts {
            pypgx: Some(vec![("CYP2C19", "*1/*2", "Intermediate Metabolizer")]),
            ..Opts::default()
        },
    );
    assert_eq!(
        r.rows.first().map(|r| r[3].as_str()),
        Some("ambiguous"),
        "ambiguous: the phenotypes row says ambiguous"
    );
    assert!(
        r.comp
            .first()
            .is_some_and(|c| c[1] == "ambiguous (2 possible diplotypes)" && c[3] == "pypgx only")
            && r.rec.contains(
                "WARNING: PharmCAT has no single result (2 possible diplotypes) for CYP2C19"
            ),
        "ambiguous: compared as ambiguous, pypgx's call warned about: {:?}",
        r.comp
    );
}

/// An unreadable pypgx summary: the CPIC report is still written, it says
/// the comparison could not be made, and no comparison file is written.
#[test]
fn an_unreadable_pypgx_summary_is_reported_not_fatal() {
    let work = Scratch::new("badpypgx");
    let rep = work.file(
        "report.json",
        &dumps(
            &genes(vec![("CYP2C19", dip("*2", "*2", "Poor Metabolizer"))]),
            &Dump::default(),
        )
        .unwrap(),
    );
    let bad = work.path("pypgx.tsv");
    fs::write(&bad, b"Gene\tDiplotype\n\xff\xfe\x00broken\n").unwrap();
    // a comparison left by an earlier run
    let comp = work.file(
        "comparison.tsv",
        "Gene\tPharmCAT_diplotype\tpypgx_diplotype\tMatch\tCalled_by\nCYP2C19\t*1/*1\t*1/*1\tYes\tboth\n",
    );
    let outdir = work.0.to_string_lossy().into_owned();
    let rc = og(&[
        "cpic-report",
        "--sample",
        "BADPYPGX",
        "--report",
        &rep,
        "--outdir",
        &outdir,
        "--pypgx",
        &bad,
        "--comparison",
        &comp,
    ]);
    let rec =
        fs::read_to_string(work.path("BADPYPGX_cpic_recommendations.txt")).unwrap_or_default();
    assert!(
        rc == 0 && rec.contains("CYP2C19 -- Poor"),
        "unreadable pypgx summary: exits 0 with the CPIC report written: {rc}"
    );
    assert!(
        rec.contains("Could not read the pypgx summary") && !Path::new(&comp).exists(),
        "unreadable pypgx summary: the report says no comparison was made, and an earlier one is removed"
    );
}

// --- 5. Indeterminate at gene level, a phenotype in the drug annotations ------------

const INDET_NOTE: &str = "[gene level Indeterminate; phenotype from PharmCAT's drug annotations]";

fn iann(phen: &str, label: &str, rec: Option<&str>) -> Value {
    o(vec![
        ("phenotypes", o(vec![("UGT1A1", s(phen))])),
        ("classification", s("Strong")),
        (
            "drugRecommendation",
            s(&rec.map_or_else(|| format!("Guidance for {phen}."), str::to_owned)),
        ),
        (
            "genotypes",
            a(vec![o(vec![(
                "diplotypes",
                a(vec![o(vec![("gene", s("UGT1A1")), ("label", s(label))])]),
            )])]),
        ),
    ])
}

fn ugt1a1_report(gene: Value, anns: Option<Vec<Value>>) -> Value {
    let g = genes(vec![("UGT1A1", gene)]);
    match anns {
        Some(anns) => with(g, "drugs", cpic_drugs("atazanavir", anns)),
        None => g,
    }
}

#[test]
fn an_indeterminate_gene_takes_one_drug_level_phenotype() {
    let work = Scratch::new("indet");
    let gene = || {
        with(
            dip("*1", "*80", "Indeterminate"),
            "relatedDrugs",
            a(vec![o(vec![("name", s("atazanavir"))])]),
        )
    };
    let r = plain(
        &work,
        "INDET",
        &ugt1a1_report(
            gene(),
            Some(vec![iann("Intermediate Metabolizer", "*1/*80", None)]),
        ),
    );
    assert_eq!(
        r.rc, 0,
        "Indeterminate gene, one drug-level phenotype: exits 0"
    );
    let note = format!("Intermediate Metabolizer {INDET_NOTE}");
    assert_eq!(
        r.rows,
        [strs(&["UGT1A1", "*1/*80", &note, "non-normal"])],
        "Indeterminate gene: the table shows the drug-level phenotype with its note, and its status"
    );
    assert!(
        r.rec
            .split("Affected Medications:")
            .next()
            .unwrap()
            .contains(&note),
        "Indeterminate gene: the gene results show the note"
    );
    assert!(
        r.rec.contains(&format!("  UGT1A1 -- {note}:"))
            && block(&r.rec, "UGT1A1")
                .contains("- atazanavir [CPIC, Strong]: Guidance for Intermediate Metabolizer."),
        "Indeterminate gene: listed with its drugs, the guidance for the called diplotype"
    );
    assert!(
        !r.rec.contains("Called Genes Without a Function Phenotype:"),
        "Indeterminate gene: not in the section of genes without a function phenotype"
    );
    let cases: [(&str, Option<Vec<Value>>); 4] = [
        ("no drug annotation", None),
        (
            "two drug annotations that disagree",
            Some(vec![
                iann("Intermediate Metabolizer", "*1/*80", None),
                iann("Poor Metabolizer", "*1/*80", None),
            ]),
        ),
        (
            "a drug annotation that is Indeterminate too",
            Some(vec![iann("Indeterminate", "*1/*80", None)]),
        ),
        (
            "one drug annotation with a phenotype and one Indeterminate",
            Some(vec![
                iann("Intermediate Metabolizer", "*1/*80", None),
                iann("Indeterminate", "*1/*80", None),
            ]),
        ),
    ];
    for (label, anns) in cases {
        let r = plain(
            &work,
            &format!("INDET{}", label.len()),
            &ugt1a1_report(gene(), anns),
        );
        assert!(
            r.rc == 0
                && r.rows == [strs(&["UGT1A1", "*1/*80", "Indeterminate", "unclassified"])]
                && !r.rec.contains(INDET_NOTE),
            "Indeterminate gene, {label}: stays Indeterminate and unclassified, no note: {:?}",
            r.rows
        );
    }
    let r = plain(
        &work,
        "INDETNORMAL",
        &ugt1a1_report(
            gene(),
            Some(vec![iann("Normal Metabolizer", "*1/*80", None)]),
        ),
    );
    assert_eq!(
        r.rows,
        [strs(&[
            "UGT1A1",
            "*1/*80",
            &format!("Normal Metabolizer {INDET_NOTE}"),
            "normal"
        ])],
        "Indeterminate gene, drug level Normal Metabolizer: normal, with the note"
    );
}

#[test]
fn an_indeterminate_gene_with_two_possible_diplotypes() {
    let work = Scratch::new("indettwo");
    let two = || {
        let mut d = dip("*1", "*80", "Indeterminate");
        let extra = match dip("*80", "*80", "Indeterminate") {
            Value::Object(mut m) => m.shift_remove("sourceDiplotypes").unwrap(),
            _ => unreachable!(),
        };
        if let (Value::Array(v), Value::Array(x)) = (at(&mut d, "sourceDiplotypes"), extra) {
            v.extend(x);
        }
        d
    };
    let undecided = strs(&[
        "UGT1A1",
        "*1/*80 (1 of 2 possible)",
        "Indeterminate",
        "unclassified",
    ]);
    // two possible diplotypes, both Indeterminate, whose drug annotations
    // disagree: the gene is undecided and is not read as its first one
    let r = plain(
        &work,
        "INDETTWO",
        &ugt1a1_report(
            two(),
            Some(vec![
                iann("Intermediate Metabolizer", "*1/*80", None),
                iann("Normal Metabolizer", "*80/*80", None),
            ]),
        ),
    );
    assert!(
        r.rc == 0 && r.rows == [undecided.clone()] && !r.rec.contains(INDET_NOTE),
        "Indeterminate gene, two possible diplotypes whose drug annotations disagree: stays Indeterminate, no note: {:?}",
        r.rows
    );
    // the same phenotype for both, but different guidance: the drug list
    // would show the first diplotype's, so the gene stays undecided
    let r = plain(
        &work,
        "INDETTWOREC",
        &ugt1a1_report(
            two(),
            Some(vec![
                iann("Intermediate Metabolizer", "*1/*80", None),
                iann(
                    "Intermediate Metabolizer",
                    "*80/*80",
                    Some("Other guidance."),
                ),
            ]),
        ),
    );
    assert!(
        r.rc == 0 && r.rows == [undecided] && !r.rec.contains(INDET_NOTE),
        "Indeterminate gene, two possible diplotypes with one phenotype but different guidance: stays Indeterminate: {:?}",
        r.rows
    );
    // the same phenotype and the same guidance for both: it takes the phenotype
    let r = plain(
        &work,
        "INDETTWOSAME",
        &ugt1a1_report(
            two(),
            Some(vec![
                iann("Intermediate Metabolizer", "*1/*80", None),
                iann("Intermediate Metabolizer", "*80/*80", None),
            ]),
        ),
    );
    assert_eq!(
        r.rows,
        [strs(&[
            "UGT1A1",
            "*1/*80 (1 of 2 possible)",
            &format!("Intermediate Metabolizer {INDET_NOTE}"),
            "non-normal"
        ])],
        "Indeterminate gene, two possible diplotypes with the same phenotype and guidance: takes it, with the note"
    );
}

// --- 6. label-only differences, long text, single-position fallback ------------------

#[test]
fn single_position_labels_compare_by_meaning() {
    let work = Scratch::new("labels");
    let vk = "rs9923231 reference (C)/rs9923231 variant (T)";
    let one = |label: &str, phen: &str| {
        o(vec![(
            "sourceDiplotypes",
            a(vec![o(vec![
                ("label", s(label)),
                ("phenotypes", a(vec![s(phen)])),
            ])]),
        )])
    };
    let data = genes(vec![
        ("VKORC1", one(vk, "-1639 GA")),
        (
            "IFNL3",
            one(
                "rs12979860 variant (T)/rs12979860 variant (T)",
                "unfavorable response genotype",
            ),
        ),
    ]);
    let r = run_report(
        &work,
        "LABELS",
        &data,
        &Opts {
            pypgx: Some(vec![
                ("VKORC1", "Reference/rs9923231", "no pypgx phenotype rule"),
                ("IFNL3", "Reference/rs12979860", "Unfavorable Response"),
            ]),
            ..Opts::default()
        },
    );
    assert!(
        r.comp.contains(&strs(&[
            "VKORC1",
            vk,
            "Reference/rs9923231",
            "Yes",
            "both",
            "-1639 GA",
            "no pypgx phenotype rule"
        ])),
        "VKORC1 'rsX reference (C)/rsX variant (T)' against pypgx 'Reference/rsX' is a Match: {:?}",
        r.comp
    );
    let ifnl3: Vec<&str> = r
        .comp
        .iter()
        .filter(|c| c[0] == "IFNL3")
        .map(|c| c[3].as_str())
        .collect();
    assert_eq!(
        ifnl3,
        ["No"],
        "IFNL3 'rsX variant (T)/rsX variant (T)' against pypgx 'Reference/rsX' is still a conflict"
    );
    assert!(
        r.comp
            .iter()
            .all(|c| !c.join("\t").contains("Indeterminate")),
        "the comparison shows PharmCAT's phenotype beside pypgx's 'no pypgx phenotype rule', no Indeterminate"
    );
    let text = fs::read_to_string(work.0.join("LABELS/comparison.tsv")).unwrap();
    let head: Vec<&str> = text.lines().next().unwrap().split('\t').collect();
    assert_eq!(
        head,
        [
            "Gene",
            "PharmCAT_diplotype",
            "pypgx_diplotype",
            "Match",
            "Called_by",
            "PharmCAT_phenotype",
            "pypgx_phenotype"
        ],
        "the comparison header ends with the two phenotype columns"
    );
}

#[test]
fn a_long_recommendation_prints_whole() {
    let work = Scratch::new("long");
    let long = format!(
        "{}If plasma concentrations are not available, consider a lower dose of a tricyclic in this patient.",
        "Initiate therapy with 50% of the starting dose followed by titration of dose based on \
         toxicity or therapeutic drug monitoring. "
            .repeat(4)
    );
    let data = with(
        genes(vec![("CYP2C19", dip("*2", "*2", "Poor Metabolizer"))]),
        "drugs",
        cpic_drugs(
            "amitriptyline",
            vec![o(vec![
                ("phenotypes", o(vec![("CYP2C19", s("Poor Metabolizer"))])),
                ("classification", s("Optional")),
                ("drugRecommendation", s(&long)),
            ])],
        ),
    );
    let r = plain(&work, "LONG", &data);
    assert!(
        long.len() > 400
            && r.rec
                .contains(&format!("- amitriptyline [CPIC, Optional]: {long}"))
            && !block(&r.rec, "CYP2C19").contains("..."),
        "a recommendation of {} characters prints whole",
        long.len()
    );
}

/// A gene PharmCAT leaves Unknown/Unknown, with the genotypes it read.
/// The call comes before the rsID: an rsID next to a genotype reads as a
/// leak.
fn variants(call: Option<&str>, rsid: &str) -> Value {
    o(vec![
        (
            "sourceDiplotypes",
            a(vec![o(vec![
                ("allele1", o(vec![("name", s("Unknown"))])),
                ("allele2", o(vec![("name", s("Unknown"))])),
                ("label", s("Unknown/Unknown")),
                ("phenotypes", a(vec![])),
            ])]),
        ),
        (
            "variants",
            a(vec![
                o(vec![
                    ("position", Value::Int(1)),
                    ("dbSnpId", s("rs0000001")),
                    ("call", s("G/A")),
                ]),
                o(vec![
                    ("position", Value::Int(2)),
                    ("dbSnpId", s(rsid)),
                    ("call", call.map_or(Value::Null, s)),
                    ("referenceAllele", s("C")),
                ]),
            ]),
        ),
    ])
}

const FB: &str = "Single-Position Fallback (not a PharmCAT call):";

#[test]
fn cyp4f2_and_vkorc1_single_position_fallback() {
    let work = Scratch::new("fallback1");
    let r = plain(
        &work,
        "CYP4F2HET",
        &genes(vec![("CYP4F2", variants(Some("C/T"), "rs2108622"))]),
    );
    let tail = r.rec.rsplit(FB).next().unwrap();
    assert!(
        r.rc == 0
            && r.rec.contains(FB)
            && tail.contains(
                "WARNING: PharmCAT has no result for CYP4F2; single-position fallback, not a PharmCAT call: \
                 rs2108622 reads C/T, one copy of the *3 variant (T, p.V433M)."
            )
            && tail.contains("Warfarin (CPIC): the CPIC warfarin dosing algorithm")
            && tail.contains("5-10%"),
        "CYP4F2 No Result with C/T at rs2108622: the fallback names the genotype and the *3 rule, labelled"
    );
    assert_eq!(
        r.rows,
        [strs(&[
            "CYP4F2",
            "Unknown/Unknown",
            "No Result",
            "not called"
        ])],
        "CYP4F2 fallback: PharmCAT's own result stays not called in the phenotypes table"
    );
    let r = plain(
        &work,
        "CYP4F2HOM",
        &genes(vec![("CYP4F2", variants(Some("T|T"), "rs2108622"))]),
    );
    assert!(
        r.rec
            .contains("rs2108622 reads T|T, two copies of the *3 variant"),
        "CYP4F2 T|T (phased): two copies"
    );
    let r = plain(
        &work,
        "CYP4F2REF",
        &genes(vec![("CYP4F2", variants(Some("C/C"), "rs2108622"))]),
    );
    assert!(
        r.rec.contains(
            "rs2108622 reads C/C, no copy of the *3 variant (T, p.V433M). The CPIC warfarin rule for this position does not apply."
        ) && !r.rec.contains("WARNING: PharmCAT has no result for CYP4F2"),
        "CYP4F2 C/C: no copy, the rule does not apply, and no WARNING"
    );
    let r = plain(
        &work,
        "CYP4F2MISS",
        &genes(vec![("CYP4F2", variants(None, "rs2108622"))]),
    );
    assert!(
        r.rec.contains(
            "CYP4F2 -- PharmCAT has no result; rs2108622 has no genotype in the VCF PharmCAT read: not callable."
        ) && !r.rec.rsplit(FB).next().unwrap().contains("WARNING"),
        "CYP4F2 with no genotype at rs2108622: not callable"
    );
    let r = plain(
        &work,
        "VKORC1HET",
        &genes(vec![("VKORC1", variants(Some("C/T"), "rs9923231"))]),
    );
    assert!(
        r.rec.contains(
            "WARNING: PharmCAT has no result for VKORC1; single-position fallback, not a PharmCAT call: \
             rs9923231 reads C/T, -1639 GA."
        ) && r.rec.contains(
            "Warfarin (CPIC): the CPIC warfarin guideline doses by a pharmacogenetic algorithm \
             that takes VKORC1 -1639G>A with CYP2C9; each A allele lowers the expected dose."
        ),
        "VKORC1 No Result with C/T at rs9923231: a WARNING reads -1639 GA, with the VKORC1 warfarin rule"
    );
    let called = genes(vec![(
        "CYP4F2",
        with(
            dip("*1", "*3", "Decreased Function"),
            "variants",
            a(vec![o(vec![
                ("dbSnpId", s("rs2108622")),
                ("call", s("C/T")),
            ])]),
        ),
    )]);
    let r = plain(&work, "CYP4F2CALLED", &called);
    assert!(
        !r.rec.contains(FB),
        "CYP4F2 called by PharmCAT: no fallback section"
    );
}

// --- 7. MT-RNR1 from step 20's chrM calls next to step 36's consensus ---------------

const CHILD: &str = "OG_PORT_TESTS_PGX_CHILD";
const MARK: &str = "\u{1e}og-child\u{1e}";

/// Not a test of its own: the child [`og_out`] starts runs `og pgx ARGV`
/// here, then exits with the verb's code. In the normal test run it
/// returns at once.
#[test]
fn og_child() {
    let Some(argv) = std::env::var_os(CHILD) else {
        return;
    };
    let argv: Vec<String> = argv
        .to_string_lossy()
        .split('\u{1f}')
        .map(str::to_owned)
        .collect();
    print!("{MARK}");
    eprint!("{MARK}");
    let code = og_pgx::main(&argv);
    let _ = std::io::stdout().flush();
    let _ = std::io::stderr().flush();
    std::process::exit(i32::from(code));
}

/// `og pgx ARGV` in a child copy of this test binary: (exit code, stdout).
fn og_out(argv: &[&str]) -> (i32, String) {
    let child = Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "og_child",
            "--nocapture",
            "--test-threads=1",
            "-q",
        ])
        .env(CHILD, argv.join("\u{1f}"))
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    let done = child.wait_with_output().unwrap();
    let code = done.status.code().unwrap_or(-1);
    let out = String::from_utf8_lossy(&done.stdout);
    match out.split_once(MARK) {
        Some((_, rest)) => (code, rest.to_owned()),
        None => panic!("the child never reached the verb: {out}"),
    }
}

#[test]
fn an_mt_rnr1_outside_call_is_not_checked_against_step36() {
    let work = Scratch::new("mtrnr1");
    let outside = |a1: &str, a2: Option<&str>, phen: &str, drugs: &[&str]| {
        let g = match a2 {
            Some(a2) => dip(a1, a2, phen),
            None => o(vec![(
                "sourceDiplotypes",
                a(vec![o(vec![
                    ("allele1", o(vec![("name", s(a1))])),
                    ("allele2", Value::Null),
                    ("label", s(a1)),
                    ("phenotypes", a(vec![s(phen)])),
                ])]),
            )]),
        };
        let g = with(g, "callSource", s("OUTSIDE"));
        with(
            g,
            "relatedDrugs",
            a(drugs.iter().map(|n| o(vec![("name", s(n))])).collect()),
        )
    };
    let mt = outside(
        "m.1555A>G",
        None,
        "increased risk of aminoglycoside-induced hearing loss",
        &["amikacin", "gentamicin"],
    );
    let mito = genes(vec![
        ("MT-RNR1", mt.clone()),
        (
            "CYP2D6",
            outside("*1", Some("*4"), "Intermediate Metabolizer", &["codeine"]),
        ),
    ]);
    let r = run_report(
        &work,
        "MTRNR1CONS",
        &mito,
        &Opts {
            consensus: Some(vec![("CYP2D6", "*1/*1", "yes")]),
            ..Opts::default()
        },
    );
    assert_eq!(
        r.rc, 0,
        "MT-RNR1 outside call, consensus without it: exits 0"
    );
    assert!(
        r.rows.contains(&strs(&[
            "MT-RNR1",
            "m.1555A>G",
            "increased risk of aminoglycoside-induced hearing loss",
            "non-normal"
        ])),
        "MT-RNR1 outside call, consensus without it: stays called and non-normal: {:?}",
        r.rows
    );
    let b = block(&r.rec, "MT-RNR1");
    assert!(
        r.rec.contains("  MT-RNR1 -- ") && b.contains("amikacin") && b.contains("gentamicin"),
        "MT-RNR1 outside call, consensus without it: listed with its drugs: {b}"
    );
    assert!(
        !r.rec.contains("outside call MT-RNR1"),
        "MT-RNR1 outside call, consensus without it: no stale warning for MT-RNR1"
    );
    assert!(
        r.rec
            .contains("WARNING: PharmCAT's report has the outside call CYP2D6 *1/*4")
            && !medications(&r.rec).contains("  CYP2D6 -- "),
        "control: the CYP2D6 outside call the consensus passes as another diplotype is still held back"
    );
    let rep = work.0.join("MTRNR1CONS/report.json");
    let (rc, out) = og_out(&["step36-outside", "--report", &rep.to_string_lossy()]);
    assert!(
        rc == 0 && out.split_whitespace().collect::<Vec<_>>() == ["CYP2D6"],
        "step36-outside lists CYP2D6 and not MT-RNR1: {rc} {out:?}"
    );
    let only = work.file(
        "MTRNR1ONLY.json",
        &dumps(&genes(vec![("MT-RNR1", mt)]), &Dump::default()).unwrap(),
    );
    let (rc, out) = og_out(&["step36-outside", "--report", &only]);
    assert!(
        rc == 0 && out.is_empty(),
        "step36-outside lists nothing for a report whose only outside call is MT-RNR1: {rc} {out:?}"
    );
}

// --- 8. phenoconversion --------------------------------------------------------------

/// A synthetic report: CYP2D6 *1/*2 normal metabolizer (activity score
/// `nm_score`, none for `None`) and CYP2C19 *1/*1 normal metabolizer.
fn pheno_report(nm_score: Option<&str>) -> Value {
    let mut d6 = dip("*1", "*2", "Normal Metabolizer");
    if let Value::Array(v) = at(&mut d6, "sourceDiplotypes")
        && let Value::Object(m) = &mut v[0]
    {
        m.insert("activityScore".to_owned(), nm_score.map_or(Value::Null, s));
    }
    o(vec![
        ("pharmcatVersion", s("3.4.0")),
        (
            "genes",
            o(vec![
                ("CYP2D6", d6),
                ("CYP2C19", dip("*1", "*1", "Normal Metabolizer")),
            ]),
        ),
    ])
}

/// The phenoconversion section of the recommendations text.
fn section(rec: &str) -> &str {
    between(
        rec,
        "Medications and Phenoconversion (CYP2D6, CYP2C19, CYP2C9, CYP3A4, CYP3A5):",
        "Affected Medications:",
    )
}

fn meds<'a>(m: &'a str) -> Opts<'a> {
    Opts {
        medications: Some(m),
        ..Opts::default()
    }
}

fn meds_with<'a>(m: &'a str, table: String) -> Opts<'a> {
    Opts {
        medications: Some(m),
        inhibitors: Some(table),
        ..Opts::default()
    }
}

/// The non-comment lines of a shipped table: its header and its rows.
fn table_rows(path: &str) -> (Vec<String>, Vec<Vec<String>>) {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.lines().filter(|l| !l.starts_with('#'));
    let header = lines
        .next()
        .unwrap()
        .split('\t')
        .map(str::to_owned)
        .collect();
    let rows = lines
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect();
    (header, rows)
}

/// The inhibitor table with `from` replaced by `to` on every line.
fn edited_inhibitors(work: &Scratch, name: &str, from: &str, to: &str) -> String {
    let text = fs::read_to_string(inhibitors()).unwrap();
    work.file(name, &text.replace(from, to))
}

#[test]
fn inhibitor_table_has_the_shape_the_reader_reads() {
    let work = Scratch::new("inhibitors");
    let (header, rows) = table_rows(&inhibitors());
    assert_eq!(
        header,
        [
            "Drug",
            "Gene",
            "Strength",
            "Activity_score_factor",
            "Class_source",
            "Rule_source"
        ],
        "inhibitors: the header"
    );
    // The strengths each gene has rows for: the FDA's table lists no strong
    // CYP2C9 inhibitor, and the CYP3A rows are its strong inhibitors only.
    let strengths: [(&str, &[&str]); 5] = [
        ("CYP2D6", &["strong", "moderate"]),
        ("CYP2C19", &["strong", "moderate"]),
        ("CYP2C9", &["moderate"]),
        ("CYP3A4", &["strong"]),
        ("CYP3A5", &["strong"]),
    ];
    let factor = |g: &str, st: &str| -> Option<&'static str> {
        let known = strengths
            .iter()
            .any(|(k, sts)| *k == g && sts.contains(&st));
        known.then_some(match (g, st) {
            ("CYP2D6", "strong") => "0",
            ("CYP2D6", _) => "0.5",
            _ => "-",
        })
    };
    let bad: Vec<&Vec<String>> = rows
        .iter()
        .filter(|r| {
            r.len() != 6
                || factor(&r[1], &r[2]) != Some(r[3].as_str())
                || r[0] != r[0].trim().to_lowercase()
                || r[0].is_empty()
                || r[4] != FDA_TABLE
                || r[5]
                    != if r[1] == "CYP2D6" {
                        "https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/"
                    } else {
                        "-"
                    }
        })
        .collect();
    assert!(
        bad.is_empty(),
        "inhibitors: every row has a gene, a strength, CPIC's factor for it and its sources: {bad:?}"
    );
    let keys: Vec<(&str, &str)> = rows
        .iter()
        .map(|r| (r[0].as_str(), r[1].as_str()))
        .collect();
    let mut unique = keys.clone();
    unique.sort_unstable();
    unique.dedup();
    assert_eq!(
        unique.len(),
        keys.len(),
        "inhibitors: no drug twice for one gene"
    );
    let mut table_genes: Vec<&str> = strengths.iter().map(|(g, _)| *g).collect();
    table_genes.sort_unstable();
    let mut parser_genes = PHENOCONVERSION_GENES.to_vec();
    parser_genes.sort_unstable();
    assert_eq!(
        table_genes, parser_genes,
        "inhibitors: the genes of the table are the parser's"
    );
    for (gene, sts) in strengths {
        for st in sts {
            assert!(
                rows.iter().any(|r| r[1] == gene && r[2] == *st),
                "inhibitors: {gene} has {st} inhibitors"
            );
        }
    }
    let drugs_of = |g: &str| -> std::collections::BTreeSet<&str> {
        rows.iter()
            .filter(|r| r[1] == g)
            .map(|r| r[0].as_str())
            .collect()
    };
    assert_eq!(
        drugs_of("CYP3A4"),
        drugs_of("CYP3A5"),
        "inhibitors: CYP3A4 and CYP3A5 list the same strong inhibitors"
    );
    let r = run_report(
        &work,
        "PHENOREAD",
        &pheno_report(Some("2.0")),
        &meds("paroxetine\n"),
    );
    let sec = section(&r.rec);
    assert!(
        !sec.contains("Could not read") && sec.contains("phenoconverted by paroxetine"),
        "inhibitors: the parser reads the table: {sec}"
    );
    let boosted = [
        "elvitegravir",
        "indinavir",
        "lopinavir",
        "paritaprevir",
        "saquinavir",
        "tipranavir",
    ];
    assert!(
        boosted.iter().all(|d| ["CYP3A4", "CYP3A5"]
            .iter()
            .all(|g| keys.contains(&(*d, *g)))),
        "inhibitors: each ritonavir-boosted partner the FDA classes as strong has a CYP3A4 and a CYP3A5 row"
    );
}

#[test]
fn a_strong_cyp2d6_inhibitor_phenoconverts_and_cyp2c19_gets_a_note() {
    let work = Scratch::new("pheno");
    let r = run_report(
        &work,
        "PHENO",
        &pheno_report(Some("2.0")),
        &meds("Paroxetine 20 mg daily\nfluconazole\nnotadrug\n# a comment\n\n"),
    );
    assert_eq!(r.rc, 0, "phenoconversion: exits 0");
    assert_eq!(
        cell(&r.rows, "CYP2D6", 4),
        "Poor Metabolizer (phenoconverted by paroxetine)",
        "phenoconversion: strong CYP2D6 inhibitor turns the NM into 'Poor Metabolizer (phenoconverted by paroxetine)'"
    );
    assert_eq!(
        row(&r.rows, "CYP2D6").map(|r| r[2..4].to_vec()),
        Some(strs(&["Normal Metabolizer", "normal"])),
        "phenoconversion: the genotype phenotype is kept beside it"
    );
    let sec = section(&r.rec);
    assert!(
        sec.contains(
            "WARNING: CYP2D6 genotype phenotype Normal Metabolizer; phenoconverted phenotype \
             Poor Metabolizer (phenoconverted by paroxetine)"
        ),
        "phenoconversion: a WARNING line with both phenotypes: {sec}"
    );
    assert!(
        sec.contains(
            "the CYP2D6 activity score is adjusted to 0 and the predicted phenotype is a poor metabolizer"
        ),
        "phenoconversion: CPIC's wording for a strong inhibitor"
    );
    assert!(
        cell(&r.rows, "CYP2C19", 4).starts_with("not computed: fluconazole")
            && sec.contains(
                "CYP2C19 genotype phenotype Normal Metabolizer; on fluconazole (strong CYP2C19 inhibitor)"
            ),
        "phenoconversion: CYP2C19 inhibitor gets a note and no converted phenotype"
    );
    assert!(
        sec.contains("Not in the inhibitor table: notadrug."),
        "phenoconversion: a name not in the table is listed, not dropped"
    );
    assert!(
        sec.contains("an aid for your clinician, not advice"),
        "phenoconversion: the docs line says it is an aid for the clinician"
    );

    // a line that starts with a separator still names its drug
    let r = run_report(
        &work,
        "PHENOSEP",
        &pheno_report(Some("2.0")),
        &meds(", paroxetine\n;\n"),
    );
    assert!(
        cell(&r.rows, "CYP2D6", 4) == "Poor Metabolizer (phenoconverted by paroxetine)"
            && !section(&r.rec).contains("Not in the inhibitor table"),
        "phenoconversion: ', paroxetine' is read as paroxetine"
    );
}

#[test]
fn a_moderate_inhibitor_halves_the_activity_score() {
    let work = Scratch::new("phenomod");
    let r = run_report(
        &work,
        "PHENOMOD",
        &pheno_report(Some("2.0")),
        &meds("duloxetine\n"),
    );
    assert!(
        cell(&r.rows, "CYP2D6", 4) == "Intermediate Metabolizer (phenoconverted by duloxetine)"
            && section(&r.rec).contains("activity score 2 x 0.5 = 1"),
        "phenoconversion: moderate inhibitor, AS 2 x 0.5 = 1 is an intermediate metabolizer"
    );
    assert_eq!(
        cell(&r.rows, "CYP2C19", 4),
        "-",
        "phenoconversion: a gene without an inhibitor on the list gets '-'"
    );
    let r = run_report(
        &work,
        "PHENONOAS",
        &pheno_report(None),
        &meds("duloxetine\n"),
    );
    assert_eq!(
        cell(&r.rows, "CYP2D6", 4),
        "not converted: the report gives no activity score",
        "phenoconversion: moderate inhibitor without an activity score converts nothing"
    );
}

#[test]
fn a_list_that_names_no_drug_changes_nothing() {
    let work = Scratch::new("phenonone");
    let r0 = plain(&work, "PHENONONE", &pheno_report(Some("2.0")));
    let r1 = run_report(
        &work,
        "PHENOEMPTY",
        &pheno_report(Some("2.0")),
        &meds("# nothing taken\n\n"),
    );
    let t0 = fs::read_to_string(work.0.join("PHENONONE/PHENONONE_phenotypes.tsv")).unwrap();
    let t1 = fs::read_to_string(work.0.join("PHENOEMPTY/PHENOEMPTY_phenotypes.tsv")).unwrap();
    assert_eq!(
        t0, t1,
        "phenoconversion: an empty list leaves the phenotypes table as without one"
    );
    assert!(
        section(&r1.rec).contains("The medication list names no drug: no phenotype is converted.")
            && !r1.rec.replace("PHENOEMPTY", "").contains("phenoconverted"),
        "phenoconversion: an empty list says so and converts nothing"
    );
    assert!(
        !r0.rec.contains("Medications and Phenoconversion"),
        "phenoconversion: no list, no section"
    );
}

/// Control: the inhibitor's factor set to 1 must not make a poor
/// metabolizer.
#[test]
fn control_a_factor_of_1_converts_nothing() {
    let work = Scratch::new("phenobroken");
    let broken = edited_inhibitors(
        &work,
        "inhibitors-factor-1.tsv",
        "paroxetine\tCYP2D6\tstrong\t0\t",
        "paroxetine\tCYP2D6\tstrong\t1\t",
    );
    let r = run_report(
        &work,
        "PHENOBROKEN",
        &pheno_report(Some("2.0")),
        &meds_with("paroxetine\n", broken),
    );
    assert!(
        !cell(&r.rows, "CYP2D6", 4).starts_with("Poor Metabolizer"),
        "control: with the factor set to 1 the NM is not converted to a poor metabolizer: {:?}",
        row(&r.rows, "CYP2D6")
    );
}

#[test]
fn a_nan_score_or_a_broken_inhibitor_row_converts_nothing() {
    let work = Scratch::new("phenonan");
    let r = run_report(
        &work,
        "PHENONAN",
        &pheno_report(Some("NaN")),
        &meds("duloxetine\n"),
    );
    assert_eq!(
        cell(&r.rows, "CYP2D6", 4),
        "not converted: the report gives no activity score",
        "phenoconversion: a NaN activity score converts nothing"
    );
    let unreadable = "Could not read the medication list or the inhibitor table";
    let short = work.file(
        "inhibitors-short-row.tsv",
        "Drug\tGene\tStrength\tActivity_score_factor\tClass_source\tRule_source\nparoxetine\n",
    );
    let r = run_report(
        &work,
        "PHENOSHORT",
        &pheno_report(Some("2.0")),
        &meds_with("paroxetine\n", short),
    );
    assert!(
        r.rc == 0 && section(&r.rec).contains(unreadable),
        "phenoconversion: a table row with no Gene cell is reported as unreadable, the report still written"
    );
    let blank = edited_inhibitors(
        &work,
        "inhibitors-blank-gene.tsv",
        "paroxetine\tCYP2D6\t",
        "paroxetine\t \t",
    );
    let r = run_report(
        &work,
        "PHENOBLANK",
        &pheno_report(Some("2.0")),
        &meds_with("paroxetine\n", blank),
    );
    assert!(
        r.rc == 0 && section(&r.rec).contains(unreadable),
        "phenoconversion: a row with a blank Gene cell is reported as unreadable, never read as a known drug"
    );
    let nofac = edited_inhibitors(
        &work,
        "inhibitors-no-factor.tsv",
        "paroxetine\tCYP2D6\tstrong\t0\t",
        "paroxetine\tCYP2D6\tstrong\t\t",
    );
    let r = run_report(
        &work,
        "PHENONOFAC",
        &pheno_report(Some("2.0")),
        &meds_with("paroxetine\n", nofac),
    );
    assert!(
        r.rc == 0 && section(&r.rec).contains(unreadable),
        "phenoconversion: a CYP2D6 row with an empty factor is reported as unreadable, the report still written"
    );
    let outdir = work.0.to_string_lossy().into_owned();
    assert_eq!(
        og(&[
            "cpic-report",
            "--sample",
            "X",
            "--report",
            "r.json",
            "--outdir",
            &outdir,
            "--medications",
            "m.txt"
        ]),
        2,
        "--medications without --inhibitors is refused"
    );
}

#[test]
fn cyp2c9_and_cyp3a5_inhibitors_give_a_note() {
    let work = Scratch::new("pheno2c9");
    let mut rep = pheno_report(Some("2.0"));
    if let Value::Object(m) = at(&mut rep, "genes") {
        m.insert("CYP2C9".to_owned(), dip("*1", "*1", "Normal Metabolizer"));
        m.insert("CYP3A5".to_owned(), dip("*3", "*3", "Poor Metabolizer"));
    }
    let r = run_report(
        &work,
        "PHENO2C9",
        &rep,
        &meds("fluconazole\nketoconazole\n"),
    );
    assert_eq!(
        cell(&r.rows, "CYP2C9", 4),
        "not computed: fluconazole inhibits CYP2C9; CPIC gives no rule that converts the CYP2C9 phenotype",
        "CYP2C9: a moderate inhibitor (fluconazole) changes the reported column to a note"
    );
    assert_eq!(
        row(&r.rows, "CYP2C9").map(|r| r[2..4].to_vec()),
        Some(strs(&["Normal Metabolizer", "normal"])),
        "CYP2C9: the genotype phenotype is kept"
    );
    let sec = section(&r.rec);
    assert!(
        sec.contains(
            "WARNING: CYP2C9 genotype phenotype Normal Metabolizer; on fluconazole (moderate CYP2C9 inhibitor)."
        ),
        "CYP2C9: a WARNING line names the inhibitor: {sec}"
    );
    assert!(
        cell(&r.rows, "CYP3A5", 4).starts_with("not computed: ketoconazole inhibits CYP3A5"),
        "CYP3A5: a strong CYP3A inhibitor (ketoconazole) changes the reported column to a note"
    );
    assert!(
        !sec.contains("can lower CYP3A5 activity") && sec.contains("can lower CYP2C9 activity"),
        "CYP3A5 *3/*3 (poor metabolizer): no line says the inhibitor can lower its activity; CYP2C9 keeps it"
    );
    let r = run_report(
        &work,
        "PHENOBOOSTED",
        &rep,
        &meds("lopinavir/ritonavir 400/100 mg\n"),
    );
    assert!(
        cell(&r.rows, "CYP3A5", 4)
            .starts_with("not computed: lopinavir, ritonavir inhibits CYP3A5")
            && !section(&r.rec).contains("Not in the inhibitor table"),
        "a boosted combination written 'lopinavir/ritonavir' is matched as both drugs"
    );
    let r = run_report(&work, "PHENO2C9NONE", &rep, &meds("paroxetine\n"));
    assert!(
        cell(&r.rows, "CYP2C9", 4) == "-" && cell(&r.rows, "CYP3A5", 4) == "-",
        "CYP2C9 and CYP3A5: without an inhibitor of theirs on the list the column stays '-'"
    );
}

// --- 9. disease-variant contraindications ------------------------------------------

const HITS_HEADER: &str =
    "chrom\tpos\tref\talt\tgenotype\tclinvar_id\tgeneinfo\tclnsig\tclnrevstat\n";

/// One synthetic row of step 06's hits table (no one's variant).
fn hit(gene: &str, clnsig: &str, gt: &str) -> String {
    format!(
        "{HITS_HEADER}chr15\t89300000\tC\tT\t{gt}\t9900001\t{gene}:1\t{clnsig}\tcriteria_provided,_single_submitter\n"
    )
}

fn contra_section(rec: &str) -> &str {
    between(
        rec,
        "Drug Contraindications From Disease Variants (ClinVar hits, step 06):",
        "Uncallable Genes:",
    )
}

fn contra_opts<'a>(table: String, hits: Option<String>, m: Option<&'a str>) -> Opts<'a> {
    Opts {
        contraindications: Some(table),
        hits,
        medications: m,
        ..Opts::default()
    }
}

#[test]
fn contraindication_table_has_the_shape_the_reader_reads() {
    let work = Scratch::new("contraread");
    let (header, rows) = table_rows(&contra());
    assert_eq!(header, CONTRA_COLUMNS, "contraindications: the header");
    let bad: Vec<&Vec<String>> = rows
        .iter()
        .filter(|r| {
            r.len() != 8
                || !r.iter().all(|c| !c.trim().is_empty())
                || !["A", "B", "A/B"].contains(&r[3].as_str())
                || r[5]
                    != format!(
                        "https://api.cpicpgx.org/v1/pair_view?genesymbol=eq.{}",
                        r[0]
                    )
                || !r[6].starts_with("https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=")
                || !["any", "biallelic"].contains(&r[7].as_str())
                || r[2]
                    .split(',')
                    .any(|n| n != n.trim().to_lowercase() || n.is_empty())
                || !r[2]
                    .split(',')
                    .any(|n| Some(n) == r[1].split_whitespace().next())
        })
        .collect();
    assert!(
        bad.is_empty(),
        "contraindications: every row has a gene, a drug, its names, a CPIC level, a note and its sources: {bad:?}"
    );
    let keys: Vec<(&str, &str)> = rows
        .iter()
        .map(|r| (r[0].as_str(), r[1].as_str()))
        .collect();
    let mut unique = keys.clone();
    unique.sort_unstable();
    unique.dedup();
    assert_eq!(
        unique.len(),
        keys.len(),
        "contraindications: no gene-drug pair twice"
    );
    for pair in [
        ("POLG", "valproic acid"),
        ("ASL", "valproic acid"),
        ("ASS1", "valproic acid"),
        ("CPS1", "valproic acid"),
        ("NAGS", "valproic acid"),
        ("OTC", "valproic acid"),
        ("HPRT1", "mycophenolic acid"),
    ] {
        assert!(
            keys.contains(&pair),
            "contraindications: POLG with valproic acid, the urea-cycle genes and HPRT1 are there: {pair:?}"
        );
    }
    assert!(
        !rows.iter().any(|r| r[0] == "SCN1A"),
        "contraindications: SCN1A is left out (no label or guideline names it)"
    );
    let r = run_report(
        &work,
        "CONTRAREAD",
        &pheno_report(Some("2.0")),
        &contra_opts(contra(), None, None),
    );
    let mut order: Vec<&str> = Vec::new();
    for r in &rows {
        if !order.contains(&r[0].as_str()) {
            order.push(&r[0]);
        }
    }
    let sec = contra_section(&r.rec);
    assert!(
        !sec.contains("Could not read")
            && sec.contains(&format!("Genes it would check: {}.", order.join(", "))),
        "contraindications: the parser reads the table: {sec}"
    );
}

#[test]
fn a_pathogenic_polg_hit_adds_the_valproate_warning() {
    let work = Scratch::new("contraplp");
    let rep = pheno_report(Some("2.0"));
    let r = run_report(
        &work,
        "CONTRAPLP",
        &rep,
        &contra_opts(contra(), Some(hit("POLG", "Pathogenic", "0/1")), None),
    );
    let sec = contra_section(&r.rec);
    assert_eq!(r.rc, 0, "contraindications: exits 0");
    assert!(
        sec.contains(
            "WARNING: valproic acid -- POLG chr15:89300000 C>T, genotype 0/1, ClinVar Pathogenic."
        ),
        "contraindications: a pathogenic POLG hit adds the valproic acid WARNING: {sec}"
    );
    assert!(
        sec.contains("mitochondrial disorders caused by POLG mutations")
            && sec.contains("(provisional CPIC level A/B: no CPIC guideline covers this pair)")
            && sec.contains("dailymed.nlm.nih.gov"),
        "contraindications: the line carries the label's note, the CPIC level and the label"
    );
    assert!(
        sec.contains("Review it with your clinician"),
        "contraindications: the section says to review it with a clinician"
    );
    let r = run_report(
        &work,
        "CONTRALP",
        &rep,
        &contra_opts(
            contra(),
            Some(hit("POLG", "Pathogenic/Likely_pathogenic", "1/1")),
            None,
        ),
    );
    assert!(
        contra_section(&r.rec).contains("WARNING: valproic acid -- POLG"),
        "contraindications: Pathogenic/Likely_pathogenic counts"
    );
    for (name, sig) in [
        ("CONTRABEN", "Benign"),
        ("CONTRACONF", "Conflicting_classifications_of_pathogenicity"),
        ("CONTRARISK", "risk_factor"),
    ] {
        let r = run_report(
            &work,
            name,
            &rep,
            &contra_opts(contra(), Some(hit("POLG", sig, "0/1")), None),
        );
        let sec = contra_section(&r.rec);
        assert!(
            !sec.contains("WARNING")
                && sec.contains("No pathogenic or likely pathogenic ClinVar hit in POLG"),
            "contraindications: a POLG hit classed {sig} adds no drug line: {sec}"
        );
    }
    let r = run_report(
        &work,
        "CONTRAOTHER",
        &rep,
        &contra_opts(contra(), Some(hit("BRCA2", "Pathogenic", "0/1")), None),
    );
    assert!(
        !contra_section(&r.rec).contains("WARNING"),
        "contraindications: a pathogenic hit in a gene outside the table adds no drug line"
    );
    let r = run_report(
        &work,
        "CONTRAMED",
        &rep,
        &contra_opts(
            contra(),
            Some(hit("POLG", "Likely_pathogenic", "0/1")),
            Some("Valproate 500 mg\n"),
        ),
    );
    assert!(
        contra_section(&r.rec).contains("ClinVar Likely_pathogenic; ON THE MEDICATION LIST."),
        "contraindications: valproate on the medication list is marked"
    );
    let r = run_report(
        &work,
        "CONTRANOTON",
        &rep,
        &contra_opts(
            contra(),
            Some(hit("POLG", "Pathogenic", "0/1")),
            Some("paroxetine\n"),
        ),
    );
    assert!(
        contra_section(&r.rec).contains("WARNING: valproic acid")
            && !r.rec.contains("ON THE MEDICATION LIST"),
        "contraindications: a drug not on the list is not marked"
    );
}

#[test]
fn no_or_an_unreadable_hits_table() {
    let work = Scratch::new("contranohits");
    let rep = pheno_report(Some("2.0"));
    let r = run_report(
        &work,
        "CONTRANOHITS",
        &rep,
        &contra_opts(contra(), None, None),
    );
    let sec = contra_section(&r.rec);
    assert!(
        sec.contains("Not checked: no ClinVar hits table") && sec.contains("POLG"),
        "contraindications: without a hits table the genes are listed as not checked"
    );
    let r = run_report(
        &work,
        "CONTRABADHITS",
        &rep,
        &contra_opts(contra(), Some("nothing\there\n".to_owned()), None),
    );
    assert!(
        r.rc == 0 && contra_section(&r.rec).contains("Could not read the ClinVar hits table"),
        "contraindications: an unreadable hits table is reported, the report still written"
    );
    let r = plain(&work, "CONTRANONE", &rep);
    assert!(
        !r.rec.contains("Drug Contraindications"),
        "contraindications: without the table there is no section"
    );
    let outdir = work.0.to_string_lossy().into_owned();
    assert_eq!(
        og(&[
            "cpic-report",
            "--sample",
            "X",
            "--report",
            "r.json",
            "--outdir",
            &outdir,
            "--clinvar-hits",
            "h.tsv"
        ]),
        2,
        "--clinvar-hits without --contraindications is refused"
    );
}

/// Control: the same hit with the POLG row taken out of the table gives no
/// valproate line.
#[test]
fn control_without_the_polg_row_no_valproate_line() {
    let work = Scratch::new("contractrl");
    let text = fs::read_to_string(contra()).unwrap();
    let kept: String = text
        .split_inclusive('\n')
        .filter(|l| !l.starts_with("POLG\t"))
        .collect();
    assert_ne!(kept, text, "the control took the POLG row out");
    let nopolg = work.file("contraindications-no-polg.tsv", &kept);
    let r = run_report(
        &work,
        "CONTRACTRL",
        &pheno_report(Some("2.0")),
        &contra_opts(nopolg, Some(hit("POLG", "Pathogenic", "0/1")), None),
    );
    assert!(
        !contra_section(&r.rec).contains("valproic acid"),
        "control: without the POLG row the same hit gives no valproic acid line"
    );
}
