//! `og pgx cpic-report` on a gene whose call rests on a position PharmCAT
//! never read (issue #411, plat-49).
//!
//! PharmCAT calls a gene from the positions it read. An allele defined at a
//! position with no call cannot be called, so on an exome, which misses
//! CYP2C19 *2's site when the capture does, the gene comes out `*1/*1`,
//! Normal Metabolizer. PharmCAT's report lists every position of the gene in
//! `variants`, a position it never read with no `call`, and with `alleles`
//! the alleles defined there. Such a gene is now `incomplete`, never
//! `normal`: the phenotypes table says so, the text lists the unread
//! positions and gives no drug guidance for it. A non-normal call keeps its
//! status and drug list, marked the same way.
//!
//! The control: the same report with every position read (a whole genome)
//! writes what it wrote before, byte for byte the report without `variants`.

use std::fs;
use std::path::PathBuf;

const FIXTURES: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../../tests/fixtures/pharmcat/"
);

struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Scratch {
        let d = std::env::temp_dir().join(format!(
            "og-pgx-i411-incomplete-{name}-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        Scratch(d)
    }

    fn read(&self, name: &str) -> String {
        fs::read_to_string(self.0.join(name)).unwrap_or_default()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// `og pgx cpic-report` on the report text: (exit code, recommendations
/// text, phenotypes table).
fn run(name: &str, report: &str) -> (u8, String, String) {
    let s = Scratch::new(name);
    let rep = s.0.join("report.json");
    fs::write(&rep, report).unwrap();
    let out = s.0.join("out");
    let args: Vec<String> = [
        "cpic-report",
        "--sample",
        "S",
        "--report",
        &rep.to_string_lossy(),
        "--outdir",
        &out.to_string_lossy(),
        "--date",
        "2026-01-01",
    ]
    .map(str::to_owned)
    .to_vec();
    let rc = og_pgx::main(&args);
    (
        rc,
        s.read("out/S_cpic_recommendations.txt"),
        s.read("out/S_phenotypes.tsv"),
    )
}

fn fixture(name: &str) -> String {
    fs::read_to_string(format!("{FIXTURES}{name}")).unwrap()
}

fn row<'a>(table: &'a str, gene: &str) -> Vec<&'a str> {
    table
        .lines()
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .find(|r| r[0] == gene)
        .unwrap_or_else(|| panic!("no {gene} row in {table}"))
}

/// One section of the recommendations: its head line and every line up to
/// the next head (a line that is not blank, not indented and not a rule).
fn section(text: &str, head: &str) -> Option<String> {
    let mut lines = text.lines().skip_while(|l| *l != head);
    let first = lines.next()?;
    let body = lines.take_while(|l| {
        l.is_empty() || l.starts_with(' ') || l.starts_with('-') || l.starts_with('=')
    });
    Some(
        std::iter::once(first)
            .chain(body)
            .collect::<Vec<_>>()
            .join("\n"),
    )
}

const UNREAD: &str = "Genes Called With Positions Not Read:";

#[test]
fn an_exome_missing_the_cyp2c19_star2_site_reads_incomplete_not_normal() {
    let (rc, text, table) = run("exome", &fixture("exome-cyp2c19-unread.json"));
    assert_eq!(rc, 0);
    assert_eq!(
        row(&table, "CYP2C19"),
        [
            "CYP2C19",
            "*1/*1",
            "Normal Metabolizer [incomplete: positions not read: rs4244285]",
            "incomplete"
        ],
        "the gene whose *2 site was never read is incomplete, not normal"
    );
    let sec = section(&text, UNREAD).expect("a section lists the unread positions");
    assert!(
        sec.contains(
            "  CYP2C19 -- *1/*1, Normal Metabolizer: incomplete, 1 defining position not read (rs4244285)."
        ),
        "{sec}"
    );
    assert!(
        sec.contains("no drug guidance"),
        "the section says an incomplete gene gets no drug guidance: {sec}"
    );
    let meds = section(&text, "Affected Medications:").unwrap();
    assert!(!meds.contains("CYP2C19"), "{meds}");
    assert!(
        text.contains("CYP2C19      *1/*1                          Normal Metabolizer [incomplete: positions not read: rs4244285]\n"),
        "the gene results line carries the mark: {text}"
    );
}

#[test]
fn a_non_normal_call_with_an_unread_position_keeps_its_drugs_and_the_mark() {
    let (_, text, table) = run("exome-nn", &fixture("exome-cyp2c19-unread.json"));
    assert_eq!(
        row(&table, "CYP2C9"),
        [
            "CYP2C9",
            "*1/*2",
            "Intermediate Metabolizer [incomplete: positions not read: rs1057910]",
            "non-normal"
        ]
    );
    let meds = section(&text, "Affected Medications:").unwrap();
    assert!(
        meds.contains(
            "  CYP2C9 -- Intermediate Metabolizer [incomplete: positions not read: rs1057910]:"
        ),
        "{meds}"
    );
    let sec = section(&text, UNREAD).unwrap();
    assert!(
        sec.contains(
            "  CYP2C9 -- *1/*2, Intermediate Metabolizer: 1 defining position not read (rs1057910)."
        ),
        "{sec}"
    );
    assert_eq!(
        row(&table, "CYP3A5"),
        ["CYP3A5", "*3/*3", "Poor Metabolizer", "non-normal"],
        "a gene with every position read is unmarked"
    );
    assert!(!sec.contains("CYP3A5"), "{sec}");
}

/// The report with every `"variants"` list removed: what the parser read
/// before it looked at positions.
fn without_variants(report: &str) -> String {
    report
        .lines()
        .filter(|l| !l.trim_start().starts_with("{\"gene\": "))
        .collect::<Vec<_>>()
        .join("\n")
        .replace(",\n   \"variants\": [\n   ]", "")
}

#[test]
fn a_whole_genome_with_every_position_read_is_unchanged() {
    let full = fixture("wgs-cyp2c19-complete.json");
    let stripped = without_variants(&full);
    assert!(
        !stripped.contains("variants") && stripped.contains("CYP3A5"),
        "the control removes the positions and keeps the genes: {stripped}"
    );
    let (rc, text, table) = run("wgs", &full);
    assert_eq!(rc, 0);
    assert_eq!(
        row(&table, "CYP2C19"),
        ["CYP2C19", "*1/*1", "Normal Metabolizer", "normal"]
    );
    assert_eq!(
        row(&table, "CYP2C9"),
        ["CYP2C9", "*1/*2", "Intermediate Metabolizer", "non-normal"]
    );
    assert!(
        !text.contains("incomplete") && !text.contains(UNREAD),
        "{text}"
    );
    let (_, text0, table0) = run("wgs-control", &stripped);
    assert_eq!(
        (text, table),
        (text0, table0),
        "every position read: the same report as without positions"
    );
}

/// A one-gene PharmCAT 3.x report with these diplotypes and positions.
fn one_gene(source: &str, dips: &[(&str, &str)], variants: &[&str]) -> String {
    let dips: Vec<String> = dips
        .iter()
        .map(|(label, phen)| {
            let (a1, a2) = label.split_once('/').unwrap();
            format!(
                r#"{{"allele1": {{"name": "{a1}"}}, "allele2": {{"name": "{a2}"}}, "label": "{label}", "phenotypes": ["{phen}"]}}"#
            )
        })
        .collect();
    format!(
        r#"{{"pharmcatVersion": "3.4.0", "genes": {{"CYP2C19": {{"callSource": "{source}", "sourceDiplotypes": [{}], "variants": [{}]}}}}}}"#,
        dips.join(", "),
        variants.join(", ")
    )
}

fn status(report: &str, name: &str) -> (String, String) {
    let (rc, _, table) = run(name, report);
    assert_eq!(rc, 0, "{report}");
    let r = row(&table, "CYP2C19");
    (r[2].to_owned(), r[3].to_owned())
}

const NM: (&str, &str) = ("*1/*1", "Normal Metabolizer");
const READ: &str = r#"{"position": 1, "dbSnpId": "rs12248560", "call": "C/C", "alleles": ["*17"]}"#;

#[test]
fn only_an_unread_position_that_defines_an_allele_counts() {
    let unread = |call: &str| {
        format!(r#"{{"position": 94781859, "dbSnpId": "rs4244285"{call}, "alleles": ["*2"]}}"#)
    };
    for (call, why) in [
        ("", "no call field (PharmCAT leaves a null out)"),
        (r#", "call": null"#, "a null call"),
        (r#", "call": "  ""#, "a blank call"),
    ] {
        assert_eq!(
            status(&one_gene("MATCHER", &[NM], &[READ, &unread(call)]), "blank").1,
            "incomplete",
            "{why}"
        );
    }
    let no_alleles =
        r#"{"position": 94781859, "dbSnpId": "rs4244285", "call": null, "alleles": []}"#;
    assert_eq!(
        status(&one_gene("MATCHER", &[NM], &[READ, no_alleles]), "noall"),
        ("Normal Metabolizer".to_owned(), "normal".to_owned()),
        "a position that defines no allele cannot hide one"
    );
    let no_rsid = r#"{"chromosome": "chr10", "position": 94781859, "alleles": ["*2"]}"#;
    assert_eq!(
        status(&one_gene("MATCHER", &[NM], &[no_rsid]), "norsid").0,
        "Normal Metabolizer [incomplete: positions not read: chr10:94781859]",
        "a position without an rsID is named by its place"
    );
    let many: Vec<String> = (1..=5)
        .map(|i| format!(r#"{{"position": {i}, "dbSnpId": "rs{i}", "alleles": ["*{i}"]}}"#))
        .collect();
    let many: Vec<&str> = many.iter().map(String::as_str).collect();
    let (rc, text, table) = run("many", &one_gene("MATCHER", &[NM], &many));
    assert_eq!(rc, 0);
    assert_eq!(
        row(&table, "CYP2C19")[2],
        "Normal Metabolizer [incomplete: positions not read: rs1, rs2, rs3 and 2 more]"
    );
    assert!(
        text.contains("5 defining positions not read (rs1, rs2, rs3, rs4, rs5)."),
        "the section names every unread position: {text}"
    );
}

#[test]
fn not_called_ambiguous_and_outside_genes_are_left_as_they_were() {
    let unread = r#"{"position": 94781859, "dbSnpId": "rs4244285", "alleles": ["*2"]}"#;
    assert_eq!(
        status(&one_gene("OUTSIDE", &[NM], &[unread]), "outside"),
        ("Normal Metabolizer".to_owned(), "normal".to_owned()),
        "an outside call reads no position (PharmCAT marks none missing)"
    );
    assert_eq!(
        status(
            &one_gene("MATCHER", &[("Unknown/Unknown", "No Result")], &[unread]),
            "nocall"
        )
        .1,
        "not called"
    );
    assert_eq!(
        status(
            &one_gene(
                "MATCHER",
                &[NM, ("*1/*2", "Intermediate Metabolizer")],
                &[unread]
            ),
            "ambig"
        )
        .1,
        "ambiguous",
        "several possible results stay ambiguous"
    );
}
