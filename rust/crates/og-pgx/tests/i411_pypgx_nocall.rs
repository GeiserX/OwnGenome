//! `og pgx cpic-report` with pypgx's `Indeterminate` and with combination
//! products on the medication list (issue #411, pgx-06, pgx-30, pgx-16):
//!
//! - pypgx writes `Indeterminate` for a gene it could not call, and step 32
//!   writes it for CYP2D6 when the depth check fails. In any case, and as
//!   one allele (`*1/Indeterminate`), it is no call: a gene PharmCAT could
//!   not call gets no comparison row and no warning, and a gene PharmCAT
//!   did call is `PharmCAT only`, never a mismatch;
//! - every drug of a combination product is checked: `olanzapine/fluoxetine`
//!   converts CYP2D6 by fluoxetine, `naltrexone/bupropion 8/90 mg` by
//!   bupropion and is not listed as unchecked, and the second drug of a
//!   combination is marked as on the medication list in the drug
//!   contraindications.
//!
//! Each case has its control: a real pypgx call still warns, and a list
//! without the inhibitor converts nothing.

use std::fs;
use std::path::PathBuf;

const ASSETS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../../assets/");

struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Scratch {
        let d =
            std::env::temp_dir().join(format!("og-pgx-i411-nocall-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        Scratch(d)
    }

    fn put(&self, name: &str, body: &str) -> String {
        let p = self.0.join(name);
        fs::write(&p, body).unwrap();
        p.to_string_lossy().into_owned()
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

/// One gene of a PharmCAT 3.x report with a single diplotype.
fn gene(name: &str, a1: &str, a2: &str, phenotype: &str, score: Option<&str>) -> String {
    let score = score.map_or_else(String::new, |s| format!(r#", "activityScore": "{s}""#));
    format!(
        r#""{name}": {{"sourceDiplotypes": [{{"allele1": {{"name": "{a1}"}}, "allele2": {{"name": "{a2}"}}, "label": "{a1}/{a2}", "phenotypes": ["{phenotype}"]{score}}}], "relatedDrugs": [{{"name": "codeine"}}]}}"#
    )
}

fn report(genes: &[String]) -> String {
    format!(
        r#"{{"pharmcatVersion": "3.4.0", "genes": {{{}}}}}"#,
        genes.join(", ")
    )
}

/// `og pgx cpic-report` on `report` with the extra options: (exit code,
/// recommendations text, phenotype rows, comparison rows).
fn run(
    s: &Scratch,
    report: &str,
    extra: &[&str],
) -> (u8, String, Vec<Vec<String>>, Vec<Vec<String>>) {
    let rep = s.put("report.json", report);
    let out = s.0.join("out").to_string_lossy().into_owned();
    let mut args: Vec<String> = [
        "cpic-report",
        "--sample",
        "S",
        "--report",
        &rep,
        "--outdir",
        &out,
        "--date",
        "2026-01-01",
    ]
    .map(str::to_owned)
    .to_vec();
    args.extend(extra.iter().map(|a| (*a).to_owned()));
    let rc = og_pgx::main(&args);
    let rows = |text: String| -> Vec<Vec<String>> {
        text.lines()
            .skip(1)
            .map(|l| l.split('\t').map(str::to_owned).collect())
            .collect()
    };
    (
        rc,
        s.read("out/S_cpic_recommendations.txt"),
        rows(s.read("out/S_phenotypes.tsv")),
        rows(s.read("comparison.tsv")),
    )
}

fn row<'a>(rows: &'a [Vec<String>], gene: &str) -> Option<&'a Vec<String>> {
    rows.iter().find(|r| r.first().is_some_and(|g| g == gene))
}

/// CYP2D6 not called by PharmCAT, CYP2C9 called; pypgx's CYP2D6 and CYP2C9
/// diplotypes as given.
fn pypgx_run(name: &str, cyp2d6: &str, cyp2c9: &str) -> (String, Vec<Vec<String>>) {
    let s = Scratch::new(name);
    let rep = report(&[
        gene("CYP2D6", "Unknown", "Unknown", "No Result", None),
        gene("CYP2C9", "*1", "*1", "Normal Metabolizer", None),
    ]);
    let pg = s.put(
        "pypgx.tsv",
        &format!(
            "Gene\tDiplotype\tPhenotype\tCNV_call\tSource\nCYP2D6\t{cyp2d6}\tPoor Metabolizer\t.\tbam\nCYP2C9\t{cyp2c9}\tIndeterminate\t.\tvcf\n"
        ),
    );
    let comp = s.0.join("comparison.tsv").to_string_lossy().into_owned();
    let (rc, rec, _, comp) = run(&s, &rep, &["--pypgx", &pg, "--comparison", &comp]);
    assert_eq!(rc, 0, "{rec}");
    (rec, comp)
}

#[test]
fn indeterminate_pypgx_row_is_no_call() {
    for indeterminate in [
        "Indeterminate",
        "indeterminate",
        "*1/Indeterminate",
        "Indeterminate/Indeterminate",
    ] {
        let (rec, comp) = pypgx_run("indet", indeterminate, indeterminate);
        assert!(
            row(&comp, "CYP2D6").is_none(),
            "{indeterminate}: no caller called CYP2D6, so no comparison row: {comp:?}"
        );
        let c9 = row(&comp, "CYP2C9").expect("PharmCAT called CYP2C9");
        assert_eq!(
            (c9[3].as_str(), c9[4].as_str()),
            ("PharmCAT only", "PharmCAT only"),
            "{indeterminate}: never a mismatch"
        );
        assert!(
            !rec.contains("pypgx (step 32) called") && !rec.contains("Drugs affected by CYP2D6"),
            "{indeterminate}: no warning reads it as a pypgx call: {rec}"
        );
        assert!(
            rec.contains("No gene that PharmCAT could not call has a pypgx call."),
            "{indeterminate}: {rec}"
        );
    }
}

#[test]
fn a_real_pypgx_call_still_warns() {
    let (rec, comp) = pypgx_run("control", "*4/*4", "*1/*2");
    assert!(
        rec.contains(
            "WARNING: PharmCAT has no result for CYP2D6, but pypgx (step 32) called *4/*4 (Poor Metabolizer)."
        ),
        "{rec}"
    );
    assert_eq!(
        row(&comp, "CYP2D6").map(|r| r[3].as_str()),
        Some("pypgx only")
    );
    assert_eq!(row(&comp, "CYP2C9").map(|r| r[3].as_str()), Some("No"));
}

/// The phenoconversion column of CYP2D6 and the recommendations text for a
/// CYP2D6 call and a medication list.
fn convert(name: &str, cyp2d6: &[&str; 4], meds: &str) -> (String, String) {
    let s = Scratch::new(name);
    let rep = report(&[gene(
        "CYP2D6",
        cyp2d6[0],
        cyp2d6[1],
        cyp2d6[2],
        Some(cyp2d6[3]),
    )]);
    let list = s.put("meds.txt", meds);
    let inh = format!("{ASSETS}pgx_inhibitors.tsv");
    let (rc, rec, rows, _) = run(&s, &rep, &["--medications", &list, "--inhibitors", &inh]);
    assert_eq!(rc, 0, "{rec}");
    let col = row(&rows, "CYP2D6")
        .and_then(|r| r.get(4).cloned())
        .unwrap_or_default();
    (col, rec)
}

const POOR: [&str; 4] = ["*4", "*4", "Poor Metabolizer", "0.0"];
const NORMAL: [&str; 4] = ["*1", "*1", "Normal Metabolizer", "2.0"];

#[test]
fn olanzapine_fluoxetine_is_checked_for_fluoxetine() {
    let (col, rec) = convert("olafluo", &POOR, "olanzapine/fluoxetine 6/25 mg\n");
    assert_eq!(
        col, "Poor Metabolizer (phenoconverted by fluoxetine)",
        "{rec}"
    );
    assert!(!rec.contains("Not in the inhibitor table"), "{rec}");
    // Control: the first drug alone inhibits nothing in the table.
    let (col, rec) = convert("ola", &POOR, "olanzapine 5 mg\n");
    assert_eq!(col, "-", "{rec}");
    assert!(
        rec.contains("Not in the inhibitor table: olanzapine 5 mg."),
        "{rec}"
    );
}

#[test]
fn naltrexone_bupropion_converts_by_bupropion() {
    let (col, rec) = convert("nalbup", &NORMAL, "naltrexone/bupropion 8/90 mg\n");
    assert_eq!(
        col, "Poor Metabolizer (phenoconverted by bupropion)",
        "{rec}"
    );
    assert!(
        rec.contains(
            "WARNING: CYP2D6 genotype phenotype Normal Metabolizer; phenoconverted phenotype Poor Metabolizer (phenoconverted by bupropion), activity score adjusted to 0."
        ),
        "{rec}"
    );
    assert!(!rec.contains("Not in the inhibitor table"), "{rec}");
    // `and`, `+`, `;` and brackets split a line too, and a trailing stop
    // is not part of the drug's name.
    for line in [
        "naltrexone and bupropion",
        "naltrexone + bupropion",
        "naltrexone; bupropion",
        "Contrave (naltrexone/bupropion) 8/90 mg",
        "naltrexone, bupropion.",
    ] {
        let (col, rec) = convert("nalbup2", &NORMAL, &format!("{line}\n"));
        assert_eq!(
            col, "Poor Metabolizer (phenoconverted by bupropion)",
            "{line}: {rec}"
        );
    }
}

/// The contraindication section for one pathogenic POLG hit and a list.
fn contra(name: &str, meds: &str) -> String {
    let s = Scratch::new(name);
    let rep = report(&[gene(
        "CYP2D6",
        "*1",
        "*1",
        "Normal Metabolizer",
        Some("2.0"),
    )]);
    let list = s.put("meds.txt", meds);
    let hits = s.put(
        "hits.tsv",
        "chrom\tpos\tref\talt\tgenotype\tclinvar_id\tgeneinfo\tclnsig\tclnrevstat\nchr15\t89300000\tC\tT\t0/1\t1\tPOLG:1\tLikely_pathogenic\tx\n",
    );
    let (inh, table) = (
        format!("{ASSETS}pgx_inhibitors.tsv"),
        format!("{ASSETS}pgx_contraindications.tsv"),
    );
    let (rc, rec, _, _) = run(
        &s,
        &rep,
        &[
            "--medications",
            &list,
            "--inhibitors",
            &inh,
            "--contraindications",
            &table,
            "--clinvar-hits",
            &hits,
        ],
    );
    assert_eq!(rc, 0, "{rec}");
    rec
}

#[test]
fn a_combinations_second_drug_is_on_the_list() {
    let want = "  WARNING: valproic acid -- POLG chr15:89300000 C>T, genotype 0/1, ClinVar Likely_pathogenic; ON THE MEDICATION LIST.";
    let rec = contra("contra", "lamotrigine/valproate 100/500 mg\n");
    assert!(rec.lines().any(|l| l == want), "{rec}");
    // A brand name with its generic in brackets.
    for line in [
        "Epilim (sodium valproate) 500 mg",
        "Depakote (divalproex sodium)",
    ] {
        let rec = contra("contra-brand", &format!("{line}\n"));
        assert!(rec.lines().any(|l| l == want), "{line}: {rec}");
    }
    // Control: without valproate the same hit warns, not as on the list.
    let rec = contra("contra-control", "lamotrigine 100 mg\n");
    assert!(
        rec.contains("WARNING: valproic acid -- POLG") && !rec.contains("ON THE MEDICATION LIST"),
        "{rec}"
    );
}
