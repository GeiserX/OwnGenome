//! What step 27 (`og pgx cpic-report --consensus`) makes of step 36's
//! consensus table (`og pgx outside-calls`):
//!
//! 8. CYP2D6 held back: the CPIC report says so, names the drugs it
//!    affects, and drops the pypgx-only warning it prints without the table;
//! 9. HLA passed on: the CPIC report lists it in the outside-call section
//!    and marks the gene PharmCAT reports with callSource OUTSIDE;
//! 10. a stale report: PharmCAT's report has outside calls the consensus
//!     does not pass on (another HLA-B, a CYP2D6 now held back, or no
//!     readable table at all); they give no drug guidance and the report
//!     says to rerun step 07.
//!
//! The consensus tables are the literal text `og pgx outside-calls` writes
//! for the named fixtures of `tests/fixtures/pgx/`; `og-steps`'s
//! `tests/pgx_outside_calls.rs` holds the verb to exactly these tables.

use std::fs;
use std::path::{Path, PathBuf};

const HEADER: &str = "Gene\tResult\tOutside_call\tReason\tEvidence\n";
const NO_HLA: &str = "HLA-A\tnot typed\tno\tHLA typing (step 08) did not run\t-\n\
                      HLA-B\tnot typed\tno\tHLA typing (step 08) did not run\t-\n";
const NO_CYP2D6: &str = "CYP2D6\tindeterminate\tno\tno caller made a call\
                         \tpypgx: not run; Cyrius: not run; depth check: not run\n";

/// `--pypgx pypgx_cyp2d6_1_4.tsv --depth-check depth_check_ok.tsv`.
fn held_back() -> String {
    format!(
        "{HEADER}{NO_HLA}CYP2D6\tindeterminate\tno\tone caller only: a second caller must agree (Cyrius: setup.sh --cyrius)\
         \tpypgx: *1/*4; Cyrius: not run; depth check: ok\n"
    )
}

/// `--hla t1k_genotype_docs_example.tsv`.
fn hla_docs_example() -> String {
    format!(
        "{HEADER}HLA-A\t*02:01/*24:02\tyes\tT1K (step 08)\tHLA-A*02:01:01 (quality 60); HLA-A*24:02:01 (quality 60)\n\
         HLA-B\t*15:02/*57:01\tyes\tT1K (step 08)\tHLA-B*15:02:01 (quality 60); HLA-B*57:01:01 (quality 48)\n\
         {NO_CYP2D6}"
    )
}

/// `--hla t1k_genotype.tsv`.
fn hla_other() -> String {
    format!(
        "{HEADER}HLA-A\t*02:01/*24:02\tyes\tT1K (step 08)\tHLA-A*02:01:01 (quality 60); HLA-A*24:02:01 (quality 60)\n\
         HLA-B\t*07:02/*44:02\tyes\tT1K (step 08)\tHLA-B*07:02:01 (quality 60); HLA-B*44:02:01 (quality 48)\n\
         {NO_CYP2D6}"
    )
}

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

/// A scratch directory of one test: emptied when made, removed when
/// dropped.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Self {
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("pgx-cons-{name}"));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

enum Table {
    None,
    Text(String),
    Missing,
}

/// The recommendations text of `og pgx cpic-report` on a PharmCAT fixture
/// report, with pypgx's fixture summary and a consensus table when given.
fn cpic(tag: &str, report: &str, pypgx: Option<&str>, table: Table) -> String {
    let tmp = Scratch::new(tag);
    let dir = tmp.0.to_string_lossy().into_owned();
    let mut args: Vec<String> = vec![
        "cpic-report".into(),
        "--sample".into(),
        "T".into(),
        "--outdir".into(),
        dir,
        "--report".into(),
        repo()
            .join("tests/fixtures/pharmcat")
            .join(report)
            .to_string_lossy()
            .into_owned(),
    ];
    if let Some(p) = pypgx {
        args.push("--pypgx".into());
        args.push(
            repo()
                .join("tests/fixtures/pgx")
                .join(p)
                .to_string_lossy()
                .into_owned(),
        );
    }
    let cons = tmp.0.join(match table {
        Table::Missing => "missing.tsv",
        _ => "consensus.tsv",
    });
    match &table {
        Table::None => {}
        Table::Text(t) => {
            fs::write(&cons, t).unwrap();
            args.push("--consensus".into());
            args.push(cons.to_string_lossy().into_owned());
        }
        Table::Missing => {
            args.push("--consensus".into());
            args.push(cons.to_string_lossy().into_owned());
        }
    }
    let rc = og_pgx::main(&args);
    if rc != 0 {
        return format!("exit {rc}");
    }
    fs::read_to_string(tmp.0.join("T_cpic_recommendations.txt")).unwrap()
}

/// The Affected Medications section of a recommendations text.
fn medications(text: &str) -> &str {
    let rest = text
        .split_once("Affected Medications:")
        .map_or(text, |(_, r)| r);
    rest.split("Calls From Other Tools")
        .next()
        .unwrap_or_default()
}

const DOCS: &str = "pharmcat-docs-example.json";

/// Case 8: the 3.4.0 report of the fixture: PharmCAT has no CYP2D6 result.
#[test]
fn case8_cyp2d6_held_back() {
    let plain = cpic(
        "8plain",
        "report-3.4.0.json",
        Some("pypgx_cyp2d6_1_4.tsv"),
        Table::None,
    );
    assert!(
        plain.contains("WARNING: PharmCAT has no result for CYP2D6, but pypgx"),
        "8 control: without the table, pypgx's call alone draws a warning"
    );
    let held = cpic(
        "8held",
        "report-3.4.0.json",
        Some("pypgx_cyp2d6_1_4.tsv"),
        Table::Text(held_back()),
    );
    assert!(
        held.contains("Calls From Other Tools (outside calls, step 36):"),
        "8 the report has the outside-call section"
    );
    assert!(
        held.contains("CYP2D6   indeterminate            not passed to PharmCAT: one caller only"),
        "8 CYP2D6 is indeterminate and why"
    );
    assert!(
        held.contains("No drug guidance is given for CYP2D6 here. Drugs affected by CYP2D6: "),
        "8 it names the drugs CYP2D6 affects"
    );
    assert!(
        !held.contains("PharmCAT has no result for CYP2D6"),
        "8 no pypgx-only warning for the held-back gene"
    );
}

/// Case 9: PharmCAT's example report has HLA-B from an outside call (and CYP2D6
/// and MT-RNR1, which this consensus does not pass on: case 10).
#[test]
fn case9_hla_passed_on() {
    let passed = cpic("9", DOCS, None, Table::Text(hla_docs_example()));
    assert!(
        passed.contains(
            "HLA-B    *15:02/*57:01            passed to PharmCAT from T1K (step 08); PharmCAT reports it as an outside call."
        ),
        "9 HLA-B is listed as passed from T1K"
    );
    assert!(
        passed.contains("[outside call]"),
        "9 the gene line marks PharmCAT's outside call"
    );
    assert!(
        passed.contains("  HLA-B -- "),
        "9 the confirmed HLA-B keeps its drug guidance"
    );
    assert!(
        !passed.contains("outside call HLA-B"),
        "9 no stale warning for the confirmed HLA-B"
    );
}

/// Case 10: outside calls the consensus does not pass on.
#[test]
fn case10_a_stale_report() {
    let plain = cpic("10plain", DOCS, None, Table::None);
    assert!(
        medications(&plain).contains("  CYP2D6 -- "),
        "10 control: without the table, the example's CYP2D6 has drug guidance"
    );
    assert!(
        medications(&plain).contains("  HLA-B -- "),
        "10 control: without the table, the example's HLA-B has drug guidance"
    );
    let other = cpic("10other", DOCS, None, Table::Text(hla_other()));
    assert!(
        !medications(&other).contains("  HLA-B -- "),
        "10 HLA-B of another call gets no drug guidance"
    );
    assert!(
        other.contains(
            "WARNING: PharmCAT's report has the outside call HLA-B *15:02/*57:01, which the consensus does not pass on"
        ),
        "10 the report says the HLA-B outside call is stale"
    );
    assert!(
        other.contains(
            "HLA-B    *07:02/*44:02            passed to PharmCAT from T1K (step 08); PharmCAT's report does not list it"
        ),
        "10 the consensus HLA-B is not shown as PharmCAT's"
    );
    let held = cpic("10held", DOCS, None, Table::Text(held_back()));
    assert!(
        !medications(&held).contains("  CYP2D6 -- "),
        "10 a held-back CYP2D6 gets no drug guidance from an older report"
    );
    assert!(
        held.contains(
            "WARNING: PharmCAT's report has the outside call CYP2D6 *1/*3, which the consensus does not pass on"
        ),
        "10 the report says the CYP2D6 outside call is stale"
    );
    assert!(
        held.contains("Rerun step 07, then this step. No drug guidance is given for CYP2D6 here."),
        "10 and says how to fix it"
    );
    let gone = cpic("10gone", DOCS, None, Table::Missing);
    assert!(
        !medications(&gone).contains("  HLA-B -- "),
        "10 with the table missing, no outside call gets drug guidance"
    );
    assert!(
        gone.contains("Run step 36 again, then steps 07 and 27."),
        "10 with the table missing, the report says to run step 36 again"
    );
}
