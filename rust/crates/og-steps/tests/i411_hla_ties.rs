//! Step 36 and step 08's concordance table on T1K's tied candidates and
//! expression suffixes, and the CYP2D6 rule on pypgx's evidence columns.
//!
//! T1K writes equally likely candidates of one allele comma-separated
//! (`HLA-B*57:01:01,HLA-B*58:01:01`), with one quality for the rank. The
//! Python this replaced cut the whole field to two fields, so PharmCAT got
//! whichever candidate T1K printed first. These cases pin:
//!
//! 1. a B*57:01/B*58:01 tie gives the same `ambiguous` answer in both print
//!    orders, in either allele, with or without an OptiType that matches one
//!    candidate, and never reaches the outside-call file;
//! 2. a tie whose candidates agree at two fields passes as before;
//! 3. a tie that names no allele PharmCAT acts on still passes the gene
//!    (a sure B*57:01 on the other copy is not lost), order-independently;
//! 4. a null allele keeps its suffix (`*02:53N`), does not agree with
//!    OptiType's `*02:53` in step 08's table, and reaches PharmCAT without
//!    it, with or without OptiType (OptiType never calls a null); a null of
//!    a drug allele is ambiguous, a null on the other copy of a sure
//!    B*57:01 withholds nothing, with or without OptiType;
//! 5. T1K's single allele keeps the caveat "one allele, read as homozygous",
//!    and two alleles that OptiType agrees with keep their reason;
//! 6. CYP2D6: when pypgx and Cyrius agree and the depth check passed, the
//!    call is passed whatever pypgx's Defining_filter and Quality columns
//!    say (the owner's rule: the reports show that flag, step 36 does not
//!    withhold on it).

use std::fs;
use std::path::{Path, PathBuf};

use og_steps::pgx_rules::{
    Decision, HLA_GENES, OPTITYPE_GENES, concordance_rows, hla_row, outside_calls_cli,
    read_optitype, read_t1k, two_fields,
};

struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Self {
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("i411-hla-ties-{name}"));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }

    fn put(&self, name: &str, text: &str) -> String {
        let p = self.0.join(name);
        fs::write(&p, text).unwrap();
        p.to_string_lossy().into_owned()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// A T1K genotype row: gene, then up to two alleles (`.` for none).
fn t1k_row(gene: &str, a1: &str, q1: u32, a2: &str, q2: u32) -> String {
    let n = if a2 == "." { 1 } else { 2 };
    format!("{gene}\t{n}\t{a1}\t10\t{q1}\t{a2}\t9\t{q2}\n")
}

/// OptiType's result with these HLA-A and HLA-B pairs.
fn optitype(a: (&str, &str), b: (&str, &str)) -> String {
    format!(
        "\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n0\t{}\t{}\t{}\t{}\tC*07:02\tC*07:02\t100\t90.0\n",
        a.0, a.1, b.0, b.1
    )
}

/// `og pgx outside-calls` on these files: (outside calls, consensus table).
fn outside(dir: &Scratch, args: &[(&str, &str)], require: bool) -> (String, String) {
    let calls = dir.0.join("calls.tsv").to_string_lossy().into_owned();
    let cons = dir.0.join("cons.tsv").to_string_lossy().into_owned();
    let mut argv = vec![
        "--calls".to_owned(),
        calls.clone(),
        "--consensus".to_owned(),
        cons.clone(),
    ];
    for (k, v) in args {
        argv.push(format!("--{k}"));
        argv.push((*v).to_owned());
    }
    if require {
        argv.push("--require-optitype".to_owned());
    }
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = outside_calls_cli(&argv, &mut out, &mut err);
    assert_eq!(code, 0, "{}", String::from_utf8_lossy(&err));
    (
        fs::read_to_string(calls).unwrap(),
        fs::read_to_string(cons).unwrap(),
    )
}

/// The consensus row of `gene`, without its Evidence cell (which prints
/// T1K's field as written, so it differs with the print order).
fn decided(cons: &str, gene: &str) -> String {
    let line = cons
        .lines()
        .find(|l| l.starts_with(&format!("{gene}\t")))
        .unwrap_or_else(|| panic!("no {gene} row in {cons}"));
    line.rsplit_once('\t').unwrap().0.to_owned()
}

const A_ROW: &str = "HLA-A\t2\tHLA-A*01:01:01\t10\t60\tHLA-A*26:01:01\t9\t60\n";
const AMBIGUOUS_B: &str = "HLA-B\tindeterminate\tno\tambiguous: B*57:01 possible, B*58:01 possible \
                           (T1K ties *57:01 and *58:01 on one allele)";

#[test]
fn a_drug_allele_tie_gives_the_same_ambiguous_answer_in_both_print_orders() {
    let d = Scratch::new("tie");
    let ot = d.put(
        "ot.tsv",
        &optitype(("A*01:01", "A*26:01"), ("B*07:02", "B*57:01")),
    );
    for (tag, b) in [
        (
            "b-second",
            t1k_row(
                "HLA-B",
                "HLA-B*07:02:01",
                60,
                "HLA-B*57:01:01,HLA-B*58:01:01",
                48,
            ),
        ),
        (
            "b-second-rev",
            t1k_row(
                "HLA-B",
                "HLA-B*07:02:01",
                60,
                "HLA-B*58:01:01,HLA-B*57:01:01",
                48,
            ),
        ),
        (
            "b-first",
            t1k_row(
                "HLA-B",
                "HLA-B*57:01:01,HLA-B*58:01:01",
                48,
                "HLA-B*07:02:01",
                60,
            ),
        ),
        (
            "b-first-rev",
            t1k_row(
                "HLA-B",
                "HLA-B*58:01:01, HLA-B*57:01:01",
                48,
                "HLA-B*07:02:01",
                60,
            ),
        ),
    ] {
        let hla = d.put(&format!("{tag}.tsv"), &format!("{A_ROW}{b}"));
        // T1K alone, T1K with an OptiType that picks one candidate, and
        // the second caller required: the tie decides in none of them.
        for (opts, require) in [
            (vec![("hla", hla.as_str())], false),
            (
                vec![("hla", hla.as_str()), ("optitype", ot.as_str())],
                false,
            ),
            (vec![("hla", hla.as_str()), ("optitype", ot.as_str())], true),
        ] {
            let (calls, cons) = outside(&d, &opts, require);
            assert!(
                !calls.contains("HLA-B"),
                "{tag}: a tied HLA-B reached PharmCAT: {calls}"
            );
            assert!(calls.contains("HLA-A\t*01:01/*26:01\n"), "{tag}: {calls}");
            assert_eq!(decided(&cons, "HLA-B"), AMBIGUOUS_B, "{tag}");
        }
        let t1k = read_t1k(&hla, &OPTITYPE_GENES).unwrap();
        let rows = concordance_rows(&t1k, &read_optitype(&ot).unwrap());
        let b = &rows[1];
        assert_eq!(b[0], "HLA-B");
        assert_eq!(b[3], "not compared", "{tag}");
        assert!(
            b[4].starts_with("ambiguous: B*57:01 possible, B*58:01 possible"),
            "{tag}: {}",
            b[4]
        );
        assert!(
            b[4].ends_with("; step 36 does not pass it to PharmCAT"),
            "{tag}"
        );
    }
}

#[test]
fn a_tie_that_agrees_at_two_fields_still_passes() {
    let d = Scratch::new("same");
    let hla = d.put(
        "t.tsv",
        &t1k_row(
            "HLA-A",
            "HLA-A*02:01:01,HLA-A*02:01:02",
            60,
            "HLA-A*31:01:02",
            60,
        ),
    );
    let t1k = read_t1k(&hla, &HLA_GENES).unwrap();
    let row = hla_row("HLA-A", Some(&t1k), None, false);
    assert_eq!(row.decision, Decision::Passed("*02:01/*31:01".to_owned()));
    assert_eq!(row.reason, "T1K (step 08)");
}

#[test]
fn a_tie_without_a_drug_allele_keeps_a_sure_one_on_the_other_copy() {
    let d = Scratch::new("other");
    let mut seen = Vec::new();
    for tie in [
        "HLA-B*07:02:01,HLA-B*08:01:01",
        "HLA-B*08:01:01,HLA-B*07:02:01",
    ] {
        let hla = d.put("t.tsv", &t1k_row("HLA-B", "HLA-B*57:01:01", 60, tie, 60));
        let t1k = read_t1k(&hla, &HLA_GENES).unwrap();
        let row = hla_row("HLA-B", Some(&t1k), None, false);
        assert_eq!(
            row.decision,
            Decision::Passed("*57:01/*07:02".to_owned()),
            "{tie}"
        );
        assert_eq!(
            row.reason,
            "T1K (step 08); T1K ties *07:02 and *08:01 on one allele, passed as *07:02: \
             PharmCAT acts on none of them"
        );
        // OptiType's pair is one of T1K's choices: agreement, its pair passed.
        let ot = read_optitype(&d.put(
            "o.tsv",
            &optitype(("A*01:01", "A*01:01"), ("B*08:01", "B*57:01")),
        ))
        .unwrap();
        let row = hla_row("HLA-B", Some(&t1k), Some(&ot), true);
        assert_eq!(
            row.decision,
            Decision::Passed("*57:01/*08:01".to_owned()),
            "{tie}"
        );
        let conc = concordance_rows(&t1k, &ot);
        assert_eq!(conc[1][1], "*57:01/*07:02,*08:01");
        assert_eq!(conc[1][3], "yes");
        seen.push(row.reason);
    }
    assert_eq!(seen[0], seen[1]);
}

#[test]
fn a_null_allele_keeps_its_suffix() {
    assert_eq!(two_fields("HLA-A*02:53:01N", "HLA-A"), "*02:53N");
    assert_eq!(two_fields("HLA-A*24:09N", "HLA-A"), "*24:09N");
    assert_eq!(two_fields("HLA-A*01:01:01:02N", "HLA-A"), "*01:01N");
    assert_eq!(two_fields("HLA-B*57:01:01:02", "HLA-B"), "*57:01");
    let d = Scratch::new("null");
    let hla = d.put(
        "t.tsv",
        &t1k_row("HLA-A", "HLA-A*02:53:01N", 60, "HLA-A*31:01:02", 60),
    );
    let t1k = read_t1k(&hla, &OPTITYPE_GENES).unwrap();
    // OptiType names no null: its *02:53 is not T1K's *02:53N.
    let ot = read_optitype(&d.put(
        "o.tsv",
        &optitype(("A*02:53", "A*31:01"), ("B*07:02", "B*08:01")),
    ))
    .unwrap();
    let conc = concordance_rows(&t1k, &ot);
    assert_eq!(conc[0][1], "*02:53N/*31:01");
    assert_eq!(conc[0][3], "no");
    assert_eq!(
        conc[0][4],
        "T1K's call has an expression suffix, which OptiType never calls; \
         step 36 passes it to PharmCAT without the suffix"
    );
    // Step 36 compares OptiType with the names PharmCAT reads: the sure
    // *31:01 on the other copy reaches PharmCAT.
    for need in [false, true] {
        let row = hla_row("HLA-A", Some(&t1k), Some(&ot), need);
        assert_eq!(row.decision, Decision::Passed("*02:53/*31:01".to_owned()));
        assert_eq!(
            row.reason,
            "T1K and OptiType (step 08), which agree at two fields; \
             T1K's *02:53N passed as *02:53: PharmCAT acts on neither"
        );
    }
    // An OptiType that disagrees on the other copy still holds it back.
    let other = read_optitype(&d.put(
        "o.tsv",
        &optitype(("A*02:53", "A*01:01"), ("B*07:02", "B*08:01")),
    ))
    .unwrap();
    let row = hla_row("HLA-A", Some(&t1k), Some(&other), true);
    assert_eq!(row.decision, Decision::Withheld("indeterminate"));
    assert_eq!(row.reason, "T1K and OptiType disagree at two fields");
    assert_eq!(
        concordance_rows(&t1k, &other)[0][4],
        "the callers disagree at two fields; step 36 does not pass it to PharmCAT"
    );
    // T1K alone: PharmCAT reads no suffix, and *02:53 is no allele it acts
    // on; the sure *31:01 on the other copy reaches it.
    let row = hla_row("HLA-A", Some(&t1k), None, false);
    assert_eq!(row.decision, Decision::Passed("*02:53/*31:01".to_owned()));
    assert_eq!(
        row.reason,
        "T1K (step 08); T1K's *02:53N passed as *02:53: PharmCAT acts on neither"
    );
}

#[test]
fn a_null_of_a_drug_allele_is_ambiguous_but_a_null_beside_one_is_not() {
    let d = Scratch::new("null-drug");
    // A hypothetical null form of B*57:01.
    let hla = d.put(
        "t.tsv",
        &t1k_row("HLA-B", "HLA-B*57:01:99N", 60, "HLA-B*07:02:01", 60),
    );
    let t1k = read_t1k(&hla, &HLA_GENES).unwrap();
    let row = hla_row("HLA-B", Some(&t1k), None, false);
    assert_eq!(row.decision, Decision::Withheld("indeterminate"));
    assert_eq!(
        row.reason,
        "ambiguous: B*57:01 possible (T1K calls *57:01N, with an expression suffix)"
    );
    // A null on one copy says nothing about the other: B*57:01 still goes.
    let hla = d.put(
        "t.tsv",
        &t1k_row("HLA-B", "HLA-B*08:67:01N", 50, "HLA-B*57:01:01", 40),
    );
    let t1k = read_t1k(&hla, &HLA_GENES).unwrap();
    let row = hla_row("HLA-B", Some(&t1k), None, false);
    assert_eq!(row.decision, Decision::Passed("*08:67/*57:01".to_owned()));
    // OptiType, the opt-in second caller, never calls the null: its
    // *08:67 is the name PharmCAT reads, so B*57:01 still goes.
    let ot = read_optitype(&d.put(
        "o.tsv",
        &optitype(("A*01:01", "A*01:01"), ("B*57:01", "B*08:67")),
    ))
    .unwrap();
    let row = hla_row("HLA-B", Some(&t1k), Some(&ot), true);
    assert_eq!(row.decision, Decision::Passed("*08:67/*57:01".to_owned()));
    assert_eq!(
        row.reason,
        "T1K and OptiType (step 08), which agree at two fields; \
         T1K's *08:67N passed as *08:67: PharmCAT acts on neither"
    );
}

#[test]
fn one_allele_keeps_the_homozygous_caveat_and_agreement_its_reason() {
    let d = Scratch::new("hom");
    let hla = d.put("t.tsv", &t1k_row("HLA-B", "HLA-B*35:01:01", 60, ".", 0));
    let t1k = read_t1k(&hla, &HLA_GENES).unwrap();
    let row = hla_row("HLA-B", Some(&t1k), None, false);
    assert_eq!(row.decision, Decision::Passed("*35:01/*35:01".to_owned()));
    assert_eq!(row.reason, "T1K (step 08): one allele, read as homozygous");
    // Two alleles that agree: the reason the e2e cases match stays as it was.
    let hla = d.put(
        "t.tsv",
        &t1k_row("HLA-B", "HLA-B*35:01:01", 60, "HLA-B*07:02:01", 60),
    );
    let t1k = read_t1k(&hla, &HLA_GENES).unwrap();
    let ot = read_optitype(&d.put(
        "o.tsv",
        &optitype(("A*01:01", "A*01:01"), ("B*07:02", "B*35:01")),
    ))
    .unwrap();
    let row = hla_row("HLA-B", Some(&t1k), Some(&ot), true);
    assert_eq!(row.decision, Decision::Passed("*35:01/*07:02".to_owned()));
    assert_eq!(
        row.reason,
        "T1K and OptiType (step 08), which agree at two fields"
    );
}

#[test]
fn cyp2d6_agreement_is_not_withheld_for_pypgx_evidence_columns() {
    let d = Scratch::new("cyp");
    let pypgx = d.put(
        "pypgx.tsv",
        "Gene\tDiplotype\tPhenotype\tCNV_call\tSource\tStatus\tDefining_variants\tMin_GQ\tWeakest_variant\tDefining_filter\tQuality\n\
         CYP2D6\t*1/*4\tIntermediate Metabolizer\tNormal\tBAM\tcalled\t1/2 in the VCF\t8\tchr22:42128945 C>T GQ 8\t\
         not PASS: chr22:42128945 C>T RefCall\tlow-confidence: check the reads\n",
    );
    // The header and the call on separate lines (scripts/ci/check-personal-data.sh).
    let head = "Sample\tGenotype\tFilter\n";
    let cyrius = d.put("cyrius.tsv", &format!("{head}S1\t*4/*1\tPASS\n"));
    let depth = d.put(
        "depth.tsv",
        "metric\tvalue\nstatus\tok\nmessage\tCYP2D6 depth is fit\n",
    );
    let (calls, cons) = outside(
        &d,
        &[
            ("pypgx", pypgx.as_str()),
            ("cyrius", cyrius.as_str()),
            ("depth-check", depth.as_str()),
        ],
        false,
    );
    assert_eq!(calls, "CYP2D6\t*1/*4\n");
    assert_eq!(
        decided(&cons, "CYP2D6"),
        "CYP2D6\t*1/*4\tyes\tpypgx and Cyrius agree; the depth check passed"
    );
}
