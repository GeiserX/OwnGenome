//! `og pgx outside-calls` with OptiType, the opt-in second HLA caller.
//!
//! With `HLA_QC=true` (bash) or `--hla_qc` (Nextflow) step 08 also runs
//! OptiType on the HLA-A, -B and -C reads, and step 36 passes HLA-A or HLA-B
//! to PharmCAT only when T1K and OptiType give the same two alleles at two
//! fields. Each case runs the verb on synthetic files:
//!
//! 1. agree: both callers agree, HLA-A and HLA-B reach PharmCAT;
//! 2. disagree: OptiType differs at HLA-B, HLA-B stays 'indeterminate', both
//!    calls are printed, it is not in the outside-call file; HLA-A still
//!    passes;
//! 3. homozygous: T1K's one allele and OptiType's A1 == A2 agree;
//! 4. no call: a file step 08 wrote instead of a result (`# OptiType not
//!    run: ...`), both genes withheld;
//! 5. asked, missing: `--require-optitype` without `--optitype`, both
//!    withheld;
//! 6. T1K alone: no OptiType asked for, T1K's call passes, as before;
//! 7. quality 0: T1K's quality-0 rule still wins over an agreement;
//! 8. concordance: the `hla-concordance` verb, HLA-A, -B, -C, yes/no;
//! 9. control: the same inputs without OptiType's result, so without the
//!    agreement rule: the disagreeing HLA-B reaches the outside-call file
//!    and case 2's checks report it, so case 2 can fail (the crate's own
//!    planted control is the unit test
//!    `optitype_agreeing_passes_and_disagreeing_withholds`).

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use og_steps::pgx_rules::{hla_concordance_cli, outside_calls_cli};

const T1K: &str = "HLA-A\t2\tHLA-A*01:01:01\t110.5\t60\tHLA-A*26:01:01\t95.2\t60\t\n\
                   HLA-B\t2\tHLA-B*57:01:01\t101.0\t60\tHLA-B*38:01:01\t88.7\t48\t\n\
                   HLA-C\t1\tHLA-C*06:02:01\t150.3\t60\t.\t0\t-1\t\n";
const T1K_HOM: &str = "HLA-A\t1\tHLA-A*02:01:01\t200.0\t60\t.\t0\t-1\t\n\
                       HLA-B\t2\tHLA-B*57:01:01\t101.0\t60\tHLA-B*38:01:01\t88.7\t48\t\n";
const T1K_Q0: &str = "HLA-A\t2\tHLA-A*01:01:01\t110.5\t60\tHLA-A*26:01:01\t95.2\t60\t\n\
                      HLA-B\t2\tHLA-B*57:01:01\t101.0\t0\tHLA-B*38:01:01\t88.7\t48\t\n";
const OT_HEAD: &str = "\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n";
const OT_AGREE: &str = "0\tA*26:01\tA*01:01\tB*38:01\tB*57:01\tC*06:02\tC*06:02\t812\t790.4\n";
const OT_DISAGREE: &str = "0\tA*01:01\tA*26:01\tB*58:01\tB*38:01\tC*06:02\tC*07:01\t812\t790.4\n";
const OT_HOM: &str = "0\tA*02:01\tA*02:01\tB*57:01\tB*38:01\t\t\t500\t480.0\n";
const OT_NONE: &str = "# OptiType not run: no reads at HLA-A, HLA-B or HLA-C in this BAM\n";
const BOTH: &str = "T1K and OptiType (step 08), which agree at two fields";

/// A scratch directory of one test: emptied when made, removed when
/// dropped.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Self {
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("pgx-hla-{name}"));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }

    fn file(&self, name: &str, text: &str) -> String {
        let p = self.path(name);
        fs::write(&p, text).unwrap();
        p
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).to_string_lossy().into_owned()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn ot(body: &str) -> String {
    format!("{OT_HEAD}{body}")
}

/// One run of the verb: the outside-call text and the consensus rows by
/// gene, or the error when it failed.
struct Run {
    calls: Option<String>,
    rows: BTreeMap<String, Vec<String>>,
}

fn run(tmp: &Scratch, t1k: &str, optitype: Option<&str>, require: bool) -> Run {
    let hla = tmp.file("t1k.tsv", t1k);
    let (calls, cons) = (tmp.path("calls.tsv"), tmp.path("cons.tsv"));
    let _ = fs::remove_file(&calls);
    let _ = fs::remove_file(&cons);
    let mut args = vec![
        "--calls".to_owned(),
        calls.clone(),
        "--consensus".to_owned(),
        cons.clone(),
        "--hla".to_owned(),
        hla,
    ];
    if let Some(text) = optitype {
        args.push("--optitype".to_owned());
        args.push(tmp.file("optitype.tsv", text));
    }
    if require {
        args.push("--require-optitype".to_owned());
    }
    let (mut out, mut err) = (Vec::new(), Vec::new());
    if outside_calls_cli(&args, &mut out, &mut err) != 0 {
        let mut rows = BTreeMap::new();
        rows.insert(
            "error".to_owned(),
            vec![String::from_utf8_lossy(&err).into_owned()],
        );
        return Run { calls: None, rows };
    }
    let rows = fs::read_to_string(&cons)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| {
            let c: Vec<String> = l.split('\t').map(str::to_owned).collect();
            (c[0].clone(), c)
        })
        .collect();
    Run {
        calls: Some(fs::read_to_string(&calls).unwrap()),
        rows,
    }
}

impl Run {
    fn passed(&self, gene: &str, dip: &str) -> bool {
        self.calls
            .as_deref()
            .is_some_and(|c| c.contains(&format!("{gene}\t{dip}\n")))
    }

    fn row(&self, gene: &str) -> Vec<String> {
        self.rows.get(gene).cloned().unwrap_or_default()
    }

    fn cell(&self, gene: &str, i: usize) -> String {
        self.row(gene).get(i).cloned().unwrap_or_default()
    }

    /// The gene is not in the outside-call file, its row is
    /// 'indeterminate', 'no', and its reason passes `reason`.
    fn withheld(&self, gene: &str, reason: impl Fn(&str) -> bool) -> bool {
        let r = self.row(gene);
        self.calls
            .as_deref()
            .is_some_and(|c| !c.lines().any(|l| l.starts_with(&format!("{gene}\t"))))
            && r.get(1..3) == Some(&["indeterminate".to_owned(), "no".to_owned()][..])
            && r.get(3).is_some_and(|s| reason(s))
    }
}

/// Case 2's checks; the names of the ones that failed.
fn disagree_checks(tmp: &Scratch, optitype: Option<&str>) -> Vec<&'static str> {
    let r = run(tmp, T1K, optitype, false);
    let mut failed = Vec::new();
    if !r.withheld("HLA-B", |s| s == "T1K and OptiType disagree at two fields") {
        failed.push("HLA-B withheld when OptiType disagrees");
    }
    if r.cell("HLA-B", 4)
        != "T1K: HLA-B*57:01:01 (quality 60); HLA-B*38:01:01 (quality 48); OptiType: *58:01/*38:01"
    {
        failed.push("both HLA-B calls are printed");
    }
    if !r.passed("HLA-A", "*01:01/*26:01") {
        failed.push("HLA-A, where they agree, still passes");
    }
    failed
}

#[test]
fn case1_both_callers_agree() {
    let tmp = Scratch::new("agree");
    let r = run(&tmp, T1K, Some(&ot(OT_AGREE)), false);
    assert!(
        r.passed("HLA-A", "*01:01/*26:01"),
        "1 HLA-A passes when both agree (either order): {:?}",
        r.calls
    );
    assert!(
        r.passed("HLA-B", "*57:01/*38:01"),
        "1 HLA-B passes when both agree: {:?}",
        r.calls
    );
    assert_eq!(r.cell("HLA-B", 3), BOTH, "1 the reason names both callers");
    assert!(
        r.calls.as_deref().is_some_and(|c| !c.contains("HLA-C")),
        "1 HLA-C is never an outside call: {:?}",
        r.calls
    );
}

#[test]
fn case2_optitype_disagrees_at_hla_b() {
    let tmp = Scratch::new("disagree");
    let failed = disagree_checks(&tmp, Some(&ot(OT_DISAGREE)));
    assert!(failed.is_empty(), "2 failed: {failed:?}");
}

#[test]
fn case3_homozygous_agreement() {
    let tmp = Scratch::new("hom");
    let r = run(&tmp, T1K_HOM, Some(&ot(OT_HOM)), false);
    assert!(
        r.passed("HLA-A", "*02:01/*02:01"),
        "3 T1K's one allele agrees with OptiType's A1 == A2: {:?}",
        r.calls
    );
    assert_eq!(
        r.cell("HLA-A", 3),
        BOTH,
        "3 the homozygous HLA-A row says both callers"
    );
}

#[test]
fn case4_optitype_made_no_call() {
    let tmp = Scratch::new("nocall");
    let r = run(&tmp, T1K, Some(OT_NONE), false);
    for g in ["HLA-A", "HLA-B"] {
        assert!(
            r.withheld(g, |s| s == "OptiType, the second caller, made no call"),
            "4 {g} withheld: OptiType made no call: {:?} {:?}",
            r.calls,
            r.row(g)
        );
    }
}

#[test]
fn case5_asked_for_without_a_result() {
    let tmp = Scratch::new("missing");
    let r = run(&tmp, T1K, None, true);
    for g in ["HLA-A", "HLA-B"] {
        assert!(
            r.withheld(g, |s| s
                .starts_with("the second caller was asked for (HLA_QC)")),
            "5 {g} withheld: OptiType asked for, no result: {:?} {:?}",
            r.calls,
            r.row(g)
        );
    }
}

#[test]
fn case6_t1k_alone_as_before() {
    let tmp = Scratch::new("alone");
    let r = run(&tmp, T1K, None, false);
    assert!(
        r.passed("HLA-B", "*57:01/*38:01"),
        "6 without OptiType, T1K's HLA-B passes: {:?}",
        r.calls
    );
    assert_eq!(r.cell("HLA-B", 3), "T1K (step 08)", "6 the reason is T1K's");
}

#[test]
fn case7_t1k_quality_0_still_withholds() {
    let tmp = Scratch::new("q0");
    let r = run(&tmp, T1K_Q0, Some(&ot(OT_AGREE)), false);
    assert!(
        r.withheld("HLA-B", |s| s.contains("quality 0 or below")),
        "7 a T1K quality-0 allele withholds HLA-B even when OptiType agrees: {:?} {:?}",
        r.calls,
        r.row("HLA-B")
    );
}

#[test]
fn case8_the_concordance_table_byte_for_byte() {
    let tmp = Scratch::new("concordance");
    let t1k = tmp.file("c_t1k.tsv", T1K);
    let ot_path = tmp.file("c_ot.tsv", &ot(OT_DISAGREE));
    let out = tmp.path("c_out.tsv");
    let args: Vec<String> = ["--t1k", &t1k, "--optitype", &ot_path, "--out", &out]
        .map(str::to_owned)
        .to_vec();
    let (mut so, mut se) = (Vec::new(), Vec::new());
    let code = hla_concordance_cli(&args, &mut so, &mut se);
    let want = "Gene\tT1K\tOptiType\tAgree\tNote\n\
                HLA-A\t*01:01/*26:01\t*01:01/*26:01\tyes\t-\n\
                HLA-B\t*57:01/*38:01\t*58:01/*38:01\tno\tthe callers disagree at two fields; step 36 does not pass it to PharmCAT\n\
                HLA-C\t*06:02/*06:02\t*06:02/*07:01\tno\tthe callers disagree at two fields\n";
    let got = if code == 0 {
        fs::read_to_string(&out).unwrap_or_default()
    } else {
        String::from_utf8_lossy(&se).into_owned()
    };
    assert_eq!(got, want, "8 the concordance table, byte for byte");
}

#[test]
fn case9_control_without_the_agreement_rule_case2_fails() {
    let tmp = Scratch::new("control");
    let failed = disagree_checks(&tmp, None);
    assert!(
        !failed.is_empty(),
        "9 with the agreement rule disabled, case 2 rejects the disagreeing HLA-B"
    );
}
