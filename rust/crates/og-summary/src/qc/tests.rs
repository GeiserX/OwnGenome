use super::*;
use crate::testutil::{Dir, run};

const HEAD: &str = "#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\toriginal_pedigree_sex\t\
gt_depth_mean\tn_hom_ref\tn_het\tn_hom_alt\tp_middling_ab\tX_depth_mean\tX_n\tX_hom_ref\t\
X_het\tX_hom_alt\tY_depth_mean\tY_n\n";

struct Case<'a> {
    sex: &'a str,
    ped: &'a str,
    x_n: &'a str,
    x_het: &'a str,
    x_hom_alt: &'a str,
    y_depth: &'a str,
    mid: &'a str,
}

const MALE: Case = Case {
    sex: "1",
    ped: "-9",
    x_n: "39",
    x_het: "0",
    x_hom_alt: "39",
    y_depth: "14.0",
    mid: "0.010",
};

fn samples(c: &Case) -> String {
    format!(
        "{HEAD}S\tS\t-9\t-9\t{}\t-9\t{}\t30.0\t60\t40\t20\t{}\t15.0\t{}\t0\t{}\t{}\t{}\t5\n",
        c.sex, c.ped, c.mid, c.x_n, c.x_het, c.x_hom_alt, c.y_depth
    )
}

/// The verdict for one somalier row, as a map.
fn verdict(c: &Case, declared: Option<Sex>) -> HashMap<&'static str, String> {
    let d = Dir::new("qc-verdict");
    let path = d.put("s.samples.tsv", samples(c));
    let i = Inputs {
        sample: "S",
        samples_tsv: &path,
        declared,
        freemix_warn: FREEMIX_WARN,
        ..Inputs::default()
    };
    table(&i).unwrap().into_iter().collect()
}

#[test]
fn somalier_pedigree_sex_is_a_call_only_when_chrx_agrees() {
    // The checks of step 33 with the declared sex as somalier's pedigree.
    let c = Case {
        ped: "male",
        x_n: "2",
        x_hom_alt: "2",
        ..MALE
    };
    let t = verdict(&c, Some(Sex::Male));
    assert_eq!(
        (t["inferred_sex"].as_str(), t["sex_check"].as_str()),
        ("unknown", "not_checked")
    );
    let c = Case {
        ped: "female",
        ..MALE
    };
    let t = verdict(&c, Some(Sex::Female));
    assert_eq!(
        (t["inferred_sex"].as_str(), t["sex_check"].as_str()),
        ("male", "mismatch")
    );
    let c = Case {
        ped: "male",
        ..MALE
    };
    assert_eq!(verdict(&c, Some(Sex::Male))["sex_check"], "ok");
    // somalier calls the sample low quality: 8% middling allele balance, or nan.
    for mid in ["0.080", "nan"] {
        let c = Case {
            ped: "male",
            mid,
            ..MALE
        };
        assert_eq!(
            verdict(&c, Some(Sex::Male))["sex_check"],
            "not_checked",
            "{mid}"
        );
    }
    // No pedigree (-9): somalier's sex column is read as before.
    let t = verdict(
        &Case {
            x_n: "2",
            x_hom_alt: "2",
            ..MALE
        },
        Some(Sex::Male),
    );
    assert_eq!(
        (t["inferred_sex"].as_str(), t["sex_check"].as_str()),
        ("male", "ok")
    );
}

#[test]
fn x_and_y_disagree() {
    // A pedigree female with chrY reads whose chrX cannot tell: unknown.
    let c = Case {
        sex: "-2",
        ped: "female",
        x_n: "4",
        x_hom_alt: "4",
        ..MALE
    };
    let t = verdict(&c, Some(Sex::Female));
    assert_eq!(t["sex_check"], "not_checked");
    assert!(t["sex_check_reason"].contains("could not tell"), "{t:?}");
    // chrX heterozygous and chrY reads: the aneuploidy line, declared or not.
    let c = Case {
        sex: "-2",
        ped: "female",
        x_n: "40",
        x_het: "20",
        x_hom_alt: "20",
        ..MALE
    };
    for declared in [Some(Sex::Female), None] {
        let t = verdict(&c, declared);
        assert_eq!(t["sex_check"], "not_checked");
        assert!(t["sex_check_reason"].contains("chrY has reads"), "{t:?}");
    }
    let c = Case { sex: "-2", ..MALE };
    assert!(verdict(&c, None)["sex_check_reason"].contains("chrY has reads"));
}

#[test]
fn x_sex_is_somaliers_own_rule() {
    let row = |n: &str, het: &str, hom: &str, mid: Option<&str>| -> Row {
        let mut r: Row = [("X_n", n), ("X_het", het), ("X_hom_alt", hom)]
            .into_iter()
            .map(|(k, v)| (k.to_owned(), v.to_owned()))
            .collect();
        if let Some(m) = mid {
            r.insert("p_middling_ab".into(), m.into());
        }
        r
    };
    assert_eq!(x_sex(&row("11", "0", "11", None)), Some(Sex::Male));
    assert_eq!(
        x_sex(&row("10", "0", "10", None)),
        None,
        "10 sites are not enough"
    );
    assert_eq!(
        x_sex(&row("40", "20", "20", Some("0.01"))),
        Some(Sex::Female)
    );
    assert_eq!(
        x_sex(&row("40", "5", "20", None)),
        None,
        "0.25 is in between"
    );
    assert_eq!(x_sex(&row("40", "1", "0", None)), Some(Sex::Female));
    assert_eq!(x_sex(&row("40", "0", "0", None)), None);
    assert_eq!(x_sex(&row("40", "0", "20", Some("0.06"))), None);
    assert_eq!(x_sex(&row("40", "0", "20", Some("x"))), None);
    assert_eq!(x_sex(&row("nan", "0", "20", None)), None);
}

#[test]
fn the_table_byte_for_byte() {
    let d = Dir::new("qc-table");
    let samples = d.put("s.samples.tsv", samples(&MALE));
    let selfsm = d.put(
        "s.selfSM",
        "#SEQ_ID\tRG\tCHIP_ID\t#SNPS\t#READS\tAVG_DP\tFREEMIX\nS\tNA\tNA\t10000\t5000\t2.1\t0.0412\n",
    );
    let pairs = d.put(
        "s.pairs.tsv",
        "#sample_a\tsample_b\trelatedness\nS\tT\t0.95\nU\tS\t0.991\nS\tV\t0.5\nS\tS\t1.0\nW\tS\tnan\n",
    );
    let out = d.path("out.tsv");
    let (code, stdout, stderr) = run(
        cli,
        &[
            "--sample",
            "S",
            "--somalier-samples",
            &samples,
            "--selfsm",
            &selfsm,
            "--somalier-pairs",
            &pairs,
            "--declared-sex",
            "male",
            "--marker-check",
            "passed",
            "--out",
            &out,
        ],
    );
    assert_eq!((code, stderr.as_str()), (0, ""));
    let want = "key\tvalue\nsample\tS\nsomalier_id\tS\ndeclared_sex\tmale\ninferred_sex\tmale\nsex_check\tok\n\
sex_check_reason\tdeclared male, somalier infers male\nsites_genotyped\t120\ndepth_at_sites\t30.0\n\
x_sites\t39\nx_het\t0\nx_hom_alt\t39\nx_depth_ratio\t0.50\ny_sites\t5\ny_depth_ratio\t0.47\n\
freemix\t0.0412\nfreemix_warn_above\t0.03\ncontamination\twarn\npanel_markers\t10000\n\
verifybamid2_marker_check\tpassed\nsame_person_as\tT (0.95), U (0.99)\n";
    assert_eq!(d.read("out.tsv"), want);
    let echoed: String = want
        .lines()
        .skip(1)
        .map(|l| format!("  {}\n", l.replacen('\t', ": ", 1)))
        .collect();
    assert_eq!(stdout, echoed);
}

#[test]
fn the_ratios_are_chrx_and_chry_over_the_autosomes() {
    // Red control for a swap of the two ratios: X 15/30, Y 14/30.
    let t = verdict(&MALE, None);
    assert_eq!(
        (t["x_depth_ratio"].as_str(), t["y_depth_ratio"].as_str()),
        ("0.50", "0.47")
    );
    let t = verdict(
        &Case {
            y_depth: "1.5",
            ..MALE
        },
        None,
    );
    assert_eq!(t["y_depth_ratio"], "0.05");
}

#[test]
fn the_threshold_and_the_somalier_name() {
    let d = Dir::new("qc-threshold");
    let samples = d.put("s.tsv", samples(&MALE).replace("S\tS\t", "S\tSM_x\t"));
    let selfsm = d.put("s.selfSM", "#SEQ_ID\tFREEMIX\nS\t0.05\n");
    let out = d.path("o.tsv");
    let args = [
        "--sample",
        "S",
        "--somalier-id",
        "SM_x",
        "--somalier-samples",
        &samples,
        "--selfsm",
        &selfsm,
        "--freemix-warn",
        "0.08",
        "--out",
        &out,
    ];
    let (code, stdout, _) = run(cli, &args);
    assert_eq!(code, 0);
    assert!(stdout.contains("  somalier_id: SM_x\n"), "{stdout}");
    assert!(
        stdout.contains("  freemix_warn_above: 0.08\n  contamination: ok\n"),
        "{stdout}"
    );
    assert!(stdout.contains("  panel_markers: \n"), "{stdout}");
    // --sample is required.
    let (code, _, stderr) = run(cli, &args[2..]);
    assert_eq!(code, 2);
    assert!(stderr.contains("required: --sample"), "{stderr}");
}

#[test]
fn unreadable_inputs_exit_1_with_a_message() {
    let d = Dir::new("qc-bad");
    let path = d.put("s.tsv", samples(&MALE));
    let out = d.path("o.tsv");
    let go = |extra: &[&str]| {
        let mut a = vec!["--sample", "Q", "--somalier-samples", &path, "--out", &out];
        a.extend_from_slice(extra);
        run(cli, &a)
    };
    let (code, _, err) = go(&[]);
    assert_eq!(
        (code, err.as_str()),
        (1, "ERROR: s.tsv has no row for sample 'Q' (it has: S)\n")
    );
    let selfsm = d.put("e.selfSM", "#SEQ_ID\tFREEMIX\nQ\tNA\n");
    let (code, _, err) = go(&["--somalier-id", "S", "--selfsm", &selfsm]);
    assert_eq!(
        (code, err.as_str()),
        (1, "ERROR: e.selfSM has no FREEMIX value\n")
    );
    // Python stopped on an OverflowError traceback here.
    let inf = d.put(
        "inf.tsv",
        samples(&MALE).replace("\t60\t40\t20\t", "\tinf\t40\t20\t"),
    );
    let a = ["--sample", "S", "--somalier-samples", &inf, "--out", &out];
    let (code, _, err) = run(cli, &a);
    assert_eq!(code, 1);
    assert!(
        err.starts_with("ERROR: inf.tsv: n_hom_ref of 'S' is not a count"),
        "{err}"
    );
    let (code, _, err) = go(&["--freemix-warn", "x"]);
    assert_eq!(code, 2);
    assert!(err.contains("invalid float value: 'x'"), "{err}");
    assert!(!std::path::Path::new(&out).exists());
}
