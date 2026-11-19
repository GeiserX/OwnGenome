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
            "--verifybamid2-markers",
            "3000",
            "--out",
            &out,
        ],
    );
    assert_eq!((code, stderr.as_str()), (0, ""));
    let want = "key\tvalue\nsample\tS\nsomalier_id\tS\ndeclared_sex\tmale\ninferred_sex\tmale\nsex_check\tok\n\
sex_check_reason\tdeclared male, somalier infers male\nsites_genotyped\t120\ndepth_at_sites\t30.0\n\
x_sites\t39\nx_het\t0\nx_hom_alt\t39\nx_depth_ratio\t0.50\ny_sites\t5\ny_depth_ratio\t0.47\n\
freemix\t0.0412\nfreemix_warn_above\t0.03\ncontamination\twarn\npanel_markers\t10000\n\
verifybamid2_marker_check\tpassed\nsame_person_as\tT (0.95), U (0.99)\n\
verifybamid2_markers\t3000\ncontamination_confidence\tuncertain\n\
contamination_confidence_reason\tFREEMIX rests on 3000 markers; a firm estimate needs 5,000\n\
identity\tuncertain\nidentity_reason\tsomalier genotyped 120 sites; telling one person from another needs 1,000\n";
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

/// The verdict for `samples` (one row, sample S) with `pairs` and the
/// VerifyBamID2 inputs, as a map.
fn evidence(
    samples: &str,
    pairs: Option<&str>,
    selfsm: bool,
    marker_check: &str,
    markers: Option<u64>,
) -> HashMap<&'static str, String> {
    let d = Dir::new("qc-evidence");
    let s = d.put("s.samples.tsv", samples);
    let p = pairs.map(|p| d.put("s.pairs.tsv", p));
    let sm = selfsm.then(|| {
        d.put(
            "s.selfSM",
            "#SEQ_ID\tRG\tCHIP_ID\t#SNPS\t#READS\tAVG_DP\tFREEMIX\nS\tNA\tNA\t100000\t5000\t2.1\t0.0050\n",
        )
    });
    let i = Inputs {
        sample: "S",
        samples_tsv: &s,
        pairs_tsv: p.as_deref(),
        selfsm: sm.as_deref(),
        freemix_warn: FREEMIX_WARN,
        marker_check,
        markers,
        ..Inputs::default()
    };
    table(&i).unwrap().into_iter().collect()
}

/// A somalier row for S with these genotype counts.
fn sites(hom_ref: u32, het: u32, hom_alt: u32) -> String {
    samples(&MALE).replace(
        "\t60\t40\t20\t",
        &format!("\t{hom_ref}\t{het}\t{hom_alt}\t"),
    )
}

const PAIRS: &str = "#sample_a\tsample_b\trelatedness\tibs0\tibs2\tn\n";

#[test]
fn identity_on_300_sites_is_uncertain() {
    // 300 informative sites: a sliced or very shallow BAM.
    let t = evidence(&sites(100, 120, 80), None, false, "", None);
    assert_eq!(t["sites_genotyped"], "300");
    assert_eq!(t["identity"], "uncertain");
    assert_eq!(
        t["identity_reason"],
        "somalier genotyped 300 sites; telling one person from another needs 1,000"
    );
    // Also with another sample on the run that shares many sites.
    let pairs = format!("{PAIRS}S\tT\t0.02\t40\t300\t300\n");
    let t = evidence(&sites(100, 120, 80), Some(&pairs), false, "", None);
    assert_eq!(t["identity"], "uncertain");
    // 999 is below the floor, 1,000 is at it.
    let t = evidence(&sites(333, 333, 333), None, false, "", None);
    assert_eq!(t["identity"], "uncertain");
    let t = evidence(&sites(334, 333, 333), None, false, "", None);
    assert_eq!(t["identity"], "not_checked");
}

#[test]
fn identity_needs_another_sample_on_enough_shared_sites() {
    let many = sites(5000, 6000, 4000);
    // One sample alone: nothing to compare with.
    let t = evidence(&many, None, false, "", None);
    assert_eq!(t["identity"], "not_checked");
    assert!(t["identity_reason"].contains("no other sample"), "{t:?}");
    let t = evidence(&many, Some(PAIRS), false, "", None);
    assert_eq!(t["identity"], "not_checked");
    // A pair that shares 12,000 sites is a comparison; one of 300 is not.
    let pairs = format!("{PAIRS}S\tT\t0.02\t40\t300\t12000\nU\tS\t0.5\t0\t300\t300\n");
    let t = evidence(&many, Some(&pairs), false, "", None);
    assert_eq!(t["identity"], "checked");
    assert_eq!(
        t["identity_reason"],
        "compared with 1 other sample of the run on 1,000 or more shared sites"
    );
    let pairs = format!("{PAIRS}U\tS\t0.97\t0\t300\t300\nS\tV\t0.1\t0\t300\t999\n");
    let t = evidence(&many, Some(&pairs), false, "", None);
    assert_eq!(t["identity"], "uncertain");
    assert_eq!(
        t["identity_reason"],
        "no other sample of the run is known to share 1,000 or more genotyped sites with this one"
    );
    // The same-person flag still shows; the identity row says it is thin.
    assert_eq!(t["same_person_as"], "U (0.97)");
    // A pair that shares no site: somalier writes nan and n 0.
    let pairs = format!("{PAIRS}S\tT\tnan\t0\t0\t0\n");
    let t = evidence(&many, Some(&pairs), false, "", None);
    assert_eq!(t["identity"], "uncertain");
    // An older pairs table without n: the shared sites are not known.
    let pairs = "#sample_a\tsample_b\trelatedness\nS\tT\t0.1\nS\tU\t0.2\nS\tS\t1.0\n";
    let t = evidence(&many, Some(pairs), false, "", None);
    assert_eq!(t["identity"], "uncertain");
    // Two pairs on enough shared sites, and a self pair that does not count.
    let pairs = format!(
        "{PAIRS}S\tT\t0.1\t9\t300\t15000\nU\tS\t0.2\t9\t300\t1000\nS\tS\t1.0\t0\t300\t15000\n"
    );
    let t = evidence(&many, Some(&pairs), false, "", None);
    assert_eq!(
        (t["identity"].as_str(), t["identity_reason"].as_str()),
        (
            "checked",
            "compared with 2 other samples of the run on 1,000 or more shared sites"
        )
    );
}

#[test]
fn contamination_on_fewer_than_5000_markers_is_uncertain() {
    let many = sites(5000, 6000, 4000);
    // The marker check passed, yet only 4,999 markers were kept.
    let t = evidence(&many, None, true, "passed", Some(4999));
    assert_eq!(t["contamination"], "ok");
    assert_eq!(t["verifybamid2_markers"], "4999");
    assert_eq!(t["contamination_confidence"], "uncertain");
    assert_eq!(
        t["contamination_confidence_reason"],
        "FREEMIX rests on 4999 markers; a firm estimate needs 5,000"
    );
    let t = evidence(&many, None, true, "passed", Some(5000));
    assert_eq!(
        (
            t["contamination_confidence"].as_str(),
            t["contamination_confidence_reason"].as_str()
        ),
        ("ok", "FREEMIX rests on 5000 markers")
    );
    // Relaxed mode: uncertain whatever the count.
    for markers in [Some(320), Some(60000), None] {
        let t = evidence(&many, None, true, "skipped", markers);
        assert_eq!(t["contamination_confidence"], "uncertain", "{markers:?}");
        assert!(
            t["contamination_confidence_reason"].contains("skipped its marker check"),
            "{t:?}"
        );
    }
    // A passed check with no count: VerifyBamID2's own floor held.
    let t = evidence(&many, None, true, "passed", None);
    assert_eq!(t["contamination_confidence"], "ok");
    assert_eq!(t["verifybamid2_markers"], "");
    // Neither the check nor the count known.
    let t = evidence(&many, None, true, "", None);
    assert_eq!(t["contamination_confidence"], "uncertain");
    // VerifyBamID2 did not run: no confidence either way.
    let t = evidence(&many, None, false, "", Some(320));
    assert_eq!(
        (
            t["contamination"].as_str(),
            t["contamination_confidence"].as_str(),
            t["contamination_confidence_reason"].as_str()
        ),
        ("not_run", "", "")
    );
}

#[test]
fn the_markers_option() {
    let d = Dir::new("qc-markers");
    let samples = d.put("s.tsv", samples(&MALE));
    let selfsm = d.put("s.selfSM", "#SEQ_ID\tFREEMIX\nS\t0.001\n");
    let out = d.path("o.tsv");
    let go = |m: &str| {
        run(
            cli,
            &[
                "--sample",
                "S",
                "--somalier-samples",
                &samples,
                "--selfsm",
                &selfsm,
                "--marker-check",
                "passed",
                "--verifybamid2-markers",
                m,
                "--out",
                &out,
            ],
        )
    };
    // Empty: the log had no count.
    let (code, stdout, _) = go("");
    assert_eq!(code, 0);
    assert!(
        stdout.contains("  verifybamid2_markers: \n  contamination_confidence: ok\n"),
        "{stdout}"
    );
    let (code, _, stderr) = go("many");
    assert_eq!(code, 2);
    assert!(
        stderr.contains("argument --verifybamid2-markers: invalid int value: 'many'"),
        "{stderr}"
    );
}

#[test]
fn thousands_has_commas() {
    assert_eq!(thousands(0), "0");
    assert_eq!(thousands(999), "999");
    assert_eq!(thousands(1000), "1,000");
    assert_eq!(thousands(17766), "17,766");
    assert_eq!(thousands(1234567), "1,234,567");
}
