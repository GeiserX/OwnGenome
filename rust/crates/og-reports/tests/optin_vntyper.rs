//! `og vntyper report` (step 52): the reading of VNtyper's MUC1 VNTR
//! result. Moved from tests/test_vntyper_report.py (its cases 1 to 5 and
//! 7); the summary's reading is in og-summary's `tests/optin_summary.rs`
//! and the report lines in og-render's `tests/optin_report.rs`.
//!
//! The run folders are made up per case, in the layout VNtyper 2.0.47
//! writes: `kestrel/kestrel_result.tsv` after its `##` banner lines, and
//! `coverage/coverage_summary.tsv` as one header and one data line.

use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

const FIELDS: [&str; 13] = [
    "status",
    "confidence",
    "variant",
    "flag",
    "alt_depth",
    "region_depth",
    "depth_score",
    "vntr_mean_depth",
    "vntr_uncovered_pct",
    "coverage_qc",
    "sensitivity",
    "message",
    "advice",
];
const ADVICE_NONE: &str = "A negative here lowers the chance of the frameshifts VNtyper looks for; it does not rule ADTKD-MUC1 out. A family history of kidney failure with gout or a bland urine sediment still deserves a clinical test (SNaPshot for dupC).";
const BANNER: &str = "##VNtyper Kestrel result\n##source=test\n";
const POS_HEADER: &str = "Motifs\tVariant\tPOS\tREF\tALT\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\tFlag\thaplo_count\tNomenclature\tMotif_sequence\n";
const NEG: &str = "Motif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\nNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNegative\n";
const COV_HEADER: &str = "mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n";

fn call(conf: &str, flag: &str) -> String {
    format!(
        "{POS_HEADER}2\tInsertion\t67\tG\tGC\t416\t7110\t0.0585\t{conf}\t{flag}\t3\t59dupC\tCCCCACC\n"
    )
}

fn hp() -> String {
    call("High_Precision*", "Not flagged")
}

struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Dir {
        let d = std::env::temp_dir().join(format!("og-optin-vntyper-{}-{tag}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        Dir(d)
    }
    fn p(&self, n: &str) -> String {
        self.0.join(n).to_str().unwrap().to_owned()
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// A run folder; `mean`, `qc` and `uncovered` go into the coverage line.
fn run_dir(
    t: &Dir,
    name: &str,
    kestrel: Option<&str>,
    mean: &str,
    qc: &str,
    uncovered: &str,
) -> String {
    let d = t.p(name);
    fs::create_dir_all(format!("{d}/kestrel")).unwrap();
    fs::create_dir_all(format!("{d}/coverage")).unwrap();
    if let Some(k) = kestrel {
        fs::write(
            format!("{d}/kestrel/kestrel_result.tsv"),
            format!("{BANNER}{k}"),
        )
        .unwrap();
    }
    fs::write(
        format!("{d}/coverage/coverage_summary.tsv"),
        format!("{COV_HEADER}{mean}\t1000.00\t300.00\t95\t1942\t4500\t0\t{uncovered}\t{qc}\n"),
    )
    .unwrap();
    d
}

fn og(args: &[&str]) -> (u8, String) {
    let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
    let (mut o, mut e) = (Vec::new(), Vec::new());
    let code = og_reports::vntyper::cli(Some("report"), &args, &mut o, &mut e);
    (code, String::from_utf8(e).unwrap())
}

type Reading = HashMap<String, String>;

/// The verb on run folder `d`: the table as a map.
fn og_read(d: &str, exit_code: i64, log_text: &str) -> Reading {
    let (out, log) = (format!("{d}.tsv"), format!("{d}.log"));
    fs::write(&log, log_text).unwrap();
    let (code, err) = og(&[
        "--run-dir",
        d,
        "--exit-code",
        &exit_code.to_string(),
        "--log",
        &log,
        "--out",
        &out,
    ]);
    assert_eq!(code, 0, "{err}");
    fs::read_to_string(&out)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| {
            let (k, v) = l.split_once('\t').unwrap();
            (k.to_owned(), v.to_owned())
        })
        .collect()
}

/// The control: a reading that calls every Kestrel call positive.
fn every_call_positive(d: &str, exit_code: i64, log_text: &str) -> Reading {
    let mut r = og_read(d, exit_code, log_text);
    if r["status"] == "possible" {
        r.insert("status".into(), "positive".into());
    }
    r
}

type Read = fn(&str, i64, &str) -> Reading;

/// (description, held, reading) of cases 1 to 4, read with `read`.
fn cases(read: Read, t: &Dir) -> Vec<(&'static str, bool, Reading)> {
    let rd = |name, k: Option<&str>, mean, qc, unc| run_dir(t, name, k, mean, qc, unc);
    let mut out = Vec::new();
    let r = read(&rd("pos", Some(&hp()), "1075.32", "PASS", "0.00"), 0, "");
    out.push((
        "1. High_Precision*, not flagged, 1,075x: positive with the name",
        r["status"] == "positive"
            && r["variant"] == "59dupC"
            && r["confidence"] == "High_Precision*"
            && r["alt_depth"] == "416",
        r.clone(),
    ));
    out.push((
        "1. the sensitivity line names the depth and the 100x band (160 of 161)",
        r["sensitivity"].contains("1075x") && r["sensitivity"].contains("160 of 161 runs (99%)"),
        r,
    ));
    let r = read(
        &rd(
            "low",
            Some(&call("Low_Precision", "Not flagged")),
            "1075.32",
            "PASS",
            "0.00",
        ),
        0,
        "",
    );
    out.push((
        "2. a Low_Precision call is possible, not positive",
        r["status"] == "possible" && r["message"].contains("low precision"),
        r,
    ));
    let r = read(
        &rd(
            "flag",
            Some(&call("High_Precision*", "False_Positive_4bp_Insertion")),
            "1075.32",
            "PASS",
            "0.00",
        ),
        0,
        "",
    );
    out.push((
        "2. a flagged High_Precision call is possible, with the flag",
        r["status"] == "possible" && r["message"].contains("False_Positive_4bp_Insertion"),
        r,
    ));
    let r = read(&rd("holes", Some(&hp()), "1075.32", "FAIL", "62.00"), 0, "");
    out.push((
        "2. a call over a VNTR 62% without reads is possible",
        r["status"] == "possible" && r["message"].contains("62% of the VNTR has no read"),
        r,
    ));
    let r = read(&rd("nounc", Some(&hp()), "1075.32", "PASS", "NA"), 0, "");
    out.push((
        "2. a call with no uncovered share from VNtyper is possible, not positive",
        r["status"] == "possible" && r["message"].contains("no share of the VNTR without reads"),
        r,
    ));
    let r = read(&rd("neg", Some(NEG), "147.94", "PASS", "0.00"), 0, "");
    out.push((
        "3. no call at 148x, coverage PASS: negative, 100x band",
        r["status"] == "negative"
            && r["variant"] == "."
            && r["sensitivity"].contains("160 of 161")
            && r["advice"] == ADVICE_NONE,
        r,
    ));
    let r = read(&rd("red", Some(NEG), "70.00", "REDUCED", "0.00"), 0, "");
    out.push((
        "3. no call at 70x, coverage REDUCED: negative, 50x-100x band (95 of 108)",
        r["status"] == "negative"
            && r["sensitivity"].contains("95 of 108 runs (88%)")
            && r["sensitivity"].contains("between 50x and 100x"),
        r,
    ));
    let r = read(&rd("thin", Some(NEG), "20.00", "FAIL", "0.00"), 0, "");
    out.push((
        "3. no call at 20x, coverage FAIL: not assessed, below-50x band (48 of 109)",
        r["status"] == "not assessed"
            && r["message"].contains("cannot be ruled out")
            && r["sensitivity"].contains("48 of 109 runs (44%)"),
        r,
    ));
    let r = read(&rd("na", Some(NEG), "NA", "NOT_EVALUATED", "0.00"), 0, "");
    out.push((
        "3. no measured depth: not assessed, sensitivity not known",
        r["status"] == "not assessed" && r["sensitivity"].starts_with("Sensitivity not known"),
        r,
    ));
    let r = read(&rd("rc", Some(&hp()), "1075.32", "PASS", "0.00"), 1, "");
    out.push((
        "4. VNtyper's exit code 1: not assessed, whatever it wrote",
        r["status"] == "not assessed" && r["message"].contains("exit code 1"),
        r,
    ));
    let r = read(
        &rd("empty", None, "1075.32", "PASS", "0.00"),
        1,
        "ValueError: FASTQ layout 'empty' cannot be consumed without dropping reads.",
    );
    out.push((
        "4. exit 1 on a BAM with no read over MUC1: not assessed, saying so",
        r["status"] == "not assessed" && r["message"].contains("no read in the MUC1 region"),
        r,
    ));
    let header_only = format!("{}\n", NEG.split('\n').next().unwrap());
    let r = read(
        &rd("hdr", Some(&header_only), "147.94", "PASS", "0.00"),
        0,
        "",
    );
    out.push((
        "4. a Kestrel file with a header and no result row: not assessed, never negative",
        r["status"] == "not assessed" && r["message"].contains("no result row"),
        r,
    ));
    let r = read(&rd("nok", None, "1075.32", "PASS", "0.00"), 0, "");
    out.push((
        "4. no kestrel_result.tsv: not assessed, naming the file",
        r["status"] == "not assessed" && r["message"].contains("kestrel/kestrel_result.tsv"),
        r,
    ));
    out
}

#[test]
fn readings_of_each_kestrel_and_coverage_result() {
    let t = Dir::new("real");
    let failed: Vec<String> = cases(og_read, &t)
        .into_iter()
        .filter(|c| !c.1)
        .map(|c| format!("{}: {:?}", c.0, c.2))
        .collect();
    assert!(failed.is_empty(), "{failed:#?}");
}

#[test]
fn control_a_reading_that_calls_every_call_positive_fails_case_2() {
    let t = Dir::new("broken");
    let failed: Vec<&str> = cases(every_call_positive, &t)
        .into_iter()
        .filter(|c| !c.1)
        .map(|c| c.0)
        .collect();
    assert!(failed.iter().any(|x| x.starts_with("2.")), "{failed:?}");
}

#[test]
fn the_cli_writes_every_metric_once_and_refuses_an_unreadable_result() {
    let t = Dir::new("cli");
    let d = run_dir(&t, "cli", Some(&hp()), "1075.32", "PASS", "0.00");
    let out = t.p("cli.tsv");
    let (code, err) = og(&["--run-dir", &d, "--out", &out]);
    assert_eq!(code, 0, "{err}");
    let text = fs::read_to_string(&out).unwrap();
    let rows: Vec<Vec<&str>> = text.lines().map(|l| l.split('\t').collect()).collect();
    assert_eq!(rows[0], ["metric", "value"]);
    let names: Vec<&str> = rows[1..].iter().map(|r| r[0]).collect();
    assert_eq!(names, FIELDS);
    assert!(rows.iter().all(|r| r.len() == 2));
    let bad = run_dir(
        &t,
        "bad",
        Some("Confidence\n\"x\n"),
        "1075.32",
        "PASS",
        "0.00",
    );
    fs::write(
        format!("{bad}/kestrel/kestrel_result.tsv"),
        b"\xff\xfe\x00bad",
    )
    .unwrap();
    let (code, err) = og(&["--run-dir", &bad, "--out", &t.p("bad.tsv")]);
    assert!(code == 2 && err.contains("cannot read"), "{code} {err}");
}
