use super::*;
use crate::testutil::{Dir, run};

/// The folder og prs format writes for one score.
#[allow(clippy::too_many_arguments)]
fn formatted(
    d: &Dir,
    dir: &str,
    pid: &str,
    trait_: &str,
    rows: usize,
    chrx: &str,
    applies: &str,
    evals: &str,
) {
    let body: String = "1\t1\tA\tG\t1\n".repeat(rows);
    d.put_gz(
        &format!("{dir}/{pid}.txt.gz"),
        &format!(
            "#pgs_id={pid}\n#trait_reported={trait_}\n#genome_build=GRCh38\n\
             chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n{body}"
        ),
    );
    d.put(
        &format!("{dir}/chrx_rows.tsv"),
        format!("pgs_id\tchrx_rows\tchrx\n{pid}\t{chrx}\n"),
    );
    d.put(
        &format!("{dir}/score_notes.tsv"),
        format!("pgs_id\tapplies_to\tambiguous_rows\tambiguous\n{pid}\t{applies}\t0\tdropped\n"),
    );
    let ev = if evals.is_empty() {
        String::new()
    } else {
        format!("{pid}\t{evals}\n")
    };
    d.put(
        &format!("{dir}/score_evaluations.tsv"),
        format!("pgs_id\tppm_id\tancestry\teffect\n{ev}"),
    );
}

/// A pgsc_calc --outdir for sampleset `sample`: PGS000001 matched 2 of 3,
/// PGS000002 dropped. The reference's rows come after the target's.
fn results(d: &Dir, dir: &str, panel: bool) {
    d.put(
        &format!("{dir}/sample/match/sample_summary.csv"),
        "dataset,accession,match_status,count\nsample,PGS000001,matched,2\n\
         sample,PGS000001,unmatched,1\nsample,PGS000002,unmatched,5\n",
    );
    if !panel {
        d.put_gz(
            &format!("{dir}/sample/score/aggregated_scores.txt.gz"),
            "sampleset\tFID\tIID\tPGS\tSUM\tDENOM\nsample\tHG\tHG\tPGS000001\t62.0\t6\n",
        );
        return;
    }
    d.put_gz(
        &format!("{dir}/sample/score/sample_pgs.txt.gz"),
        "sampleset\tFID\tIID\tPGS\tSUM\tpercentile_MostSimilarPop\n\
         sample\tHG\tHG\tPGS000001\t1.25\t94.44\nreference\tR1\tR1\tPGS000001\t10.0\t12.3\n",
    );
    d.put_gz(
        &format!("{dir}/sample/score/sample_popsimilarity.txt.gz"),
        "sampleset\tFID\tIID\tPC1\tPC2\tPC10\tRF_P_AFR\tRF_P_EUR\tMostSimilarPop\tMostSimilarPop_LowConfidence\tREFERENCE\n\
         reference\tR1\tR1\t-1\t2\t0\t0.9\t0.1\tAFR\tFalse\tTrue\n\
         sample\tHG\tHG\t0.012\t-0.034\t0.5\t0.05\t0.95\tEUR\tFalse\tFalse\n",
    );
}

/// The scores folder of two scores: PGS000001 with one chrX row left out
/// and an EUR evaluation, PGS000002 with none of the side files' rows.
fn two_scores(d: &Dir) -> String {
    formatted(
        d,
        "pgs",
        "PGS000001",
        "Coronary artery disease",
        3,
        "1\tleft out",
        "",
        "PPM000591;PPM000594\tEUR;AFR\tHR 1.53;HR 1.27",
    );
    d.put_gz(
        "pgs/PGS000002.txt.gz",
        "#pgs_id=PGS000002\n#genome_build=GRCh38\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n\
         1\t1\tA\t\t1\n1\t2\tA\t\t1\n1\t3\tA\t\t1\n1\t4\tA\t\t1\n1\t5\tA\t\t1\n",
    );
    d.path("pgs")
}

fn go(d: &Dir, args: &[&str]) -> (u8, String, String, Vec<Vec<String>>) {
    let out = d.path("out.tsv");
    let mut a = vec![
        "--sample",
        "S",
        "--sampleset",
        "sample",
        "--input-kind",
        "gvcf",
        "--out",
        &out,
    ];
    a.extend_from_slice(args);
    let (code, stdout, stderr) = run(cli, &a);
    let rows = std::fs::read_to_string(&out)
        .map(|t| {
            t.lines()
                .map(|l| l.split('\t').map(str::to_owned).collect())
                .collect()
        })
        .unwrap_or_default();
    (code, stdout, stderr, rows)
}

#[test]
fn no_panel_raw_sums() {
    let d = Dir::new("prs-raw");
    let pgs = two_scores(&d);
    results(&d, "raw", false);
    let res = d.path("raw");
    let (code, stdout, stderr, rows) = go(&d, &["--results", &res, "--scores", &pgs]);
    assert_eq!((code, stderr.as_str()), (0, ""));
    let cols: Vec<&str> = COLUMNS[..12]
        .iter()
        .chain(&EFFECT_COLUMNS)
        .chain(&COLUMNS[12..])
        .copied()
        .collect();
    assert_eq!(rows[0], cols);
    assert_eq!(
        rows[1],
        [
            "Coronary artery disease",
            "PGS000001",
            "62",
            "2",
            "3",
            "66.7",
            "NA",
            "NA",
            "1 left out",
            "0",
            "",
            "",
            "EUR;AFR",
            "",
            "gvcf"
        ]
    );
    assert_eq!(
        rows[2],
        [
            "PGS000002",
            "PGS000002",
            "NA",
            "0",
            "5",
            "0.0",
            "NA",
            "NA",
            "",
            "",
            "",
            "",
            "",
            "",
            "gvcf"
        ]
    );
    assert_eq!(
        stdout,
        "  Coronary artery disease (PGS000001): sum 62, 2 of 3 variants matched (66.7%)\n  \
         PGS000002 (PGS000002): sum NA, 0 of 5 variants matched (0.0%)\n"
    );
}

#[test]
fn panel_percentile_group_effect_and_ancestry() {
    let d = Dir::new("prs-panel");
    let pgs = two_scores(&d);
    results(&d, "adj", true);
    let (res, anc) = (d.path("adj"), d.path("anc.tsv"));
    let (code, stdout, _, rows) = go(
        &d,
        &[
            "--results",
            &res,
            "--scores",
            &pgs,
            "--panel",
            "pgsc_1000G_v1",
            "--ancestry-out",
            &anc,
        ],
    );
    assert_eq!(code, 0);
    // The target's values, not the reference sample's.
    assert_eq!(rows[1][2..8], ["1.25", "2", "3", "66.7", "94.4", "EUR"]);
    assert_eq!(rows[1][13], "HR 1.53 (EUR, PPM000591)");
    assert!(
        stdout
            .contains("percentile 94.4 among EUR, effect in that group HR 1.53 (EUR, PPM000591)\n"),
        "{stdout}"
    );
    assert_eq!(
        d.read("anc.tsv"),
        "key\tvalue\nsample\tS\nreference_panel\tpgsc_1000G_v1\npopulation\tEUR\npopulation_low_confidence\tFalse\n\
         probability_AFR\t0.05\nprobability_EUR\t0.95\nPC1\t0.012\nPC2\t-0.034\nPC10\t0.5\n"
    );
}

#[test]
fn a_score_for_the_other_sex_keeps_its_sum_and_loses_its_percentile() {
    let d = Dir::new("prs-sex");
    formatted(
        &d,
        "pgs",
        "PGS000001",
        "Prostate cancer",
        3,
        "0\tleft out",
        "male",
        "PPM000591\tEUR\tHR 1.53",
    );
    results(&d, "adj", true);
    let (res, pgs) = (d.path("adj"), d.path("pgs"));
    let base = [
        "--results",
        res.as_str(),
        "--scores",
        pgs.as_str(),
        "--panel",
        "p",
    ];
    let row = |sex: &str| {
        let mut a = base.to_vec();
        a.extend(["--sex", sex]);
        let (code, stdout, _, rows) = go(&d, &a);
        assert_eq!(code, 0);
        (rows[1].clone(), stdout)
    };
    let (f, out) = row("female");
    assert_eq!(
        (f[2].as_str(), f[6].as_str(), f[11].as_str(), f[13].as_str()),
        ("1.25", "NA", SEX_NOT_APPLICABLE, "")
    );
    assert!(
        out.contains(", male-specific: not applicable to the declared sex, no percentile\n"),
        "{out}"
    );
    let (m, _) = row("male");
    assert_eq!(
        (m[6].as_str(), m[11].as_str(), m[13].as_str()),
        ("94.4", SEX_APPLIES, "HR 1.53 (EUR, PPM000591)")
    );
    let (n, out) = row("");
    assert_eq!((n[6].as_str(), n[11].as_str()), ("94.4", SEX_NOT_DECLARED));
    assert!(
        out.contains(", male-specific score: no sex declared"),
        "{out}"
    );
}

#[test]
fn zero_matches_and_below_threshold_are_rows() {
    let d = Dir::new("prs-none");
    let pgs = two_scores(&d);
    let none = d.path("none");
    let (code, _, _, rows) = go(
        &d,
        &["--results", &none, "--scores", &pgs, "--zero-matches"],
    );
    assert_eq!(code, 0);
    assert_eq!(rows[1][2..6], ["NA", "0", "3", "0.0"]);
    let log = d.put(
        "below.log",
        "ERROR Score PGS000001 fails minimum matching threshold (33.33% variants match)\n\
         Score  fails minimum matching threshold (1% variants match)\nScore PGS000002 fails minimum matching threshold (x% variants match)\n",
    );
    let (code, stdout, _, rows) = go(
        &d,
        &[
            "--results",
            &none,
            "--scores",
            &pgs,
            "--below-threshold",
            &log,
        ],
    );
    assert_eq!(code, 0);
    assert_eq!(rows[1][2..6], ["NA", "NA", "3", "33.3"]);
    assert_eq!(rows[2][2..6], ["NA", "NA", "5", "NA"]);
    assert!(stdout.starts_with(
        "  Coronary artery disease (PGS000001): no sum, under pgsc_calc's minimum overlap (33.3% of 3 variants matched)\n"
    ));
    // Neither flag and no pgsc_calc output: an error, and no table.
    std::fs::remove_file(d.path("out.tsv")).unwrap();
    let (code, _, err, rows) = go(&d, &["--results", &none, "--scores", &pgs]);
    assert_eq!(code, 1);
    assert!(
        err.starts_with("ERROR: pgsc_calc wrote no match summary ("),
        "{err}"
    );
    assert!(rows.is_empty());
    let (code, _, err, _) = go(
        &d,
        &["--results", &none, "--scores", &none, "--zero-matches"],
    );
    assert_eq!(
        (code, err),
        (1, format!("ERROR: no formatted scoring file in {none}\n"))
    );
}

#[test]
fn a_panel_run_names_the_chrx_rows_its_matches_did_not_keep() {
    let d = Dir::new("prs-chrx");
    formatted(
        &d,
        "pgs",
        "PGS000001",
        "Coronary artery disease",
        4,
        "1\tscored",
        "",
        "",
    );
    results(&d, "adj", true);
    let (res, pgs) = (d.path("adj"), d.path("pgs"));
    let head =
        "row_nr,accession,chr_name,chr_position,effect_allele,other_allele,match_status,dataset\n";
    let chrx = |log: &str, n: &str| {
        d.put_gz(
            "adj/sample/match/sample_log.csv.gz",
            &format!("{head}0,PGS000001,1,1100,A,G,matched,sample\n{log}"),
        );
        d.put(
            "pgs/chrx_rows.tsv",
            format!("pgs_id\tchrx_rows\tchrx\nPGS000001\t{n}\tscored\n"),
        );
        let (code, _, _, rows) = go(&d, &["--results", &res, "--scores", &pgs, "--panel", "p"]);
        assert_eq!(code, 0);
        rows[1][8].clone()
    };
    let x = "1,PGS000001_hmPOS_GRCh38,X,51000000,A,G,";
    assert_eq!(
        chrx(&format!("{x}excluded,sample\n"), "1"),
        "1 left out (not matched in the panel run)"
    );
    assert_eq!(chrx(&format!("{x}matched,sample\n"), "1"), "1 scored");
    assert_eq!(
        chrx(&format!("{x}matched,sample\n"), "3"),
        "1 scored, 2 left out (not matched in the panel run)"
    );
    // Without --panel the log is not read.
    let (_, _, _, rows) = go(&d, &["--results", &res, "--scores", &pgs]);
    assert_eq!(rows[1][8], "3 scored");
}

#[test]
fn a_list_without_evaluations_keeps_the_older_columns() {
    let d = Dir::new("prs-plain");
    formatted(
        &d,
        "pgs",
        "PGS000001",
        "Coronary artery disease",
        3,
        "0\tleft out",
        "",
        "",
    );
    results(&d, "raw", false);
    let (res, pgs) = (d.path("raw"), d.path("pgs"));
    let (code, _, _, rows) = go(&d, &["--results", &res, "--scores", &pgs]);
    assert_eq!(code, 0);
    assert_eq!(rows[0], COLUMNS);
    assert_eq!(rows[1].len(), 13);
}

#[test]
fn the_shipped_score_list_reads() {
    // assets/pgs_scores.tsv, as og prs format copies its evaluations.
    let assets = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../../assets/pgs_scores.tsv"
    );
    let l = Labels::read(Path::new(assets)).unwrap();
    let ev: Vec<(&str, &str, &str)> = l
        .evals("PGS000018")
        .iter()
        .map(|e| (e.ancestry.as_str(), e.effect.as_str(), e.ppm.as_str()))
        .collect();
    assert_eq!(
        ev,
        [
            ("EUR", "HR 1.53", "PPM000591"),
            ("AFR", "HR 1.27", "PPM000594"),
            ("AMR", "HR 1.53", "PPM000597")
        ]
    );
    let groups: Vec<&str> = l
        .evals("PGS002727")
        .iter()
        .map(|e| e.ancestry.as_str())
        .collect();
    assert_eq!(groups, ["EUR", "AFR", "AMR", "EAS", "SAS"]);
    assert!(l.evals("PGS002802").is_empty() && l.evals("PGS000065").is_empty());
}

#[test]
fn usage_errors_exit_2() {
    let d = Dir::new("prs-usage");
    let (code, _, err, _) = go(&d, &["--results", "r", "--scores", "s", "--sex", "other"]);
    assert_eq!(code, 2);
    assert!(
        err.contains("argument --sex: invalid choice: 'other'"),
        "{err}"
    );
}
