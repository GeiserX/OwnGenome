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

/// A one-sample VCF with these GT cells, one record each; gzipped as
/// step 25 and PRS_SCORE_SITES write it.
fn target(d: &Dir, name: &str, gts: &[&str]) -> String {
    let mut body = String::from(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
    );
    for (k, gt) in gts.iter().enumerate() {
        body.push_str(&format!(
            "1\t{}\t.\tA\tG\t.\tPASS\t.\tGT:GQ\t{gt}:30\n",
            k + 1
        ));
    }
    d.put_gz(name, &body)
}

/// og prs table on the panel results of [`results`] and the scores of
/// [`two_scores`], for this input kind and genotypes.
fn with_target(d: &Dir, kind: &str, gts: Option<&[&str]>) -> (u8, String, Vec<Vec<String>>) {
    let pgs = two_scores(d);
    results(d, "adj", true);
    let (res, out) = (d.path("adj"), d.path("out.tsv"));
    let t = gts.map(|g| target(d, "target.vcf.gz", g));
    let mut a = vec![
        "--sample",
        "S",
        "--sampleset",
        "sample",
        "--input-kind",
        kind,
        "--results",
        &res,
        "--scores",
        &pgs,
        "--panel",
        "p",
        "--out",
        &out,
    ];
    if let Some(t) = &t {
        a.extend(["--target", t.as_str()]);
    }
    let (code, stdout, stderr) = run(cli, &a);
    assert_eq!(stderr, "");
    let rows = std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect();
    (code, stdout, rows)
}

/// The cell of `col` in row `k` (0 is the header).
fn cell<'a>(rows: &'a [Vec<String>], k: usize, col: &str) -> &'a str {
    let i = rows[0].iter().position(|c| c == col).unwrap();
    &rows[k][i]
}

#[test]
fn a_variant_only_input_keeps_its_sum_and_loses_its_percentile() {
    // What a variant-only VCF gives pgsc_calc: het and hom-alt calls only.
    let d = Dir::new("prs-variants-only");
    let (code, stdout, rows) = with_target(&d, "vcf", Some(&["0/1", "1/1", "0|1", "1/1"]));
    assert_eq!(code, 0);
    let mut cols: Vec<&str> = COLUMNS[..12].to_vec();
    cols.extend(EFFECT_COLUMNS);
    cols.extend([NOTE_COLUMN, "Input"]);
    assert_eq!(rows[0], cols);
    assert_eq!(cell(&rows, 1, "Score_SUM"), "1.25");
    assert_eq!(cell(&rows, 1, "Percentile"), "NA");
    assert_eq!(cell(&rows, 1, "Effect"), "");
    assert_eq!(cell(&rows, 1, NOTE_COLUMN), WITHHELD);
    assert_eq!(cell(&rows, 1, "Input"), "vcf");
    // PGS000002 had no percentile to withhold.
    assert_eq!(cell(&rows, 2, NOTE_COLUMN), "");
    assert!(
        stdout.contains(
            "sum 1.25, 2 of 3 variants matched (66.7%), percentile withheld: the input lists variants only\n"
        ),
        "{stdout}"
    );
    assert!(!stdout.contains("percentile 94.4"), "{stdout}");
    assert!(
        stdout.ends_with(
            "NOTE: hom-ref sites are absent from this input (0 of 4 genotyped positions are \
             reference calls), so a missing site adds nothing to the sum and the sums are \
             biased; not comparable to published distributions, and no percentile is given\n"
        ),
        "{stdout}"
    );
}

#[test]
fn a_chip_shaped_input_keeps_its_percentile() {
    // scripts/chip-to-vcf.sh writes every probe: most are 0/0.
    let d = Dir::new("prs-chip");
    let gts = ["0/0", "0/0", "0/1", "0/0", "1/1", "0/0", "./.", "0/0"];
    let (code, stdout, rows) = with_target(&d, "vcf", Some(&gts));
    assert_eq!(code, 0);
    assert_eq!(cell(&rows, 1, "Percentile"), "94.4");
    assert_eq!(cell(&rows, 1, "Effect"), "HR 1.53 (EUR, PPM000591)");
    assert_eq!(cell(&rows, 1, NOTE_COLUMN), "");
    assert_eq!(cell(&rows, 1, "Input"), "vcf");
    assert!(stdout.contains("percentile 94.4 among EUR"), "{stdout}");
    assert!(
        !stdout.contains("NOTE:") && !stdout.contains("withheld"),
        "{stdout}"
    );
}

#[test]
fn the_genotypes_decide_never_the_input_kind() {
    // A gVCF whose reference blocks miss every site is as blind as a
    // variant-only VCF.
    let d = Dir::new("prs-kind");
    let (_, _, rows) = with_target(&d, "gvcf", Some(&["0/1", "1/1"]));
    assert_eq!(cell(&rows, 1, "Percentile"), "NA");
    assert_eq!(cell(&rows, 1, NOTE_COLUMN), WITHHELD);
    let (_, _, rows) = with_target(&d, "gvcf", Some(&["0/0", "0/1"]));
    assert_eq!(cell(&rows, 1, "Percentile"), "94.4");
    assert_eq!(cell(&rows, 1, NOTE_COLUMN), "");
    // Without --target there is no evidence and no note column: the table
    // of older callers.
    let (_, stdout, rows) = with_target(&d, "vcf", None);
    assert_eq!(rows[0].len(), 15);
    assert!(!rows[0].iter().any(|c| c == NOTE_COLUMN));
    assert_eq!(cell(&rows, 1, "Percentile"), "94.4");
    assert!(!stdout.contains("NOTE:"), "{stdout}");
}

#[test]
fn a_withheld_percentile_leaves_other_rows_alone() {
    // No panel: nothing to withhold, but the note on the input stays.
    let d = Dir::new("prs-raw-target");
    let pgs = two_scores(&d);
    results(&d, "raw", false);
    let res = d.path("raw");
    let t = target(&d, "target.vcf.gz", &["0/1"]);
    let (code, stdout, _, rows) = go(&d, &["--results", &res, "--scores", &pgs, "--target", &t]);
    assert_eq!(code, 0);
    assert_eq!(cell(&rows, 1, "Percentile"), "NA");
    assert_eq!(cell(&rows, 1, NOTE_COLUMN), "");
    assert!(
        stdout.contains("NOTE: hom-ref sites are absent"),
        "{stdout}"
    );
    // A score of the other sex says that, not the input's note.
    let d = Dir::new("prs-sex-target");
    formatted(
        &d,
        "pgs",
        "PGS000001",
        "Prostate cancer",
        3,
        "0\tleft out",
        "male",
        "",
    );
    results(&d, "adj", true);
    let (res, pgs) = (d.path("adj"), d.path("pgs"));
    let t = target(&d, "target.vcf.gz", &["1/1"]);
    let (_, _, _, rows) = go(
        &d,
        &[
            "--results",
            &res,
            "--scores",
            &pgs,
            "--panel",
            "p",
            "--sex",
            "female",
            "--target",
            &t,
        ],
    );
    assert_eq!(cell(&rows, 1, "Sex_Check"), SEX_NOT_APPLICABLE);
    assert_eq!(cell(&rows, 1, NOTE_COLUMN), "");
    // An unreadable target is an error, and no table.
    std::fs::remove_file(d.path("out.tsv")).unwrap();
    let (code, _, err, rows) = go(
        &d,
        &[
            "--results",
            &res,
            "--scores",
            &pgs,
            "--target",
            &d.path("none.vcf.gz"),
        ],
    );
    assert_eq!(code, 1);
    assert!(err.contains("No such file or directory"), "{err}");
    assert!(rows.is_empty());
}

#[test]
fn genotypes_count_called_and_reference_calls() {
    let d = Dir::new("prs-genotypes");
    let p = d.put(
        "t.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
         1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\n\
         1\t2\t.\tA\tG\t.\tPASS\t.\tGQ:GT\t30:0|0\n\
         X\t3\t.\tA\tG\t.\tPASS\t.\tGT\t0\n\
         1\t4\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n\
         1\t5\t.\tA\tG\t.\tPASS\t.\tGT\t./.\n\
         1\t6\t.\tA\tG\t.\tPASS\t.\tGT\t0/.\n\
         1\t7\t.\tA\tG\t.\tPASS\t.\tGQ\t30\n\
         1\t8\t.\tA\tG\t.\tPASS\t.\n\
         1\t9\t.\tA\tG\t.\tPASS\t.\tGT\t1/2\r\n",
    );
    assert_eq!(
        Genotypes::read(&p).unwrap(),
        Genotypes {
            called: 5,
            hom_ref: 3
        }
    );
    // One called site in five as 0/0 is enough; fewer, or none called, is not.
    let g = |called, hom_ref| Genotypes { called, hom_ref }.lists_reference();
    assert!(g(5, 1) && g(1, 1) && g(430_000, 300_000));
    assert!(!g(6, 1) && !g(4, 0) && !g(0, 0));
}

/// The panel run of [`results`], with the target's random-forest
/// probabilities and pgsc_calc's placement.
fn placed(d: &Dir, dir: &str, afr: &str, eur: &str, top: &str) {
    results(d, dir, true);
    d.put_gz(
        &format!("{dir}/sample/score/sample_popsimilarity.txt.gz"),
        &format!(
            "sampleset\tFID\tIID\tPC1\tRF_P_AFR\tRF_P_EUR\tMostSimilarPop\tMostSimilarPop_LowConfidence\tREFERENCE\n\
             reference\tR1\tR1\t-1\t0.9\t0.1\tAFR\tFalse\tTrue\n\
             sample\tHG\tHG\t0.012\t{afr}\t{eur}\t{top}\tFalse\tFalse\n"
        ),
    );
}

#[test]
fn a_sample_between_two_groups_is_told_the_percentile_ranks_it_among_one() {
    // pgsc_calc calls none of these low-confidence: its top probability is
    // 0.5 or more. Two groups at 0.2 or more each read as mixed.
    for (afr, eur, want) in [
        ("0.5", "0.5", true),
        ("0.45", "0.55", true),
        ("0.2", "0.8", true),
        ("0.19", "0.81", false),
        ("0.05", "0.95", false),
    ] {
        let d = Dir::new("prs-mixed");
        let pgs = two_scores(&d);
        placed(&d, "adj", afr, eur, "EUR");
        let (res, anc) = (d.path("adj"), d.path("anc.tsv"));
        let (code, stdout, _, rows) = go(
            &d,
            &[
                "--results",
                &res,
                "--scores",
                &pgs,
                "--panel",
                "p",
                "--ancestry-out",
                &anc,
            ],
        );
        assert_eq!(code, 0);
        // The ancestry table keeps pgsc_calc's own flag and probabilities.
        assert!(
            d.read("anc.tsv").contains(&format!(
                "population_low_confidence\tFalse\nprobability_AFR\t{afr}\nprobability_EUR\t{eur}\n"
            )),
            "{afr}/{eur}: {}",
            d.read("anc.tsv")
        );
        // The percentile itself stays: the note says what it is among.
        assert_eq!(rows[1][6..8], ["94.4", "EUR"]);
        let note = format!(
            "NOTE: the genetic ancestry sits between EUR (probability {eur}) and AFR (probability \
             {afr}), each at 0.2 or more, and each percentile ranks the sample among EUR samples \
             only, so it may not fit\n"
        );
        assert_eq!(stdout.contains(&note), want, "{afr}/{eur}: {stdout}");
        assert_eq!(stdout.contains("NOTE: the genetic ancestry"), want);
    }
}

#[test]
fn no_percentile_no_ancestry_note() {
    // A score for the other sex has no percentile, so nothing is ranked.
    let d = Dir::new("prs-mixed-none");
    formatted(
        &d,
        "pgs",
        "PGS000001",
        "Prostate cancer",
        3,
        "0\tleft out",
        "male",
        "",
    );
    placed(&d, "adj", "0.5", "0.5", "EUR");
    let (res, pgs, anc) = (d.path("adj"), d.path("pgs"), d.path("anc.tsv"));
    let (code, stdout, _, rows) = go(
        &d,
        &[
            "--results",
            &res,
            "--scores",
            &pgs,
            "--panel",
            "p",
            "--sex",
            "female",
            "--ancestry-out",
            &anc,
        ],
    );
    assert_eq!((code, rows[1][6].as_str()), (0, "NA"));
    assert!(!stdout.contains("NOTE"), "{stdout}");
}

#[test]
fn mixed_reads_the_two_most_probable_groups() {
    let p = |v: &[(&str, f64)]| -> Vec<(String, f64)> {
        v.iter().map(|(g, x)| ((*g).to_owned(), *x)).collect()
    };
    let three = p(&[("AFR", 0.25), ("EUR", 0.5), ("SAS", 0.25)]);
    // Ties go by name; the second group is never the placement.
    assert_eq!(mixed(&three, "EUR"), Some((true, "AFR".to_owned())));
    assert_eq!(mixed(&three, "AFR"), Some((true, "EUR".to_owned())));
    assert_eq!(
        mixed(&p(&[("AFR", 0.1), ("EUR", 0.8), ("SAS", 0.1)]), "EUR"),
        Some((false, "AFR".to_owned()))
    );
    // A cell that is not a number is not a group.
    assert_eq!(mixed(&p(&[("AFR", f64::NAN), ("EUR", 1.0)]), "EUR"), None);
    assert_eq!(mixed(&p(&[("EUR", 1.0)]), "EUR"), None);
    // The ancestry table's rows: only probability_ rows with a number count.
    assert_eq!(
        probabilities([
            ("population", "EUR"),
            ("probability_AFR", "0.45"),
            ("probability_EUR", "0.55"),
            ("probability_SAS", "NA"),
            ("PC1", "0.1"),
        ]),
        p(&[("AFR", 0.45), ("EUR", 0.55)])
    );
}

// --- percentile calibration ----------------------------------------------------
//
// og prs table reports the percentile pgsc_calc gave: the rank of the
// sample's sum among the sums of the reference group it was placed in. For
// people drawn from that group's own allele frequencies the reported
// percentiles must then be uniform from 0 to 100. These tests simulate the
// genotypes, compute the percentile the way pgsc_calc's ancestry step does
// (scipy's percentileofscore among the group's reference sums, pgscatalog-utils
// _ancestry/tools.py), write it into a pgsc_calc output folder with the panel's
// own rows beside the sample's, and read back what og prs table reports.

/// splitmix64: the same simulated people on every run.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in [0, 1).
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Score sites in the simulation.
const SITES: usize = 200;

/// One reference group's effect-allele frequencies, and the score's weights.
struct Group {
    freq: Vec<f64>,
    weight: Vec<f64>,
}

impl Group {
    fn new(rng: &mut Rng) -> Group {
        Group {
            freq: (0..SITES).map(|_| 0.05 + 0.9 * rng.unit()).collect(),
            weight: (0..SITES).map(|_| rng.unit() - 0.5).collect(),
        }
    }

    /// The sum of one person of the group: two alleles at each site, each the
    /// effect allele at the group's frequency.
    fn person(&self, rng: &mut Rng) -> f64 {
        self.freq
            .iter()
            .zip(&self.weight)
            .map(|(p, w)| {
                let dose = u8::from(rng.unit() < *p) + u8::from(rng.unit() < *p);
                w * f64::from(dose)
            })
            .sum()
    }
}

/// scipy's `percentileofscore(reference, x)` (kind "rank").
fn percentile_of(reference: &[f64], x: f64) -> f64 {
    let left = reference.iter().filter(|r| **r < x).count();
    let right = reference.iter().filter(|r| **r <= x).count();
    let plus1 = usize::from(right > left);
    (left + right + plus1) as f64 * 50.0 / reference.len() as f64
}

/// The p-value of a Kolmogorov-Smirnov test of `u` against the uniform
/// distribution on [0, 1] (Kolmogorov's limit with Stephens' small-sample
/// correction). A percentile is a rank among `reference` people, a sample
/// of the group itself, so the test is the two-sample one: its size is
/// n * reference / (n + reference). `None` for a known uniform.
fn ks_uniform(u: &mut [f64], reference: Option<usize>) -> f64 {
    u.sort_by(f64::total_cmp);
    let n = u.len() as f64;
    let d = u
        .iter()
        .enumerate()
        .map(|(i, x)| ((i + 1) as f64 / n - x).max(x - i as f64 / n))
        .fold(0.0, f64::max);
    let size = reference.map_or(n, |r| n * r as f64 / (n + r as f64));
    let l = (size.sqrt() + 0.12 + 0.11 / size.sqrt()) * d;
    if l < 0.2 {
        return 1.0; // the series is 1 to many digits there, and converges slowly
    }
    let p: f64 = (1..=100i32)
        .map(|k| {
            let sign = if k % 2 == 1 { 2.0 } else { -2.0 };
            let k = f64::from(k);
            sign * (-2.0 * k * k * l * l).exp()
        })
        .sum();
    p.clamp(0.0, 1.0)
}

/// What og prs table reports for each person: pgsc_calc's panel run is
/// written with the person's sum and percentile among `reference`, the
/// panel's own rows after it, and read back.
fn reported(tag: &str, people: &[f64], reference: &[f64]) -> Vec<f64> {
    let d = Dir::new(tag);
    formatted(
        &d,
        "pgs",
        "PGS000001",
        "Simulated trait",
        SITES,
        "0\tleft out",
        "",
        "",
    );
    d.put(
        "res/sample/match/sample_summary.csv",
        format!("dataset,accession,match_status,count\nsample,PGS000001,matched,{SITES}\n"),
    );
    d.put(
        "res/sample/score/sample_popsimilarity.txt.gz",
        "sampleset\tFID\tIID\tRF_P_AFR\tRF_P_EUR\tMostSimilarPop\tMostSimilarPop_LowConfidence\tREFERENCE\n\
         sample\tS\tS\t0.02\t0.98\tEUR\tFalse\tFalse\n",
    );
    let head = "sampleset\tFID\tIID\tPGS\tSUM\tZ_MostSimilarPop\tpercentile_MostSimilarPop\tZ_norm1\tZ_norm2\n";
    let row = |set: &str, id: &str, x: f64| {
        format!(
            "{set}\t{id}\t{id}\tPGS000001\t{x:.6}\t0\t{}\t0\t0\n",
            percentile_of(reference, x)
        )
    };
    let panel: String = reference
        .iter()
        .take(20)
        .enumerate()
        .map(|(k, x)| row("reference", &format!("R{k}"), *x))
        .collect();
    let (res, scores) = (d.path("res"), d.path("pgs"));
    people
        .iter()
        .map(|x| {
            d.put(
                "res/sample/score/sample_pgs.txt.gz",
                format!("{head}{}{panel}", row("sample", "S", *x)),
            );
            let inputs = Inputs {
                results: &res,
                sampleset: "sample",
                scores: &scores,
                input_kind: "gvcf",
                outcome: Outcome::Scored,
                sex: "",
                panel: true,
                genotypes: None,
            };
            let (rows, _) = table(&inputs).unwrap();
            assert_eq!(rows[0][7], "EUR");
            rows[0][6].parse::<f64>().unwrap()
        })
        .collect()
}

/// The KS p-value of what og prs table reports for 1,000 people of one
/// group, scored against 2,000 reference people of the group whose sums are
/// moved by `shift` standard deviations.
fn calibration(tag: &str, shift: f64) -> f64 {
    let mut rng = Rng(0x2026_1120);
    let group = Group::new(&mut rng);
    let mut reference: Vec<f64> = (0..2000).map(|_| group.person(&mut rng)).collect();
    let people: Vec<f64> = (0..1000).map(|_| group.person(&mut rng)).collect();
    let n = reference.len() as f64;
    let mean = reference.iter().sum::<f64>() / n;
    let sd = (reference.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
    for x in &mut reference {
        *x += shift * sd;
    }
    let mut u: Vec<f64> = reported(tag, &people, &reference)
        .iter()
        .map(|p| p / 100.0)
        .collect();
    ks_uniform(&mut u, Some(reference.len()))
}

#[test]
fn percentiles_of_people_from_the_reference_group_are_uniform() {
    let p = calibration("prs-calibrated", 0.0);
    assert!(p > 0.01, "KS p = {p}");
}

#[test]
fn a_reference_group_half_a_standard_deviation_off_fails_the_check() {
    // The same people against a reference whose mean sits 0.5 SD higher, as
    // when a person is ranked among a group their ancestry does not match:
    // the check must see it.
    let p = calibration("prs-shifted", 0.5);
    assert!(p < 0.01, "KS p = {p}");
}

#[test]
fn the_ks_test_reads_known_samples() {
    // An evenly spread sample is uniform; one squeezed into the lower half is not.
    let mut even: Vec<f64> = (0..1000).map(|i| (f64::from(i) + 0.5) / 1000.0).collect();
    assert!((ks_uniform(&mut even, None) - 1.0).abs() < 1e-9);
    let mut low: Vec<f64> = (0..1000).map(|i| f64::from(i) / 2000.0).collect();
    assert!(ks_uniform(&mut low, None) < 1e-6);
    // D = 0.05 over 1,000 values: p 0.013 against a known uniform, 0.15
    // against a reference of 1,000 people (the two-sample test).
    let mut off: Vec<f64> = (0..1000)
        .map(|i| (f64::from(i) + 0.5) / 1000.0 * 0.95)
        .collect();
    let (known, sampled) = (
        ks_uniform(&mut off.clone(), None),
        ks_uniform(&mut off, Some(1000)),
    );
    assert!((0.01..0.016).contains(&known), "{known}");
    assert!((0.13..0.18).contains(&sampled), "{sampled}");
}
