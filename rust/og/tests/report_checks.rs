//! The checks of the old `tests/test_render_report.py`, moved here when the
//! pipeline's reports moved to `og report render`: the PRS path around
//! pgsc_calc (`og prs format`, `og prs table`, then og summary collect), and
//! the cards of both reports on small synthetic sample folders.
//!
//!  1. og prs format writes a harmonised file as pgsc_calc's custom GRCh38
//!     file: chr_name/chr_position from hm_chr/hm_pos without "chr", the
//!     other allele from other_allele or hm_inferOtherAllele (never "A/C"),
//!     rows without a GRCh38 position dropped, chrX and chrY rows dropped
//!     from the score and the allele list, the label from the list; it
//!     refuses a GRCh37 file, a non-additive one and a listed id with no
//!     file; with --keep-x the chrX rows stay and chrY still goes.
//!  2. prs-table without a panel; 3. with a panel (the target sample's row
//!     only, and the ancestry table); 4. --zero-matches and
//!     --below-threshold; 4c. a sex-specific score from the real
//!     assets/pgs_scores.tsv, in the table and both reports, with the red
//!     control of a list without its applies_to column; 4d. strand-ambiguous
//!     rows; 4e. a chrX row and the panel's match log.
//!  5. both reports' PRS card: the percentile with its group, the raw-score
//!     line, the withheld percentiles of a variant-only input, the chrX
//!     note and the ancestry-transfer caveat; 5b. the fine-grained ancestry
//!     section beside it.
//!  6. the ClinVar card: no green badge, the frequency request, and step
//!     06's cDNA read check.
//!  7. the telomere card; 8. the mosaic card; 9. the paralog card;
//!  10. the blood group card; 12. the raw-data export; 13. the couple card.
//!      (11, the traits card, is genome_traits_report.rs.)
//!  14. the chapters, on DEMO-001 (og dev demo-sample): no off-machine
//!      address but the marked links, no green in the ClinVar section, each
//!      text chapter holds its sections, the card titles are the golden list,
//!      the "New since" badge only where a result changed, each caveat once
//!      under its card title. Each detector also runs on a planted bad input.
//!      The frame's own rules are in og-render's port_report_core.rs.

use std::collections::{BTreeMap, HashMap};
use std::io::{Read, Write};
use std::num::NonZero;
use std::path::PathBuf;
use std::process::{Command, Output};
use std::time::{Duration, UNIX_EPOCH};

use og_render::escape::e;
use og_render::notes::{EH_NOTE, MOSAIC_NOTE, NOTES, PRS_ANCESTRY_NOTE};
use og_summary::collect::{Group, Input};
use pycompat::json::{self, Dump, Value};

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");

const PRS_NOT_ADJUSTED: &str =
    "Polygenic score percentiles: no ancestry reference panel was used, so the scores are raw sums";
const PRS_WITHHELD: &str = "Polygenic score percentiles: the input lists variants only, so the percentiles were \
                            withheld and the scores are raw sums";

/// The columns of og prs table without evaluations: Input last.
const PRS_COLUMNS: [&str; 13] = [
    "Condition",
    "PGS_ID",
    "Score_SUM",
    "Variants_Matched",
    "Variants_Total",
    "Matched_Pct",
    "Percentile",
    "Ancestry_Group",
    "ChrX",
    "Ambiguous",
    "Applies_To",
    "Sex_Check",
    "Input",
];

const HM: &str = "#pgs_id=PGS000001\n#trait_reported=Something else\n#HmPOS_build=GRCh38\n\
rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\thm_chr\thm_pos\thm_inferOtherAllele\n\
rs1\t1\t100\tA\tG\t0.5\t1\t1100\t\n\
rs2\t1\t200\tc\t\t-0.25\tchr2\t2200\tT\n\
rs3\t1\t300\tG\t\t1\t3\t3300\tA/C\n\
rs4\t1\t400\tT\tC\t2\t\t\t\n\
rs5\tX\t500\tA\tG\t1\tX\t51000000\t\n\
rs6\tY\t600\tC\tT\t1\tchrY\t2800000\t\n";

const ROWS1: [&str; 3] = [
    "1\t1100\tA\tG\t0.5",
    "2\t2200\tC\tT\t-0.25",
    "3\t3300\tG\t\t1",
];

fn repo(rel: &str) -> String {
    format!("{REPO}/{rel}")
}

struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Scratch {
        let d =
            std::env::temp_dir().join(format!("og-report-checks-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        Scratch(d)
    }

    /// The absolute path of `rel` under the scratch folder.
    fn p(&self, rel: &str) -> String {
        self.0.join(rel).to_string_lossy().into_owned()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn parent(path: &str) {
    std::fs::create_dir_all(std::path::Path::new(path).parent().unwrap()).unwrap();
}

fn put(path: &str, text: &str) {
    parent(path);
    std::fs::write(path, text).unwrap();
}

fn put_gz(path: &str, text: &str) {
    parent(path);
    let f = std::fs::File::create(path).unwrap();
    let mut w = og_io::bgzf::Writer::new(f, NonZero::<usize>::MIN).unwrap();
    w.write_all(text.as_bytes()).unwrap();
    w.finish().unwrap();
}

/// A file as text, gzip read through.
fn read(path: &str) -> String {
    let mut s = String::new();
    og_io::open(path)
        .unwrap_or_else(|e| panic!("{path}: {e}"))
        .read_to_string(&mut s)
        .unwrap();
    s
}

fn og(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_og"))
        .args(args)
        .output()
        .expect("og runs")
}

fn rc(args: &[&str]) -> i32 {
    og(args).status.code().unwrap_or(-1)
}

fn prs_table(argv: &[&str]) -> i32 {
    let mut a = vec!["prs", "table"];
    a.extend(argv);
    rc(&a)
}

fn og_format(argv: &[&str]) -> i32 {
    let mut a = vec!["prs", "format"];
    a.extend(argv);
    rc(&a)
}

/// The score folder og prs format writes for one score.
#[allow(clippy::too_many_arguments)]
fn formatted(
    out: &str,
    pid: &str,
    trait_: &str,
    rows: &[&str],
    nx: u32,
    chrx: &str,
    applies: &str,
    amb: &str,
    evals: &str,
) {
    let body: String = rows.iter().map(|r| format!("{r}\n")).collect();
    put_gz(
        &format!("{out}/{pid}.txt.gz"),
        &format!(
            "#pgs_id={pid}\n#pgs_name={pid}\n#trait_reported={trait_}\n#genome_build=GRCh38\n\
             chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n{body}"
        ),
    );
    put(
        &format!("{out}/chrx_rows.tsv"),
        &format!("pgs_id\tchrx_rows\tchrx\n{pid}\t{nx}\t{chrx}\n"),
    );
    put(
        &format!("{out}/score_notes.tsv"),
        &format!("pgs_id\tapplies_to\tambiguous_rows\tambiguous\n{pid}\t{applies}\t{amb}\n"),
    );
    let ev = if evals.is_empty() {
        String::new()
    } else {
        format!("{pid}\t{evals}\n")
    };
    put(
        &format!("{out}/score_evaluations.tsv"),
        &format!("pgs_id\tppm_id\tancestry\teffect\n{ev}"),
    );
}

/// Each file of a score folder by name, gzip files read as text.
fn folder_text(d: &str) -> BTreeMap<String, String> {
    std::fs::read_dir(d)
        .unwrap_or_else(|e| panic!("{d}: {e}"))
        .map(|x| {
            let n = x.unwrap().file_name().to_string_lossy().into_owned();
            let t = read(&format!("{d}/{n}"));
            (n, t)
        })
        .collect()
}

/// og prs format ARGS must write WANT_DIR's files, into WANT_DIR + "_og".
fn og_writes(desc: &str, args: &[&str], want_dir: &str) {
    let out = format!("{want_dir}_og");
    let mut a = args.to_vec();
    a.extend(["--out", &out]);
    assert_eq!(og_format(&a), 0, "og prs format: {desc}: exit code");
    assert_eq!(
        folder_text(&out),
        folder_text(want_dir),
        "og prs format: {desc}"
    );
}

fn read_rows(path: &str) -> (Vec<String>, Vec<Vec<String>>) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut lines = text
        .lines()
        .map(|l| l.split('\t').map(str::to_owned).collect::<Vec<_>>());
    let head = lines.next().unwrap_or_default();
    (head, lines.collect())
}

fn first_row(path: &str) -> HashMap<String, String> {
    let (h, r) = read_rows(path);
    h.into_iter()
        .zip(r.into_iter().next().unwrap_or_default())
        .collect()
}

/// A pgsc_calc --outdir for sampleset 'sample': two scores, ACC passed,
/// PGS000002 was dropped. The reference samples' rows come after the
/// target's, as a reader that forgot to skip them would take their values.
fn pgsc_results(d: &str, adjusted: bool, acc: &str) {
    put(
        &format!("{d}/sample/match/sample_summary.csv"),
        &format!(
            "dataset,accession,ambiguous,is_multiallelic,match_flipped,duplicate_best_match,duplicate_ID,match_IDs,match_status,count,score_pass,match_rate\n\
             sample,{acc},false,false,false,false,false,NA,matched,2,true,0.67\n\
             sample,{acc},false,false,false,false,false,NA,unmatched,1,true,0.67\n\
             sample,PGS000002,false,false,false,false,false,NA,unmatched,5,false,0.0\n"
        ),
    );
    if !adjusted {
        put_gz(
            &format!("{d}/sample/score/aggregated_scores.txt.gz"),
            &format!("sampleset\tFID\tIID\tPGS\tSUM\tDENOM\nsample\tHG\tHG\t{acc}\t62.0\t6\n"),
        );
        return;
    }
    put_gz(
        &format!("{d}/sample/score/sample_pgs.txt.gz"),
        &format!(
            "sampleset\tFID\tIID\tPGS\tSUM\tZ_MostSimilarPop\tZ_norm1\tZ_norm2\tpercentile_MostSimilarPop\n\
             sample\tHG\tHG\t{acc}\t1.25\t1.6\t1.5\t1.4\t94.44\n\
             reference\tR1\tR1\t{acc}\t10.0\t0.1\t0.1\t0.1\t12.3\n"
        ),
    );
    put_gz(
        &format!("{d}/sample/score/sample_popsimilarity.txt.gz"),
        "sampleset\tFID\tIID\tPC1\tPC2\tPC10\tSuperPop\tUnrelated\tRF_P_AFR\tRF_P_EUR\tMostSimilarPop\tMostSimilarPop_LowConfidence\tREFERENCE\n\
         reference\tR1\tR1\t-1\t2\t0\tAFR\tTrue\t0.9\t0.1\tAFR\tFalse\tTrue\n\
         sample\tHG\tHG\t0.012\t-0.034\t0.5\t\t\t0.05\t0.95\tEUR\tFalse\tFalse\n",
    );
}

/// The summary og summary collect writes for SAMPLE in D.
fn summary_of(sample: &str, d: &str) -> Value {
    og_summary::collect::post::collect(&Input {
        sample: sample.into(),
        sample_dir: d.into(),
        declared_sex: None,
        assets: repo("assets"),
        group: Group::All,
        tables: repo("bin"),
    })
    .unwrap_or_else(|e| panic!("collect {d}: {e}"))
}

/// (text report, HTML report) of a summary.
fn reports(v: &Value) -> (String, String) {
    let text = json::dumps(v, &Dump::default()).unwrap();
    let sm = og_render::model::Summary::parse(&text).unwrap();
    (
        og_render::text::report(&sm).unwrap(),
        og_render::html::report(&sm).unwrap(),
    )
}

/// (summary, text report, HTML report) of sample S in D.
fn render(d: &str) -> (Value, String, String) {
    let v = summary_of("S", d);
    let (t, h) = reports(&v);
    (v, t, h)
}

fn at<'a>(v: &'a Value, path: &[&str]) -> Option<&'a Value> {
    path.iter().try_fold(v, |v, k| match v {
        Value::Object(m) => m.get(*k),
        _ => None,
    })
}

fn set(v: &mut Value, path: &[&str], x: Value) {
    let (last, head) = path.split_last().unwrap();
    let mut cur = v;
    for k in head {
        cur = match cur {
            Value::Object(m) => m.get_mut(*k).unwrap_or_else(|| panic!("no {k}")),
            _ => panic!("{k}: not an object"),
        };
    }
    match cur {
        Value::Object(m) => {
            m.insert((*last).to_owned(), x);
        }
        _ => panic!("{last}: not an object"),
    }
}

fn j(text: &str) -> Value {
    json::loads(text).unwrap()
}

fn s(x: &str) -> Value {
    Value::Str(x.to_owned())
}

fn not_assessed(v: &Value) -> Vec<String> {
    match at(v, &["not_assessed"]) {
        Some(Value::Array(a)) => a
            .iter()
            .filter_map(|x| match x {
                Value::Str(s) => Some(s.clone()),
                _ => None,
            })
            .collect(),
        _ => Vec::new(),
    }
}

/// The card of an HTML report titled TITLE, up to the next card title.
fn card<'a>(html: &'a str, title: &str) -> &'a str {
    html.split_once(&format!("<h2>{title}</h2>"))
        .map_or("", |x| x.1.split("<h2>").next().unwrap())
}

#[test]
fn prs_format_writes_the_folder_pgsc_calc_reads() {
    let w = Scratch::new("format");
    let wp = |r: &str| w.p(r);
    put_gz(&wp("in/PGS000001_hmPOS_GRCh38.txt.gz"), HM);
    put(
        &wp("labels.tsv"),
        "# comment\npgs_id\ttrait_reported\nPGS000001\tCoronary artery disease\n",
    );
    formatted(
        &wp("pgs"),
        "PGS000001",
        "Coronary artery disease",
        &ROWS1,
        1,
        "left out",
        "",
        "0\tdropped",
        "",
    );
    og_writes(
        "the header pgsc_calc reads, the rows of the autosomes, chrX left out",
        &[
            "--scores",
            &wp("in"),
            "--labels",
            &wp("labels.tsv"),
            "--alleles",
            &wp("alleles.tsv"),
        ],
        &wp("pgs"),
    );
    let al = read(&wp("alleles.tsv"));
    assert_eq!(
        al.lines().collect::<Vec<_>>(),
        [
            "chr1\t1100\tA",
            "chr1\t1100\tG",
            "chr2\t2200\tC",
            "chr2\t2200\tT",
            "chr3\t3300\tG"
        ],
        "every effect and other allele of the autosomes, chr-prefixed and sorted"
    );
    let mut rx = ROWS1.to_vec();
    rx.push("X\t51000000\tA\tG\t1");
    formatted(
        &wp("pgsx"),
        "PGS000001",
        "Coronary artery disease",
        &rx,
        1,
        "scored",
        "",
        "0\tdropped",
        "",
    );
    og_writes(
        "--keep-x: the chrX row stays, the chrY row still goes, chrx_rows.tsv says scored",
        &[
            "--scores",
            &wp("in"),
            "--labels",
            &wp("labels.tsv"),
            "--alleles",
            &wp("allelesx.tsv"),
            "--keep-x",
        ],
        &wp("pgsx"),
    );
    let alx = read(&wp("allelesx.tsv"));
    let alx: Vec<&str> = alx.lines().collect();
    assert_eq!(
        alx[alx.len() - 2..],
        ["chrX\t51000000\tA", "chrX\t51000000\tG"],
        "--keep-x: chrX alleles last"
    );
    assert!(
        !alx.iter().any(|a| a.starts_with("chrY")),
        "--keep-x: no chrY allele"
    );
    put_gz(
        &wp("bad37/PGS000009.txt.gz"),
        &HM.replace("HmPOS_build=GRCh38", "HmPOS_build=GRCh37"),
    );
    assert_eq!(
        og_format(&[
            "--scores",
            &wp("bad37"),
            "--out",
            &wp("o37"),
            "--alleles",
            &wp("a37")
        ]),
        1,
        "refuses a file harmonised to GRCh37"
    );
    put_gz(
        &wp("badnonadd/PGS000010.txt.gz"),
        "#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\teffect_weight\tis_recessive\n1\t5\tA\t1\tTrue\n",
    );
    assert_eq!(
        og_format(&[
            "--scores",
            &wp("badnonadd"),
            "--out",
            &wp("onon"),
            "--alleles",
            &wp("anon")
        ]),
        1,
        "refuses a non-additive score"
    );
    assert_eq!(
        og_format(&[
            "--scores",
            &wp("in"),
            "--ids",
            "PGS000001,PGS000777",
            "--out",
            &wp("omiss"),
            "--alleles",
            &wp("amiss")
        ]),
        1,
        "an id of the list with no file is an error"
    );

    // 4d. strand-ambiguous rows against a small reference: chr1 is GAAGCCAT
    put(&wp("ref/ref.fa"), ">chr1\nGAAG\nCCAT\n");
    put(&wp("ref/ref.fa.fai"), "chr1\t8\t6\t4\t5\n");
    let amb = "#pgs_id=PGS000003\n#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n\
               1\t2\tA\tG\t1\n1\t3\tA\tT\t1\n1\t4\tC\tG\t1\n1\t5\tT\tC\t1\n";
    put_gz(&wp("ambin/PGS000003.txt.gz"), amb);
    // chr1:7 is A: T/C agrees only on the other strand (A/G)
    put_gz(
        &wp("flipin/PGS000003.txt.gz"),
        &format!("{amb}1\t7\tT\tC\t1\n"),
    );
    let amb_run = |src: &str, tag: &str, extra: &[&str]| -> (i32, Vec<String>, Vec<String>) {
        let out = wp(&format!("amb_{tag}"));
        let alleles = format!("{out}.alleles");
        let src = wp(src);
        let mut a = vec!["--scores", &src, "--out", &out, "--alleles", &alleles];
        a.extend(extra);
        let code = og_format(&a);
        if code != 0 {
            return (code, Vec::new(), Vec::new());
        }
        let body = read(&format!("{out}/PGS000003.txt.gz"))
            .lines()
            .skip(5)
            .map(|l| l.split('\t').nth(1).unwrap_or("").to_owned())
            .collect();
        let notes = read(&format!("{out}/score_notes.tsv"))
            .lines()
            .skip(1)
            .map(str::to_owned)
            .collect();
        (code, body, notes)
    };
    let fa = wp("ref/ref.fa");
    let refargs = ["--keep-ambiguous", "--reference", fa.as_str()];
    assert_eq!(
        amb_run("ambin", "off", &[]),
        (
            0,
            vec!["2".into(), "3".into(), "4".into(), "5".into()],
            vec!["PGS000003\t\t2\tdropped".into()]
        ),
        "4d: default: the 2 ambiguous rows are written and counted as dropped"
    );
    assert_eq!(
        amb_run("ambin", "on", &refargs),
        (
            0,
            vec!["2".into(), "3".into(), "4".into(), "5".into()],
            vec!["PGS000003\t\t2\tkept".into()]
        ),
        "4d: --keep-ambiguous, every other row agrees with the reference: kept"
    );
    assert_eq!(
        amb_run("flipin", "flip", &refargs),
        (
            0,
            vec!["2".into(), "5".into(), "7".into()],
            vec!["PGS000003\t\t2\tdropped, strand check failed".into()]
        ),
        "4d: --keep-ambiguous, a row that agrees only on the other strand: the ambiguous rows are left out"
    );
    assert_eq!(
        amb_run("ambin", "noref", &["--keep-ambiguous"]).0,
        1,
        "4d: --keep-ambiguous without --reference"
    );
    formatted(
        &wp("amb_off"),
        "PGS000003",
        "PGS000003",
        &[
            "1\t2\tA\tG\t1",
            "1\t3\tA\tT\t1",
            "1\t4\tC\tG\t1",
            "1\t5\tT\tC\t1",
        ],
        0,
        "left out",
        "",
        "2\tdropped",
        "",
    );
    og_writes(
        "4d: the folder of the default run",
        &[
            "--scores",
            &wp("ambin"),
            "--alleles",
            &wp("amb_off.alleles2"),
        ],
        &wp("amb_off"),
    );
}

/// The score folder the prs-table checks read, as og prs format writes it,
/// with a second score.
fn score_folder(dir: &str) {
    formatted(
        dir,
        "PGS000001",
        "Coronary artery disease",
        &ROWS1,
        1,
        "left out",
        "",
        "0\tdropped",
        "",
    );
    put_gz(
        &format!("{dir}/PGS000002.txt.gz"),
        &format!(
            "#pgs_id=PGS000002\n#trait_reported=Trait two\n#genome_build=GRCh38\n\
             chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n{}",
            "1\t1\tA\t\t1\n".repeat(5)
        ),
    );
}

#[test]
fn prs_table_reads_pgsc_calc_and_the_reports_show_sex_specific_scores() {
    let w = Scratch::new("table");
    let wp = |r: &str| w.p(r);
    let pgs = wp("pgs");
    score_folder(&pgs);
    let out = wp("S_prs_summary.tsv");

    // 2. no panel
    pgsc_results(&wp("raw"), false, "PGS000001");
    let base = |res: &str, kind: &'static str| -> Vec<String> {
        [
            "--sample",
            "S",
            "--results",
            res,
            "--sampleset",
            "sample",
            "--scores",
            &pgs,
            "--input-kind",
            kind,
        ]
        .iter()
        .map(|x| (*x).to_owned())
        .collect()
    };
    let run = |mut a: Vec<String>, extra: &[&str]| -> i32 {
        a.extend(extra.iter().map(|x| (*x).to_owned()));
        prs_table(&a.iter().map(String::as_str).collect::<Vec<_>>())
    };
    assert_eq!(
        run(base(&wp("raw"), "gvcf"), &["--out", &out]),
        0,
        "prs-table exits 0"
    );
    let (hdr, rows) = read_rows(&out);
    assert_eq!(hdr, PRS_COLUMNS, "the columns, Input last");
    assert_eq!(
        rows[0],
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
            "gvcf"
        ],
        "sum 62 (not 62.0), 2 of 3 matched, no percentile, its chrX row left out"
    );
    assert_eq!(
        rows[1],
        [
            "Trait two",
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
            "gvcf"
        ],
        "a score pgsc_calc dropped has no sum; no chrX entry for a score og prs format did not write"
    );

    // 3. with a panel
    pgsc_results(&wp("adj"), true, "PGS000001");
    let anc = wp("S_ancestry.tsv");
    assert_eq!(
        run(
            base(&wp("adj"), "gvcf"),
            &[
                "--panel",
                "pgsc_1000G_v1",
                "--ancestry-out",
                &anc,
                "--out",
                &out
            ]
        ),
        0,
        "prs-table with a panel exits 0"
    );
    let (_, rows) = read_rows(&out);
    assert_eq!(
        rows[0][2..8],
        ["1.25", "2", "3", "66.7", "94.4", "EUR"],
        "the target's row, not the reference's"
    );
    let kv: Vec<(String, String)> = std::fs::read_to_string(&anc)
        .unwrap()
        .lines()
        .filter_map(|l| {
            l.split_once('\t')
                .map(|(a, b)| (a.to_owned(), b.to_owned()))
        })
        .collect();
    let get = |k: &str| kv.iter().find(|x| x.0 == k).map(|x| x.1.as_str());
    assert_eq!(
        (
            get("population"),
            get("reference_panel"),
            get("probability_EUR"),
            get("PC1"),
            get("PC2")
        ),
        (
            Some("EUR"),
            Some("pgsc_1000G_v1"),
            Some("0.95"),
            Some("0.012"),
            Some("-0.034")
        ),
        "ancestry table: {kv:?}"
    );
    let pcs: Vec<&str> = kv
        .iter()
        .filter(|x| x.0.starts_with("PC"))
        .map(|x| x.0.as_str())
        .collect();
    assert_eq!(pcs, ["PC1", "PC2", "PC10"], "the PCs in order");

    // 4. zero matches
    assert_eq!(
        run(base(&wp("none"), "vcf"), &["--zero-matches", "--out", &out]),
        0
    );
    let (_, rows) = read_rows(&out);
    let got: Vec<Vec<&str>> = rows
        .iter()
        .map(|r| r[2..6].iter().chain(&r[12..]).map(String::as_str).collect())
        .collect();
    assert_eq!(
        got,
        [
            ["NA", "0", "3", "0.0", "vcf"],
            ["NA", "0", "5", "0.0", "vcf"]
        ],
        "--zero-matches"
    );
    assert_eq!(
        run(base(&wp("none"), "vcf"), &["--out", &out]),
        1,
        "no pgsc_calc output, no --zero-matches"
    );

    // 4b. every score under pgsc_calc's minimum overlap
    put(
        &wp("below.log"),
        "ERROR Score PGS000001 fails minimum matching threshold (33.33% variants match)\n\
         ERROR pgscatalog.core.lib.pgsexceptions.ZeroMatchesError: All scores fail to meet match threshold 0.75\n",
    );
    assert_eq!(
        run(
            base(&wp("none"), "gvcf"),
            &["--below-threshold", &wp("below.log"), "--out", &out]
        ),
        0
    );
    let (_, rows) = read_rows(&out);
    let got: Vec<&[String]> = rows.iter().map(|r| &r[2..6]).collect();
    assert_eq!(
        got,
        [["NA", "NA", "3", "33.3"], ["NA", "NA", "5", "NA"]],
        "--below-threshold"
    );

    // 4c. sex-specific scores, from the real score list
    let assets = repo("assets/pgs_scores.tsv");
    put_gz(
        &wp("sexin/PGS000662_hmPOS_GRCh38.txt.gz"),
        &HM.replace("PGS000001", "PGS000662"),
    );
    let sex_rows = |labels: &str, sex: &str, tag: &str, applies: &str, evals: &str| {
        let pgs = wp(&format!("sex_{tag}/pgs"));
        let res = wp(&format!("sex_{tag}/res"));
        formatted(
            &pgs,
            "PGS000662",
            "Prostate cancer",
            &ROWS1,
            1,
            "left out",
            applies,
            "0\tdropped",
            evals,
        );
        og_writes(
            &format!("4c: PGS000662 labelled from {labels}"),
            &[
                "--scores",
                &wp("sexin"),
                "--labels",
                labels,
                "--alleles",
                &wp(&format!("sex_{tag}/alleles.tsv")),
            ],
            &pgs,
        );
        pgsc_results(&res, true, "PGS000662");
        let o = wp(&format!("sex_{tag}/S_prs_summary.tsv"));
        let mut a = vec![
            "--sample",
            "S",
            "--results",
            &res,
            "--sampleset",
            "sample",
            "--scores",
            &pgs,
            "--input-kind",
            "gvcf",
            "--panel",
            "pgsc_1000G_v1",
            "--out",
            &o,
        ];
        if !sex.is_empty() {
            a.extend(["--sex", sex]);
        }
        let code = prs_table(&a);
        (code, first_row(&o), o)
    };
    let (male_applies, male_evals) = ("male", "PPM022982;PPM015482\tEUR;AFR\tHR 1.72;OR 1.86");
    let (code, female, f_out) = sex_rows(&assets, "female", "f", male_applies, male_evals);
    assert_eq!(code, 0, "4c: prs-table exits 0 on PGS000662");
    let g = |m: &HashMap<String, String>, k: &str| m.get(k).cloned().unwrap_or_default();
    assert_eq!(
        [
            g(&female, "Score_SUM"),
            g(&female, "Percentile"),
            g(&female, "Applies_To"),
            g(&female, "Sex_Check")
        ],
        ["1.25", "NA", "male", "not applicable to the declared sex"],
        "4c: declared female, PGS000662 (male): raw sum kept, no percentile"
    );
    let (_, male, _) = sex_rows(&assets, "male", "m", male_applies, male_evals);
    assert_eq!(
        [g(&male, "Percentile"), g(&male, "Sex_Check")],
        ["94.4", "applies"],
        "4c: declared male"
    );
    let (_, nosex, n_out) = sex_rows(&assets, "", "n", male_applies, male_evals);
    assert_eq!(
        [g(&nosex, "Percentile"), g(&nosex, "Sex_Check")],
        ["94.4", "no sex declared"],
        "4c: no sex"
    );
    for (tag, src) in [("f", &f_out), ("n", &n_out)] {
        let dd = wp(&format!("sex_sample_{tag}/S"));
        put(
            &format!("{dd}/prs/S_prs_summary.tsv"),
            &std::fs::read_to_string(src).unwrap(),
        );
        put(
            &format!("{dd}/ancestry/S_ancestry.tsv"),
            "key\tvalue\nsample\tS\nreference_panel\tpgsc_1000G_v1\npopulation\tEUR\n",
        );
        let (_, tx, ht) = render(&dd);
        if tag == "f" {
            assert!(
                tx.contains("PGS000662  not applicable to the declared sex (male-specific score)")
                    && ht.contains(
                        "<td>not applicable to the declared sex (male-specific score)</td>"
                    )
                    && !tx.contains("percentile 94.4")
                    && !ht.contains("<td>94.4 (EUR)</td>"),
                "4c: declared female: both reports print not applicable, no percentile:\n{tx}"
            );
        } else {
            assert!(
                tx.contains("percentile 94.4 (EUR)")
                    && ht.contains("<td>94.4 (EUR)</td>")
                    && tx.contains(
                        "No sex was declared, so these sex-specific scores are shown for any sex: PGS000662 (male)."
                    )
                    && ht.contains("PGS000662 (male)."),
                "4c: no sex: both reports keep the percentile and name the score in a note:\n{tx}"
            );
        }
    }
    // Red control: the score list without its applies_to column.
    let two: String = std::fs::read_to_string(&assets)
        .unwrap()
        .lines()
        .map(|l| {
            if l.starts_with('#') {
                format!("{l}\n")
            } else {
                format!("{}\n", l.split('\t').take(2).collect::<Vec<_>>().join("\t"))
            }
        })
        .collect();
    put(&wp("two_columns.tsv"), &two);
    let (_, ctl, _) = sex_rows(&wp("two_columns.tsv"), "female", "ctl", "", "");
    assert_eq!(
        [g(&ctl, "Percentile"), g(&ctl, "Sex_Check")],
        ["94.4", ""],
        "4c: control: without applies_to the female sample would get the percentile"
    );

    // 4d. the Ambiguous column
    formatted(
        &wp("amb_off"),
        "PGS000003",
        "PGS000003",
        &[
            "1\t2\tA\tG\t1",
            "1\t3\tA\tT\t1",
            "1\t4\tC\tG\t1",
            "1\t5\tT\tC\t1",
        ],
        0,
        "left out",
        "",
        "2\tdropped",
        "",
    );
    put(
        &wp("amb_res/sample/match/sample_summary.csv"),
        "dataset,accession,ambiguous,is_multiallelic,match_flipped,duplicate_best_match,duplicate_ID,match_IDs,match_status,count,score_pass,match_rate\n\
         sample,PGS000003,false,false,false,false,false,NA,matched,2,true,0.5\n\
         sample,PGS000003,true,false,false,false,false,NA,excluded,2,true,0.5\n",
    );
    put_gz(
        &wp("amb_res/sample/score/aggregated_scores.txt.gz"),
        "sampleset\tFID\tIID\tPGS\tSUM\tDENOM\nsample\tHG\tHG\tPGS000003\t2\t2\n",
    );
    let (ar, ao, at_) = (wp("amb_res"), wp("amb_off"), wp("amb.tsv"));
    let code = prs_table(&[
        "--sample",
        "S",
        "--results",
        &ar,
        "--sampleset",
        "sample",
        "--scores",
        &ao,
        "--input-kind",
        "gvcf",
        "--out",
        &at_,
    ]);
    assert_eq!(
        (code, g(&first_row(&at_), "Ambiguous")),
        (0, "2 dropped".to_owned()),
        "4d: the Ambiguous column"
    );

    // 4e. a chrX row and the panel's intersection: pgsc_calc's match log
    let mut rx = ROWS1.to_vec();
    rx.push("X\t51000000\tA\tG\t1");
    formatted(
        &wp("pgsx_panel"),
        "PGS000001",
        "Coronary artery disease",
        &rx,
        1,
        "scored",
        "",
        "0\tdropped",
        "",
    );
    pgsc_results(&wp("xpanel"), true, "PGS000001");
    let log =
        "row_nr,accession,chr_name,chr_position,effect_allele,other_allele,match_status,dataset\n";
    let chrx_of = || {
        let (res, sc, o) = (wp("xpanel"), wp("pgsx_panel"), wp("xp.tsv"));
        let code = prs_table(&[
            "--sample",
            "S",
            "--results",
            &res,
            "--sampleset",
            "sample",
            "--scores",
            &sc,
            "--input-kind",
            "gvcf",
            "--panel",
            "pgsc_1000G_v1",
            "--out",
            &o,
        ]);
        assert_eq!(code, 0, "og prs table on the chrX panel");
        g(&first_row(&o), "ChrX")
    };
    put_gz(
        &wp("xpanel/sample/match/sample_log.csv.gz"),
        &format!(
            "{log}0,PGS000001,1,1100,A,G,matched,sample\n1,PGS000001,X,51000000,A,G,excluded,sample\n"
        ),
    );
    assert_eq!(
        chrx_of(),
        "1 left out (not matched in the panel run)",
        "4e: the chrX row excluded"
    );
    put_gz(
        &wp("xpanel/sample/match/sample_log.csv.gz"),
        &format!(
            "{log}0,PGS000001,1,1100,A,G,matched,sample\n1,PGS000001,X,51000000,A,G,matched,sample\n"
        ),
    );
    assert_eq!(chrx_of(), "1 scored", "4e: the chrX row matched");
    put(
        &wp("pgsx_panel/chrx_rows.tsv"),
        "pgs_id\tchrx_rows\tchrx\nPGS000001\t3\tscored\n",
    );
    assert_eq!(
        chrx_of(),
        "1 scored, 2 left out (not matched in the panel run)",
        "4e: 1 of 3 chrX rows matched"
    );
}

#[test]
fn the_prs_card_and_the_fine_ancestry_section_in_both_reports() {
    let w = Scratch::new("prs-reports");
    let cols = PRS_COLUMNS.join("\t");
    let d = w.p("sample_raw/S");
    put(
        &format!("{d}/prs/S_prs_summary.tsv"),
        &format!(
            "{cols}\nCoronary artery disease\tPGS000001\t62\t2\t3\t66.7\tNA\tNA\t8 left out\tgvcf\n\
             Trait two\tPGS000002\t1\t5\t5\t100.0\tNA\tNA\t3 scored\tgvcf\n"
        ),
    );
    let (summ, txt, html) = render(&d);
    assert!(
        txt.contains("Raw score only") && html.contains("Raw score only"),
        "no panel: the raw-score line"
    );
    assert!(
        txt.contains("Rows on chrX were left out, because the sample's sex was not given: PGS000001 (8).")
            && html.contains("Rows on chrX were left out, because the sample&#x27;s sex was not given: PGS000001 (8).")
            && !txt.contains("PGS000002 (3)"),
        "chrX left out: both reports name the score and its row count, not the scored one:\n{txt}"
    );
    assert!(
        html.contains("<td>raw score only</td>"),
        "no panel: the HTML row says raw score only"
    );
    let na = not_assessed(&summ);
    assert!(
        na.iter().any(|x| x == PRS_NOT_ADJUSTED) && !na.iter().any(|x| x == PRS_WITHHELD),
        "no panel: not-assessed names the missing percentiles: {na:?}"
    );
    let caveat = PRS_ANCESTRY_NOTE;
    assert!(
        txt.contains(caveat) && html.contains(caveat),
        "no panel: the ancestry-transfer caveat in both reports"
    );

    // A panel run whose input lists variants only: the percentiles withheld.
    let d = w.p("sample_withheld/S");
    let withheld = "percentile withheld: the input lists variants only";
    put(
        &format!("{d}/prs/S_prs_summary.tsv"),
        &format!(
            "{}\tPercentile_Note\tInput\nCoronary artery disease\tPGS000001\t62\t2\t3\t66.7\tNA\tEUR\t0\t0\t\t\t{withheld}\tvcf\n",
            PRS_COLUMNS[..12].join("\t")
        ),
    );
    let (summ, txt, html) = render(&d);
    assert!(
        txt.contains("Percentile withheld: the input lists variants only")
            && html.contains("Percentile withheld: the input lists variants only"),
        "withheld: both reports say why:\n{txt}"
    );
    let both = format!("{txt}{html}");
    assert!(
        !both.contains("no ancestry reference panel was installed")
            && !both.contains("--ancestry-panel adds it"),
        "withheld: no report says no panel was installed, nor asks for one:\n{txt}"
    );
    let na = not_assessed(&summ);
    assert!(
        na.iter().any(|x| x == PRS_WITHHELD) && !na.iter().any(|x| x == PRS_NOT_ADJUSTED),
        "withheld: not-assessed gives the true reason: {na:?}"
    );

    let d = w.p("sample_adj/S");
    put(
        &format!("{d}/prs/S_prs_summary.tsv"),
        &format!(
            "{cols}\nCoronary artery disease\tPGS000001\t1.25\t2\t3\t66.7\t94.4\tEUR\t0\tgvcf\n"
        ),
    );
    put(
        &format!("{d}/ancestry/S_ancestry.tsv"),
        "key\tvalue\nsample\tS\nreference_panel\tpgsc_1000G_v1\npopulation\tEUR\n",
    );
    let (summ, txt, html) = render(&d);
    assert!(
        txt.contains("percentile 94.4 (EUR)"),
        "panel: the percentile with its group in the text:\n{txt}"
    );
    assert!(
        html.contains("<td>94.4 (EUR)</td>"),
        "panel: the percentile with its group in the HTML"
    );
    assert!(
        txt.contains("among the EUR samples of the pgsc_1000G_v1 reference panel")
            && !txt.contains("Raw score only")
            && !html.contains("Raw score only")
            && !txt.contains("Rows on chrX"),
        "panel: the note names the group and the panel, no raw-score line, no chrX line:\n{txt}"
    );
    assert!(
        !not_assessed(&summ).iter().any(|x| x == PRS_NOT_ADJUSTED),
        "panel: no not-assessed line"
    );
    assert!(
        txt.contains(caveat) && html.contains(caveat),
        "panel: the ancestry-transfer caveat in both reports"
    );

    // 5b. the fine-grained ancestry section (step 26b), beside the PRS card
    put(
        &format!("{d}/ancestry_fine/S_ancestry_fine.tsv"),
        "level\tcode\tgroup\tfraction\tse\tlow\thigh\n\
         continent\tWEA\tWest Eurasia (Europe, Middle East, North Africa)\t0.8200\t0.0300\t0.7600\t0.8800\n\
         continent\tCSA\tCentral and South Asia\t0.1500\t0.0300\t0.0900\t0.2100\n\
         continent\tAFR\tAfrica\t0.0300\t0.0100\t0.0100\t0.0500\n\
         region\tLEV\tLevant\t0.4600\t0.1000\t0.2600\t0.6600\n\
         region\tSAR\tSardinia\t0.0050\t0.0200\t0.0000\t0.0600\n\
         region\tOCE\tOceania\t0.0010\t0.0005\t0.0000\t0.0020\n\
         painting\tWEA\tWest Eurasia (Europe, Middle East, North Africa)\t0.6400\t.\t.\t.\n\
         painting\tAFR\tAfrica\t0.3600\t.\t.\t.\n",
    );
    let mut png = b"\x89PNG\r\n\x1a\n".to_vec();
    png.extend([b'0'; 64]);
    std::fs::write(format!("{d}/ancestry_fine/S_painting.png"), png).unwrap();
    let (_, txt, html) = render(&d);
    assert!(
        txt.contains("percentile 94.4 (EUR)") && html.contains("<td>94.4 (EUR)</td>"),
        "26b: both reports still render with the PRS card beside this one"
    );
    assert!(
        txt.contains("## Fine-Grained Ancestry\n")
            && txt
                .contains("West Eurasia (Europe, Middle East, North Africa)     82% (76% to 88%)"),
        "26b: the text section has the group, its fraction and its bootstrap range:\n{txt}"
    );
    let fine = txt
        .split("## Fine-Grained")
        .nth(1)
        .unwrap_or("")
        .split("##")
        .next()
        .unwrap();
    assert!(
        txt.contains("Sardinia") && !fine.contains("Oceania"),
        "26b: Sardinia shown, Oceania not:\n{fine}"
    );
    assert!(
        html.contains("<h2>Fine-Grained Ancestry</h2>")
            && html.contains("West Eurasia (Europe, Middle East, North Africa) 64%")
            && html.contains("src=\"data:image/png;base64,"),
        "26b: the painted shares and the picture are in the HTML card"
    );
    assert!(
        txt.contains("not a nationality") && html.contains("not a nationality"),
        "26b: both reports say a group name is not a nationality"
    );
}

const VCF_HEAD: &str =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

#[test]
fn the_clinvar_card_never_shows_green_and_carries_the_cdna_check() {
    let w = Scratch::new("clinvar");
    let d = w.p("clinvar_hom/S");
    put(
        &format!("{d}/clinvar/S_clinvar_hits.vcf"),
        &format!(
            "{VCF_HEAD}chr1\t20\t2\tT\tG\t50\tPASS\tGENEINFO=GENEB:22;CLNSIG=Pathogenic;\
             CLNREVSTAT=reviewed_by_expert_panel\tGT\t1/1\n"
        ),
    );
    let (_, txt, html) = render(&d);
    let badge_line = html
        .lines()
        .find(|l| l.contains("ClinVar matches"))
        .unwrap_or("");
    assert!(
        !badge_line.contains("badge-green") && badge_line.contains("badge-yellow\">1<"),
        "a homozygous 3-star P/LP hit gets no green badge: {badge_line}"
    );
    let freq = "a heterozygous hit in a dominant gene as well as a homozygous one: a hit above 1% in gnomAD";
    assert!(
        txt.contains(freq)
            && html.contains(freq)
            && !format!("{txt}{html}").contains("Look up a homozygous P/LP hit"),
        "both reports ask for every hit's population frequency:\n{txt}"
    );
    assert!(
        txt.contains("frequency: gnomAD AF unknown") && html.contains("<td>gnomAD AF unknown</td>"),
        "without step 13 the frequency is unknown:\n{txt}"
    );
    put(&format!("{d}/clinvar/S_clinvar_hits.vcf"), VCF_HEAD);
    let (_, _, html0) = render(&d);
    assert!(
        html0.contains(
            "ClinVar matches</span><span class=\"value\"><span class=\"badge badge-gray\">0<"
        ),
        "no hit gets a gray badge, not a green all-clear"
    );

    // 6b. step 06's cDNA read check beside each hit
    let d = w.p("clinvar_cdna/S");
    put(
        &format!("{d}/clinvar/S_clinvar_hits.vcf"),
        &format!(
            "{VCF_HEAD}chr1\t20\t2\tT\tG\t50\tPASS\tGENEINFO=GENEB:22;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel\tGT\t0/1\n\
             chr2\t30\t3\tC\tA\t50\tPASS\tGENEINFO=GENEC:33;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel\tGT\t0/1\n"
        ),
    );
    let (summ, _, html) = render(&d);
    assert!(
        !html.contains("cDNA read check")
            && at(&summ, &["sections", "clinvar", "data", "cdna_flagged"]).is_none(),
        "no table, no column (a Nextflow folder)"
    );
    let flag = "possible cDNA/transgene contamination: check the reads";
    put(
        &format!("{d}/clinvar/S_clinvar_cdna_check.tsv"),
        &format!(
            "chrom\tpos\tref\talt\talt_reads\talt_clipped\tref_reads\tref_clipped\tcdna_check\tcdna_note\n\
             chr1\t20\tT\tG\t9\t8\t11\t0\tflagged\t{flag}\nchr2\t30\tC\tA\t12\t0\t10\t0\tclear\t\n"
        ),
    );
    let (summ, txt, html) = render(&d);
    let shown = format!("{flag} (8 of 9 alt reads and 0 of 11 ref reads clipped at an exon edge)");
    assert!(
        txt.contains("chr1:20 T>G")
            && txt.contains(&format!("reads: {shown}"))
            && txt.contains("reads: clear"),
        "the text report prints the flag under the flagged hit:\n{txt}"
    );
    assert!(
        html.contains("<th>cDNA read check</th>")
            && html.contains(&format!("<td>{shown}</td>"))
            && html.contains("<td>clear</td>"),
        "the HTML hits table has the column, the flag and clear"
    );
    assert!(
        html.contains("Flagged by the cDNA read check</span><span class=\"value\"><span class=\"badge badge-red\">1<"),
        "the card counts one flagged hit"
    );
    let p = w.p("summary.json");
    put(&p, &json::dumps(&summ, &Dump::default()).unwrap());
    assert_eq!(
        rc(&[
            "dev",
            "schema",
            &repo("tests/schema/summary.schema.json"),
            &p
        ]),
        0,
        "the summary validates against the schema"
    );
}

#[test]
fn the_telomere_mosaic_paralog_and_export_cards() {
    let w = Scratch::new("cards");
    // 7. telomere card
    let d = w.p("telomere/S");
    put(
        &format!("{d}/telomere/S/S_summary.tsv"),
        "PID\tsample\tintratelomeric_reads\ttotal_reads_with_tel_gc\ttel_content\nS\ttumor\t61873\t185204533\t334.08\n",
    );
    let (_, txt, html) = render(&d);
    assert!(
        txt.contains("## Telomere content (relative, TelomereHunter)")
            && txt.contains(
                "Telomere content: 334.08 intratelomeric reads per million reads with 48-52% GC"
            ),
        "telomere: the text heading and the value with its unit:\n{txt}"
    );
    assert!(
        html.contains("<h2>Telomere content (relative)</h2>")
            && html.contains("intratelomeric reads per million reads with 48-52% GC")
            && html.contains("not a telomere length"),
        "telomere: the HTML card title and the unit"
    );
    let mut said = Vec::new();
    for dir in [
        repo("bin"),
        repo("rust/crates/og-render/src"),
        repo("rust/crates/og-render/src/sections"),
    ] {
        for f in std::fs::read_dir(&dir).unwrap() {
            let p = f.unwrap().path();
            if p.is_file()
                && std::fs::read_to_string(&p)
                    .unwrap_or_default()
                    .contains("Telomere Length")
            {
                said.push(p);
            }
        }
    }
    assert!(
        said.is_empty(),
        "no report source says 'Telomere Length': {said:?}"
    );

    // 8. mosaic loss of chrY or chrX
    let d = w.p("mosaic/S");
    let mosaic_table = |rows: &[(&str, &str)]| {
        let body: String = rows
            .iter()
            .map(|(k, v)| format!("mosaic_{k}\t{v}\n"))
            .collect();
        put(
            &format!("{d}/qc/S_sample_qc.tsv"),
            &format!("key\tvalue\nsample\tS\ninferred_sex\tmale\n{body}"),
        );
        let (_, t, h) = render(&d);
        (t, h)
    };
    let male = [
        ("sex", "male"),
        ("chrom", "chrY"),
        ("depth_ratio", "0.351"),
        ("depth_ratio_expected", "0.5"),
        ("loss_depth", "0.298"),
        ("loss_depth_low", "0.293"),
        ("loss_depth_high", "0.303"),
        ("floor", "0.100"),
        ("call", "detected"),
        ("reason", "the low end of the 95% range is above 10%"),
    ];
    let with = |over: &[(&'static str, &'static str)]| -> Vec<(&'static str, &'static str)> {
        let mut v: Vec<(&str, &str)> = male.to_vec();
        for (k, x) in over {
            match v.iter_mut().find(|y| y.0 == *k) {
                Some(y) => y.1 = x,
                None => v.push((k, x)),
            }
        }
        v
    };
    let (txt, html) = mosaic_table(&male);
    let line = "From depth: estimated share of the cells that lost the Y 29.8% (95% range 29.3% to 30.3%); \
                chrY to autosome depth ratio 0.351";
    let card_m = card(&html, "Mosaic Loss of Chromosome Y or X");
    assert!(
        txt.contains("## Mosaic loss of chromosome Y or X")
            && txt.contains("chrY: detected (the low end")
            && txt.contains(line),
        "mosaic: the text report heading, call and estimate:\n{txt}"
    );
    assert!(
        card_m.contains("badge-yellow\">detected<") && card_m.contains(line),
        "mosaic: the HTML card:\n{card_m}"
    );
    for t in [txt.as_str(), card_m] {
        assert!(
            t.contains("population-level association, not a diagnosis")
                && t.contains("Below about 10% of the cells"),
            "mosaic: population-level, not a diagnosis, and the 10% floor"
        );
    }
    assert!(
        !txt.contains("allele fractions at"),
        "mosaic: a man has no allele-fraction line"
    );
    let (txt, html) = mosaic_table(&with(&[
        ("sex", "female"),
        ("chrom", "chrX"),
        ("depth_ratio", "0.849"),
        ("depth_ratio_expected", "1.0"),
        ("x_het_sites", "600"),
        ("loss_af", "0.316"),
        ("loss_af_low", "0.288"),
        ("loss_af_high", "0.340"),
    ]));
    assert!(
        txt.contains("of the cells that lost one X")
            && txt.contains("From chrX allele fractions at 600 heterozygous sites: 31.6% (95% range 28.8% to 34.0%)")
            && html.contains("31.6%"),
        "mosaic: a woman shows one X and the allele-fraction line:\n{txt}"
    );
    let (txt, html) = mosaic_table(&with(&[
        ("depth_ratio", ""),
        ("loss_depth", "0.190"),
        ("loss_depth_low", "-0.002"),
        ("loss_depth_high", "0.382"),
        ("call", "not_detected"),
    ]));
    let card_m = card(&html, "Mosaic Loss of Chromosome Y or X");
    assert!(
        txt.contains("chrY: not detected above 10%")
            && !txt.contains("below 10%")
            && txt.contains("depth ratio ?,")
            && card_m.contains("badge-gray\">not detected above 10%<"),
        "mosaic: not detected says 'above 10%', a missing depth ratio prints '?':\n{txt}"
    );
    put(
        &format!("{d}/qc/S_sample_qc.tsv"),
        "key\tvalue\nsample\tS\ninferred_sex\tmale\n",
    );
    let (_, txt, html) = render(&d);
    assert!(
        card(&html, "Mosaic Loss of Chromosome Y or X").contains("badge-gray\">Not run<")
            && !txt.contains("## Mosaic loss"),
        "mosaic: no estimate, the card says Not run and the text has no section"
    );

    // 9. paralog card (step 35b)
    assert!(
        !html.contains("Paralog Genes (Paraphase)"),
        "paralogs: no card when step 35b did not run"
    );
    let d = w.p("paralogs/S");
    let f = format!("{d}/paralogs/S_paraphase_copy_number.tsv");
    put(
        &f,
        "region\tgene\tcopy_number\tstatus\tdetail\ttotal_cn\tregion_depth\tgenome_depth\npms2\tPMS2\t2\tcalled\t.\t4\t60\t30\n",
    );
    std::fs::File::options()
        .write(true)
        .open(&f)
        .unwrap()
        .set_modified(UNIX_EPOCH + Duration::from_secs(1))
        .unwrap();
    put(
        &format!("{d}/logs/run_status.tsv"),
        "meta\tstarted_epoch\t1000\nstep\t35b\tfailed\n",
    );
    let (summ, _, html) = render(&d);
    let card_p = card(&html, "Paralog Genes (Paraphase)");
    assert!(
        at(&summ, &["sections", "paralogs", "state"]) == Some(&s("stale"))
            && card_p.contains("badge-yellow\">Stale<")
            && card_p.contains("<td>PMS2</td><td>2</td>"),
        "paralogs: a result older than a run where step 35b failed is stale, with its PMS2 row:\n{card_p}"
    );

    // 12. raw-data export card
    let d = w.p("consumer/S");
    std::fs::create_dir_all(&d).unwrap();
    let (_, txt, html) = render(&d);
    assert!(
        !txt.contains("Raw-Data Export") && !html.contains("Raw-Data Export"),
        "raw-data export: no card and no text section without the step's table"
    );
    put(
        &format!("{d}/consumer_export/S_consumer_export.tsv"),
        "key\tvalue\nsites\t200\ncalled\t150\nno_call\t50\nfrom_reference_blocks\t120\nfile_23andme\tS_23andme_v5.zip\n\
         file_ancestrydna\tS_ancestrydna_layout.zip\nno_call: no gVCF record\t30\nno_call: filtered call\t20\n",
    );
    let (_, txt, html) = render(&d);
    let want = [
        "Sites with a genotype: 150 of 200 (75.0%); '--' at 50",
        "Homozygous reference from the gVCF's reference blocks: 120",
        "consumer_export/S_23andme_v5.zip and consumer_export/S_ancestrydna_layout.zip",
        "'--' because: no gVCF record 30; filtered call 20",
        "Uploading one to a service shares your genotypes",
    ];
    assert!(
        txt.contains("## Raw-Data Export (23andMe v5 layout)")
            && want.iter().all(|x| txt.contains(x)),
        "raw-data export: the text section has every line:\n{txt}"
    );
    assert!(
        html.contains("<h2>Raw-Data Export (23andMe v5 layout)</h2>")
            && want.iter().all(|x| html.contains(&*e(x))),
        "raw-data export: the HTML card has every line"
    );
}

#[test]
fn the_blood_group_and_couple_cards() {
    let w = Scratch::new("blood-couple");
    // 10. blood groups
    let d = w.p("blood/S");
    put(
        &format!("{d}/bloodtype/S_bloodtype.tsv"),
        "system\tphenotype\tgenotype\tconfidence\tnote\nABO\tO\tO.01/O.01\thigh\t\n\
         RhD\tD+\tRHD/RHD (RHD copies by read depth 1.78)\thigh\t\nDuffy\tFy(a+b-)\tFY*A/FY*B GATA\thigh\tc.-67T>C\n\
         HPA-1\tHPA-1(a+b+)\t1a/1b\thigh\t\n",
    );
    let (_, txt, html) = render(&d);
    assert!(
        txt.contains("## Blood Groups (predicted from the genome)")
            && [
                "ABO       O ",
                "RhD       D+ ",
                "Duffy     Fy(a+b-) ",
                "HPA-1     HPA-1(a+b+) "
            ]
            .iter()
            .all(|x| txt.contains(x)),
        "blood groups: the text section has every system:\n{txt}"
    );
    assert!(
        html.contains("<h2>Blood Groups (predicted)</h2>")
            && html.contains("<td>Fy(a+b-)</td>")
            && html.contains("<td>FY*A/FY*B GATA</td>")
            && html.contains(">D+</span>"),
        "blood groups: the HTML card has ABO, RhD and the table"
    );
    assert_eq!(
        (
            txt.matches("not typed from blood").count(),
            html.matches("not typed from blood").count()
        ),
        (1, 1),
        "blood groups: both reports say it is predicted, not typed from blood"
    );
    let empty = w.p("noblood/S");
    std::fs::create_dir_all(&empty).unwrap();
    assert!(
        !render(&empty).2.contains("Blood Groups"),
        "blood groups: no card without step 41"
    );
    let words = |t: &str| {
        let l = t.to_lowercase();
        ["transfus", "donat", "donor"].iter().any(|x| l.contains(x))
    };
    let src = repo("rust/crates/og-bloodtype/src");
    let mut names: Vec<PathBuf> = std::fs::read_dir(&src)
        .unwrap()
        .map(|f| f.unwrap().path())
        .collect();
    names.sort();
    let code: String = names
        .iter()
        .filter(|p| p.extension().is_some_and(|x| x == "rs"))
        .map(|p| std::fs::read_to_string(p).unwrap())
        .collect();
    let said: Vec<&str> = [
        ("text report", &txt),
        ("HTML report", &html),
        ("og bloodtype", &code),
    ]
    .iter()
    .filter(|(_, t)| words(t))
    .map(|(n, _)| *n)
    .collect();
    assert!(
        said.is_empty(),
        "blood groups: nothing about transfusion or donation: {said:?}"
    );
    assert!(
        words("for a Transfusion"),
        "control: the wording check finds the word"
    );

    // 13. couple card (step 39)
    let genes = std::fs::read_to_string(repo("assets/carrier_genes.tsv")).unwrap();
    let g: Vec<&str> = genes
        .lines()
        .find(|x| x.starts_with("NAGA\t"))
        .unwrap()
        .split('\t')
        .collect();
    let pos = (g[2].parse::<i64>().unwrap() + g[3].parse::<i64>().unwrap()) / 2;
    let cv = |p: i64| {
        format!(
            "{VCF_HEAD}chr22\t{p}\t1\tA\tG\t.\tPASS\tGENEINFO=NAGA:4668;CLNSIG=Pathogenic;\
             CLNREVSTAT=reviewed_by_expert_panel\tGT\t0/1\n"
        )
    };
    let c = w.p("couple");
    put(
        &format!("{c}/ped.tsv"),
        "sample\tsex\tfather\tmother\tpartner\nP1\tfemale\t.\t.\tP2\nP2\tmale\t.\t.\tP1\n",
    );
    put(&format!("{c}/P1.vcf"), &cv(pos));
    put(&format!("{c}/P2.vcf"), &cv(pos + 5));
    let out = og(&[
        "couple",
        "screen",
        "--pedigree",
        &format!("{c}/ped.tsv"),
        "--genes",
        &repo("assets/carrier_genes.tsv"),
        "--clinvar",
        &format!("P1={c}/P1.vcf"),
        "--clinvar",
        &format!("P2={c}/P2.vcf"),
        "--out-dir",
        &format!("{c}/out"),
    ]);
    assert!(
        out.status.success(),
        "og couple screen: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let d = format!("{c}/P2");
    let table = format!("{d}/couple/P2_couple_carrier.tsv");
    put(
        &table,
        &std::fs::read_to_string(format!("{c}/out/P2_couple_carrier.tsv")).unwrap(),
    );
    let (txt, html) = reports(&summary_of("P2", &d));
    let card_c = card(&html, "Couple Carrier Screen");
    let want = [
        "Partner: P1".to_owned(),
        "Shared: NAGA (Schindler disease): 1 in 4 per pregnancy".to_owned(),
        format!("This sample: chr22:{} A&gt;G het", pos + 5),
        format!("P1: chr22:{pos} A&gt;G het, Pathogenic, 3 stars"),
        "genetic counsellor".to_owned(),
    ];
    assert!(
        want.iter().all(|x| card_c.contains(x.as_str())),
        "couple: the card:\n{card_c}"
    );
    assert!(
        want.iter().all(|x| txt.contains(&x.replace("&gt;", ">")))
            && txt.contains("## Couple Carrier Screen (step 39)"),
        "couple: the text report says the same:\n{txt}"
    );
    std::fs::remove_file(&table).unwrap();
    let (txt, html) = reports(&summary_of("P2", &d));
    assert!(
        !html.contains("Couple Carrier Screen") && !txt.contains("Couple Carrier Screen"),
        "couple: without the table there is no card and no text section"
    );
}

/// The card titles of DEMO-001's HTML report: the chapters move cards, they
/// never add, drop or rename one.
const GOLDEN_TITLES: [&str; 48] = [
    "Quality Control",
    "Mosaic Loss of Chromosome Y or X",
    "Variant Calling",
    "ClinVar Screening",
    "Conflicting classifications with a P/LP submission: review, not counted (1 of 1)",
    "Carrier List (this sample)",
    "Couple Carrier Screen",
    "Pharmacogenomics",
    "pypgx Pharmacogenomics",
    "PGx Copy Number",
    "CPIC Drug Recommendations",
    "Genes With a Non-Normal or Unresolved Phenotype (CPIC)",
    "CYP2D6 Across Callers",
    "HLA Typing",
    "Structural Variants",
    "Structural Variants, Class 4 and 5 (AnnotSV)",
    "Cancer Predisposition",
    "Repeat Expansions",
    "Runs of Homozygosity",
    "Mitochondrial Haplogroup",
    "Paralog Genes (Paraphase)",
    "SMN1/SMN2 Copy Number (Parascopy)",
    "Y-Chromosome Haplogroup",
    "Y Placement on the Current Public Tree",
    "Y-STR Markers (HipSTR)",
    "Raw-Data Export (23andMe v5 layout)",
    "Blood Groups (predicted)",
    "Archaic Ancestry (Neanderthal, Denisovan)",
    "Lp(a): LPA KIV-2 Repeat Count",
    "MUC1 VNTR Frameshifts (VNtyper)",
    "Locityper Second Opinion at Hard Loci",
    "Viral DNA in the Unmapped Reads",
    "Clonal Haematopoiesis (CHIP)",
    "Mosaic Chromosomal Alterations (MoChA)",
    "Mobile Element Insertions",
    "Telomere content (relative)",
    "Mitochondrial Analysis",
    "Clinical Variant Filter",
    "Variant Prioritization (slivar)",
    "Symptom-Driven Ranking (Exomiser)",
    "Polygenic Risk Scores",
    "Fine-Grained Ancestry",
    "Traits",
    "Ancient Populations",
    "Genes This Genome Could Not Read Reliably",
    "ClinVar Hits (best-reviewed first, 3 of 3)",
    "Steps Not Run",
    "Not Assessed by This Pipeline",
];

/// Every card title of an HTML page, sorted: each `<h2>`, whatever its
/// attributes, except the chapter and contents headings.
fn titles(html: &str) -> Vec<String> {
    let mut t: Vec<String> = html
        .split("<h2")
        .skip(1)
        .filter_map(|x| x.split_once('>'))
        .filter(|(attrs, _)| !attrs.contains("chapter-title") && !attrs.contains("toc-title"))
        .filter_map(|(_, x)| x.split_once("</h2>").map(|y| y.0))
        .filter(|x| !x.contains('<'))
        .map(str::to_owned)
        .collect();
    t.sort();
    t
}

fn golden() -> Vec<String> {
    let mut g: Vec<String> = GOLDEN_TITLES.iter().map(|x| (*x).to_owned()).collect();
    g.sort();
    g
}

/// Every http(s) src or href that is not the loopback viewer and not inside
/// a link marked as opening a website.
fn external_problems(html: &str) -> Vec<String> {
    const OPEN: &str = "<span class=\"offsite\"><a href=\"";
    let mut kept = String::new();
    let mut rest = html;
    while let Some(i) = rest.find(OPEN) {
        kept.push_str(&rest[..i]);
        let tail = &rest[i..];
        match tail.find("</span></span>") {
            Some(end)
                if tail[..end].contains("</a> <span class=\"offsite-note\">(opens a website") =>
            {
                rest = &tail[end + "</span></span>".len()..];
            }
            _ => {
                kept.push_str(OPEN);
                rest = &tail[OPEN.len()..];
            }
        }
    }
    kept.push_str(rest);
    // ASCII only, so a byte offset in `low` is the same offset in `kept`.
    let low = kept.to_ascii_lowercase();
    let mut found = Vec::new();
    for attr in ["src", "href"] {
        let mut from = 0;
        while let Some(i) = low[from..].find(attr) {
            let mut k = from + i + attr.len();
            from = k;
            let b = low.as_bytes();
            while k < b.len() && b[k].is_ascii_whitespace() {
                k += 1;
            }
            if b.get(k) != Some(&b'=') {
                continue;
            }
            k += 1;
            while k < b.len() && b[k].is_ascii_whitespace() {
                k += 1;
            }
            if matches!(b.get(k), Some(b'"' | b'\'')) {
                k += 1;
            }
            if low[k..].starts_with("http:") || low[k..].starts_with("https:") {
                let end = kept[k..]
                    .find(['"', '\'', ' ', '>'])
                    .map_or(kept.len(), |x| k + x);
                let url = &kept[k..end];
                if !url.starts_with(og_render::html::VIEWER_URL) {
                    found.push(url.to_owned());
                }
            }
        }
    }
    found
}

fn frame_of<'a>(html: &'a str, sid: &str) -> &'a str {
    html.split_once(&format!("id=\"{sid}\">"))
        .map_or("", |x| x.1.split("</section>").next().unwrap())
}

/// The text report's chapter titles, in order, each with the text under it.
fn text_chapters(txt: &str) -> Vec<(String, String)> {
    let mut out: Vec<(String, String)> = Vec::new();
    for l in txt.split_inclusive('\n') {
        let head = l
            .strip_prefix("#  ")
            .and_then(|r| {
                r.split_once(". ")
                    .filter(|(n, _)| !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()))
            })
            .and_then(|(_, t)| t.split_once(": ").map(|x| x.0.to_owned()));
        match head {
            Some(t) => out.push((t, String::new())),
            None => {
                if let Some(last) = out.last_mut() {
                    last.1.push_str(l);
                }
            }
        }
    }
    out
}

fn under<'a>(txt: &'a str, title: &str) -> &'a str {
    txt.find(title).map_or("", |i| {
        let end = (i + 120).min(txt.len());
        let end = (end..=txt.len())
            .find(|&k| txt.is_char_boundary(k))
            .unwrap();
        &txt[i..end]
    })
}

#[test]
fn the_chapters_on_the_demo_sample() {
    let w = Scratch::new("chapters");
    let g = w.p("demo");
    std::fs::create_dir_all(&g).unwrap();
    let out = og(&[
        "dev",
        "demo-sample",
        "--genome-dir",
        &g,
        "--records",
        "20000",
        "--repo",
        REPO,
    ]);
    assert!(
        out.status.success(),
        "og dev demo-sample: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let sd = format!("{g}/DEMO-001");
    let summ = summary_of("DEMO-001", &sd);
    let (txt, html) = reports(&summ);

    // The haplogroup links: FamilyTreeDNA and YFull say the address carries
    // the haplogroup, the docs link does not.
    for (kind, hg) in [("y", "R-M269*(xR-L23)"), ("mt", "H1a+16189")] {
        let hd = w.p(&format!("hg-{kind}/S"));
        let title = if kind == "y" {
            put(
                &format!("{hd}/y_haplogroup/S_y_haplogroup.txt"),
                &format!("Sample_name\tHg\tValid_markers\tQC-score\nS\t{hg}\t9\t1\n"),
            );
            "Y-Chromosome Haplogroup"
        } else {
            put(
                &format!("{hd}/mito/S_haplogroup.txt"),
                &format!("SampleID\tHaplogroup\nS\t{hg}\n"),
            );
            "Mitochondrial Haplogroup"
        };
        let (_, _, h) = render(&hd);
        let c = card(&h, title);
        let marks: Vec<&str> = c
            .split("<span class=\"offsite-note\">(")
            .skip(1)
            .map(|x| x.split(")</span>").next().unwrap())
            .collect();
        let mut want =
            vec![og_render::sections::ancestry::OFFSITE_HAPLOGROUP; marks.len().saturating_sub(1)];
        want.push(og_render::html::OFFSITE);
        assert!(
            marks.len() >= 3 && marks == want,
            "{hg}: the offsite notes: {marks:?}"
        );
    }

    assert_eq!(
        external_problems(&html),
        Vec::<String>::new(),
        "no off-machine address but the marked links"
    );
    for u in [
        "https://discover.familytreedna.com/",
        "https://www.yfull.com/",
        "https://github.com/",
    ] {
        assert!(
            html.contains(&format!("<span class=\"offsite\"><a href=\"{u}")),
            "the marked link {u} is there"
        );
    }
    for bad in [
        "<a href=\"https://example.org/x\">x</a>",
        "<img src=\"https://example.org/i.png\">",
        "<a class=\"offsite\" href=\"https://example.org/y\">y</a>",
    ] {
        assert!(
            !external_problems(&format!("{html}{bad}")).is_empty(),
            "red control: {bad} is caught"
        );
    }

    let sids: Vec<&str> = og_render::chapters::CHAPTERS
        .iter()
        .flat_map(|c| c.sections.iter().map(|x| x.id))
        .collect();
    let cv = frame_of(&html, "clinvar");
    assert!(
        cv.contains("<h2>ClinVar Screening</h2>")
            && cv.contains("<h2>ClinVar Hits")
            && !cv.contains("badge-green"),
        "the ClinVar section has its two cards and no green badge"
    );
    assert!(
        frame_of(&html, "qc").contains("badge-green"),
        "red control: the detector sees green in the QC block"
    );
    let mut nohit = summ.clone();
    set(
        &mut nohit,
        &["sections", "clinvar", "data", "count"],
        Value::Int(0),
    );
    set(
        &mut nohit,
        &["sections", "clinvar", "data", "hits"],
        j("[]"),
    );
    let nohit_html = reports(&nohit).1;
    let nv = frame_of(&nohit_html, "clinvar");
    assert!(
        nv.contains("<h2>ClinVar Screening</h2>")
            && nv.contains("<span class=\"badge badge-gray\">0</span>")
            && !nv.contains("badge-green"),
        "no ClinVar hit: the screen card with a gray 0, still no green:\n{nv}"
    );
    assert!(
        frame_of(&html, "cpic").contains("CPIC Drug Recommendations"),
        "the CPIC table sits in the CPIC section"
    );

    let heads = [
        ("QUALITY", "## Quality control"),
        ("CLINICAL", "## ClinVar Pathogenic Screen"),
        ("CARRIER", "## Carrier List (this sample)"),
        ("MEDICINES", "## Pharmacogenomics (PharmCAT)"),
        ("ANCESTRY", "## Mitochondrial Haplogroup"),
        ("TRAITS AND SCORES", "## Polygenic Risk Scores"),
        (
            "OTHER AND EXPORTS",
            "## Raw-Data Export (23andMe v5 layout)",
        ),
    ];
    let chunks = text_chapters(&txt);
    for (c, h) in heads {
        let body = chunks
            .iter()
            .find(|x| x.0 == c)
            .map_or("", |x| x.1.as_str());
        assert!(body.contains(h), "text chapter {c} holds {h}");
    }

    assert_eq!(
        titles(&html),
        golden(),
        "the card titles are the golden list"
    );
    assert_ne!(
        titles(&html.replace("<h2>Traits</h2>", "")),
        golden(),
        "red control: a dropped card is caught"
    );
    assert_ne!(
        titles(&html.replace(
            "<h2>Traits</h2>",
            "<h2>Traits</h2><h2 class=\"card-title\">Extra</h2>"
        )),
        golden(),
        "red control: an extra card with attributes on its heading is caught"
    );

    // "New since" and the database re-check
    let badge = "New since 2026-01-01T00:00:00Z";
    let mut since = summ.clone();
    set(
        &mut since,
        &["new_since"],
        j(r#"{"previous_utc": "2026-01-01T00:00:00Z", "sections": ["clinvar"]}"#),
    );
    set(
        &mut since,
        &["recheck"],
        j(
            r#"{"checked_utc": "2026-02-02T00:00:00Z", "affected_steps": ["06", "24"], "rerun": null,
              "items": [{"item": "clinvar_file_date", "label": "ClinVar", "class": "refresh-available"},
                        {"item": "VEP_CACHE_RELEASE", "label": "VEP", "class": "no-newer-release"}]}"#,
        ),
    );
    let (t2, h2) = reports(&since);
    let badge_problems = |page: &str, want: &[&str]| -> Vec<String> {
        let got: Vec<&str> = sids
            .iter()
            .copied()
            .filter(|sid| frame_of(page, sid).contains(badge))
            .collect();
        let n = page.matches(badge).count();
        if got == want && n == want.len() {
            Vec::new()
        } else {
            vec![format!("the badge in {got:?}, {n} on the page")]
        }
    };
    assert_ne!(
        badge_problems(&h2.replace(badge, ""), &["clinvar"]),
        Vec::<String>::new(),
        "red control: a missing badge is caught"
    );
    assert_eq!(
        badge_problems(&h2, &["clinvar"]),
        Vec::<String>::new(),
        "the badge only where the result changed"
    );
    assert!(
        !html.contains("New since") && !txt.contains("New since"),
        "no badge without a previous run"
    );
    let mut every = since.clone();
    let keys: Vec<Value> = og_render::chapters::CHAPTERS
        .iter()
        .flat_map(|c| c.sections.iter().flat_map(|x| x.keys.iter().map(|k| s(k))))
        .collect();
    set(&mut every, &["new_since", "sections"], Value::Array(keys));
    assert!(
        !badge_problems(&reports(&every).1, &["clinvar"]).is_empty(),
        "red control: a badge on every section"
    );
    let text_with = |sections: &[&str], missing: Option<&str>| -> (String, String) {
        let mut x = since.clone();
        if let Some(m) = missing {
            set(&mut x, &["sections", m, "state"], s("missing"));
        }
        set(
            &mut x,
            &["new_since", "sections"],
            Value::Array(sections.iter().map(|k| s(k)).collect()),
        );
        let (t, h) = reports(&x);
        (h, t)
    };
    let cyp = |h: &str, t: &str| {
        let _ = h;
        under(t, "## CYP2D6 across callers").contains(&format!("[{badge}]"))
    };
    let (h, t) = text_with(&["pypgx"], None);
    assert!(
        cyp(&h, &t),
        "the text CYP2D6 block carries the badge when pypgx changed:\n{}",
        under(&t, "## CYP2D6")
    );
    let (h, t) = text_with(&["clinvar"], None);
    assert!(
        !cyp(&h, &t),
        "red control: the CYP2D6 check fails when only ClinVar changed"
    );
    let (h, t) = text_with(&["roh"], Some("roh"));
    assert!(
        frame_of(&h, "roh").contains(badge)
            && frame_of(&h, "roh").contains("Not run")
            && !t.contains("## Runs of Homozygosity")
            && !t.contains(badge),
        "a section gone missing keeps its badge on its HTML card, the text leaves the section out"
    );
    assert_eq!(
        titles(&h2),
        golden(),
        "new since: the card titles stay the golden list"
    );

    let hist = format!("{sd}/history/2026-01-01T000000Z");
    std::fs::create_dir_all(&hist).unwrap();
    let mut old = summ.clone();
    set(&mut old, &["generated_utc"], s("2026-01-01T00:00:00Z"));
    let has_hits = matches!(at(&summ, &["sections", "clinvar", "data", "hits"]), Some(Value::Array(a)) if !a.is_empty());
    if let Some(Value::Object(sec)) = match &mut old {
        Value::Object(m) => m.get_mut("sections"),
        _ => None,
    } && let Some(Value::Object(cvs)) = sec.get_mut("clinvar")
        && let Some(Value::Object(data)) = cvs.get_mut("data")
        && let Some(Value::Array(hits)) = data.get_mut("hits")
        && let Some(Value::Object(h0)) = hits.first_mut()
    {
        h0.insert("significance".to_owned(), s("Uncertain significance"));
    }
    put(
        &format!("{hist}/summary.json"),
        &json::dumps(&old, &Dump::default()).unwrap(),
    );
    if std::path::Path::new(&format!("{sd}/run_manifest.tsv")).is_file() {
        std::fs::copy(
            format!("{sd}/run_manifest.tsv"),
            format!("{hist}/run_manifest.tsv"),
        )
        .unwrap();
    }
    let again = summary_of("DEMO-001", &sd);
    std::fs::remove_dir_all(format!("{sd}/history")).unwrap();
    assert!(
        has_hits
            && at(&summ, &["new_since"]).is_none()
            && at(&again, &["new_since"])
                == Some(&j(
                    r#"{"previous_utc": "2026-01-01T00:00:00Z", "sections": ["clinvar"]}"#
                )),
        "new since: og summary collect compares with the newest history/ folder: {:?}",
        at(&again, &["new_since"])
    );
    assert!(
        t2.matches(&format!("[{badge}]")).count() == 1
            && t2.contains(&format!("## ClinVar Pathogenic Screen\n---\n  [{badge}]")),
        "new since: the text marks the ClinVar section once"
    );

    let shown: Vec<&str> = NOTES
        .iter()
        .map(|x| x.1)
        .filter(|n| html.contains(&*e(n)))
        .collect();
    // Every caveat shown is its own caveat paragraph, at most once a card
    // (one note may sit on several cards: notes.rs lists where).
    let loose_of = |page: &str| -> Vec<String> {
        shown
            .iter()
            .copied()
            .filter(|n| {
                let wrapped = format!("<p class=\"caveat\">{}</p>", e(n));
                *n != PRS_ANCESTRY_NOTE
                    && (page.matches(&*e(n)).count() != page.matches(&wrapped).count()
                        || page
                            .split("<h2")
                            .any(|card| card.matches(&wrapped).count() > 1))
            })
            .map(|n| n[..40].to_owned())
            .collect()
    };
    let loose = loose_of(&html);
    assert!(
        shown.len() > 20 && loose.is_empty(),
        "every caveat shown is under a card title, once a card: {loose:?} {}",
        shown.len()
    );
    let eh = format!("<p class=\"caveat\">{}", e(EH_NOTE));
    assert!(html.contains(&eh), "the ExpansionHunter caveat is shown");
    let moved = html.replace(&eh, &format!("<p>{}", e(EH_NOTE)));
    assert_eq!(
        loose_of(&moved),
        vec![EH_NOTE[..40].to_owned()],
        "red control: a caveat left in a body is caught"
    );
    let twice = html.replacen(&eh, &format!("{eh}</p>\n{eh}"), 1);
    assert_eq!(
        loose_of(&twice),
        vec![EH_NOTE[..40].to_owned()],
        "red control: a caveat shown twice is caught"
    );
    let mo = html
        .split_once("<h2>Mosaic Loss of Chromosome Y or X</h2>")
        .map_or("", |x| x.1)
        .trim_start();
    let after_meta = mo
        .strip_prefix("<p class=\"sec-meta\">Step 33 ")
        .and_then(|r| r.split_once("</p>"))
        .filter(|(m, _)| !m.contains('\n'))
        .map(|x| x.1.trim_start());
    assert!(
        after_meta
            .is_some_and(|r| r.starts_with(&format!("<p class=\"caveat\">{}", e(MOSAIC_NOTE)))),
        "the caveat comes right after the title and the step line:\n{}",
        &mo[..mo.len().min(400)]
    );
}
