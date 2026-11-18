//! `og locityper report` (step 53): the table from what Locityper genotype
//! wrote. Moved from tests/test_locityper_report.py (its cases 1 to 4 and
//! 6); the summary's reading is in og-summary's `tests/optin_summary.rs`,
//! the report lines in og-render's `tests/optin_report.rs`, and run-all.sh's
//! plan and the DONE guards (cases 7 to 9) in the fake-docker case
//! `tests/fake-docker/port-report-optin-locityper.sh`.
//!
//! The output folders are made up per case, in the layout Locityper 1.7.14
//! writes (`loci/<locus>/res.json.gz`, a `success` file beside each
//! finished locus).

use std::collections::HashMap;
use std::fs;
use std::io::Write as _;
use std::num::NonZero;
use std::path::PathBuf;

const COLUMNS: &str = "locus\tregion\tstatus\tgenotype\thaplotype_1\thaplotype_2\tquality\ttotal_reads\tunexplained_reads\tweighted_distance\twarnings\tfilter\treason";
const BED: &str = "# a comment line
chr6\t32652603\t32669077\tHLA-DQB1
chr22\t42122284\t42152960\tCYP2D6
chr1\t161506761\t161681047\tFCGR
";

/// A gzip file, as Locityper writes `res.json.gz`.
fn gz(data: &[u8]) -> Vec<u8> {
    let mut w = og_io::bgzf::Writer::new(Vec::new(), NonZero::new(1).unwrap()).unwrap();
    w.write_all(data).unwrap();
    w.finish().unwrap()
}

/// One `res.json.gz` body; `None` leaves a key out.
struct Res {
    genotype: &'static str,
    total: Option<&'static str>,
    unexpl: Option<&'static str>,
    wd: Option<&'static str>,
    warnings: Option<&'static str>,
}

fn res() -> Res {
    Res {
        genotype: "HG00733.1,NA19240.2",
        total: Some("900"),
        unexpl: Some("12"),
        wd: Some("0.2"),
        warnings: None,
    }
}

impl Res {
    fn json(&self, pad: bool) -> String {
        let mut f = vec![
            format!("\"genotype\": \"{}\"", self.genotype),
            "\"quality\": 27.8".into(),
        ];
        for (k, v) in [
            ("total_reads", self.total),
            ("unexpl_reads", self.unexpl),
            ("weight_dist", self.wd),
            ("warnings", self.warnings),
        ] {
            if let Some(v) = v {
                f.push(format!("\"{k}\": {v}"));
            }
        }
        f.push("\"options\": 120, \"dist_type\": \"edit\"".into());
        if pad {
            f.push(format!("\"pad\": \"{}\"", "x".repeat(4000)));
        }
        format!("{{{}}}", f.join(", "))
    }
}

/// What `res.json.gz` holds.
enum Body {
    None,
    Res(Res),
    Raw(&'static [u8]),
    /// Cut in half, as a full disk leaves it.
    Truncated(Res),
    /// Its deflate body corrupt.
    Damaged(Res),
}

struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Dir {
        let d =
            std::env::temp_dir().join(format!("og-optin-locityper-{}-{tag}", std::process::id()));
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

fn gt_dir(root: &str, locus: &str, body: Body, success: bool) {
    let d = PathBuf::from(root).join("loci").join(locus);
    fs::create_dir_all(&d).unwrap();
    let bytes = match body {
        Body::None => None,
        Body::Res(r) => Some(gz(r.json(false).as_bytes())),
        Body::Raw(b) => Some(gz(b)),
        Body::Truncated(r) => {
            let w = gz(r.json(false).as_bytes());
            Some(w[..w.len() / 2].to_vec())
        }
        Body::Damaged(r) => {
            let mut w = gz(r.json(true).as_bytes());
            // Past the 18-byte BGZF header: the deflate data itself.
            for b in &mut w[20..40] {
                *b ^= 0xFF;
            }
            Some(w)
        }
    };
    if let Some(b) = bytes {
        fs::write(d.join("res.json.gz"), b).unwrap();
    }
    if success {
        fs::write(d.join("success"), "").unwrap();
    }
}

fn og(args: &[&str]) -> (u8, String, String) {
    let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
    let (mut o, mut e) = (Vec::new(), Vec::new());
    let code = og_reports::locityper::cli(Some("report"), &args, &mut o, &mut e);
    (
        code,
        String::from_utf8(o).unwrap(),
        String::from_utf8(e).unwrap(),
    )
}

type Row = HashMap<String, String>;

fn table(text: &str) -> Vec<Row> {
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    lines
        .map(|l| {
            head.iter()
                .zip(l.split('\t'))
                .map(|(k, v)| ((*k).to_owned(), v.to_owned()))
                .collect()
        })
        .collect()
}

/// The verb on CYP2D6 alone: its row.
fn og_row(d: &Dir, n: usize, body: Body, success: bool, exit_code: i64) -> Row {
    let root = d.p(&format!("c{n}"));
    gt_dir(&root, "CYP2D6", body, success);
    let bed = d.p(&format!("c{n}.bed"));
    fs::write(&bed, "chr22\t42122284\t42152960\tCYP2D6\n").unwrap();
    let out = d.p(&format!("c{n}.tsv"));
    let (code, _, err) = og(&[
        "--gt-dir",
        &root,
        "--loci",
        &bed,
        "--exit-code",
        &exit_code.to_string(),
        "--out",
        &out,
    ]);
    assert_eq!(code, 0, "{err}");
    table(&fs::read_to_string(&out).unwrap()).remove(0)
}

/// The control: a reading whose filter passes every called locus.
fn pass_everything(d: &Dir, n: usize, body: Body, success: bool, exit_code: i64) -> Row {
    let mut r = og_row(d, n, body, success, exit_code);
    if r["status"] == "called" {
        r.insert("filter".into(), "pass".into());
        r.insert("reason".into(), ".".into());
    }
    r
}

type Read = fn(&Dir, usize, Body, bool, i64) -> Row;

/// (description, held) of cases 1 to 3, read with `read`.
fn cases(read: Read, d: &Dir) -> Vec<(&'static str, bool, Row)> {
    let mut n = 0;
    let mut row = |body: Body, success: bool, exit_code: i64| {
        n += 1;
        read(d, n, body, success, exit_code)
    };
    let mut out = Vec::new();
    let r = row(Body::Res(res()), true, 0);
    out.push((
        "1. a clean call is called and passes, both haplotypes in their columns",
        r["status"] == "called"
            && r["filter"] == "pass"
            && r["haplotype_1"] == "HG00733.1"
            && r["haplotype_2"] == "NA19240.2"
            && r["quality"] == "27.8"
            && r["weighted_distance"] == "0.200"
            && r["reason"] == ".",
        r,
    ));
    let r = row(
        Body::Res(Res {
            wd: Some("31.5"),
            ..res()
        }),
        true,
        0,
    );
    out.push((
        "2. a weighted distance above 30 fails, saying so",
        r["filter"] == "fail" && r["reason"].contains("weighted distance 31.50"),
        r,
    ));
    let r = row(
        Body::Res(Res {
            total: Some("3000"),
            unexpl: Some("1200"),
            ..res()
        }),
        true,
        0,
    );
    out.push((
        "2. 1200 of 3000 unexplained reads fails",
        r["filter"] == "fail" && r["reason"].contains("1200 of 3000"),
        r,
    ));
    let r = row(
        Body::Res(Res {
            total: Some("10000"),
            unexpl: Some("1500"),
            ..res()
        }),
        true,
        0,
    );
    out.push((
        "2. 1500 of 10000 unexplained reads passes (under 20%)",
        r["filter"] == "pass",
        r,
    ));
    let r = row(
        Body::Res(Res {
            warnings: Some("[\"UnexplainedReads\"]"),
            ..res()
        }),
        true,
        0,
    );
    out.push((
        "2. a Locityper warning fails, naming it",
        r["filter"] == "fail"
            && r["warnings"] == "UnexplainedReads"
            && r["reason"].contains("UnexplainedReads"),
        r,
    ));
    let r = row(Body::Res(Res { wd: None, ..res() }), true, 0);
    out.push((
        "2. no weighted distance fails, pointing at locityper augment",
        r["filter"] == "fail" && r["reason"].contains("augment") && r["weighted_distance"] == ".",
        r,
    ));
    let r = row(
        Body::Res(Res {
            total: None,
            unexpl: None,
            ..res()
        }),
        true,
        0,
    );
    out.push((
        "2. no read counts fails, saying so, and the counts read '.' (not 0 of 0)",
        r["filter"] == "fail"
            && r["reason"].contains("no read counts")
            && r["total_reads"] == "."
            && r["unexplained_reads"] == ".",
        r,
    ));
    let r = row(
        Body::Res(Res {
            total: Some("null"),
            ..res()
        }),
        true,
        0,
    );
    out.push((
        "2. a null total_reads fails too",
        r["filter"] == "fail" && r["reason"].contains("no read counts"),
        r,
    ));
    let r = row(Body::Res(res()), false, 101);
    out.push((
        "3. no success file: not called, with Locityper's exit code",
        r["status"] == "not called"
            && r["reason"].contains("exit code 101")
            && r["genotype"] == "."
            && r["filter"] == ".",
        r,
    ));
    let r = row(Body::None, true, 0);
    out.push((
        "3. no res.json.gz: not called",
        r["status"] == "not called" && r["reason"].contains("no result"),
        r,
    ));
    let r = row(
        Body::Res(Res {
            genotype: "HG00733.1",
            ..res()
        }),
        true,
        0,
    );
    out.push((
        "3. one haplotype only: not called",
        r["status"] == "not called" && r["reason"].contains("no diploid"),
        r,
    ));
    for (what, body) in [
        (
            "3. an unreadable res.json.gz: not called",
            Body::Raw(b"{not json"),
        ),
        (
            "3. a truncated res.json.gz: not called, the table still written",
            Body::Truncated(res()),
        ),
        (
            "3. a corrupt deflate body: not called, the table still written",
            Body::Damaged(res()),
        ),
    ] {
        let r = row(body, true, 0);
        out.push((
            what,
            r["status"] == "not called" && r["reason"].contains("cannot be read"),
            r,
        ));
    }
    out
}

#[test]
fn the_filter_and_what_is_not_called() {
    let d = Dir::new("cases");
    let failed: Vec<String> = cases(og_row, &d)
        .into_iter()
        .filter(|c| !c.1)
        .map(|c| format!("{}: {:?}", c.0, c.2))
        .collect();
    assert!(failed.is_empty(), "{failed:#?}");
}

#[test]
fn control_a_reading_that_passes_everything_fails_the_filter_cases() {
    let d = Dir::new("control");
    let failed: Vec<&str> = cases(pass_everything, &d)
        .into_iter()
        .filter(|c| !c.1)
        .map(|c| c.0)
        .collect();
    assert!(failed.iter().any(|x| x.starts_with("2.")), "{failed:?}");
}

#[test]
fn the_cli_writes_a_row_per_bed_locus_in_order() {
    let d = Dir::new("cli");
    let bed = d.p("loci.bed");
    fs::write(&bed, BED).unwrap();
    let root = d.p("gt");
    gt_dir(&root, "CYP2D6", Body::Res(res()), true);
    gt_dir(
        &root,
        "HLA-DQB1",
        Body::Res(Res {
            wd: Some("50"),
            ..res()
        }),
        true,
    );
    let out = d.p("cli.tsv");
    let (code, stdout, err) = og(&[
        "--gt-dir",
        &root,
        "--loci",
        &bed,
        "--exit-code",
        "0",
        "--out",
        &out,
    ]);
    assert_eq!(code, 0, "{err}");
    let text = fs::read_to_string(&out).unwrap();
    let lines: Vec<Vec<&str>> = text.lines().map(|l| l.split('\t').collect()).collect();
    assert_eq!(lines[0].join("\t"), COLUMNS);
    let names: Vec<&str> = lines[1..].iter().map(|r| r[0]).collect();
    assert_eq!(names, ["HLA-DQB1", "CYP2D6", "FCGR"]);
    assert!(lines.iter().all(|r| r.len() == 13));
    assert_eq!(lines[1][1], "chr6:32652604-32669077");
    let st: Vec<&str> = lines[1..].iter().map(|r| r[2]).collect();
    assert_eq!(st, ["called", "called", "not called"]);
    assert!(stdout.contains("2 of 3 loci called, 1 pass"), "{stdout}");
    let bad = d.p("bad.bed");
    fs::write(&bad, "chr6\tx\t10\tA\n").unwrap();
    let (code, _, err) = og(&["--gt-dir", &root, "--loci", &bad, "--out", &d.p("bad.tsv")]);
    assert!(code == 1 && err.contains("ERROR"), "{code} {err}");
}
