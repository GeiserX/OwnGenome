//! `og locityper report`: the authors' filter reads pass or fail only at
//! the HLA and KIR loci they calibrated it on. At CYP2D6 and FCGR it reads
//! "uncalibrated", never a green pass, and the reason says what the rule
//! alone would give there.
//!
//! Red first: before the fix a clean CYP2D6 call read "pass" and a far
//! one "fail", the same as at HLA-DQB1.

use std::fs;
use std::io::Write as _;
use std::num::NonZero;
use std::path::PathBuf;

/// The five loci of `assets/locityper_loci.bed`.
const BED: &str = "# the shipped loci
chr1\t161506761\t161681047\tFCGR
chr6\t32484347\t32593538\tHLA-DRB1
chr6\t32652603\t32669077\tHLA-DQB1
chr19\t54723153\t54868521\tKIR
chr22\t42122284\t42152960\tCYP2D6
";

/// A clean call: no warning, weighted distance 0.2, 12 of 900 reads
/// unexplained.
const CLEAN: &str = r#"{"genotype": "HG00733.1,NA19240.2", "quality": 27.8, "total_reads": 900, "unexpl_reads": 12, "weight_dist": 0.2, "options": 120, "dist_type": "edit"}"#;
/// The same call with a weighted distance of 31.5, above the authors' 30.
const FAR: &str = r#"{"genotype": "HG00733.1,NA19240.2", "quality": 27.8, "total_reads": 900, "unexpl_reads": 12, "weight_dist": 31.5, "options": 120, "dist_type": "edit"}"#;

struct Dir(PathBuf);

impl Dir {
    fn new(name: &str) -> Dir {
        let p =
            std::env::temp_dir().join(format!("og-i411-locityper-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&p);
        fs::create_dir_all(&p).unwrap();
        Dir(p)
    }

    fn p(&self, s: &str) -> String {
        self.0.join(s).to_str().unwrap().to_owned()
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// `res.json.gz` and `success` of one locus, as Locityper 1.7.14 writes
/// them.
fn locus(root: &str, name: &str, json: &str) {
    let d = PathBuf::from(root).join("loci").join(name);
    fs::create_dir_all(&d).unwrap();
    let mut w = og_io::bgzf::Writer::new(Vec::new(), NonZero::new(1).unwrap()).unwrap();
    w.write_all(json.as_bytes()).unwrap();
    fs::write(d.join("res.json.gz"), w.finish().unwrap()).unwrap();
    fs::write(d.join("success"), "").unwrap();
}

/// The verb on the five loci, each with `json(locus)`: (stdout, rows as
/// (locus, filter, reason)).
fn report(d: &Dir, json: fn(&str) -> &'static str) -> (String, Vec<[String; 3]>) {
    let bed = d.p("loci.bed");
    fs::write(&bed, BED).unwrap();
    let root = d.p("gt");
    for name in ["FCGR", "HLA-DRB1", "HLA-DQB1", "KIR", "CYP2D6"] {
        locus(&root, name, json(name));
    }
    let out = d.p("out.tsv");
    let args: Vec<String> = ["--gt-dir", &root, "--loci", &bed, "--out", &out]
        .iter()
        .map(|s| (*s).to_owned())
        .collect();
    let (mut o, mut e) = (Vec::new(), Vec::new());
    let code = og_reports::locityper::cli(Some("report"), &args, &mut o, &mut e);
    assert_eq!(code, 0, "{}", String::from_utf8_lossy(&e));
    let text = fs::read_to_string(&out).unwrap();
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |n: &str| head.iter().position(|h| *h == n).unwrap();
    let (f, r) = (col("filter"), col("reason"));
    let rows = lines
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            assert_eq!(c.len(), 13, "{l}");
            [c[0].to_owned(), c[f].to_owned(), c[r].to_owned()]
        })
        .collect();
    (String::from_utf8(o).unwrap(), rows)
}

fn filters(rows: &[[String; 3]]) -> Vec<(&str, &str)> {
    rows.iter()
        .map(|r| (r[0].as_str(), r[1].as_str()))
        .collect()
}

#[test]
fn red_first_a_clean_cyp2d6_call_reads_uncalibrated_not_pass() {
    let d = Dir::new("clean");
    let (stdout, rows) = report(&d, |_| CLEAN);
    assert_eq!(
        filters(&rows),
        [
            ("FCGR", "uncalibrated"),
            ("HLA-DRB1", "pass"),
            ("HLA-DQB1", "pass"),
            ("KIR", "pass"),
            ("CYP2D6", "uncalibrated"),
        ]
    );
    for r in rows.iter().filter(|r| r[1] == "pass") {
        assert_eq!(r[2], ".", "{r:?}");
    }
    for r in rows.iter().filter(|r| r[1] == "uncalibrated") {
        assert_eq!(
            r[2],
            "the authors calibrated this filter on HLA and KIR loci only; by its rule alone it passes",
            "{r:?}"
        );
    }
    assert_eq!(
        stdout,
        "5 of 5 loci called, 3 pass the authors' filter, 2 uncalibrated (outside HLA and KIR)\n"
    );
}

#[test]
fn a_far_call_fails_at_hla_and_kir_and_reads_uncalibrated_elsewhere() {
    let d = Dir::new("far");
    let (stdout, rows) = report(&d, |_| FAR);
    assert_eq!(
        filters(&rows),
        [
            ("FCGR", "uncalibrated"),
            ("HLA-DRB1", "fail"),
            ("HLA-DQB1", "fail"),
            ("KIR", "fail"),
            ("CYP2D6", "uncalibrated"),
        ]
    );
    for r in &rows {
        let why = "weighted distance 31.50 is above 30";
        if r[1] == "fail" {
            assert_eq!(r[2], why, "{r:?}");
        } else {
            assert_eq!(
                r[2],
                format!(
                    "the authors calibrated this filter on HLA and KIR loci only; by its rule alone it fails: {why}"
                ),
                "{r:?}"
            );
        }
    }
    assert_eq!(
        stdout,
        "5 of 5 loci called, 0 pass the authors' filter, 2 uncalibrated (outside HLA and KIR)\n"
    );
}

#[test]
fn hla_and_kir_loci_by_name() {
    use og_reports::locityper::calibrated;
    for l in ["HLA-DRB1", "HLA-DQB1", "HLA-A", "KIR", "KIR2DL1"] {
        assert!(calibrated(l), "{l}");
    }
    for l in ["CYP2D6", "FCGR", "FCGR3A", "C4", "LPA", "hla-a", "HLA", ""] {
        assert!(!calibrated(l), "{l}");
    }
}

#[test]
fn only_hla_and_kir_calls_count_as_passing() {
    // With every locus clean, HLA-DQB1 alone passes: the CYP2D6 row adds
    // nothing to the pass count.
    let d = Dir::new("mixed");
    let (stdout, rows) = report(&d, |l| if l == "HLA-DQB1" { CLEAN } else { FAR });
    assert_eq!(
        filters(&rows),
        [
            ("FCGR", "uncalibrated"),
            ("HLA-DRB1", "fail"),
            ("HLA-DQB1", "pass"),
            ("KIR", "fail"),
            ("CYP2D6", "uncalibrated"),
        ]
    );
    assert!(
        stdout.starts_with("5 of 5 loci called, 1 pass the authors' filter, 2 uncalibrated"),
        "{stdout}"
    );
}
