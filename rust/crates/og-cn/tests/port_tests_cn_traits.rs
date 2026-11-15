//! `og cn call --locus HBA`: step 35's HBA1/HBA2 count from short-read
//! depth, through the `og cn` command line, on `samtools bedcov` tables made
//! up for each case (summed depth = mean depth x length).
//!
//! At 30x in the flanks, with a quarter of the gene-body depth from reads of
//! MAPQ 0 (as BWA-MEM and minimap2 give on HG002) unless a case says
//! otherwise:
//!   1. no deletion: 30x over HBA2 and HBA1 is four alpha-globin genes;
//!   2. a planted -a3.7-like loss (HBA2 at half depth) is three, a silent
//!      carrier, marked unconfirmed;
//!   3. both genes at half depth (--SEA-like) is two, with the warning that
//!      depth cannot tell --/aa from -a/-a;
//!   4. a triplication (HBA2 at 1.5x) is five;
//!   5. an estimate between two whole numbers, flanks under 8x, and flanks
//!      that disagree are each "not assessed", with the reason;
//!   6. no read of MAPQ 0 over the genes (GIAB's novoalign BAM of HG002) is
//!      "not assessed", naming the cause; 1% is too, 3% is called; no gene at
//!      all (no reads to judge by) is still 0;
//!   7. the command line: `og cn bed --locus HBA` prints the four regions,
//!      `og cn call --table` writes the metric table and exits 0; a table
//!      with other regions, or a short row, in either input exits 2, and so
//!      does a missing `--bedcov-mapq1`;
//!   8. the demo sample's HBA table (`tests/demo/make_demo_sample.py`) is the
//!      one `og cn call` writes for the demo's depths.
//!
//! This was `tests/test_hba_copy_number.py`.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

/// Step 35's HBA table, in its order.
const ORDER: [&str; 9] = [
    "status",
    "alpha_globin_genes",
    "estimate",
    "hba2_depth",
    "hba1_depth",
    "flank_left_depth",
    "flank_right_depth",
    "mapq0_share",
    "message",
];

fn catalog() -> String {
    format!("{REPO}/bin/cn_loci.tsv")
}

/// `og cn ARGS`: (exit code, stdout, stderr).
fn og_cn(args: &[&str]) -> (u8, String, String) {
    let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = og_cn::cli::main(&args, &mut out, &mut err);
    (
        code,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

/// The HBA row's regions as `og cn bed --locus HBA` prints them:
/// (chrom, start, end, name).
fn regions() -> Vec<(String, i64, i64, String)> {
    let (code, out, err) = og_cn(&["bed", "--catalog", &catalog(), "--locus", "HBA"]);
    assert_eq!(code, 0, "og cn bed: {err}");
    out.lines()
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (
                c[0].to_owned(),
                c[1].parse().unwrap(),
                c[2].parse().unwrap(),
                c[3].to_owned(),
            )
        })
        .collect()
}

/// A fresh folder for one case.
fn dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-cn-hba-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

type Depth = BTreeMap<&'static str, f64>;

fn depths(hba2: f64, hba1: f64, left: f64, right: f64) -> Depth {
    BTreeMap::from([
        ("flank_left", left),
        ("HBA2", hba2),
        ("HBA1", hba1),
        ("flank_right", right),
    ])
}

fn at(hba2: f64, hba1: f64) -> Depth {
    depths(hba2, hba1, 30.0, 30.0)
}

/// The depths of reads of MAPQ 1 or more: the flanks whole, the genes
/// without the given share.
fn mapq1(depth: &Depth, share: f64) -> Depth {
    depth
        .iter()
        .map(|(n, v)| {
            let v = if matches!(*n, "HBA2" | "HBA1") {
                v * (1.0 - share)
            } else {
                *v
            };
            (*n, v)
        })
        .collect()
}

/// A bedcov table with the given mean depth per region name.
fn bedcov(path: &Path, depth: &Depth) {
    let mut text = String::new();
    for (c, s, e, n) in regions() {
        // Truncated toward zero, as samtools' integer sum.
        let sum = (depth[n.as_str()] * (e - s) as f64) as i64;
        text.push_str(&format!("{c}\t{s}\t{e}\t{n}\t{sum}\n"));
    }
    std::fs::write(path, text).unwrap();
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

/// The command line of one call on the tables in `d`.
fn call_args(d: &Path, with_mapq1: bool) -> Vec<String> {
    let mut a = vec![
        "call".to_owned(),
        "--catalog".to_owned(),
        catalog(),
        "--locus".to_owned(),
        "HBA".to_owned(),
        "--bedcov".to_owned(),
        s(&d.join("cov.tsv")),
    ];
    if with_mapq1 {
        a.extend(["--bedcov-mapq1".to_owned(), s(&d.join("cov_q1.tsv"))]);
    }
    a.extend(["--table".to_owned(), s(&d.join("hba.tsv"))]);
    a
}

fn call(args: &[String]) -> (u8, String, String) {
    og_cn(&args.iter().map(String::as_str).collect::<Vec<_>>())
}

/// The metric table `og cn call` wrote, in its order; empty when it wrote none.
fn table(path: &Path) -> Vec<(String, String)> {
    let Ok(text) = std::fs::read_to_string(path) else {
        return Vec::new();
    };
    text.lines()
        .skip(1)
        .map(|l| {
            let (k, v) = l.split_once('\t').unwrap();
            (k.to_owned(), v.to_owned())
        })
        .collect()
}

/// Step 35's HBA table for these depths, as `og cn call` writes it.
fn run_share(name: &str, depth: &Depth, share: f64) -> BTreeMap<String, String> {
    let d = dir(name);
    bedcov(&d.join("cov.tsv"), depth);
    bedcov(&d.join("cov_q1.tsv"), &mapq1(depth, share));
    let (code, _, err) = call(&call_args(&d, true));
    assert_eq!(code, 0, "{name}: og cn call: {err}");
    let t = table(&d.join("hba.tsv")).into_iter().collect();
    let _ = std::fs::remove_dir_all(&d);
    t
}

fn run(name: &str, depth: &Depth) -> BTreeMap<String, String> {
    run_share(name, depth, 0.25)
}

fn get<'a>(t: &'a BTreeMap<String, String>, k: &str) -> &'a str {
    t.get(k).map_or("<missing>", String::as_str)
}

#[test]
fn no_deletion_is_four_genes_called() {
    let r = run("four", &at(30.0, 30.0));
    assert_eq!(
        (
            get(&r, "status"),
            get(&r, "alpha_globin_genes"),
            get(&r, "estimate")
        ),
        ("called", "4", "4.00"),
        "{r:?}"
    );
    assert!(get(&r, "message").contains("no deletion seen"), "{r:?}");
}

#[test]
fn hba2_at_half_depth_is_three_genes_a_silent_carrier() {
    let r = run("three", &at(15.0, 30.0));
    assert_eq!(
        (get(&r, "status"), get(&r, "alpha_globin_genes")),
        ("called", "3"),
        "{r:?}"
    );
    assert!(get(&r, "message").contains("silent carrier"), "{r:?}");
}

#[test]
fn a_count_under_four_says_it_is_unconfirmed_and_four_does_not() {
    let r = run("unconfirmed", &at(15.0, 30.0));
    let m = get(&r, "message");
    assert!(
        m.contains("unconfirmed") && m.contains("needs a clinical test"),
        "{r:?}"
    );
    for (name, d) in [("four-u", at(30.0, 30.0)), ("five-u", at(45.0, 30.0))] {
        let r = run(name, &d);
        assert!(!get(&r, "message").contains("unconfirmed"), "{name}: {r:?}");
    }
}

#[test]
fn both_genes_at_half_depth_is_two_and_depth_cannot_tell_which() {
    let r = run("two", &at(15.0, 15.0));
    assert_eq!(
        (get(&r, "status"), get(&r, "alpha_globin_genes")),
        ("called", "2"),
        "{r:?}"
    );
    assert!(
        get(&r, "message").contains("cannot tell these apart"),
        "{r:?}"
    );
}

#[test]
fn hba2_at_one_and_a_half_times_is_five_a_triplication() {
    let r = run("five", &at(45.0, 30.0));
    assert_eq!(get(&r, "alpha_globin_genes"), "5", "{r:?}");
    assert!(get(&r, "message").contains("triplication"), "{r:?}");
}

#[test]
fn an_estimate_of_3_50_is_not_assessed_naming_both_whole_numbers() {
    let r = run("between", &at(22.5, 30.0));
    assert_eq!(
        (get(&r, "status"), get(&r, "estimate")),
        ("not assessed", "3.50"),
        "{r:?}"
    );
    assert!(get(&r, "message").contains("between 3 and 4"), "{r:?}");
}

#[test]
fn flanks_under_8x_are_not_assessed_too_low() {
    let r = run("low", &depths(5.0, 5.0, 6.0, 6.0));
    assert_eq!(
        (get(&r, "status"), get(&r, "alpha_globin_genes")),
        ("not assessed", "."),
        "{r:?}"
    );
    assert!(get(&r, "message").contains("too low"), "{r:?}");
}

#[test]
fn flanks_that_disagree_are_not_assessed() {
    let r = run("disagree", &depths(30.0, 30.0, 30.0, 15.0));
    assert_eq!(get(&r, "status"), "not assessed", "{r:?}");
    assert!(get(&r, "message").contains("disagree"), "{r:?}");
}

#[test]
fn no_mapq0_read_over_the_genes_is_not_assessed_naming_the_aligner() {
    // The depths of GIAB's novoalign BAM of HG002, which gives no read MAPQ 0.
    let r = run_share("novo", &depths(15.5, 25.68, 30.11, 31.25), 0.0);
    assert_eq!(
        (
            get(&r, "status"),
            get(&r, "alpha_globin_genes"),
            get(&r, "mapq0_share")
        ),
        ("not assessed", ".", "0.000"),
        "{r:?}"
    );
    let m = get(&r, "message");
    assert!(
        m.contains("mapping quality 0") && m.contains("novoalign"),
        "{r:?}"
    );
}

#[test]
fn no_mapq0_read_with_an_estimate_of_3_is_not_assessed_never_three() {
    let r = run_share("novo-three", &at(15.0, 30.0), 0.0);
    assert_eq!(
        (get(&r, "status"), get(&r, "estimate")),
        ("not assessed", "3.00"),
        "{r:?}"
    );
    assert!(get(&r, "message").contains("mapping quality 0"), "{r:?}");
}

#[test]
fn one_percent_at_mapq0_is_not_assessed_three_percent_is_called() {
    let r = run_share("one-pct", &at(30.0, 30.0), 0.01);
    assert_eq!(get(&r, "status"), "not assessed", "{r:?}");
    assert!(get(&r, "message").contains("1.0%"), "{r:?}");
    let r = run_share("three-pct", &at(30.0, 30.0), 0.03);
    assert_eq!(
        (
            get(&r, "status"),
            get(&r, "alpha_globin_genes"),
            get(&r, "mapq0_share")
        ),
        ("called", "4", "0.030"),
        "{r:?}"
    );
}

#[test]
fn no_gene_at_all_is_still_called_0() {
    let r = run_share("none", &at(0.0, 0.0), 0.0);
    assert_eq!(
        (get(&r, "status"), get(&r, "alpha_globin_genes")),
        ("called", "0"),
        "{r:?}"
    );
    assert!(get(&r, "message").contains("Bart"), "{r:?}");
}

#[test]
fn bed_prints_four_regions_the_genes_between_the_flanks() {
    let names: Vec<String> = regions().into_iter().map(|r| r.3).collect();
    assert_eq!(names, ["flank_left", "HBA2", "HBA1", "flank_right"]);
}

#[test]
fn call_exits_0_and_writes_the_table_in_its_order() {
    let d = dir("cli");
    bedcov(&d.join("cov.tsv"), &at(15.0, 30.0));
    bedcov(&d.join("cov_q1.tsv"), &mapq1(&at(15.0, 30.0), 0.25));
    let (code, out, err) = call(&call_args(&d, true));
    let t = table(&d.join("hba.tsv"));
    assert_eq!(code, 0, "{err}");
    let keys: Vec<&str> = t.iter().map(|(k, _)| k.as_str()).collect();
    assert_eq!(keys, ORDER);
    let m: BTreeMap<_, _> = t.into_iter().collect();
    assert_eq!(
        (
            get(&m, "status"),
            get(&m, "alpha_globin_genes"),
            get(&m, "mapq0_share")
        ),
        ("called", "3", "0.250")
    );
    assert!(out.contains("HBA1/HBA2: three"), "{out}");
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn call_without_bedcov_mapq1_exits_2() {
    let d = dir("no-q1");
    bedcov(&d.join("cov.tsv"), &at(15.0, 30.0));
    let (code, _, err) = call(&call_args(&d, false));
    assert_eq!(code, 2, "{err}");
    assert!(err.contains("--bedcov-mapq1"), "{err}");
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn call_refuses_other_regions_or_a_short_row_in_either_table() {
    for (name, file) in [("all reads", "cov.tsv"), ("MAPQ 1 or more", "cov_q1.tsv")] {
        let d = dir(&format!("broken-{file}"));
        let rewrite = |d: &Path| {
            bedcov(&d.join("cov.tsv"), &at(15.0, 30.0));
            bedcov(&d.join("cov_q1.tsv"), &mapq1(&at(15.0, 30.0), 0.25));
        };
        rewrite(&d);
        let p = d.join(file);
        let mut text = std::fs::read_to_string(&p).unwrap();
        text.push_str("chr16\t1\t100\tother\t500\n");
        std::fs::write(&p, text).unwrap();
        let (code, _, err) = call(&call_args(&d, true));
        assert_eq!(code, 2, "a table ({name}) with other regions: {err}");
        assert!(err.contains("differ"), "{name}: {err}");
        rewrite(&d);
        std::fs::write(&p, "chr16\t100000\t140000\n").unwrap();
        let (code, _, err) = call(&call_args(&d, true));
        assert_eq!(code, 2, "a short row ({name}): {err}");
        assert!(err.contains("expected chrom"), "{name}: {err}");
        let _ = std::fs::remove_dir_all(&d);
    }
}

/// The `{"name": number, ...}` literal assigned to `var` in `src`.
fn py_dict(src: &str, var: &str) -> Depth {
    let line = src
        .lines()
        .find(|l| l.starts_with(&format!("{var} = {{")))
        .unwrap_or_else(|| panic!("the demo builder has no {var} line"));
    let body = line.split_once('{').unwrap().1.rsplit_once('}').unwrap().0;
    body.split(',')
        .map(|kv| {
            let (k, v) = kv.split_once(':').unwrap();
            let k = k.trim().trim_matches('"');
            let name = ["flank_left", "HBA2", "HBA1", "flank_right"]
                .into_iter()
                .find(|n| *n == k)
                .unwrap_or_else(|| panic!("{var}: unknown region {k}"));
            (name, v.trim().parse().unwrap())
        })
        .collect()
}

/// The string the parenthesised literals assigned to `var` in `src` join to.
fn py_str(src: &str, var: &str) -> String {
    let start = src
        .find(&format!("\n{var} = ("))
        .unwrap_or_else(|| panic!("the demo builder has no {var} literal"));
    let rest = &src[start..];
    let body = &rest[rest.find('(').unwrap() + 1..rest.find(")\n").unwrap()];
    let mut out = String::new();
    for piece in body.split('"').skip(1).step_by(2) {
        out.push_str(&piece.replace("\\t", "\t").replace("\\n", "\n"));
    }
    out
}

#[test]
fn the_demo_hba_table_is_what_og_cn_call_writes_for_its_depths() {
    let src = std::fs::read_to_string(format!("{REPO}/tests/demo/make_demo_sample.py")).unwrap();
    let depth = py_dict(&src, "HBA_DEPTH");
    let q1 = py_dict(&src, "HBA_Q1");
    let want = py_str(&src, "HBA_TABLE");
    assert!(want.starts_with("metric\tvalue\nstatus\t"), "{want:?}");
    let d = dir("demo");
    bedcov(&d.join("cov.tsv"), &depth);
    bedcov(
        &d.join("cov_q1.tsv"),
        &depth.iter().map(|(n, v)| (*n, v * q1[n])).collect(),
    );
    let (code, _, err) = call(&call_args(&d, true));
    assert_eq!(code, 0, "{err}");
    assert_eq!(std::fs::read_to_string(d.join("hba.tsv")).unwrap(), want);
    let _ = std::fs::remove_dir_all(&d);
}
