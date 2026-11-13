use std::path::PathBuf;

use super::add_callers::{self, Opts};
use super::postfilter;
use super::*;

/// A scratch folder of its own for one test.
fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-sv-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn put(dir: &std::path::Path, name: &str, text: &str) -> String {
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p.to_string_lossy().into_owned()
}

fn opts(base: &str) -> Opts {
    Opts {
        base: base.to_owned(),
        base_name: "manta".to_owned(),
        delly: None,
        sex: Sex::Unknown,
        par_bed: None,
        cnvpytor: None,
        max_eval: 0.01,
        max_q0: 0.5,
        min_overlap: 0.5,
    }
}

const HEAD: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##contig=<ID=chrX>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n";

fn ids(vcf: &str) -> Vec<String> {
    vcf.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').nth(2).unwrap().to_owned())
        .collect()
}

#[test]
fn the_depth_rule_by_sex_and_par() {
    let del =
        |cn, sex, par| depth_contradicts(SvKind::Del, "chrX", 5000, Rdcn::Copies(cn), sex, par);
    let dup =
        |cn, sex, par| depth_contradicts(SvKind::Dup, "chrX", 5000, Rdcn::Copies(cn), sex, par);
    // Red-first control: one copy is a male's normal on chrX outside the PARs.
    assert!(
        del(1, Sex::Male, false),
        "a male chrX DEL with RDCN 1 must be set aside"
    );
    assert!(!del(0, Sex::Male, false));
    assert!(!del(1, Sex::Female, false));
    assert!(!del(1, Sex::Unknown, false));
    assert!(del(2, Sex::Unknown, false));
    assert!(
        !del(1, Sex::Male, true),
        "inside the PAR, 2 copies are normal"
    );
    assert!(!dup(2, Sex::Male, false));
    assert!(dup(2, Sex::Female, false));
    assert!(!dup(2, Sex::Unknown, false));
    assert!(dup(1, Sex::Unknown, false));
    // Autosomes, short calls, no field, no evidence.
    assert!(depth_contradicts(
        SvKind::Del,
        "chr1",
        5000,
        Rdcn::Copies(2),
        Sex::Male,
        false
    ));
    assert!(!depth_contradicts(
        SvKind::Del,
        "chr1",
        5000,
        Rdcn::Copies(1),
        Sex::Male,
        false
    ));
    assert!(!depth_contradicts(
        SvKind::Del,
        "chr1",
        1000,
        Rdcn::Copies(2),
        Sex::Male,
        false
    ));
    assert!(!depth_contradicts(
        SvKind::Del,
        "chr1",
        5000,
        Rdcn::Absent,
        Sex::Male,
        false
    ));
    assert!(depth_contradicts(
        SvKind::Del,
        "chr1",
        5000,
        Rdcn::NoEvidence,
        Sex::Male,
        false
    ));
}

#[test]
fn rdcn_is_typed() {
    assert_eq!(Rdcn::from_format(&["GT"], &["0/1"]), Rdcn::Absent);
    assert_eq!(
        Rdcn::from_format(&["GT", "RDCN"], &["0/1"]),
        Rdcn::NoEvidence
    );
    assert_eq!(
        Rdcn::from_format(&["GT", "RDCN"], &["0/1", "-1"]),
        Rdcn::NoEvidence
    );
    assert_eq!(
        Rdcn::from_format(&["GT", "RDCN"], &["0/1", "1.5"]),
        Rdcn::NoEvidence
    );
    assert_eq!(
        Rdcn::from_format(&["GT", "RDCN"], &["0/1", "3"]),
        Rdcn::Copies(3)
    );
    assert_eq!(
        Rdcn::from_format(&["RDCN"], &["99999999999999999999999"]),
        Rdcn::Copies(u64::MAX)
    );
    assert_eq!(Rdcn::from_number(2.0), Rdcn::Copies(2));
    assert_eq!(Rdcn::from_number(1.5), Rdcn::NoEvidence);
    assert_eq!(Rdcn::from_number(-1.0), Rdcn::NoEvidence);
    assert_eq!(Rdcn::from_number(f64::INFINITY), Rdcn::NoEvidence);
}

#[test]
fn calls_the_sample_carries() {
    let rec = |filt: &'static str, gt: &'static str| {
        vec!["c", "1", ".", "N", "<DEL>", ".", filt, ".", "GT", gt]
    };
    assert!(call_holds(&rec("PASS", "0/1"), 1));
    assert!(call_holds(&rec(".", "0|1"), 1));
    assert!(!call_holds(&rec("MinGQ", "0/1"), 1));
    assert!(!call_holds(&rec("PASS", "0/0"), 1));
    assert!(!call_holds(&rec("PASS", "./."), 1));
    assert!(
        call_holds(&rec("PASS", "0/0")[..8], 0),
        "sites-only: no genotype"
    );
}

#[test]
fn reciprocal_overlap_is_measured_on_the_longer_call() {
    let mut ix = OverlapIndex::default();
    ix.add("chr1", SvKind::Del, 10_000, 20_000);
    // Red-first control: 4 kb inside 10 kb is 40% of the longer call.
    assert!(!ix.covers("chr1", SvKind::Del, 10_000, 14_000, 0.5));
    assert!(ix.covers("chr1", SvKind::Del, 10_000, 16_000, 0.5));
    assert!(!ix.covers("chr1", SvKind::Dup, 10_000, 20_000, 0.5));
    assert!(!ix.covers("chr2", SvKind::Del, 10_000, 20_000, 0.5));
    // A long call far left still reaches a query through max_len.
    ix.add("chr1", SvKind::Del, 1_000, 100_000);
    assert!(ix.covers("chr1", SvKind::Del, 40_000, 100_000, 0.5));
}

#[test]
fn chrom_keys_and_info_values() {
    assert_eq!(chrom_key("chr20"), "20");
    assert_eq!(chrom_key("CHRX"), "X");
    assert_eq!(chrom_key("20"), "20");
    assert_eq!(chrom_key("ch"), "ch");
    assert_eq!(
        info_value("A=1;SVCALLER=delly;SVCALLER=x", "SVCALLER"),
        "delly"
    );
    assert_eq!(info_value("SVCALLER", "SVCALLER"), "");
    assert_eq!(info_value("", "SVCALLER"), "");
}

#[test]
fn add_callers_adds_what_no_call_covers() {
    let d = scratch("add");
    let base = put(
        &d,
        "base.vcf",
        &format!(
            "{HEAD}chr1\t10000\tM1\tN\t<DEL>\t50\tPASS\tEND=20000;SVTYPE=DEL\tGT\t0/1\nchr1\t500\tM0\tN\t<INV>\t50\tPASS\t.\tGT\t0/1\n"
        ),
    );
    let delly = put(
        &d,
        "delly.vcf",
        &format!(
            "{HEAD}chr1\t10000\tD_inside\tN\t<DEL>\t60\tPASS\tIMPRECISE;SVTYPE=DEL;END=14000\tGT:RDCN\t0/1:1\n\
             chr1\t10100\tD_same\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=19900\tGT:RDCN\t0/1:1\n\
             chrX\t50000000\tD_x1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=51000000\tGT:RDCN\t1/1:1\n\
             chr1\t1\tD_low\tN\t<DEL>\t60\tLowQual\tSVTYPE=DEL;END=9000\tGT:RDCN\t0/1:1\n"
        ),
    );
    let mut o = opts(&base);
    o.delly = Some(delly.clone());
    o.sex = Sex::Male;
    let (vcf, note) = add_callers::run_text(&o).unwrap();
    assert_eq!(ids(&vcf), ["M0", "M1", "D_inside"], "{vcf}");
    assert!(
        vcf.contains("\tEND=14000;SVTYPE=DEL;IMPRECISE;SVCALLER=delly\tGT\t0/1\n"),
        "{vcf}"
    );
    assert!(vcf.contains("\tSVCALLER=manta\tGT\t0/1\n"), "{vcf}");
    assert!(vcf.contains("##INFO=<ID=SVCALLER,"));
    assert_eq!(
        note,
        "Callers beside manta: delly: 2 PASS DEL/DUP calls, 1 added, and 1 PASS calls longer than 1000 bp set aside because read depth (RDCN) does not show the change"
    );
    o.sex = Sex::Female;
    let (vcf, _) = add_callers::run_text(&o).unwrap();
    assert!(ids(&vcf).contains(&"D_x1".to_owned()), "{vcf}");
    std::fs::remove_dir_all(d).unwrap();
}

#[test]
fn add_callers_input_errors_are_messages() {
    let d = scratch("adderr");
    let base = put(
        &d,
        "base.vcf",
        &format!("{HEAD}chr1\tx\tM1\tN\t<DEL>\t50\tPASS\t.\tGT\t0/1\n"),
    );
    let e = add_callers::run_text(&opts(&base)).unwrap_err().to_string();
    assert!(
        e.ends_with("record 1: POS 'x' is not a whole number"),
        "{e}"
    );
    let short = put(&d, "short.vcf", &format!("{HEAD}chr1\t5\tM1\n"));
    let e = add_callers::run_text(&opts(&short))
        .unwrap_err()
        .to_string();
    assert!(
        e.ends_with("record 1 has fewer than the 8 columns of a VCF record"),
        "{e}"
    );
    let gone = d.join("gone.vcf").to_string_lossy().into_owned();
    let e = add_callers::run_text(&opts(&gone)).unwrap_err().to_string();
    assert!(e.contains("No such file or directory"), "{e}");
    std::fs::remove_dir_all(d).unwrap();
}

const C45_HEAD: &str = "AnnotSV_ID\tSV_chrom\tSV_start\tSV_end\tSV_length\tSV_type\tID\tINFO\tFORMAT\tS1\tAnnotation_mode\tGene_name\tB_loss_AFmax\tACMG_class\tAnnotSV_ranking_score\n";

fn c45_row(id: &str, start: i64, end: i64, t: &str, caller: &str, fmt: &str, s: &str) -> String {
    format!(
        "a\t20\t{start}\t{end}\t{}\t{t}\t{id}\tSVTYPE={t};SVCALLER={caller}\t{fmt}\t{s}\tfull\tG\t\t4\t0.9\n",
        start - end
    )
}

fn postfilter_rows(
    d: &std::path::Path,
    body: &str,
) -> Vec<std::collections::HashMap<String, String>> {
    let c45 = put(d, "c45.tsv", &format!("{C45_HEAD}{body}"));
    let out = d.join("checked.tsv").to_string_lossy().into_owned();
    let o = postfilter::Opts {
        class45: c45,
        genes: None,
        rules: None,
        consensus: None,
        sex: Sex::Unknown,
        par_bed: None,
        out: out.clone(),
    };
    postfilter::run(&o).unwrap();
    let head = postfilter::header_of(&out);
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| {
            head.iter()
                .cloned()
                .zip(l.split('\t').map(str::to_owned))
                .collect()
        })
        .collect()
}

#[test]
fn duphold_thresholds_move_the_class() {
    let d = scratch("duphold");
    let rows = postfilter_rows(
        &d,
        &[
            c45_row("D075", 1000, 50000, "DEL", "manta", "GT:DHFFC", "0/1:0.75"),
            c45_row("D065", 1000, 50000, "DEL", "manta", "GT:DHFFC", "0/1:0.65"),
            c45_row("U125", 1000, 50000, "DUP", "manta", "GT:DHBFC", "0/1:1.25"),
            c45_row("R15", 1000, 50000, "DEL", "delly", "GT:RDCN", "0/1:1.5"),
        ]
        .concat(),
    );
    let get = |i: usize, c: &str| rows[i][c].clone();
    // Red-first control: DHFFC 0.75 is not below 0.7.
    assert_eq!(
        (get(0, "Depth_support"), get(0, "ACMG_class")),
        ("no".into(), "3".into())
    );
    assert_eq!(get(0, "Rule"), "no_depth_support");
    assert_eq!(
        (get(1, "Depth_support"), get(1, "ACMG_class")),
        ("yes".into(), "4".into())
    );
    assert_eq!(get(2, "Depth_support"), "no");
    // A fractional RDCN is no evidence, not "no".
    assert_eq!(
        (get(3, "Depth_support"), get(3, "ACMG_class")),
        ("not checked".into(), "4".into())
    );
    std::fs::remove_dir_all(d).unwrap();
}

#[test]
fn a_missing_table_is_a_message() {
    let d = scratch("pfmissing");
    let o = postfilter::Opts {
        class45: d.join("none.tsv").to_string_lossy().into_owned(),
        genes: None,
        rules: None,
        consensus: None,
        sex: Sex::Unknown,
        par_bed: None,
        out: d.join("o.tsv").to_string_lossy().into_owned(),
    };
    let e = postfilter::run(&o).unwrap_err().to_string();
    assert!(
        e.starts_with("cannot read ") && e.contains("No such file"),
        "{e}"
    );
    assert!(!d.join("o.tsv").exists());
    std::fs::remove_dir_all(d).unwrap();
}
