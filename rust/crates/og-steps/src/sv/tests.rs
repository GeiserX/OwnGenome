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
    assert_eq!(chrom_key("éé"), "éé");
    assert_eq!(chrom_key("é"), "é");
    assert_eq!(chrom_key("☃1"), "☃1");
    assert_eq!(chrom_key("chr☃"), "☃");
    assert_eq!(chrom_key("cHrX"), "X");
    assert_eq!(chrom_key("chr"), "");
    assert_eq!(chrom_key(""), "");
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
        vcf.contains("\tEND=14000;SVTYPE=DEL;IMPRECISE;SVCALLER=delly\tGT:RDCN\t0/1:1\n"),
        "a Delly call keeps its RDCN beside GT: {vcf}"
    );
    assert!(vcf.contains("\tSVCALLER=manta\tGT\t0/1\n"), "{vcf}");
    assert!(vcf.contains("##INFO=<ID=SVCALLER,"));
    assert!(
        vcf.contains("\n##FORMAT=<ID=RDCN,Number=1,Type=Integer,"),
        "{vcf}"
    );
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
fn add_callers_keeps_rdcn_as_an_integer_field() {
    let d = scratch("add-rdcn-cell");
    let base = put(&d, "base.vcf", HEAD);
    // Calls of 1,000 bp or less pass add-callers whatever their RDCN says.
    let delly = put(
        &d,
        "delly.vcf",
        &format!(
            "{HEAD}chr1\t1000\tFRAC\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=1500\tGT:RDCN\t0/1:1.5\n\
             chr1\t3000\tNEG\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=3500\tGT:RDCN\t0/1:-1\n\
             chr1\t5000\tCUT\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=5500\tGT:GQ:RDCN\t0/1:40\n"
        ),
    );
    let mut o = opts(&base);
    o.delly = Some(delly);
    let (vcf, _) = add_callers::run_text(&o).unwrap();
    let format_of = |id: &str| {
        vcf.lines()
            .find(|l| l.split('\t').nth(2) == Some(id))
            .map(|l| l.split('\t').skip(8).collect::<Vec<_>>().join("\t"))
    };
    assert_eq!(
        (format_of("FRAC"), format_of("NEG"), format_of("CUT")),
        (
            Some("GT:RDCN\t0/1:.".to_owned()),
            Some("GT:RDCN\t0/1:-1".to_owned()),
            Some("GT:RDCN\t0/1:.".to_owned())
        ),
        "RDCN goes along as a whole number, else `.` (the header says Integer): {vcf}"
    );
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
    consensus: Option<String>,
) -> Vec<std::collections::HashMap<String, String>> {
    let c45 = put(d, "c45.tsv", &format!("{C45_HEAD}{body}"));
    let out = d.join("checked.tsv").to_string_lossy().into_owned();
    let o = postfilter::Opts {
        class45: c45,
        genes: None,
        rules: None,
        consensus,
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
        None,
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
fn a_delly_call_without_rdcn_has_unknown_depth_support() {
    let d = scratch("missing-rdcn");
    let rows = postfilter_rows(
        &d,
        &[
            c45_row("ABSENT", 1000, 50000, "DEL", "delly", "GT", "0/1"),
            c45_row("MISSING", 1000, 50000, "DEL", "delly", "GT:RDCN", "0/1:."),
        ]
        .concat(),
        None,
    );
    for row in &rows {
        assert_eq!(row["Depth_support"], "unknown", "{}", row["ID"]);
        assert_eq!(row["ACMG_class"], "4", "{}", row["ID"]);
    }
    std::fs::remove_dir_all(d).unwrap();
}

const CONSENSUS_HEAD: &str =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmanta\tdelly\n";

#[test]
fn support_rejects_a_dup_at_the_same_position_as_a_del() {
    let d = scratch("consensus-type");
    let consensus = put(
        &d,
        "consensus.vcf",
        &format!(
            "{CONSENSUS_HEAD}chr20\t1000\tCLUSTER\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=5000;SUPP=2;SUPP_VEC=11\tGT:TY:ID\t0/1:DUP:SHARED\t0/1:DUP:D1\nchr20\t1000\tDECOY\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000;SUPP=2;SUPP_VEC=11\tGT:TY:ID\t0/1:DEL:M2\t0/1:DEL:D2\n"
        ),
    );
    let body = [
        c45_row("SHARED", 1000, 5000, "DEL", "manta", "GT", "0/1"),
        c45_row("SHARED", 1000, 5000, "DUP", "manta", "GT", "0/1"),
    ]
    .concat();
    let rows = postfilter_rows(&d, &body, Some(consensus));
    assert_eq!(rows[0]["Support"], "1");
    assert_eq!(rows[1]["Support"], "2");
    assert_eq!(rows[1]["Support_callers"], "manta,delly");
    std::fs::remove_dir_all(d).unwrap();
}

#[test]
fn support_rejects_filtered_records_and_non_voting_members() {
    let d = scratch("consensus-filtered");
    let consensus = put(
        &d,
        "consensus.vcf",
        &format!(
            "{CONSENSUS_HEAD}chr20\t1000\tFILTERED\tN\t<DEL>\t.\tLowQual\tSVTYPE=DEL;END=5000;SUPP=2;SUPP_VEC=11\tGT:TY:ID\t0/1:DEL:M1\t0/1:DEL:D1\n\
             chr20\t10000\tNOT_VOTED\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=15000;SUPP=1;SUPP_VEC=10\tGT:TY:ID\t0/1:DEL:M2\t0/0:DEL:D2\n"
        ),
    );
    let body = [
        c45_row("M1", 1000, 5000, "DEL", "manta", "GT", "0/1"),
        c45_row("D2", 10000, 15000, "DEL", "delly", "GT", "0/1"),
    ]
    .concat();
    let rows = postfilter_rows(&d, &body, Some(consensus));
    assert_eq!(rows[0]["Support"], "1");
    assert_eq!(rows[1]["Support_callers"], "delly");
    std::fs::remove_dir_all(d).unwrap();
}

#[test]
fn support_requires_the_original_id_even_at_the_old_distance_limit() {
    let d = scratch("consensus-distance");
    let consensus = put(
        &d,
        "consensus.vcf",
        &format!(
            "{CONSENSUS_HEAD}chr20\t2000\tCLUSTER\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=6000;SUPP=2;SUPP_VEC=11\tGT:TY:ID\t0/1:DEL:M1\t0/1:DEL:D1\n"
        ),
    );
    let body = [
        c45_row("M2", 1000, 5000, "DEL", "manta", "GT", "0/1"),
        c45_row("M1", 1000, 5000, "DEL", "manta", "GT", "0/1"),
        c45_row("D1", 2000, 6000, "DEL", "manta", "GT", "0/1"),
    ]
    .concat();
    let rows = postfilter_rows(&d, &body, Some(consensus));
    assert_eq!(rows[0]["Support"], "1");
    assert_eq!(rows[1]["Support"], "2");
    assert_eq!(rows[1]["Support_callers"], "manta,delly");
    assert_eq!(rows[2]["Support"], "1");
    std::fs::remove_dir_all(d).unwrap();
}

/// AnnotSV's class 4 and 5 rows for the records of a VCF, with the cells
/// AnnotSV copies from each record: ID, INFO, FORMAT and the sample.
fn annotsv_rows(vcf: &str) -> String {
    vcf.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let start: i64 = f[1].parse().unwrap();
            let end: i64 = info_value(f[7], "END").parse().unwrap();
            let t = info_value(f[7], "SVTYPE");
            format!(
                "a\t{}\t{start}\t{end}\t{}\t{t}\t{}\t{}\t{}\t{}\tfull\tG\t\t4\t0.9\n",
                chrom_key(f[0]),
                start - end,
                f[2],
                f[7],
                f[8],
                f[9]
            )
        })
        .collect()
}

#[test]
fn a_delly_call_step_05_adds_keeps_its_rdcn_for_the_depth_check() {
    let d = scratch("chain-rdcn");
    let base = put(
        &d,
        "base.vcf",
        &format!("{HEAD}chr1\t500000\tM1\tN\t<DEL>\t50\tPASS\tEND=520000;SVTYPE=DEL\tGT\t0/1\n"),
    );
    // KEEP: Delly 2.x FORMAT, RDCN 1 on a 100 kb deletion, which passes
    // add-callers' depth gate. NORDCN: a Delly VCF with no RDCN field.
    let delly = put(
        &d,
        "delly.vcf",
        &format!(
            "{HEAD}chr1\t1000000\tKEEP\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=1100000\tGT:GQ:FT:RDCN:DV\t0/1:40:PASS:1:10\n\
             chr1\t2000000\tNORDCN\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=2100000\tGT:GQ\t0/1:40\n"
        ),
    );
    let mut o = opts(&base);
    o.delly = Some(delly);
    let (vcf, _) = add_callers::run_text(&o).unwrap();
    assert_eq!(ids(&vcf), ["M1", "KEEP", "NORDCN"], "{vcf}");
    let rows = postfilter_rows(&d, &annotsv_rows(&vcf), None);
    let depth = |id: &str| {
        rows.iter()
            .find(|r| r["ID"] == id)
            .map(|r| (r["Depth_support"].as_str(), r["ACMG_class"].as_str()))
    };
    assert_eq!(
        (depth("KEEP"), depth("NORDCN")),
        (Some(("yes", "4")), Some(("unknown", "4"))),
        "the RDCN add-callers checked reaches the depth check; only a call with no RDCN is unknown: {vcf}"
    );
    assert!(
        vcf.contains(
            "\tKEEP\tN\t<DEL>\t60\tPASS\tEND=1100000;SVTYPE=DEL;SVCALLER=delly\tGT:RDCN\t0/1:1\n"
        ) && vcf.contains(
            "\tNORDCN\tN\t<DEL>\t60\tPASS\tEND=2100000;SVTYPE=DEL;SVCALLER=delly\tGT\t0/1\n"
        ),
        "{vcf}"
    );
    std::fs::remove_dir_all(d).unwrap();
}

#[test]
fn a_source_id_in_two_clusters_has_unknown_support() {
    let d = scratch("consensus-twice");
    let consensus = put(
        &d,
        "consensus.vcf",
        &format!(
            "{CONSENSUS_HEAD}chr20\t1000\tA\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000;SUPP=2;SUPP_VEC=11\tGT:TY:ID\t0/1:DEL:TWICE\t0/1:DEL:D1\n\
             chr20\t9000\tB\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=12000;SUPP=1;SUPP_VEC=10\tGT:TY:ID\t0/1:DEL:TWICE\t./.:NaN:NaN\n\
             chr20\t20000\tC\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=25000;SUPP=2;SUPP_VEC=11\tGT:TY:ID\t0/1:DEL:ONCE\t0/1:DEL:D3\n"
        ),
    );
    let body = [
        c45_row("TWICE", 1000, 5000, "DEL", "manta", "GT", "0/1"),
        c45_row("ONCE", 20000, 25000, "DEL", "manta", "GT", "0/1"),
    ]
    .concat();
    let rows = postfilter_rows(&d, &body, Some(consensus));
    let support = |i: usize| {
        (
            rows[i]["Support"].as_str(),
            rows[i]["Support_callers"].as_str(),
        )
    };
    assert_eq!(
        support(0),
        ("", ""),
        "a source ID in two consensus records of one contig and type has no unique event: unknown"
    );
    assert_eq!(support(1), ("2", "manta,delly"), "control: a unique member");
    std::fs::remove_dir_all(d).unwrap();
}

#[test]
fn support_is_unknown_when_the_consensus_cannot_tell() {
    let d = scratch("consensus-unknown");
    // Three callers. Manta and Delly always vote with their IDs; CNVpytor
    // votes once with no ID, as for a VCF whose records have ID `.`.
    let head = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmanta\tdelly\tcnvpytor\n";
    let consensus = put(
        &d,
        "consensus.vcf",
        &format!(
            "{head}chr20\t1000\tA\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000;SUPP=2;SUPP_VEC=110\tGT:TY:ID\t0/1:DEL:M1\t0/1:DEL:D1\t./.:NaN:NaN\n\
             chr20\t9000\tB\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=12000;SUPP=2;SUPP_VEC=101\tGT:TY:ID\t0/1:DEL:M2\t./.:NaN:NaN\t0/1:DEL:.\n"
        ),
    );
    let body = [
        c45_row(".", 1000, 5000, "DEL", "manta", "GT", "0/1"),
        c45_row("", 1000, 5000, "DEL", "manta", "GT", "0/1"),
        c45_row("M1", 1000, 5000, "DEL", "sv_vcf", "GT", "0/1"),
        c45_row("CNV9", 9000, 12000, "DEL", "cnvpytor", "GT", "0/1"),
        c45_row("M9", 30000, 35000, "DEL", "manta", "GT", "0/1"),
        c45_row("D1", 1000, 5000, "DEL", "delly", "GT", "0/1"),
    ]
    .concat();
    let rows = postfilter_rows(&d, &body, Some(consensus));
    let support = |i: usize| {
        (
            rows[i]["Support"].as_str(),
            rows[i]["Support_callers"].as_str(),
        )
    };
    assert_eq!(support(0), ("", ""), "a row with ID `.`: unknown");
    assert_eq!(support(1), ("", ""), "a row with no ID: unknown");
    assert_eq!(
        support(2),
        ("", ""),
        "a caller that is no consensus column (SV_VCF's sv_vcf): unknown"
    );
    assert_eq!(
        support(3),
        ("", ""),
        "a CNVpytor row when CNVpytor voted with no ID: unknown"
    );
    assert_eq!(
        support(4),
        ("1", "manta"),
        "control: a Manta row no event holds, Manta's members all named: single caller"
    );
    assert_eq!(support(5), ("2", "manta,delly"), "control: a member");

    // A consensus whose records carry no FORMAT/ID cannot name any member.
    let no_ids = put(
        &d,
        "no_ids.vcf",
        &format!(
            "{CONSENSUS_HEAD}chr20\t1000\tA\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000;SUPP=2;SUPP_VEC=11\tGT\t0/1\t0/1\n"
        ),
    );
    let rows = postfilter_rows(
        &d,
        &c45_row("M1", 1000, 5000, "DEL", "manta", "GT", "0/1"),
        Some(no_ids),
    );
    assert_eq!(
        (
            rows[0]["Support"].as_str(),
            rows[0]["Support_callers"].as_str()
        ),
        ("", ""),
        "a consensus without FORMAT/ID: unknown, not single caller"
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

// og sv survivor-prep: which records vote in step 22's SURVIVOR merge.

use super::survivor_prep::{self, NoVote};

fn prep_text(input: &str, caller: &str) -> (String, survivor_prep::Counts) {
    let mut out = Vec::new();
    let c = survivor_prep::prep(input.as_bytes(), "in.vcf", caller, &mut out).unwrap();
    (String::from_utf8(out).unwrap(), c)
}

const SV_HEAD: &str = "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n";

fn sv_rec(id: &str, filter: &str, info: &str, format: &str, sample: &str) -> String {
    format!(
        "chr20\t1000\t{id}\tN\t<DEL>\t60\t{filter}\tSVTYPE=DEL;END=5000;{info}\t{format}\t{sample}\n"
    )
}

#[test]
fn a_reference_genotype_gives_no_vote() {
    for gt in ["0/0", "0|0", "0", "0/.", ".|0", "00"] {
        assert!(survivor_prep::reference_only(gt), "{gt} calls no ALT");
    }
    for gt in ["0/1", "1/1", "1", "./1", "0|2", "./.", ".", ""] {
        assert!(
            !survivor_prep::reference_only(gt),
            "{gt} is not reference only"
        );
    }
    let input = [
        SV_HEAD.to_owned(),
        sv_rec("REF00", "PASS", "CT=3to5", "GT:FT", "0/0:PASS"),
        sv_rec("REF0", "PASS", "CT=3to5", "GT", "0"),
        sv_rec("HALF", "PASS", "CT=3to5", "GT", "0/."),
        sv_rec("HET", "PASS", "CT=3to5", "GT:FT", "0/1:PASS"),
        sv_rec("NOCALL", ".", "CT=3to5", "GT", "./."),
    ]
    .concat();
    let (out, c) = prep_text(&input, "delly");
    // Red-first control: before this verb, step 22 kept every PASS record,
    // so REF00 made a Manta-only call read "2 callers".
    assert_eq!(ids(&out), ["HET", "NOCALL"]);
    assert_eq!((c.records, c.kept, c.reference), (5, 2, 3));
}

#[test]
fn a_sample_filter_or_filter_gives_no_vote() {
    let input = [
        SV_HEAD.to_owned(),
        sv_rec("LOWQ", "PASS", "CT=3to5", "GT:GQ:FT", "0/1:30:LowQual"),
        sv_rec("FTDOT", "PASS", "CT=3to5", "GT:GQ:FT", "0/1:30:."),
        sv_rec("SHORTCOL", "PASS", "CT=3to5", "GT:GQ:FT", "0/1"),
        sv_rec("SITELOW", "LowQual", "CT=3to5", "GT:FT", "0/1:PASS"),
    ]
    .concat();
    let (out, c) = prep_text(&input, "delly");
    assert_eq!(ids(&out), ["FTDOT", "SHORTCOL"]);
    assert_eq!((c.sample_filter, c.filter), (1, 1));
    let line = sv_rec("X", "PASS", "", "GT:FT", "1/1:LowQual");
    let v: Vec<&str> = line.trim_end().split('\t').collect();
    assert_eq!(
        survivor_prep::verdict(&v, false),
        Some(NoVote::SampleFilter)
    );
}

#[test]
fn cnvpytor_votes_by_its_quality_rule() {
    let cnv = |id: &str, p1: &str, q0: &str, gt: &str| {
        sv_rec(
            id,
            "PASS",
            &format!("pytorP1={p1};pytorQ0={q0};pytorPN=0"),
            "GT",
            gt,
        )
    };
    let input = [
        SV_HEAD.to_owned(),
        cnv("GOOD", "1e-08", "0.01", "./."),
        cnv("NOISY", "0.5", "0", "./."),
        cnv("MAPQ0", "1e-08", "0.9", "./."),
        cnv("EDGE", "0.01", "0", "./."),
        cnv("NAN", "nan", "0", "./."),
        sv_rec("NOTAGS", "PASS", "SVLEN=-4000", "GT", "./."),
    ]
    .concat();
    let (out, c) = prep_text(&input, "cnvpytor");
    // ./. is how the table-made VCF writes every call: it keeps the vote.
    assert_eq!(ids(&out), ["GOOD"]);
    assert_eq!((c.kept, c.quality), (1, 5));
    assert!(
        c.to_string().ends_with("5 below CNVpytor's quality rule"),
        "{c}"
    );
    // The rule is CNVpytor's own: another caller's records keep their vote.
    let (out, c) = prep_text(&input, "delly");
    assert_eq!(c.kept, 6, "{out}");
    assert_eq!(
        survivor_prep::verdict(&["c", "1", ".", "N", "<DEL>", ".", "PASS", "."], true),
        Some(NoVote::Quality)
    );
}

#[test]
fn a_header_without_gt_or_samples_does_not_stop_it() {
    // Sites only and no GT declared: SURVIVOR still gets a GT column.
    let sites = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr20\t1000\tS1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000\n";
    let (out, c) = prep_text(sites, "tiddit");
    assert_eq!(c.kept, 1);
    assert_eq!(
        out,
        "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttiddit\nchr20\t1000\tS1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000\tGT\t./.\n"
    );
    // A sample column whose FORMAT has no GT keeps its vote, and the
    // header gets no second GT line.
    let no_gt = "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\nchr20\t1000\tC1\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=5000\tCN\t3\t2\n";
    let (out, c) = prep_text(no_gt, "manta");
    assert_eq!(c.kept, 1);
    assert_eq!(out.matches("##FORMAT=<ID=GT").count(), 1);
    // One sample column, named after the caller.
    assert!(out.contains("\tFORMAT\tmanta\n"), "{out}");
    assert!(out.ends_with("\tCN\t3\n"), "{out}");
}

#[test]
fn survivor_prep_cli_writes_the_file_and_refuses_bad_input() {
    let d = scratch("svprep");
    let input = put(
        &d,
        "delly.vcf",
        &[
            SV_HEAD.to_owned(),
            sv_rec("REF00", "PASS", "", "GT", "0/0"),
            sv_rec("HET", "PASS", "", "GT", "0/1"),
        ]
        .concat(),
    );
    let out = d.join("delly.out.vcf").to_string_lossy().into_owned();
    let run = |a: &[&str]| {
        let argv: Vec<String> = a.iter().map(|s| (*s).to_owned()).collect();
        let (mut o, mut e) = (Vec::new(), Vec::new());
        let code = survivor_prep::cli(&argv, &mut o, &mut e);
        (
            code,
            String::from_utf8(o).unwrap(),
            String::from_utf8(e).unwrap(),
        )
    };
    let (code, _, err) = run(&["--caller", "delly", &input, &out]);
    assert_eq!(code, 0, "{err}");
    assert_eq!(ids(&std::fs::read_to_string(&out).unwrap()), ["HET"]);
    assert!(
        err.contains("og sv survivor-prep delly: 2 records, 1 vote;")
            && err.contains("1 genotyped reference"),
        "{err}"
    );
    let (code, stdout, _) = run(&["--caller", "delly", &input, "-"]);
    assert_eq!((code, ids(&stdout)), (0, vec!["HET".to_owned()]));
    let (code, _, err) = run(&[
        "--caller",
        "delly",
        &d.join("none.vcf").to_string_lossy(),
        &out,
    ]);
    assert!(code == 1 && err.starts_with("ERROR: cannot read "), "{err}");
    let short = put(&d, "short.vcf", &format!("{SV_HEAD}chr20\t1000\tX\n"));
    let bad_out = d.join("bad.vcf").to_string_lossy().into_owned();
    let (code, _, err) = run(&["--caller", "delly", &short, &bad_out]);
    assert!(
        code == 1 && err.contains("line 4: a record needs at least 8 columns"),
        "{err}"
    );
    assert!(!d.join("bad.vcf").exists() && !d.join("bad.vcf.part").exists());
    assert_eq!(run(&["--caller", "a b", &input, &out]).0, 2);
    assert_eq!(run(&["--caller", "delly", &input]).0, 2);
    assert_eq!(run(&[&input, &out]).0, 2);
    std::fs::remove_dir_all(d).unwrap();
}

// og sv consensus-check: a SURVIVOR record is one event two or more callers
// agree on. The records below are what the pinned SURVIVOR 1.0.7 wrote for
// made-up calls (`merge LIST 1000 2 1 1 0 50`, columns manta, delly,
// cnvpytor), so each red case is a merge SURVIVOR really makes.

use super::consensus_check::{self, Apart, Member};

const CC_HEAD: &str = "##fileformat=VCFv4.1\n##source=SURVIVOR\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmanta\tdelly\tcnvpytor\n";
const CC_FORMAT: &str = "GT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO";
const CC_EMPTY: &str = "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN";

/// Two deletions of 100 bp that start 900 bp apart and share no base.
const CC_DISJOINT: &str = "chr20\t1000900\td1\tN\t<DEL>\t60\tPASS\tSUPP=2;SUPP_VEC=110;SVLEN=-100;SVTYPE=DEL;SVMETHOD=SURVIVOR1.0.7;CHR2=chr20;END=1001000;CIPOS=-900,0;CIEND=-900,0;STRANDS=+-\tGT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO\t0/1:NA:100:0,0:+-:60:DEL:m1:NA:NA:chr20_1000000-chr20_1000100\t0/1:NA:100:0,0:+-:60:DEL:d1:NA:NA:chr20_1000900-chr20_1001000\t./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN";
/// A Manta deletion and a Delly duplication at the same place: SURVIVOR
/// wrote one DUP with two callers, type and strand flags set.
const CC_DEL_DUP: &str = "chr20\t7000000\td7\tN\t<DUP>\t60\tPASS\tSUPP=2;SUPP_VEC=110;SVLEN=-3000;SVTYPE=DUP;SVMETHOD=SURVIVOR1.0.7;CHR2=chr20;END=7003000;CIPOS=0,0;CIEND=0,0;STRANDS=-+\tGT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO\t0/1:NA:3000:0,0:+-:60:DEL:m7:NA:NA:chr20_7000000-chr20_7003000\t0/1:NA:3000:0,0:-+:60:DUP:d7:NA:NA:chr20_7000000-chr20_7003000\t./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN";
/// The same 500 bp deletion from all three callers.
const CC_AGREE: &str = "chr20\t2000100\td2\tN\t<DEL>\t60\tPASS\tSUPP=3;SUPP_VEC=111;SVLEN=-500;SVTYPE=DEL;SVMETHOD=SURVIVOR1.0.7;CHR2=chr20;END=2000600;CIPOS=-100,0;CIEND=-100,0;STRANDS=+-\tGT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO\t0/1:NA:500:0,0:+-:60:DEL:m2:NA:NA:chr20_2000000-chr20_2000500\t0/1:NA:500:0,0:+-:60:DEL:d2:NA:NA:chr20_2000100-chr20_2000600\t./.:NA:500:0,0:+-:.:DEL:c2:NA:NA:chr20_2000050-chr20_2000550";
/// A Manta inversion and a Delly breakend pair (BND, 3to3) of it.
const CC_BND_INV: &str = "chr20\t3000100\td3\tN\tN]chr20:3005100]\t60\tPASS\tSUPP=2;SUPP_VEC=110;SVLEN=5000;SVTYPE=INV;SVMETHOD=SURVIVOR1.0.7;CHR2=chr20;END=3005100;CIPOS=-100,0;CIEND=-100,0;STRANDS=++\tGT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO\t0/1:NA:5000:0,0:++:60:INV:m3:NA:NA:chr20_3000000-chr20_3005000\t0/1:NA:5000:0,0:++:60:INV:d3:NA:NA:chr20_3000100-chr20_3005100\t./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN";
/// One translocation, chr20 to chr21, from two callers.
const CC_TRA: &str = "chr20\t5000300\tb5\tN\tN[chr21:600000[\t60\tPASS\tSUPP=2;SUPP_VEC=110;SVLEN=0;SVTYPE=TRA;SVMETHOD=SURVIVOR1.0.7;CHR2=chr21;END=600400;CIPOS=-300,0;CIEND=-400,0;STRANDS=+-\tGT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO\t0/1:NA:4400000:0,0:+-:60:TRA:a5:NA:NA:chr20_5000000-chr21_600000\t0/1:NA:4399900:0,0:+-:60:TRA:b5:NA:NA:chr20_5000300-chr21_600400\t./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN";
/// Delly gave two calls; both overlap Manta's by half or more.
const CC_TWO_FROM_ONE: &str = "chr20\t6000050\td6b\tN\t<DEL>\t60\tPASS\tSUPP=2;SUPP_VEC=110;SVLEN=-2000;SVTYPE=DEL;SVMETHOD=SURVIVOR1.0.7;CHR2=chr20;END=6002050;CIPOS=-50,850;CIEND=-150,0;STRANDS=+-\tGT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO\t0/1:NA:2000:0,0:+-:60:DEL:m6:NA:NA:chr20_6000000-chr20_6002000\t0/1:NA:2000:0,0:+-:60,60:DEL,DEL:d6b:NA:NA:chr20_6000050-chr20_6002050,chr20_6000900-chr20_6001900\t./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN";
/// Three 600 bp deletions 300 and 400 bp apart: the first and the last
/// share no base.
const CC_CHAIN: &str = "chr20\t6000300\tb6\tN\t<DEL>\t60\tPASS\tSUPP=3;SUPP_VEC=111;SVLEN=-600;SVTYPE=DEL;SVMETHOD=SURVIVOR1.0.7;CHR2=chr20;END=6000900;CIPOS=-300,400;CIEND=-300,400;STRANDS=+-\tGT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO\t0/1:NA:600:0,0:+-:60:DEL:a6:NA:NA:chr20_6000000-chr20_6000600\t0/1:NA:600:0,0:+-:60:DEL:b6:NA:NA:chr20_6000300-chr20_6000900\t0/1:NA:600:0,0:+-:60:DEL:c6:NA:NA:chr20_6000700-chr20_6001300";

/// A Delly-shaped insertion pair (`<INS>`, END = POS + 1, SVLEN, CT=NtoN
/// in the input) starting 200 bp apart, of `a` and `b` bp, as SURVIVOR
/// merged it. A Manta insertion has no SVLEN SURVIVOR reads, so it never
/// merges; Delly's does.
fn cc_ins(a: i64, b: i64) -> String {
    let mean = ((a + b) as f64 / 2.0).round() as i64;
    format!(
        "chr20\t4000200\tb4\tN\t<INS>\t60\tPASS\tSUPP=2;SUPP_VEC=110;SVLEN={mean};SVTYPE=INS;SVMETHOD=SURVIVOR1.0.7;CHR2=chr20;END=4000201;CIPOS=-200,0;CIEND=0,0;STRANDS=++\t{CC_FORMAT}\t0/1:NA:{a}:0,0:++:60:INS:a4:NA:NA:chr20_4000000-chr20_4000001\t0/1:NA:{b}:0,0:++:60:INS:b4:NA:NA:chr20_4000200-chr20_4000201\t{CC_EMPTY}"
    )
}

fn cc_text(records: &[&str]) -> (String, consensus_check::Counts) {
    let input = format!("{CC_HEAD}{}\n", records.join("\n"));
    let mut out = Vec::new();
    let c = consensus_check::check(input.as_bytes(), "merged.vcf", &mut out).unwrap();
    (String::from_utf8(out).unwrap(), c)
}

fn cc_records(out: &str) -> Vec<&str> {
    out.lines().filter(|l| !l.starts_with('#')).collect()
}

fn member(kind: &str, strands: &str, start: i64, end: i64, len: i64) -> Member {
    Member {
        caller: 0,
        index: 0,
        kind: kind.to_owned(),
        strands: strands.to_owned(),
        chrom2: "chr20".to_owned(),
        start,
        end,
        len,
        other_contig: false,
    }
}

#[test]
fn two_deletions_900_bp_apart_are_not_one_event() {
    let (out, c) = cc_text(&[CC_DISJOINT]);
    // Red-first: SURVIVOR's record reads SUPP=2 for two events of one caller each.
    assert!(cc_records(&out).is_empty(), "{out}");
    assert_eq!((c.records, c.split, c.dropped, c.written()), (1, 1, 1, 0));
    let a = member("DEL", "+-", 1_000_000, 1_000_100, 100);
    let b = member("DEL", "+-", 1_000_900, 1_001_000, 100);
    assert_eq!(consensus_check::apart(&a, &b), Some(Apart::Overlap));
    // Half the longer call is enough; a base less is not.
    let c50 = member("DEL", "+-", 1_000_050, 1_000_150, 100);
    assert_eq!(consensus_check::apart(&a, &c50), None);
    let c49 = member("DEL", "+-", 1_000_051, 1_000_151, 100);
    assert_eq!(consensus_check::apart(&a, &c49), Some(Apart::Overlap));
}

#[test]
fn a_deletion_and_a_duplication_are_not_one_event() {
    let (out, c) = cc_text(&[CC_DEL_DUP]);
    assert!(cc_records(&out).is_empty(), "{out}");
    assert_eq!((c.split, c.dropped), (1, 1));
    let del = member("DEL", "+-", 7_000_000, 7_003_000, 3000);
    let dup = member("DUP", "-+", 7_000_000, 7_003_000, 3000);
    assert_eq!(consensus_check::apart(&del, &dup), Some(Apart::Type));
    let del_other_strands = member("DEL", "-+", 7_000_000, 7_003_000, 3000);
    assert_eq!(
        consensus_check::apart(&del, &del_other_strands),
        Some(Apart::Strands)
    );
    // DUP:TANDEM is a DUP.
    let tandem = member("DUP:TANDEM", "-+", 7_000_000, 7_003_000, 3000);
    assert_eq!(consensus_check::apart(&dup, &tandem), None);
}

#[test]
fn delly_shaped_insertions_need_a_like_size() {
    // 1,000 bp and 300 bp: SURVIVOR merged them (ends within 1 kb).
    let unlike = cc_ins(1000, 300);
    let like = cc_ins(1000, 600);
    let (out, c) = cc_text(&[&unlike, &like]);
    let recs = cc_records(&out);
    assert_eq!(recs, [like.as_str()], "{out}");
    assert_eq!((c.records, c.kept, c.split, c.dropped), (2, 1, 1, 1));
    let a = member("INS", "++", 4_000_000, 4_000_001, 1000);
    let half = member("INS", "++", 4_000_200, 4_000_201, 500);
    assert_eq!(consensus_check::apart(&a, &half), None);
    let under = member("INS", "++", 4_000_200, 4_000_201, 499);
    assert_eq!(consensus_check::apart(&a, &under), Some(Apart::Size));
    let far = member("INS", "++", 4_001_001, 4_001_002, 1000);
    assert_eq!(consensus_check::apart(&a, &far), Some(Apart::Distance));
    let near = member("INS", "++", 4_001_000, 4_001_001, 1000);
    assert_eq!(consensus_check::apart(&a, &near), None);
}

#[test]
fn a_breakend_pair_of_an_inversion_and_a_translocation_stay() {
    // SURVIVOR's own rule counts BND, INV and TRA as one type.
    let (out, c) = cc_text(&[CC_BND_INV, CC_TRA]);
    assert_eq!(cc_records(&out), [CC_BND_INV, CC_TRA]);
    assert_eq!((c.kept, c.split), (2, 0));
    let inv = member("INV", "++", 3_000_000, 3_005_000, 5000);
    let bnd = member("BND", "++", 3_000_100, 3_005_100, 5000);
    assert_eq!(consensus_check::apart(&inv, &bnd), None);
    let tra = |start: i64, chrom2: &str, end: i64| Member {
        chrom2: chrom2.to_owned(),
        other_contig: true,
        len: 0,
        ..member("TRA", "+-", start, end, 0)
    };
    assert_eq!(
        consensus_check::apart(
            &tra(5_000_000, "chr21", 600_000),
            &tra(5_000_300, "chr21", 601_000)
        ),
        None
    );
    assert_eq!(
        consensus_check::apart(
            &tra(5_000_000, "chr21", 600_000),
            &tra(5_000_300, "chr21", 601_001)
        ),
        Some(Apart::Distance)
    );
    assert_eq!(
        consensus_check::apart(
            &tra(5_000_000, "chr21", 600_000),
            &tra(5_000_000, "chr22", 600_000)
        ),
        Some(Apart::Distance)
    );
}

#[test]
fn a_merge_of_one_event_is_written_as_survivor_wrote_it() {
    let (out, c) = cc_text(&[CC_AGREE, CC_TWO_FROM_ONE]);
    assert_eq!(cc_records(&out), [CC_AGREE, CC_TWO_FROM_ONE]);
    assert_eq!((c.records, c.kept, c.written()), (2, 2, 2));
    // The header gains one line, just above #CHROM.
    let head: Vec<&str> = out.lines().filter(|l| l.starts_with('#')).collect();
    assert_eq!(head.len(), 4, "{out}");
    assert_eq!(head[2], consensus_check::HEADER);
    assert!(head[3].starts_with("#CHROM\t"));
    assert_eq!(
        c.to_string(),
        "2 records, 2 written; 2 one event as merged, 0 not one event: split into 0 records of 2+ callers, 0 dropped with fewer than 2 callers left"
    );
}

#[test]
fn a_chain_keeps_the_first_event_and_recomputes_its_coordinates() {
    let (out, c) = cc_text(&[CC_CHAIN]);
    // a6 and b6 overlap by half; c6 overlaps b6 but shares no base with a6,
    // so it is an event of its own, of one caller, and is dropped.
    let want = format!(
        "chr20\t6000000\ta6\tN\t<DEL>\t60\tPASS\tSUPP=2;SUPP_VEC=110;SVLEN=-600;SVTYPE=DEL;SVMETHOD=SURVIVOR1.0.7;CHR2=chr20;END=6000600;CIPOS=0,300;CIEND=0,300;STRANDS=+-\t{CC_FORMAT}\t0/1:NA:600:0,0:+-:60:DEL:a6:NA:NA:chr20_6000000-chr20_6000600\t0/1:NA:600:0,0:+-:60:DEL:b6:NA:NA:chr20_6000300-chr20_6000900\t{CC_EMPTY}"
    );
    assert_eq!(cc_records(&out), [want.as_str()]);
    assert_eq!((c.split, c.from_split, c.dropped), (1, 1, 0));
}

#[test]
fn a_caller_keeps_only_its_calls_of_the_event() {
    // Delly's second call (100 bp inside Manta's 2 kb) is not the event.
    let rec = CC_TWO_FROM_ONE.replace("chr20_6000900-chr20_6001900", "chr20_6000900-chr20_6001000");
    let (out, c) = cc_text(&[&rec]);
    let want = format!(
        "chr20\t6000000\tm6\tN\t<DEL>\t60\tPASS\tSUPP=2;SUPP_VEC=110;SVLEN=-2000;SVTYPE=DEL;SVMETHOD=SURVIVOR1.0.7;CHR2=chr20;END=6002000;CIPOS=0,50;CIEND=0,50;STRANDS=+-\t{CC_FORMAT}\t0/1:NA:2000:0,0:+-:60:DEL:m6:NA:NA:chr20_6000000-chr20_6002000\t0/1:NA:2000:0,0:+-:60:DEL:.:NA:NA:chr20_6000050-chr20_6002050\t{CC_EMPTY}"
    );
    // Delly's ID named its second call, so the kept column says "." there.
    assert_eq!(cc_records(&out), [want.as_str()]);
    assert_eq!((c.split, c.from_split), (1, 1));
}

#[test]
fn two_events_in_one_record_become_two_records() {
    let rec = format!(
        "chr20\t1000\tdB\tN\t<DEL>\t60\tPASS\tSUPP=2;SUPP_VEC=110;SVLEN=-1000;SVTYPE=DEL;SVMETHOD=SURVIVOR1.0.7;CHR2=chr20;END=2000;CIPOS=0,4000;CIEND=0,4000;STRANDS=+-\t{CC_FORMAT}\t0/1:NA:1000:0,0:+-:60,60:DEL,DEL:mB:NA:NA:chr20_1000-chr20_2000,chr20_5000-chr20_6000\t0/1:NA:1000:0,0:+-:60,60:DEL,DEL:dB:NA:NA:chr20_1100-chr20_2100,chr20_5100-chr20_6100\t{CC_EMPTY}"
    );
    let (out, c) = cc_text(&[&rec]);
    let first = format!(
        "chr20\t1000\t.\tN\t<DEL>\t60\tPASS\tSUPP=2;SUPP_VEC=110;SVLEN=-1000;SVTYPE=DEL;SVMETHOD=SURVIVOR1.0.7;CHR2=chr20;END=2000;CIPOS=0,100;CIEND=0,100;STRANDS=+-\t{CC_FORMAT}\t0/1:NA:1000:0,0:+-:60:DEL:.:NA:NA:chr20_1000-chr20_2000\t0/1:NA:1000:0,0:+-:60:DEL:.:NA:NA:chr20_1100-chr20_2100\t{CC_EMPTY}"
    );
    let second = format!(
        "chr20\t5000\tmB\tN\t<DEL>\t60\tPASS\tSUPP=2;SUPP_VEC=110;SVLEN=-1000;SVTYPE=DEL;SVMETHOD=SURVIVOR1.0.7;CHR2=chr20;END=6000;CIPOS=0,100;CIEND=0,100;STRANDS=+-\t{CC_FORMAT}\t0/1:NA:1000:0,0:+-:60:DEL:mB:NA:NA:chr20_5000-chr20_6000\t0/1:NA:1000:0,0:+-:60:DEL:dB:NA:NA:chr20_5100-chr20_6100\t{CC_EMPTY}"
    );
    assert_eq!(cc_records(&out), [first.as_str(), second.as_str()]);
    assert_eq!((c.split, c.from_split, c.written()), (1, 2, 2));
}

#[test]
fn a_record_it_cannot_read_is_written_as_it_was() {
    let no_co = CC_DISJOINT.replace(":AAL:CO\t", ":AAL\t");
    let short_vec = CC_DISJOINT.replace("SUPP_VEC=110", "SUPP_VEC=11");
    let other_contig =
        CC_DISJOINT.replace("chr20_1000000-chr20_1000100", "chr2_1000000-chr2_1000100");
    let (out, c) = cc_text(&[&no_co, &short_vec, &other_contig]);
    assert_eq!(
        cc_records(&out),
        [no_co.as_str(), short_vec.as_str(), other_contig.as_str()]
    );
    assert_eq!((c.unread, c.written()), (3, 3));
    assert!(
        c.to_string()
            .ends_with("; 3 could not be read and are written as SURVIVOR wrote them")
    );
}

#[test]
fn consensus_check_cli_writes_the_file_and_refuses_bad_input() {
    let d = scratch("svcc");
    let input = put(
        &d,
        "merged.vcf",
        &format!("{CC_HEAD}{CC_DISJOINT}\n{CC_AGREE}\n"),
    );
    let out = d.join("checked.vcf").to_string_lossy().into_owned();
    let run = |a: &[&str]| {
        let argv: Vec<String> = a.iter().map(|s| (*s).to_owned()).collect();
        let (mut o, mut e) = (Vec::new(), Vec::new());
        let code = consensus_check::cli(&argv, &mut o, &mut e);
        (
            code,
            String::from_utf8(o).unwrap(),
            String::from_utf8(e).unwrap(),
        )
    };
    let (code, _, err) = run(&[&input, &out]);
    assert_eq!(code, 0, "{err}");
    assert_eq!(ids(&std::fs::read_to_string(&out).unwrap()), ["d2"]);
    assert_eq!(
        err,
        "og sv consensus-check: 2 records, 1 written; 1 one event as merged, 1 not one event: split into 0 records of 2+ callers, 1 dropped with fewer than 2 callers left\n"
    );
    let (code, stdout, _) = run(&[&input, "-"]);
    assert_eq!((code, ids(&stdout)), (0, vec!["d2".to_owned()]));
    let (code, _, err) = run(&[&d.join("none.vcf").to_string_lossy(), &out]);
    assert!(code == 1 && err.starts_with("ERROR: cannot read "), "{err}");
    let (code, stdout, _) = run(&["--help"]);
    assert!(code == 0 && stdout.starts_with("usage: og sv consensus-check"));
    assert_eq!(run(&[&input]).0, 2);
    assert_eq!(run(&["--og-scaffold-probe"]).0, 2);
    std::fs::remove_dir_all(d).unwrap();
}
