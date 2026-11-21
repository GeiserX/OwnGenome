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
