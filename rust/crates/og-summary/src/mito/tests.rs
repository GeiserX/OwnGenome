use super::*;
use crate::testutil::{Dir, run};

fn cli_of(verb: &'static str) -> impl FnOnce(&[String], &mut Vec<u8>, &mut Vec<u8>) -> u8 {
    move |a, o, e| cli(Some(verb), a, o, e)
}

#[test]
fn min_allele_trims_the_shared_bases() {
    let k = |p, r: &str, a: &str| min_allele(p, r, a);
    assert_eq!(k(302, "AC", "A"), (302, "AC".into(), "A".into()));
    assert_eq!(k(302, "AC", "ACC"), (302, "A".into(), "AC".into()));
    assert_eq!(k(10, "GAT", "GCT"), (11, "A".into(), "C".into()));
    assert_eq!(k(10, "A", "G"), (10, "A".into(), "G".into()));
    assert_eq!(k(10, "AA", "AA"), (10, "A".into(), "A".into()));
}

#[test]
fn left_allele_moves_an_indel_to_the_left_end_of_its_repeat() {
    let k = |p, r: &str, a: &str| left_allele(p, r, a);
    let al = |p, r: &str, a: &str| (p, r.to_owned(), a.to_owned());
    // One more C in m.311-315, spelt three ways.
    assert_eq!(k(315, "C", "CC"), al(310, "T", "TC"));
    assert_eq!(k(312, "CC", "CCC"), al(310, "T", "TC"));
    assert_eq!(k(310, "T", "TC"), al(310, "T", "TC"));
    // One C fewer in m.303-309, and one more in m.16190-16193.
    assert_eq!(k(308, "CC", "C"), al(302, "AC", "A"));
    assert_eq!(k(16193, "C", "CC"), al(16189, "T", "TC"));
    // One CA fewer in the m.514-523 CA repeat.
    assert_eq!(k(521, "ACA", "A"), al(513, "GCA", "G"));
    // An SNV, a REF that is not the rCRS, and nothing before m.1 stay.
    assert_eq!(k(3243, "A", "G"), al(3243, "A", "G"));
    assert_eq!(k(315, "G", "GG"), al(315, "G", "GG"));
    assert_eq!(k(1, "G", "GG"), al(1, "G", "GG"));
    assert_eq!(k(16570, "A", "AA"), al(16570, "A", "AA"));
}

fn idxstats(chrm: u64, per_chrom: u64, m: &str) -> String {
    let mut s: String = (1..=22)
        .map(|i| format!("chr{i}\t50000000\t{per_chrom}\t0\n"))
        .collect();
    s.push_str(&format!(
        "chrX\t156040895\t{per_chrom}\t0\n{m}\t16569\t{chrm}\t0\n*\t0\t0\t500\n"
    ));
    s
}

#[test]
fn copies_from_idxstats() {
    // 200 reads per kb of autosome is about 30x; chrM at 100 times that depth.
    let per = 200 * 50_000;
    let chrm = (100.0 * 16569.0 * per as f64 / 50_000_000.0).round() as u64;
    let c = copies(&idxstats(chrm, per, "chrM")).unwrap();
    let (ratio, n) = c.call.unwrap();
    assert!((n - 200.0).abs() <= 10.0 && ratio == 100.0, "{c:?}");
    let half = copies(&idxstats(chrm / 2, per, "chrM"))
        .unwrap()
        .call
        .unwrap()
        .1;
    assert!((half - 100.0).abs() <= 5.0, "{half}");
    assert_eq!(
        copies(&idxstats(chrm, 500, "chrM")).unwrap().call,
        Err(CnReason::LowAutosomalDepth)
    );
    assert_eq!(
        copies(&idxstats(0, per, "chrM")).unwrap().call,
        Err(CnReason::NoChrmReads)
    );
    let plain = idxstats(chrm, per, "chrM")
        .replace("chrM", "MT")
        .replace("chr", "");
    assert!(copies(&plain).unwrap().call.is_ok());
    // Python knew only chrM and MT; chrMT and M are chrM too now.
    for m in ["chrMT", "M"] {
        assert!(copies(&idxstats(chrm, per, m)).unwrap().call.is_ok(), "{m}");
    }
    assert_eq!(
        copies("chr1\t10\tx\t0\n").unwrap_err(),
        "invalid literal for int() with base 10: 'x'"
    );
}

#[test]
fn copies_table_and_line() {
    let d = Dir::new("mito-copies");
    let ip = d.put("s.idxstats", idxstats(331_380, 10_000_000, "chrM"));
    let out = d.path("cn.tsv");
    let (code, stdout, _) = run(cli_of("copies"), &["--idxstats", &ip, "--out", &out]);
    assert_eq!(code, 0);
    assert_eq!(
        d.read("cn.tsv"),
        "metric\tvalue\nstatus\tcalled\nreason\t.\nchrm_reads\t331380\nchrm_length\t16569\n\
         autosomal_reads\t220000000\nautosomal_length\t1100000000\ndepth_ratio\t100.0\nmtdna_copies\t200.0\n"
    );
    assert_eq!(stdout, "  mtDNA copy number: 200.0\n");
    let ip = d.put("t.idxstats", "chr1\t1000\t5\t0\n");
    let (code, stdout, _) = run(cli_of("copies"), &["--idxstats", &ip, "--out", &out]);
    assert_eq!(code, 0);
    assert_eq!(
        stdout,
        "  mtDNA copy number: not assessed (no_chrM_reads)\n"
    );
    assert_eq!(
        d.read("cn.tsv"),
        "metric\tvalue\nstatus\tnot assessed\nreason\tno_chrM_reads\nchrm_reads\t.\nchrm_length\t.\n\
         autosomal_reads\t5\nautosomal_length\t1000\ndepth_ratio\t.\nmtdna_copies\t.\n"
    );
}

const EK_HEAD: &str = "\"Title\";\"5' breakpoint\";\"3' breakpoint\";\"Freq\";\"Freq For\";\"Freq Rev\";\"5' Blast\";\
\"3' Blast\";\"5' Depth\";\"3' Depth\";\"Repetition\"\n";

fn eklipse(rows: &[(i64, i64, &str)]) -> String {
    let mut s = EK_HEAD.to_owned();
    for (b5, b3, v) in rows {
        s.push_str(&format!("\"S\";\"{b5}\";\"{b3}\";\"{v}\";\"{v}\";\"{v}\";\"12\";\"14\";\"390.0\";\"402.0\";\"\"\n"));
    }
    s
}

#[test]
fn deletions_from_eklipse() {
    let d = Dir::new("mito-del");
    let p = d.put(
        "e.csv",
        eklipse(&[
            (8469, 13447, "10,2"),
            (16000, 500, "2,0"),
            (5000, 7000, "0,4567"),
        ]),
    );
    let rows = deletions(&p).unwrap();
    assert_eq!(
        rows.iter().map(|r| r.level_pct).collect::<Vec<_>>(),
        [10.2, 2.0, 0.457]
    );
    assert_eq!(rows[0].length, 4977);
    assert_eq!(rows[1].length, 1068, "across the origin");
    assert_eq!((rows[0].depth_5, rows[0].repeat.as_str()), (390, "none"));
    let out = d.path("del.tsv");
    let (code, stdout, _) = run(cli_of("deletions"), &["--eklipse", &p, "--out", &out]);
    assert_eq!(
        (code, stdout.as_str()),
        (0, "  3 large deletions (eKLIPse)\n")
    );
    assert_eq!(
        d.read("del.tsv").lines().nth(1).unwrap(),
        "8469\t13447\t4977\t10.2\t10.2\t10.2\t12\t14\t390\t402\tnone"
    );
    // CRLF line endings read the same.
    let crlf = d.put(
        "crlf.csv",
        eklipse(&[(8469, 13447, "10,2")]).replace('\n', "\r\n"),
    );
    assert_eq!(deletions(&crlf).unwrap()[0], rows[0]);
    let empty = d.put("none.csv", EK_HEAD);
    assert!(deletions(&empty).unwrap().is_empty());
}

#[test]
fn deletions_refuse_what_is_not_eklipse() {
    let d = Dir::new("mito-del-bad");
    let bad = d.put("bad.csv", "pos\tref\n");
    assert_eq!(
        deletions(&bad).unwrap_err(),
        format!("{bad}: not eKLIPse's deletions table (header ['pos\\tref'])")
    );
    let short = d.put("short.csv", format!("{EK_HEAD}\"S\";\"1\";\"2\"\n"));
    assert!(
        deletions(&short)
            .unwrap_err()
            .ends_with("a row with 3 fields, want 11: S;1;2")
    );
    let inf = d.put("inf.csv", eklipse(&[(1, 2, "inf")]));
    assert_eq!(deletions(&inf).unwrap_err(), "'inf' is not a number");
    let out = d.path("o.tsv");
    let (code, _, err) = run(
        cli_of("deletions"),
        &["--eklipse", &d.path("missing.csv"), "--out", &out],
    );
    assert_eq!(code, 1);
    assert!(
        err.starts_with("ERROR: [Errno 2] No such file or directory: '"),
        "{err}"
    );
}

const MITO_VCF: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
chrM\t73\t.\tA\tG\t.\tPASS\t.\tGT:AF\t1:1.00\n\
chrM\t302\t.\tAC\tA,ACC\t.\tPASS\t.\tGT:AF\t0/1/2:0.3,0.4\n\
chrM\t3243\t.\tA\tG\t.\tPASS\t.\tGT:AF\t0/1:0.3\n\
chrM\t16189\t.\tT\tC\t.\tPASS\t.\tGT:AF\t0/1:0.5\n";

const GNOMAD: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chrM\t302\t.\tA\tAC\t.\tPASS\tAN=56434;AF_hom=0.01;AF_het=0.2;max_hl=0.99\n\
chrM\t3243\t.\tA\tG\t.\tPASS\tAN=56434;AF_hom=0.00002;AF_het=0.0003;max_hl=0.9\n";

#[test]
fn annotate_without_clinvar() {
    let d = Dir::new("mito-annot");
    let vcf = d.put_gz("S_chrM_filtered.vcf.gz", MITO_VCF);
    let gn = d.put("gnomad.vcf", GNOMAD);
    let (head, rows) = annotate(&vcf, Some(&d.path("cv.vcf.gz")), Some(&gn), false).unwrap();
    assert_eq!(
        head,
        [
            format!(
                "# clinvar\tnot annotated\tno {} with its .tbi",
                d.path("cv.vcf.gz")
            ),
            "# gnomad\tgnomad.vcf\t-".to_owned()
        ]
    );
    let keys: Vec<(&str, &str, &str, &str)> = rows
        .iter()
        .map(|r| (r[1].as_str(), r[2].as_str(), r[3].as_str(), r[5].as_str()))
        .collect();
    assert_eq!(
        keys,
        [
            ("73", "A", "G", "1.00"),
            ("302", "AC", "A", "0.3"),
            ("302", "A", "AC", "0.4"),
            ("3243", "A", "G", "0.3"),
            ("16189", "T", "C", "0.5")
        ]
    );
    assert_eq!(rows[2][10..], ["0.01", "0.2", "0.99"]);
    assert_eq!(rows[3][6..10], [".", ".", ".", "."]);
    let out = d.path("a.tsv");
    let (code, stdout, _) = run(
        cli_of("annotate"),
        &["--vcf", &vcf, "--gnomad", &gn, "--out", &out],
    );
    assert_eq!(code, 0);
    assert_eq!(
        stdout,
        "  clinvar: not annotated (no ClinVar file given)\n  gnomad: gnomad.vcf\n  \
         5 chrM alleles: 0 with a ClinVar record, 2 in gnomAD\n"
    );
    assert!(d.read("a.tsv").starts_with(
        "# clinvar\tnot annotated\tno ClinVar file given\n# gnomad\tgnomad.vcf\t-\nchrom\tpos\tref\talt\tfilter\taf\t"
    ));
}

#[test]
fn verbs_and_usage() {
    let (code, _, err) = run(|a, o, e| cli(Some("nope"), a, o, e), &[]);
    assert_eq!(code, 2);
    assert!(err.contains("annotate, deletions or copies"), "{err}");
    let (code, _, err) = run(cli_of("copies"), &["--out", "x"]);
    assert_eq!(code, 2);
    assert!(err.contains("required: --idxstats"), "{err}");
}

/// Writes `text` as a BGZF VCF with its `.tbi` and returns its path.
fn put_indexed(d: &Dir, name: &str, text: &str) -> String {
    use std::fs::File;
    use std::num::NonZero;
    let p = d.path(name);
    let mut w = og_io::bgzf::Writer::new(File::create(&p).unwrap(), NonZero::<usize>::MIN).unwrap();
    w.write_all(text.as_bytes()).unwrap();
    w.finish().unwrap();
    og_io::tabix_write::build(File::open(&p).unwrap())
        .unwrap()
        .write_tbi(std::path::Path::new(&format!("{p}.tbi")), true)
        .unwrap();
    p
}

/// Indels in the D-loop poly-C tracts, each spelled at the other end of its
/// repeat from the known record: the rCRS has A at m.302, C at m.303-309,
/// T at m.310, C at m.311-315 and G at m.316.
const SHIFTED_VCF: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
chrM\t302\t.\tAC\tACC\t.\tPASS\t.\tGT:AF\t0/1:0.2\n\
chrM\t308\tdel\tCC\tC\t.\tPASS\t.\tGT:AF\t0/1:0.3\n\
chrM\t315\tins\tC\tCC\t.\tPASS\t.\tGT:AF\t0/1:0.4\n\
chrM\t16193\t.\tC\tCC\t.\tPASS\t.\tGT:AF\t0/1:0.5\n";

/// ClinVar and gnomAD spell them leftmost; ClinVar also spells the m.303
/// insertion a base to the right of its exact match at m.302.
const SHIFTED_CLINVAR: &str = "##fileformat=VCFv4.2\n##fileDate=2026-11-01\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chrM\t302\t8\tA\tAC\t.\t.\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter;CLNDN=exact\n\
chrM\t302\t11\tAC\tA\t.\t.\tCLNSIG=Benign;CLNREVSTAT=no_assertion_criteria_provided;CLNDN=deletion\n\
chrM\t303\t10\tC\tCC\t.\t.\tCLNSIG=Likely_benign;CLNDN=shifted\n\
chrM\t310\t12\tT\tTC\t.\t.\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter;CLNDN=poly-C\n";

const SHIFTED_GNOMAD: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chrM\t302\t.\tAC\tA\t.\tPASS\tAF_hom=0.03;AF_het=0.1;max_hl=0.5\n\
chrM\t310\t.\tT\tTC\t.\tPASS\tAF_hom=0.4;AF_het=0.3;max_hl=0.99\n\
chrM\t16189\t.\tT\tTC\t.\tPASS\tAF_hom=0.05;AF_het=0.2;max_hl=0.98\n";

#[test]
fn annotate_matches_an_indel_anywhere_in_its_repeat() {
    let d = Dir::new("mito-shifted");
    let vcf = d.put_gz("S_chrM_filtered.vcf.gz", SHIFTED_VCF);
    let cv = put_indexed(&d, "clinvar_chr.vcf.gz", SHIFTED_CLINVAR);
    let gn = d.put("gnomad.vcf", SHIFTED_GNOMAD);
    let out = d.path("a.tsv");
    let (code, stdout, err) = run(
        cli_of("annotate"),
        &[
            "--vcf",
            &vcf,
            "--clinvar",
            &cv,
            "--gnomad",
            &gn,
            "--out",
            &out,
        ],
    );
    assert_eq!(code, 0, "{err}");
    let rows: Vec<Vec<String>> = d
        .read("a.tsv")
        .lines()
        .filter(|l| !l.starts_with('#') && !l.starts_with("chrom\t"))
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect();
    let cols = |r: &Vec<String>| -> Vec<String> {
        [1, 2, 3, 6, 9, 10].iter().map(|&k| r[k].clone()).collect()
    };
    let want: [[&str; 6]; 4] = [
        // An exact spelling wins over a shifted one (m.303 C>CC).
        ["302", "A", "AC", "8", "exact", "."],
        // m.308 CC>C is m.302 AC>A, the deletion of one C of m.303-309.
        ["308", "CC", "C", "11", "deletion", "0.03"],
        // m.315 C>CC is m.310 T>TC, the m.310 poly-C insertion.
        ["315", "C", "CC", "12", "poly-C", "0.4"],
        // m.16193 C>CC is m.16189 T>TC, in the m.16184-16193 tract.
        ["16193", "C", "CC", ".", ".", "0.05"],
    ];
    assert_eq!(rows.iter().map(cols).collect::<Vec<_>>(), want, "{stdout}");
    assert!(
        stdout.ends_with("  4 chrM alleles: 3 with a ClinVar record, 3 in gnomAD\n"),
        "{stdout}"
    );
}

#[test]
fn annotate_keeps_an_allele_off_the_rcrs_as_it_is() {
    // A REF that is not the rCRS there (another chrM) is not moved: m.315
    // is C on the rCRS, so G>GG stays where it is. Moved over the rCRS it
    // would read m.314 C>CG and meet that record.
    let d = Dir::new("mito-off-rcrs");
    let vcf = d.put_gz(
        "S_chrM_filtered.vcf.gz",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
         chrM\t315\t.\tG\tGG\t.\tPASS\t.\tGT:AF\t0/1:0.4\n",
    );
    let gn = d.put(
        "gnomad.vcf",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         chrM\t314\t.\tC\tCG\t.\tPASS\tAF_hom=0.7;AF_het=0.1;max_hl=0.9\n",
    );
    let out = d.path("a.tsv");
    let (code, _, err) = run(
        cli_of("annotate"),
        &["--vcf", &vcf, "--gnomad", &gn, "--out", &out],
    );
    assert_eq!(code, 0, "{err}");
    let row = d.read("a.tsv").lines().last().unwrap().to_owned();
    assert!(row.starts_with("chrM\t315\tG\tGG\tPASS\t0.4\t"), "{row}");
    assert!(row.ends_with("\t.\t.\t."), "{row}");
}

#[test]
fn annotate_says_when_the_numt_filter_had_no_depth() {
    let d = Dir::new("mito-numt");
    let vcf = d.put_gz("S_chrM_filtered.vcf.gz", MITO_VCF);
    let out = d.path("a.tsv");
    let (code, stdout, err) = run(
        cli_of("annotate"),
        &["--vcf", &vcf, "--numt-default", "--out", &out],
    );
    assert_eq!(code, 0, "{err}");
    let head: Vec<String> = d
        .read("a.tsv")
        .lines()
        .take_while(|l| l.starts_with('#'))
        .map(str::to_owned)
        .collect();
    assert_eq!(
        head[2],
        "# numt\tNuMT filter used a default depth\tno median autosomal depth was given, \
         so it ran at depth 0 and marked no call as a possible NuMT"
    );
    assert!(
        stdout.contains("\n  numt: NuMT filter used a default depth (no median autosomal depth"),
        "{stdout}"
    );
    // Without the flag the table has no numt line, as before.
    let (code, _, _) = run(cli_of("annotate"), &["--vcf", &vcf, "--out", &out]);
    assert_eq!(code, 0);
    assert!(!d.read("a.tsv").contains("# numt"));
}
