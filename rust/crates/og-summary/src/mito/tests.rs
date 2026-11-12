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
    let (head, rows) = annotate(&vcf, Some(&d.path("cv.vcf.gz")), Some(&gn)).unwrap();
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
