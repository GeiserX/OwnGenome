//! `og family` (step 38 and the Nextflow FAMILY workflow) on synthetic
//! pedigrees, somalier tables and joint VCFs. These are the checks of the
//! retired `tests/test_family_report.py`, one for one:
//!
//! 1. prepare: a trio with a grandparent named but not sequenced writes
//!    members.tsv, somalier's pedigree (sequenced members only) and
//!    WhatsHap's trios; a pedigree problem stops it with exit 2 and names
//!    the problem (no header, an unknown column, a female father, a person
//!    who is their own parent, a loop, a second row, a sex the samplesheet
//!    contradicts, a gVCF whose sample is another, fewer than two
//!    sequenced samples);
//! 2. (the expected relatedness of each relation is og-pedigree's test
//!    `expected_relatedness`);
//! 3. check: somalier numbers that fit the trio pass; the same numbers with
//!    two names swapped in the pedigree exit 1 and name the pair; too few
//!    sites is not_checked, not a pass and not a stop; a first-degree pair
//!    with many opposite homozygotes is not parent and child;
//! 4. report: the allele from each parent where the genotypes tell, an
//!    all-heterozygous site left ambiguous, a planted new allele listed as
//!    a de novo candidate (and an opposite-homozygote site not), a son's
//!    chrX outside the PARs taken from the mother only, inside a PAR from
//!    both, the allele of an unsequenced father inferred from a duo,
//!    compound heterozygotes in trans and in cis, WhatsHap's phase compared
//!    with the genotype origin, and the same bytes on a second run; an
//!    unsorted joint VCF stops it with exit 2 and writes nothing;
//! 5. a daughter's X, the families file (`--families`, `--family`,
//!    `--which`, `--write-pedigree`, read back by `og live family state`)
//!    and `og family plan`.
//!
//! `og family` prints its messages itself, so each call runs in a child
//! process: this test binary, with [`CHILD`] set (see [`child_entry`]).

use std::collections::BTreeMap;
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use pycompat::json::{self, Value};

const CHILD: &str = "OG_TEST_FAMILY_REPORT_CHILD";
const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

/// In the child process only: runs `og family` with the arguments in
/// [`CHILD`] (separated by U+001F) and exits with its code. In a normal
/// test run it does nothing.
#[test]
fn child_entry() {
    let Ok(joined) = std::env::var(CHILD) else {
        return;
    };
    let args: Vec<String> = joined.split('\u{1f}').map(str::to_owned).collect();
    let code = og_steps::family::cli(&args);
    let _ = std::io::stdout().flush();
    std::process::exit(i32::from(code));
}

/// `og family ARGS` in a child process: (exit code, stdout then stderr).
fn run(args: &[&str]) -> (i32, String) {
    let res = Command::new(std::env::current_exe().unwrap())
        .args([
            "child_entry",
            "--exact",
            "--nocapture",
            "--quiet",
            "--test-threads=1",
        ])
        .env(CHILD, args.join("\u{1f}"))
        .output()
        .unwrap();
    // The test harness's own lines (quiet: only these) come before the
    // child's output, which exits before the harness writes more.
    let stdout = String::from_utf8_lossy(&res.stdout);
    let stdout = stdout
        .strip_prefix("\nrunning 1 test\n")
        .or_else(|| stdout.strip_prefix("running 1 test\n"))
        .unwrap_or(&stdout);
    (
        res.status.code().unwrap_or(-1),
        format!("{stdout}{}", String::from_utf8_lossy(&res.stderr)),
    )
}

fn owned(v: &[String]) -> Vec<&str> {
    v.iter().map(String::as_str).collect()
}

struct Tmp(PathBuf);

impl Tmp {
    fn new(name: &str) -> Tmp {
        let d =
            std::env::temp_dir().join(format!("og-family-report-{}-{name}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        Tmp(d)
    }
    fn j(&self, rel: &str) -> String {
        self.0.join(rel).display().to_string()
    }
}

impl Drop for Tmp {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn write(path: &str, text: &str) -> String {
    if let Some(parent) = Path::new(path).parent() {
        fs::create_dir_all(parent).unwrap();
    }
    fs::write(path, text).unwrap();
    path.to_owned()
}

fn read(path: &str) -> String {
    fs::read_to_string(path).unwrap_or_default()
}

fn gz(path: &str, text: &str) {
    let mut e = GzEncoder::new(fs::File::create(path).unwrap(), Compression::default());
    e.write_all(text.as_bytes()).unwrap();
    e.finish().unwrap();
}

fn gvcf(d: &Tmp, name: &str, sample: &str) -> String {
    let path = d.j(name);
    gz(
        &path,
        &format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n\
             chr1\t1\t.\tA\t<*>\t0\t.\tEND=10\tGT:GQ\t0/0:50\n"
        ),
    );
    path
}

type Row = BTreeMap<String, String>;

fn tsv(path: &str) -> Vec<Row> {
    let text = read(path);
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    lines
        .map(|l| {
            head.iter()
                .zip(l.split('\t'))
                .map(|(k, v)| ((*k).to_owned(), v.to_owned()))
                .collect()
        })
        .collect()
}

/// The header of a TSV, in order.
fn columns(path: &str) -> Vec<String> {
    read(path)
        .lines()
        .next()
        .unwrap_or("")
        .split('\t')
        .map(str::to_owned)
        .collect()
}

fn get<'a>(r: &'a Row, k: &str) -> &'a str {
    r.get(k).map_or("", String::as_str)
}

const TRIO: &str = "sample\tsex\tfather\tmother\tpartner\n\
                    kid\tmale\tdad\tmum\t.\n\
                    dad\tmale\tgrandpa\t.\tmum\n\
                    mum\tfemale\t.\t.\tdad\n";

/// `og family prepare` on this pedigree, with a gVCF for each sample.
fn prepare(d: &Tmp, ped_text: &str, samples: &[&str], extra: &[&str]) -> (i32, String) {
    let ped = write(&d.j("ped.tsv"), ped_text);
    let mut args = vec![
        "prepare".to_owned(),
        "--pedigree".into(),
        ped,
        "--out-dir".into(),
        d.j("out"),
    ];
    for s in samples {
        args.push("--gvcf".into());
        args.push(format!("{s}={}", gvcf(d, &format!("{s}.g.vcf.gz"), s)));
    }
    args.extend(extra.iter().map(|e| (*e).to_owned()));
    run(&owned(&args))
}

const KDM: &[&str] = &["kid", "dad", "mum"];

#[test]
fn prepare_writes_the_pedigrees_and_stops_on_a_bad_one() {
    {
        let d = Tmp::new("prepare");
        let (rc, out) = prepare(&d, TRIO, KDM, &[]);
        assert!(rc == 0, "prepare: a trio passes (exit {rc}): {out}");
        let m: BTreeMap<String, Row> = tsv(&d.j("out/members.tsv"))
            .into_iter()
            .map(|r| (get(&r, "sample").to_owned(), r))
            .collect();
        let g = m.get("grandpa").cloned().unwrap_or_default();
        assert!(
            get(&g, "sequenced") == "no"
                && get(&g, "sex") == "male"
                && get(&g, "in_pedigree") == "no",
            "members.tsv names the grandparent, not sequenced, male: {g:?}"
        );
        assert!(
            KDM.iter().all(|s| get(&m[*s], "sequenced") == "yes"),
            "members.tsv: the three are sequenced: {m:?}"
        );
        let mut ped: Vec<String> = read(&d.j("out/family.ped"))
            .lines()
            .map(str::to_owned)
            .collect();
        ped.sort();
        assert_eq!(
            ped,
            [
                "family\tdad\t-9\t-9\t1\t-9",
                "family\tkid\tdad\tmum\t1\t-9",
                "family\tmum\t-9\t-9\t2\t-9"
            ],
            "somalier's pedigree: three rows, the kid with both parents, the unsequenced grandpa as -9"
        );
        let trios: Vec<String> = read(&d.j("out/trios.ped"))
            .lines()
            .map(str::to_owned)
            .collect();
        assert_eq!(
            trios,
            ["family\tkid\tdad\tmum\t1\t0"],
            "WhatsHap's trios: the kid"
        );
        assert!(
            out.contains("kid and mum: declared parent-child"),
            "prepare says what it will check: {out}"
        );

        // a duo: no trio for WhatsHap
        let (rc, out) = prepare(&d, TRIO, &["kid", "mum"], &[]);
        assert!(rc == 0, "prepare: a duo passes: {out}");
        assert_eq!(read(&d.j("out/trios.ped")), "", "a duo has no trio: {out}");
    }

    let bad: Vec<(&str, String, &str)> = vec![
        (
            "no header",
            "kid\tmale\tdad\tmum\t.\n".into(),
            "header with a 'sample' column",
        ),
        (
            "an unknown column",
            "sample\tsex\tdad\nkid\tmale\tx\n".into(),
            "unknown column dad",
        ),
        (
            "a female father",
            TRIO.replace("dad\tmale", "dad\tfemale"),
            "'dad' is the father of 'kid' but its sex is female",
        ),
        (
            "own parent",
            TRIO.replace("kid\tmale\tdad", "kid\tmale\tkid"),
            "'kid' is named as its own father",
        ),
        (
            "a loop",
            TRIO.replace("dad\tmale\tgrandpa", "dad\tmale\tkid"),
            "among its own ancestors",
        ),
        (
            "a second row",
            format!("{TRIO}kid\tmale\t.\t.\t.\n"),
            "sample 'kid' has a second row",
        ),
        (
            "same father and mother",
            TRIO.replace("kid\tmale\tdad\tmum", "kid\tmale\tdad\tdad"),
            "same person, 'dad', as father and mother",
        ),
        (
            "a bad sex",
            TRIO.replace("kid\tmale", "kid\tboy"),
            "sex 'boy'",
        ),
        (
            "partners that disagree",
            TRIO.replace("mum\tfemale\t.\t.\tdad", "mum\tfemale\t.\t.\tkid"),
            "'dad' names 'mum' as partner, but 'mum' names 'kid'",
        ),
    ];
    for (i, (desc, text, msg)) in bad.iter().enumerate() {
        let d = Tmp::new(&format!("prepare-bad-{i}"));
        let (rc, out) = prepare(&d, text, KDM, &[]);
        assert!(
            rc == 2 && out.contains(msg),
            "prepare stops on {desc} (exit {rc}) and says so: {out}"
        );
    }
    let d = Tmp::new("prepare-stops");
    let (rc, out) = prepare(&d, TRIO, KDM, &["--declared-sex", "mum=male"]);
    assert!(
        rc == 2 && out.contains("the pedigree says 'mum' is female, the samplesheet says male"),
        "prepare stops when the samplesheet's sex contradicts the pedigree's: {out}"
    );
    let (rc, out) = prepare(&d, TRIO, &["kid"], &[]);
    assert!(
        rc == 2 && out.contains("needs two or more"),
        "prepare stops with one sequenced sample: {out}"
    );
    let ped = write(&d.j("ped.tsv"), TRIO);
    let kid = format!("kid={}", gvcf(&d, "a.g.vcf.gz", "kid"));
    let dad = format!("dad={}", gvcf(&d, "b.g.vcf.gz", "HG003"));
    let o2 = d.j("o2");
    let (rc, out) = run(&[
        "prepare",
        "--pedigree",
        &ped,
        "--out-dir",
        &o2,
        "--gvcf",
        &kid,
        "--gvcf",
        &dad,
    ]);
    assert!(
        rc == 2
            && out.contains("holds the sample HG003, not 'dad'")
            && out.contains("bcftools reheader"),
        "prepare stops on a gVCF whose sample is another, with the rename command: {out}"
    );
}

const PAIRS_HEAD: &str = "#sample_a\tsample_b\trelatedness\tibs0\tibs2\tconcordance\thets_a\thets_b\thets_ab\
                          \tshared_hets\thom_alts_a\thom_alts_b\tshared_hom_alts\tn\tx_ibs0\tx_ibs2\
                          \texpected_relatedness\n";

/// somalier's pairs table: (a, b, relatedness, ibs0, n) per row.
fn pairs_file(path: &str, rows: &[(&str, &str, &str, u32, u32)]) -> String {
    let mut text = String::from(PAIRS_HEAD);
    for (a, b, rel, ibs0, n) in rows {
        text += &format!("{a}\t{b}\t{rel}\t{ibs0}\t100\t0.9\t1\t1\t1\t1\t1\t1\t1\t{n}\t0\t0\t-1\n");
    }
    write(path, &text)
}

const TRIO_PAIRS: &[(&str, &str, &str, u32, u32)] = &[
    ("dad", "kid", "0.49", 2, 5000),
    ("kid", "mum", "0.51", 1, 5000),
    ("dad", "mum", "0.01", 600, 5000),
];

#[test]
fn check_against_somalier() {
    let d = Tmp::new("check");
    prepare(&d, TRIO, KDM, &[]);
    let members = d.j("out/members.tsv");
    let pairs = pairs_file(&d.j("pairs.tsv"), TRIO_PAIRS);
    let out_tsv = d.j("rel.tsv");
    let (rc, out) = run(&[
        "check",
        "--members",
        &members,
        "--pairs",
        &pairs,
        "--out",
        &out_tsv,
    ]);
    let rows = tsv(&out_tsv);
    assert!(
        rc == 0 && rows.iter().all(|r| get(r, "verdict") == "ok"),
        "check: the trio's numbers pass (exit {rc}): {out}"
    );
    assert_eq!(rows.len(), 3, "check writes the three pairs: {rows:?}");

    // The same data, the pedigree with kid and dad swapped: dad is now the
    // child of kid and mum, and dad-mum (unrelated) is declared parent-child.
    let swapped = TRIO
        .replace("kid\tmale\tdad\tmum", "@")
        .replace("dad\tmale\tgrandpa\t.\tmum", "kid\tmale\tgrandpa\t.\tmum")
        .replace('@', "dad\tmale\tkid\tmum")
        .replace("mum\tfemale\t.\t.\tdad", "mum\tfemale\t.\t.\tkid");
    let d2 = Tmp(d.0.join("swap"));
    fs::create_dir_all(&d2.0).unwrap();
    let (rc0, out0) = prepare(&d2, &swapped, KDM, &[]);
    assert!(rc0 == 0, "the swapped pedigree is itself valid: {out0}");
    let (rc, out) = run(&[
        "check",
        "--members",
        &d2.j("out/members.tsv"),
        "--pairs",
        &pairs,
        "--out",
        &d2.j("rel.tsv"),
    ]);
    assert!(
        rc == 1,
        "check: two names swapped in the pedigree stop it (exit {rc}): {out}"
    );
    assert!(
        out.contains(
            "dad and mum: declared parent-child (first degree), but somalier measures unrelated"
        ),
        "and the message names the pair that disagrees: {out}"
    );
    assert!(
        out.contains("sample swap"),
        "and says what it may be: {out}"
    );

    let few_rows: Vec<_> = TRIO_PAIRS
        .iter()
        .map(|&(a, b, r, i, _)| (a, b, r, i, 300))
        .collect();
    let few = pairs_file(&d.j("few.tsv"), &few_rows);
    let (rc, out) = run(&[
        "check",
        "--members",
        &members,
        "--pairs",
        &few,
        "--out",
        &out_tsv,
    ]);
    let rows = tsv(&out_tsv);
    assert!(
        rc == 0 && rows.iter().all(|r| get(r, "verdict") == "not_checked"),
        "check: 300 sites is not_checked and does not stop (exit {rc}): {out}"
    );
    let mut sib_rows = vec![("dad", "kid", "0.5", 400, 5000)];
    sib_rows.extend_from_slice(&TRIO_PAIRS[1..]);
    let sib = pairs_file(&d.j("sib.tsv"), &sib_rows);
    let (rc, out) = run(&[
        "check",
        "--members",
        &members,
        "--pairs",
        &sib,
        "--out",
        &out_tsv,
    ]);
    assert!(
        rc == 1 && out.contains("siblings do not"),
        "check: first degree with 8% opposite homozygotes is not parent and child (exit {rc}): {out}"
    );
    let lost = pairs_file(&d.j("lost.tsv"), &TRIO_PAIRS[..2]);
    let (rc, out) = run(&[
        "check",
        "--members",
        &members,
        "--pairs",
        &lost,
        "--out",
        &out_tsv,
    ]);
    assert!(
        rc == 2 && out.contains("no row for dad and mum"),
        "check: a pair somalier did not compare stops it: {out}"
    );
}

const HEAD: &str = "##fileformat=VCFv4.2\n\
                    ##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. \
                    Format: Allele|Consequence|IMPACT|SYMBOL\">\n";

type Site<'a> = (&'a str, u64, &'a str, &'a str, Vec<&'a str>);

/// A gzipped joint VCF with these samples; a genotype without `:` gets
/// DP 30, AD 15,15 and GQ 50.
fn joint(path: &str, rows: &[Site], samples: &[&str], phased: bool) -> String {
    let mut text = format!(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n",
        samples.join("\t")
    );
    for (chrom, pos, r, a, gts) in rows {
        let fmt = if phased {
            "GT:DP:AD:GQ:PS"
        } else {
            "GT:DP:AD:GQ"
        };
        let cells: Vec<String> = gts
            .iter()
            .map(|g| {
                let cell = if g.contains(':') {
                    (*g).to_owned()
                } else {
                    format!("{g}:30:15,15:50")
                };
                if phased { cell + ":1" } else { cell }
            })
            .collect();
        text += &format!(
            "{chrom}\t{pos}\t.\t{r}\t{a}\t50\t.\t.\t{fmt}\t{}\n",
            cells.join("\t")
        );
    }
    gz(path, &text);
    path.to_owned()
}

fn read_inh(path: &str) -> BTreeMap<(String, u64), Row> {
    let mut text = String::new();
    if let Ok(f) = fs::File::open(path) {
        GzDecoder::new(f).read_to_string(&mut text).unwrap();
    }
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    lines
        .map(|l| {
            let r: Vec<&str> = l.split('\t').collect();
            let row: Row = head
                .iter()
                .zip(&r)
                .map(|(k, v)| ((*k).to_owned(), (*v).to_owned()))
                .collect();
            ((r[0].to_owned(), r[1].parse().unwrap()), row)
        })
        .collect()
}

fn par() -> String {
    format!("{REPO}/assets/par_grch38.bed")
}

fn got(inh: &BTreeMap<(String, u64), Row>, chrom: &str, pos: u64) -> (String, String, String) {
    let r = inh
        .get(&(chrom.to_owned(), pos))
        .cloned()
        .unwrap_or_default();
    (
        get(&r, "from_father").to_owned(),
        get(&r, "from_mother").to_owned(),
        get(&r, "how").to_owned(),
    )
}

fn t3(a: &str, b: &str, c: &str) -> (String, String, String) {
    (a.to_owned(), b.to_owned(), c.to_owned())
}

#[test]
fn report_alleles_de_novo_compound_het_phase() {
    //             dad      kid      mum
    let rows: Vec<Site> = vec![
        ("chr1", 100, "A", "G", vec!["0/0", "0/1", "0/1"]), // from mum
        ("chr1", 200, "C", "T", vec!["0/1", "0/1", "0/1"]), // ambiguous
        (
            "chr1",
            300,
            "G",
            "A",
            vec!["0/0:30:30,0:50", "0/1:30:14,16:50", "0/0:30:30,0:50"],
        ), // planted de novo
        ("chr1", 400, "T", "C", vec!["1/1", "0/0", "0/0"]), // opposite homozygote: inconsistent, not de novo
        ("chr1", 500, "A", "C", vec!["0/1", "1/1", "0/1"]), // both alt
        ("chr2", 100, "A", "G", vec!["0/1", "0/1", "0/0"]), // gene G1, from dad
        ("chr2", 200, "C", "T", vec!["0/0", "0/1", "0/1"]), // gene G1, from mum -> in trans
        ("chr3", 100, "A", "G", vec!["0/0", "0/1", "0/1"]), // gene G2, from mum
        ("chr3", 200, "C", "T", vec!["0/0", "0/1", "1/1"]), // gene G2, from mum -> in cis
        ("chrX", 500000, "C", "T", vec!["0/1", "0/1", "0/0"]), // PAR1: autosomal, from dad
        ("chrX", 74000000, "G", "A", vec!["0/0", "1/1", "0/1"]), // son's X: from mum only
        ("chrX", 74000100, "G", "A", vec!["1/1", "0/0", "0/1"]), // son's X: dad's allele does not count
        ("chrM", 100, "A", "G", vec!["0/0", "0/1", "0/0"]),      // left out
    ];
    let samples = ["dad", "kid", "mum"];
    let par = par();
    {
        let d = Tmp::new("report");
        prepare(&d, TRIO, KDM, &[]);
        let out = d.j("out");
        let members = d.j("out/members.tsv");
        let pairs = pairs_file(&d.j("pairs.tsv"), TRIO_PAIRS);
        let rel = d.j("out/relationships.tsv");
        run(&[
            "check",
            "--members",
            &members,
            "--pairs",
            &pairs,
            "--out",
            &rel,
        ]);
        let vcf = joint(&d.j("family.vcf.gz"), &rows, &samples, false);
        let cand = write(
            &d.j("kid_clinical.vcf"),
            &format!(
                "{HEAD}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tkid\n\
                 chr2\t100\t.\tA\tG\t50\tPASS\tCSQ=G|missense_variant|MODERATE|G1\tGT\t0/1\n\
                 chr2\t200\t.\tC\tT\t50\tPASS\tCSQ=T|stop_gained|HIGH|G1\tGT\t0/1\n\
                 chr3\t100\t.\tA\tG\t50\tPASS\tCSQ=G|missense_variant|MODERATE|G2\tGT\t0/1\n\
                 chr3\t200\t.\tC\tT\t50\tPASS\tCSQ=T|missense_variant|MODERATE|G2\tGT\t0/1\n"
            ),
        );
        let kid_cand = format!("kid={cand}");
        let args = [
            "report",
            "--members",
            &members,
            "--relationships",
            &rel,
            "--vcf",
            &vcf,
            "--candidates",
            &kid_cand,
            "--par-bed",
            &par,
            "--out-dir",
            &out,
        ];
        let (rc, log) = run(&args);
        assert!(rc == 0, "report runs (exit {rc}): {log}");
        let inh = read_inh(&d.j("out/kid_inheritance.tsv.gz"));
        assert_eq!(
            got(&inh, "chr1", 100),
            t3("A", "G", "genotype"),
            "chr1:100 the G came from the mother"
        );
        assert_eq!(
            got(&inh, "chr1", 200),
            t3(".", ".", "ambiguous"),
            "chr1:200 all heterozygous: ambiguous"
        );
        assert_eq!(
            got(&inh, "chr1", 400).2,
            "inconsistent",
            "chr1:400 opposite homozygote: inconsistent"
        );
        assert_eq!(
            got(&inh, "chr1", 500),
            t3("C", "C", "genotype"),
            "chr1:500 the alt from both"
        );
        assert_eq!(
            got(&inh, "chrX", 74000000),
            t3(".", "A", "genotype"),
            "chrX:74000000 the son's X from the mother only"
        );
        assert_eq!(
            got(&inh, "chrX", 74000100),
            t3(".", "G", "genotype"),
            "chrX:74000100 the father's X allele does not count for a son"
        );
        assert_eq!(
            got(&inh, "chrX", 500000),
            t3("T", "C", "genotype"),
            "chrX:500000 (PAR1) from the father, as an autosome"
        );
        assert!(
            !inh.contains_key(&("chrM".to_owned(), 100)),
            "chrM is left out"
        );
        let dn = tsv(&d.j("out/de_novo.tsv"));
        let dn_cells: Vec<[&str; 4]> = dn
            .iter()
            .map(|r| {
                [
                    get(r, "chrom"),
                    get(r, "pos"),
                    get(r, "alt"),
                    get(r, "quality"),
                ]
            })
            .collect();
        assert_eq!(
            dn_cells,
            [["chr1", "300", "A", "pass"]],
            "the planted new allele is the one de novo candidate, quality pass"
        );
        let ch = tsv(&d.j("out/compound_het.tsv"));
        let rel_of: BTreeMap<&str, [&str; 3]> = ch
            .iter()
            .map(|r| {
                (
                    get(r, "gene"),
                    [get(r, "a_from"), get(r, "b_from"), get(r, "phase")],
                )
            })
            .collect();
        assert_eq!(
            rel_of.get("G1"),
            Some(&["father", "mother", "in trans"]),
            "G1: one from each parent, in trans: {ch:?}"
        );
        assert_eq!(
            rel_of.get("G2"),
            Some(&["mother", "mother", "in cis"]),
            "G2: both from the mother, in cis: {ch:?}"
        );
        let md = read(&d.j("out/family_report.md"));
        assert!(
            md.contains("De novo candidates (an allele neither parent has): 1, 1")
                && md.contains(
                    "in trans (one from each parent, both copies of the gene hit) 1, in cis (both from one parent) 1"
                ),
            "the report states the de novo count and the pairs: {md}"
        );

        let first = fs::read(d.j("out/kid_inheritance.tsv.gz")).unwrap();
        run(&args);
        assert!(
            fs::read(d.j("out/kid_inheritance.tsv.gz")).unwrap() == first,
            "a second run writes the same bytes"
        );

        // A low-quality copy of the planted site: still listed, quality low
        let low_rows: Vec<Site> = vec![(
            "chr1",
            300,
            "G",
            "A",
            vec!["0/0:5:5,0:10", "0/1:30:14,16:50", "0/0:30:30,0:50"],
        )];
        let low = joint(&d.j("low.vcf.gz"), &low_rows, &samples, false);
        let low_out = d.j("low");
        run(&[
            "report",
            "--members",
            &members,
            "--relationships",
            &rel,
            "--vcf",
            &low,
            "--par-bed",
            &par,
            "--out-dir",
            &low_out,
        ]);
        let dn = tsv(&d.j("low/de_novo.tsv"));
        assert_eq!(
            dn.iter().map(|r| get(r, "quality")).collect::<Vec<_>>(),
            ["low"],
            "a candidate with a parent at depth 5 is listed as low: {dn:?}"
        );

        // WhatsHap's phase: kid's first haplotype holds the father's allele at chr2:100
        let ph_rows: Vec<Site> = vec![
            ("chr2", 100, "A", "G", vec!["0|1", "1|0", "0|0"]),
            ("chr2", 200, "C", "T", vec!["0|0", "0|1", "0|1"]),
        ];
        let ph = joint(&d.j("phased.vcf.gz"), &ph_rows, &samples, true);
        let ph_out = d.j("ph");
        run(&[
            "report",
            "--members",
            &members,
            "--relationships",
            &rel,
            "--vcf",
            &vcf,
            "--phased",
            &ph,
            "--par-bed",
            &par,
            "--out-dir",
            &ph_out,
        ]);
        let md = read(&d.j("ph/family_report.md"));
        assert!(
            md.contains(
                "WhatsHap phased 2 heterozygous sites of kid by inheritance. At 2 of them the genotypes also give the \
                 origin, and the phase agrees with it at 2 (100.0%)"
            ),
            "WhatsHap's phase agrees with the genotype origin at both sites: {md}"
        );

        // An unsorted joint VCF would split the chromosome map wrongly: it stops, writing nothing
        let bad = joint(
            &d.j("unsorted.vcf.gz"),
            &[rows[1].clone(), rows[0].clone()],
            &samples,
            false,
        );
        let uns = d.j("unsorted");
        let (rc, log) = run(&[
            "report",
            "--members",
            &members,
            "--relationships",
            &rel,
            "--vcf",
            &bad,
            "--out-dir",
            &uns,
        ]);
        assert!(
            rc == 2
                && log.contains("is not sorted: chr1:100 after chr1:200")
                && !Path::new(&d.j("unsorted/kid_inheritance.tsv.gz")).exists(),
            "an unsorted joint VCF stops the report (exit 2), naming the two records: exit {rc}: {log}"
        );
    }

    // A duo: the father is not sequenced; his allele is inferred
    let d = Tmp::new("report-duo");
    prepare(&d, TRIO, &["kid", "mum"], &[]);
    let out = d.j("out");
    let members = d.j("out/members.tsv");
    let pairs = pairs_file(&d.j("pairs.tsv"), &[("kid", "mum", "0.5", 1, 5000)]);
    let rel = d.j("out/relationships.tsv");
    let (rc, log) = run(&[
        "check",
        "--members",
        &members,
        "--pairs",
        &pairs,
        "--out",
        &rel,
    ]);
    assert!(rc == 0, "a duo passes the check: {log}");
    let duo_rows: Vec<Site> = vec![
        ("chr1", 100, "A", "G", vec!["0/1", "0/0"]),
        ("chr1", 200, "A", "G", vec!["0/1", "0/1"]),
        ("chr1", 300, "A", "G", vec!["0/1", "1/1"]),
    ];
    let vcf = joint(&d.j("duo.vcf.gz"), &duo_rows, &["kid", "mum"], false);
    run(&[
        "report",
        "--members",
        &members,
        "--relationships",
        &rel,
        "--vcf",
        &vcf,
        "--par-bed",
        &par,
        "--out-dir",
        &out,
    ]);
    let inh = read_inh(&d.j("out/kid_inheritance.tsv.gz"));
    let r100 = &inh[&("chr1".to_owned(), 100)];
    assert_eq!(
        (get(r100, "from_father"), get(r100, "father_gt")),
        ("G", "."),
        "duo: the G the mother lacks came from the unsequenced father: {r100:?}"
    );
    assert_eq!(
        get(&inh[&("chr1".to_owned(), 200)], "how"),
        "ambiguous",
        "duo: child and mother both heterozygous is ambiguous"
    );
    let md = read(&d.j("out/family_report.md"));
    assert!(
        md.contains("The father is not sequenced. At 2 of these sites (66.7%)")
            && md.contains("needs both parents sequenced"),
        "duo: the report states the inferred share of the father, and calls no de novo: {md}"
    );
}

/// A daughter's X comes from both parents: an allele the father carries is
/// not de novo.
#[test]
fn daughter_x_from_both_parents() {
    let daughter = TRIO.replace("kid\tmale", "kid\tfemale");
    //             dad    kid    mum
    let rows: Vec<Site> = vec![
        ("chrX", 5000000, "A", "G", vec!["1/1", "1/1", "0/0"]), // dad has G: mother's dropout, not de novo
        (
            "chrX",
            5000100,
            "A",
            "G",
            vec!["0/0", "0/1:30:14,16:50", "0/0"],
        ), // neither parent has G: de novo
        ("chrX", 5000200, "A", "G", vec!["1/1", "0/1", "0/0"]), // one from each parent
    ];
    let d = Tmp::new("daughter");
    prepare(&d, &daughter, KDM, &[]);
    let out = d.j("out");
    let members = d.j("out/members.tsv");
    let pairs = pairs_file(&d.j("pairs.tsv"), TRIO_PAIRS);
    let rel = d.j("out/relationships.tsv");
    run(&[
        "check",
        "--members",
        &members,
        "--pairs",
        &pairs,
        "--out",
        &rel,
    ]);
    let vcf = joint(&d.j("family.vcf.gz"), &rows, &["dad", "kid", "mum"], false);
    let par = par();
    let (rc, log) = run(&[
        "report",
        "--members",
        &members,
        "--relationships",
        &rel,
        "--vcf",
        &vcf,
        "--par-bed",
        &par,
        "--out-dir",
        &out,
    ]);
    assert!(rc == 0, "daughter: report runs (exit {rc}): {log}");
    let inh = read_inh(&d.j("out/kid_inheritance.tsv.gz"));
    assert_eq!(
        got(&inh, "chrX", 5000000).2,
        "inconsistent",
        "daughter chrX:5000000 (dad 1/1, mum 0/0, kid 1/1) is inconsistent"
    );
    let (f, m, _) = got(&inh, "chrX", 5000200);
    assert_eq!(
        (f.as_str(), m.as_str()),
        ("G", "A"),
        "daughter chrX:5000200 one allele from each parent"
    );
    let dn = tsv(&d.j("out/de_novo.tsv"));
    assert_eq!(
        dn.iter()
            .map(|r| (get(r, "chrom"), get(r, "pos")))
            .collect::<Vec<_>>(),
        [("chrX", "5000100")],
        "daughter: only the allele neither parent has is a de novo candidate, not the father's: {dn:?}"
    );
}

fn two() -> String {
    format!("{REPO}/tests/fixtures/families/two_families.tsv")
}

fn plan(args: &[&str]) -> (i32, String) {
    let mut a = vec!["plan"];
    a.extend_from_slice(args);
    run(&a)
}

/// The plan table by analysis.
fn by_analysis(path: &str) -> BTreeMap<String, Vec<Row>> {
    let mut out: BTreeMap<String, Vec<Row>> = BTreeMap::new();
    for r in tsv(path) {
        out.entry(get(&r, "analysis").to_owned())
            .or_default()
            .push(r);
    }
    out
}

fn first<'a>(p: &'a BTreeMap<String, Vec<Row>>, analysis: &str, col: &str) -> &'a str {
    p.get(analysis)
        .and_then(|v| v.first())
        .map_or("", |r| get(r, col))
}

/// `og live family state DIR`: the JSON it prints.
fn live_state(dir: &str) -> (u8, String, Value) {
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = og_live::family::cli(&["state".to_owned(), dir.to_owned()], &mut out, &mut err);
    let text = String::from_utf8_lossy(&out).into_owned() + &String::from_utf8_lossy(&err);
    let v = if code == 0 {
        json::loads(&String::from_utf8_lossy(&out)).unwrap_or(Value::Null)
    } else {
        Value::Null
    };
    (code, text, v)
}

fn field<'a>(v: &'a Value, k: &str) -> Option<&'a Value> {
    match v {
        Value::Object(o) => o.get(k),
        _ => None,
    }
}

fn text_of(v: Option<&Value>) -> Option<&str> {
    match v {
        Some(Value::Str(s)) => Some(s),
        _ => None,
    }
}

#[test]
fn families_file() {
    let two = two();
    let fix = format!("{REPO}/tests/fixtures/families");
    {
        let d = Tmp::new("families");
        // the two-family file splits into two pedigrees
        for (fam, want) in [
            ("family_a", vec!["child_a", "father_a", "mother_a"]),
            ("family_b", vec!["partner_x", "partner_y"]),
        ] {
            let dir = d.0.display().to_string();
            let (rc, out) = plan(&[
                "--families",
                &two,
                "--family",
                fam,
                "--write-pedigree",
                &dir,
            ]);
            let ped = d.j(&format!("{fam}/pedigree.tsv"));
            let rows = tsv(&ped);
            assert!(
                rc == 0 && rows.iter().map(|r| get(r, "sample")).collect::<Vec<_>>() == want,
                "plan --write-pedigree: {fam} gets its own pedigree (exit {rc}): {out}{rows:?}"
            );
            assert!(
                !rows.is_empty()
                    && columns(&ped) == ["sample", "sex", "father", "mother", "partner"],
                "{fam}'s pedigree has the pedigree columns, no family column: {:?}",
                rows.first()
            );
        }
        let (code, text, st) = live_state(&d.j("family_a"));
        let kid = match field(&st, "members") {
            Some(Value::Array(m)) => m
                .iter()
                .find(|x| text_of(field(x, "sample")) == Some("child_a")),
            _ => None,
        };
        assert!(
            code == 0
                && field(&st, "error").is_none()
                && kid.is_some_and(|k| {
                    (text_of(field(k, "father")), text_of(field(k, "mother")))
                        == (Some("father_a"), Some("mother_a"))
                }),
            "the written pedigree reads back (og live family state): child_a's parents are father_a and mother_a: {text}"
        );

        // a file without the column gives the same members.tsv as before the column existed
        let example = format!("{REPO}/assets/pedigree.example.tsv");
        let gvcfs: Vec<String> = ["child_a", "father_a", "mother_a"]
            .iter()
            .flat_map(|s| {
                [
                    "--gvcf".to_owned(),
                    format!("{s}={}", gvcf(&d, &format!("{s}.g.vcf.gz"), s)),
                ]
            })
            .collect();
        let mut args = vec![
            "prepare".to_owned(),
            "--pedigree".into(),
            example,
            "--out-dir".into(),
            d.j("old"),
        ];
        args.extend(gvcfs.iter().cloned());
        let (rc, out) = run(&owned(&args));
        let got_members = if rc == 0 {
            fs::read(d.j("old/members.tsv")).unwrap_or_default()
        } else {
            Vec::new()
        };
        assert!(
            got_members == fs::read(format!("{fix}/members.golden.tsv")).unwrap(),
            "a pedigree without the family column: members.tsv is byte for byte the golden file: {out}{:?}",
            String::from_utf8_lossy(&got_members)
        );
        // the same family from the families file: the same members.tsv
        args[2] = two.clone();
        args[4] = d.j("new");
        let mut with_family = args.clone();
        with_family.extend(["--family".to_owned(), "family_a".to_owned()]);
        let (rc, out) = run(&owned(&with_family));
        assert!(
            rc == 0 && fs::read(d.j("new/members.tsv")).unwrap_or_default() == got_members,
            "the same family from the families file gives the same members.tsv: {out}"
        );
        let (rc, out) = run(&owned(&args));
        assert!(
            rc == 2 && out.contains("holds 2 families (family_a, family_b)"),
            "prepare on a two-family file without --family stops and lists the families: {out}"
        );
        let (rc, out) = plan(&["--families", &two]);
        assert!(
            rc == 2 && out.contains("holds 2 families (family_a, family_b)"),
            "plan on a two-family file without --family or --sample stops and lists the families: {out}"
        );
        let (rc, out) = plan(&["--families", &two, "--family", "family_c"]);
        assert!(
            rc == 2 && out.contains("no family 'family_c'; its families are: family_a, family_b"),
            "an unknown --family is named, with the file's families: {out}"
        );
        let (rc, out) = plan(&["--families", &two, "--sample", "partner_y", "--which"]);
        assert!(
            rc == 0 && out.trim().lines().last() == Some("family_b"),
            "--which prints the family of the sample: {out}"
        );
        let (rc, out) = plan(&["--families", &two, "--sample", "nobody", "--which"]);
        assert!(
            rc == 0 && !out.contains("family_"),
            "--which: a sample in no family prints no family and exits 0: {out}"
        );
    }

    let head = "family\tsample\tsex\tfather\tmother\tpartner\n";
    let bad = [
        (
            "a sample in two families",
            format!(
                "{head}fa\tkid\tmale\t.\t.\t.\nfb\tkid\tmale\t.\t.\t.\nfb\tmum\tfemale\t.\t.\t.\n"
            ),
            "line 3: sample 'kid' is in family 'fb' here and in family 'fa' on line 2",
        ),
        (
            "a family name breaking the sample-name rule",
            format!("{head}f/a\tkid\tmale\t.\t.\t.\n"),
            "line 2: family name 'f/a' may hold only letters",
        ),
        (
            "a family name equal to a sample",
            format!("{head}kid\tkid\tmale\t.\t.\t.\nkid\tmum\tfemale\t.\t.\t.\n"),
            "the family name 'kid' is also a person's name (line 2)",
        ),
        (
            "a parent in another family",
            format!("{head}fa\tkid\tmale\tdad\t.\t.\nfb\tdad\tmale\t.\t.\t.\n"),
            "line 2: the father of 'kid', 'dad', is in family 'fb' (line 3)",
        ),
        (
            "a partner in another family",
            format!("{head}fa\tdad\tmale\t.\t.\tmum\nfb\tmum\tfemale\t.\t.\tdad\n"),
            "line 2: the partner of 'dad', 'mum', is in family 'fb' (line 3)",
        ),
        (
            "the column filled on some rows only",
            format!("{head}fa\tkid\tmale\t.\t.\t.\n.\tmum\tfemale\t.\t.\t.\n"),
            "line 3: no family name, but other rows have one (line 2: 'fa')",
        ),
    ];
    for (i, (desc, text, msg)) in bad.iter().enumerate() {
        let d = Tmp::new(&format!("families-bad-{i}"));
        let f = write(&d.j("f.tsv"), text);
        let (rc, out) = plan(&["--families", &f, "--family", "fa"]);
        assert!(
            rc == 2 && out.contains(msg) && out.contains("f.tsv"),
            "plan stops on {desc} (exit {rc}), naming the file and line: {out}"
        );
    }
    let (rc, out) = plan(&[
        "--families",
        &two,
        "--sample",
        "partner_x",
        "--sex",
        "female",
    ]);
    assert!(
        rc == 2
            && out.contains("says 'partner_x' is male, the sex given to run-all.sh says female"),
        "a sex that contradicts the file stops plan (exit {rc}): {out}"
    );
    let (rc, out) = plan(&["--families", &two, "--sample", "partner_x", "--sex", "male"]);
    assert!(rc == 0, "the sex the file gives passes: {out}");
}

#[test]
fn plan_runs_skips_and_waits() {
    let two = two();
    let d = Tmp::new("plan");
    let out = d.j("plan.tsv");
    // one member with a gVCF: everything waits, naming the missing one
    let (rc, log) = plan(&[
        "--families",
        &two,
        "--sample",
        "partner_x",
        "--have-gvcf",
        "partner_x",
        "--out",
        &out,
    ]);
    let p = if rc == 0 {
        by_analysis(&out)
    } else {
        BTreeMap::new()
    };
    assert!(
        rc == 0 && columns(&out) == ["family", "analysis", "members", "decision", "reason"],
        "plan writes the columns family, analysis, members, decision, reason: {log}"
    );
    assert!(
        first(&p, "joint_call", "decision") == "waiting"
            && first(&p, "joint_call", "reason")
                == "partner_y has no gVCF yet: run run-all.sh partner_y first",
        "one member with a gVCF: the joint call waits for partner_y: {p:?}"
    );
    assert_eq!(
        first(&p, "couple_carrier", "decision"),
        "waiting",
        "and the couple screen waits too: {p:?}"
    );
    assert!(
        p.values().flatten().all(|r| get(r, "family") == "family_b"),
        "every row is family_b's: {p:?}"
    );

    // the trio: phasing runs
    plan(&[
        "--families",
        &two,
        "--family",
        "family_a",
        "--have-gvcf",
        "child_a",
        "father_a",
        "mother_a",
        "--have-clinvar",
        "child_a",
        "father_a",
        "mother_a",
        "--ibd-panel",
        "yes",
        "--out",
        &out,
    ]);
    let p = by_analysis(&out);
    assert!(
        first(&p, "joint_call", "decision") == "runs"
            && first(&p, "joint_call", "members") == "child_a,father_a,mother_a",
        "a trio: the joint call runs on the three: {p:?}"
    );
    assert_eq!(
        (
            first(&p, "trio_phasing", "decision"),
            first(&p, "trio_phasing", "members")
        ),
        ("runs", "child_a"),
        "a trio: phasing by inheritance runs for child_a: {p:?}"
    );
    assert_eq!(
        first(&p, "inheritance_map", "decision"),
        "runs",
        "a trio: the chromosome map runs: {p:?}"
    );
    assert_eq!(
        first(&p, "shared_segments", "decision"),
        "runs",
        "with the panel the shared segments run: {p:?}"
    );
    assert_eq!(
        (
            first(&p, "couple_carrier", "decision"),
            first(&p, "couple_carrier", "members")
        ),
        ("runs", "father_a,mother_a"),
        "the parents, partners with ClinVar hits, get the couple screen: {p:?}"
    );

    // a duo: phasing skipped, the map runs; no panel: segments skipped; a partner without ClinVar: couple skipped
    plan(&[
        "--families",
        &two,
        "--family",
        "family_a",
        "--have-gvcf",
        "child_a",
        "mother_a",
        "--have-clinvar",
        "mother_a",
        "--ibd-panel",
        "no",
        "--out",
        &out,
    ]);
    let p = by_analysis(&out);
    assert!(
        first(&p, "joint_call", "decision") == "runs"
            && first(&p, "joint_call", "members") == "child_a,mother_a"
            && first(&p, "joint_call", "reason").contains("left out: father_a has no gVCF yet"),
        "a duo: the joint call runs on two and names the one left out: {p:?}"
    );
    assert!(
        first(&p, "trio_phasing", "decision") == "skipped"
            && first(&p, "trio_phasing", "reason").contains("no child has both parents sequenced"),
        "a duo: phasing by inheritance is skipped: {p:?}"
    );
    assert_eq!(
        first(&p, "inheritance_map", "decision"),
        "runs",
        "a duo: the chromosome map runs: {p:?}"
    );
    assert!(
        first(&p, "shared_segments", "decision") == "skipped"
            && first(&p, "shared_segments", "reason").contains("setup.sh --ancestry-fine-panel"),
        "--ibd-panel no: the shared segments are skipped: {p:?}"
    );
    assert!(
        first(&p, "couple_carrier", "decision") == "waiting"
            && first(&p, "couple_carrier", "reason")
                == "father_a has no ClinVar hits yet (step 06 writes clinvar/father_a_clinvar_hits.vcf): \
                    run run-all.sh father_a first",
        "a partner without ClinVar hits: the couple screen waits for him: {p:?}"
    );
    plan(&[
        "--families",
        &two,
        "--family",
        "family_b",
        "--have-gvcf",
        "partner_x",
        "partner_y",
        "--have-clinvar",
        "partner_x",
        "--out",
        &out,
    ]);
    let p = by_analysis(&out);
    assert!(
        first(&p, "couple_carrier", "decision") == "waiting"
            && first(&p, "couple_carrier", "reason").contains("partner_y has no ClinVar hits yet"),
        "a partner with a gVCF but no ClinVar hits: the couple screen waits, naming her: {p:?}"
    );
    assert!(
        first(&p, "trio_phasing", "decision") == "skipped"
            && first(&p, "inheritance_map", "decision") == "skipped",
        "partners only: no child, so no phasing and no map: {p:?}"
    );
    // step 39 never reads a gVCF: a partner run from a VCF alone, with ClinVar hits, still gets the screen
    plan(&[
        "--families",
        &two,
        "--family",
        "family_b",
        "--have-gvcf",
        "partner_x",
        "--have-clinvar",
        "partner_x",
        "partner_y",
        "--out",
        &out,
    ]);
    let p = by_analysis(&out);
    assert_eq!(
        (
            first(&p, "couple_carrier", "decision"),
            first(&p, "couple_carrier", "members")
        ),
        ("runs", "partner_x,partner_y"),
        "both partners with ClinVar hits, one without a gVCF: the couple screen runs: {p:?}"
    );
    assert_eq!(
        first(&p, "joint_call", "decision"),
        "waiting",
        "while the joint call waits for partner_y's gVCF: {p:?}"
    );
    plan(&[
        "--families",
        &two,
        "--family",
        "family_b",
        "--have-gvcf",
        "partner_x",
        "partner_y",
        "--have-clinvar",
        "partner_x",
        "--out",
        &out,
    ]);
    let p = by_analysis(&out);
    assert!(
        first(&p, "trio_phasing", "decision") == "skipped"
            && first(&p, "inheritance_map", "decision") == "skipped",
        "partners only, again: no child, so no phasing and no map: {p:?}"
    );

    // prepare writes the same table as plan, for the same members
    let (rc, log) = plan(&[
        "--families",
        &two,
        "--family",
        "family_a",
        "--have-gvcf",
        "child_a",
        "father_a",
        "mother_a",
        "--have-clinvar",
        "father_a",
        "--ibd-panel",
        "yes",
        "--out",
        &out,
    ]);
    let mut args = vec![
        "prepare".to_owned(),
        "--pedigree".into(),
        two.clone(),
        "--family".into(),
        "family_a".into(),
        "--have-clinvar".into(),
        "father_a".into(),
        "--ibd-panel".into(),
        "yes".into(),
        "--pedigree-out".into(),
        d.j("ped.tsv"),
        "--out-dir".into(),
        d.j("prep"),
    ];
    for s in ["child_a", "father_a", "mother_a"] {
        args.push("--gvcf".into());
        args.push(format!("{s}={}", gvcf(&d, &format!("{s}.g.vcf.gz"), s)));
    }
    let (rc2, log2) = run(&owned(&args));
    let pp = d.j("prep/family_plan.tsv");
    assert!(
        rc == 0
            && rc2 == 0
            && Path::new(&pp).exists()
            && fs::read(&pp).unwrap() == fs::read(&out).unwrap(),
        "prepare's family_plan.tsv is plan's, byte for byte: {log}{log2}"
    );
    assert!(
        rc2 == 0
            && tsv(&d.j("ped.tsv"))
                .iter()
                .map(|r| get(r, "sample"))
                .collect::<Vec<_>>()
                == ["child_a", "father_a", "mother_a"],
        "prepare --pedigree-out writes the family's own pedigree: {log2}"
    );
}
