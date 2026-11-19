//! `og archaic` (step 43) on a synthetic archaic panel, then the summary and
//! the report lines; on amd64 in CI also with the real IBDmix image. It
//! replaces `tests/test_archaic_report.py`.
//!
//! 1. table: a 1000 Genomes site keeps its frequencies and gets each archaic
//!    genome's ALT count, 10 plus it outside that genome's mask or inside an
//!    excluded region; an archaic record with another ALT or another REF is
//!    unknown (9); a site only an archaic genome varies at gets frequency 0;
//!    a site no archaic genome was called at is left out; a position given
//!    twice by 1000 Genomes is left out;
//! 2. strict: a lifted row outside the GRCh38 strict mask is kept once per
//!    position with every archaic call marked masked (9 stays 9);
//! 3. genotypes: a masked archaic call goes to IBDmix as its ALT count and
//!    its position into that genome's mask file; the stand-in columns carry
//!    the population's ALT count; your column is your ALT count, 9 for a
//!    non-PASS record, and where your VCF has no record 0 inside a reference
//!    block of the gVCF and 9 outside every one; a site where your VCF has
//!    another ALT is left out; a chromosome your VCF has no record on is not
//!    analysed; the split-VCF rules; unsorted input stops and keeps the
//!    previous outputs;
//! 4. report: Altai and Vindija segments merge into Neanderthal ones, a
//!    Denisova segment that overlaps one is not Denisovan, segments under
//!    50 kb or LOD 4 are dropped, the share is of the analysed span, and the
//!    genes of the longest segments are named; og summary collect reads the
//!    summary and marks a damaged one unreadable, and the report lines give
//!    both shares;
//! 5. with IBDmix itself (`IBDMIX_IMAGE`, docker; amd64 only, run in CI or
//!    with `OG_DOCKER_TESTS=1`): a planted Neanderthal and a planted
//!    Denisovan stretch come out as one segment each, the same person
//!    without them gets none; masked mismatches stop a segment where the
//!    same calls written as unknown let it run through; a stretch with no
//!    call gives a segment only when a gVCF reference block says it was read.

use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use og_ancestry::archaic::cli;
use pycompat::json::Value;

const TABLE_HEADER: &str = "chrom\tpos\tref\talt\tAltai\tVindija\tDenisova\tAF_ALL\tAF_AFR\tAF_AMR\tAF_EAS\tAF_EUR\tAF_SAS";
const VCF_HEAD: &str =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";
const ARCHAICS: [&str; 3] = ["Altai", "Vindija", "Denisova"];

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

/// A scratch folder that is removed when the test ends, also when it fails.
struct Tmp(PathBuf);

impl std::ops::Deref for Tmp {
    type Target = Path;
    fn deref(&self) -> &Path {
        &self.0
    }
}

impl AsRef<Path> for Tmp {
    fn as_ref(&self) -> &Path {
        &self.0
    }
}

impl AsRef<std::ffi::OsStr> for Tmp {
    fn as_ref(&self) -> &std::ffi::OsStr {
        self.0.as_os_str()
    }
}

impl Tmp {
    /// A new folder that only this test owns. `create_dir` fails when the
    /// path exists, so a run never shares or deletes another run's folder.
    fn new(kind: &str, tag: &str) -> Tmp {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let d = std::env::temp_dir().join(format!(
            "og-ancestry-{kind}-{}-{tag}-{nanos}",
            std::process::id()
        ));
        std::fs::create_dir(&d).unwrap();
        Tmp(d)
    }
}

impl Drop for Tmp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn tmp(tag: &str) -> Tmp {
    Tmp::new("archaic", tag)
}

fn p(d: &Path, name: &str) -> String {
    d.join(name).display().to_string()
}

fn write(path: &str, text: &str) {
    std::fs::write(path, text).unwrap();
}

fn gzip(path: &str, text: &str) {
    let mut z = GzEncoder::new(std::fs::File::create(path).unwrap(), Compression::default());
    z.write_all(text.as_bytes()).unwrap();
    z.finish().unwrap();
}

fn gunzip(path: &str) -> String {
    let mut s = String::new();
    GzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut s)
        .unwrap();
    s
}

fn og(args: &[&str]) -> (u8, String) {
    let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = cli(&args, &mut out, &mut err);
    (code, String::from_utf8(err).unwrap())
}

fn ref_block(chrom: &str, start: u32, end: u32) -> String {
    format!("{chrom}\t{start}\t.\tN\t<*>\t0\t.\tEND={end}\tGT:GQ:MIN_DP\t0/0:50:30\n")
}

/// The records of a table VCF by position, and each one's INFO.
type Rows = BTreeMap<u32, Vec<String>>;
type Info = BTreeMap<u32, BTreeMap<String, String>>;

fn vcf_rows(path: &str) -> (Rows, Info) {
    let rows: Rows = std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let c: Vec<String> = l.split('\t').map(str::to_owned).collect();
            (c[1].parse().unwrap(), c)
        })
        .collect();
    let info = rows
        .iter()
        .map(|(k, c)| {
            (
                *k,
                c[7].split(';')
                    .map(|kv| {
                        let (a, b) = kv.split_once('=').unwrap();
                        (a.to_owned(), b.to_owned())
                    })
                    .collect(),
            )
        })
        .collect();
    (rows, info)
}

fn codes(i: &Info, pos: u32) -> [&str; 3] {
    ARCHAICS.map(|g| i[&pos][g].as_str())
}

/// The inputs of part 1 in `d`.
fn table_inputs(d: &Path) -> Vec<String> {
    let afs = "0.1\t0.2\t0.05\t0\t0.15\t0.1";
    write(
        &p(d, "modern.tsv"),
        &format!(
            // every archaic genome called; Vindija has another ALT and
            // Denisova another REF; no archaic record; a split multi-allelic
            // site.
            "100\tA\tG\t{afs}\n200\tC\tT\t{afs}\n300\tG\tA\t{afs}\n400\tT\tC\t{afs}\n400\tT\tG\t{afs}\n"
        ),
    );
    write(
        &p(d, "altai.tsv"),
        "100\tA\tG\t1/1\n200\tC\t.\t0/0\n500\tG\tC\t0/1\n600\tA\t.\t0/0\n",
    );
    write(
        &p(d, "vindija.tsv"),
        "100\tA\t.\t0/0\n200\tC\tA\t0/1\n500\tG\t.\t0/0\n600\tA\t.\t0/0\n",
    );
    write(
        &p(d, "denisova.tsv"),
        "100\tA\tG\t0|1\n200\tT\tC\t0/1\n400\tT\tC\t1/1\n500\tG\tT\t1/1\n",
    );
    [
        "table".to_owned(),
        "--chrom".to_owned(),
        "22".to_owned(),
        "--modern".to_owned(),
        p(d, "modern.tsv"),
        "--archaic".to_owned(),
        format!("Altai={}", p(d, "altai.tsv")),
        "--archaic".to_owned(),
        format!("Vindija={}", p(d, "vindija.tsv")),
        "--archaic".to_owned(),
        format!("Denisova={}", p(d, "denisova.tsv")),
    ]
    .to_vec()
}

fn table(common: &[String], extra: &[&str], out: &str) -> (u8, String) {
    let mut a: Vec<&str> = common.iter().map(String::as_str).collect();
    a.extend(extra);
    a.extend(["--out", out]);
    og(&a)
}

#[test]
fn the_table_codes_each_archaic_call() {
    let d = tmp("table");
    let common = table_inputs(&d);
    let (code, err) = table(&common, &[], &p(&d, "t.vcf"));
    assert_eq!(code, 0, "{err}");
    let (recs, info) = vcf_rows(&p(&d, "t.vcf"));
    assert_eq!(
        recs.keys().copied().collect::<Vec<_>>(),
        [100, 200, 500],
        "the modern sites with an archaic call and the archaic-only site"
    );
    assert_eq!(
        codes(&info, 100),
        ["2", "0", "1"],
        "each genome's ALT count"
    );
    assert_eq!(
        [&info[&100]["AF_ALL"], &info[&100]["AF_EUR"]],
        ["0.1", "0.15"]
    );
    assert_eq!(
        codes(&info, 200),
        ["0", "9", "9"],
        "another ALT or another REF is unknown"
    );
    assert_eq!(recs[&500][4], "C", "the first non-reference ALT");
    assert_eq!(codes(&info, 500), ["1", "0", "9"]);
    assert_eq!(info[&500]["AF_ALL"], "0");
    assert_eq!(recs[&100][0], "chr22", "the chr prefix for the liftover");
}

#[test]
fn a_call_outside_its_mask_or_inside_an_excluded_region_is_masked() {
    let d = tmp("mask");
    let common = table_inputs(&d);
    // 1-based 51-150 and 451-550: 100 and 500, not 200.
    write(&p(&d, "altai.bed"), "22\t50\t150\n22\t450\t550\n");
    write(&p(&d, "vindija.bed"), "22\t150\t250\n"); // holds 200, not 100
    gzip(&p(&d, "denisova.bed.gz"), "chr22\t99\t100\n"); // 100 only, named chr22

    // Control: without --mask the same inputs keep every call.
    let (code, err) = table(&common, &[], &p(&d, "nomask.vcf"));
    assert_eq!(code, 0, "{err}");
    let (nm_rows, nm) = vcf_rows(&p(&d, "nomask.vcf"));
    assert_eq!(nm[&100]["Altai"], "2");
    assert_eq!(nm_rows.keys().copied().collect::<Vec<_>>(), [100, 200, 500]);

    let masks = [
        format!("Altai={}", p(&d, "altai.bed")),
        format!("Vindija={}", p(&d, "vindija.bed")),
        format!("Denisova={}", p(&d, "denisova.bed.gz")),
    ];
    let (code, err) = table(
        &common,
        &[
            "--mask", &masks[0], "--mask", &masks[1], "--mask", &masks[2],
        ],
        &p(&d, "masked.vcf"),
    );
    assert_eq!(code, 0, "{err}");
    let (m_rows, mi) = vcf_rows(&p(&d, "masked.vcf"));
    assert_eq!(
        [&mi[&100]["Altai"], &mi[&100]["Denisova"]],
        ["2", "1"],
        "inside its mask a genome keeps its call; a mask is read with and without chr"
    );
    assert_eq!(
        [&mi[&100]["Vindija"], &nm[&100]["Vindija"]],
        ["10", "0"],
        "outside its mask a call is masked (10 for 0/0)"
    );
    assert_eq!(m_rows.keys().copied().collect::<Vec<_>>(), [100, 200, 500]);
    assert_eq!(
        mi[&200]["Altai"], "10",
        "a row whose only call is masked stays"
    );
    assert_eq!(m_rows[&500][4], "C");
    assert_eq!([&mi[&500]["Altai"], &mi[&500]["Vindija"]], ["1", "10"]);
    assert_eq!(
        [&mi[&200]["Vindija"], &mi[&500]["Denisova"]],
        ["9", "9"],
        "unknown stays 9"
    );

    // The excluded regions mask every call at a site inside one.
    write(&p(&d, "excl.bed"), "22\t90\t120\tsegdup\n");
    let excl = p(&d, "excl.bed");
    let (code, err) = table(&common, &["--exclude", &excl], &p(&d, "excl.vcf"));
    assert_eq!(code, 0, "{err}");
    let (ex_rows, ex) = vcf_rows(&p(&d, "excl.vcf"));
    assert_eq!(codes(&ex, 100), ["12", "10", "11"]);
    assert_eq!(codes(&nm, 100), ["2", "0", "1"]);
    assert_eq!(ex_rows.keys().copied().collect::<Vec<_>>(), [100, 200, 500]);
    assert_eq!([&ex[&200]["Altai"], &ex[&500]["Altai"]], ["0", "1"]);

    write(&p(&d, "excl_other.bed"), "7\t90\t120\tsegdup\n");
    let other = p(&d, "excl_other.bed");
    let (code, err) = table(&common, &["--exclude", &other], &p(&d, "x.vcf"));
    assert!(code == 1 && err.contains("no interval"), "{code} {err}");
    let (code, err) = table(
        &common,
        &["--mask", &masks[0], "--mask", &masks[0]],
        &p(&d, "x.vcf"),
    );
    assert!(code == 2 && err.contains("twice"), "{code} {err}");
    let nobody = format!("Nobody={}", p(&d, "altai.bed"));
    let (code, err) = table(&common, &["--mask", &nobody], &p(&d, "x.vcf"));
    assert_eq!(code, 2, "{err}");
    write(&p(&d, "other.bed"), "7\t50\t150\n");
    let elsewhere = format!("Altai={}", p(&d, "other.bed"));
    let (code, err) = table(&common, &["--mask", &elsewhere], &p(&d, "x.vcf"));
    assert!(
        code == 1 && err.contains("no interval"),
        "a mask with no interval on the chromosome stops: {code} {err}"
    );
}

#[test]
fn strict_marks_rows_outside_the_mask_and_keeps_each_position_once() {
    let d = tmp("strict");
    // Strict mask 1-based 101-300 on chr22.
    let rows = [
        "chr22\t100\tA\tG\t2\t0\t9\t0.5\t0\t0\t0\t0\t0",
        "chr22\t150\tA\tG\t1\t1\t1\t0.5\t0\t0\t0\t0\t0",
        "chr22\t300\tG\tA\t0\t0\t0\t0.5\t0\t0\t0\t0\t0",
        "chr22\t301\tC\tT\t0\t10\t12\t0.5\t0\t0\t0\t0\t0",
        "chr22\t301\tC\tT\t2\t2\t2\t0.5\t0\t0\t0\t0\t0",
    ];
    write(&p(&d, "lifted.tsv"), &(rows.join("\n") + "\n"));
    write(&p(&d, "strict.bed"), "chr22\t100\t300\n");
    let strict = |mask: &str, rows: &str, out: &str| {
        og(&[
            "strict", "--chrom", "22", "--strict", mask, "--rows", rows, "--out", out,
        ])
    };
    let (code, err) = strict(
        &p(&d, "strict.bed"),
        &p(&d, "lifted.tsv"),
        &p(&d, "marked.tsv"),
    );
    assert_eq!(code, 0, "{err}");
    let got: Vec<Vec<String>> = std::fs::read_to_string(p(&d, "marked.tsv"))
        .unwrap()
        .lines()
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect();
    assert_eq!(
        got.iter().map(|c| c[1].as_str()).collect::<Vec<_>>(),
        ["100", "150", "300", "301"]
    );
    let by = |pos: &str| -> Vec<&str> {
        got.iter().find(|c| c[1] == pos).unwrap()[4..7]
            .iter()
            .map(String::as_str)
            .collect()
    };
    assert_eq!(by("100"), ["12", "10", "9"], "outside: 2, 0, unknown");
    assert_eq!(by("150"), ["1", "1", "1"], "inside: untouched");
    assert_eq!(by("300"), ["0", "0", "0"], "the last base inside too");
    assert_eq!(
        by("301"),
        ["10", "10", "12"],
        "a masked call stays, the first of a duplicate wins"
    );
    write(
        &p(&d, "unsorted.tsv"),
        &([rows[0], rows[3], rows[2]].join("\n") + "\n"),
    );
    let (code, err) = strict(
        &p(&d, "strict.bed"),
        &p(&d, "unsorted.tsv"),
        &p(&d, "marked3.tsv"),
    );
    assert!(code == 2 && err.contains("sorted"), "{code} {err}");
    assert!(!d.join("marked3.tsv").exists());
    write(&p(&d, "strict_other.bed"), "chr21\t0\t10\n");
    let (code, err) = strict(
        &p(&d, "strict_other.bed"),
        &p(&d, "lifted.tsv"),
        &p(&d, "marked2.tsv"),
    );
    assert!(code != 0 && err.contains("no interval"), "{code} {err}");
}

/// A genotype file by position.
fn gt_rows(path: &str) -> BTreeMap<u32, Vec<String>> {
    let mut rows = BTreeMap::new();
    for l in gunzip(path).lines().skip(1) {
        let c: Vec<String> = l.split('\t').map(str::to_owned).collect();
        let pos: u32 = c[1].parse().unwrap();
        assert!(
            rows.insert(pos, c).is_none(),
            "two genotype rows for position {pos} in {path}"
        );
    }
    rows
}

const ROWS: [&str; 5] = [
    "chr22\t1000\tA\tG\t2\t2\t0\t0.5\t0.5\t0.5\t0.5\t0.5\t0.5",
    "chr22\t2000\tC\tT\t0\t0\t0\t0.01\t0\t0\t0\t0.02\t0",
    "chr22\t3000\tG\tA\t1\t1\t2\t0\t0\t0\t0\t0\t0",
    "chr22\t4000\tT\tC\t2\t9\t0\t0.2\t0\t0\t0\t0\t0",
    "chr22\t5000\tA\tC\t0\t0\t0\t0.2\t0\t0\t0\t0\t0",
];

/// The table of part 2 and the person's VCF.
fn genotype_inputs(d: &Path) -> String {
    let tdir = p(d, "table");
    std::fs::create_dir_all(&tdir).unwrap();
    gzip(
        &format!("{tdir}/chr22.tsv.gz"),
        &format!("{TABLE_HEADER}\n{}\n", ROWS.join("\n")),
    );
    gzip(
        &format!("{tdir}/chr21.tsv.gz"),
        &format!("{TABLE_HEADER}\nchr21\t1000\tA\tG\t2\t2\t0\t0.5\t0.5\t0.5\t0.5\t0.5\t0.5\n"),
    );
    write(
        &p(d, "s.vcf"),
        &format!(
            "{VCF_HEAD}chr22\t1000\t.\tA\tG\t50\tPASS\t.\tGT:GQ\t0/1:30\n\
             chr22\t3000\t.\tG\tA\t50\tRefCall\t.\tGT\t1/1\n\
             chr22\t5000\t.\tA\tT\t50\tPASS\t.\tGT\t1/1\n"
        ),
    );
    tdir
}

#[test]
fn genotypes_code_your_calls_and_the_stand_ins() {
    let d = tmp("genotypes");
    let tdir = genotype_inputs(&d);
    // A reference block over 1-2500: 2000 was read; 4000 was not.
    write(
        &p(&d, "s.g.vcf"),
        &format!("{VCF_HEAD}{}", ref_block("chr22", 1, 2500)),
    );
    let (code, err) = og(&[
        "genotypes",
        "--table-dir",
        &tdir,
        "--vcf",
        &p(&d, "s.vcf"),
        "--gvcf",
        &p(&d, "s.g.vcf"),
        "--pop",
        "EUR",
        "--pseudo",
        "10",
        "--out-dir",
        &p(&d, "gt"),
    ]);
    assert_eq!(code, 0, "{err}");
    let text = gunzip(&p(&d, "gt/chr22.gt.gz"));
    let head: Vec<&str> = text.lines().next().unwrap().split('\t').collect();
    assert_eq!(
        head[..8],
        [
            "chrom", "pos", "ref", "alt", "Altai", "Vindija", "Denisova", "SAMPLE"
        ]
    );
    assert_eq!(head.len(), 18, "ten stand-ins");
    let g = gt_rows(&p(&d, "gt/chr22.gt.gz"));
    assert_eq!(g[&1000][0], "22", "an integer chromosome for IBDmix");
    assert_eq!(g[&1000][7], "1", "a heterozygous call");
    assert_eq!(g[&2000][7], "0", "no record inside a reference block");
    assert_eq!(g[&4000][7], "9", "no record outside every block");
    assert_eq!(g[&3000][7], "9", "a non-PASS record");
    assert!(!g.contains_key(&5000), "another ALT is left out");
    let sum = |pos: u32| -> u32 { g[&pos][8..].iter().map(|x| x.parse::<u32>().unwrap()).sum() };
    assert_eq!(sum(1000), 10, "0.5 of 20 alleles");
    assert_eq!(sum(2000), 0, "0.02 of 20 rounds to 0");
    assert!(g[&1000][4..].iter().all(|x| x.len() == 1));
    assert!(
        !d.join("gt/chr21.gt.gz").exists(),
        "a chromosome the VCF has no record on is not analysed"
    );
    let spans = std::fs::read_to_string(p(&d, "gt/spans.tsv")).unwrap();
    let span: Vec<&str> = spans.lines().nth(1).unwrap().split('\t').collect();
    assert_eq!(span[..4], ["chr22", "1000", "2000", "2"]);
    let samples = std::fs::read_to_string(p(&d, "gt/samples.txt")).unwrap();
    let samples: Vec<&str> = samples.split_whitespace().collect();
    assert_eq!((samples[0], samples.len()), ("SAMPLE", 11));
    assert_eq!(
        std::fs::read_to_string(p(&d, "gt/chr22.Altai.mask.bed")).unwrap(),
        "",
        "an unmasked table writes empty masks"
    );

    write(
        &p(&d, "two.vcf"),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n",
    );
    let (code, err) = og(&[
        "genotypes",
        "--table-dir",
        &tdir,
        "--vcf",
        &p(&d, "two.vcf"),
        "--out-dir",
        &p(&d, "gt2"),
    ]);
    assert!(code != 0 && err.contains("one-sample"), "{code} {err}");
}

#[test]
fn a_split_vcf_keeps_each_snv_call() {
    let d = tmp("split");
    let sdir = p(&d, "table_split");
    std::fs::create_dir_all(&sdir).unwrap();
    let rows: String = [
        (1000, "A\tG"),
        (2000, "C\tT"),
        (3000, "G\tA"),
        (4000, "T\tC"),
        (5000, "A\tG"),
        (6000, "C\tT"),
        (7000, "G\tA"),
        (8000, "T\tC"),
        (9000, "A\tC"),
    ]
    .iter()
    .map(|(pos, ra)| format!("chr22\t{pos}\t{ra}\t0\t0\t0\t0.2\t0\t0\t0\t0\t0\n"))
    .collect();
    gzip(
        &format!("{sdir}/chr22.tsv.gz"),
        &format!("{TABLE_HEADER}\n{rows}"),
    );
    let recs = [
        "1000\tA\tAT\tPASS\t0/1", // an insertion, then the SNV
        "1000\tA\tG\tPASS\t0/1",
        "2000\tC\tT\tPASS\t1/1", // the SNV, then an insertion
        "2000\tC\tCA\tPASS\t0/1",
        "3000\tG\tA,GT\tPASS\t1/2", // a multi-allelic record, then the SNV
        "3000\tG\tA\tPASS\t0/1",
        "4000\tT\tC\tPASS\t0/1", // two SNVs, both carried
        "4000\tT\tG\tPASS\t0/1",
        "4000\tT\tC\tPASS\t0/1",
        "5000\tA\tT\tPASS\t0/0", // a split 0/1 at A>G,T: 0/0, then 0/1
        "5000\tA\tG\tPASS\t0/1",
        "6000\tC\tT\tPASS\t0/1", // 0/1, then the 0/0 half
        "6000\tC\tG\tPASS\t0/0",
        "7000\tG\tA\tLowQual\t0/1", // filtered, then usable
        "7000\tG\tA\tPASS\t1/1",
        "8000\tT\tA\tPASS\t0/0", // two 0/0 halves: the reference
        "8000\tT\tG\tPASS\t0/0",
    ];
    let body: String = recs
        .iter()
        .map(|r| {
            let c: Vec<&str> = r.split('\t').collect();
            format!(
                "chr22\t{}\t.\t{}\t{}\t50\t{}\t.\tGT\t{}\n",
                c[0], c[1], c[2], c[3], c[4]
            )
        })
        .collect();
    write(&p(&d, "split.vcf"), &format!("{VCF_HEAD}{body}"));
    let (code, err) = og(&[
        "genotypes",
        "--table-dir",
        &sdir,
        "--vcf",
        &p(&d, "split.vcf"),
        "--pseudo",
        "10",
        "--out-dir",
        &p(&d, "gt3"),
    ]);
    assert_eq!(code, 0, "{err}");
    let g = gt_rows(&p(&d, "gt3/chr22.gt.gz"));
    let code = |pos: u32| g.get(&pos).map_or("left out", |c| c[7].as_str());
    assert_eq!(code(1000), "1", "an insertion before the SNV");
    assert_eq!(code(2000), "2", "an insertion after the SNV");
    assert_eq!(code(3000), "1", "a multi-allelic record beside the SNV");
    assert_eq!(code(4000), "9", "two carried SNVs at one position");
    assert_eq!(code(5000), "1", "a 0/0 half before the carried SNV");
    assert_eq!(code(6000), "1", "a 0/0 half after the carried SNV");
    assert_eq!(code(7000), "2", "a usable record wins over a filtered one");
    assert_eq!(code(8000), "0", "two 0/0 halves are the reference");
}

#[test]
fn masked_calls_go_to_ibdmix_with_a_mask_per_genome() {
    let d = tmp("masked");
    let mdir = p(&d, "table_masked");
    std::fs::create_dir_all(&mdir).unwrap();
    let rows: String = [
        (100, "A\tG", "12\t0\t9"),
        (200, "C\tT", "10\t11\t0"),
        (300, "G\tA", "2\t0\t0"),
        (400, "T\tC", "10\t0\t0"),
        (500, "A\tC", "0\t0\t0"),
        (600, "C\tT", "11\t0\t0"),
    ]
    .iter()
    .map(|(pos, ra, c)| format!("chr22\t{pos}\t{ra}\t{c}\t0.2\t0\t0\t0\t0\t0\n"))
    .collect();
    gzip(
        &format!("{mdir}/chr22.tsv.gz"),
        &format!("{TABLE_HEADER}\n{rows}"),
    );
    // Another ALT: row 500 is left out.
    write(
        &p(&d, "m.vcf"),
        &format!("{VCF_HEAD}chr22\t500\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\n"),
    );
    write(
        &p(&d, "m.g.vcf"),
        &format!("{VCF_HEAD}{}", ref_block("chr22", 1, 1000)),
    );
    let (code, err) = og(&[
        "genotypes",
        "--table-dir",
        &mdir,
        "--vcf",
        &p(&d, "m.vcf"),
        "--gvcf",
        &p(&d, "m.g.vcf"),
        "--pseudo",
        "10",
        "--out-dir",
        &p(&d, "gtm"),
    ]);
    assert_eq!(code, 0, "{err}");
    let g = gt_rows(&p(&d, "gtm/chr22.gt.gz"));
    let calls = |pos: u32| g[&pos][4..7].to_vec();
    assert_eq!(calls(100), ["2", "0", "9"], "12 is 2");
    assert_eq!(calls(200), ["0", "1", "0"], "10 and 11 are 0 and 1");
    assert_eq!(calls(600), ["1", "0", "0"]);
    let mask =
        |a: &str| std::fs::read_to_string(p(&d, &format!("gtm/chr22.{a}.mask.bed"))).unwrap();
    assert_eq!(
        mask("Altai"),
        "22\t99\t200\n22\t399\t600\n",
        "the left-out row 500 does not close the run"
    );
    assert_eq!(mask("Vindija"), "22\t199\t200\n");
    assert_eq!(mask("Denisova"), "", "written and empty");

    let bdir = p(&d, "table_bad");
    std::fs::create_dir_all(&bdir).unwrap();
    gzip(
        &format!("{bdir}/chr22.tsv.gz"),
        &format!("{TABLE_HEADER}\nchr22\t100\tA\tG\t7\t0\t0\t0.2\t0\t0\t0\t0\t0\n"),
    );
    let (code, err) = og(&[
        "genotypes",
        "--table-dir",
        &bdir,
        "--vcf",
        &p(&d, "m.vcf"),
        "--out-dir",
        &p(&d, "gtb"),
    ]);
    assert!(code == 1 && err.contains("code '7'"), "{code} {err}");
}

#[test]
fn unsorted_input_stops_and_keeps_the_previous_outputs() {
    let d = tmp("order");
    let tdir = genotype_inputs(&d);
    for (label, records) in [
        (
            "within-table",
            [(3000, "G", "A", "GT"), (1000, "A", "G", "GT")],
        ),
        (
            "past-table",
            [(9000, "A", "G", "GT"), (8000, "A", "G", "GT")],
        ),
        (
            "ignored-call",
            [(3000, "G", "A", "GT"), (1000, "A", ".", "GQ")],
        ),
    ] {
        let bad = p(&d, &format!("{label}.vcf"));
        let body: String = records
            .iter()
            .map(|(pos, r, a, f)| format!("chr22\t{pos}\t.\t{r}\t{a}\t.\tPASS\t.\t{f}\t0/1\n"))
            .collect();
        write(&bad, &format!("{VCF_HEAD}{body}"));
        let out = p(&d, &format!("order-{label}"));
        std::fs::create_dir_all(&out).unwrap();
        let protected: Vec<String> = std::iter::once("chr22.gt.gz".to_owned())
            .chain(ARCHAICS.iter().map(|g| format!("chr22.{g}.mask.bed")))
            .collect();
        for name in &protected {
            write(&format!("{out}/{name}"), "previous output\n");
        }
        let (code, err) = og(&[
            "genotypes",
            "--table-dir",
            &tdir,
            "--vcf",
            &bad,
            "--out-dir",
            &out,
        ]);
        assert!(
            code == 2 && err.contains("not sorted"),
            "{label}: {code} {err}"
        );
        for name in &protected {
            assert_eq!(
                std::fs::read(format!("{out}/{name}")).unwrap(),
                b"previous output\n",
                "{label} keeps {name}"
            );
        }
    }
    let udir = p(&d, "table_unordered");
    std::fs::create_dir_all(&udir).unwrap();
    gzip(
        &format!("{udir}/chr22.tsv.gz"),
        &format!(
            "{TABLE_HEADER}\n{}\n",
            [ROWS[1], ROWS[0], ROWS[2], ROWS[3], ROWS[4]].join("\n")
        ),
    );
    let out = p(&d, "gt_unordered");
    let (code, err) = og(&[
        "genotypes",
        "--table-dir",
        &udir,
        "--vcf",
        &p(&d, "s.vcf"),
        "--out-dir",
        &out,
    ]);
    assert!(code == 2 && err.contains("not sorted"), "{code} {err}");
    let left: Vec<String> = match std::fs::read_dir(&out) {
        Ok(r) => r
            .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
            .collect(),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Vec::new(),
        Err(e) => panic!("cannot read {out}: {e}"),
    };
    assert!(
        !left
            .iter()
            .any(|n| n.ends_with(".gt.gz") || n.ends_with(".mask.bed")),
        "no partial genotype or mask: {left:?}"
    );
}

/// The archaic section og summary collect writes for S1 in `dir`.
fn archaic_section(dir: &Path) -> Value {
    let input = og_summary::collect::Input {
        sample: "S1".to_owned(),
        sample_dir: dir.display().to_string(),
        declared_sex: None,
        assets: repo().join("assets").display().to_string(),
        group: og_summary::collect::Group::All,
        tables: repo().join("bin").display().to_string(),
    };
    let v = og_summary::collect::post::collect(&input).unwrap();
    let get = |v: &Value, k: &str| pycompat::value::get(v, k).unwrap().unwrap().clone();
    get(&get(&v, "sections"), "archaic")
}

/// Both reports of the summary og summary collect writes for S1.
fn reports(dir: &Path) -> (String, String) {
    let input = og_summary::collect::Input {
        sample: "S1".to_owned(),
        sample_dir: dir.display().to_string(),
        declared_sex: None,
        assets: repo().join("assets").display().to_string(),
        group: og_summary::collect::Group::All,
        tables: repo().join("bin").display().to_string(),
    };
    let v = og_summary::collect::post::collect(&input).unwrap();
    let s = og_render::model::Summary::parse(
        &pycompat::json::dumps(&v, &pycompat::json::Dump::default()).unwrap(),
    )
    .unwrap();
    (
        og_render::text::report(&s).unwrap(),
        og_render::html::report(&s).unwrap(),
    )
}

fn field<'a>(v: &'a Value, k: &str) -> &'a Value {
    pycompat::value::get(v, k).unwrap().unwrap()
}

#[test]
fn the_report_merges_segments_and_both_reports_read_it() {
    let d = tmp("report");
    let idir = p(&d, "ibd");
    std::fs::create_dir_all(&idir).unwrap();
    let head = "ID\tchrom\tstart\tend\tslod\n";
    write(
        &format!("{idir}/chr22.Altai.txt"),
        &format!(
            "{head}SAMPLE\t22\t1000000\t1300000\t40\nP001\t22\t5000000\t9000000\t99\n\
             SAMPLE\t22\t2000000\t2040000\t30\nSAMPLE\t22\t3000000\t3100000\t3.5\n"
        ),
    );
    write(
        &format!("{idir}/chr22.Vindija.txt"),
        &format!("{head}SAMPLE\t22\t1200000\t1500000\t50\n"),
    );
    write(
        &format!("{idir}/chr22.Denisova.txt"),
        &format!("{head}SAMPLE\t22\t1400000\t1450000\t12\nSAMPLE\t22\t6000000\t6200000\t20\n"),
    );
    write(
        &p(&d, "spans.tsv"),
        "chrom\tfirst\tlast\tsites\nchr22\t1\t10000000\t5000\n",
    );
    write(
        &p(&d, "genes.gtf"),
        "chr22\tHAVANA\tgene\t1100000\t1110000\t.\t+\t.\tgene_id \"G1\"; gene_type \"protein_coding\"; gene_name \"GENEA\";\n\
         chr22\tHAVANA\tgene\t1120000\t1130000\t.\t+\t.\tgene_id \"G2\"; gene_type \"lncRNA\"; gene_name \"LNC1\";\n\
         chr22\tHAVANA\tgene\t6100000\t6300000\t.\t-\t.\tgene_id \"G3\"; gene_type \"protein_coding\"; gene_name \"GENEB\";\n",
    );
    let (code, err) = og(&[
        "report",
        "--sample",
        "S1",
        "--ibdmix-dir",
        &idir,
        "--spans",
        &p(&d, "spans.tsv"),
        "--genes",
        &p(&d, "genes.gtf"),
        "--out-dir",
        &p(&d, "out"),
    ]);
    assert_eq!(code, 0, "{err}");
    let segs = std::fs::read_to_string(p(&d, "out/S1_archaic_segments.tsv")).unwrap();
    let segs: Vec<Vec<&str>> = segs.lines().map(|l| l.split('\t').collect()).collect();
    assert_eq!(
        segs[0],
        [
            "chrom",
            "start",
            "end",
            "length_bp",
            "ancestry",
            "matched",
            "lod",
            "genes"
        ]
    );
    let body = &segs[1..];
    assert_eq!(body.len(), 2, "{body:?}");
    assert_eq!(
        body[0][..6],
        [
            "chr22",
            "1000000",
            "1500000",
            "500000",
            "Neanderthal",
            "Altai,Vindija,Denisova"
        ],
        "Altai 1.0-1.3 and Vindija 1.2-1.5 Mb merge"
    );
    assert_eq!(body[0][6], "50.0", "the higher LOD");
    assert_eq!(
        body[1][..5],
        ["chr22", "6000000", "6200000", "200000", "Denisovan"],
        "a Denisova stretch inside a Neanderthal one is not Denisovan; one apart is"
    );
    assert_eq!(body[0][7], "GENEA", "protein-coding genes only");
    let summ: BTreeMap<String, String> = std::fs::read_to_string(p(&d, "out/S1_archaic.tsv"))
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| {
            let (k, v) = l.split_once('\t').unwrap();
            (k.to_owned(), v.to_owned())
        })
        .collect();
    assert_eq!(summ["neanderthal_percent"], "5.00");
    assert_eq!(summ["denisovan_percent"], "2.00");
    assert_eq!(
        [&summ["neanderthal_segments"], &summ["denisovan_segments"]],
        ["1", "1"]
    );
    assert!(summ["longest_1"].starts_with("chr22:1000000-1500000 Neanderthal"));
    assert!(summ["longest_1"].contains("GENEA"));
    let txt = std::fs::read_to_string(p(&d, "out/S1_archaic.txt")).unwrap();
    assert!(txt.contains("An estimate, not a measurement."), "{txt}");
    let (code, err) = og(&[
        "report",
        "--sample",
        "S1",
        "--ibdmix-dir",
        &p(&d, "nothing"),
        "--spans",
        &p(&d, "spans.tsv"),
        "--out-dir",
        &p(&d, "out2"),
    ]);
    assert!(
        code != 0 && err.contains("no IBDmix output"),
        "{code} {err}"
    );

    // The summary and the report lines read it.
    let sd = d.join("sample");
    std::fs::create_dir_all(sd.join("archaic")).unwrap();
    let tsv = sd.join("archaic/S1_archaic.tsv");
    std::fs::copy(p(&d, "out/S1_archaic.tsv"), &tsv).unwrap();
    let sec = archaic_section(&sd);
    assert_eq!(field(&sec, "state"), &Value::Str("ok".into()));
    assert_eq!(field(&sec, "step"), &Value::Str("43".into()));
    let data = field(&sec, "data");
    assert_eq!(
        field(data, "neanderthal_percent"),
        &Value::Str("5.00".into())
    );
    assert_eq!(field(data, "denisovan_segments"), &Value::Int(1));
    let Value::Array(longest) = field(data, "longest") else {
        panic!("longest is not a list")
    };
    assert_eq!(longest.len(), 2);
    let lines = og_render::sections::archaic::lines(data).unwrap();
    assert!(
        lines[0].starts_with("Neanderthal: 5.00% (0.5 Mb in 1 segments)"),
        "{lines:?}"
    );
    assert!(lines[1].starts_with("Denisovan: 2.00%"), "{lines:?}");
    let (txt, page) = reports(&sd);
    let shares = "Neanderthal: 5.00% (0.5 Mb in 1 segments)";
    assert!(txt.contains(shares), "the text report: {txt}");
    assert!(
        page.contains("<h2>Archaic Ancestry (Neanderthal, Denisovan)</h2>")
            && page.contains(shares),
        "the HTML report"
    );
    for damaged in [
        "key\tvalue\nneanderthal_percent\tlots\n",
        "key\tvalue\nneanderthal_segments\n",
    ] {
        std::fs::write(&tsv, damaged).unwrap();
        assert_eq!(
            field(&archaic_section(&sd), "state"),
            &Value::Str("unreadable".into()),
            "{damaged:?}"
        );
    }
    let nowhere = archaic_section(&d.join("nowhere"));
    assert_eq!(field(&nowhere, "state"), &Value::Str("missing".into()));
    assert_eq!(field(&nowhere, "data"), &Value::Object(Default::default()));
}

/// Python's `random.Random(seed)` for a small seed: MT19937 seeded by
/// `init_by_array([seed])`, and `random()` from two draws. The IBDmix
/// fixture below is the one the Python test drew.
struct PyRandom {
    mt: [u32; 624],
    i: usize,
}

impl PyRandom {
    fn new(seed: u32) -> PyRandom {
        let mut mt = [0u32; 624];
        mt[0] = 19_650_218;
        for i in 1..624 {
            mt[i] = 1_812_433_253u32
                .wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        let key = [seed];
        let (mut i, mut j) = (1usize, 0usize);
        for _ in 0..624 {
            mt[i] = (mt[i] ^ (mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1_664_525))
                .wrapping_add(key[j])
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= 624 {
                mt[0] = mt[623];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
        }
        for _ in 0..623 {
            mt[i] = (mt[i] ^ (mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1_566_083_941))
                .wrapping_sub(i as u32);
            i += 1;
            if i >= 624 {
                mt[0] = mt[623];
                i = 1;
            }
        }
        mt[0] = 0x8000_0000;
        PyRandom { mt, i: 624 }
    }

    fn next_u32(&mut self) -> u32 {
        if self.i >= 624 {
            for k in 0..624 {
                let y = (self.mt[k] & 0x8000_0000) | (self.mt[(k + 1) % 624] & 0x7fff_ffff);
                let mut v = self.mt[(k + 397) % 624] ^ (y >> 1);
                if y & 1 == 1 {
                    v ^= 0x9908_b0df;
                }
                self.mt[k] = v;
            }
            self.i = 0;
        }
        let mut y = self.mt[self.i];
        self.i += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^ (y >> 18)
    }

    fn random(&mut self) -> f64 {
        let a = f64::from(self.next_u32() >> 5);
        let b = f64::from(self.next_u32() >> 6);
        (a * 67_108_864.0 + b) / 9_007_199_254_740_992.0
    }

    fn uniform(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.random()
    }
}

#[test]
fn the_fixture_generator_draws_what_python_drew() {
    // Python 3: random.Random(11), three random() then uniform(0.01, 0.04).
    let mut r = PyRandom::new(11);
    let got = [r.random(), r.random(), r.random(), r.uniform(0.01, 0.04)];
    assert_eq!(
        got,
        [
            0.452_379_553_509_818_6,
            0.559_772_386_080_496,
            0.924_210_584_023_729_4,
            0.023_969_502_102_993_2
        ]
    );
}

#[cfg(target_arch = "x86_64")]
mod ibdmix {
    //! Part 4: the chain with IBDmix itself.

    use super::*;
    use std::process::Command;

    const NEA: (u32, u32) = (2_000_000, 2_600_000);
    const DEN: (u32, u32) = (3_200_000, 3_600_000);
    const MID: (u32, u32) = (2_200_000, 2_400_000);
    const MASKED: u32 = 10;

    /// Whether the docker tests run here: always in CI.
    fn docker_tests(test: &str) -> bool {
        let on = std::env::var_os("CI").is_some_and(|v| !v.is_empty())
            || std::env::var("OG_DOCKER_TESTS").as_deref() == Ok("1");
        if !on {
            eprintln!(
                "skip: {test} runs IBDmix in docker; set OG_DOCKER_TESTS=1 (CI always runs it)"
            );
        }
        on
    }

    fn image() -> String {
        let text = std::fs::read_to_string(repo().join("versions.env")).unwrap();
        let line = text
            .lines()
            .find(|l| l.starts_with("IBDMIX_IMAGE="))
            .expect("IBDMIX_IMAGE in versions.env");
        let v = &line["IBDMIX_IMAGE=".len()..];
        v.split('#')
            .next()
            .unwrap()
            .trim()
            .trim_matches('"')
            .to_owned()
    }

    fn id(flag: &str) -> String {
        let o = Command::new("id").arg(flag).output().unwrap();
        String::from_utf8(o.stdout).unwrap().trim().to_owned()
    }

    fn open_dir(name: &str) -> Tmp {
        let d = tmp(name);
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&d, std::fs::Permissions::from_mode(0o777)).unwrap();
        d
    }

    /// A 4 Mb chr22 table and one person, drawn as the Python test drew
    /// them. A third of the sites carry an allele the Neanderthals share
    /// and modern people rarely do, a third one only the Denisovan carries,
    /// a third are common. Inside a planted stretch one copy carries the
    /// archaic allele. `mid`: the middle of the planted Neanderthal stretch
    /// holds sites the person contradicts, with the Neanderthals' calls
    /// written masked (10 plus the call) or unknown (9).
    fn synthetic(d: &Path, plant: bool, mid: Option<&str>) {
        let mut rnd = PyRandom::new(11);
        std::fs::create_dir_all(d.join("table")).unwrap();
        let (mut rows, mut calls) = (Vec::new(), Vec::new());
        for (i, pos) in (1_000_000u32..5_000_000).step_by(1000).enumerate() {
            let (f, al, vi, de) = match i % 3 {
                0 => (rnd.uniform(0.01, 0.04), 2, 2, 0),
                1 => (rnd.uniform(0.01, 0.04), 0, 0, 2),
                _ => {
                    let f = rnd.uniform(0.1, 0.9);
                    let n = u32::from(rnd.random() < f) + u32::from(rnd.random() < f);
                    (f, n, n, n)
                }
            };
            let (mut ta, mut tv) = (al, vi);
            let mut g = [rnd.random() < f, rnd.random() < f];
            if plant && NEA.0 <= pos && pos < NEA.1 {
                g[0] = al == 2;
            }
            if let Some(m) = mid
                && MID.0 <= pos
                && pos < MID.1
            {
                if i.is_multiple_of(3) {
                    g = [false, false];
                }
                (ta, tv) = if m == "masked" {
                    (al + MASKED, vi + MASKED)
                } else {
                    (9, 9)
                };
            }
            let fs = format!("{f:.4}");
            rows.push(format!(
                "chr22\t{pos}\tA\tG\t{ta}\t{tv}\t{de}\t{}",
                [fs.as_str(); 6].join("\t")
            ));
            if plant && DEN.0 <= pos && pos < DEN.1 {
                g[0] = de == 2;
            }
            if g.iter().any(|x| *x) {
                let gt = if g.iter().all(|x| *x) { "1/1" } else { "0/1" };
                calls.push(format!("chr22\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT\t{gt}"));
            }
        }
        gzip(
            &p(d, "table/chr22.tsv.gz"),
            &format!("{TABLE_HEADER}\n{}\n", rows.join("\n")),
        );
        write(&p(d, "s.vcf"), &format!("{VCF_HEAD}{}\n", calls.join("\n")));
    }

    /// genotypes, IBDmix with the flags of scripts/43-archaic.sh, report;
    /// the segments as (start, end, length, ancestry).
    fn chain(d: &Path, gvcf: Option<String>) -> Vec<(u32, u32, u32, String)> {
        let gvcf =
            gvcf.unwrap_or_else(|| format!("{VCF_HEAD}{}", ref_block("chr22", 1, 6_000_000)));
        write(&p(d, "s.g.vcf"), &gvcf);
        let (code, err) = og(&[
            "genotypes",
            "--table-dir",
            &p(d, "table"),
            "--vcf",
            &p(d, "s.vcf"),
            "--gvcf",
            &p(d, "s.g.vcf"),
            "--out-dir",
            &p(d, "gt"),
        ]);
        assert_eq!(code, 0, "{err}");
        std::fs::create_dir_all(d.join("ibd")).unwrap();
        let sh = "set -e; gzip -dc /w/gt/chr22.gt.gz > /w/gt/chr22.gt; for a in Altai Vindija Denisova; do \
                  ibdmix -g /w/gt/chr22.gt -s /w/gt/samples.txt -n $a -r /w/gt/chr22.$a.mask.bed -o /w/ibd/chr22.$a.txt \
                  -d 4 -m 1 -a 0.01 -e 0.002 -c 2; done";
        let o = Command::new("docker")
            .args(["run", "--rm", "--network", "none", "-u"])
            .arg(format!("{}:{}", id("-u"), id("-g")))
            .arg("-v")
            .arg(format!("{}:/w", d.display()))
            .arg(image())
            .args(["sh", "-c", sh])
            .output()
            .unwrap_or_else(|e| panic!("docker: {e}"));
        assert!(
            o.status.success(),
            "IBDmix: {}{}",
            String::from_utf8_lossy(&o.stderr),
            String::from_utf8_lossy(&o.stdout)
        );
        let (code, err) = og(&[
            "report",
            "--sample",
            "S",
            "--ibdmix-dir",
            &p(d, "ibd"),
            "--spans",
            &p(d, "gt/spans.tsv"),
            "--out-dir",
            &p(d, "out"),
        ]);
        assert_eq!(code, 0, "{err}");
        let body: Vec<(u32, u32, u32, String)> =
            std::fs::read_to_string(p(d, "out/S_archaic_segments.tsv"))
                .unwrap()
                .lines()
                .skip(1)
                .map(|l| {
                    let c: Vec<&str> = l.split('\t').collect();
                    (
                        c[1].parse().unwrap(),
                        c[2].parse().unwrap(),
                        c[3].parse().unwrap(),
                        c[4].to_owned(),
                    )
                })
                .collect();
        eprintln!("segments in {}: {body:?}", d.display());
        body
    }

    fn cover(segs: &[(u32, u32, u32, String)], kind: &str, (lo, hi): (u32, u32)) -> u32 {
        segs.iter()
            .filter(|s| s.3 == kind)
            .map(|s| s.1.min(hi).saturating_sub(s.0.max(lo)))
            .sum()
    }

    #[test]
    fn planted_stretches_come_out_as_one_segment_each() {
        if !docker_tests("planted_stretches_come_out_as_one_segment_each") {
            return;
        }
        let d = open_dir("real_planted");
        synthetic(&d, true, None);
        let body = chain(&d, None);
        let nea: Vec<_> = body
            .iter()
            .filter(|s| s.3 == "Neanderthal")
            .cloned()
            .collect();
        let den: Vec<_> = body
            .iter()
            .filter(|s| s.3 == "Denisovan")
            .cloned()
            .collect();
        assert_eq!(nea.len(), 1, "{body:?}");
        assert!(
            cover(&nea, "Neanderthal", NEA) * 10 >= 9 * (NEA.1 - NEA.0),
            "{nea:?}"
        );
        assert!(
            nea[0].2 <= NEA.1 - NEA.0 + 100_000,
            "no more than 100 kb beyond it"
        );
        assert_eq!(den.len(), 1, "{body:?}");
        assert!(
            cover(&den, "Denisovan", DEN) * 10 >= 9 * (DEN.1 - DEN.0),
            "{den:?}"
        );

        let d = open_dir("real_control");
        synthetic(&d, false, None);
        assert!(chain(&d, None).is_empty(), "without the planted stretches");
    }

    #[test]
    fn masked_mismatches_stop_a_segment_that_unknown_calls_let_through() {
        if !docker_tests("masked_mismatches_stop_a_segment_that_unknown_calls_let_through") {
            return;
        }
        let span = MID.1 - MID.0;
        let d = open_dir("real_mid_unknown");
        synthetic(&d, true, Some("unknown"));
        let unknown = cover(&chain(&d, None), "Neanderthal", MID);
        assert!(
            unknown * 10 >= 9 * span,
            "control: written as unknown the segment runs through ({unknown})"
        );
        let d = open_dir("real_mid_masked");
        synthetic(&d, true, Some("masked"));
        let masked = cover(&chain(&d, None), "Neanderthal", MID);
        assert!(
            masked * 10 <= span,
            "written as masked no segment covers the middle ({masked})"
        );
    }

    #[test]
    fn a_stretch_with_no_call_needs_a_reference_block() {
        if !docker_tests("a_stretch_with_no_call_needs_a_reference_block") {
            return;
        }
        let (lo, hi) = (200_000u32, 349_000u32);
        let mut got = Vec::new();
        for covered in [true, false] {
            let d = open_dir(&format!("real_gap_{covered}"));
            std::fs::create_dir_all(d.join("table")).unwrap();
            let rows: String = (lo..=hi)
                .step_by(1000)
                .map(|pos| format!("chr22\t{pos}\tC\tT\t0\t0\t0\t{}\n", ["0.9"; 6].join("\t")))
                .collect();
            gzip(
                &p(&d, "table/chr22.tsv.gz"),
                &format!("{TABLE_HEADER}\n{rows}"),
            );
            write(
                &p(&d, "s.vcf"),
                &format!("{VCF_HEAD}chr22\t10000000\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"),
            );
            // Without the block the person has one read site left, at the start.
            let block = if covered {
                ref_block("chr22", 1, 400_000)
            } else {
                ref_block("chr22", 1, lo)
            };
            let body = chain(&d, Some(format!("{VCF_HEAD}{block}")));
            got.push(
                body.iter()
                    .filter(|s| s.3 == "Neanderthal")
                    .map(|s| s.2)
                    .sum::<u32>(),
            );
        }
        assert!(
            got[0] >= 100_000,
            "control: a stretch read as the reference gives a segment ({got:?})"
        );
        assert_eq!(got[1], 0, "no call and no reference block: no segment");
    }
}
