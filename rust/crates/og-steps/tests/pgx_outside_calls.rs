//! `og pgx outside-calls`: what reaches PharmCAT as an outside call, and
//! `og pgx pypgx-summary`, step 32's summary.
//!
//! CYP2D6 reaches PharmCAT only when pypgx and a second caller (Cyrius) give
//! the same diplotype and the depth check passed; HLA-A and HLA-B come from
//! T1K. Each case runs the verb on synthetic inputs (`tests/fixtures/pgx/`)
//! and compares both files it writes, byte for byte:
//!
//! 1. disagree: pypgx and Cyrius differ, no outside call, 'indeterminate';
//! 2. one missing: pypgx only (Cyrius not run), no outside call;
//! 3. agree: both give the same diplotype, depth ok, the call is passed;
//! 4. depth: both agree, but the depth check found multi-mapped reads (4b:
//!    no depth check at all);
//! 5. no call: Cyrius ran but made no call, one caller only;
//! 6. HLA: T1K's HLA-A and HLA-B, truncated to two fields;
//! 7. HLA edges: one allele (read as homozygous); an allele of quality 0.
//!
//! What step 27 (`og pgx cpic-report --consensus`) makes of these tables
//! (cases 8 to 10) is tested in `og-pgx`
//! (`tests/port_tests_pgx_outside.rs`), on the tables
//! [`consensus_tables_the_cpic_report_tests_read`] holds to this verb.
//!
//! Step 32's summary and the probable CYP2D6 call:
//!
//! 11. the star alleles of a diplotype and the defining-variant evidence of
//!     one call are the crate's own unit tests (`pypgx_summary::tests`);
//! 13. pypgx alone with the depth check ok and PASS defining variants: the
//!     outside-call file stays empty and step 36's row stays indeterminate
//!     (the reports' probable call itself is `og summary collect`'s rule,
//!     tested in `og-summary`'s `collect/pgx.rs`);
//! 14. the whole verb on pypgx results archives, pypgx's own tables (rows
//!     of pypgx 0.26.0's gene and allele tables,
//!     `tests/fixtures/pgx/summary/`) and an indexed VCF: GSTT1 not
//!     callable on this reference (not FAILED), a gene without results
//!     FAILED, pypgx's CNV call for a BAM gene only, the evidence columns
//!     per gene, byte for byte what the Python it replaced wrote with the
//!     real pypgx and pysam; the tables found beside the pypgx on PATH;
//! 15. both files or neither: a consensus table that cannot be written
//!     leaves no outside-call file behind.
//!
//! The PATH cases run the verb in a child copy of this test binary
//! ([`og_child`]), so the child's PATH is its own.

use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use og_steps::pgx_rules::outside_calls_cli;

const HEADER: &str = "Gene\tResult\tOutside_call\tReason\tEvidence\n";
const NO_HLA: &str = "HLA-A\tnot typed\tno\tHLA typing (step 08) did not run\t-\n\
                      HLA-B\tnot typed\tno\tHLA typing (step 08) did not run\t-\n";
const MULTIMAPPED: &str = "CYP2D6 copy number unreliable: reads are multi-mapped \
                           (was this BAM aligned to a reference with ALT contigs?)";
const NO_CYP2D6: &str = "CYP2D6\tindeterminate\tno\tno caller made a call\
                         \tpypgx: not run; Cyrius: not run; depth check: not run\n";

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../tests/fixtures/pgx")
}

fn fx(name: &str) -> String {
    fixtures().join(name).to_string_lossy().into_owned()
}

/// A scratch directory of one test: emptied when made, removed when
/// dropped.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Self {
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("pgx-outside-{name}"));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).to_string_lossy().into_owned()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// The two files `og pgx outside-calls` writes for these inputs (option,
/// fixture name), or the exit code and stderr when it failed.
fn run(tag: &str, inputs: &[(&str, &str)]) -> (String, String) {
    let tmp = Scratch::new(tag);
    let (calls, cons) = (tmp.path("calls.tsv"), tmp.path("consensus.tsv"));
    let mut args = vec![
        "--calls".to_owned(),
        calls.clone(),
        "--consensus".to_owned(),
        cons.clone(),
    ];
    for (opt, name) in inputs {
        args.push(format!("--{opt}"));
        args.push(fx(name));
    }
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = outside_calls_cli(&args, &mut out, &mut err);
    if code != 0 {
        let e = String::from_utf8_lossy(&err);
        return (format!("exit {code}: {}", e.trim()), String::new());
    }
    (
        fs::read_to_string(calls).unwrap(),
        fs::read_to_string(cons).unwrap(),
    )
}

fn case(desc: &str, want_calls: &str, want_consensus: &str, inputs: &[(&str, &str)]) {
    let (calls, cons) = run(&desc.replace(' ', "_"), inputs);
    assert_eq!(calls, want_calls, "{desc}: outside calls");
    assert_eq!(cons, want_consensus, "{desc}: consensus");
}

const PYPGX_12: (&str, &str) = ("pypgx", "pypgx_cyp2d6_1_2.tsv");
const CYRIUS_12: (&str, &str) = ("cyrius", "cyrius_1_2.tsv");
const DEPTH_OK: (&str, &str) = ("depth-check", "depth_check_ok.tsv");

/// Written first: a naive "take pypgx" passes pypgx's call on and fails
/// here.
#[test]
fn case1_pypgx_and_cyrius_disagree() {
    case(
        "1 disagree",
        "",
        &format!(
            "{HEADER}{NO_HLA}CYP2D6\tindeterminate\tno\tpypgx and Cyrius disagree\
             \tpypgx: *1/*4; Cyrius: *1/*2; depth check: ok\n"
        ),
        &[("pypgx", "pypgx_cyp2d6_1_4.tsv"), CYRIUS_12, DEPTH_OK],
    );
}

#[test]
fn case2_one_caller_missing() {
    case(
        "2 one missing",
        "",
        &format!(
            "{HEADER}{NO_HLA}CYP2D6\tindeterminate\tno\tone caller only: a second caller must agree (Cyrius: setup.sh --cyrius)\
             \tpypgx: *1/*2; Cyrius: not run; depth check: ok\n"
        ),
        &[PYPGX_12, DEPTH_OK],
    );
}

#[test]
fn case3_both_agree_and_the_depth_check_passed() {
    case(
        "3 agree",
        "CYP2D6\t*1/*2\n",
        &format!(
            "{HEADER}{NO_HLA}CYP2D6\t*1/*2\tyes\tpypgx and Cyrius agree; the depth check passed\
             \tpypgx: *1/*2; Cyrius: *1/*2; depth check: ok\n"
        ),
        &[PYPGX_12, CYRIUS_12, DEPTH_OK],
    );
}

#[test]
fn case4_depth_unreliable_or_not_run() {
    case(
        "4 depth unreliable",
        "",
        &format!(
            "{HEADER}{NO_HLA}CYP2D6\tindeterminate\tno\t{MULTIMAPPED}\
             \tpypgx: *1/*2; Cyrius: *1/*2; depth check: unreliable\n"
        ),
        &[
            PYPGX_12,
            CYRIUS_12,
            ("depth-check", "depth_check_multimapped.tsv"),
        ],
    );
    case(
        "4b no depth check",
        "",
        &format!(
            "{HEADER}{NO_HLA}CYP2D6\tindeterminate\tno\tthe CYP2D6 depth check did not run\
             \tpypgx: *1/*2; Cyrius: *1/*2; depth check: not run\n"
        ),
        &[PYPGX_12, CYRIUS_12],
    );
}

#[test]
fn case5_cyrius_made_no_call() {
    case(
        "5 Cyrius made no call",
        "",
        &format!(
            "{HEADER}{NO_HLA}CYP2D6\tindeterminate\tno\tone caller only: a second caller must agree (Cyrius: setup.sh --cyrius)\
             \tpypgx: *1/*2; Cyrius: no call (None, Not_assigned_to_haplotypes); depth check: ok\n"
        ),
        &[PYPGX_12, ("cyrius", "cyrius_none.tsv"), DEPTH_OK],
    );
}

#[test]
fn case6_hla_from_t1k() {
    case(
        "6 HLA from T1K",
        "HLA-A\t*02:01/*24:02\nHLA-B\t*07:02/*44:02\n",
        &format!(
            "{HEADER}HLA-A\t*02:01/*24:02\tyes\tT1K (step 08)\tHLA-A*02:01:01 (quality 60); HLA-A*24:02:01 (quality 60)\n\
             HLA-B\t*07:02/*44:02\tyes\tT1K (step 08)\tHLA-B*07:02:01 (quality 60); HLA-B*44:02:01 (quality 48)\n\
             {NO_CYP2D6}"
        ),
        &[("hla", "t1k_genotype.tsv")],
    );
}

#[test]
fn case7_hla_one_allele_and_an_allele_of_quality_0() {
    case(
        "7 HLA one allele, and an allele of quality 0",
        "HLA-A\t*02:01/*02:01\n",
        &format!(
            "{HEADER}HLA-A\t*02:01/*02:01\tyes\tT1K (step 08): one allele, read as homozygous\tHLA-A*02:01:01 (quality 60)\n\
             HLA-B\tindeterminate\tno\tT1K gives an allele quality 0 or below (T1K says to ignore it)\
             \tHLA-B*07:02:01 (quality 60); HLA-B*44:02:01 (quality 0)\n\
             {NO_CYP2D6}"
        ),
        &[("hla", "t1k_genotype_homozygous_lowqual.tsv")],
    );
}

/// The consensus tables `og-pgx`'s `tests/port_tests_pgx_outside.rs` hands
/// to `og pgx cpic-report` (cases 8 to 10), as literal text there: they
/// must be what this verb writes for those inputs.
#[test]
fn consensus_tables_the_cpic_report_tests_read() {
    let held = format!(
        "{HEADER}{NO_HLA}CYP2D6\tindeterminate\tno\tone caller only: a second caller must agree (Cyrius: setup.sh --cyrius)\
         \tpypgx: *1/*4; Cyrius: not run; depth check: ok\n"
    );
    case(
        "8 CYP2D6 held back",
        "",
        &held,
        &[("pypgx", "pypgx_cyp2d6_1_4.tsv"), DEPTH_OK],
    );
    case(
        "9 HLA of PharmCAT's example",
        "HLA-A\t*02:01/*24:02\nHLA-B\t*15:02/*57:01\n",
        &format!(
            "{HEADER}HLA-A\t*02:01/*24:02\tyes\tT1K (step 08)\tHLA-A*02:01:01 (quality 60); HLA-A*24:02:01 (quality 60)\n\
             HLA-B\t*15:02/*57:01\tyes\tT1K (step 08)\tHLA-B*15:02:01 (quality 60); HLA-B*57:01:01 (quality 48)\n\
             {NO_CYP2D6}"
        ),
        &[("hla", "t1k_genotype_docs_example.tsv")],
    );
}

const SUMMARY_HEADER: [&str; 11] = [
    "Gene",
    "Diplotype",
    "Phenotype",
    "CNV_call",
    "Source",
    "Status",
    "Defining_variants",
    "Min_GQ",
    "Weakest_variant",
    "Defining_filter",
    "Quality",
];

#[test]
fn case13_pypgx_alone_passes_nothing() {
    let tmp = Scratch::new("13");
    let summ = tmp.path("pypgx.tsv");
    fs::write(
        &summ,
        format!(
            "{}\nCYP2D6\t*1/*2\tNormal Metabolizer\tNormal\tBAM\tcalled\t1/1 in the VCF\t50\
             \tchr22:42126611 C>G GQ 50 AD 20,22\tPASS\tok\n",
            SUMMARY_HEADER.join("\t")
        ),
    )
    .unwrap();
    let (calls, cons) = (tmp.path("c.tsv"), tmp.path("k.tsv"));
    let args: Vec<String> = [
        "--calls",
        &calls,
        "--consensus",
        &cons,
        "--pypgx",
        &summ,
        "--depth-check",
        &fx("depth_check_ok.tsv"),
    ]
    .map(str::to_owned)
    .to_vec();
    assert_eq!(
        outside_calls_cli(&args, &mut Vec::new(), &mut Vec::new()),
        0
    );
    assert_eq!(
        fs::read_to_string(&calls).unwrap(),
        "",
        "13 the outside-call file stays empty for a probable call"
    );
    let cons = fs::read_to_string(&cons).unwrap();
    let last: Vec<&str> = cons.lines().last().unwrap().split('\t').take(4).collect();
    assert_eq!(
        last,
        [
            "CYP2D6",
            "indeterminate",
            "no",
            "one caller only: a second caller must agree (Cyrius: setup.sh --cyrius)"
        ],
        "13 step 36's row stays indeterminate, one caller only"
    );
}

// --- 14: the verb on pypgx's results archives and tables and an indexed VCF

/// `sample.vcf.gz` (bgzip and tabix of bcftools 1.24) holds these records.
const SUMMARY_VCF: &str = "##fileformat=VCFv4.2
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr22,length=50818468>
##FILTER=<ID=RefCall,Description=\"ref\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"q\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
chr10\t94761900\t.\tC\tT\t50\tPASS\t.\tGT:GQ:AD\t0/1:60:0,31
chr10\t94775360\t.\tAGGGGGGGG\tA\t50\tPASS\t.\tGT:GQ:AD\t0/1:40:10,10
chr10\t94775367\t.\tA\tG\t50\tPASS\t.\tGT:GQ:AD\t0/1:48:20,19
chr10\t94781859\t.\tG\tA,C\t50\tPASS\t.\tGT:GQ:AD\t0/1:12:26,3,.
chr22\t42126611\t.\tC\tG\t50\tRefCall\t.\tGT:GQ:AD\t0/0:50:20,22
";

/// What `bin/pgx_outside_calls.py pypgx-summary` wrote for these inputs in
/// the pypgx 0.26.0 image, with the real pypgx and pysam (the parent of the
/// commit that deleted it).
fn summary_want() -> Vec<Vec<String>> {
    let rows: [&[&str]; 7] = [
        &[
            "GSTM1",
            "*0/*0",
            "no pypgx phenotype rule",
            "WholeDel2",
            "BAM",
            "called",
            "none",
            "N/A",
            "N/A",
            "none",
            "N/A",
        ],
        &[
            "GSTT1",
            "N/A",
            "N/A",
            "N/A",
            "BAM",
            "not callable on this reference",
            "N/A",
            "N/A",
            "N/A",
            "N/A",
            "N/A",
        ],
        &[
            "CYP2D6",
            "*1/*2x2",
            "Normal Metabolizer",
            "Duplication",
            "BAM",
            "called",
            "1/2 in the VCF",
            "50",
            "chr22:42126611 C>G GQ 50 AD 20,22",
            "not PASS: chr22:42126611 C>G RefCall; not in the VCF: chr22:42127941 G>A",
            "ok",
        ],
        &[
            "CYP2B6", "FAILED", "N/A", "N/A", "BAM", "failed", "N/A", "N/A", "N/A", "N/A", "N/A",
        ],
        &[
            "CYP2C19",
            "*2/*17",
            "Intermediate Metabolizer",
            "N/A",
            "VCF",
            "called",
            "3/3 in the VCF",
            "12",
            "chr10:94781859 G>A GQ 12 AD 26,3,.",
            "PASS",
            "low-confidence: check the reads",
        ],
        &[
            "CYP2C9",
            "Indeterminate",
            "Indeterminate",
            "N/A",
            "VCF",
            "no call",
            "N/A",
            "N/A",
            "N/A",
            "N/A",
            "N/A",
        ],
        &[
            "NUDT15",
            "*1/*99",
            "N/A",
            "N/A",
            "VCF",
            "called",
            "unknown",
            "N/A",
            "N/A",
            "not in pypgx's table: *99",
            "N/A",
        ],
    ];
    std::iter::once(&SUMMARY_HEADER[..])
        .chain(rows)
        .map(|r| r.iter().map(|s| (*s).to_owned()).collect())
        .collect()
}

const SUMMARY_LOG: &str = "Genes called: 4/7\nWARNING: star alleles resting on a variant below GQ 20, check the reads: \
                           CYP2C19 *2/*17 (chr10:94781859 G>A GQ 12 AD 26,3,.)\n";

/// A zip file of `(name, content)` members, deflated, as Python's
/// `zipfile.ZipFile(..., ZIP_DEFLATED)` writes it.
fn write_zip(path: &Path, members: &[(&str, &str)]) {
    use flate2::Compression;
    use flate2::write::DeflateEncoder;
    let mut out = Vec::new();
    let mut central = Vec::new();
    for (name, body) in members {
        let mut crc = flate2::Crc::new();
        crc.update(body.as_bytes());
        let mut enc = DeflateEncoder::new(Vec::new(), Compression::default());
        enc.write_all(body.as_bytes()).unwrap();
        let data = enc.finish().unwrap();
        let offset = u32::try_from(out.len()).unwrap();
        let sizes = [
            crc.sum(),
            u32::try_from(data.len()).unwrap(),
            u32::try_from(body.len()).unwrap(),
        ];
        let name_len = u16::try_from(name.len()).unwrap();
        // Local header: version 20, no flags, deflated, 1980-01-01 00:00.
        out.extend_from_slice(&0x0403_4b50_u32.to_le_bytes());
        for v in [20_u16, 0, 8, 0, 0x21] {
            out.extend_from_slice(&v.to_le_bytes());
        }
        for v in sizes {
            out.extend_from_slice(&v.to_le_bytes());
        }
        out.extend_from_slice(&name_len.to_le_bytes());
        out.extend_from_slice(&0_u16.to_le_bytes());
        out.extend_from_slice(name.as_bytes());
        out.extend_from_slice(&data);
        central.extend_from_slice(&0x0201_4b50_u32.to_le_bytes());
        for v in [20_u16, 20, 0, 8, 0, 0x21] {
            central.extend_from_slice(&v.to_le_bytes());
        }
        for v in sizes {
            central.extend_from_slice(&v.to_le_bytes());
        }
        for v in [name_len, 0, 0, 0, 0] {
            central.extend_from_slice(&v.to_le_bytes());
        }
        central.extend_from_slice(&0_u32.to_le_bytes());
        central.extend_from_slice(&offset.to_le_bytes());
        central.extend_from_slice(name.as_bytes());
    }
    let cd_offset = u32::try_from(out.len()).unwrap();
    let count = u16::try_from(members.len()).unwrap();
    let cd_len = u32::try_from(central.len()).unwrap();
    out.extend_from_slice(&central);
    out.extend_from_slice(&0x0605_4b50_u32.to_le_bytes());
    for v in [0_u16, 0, count, count] {
        out.extend_from_slice(&v.to_le_bytes());
    }
    out.extend_from_slice(&cd_len.to_le_bytes());
    out.extend_from_slice(&cd_offset.to_le_bytes());
    out.extend_from_slice(&0_u16.to_le_bytes());
    fs::create_dir_all(path.parent().unwrap()).unwrap();
    fs::write(path, out).unwrap();
}

/// A pypgx results archive: `metadata.txt` and the sample table `data.tsv`
/// in one folder, as pypgx's `Archive.to_file` writes it.
fn results_zip(path: &Path, genotype: &str, phenotype: &str, cnv: &str) {
    write_zip(
        path,
        &[
            (
                "abc123/metadata.txt",
                "Gene=X\nAssembly=GRCh38\nSemanticType=SampleTable[Results]\n",
            ),
            (
                "abc123/data.tsv",
                &format!(
                    "\tVariantData\tCNV\tCandidates\tGenotype\tPhenotype\n\
                     S\tx\t{cnv}\t*1\t{genotype}\t{phenotype}\n"
                ),
            ),
        ],
    );
}

const CHILD: &str = "OG_PORT_TESTS_PGX_CHILD";
const MARK: &str = "\u{1e}og-child\u{1e}";

/// Not a test of its own: the child [`og`] starts runs
/// `og pgx pypgx-summary` here with the argv it was handed, then exits with
/// the verb's code. In the normal test run it returns at once.
#[test]
fn og_child() {
    let Some(argv) = std::env::var_os(CHILD) else {
        return;
    };
    let argv: Vec<String> = argv
        .to_string_lossy()
        .split('\u{1f}')
        .map(str::to_owned)
        .collect();
    let (mut out, mut err) = (std::io::stdout().lock(), std::io::stderr().lock());
    let _ = write!(out, "{MARK}");
    let _ = write!(err, "{MARK}");
    let code = og_steps::pypgx_summary::cli(&argv, &mut out, &mut err);
    let _ = out.flush();
    let _ = err.flush();
    std::process::exit(i32::from(code));
}

/// `og pgx pypgx-summary ARGV` in a child copy of this test binary, with
/// `path` as its PATH: (exit code, stdout, stderr).
fn og(argv: &[String], path: &str) -> (i32, String, String) {
    let child = Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "og_child",
            "--nocapture",
            "--test-threads=1",
            "-q",
        ])
        .env(CHILD, argv.join("\u{1f}"))
        .env("PATH", path)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let out = child.wait_with_output().unwrap();
    let code = out.status.code().unwrap_or(-1);
    let (so, se) = (
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr),
    );
    let after = |s: &str| match s.split_once(MARK) {
        Some((_, rest)) => rest.to_owned(),
        None => panic!("the child never reached the verb: {s}"),
    };
    (code, after(&so), after(&se))
}

/// The verb on the five results archives: (exit code, stdout, stderr,
/// rows of the summary).
fn summary(
    tmp: &Scratch,
    extra: &[&str],
    path: Option<&str>,
) -> (i32, String, String, Vec<Vec<String>>) {
    let outbase = tmp.0.join("out");
    let out = tmp.path("summary.tsv");
    for (gene, g, p, c) in [
        ("CYP2C19", "*2/*17", "Intermediate Metabolizer", "Normal"),
        ("GSTM1", "*0/*0", "Indeterminate", "WholeDel2"),
        ("CYP2D6", "*1/*2x2", "Normal Metabolizer", "Duplication"),
        ("CYP2C9", "Indeterminate", "Indeterminate", "NA"),
        ("NUDT15", "*1/*99", "", "Normal"),
    ] {
        results_zip(&outbase.join(gene).join("results.zip"), g, p, c);
    }
    let mut argv: Vec<String> = [
        "--outbase",
        &outbase.to_string_lossy(),
        "--bam-genes",
        "GSTM1 GSTT1 CYP2D6 CYP2B6",
        "--vcf-genes",
        "CYP2C19 CYP2C9 NUDT15",
        "--vcf",
        &fx("summary/sample.vcf.gz"),
        "--not-callable",
        "GSTT1",
        "--out",
        &out,
    ]
    .map(str::to_owned)
    .to_vec();
    argv.extend(extra.iter().map(|s| (*s).to_owned()));
    let (code, so, se) = match path {
        Some(p) => og(&argv, p),
        None => {
            let (mut so, mut se) = (Vec::new(), Vec::new());
            let code = og_steps::pypgx_summary::cli(&argv, &mut so, &mut se);
            (
                i32::from(code),
                String::from_utf8(so).unwrap(),
                String::from_utf8(se).unwrap(),
            )
        }
    };
    let rows = fs::read_to_string(&out)
        .map(|t| {
            t.lines()
                .map(|l| l.split('\t').map(str::to_owned).collect())
                .collect()
        })
        .unwrap_or_default();
    (code, so, se, rows)
}

#[test]
fn case14_the_fixture_vcf_holds_the_records_the_test_names() {
    let raw = fs::read(fx("summary/sample.vcf.gz")).unwrap();
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(&raw[..])
        .read_to_string(&mut text)
        .unwrap();
    assert_eq!(
        text, SUMMARY_VCF,
        "14 the fixture VCF holds the records the test names"
    );
}

#[test]
fn case14_the_summary_row_by_row_as_the_python_wrote_it() {
    let tmp = Scratch::new("14a");
    let data = fx("summary");
    let (code, so, se, rows) = summary(&tmp, &["--pypgx-data", &data], None);
    assert_eq!(
        (code, se.as_str()),
        (0, ""),
        "14 the subcommand exits 0, without a warning"
    );
    assert_eq!(
        rows,
        summary_want(),
        "14 the summary, row by row, as the Python wrote it"
    );
    assert_eq!(so, SUMMARY_LOG, "14 the log names the low-confidence call");
}

#[test]
fn case14_pypgx_tables_beside_the_pypgx_on_path() {
    // pypgx's tables found beside the pypgx on PATH, for the Python its #!
    // line names (the image also has a python3.1 link to python3.12).
    let tmp = Scratch::new("14b");
    let prefix = tmp.0.join("prefix");
    let data = prefix.join("lib/python3.12/site-packages/pypgx/api/data");
    fs::create_dir_all(&data).unwrap();
    for name in ["gene-table.csv", "allele-table.csv"] {
        fs::copy(fixtures().join("summary").join(name), data.join(name)).unwrap();
    }
    std::os::unix::fs::symlink("python3.12", prefix.join("lib/python3.1")).unwrap();
    fs::create_dir_all(prefix.join("bin")).unwrap();
    let pypgx = prefix.join("bin/pypgx");
    fs::write(&pypgx, "#!/usr/local/bin/python3.12\n").unwrap();
    {
        use std::os::unix::fs::PermissionsExt;
        fs::set_permissions(&pypgx, fs::Permissions::from_mode(0o755)).unwrap();
    }
    let path = format!(
        "{}:{}",
        prefix.join("bin").display(),
        std::env::var("PATH").unwrap_or_default()
    );
    let (code, _, _, rows) = summary(&tmp, &[], Some(&path));
    assert_eq!(
        (code, rows),
        (0, summary_want()),
        "14 pypgx's tables are found beside the pypgx on PATH"
    );
}

#[test]
fn case14_without_pypgx_tables_the_evidence_says_not_looked_up() {
    let tmp = Scratch::new("14c");
    let (code, _, se, rows) = summary(&tmp, &[], Some("/usr/bin:/bin"));
    let evidence = rows.get(3).map(|r| r[6..].to_vec());
    assert_eq!(
        (code, evidence, se.contains("WARNING: no pypgx tables")),
        (
            0,
            Some(
                ["not looked up", "N/A", "N/A", "N/A", "N/A"]
                    .map(str::to_owned)
                    .to_vec()
            ),
            true
        ),
        "14 without pypgx's tables the evidence says not looked up, and why"
    );
}

#[test]
fn case15_both_files_or_neither() {
    let tmp = Scratch::new("15");
    let calls = tmp.path("calls.tsv");
    let args: Vec<String> = [
        "--calls",
        &calls,
        "--consensus",
        &tmp.path("no-such-dir/consensus.tsv"),
        "--hla",
        &fx("t1k_genotype.tsv"),
    ]
    .map(str::to_owned)
    .to_vec();
    let code = outside_calls_cli(&args, &mut Vec::new(), &mut Vec::new());
    let left: Vec<_> = fs::read_dir(&tmp.0).unwrap().collect();
    assert_eq!(
        (code, Path::new(&calls).exists(), left.len()),
        (1, false, 0),
        "15 the consensus cannot be written: exit 1, and no outside-call file"
    );
}
