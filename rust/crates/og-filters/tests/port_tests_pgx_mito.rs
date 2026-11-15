//! `og pgx mito-merge`: the MT-RNR1 outside call PharmCAT gets from the
//! chrM calls of step 20 / MITO_VARIANTS (step 07 and PHARMCAT run it).
//!
//! The test builds a stand-in for `pharmcat.jar`, a zip holding
//! `org/pharmgkb/pharmcat/phenotype/MT_RNR1.json` with some of PharmCAT
//! 3.4.0's MT-RNR1 alleles and their function values, writes invented
//! Mutect2-style chrM VCFs (GT 0/1 at every AF, as Mutect2 writes them), runs
//! the verb in a child copy of this test binary ([`og_child`]) so its log on
//! stdout and stderr is read as the `og` binary writes it, and checks:
//!
//! 1. a homoplasmic PASS m.1555A>G (AF 0.99) is the call; step 36's lines
//!    are copied before it and an MT-RNR1 line there is dropped;
//! 2. AF 0.95 counts as homoplasmic, 0.949 does not; `--min-af` moves the
//!    bound;
//! 3. a heteroplasmic, a filtered or an AF-less m.1555A>G passes nothing and
//!    says why; so do a homoplasmic PASS change at a named position that
//!    PharmCAT does not name and a PASS deletion near m.961; a filtered one,
//!    or a low-AF unnamed change, does not;
//! 4. PASS chrM records with no named allele give Reference when samtools
//!    depth shows 20 reads or more at every named position; without a depth
//!    file, with one position below 20 (`--min-depth` moves the bound) or
//!    absent from it, no call; a homoplasmic named allele needs no depth; a
//!    file with no PASS chrM record gives no call; records on other contigs
//!    are ignored;
//! 5. two homoplasmic named alleles: the increased-risk one is the call, the
//!    other is listed;
//! 6. INFO AF when FORMAT has none, the AF of the second ALT of a
//!    multi-allelic record, a bgzip/gzip input and the contig name MT are
//!    read;
//! 7. a jar without the MT-RNR1 definitions, or an unreadable VCF, exits 2.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const INC: &str = "Increased risk of aminoglycoside-induced hearing loss";
const UNC: &str = "Uncertain risk of aminoglycoside-induced hearing loss";
const NOR: &str = "Normal risk of aminoglycoside-induced hearing loss";
const ALLELES: [(&str, &str); 8] = [
    ("Reference", NOR),
    ("m.827A>G", NOR),
    ("m.960C>del", UNC),
    ("m.961T>del", UNC),
    ("m.961T>del+Cn", UNC),
    ("m.1095T>C", INC),
    ("m.1189T>C", UNC),
    ("m.1555A>G", INC),
];
const HEADER: &str = "##fileformat=VCFv4.2\n##contig=<ID=chrM,length=16569>\n\
                      #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

const CHILD: &str = "OG_PORT_TESTS_PGX_CHILD";
const MARK: &str = "\u{1e}og-child\u{1e}";

/// Not a test of its own: the child [`og`] starts runs `og pgx mito-merge`
/// here with the argv it was handed, then exits with the verb's code. In
/// the normal test run it returns at once.
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
    print!("{MARK}");
    eprint!("{MARK}");
    let code = og_filters::mito::cli(&argv);
    let _ = std::io::stdout().flush();
    let _ = std::io::stderr().flush();
    std::process::exit(i32::from(code));
}

/// `og pgx mito-merge ARGV` in a child copy of this test binary: (exit
/// code, stdout and stderr).
fn og(argv: &[String]) -> (i32, String) {
    let child = Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "og_child",
            "--nocapture",
            "--test-threads=1",
            "-q",
        ])
        .env(CHILD, argv.join("\u{1f}"))
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let out = child.wait_with_output().unwrap();
    let after = |b: &[u8]| {
        let s = String::from_utf8_lossy(b);
        match s.split_once(MARK) {
            Some((_, rest)) => rest.to_owned(),
            None => panic!("the child never reached the verb: {s}"),
        }
    };
    (
        out.status.code().unwrap_or(-1),
        after(&out.stdout) + &after(&out.stderr),
    )
}

/// A zip file of stored members, as Python's `zipfile.ZipFile(path, "w")`
/// writes it.
fn write_zip(path: &Path, members: &[(&str, &str)]) {
    let mut out = Vec::new();
    let mut central = Vec::new();
    for (name, body) in members {
        let mut crc = flate2::Crc::new();
        crc.update(body.as_bytes());
        let offset = u32::try_from(out.len()).unwrap();
        let size = u32::try_from(body.len()).unwrap();
        let name_len = u16::try_from(name.len()).unwrap();
        // Local header: version 20, no flags, stored, 1980-01-01 00:00.
        out.extend_from_slice(&0x0403_4b50_u32.to_le_bytes());
        for v in [20_u16, 0, 0, 0, 0x21] {
            out.extend_from_slice(&v.to_le_bytes());
        }
        for v in [crc.sum(), size, size] {
            out.extend_from_slice(&v.to_le_bytes());
        }
        out.extend_from_slice(&name_len.to_le_bytes());
        out.extend_from_slice(&0_u16.to_le_bytes());
        out.extend_from_slice(name.as_bytes());
        out.extend_from_slice(body.as_bytes());
        central.extend_from_slice(&0x0201_4b50_u32.to_le_bytes());
        for v in [20_u16, 20, 0, 0, 0, 0x21] {
            central.extend_from_slice(&v.to_le_bytes());
        }
        for v in [crc.sum(), size, size] {
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
    fs::write(path, out).unwrap();
}

/// One Mutect2-style chrM record.
fn rec(pos: u32, r: &str, alt: &str, af: Option<f64>, filt: &str, chrom: &str) -> String {
    let (fmt, val) = match af {
        Some(af) => ("GT:AF:DP", format!("0/1:{af}:3000")),
        None => ("GT:DP", "0/1:3000".to_owned()),
    };
    format!("{chrom}\t{pos}\t.\t{r}\t{alt}\t.\t{filt}\t.\t{fmt}\t{val}\n")
}

fn pass(pos: u32, r: &str, alt: &str, af: f64) -> String {
    rec(pos, r, alt, Some(af), "PASS", "chrM")
}

fn background() -> String {
    pass(263, "A", "G", 0.998) + &pass(750, "A", "G", 0.999)
}

/// The scratch folder of one test, with its stand-in jar.
struct Case {
    dir: PathBuf,
    jar: String,
    n: usize,
}

/// What one run gave: exit code, the log, the lines of the output file
/// (`None` when it was not written).
struct Out {
    rc: i32,
    log: String,
    lines: Option<Vec<String>>,
}

impl Out {
    fn mt(&self) -> Vec<&str> {
        self.lines
            .iter()
            .flatten()
            .filter_map(|l| l.strip_prefix("MT-RNR1\t"))
            .collect()
    }
}

#[derive(Default)]
struct Run<'a> {
    outside: Option<&'a str>,
    gz: bool,
    extra: &'a [&'a str],
    jar: Option<String>,
    /// Depth at the positions 648-1601: 500, or the value given here
    /// (`None`: the position is left out).
    depth: Option<&'a [(u32, Option<u32>)]>,
}

impl Case {
    fn new(name: &str) -> Self {
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("pgx-mito-{name}"));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let named: Vec<String> = ALLELES
            .iter()
            .map(|(n, f)| {
                format!(
                    "{{\"name\": \"{n}\", \"functionValue\": \"{f}\", \"activityValue\": null}}"
                )
            })
            .collect();
        let json = format!(
            "{{\"gene\": \"MT-RNR1\", \"namedAlleles\": [{}]}}",
            named.join(", ")
        );
        let jar = dir.join("pharmcat.jar");
        write_zip(
            &jar,
            &[("org/pharmgkb/pharmcat/phenotype/MT_RNR1.json", &json)],
        );
        Self {
            jar: jar.to_string_lossy().into_owned(),
            dir,
            n: 0,
        }
    }

    fn path(&self, name: &str) -> String {
        self.dir.join(name).to_string_lossy().into_owned()
    }

    fn run(&mut self, records: &str, o: &Run) -> Out {
        self.n += 1;
        let n = self.n;
        let mut args: Vec<String> = Vec::new();
        let vcf = self.path(&format!("m{n}.vcf{}", if o.gz { ".gz" } else { "" }));
        let text = format!("{HEADER}{records}");
        if o.gz {
            let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
            enc.write_all(text.as_bytes()).unwrap();
            fs::write(&vcf, enc.finish().unwrap()).unwrap();
        } else {
            fs::write(&vcf, text).unwrap();
        }
        let out = self.path(&format!("o{n}.tsv"));
        args.extend(
            [
                "--jar",
                o.jar.as_deref().unwrap_or(&self.jar),
                "--mito",
                &vcf,
                "--out",
                &out,
            ]
            .map(str::to_owned),
        );
        args.extend(o.extra.iter().map(|s| (*s).to_owned()));
        if let Some(depth) = o.depth {
            let path = self.path(&format!("d{n}.depth"));
            let mut t = String::new();
            for p in 648..1602 {
                let d = depth
                    .iter()
                    .find(|(q, _)| *q == p)
                    .map_or(Some(500), |(_, d)| *d);
                if let Some(d) = d {
                    t.push_str(&format!("chrM\t{p}\t{d}\n"));
                }
            }
            fs::write(&path, t).unwrap();
            args.push("--depth".to_owned());
            args.push(path);
        }
        if let Some(text) = o.outside {
            let path = self.path(&format!("c{n}.tsv"));
            fs::write(&path, text).unwrap();
            args.push("--outside".to_owned());
            args.push(path);
        }
        let (rc, log) = og(&args);
        let lines = fs::read_to_string(&out)
            .ok()
            .map(|t| t.lines().map(str::to_owned).collect());
        Out { rc, log, lines }
    }
}

impl Drop for Case {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.dir);
    }
}

#[test]
fn case1_a_homoplasmic_m1555ag_is_the_call() {
    let mut c = Case::new("1");
    let r = c.run(
        &(background() + &pass(1555, "A", "G", 0.99)),
        &Run {
            outside: Some("HLA-A\t*01:01/*26:01\nMT-RNR1\tm.1095T>C\n"),
            ..Run::default()
        },
    );
    assert_eq!(r.rc, 0, "1: exit 0: {}", r.log);
    assert_eq!(r.mt(), ["m.1555A>G"], "1: MT-RNR1 m.1555A>G is the call");
    assert_eq!(
        r.lines.as_ref().and_then(|l| l.first()).map(String::as_str),
        Some("HLA-A\t*01:01/*26:01"),
        "1: step 36's HLA-A line is copied first"
    );
    assert_eq!(
        r.mt().len(),
        1,
        "1: the MT-RNR1 line of step 36's file is dropped"
    );
    assert!(
        r.log.contains(&format!(
            "MT-RNR1 m.1555A>G passed to PharmCAT as an outside call (homoplasmic, chrM:1555 A>G AF 0.990; {INC}"
        )),
        "1: the log says it passed m.1555A>G, homoplasmic, with the risk: {}",
        r.log
    );
    assert!(
        r.log.contains("dropped the MT-RNR1 line"),
        "1: the log names the dropped line: {}",
        r.log
    );
}

#[test]
fn case2_the_af_bound() {
    let mut c = Case::new("2");
    let r = c.run(
        &(background() + &pass(1555, "A", "G", 0.95)),
        &Run::default(),
    );
    assert_eq!(
        r.mt(),
        ["m.1555A>G"],
        "2: AF 0.95 is homoplasmic: {}",
        r.log
    );
    let r = c.run(
        &(background() + &pass(1555, "A", "G", 0.949)),
        &Run::default(),
    );
    assert!(r.mt().is_empty(), "2: AF 0.949 passes nothing: {}", r.log);
    let r = c.run(
        &(background() + &pass(1555, "A", "G", 0.9)),
        &Run {
            extra: &["--min-af", "0.85"],
            ..Run::default()
        },
    );
    assert_eq!(
        r.mt(),
        ["m.1555A>G"],
        "2: --min-af 0.85 takes AF 0.9: {}",
        r.log
    );
}

#[test]
fn case3_what_passes_nothing() {
    let mut c = Case::new("3");
    let cases = [
        (
            "heteroplasmic",
            pass(1555, "A", "G", 0.4),
            "is heteroplasmic (AF below 0.95)",
        ),
        (
            "filtered",
            rec(1555, "A", "G", Some(0.99), "possible_numt", "chrM"),
            "is filtered (possible_numt)",
        ),
        (
            "without an allele fraction",
            rec(1555, "A", "G", None, "PASS", "chrM"),
            "has no allele fraction",
        ),
        (
            "a change PharmCAT does not name",
            pass(1555, "A", "C", 0.99),
            "a change at an MT-RNR1 position that PharmCAT does not name",
        ),
        (
            "a deletion near m.961",
            pass(955, "AC", "A", 0.99),
            "at or near m.960",
        ),
        (
            "an insertion at m.961",
            pass(961, "T", "TC", 0.99),
            "(chrM:961 T>TC AF 0.990) at or near m.960",
        ),
    ];
    for (desc, records, want) in cases {
        let r = c.run(
            &(background() + &records),
            &Run {
                outside: Some("HLA-B\t*35:01/*44:02\n"),
                ..Run::default()
            },
        );
        assert_eq!(r.rc, 0, "3: {desc}: exit 0: {}", r.log);
        assert!(
            r.mt().is_empty(),
            "3: {desc}: no MT-RNR1 line: {:?}",
            r.lines
        );
        assert_eq!(
            r.lines,
            Some(vec!["HLA-B\t*35:01/*44:02".to_owned()]),
            "3: {desc}: step 36's line is still copied"
        );
        assert!(
            r.log.contains(want),
            "3: {desc}: the log says why ({want}): {}",
            r.log
        );
        assert!(
            r.log.contains("it stays No Result"),
            "3: {desc}: the log says the gene stays No Result: {}",
            r.log
        );
    }
    // m.960 is named only as a deletion: an unnamed SNV there still holds
    // Reference back.
    let r = c.run(
        &(background() + &pass(960, "C", "T", 0.99)),
        &Run {
            depth: Some(&[]),
            ..Run::default()
        },
    );
    assert!(
        r.mt().is_empty(),
        "3: an unnamed SNV at m.960 (named only as a deletion): no MT-RNR1 line: {}",
        r.log
    );
    assert!(
        r.log
            .contains("a change at an MT-RNR1 position that PharmCAT does not name"),
        "3: an unnamed SNV at m.960: the log says why: {}",
        r.log
    );
    let r = c.run(
        &(background()
            + &rec(1555, "A", "C", Some(0.3), "weak_evidence", "chrM")
            + &rec(955, "AC", "A", Some(0.2), "weak_evidence", "chrM")),
        &Run {
            depth: Some(&[]),
            ..Run::default()
        },
    );
    assert_eq!(
        r.mt(),
        ["Reference"],
        "3: a filtered unnamed change and a filtered deletion near m.961 hold nothing back: {}",
        r.log
    );
}

#[test]
fn case4_reference_and_no_pass_record() {
    let mut c = Case::new("4");
    let r = c.run(
        &(background()
            + &pass(2706, "A", "G", 0.999)
            + &rec(1555, "A", "G", Some(0.99), "PASS", "chr1")),
        &Run {
            depth: Some(&[]),
            ..Run::default()
        },
    );
    assert_eq!(
        r.mt(),
        ["Reference"],
        "4: no named allele on chrM, depth 500: Reference: {}",
        r.log
    );
    assert!(
        r.log.contains(
            "all 6 named positions have depth 20 or more (lowest 500); 3 PASS chrM records read"
        ),
        "4: the log says all 6 named positions are covered and counts the PASS records: {}",
        r.log
    );
    let r = c.run(
        &background(),
        &Run {
            depth: Some(&[(1555, Some(20)), (961, Some(20))]),
            ..Run::default()
        },
    );
    assert_eq!(
        r.mt(),
        ["Reference"],
        "4: depth exactly 20 counts: {}",
        r.log
    );
    let r = c.run(&background(), &Run::default());
    assert!(r.mt().is_empty(), "4: no depth file: no call: {}", r.log);
    assert!(
        r.log.contains("there is no read depth") && r.log.contains("stays No Result"),
        "4: no depth file: the log says so: {}",
        r.log
    );
    let r = c.run(
        &background(),
        &Run {
            depth: Some(&[(1555, Some(5))]),
            ..Run::default()
        },
    );
    assert!(
        r.mt().is_empty(),
        "4: chrM:1555 at depth 5: no call: {}",
        r.log
    );
    assert!(
        r.log
            .contains("1 of its 6 positions have depth below 20 (chrM:1555 depth 5)"),
        "4: chrM:1555 at depth 5: the log names it: {}",
        r.log
    );
    let r = c.run(
        &background(),
        &Run {
            depth: Some(&[(1555, Some(5))]),
            extra: &["--min-depth", "5"],
            ..Run::default()
        },
    );
    assert_eq!(
        r.mt(),
        ["Reference"],
        "4: --min-depth 5 takes depth 5: {}",
        r.log
    );
    let r = c.run(
        &background(),
        &Run {
            depth: Some(&[(960, None)]),
            ..Run::default()
        },
    );
    assert!(
        r.mt().is_empty() && r.log.contains("chrM:960 depth 0"),
        "4: a named position absent from the depth file counts as depth 0: {}",
        r.log
    );
    let r = c.run(
        &(background() + &pass(1555, "A", "G", 0.99)),
        &Run::default(),
    );
    assert_eq!(
        r.mt(),
        ["m.1555A>G"],
        "4: a homoplasmic named allele needs no depth file: {}",
        r.log
    );
    let r = c.run(
        &rec(263, "A", "G", Some(0.99), "weak_evidence", "chrM"),
        &Run::default(),
    );
    assert_eq!(r.rc, 0, "4: no PASS chrM record: exit 0: {}", r.log);
    assert!(r.mt().is_empty(), "4: no PASS chrM record: no call");
    assert!(
        r.log.contains("the chrM call set has no PASS record"),
        "4: no PASS chrM record: the log says so: {}",
        r.log
    );
    assert_eq!(
        r.lines,
        Some(Vec::new()),
        "4: an empty output file is still written"
    );
}

#[test]
fn case5_two_homoplasmic_named_alleles() {
    let mut c = Case::new("5");
    let r = c.run(
        &(background() + &pass(1189, "T", "C", 0.99) + &pass(1555, "A", "G", 0.98)),
        &Run::default(),
    );
    assert_eq!(
        r.mt(),
        ["m.1555A>G"],
        "5: the increased-risk allele is the call: {}",
        r.log
    );
    assert!(
        r.log.contains("also homoplasmic, not passed") && r.log.contains("m.1189T>C"),
        "5: the other is listed: {}",
        r.log
    );
}

#[test]
fn case6_af_sources_multi_allelic_gzip_and_mt() {
    let mut c = Case::new("6");
    let r = c.run(
        &(background() + "chrM\t1555\t.\tA\tG\t.\tPASS\tAF=0.97\tGT\t0/1\n"),
        &Run::default(),
    );
    assert_eq!(
        r.mt(),
        ["m.1555A>G"],
        "6: INFO AF is read when FORMAT has none: {}",
        r.log
    );
    let r = c.run(
        &(background() + "chrM\t1555\t.\tA\tC,G\t.\tPASS\t.\tGT:AF\t0/1/2:0.01,0.98\n"),
        &Run::default(),
    );
    assert_eq!(
        r.mt(),
        ["m.1555A>G"],
        "6: the AF of the second ALT: {}",
        r.log
    );
    let r = c.run(
        &(background() + &pass(1555, "A", "G", 0.99)),
        &Run {
            gz: true,
            ..Run::default()
        },
    );
    assert_eq!(r.mt(), ["m.1555A>G"], "6: a gzip VCF: {}", r.log);
    let r = c.run(
        &(rec(263, "A", "G", Some(0.99), "PASS", "MT")
            + &rec(1555, "A", "G", Some(0.99), "PASS", "MT")),
        &Run::default(),
    );
    assert_eq!(r.mt(), ["m.1555A>G"], "6: contig MT: {}", r.log);
}

#[test]
fn case7_errors_exit_2() {
    let mut c = Case::new("7");
    let bad = c.dir.join("empty.jar");
    write_zip(&bad, &[("logback.xml", "<configuration/>")]);
    let r = c.run(
        &background(),
        &Run {
            jar: Some(bad.to_string_lossy().into_owned()),
            ..Run::default()
        },
    );
    assert_eq!(r.rc, 2, "7: a jar without MT_RNR1.json exits 2: {}", r.log);
    assert!(
        r.log
            .contains("org/pharmgkb/pharmcat/phenotype/MT_RNR1.json"),
        "7: and names the resource: {}",
        r.log
    );
    let r = c.run("chrM\t1555\tA\n", &Run::default());
    assert_eq!(
        r.rc, 2,
        "7: a record with too few columns exits 2: {}",
        r.log
    );
}
