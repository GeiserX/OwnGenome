//! `og viral report` (step 47): the viral DNA in the reads that do not map
//! to the human genome. Moved from tests/test_viral_report.py (its cases 1
//! to 7); the summary's reading is in og-summary's `tests/optin_summary.rs`
//! and the report lines in og-render's `tests/optin_report.rs`. The read
//! selection itself is in `i411_viral_calls.rs`.
//!
//! Each case writes a viral FASTA, the `samtools coverage` and `depth -aa`
//! files with chosen depths, the alignments that cover the same bases, a
//! human coverage file and `samtools idxstats`, runs the verb and reads the
//! rows it writes.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

const HEADER: &str = "virus\taccession\tlength\treads\tbreadth\tmean_depth\tmode_depth\thuman_depth\trelative_depth\tcopies_per_cell\tcall\tnote";
const ICI: &str = "inherited integrated HHV-6 likely";
const CONTROL: &str = "lab control (phiX)";
const FOUND: &str = "viral DNA found";
const FEW: &str = "a few reads";
const NONE: &str = "not found";
const BELOW: &str = "below threshold";
const COV_HEAD: &str =
    "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n";

const A: &str = "NC_001664.4";
const B: &str = "NC_000898.1";
const EBV: &str = "NC_007605.1";
const HCV: &str = "NC_004102.1";
const PHIX: &str = "NC_001422.1";
const GENOMES: [(&str, &str, usize); 5] = [
    (
        A,
        "Human betaherpesvirus 6A, variant A DNA, complete virion genome, isolate U1102",
        2000,
    ),
    (B, "Human herpesvirus 6B, complete genome", 2000),
    (EBV, "Human gammaherpesvirus 4, complete genome", 1500),
    (HCV, "Hepatitis C virus genotype 1, complete genome", 900),
    (PHIX, "Escherichia phage phiX174, complete genome", 500),
];

type Depth = Box<dyn Fn(usize) -> u32>;

fn flat(x: u32) -> Depth {
    Box::new(move |_| x)
}

/// A fixed pseudo-random sequence: not low complexity.
fn sequence(len: usize, seed: &mut u32) -> String {
    (0..len)
        .map(|_| {
            *seed = seed.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            char::from(b"ACGT"[(*seed >> 16) as usize % 4])
        })
        .collect()
}

/// `n` SAM records (MAPQ 60) over the covered positions: the positions in
/// `n` chunks in order, each chunk one record whose CIGAR skips (N) the
/// gaps; past one read per position, the rest on the reverse strand.
fn hits(acc: &str, covered: &[usize], n: usize, seed: &mut u32) -> String {
    let mut out = String::new();
    let k = covered.len();
    if k == 0 || n == 0 {
        return out;
    }
    let chunks: Vec<Vec<usize>> = if n <= k {
        (0..n)
            .map(|i| covered[i * k / n..(i + 1) * k / n].to_vec())
            .collect()
    } else {
        let mut c: Vec<Vec<usize>> = Vec::new();
        for _ in 0..n / k {
            c.extend(covered.iter().map(|p| vec![*p]));
        }
        c.extend(covered[..n % k].iter().map(|p| vec![*p]));
        c
    };
    for (i, c) in chunks.iter().enumerate() {
        let (mut cigar, mut run, mut prev) = (String::new(), 0, None::<usize>);
        for &p in c {
            if let Some(q) = prev
                && p != q + 1
            {
                let _ = write!(cigar, "{run}M{}N", p - q - 1);
                run = 0;
            }
            run += 1;
            prev = Some(p);
        }
        let _ = write!(cigar, "{run}M");
        let flag = if i >= k { 16 } else { 0 };
        let _ = writeln!(
            out,
            "r{i}\t{flag}\t{acc}\t{}\t60\t{cigar}\t*\t0\t0\t{}\t*",
            c[0],
            sequence(c.len(), seed)
        );
    }
    out
}

/// One case's inputs.
#[derive(Default)]
struct Case {
    depths: HashMap<&'static str, Depth>,
    reads: HashMap<&'static str, usize>,
    /// (rname, length, covbases, meandepth) of the human coverage; chr1
    /// and chr22 at 30x when `None`.
    human: Option<Vec<(&'static str, u64, u64, f64)>>,
}

fn case(depths: Vec<(&'static str, Depth)>, reads: &[(&'static str, usize)]) -> Case {
    Case {
        depths: depths.into_iter().collect(),
        reads: reads.iter().copied().collect(),
        human: None,
    }
}

struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Dir {
        let d = std::env::temp_dir().join(format!("og-optin-viral-{}-{tag}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        Dir(d)
    }
    fn p(&self, n: &str) -> String {
        self.0.join(n).to_str().unwrap().to_owned()
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn write_inputs(d: &Dir, c: &Case) {
    let mut seed = 47u32;
    let (mut fa, mut cov, mut dep, mut sam) = (
        String::new(),
        COV_HEAD.to_owned(),
        String::new(),
        String::new(),
    );
    for (acc, title, len) in GENOMES {
        let _ = write!(fa, ">{acc} {title}\n{}\n\n", "A".repeat(len));
        let depths: Vec<u32> = (1..=len)
            .map(|p| c.depths.get(acc).map_or(0, |f| f(p)))
            .collect();
        let covered: Vec<usize> = (1..=len).filter(|p| depths[p - 1] > 0).collect();
        let n = c.reads.get(acc).copied().unwrap_or(0);
        sam += &hits(acc, &covered, n, &mut seed);
        let mean = depths.iter().map(|x| f64::from(*x)).sum::<f64>() / len as f64;
        let _ = writeln!(
            cov,
            "{acc}\t1\t{len}\t{n}\t{}\t{}\t{mean}\t30\t60",
            covered.len(),
            100.0 * covered.len() as f64 / len as f64
        );
        for (i, x) in depths.iter().enumerate() {
            let _ = writeln!(dep, "{acc}\t{}\t{x}", i + 1);
        }
    }
    let human = c.human.clone().unwrap_or_else(|| {
        vec![
            ("chr1", 100_000, 100_000, 30.0),
            ("chr22", 50_000, 50_000, 30.0),
        ]
    });
    let mut hcov = COV_HEAD.to_owned();
    for (r, len, covbases, mean) in human {
        let _ = writeln!(
            hcov,
            "{r}\t1\t{len}\t1000\t{covbases}\t{}\t{mean}\t30\t60",
            100.0 * covbases as f64 / len as f64
        );
    }
    for (n, t) in [
        ("v.fa", fa),
        ("vcov.tsv", cov),
        ("vdepth.tsv", dep),
        ("hits.sam", sam),
        ("hcov.tsv", hcov),
        (
            "idx.tsv",
            "chr1\t100000\t1000\t2\n*\t0\t0\t1000\n".to_owned(),
        ),
    ] {
        fs::write(d.p(n), t).unwrap();
    }
}

/// `og viral report` on the inputs in `d`, into `out`: exit code and stderr.
fn og(d: &Dir, out: &str) -> (u8, String) {
    let args: Vec<String> = [
        "--sample",
        "S",
        "--viral-fasta",
        &d.p("v.fa"),
        "--viral-coverage",
        &d.p("vcov.tsv"),
        "--viral-depth",
        &d.p("vdepth.tsv"),
        "--human-coverage",
        &d.p("hcov.tsv"),
        "--viral-hits",
        &d.p("hits.sam"),
        "--idxstats",
        &d.p("idx.tsv"),
        "--out",
        out,
    ]
    .iter()
    .map(|s| (*s).to_owned())
    .collect();
    let (mut o, mut e) = (Vec::new(), Vec::new());
    let code = og_reports::viral::cli(Some("report"), &args, &mut o, &mut e);
    (code, String::from_utf8(e).unwrap())
}

type Rows = Vec<(String, HashMap<String, String>)>;

/// The rows by accession, in the table's order, or og's stderr.
fn run_raw(tag: &str, c: &Case, damage: impl Fn(&Dir)) -> Result<Rows, String> {
    let d = Dir::new(tag);
    write_inputs(&d, c);
    damage(&d);
    let out = d.p("out.tsv");
    let (code, err) = og(&d, &out);
    if code != 0 {
        return Err(err);
    }
    let text = fs::read_to_string(&out).unwrap();
    let mut lines = text.lines();
    let head = lines.next().unwrap();
    assert_eq!(head, HEADER, "the header is the documented one");
    let head: Vec<&str> = head.split('\t').collect();
    Ok(lines
        .map(|l| {
            let r: HashMap<String, String> = head
                .iter()
                .zip(l.split('\t'))
                .map(|(k, v)| ((*k).to_owned(), v.to_owned()))
                .collect();
            (r["accession"].clone(), r)
        })
        .collect())
}

fn run(tag: &str, c: &Case) -> HashMap<String, HashMap<String, String>> {
    run_raw(tag, c, |_| {}).unwrap().into_iter().collect()
}

#[test]
fn integration_at_half_the_human_depth_and_at_the_threshold() {
    let ici = run_raw("ici", &case(vec![(B, flat(15))], &[(B, 400)]), |_| {}).unwrap();
    let order: Vec<&str> = ici.iter().map(|r| r.0.as_str()).collect();
    assert_eq!(order, GENOMES.map(|g| g.0), "rows in FASTA order");
    let ici: HashMap<_, _> = ici.into_iter().collect();
    let b = &ici[B];
    assert_eq!(
        [
            b["call"].as_str(),
            &b["relative_depth"],
            &b["copies_per_cell"],
            &b["breadth"],
            &b["virus"]
        ],
        [ICI, "0.500", "1.0000", "1.0000", "HHV-6B"]
    );
    assert_eq!(ici[A]["call"], NONE);
    assert!(b["note"].contains("Liu et al. 2020"));
    let above = run("above", &case(vec![(B, flat(9))], &[(B, 300)]));
    assert_eq!(
        (
            above[B]["relative_depth"].as_str(),
            above[B]["call"].as_str()
        ),
        ("0.300", ICI)
    );
    let below = run("below", &case(vec![(B, flat(8))], &[(B, 300)]));
    assert_eq!(
        (
            below[B]["relative_depth"].as_str(),
            below[B]["call"].as_str()
        ),
        ("0.267", FOUND)
    );
    assert!(below[B]["note"].contains("Below the integration threshold"));
    // The mode, not the mean: a pile-up on a short stretch is no integration.
    let pile = run(
        "pile",
        &case(
            vec![(B, Box::new(|p| if p <= 100 { 400 } else { 0 }))],
            &[(B, 300)],
        ),
    );
    assert_eq!(
        (pile[B]["mode_depth"].as_str(), pile[B]["call"].as_str()),
        ("0", FEW)
    );
    assert_eq!(pile[B]["mean_depth"], "20.000");
}

#[test]
fn low_level_few_reads_and_none() {
    let low = run(
        "low",
        &case(
            vec![(EBV, Box::new(|p| u32::from(p % 2 == 1)))],
            &[(EBV, 40)],
        ),
    );
    assert_eq!(low[EBV]["call"], FOUND);
    assert!(low[EBV]["relative_depth"].parse::<f64>().unwrap() < 0.1);
    let few = run(
        "few",
        &case(vec![(EBV, Box::new(|p| u32::from(p <= 100)))], &[(EBV, 12)]),
    );
    assert_eq!(few[EBV]["call"], FEW);
    let nine = run(
        "nine",
        &case(vec![(EBV, Box::new(|p| u32::from(p <= 1000)))], &[(EBV, 9)]),
    );
    assert_eq!(nine[EBV]["call"], BELOW);
    let hcv = run("hcv", &case(vec![(HCV, flat(2))], &[(HCV, 50)]));
    assert!(hcv[HCV]["note"].contains("RNA virus with no DNA stage"));
    let ebv = run("ebv15", &case(vec![(EBV, flat(15))], &[(EBV, 500)]));
    assert_eq!(
        ebv[EBV]["call"], FOUND,
        "the HHV-6 rule never applies to EBV"
    );
}

#[test]
fn phix_is_the_lab_control() {
    let px = run("phix", &case(vec![(PHIX, flat(30))], &[(PHIX, 200)]));
    assert_eq!(px[PHIX]["call"], CONTROL);
    let one = run(
        "phix1",
        &case(vec![(PHIX, Box::new(|p| u32::from(p < 5)))], &[(PHIX, 1)]),
    );
    assert_eq!(one[PHIX]["call"], CONTROL);
    assert_eq!(run("phix0", &Case::default())[PHIX]["call"], NONE);
}

#[test]
fn both_species_over_the_threshold_are_both_called() {
    let both = run(
        "both",
        &case(vec![(A, flat(10)), (B, flat(15))], &[(A, 300), (B, 400)]),
    );
    assert_eq!(both[B]["call"], ICI);
    assert!(both[A]["call"] == ICI && both[A]["note"].contains("both pass"));
}

#[test]
fn the_human_depth() {
    let mut c = case(vec![(B, flat(15))], &[(B, 400)]);
    c.human = Some(vec![("chrX", 1000, 1000, 30.0)]);
    let nohum = run("nohuman", &c);
    assert_eq!(nohum[B]["human_depth"], ".");
    assert!(nohum[B]["call"] == FOUND && nohum[B]["note"].contains("integration is not assessed"));
    let c = Case {
        human: Some(vec![
            ("chr1", 1000, 500, 10.0),
            ("22", 1000, 1000, 40.0),
            ("chrX", 1000, 1000, 99.0),
            ("chrUn_KI270302v1", 1000, 1000, 99.0),
            ("chrEBV", 1000, 1000, 99.0),
        ]),
        ..Case::default()
    };
    assert_eq!(
        run("human", &c)[B]["human_depth"],
        format!("{:.2}", 50_000.0 / 1500.0),
        "over the autosomal bases with reads"
    );
}

#[test]
fn damaged_input() {
    let drop_hcv = |d: &Dir| {
        let t = fs::read_to_string(d.p("vdepth.tsv")).unwrap();
        let kept: String = t
            .lines()
            .filter(|l| !l.starts_with(&format!("{HCV}\t")))
            .map(|l| format!("{l}\n"))
            .collect();
        fs::write(d.p("vdepth.tsv"), kept).unwrap();
    };
    let gone: HashMap<_, _> = run_raw("gone", &case(vec![(B, flat(15))], &[(B, 400)]), drop_hcv)
        .unwrap()
        .into_iter()
        .collect();
    assert_eq!(
        [
            gone[HCV]["call"].as_str(),
            &gone[HCV]["breadth"],
            &gone[B]["call"]
        ],
        [NONE, "0.0000", ICI]
    );
    let cut = |d: &Dir| {
        let t = fs::read_to_string(d.p("vdepth.tsv")).unwrap();
        let l: Vec<&str> = t.lines().collect();
        let kept: String = l[..l.len() - 10].iter().map(|l| format!("{l}\n")).collect();
        fs::write(d.p("vdepth.tsv"), kept).unwrap();
    };
    let err = run_raw("cut", &Case::default(), cut).unwrap_err();
    assert!(err.contains("samtools depth -aa"), "{err}");
    let bad = |d: &Dir| fs::write(d.p("hcov.tsv"), "chr1\t1\t2\n").unwrap();
    let err = run_raw("badcov", &Case::default(), bad).unwrap_err();
    assert!(err.contains("not samtools coverage output"), "{err}");
}
