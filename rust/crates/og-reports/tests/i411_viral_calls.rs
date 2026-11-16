//! `og viral report` on synthetic inputs (#411 trait-51 to trait-56): which
//! reads count as viral evidence, and what a call says.
//!
//! Each case writes a viral FASTA of small genomes with real accessions, the
//! `samtools coverage`, `depth -aa` and `view -F 4` files of the viral BAM,
//! the human coverage at 30x and `samtools idxstats` of the sample's BAM,
//! then reads the table `report` writes.

use og_reports::viral::{Files, report};
use std::collections::HashMap;
use std::fmt::Write as _;
use std::path::PathBuf;

const A: &str = "NC_001664.4";
const B: &str = "NC_000898.1";
const EBV: &str = "NC_007605.1";
const HBV: &str = "NC_003977.2";
const VZV: &str = "NC_001348.1";
const PHIX: &str = "NC_001422.1";
const GENOMES: [(&str, &str, usize); 6] = [
    (A, "Human betaherpesvirus 6A, complete genome", 2000),
    (B, "Human herpesvirus 6B, complete genome", 2000),
    (EBV, "Human gammaherpesvirus 4, complete genome", 1500),
    (HBV, "Hepatitis B virus, complete genome", 3200),
    (VZV, "Human alphaherpesvirus 3, complete genome", 2000),
    (PHIX, "Escherichia phage phiX174, complete genome", 500),
];
const COV_HEAD: &str =
    "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n";

/// A fixed pseudo-random sequence: not low complexity.
fn sequence(len: usize, seed: u32) -> String {
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            char::from(b"ACGT"[(x >> 16) as usize % 4])
        })
        .collect()
}

/// One SAM record on `acc` at 1-based `pos`, `len` bases, with `seq`.
fn sam(acc: &str, flag: u16, pos: usize, mapq: u8, seq: &str) -> String {
    format!(
        "r\t{flag}\t{acc}\t{pos}\t{mapq}\t{}M\t*\t0\t0\t{seq}\t*\n",
        seq.len()
    )
}

/// `len` bases of genome `acc` from 1-based `pos`.
fn slice(acc: &str, pos: usize, len: usize) -> String {
    let (i, &(_, _, length)) = GENOMES
        .iter()
        .enumerate()
        .find(|(_, g)| g.0 == acc)
        .unwrap();
    sequence(length, i as u32 + 1)[pos - 1..pos - 1 + len].to_owned()
}

/// A genome's slice, as a read placed at 1-based `pos`, MAPQ 60.
fn read(acc: &str, flag: u16, pos: usize, len: usize) -> String {
    sam(acc, flag, pos, 60, &slice(acc, pos, len))
}

/// A per-base depth function.
type Depth = Box<dyn Fn(usize) -> u32>;

fn depth(f: impl Fn(usize) -> u32 + 'static) -> Depth {
    Box::new(f)
}

/// The inputs of one case.
#[derive(Default)]
struct Case {
    /// Per-base depth by genome; zero elsewhere.
    depth: HashMap<&'static str, Depth>,
    /// `numreads` of `samtools coverage` by genome.
    aligned: HashMap<&'static str, usize>,
    /// The `samtools view -F 4` lines.
    hits: String,
    /// The unplaced unmapped reads in the sample's BAM.
    unplaced: Option<u64>,
}

struct Out {
    rows: HashMap<String, HashMap<String, String>>,
    stdout: String,
}

fn run(tag: &str, c: Case) -> Out {
    let d: PathBuf =
        std::env::temp_dir().join(format!("og-i411-viral-{tag}-{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    let p = |n: &str| d.join(n).to_string_lossy().into_owned();
    let (mut fa, mut cov, mut dep) = (String::new(), COV_HEAD.to_owned(), String::new());
    for (i, (acc, title, len)) in GENOMES.iter().enumerate() {
        let _ = write!(fa, ">{acc} {title}\n{}\n", sequence(*len, i as u32 + 1));
        let depths: Vec<u32> = (1..=*len)
            .map(|q| c.depth.get(acc).map_or(0, |f| f(q)))
            .collect();
        let covered = depths.iter().filter(|x| **x > 0).count();
        let mean = depths.iter().map(|x| f64::from(*x)).sum::<f64>() / *len as f64;
        let _ = writeln!(
            cov,
            "{acc}\t1\t{len}\t{}\t{covered}\t0\t{mean}\t30\t60",
            c.aligned.get(acc).copied().unwrap_or(0)
        );
        for (q, x) in depths.iter().enumerate() {
            let _ = writeln!(dep, "{acc}\t{}\t{x}", q + 1);
        }
    }
    let files = [
        ("v.fa", fa),
        ("vcov.tsv", cov),
        ("vdepth.tsv", dep),
        (
            "hcov.tsv",
            format!("{COV_HEAD}chr1\t1\t100000\t1000\t100000\t100\t30\t30\t60\n"),
        ),
        ("hits.sam", format!("@HD\tVN:1.6\n{}", c.hits)),
        (
            "idx.tsv",
            format!(
                "chr1\t100000\t1000\t2\n*\t0\t0\t{}\n",
                c.unplaced.unwrap_or(5000)
            ),
        ),
    ];
    for (n, text) in &files {
        std::fs::write(d.join(n), text).unwrap();
    }
    let stdout = report(
        "S",
        &Files {
            viral_fasta: &p("v.fa"),
            viral_coverage: &p("vcov.tsv"),
            viral_depth: &p("vdepth.tsv"),
            human_coverage: &p("hcov.tsv"),
            viral_hits: &p("hits.sam"),
            idxstats: &p("idx.tsv"),
            out: &p("out.tsv"),
        },
    )
    .unwrap();
    let table = std::fs::read_to_string(d.join("out.tsv")).unwrap();
    std::fs::remove_dir_all(&d).unwrap();
    let mut lines = table.lines().map(|l| l.split('\t').collect::<Vec<_>>());
    let head: Vec<String> = lines
        .next()
        .unwrap()
        .iter()
        .map(|s| (*s).to_owned())
        .collect();
    let rows = lines
        .map(|r| {
            let row: HashMap<String, String> = head
                .iter()
                .cloned()
                .zip(r.iter().map(|s| (*s).to_owned()))
                .collect();
            (row["accession"].clone(), row)
        })
        .collect();
    Out { rows, stdout }
}

fn cell<'a>(o: &'a Out, acc: &str, k: &str) -> &'a str {
    &o.rows[acc][k]
}

/// Red first: 4 reads on EBV read "below threshold", with the count, and the
/// stdout line lists it instead of saying nothing was found.
#[test]
fn four_reads_on_ebv_read_below_threshold() {
    let hits: String = (0..4).map(|i| read(EBV, 0, 1 + i * 300, 150)).collect();
    let o = run(
        "ebv4",
        Case {
            depth: HashMap::from([(
                EBV,
                depth(|q: usize| u32::from((q - 1) % 300 < 150 && q <= 1050)),
            )]),
            aligned: HashMap::from([(EBV, 4)]),
            hits,
            ..Case::default()
        },
    );
    assert_eq!(cell(&o, EBV, "call"), "below threshold");
    assert_eq!(cell(&o, EBV, "reads"), "4");
    assert_eq!(cell(&o, EBV, "breadth"), "0.4000");
    assert_eq!(cell(&o, EBV, "copies_per_cell"), ".");
    assert!(cell(&o, EBV, "note").starts_with("4 counted reads, under the 10"));
    assert_eq!(cell(&o, HBV, "call"), "not found");
    assert!(
        o.stdout
            .contains("Human gammaherpesvirus 4 (below threshold, 4 reads)"),
        "{}",
        o.stdout
    );
}

/// Nothing counted anywhere: the stdout line says what a negative means.
#[test]
fn a_negative_says_no_reads_above_the_threshold() {
    let o = run("none", Case::default());
    assert_eq!(
        o.stdout,
        "S: 6 viral genomes; no reads above the threshold at this depth\n"
    );
    assert!(o.rows.values().all(|r| r["call"] == "not found"));
}

/// An HHV-6 genome at half the autosomal depth is integrated HHV-6, read
/// from the mode; a 400x pile-up on 100 bases does not change its one copy
/// per cell (red first: the mean gave 2.2833).
#[test]
fn hhv6_at_half_the_autosomal_depth_reads_integrated() {
    let hits: String = (0..13).map(|i| read(B, 0, 1 + i * 150, 150)).collect();
    let o = run(
        "ici",
        Case {
            depth: HashMap::from([(B, depth(|q: usize| if q <= 100 { 400 } else { 15 }))]),
            aligned: HashMap::from([(B, 400)]),
            hits,
            ..Case::default()
        },
    );
    assert_eq!(cell(&o, B, "call"), "inherited integrated HHV-6 likely");
    assert_eq!(cell(&o, B, "relative_depth"), "0.500");
    assert_eq!(cell(&o, B, "copies_per_cell"), "1.0000");
    assert_eq!(cell(&o, A, "call"), "not found");
}

/// A 400x pile-up on 100 bases and nothing else is no copy per cell.
#[test]
fn a_pile_up_has_no_copies_per_cell() {
    let hits: String = (0..12).map(|i| read(B, 0, 1 + i * 5, 40)).collect();
    let o = run(
        "pile",
        Case {
            depth: HashMap::from([(B, depth(|q: usize| if q <= 100 { 400 } else { 0 }))]),
            aligned: HashMap::from([(B, 300)]),
            hits,
            ..Case::default()
        },
    );
    assert_eq!(cell(&o, B, "call"), "a few reads");
    assert_eq!(cell(&o, B, "mean_depth"), "20.000");
    assert_eq!(cell(&o, B, "copies_per_cell"), ".");
}

/// Both species pass: both are called, and neither note claims the other
/// has the higher depth.
#[test]
fn both_hhv6_species_pass_and_both_are_called() {
    let o = run(
        "both",
        Case {
            depth: HashMap::from([(A, depth(|_| 15)), (B, depth(|_| 15))]),
            aligned: HashMap::from([(A, 400), (B, 400)]),
            ..Case::default()
        },
    );
    for acc in [A, B] {
        assert_eq!(cell(&o, acc, "call"), "inherited integrated HHV-6 likely");
        assert!(cell(&o, acc, "note").contains("both pass"));
        assert!(!cell(&o, acc, "note").contains("higher depth"));
    }
}

/// Red first: 3 HBV molecules read 4 times each (the same genome, strand and
/// 5' end) are 3 counted reads, below threshold, not viral DNA found; 12
/// distinct reads at the same places are found.
#[test]
fn duplicates_do_not_make_hbv_found() {
    let dups: String = (0..3)
        .flat_map(|i| (0..4).map(move |_| read(HBV, 0, 1 + i * 1000, 150)))
        .collect();
    let o = run(
        "dup",
        Case {
            aligned: HashMap::from([(HBV, 12)]),
            hits: dups,
            ..Case::default()
        },
    );
    assert_eq!(cell(&o, HBV, "reads"), "3");
    assert_eq!(cell(&o, HBV, "call"), "below threshold");
    let distinct: String = (0..3)
        .flat_map(|i| (0..4).map(move |j| read(HBV, 0, 1 + i * 1000 + j * 20, 150)))
        .collect();
    let o = run(
        "dup-control",
        Case {
            aligned: HashMap::from([(HBV, 12)]),
            hits: distinct,
            ..Case::default()
        },
    );
    assert_eq!(cell(&o, HBV, "reads"), "12");
    assert_eq!(cell(&o, HBV, "call"), "viral DNA found");
}

/// Red first: low-complexity reads at MAPQ 60 and reads at MAPQ 0 are not
/// counted, so 60 aligned records on VZV read "not found", with a note.
#[test]
fn low_complexity_and_mapq_0_reads_are_not_counted() {
    let mut hits = String::new();
    for i in 0..20 {
        hits += &sam(VZV, 0, 1 + i * 50, 60, &"A".repeat(150));
        hits += &sam(VZV, 16, 1 + i * 50, 60, &"CA".repeat(75));
        hits += &sam(VZV, 0, 1 + i * 90, 0, &slice(VZV, 1 + i * 90, 150));
    }
    let o = run(
        "lowcx",
        Case {
            depth: HashMap::from([(VZV, depth(|_| 4))]),
            aligned: HashMap::from([(VZV, 60)]),
            hits,
            ..Case::default()
        },
    );
    assert_eq!(cell(&o, VZV, "reads"), "0");
    assert_eq!(cell(&o, VZV, "breadth"), "0.0000");
    assert_eq!(cell(&o, VZV, "call"), "not found");
    assert!(cell(&o, VZV, "note").starts_with("60 aligned records, none counted"));
}

/// A supplementary record adds breadth, not a read.
#[test]
fn supplementary_records_add_breadth_not_reads() {
    let mut hits: String = (0..10).map(|i| read(PHIX, 0, 1 + i * 10, 10)).collect();
    hits += &read(EBV, 0, 1, 100);
    hits += &read(EBV, 2048, 1001, 100);
    let o = run(
        "supp",
        Case {
            aligned: HashMap::from([(EBV, 2), (PHIX, 10)]),
            hits,
            ..Case::default()
        },
    );
    assert_eq!(cell(&o, EBV, "reads"), "1");
    assert_eq!(cell(&o, EBV, "breadth"), "0.1333");
    assert_eq!(cell(&o, PHIX, "call"), "lab control (phiX)");
}

/// Red first (#411 trait-53): a BAM with no unplaced unmapped read was cut
/// or filtered, so every genome reads "not assessed", not "not found".
#[test]
fn no_unplaced_unmapped_read_is_not_assessed() {
    let o = run(
        "cut",
        Case {
            unplaced: Some(0),
            ..Case::default()
        },
    );
    assert!(o.rows.values().all(|r| r["call"] == "not assessed"));
    assert!(o.rows[B]["note"].contains("no unplaced unmapped read"));
    assert_eq!(
        o.stdout,
        "S: 6 viral genomes; not assessed: the BAM has no unplaced unmapped read\n"
    );
}
