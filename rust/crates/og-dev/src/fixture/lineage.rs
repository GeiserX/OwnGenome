//! `og dev fixture <kind>`: the inputs the ancestry, family, sample QC, PRS,
//! Y, CHIP, MoChA, viral and live e2e cases make for themselves, and the few
//! table checks they make on what a step wrote.
//!
//! Each kind replaces a Python heredoc of a case under `tests/e2e/`. It takes
//! the same arguments in the same order and writes the same bytes, so a case
//! keeps every count and threshold it asserts. The random ones port CPython's
//! `random.Random` on the same seed ([`PyRandom`]).
//!
//! Exit 0 when written, 1 when an input cannot be read or has the wrong
//! shape, 2 on a usage error.

use super::Mt19937;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Read, Seek, SeekFrom, Write};

mod ancient;
mod chip;
mod family;
mod mocha;
mod prs;
mod qc;
mod truth;
mod viral;
mod ystr;

/// One kind: its name, its usage line and the function that runs it.
type Run = fn(&[String], &mut dyn Write, &mut dyn Write) -> Result<(), Fail>;

const TABLE: &[(&str, &str, Run)] = &[
    (
        "prs-sites",
        "BED TRUTH.vcf.gz FASTA HOMREF [SNVS [--no-ambiguous]]",
        prs::sites,
    ),
    (
        "x-ref-sites",
        "BLOCKS.tsv VARIANTS.txt FASTA CHROM",
        prs::x_ref_sites,
    ),
    (
        "csv-count",
        "FILE.csv[.gz] KEY=V1,V2|KEY!=V1,V2...",
        prs::csv_count,
    ),
    (
        "truth-genotypes",
        "SITES.tsv FAI [--chrom C] NAME=CALLS.vcf,CONFIDENT.bed...",
        truth::genotypes,
    ),
    (
        "ancient-panel",
        "REGIONS.bed TRUTH.vcf.gz CONFIDENT.bed FASTA OUT_DIR CHAIN",
        ancient::panel,
    ),
    ("plant-snv", "EVERY FAI READS.sam POS", chip::plant_snv),
    ("base-count", "POS BASE [READS.sam]", chip::base_count),
    ("sites-fasta", "SITES.tsv.gz OUT.fasta", chip::sites_fasta),
    (
        "mt-roundtrip",
        "LIFTED.tsv CALLS.tsv SITES.tsv.gz CHIP.txt",
        chip::mt_roundtrip,
    ),
    ("sites-vcf", "FAI SNVS.tsv OUT.vcf", family::sites_vcf),
    (
        "segments-check",
        "SEGMENTS.tsv DNAPAINTER.csv",
        family::segments_check,
    ),
    (
        "contig-cut",
        "FAI SITES.vcf.gz OUT_SITES.vcf PANEL_PREFIX OUT_PREFIX",
        qc::contig_cut,
    ),
    (
        "add-sites",
        "SITES.vcf EXTRA.tsv FAI OUT.vcf",
        qc::add_sites,
    ),
    (
        "sex-loss-woman",
        "WINDOWS.bed.gz OUT_DIR NAME LOSS",
        qc::sex_loss_woman,
    ),
    ("mocha-vcf", "FASTA OUT_SHIFT.vcf OUT_FLAT.vcf", mocha::vcf),
    ("viral-spike", "VIRAL.fa HUMAN_DEPTH FRACTION", viral::spike),
    (
        "y-seed",
        "CALLS.vcf.gz MASK.fasta.gz PRIVATE.tsv OUT.vcf TABLE...",
        ystr::y_seed,
    ),
    ("ystr-periodic", "REGIONS.bed REGIONS.fa", ystr::periodic),
    (
        "ystr-readcount",
        "PANEL.tsv CALLS.tsv READS.sam",
        ystr::readcount,
    ),
    ("ystr-cut", "START END READS.sam", ystr::cut),
    ("free-port", "", free_port),
];

/// The usage lines of every kind, one per line.
pub fn usage() -> String {
    TABLE
        .iter()
        .map(|(name, args, _)| format!("       og dev fixture {name} {args}\n"))
        .collect()
}

/// Runs `og dev fixture KIND ARGS` when KIND is one of these kinds; `None`
/// for a kind this module does not have.
pub fn cli(kind: &str, args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> Option<u8> {
    let (name, line, run) = TABLE.iter().find(|(name, _, _)| *name == kind)?;
    if matches!(args, [a] if a == "-h" || a == "--help") {
        return Some(match writeln!(out, "usage: og dev fixture {name} {line}") {
            Ok(()) => 0,
            Err(e) => {
                let _ = writeln!(err, "og dev fixture {name}: {e}");
                1
            }
        });
    }
    Some(match run(args, out, err) {
        Ok(()) => 0,
        Err(Fail::Usage(msg)) => {
            let _ = writeln!(
                err,
                "usage: og dev fixture {name} {line}\nog dev fixture {name}: error: {msg}"
            );
            2
        }
        Err(Fail::Run(msg)) => {
            let _ = writeln!(err, "og dev fixture {name}: {msg}");
            1
        }
    })
}

/// Why a kind stopped.
#[derive(Debug)]
pub enum Fail {
    /// The arguments are wrong (exit 2).
    Usage(String),
    /// An input is missing or has the wrong shape (exit 1).
    Run(String),
}

impl From<io::Error> for Fail {
    fn from(e: io::Error) -> Self {
        Fail::Run(e.to_string())
    }
}

/// A [`Fail::Run`] from any message.
fn bad(msg: impl Into<String>) -> Fail {
    Fail::Run(msg.into())
}

/// Exactly `n` positional arguments, or a usage error.
fn want(args: &[String], n: usize) -> Result<(), Fail> {
    if args.len() == n {
        Ok(())
    } else {
        Err(Fail::Usage(format!(
            "takes {n} arguments, got {}",
            args.len()
        )))
    }
}

/// A whole number argument.
fn int<T: std::str::FromStr>(s: &str, what: &str) -> Result<T, Fail> {
    s.trim()
        .parse()
        .map_err(|_| Fail::Usage(format!("{what} '{s}' is not a number")))
}

/// `int(s)` of a table cell, or an error naming the cell.
fn cell_int(s: &str) -> Result<i64, Fail> {
    s.trim()
        .parse()
        .map_err(|_| bad(format!("'{s}' is not a whole number")))
}

/// Opens a file, plain or gzip, with the path in the error.
fn open(path: &str) -> Result<Box<dyn BufRead>, Fail> {
    og_io::open(path).map_err(|e| bad(format!("{path}: {e}")))
}

/// Calls `f` on each line of `path` (plain or gzip) without its `\n`, as
/// `for l in open(path): l.rstrip("\n")` reads it.
fn each_line(path: &str, mut f: impl FnMut(&str) -> Result<(), Fail>) -> Result<(), Fail> {
    let mut input = open(path)?;
    let mut buf = Vec::new();
    loop {
        buf.clear();
        if input
            .read_until(b'\n', &mut buf)
            .map_err(|e| bad(format!("{path}: {e}")))?
            == 0
        {
            return Ok(());
        }
        if buf.last() == Some(&b'\n') {
            buf.pop();
        }
        let line = std::str::from_utf8(&buf).map_err(|_| bad(format!("{path}: not UTF-8")))?;
        f(line)?;
    }
}

/// Every line of `path` without its `\n`.
fn lines(path: &str) -> Result<Vec<String>, Fail> {
    let mut v = Vec::new();
    each_line(path, |l| {
        v.push(l.to_string());
        Ok(())
    })?;
    Ok(v)
}

/// A file to write, buffered, with the path in the error.
fn create(path: &str) -> Result<BufWriter<File>, Fail> {
    File::create(path)
        .map(BufWriter::new)
        .map_err(|e| bad(format!("{path}: {e}")))
}

/// CPython's `random.Random(seed)`, as far as the cases use it.
pub struct PyRandom(Mt19937);

impl PyRandom {
    /// `random.Random(seed)` and `random.seed(seed)`.
    pub fn new(seed: u32) -> Self {
        PyRandom(Mt19937::new(seed))
    }

    /// `random()`: 53 random bits in [0, 1).
    pub fn random(&mut self) -> f64 {
        let a = f64::from(self.0.next_u32() >> 5);
        let b = f64::from(self.0.next_u32() >> 6);
        (a * 67_108_864.0 + b) * (1.0 / 9_007_199_254_740_992.0)
    }

    /// `uniform(a, b)`.
    pub fn uniform(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.random()
    }

    /// `randrange(0, n)` for `0 < n < 2^32`.
    pub fn below(&mut self, n: u64) -> u64 {
        let n = u32::try_from(n).expect("a range below 2^32");
        u64::from(self.0.below(n))
    }

    /// `randint(a, b)`.
    pub fn randint(&mut self, a: i64, b: i64) -> i64 {
        let width = u64::try_from(b - a + 1).expect("b >= a");
        a + i64::try_from(self.below(width)).expect("small")
    }
}

/// Random access to the bases of an indexed FASTA (`.fai` beside it), read
/// as Python's `seek` and `read(1)` on the file read them, through a cache
/// of one 64 KiB window.
pub struct Fasta {
    file: File,
    fai: HashMap<String, [u64; 4]>,
    window: Vec<u8>,
    start: u64,
}

impl Fasta {
    /// Opens `path` and reads `path.fai`.
    pub fn open(path: &str) -> Result<Self, Fail> {
        let fai = read_fai(&format!("{path}.fai"))?;
        let file = File::open(path).map_err(|e| bad(format!("{path}: {e}")))?;
        Ok(Fasta {
            file,
            fai,
            window: Vec::new(),
            start: 0,
        })
    }

    /// The length of a contig, if the index has it.
    pub fn len(&self, chrom: &str) -> Option<u64> {
        self.fai.get(chrom).map(|f| f[0])
    }

    /// The upper-case base at 1-based `pos` of `chrom`; `None` past the end
    /// of the file or for a contig the index lacks.
    pub fn base(&mut self, chrom: &str, pos: u64) -> Result<Option<u8>, Fail> {
        let Some(&[_, off, lb, lw]) = self.fai.get(chrom) else {
            return Ok(None);
        };
        if pos == 0 || lb == 0 {
            return Ok(None);
        }
        let at = off + (pos - 1) / lb * lw + (pos - 1) % lb;
        if at < self.start || at >= self.start + self.window.len() as u64 {
            self.start = at;
            self.window.clear();
            self.file.seek(SeekFrom::Start(at))?;
            (&mut self.file)
                .take(1 << 16)
                .read_to_end(&mut self.window)?;
        }
        Ok(self
            .window
            .get(usize::try_from(at - self.start).unwrap_or(usize::MAX))
            .map(u8::to_ascii_uppercase))
    }
}

/// A `.fai`: contig to `[length, offset, line bases, line width]`.
fn read_fai(path: &str) -> Result<HashMap<String, [u64; 4]>, Fail> {
    let mut fai = HashMap::new();
    each_line(path, |l| {
        let c: Vec<&str> = l.split_whitespace().collect();
        if c.len() >= 5 {
            let mut v = [0u64; 4];
            for (slot, x) in v.iter_mut().zip(&c[1..5]) {
                *slot = x
                    .parse()
                    .map_err(|_| bad(format!("{path}: '{x}' is not a number")))?;
            }
            fai.insert(c[0].to_string(), v);
        }
        Ok(())
    })?;
    Ok(fai)
}

/// The contigs of a `.fai` in its order, with their length as written.
fn fai_contigs(path: &str) -> Result<Vec<(String, String)>, Fail> {
    let mut v = Vec::new();
    each_line(path, |l| {
        let mut c = l.split('\t');
        let name = c.next().unwrap_or("").to_string();
        let len = c.next().unwrap_or("").to_string();
        v.push((name, len));
        Ok(())
    })?;
    Ok(v)
}

/// The intervals of a BED on one contig, for `bisect` lookups: Python's
/// `confident(c, p)` of the cases (the last interval starting before `p`
/// must reach `p`).
#[derive(Default)]
pub struct Confident {
    starts: Vec<i64>,
    ends: Vec<i64>,
}

impl Confident {
    /// From `(start, end)` pairs, sorted here.
    pub fn new(mut iv: Vec<(i64, i64)>) -> Self {
        iv.sort_unstable();
        Confident {
            starts: iv.iter().map(|x| x.0).collect(),
            ends: iv.iter().map(|x| x.1).collect(),
        }
    }

    /// Whether 1-based `p` lies in the last interval that starts before it.
    pub fn has(&self, p: i64) -> bool {
        let i = self.starts.partition_point(|&s| s < p);
        i > 0 && self.ends[i - 1] >= p
    }
}

/// A BED per contig, from its first three columns (lines with fewer are
/// skipped).
fn read_bed(path: &str) -> Result<HashMap<String, Confident>, Fail> {
    let mut iv: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
    each_line(path, |l| {
        let c: Vec<&str> = l.split_whitespace().collect();
        if c.len() >= 3 {
            iv.entry(c[0].to_string())
                .or_default()
                .push((cell_int(c[1])?, cell_int(c[2])?));
        }
        Ok(())
    })?;
    Ok(iv
        .into_iter()
        .map(|(k, v)| (k, Confident::new(v)))
        .collect())
}

/// One CIGAR operation: its length and its letter.
type Op = (i64, u8);

/// A CIGAR string as `(length, op)` pairs, the way `re.findall(r"(\d+)([MIDNSHP=X])")`
/// reads it: anything else is skipped.
fn cigar(s: &str) -> Vec<Op> {
    let mut ops = Vec::new();
    let mut n: Option<i64> = None;
    for b in s.bytes() {
        if b.is_ascii_digit() {
            n = Some(n.unwrap_or(0) * 10 + i64::from(b - b'0'));
        } else if let Some(len) = n.take()
            && b"MIDNSHP=X".contains(&b)
        {
            ops.push((len, b));
        }
    }
    ops
}

/// `og dev fixture free-port`: a TCP port on 127.0.0.1 that was free a
/// moment ago (the OS picks it), for a case that starts a server.
fn free_port(args: &[String], out: &mut dyn Write, _err: &mut dyn Write) -> Result<(), Fail> {
    want(args, 0)?;
    let port = std::net::TcpListener::bind(("127.0.0.1", 0))?
        .local_addr()?
        .port();
    writeln!(out, "{port}")?;
    Ok(())
}

#[cfg(test)]
pub(crate) mod testutil {
    use std::path::PathBuf;

    /// A fresh scratch folder for one test.
    pub fn dir(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("og-lineage-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).expect("scratch");
        d
    }

    /// Runs a kind; returns its exit code, stdout and stderr.
    pub fn run(kind: &str, args: &[&str]) -> (u8, String, String) {
        let args: Vec<String> = args.iter().map(|s| s.to_string()).collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = super::cli(kind, &args, &mut out, &mut err).expect("a kind of this module");
        (
            code,
            String::from_utf8(out).expect("utf8"),
            String::from_utf8(err).expect("utf8"),
        )
    }

    /// 64-bit FNV-1a, to pin bytes in a test.
    pub fn fnv(bytes: &[u8]) -> u64 {
        bytes.iter().fold(0xcbf2_9ce4_8422_2325, |h, &b| {
            (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3)
        })
    }

    /// A FASTA with its `.fai`: each contig a seeded random sequence of
    /// `ACGT` (and an `N` run where asked), 60 bases a line.
    pub fn fasta(path: &std::path::Path, contigs: &[(&str, usize)], seed: u32) {
        let mut rng = super::PyRandom::new(seed);
        let (mut fa, mut fai, mut off) = (String::new(), String::new(), 0usize);
        for (name, n) in contigs {
            let head = format!(">{name}\n");
            off += head.len();
            fa.push_str(&head);
            fai.push_str(&format!("{name}\t{n}\t{off}\t60\t61\n"));
            let seq: Vec<u8> = (0..*n)
                .map(|_| b"ACGTacgt"[rng.below(8) as usize])
                .collect();
            for chunk in seq.chunks(60) {
                fa.push_str(std::str::from_utf8(chunk).expect("ascii"));
                fa.push('\n');
            }
            off += n + n.div_ceil(60);
        }
        std::fs::write(path, fa).expect("fasta");
        std::fs::write(format!("{}.fai", path.display()), fai).expect("fai");
    }
}

#[cfg(test)]
mod tests {
    use super::testutil::run;
    use super::*;

    #[test]
    fn random_is_cpythons() {
        // python3 -c 'import random; r = random.Random(44); print([r.random() for _ in range(2)], r.randint(800, 1200), r.uniform(0.05, 0.95))'
        let mut r = PyRandom::new(44);
        assert_eq!(r.random(), 0.408_535_879_254_493_75);
        assert_eq!(r.random(), 0.541_972_009_104_596_2);
        assert_eq!(r.randint(800, 1200), 859);
        assert_eq!(r.uniform(0.05, 0.95), 0.208_999_079_775_899_43);
    }

    #[test]
    fn cigar_reads_as_the_regex_does() {
        assert_eq!(
            cigar("5S10M2I3D1N4="),
            vec![
                (5, b'S'),
                (10, b'M'),
                (2, b'I'),
                (3, b'D'),
                (1, b'N'),
                (4, b'=')
            ]
        );
        assert_eq!(cigar("*"), vec![]);
    }

    #[test]
    fn confident_is_bisects() {
        let c = Confident::new(vec![(10, 20), (30, 40)]);
        assert!(!c.has(10) && c.has(11) && c.has(20) && !c.has(21) && c.has(31) && !c.has(41));
    }

    #[test]
    fn unknown_kind_is_none_and_help_prints_usage() {
        let (mut out, mut err) = (Vec::new(), Vec::new());
        assert!(cli("no-such", &[], &mut out, &mut err).is_none());
        let (code, out, _) = run("ystr-cut", &["--help"]);
        assert_eq!(
            (code, out.as_str()),
            (0, "usage: og dev fixture ystr-cut START END READS.sam\n")
        );
        let (code, _, err) = run("ystr-cut", &["1"]);
        assert_eq!(code, 2, "{err}");
        assert!(usage().lines().count() == TABLE.len());
    }

    #[test]
    fn free_port_prints_a_port() {
        let (code, out, _) = run("free-port", &[]);
        assert_eq!(code, 0);
        assert!(out.trim().parse::<u16>().expect("a port") > 0);
    }
}
