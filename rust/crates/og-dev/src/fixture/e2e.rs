//! The inputs the genotype, repeat and copy-number e2e cases make for
//! themselves: made-up references, reads and records with fixed seeds.
//! Each kind writes the bytes the inline Python of its case wrote, so a
//! case keeps what it asserts.
//!
//! - `append-n SRC OUT CONTIG LENGTH`: SRC, then CONTIG as LENGTH `N`s, 60
//!   a line (a stand-in contig a reference lacks).
//! - `rfc1 KIND OUT.sam [--samd12 UNIT]`: made-up reads at the GRCh38 RFC1
//!   repeat ([`rfc1`]); `--samd12` adds reads of UNIT and reference reads at
//!   SAMD12 on chr8.
//! - `str-extra REF.fa HET_NAME HET.sam HOM_NAME HOM.sam`: a reference that
//!   is N except a window around the CFTR TG/poly-T and UGT1A1 TA tracts,
//!   and read pairs of two samples (`e2e/str_extra.rs`).
//! - `tas2r38 REF_DIR GENOME_DIR`: a chr7 that is N except a window around
//!   TAS2R38, with gVCFs, VCFs and reads of four samples (`e2e/tas2r38.rs`).
//! - `mito-sam --fasta FA --header SAM --sample S --out OUT [--depth X]
//!   [--deletion FIRST-LAST] [--level F]`: single-end chrM reads, some across
//!   a deletion's junction (`e2e/mito_sam.rs`).
//! - `hba-tiles FASTA HEADER SAMPLE OUT`: reads tiled every 5 bases over a
//!   chr16 window, MAPQ 0 where a read fits twice; a sample whose name ends
//!   in `3` loses every other read over HBA2 (`e2e/hba_tiles.rs`).
//! - `cut-haplotype START END FQ NAMES`: SAM on stdin; the haplotype-1 reads
//!   lose bases START..END-1 and go to FQ as FASTQ (`e2e/cut_haplotype.rs`).
//! - `contains FILE PART`: exit 0 when FILE holds PART's bytes, 1 when not.

mod append_n;
mod contains;
mod cut_haplotype;
mod hba_tiles;
mod mito_sam;
pub mod rfc1;
mod str_extra;
mod tas2r38;

use super::Mt19937;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

/// The usage lines of these kinds, one a line.
pub const USAGE: &str = "       og dev fixture append-n SRC OUT CONTIG LENGTH\n       og dev fixture rfc1 KIND OUT.sam [--samd12 UNIT]\n       og dev fixture str-extra REF.fa HET_NAME HET.sam HOM_NAME HOM.sam\n       og dev fixture tas2r38 REF_DIR GENOME_DIR\n       og dev fixture mito-sam --fasta FA --header SAM --sample S --out OUT [--depth X] [--deletion FIRST-LAST] [--level F]\n       og dev fixture hba-tiles FASTA HEADER SAMPLE OUT\n       og dev fixture cut-haplotype START END FQ NAMES < SAM\n       og dev fixture contains FILE PART\n";

/// What a kind ended in: written, a usage error, or a failure.
pub enum Fail {
    Usage(String),
    Run(String),
}

impl From<io::Error> for Fail {
    fn from(e: io::Error) -> Self {
        Fail::Run(e.to_string())
    }
}

/// Runs kind `kind` with `args`, or `None` when no kind here has that name.
pub fn cli(kind: &str, args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> Option<u8> {
    let result = match kind {
        "append-n" => append_n::run(args),
        "rfc1" => rfc1::run(args),
        "str-extra" => str_extra::run(args),
        "tas2r38" => tas2r38::run(args),
        "mito-sam" => mito_sam::run(args, out),
        "hba-tiles" => hba_tiles::run(args),
        "cut-haplotype" => cut_haplotype::run(args, &mut io::stdin().lock(), out, err),
        "contains" => return Some(contains::run(args, err)),
        _ => return None,
    };
    Some(
        match result.and_then(|()| out.flush().map_err(Fail::from)) {
            Ok(()) => 0,
            Err(Fail::Usage(msg)) => {
                let _ = writeln!(err, "usage:\n{USAGE}og dev fixture: error: {msg}");
                2
            }
            Err(Fail::Run(msg)) => {
                let _ = writeln!(err, "og dev fixture {kind}: {msg}");
                1
            }
        },
    )
}

/// `n` arguments exactly, or a usage error naming them.
fn exactly<'a>(args: &'a [String], n: usize, names: &str) -> Result<&'a [String], Fail> {
    if args.len() == n {
        Ok(args)
    } else {
        Err(Fail::Usage(format!("takes {names}")))
    }
}

/// A whole number argument.
fn number(s: &str, what: &str) -> Result<usize, Fail> {
    s.parse()
        .map_err(|_| Fail::Usage(format!("{what} '{s}' is not a whole number")))
}

/// A file to write, buffered; the error names the path.
fn create(path: &Path) -> Result<BufWriter<File>, Fail> {
    File::create(path)
        .map(BufWriter::new)
        .map_err(|e| Fail::Run(format!("{}: {e}", path.display())))
}

/// A file to read, whole; the error names the path.
fn read(path: &Path) -> Result<Vec<u8>, Fail> {
    std::fs::read(path).map_err(|e| Fail::Run(format!("{}: {e}", path.display())))
}

/// CPython's `random.Random(seed)`, as far as the cases use it.
pub struct PyRandom(Mt19937);

impl PyRandom {
    /// `random.Random(seed)` (and `random.seed(seed)`).
    pub fn new(seed: u32) -> Self {
        PyRandom(Mt19937::new(seed))
    }

    /// `randrange(0, n)`, for `0 < n < 2^32`.
    pub fn below(&mut self, n: usize) -> usize {
        let n = u32::try_from(n).expect("a range below 2^32");
        self.0.below(n) as usize
    }

    /// `choice(seq)`.
    pub fn choice<T: Copy>(&mut self, seq: &[T]) -> T {
        seq[self.below(seq.len())]
    }

    /// `random()`: 53 random bits in [0, 1).
    pub fn random(&mut self) -> f64 {
        let a = f64::from(self.0.next_u32() >> 5);
        let b = f64::from(self.0.next_u32() >> 6);
        (a * 67_108_864.0 + b) * (1.0 / 9_007_199_254_740_992.0)
    }

    /// `"".join(choice("ACGT") for _ in range(n))`.
    pub fn bases(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| self.choice(b"ACGT")).collect()
    }
}

/// The reverse complement over `ACGTN`; other bytes stay as they are.
pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            x => x,
        })
        .collect()
}

/// Writes `n` `N`s as FASTA lines of 60, continuing a line that already
/// holds `col` bases.
fn n_lines(w: &mut dyn Write, mut n: usize, col: usize) -> io::Result<usize> {
    const NS: [u8; 60] = [b'N'; 60];
    let mut col = col;
    while n > 0 {
        let k = n.min(60 - col);
        w.write_all(&NS[..k])?;
        col += k;
        n -= k;
        if col == 60 {
            w.write_all(b"\n")?;
            col = 0;
        }
    }
    Ok(col)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// FNV-1a, 64 bits: a pin of the bytes a kind writes.
    pub fn fnv(bytes: &[u8]) -> u64 {
        bytes.iter().fold(0xcbf2_9ce4_8422_2325, |h, &b| {
            (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3)
        })
    }

    /// A writer that keeps only an FNV-1a hash, the line count and the
    /// longest line of what it is given.
    pub struct Hashed {
        pub hash: u64,
        pub lines: usize,
        pub longest: usize,
        col: usize,
    }

    impl Default for Hashed {
        fn default() -> Self {
            Hashed {
                hash: fnv(b""),
                lines: 0,
                longest: 0,
                col: 0,
            }
        }
    }

    impl Write for Hashed {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            for &b in buf {
                self.hash = (self.hash ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3);
                if b == b'\n' {
                    self.lines += 1;
                    self.longest = self.longest.max(self.col);
                    self.col = 0;
                } else {
                    self.col += 1;
                }
            }
            Ok(buf.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    pub fn scratch(name: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("og-e2e-fixture-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    pub fn run(kind: &str, args: &[&str]) -> (Option<u8>, String, String) {
        let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = cli(kind, &args, &mut out, &mut err);
        (
            code,
            String::from_utf8(out).unwrap(),
            String::from_utf8(err).unwrap(),
        )
    }

    #[test]
    fn random_is_cpythons() {
        // python3 -c 'import random; r = random.Random(1);
        //   print(r.random(), r.randrange(0, 16569), r.choice((0, 16)))'
        let mut r = PyRandom::new(1);
        assert_eq!(r.random(), 0.134_364_244_112_401_22);
        assert_eq!(r.below(16569), 2067);
        assert_eq!(r.choice(&[0, 16]), 16);
        // r = random.Random(7); "".join(r.choice("ACGT") for _ in range(12))
        assert_eq!(PyRandom::new(7).bases(12), b"GCTAAAGACAAT");
    }

    #[test]
    fn n_lines_fill_60_a_line_and_carry_a_partial_line() {
        let mut v = Vec::new();
        assert_eq!(n_lines(&mut v, 130, 0).unwrap(), 10);
        assert_eq!(v.len(), 132);
        assert_eq!(n_lines(&mut v, 50, 10).unwrap(), 0);
        assert_eq!(v.iter().filter(|&&b| b == b'\n').count(), 3);
    }

    #[test]
    fn unknown_kinds_are_left_to_the_caller_and_bad_arguments_exit_2() {
        assert_eq!(run("no-such-kind", &[]).0, None);
        let (code, _, err) = run("append-n", &["a", "b"]);
        assert_eq!(code, Some(2));
        assert!(
            err.contains("og dev fixture: error: takes SRC OUT CONTIG LENGTH"),
            "{err}"
        );
        let (code, _, err) = run(
            "append-n",
            &["/no/such/file", "/no/such/out", "chr11", "60"],
        );
        assert_eq!(code, Some(1), "{err}");
    }
}
