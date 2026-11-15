//! `rfc1 KIND OUT.sam [--samd12 UNIT]`: made-up reads at the GRCh38 RFC1
//! repeat (chr4:39348425-39348479, the reference (AAAAG)11) with random
//! flanks, the read sets of `tests/test_rfc1_motif.py`, as unsorted SAM:
//!
//! - `benign`: reads across an expanded AAAAG allele, and reads whose AAAAG
//!   part is soft-clipped;
//! - `pathogenic`: AAGGG reads soft-clipped next to the repeat, unmapped and
//!   placed at their mate, mates mapped to chr1, and reverse-strand reads
//!   across the repeat; three reads of an expanded AAAAG allele; eight reads
//!   across the unexpanded allele; a duplicate, a secondary, a QC-failed
//!   alignment, a homopolymer and a run of three units;
//! - `mixed` (AAAAG and AAAGG), `acagg`, `unknown` (AAGAG) and `few` (two
//!   AAGGG reads).
//!
//! `--samd12 UNIT` adds chr8 and, at SAMD12, five reads of UNIT x30 and six
//! reads of the reference allele (TAAAA)20 TGAAA with a random flank.

use super::{Fail, PyRandom, create, revcomp};
use std::io::Write;
use std::path::Path;

pub const CHROM: &str = "chr4";
pub const OTHER: &str = "chr1";
/// 0-based start and end (exclusive) of the repeat: 55 bases.
pub const START: usize = 39_348_424;
pub const END: usize = 39_348_479;
const KINDS: [&str; 6] = ["benign", "pathogenic", "mixed", "acagg", "unknown", "few"];

/// One alignment, as the SAM line needs it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Read {
    pub name: String,
    pub contig: &'static str,
    /// 0-based.
    pub start: usize,
    /// Empty for an unmapped read.
    pub cigar: String,
    pub seq: Vec<u8>,
    pub read1: bool,
    pub unmapped: bool,
    pub reverse: bool,
    pub dup: bool,
    pub secondary: bool,
    pub qcfail: bool,
    pub mate: (&'static str, usize),
}

impl Read {
    fn new(name: String, contig: &'static str, start: usize, cigar: &str, seq: Vec<u8>) -> Self {
        Read {
            name,
            contig,
            start,
            cigar: cigar.to_owned(),
            seq,
            read1: true,
            unmapped: false,
            reverse: false,
            dup: false,
            secondary: false,
            qcfail: false,
            mate: (contig, start),
        }
    }

    pub fn flag(&self) -> u16 {
        1 | if self.read1 { 64 } else { 128 }
            | if self.unmapped { 4 } else { 0 }
            | if self.reverse { 16 } else { 0 }
            | if self.secondary { 256 } else { 0 }
            | if self.qcfail { 512 } else { 0 }
            | if self.dup { 1024 } else { 0 }
    }

    /// The SAM line, with no tags and a newline.
    pub fn sam(&self) -> String {
        let mate = if self.mate.0 == self.contig {
            "="
        } else {
            self.mate.0
        };
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{mate}\t{}\t0\t{}\t{}\n",
            self.name,
            self.flag(),
            self.contig,
            self.start + 1,
            if self.unmapped { 0 } else { 60 },
            if self.cigar.is_empty() {
                "*"
            } else {
                &self.cigar
            },
            self.mate.1 + 1,
            String::from_utf8_lossy(&self.seq),
            "I".repeat(self.seq.len()),
        )
    }
}

/// The 300 random bases either side of the repeat: `random.Random(7)`.
pub struct Flanks {
    pub left: Vec<u8>,
    pub right: Vec<u8>,
}

impl Flanks {
    pub fn new() -> Self {
        let mut rng = PyRandom::new(7);
        let left = rng.bases(300);
        let right = rng.bases(300);
        Flanks { left, right }
    }

    fn left_tail(&self, n: usize) -> &[u8] {
        &self.left[self.left.len() - n..]
    }

    /// A 150-base read aligned across the repeat: `offset` flank bases, the
    /// repeat (55 reference bases plus an insertion of the rest), flank.
    fn spanning(&self, name: String, unit: &str, units: usize) -> Read {
        let offset = 40;
        let rep = unit.repeat(units);
        let right = 150 - offset - rep.len();
        let ins = rep.len() - (END - START);
        let mut cigar = format!("{offset}M{}M", END - START);
        if ins > 0 {
            cigar += &format!("{ins}I");
        }
        cigar += &format!("{right}M");
        let seq = [self.left_tail(offset), rep.as_bytes(), &self.right[..right]].concat();
        Read::new(name, CHROM, START - offset, &cigar, seq)
    }

    /// 60 aligned bases that end at the repeat start, 90 soft-clipped bases.
    fn clipped(&self, name: String, unit: &str) -> Read {
        self.clip(name, unit.repeat(18).as_bytes())
    }

    fn clip(&self, name: String, tail: &[u8]) -> Read {
        let seq = [self.left_tail(60), tail].concat();
        Read::new(name, CHROM, START - 60, "60M90S", seq)
    }

    /// (reads at the locus, mates elsewhere) of one kind.
    pub fn case(&self, kind: &str) -> (Vec<Read>, Vec<Read>) {
        let (mut at, mut elsewhere) = (Vec::new(), Vec::new());
        let n = |p: &str, i: usize| format!("{p}{i}");
        match kind {
            "benign" => {
                at.extend((0..6).map(|i| self.spanning(n("b", i), "AAAAG", 20)));
                at.extend((0..4).map(|i| self.clipped(n("bc", i), "AAAAG")));
            }
            "pathogenic" => {
                at.extend((0..6).map(|i| self.clipped(n("pc", i), "AAGGG")));
                for i in 0..2 {
                    let mut r = Read::new(
                        n("pu", i),
                        CHROM,
                        START - 100,
                        "",
                        "AAGGG".repeat(30).into_bytes(),
                    );
                    r.read1 = false;
                    r.unmapped = true;
                    at.push(r);
                }
                for i in 0..3 {
                    let mut a = Read::new(
                        n("pm", i),
                        CHROM,
                        START - 120,
                        "100M",
                        self.left_tail(120)[..100].to_vec(),
                    );
                    a.mate = (OTHER, 5000);
                    let mut b = Read::new(
                        n("pm", i),
                        OTHER,
                        5000,
                        "150M",
                        "AAGGG".repeat(30).into_bytes(),
                    );
                    b.read1 = false;
                    b.mate = (CHROM, START - 120);
                    at.push(a);
                    elsewhere.push(b);
                }
                for i in 0..2 {
                    let fwd = [self.left_tail(40), "AAGGG".repeat(22).as_bytes()].concat();
                    let mut r = Read::new(n("pr", i), CHROM, START - 40, "150M", revcomp(&fwd));
                    r.reverse = true;
                    at.push(r);
                }
                at.extend((0..3).map(|i| self.spanning(n("aaaag", i), "AAAAG", 20)));
                at.extend((0..8).map(|i| self.spanning(n("ref", i), "AAAAG", 11)));
                at.push(self.clipped("pc0".into(), "AAGGG"));
                let tail = "AAGGG".repeat(18);
                let mut dup = self.clip("dup".into(), tail.as_bytes());
                dup.dup = true;
                let mut sec = self.clip("sec".into(), tail.as_bytes());
                sec.secondary = true;
                let mut qc = self.clip("qc".into(), tail.as_bytes());
                qc.qcfail = true;
                at.extend([dup, sec, qc]);
                at.push(self.clip("poly".into(), &[b'A'; 90]));
                let short = ["AAGGG".repeat(3).as_bytes(), &self.right[..75]].concat();
                at.push(self.clip("short".into(), &short));
            }
            "mixed" => {
                at.extend((0..6).map(|i| self.spanning(n("m", i), "AAAAG", 20)));
                at.extend((0..4).map(|i| self.spanning(n("n", i), "AAAGG", 20)));
            }
            "acagg" => at.extend((0..5).map(|i| self.clipped(n("a", i), "ACAGG"))),
            "unknown" => at.extend((0..5).map(|i| self.clipped(n("u", i), "AAGAG"))),
            "few" => at.extend((0..2).map(|i| self.clipped(n("f", i), "AAGGG"))),
            _ => {}
        }
        (at, elsewhere)
    }
}

impl Default for Flanks {
    fn default() -> Self {
        Self::new()
    }
}

/// The SAMD12 reads: five of `unit` x30 and six of the reference allele,
/// with a 60-base flank of `random.Random(3)`.
fn samd12(unit: &str) -> String {
    let start = 118_366_812;
    let mut s = String::new();
    let seq = unit.repeat(30);
    for i in 0..5 {
        s += &format!(
            "f{i}\t0\tchr8\t{}\t60\t150M\t*\t0\t0\t{seq}\t{}\n",
            start + 1,
            "I".repeat(seq.len())
        );
    }
    let left = String::from_utf8(PyRandom::new(3).bases(60)).expect("ACGT");
    let allele = format!("{left}{}TGAAA{}", "TAAAA".repeat(20), &left[..30]);
    for i in 0..6 {
        s += &format!(
            "fr{i}\t0\tchr8\t{}\t60\t{}M\t*\t0\t0\t{allele}\t{}\n",
            start - 60 + 1,
            allele.len(),
            "I".repeat(allele.len())
        );
    }
    s
}

/// The whole SAM of one kind.
pub fn sam(kind: &str, samd12_unit: Option<&str>) -> String {
    let mut s = format!(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:{CHROM}\tLN:190214555\n@SQ\tSN:{OTHER}\tLN:248956422\n"
    );
    if samd12_unit.is_some() {
        s += "@SQ\tSN:chr8\tLN:145138636\n";
    }
    let (at, elsewhere) = Flanks::new().case(kind);
    for r in at.iter().chain(&elsewhere) {
        s += &r.sam();
    }
    if let Some(unit) = samd12_unit {
        s += &samd12(unit);
    }
    s
}

pub fn run(args: &[String]) -> Result<(), Fail> {
    let (unit, rest) = match args {
        [k, o, flag, u] if flag == "--samd12" => (Some(u.as_str()), [k, o]),
        [k, o] => (None, [k, o]),
        _ => {
            return Err(Fail::Usage(
                "rfc1 takes KIND OUT.sam [--samd12 UNIT]".into(),
            ));
        }
    };
    if !KINDS.contains(&rest[0].as_str()) {
        return Err(Fail::Usage(format!(
            "rfc1 KIND '{}' is not one of {}",
            rest[0],
            KINDS.join(", ")
        )));
    }
    if unit.is_some_and(|u| u.is_empty() || !u.bytes().all(|b| b"ACGT".contains(&b))) {
        return Err(Fail::Usage("--samd12 UNIT is bases of ACGT".into()));
    }
    let mut out = create(Path::new(rest[1]))?;
    out.write_all(sam(rest[0], unit).as_bytes())?;
    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixture::e2e::tests::fnv;

    // The pins are of the bytes the inline Python of
    // tests/e2e/expansions-and-rfc1-1-motif.sh and str-catalog-1-strchive.sh
    // wrote, from tests/test_rfc1_motif.py's case().
    #[test]
    fn writes_the_pythons_bytes() {
        assert_eq!(fnv(sam("pathogenic", None).as_bytes()), PIN_PATHOGENIC);
        assert_eq!(fnv(sam("benign", None).as_bytes()), PIN_BENIGN);
        assert_eq!(fnv(sam("pathogenic", Some("TTTCA")).as_bytes()), PIN_SAMD12);
    }

    #[test]
    fn the_pathogenic_case_has_what_the_e2e_case_counts() {
        let (at, elsewhere) = Flanks::new().case("pathogenic");
        let clipped = at
            .iter()
            .filter(|r| r.cigar == "60M90S" && r.seq.ends_with(b"AAGGGAAGGG"))
            .count();
        // 6 soft-clipped plus the same read again, the duplicate, the secondary and the QC-failed one.
        assert_eq!(clipped, 10);
        assert_eq!(at.iter().filter(|r| r.unmapped).count(), 2);
        assert_eq!(elsewhere.len(), 3);
        assert!(
            elsewhere
                .iter()
                .all(|r| r.contig == OTHER && r.mate == (CHROM, START - 120))
        );
        assert_eq!(at.iter().filter(|r| r.reverse).count(), 2);
        assert!(at.iter().all(|r| r.seq.len() == 150 || r.cigar == "100M"));
        let line = at.iter().find(|r| r.name == "pu0").unwrap().sam();
        assert!(
            line.starts_with("pu0\t133\tchr4\t39348325\t0\t*\t=\t39348325\t0\t"),
            "{line}"
        );
        let line = at.iter().find(|r| r.name == "pm0").unwrap().sam();
        assert!(
            line.starts_with("pm0\t65\tchr4\t39348305\t60\t100M\tchr1\t5001\t0\t"),
            "{line}"
        );
    }

    #[test]
    fn spanning_reads_insert_what_the_reference_lacks() {
        let f = Flanks::new();
        assert_eq!(f.spanning("x".into(), "AAAAG", 20).cigar, "40M55M45I10M");
        assert_eq!(f.spanning("x".into(), "AAAAG", 11).cigar, "40M55M55M");
    }

    #[test]
    fn bad_arguments() {
        use crate::fixture::e2e::tests::run;
        assert_eq!(run("rfc1", &["nope", "/dev/null"]).0, Some(2));
        assert_eq!(
            run("rfc1", &["benign", "/dev/null", "--samd12", "XY"]).0,
            Some(2)
        );
        assert_eq!(run("rfc1", &["benign"]).0, Some(2));
    }

    const PIN_PATHOGENIC: u64 = 7_022_458_020_176_947_217;
    const PIN_BENIGN: u64 = 5_018_501_953_293_861_024;
    const PIN_SAMD12: u64 = 5_496_064_447_239_138_924;
}
