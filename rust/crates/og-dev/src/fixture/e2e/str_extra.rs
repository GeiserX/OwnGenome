//! `str-extra REF.fa HET_NAME HET.sam HOM_NAME HOM.sam`: the made-up
//! reference and reads of the two cases on the repeats STRchive leaves out
//! (`assets/str_extra_loci.json`).
//!
//! REF.fa holds chr7 and chr2 of the given lengths, N everywhere except 1500
//! random bases (`random.Random(11)`) either side of the GRCh38 tract at its
//! GRCh38 place: CFTR TG11 T7 and UGT1A1 A(TA)7. Each sample's SAM holds
//! read pairs every 4 bases from each of its two haplotypes, the CIGAR
//! putting the tract's length change at the tract's end:
//!
//! - HET: CFTR TG12-T5 and TG11-T7, UGT1A1 TA7 and TA8;
//! - HOM: CFTR TG11-T7 and UGT1A1 TA7 on both copies.

use super::{Fail, PyRandom, create, exactly};
use std::fmt::Write as _;
use std::io::{self, Write};
use std::path::Path;

const W: usize = 1500;

/// A locus: contig, contig length, 0-based tract start, the bases before
/// and after the tract, and the reference tract.
struct Locus {
    contig: &'static str,
    length: usize,
    start: usize,
    left: &'static str,
    right: &'static str,
    tract: &'static [(&'static str, usize)],
}

const LOCI: [Locus; 2] = [
    Locus {
        contig: "chr7",
        length: 117_551_000,
        start: 117_548_606,
        left: "TTTTGA",
        right: "AACAGG",
        tract: &[("TG", 11), ("T", 7)],
    },
    Locus {
        contig: "chr2",
        length: 233_763_000,
        start: 233_760_234,
        left: "TTGCCA",
        right: "AGTAGG",
        tract: &[("TA", 7)],
    },
];

type Hap = &'static [(&'static str, usize)];
const HET: [[Hap; 2]; 2] = [
    [&[("TG", 12), ("T", 5)], &[("TG", 11), ("T", 7)]],
    [&[("TA", 7)], &[("TA", 8)]],
];
const HOM: [[Hap; 2]; 2] = [
    [&[("TG", 11), ("T", 7)], &[("TG", 11), ("T", 7)]],
    [&[("TA", 7)], &[("TA", 7)]],
];

fn tract(units: Hap) -> String {
    units.iter().map(|(m, n)| m.repeat(*n)).collect()
}

/// The random flanks of each locus, in order.
fn flanks() -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut rng = PyRandom::new(11);
    LOCI.iter()
        .map(|l| {
            let mut lf = rng.bases(W - l.left.len());
            lf.extend_from_slice(l.left.as_bytes());
            let mut rf = l.right.as_bytes().to_vec();
            rf.extend(rng.bases(W - l.right.len()));
            (lf, rf)
        })
        .collect()
}

/// Writes FASTA lines of 60, a run of `N`s or bases at a time.
pub(crate) struct Lines<'a> {
    w: &'a mut dyn Write,
    col: usize,
}

impl<'a> Lines<'a> {
    pub(crate) fn new(w: &'a mut dyn Write) -> Self {
        Lines { w, col: 0 }
    }

    pub(crate) fn ns(&mut self, n: usize) -> io::Result<()> {
        self.col = super::n_lines(self.w, n, self.col)?;
        Ok(())
    }

    pub(crate) fn bases(&mut self, mut s: &[u8]) -> io::Result<()> {
        while !s.is_empty() {
            let k = s.len().min(60 - self.col);
            self.w.write_all(&s[..k])?;
            self.col += k;
            s = &s[k..];
            if self.col == 60 {
                self.w.write_all(b"\n")?;
                self.col = 0;
            }
        }
        Ok(())
    }

    /// Ends a partial last line.
    pub(crate) fn end(&mut self) -> io::Result<()> {
        if self.col > 0 {
            self.w.write_all(b"\n")?;
            self.col = 0;
        }
        Ok(())
    }
}

fn reference(w: &mut dyn Write, flanks: &[(Vec<u8>, Vec<u8>)]) -> io::Result<()> {
    for (l, (lf, rf)) in LOCI.iter().zip(flanks) {
        writeln!(w, ">{}", l.contig)?;
        let t = tract(l.tract);
        let mut lines = Lines::new(w);
        lines.ns(l.start - W)?;
        lines.bases(lf)?;
        lines.bases(t.as_bytes())?;
        lines.bases(rf)?;
        lines.ns(l.length - (l.start - W + lf.len() + t.len() + rf.len()))?;
        lines.end()?;
    }
    Ok(())
}

/// The 1-based POS and CIGAR of haplotype bases [a, a+150), or `None` when
/// the read starts or ends inside the tract.
fn place(
    a: usize,
    w0: usize,
    t0: usize,
    t1: usize,
    hap_len: usize,
    ref_len: usize,
) -> Option<(usize, String)> {
    const N: usize = 150;
    if a + N <= t0 {
        return Some((w0 + a + 1, format!("{N}M")));
    }
    if a >= t1 {
        return Some(((w0 + a + 1 + ref_len) - hap_len, format!("{N}M")));
    }
    if a >= t0 || a + N <= t1 {
        return None;
    }
    let pre = t0 - a + hap_len.min(ref_len);
    let post = N - (t0 - a) - hap_len;
    let cigar = match hap_len.cmp(&ref_len) {
        std::cmp::Ordering::Equal => format!("{N}M"),
        std::cmp::Ordering::Greater => format!("{pre}M{}I{post}M", hap_len - ref_len),
        std::cmp::Ordering::Less => format!("{pre}M{}D{post}M", ref_len - hap_len),
    };
    Some((w0 + a + 1, cigar))
}

fn reads(out: &mut String, l: &Locus, (lf, rf): &(Vec<u8>, Vec<u8>), haps: &[Hap; 2], tag: &str) {
    let w0 = l.start - W;
    let ref_len = tract(l.tract).len();
    let q = "I".repeat(150);
    for (h, units) in haps.iter().enumerate() {
        let t = tract(units);
        let seq = [lf.as_slice(), t.as_bytes(), rf.as_slice()].concat();
        let (t0, t1) = (lf.len(), lf.len() + t.len());
        let mut k = 0;
        for a in (0..seq.len() - 400).step_by(4) {
            let b = a + 300 + (a % 7) * 10 - 150;
            let (Some((p1, c1)), Some((p2, c2))) = (
                place(a, w0, t0, t1, t.len(), ref_len),
                place(b, w0, t0, t1, t.len(), ref_len),
            ) else {
                continue;
            };
            k += 1;
            let name = format!("{tag}{}h{h}r{k}", l.contig);
            let (r1, r2) = (
                String::from_utf8_lossy(&seq[a..a + 150]),
                String::from_utf8_lossy(&seq[b..b + 150]),
            );
            let tlen = p2 + 150 - p1;
            let c = l.contig;
            let _ = writeln!(
                out,
                "{name}\t99\t{c}\t{p1}\t60\t{c1}\t=\t{p2}\t{tlen}\t{r1}\t{q}\tRG:Z:{tag}"
            );
            let _ = writeln!(
                out,
                "{name}\t147\t{c}\t{p2}\t60\t{c2}\t=\t{p1}\t-{tlen}\t{r2}\t{q}\tRG:Z:{tag}"
            );
        }
    }
}

/// The SAM of one sample.
fn sam(name: &str, haps: &[[Hap; 2]; 2], flanks: &[(Vec<u8>, Vec<u8>)]) -> String {
    let mut s = String::from("@HD\tVN:1.6\tSO:unsorted\n");
    for l in &LOCI {
        let _ = writeln!(s, "@SQ\tSN:{}\tLN:{}", l.contig, l.length);
    }
    let _ = writeln!(s, "@RG\tID:{name}\tSM:{name}");
    for ((l, f), h) in LOCI.iter().zip(flanks).zip(haps) {
        reads(&mut s, l, f, h, name);
    }
    s
}

pub fn run(args: &[String]) -> Result<(), Fail> {
    let a = exactly(args, 5, "REF.fa HET_NAME HET.sam HOM_NAME HOM.sam")?;
    let f = flanks();
    let mut out = create(Path::new(&a[0]))?;
    reference(&mut out, &f)?;
    out.flush()?;
    for (name, path, haps) in [(&a[1], &a[2], &HET), (&a[3], &a[4], &HOM)] {
        let mut out = create(Path::new(path))?;
        out.write_all(sam(name, haps, &f).as_bytes())?;
        out.flush()?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixture::e2e::tests::{Hashed, fnv};

    // Pins of the bytes the inline Python of
    // tests/e2e/genotype-rule-rows-1-str-extra.sh wrote (samples EXA, EXB).
    #[test]
    fn writes_the_pythons_reads() {
        let f = flanks();
        assert_eq!(fnv(sam("EXA", &HET, &f).as_bytes()), PIN_EXA);
        assert_eq!(fnv(sam("EXB", &HOM, &f).as_bytes()), PIN_EXB);
    }

    #[test]
    fn the_window_holds_the_grch38_tract_at_its_place() {
        let f = flanks();
        let (lf, rf) = &f[0];
        let win = [
            lf.as_slice(),
            tract(LOCI[0].tract).as_bytes(),
            rf.as_slice(),
        ]
        .concat();
        // chr7:117548601-117548641, 1-based, starts 5 bases before the tract.
        let at = 117_548_600 - (LOCI[0].start - W);
        assert_eq!(
            &win[at..at + 41],
            b"TTTTGATGTGTGTGTGTGTGTGTGTGTGTTTTTTTAACAGG"
        );
    }

    #[test]
    fn the_reference_is_n_outside_the_windows_and_60_a_line() {
        // Both contigs are 351 MB: hash them on the way, line by line.
        let mut h = Hashed::default();
        reference(&mut h, &flanks()).unwrap();
        assert_eq!(
            h.lines,
            2 + 117_551_000_usize.div_ceil(60) + 233_763_000_usize.div_ceil(60)
        );
        assert_eq!(h.longest, 60);
        assert_eq!(h.hash, PIN_REF);
    }

    #[test]
    fn a_read_across_a_longer_tract_carries_the_insertion() {
        // TG12-T5 is 29 bases against the reference's 29: no change; TA8 is 2 longer.
        assert_eq!(
            place(1400, 0, 1500, 1516, 16, 14),
            Some((1401, "114M2I34M".into()))
        );
        assert_eq!(
            place(1400, 0, 1500, 1512, 12, 14),
            Some((1401, "112M2D38M".into()))
        );
        assert_eq!(
            place(1520, 0, 1500, 1516, 16, 14),
            Some((1519, "150M".into()))
        );
        assert_eq!(place(1505, 0, 1500, 1516, 16, 14), None);
    }

    const PIN_EXA: u64 = 12_879_245_652_337_382_645;
    const PIN_EXB: u64 = 10_035_011_117_999_463_469;
    const PIN_REF: u64 = 8_076_056_084_048_732_156;
}
