//! Reads at a variant: what one alignment shows at a site, where it is
//! clipped, and the pieces its SA tag names.
//!
//! The three read tools (`og cdna check`, `og cdna scan`, `og comphet
//! phase`) look at a handful of bases of each read. [`allele_at`] walks the
//! CIGAR once, keeps only the bases of the site's window and stops past it,
//! so a read costs no map of its bases. [`ReadView`] is what the tools read
//! of an alignment: a BAM [`Record`], or a [`SamRead`] parsed from SAM text,
//! so the unit tests feed reads without a BAM.

use pycompat::num::parse_int;
use pycompat::text::{is_digit, isalpha, isdigit, upper};
use rust_htslib::bam::Record;
use rust_htslib::bam::record::Aux;

/// CIGAR operation codes, as the BAM format numbers them.
pub mod op {
    pub const M: u8 = 0;
    pub const INS: u8 = 1;
    pub const DEL: u8 = 2;
    pub const SKIP: u8 = 3;
    pub const SOFT: u8 = 4;
    pub const HARD: u8 = 5;
    pub const PAD: u8 = 6;
    pub const EQ: u8 = 7;
    pub const DIFF: u8 = 8;
}

/// SAM flag bits the tools read.
pub mod flag {
    pub const UNMAPPED: u16 = 0x4;
    pub const REVERSE: u16 = 0x10;
    pub const READ1: u16 = 0x40;
    pub const READ2: u16 = 0x80;
    pub const SECONDARY: u16 = 0x100;
    pub const QCFAIL: u16 = 0x200;
    pub const DUPLICATE: u16 = 0x400;
    pub const SUPPLEMENTARY: u16 = 0x800;
}

const OP_CHARS: &[u8; 9] = b"MIDNSHP=X";

/// One CIGAR operation: its code ([`op`]) and length.
pub type Op = (u8, i64);

/// The operations of a raw BAM CIGAR.
pub fn ops(raw: &[u32]) -> impl Iterator<Item = Op> + Clone + '_ {
    raw.iter().map(|&c| ((c & 0xf) as u8, i64::from(c >> 4)))
}

/// Whether `code` moves along the reference (M, D, N, = and X).
pub fn on_reference(code: u8) -> bool {
    matches!(code, op::M | op::DEL | op::SKIP | op::EQ | op::DIFF)
}

/// What the read tools read of one alignment.
pub trait ReadView {
    fn flags(&self) -> u16;
    fn mapq(&self) -> u8;
    /// The 0-based leftmost reference position.
    fn pos(&self) -> i64;
    /// The CIGAR as BAM stores it: `len << 4 | op`.
    fn raw_cigar(&self) -> &[u32];
    fn seq_len(&self) -> usize;
    /// The read base at `i`, `None` past the sequence's end.
    fn base(&self, i: usize) -> Option<u8>;
    fn qname(&self) -> &[u8];
    /// The SA tag's text; `None` when the read has none, or one that is not text.
    fn sa(&self) -> Option<&str>;

    fn has(&self, bits: u16) -> bool {
        self.flags() & bits != 0
    }
}

impl ReadView for Record {
    fn flags(&self) -> u16 {
        Record::flags(self)
    }
    fn mapq(&self) -> u8 {
        Record::mapq(self)
    }
    fn pos(&self) -> i64 {
        Record::pos(self)
    }
    fn raw_cigar(&self) -> &[u32] {
        Record::raw_cigar(self)
    }
    fn seq_len(&self) -> usize {
        Record::seq_len(self)
    }
    fn base(&self, i: usize) -> Option<u8> {
        (i < Record::seq_len(self)).then(|| self.seq()[i])
    }
    fn qname(&self) -> &[u8] {
        Record::qname(self)
    }
    fn sa(&self) -> Option<&str> {
        match self.aux(b"SA") {
            Ok(Aux::String(s)) => Some(s),
            _ => None,
        }
    }
}

/// Which alignments a tool counts. Unmapped, secondary, duplicate and
/// QC-failed alignments and those without a CIGAR never count.
#[derive(Debug, Clone, Copy)]
pub struct Filter {
    /// comphet skips supplementary alignments; the cDNA tools keep them,
    /// since a construct's junction piece is often one.
    pub skip_supplementary: bool,
    pub min_mapq: u8,
}

impl Filter {
    /// `og cdna check` and `og cdna scan`: supplementary alignments count, any MAPQ.
    pub const CDNA: Filter = Filter {
        skip_supplementary: false,
        min_mapq: 0,
    };
    /// `og comphet phase`: primary alignments of MAPQ 20 or more.
    pub const PHASE: Filter = Filter {
        skip_supplementary: true,
        min_mapq: 20,
    };

    pub fn keeps(&self, r: &impl ReadView) -> bool {
        let never = flag::UNMAPPED | flag::SECONDARY | flag::DUPLICATE | flag::QCFAIL;
        !(r.has(never)
            || r.raw_cigar().is_empty()
            || (self.skip_supplementary && r.has(flag::SUPPLEMENTARY))
            || r.mapq() < self.min_mapq)
    }
}

/// The allele a read shows at a site.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Seen {
    Ref,
    Alt,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum Shape {
    /// SNV or MNV: the read's bases over REF.
    Same,
    /// REF longer: the bases after the anchor are deleted.
    Del,
    /// ALT longer: these bases (upper case) are inserted after the anchor.
    Ins(String),
    /// The anchors differ: no read shows either allele.
    Never,
}

/// A VCF allele pair at a 0-based position, ready to compare with reads.
#[derive(Debug, Clone)]
pub struct Site {
    pub pos0: i64,
    ref_len: i64,
    /// The bases before the gap: min(len(REF), len(ALT)).
    anchor: i64,
    ref_up: String,
    alt_up: String,
    shape: Shape,
}

impl Site {
    /// Lengths count code points, as Python's `len` does.
    pub fn new(pos0: i64, reference: &str, alt: &str) -> Site {
        let (r, a) = (reference.chars().count(), alt.chars().count());
        let n = r.min(a);
        let head = |s: &str| upper(&s.chars().take(n).collect::<String>());
        let shape = if r == a {
            Shape::Same
        } else if head(reference) != head(alt) {
            Shape::Never
        } else if r > a {
            Shape::Del
        } else {
            Shape::Ins(upper(&alt.chars().skip(n).collect::<String>()))
        };
        Site {
            pos0,
            ref_len: r as i64,
            anchor: n as i64,
            ref_up: upper(reference),
            alt_up: upper(alt),
            shape,
        }
    }

    /// The reference window a read must reach: REF and the base after it.
    pub fn window(&self) -> (i64, i64) {
        (self.pos0, self.pos0 + self.ref_len.max(1) + 1)
    }
}

/// Not a SNV, MNV or an indel that shares its first bases: such an allele is
/// not looked for in the reads.
pub fn complex_allele(reference: &str, alt: &str) -> bool {
    let (r, a) = (reference.chars().count(), alt.chars().count());
    let n = r.min(a);
    let head = |s: &str| upper(&s.chars().take(n).collect::<String>());
    (r != a && head(reference) != head(alt)) || !isalpha(reference) || !isalpha(alt)
}

/// One reference position of a site's window as a read covers it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Cell {
    Absent,
    Deleted,
    Base(u8),
}

/// What a read shows at `site`: [`Seen::Alt`], [`Seen::Ref`], or `None`
/// when it shows neither or does not span the site. SNV and MNV compare
/// the bases; a deletion is seen when every gap base is deleted and the
/// base after it aligned with no insertion after the anchor, its reference
/// when every gap base is aligned; an insertion by the bases inserted after
/// the last anchor base, its reference when none are.
pub fn allele_at(r: &impl ReadView, site: &Site) -> Option<Seen> {
    if r.seq_len() == 0 || site.shape == Shape::Never {
        return None;
    }
    let lo = site.pos0;
    let hi = site.pos0 + site.ref_len; // inclusive: the base after REF
    let width = usize::try_from(hi - lo + 1).ok()?;
    let mut cells = vec![Cell::Absent; width];
    let mut ins: Vec<Option<Vec<u8>>> = vec![None; width];
    let at = |p: i64| usize::try_from(p - lo).ok().filter(|&i| i < width);
    let (mut rp, mut q) = (r.pos(), 0usize);
    for (code, len) in ops(r.raw_cigar()) {
        if rp > hi + 1 {
            break;
        }
        let n = usize::try_from(len).ok()?;
        match code {
            op::M | op::EQ | op::DIFF => {
                for p in rp.max(lo)..(rp + len).min(hi + 1) {
                    let i = usize::try_from(p - rp).ok()?;
                    cells[at(p)?] = Cell::Base(r.base(q + i)?);
                }
                rp += len;
                q += n;
            }
            op::INS => {
                if let Some(k) = at(rp - 1) {
                    let got = ins[k].get_or_insert_with(Vec::new);
                    for i in 0..n {
                        got.push(r.base(q + i)?);
                    }
                }
                q += n;
            }
            op::DEL | op::SKIP => {
                for p in rp.max(lo)..(rp + len).min(hi + 1) {
                    cells[at(p)?] = Cell::Deleted;
                }
                rp += len;
            }
            op::SOFT => q += n,
            _ => {}
        }
    }
    let cell = |p: i64| at(p).map_or(Cell::Absent, |i| cells[i]);
    let base = |p: i64| match cell(p) {
        Cell::Base(b) => Some(b),
        _ => None,
    };
    let bases = |from: i64, to: i64| -> Option<String> {
        (from..to)
            .map(|p| base(p).map(|b| char::from(b.to_ascii_uppercase())))
            .collect()
    };
    match &site.shape {
        Shape::Never => None,
        Shape::Same => {
            let got = bases(lo, lo + site.ref_len)?;
            if got == site.alt_up {
                Some(Seen::Alt)
            } else if got == site.ref_up {
                Some(Seen::Ref)
            } else {
                None
            }
        }
        shape => {
            let last = lo + site.anchor - 1;
            bases(lo, lo + site.anchor)?;
            if let Shape::Ins(inserted) = shape {
                base(last + 1)?;
                let got: String = ins[at(last)?]
                    .as_deref()
                    .unwrap_or_default()
                    .iter()
                    .map(|b| char::from(b.to_ascii_uppercase()))
                    .collect();
                return if &got == inserted {
                    Some(Seen::Alt)
                } else if got.is_empty() {
                    Some(Seen::Ref)
                } else {
                    None
                };
            }
            if ins[at(last)?].is_some() || base(hi).is_none() {
                return None;
            }
            let gap = last + 1..hi;
            if gap.clone().all(|p| cell(p) == Cell::Deleted) {
                Some(Seen::Alt)
            } else if gap.clone().all(|p| base(p).is_some()) {
                Some(Seen::Ref)
            } else {
                None
            }
        }
    }
}

/// One alignment's reference span and its clips: 0-based, end exclusive,
/// the clips counting soft and hard clipped bases.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Segment {
    pub start: i64,
    pub end: i64,
    pub left_clip: i64,
    pub right_clip: i64,
}

/// The [`Segment`] of an alignment at `start` with these operations.
pub fn segment(start: i64, cigar: impl Iterator<Item = Op> + Clone) -> Segment {
    let clip = |&(code, _): &Op| matches!(code, op::SOFT | op::HARD);
    let ref_len: i64 = cigar
        .clone()
        .filter(|&(c, _)| on_reference(c))
        .map(|(_, n)| n)
        .sum();
    let left = cigar.clone().take_while(clip).map(|(_, n)| n).sum();
    let ops: Vec<Op> = cigar.collect();
    let right = ops
        .iter()
        .rev()
        .take_while(|o| clip(o))
        .map(|(_, n)| n)
        .sum();
    Segment {
        start,
        end: start + ref_len,
        left_clip: left,
        right_clip: right,
    }
}

/// A CIGAR string as text: `75M25S` is `[(M, 75), (S, 25)]`. Trailing
/// digits without an operation are ignored; an unknown operation or one
/// without a length is an error.
pub fn parse_cigar(text: &str) -> Result<Vec<Op>, String> {
    let mut out = Vec::new();
    let mut n = String::new();
    for ch in text.chars() {
        if is_digit(ch) {
            n.push(ch);
            continue;
        }
        let code = u8::try_from(ch)
            .ok()
            .and_then(|c| OP_CHARS.iter().position(|&o| o == c))
            .ok_or_else(|| format!("unknown CIGAR operation {ch:?} in {text:?}"))?;
        let len = parse_int(&n)
            .ok()
            .and_then(|v| i64::try_from(v).ok())
            .ok_or_else(|| format!("CIGAR operation without a length in {text:?}"))?;
        out.push((code as u8, len));
        n.clear();
    }
    Ok(out)
}

/// One alignment an SA tag names.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SaPiece<'a> {
    pub contig: &'a str,
    /// 0-based.
    pub pos0: i64,
    pub reverse: bool,
    pub cigar: Vec<Op>,
}

/// The pieces of an SA tag `chr,pos,strand,CIGAR,mapq,nm;...`. A part with
/// fewer than four fields, or a position that is not digits, is skipped.
pub fn sa_pieces(sa: &str) -> impl Iterator<Item = Result<SaPiece<'_>, String>> {
    sa.split(';').filter_map(|part| {
        let f: Vec<&str> = part.split(',').collect();
        if f.len() < 4 || !isdigit(f[1]) {
            return None;
        }
        let piece = parse_int(f[1])
            .ok()
            .and_then(|v| i64::try_from(v - 1).ok())
            .ok_or_else(|| format!("SA position {:?} out of range", f[1]))
            .and_then(|pos0| {
                Ok(SaPiece {
                    contig: f[0],
                    pos0,
                    reverse: f[2] == "-",
                    cigar: parse_cigar(f[3])?,
                })
            });
        Some(piece)
    })
}

/// A read parsed from one SAM line (no header): what the unit tests feed
/// the tools instead of a BAM.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SamRead {
    pub qname: Vec<u8>,
    pub flags: u16,
    pub contig: String,
    pub pos: i64,
    pub mapq: u8,
    pub cigar: Vec<u32>,
    pub seq: Vec<u8>,
    pub sa: Option<String>,
}

impl SamRead {
    /// Reads `QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
    /// [TAG...]`; POS is 1-based, `*` is an empty CIGAR or sequence.
    pub fn parse(line: &str) -> Result<SamRead, String> {
        let f: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        if f.len() < 11 {
            return Err(format!("not a SAM line: {line:?}"));
        }
        let num = |s: &str| s.parse::<i64>().map_err(|_| format!("not a number: {s:?}"));
        let cigar = if f[5] == "*" {
            Vec::new()
        } else {
            parse_cigar(f[5])?
                .into_iter()
                .map(|(c, n)| u32::try_from(n).map(|n| n << 4 | u32::from(c)))
                .collect::<Result<_, _>>()
                .map_err(|_| format!("CIGAR too long: {}", f[5]))?
        };
        Ok(SamRead {
            qname: f[0].as_bytes().to_vec(),
            flags: u16::try_from(num(f[1])?).map_err(|_| format!("bad flag {}", f[1]))?,
            contig: f[2].to_owned(),
            pos: num(f[3])? - 1,
            mapq: u8::try_from(num(f[4])?).map_err(|_| format!("bad MAPQ {}", f[4]))?,
            cigar,
            seq: if f[9] == "*" {
                Vec::new()
            } else {
                f[9].as_bytes().to_vec()
            },
            sa: f[11..]
                .iter()
                .find_map(|t| t.strip_prefix("SA:Z:"))
                .map(str::to_owned),
        })
    }

    /// One past the last reference base, as htslib's `bam_endpos` (a span
    /// of no reference bases counts as 1).
    pub fn end(&self) -> i64 {
        let len: i64 = ops(&self.cigar)
            .filter(|&(c, _)| on_reference(c))
            .map(|(_, n)| n)
            .sum();
        self.pos + len.max(1)
    }
}

impl ReadView for SamRead {
    fn flags(&self) -> u16 {
        self.flags
    }
    fn mapq(&self) -> u8 {
        self.mapq
    }
    fn pos(&self) -> i64 {
        self.pos
    }
    fn raw_cigar(&self) -> &[u32] {
        &self.cigar
    }
    fn seq_len(&self) -> usize {
        self.seq.len()
    }
    fn base(&self, i: usize) -> Option<u8> {
        self.seq.get(i).copied()
    }
    fn qname(&self) -> &[u8] {
        &self.qname
    }
    fn sa(&self) -> Option<&str> {
        self.sa.as_deref()
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// A made-up contig: the same bases every run (the clip tests' genome).
    pub(crate) fn genome(n: usize) -> Vec<u8> {
        let mut s: u64 = 7;
        (0..n)
            .map(|_| {
                s = (s * 75 + 74) % 65537;
                b"ACGT"[(s % 4) as usize]
            })
            .collect()
    }

    /// A mapped read at 0-based `pos` with a CIGAR and bases.
    pub(crate) fn read(name: &str, pos: i64, cigar: &str, seq: &[u8]) -> SamRead {
        let line = format!(
            "{name}\t0\tchrT\t{}\t60\t{cigar}\t*\t0\t0\t{}\t*",
            pos + 1,
            if seq.is_empty() {
                "*".to_owned()
            } else {
                String::from_utf8_lossy(seq).into_owned()
            }
        );
        SamRead::parse(&line).unwrap()
    }

    fn s(b: &[u8]) -> String {
        String::from_utf8_lossy(b).into_owned()
    }

    #[test]
    fn a_cigar_parses_and_its_segment_has_both_clips() {
        let c = parse_cigar("5S70M2I3D20M10H").unwrap();
        assert_eq!(
            c,
            vec![
                (op::SOFT, 5),
                (op::M, 70),
                (op::INS, 2),
                (op::DEL, 3),
                (op::M, 20),
                (op::HARD, 10)
            ]
        );
        assert_eq!(
            segment(100, c.into_iter()),
            Segment {
                start: 100,
                end: 193,
                left_clip: 5,
                right_clip: 10
            }
        );
        assert!(parse_cigar("5Q").is_err());
        assert!(parse_cigar("M").is_err());
        assert_eq!(parse_cigar("10M5").unwrap(), vec![(op::M, 10)]);
    }

    #[test]
    fn sa_positions_turn_0_based_and_odd_parts_are_skipped() {
        let got: Vec<_> = sa_pieces("chr2,3001,-,80S20M,60,0;bad;chr3,x,+,5M,1,0;")
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(
            got,
            vec![SaPiece {
                contig: "chr2",
                pos0: 3000,
                reverse: true,
                cigar: vec![(op::SOFT, 80), (op::M, 20)]
            }]
        );
    }

    #[test]
    fn snv_alleles() {
        let g = genome(10_000);
        let pos0 = 1150;
        let (r, a) = (
            g[pos0],
            b"ACGT"[(b"ACGT".iter().position(|&c| c == g[pos0]).unwrap() + 1) % 4],
        );
        let site = Site::new(pos0 as i64, &s(&[r]), &s(&[a]));
        let seq = g[1100..1200].to_vec();
        assert_eq!(
            allele_at(&read("r", 1100, "100M", &seq), &site),
            Some(Seen::Ref)
        );
        let mut alt = seq.clone();
        alt[50] = a;
        assert_eq!(
            allele_at(&read("a", 1100, "100M", &alt), &site),
            Some(Seen::Alt)
        );
        assert_eq!(allele_at(&read("s", 1100, "40M", &seq[..40]), &site), None);
        // Lower case bases read as upper case.
        assert_eq!(
            allele_at(&read("l", 1100, "100M", &alt.to_ascii_lowercase()), &site),
            Some(Seen::Alt)
        );
        // No stored sequence: nothing seen.
        assert_eq!(allele_at(&read("n", 1100, "100M", b""), &site), None);
    }

    #[test]
    fn deletion_alleles() {
        let g = genome(10_000);
        let p = 1150;
        let site = Site::new(p as i64, &s(&g[p..p + 4]), &s(&g[p..=p]));
        let del = [&g[1100..=p], &g[p + 4..1200]].concat();
        assert_eq!(
            allele_at(&read("d", 1100, "51M3D46M", &del), &site),
            Some(Seen::Alt)
        );
        assert_eq!(
            allele_at(&read("r", 1100, "100M", &g[1100..1200]), &site),
            Some(Seen::Ref)
        );
        // A read that stops inside the gap shows neither.
        assert_eq!(
            allele_at(&read("e", 1100, "52M", &g[1100..1152]), &site),
            None
        );
        // A deletion that does not cover the whole gap shows neither.
        let part = [&g[1100..=p], &g[p + 3..1200]].concat();
        assert_eq!(allele_at(&read("p", 1100, "51M2D47M", &part), &site), None);
        // A deletion with the base after it unaligned shows neither.
        assert_eq!(
            allele_at(&read("t", 1100, "51M3D", &g[1100..=p]), &site),
            None
        );
        // An insertion right after the anchor makes it unreadable.
        let ins = [&g[1100..=p], b"TT", &g[p + 4..1200]].concat();
        assert_eq!(allele_at(&read("i", 1100, "51M2I3D46M", &ins), &site), None);
    }

    #[test]
    fn insertion_alleles() {
        let g = genome(10_000);
        let p = 1150;
        let r = s(&g[p..=p]);
        let site = Site::new(p as i64, &r, &format!("{r}TT"));
        let with = [&g[1100..=p], b"TT", &g[p + 1..1149 + 51]].concat();
        assert_eq!(
            allele_at(&read("i", 1100, "51M2I49M", &with), &site),
            Some(Seen::Alt)
        );
        assert_eq!(
            allele_at(&read("r", 1100, "100M", &g[1100..1200]), &site),
            Some(Seen::Ref)
        );
        let other = [&g[1100..=p], b"GA", &g[p + 1..1149 + 51]].concat();
        assert_eq!(allele_at(&read("o", 1100, "51M2I49M", &other), &site), None);
        // A read that ends at the anchor does not show what follows it.
        assert_eq!(
            allele_at(&read("e", 1100, "51M", &g[1100..=p]), &site),
            None
        );
    }

    #[test]
    fn complex_alleles_are_not_looked_for() {
        assert!(complex_allele("AC", "GTT"));
        assert!(!complex_allele("A", "G"));
        assert!(!complex_allele("ACG", "A"));
        assert!(complex_allele("A", "<DEL>"));
        assert!(complex_allele("", "A"));
        assert!(complex_allele("ac", "AGT"));
        assert!(!complex_allele("ac", "ACT"));
    }

    #[test]
    fn the_filters() {
        let mut r = read("r", 10, "10M", b"AAAAAAAAAA");
        assert!(Filter::CDNA.keeps(&r) && Filter::PHASE.keeps(&r));
        r.flags = flag::SUPPLEMENTARY;
        assert!(Filter::CDNA.keeps(&r) && !Filter::PHASE.keeps(&r));
        for f in [
            flag::UNMAPPED,
            flag::SECONDARY,
            flag::DUPLICATE,
            flag::QCFAIL,
        ] {
            r.flags = f;
            assert!(!Filter::CDNA.keeps(&r));
        }
        r.flags = 0;
        r.mapq = 19;
        assert!(Filter::CDNA.keeps(&r) && !Filter::PHASE.keeps(&r));
        r.cigar.clear();
        assert!(!Filter::CDNA.keeps(&r));
    }

    #[test]
    fn a_sam_line_reads_back() {
        let r = SamRead::parse("q\t2131\tchrT\t101\t60\t60M40H\t*\t0\t0\tACGT\t*\tNM:i:0\tSA:Z:chrT,3001,+,60S40M,60,0;")
            .unwrap();
        assert_eq!((r.pos, r.end(), r.flags), (100, 160, 2131));
        assert_eq!(r.sa(), Some("chrT,3001,+,60S40M,60,0;"));
        assert!(r.has(flag::SUPPLEMENTARY) && r.has(flag::REVERSE));
    }
}
