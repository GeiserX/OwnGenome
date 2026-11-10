//! BAM reading for og's ports of the pysam tools: htslib, linked
//! statically, read the way pysam reads it.
//!
//! `fetch` is htslib's own region query (`sam_itr_queryi`), the call pysam
//! makes, so which reads a region returns is htslib's answer on both sides.
//! What pysam adds on top of a record, `reference_end`, `cigarstring` and
//! `get_blocks`, is rebuilt here from the raw CIGAR.

use std::fmt::Write as _;
use std::io::{self, BufRead, Write};
use std::path::Path;

use rust_htslib::bam::{self, FetchDefinition, Read as _};

mod args;
pub mod cyp2d6;
pub mod motif;

pub use rust_htslib::bam::Record;
pub use rust_htslib::errors::Error;

/// A BAM opened with its index, as `pysam.AlignmentFile(path)`.
pub struct Bam {
    reader: bam::IndexedReader,
}

impl Bam {
    pub fn open(path: &Path) -> Result<Self, Error> {
        Ok(Self {
            reader: bam::IndexedReader::from_path(path)?,
        })
    }

    /// `fetch(contig, start, stop)`: every read that overlaps the 0-based,
    /// half-open `[start, stop)`, in file order. A placed unmapped read
    /// covers one base, its position.
    pub fn fetch(
        &mut self,
        contig: &str,
        start: i64,
        stop: i64,
    ) -> Result<impl Iterator<Item = Result<Record, Error>> + '_, Error> {
        self.reader.fetch(FetchDefinition::RegionString(
            contig.as_bytes(),
            start,
            stop,
        ))?;
        Ok(self.reader.records())
    }

    /// The index of a contig of the header, as `bam.references` holds it.
    pub fn tid(&self, contig: &str) -> Option<u32> {
        self.reader.header().tid(contig.as_bytes())
    }

    /// `fetch(start, stop, tid=tid)`: the same query by contig index.
    pub fn fetch_tid(
        &mut self,
        tid: i32,
        start: i64,
        stop: i64,
    ) -> Result<impl Iterator<Item = Result<Record, Error>> + '_, Error> {
        self.reader
            .fetch(FetchDefinition::Region(tid, start, stop))?;
        Ok(self.reader.records())
    }
}

const CIGAR_OPS: &[u8; 10] = b"MIDNSHP=XB";

fn ops(r: &Record) -> impl Iterator<Item = (u8, i64)> + '_ {
    r.raw_cigar().iter().map(|&c| {
        (
            CIGAR_OPS.get((c & 0xf) as usize).copied().unwrap_or(b'?'),
            i64::from(c >> 4),
        )
    })
}

/// pysam's `reference_end`: one past the last reference base the alignment
/// covers; `None` for an unmapped read or one without a CIGAR. Like htslib's
/// `bam_endpos`, which pysam returns, a span of no reference bases counts as 1.
pub fn reference_end(r: &Record) -> Option<i64> {
    if r.is_unmapped() || r.raw_cigar().is_empty() {
        return None;
    }
    let len: i64 = ops(r)
        .filter(|(op, _)| matches!(op, b'M' | b'D' | b'N' | b'=' | b'X'))
        .map(|(_, l)| l)
        .sum();
    Some(r.pos() + len.max(1))
}

/// pysam's `cigarstring`: `None` when the read has no CIGAR.
pub fn cigar_string(r: &Record) -> Option<String> {
    if r.raw_cigar().is_empty() {
        return None;
    }
    let mut s = String::new();
    for (op, l) in ops(r) {
        let _ = write!(s, "{l}{}", char::from(op));
    }
    Some(s)
}

/// pysam's `get_blocks`: the reference stretches of the aligned bases
/// (M, = and X), 0-based and half-open; D and N move past a gap.
pub fn blocks(r: &Record) -> Vec<(i64, i64)> {
    let mut pos = r.pos();
    let mut out = Vec::new();
    for (op, l) in ops(r) {
        match op {
            b'M' | b'=' | b'X' => {
                out.push((pos, pos + l));
                pos += l;
            }
            b'D' | b'N' => pos += l,
            _ => {}
        }
    }
    out
}

fn none_or<T: std::fmt::Display>(v: Option<T>) -> String {
    v.map_or_else(|| "None".to_owned(), |v| v.to_string())
}

/// `og bam-fetch-probe BAM REGIONS`: for each `contig<TAB>start<TAB>stop`
/// line of REGIONS (`#` lines skipped), one line per fetched read in the
/// format of `tests/pysam_fetch.py`, so the two can be compared byte for byte.
pub fn fetch_probe(bam: &Path, regions: impl BufRead, mut out: impl Write) -> io::Result<()> {
    let mut bam = Bam::open(bam).map_err(io::Error::other)?;
    for line in regions.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let bad = || {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("not contig<TAB>start<TAB>stop: {line}"),
            )
        };
        let mut f = line.split('\t');
        let (Some(contig), Some(start), Some(stop), None) =
            (f.next(), f.next(), f.next(), f.next())
        else {
            return Err(bad());
        };
        let (start, stop): (i64, i64) = (
            start.parse().map_err(|_| bad())?,
            stop.parse().map_err(|_| bad())?,
        );
        for r in bam.fetch(contig, start, stop).map_err(io::Error::other)? {
            let r = r.map_err(io::Error::other)?;
            let blocks: Vec<String> = blocks(&r).iter().map(|(s, e)| format!("{s}-{e}")).collect();
            writeln!(
                out,
                "{contig}:{start}-{stop}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                String::from_utf8_lossy(r.qname()),
                r.flags(),
                r.pos(),
                none_or(reference_end(&r)),
                r.mapq(),
                none_or(cigar_string(&r)),
                blocks.join(","),
            )?;
        }
    }
    out.flush()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString};

    fn mapped(cigar: Vec<Cigar>, len: usize) -> Record {
        let mut r = Record::new();
        r.set(
            b"q",
            Some(&CigarString(cigar)),
            &vec![b'A'; len],
            &vec![30; len],
        );
        r.set_pos(100);
        r.unset_unmapped();
        r
    }

    #[test]
    fn reference_end_spans_the_reference_ops() {
        let r = mapped(vec![Cigar::SoftClip(2), Cigar::Match(5), Cigar::Del(3)], 7);
        assert_eq!(reference_end(&r), Some(108));
    }

    #[test]
    fn reference_end_of_a_zero_length_span_is_pos_plus_one() {
        // htslib's bam_endpos, which pysam returns, floors the span at 1.
        let r = mapped(vec![Cigar::Ins(10)], 10);
        assert_eq!(reference_end(&r), Some(101));
    }

    #[test]
    fn reference_end_of_an_unmapped_read_is_none() {
        let mut r = mapped(vec![Cigar::Match(5)], 5);
        r.set_unmapped();
        assert_eq!(reference_end(&r), None);
    }
}
