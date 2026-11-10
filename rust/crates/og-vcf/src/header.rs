use std::io::{self, BufRead, Write};

use crate::CsqLayout;

/// Lines of a VCF one at a time, without their line ending (`\n` or
/// `\r\n`), in one reused buffer.
pub struct Lines<R> {
    inner: R,
    buf: Vec<u8>,
    number: u64,
}

impl<R: BufRead> Lines<R> {
    pub fn new(inner: R) -> Self {
        Lines {
            inner,
            buf: Vec::with_capacity(1 << 12),
            number: 0,
        }
    }

    /// The next line, or `None` at the end.
    pub fn next_line(&mut self) -> io::Result<Option<&[u8]>> {
        self.buf.clear();
        if self.inner.read_until(b'\n', &mut self.buf)? == 0 {
            return Ok(None);
        }
        self.number += 1;
        let mut line = &self.buf[..];
        if let Some(l) = line.strip_suffix(b"\n") {
            line = l.strip_suffix(b"\r").unwrap_or(l);
        }
        Ok(Some(line))
    }

    /// The 1-based number of the line [`next_line`](Lines::next_line) last
    /// returned.
    pub fn number(&self) -> u64 {
        self.number
    }
}

/// The header of a VCF, line for line as it was read.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Header {
    /// Every line before `#CHROM`, without its line ending.
    pub meta: Vec<Vec<u8>>,
    /// The `#CHROM` line.
    pub columns: Vec<u8>,
}

impl Header {
    /// Reads up to and including the `#CHROM` line. A record before it, or
    /// no `#CHROM` line at all, is an error.
    pub fn read<R: BufRead>(lines: &mut Lines<R>) -> io::Result<Header> {
        let mut meta = Vec::new();
        while let Some(line) = lines.next_line()? {
            if line.starts_with(b"#CHROM") {
                return Ok(Header {
                    meta,
                    columns: line.to_vec(),
                });
            }
            if !line.starts_with(b"#") {
                return Err(invalid(format!(
                    "line {}: a record before the #CHROM header line",
                    lines.number()
                )));
            }
            meta.push(line.to_vec());
        }
        Err(invalid("no #CHROM header line"))
    }

    /// The sample names: the columns after FORMAT.
    pub fn samples(&self) -> impl Iterator<Item = &[u8]> {
        self.columns.split(|&c| c == b'\t').skip(9)
    }

    /// The `##contig` IDs in header order.
    pub fn contigs(&self) -> impl Iterator<Item = &[u8]> {
        self.meta.iter().filter_map(|l| {
            let id = l.strip_prefix(b"##contig=<ID=")?;
            let end = id
                .iter()
                .position(|&c| c == b',' || c == b'>')
                .unwrap_or(id.len());
            Some(&id[..end])
        })
    }

    /// The layout of the CSQ annotation, from the last CSQ header line.
    pub fn csq_layout(&self) -> Option<CsqLayout> {
        self.meta
            .iter()
            .rev()
            .find_map(|l| CsqLayout::from_header_line(l))
    }

    /// Writes the header back: the meta lines, then `extra` (lines a tool
    /// adds, without line endings), then `#CHROM`, each ending in `\n`.
    pub fn write(&self, out: &mut impl Write, extra: &[&[u8]]) -> io::Result<()> {
        for line in self
            .meta
            .iter()
            .map(Vec::as_slice)
            .chain(extra.iter().copied())
        {
            out.write_all(line)?;
            out.write_all(b"\n")?;
        }
        out.write_all(&self.columns)?;
        out.write_all(b"\n")
    }
}

fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    const VCF: &[u8] = b"##fileformat=VCFv4.2\r\n##contig=<ID=chr1,length=10>\n##contig=<ID=chrM>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\nchr1\t1\t.\tA\tC\t.\t.\t.\n";

    #[test]
    fn reads_and_writes_back() {
        let mut lines = Lines::new(VCF);
        let h = Header::read(&mut lines).unwrap();
        assert_eq!(h.meta.len(), 3);
        assert_eq!(h.samples().collect::<Vec<_>>(), [&b"s1"[..], b"s2"]);
        assert_eq!(h.contigs().collect::<Vec<_>>(), [&b"chr1"[..], b"chrM"]);
        assert_eq!(h.csq_layout(), None);
        assert_eq!(lines.number(), 4);
        assert_eq!(
            lines.next_line().unwrap(),
            Some(&b"chr1\t1\t.\tA\tC\t.\t.\t."[..])
        );
        assert_eq!(lines.number(), 5);
        assert_eq!(lines.next_line().unwrap(), None);
        let mut out = Vec::new();
        h.write(&mut out, &[b"##added=1"]).unwrap();
        let want = b"##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10>\n##contig=<ID=chrM>\n##added=1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n";
        assert_eq!(out, want);
    }

    #[test]
    fn refuses_a_file_without_a_header() {
        for vcf in [
            &b""[..],
            b"##fileformat=VCFv4.2\n",
            b"##x\nchr1\t1\t.\tA\tC\t.\t.\t.\n#CHROM\n",
        ] {
            assert!(Header::read(&mut Lines::new(vcf)).is_err());
        }
    }
}
