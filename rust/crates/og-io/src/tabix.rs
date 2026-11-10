//! Random access into a bgzipped, tabix-indexed text file: the lines of a
//! region, as `tabix FILE seq:beg-end` prints them.
//!
//! The index is a `.tbi` or a `.csi` beside the file. A query reads only
//! the BGZF blocks the index points at, so a lookup in a 100 MB VCF reads a
//! few kilobytes.
//!
//! Coordinates are 0-based and half-open: the region `chr1:100-200` of the
//! `tabix` command line (1-based, both ends included) is
//! `query("chr1", 99, 200)`. A line is returned when its own interval
//! overlaps the region, with the interval rules of htslib:
//!
//! - VCF: from `POS` to `POS + len(REF)`, or to `INFO/END` when that is
//!   after `POS`.
//! - Other files: the begin and end columns the index names; begin is
//!   1-based unless the index says 0-based (`tabix -p bed`, `-0`), and with
//!   no end column the line covers one base.
//!
//! SAM text (`tabix -p sam`) is refused: no OwnGenome file is indexed that way.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use flate2::{Crc, Decompress, FlushDecompress, Status};

/// What a line's columns mean, as the index records it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Conf {
    /// The kind of file: VCF or generic columns.
    pub preset: Preset,
    /// Begin is 0-based (BED); otherwise 1-based.
    pub zero_based: bool,
    /// 1-based column of the sequence name.
    pub col_seq: usize,
    /// 1-based column of the begin.
    pub col_beg: usize,
    /// 1-based column of the end; 0 when the file has none.
    pub col_end: usize,
    /// Lines starting with this byte are headers.
    pub meta: u8,
    /// Lines skipped at the top of the file.
    pub skip: u32,
}

/// The kinds of file a tabix index describes, of those this reader reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Preset {
    /// Begin and end from the columns of [`Conf`].
    Generic,
    /// `POS`, `REF` and `INFO/END`.
    Vcf,
}

/// One chunk of a bin: the virtual offsets of its first line and of the
/// byte after its last.
type Chunk = (u64, u64);

#[derive(Debug, Default)]
struct Bin {
    /// CSI only: no line of this bin's span starts before this offset.
    loff: u64,
    chunks: Vec<Chunk>,
}

#[derive(Debug, Default)]
struct Seq {
    bins: HashMap<u32, Bin>,
    /// TBI only: per window of `1 << min_shift` bases, the offset before
    /// which no line reaching that window starts.
    linear: Vec<u64>,
}

/// A parsed `.tbi` or `.csi` index.
#[derive(Debug)]
pub struct Index {
    conf: Conf,
    min_shift: u32,
    depth: u32,
    csi: bool,
    names: Vec<String>,
    ids: HashMap<String, usize>,
    seqs: Vec<Seq>,
}

fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

/// A little-endian cursor over the inflated index.
struct Bytes<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Bytes<'a> {
    fn take(&mut self, n: usize) -> io::Result<&'a [u8]> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&e| e <= self.data.len())
            .ok_or_else(|| invalid("the index ends too early"))?;
        let out = &self.data[self.pos..end];
        self.pos = end;
        Ok(out)
    }

    fn u32(&mut self) -> io::Result<u32> {
        Ok(u32::from_le_bytes(
            self.take(4)?.try_into().expect("4 bytes"),
        ))
    }

    fn u64(&mut self) -> io::Result<u64> {
        Ok(u64::from_le_bytes(
            self.take(8)?.try_into().expect("8 bytes"),
        ))
    }

    /// An `int32_t` count or field, which the format never lets be negative.
    fn count(&mut self, what: &str) -> io::Result<usize> {
        let v = self.u32()? as i32;
        usize::try_from(v).map_err(|_| invalid(format!("the index has a negative {what}")))
    }
}

/// The tabix header: inside a `.tbi` after the magic, inside a `.csi` as
/// its auxiliary data.
fn read_conf(b: &mut Bytes) -> io::Result<(Conf, Vec<String>)> {
    let format = b.u32()?;
    let preset = match format & 0xffff {
        0 => Preset::Generic,
        2 => Preset::Vcf,
        1 => {
            return Err(invalid(
                "the index is for SAM text, which this reader does not read",
            ));
        }
        p => return Err(invalid(format!("the index has an unknown preset {p}"))),
    };
    let conf = Conf {
        preset,
        zero_based: format & 0x10000 != 0,
        col_seq: b.count("column")?,
        col_beg: b.count("column")?,
        col_end: b.count("column")?,
        meta: b.u32()? as u8,
        skip: b.u32()?,
    };
    if conf.col_seq == 0 || conf.col_beg == 0 {
        return Err(invalid("the index names no sequence or begin column"));
    }
    let l_nm = b.count("name length")?;
    let names = b
        .take(l_nm)?
        .split(|&c| c == 0)
        .filter(|n| !n.is_empty())
        .map(|n| String::from_utf8_lossy(n).into_owned())
        .collect();
    Ok((conf, names))
}

impl Index {
    /// Reads a `.tbi` or `.csi` (by its magic, either compressed or not).
    pub fn read(path: impl AsRef<Path>) -> io::Result<Index> {
        let path = path.as_ref();
        let mut data = Vec::new();
        crate::open(path)?.read_to_end(&mut data)?;
        Index::parse(&data).map_err(|e| invalid(format!("{}: {e}", path.display())))
    }

    /// Parses an inflated `.tbi` or `.csi`.
    pub fn parse(data: &[u8]) -> io::Result<Index> {
        let mut b = Bytes { data, pos: 0 };
        let magic = b.take(4)?;
        let csi = match magic {
            b"TBI\x01" => false,
            b"CSI\x01" => true,
            _ => return Err(invalid("not a tabix (.tbi) or CSI (.csi) index")),
        };
        let (min_shift, depth, conf, names, n_ref);
        if csi {
            min_shift = b.u32()?;
            depth = b.u32()?;
            if !(1..=30).contains(&min_shift) || depth > 10 || min_shift + 3 * depth > 62 {
                return Err(invalid(format!(
                    "the index has an impossible min_shift {min_shift} and depth {depth}"
                )));
            }
            let l_aux = b.count("auxiliary length")?;
            if l_aux < 28 {
                return Err(invalid(
                    "the CSI index carries no tabix header, so its sequence names are unknown",
                ));
            }
            let mut aux = Bytes {
                data: b.take(l_aux)?,
                pos: 0,
            };
            (conf, names) = read_conf(&mut aux)?;
            n_ref = b.count("sequence count")?;
        } else {
            (min_shift, depth) = (14, 5);
            n_ref = b.count("sequence count")?;
            (conf, names) = read_conf(&mut b)?;
        }
        if names.len() != n_ref {
            return Err(invalid(format!(
                "the index names {} sequences but describes {n_ref}",
                names.len()
            )));
        }
        // The pseudo-bin holds counts, not lines.
        let pseudo = bin_first(depth + 1) + 1;
        let mut seqs = Vec::with_capacity(n_ref);
        for _ in 0..n_ref {
            let mut seq = Seq::default();
            for _ in 0..b.count("bin count")? {
                let id = b.u32()?;
                let loff = if csi { b.u64()? } else { 0 };
                let n = b.count("chunk count")?;
                let mut chunks = Vec::with_capacity(n.min(1 << 16));
                for _ in 0..n {
                    chunks.push((b.u64()?, b.u64()?));
                }
                if id != pseudo {
                    seq.bins.insert(id, Bin { loff, chunks });
                }
            }
            if !csi {
                let n = b.count("window count")?;
                seq.linear = (0..n).map(|_| b.u64()).collect::<io::Result<_>>()?;
            }
            seqs.push(seq);
        }
        let ids = names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.clone(), i))
            .collect();
        Ok(Index {
            conf,
            min_shift,
            depth,
            csi,
            names,
            ids,
            seqs,
        })
    }

    /// What the columns mean.
    pub fn conf(&self) -> &Conf {
        &self.conf
    }

    /// The sequences the file holds, in file order.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// One past the last position the index can address.
    fn max_pos(&self) -> u64 {
        1 << (self.min_shift + 3 * self.depth)
    }

    /// The chunks that may hold a line overlapping `[beg, end)` of sequence
    /// `id`, sorted and merged, as htslib's `hts_itr_query` picks them.
    fn chunks(&self, id: usize, beg: u64, end: u64) -> Vec<Chunk> {
        let seq = &self.seqs[id];
        let end = end.min(self.max_pos());
        if beg >= end {
            return Vec::new();
        }
        let min_off = if self.csi {
            // The deepest existing bin holding beg: its lines start no earlier.
            let mut bin = bin_first(self.depth) as u64 + (beg >> self.min_shift);
            let mut level = self.depth;
            loop {
                if let Some(b) = seq.bins.get(&(bin as u32)) {
                    break b.loff;
                }
                if level == 0 {
                    break 0;
                }
                bin = (bin - 1) >> 3;
                level -= 1;
            }
        } else {
            // A window with no offset (an old tabix wrote 0) takes the last
            // window before it that has one.
            let w = (beg >> self.min_shift) as usize;
            let w = w.min(seq.linear.len());
            seq.linear[..w]
                .iter()
                .chain(seq.linear.get(w))
                .rev()
                .copied()
                .find(|&o| o != 0)
                .unwrap_or(0)
        };
        let mut chunks: Vec<Chunk> = Vec::new();
        for level in 0..=self.depth {
            let shift = self.min_shift + 3 * (self.depth - level);
            let first = bin_first(level) as u64;
            for bin in first + (beg >> shift)..=first + ((end - 1) >> shift) {
                if let Some(b) = seq.bins.get(&(bin as u32)) {
                    chunks.extend(b.chunks.iter().filter(|c| c.1 > min_off));
                }
            }
        }
        chunks.sort_unstable();
        let mut merged: Vec<Chunk> = Vec::with_capacity(chunks.len());
        for c in chunks {
            match merged.last_mut() {
                Some(last) if c.0 <= last.1 => last.1 = last.1.max(c.1),
                _ => merged.push(c),
            }
        }
        merged
    }
}

/// The number of the first bin of `level` (0 is the root).
fn bin_first(level: u32) -> u32 {
    ((1u64 << (3 * level)) - 1) as u32 / 7
}

/// The `[beg, end)` of a data line, by the rules of htslib's `tbx_parse1`;
/// `None` for a line it cannot place.
pub fn interval<'a>(conf: &Conf, line: &'a [u8]) -> Option<(&'a [u8], u64, u64)> {
    let mut seq = None;
    let mut beg = None;
    let mut end = None;
    for (i, field) in line.split(|&c| c == b'\t').enumerate() {
        let col = i + 1;
        if col == conf.col_seq {
            seq = Some(field);
        }
        if col == conf.col_beg {
            let b = int(field)?;
            let b = if conf.zero_based { b } else { b - 1 };
            let b = b.max(0) as u64;
            beg = Some(b);
            end = Some(b + 1);
        } else if let Some(b) = beg {
            match conf.preset {
                Preset::Generic if col == conf.col_end => end = Some(int(field)?.max(0) as u64),
                Preset::Vcf if col == 4 && !field.is_empty() => end = Some(b + field.len() as u64),
                Preset::Vcf if col == 8 => {
                    if let Some(e) = info_end(field)
                        && e > b
                    {
                        end = Some(e);
                    }
                }
                _ => {}
            }
        }
    }
    let (beg, mut end) = (beg?, end?);
    if end <= beg {
        end = beg + 1;
    }
    Some((seq?, beg, end))
}

/// `INFO/END` as htslib finds it: the key at the start or after a `;`, its
/// leading digits read; `.` or no digits is no value.
fn info_end(info: &[u8]) -> Option<u64> {
    let at = if info.starts_with(b"END=") {
        4
    } else {
        info.windows(5).position(|w| w == b";END=")? + 5
    };
    let digits = info[at..].iter().take_while(|c| c.is_ascii_digit()).count();
    std::str::from_utf8(&info[at..at + digits])
        .ok()?
        .parse()
        .ok()
}

/// A whole number field: optional sign, then digits up to the end.
fn int(field: &[u8]) -> Option<i64> {
    std::str::from_utf8(field).ok()?.parse().ok()
}

/// A bgzipped file opened with its index, for region queries.
pub struct Reader {
    file: File,
    index: Index,
    /// The block in memory: its file offset, the offset of the next block,
    /// its bytes, and how many of them have been read.
    block: Option<(u64, u64)>,
    data: Vec<u8>,
    pos: usize,
    cdata: Vec<u8>,
}

impl Reader {
    /// Opens `path` with `path.tbi`, or `path.csi` when there is no `.tbi`.
    pub fn open(path: impl AsRef<Path>) -> io::Result<Reader> {
        let path = path.as_ref();
        let beside = |ext: &str| {
            let mut p = path.as_os_str().to_owned();
            p.push(ext);
            PathBuf::from(p)
        };
        let index = [beside(".tbi"), beside(".csi")]
            .into_iter()
            .find(|p| p.is_file())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("{} has no .tbi or .csi index beside it", path.display()),
                )
            })?;
        Reader::with_index(path, index)
    }

    /// Opens `path` with the index at `index`.
    pub fn with_index(path: impl AsRef<Path>, index: impl AsRef<Path>) -> io::Result<Reader> {
        Ok(Reader {
            file: File::open(path)?,
            index: Index::read(index)?,
            block: None,
            data: Vec::new(),
            pos: 0,
            cdata: Vec::new(),
        })
    }

    /// The index.
    pub fn index(&self) -> &Index {
        &self.index
    }

    /// The lines whose interval overlaps `[beg, end)` of `seq`, in file
    /// order, without their line ending. A sequence the index does not name
    /// has none.
    pub fn query(&mut self, seq: &str, beg: u64, end: u64) -> Query<'_> {
        let chunks = match self.index.ids.get(seq) {
            Some(&id) => self.index.chunks(id, beg, end),
            None => Vec::new(),
        };
        Query {
            reader: self,
            seq: seq.as_bytes().to_vec(),
            beg,
            end,
            chunks: chunks.into_iter(),
            chunk_end: 0,
            done: false,
        }
    }

    /// The virtual offset of the next unread byte, a block's end read as
    /// the start of the next block.
    fn tell(&self) -> u64 {
        match self.block {
            Some((_, next)) if self.pos == self.data.len() => next << 16,
            Some((at, _)) => (at << 16) | self.pos as u64,
            None => 0,
        }
    }

    fn seek(&mut self, voff: u64) -> io::Result<()> {
        let within = (voff & 0xffff) as usize;
        self.load(voff >> 16)?;
        if within > self.data.len() {
            return Err(invalid(format!(
                "the index points at byte {within} of a {}-byte block",
                self.data.len()
            )));
        }
        self.pos = within;
        Ok(())
    }

    /// Reads and inflates the BGZF block at file offset `at`; past the end
    /// of the file it is empty.
    fn load(&mut self, at: u64) -> io::Result<()> {
        if self.block.is_some_and(|(a, _)| a == at) {
            return Ok(());
        }
        self.block = None;
        self.data.clear();
        self.pos = 0;
        self.file.seek(SeekFrom::Start(at))?;
        let mut head = [0u8; 12];
        match self.file.read(&mut head[..1])? {
            0 => {
                self.block = Some((at, at));
                return Ok(());
            }
            _ => self.file.read_exact(&mut head[1..])?,
        }
        let bad = || invalid(format!("no BGZF block at byte {at}"));
        if head[..4] != [0x1f, 0x8b, 8, 4] {
            return Err(bad());
        }
        let xlen = u16::from_le_bytes([head[10], head[11]]) as usize;
        let mut extra = vec![0u8; xlen];
        self.file.read_exact(&mut extra)?;
        let mut bsize = None;
        let mut x = &extra[..];
        while x.len() >= 4 {
            let len = u16::from_le_bytes([x[2], x[3]]) as usize;
            if x[..2] == *b"BC" && len == 2 && x.len() >= 6 {
                bsize = Some(u16::from_le_bytes([x[4], x[5]]) as usize + 1);
            }
            x = x.get(4 + len..).unwrap_or_default();
        }
        let total = bsize.ok_or_else(bad)?;
        let rest = total.checked_sub(12 + xlen + 8).ok_or_else(bad)?;
        self.cdata.resize(rest + 8, 0);
        self.file.read_exact(&mut self.cdata)?;
        let (body, tail) = self.cdata.split_at(rest);
        let crc = u32::from_le_bytes(tail[..4].try_into().expect("4 bytes"));
        let isize = u32::from_le_bytes(tail[4..].try_into().expect("4 bytes")) as usize;
        if isize > 1 << 16 {
            return Err(bad());
        }
        self.data.resize(isize, 0);
        let mut d = Decompress::new(false);
        let status = d
            .decompress(body, &mut self.data, FlushDecompress::Finish)
            .map_err(|e| invalid(format!("BGZF block at byte {at}: {e}")))?;
        let mut sum = Crc::new();
        sum.update(&self.data);
        if status != Status::StreamEnd || d.total_out() as usize != isize || sum.sum() != crc {
            return Err(invalid(format!("the BGZF block at byte {at} is damaged")));
        }
        self.block = Some((at, at + total as u64));
        Ok(())
    }

    /// The next line into `out` (without `\n`); `false` at the end of the file.
    fn read_line(&mut self, out: &mut Vec<u8>) -> io::Result<bool> {
        out.clear();
        loop {
            let Some((at, next)) = self.block else {
                return Ok(false);
            };
            let rest = &self.data[self.pos..];
            if let Some(k) = rest.iter().position(|&c| c == b'\n') {
                out.extend_from_slice(&rest[..k]);
                self.pos += k + 1;
                return Ok(true);
            }
            out.extend_from_slice(rest);
            self.pos = self.data.len();
            if next == at {
                return Ok(!out.is_empty());
            }
            self.load(next)?;
        }
    }
}

/// The lines of one region: see [`Reader::query`].
pub struct Query<'a> {
    reader: &'a mut Reader,
    seq: Vec<u8>,
    beg: u64,
    end: u64,
    chunks: std::vec::IntoIter<Chunk>,
    chunk_end: u64,
    done: bool,
}

impl Query<'_> {
    fn step(&mut self) -> io::Result<Option<Vec<u8>>> {
        let mut line = Vec::new();
        loop {
            if self.reader.tell() >= self.chunk_end || self.reader.block.is_none() {
                let Some((start, end)) = self.chunks.next() else {
                    return Ok(None);
                };
                if self.reader.tell() != start || self.reader.block.is_none() {
                    self.reader.seek(start)?;
                }
                self.chunk_end = end;
                continue;
            }
            if !self.reader.read_line(&mut line)? {
                return Ok(None);
            }
            let conf = &self.reader.index.conf;
            if line.first() == Some(&conf.meta) {
                continue;
            }
            let (seq, beg, end) = interval(conf, &line).ok_or_else(|| {
                invalid(format!(
                    "cannot read the position of the line {:?}",
                    String::from_utf8_lossy(&line)
                ))
            })?;
            if seq != self.seq.as_slice() || beg >= self.end {
                return Ok(None);
            }
            if end >= self.beg {
                return Ok(Some(line));
            }
        }
    }
}

impl Iterator for Query<'_> {
    type Item = io::Result<Vec<u8>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let r = self.step().transpose();
        if !matches!(r, Some(Ok(_))) {
            self.done = true;
        }
        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bgzf::Writer;
    use std::collections::BTreeMap;
    use std::io::Write;
    use std::num::NonZero;

    const BLOCK_INPUT: usize = 0xff00;

    fn vcf_conf() -> Conf {
        Conf {
            preset: Preset::Vcf,
            zero_based: false,
            col_seq: 1,
            col_beg: 2,
            col_end: 0,
            meta: b'#',
            skip: 0,
        }
    }

    fn bed_conf() -> Conf {
        Conf {
            preset: Preset::Generic,
            zero_based: true,
            col_seq: 1,
            col_beg: 2,
            col_end: 3,
            meta: b'#',
            skip: 0,
        }
    }

    fn temp(name: &str, data: &[u8]) -> PathBuf {
        let p = std::env::temp_dir().join(format!("og-tabix-{}-{name}", std::process::id()));
        std::fs::write(&p, data).unwrap();
        p
    }

    /// xorshift64*: the same numbers on every run.
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    /// The text compressed by og-io's writer, and the file offset of each
    /// block (each holds `BLOCK_INPUT` bytes but the last).
    fn bgzip(text: &[u8]) -> (Vec<u8>, Vec<u64>) {
        let mut w = Writer::new(Vec::new(), NonZero::<usize>::MIN).unwrap();
        w.write_all(text).unwrap();
        let data = w.finish().unwrap();
        let mut starts = Vec::new();
        let mut at = 0;
        while at < data.len() {
            starts.push(at as u64);
            at += u16::from_le_bytes([data[at + 16], data[at + 17]]) as usize + 1;
        }
        (data, starts)
    }

    /// The virtual offset of byte `u` of the text. `end` writes a block's
    /// end as htslib's `bgzf_tell` does after a line that fills the block:
    /// `(block << 16) | 65280`, not the start of the next block.
    fn voff(starts: &[u64], u: usize, end: bool) -> u64 {
        let (k, w) = (u / BLOCK_INPUT, u % BLOCK_INPUT);
        if end && w == 0 && k > 0 {
            (starts[k - 1] << 16) | BLOCK_INPUT as u64
        } else {
            (starts[k] << 16) | w as u64
        }
    }

    /// One sequence of a built index: chunks by bin, and the linear index.
    type Built = (BTreeMap<u32, Vec<Chunk>>, Vec<u64>);

    fn reg2bin(beg: u64, end: u64, min_shift: u32, depth: u32) -> u32 {
        let end = end - 1;
        let mut s = min_shift;
        for level in (1..=depth).rev() {
            if beg >> s == end >> s {
                return bin_first(level) + (beg >> s) as u32;
            }
            s += 3;
        }
        0
    }

    /// An index of `text` built the way htslib's `hts_idx_push` builds one:
    /// TBI when `csi` is `None`, else CSI with that min_shift and depth.
    fn build_index(text: &[u8], starts: &[u64], conf: &Conf, csi: Option<(u32, u32)>) -> Vec<u8> {
        let (min_shift, depth) = csi.unwrap_or((14, 5));
        let mut names: Vec<Vec<u8>> = Vec::new();
        let mut seqs: Vec<Built> = Vec::new();
        let mut u = 0;
        for line in text.split_inclusive(|&c| c == b'\n') {
            let (from, to) = (u, u + line.len());
            u = to;
            let body = line.strip_suffix(b"\n").unwrap_or(line);
            if body.first() == Some(&conf.meta) {
                continue;
            }
            let (seq, beg, end) = interval(conf, body).unwrap();
            if names.last().map(Vec::as_slice) != Some(seq) {
                names.push(seq.to_vec());
                seqs.push(Default::default());
            }
            let (bins, linear) = seqs.last_mut().unwrap();
            let (a, b) = (voff(starts, from, false), voff(starts, to, true));
            let chunks = bins.entry(reg2bin(beg, end, min_shift, depth)).or_default();
            match chunks.last_mut() {
                Some(last) if last.1 == a => last.1 = b,
                _ => chunks.push((a, b)),
            }
            for w in (beg >> min_shift)..=((end - 1) >> min_shift) {
                let w = w as usize;
                if linear.len() <= w {
                    linear.resize(w + 1, 0);
                }
                if linear[w] == 0 {
                    linear[w] = a;
                }
            }
        }
        let mut head = Vec::new();
        let format = match conf.preset {
            Preset::Generic => 0,
            Preset::Vcf => 2,
        } | if conf.zero_based { 0x10000 } else { 0 };
        for v in [
            format,
            conf.col_seq as u32,
            conf.col_beg as u32,
            conf.col_end as u32,
            conf.meta as u32,
            conf.skip,
        ] {
            head.extend(v.to_le_bytes());
        }
        let nm: Vec<u8> = names
            .iter()
            .flat_map(|n| n.iter().copied().chain([0]))
            .collect();
        head.extend((nm.len() as u32).to_le_bytes());
        head.extend(&nm);
        let mut out = Vec::new();
        match csi {
            None => {
                out.extend(b"TBI\x01");
                out.extend((names.len() as u32).to_le_bytes());
                out.extend(&head);
            }
            Some(_) => {
                out.extend(b"CSI\x01");
                out.extend(min_shift.to_le_bytes());
                out.extend(depth.to_le_bytes());
                out.extend((head.len() as u32).to_le_bytes());
                out.extend(&head);
                out.extend((names.len() as u32).to_le_bytes());
            }
        }
        for (bins, mut linear) in seqs {
            // Empty windows take the offset before them, as hts_idx_finish does.
            for w in 1..linear.len() {
                if linear[w] == 0 {
                    linear[w] = linear[w - 1];
                }
            }
            out.extend((bins.len() as u32).to_le_bytes());
            for (bin, chunks) in bins {
                out.extend(bin.to_le_bytes());
                if csi.is_some() {
                    // The offset of the bin's first window (htslib's update_loff).
                    let mut level = 0;
                    while bin_first(level + 1) <= bin {
                        level += 1;
                    }
                    let first = ((bin - bin_first(level)) as u64)
                        << (min_shift + 3 * (depth - level))
                        >> min_shift;
                    let loff = linear
                        .get(first as usize)
                        .or(linear.last())
                        .copied()
                        .unwrap_or(0);
                    out.extend(loff.to_le_bytes());
                }
                out.extend((chunks.len() as u32).to_le_bytes());
                for (a, b) in chunks {
                    out.extend(a.to_le_bytes());
                    out.extend(b.to_le_bytes());
                }
            }
            if csi.is_none() {
                // Windows reached only through the backfill: keep them 0, as
                // a tabix before 0.1.4 left them, so the backward walk runs.
                out.extend((linear.len() as u32).to_le_bytes());
                for o in linear {
                    out.extend(o.to_le_bytes());
                }
            }
        }
        out
    }

    /// Every data line with its interval, parsed once.
    type Parsed = Vec<(Vec<u8>, (Vec<u8>, u64, u64))>;

    fn parse_all(text: &[u8], conf: &Conf) -> Parsed {
        text.split(|&c| c == b'\n')
            .filter(|l| !l.is_empty() && l.first() != Some(&conf.meta))
            .map(|l| {
                let (s, b, e) = interval(conf, l).unwrap();
                (l.to_vec(), (s.to_vec(), b, e))
            })
            .collect()
    }

    /// What a query must return: every data line whose interval overlaps.
    fn scan(all: &Parsed, seq: &str, beg: u64, end: u64) -> Vec<Vec<u8>> {
        all.iter()
            .filter(|(_, (s, b, e))| s == seq.as_bytes() && *b < end && *e > beg)
            .map(|(l, _)| l.clone())
            .collect()
    }

    fn query(r: &mut Reader, seq: &str, beg: u64, end: u64) -> Vec<Vec<u8>> {
        r.query(seq, beg, end).collect::<io::Result<_>>().unwrap()
    }

    /// A VCF of three sequences over many 16 kb windows and six BGZF
    /// blocks: long REFs, INFO/END far ahead, an END before POS, and keys
    /// that only end in `END=`.
    fn vcf(seed: u64) -> Vec<u8> {
        let mut rng = Rng(seed);
        let mut t =
            b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_vec();
        for seq in ["chr1", "chr2", "chrX"] {
            let mut pos = 1 + rng.below(50);
            for i in 0..1500 {
                let reflen = match rng.below(400) {
                    0 => 20_000,
                    1..=50 => 1 + rng.below(30),
                    _ => 1,
                } as usize;
                let info = match rng.below(25) {
                    0 => format!("SVTYPE=DEL;END={}", pos + 1 + rng.below(70_000)),
                    1 => format!("END={}", pos + rng.below(3)),
                    2 => format!("SVEND={};CIEND=-5,5", pos + 90_000),
                    3 => format!("END={}", pos.saturating_sub(5)),
                    _ => ".".to_owned(),
                };
                let r = "A".repeat(reflen);
                t.extend(format!("{seq}\t{pos}\tid{i}\t{r}\tC\t50\tPASS\t{info}\n").as_bytes());
                pos += rng.below(400);
            }
        }
        t
    }

    /// Regions at every edge that matters: each line's own begin and end,
    /// one base either side, and the 16 kb window edges.
    fn regions(text: &[u8], conf: &Conf) -> Vec<(String, u64, u64)> {
        let mut out = Vec::new();
        for (i, l) in text.split(|&c| c == b'\n').enumerate() {
            if l.is_empty() || l[0] == conf.meta || i % 11 != 0 {
                continue;
            }
            let (s, b, e) = interval(conf, l).unwrap();
            let s = String::from_utf8(s.to_vec()).unwrap();
            for (x, y) in [
                (b, b + 1),
                (b.saturating_sub(1), b),
                (e - 1, e),
                (e, e + 1),
                (b, e),
                (b.saturating_sub(300), e + 300),
            ] {
                out.push((s.clone(), x, y));
            }
            let w = (b >> 14) << 14;
            out.push((s.clone(), w.saturating_sub(1), w));
            out.push((s.clone(), w, w + 1));
            out.push((s.clone(), w + 16383, w + 16385));
        }
        for s in ["chr1", "chr2", "chrX"] {
            out.push((s.to_owned(), 0, 1 << 29));
            out.push((s.to_owned(), 5_000_000, 6_000_000));
        }
        out.push(("chr3".to_owned(), 0, 1 << 29));
        out
    }

    fn check(text: &[u8], conf: &Conf, csi: Option<(u32, u32)>, name: &str) {
        let (data, starts) = bgzip(text);
        assert!(starts.len() >= 4, "only {} blocks", starts.len());
        let file = temp(&format!("{name}.gz"), &data);
        let idx = temp(
            &format!("{name}.idx"),
            &build_index(text, &starts, conf, csi),
        );
        let mut r = Reader::with_index(&file, &idx).unwrap();
        assert_eq!(r.index().conf(), conf);
        let all = parse_all(text, conf);
        let mut hits = 0;
        for (s, b, e) in regions(text, conf) {
            let want = scan(&all, &s, b, e);
            hits += want.len();
            assert_eq!(query(&mut r, &s, b, e), want, "{name} {s}:[{b},{e})");
        }
        assert!(hits > 10_000, "{name}: the regions found only {hits} lines");
        let _ = std::fs::remove_file(file);
        let _ = std::fs::remove_file(idx);
    }

    #[test]
    fn vcf_queries_equal_a_full_scan_tbi_and_csi() {
        let text = vcf(7);
        check(&text, &vcf_conf(), None, "vcf-tbi");
        check(&text, &vcf_conf(), Some((14, 5)), "vcf-csi");
        check(&text, &vcf_conf(), Some((10, 7)), "vcf-csi-deep");
    }

    #[test]
    fn bed_queries_equal_a_full_scan() {
        let mut rng = Rng(11);
        let mut text = b"#chrom\tstart\tend\n".to_vec();
        for seq in ["chr1", "chr2"] {
            let mut beg = rng.below(100);
            for _ in 0..4000 {
                let len = if rng.below(30) == 0 {
                    40_000
                } else {
                    rng.below(200)
                };
                text.extend(format!("{seq}\t{beg}\t{}\tx\n", beg + len).as_bytes());
                beg += rng.below(120);
            }
        }
        check(&text, &bed_conf(), None, "bed-tbi");
        check(&text, &bed_conf(), Some((14, 5)), "bed-csi");
    }

    #[test]
    fn reads_the_stub_vcf_tabix_wrote() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../../assets/stub/sample.vcf.gz"
        );
        let mut r = Reader::open(path).unwrap();
        assert_eq!(r.index().names(), ["chr1"]);
        assert_eq!(r.index().conf(), &vcf_conf());
        let rec = b"chr1\t5\t.\tN\tA\t50\tPASS\t.\tGT\t0/1".to_vec();
        // POS 5 is the half-open [4, 5).
        assert_eq!(query(&mut r, "chr1", 4, 5), std::slice::from_ref(&rec));
        assert_eq!(query(&mut r, "chr1", 0, 1 << 29), [rec]);
        assert!(query(&mut r, "chr1", 5, 6).is_empty());
        assert!(query(&mut r, "chr1", 3, 4).is_empty());
        assert!(query(&mut r, "chr2", 0, 100).is_empty());
    }

    #[test]
    fn intervals_follow_htslib() {
        let v = vcf_conf();
        let iv = |c: &Conf, l: &str| interval(c, l.as_bytes()).map(|(s, b, e)| (s.to_vec(), b, e));
        let s = b"c".to_vec();
        assert_eq!(
            iv(&v, "c\t10\t.\tACG\tA\t.\t.\t."),
            Some((s.clone(), 9, 12))
        );
        assert_eq!(
            iv(&v, "c\t10\t.\tA\t<DEL>\t.\t.\tEND=500"),
            Some((s.clone(), 9, 500))
        );
        assert_eq!(
            iv(&v, "c\t10\t.\tA\t<DEL>\t.\t.\tX=1;END=20;Y"),
            Some((s.clone(), 9, 20))
        );
        // CIEND= and SVEND= are not END=; an END at or before POS is ignored.
        assert_eq!(
            iv(&v, "c\t10\t.\tA\tT\t.\t.\tCIEND=5;SVEND=900"),
            Some((s.clone(), 9, 10))
        );
        assert_eq!(
            iv(&v, "c\t10\t.\tAA\tT\t.\t.\tEND=9"),
            Some((s.clone(), 9, 11))
        );
        assert_eq!(
            iv(&v, "c\t10\t.\tA\tT\t.\t.\tEND=."),
            Some((s.clone(), 9, 10))
        );
        assert_eq!(iv(&v, "c\tx\t.\tA\tT"), None);
        let b = bed_conf();
        assert_eq!(iv(&b, "c\t0\t10"), Some((s.clone(), 0, 10)));
        // A row that covers no base is placed on one, as tabix places it.
        assert_eq!(iv(&b, "c\t7\t7"), Some((s.clone(), 7, 8)));
        let one = Conf {
            zero_based: false,
            col_end: 0,
            ..bed_conf()
        };
        assert_eq!(iv(&one, "c\t1\tq"), Some((s, 0, 1)));
    }

    #[test]
    fn refuses_what_it_cannot_read() {
        let err = |d: &[u8]| Index::parse(d).unwrap_err().to_string();
        assert!(err(b"BAI\x01").contains("not a tabix"));
        assert!(err(b"TBI\x01\x01\x00").contains("ends too early"));
        let mut sam = b"TBI\x01".to_vec();
        for v in [0u32, 1, 3, 4, 0, b'@' as u32, 0, 0] {
            sam.extend(v.to_le_bytes());
        }
        assert!(err(&sam).contains("SAM"));
        let mut csi = b"CSI\x01".to_vec();
        for v in [14u32, 5, 0, 0] {
            csi.extend(v.to_le_bytes());
        }
        assert!(err(&csi).contains("no tabix header"));
        let missing = Reader::open("/nonexistent/x.vcf.gz").err().unwrap();
        assert!(missing.to_string().contains("no .tbi or .csi"));
    }

    #[test]
    fn a_damaged_block_is_an_error() {
        let text = vcf(3);
        let (mut data, starts) = bgzip(&text);
        let idx = build_index(&text, &starts, &vcf_conf(), None);
        let at = starts[1] as usize + 40;
        data[at] ^= 0xff;
        let file = temp("damaged.gz", &data);
        let idx = temp("damaged.tbi", &idx);
        let mut r = Reader::with_index(&file, &idx).unwrap();
        let all: Vec<_> = r.query("chr1", 0, 1 << 29).collect();
        assert!(all.last().unwrap().is_err(), "a damaged block read cleanly");
        let _ = std::fs::remove_file(file);
        let _ = std::fs::remove_file(idx);
    }
}
