//! `og vcf index`: the `.tbi` index of a BGZF VCF, as `tabix -p vcf` and
//! `pysam.tabix_index(preset="vcf")` write it.
//!
//! The index is built the way htslib 1.24 builds one (`tbx_index`,
//! `hts_idx_push`, `hts_idx_finish`; the htslib of the pipeline's pysam
//! and bcftools), so its bytes, once inflated, are htslib's: the same bins
//! in the same order (the order of htslib's hash table), the same chunks
//! and the same linear index. Only the deflate stream of the `.tbi` can
//! differ, as two deflate libraries differ; every reader inflates both to
//! the same bytes.
//!
//! What htslib does, and so this does too:
//!
//! - A line that starts with `#` is a header.
//! - A record starts at `POS - 1` and ends at the largest of `INFO/END`
//!   (when after the start), the REF length, the `|SVLEN|` of a `<DEL>`,
//!   `<DUP>`, `<INV>` or `<CNV>` ALT, and a sample's `LEN` when an ALT is
//!   `<*>` or `<NON_REF>`. Numbers are read as C reads them (`strtoll`,
//!   `atoll`); a POS of 0 or below is placed at the first base.
//! - Refused, and no index written: a file that is not BGZF, a line with
//!   no position to read (fewer than two columns, or a POS that is not a
//!   number), a sequence that starts again after another one, positions
//!   that go back on one sequence, and a record past 2^29 (a `.tbi` cannot
//!   hold it).
//!
//! [`build_positions`] indexes a tab-separated file with one base per
//! line, as `tabix -s1 -b2 -e2` does: the sequence in column 1, the 1-based
//! position in column 2 (read as C's `strtoll` with base 0), each line over
//! that one base. A POS of 0 or below is placed at the first base, lines
//! that start with `#` are headers, and the same files are refused.
//!
//! It reads the file once, one BGZF block at a time: memory is one block,
//! one line and the index itself.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, BufReader, Read, Write};
use std::mem;
use std::num::NonZero;
use std::path::{Path, PathBuf};

use flate2::{Crc, Decompress, FlushDecompress, Status};

use crate::AtomicFile;
use crate::bgzf;

/// The `.tbi` scheme: 16 kb windows, five levels of bins below the root.
const MIN_SHIFT: u32 = 14;
const N_LVLS: u32 = 5;
/// The bins of the scheme: the root, then 8, 64, ... 32768 of them.
const N_BINS: u32 = ((1 << (3 * N_LVLS + 3)) - 1) / 7;
/// The pseudo-bin holding each sequence's offsets and record count.
const META_BIN: u32 = N_BINS + 1;
/// The last position a `.tbi` can hold.
const MAX_POS: i64 = 1 << (MIN_SHIFT + 3 * N_LVLS);
/// htslib merges a bin into its parent when its chunks span fewer
/// compressed bytes than this.
const MIN_MARKER_DIST: u64 = 0x10000;
/// A window of the linear index that no record has reached yet.
const UNSET: u64 = u64::MAX;
/// The tabix header of the `vcf` preset: VCF, sequence in column 1, begin
/// in column 2, no end column, `#` for headers, no lines skipped.
const VCF_CONF: [i32; 6] = [2, 1, 2, 0, b'#' as i32, 0];
/// The tabix header of `tabix -s1 -b2 -e2`: a generic file, sequence in
/// column 1, begin and end in column 2, `#` for headers, no lines skipped.
const POSITIONS_CONF: [i32; 6] = [0, 1, 2, 2, b'#' as i32, 0];

/// A chunk: the virtual offsets of its first byte and of the byte after it.
type Chunk = (u64, u64);

fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

/// Reads into `buf` until it is full or the input ends; the bytes read.
fn read_full_or_end(r: &mut impl Read, buf: &mut [u8]) -> io::Result<usize> {
    let mut n = 0;
    while n < buf.len() {
        match r.read(&mut buf[n..]) {
            Ok(0) => break,
            Ok(k) => n += k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(n)
}

/// The BGZF blocks of a file, inflated one at a time.
struct Blocks<R> {
    inner: R,
    /// The file offset after the last block read.
    pos: u64,
    cdata: Vec<u8>,
    inflater: Decompress,
}

impl<R: Read> Blocks<R> {
    /// Inflates blocks into `data` until one holds bytes, as htslib skips
    /// empty blocks; that block's file offset, or `None` at the end.
    fn next_data(&mut self, data: &mut Vec<u8>) -> io::Result<Option<u64>> {
        loop {
            let at = self.pos;
            let mut head = [0u8; 12];
            let n = read_full_or_end(&mut self.inner, &mut head)?;
            let not_bgzf = || {
                invalid(if at == 0 {
                    "not BGZF-compressed (compress it with bgzip)".to_owned()
                } else {
                    format!("no BGZF block at byte {at}")
                })
            };
            if n == 0 && at > 0 {
                return Ok(None);
            }
            // An empty file is not BGZF either: even an empty VCF.gz holds
            // the empty end-of-file block.
            if n < head.len() || head[..4] != [0x1f, 0x8b, 8, 4] {
                return Err(not_bgzf());
            }
            let xlen = u16::from_le_bytes([head[10], head[11]]) as usize;
            let mut extra = vec![0u8; xlen];
            self.inner.read_exact(&mut extra)?;
            let mut total = None;
            let mut x = &extra[..];
            while x.len() >= 4 {
                let len = u16::from_le_bytes([x[2], x[3]]) as usize;
                if x[..2] == *b"BC" && len == 2 && x.len() >= 6 {
                    total = Some(u16::from_le_bytes([x[4], x[5]]) as usize + 1);
                }
                x = x.get(4 + len..).unwrap_or_default();
            }
            let total = total.ok_or_else(not_bgzf)?;
            let body = total.checked_sub(12 + xlen + 8).ok_or_else(not_bgzf)?;
            self.cdata.resize(body + 8, 0);
            self.inner.read_exact(&mut self.cdata)?;
            let (deflated, tail) = self.cdata.split_at(body);
            let crc = u32::from_le_bytes(tail[..4].try_into().expect("4 bytes"));
            let isize = u32::from_le_bytes(tail[4..].try_into().expect("4 bytes")) as usize;
            let damaged = || invalid(format!("the BGZF block at byte {at} is damaged"));
            if isize > 1 << 16 {
                return Err(damaged());
            }
            data.clear();
            data.resize(isize, 0);
            self.inflater.reset(false);
            let status = self
                .inflater
                .decompress(deflated, data, FlushDecompress::Finish)
                .map_err(|_| damaged())?;
            let mut sum = Crc::new();
            sum.update(data);
            if status != Status::StreamEnd
                || self.inflater.total_out() as usize != isize
                || sum.sum() != crc
            {
                return Err(damaged());
            }
            self.pos = at + total as u64;
            if !data.is_empty() {
                return Ok(Some(at));
            }
        }
    }
}

/// Lines of a BGZF file with htslib's virtual offsets: `bgzf_getline` and
/// `bgzf_tell`.
struct Lines<R> {
    blocks: Blocks<R>,
    data: Vec<u8>,
    off: usize,
    /// htslib's `block_address`: the block being read, or, once a block is
    /// used up, the file offset after it. The end of the file leaves it
    /// where it was.
    addr: u64,
}

impl<R: Read> Lines<R> {
    fn new(inner: R) -> Self {
        Lines {
            blocks: Blocks {
                inner,
                pos: 0,
                cdata: Vec::new(),
                inflater: Decompress::new(false),
            },
            data: Vec::new(),
            off: 0,
            addr: 0,
        }
    }

    /// The virtual offset of the next unread byte.
    fn tell(&self) -> u64 {
        (self.addr << 16) | self.off as u64
    }

    /// The next line into `line`, without its `\n` and a `\r` before it;
    /// `false` at the end.
    fn read_line(&mut self, line: &mut Vec<u8>) -> io::Result<bool> {
        line.clear();
        let mut ended = false;
        while !ended {
            if self.off >= self.data.len() {
                match self.blocks.next_data(&mut self.data)? {
                    Some(at) => {
                        self.addr = at;
                        self.off = 0;
                    }
                    None => {
                        self.data.clear();
                        self.off = 0;
                        break;
                    }
                }
            }
            let rest = &self.data[self.off..];
            match rest.iter().position(|&c| c == b'\n') {
                Some(k) => {
                    line.extend_from_slice(&rest[..k]);
                    self.off += k + 1;
                    ended = true;
                }
                None => {
                    line.extend_from_slice(rest);
                    self.off = self.data.len();
                }
            }
            if self.off >= self.data.len() {
                self.addr = self.blocks.pos;
                self.data.clear();
                self.off = 0;
            }
        }
        if !ended && line.is_empty() {
            return Ok(false);
        }
        if line.last() == Some(&b'\r') {
            line.pop();
        }
        Ok(true)
    }
}

/// C's `strtoll(s, &end, base)` for base 0 (decimal, octal after a `0`,
/// hex after `0x`) or 10: leading space, a sign, then digits up to the
/// first other byte, clamped on overflow. `None` when no digit was read.
fn strtoll(s: &[u8], auto: bool) -> Option<i64> {
    let mut i = s
        .iter()
        .position(|c| !matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r'))
        .unwrap_or(s.len());
    let neg = match s.get(i) {
        Some(b'-') => {
            i += 1;
            true
        }
        Some(b'+') => {
            i += 1;
            false
        }
        _ => false,
    };
    let hex = auto
        && s.get(i) == Some(&b'0')
        && matches!(s.get(i + 1), Some(b'x' | b'X'))
        && s.get(i + 2).is_some_and(u8::is_ascii_hexdigit);
    let (base, from) = match (hex, auto && s.get(i) == Some(&b'0')) {
        (true, _) => (16, i + 2),
        (false, true) => (8, i),
        _ => (10, i),
    };
    let mut value: u64 = 0;
    let mut any = false;
    for &c in &s[from..] {
        let Some(d) = char::from(c).to_digit(base) else {
            break;
        };
        any = true;
        value = value
            .saturating_mul(u64::from(base))
            .saturating_add(u64::from(d));
    }
    if !any {
        return None;
    }
    Some(match neg {
        false => i64::try_from(value).unwrap_or(i64::MAX),
        true if value > 1 << 63 => i64::MIN,
        true => 0i64.wrapping_sub_unsigned(value),
    })
}

/// C's `atoll`: base 10, 0 when there is no digit.
fn atoll(s: &[u8]) -> i64 {
    strtoll(s, false).unwrap_or(0)
}

/// Where the value of INFO key `key` (with its `=`) starts, as htslib
/// finds it: at the start of the field, or after a `;` when the key's
/// first appearance is not at the start.
fn info_value<'a>(info: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    let n = key.len();
    if info.starts_with(key) {
        return Some(&info[n..]);
    }
    info.windows(n).position(|w| w == key)?;
    let at = info
        .windows(n + 1)
        .position(|w| w[0] == b';' && &w[1..] == key)?;
    Some(&info[at + n + 1..])
}

/// A symbolic ALT whose length is SVLEN (htslib's
/// `svlen_on_ref_for_vcf_alt`): `<DEL>`, `<DUP>`, `<INV>`, `<CNV>`, or
/// one of them with a `:` subtype.
fn svlen_on_ref(alt: &[u8]) -> bool {
    alt.len() >= 5
        && alt[0] == b'<'
        && matches!(alt[4], b'>' | b':')
        && matches!(&alt[..4], b"<CNV" | b"<DEL" | b"<DUP" | b"<INV")
        && alt[alt.len() - 1] == b'>'
}

/// A VCF line's sequence name and 0-based `[beg, end)`, by the rules of
/// htslib 1.24's `tbx_parse1` (the htslib of the pipeline's pysam and
/// bcftools); `None` for a line it cannot place. A NUL byte splits columns
/// as a tab does, as in htslib. The flag says the POS was 0 or below and
/// moved to the first base. `sv` is scratch space.
///
/// The end is the largest of: `INFO/END` when after POS, the REF length,
/// the largest `|SVLEN|` of a `<DEL>`, `<DUP>`, `<INV>` or `<CNV>` ALT, and
/// the largest sample `LEN` when an ALT is `<*>` or `<NON_REF>`.
fn parse<'a>(line: &'a [u8], sv: &mut Vec<bool>) -> Option<(&'a [u8], i64, i64, bool)> {
    let (mut name, mut beg, mut end, mut clamped) = (None, -1i64, -1i64, false);
    let (mut reflen, mut svlen, mut fmtlen) = (0i64, 0i64, 0i64);
    // Per allele, REF first: an ALT whose SVLEN counts. The buffer is the
    // caller's, so a line costs no allocation.
    sv.clear();
    let (mut alleles, mut use_svlen, mut getlen) = (0usize, false, false);
    let mut lenpos = None;
    let len = |f: &[u8]| i64::try_from(f.len()).unwrap_or(i64::MAX);
    for (i, field) in line.split(|&c| c == b'\t' || c == 0).enumerate() {
        match i + 1 {
            1 => name = Some(field),
            2 => {
                beg = strtoll(field, true)?.saturating_sub(1);
                if beg < 0 {
                    beg = 0;
                    clamped = true;
                }
                end = end.max(1);
            }
            4 => {
                if !field.is_empty() {
                    end = beg.saturating_add(len(field));
                }
                alleles += 1;
                sv.push(false);
                reflen = len(field);
            }
            5 => {
                for alt in field.split(|&c| c == b',') {
                    if alleles >= 65536 {
                        break;
                    }
                    alleles += 1;
                    let on_ref = svlen_on_ref(alt);
                    sv.push(on_ref);
                    use_svlen |= on_ref;
                    getlen |= !on_ref && (alt == b"<*>" || alt == b"<NON_REF>");
                }
            }
            8 => {
                if let Some(v) = info_value(field, b"END=")
                    && v.first() != Some(&b'.')
                {
                    let e = strtoll(v, true).unwrap_or(0);
                    if e > beg {
                        end = e;
                    }
                }
                // SVLEN values pair with ALTs 1, 2, ...; the next value
                // starts after the next comma anywhere in the field.
                let mut rest = info_value(field, b"SVLEN=");
                let mut d = 1;
                while let Some(v) = rest
                    && d < alleles
                {
                    let tmp = if use_svlen && sv.get(d).copied().unwrap_or(false) {
                        atoll(v).saturating_abs()
                    } else {
                        1
                    };
                    svlen = svlen.max(tmp);
                    rest = v.iter().position(|&c| c == b',').map(|k| &v[k + 1..]);
                    d += 1;
                }
            }
            9 if getlen => match field.split(|&c| c == b':').position(|k| k == b"LEN") {
                Some(p) => lenpos = Some(p),
                None => break,
            },
            n if n > 9 && getlen => {
                if let Some(p) = lenpos {
                    let tmp = field.split(|&c| c == b':').nth(p).map_or(0, atoll);
                    fmtlen = fmtlen.max(tmp);
                }
            }
            _ => {}
        }
    }
    if beg < 0 {
        return None;
    }
    end = end.max(beg.saturating_add(reflen.max(svlen).max(fmtlen)));
    Some((name?, beg, end, clamped))
}

/// The bin of `[beg, end)`: the smallest one that holds it (htslib's
/// `hts_reg2bin`).
fn reg2bin(beg: i64, end: i64) -> u32 {
    let end = end - 1;
    let mut shift = MIN_SHIFT;
    let mut first = bin_first(N_LVLS);
    for level in (1..=N_LVLS).rev() {
        if beg >> shift == end >> shift {
            return first + (beg >> shift) as u32;
        }
        shift += 3;
        first -= 1 << (3 * (level - 1));
    }
    0
}

/// The number of the first bin of `level` (0 is the root).
fn bin_first(level: u32) -> u32 {
    ((1u32 << (3 * level)) - 1) / 7
}

/// A slot of the hash table.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Slot {
    Empty,
    Deleted,
    Live,
}

/// htslib's `khash` map from bin number to chunk list, emulated so its
/// slot order, the order a `.tbi` lists bins in, is htslib's: open
/// addressing with the number as its own hash, quadratic probing, a table
/// of a power of two doubled past 77% full, and the kick-out rehash.
#[derive(Debug, Default)]
struct BinTable {
    keys: Vec<u32>,
    vals: Vec<Vec<Chunk>>,
    slots: Vec<Slot>,
    size: usize,
    occupied: usize,
    upper: usize,
}

/// The fill at which a table of `n` slots grows (`__ac_HASH_UPPER`).
fn upper_bound(n: usize) -> usize {
    (n as f64 * 0.77 + 0.5) as usize
}

impl BinTable {
    /// `kh_resize`: a table of `want` slots (rounded up to a power of two,
    /// at least 4), rehashed in place.
    fn resize(&mut self, want: usize) {
        let new_n = want.next_power_of_two().max(4);
        if self.size >= upper_bound(new_n) {
            return;
        }
        let old_n = self.slots.len();
        let mut fresh = vec![Slot::Empty; new_n];
        if old_n < new_n {
            self.keys.resize(new_n, 0);
            self.vals.resize_with(new_n, Vec::new);
        }
        let mask = new_n - 1;
        for j in 0..old_n {
            if self.slots[j] != Slot::Live {
                continue;
            }
            let mut key = self.keys[j];
            let mut val = mem::take(&mut self.vals[j]);
            self.slots[j] = Slot::Deleted;
            loop {
                let mut i = key as usize & mask;
                let mut step = 0;
                while fresh[i] != Slot::Empty {
                    step += 1;
                    i = (i + step) & mask;
                }
                fresh[i] = Slot::Live;
                if i < old_n && self.slots[i] == Slot::Live {
                    mem::swap(&mut key, &mut self.keys[i]);
                    mem::swap(&mut val, &mut self.vals[i]);
                    self.slots[i] = Slot::Deleted;
                } else {
                    self.keys[i] = key;
                    self.vals[i] = val;
                    break;
                }
            }
        }
        self.keys.truncate(new_n);
        self.vals.truncate(new_n);
        self.slots = fresh;
        self.occupied = self.size;
        self.upper = upper_bound(new_n);
    }

    /// `kh_put`: the slot of `key`, added when absent.
    fn put(&mut self, key: u32) -> usize {
        if self.occupied >= self.upper {
            let n = self.slots.len();
            self.resize(if n > self.size << 1 { n - 1 } else { n + 1 });
        }
        let n = self.slots.len();
        let mask = n - 1;
        let mut i = key as usize & mask;
        let x = if self.slots[i] == Slot::Empty {
            i
        } else {
            let (last, mut site, mut step, mut x) = (i, n, 0, n);
            while self.slots[i] != Slot::Empty
                && (self.slots[i] == Slot::Deleted || self.keys[i] != key)
            {
                if self.slots[i] == Slot::Deleted {
                    site = i;
                }
                step += 1;
                i = (i + step) & mask;
                if i == last {
                    x = site;
                    break;
                }
            }
            if x == n {
                x = if self.slots[i] == Slot::Empty && site != n {
                    site
                } else {
                    i
                };
            }
            x
        };
        match self.slots[x] {
            Slot::Empty => {
                self.occupied += 1;
                self.size += 1;
            }
            Slot::Deleted => self.size += 1,
            Slot::Live => return x,
        }
        self.keys[x] = key;
        self.vals[x].clear();
        self.slots[x] = Slot::Live;
        x
    }

    /// `kh_get`: the slot holding `key`.
    fn get(&self, key: u32) -> Option<usize> {
        let n = self.slots.len();
        if n == 0 {
            return None;
        }
        let mask = n - 1;
        let (mut i, mut step) = (key as usize & mask, 0);
        let last = i;
        while self.slots[i] != Slot::Empty
            && (self.slots[i] == Slot::Deleted || self.keys[i] != key)
        {
            step += 1;
            i = (i + step) & mask;
            if i == last {
                return None;
            }
        }
        (self.slots[i] == Slot::Live).then_some(i)
    }

    /// htslib's `insert_to_b`: one more chunk for `bin`.
    fn insert(&mut self, bin: u32, chunk: Chunk) {
        let k = self.put(bin);
        self.vals[k].push(chunk);
    }

    /// The live slots, in table order.
    fn live(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.slots.len()).filter(|&k| self.slots[k] == Slot::Live)
    }

    /// htslib's `compress_binning`: a bin whose chunks span less than
    /// [`MIN_MARKER_DIST`] moves into its parent, when the parent exists;
    /// then chunks that start in the block another one ends in are joined.
    fn compress(&mut self) {
        for level in (1..=N_LVLS).rev() {
            let start = bin_first(level);
            for k in 0..self.slots.len() {
                let key = self.keys[k];
                if self.slots[k] != Slot::Live || key >= N_BINS || key < start {
                    continue;
                }
                if level < N_LVLS {
                    self.vals[k].sort_unstable_by_key(|c| c.0);
                }
                let p = &self.vals[k];
                let span = (p[p.len() - 1].1 >> 16).wrapping_sub(p[0].0 >> 16);
                if span >= MIN_MARKER_DIST {
                    continue;
                }
                let Some(parent) = self.get((key - 1) >> 3) else {
                    continue;
                };
                let moved = mem::take(&mut self.vals[k]);
                self.vals[parent].extend(moved);
                self.slots[k] = Slot::Deleted;
                self.size -= 1;
            }
        }
        if let Some(root) = self.get(0) {
            self.vals[root].sort_unstable_by_key(|c| c.0);
        }
        for k in 0..self.slots.len() {
            if self.slots[k] != Slot::Live || self.keys[k] >= N_BINS {
                continue;
            }
            let p = &mut self.vals[k];
            let mut m = 0;
            for l in 1..p.len() {
                if p[m].1 >> 16 >= p[l].0 >> 16 {
                    p[m].1 = p[m].1.max(p[l].1);
                } else {
                    m += 1;
                    p[m] = p[l];
                }
            }
            p.truncate(m + 1);
        }
    }
}

/// One sequence of the index.
#[derive(Debug, Default)]
struct Seq {
    bins: BinTable,
    /// Per 16 kb window, the offset of the first record that reaches it.
    linear: Vec<u64>,
}

/// htslib's index under construction (`hts_idx_t` and its `z` state).
#[derive(Debug, Default)]
struct Builder {
    names: Vec<Vec<u8>>,
    ids: HashMap<Vec<u8>, usize>,
    seqs: Vec<Seq>,
    save_bin: Option<u32>,
    save_tid: usize,
    save_off: u64,
    last_tid: Option<usize>,
    last_bin: Option<u32>,
    last_off: u64,
    off_beg: u64,
    last_coor: i64,
    n_mapped: u64,
    /// The tabix header the index is written with.
    conf: [i32; 6],
}

impl Builder {
    /// `hts_idx_init`: the first record starts at `offset`; `conf` is the
    /// tabix header.
    fn new(offset: u64, conf: [i32; 6]) -> Self {
        Builder {
            conf,
            save_off: offset,
            last_off: offset,
            off_beg: offset,
            ..Builder::default()
        }
    }

    /// The sequence id of `name`, added when new.
    fn tid(&mut self, name: &[u8]) -> usize {
        if let Some(id) = self.last_tid
            && self.names[id] == name
        {
            return id;
        }
        if let Some(&id) = self.ids.get(name) {
            return id;
        }
        let id = self.names.len();
        self.names.push(name.to_vec());
        self.ids.insert(name.to_vec(), id);
        id
    }

    /// `hts_idx_push`: one record of `tid` over `[beg, end)`, which ends
    /// at virtual offset `offset`.
    fn push(&mut self, tid: usize, beg: i64, end: i64, offset: u64) -> Result<(), String> {
        let name = || String::from_utf8_lossy(&self.names[tid]).into_owned();
        if beg > MAX_POS || end > MAX_POS {
            return Err(format!(
                "{n}:{beg}-{end} is past position {MAX_POS}, the last a .tbi index can hold",
                n = name()
            ));
        }
        if self.last_tid != Some(tid) {
            if tid < self.seqs.len() {
                return Err(format!(
                    "{n} starts again after another sequence: the records are not sorted",
                    n = name()
                ));
            }
            self.last_tid = Some(tid);
            self.last_bin = None;
        } else if self.last_coor > beg {
            return Err(format!(
                "positions go back on {n}: {} is followed by {}; the records are not sorted",
                self.last_coor + 1,
                beg + 1,
                n = name()
            ));
        }
        if end < beg {
            return Err(format!(
                "the record at {n}:{} ends at {end}, before it begins",
                beg + 1,
                n = name()
            ));
        }
        if tid == self.seqs.len() {
            self.seqs.push(Seq::default());
        }
        let end = end.max(1);
        let linear = &mut self.seqs[tid].linear;
        let (first, last) = (
            (beg >> MIN_SHIFT) as usize,
            ((end - 1) >> MIN_SHIFT) as usize,
        );
        if linear.len() <= last {
            linear.resize(last + 1, UNSET);
        }
        for o in linear.iter_mut().take(last + 1).skip(first) {
            if *o == UNSET {
                *o = self.last_off;
            }
        }
        let bin = reg2bin(beg, end);
        if self.last_bin != Some(bin) {
            if let Some(saved) = self.save_bin {
                self.seqs[self.save_tid]
                    .bins
                    .insert(saved, (self.save_off, self.last_off));
                if self.last_bin.is_none() {
                    let bins = &mut self.seqs[self.save_tid].bins;
                    bins.insert(META_BIN, (self.off_beg, self.last_off));
                    bins.insert(META_BIN, (self.n_mapped, 0));
                    self.n_mapped = 0;
                    self.off_beg = self.last_off;
                }
            }
            self.save_off = self.last_off;
            self.save_bin = Some(bin);
            self.last_bin = Some(bin);
            self.save_tid = tid;
        }
        self.n_mapped += 1;
        self.last_off = offset;
        self.last_coor = beg;
        Ok(())
    }

    /// `hts_idx_finish`, then the `.tbi` as `hts_idx_save` writes it,
    /// before compression, into `out`. `end` is the virtual offset of the
    /// file's end.
    fn finish(mut self, end: u64, out: &mut impl Write) -> io::Result<()> {
        if let Some(saved) = self.save_bin {
            let bins = &mut self.seqs[self.save_tid].bins;
            bins.insert(saved, (self.save_off, end));
            bins.insert(META_BIN, (self.off_beg, end));
            bins.insert(META_BIN, (self.n_mapped, 0));
        }
        let int = |out: &mut dyn Write, v: usize| {
            out.write_all(&i32::try_from(v).unwrap_or(i32::MAX).to_le_bytes())
        };
        out.write_all(b"TBI\x01")?;
        int(out, self.seqs.len())?;
        for v in self.conf {
            out.write_all(&v.to_le_bytes())?;
        }
        int(out, self.names.iter().map(|n| n.len() + 1).sum())?;
        for n in &self.names {
            out.write_all(n)?;
            out.write_all(&[0])?;
        }
        for seq in &mut self.seqs {
            // An empty window takes the offset of the next one: the last
            // window always has one.
            for w in (0..seq.linear.len().saturating_sub(1)).rev() {
                if seq.linear[w] == UNSET {
                    seq.linear[w] = seq.linear[w + 1];
                }
            }
            seq.bins.compress();
            let bins = &seq.bins;
            int(out, bins.size)?;
            for k in bins.live() {
                out.write_all(&bins.keys[k].to_le_bytes())?;
                int(out, bins.vals[k].len())?;
                for &(a, b) in &bins.vals[k] {
                    out.write_all(&a.to_le_bytes())?;
                    out.write_all(&b.to_le_bytes())?;
                }
            }
            int(out, seq.linear.len())?;
            for o in &seq.linear {
                out.write_all(&o.to_le_bytes())?;
            }
        }
        // No record without coordinates: htslib writes the count anyway.
        out.write_all(&0u64.to_le_bytes())
    }
}

/// An index read from a VCF or a positions file, not yet written.
#[derive(Debug)]
pub struct Index {
    builder: Builder,
    end: u64,
    /// How many records had a POS of 0 or below.
    pub clamped: u64,
}

impl Index {
    /// Writes the `.tbi`, inflated, into `out`.
    pub fn write_to(self, out: &mut impl Write) -> io::Result<()> {
        self.builder.finish(self.end, out)
    }

    /// The `.tbi`, inflated.
    pub fn into_bytes(self) -> Vec<u8> {
        let mut out = Vec::new();
        self.write_to(&mut out).expect("a Vec takes every write");
        out
    }

    /// Writes the `.tbi` to `path` as BGZF, atomically: a reader sees the
    /// whole index or none. With `force` it replaces an index there; without
    /// it, an index that appears while this one is built stays, and the
    /// error is `AlreadyExists`.
    pub fn write_tbi(self, path: &Path, force: bool) -> io::Result<()> {
        if force {
            let mut w = bgzf::Writer::new(AtomicFile::create(path)?, NonZero::<usize>::MIN)?;
            self.write_to(&mut w)?;
            return w.finish()?.commit();
        }
        // A hard link never replaces its target, so it publishes the
        // finished temporary file without the race of a check and a rename.
        let mut tmp = path.as_os_str().to_owned();
        tmp.push(format!(".tmp-{}", std::process::id()));
        let tmp = PathBuf::from(tmp);
        let done = File::create(&tmp)
            .and_then(|f| {
                let mut w = bgzf::Writer::new(io::BufWriter::new(f), NonZero::<usize>::MIN)?;
                self.write_to(&mut w)?;
                // finish flushes the BufWriter, so a write error shows here.
                w.finish().map(drop)
            })
            .and_then(|()| fs::hard_link(&tmp, path));
        let _ = fs::remove_file(&tmp);
        done
    }
}

/// Reads the BGZF VCF `input` and builds its index.
pub fn build(input: impl Read) -> io::Result<Index> {
    build_as(input, VCF_CONF)
}

/// Reads the BGZF tab-separated `input`, one base per line (sequence in
/// column 1, 1-based position in column 2), and builds the index
/// `tabix -s1 -b2 -e2` writes for it.
pub fn build_positions(input: impl Read) -> io::Result<Index> {
    build_as(input, POSITIONS_CONF)
}

/// A line's sequence name and 0-based `[beg, end)` under `-s1 -b2 -e2`, by
/// htslib's `tbx_parse1`: with the begin and end in one column, the end is
/// POS and the begin one before it; a begin below 0 moves to 0 (the flag),
/// an end below 1 to 1. `None` for a line it cannot place.
fn parse_position(line: &[u8]) -> Option<(&[u8], i64, i64, bool)> {
    let mut cols = line.split(|&c| c == b'\t' || c == 0);
    let name = cols.next()?;
    let pos = strtoll(cols.next()?, true)?;
    let beg = pos.saturating_sub(1);
    Some((name, beg.max(0), pos.max(1), beg < 0))
}

/// [`build`] with the tabix header `conf`, which also picks the parser.
fn build_as(input: impl Read, conf: [i32; 6]) -> io::Result<Index> {
    let mut lines = Lines::new(input);
    let mut line = Vec::new();
    let mut sv = Vec::new();
    let mut builder: Option<Builder> = None;
    let mut clamped = 0;
    let mut header_end = 0;
    let mut lineno = 0u64;
    while lines.read_line(&mut line)? {
        lineno += 1;
        if line.first() == Some(&b'#') {
            header_end = lines.tell();
            continue;
        }
        let b = builder.get_or_insert_with(|| Builder::new(header_end, conf));
        let parsed = if conf == VCF_CONF {
            parse(&line, &mut sv)
        } else {
            parse_position(&line)
        };
        let Some((name, beg, end, low)) = parsed else {
            let shown = String::from_utf8_lossy(&line[..line.len().min(60)]).into_owned();
            return Err(invalid(format!(
                "line {lineno} has no position to read: {shown:?}"
            )));
        };
        clamped += u64::from(low);
        let tid = b.tid(name);
        b.push(tid, beg, end, lines.tell())
            .map_err(|e| invalid(format!("line {lineno}: {e}")))?;
    }
    Ok(Index {
        builder: builder.unwrap_or_else(|| Builder::new(header_end, conf)),
        end: lines.tell(),
        clamped,
    })
}

const USAGE: &str = "usage: og vcf index [-f] VCF.gz\n";

/// `og vcf index [-f] VCF.gz`: writes `VCF.gz.tbi`. Exit 0 done, 1 an
/// error (nothing written), 2 usage.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let mut usage = |msg: &str| {
        let _ = writeln!(err, "{USAGE}og vcf index: error: {msg}");
        2
    };
    let (mut force, mut vcf) = (false, None);
    for a in args {
        match a.as_str() {
            "-h" | "--help" => {
                let _ = write!(out, "{USAGE}");
                return 0;
            }
            "-f" | "--force" => force = true,
            s if s.starts_with('-') => return usage(&format!("unrecognized arguments: {s}")),
            s if vcf.is_none() => vcf = Some(PathBuf::from(s)),
            s => return usage(&format!("unrecognized arguments: {s}")),
        }
    }
    let Some(vcf) = vcf else {
        return usage("the following arguments are required: VCF.gz");
    };
    let mut tbi = vcf.clone().into_os_string();
    tbi.push(".tbi");
    let tbi = PathBuf::from(tbi);
    if tbi.exists() && !force {
        let _ = writeln!(
            err,
            "og vcf index: {} exists; pass -f to write over it",
            tbi.display()
        );
        return 1;
    }
    let index = match File::open(&vcf).and_then(|f| build(BufReader::with_capacity(1 << 16, f))) {
        Ok(i) => i,
        Err(e) => {
            let _ = writeln!(err, "og vcf index: {}: {e}", vcf.display());
            return 1;
        }
    };
    if index.clamped > 0 {
        let _ = writeln!(
            err,
            "og vcf index: warning: {} records with a POS of 0 or below are indexed at the first base, as tabix does",
            index.clamped
        );
    }
    if let Err(e) = index.write_tbi(&tbi, force) {
        let _ = writeln!(err, "og vcf index: {}: {e}", tbi.display());
        return 1;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tabix::Reader;
    use std::process::{Command, Stdio};

    /// The repository root.
    fn repo() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
    }

    fn temp(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("og-tabix-write-{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        d.join(name)
    }

    /// Inflates a BGZF or gzip file.
    fn inflate(path: &Path) -> Vec<u8> {
        let mut out = Vec::new();
        crate::open(path).unwrap().read_to_end(&mut out).unwrap();
        out
    }

    /// The text as BGZF, one file per piece concatenated: each piece ends
    /// in an empty block, so empty blocks sit inside the file, and a
    /// piece boundary can fall at a line end or inside a line.
    fn bgzip(pieces: &[&[u8]]) -> Vec<u8> {
        let mut out = Vec::new();
        for p in pieces {
            let mut w = bgzf::Writer::new(Vec::new(), NonZero::<usize>::MIN).unwrap();
            w.write_all(p).unwrap();
            out.extend(w.finish().unwrap());
        }
        out
    }

    /// xorshift64*: the same numbers on every run.
    struct Rng(u64);
    impl Rng {
        fn below(&mut self, n: u64) -> u64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            self.0.wrapping_mul(0x2545_f491_4f6c_dd1d) % n
        }
    }

    /// A VCF over three sequences that crosses every level of bins: SNVs
    /// a few hundred bases apart, REFs up to 20 kb, END up to 200 Mb
    /// ahead (bins down to the root), an END before POS, keys that only
    /// end in `END=`, `<DEL>` and `<DUP:TANDEM>` with SVLEN, gVCF blocks
    /// with a sample LEN, and some lines ending in `\r\n`. Split at odd
    /// places, with empty blocks between the pieces.
    fn synthetic(seed: u64) -> (Vec<u8>, Vec<u8>) {
        let mut rng = Rng(seed);
        let mut t = b"##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tkid\n".to_vec();
        for seq in ["chr1", "chr2", "chrX"] {
            let mut pos = 1 + rng.below(50);
            for i in 0..4000 {
                let reflen = match rng.below(400) {
                    0 => 20_000,
                    1..=40 => 1 + rng.below(30),
                    _ => 1,
                } as usize;
                let (alt, info, sample) = match rng.below(40) {
                    0 => (
                        "C",
                        format!("SVTYPE=DEL;END={}", pos + 1 + rng.below(70_000)),
                        "0|1".to_owned(),
                    ),
                    1 => (
                        "C",
                        format!("SVTYPE=DUP;END={}", pos + 1 + rng.below(200_000_000)),
                        "0|1".to_owned(),
                    ),
                    2 => ("C", format!("END={}", pos + rng.below(3)), "0|1".to_owned()),
                    3 => (
                        "C",
                        format!("SVEND={};CIEND=-5,5", pos + 90_000),
                        "0|1".to_owned(),
                    ),
                    4 => (
                        "C",
                        format!("END={}", pos.saturating_sub(5)),
                        "0|1".to_owned(),
                    ),
                    5 => (
                        "<DEL>",
                        format!("SVLEN=-{}", 1 + rng.below(90_000)),
                        "0|1".to_owned(),
                    ),
                    6 => (
                        "C,<DUP:TANDEM>",
                        format!("SVLEN=.,{};CIPOS=-3,3", rng.below(40_000)),
                        "1|2".to_owned(),
                    ),
                    7 => (
                        "<*>",
                        "DP=30".to_owned(),
                        format!("0/0:{}", rng.below(30_000)),
                    ),
                    _ => ("C", "DP=30".to_owned(), "0|1".to_owned()),
                };
                let eol = if rng.below(15) == 0 { "\r\n" } else { "\n" };
                let r = "A".repeat(reflen);
                t.extend(
                    format!(
                        "{seq}\t{pos}\tid{i}\t{r}\t{alt}\t50\tPASS\t{info}\tGT:LEN\t{sample}{eol}"
                    )
                    .as_bytes(),
                );
                pos += rng.below(if i % 500 == 0 { 300_000 } else { 400 });
            }
        }
        let a = t.iter().position(|&c| c == b'#').unwrap() + 30;
        let b = t.len() / 3;
        let b = b + t[b..].iter().position(|&c| c == b'\n').unwrap() + 1;
        let c = 2 * t.len() / 3 + 7;
        let gz = bgzip(&[&t[..a], &t[a..b], &[], &t[b..c], &t[c..]]);
        (t, gz)
    }

    /// The data lines with their sequence and the 0-based interval the
    /// reader of [`crate::tabix`] gives them (it reads END and REF only, so
    /// its interval never ends after the one the index used).
    fn records(text: &[u8]) -> Vec<(Vec<u8>, Vec<u8>, i64, i64)> {
        let conf = crate::tabix::Conf {
            preset: crate::tabix::Preset::Vcf,
            zero_based: false,
            col_seq: 1,
            col_beg: 2,
            col_end: 0,
            meta: b'#',
            skip: 0,
        };
        text.split(|&c| c == b'\n')
            .map(|l| l.strip_suffix(b"\r").unwrap_or(l))
            .filter(|l| !l.is_empty() && l[0] != b'#')
            .map(|l| {
                let (s, b, e) = crate::tabix::interval(&conf, l).unwrap();
                let (_, ib, ie, _) = parse(l, &mut Vec::new()).unwrap();
                assert!(ib as u64 == b && ie as u64 >= e, "{l:?}");
                (l.to_vec(), s.to_vec(), b as i64, e as i64)
            })
            .collect()
    }

    /// Regions at the edges that matter: records' own ends, one base either
    /// side, the 16 kb windows and whole sequences.
    fn regions(recs: &[(Vec<u8>, Vec<u8>, i64, i64)]) -> Vec<(String, u64, u64)> {
        let mut out = Vec::new();
        for (i, (_, s, b, e)) in recs.iter().enumerate() {
            if i % 13 != 0 {
                continue;
            }
            let s = String::from_utf8(s.clone()).unwrap();
            let (b, e) = (*b as u64, *e as u64);
            let w = (b >> 14) << 14;
            for (x, y) in [
                (b, b + 1),
                (b.saturating_sub(1), b),
                (e.saturating_sub(1), e),
                (e, e + 1),
                (b.saturating_sub(300), e + 300),
                (w, w + 1),
                (w + 16383, w + 16385),
            ] {
                out.push((s.clone(), x, y));
            }
        }
        for s in ["chr1", "chr2", "chrX", "chr3"] {
            out.push((s.to_owned(), 0, 1 << 29));
            out.push((s.to_owned(), 5_000_000, 6_000_000));
        }
        out
    }

    /// How many regions a reader over `tbi` answers differently from a
    /// scan of the text, and how many lines the scan found.
    fn wrong_regions(text: &[u8], gz: &Path, tbi: &[u8], name: &str) -> (usize, usize) {
        let idx = temp(&format!("{name}.tbi"));
        let mut w = bgzf::Writer::new(File::create(&idx).unwrap(), NonZero::<usize>::MIN).unwrap();
        w.write_all(tbi).unwrap();
        w.finish().unwrap();
        let mut r = Reader::with_index(gz, &idx).unwrap();
        let recs = records(text);
        let (mut wrong, mut hits) = (0, 0);
        for (s, b, e) in regions(&recs) {
            let want: Vec<Vec<u8>> = recs
                .iter()
                .filter(|(_, rs, rb, re)| {
                    rs == s.as_bytes() && (*rb as u64) < e && (*re as u64) > b
                })
                .map(|(l, ..)| l.clone())
                .collect();
            hits += want.len();
            let got: Vec<Vec<u8>> = r.query(&s, b, e).collect::<io::Result<_>>().unwrap();
            wrong += usize::from(got != want);
        }
        (wrong, hits)
    }

    /// The byte offset, in an inflated `.tbi`, of the start of the first
    /// sequence's longest chunk (of a real bin, not the pseudo-bin), and
    /// how many chunks that sequence has.
    fn longest_chunk(tbi: &[u8]) -> (usize, usize) {
        let u32_at = |p: usize| u32::from_le_bytes(tbi[p..p + 4].try_into().unwrap()) as usize;
        let u64_at = |p: usize| u64::from_le_bytes(tbi[p..p + 8].try_into().unwrap());
        let mut p = 36 + u32_at(32);
        let n_bin = u32_at(p);
        p += 4;
        let (mut best, mut span, mut count) = (0, 0, 0);
        for _ in 0..n_bin {
            let (bin, n) = (u32_at(p), u32_at(p + 4));
            p += 8;
            for c in 0..n {
                let at = p + 16 * c;
                if bin < N_BINS as usize {
                    let len = u64_at(at + 8) - u64_at(at);
                    count += 1;
                    if len > span {
                        (best, span) = (at, len);
                    }
                }
            }
            p += 16 * n;
        }
        (best, count)
    }

    #[test]
    fn strtoll_reads_as_c_does() {
        for (s, want) in [
            ("12", Some(12)),
            ("  +7x", Some(7)),
            ("-3", Some(-3)),
            ("0x1f", Some(31)),
            ("0X", Some(0)),
            ("010", Some(8)),
            ("08", Some(0)),
            ("99999999999999999999", Some(i64::MAX)),
            ("-99999999999999999999", Some(i64::MIN)),
            ("-9223372036854775808", Some(i64::MIN)),
            ("", None),
            ("x1", None),
            ("-", None),
            (".", None),
        ] {
            assert_eq!(strtoll(s.as_bytes(), true), want, "{s:?}");
        }
        // atoll is base 10 and reads 0 when there is no digit.
        assert_eq!(atoll(b"0x1f"), 0);
        assert_eq!(atoll(b"010"), 10);
        assert_eq!(atoll(b"-500;CIPOS"), -500);
        assert_eq!(atoll(b"."), 0);
    }

    #[test]
    fn intervals_follow_tbx_parse1() {
        let p = |l: &str| {
            parse(l.as_bytes(), &mut Vec::new()).map(|(s, b, e, c)| (s.to_vec(), b, e, c))
        };
        let c = b"c".to_vec();
        assert_eq!(
            p("c\t10\t.\tACG\tA\t.\t.\t."),
            Some((c.clone(), 9, 12, false))
        );
        assert_eq!(
            p("c\t10\t.\tA\t<DEL>\t.\t.\tEND=500"),
            Some((c.clone(), 9, 500, false))
        );
        // END after POS counts, but never shortens the REF span.
        assert_eq!(
            p("c\t10\t.\tAAAA\tA\t.\t.\tX;END=11"),
            Some((c.clone(), 9, 13, false))
        );
        assert_eq!(
            p("c\t10\t.\tAA\tT\t.\t.\tEND=9"),
            Some((c.clone(), 9, 11, false))
        );
        assert_eq!(
            p("c\t10\t.\tA\tT\t.\t.\tSVEND=900"),
            Some((c.clone(), 9, 10, false))
        );
        assert_eq!(
            p("c\t10\t.\tA\tT\t.\t.\tEND=0x20"),
            Some((c.clone(), 9, 32, false))
        );
        assert_eq!(
            p("c\t10\t.\tA\tT\t.\t.\tEND=."),
            Some((c.clone(), 9, 10, false))
        );
        // SVLEN counts for <DEL>, <DUP>, <INV> and <CNV>, with a subtype
        // too, as its absolute value; not for <INS> or a plain ALT.
        assert_eq!(
            p("c\t10\t.\tA\t<DEL>\t.\t.\tSVLEN=-500"),
            Some((c.clone(), 9, 509, false))
        );
        assert_eq!(
            p("c\t10\t.\tA\t<DUP:TANDEM>\t.\t.\tX=1;SVLEN=70"),
            Some((c.clone(), 9, 79, false))
        );
        assert_eq!(
            p("c\t10\t.\tA\t<INS>\t.\t.\tSVLEN=900"),
            Some((c.clone(), 9, 10, false))
        );
        assert_eq!(
            p("c\t10\t.\tA\tT\t.\t.\tSVLEN=900"),
            Some((c.clone(), 9, 10, false))
        );
        assert_eq!(
            p("c\t10\t.\tA\t<DEL\t.\t.\tSVLEN=900"),
            Some((c.clone(), 9, 10, false))
        );
        // One SVLEN per ALT, in order; the next value starts after the
        // next comma anywhere in INFO, as htslib reads it.
        assert_eq!(
            p("c\t10\t.\tA\tT,<DEL>\t.\t.\tSVLEN=.,300"),
            Some((c.clone(), 9, 309, false))
        );
        assert_eq!(
            p("c\t10\t.\tA\t<DEL>,<INV>\t.\t.\tSVLEN=-5;CIPOS=-1,800"),
            Some((c.clone(), 9, 809, false))
        );
        assert_eq!(
            p("c\t10\t.\tA\t<DEL>\t.\t.\tEND=600;SVLEN=-50"),
            Some((c.clone(), 9, 600, false))
        );
        // A gVCF block: the largest sample LEN, when FORMAT has LEN.
        assert_eq!(
            p("c\t10\t.\tA\t<*>\t.\t.\t.\tGT:LEN\t0/0:40\t0/0:7"),
            Some((c.clone(), 9, 49, false))
        );
        assert_eq!(
            p("c\t10\t.\tA\t<NON_REF>\t.\t.\t.\tGT\t0/0:40"),
            Some((c.clone(), 9, 10, false))
        );
        assert_eq!(
            p("c\t10\t.\tA\tT\t.\t.\t.\tGT:LEN\t0/0:40"),
            Some((c.clone(), 9, 10, false))
        );
        // POS 0 goes to the first base; no REF ends where it begins.
        assert_eq!(p("c\t0\t.\tA"), Some((c.clone(), 0, 1, true)));
        assert_eq!(p("c\t7"), Some((c.clone(), 6, 6, false)));
        assert_eq!(p("c\t5x\t.\tA"), Some((c.clone(), 4, 5, false)));
        assert_eq!(p("c\tx\t.\tA"), None);
        assert_eq!(p("c"), None);
        assert_eq!(p(""), None);
    }

    #[test]
    fn bins_are_htslibs() {
        assert_eq!(reg2bin(0, 1), 4681);
        assert_eq!(reg2bin(16383, 16385), 585);
        assert_eq!(reg2bin(0, 1 << 29), 0);
        assert_eq!(reg2bin((1 << 29) - 1, 1 << 29), 4681 + 32767);
        // An empty interval at 1 is binned on the base before it.
        assert_eq!(reg2bin(1, 1), 4681);
    }

    /// Every VCF in the repository whose committed `.tbi` htslib wrote:
    /// og writes the same bytes, once inflated.
    #[test]
    fn rebuilds_the_committed_htslib_indexes_byte_for_byte() {
        let pairs = [
            "assets/stub/sample.vcf.gz",
            "tests/fixtures/pgx/summary/no-contig.vcf.gz",
            "tests/fixtures/pgx/summary/sample.vcf.gz",
            "tests/parity/inputs/bloodtype/abo.vcf.gz",
            "tests/parity/inputs/bloodtype/blocks.vcf.gz",
            "tests/parity/inputs/bloodtype/bombay.vcf.gz",
            "tests/parity/inputs/bloodtype/gvcf.vcf.gz",
            "tests/parity/inputs/bloodtype/nochr.vcf.gz",
            "tests/parity/inputs/bloodtype/vcf-quirks.vcf.gz",
            "tests/parity/inputs/lookup/ann.vcf.gz",
            "tests/parity/inputs/lookup/g.vcf.gz",
        ];
        for rel in pairs {
            let vcf = repo().join(rel);
            let mut tbi = vcf.clone().into_os_string();
            tbi.push(".tbi");
            let want = inflate(Path::new(&tbi));
            let got = build(File::open(&vcf).unwrap()).unwrap().into_bytes();
            assert_eq!(got, want, "{rel}: the index differs from htslib's");
        }
    }

    /// A synthetic VCF across bin levels and odd block boundaries: every
    /// region query through og's index returns what a scan returns. The
    /// planted control moves the start of the longest chunk to its end,
    /// and the same check must then find a wrong region.
    #[test]
    fn region_queries_through_the_index_equal_a_scan() {
        let (text, gz) = synthetic(7);
        let path = temp("synthetic.vcf.gz");
        std::fs::write(&path, &gz).unwrap();
        let tbi = build(&gz[..]).unwrap().into_bytes();
        let (wrong, hits) = wrong_regions(&text, &path, &tbi, "good");
        assert!(hits > 10_000, "the regions found only {hits} lines");
        assert_eq!(wrong, 0, "{wrong} regions differ from a scan");

        let mut bad = tbi.clone();
        let (p, chunks) = longest_chunk(&bad);
        assert!(chunks > 0, "no chunk to move");
        let end: [u8; 8] = bad[p + 8..p + 16].try_into().unwrap();
        bad[p..p + 8].copy_from_slice(&end);
        let (wrong, _) = wrong_regions(&text, &path, &bad, "planted");
        assert!(
            wrong > 0,
            "a moved chunk start went unseen by the region check"
        );
    }

    #[test]
    fn line_offsets_follow_bgzf_tell() {
        // A line that ends a block is followed by the next block's start,
        // even when empty blocks come between; the end of the file leaves
        // the offset after the last block that held bytes.
        let gz = bgzip(&[b"a\n", b"", b"bc", b"d\ne"]);
        let mut l = Lines::new(&gz[..]);
        let mut line = Vec::new();
        let starts: Vec<u64> = {
            let mut s = vec![0];
            let mut at = 0;
            while at < gz.len() {
                at += u16::from_le_bytes([gz[at + 16], gz[at + 17]]) as usize + 1;
                s.push(at as u64);
            }
            s
        };
        // Blocks: "a\n", EOF, EOF, "bc", EOF, "d\ne", EOF.
        assert!(l.read_line(&mut line).unwrap());
        assert_eq!((line.as_slice(), l.tell()), (&b"a"[..], starts[1] << 16));
        assert!(l.read_line(&mut line).unwrap());
        assert_eq!(
            (line.as_slice(), l.tell()),
            (&b"bcd"[..], (starts[5] << 16) | 2)
        );
        assert!(l.read_line(&mut line).unwrap());
        assert_eq!((line.as_slice(), l.tell()), (&b"e"[..], starts[6] << 16));
        assert!(!l.read_line(&mut line).unwrap());
        assert_eq!(l.tell(), starts[6] << 16);
    }

    #[test]
    fn refuses_what_htslib_refuses() {
        let head = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let err = |body: &str| {
            let gz = bgzip(&[format!("{head}{body}").as_bytes()]);
            build(&gz[..]).unwrap_err().to_string()
        };
        assert!(
            err("c\t5\t.\tA\nc\t4\t.\tA\n").contains("positions go back on c: 5 is followed by 4")
        );
        assert!(err("a\t5\t.\tA\nb\t4\t.\tA\na\t6\t.\tA\n").contains("a starts again"));
        assert!(err("c\t1\t.\tA\t.\t.\t.\tEND=536870913\n").contains("past position"));
        assert!(err("c\t536870913\t.\tA\n").contains("past position"));
        assert!(err("c\t5\t.\tA\n\nc\t9\t.\tA\n").contains("line 3 has no position"));
        assert!(err("c\t5\t.\tA\nc\tx\n").contains("line 3 has no position"));
        // The last base a .tbi holds is fine.
        let gz = bgzip(&[format!("{head}c\t536870912\t.\tA\n").as_bytes()]);
        assert!(build(&gz[..]).is_ok());
        // Plain gzip and plain text are not BGZF.
        let mut plain = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        plain.write_all(head.as_bytes()).unwrap();
        let plain = plain.finish().unwrap();
        for data in [plain, head.as_bytes().to_vec()] {
            let e = build(&data[..]).unwrap_err().to_string();
            assert!(e.contains("not BGZF"), "{e}");
        }
    }

    #[test]
    fn an_empty_or_header_only_file_indexes_no_sequence() {
        let head = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        for data in [bgzip(&[b""]), bgzip(&[head.as_bytes()])] {
            let tbi = build(&data[..]).unwrap().into_bytes();
            assert_eq!(&tbi[4..8], &0i32.to_le_bytes());
            assert_eq!(tbi.len(), 4 + 4 + 24 + 4 + 8);
        }
    }

    fn run(args: &[&str]) -> (u8, String, String) {
        let args: Vec<String> = args.iter().map(|s| s.to_string()).collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = cli(&args, &mut out, &mut err);
        (
            code,
            String::from_utf8(out).unwrap(),
            String::from_utf8(err).unwrap(),
        )
    }

    #[test]
    fn the_verb_writes_the_index_beside_the_vcf() {
        let src = repo().join("tests/parity/inputs/bloodtype/abo.vcf.gz");
        let vcf = temp("abo.vcf.gz");
        std::fs::copy(&src, &vcf).unwrap();
        let tbi = temp("abo.vcf.gz.tbi");
        let _ = std::fs::remove_file(&tbi);
        let v = vcf.to_str().unwrap();
        assert_eq!(run(&[v]).0, 0);
        let mut want = src.into_os_string();
        want.push(".tbi");
        assert_eq!(inflate(&tbi), inflate(Path::new(&want)));
        // It does not write over an index unless told to.
        let (code, _, err) = run(&[v]);
        assert_eq!(code, 1);
        assert!(err.contains("pass -f"), "{err}");
        assert_eq!(run(&["-f", v]).0, 0);
        // Usage errors come before any file is read.
        for bad in [&["--og-scaffold-probe", v][..], &[v, v], &[]] {
            let (code, _, err) = run(bad);
            assert_eq!(code, 2, "{bad:?}");
            assert!(err.starts_with("usage: og vcf index"), "{err}");
        }
        assert_eq!(run(&["--help"]).1, USAGE);
        // A failed build writes nothing.
        let broken = temp("broken.vcf.gz");
        std::fs::write(&broken, b"not bgzf").unwrap();
        let (code, _, err) = run(&[broken.to_str().unwrap()]);
        assert_eq!(code, 1);
        assert!(err.contains("not BGZF"), "{err}");
        assert!(!temp("broken.vcf.gz.tbi").exists());
        // An empty file is not BGZF; an empty BGZF stream is a VCF with no
        // lines.
        let empty = temp("empty.vcf.gz");
        std::fs::write(&empty, b"").unwrap();
        let (code, _, err) = run(&[empty.to_str().unwrap()]);
        assert_eq!(code, 1);
        assert!(err.contains("not BGZF"), "{err}");
        assert!(!temp("empty.vcf.gz.tbi").exists());
        let eof = temp("eof.vcf.gz");
        std::fs::write(&eof, bgzf::EOF_BLOCK).unwrap();
        assert_eq!(run(&[eof.to_str().unwrap()]).0, 0);
        assert!(temp("eof.vcf.gz.tbi").exists());
    }

    #[test]
    fn without_force_an_index_that_appeared_meanwhile_stays() {
        let src = repo().join("tests/parity/inputs/bloodtype/abo.vcf.gz");
        let index = build(BufReader::new(File::open(&src).unwrap())).unwrap();
        let tbi = temp("raced.vcf.gz.tbi");
        std::fs::write(&tbi, b"another writer").unwrap();
        let e = index.write_tbi(&tbi, false).unwrap_err();
        assert_eq!(e.kind(), io::ErrorKind::AlreadyExists);
        assert_eq!(std::fs::read(&tbi).unwrap(), b"another writer");
        let index = build(BufReader::new(File::open(&src).unwrap())).unwrap();
        index.write_tbi(&tbi, true).unwrap();
        assert_ne!(std::fs::read(&tbi).unwrap(), b"another writer");
    }

    /// A positions file like step 11's allele-frequency files: a header,
    /// a POS of 0, repeated positions, gaps up to 40 kb, three sequences,
    /// split at odd places.
    fn positions(seed: u64) -> Vec<u8> {
        let mut rng = Rng(seed);
        let mut t = b"#CHROM\tPOS\tREF,ALT\tAF\nchr1\t0\tA,C\t0.5\n".to_vec();
        for seq in ["chr1", "chr2", "chr22"] {
            let mut pos = 1;
            for _ in 0..20_000 {
                pos += match rng.below(100) {
                    0 => 40_000,
                    1..=9 => 0,
                    _ => rng.below(3000),
                };
                t.extend(format!("{seq}\t{pos}\tA,G\t0.{}\n", rng.below(99)).bytes());
            }
        }
        let (a, b) = (t.len() / 3, t.len() / 3 + 12_345);
        bgzip(&[&t[..a], &t[a..b], &t[b..]])
    }

    #[test]
    fn a_position_line_covers_its_one_base() {
        let c = b"c".as_slice();
        assert_eq!(parse_position(b"c\t10\tA,G\t0.5"), Some((c, 9, 10, false)));
        assert_eq!(parse_position(b"c\t0\tA,G"), Some((c, 0, 1, true)));
        assert_eq!(parse_position(b"c\t0x10"), Some((c, 15, 16, false)));
        assert_eq!(parse_position(b"c\tx\tA,G"), None);
        assert_eq!(parse_position(b"c"), None);
        // The header says generic, sequence 1, begin 2, end 2, `#`, no skip.
        let gz = bgzip(&[b"#h\nc\t5\tA,G\t0.5\n"]);
        let tbi = build_positions(&gz[..]).unwrap().into_bytes();
        let conf: Vec<u8> = [0i32, 1, 2, 2, 35, 0]
            .iter()
            .flat_map(|v| v.to_le_bytes())
            .collect();
        assert_eq!(&tbi[8..32], &conf[..]);
    }

    /// htslib's own `tabix -s1 -b2 -e2`, in the pipeline's bcftools image,
    /// writes the index [`build_positions`] writes (inflated). In CI this
    /// must run; elsewhere, with no Docker, it says so and passes.
    #[test]
    fn htslib_tabix_writes_the_same_positions_index() {
        let docker = Command::new("docker")
            .arg("version")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok_and(|s| s.success());
        if !docker {
            assert!(
                std::env::var_os("CI").is_none(),
                "CI must run the htslib comparison, and Docker is missing"
            );
            eprintln!("no Docker here: the htslib comparison runs in CI");
            return;
        }
        let gz = positions(5);
        let ours = build_positions(&gz[..]).unwrap().into_bytes();
        let htslib = in_image(
            &[("x.tsv.gz", &gz[..])],
            "tabix -s1 -b2 -e2 x.tsv.gz && cat x.tsv.gz.tbi",
        );
        let mut want = Vec::new();
        crate::reader(io::Cursor::new(htslib))
            .unwrap()
            .read_to_end(&mut want)
            .unwrap();
        assert!(ours.len() > 1000, "a {} byte index", ours.len());
        assert!(ours == want, "tabix -s1 -b2 -e2 writes another index");
    }

    /// The image the pipeline runs bcftools and tabix from.
    fn bcftools_image() -> String {
        let env = std::fs::read_to_string(repo().join("versions.env")).unwrap();
        env.lines()
            .find_map(|l| l.strip_prefix("BCFTOOLS_IMAGE=\""))
            .and_then(|l| l.split('"').next())
            .expect("BCFTOOLS_IMAGE in versions.env")
            .to_owned()
    }

    /// Runs `sh -c script` in the bcftools image, in a work folder holding
    /// `files`, which go in on stdin (a bind mount fails where Docker runs
    /// in a VM); its stdout.
    fn in_image(files: &[(&str, &[u8])], script: &str) -> Vec<u8> {
        let mut prelude = "mkdir -p /tmp/w && cd /tmp/w".to_owned();
        for (name, data) in files {
            prelude.push_str(&format!(
                " && dd bs={} count=1 iflag=fullblock of={name} 2>/dev/null",
                data.len()
            ));
        }
        let mut child = Command::new("docker")
            .args(["run", "-i", "--rm", "--platform", "linux/amd64"])
            .args([&bcftools_image(), "sh", "-c"])
            .arg(format!("{prelude} && {script}"))
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .expect("docker runs");
        let mut stdin = child.stdin.take().unwrap();
        let input: Vec<u8> = files.iter().flat_map(|(_, d)| d.iter().copied()).collect();
        let feeder = std::thread::spawn(move || stdin.write_all(&input));
        let out = child.wait_with_output().unwrap();
        feeder.join().unwrap().unwrap();
        assert!(out.status.success(), "the image failed: {script}");
        out.stdout
    }

    /// htslib itself, in the pipeline's bcftools image: `bcftools index -t`
    /// (htslib 1.24, the htslib of the pysam the pipeline ran before) writes
    /// the index og writes (inflated), and `bcftools view -r` over og's
    /// index returns what it returns over htslib's. The image's `tabix`
    /// binary is htslib 1.19, whose rules predate SVLEN and sample LEN, so
    /// it is not the reference. In CI this must run; elsewhere, with no
    /// Docker, it says so and passes.
    #[test]
    fn htslib_writes_the_same_index_and_reads_it_the_same() {
        let docker = Command::new("docker")
            .arg("version")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok_and(|s| s.success());
        if !docker {
            assert!(
                std::env::var_os("CI").is_none(),
                "CI must run the htslib comparison, and Docker is missing"
            );
            eprintln!("no Docker here: the htslib comparison runs in CI");
            return;
        }
        let mut inputs: Vec<(String, Vec<u8>)> = vec![("synthetic.vcf.gz".into(), synthetic(11).1)];
        let quirks = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                      c\t0\t.\tA\tC\t.\t.\t.\r\nc\t3\t.\tA\tC\t.\t.\tEND=0x50000\n\
                      #late header\nc\t16384\t.\tAC\tA\t.\t.\t.\nc\t16390\t.\tAC\t<DEL>,<INV>\t.\t.\tSVLEN=-5;CIPOS=-1,80000\n\
                      c\t16400\t.\tA\t<*>\t.\t.\t.\tGT:LEN\t0/0:3\t0/0:40000\nd\t2\t.\tA\tC\t.\t.\t.";
        inputs.push(("quirks.vcf.gz".into(), bgzip(&[quirks.as_bytes()])));
        for f in [
            "hg002_benchmark.vcf.gz",
            "hg002_chr20.g.vcf.gz",
            "hg002_smvar_trf.vcf.gz",
            "hg002_vep.vcf.gz",
        ] {
            let p = repo().join("rust/crates/og-vcf/tests/data").join(f);
            inputs.push((f.into(), std::fs::read(p).unwrap()));
        }
        let inflate_bytes = |b: &[u8]| {
            let mut out = Vec::new();
            crate::reader(io::Cursor::new(b.to_vec()))
                .unwrap()
                .read_to_end(&mut out)
                .unwrap();
            out
        };
        for (name, gz) in &inputs {
            let ours = build(&gz[..]).unwrap().into_bytes();
            let vcf = [("x.vcf.gz", &gz[..])];
            let htslib = in_image(&vcf, "bcftools index -t x.vcf.gz && cat x.vcf.gz.tbi");
            assert_eq!(
                ours,
                inflate_bytes(&htslib),
                "{name}: bcftools index -t writes another index"
            );
            if name == "quirks.vcf.gz" {
                // A header line among the records is no VCF: bcftools view
                // refuses it.
                continue;
            }
            // Every sequence whole, then 2 kb regions 99,991 bases apart
            // over its first 3 Mb.
            let view = |tbi: &[u8]| {
                in_image(
                    &[("x.vcf.gz", &gz[..]), ("x.vcf.gz.tbi", tbi)],
                    "set -e; for s in $(tabix -l x.vcf.gz); do \
                     bcftools view -H -r \"$s\" x.vcf.gz 2>>/tmp/log; \
                     i=1; while [ $i -lt 3000000 ]; do echo \"== $s:$i\"; \
                     bcftools view -H -r \"$s:$i-$((i + 1999))\" x.vcf.gz 2>>/tmp/log > /tmp/r; \
                     cut -f1-2 /tmp/r; i=$((i + 99991)); done; done",
                )
            };
            let mut og_tbi = Vec::new();
            {
                let mut w = bgzf::Writer::new(&mut og_tbi, NonZero::<usize>::MIN).unwrap();
                w.write_all(&ours).unwrap();
                w.finish().unwrap();
            }
            let (with_og, with_htslib) = (view(&og_tbi), view(&htslib));
            let records = with_og
                .split(|&c| c == b'\n')
                .filter(|l| !l.is_empty() && !l.starts_with(b"=="))
                .count();
            assert!(
                records >= 10,
                "{name}: bcftools view returned {records} records"
            );
            assert!(
                with_og == with_htslib,
                "{name}: bcftools view -r differs over og's index"
            );
        }
    }
}
