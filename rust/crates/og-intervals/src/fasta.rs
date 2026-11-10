//! Random access into an uncompressed FASTA through its samtools `.fai`.
//!
//! [`Indexed`] reads the bases of a region as `samtools faidx` prints them,
//! byte for byte (case kept). It reads the file in 64 KiB pages and keeps
//! the 64 most recently used, so lookups sorted by position, as a score
//! file's rows mostly are, read each page once.
//!
//! A gzip or BGZF FASTA is refused when it is opened: the `.fai` offsets
//! of a compressed file count inflated bytes, so seeking by them would read
//! garbage.

use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

/// Bytes per page.
const PAGE: u64 = 1 << 16;
/// Pages kept.
const PAGES: usize = 64;

/// One line of a `.fai`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Entry {
    /// Bases in the sequence.
    pub len: u64,
    /// Byte offset of its first base.
    pub offset: u64,
    /// Bases per full line.
    pub line_bases: u64,
    /// Bytes per full line, line ending included.
    pub line_width: u64,
}

/// Why a FASTA could not be opened or read.
#[derive(Debug)]
pub enum FastaError {
    Io(io::Error),
    /// The FASTA is gzip or BGZF.
    Compressed(PathBuf),
    /// A `.fai` line that cannot be used: its line (1-based) and why.
    Index {
        line: usize,
        why: String,
    },
    /// The bytes where the `.fai` puts a base are not one: the index is
    /// not this file's.
    Stale {
        name: String,
        pos: u64,
    },
}

impl fmt::Display for FastaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FastaError::Io(e) => write!(f, "{e}"),
            FastaError::Compressed(p) => write!(
                f,
                "{} is compressed; random access needs an uncompressed FASTA",
                p.display()
            ),
            FastaError::Index { line, why } => write!(f, ".fai line {line}: {why}"),
            FastaError::Stale { name, pos } => write!(
                f,
                "the .fai does not match the FASTA: {name} position {} is not a base",
                pos + 1
            ),
        }
    }
}

impl std::error::Error for FastaError {}

impl From<io::Error> for FastaError {
    fn from(e: io::Error) -> Self {
        FastaError::Io(e)
    }
}

/// A FASTA opened with its `.fai`.
pub struct Indexed {
    file: File,
    entries: HashMap<String, Entry>,
    names: Vec<String>,
    /// Page number to its slot in `pages`.
    slots: HashMap<u64, usize>,
    /// Cached pages: number, last use, bytes.
    pages: Vec<(u64, u64, Vec<u8>)>,
    clock: u64,
}

/// Parses a `.fai`: `name length offset line_bases line_width`, tab
/// separated; a sixth column (a FASTQ index) is ignored.
pub fn read_fai(input: impl BufRead) -> Result<Vec<(String, Entry)>, FastaError> {
    let mut out = Vec::new();
    for (i, line) in input.lines().enumerate() {
        let line = line?;
        let n = i + 1;
        let bad = |why: String| FastaError::Index { line: n, why };
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 5 {
            return Err(bad(format!("{} columns, want 5", c.len())));
        }
        let num = |k: usize| {
            c[k].parse::<u64>()
                .map_err(|_| bad(format!("column {} '{}' is not a whole number", k + 1, c[k])))
        };
        let e = Entry {
            len: num(1)?,
            offset: num(2)?,
            line_bases: num(3)?,
            line_width: num(4)?,
        };
        if e.len > 0 && (e.line_bases == 0 || e.line_width < e.line_bases) {
            return Err(bad(format!(
                "{} bases per line in {} bytes",
                e.line_bases, e.line_width
            )));
        }
        out.push((c[0].to_owned(), e));
    }
    Ok(out)
}

impl Indexed {
    /// Opens `path` with `path.fai`.
    pub fn open(path: impl AsRef<Path>) -> Result<Indexed, FastaError> {
        let path = path.as_ref();
        let mut fai = path.as_os_str().to_owned();
        fai.push(".fai");
        Indexed::with_index(path, PathBuf::from(fai))
    }

    /// Opens `path` with the `.fai` at `fai`.
    pub fn with_index(
        path: impl AsRef<Path>,
        fai: impl AsRef<Path>,
    ) -> Result<Indexed, FastaError> {
        let path = path.as_ref();
        let mut file = File::open(path)?;
        let mut magic = [0u8; 2];
        let n = file.read(&mut magic)?;
        if n == 2 && magic == [0x1f, 0x8b] {
            return Err(FastaError::Compressed(path.to_owned()));
        }
        let rows = read_fai(BufReader::new(File::open(fai)?))?;
        let names = rows.iter().map(|(n, _)| n.clone()).collect();
        Ok(Indexed {
            file,
            entries: rows.into_iter().collect(),
            names,
            slots: HashMap::new(),
            pages: Vec::with_capacity(PAGES),
            clock: 0,
        })
    }

    /// The sequences of the `.fai`, in its order.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// The `.fai` line of `name`.
    pub fn entry(&self, name: &str) -> Option<Entry> {
        self.entries.get(name).copied()
    }

    /// The bases of `[beg, end)` of `name` (0-based), as stored (case
    /// kept); the region is cut at the sequence's end. `None` when the
    /// `.fai` has no such sequence.
    pub fn fetch(&mut self, name: &str, beg: u64, end: u64) -> Result<Option<Vec<u8>>, FastaError> {
        let Some(e) = self.entry(name) else {
            return Ok(None);
        };
        let end = end.min(e.len);
        let mut out = Vec::with_capacity(end.saturating_sub(beg) as usize);
        let mut pos = beg;
        while pos < end {
            // The rest of pos's line, up to end.
            let line_end = end.min((pos / e.line_bases + 1) * e.line_bases + 1);
            let at = e.offset + pos / e.line_bases * e.line_width + pos % e.line_bases;
            let start = out.len();
            self.read_into(at, line_end - pos, &mut out)?;
            if let Some(k) = out[start..].iter().position(|c| !c.is_ascii_graphic()) {
                return Err(FastaError::Stale {
                    name: name.to_owned(),
                    pos: pos + k as u64,
                });
            }
            pos = line_end;
        }
        Ok(Some(out))
    }

    /// The base at 0-based `pos` of `name`; `None` for an unknown sequence
    /// or a position past its end.
    pub fn base(&mut self, name: &str, pos: u64) -> Result<Option<u8>, FastaError> {
        Ok(self
            .fetch(name, pos, pos.saturating_add(1))?
            .and_then(|b| b.first().copied()))
    }

    /// Appends `n` bytes of the file from offset `at`.
    fn read_into(&mut self, mut at: u64, mut n: u64, out: &mut Vec<u8>) -> Result<(), FastaError> {
        while n > 0 {
            let page = self.page(at / PAGE)?;
            let from = (at % PAGE) as usize;
            let take = (n as usize).min(page.len().saturating_sub(from));
            if take == 0 {
                return Err(FastaError::Io(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "the .fai points past the end of the FASTA",
                )));
            }
            out.extend_from_slice(&page[from..from + take]);
            at += take as u64;
            n -= take as u64;
        }
        Ok(())
    }

    /// Page `no`, read on a miss over the least recently used one.
    fn page(&mut self, no: u64) -> io::Result<&[u8]> {
        self.clock += 1;
        let slot = match self.slots.get(&no) {
            Some(&s) => s,
            None => {
                let s = if self.pages.len() < PAGES {
                    self.pages.push((no, 0, Vec::new()));
                    self.pages.len() - 1
                } else {
                    let s = (0..PAGES)
                        .min_by_key(|&s| self.pages[s].1)
                        .expect("a full cache");
                    self.slots.remove(&self.pages[s].0);
                    s
                };
                let buf = &mut self.pages[s].2;
                buf.clear();
                self.file.seek(SeekFrom::Start(no * PAGE))?;
                (&mut self.file).take(PAGE).read_to_end(buf)?;
                self.pages[s].0 = no;
                self.slots.insert(no, s);
                s
            }
        };
        self.pages[slot].1 = self.clock;
        Ok(&self.pages[slot].2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp(name: &str, data: &[u8]) -> PathBuf {
        let p = std::env::temp_dir().join(format!("og-fasta-{}-{name}", std::process::id()));
        std::fs::write(&p, data).unwrap();
        p
    }

    /// A FASTA of the given records with `width` bases per line and `eol`
    /// line endings, and its `.fai` as samtools writes it.
    fn fasta(records: &[(&str, Vec<u8>)], width: usize, eol: &[u8]) -> (Vec<u8>, String) {
        let (mut fa, mut fai) = (Vec::new(), String::new());
        for (name, seq) in records {
            fa.extend(format!(">{name} some description").as_bytes());
            fa.extend(eol);
            let offset = fa.len();
            for line in seq.chunks(width) {
                fa.extend(line);
                fa.extend(eol);
            }
            fai += &format!(
                "{name}\t{}\t{offset}\t{width}\t{}\n",
                seq.len(),
                width + eol.len()
            );
        }
        (fa, fai)
    }

    fn bases(seed: u64, n: usize) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGTNacgtn"[(x % 10) as usize]
            })
            .collect()
    }

    #[test]
    fn regions_equal_the_sequence_at_line_and_page_edges() {
        let recs = vec![
            ("chr1", bases(1, 300_000)),
            ("chr2", bases(2, 70_001)),
            ("chrM", bases(3, 7)),
        ];
        for (width, eol) in [(60, &b"\n"[..]), (61, b"\r\n"), (80, b"\n"), (1, b"\n")] {
            let (fa, fai) = fasta(&recs, width, eol);
            let p = temp(&format!("w{width}.fa"), &fa);
            let i = temp(&format!("w{width}.fa.fai"), fai.as_bytes());
            let mut f = Indexed::with_index(&p, &i).unwrap();
            for (name, seq) in &recs {
                let len = seq.len() as u64;
                let mut edges: Vec<u64> = vec![0, 1, len - 1, len];
                for k in 1..6 {
                    edges.extend([k * width as u64 - 1, k * width as u64, k * width as u64 + 1]);
                }
                // Bytes at page edges: positions whose byte lands on k * 64 KiB.
                for k in 1..=(fa.len() as u64 / PAGE) {
                    let byte = k * PAGE;
                    let e = f.entry(name).unwrap();
                    if byte >= e.offset {
                        let rel = byte - e.offset;
                        let p = rel / e.line_width * e.line_bases
                            + (rel % e.line_width).min(e.line_bases);
                        edges.extend([p.saturating_sub(1), p, p + 1]);
                    }
                }
                for &a in &edges {
                    // Each edge to the next few, across a line, across a page.
                    let w = width as u64;
                    for b in [
                        a,
                        a + 1,
                        a + 2,
                        a + w,
                        a + w + 1,
                        a + 3 * w - 1,
                        a + PAGE + 3,
                        len,
                    ] {
                        if a > len || b < a {
                            continue;
                        }
                        let want = seq[a as usize..(b.min(len)) as usize].to_vec();
                        assert_eq!(
                            f.fetch(name, a, b).unwrap().unwrap(),
                            want,
                            "{name} [{a},{b}) w{width}"
                        );
                    }
                    if a < len {
                        assert_eq!(f.base(name, a).unwrap(), Some(seq[a as usize]));
                    }
                }
                // Past the end: nothing, as samtools faidx prints nothing.
                assert_eq!(f.fetch(name, len, len + 10).unwrap().unwrap(), b"");
                assert_eq!(f.base(name, len).unwrap(), None);
                assert_eq!(f.fetch(name, 0, u64::MAX).unwrap().unwrap(), *seq);
            }
            assert_eq!(f.fetch("chr9", 0, 10).unwrap(), None);
            assert_eq!(f.names(), ["chr1", "chr2", "chrM"]);
            let _ = std::fs::remove_file(p);
            let _ = std::fs::remove_file(i);
        }
    }

    #[test]
    fn many_pages_through_a_small_cache() {
        // 5 MB of sequence: 80 pages, more than the cache holds, read
        // forwards, backwards and at random.
        let seq = bases(9, 5_000_000);
        let (fa, fai) = fasta(&[("s", seq.clone())], 70, b"\n");
        let p = temp("big.fa", &fa);
        let i = temp("big.fa.fai", fai.as_bytes());
        let mut f = Indexed::with_index(&p, &i).unwrap();
        let mut x = 77u64;
        let mut order: Vec<u64> = (0..5_000_000).step_by(9_973).collect();
        order.extend(order.clone().iter().rev());
        for _ in 0..2_000 {
            x = x.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            order.push((x >> 33) % 5_000_000);
        }
        for p in order {
            assert_eq!(f.base("s", p).unwrap(), Some(seq[p as usize]), "{p}");
        }
        assert!(f.pages.len() <= PAGES);
        let _ = std::fs::remove_file(p);
        let _ = std::fs::remove_file(i);
    }

    #[test]
    fn refuses_compressed_and_bad_or_stale_indexes() {
        let gz = temp("x.fa.gz", &[0x1f, 0x8b, 8, 0]);
        let fai = temp("x.fa.gz.fai", b"s\t4\t3\t4\t5\n");
        let err = Indexed::with_index(&gz, &fai).err().unwrap().to_string();
        assert!(
            err.ends_with("is compressed; random access needs an uncompressed FASTA"),
            "{err}"
        );
        for (fai, want) in [
            ("s\t4\t3\n", ".fai line 1: 3 columns, want 5"),
            (
                "s\t4\t3\t4\t5\nt\tx\t3\t4\t5\n",
                ".fai line 2: column 2 'x' is not a whole number",
            ),
            (
                "s\t4\t3\t0\t5\n",
                ".fai line 1: 0 bases per line in 5 bytes",
            ),
            (
                "s\t4\t3\t5\t4\n",
                ".fai line 1: 5 bases per line in 4 bytes",
            ),
        ] {
            let e = read_fai(fai.as_bytes()).err().unwrap().to_string();
            assert_eq!(e, want);
        }
        // An index whose line width is off by one: the newline lands
        // inside the region.
        let fa = temp("s.fa", b">s\nACGT\nACGT\n");
        let stale = temp("s.fa.fai", b"s\t8\t3\t5\t5\n");
        let mut f = Indexed::with_index(&fa, &stale).unwrap();
        let err = f.fetch("s", 0, 8).err().unwrap().to_string();
        assert_eq!(
            err,
            "the .fai does not match the FASTA: s position 5 is not a base"
        );
        let past = temp("p.fa.fai", b"s\t800\t3\t4\t5\n");
        let mut f = Indexed::with_index(&fa, &past).unwrap();
        assert!(f.fetch("s", 700, 800).is_err());
        for p in [gz, fai, fa, stale, past] {
            let _ = std::fs::remove_file(p);
        }
    }
}
