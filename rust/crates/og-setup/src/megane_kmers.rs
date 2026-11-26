//! `og setup megane-kmers --fa FASTA --out DIR [--prefix P] [--mem GB]
//! [--threads N]`: MEGAnE's k-mer set of a reference for step 50, the files
//! `DIR/P.mk` and `DIR/P.mi`, byte for byte what MEGAnE's own
//! `build_kmerset` writes for the same FASTA and `.fai`, in a few GB of
//! memory, where MEGAnE's own builder peaked at about 67 GB (RAM plus
//! swap) for GRCh38.
//!
//! The format (MEGAnE, RIKEN, MIT licence; this is a new implementation of
//! its file format and copies no code):
//!
//! - `.mk` has no header. It is an array of little-endian `u64` values in
//!   strictly increasing order. Each value is one 32-mer at 2 bits per
//!   base, the first base in bits 63 and 62, coded A=0 T=1 G=2 C=3.
//! - Each 32-base window of each contig adds its forward value and the
//!   value of its reverse complement, or one value when the two are equal
//!   (a palindrome). A value is kept when it was added 2 or more times.
//!   There are no canonical k-mers: the set holds both strands.
//! - `.mi` is the number of values in `.mk` in decimal, then `\n`.
//!
//! The input rules, each copied from MEGAnE's builder:
//!
//! - The contigs are the lines of `FASTA.fai`, in that order: a contig
//!   listed twice is counted twice. Each one is read from its offset
//!   (column 3) in chunks of `bytes per line` (column 5), of which the
//!   first `bytes per line - 1` are bases, until there are `length`
//!   (column 2) bases. Column 4 is not read.
//! - Lowercase is uppercase. `N` and `n` code as 0 and drop every window
//!   that holds them, with one exception: when the only N of a contig's
//!   first 32 bases is its first base, that first window is kept. Any other
//!   letter (IUPAC R, Y, K..., or a `-`) does not drop a window and codes
//!   as 0 in both values, so it reads as A on the forward strand and as T
//!   in the reverse complement.
//! - No window spans two contigs. A contig shorter than 32 bases adds
//!   nothing.
//!
//! Input MEGAnE's builder reads in an undefined way, or that is not the
//! FASTA its `.fai` describes, is refused with exit 2 and nothing written:
//! a carriage return in the sequence (Windows line ends, which MEGAnE would
//! read as bases), a byte of 0x80 or more (MEGAnE reads outside its table),
//! a line end in the sequence or a FASTA shorter than the `.fai` says (a
//! stale `.fai`), a gzip-compressed FASTA, and a `.fai` line that is not 5
//! whole numbers after the name or whose bytes per line is under 2.
//!
//! How it fits in a few GB. The genome is held at 4 bits per base (1.55 GB
//! for GRCh38). One pass counts the values by their top 20 bits. The 2^20
//! prefixes are then cut into ranges in order, each holding at most the
//! values one pass may hold under `--mem`. Each pass generates every value
//! again, keeps those in its range, sorts them and keeps the ones seen
//! twice or more. The ranges ascend, so the passes are written one after
//! the other and the file is sorted. `--threads` passes run at once. The
//! output is a function of the set of values alone, so it does not depend
//! on the thread count or on how the ranges were cut. A single 20-bit
//! prefix is never split: one that holds more values than a pass may hold
//! gets a pass of its own, over the budget.
//!
//! Both files are written under temporary names and renamed into place at
//! the end: a kill or an error leaves the old files or none.

use og_io::AtomicFile;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Barrier, Mutex};

/// The k-mer length; a value is one window at 2 bits per base.
pub const K: usize = 32;
/// The bits of a value the passes are cut on.
const PREFIX_BITS: u32 = 20;
const SHIFT: u32 = 64 - PREFIX_BITS;

/// A base's class, held at 4 bits: its forward code for A T G C, then N
/// (and n) and every other letter.
const A: u8 = 0;
const T: u8 = 1;
const G: u8 = 2;
const C: u8 = 3;
const N: u8 = 4;
const OTHER: u8 = 5;
/// The bytes that are refused, in the byte table only.
const BAD_CR: u8 = 0xFD;
const BAD_LF: u8 = 0xFE;
const BAD_HIGH: u8 = 0xFF;

/// The forward code and the reverse-complement code of each class: N and
/// other letters are 0 in both, as in MEGAnE's two tables.
const FWD: [u64; 6] = [0, 1, 2, 3, 0, 0];
const REV: [u64; 6] = [1, 0, 3, 2, 0, 0];

const fn class_table() -> [u8; 256] {
    let mut t = [OTHER; 256];
    let mut b = 0x80;
    while b < 256 {
        t[b] = BAD_HIGH;
        b += 1;
    }
    t[b'A' as usize] = A;
    t[b'a' as usize] = A;
    t[b'T' as usize] = T;
    t[b't' as usize] = T;
    t[b'G' as usize] = G;
    t[b'g' as usize] = G;
    t[b'C' as usize] = C;
    t[b'c' as usize] = C;
    t[b'N' as usize] = N;
    t[b'n' as usize] = N;
    t[b'\r' as usize] = BAD_CR;
    t[b'\n' as usize] = BAD_LF;
    t
}
const CLASS: [u8; 256] = class_table();

/// Why a build stopped: exit 1 for I/O, exit 2 for input or usage.
#[derive(Debug, PartialEq, Eq)]
pub enum Error {
    Io(String),
    Input(String),
}

impl Error {
    pub fn code(&self) -> u8 {
        match self {
            Error::Io(_) => 1,
            Error::Input(_) => 2,
        }
    }
    pub fn message(&self) -> &str {
        match self {
            Error::Io(m) | Error::Input(m) => m,
        }
    }
}

/// One line of the `.fai`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaiLine {
    pub name: String,
    pub len: u64,
    pub offset: u64,
    pub line_bytes: u64,
}

/// The `.fai` as MEGAnE's builder reads it: every line, in order, five
/// tab-separated fields. Columns 2, 3 and 5 must be whole numbers; column 4
/// is not used, but must be one too.
pub fn parse_fai(text: &str, path: &str) -> Result<Vec<FaiLine>, Error> {
    let mut out = Vec::new();
    for (n, line) in text.split_terminator('\n').enumerate() {
        let bad = |what: &str| Error::Input(format!("{path}: line {}: {what}", n + 1));
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 5 {
            return Err(bad(&format!(
                "{} tab-separated fields, a FASTA index line has 5 (index the FASTA again with samtools faidx)",
                f.len()
            )));
        }
        let num = |i: usize, what: &str| -> Result<u64, Error> {
            let s = f[i];
            if s.is_empty() || !s.bytes().all(|b| b.is_ascii_digit()) {
                return Err(bad(&format!("{what} {s:?} is not a whole number")));
            }
            s.parse::<u64>()
                .map_err(|_| bad(&format!("{what} {s:?} is too large")))
        };
        let len = num(1, "the length")?;
        let offset = num(2, "the offset")?;
        num(3, "the bases per line")?;
        let line_bytes = num(4, "the bytes per line")?;
        if line_bytes < 2 {
            return Err(bad(&format!(
                "{line_bytes} bytes per line: a line holds at least one base and its line end"
            )));
        }
        out.push(FaiLine {
            name: f[0].to_owned(),
            len,
            offset,
            line_bytes,
        });
    }
    Ok(out)
}

/// One contig at 4 bits per base, the first base in the low half of a byte.
pub struct Contig {
    pub name: String,
    pub len: usize,
    packed: Vec<u8>,
}

/// Every contig of the `.fai`, in its order.
pub struct Genome {
    pub contigs: Vec<Contig>,
}

impl Genome {
    /// The bytes the genome holds in memory for these `.fai` lines.
    pub fn bytes_for(fai: &[FaiLine]) -> u64 {
        fai.iter().map(|l| l.len.div_ceil(2)).sum()
    }

    /// Reads FASTA as MEGAnE's builder reads it through its `.fai`.
    pub fn load(fasta: &Path, fai: &[FaiLine]) -> Result<Genome, Error> {
        let fa = fasta.display().to_string();
        let io_err = |e: io::Error| Error::Io(format!("{fa}: {e}"));
        let mut file = File::open(fasta).map_err(io_err)?;
        let size = file.metadata().map_err(io_err)?.len();
        let mut magic = [0u8; 2];
        if size >= 2 {
            file.read_exact(&mut magic).map_err(io_err)?;
            if magic == [0x1f, 0x8b] {
                return Err(Error::Input(format!(
                    "{fa} is gzip-compressed: MEGAnE's k-mer set is built from a plain FASTA, unpack it first"
                )));
            }
        }
        let mut contigs = Vec::with_capacity(fai.len());
        let mut buf = vec![0u8; 1 << 22];
        for l in fai {
            contigs.push(read_contig(&mut file, size, l, &fa, &mut buf)?);
        }
        Ok(Genome { contigs })
    }
}

fn read_contig(
    file: &mut File,
    size: u64,
    l: &FaiLine,
    fa: &str,
    buf: &mut [u8],
) -> Result<Contig, Error> {
    let len = usize::try_from(l.len)
        .map_err(|_| Error::Input(format!("{fa}: contig {} is too long", l.name)))?;
    let mut packed = vec![0u8; len.div_ceil(2)];
    if len == 0 {
        return Ok(Contig {
            name: l.name.clone(),
            len,
            packed,
        });
    }
    let per = l.line_bytes - 1;
    let (full, rem) = (l.len / per, l.len % per);
    // The offset, from the contig's start, of the byte after its last base.
    let span = if rem > 0 {
        full * l.line_bytes + rem
    } else {
        (full - 1) * l.line_bytes + per
    };
    let stale = || {
        Error::Input(format!(
            "{fa}.fai does not describe {fa}: contig {} needs bytes up to offset {} and the file has {size}; index the FASTA again with samtools faidx",
            l.name,
            l.offset.saturating_add(span)
        ))
    };
    if l.offset.checked_add(span).is_none_or(|end| end > size) {
        return Err(stale());
    }
    file.seek(SeekFrom::Start(l.offset))
        .map_err(|e| Error::Io(format!("{fa}: {e}")))?;
    let (mut left, mut col, mut i) = (span, 0u64, 0usize);
    while left > 0 {
        let take = usize::try_from(left.min(buf.len() as u64)).unwrap_or(buf.len());
        file.read_exact(&mut buf[..take])
            .map_err(|e| Error::Io(format!("{fa}: {e}")))?;
        for &b in &buf[..take] {
            if col < per {
                let k = CLASS[b as usize];
                if k > OTHER {
                    return Err(refused(fa, &l.name, i, b, k));
                }
                packed[i >> 1] |= k << ((i & 1) << 2);
                i += 1;
            }
            col += 1;
            if col == l.line_bytes {
                col = 0;
            }
        }
        left -= take as u64;
    }
    debug_assert_eq!(i, len);
    Ok(Contig {
        name: l.name.clone(),
        len,
        packed,
    })
}

fn refused(fa: &str, name: &str, i: usize, b: u8, k: u8) -> Error {
    let at = format!("{fa}: contig {name}, base {}", i + 1);
    Error::Input(match k {
        BAD_CR => format!(
            "{at} is a carriage return: the FASTA has Windows (CRLF) line ends, which MEGAnE's builder reads as bases. Convert it to LF line ends and index it again with samtools faidx"
        ),
        BAD_LF => format!(
            "{at} is a line end: {fa}.fai does not describe this FASTA; index it again with samtools faidx"
        ),
        _ => format!(
            "{at} is the byte 0x{b:02x}, outside ASCII: MEGAnE's builder reads such a byte outside its table, so its k-mer set is not defined for this FASTA"
        ),
    })
}

/// Calls `emit` with every value contig `c` adds, in MEGAnE's order.
#[inline]
fn each_value(c: &Contig, mut emit: impl FnMut(u64)) {
    if c.len < K {
        return;
    }
    let (mut f, mut r, mut nn) = (0u64, 0u64, 0u32);
    let mut step = |i: usize, k: u8| {
        f = (f << 2) | FWD[k as usize];
        r = (r >> 2) | (REV[k as usize] << 62);
        let keep = if k == N {
            nn = K as u32 - 1;
            false
        } else if nn > 0 {
            nn -= 1;
            // In the first window the count of bases since the last N is
            // checked after this base: an N at base 0 alone keeps it.
            i == K - 1 && nn == 0
        } else {
            i >= K - 1
        };
        if keep {
            emit(f);
            if r != f {
                emit(r);
            }
        }
    };
    let pairs = c.len / 2;
    for (j, &b) in c.packed[..pairs].iter().enumerate() {
        step(2 * j, b & 15);
        step(2 * j + 1, b >> 4);
    }
    if c.len % 2 == 1 {
        step(c.len - 1, c.packed[pairs] & 15);
    }
}

/// How the work is split: threads, and the most values one pass holds.
#[derive(Debug, Clone, Copy)]
pub struct Plan {
    pub threads: usize,
    pub pass_values: u64,
}

/// What a build did.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Summary {
    /// Values added (forward and reverse complement of every window kept).
    pub added: u64,
    /// Values written to `.mk`.
    pub kept: u64,
    pub passes: usize,
}

/// Runs `job(t)` for each thread index t on `threads` threads and
/// collects the answers in index order.
fn on_threads<R: Send>(threads: usize, job: impl Fn(usize) -> R + Sync) -> Vec<R> {
    std::thread::scope(|s| {
        let hs: Vec<_> = (0..threads)
            .map(|t| {
                s.spawn({
                    let job = &job;
                    move || job(t)
                })
            })
            .collect();
        hs.into_iter()
            .map(|h| h.join().expect("a k-mer thread panicked"))
            .collect()
    })
}

/// The number of values under each 20-bit prefix.
fn histogram(g: &Genome, threads: usize) -> Vec<u64> {
    let mut order: Vec<usize> = (0..g.contigs.len()).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(g.contigs[i].len));
    let next = AtomicUsize::new(0);
    let parts = on_threads(threads, |_| {
        let mut h = vec![0u64; 1 << PREFIX_BITS];
        loop {
            let n = next.fetch_add(1, Ordering::Relaxed);
            let Some(&ci) = order.get(n) else { break };
            each_value(&g.contigs[ci], |v| h[(v >> SHIFT) as usize] += 1);
        }
        h
    });
    let mut total = vec![0u64; 1 << PREFIX_BITS];
    for h in parts {
        for (a, b) in total.iter_mut().zip(h) {
            *a += b;
        }
    }
    total
}

/// The prefix ranges `[lo, hi)` of each pass with their value counts, in
/// ascending order; prefixes that hold no value are in no pass.
fn cut(hist: &[u64], cap: u64) -> Vec<(u64, u64, u64)> {
    let mut out = Vec::new();
    let mut cur: Option<(u64, u64, u64)> = None;
    for (p, &n) in hist.iter().enumerate() {
        if n == 0 {
            continue;
        }
        let p = p as u64;
        cur = match cur {
            Some((lo, _, m)) if m + n <= cap => Some((lo, p + 1, m + n)),
            Some(done) => {
                out.push(done);
                Some((p, p + 1, n))
            }
            None => Some((p, p + 1, n)),
        };
    }
    out.extend(cur);
    out
}

/// One pass into `v`: every value under a prefix in `[lo, hi)`, sorted,
/// each value seen twice or more once. False when the pass found another
/// number of values than the histogram counted, which is a bug.
fn pass(g: &Genome, lo: u64, hi: u64, count: u64, v: &mut Vec<u64>) -> bool {
    v.clear();
    let width = hi - lo;
    for c in &g.contigs {
        each_value(c, |x| {
            if (x >> SHIFT).wrapping_sub(lo) < width {
                v.push(x);
            }
        });
    }
    if v.len() as u64 != count {
        return false;
    }
    v.sort_unstable();
    // Keep the first of each run of 2 or more, in place.
    let (mut w, mut i) = (0, 0);
    while i < v.len() {
        let mut j = i + 1;
        while j < v.len() && v[j] == v[i] {
            j += 1;
        }
        if j - i >= 2 {
            v[w] = v[i];
            w += 1;
        }
        i = j;
    }
    v.truncate(w);
    true
}

/// Writes the sorted repeated values of `g` to `mk` and returns the summary.
///
/// Thread t runs passes t, t + threads, and so on, in rounds of one pass a
/// thread; after each round the values are written in pass order. Each
/// thread has one buffer, sized once for its largest pass and used again
/// for every one, so the memory held is the plan's, whatever the allocator
/// does with freed blocks.
pub fn write_kmers(g: &Genome, plan: Plan, mk: &mut impl Write) -> io::Result<Summary> {
    let threads = plan.threads.max(1);
    let hist = histogram(g, threads);
    let added = hist.iter().sum();
    let passes = cut(&hist, plan.pass_values.max(1));
    drop(hist);
    let rounds = passes.len().div_ceil(threads);
    let bufs: Vec<Mutex<Vec<u64>>> = (0..threads).map(|_| Mutex::new(Vec::new())).collect();
    let (go, done) = (Barrier::new(threads + 1), Barrier::new(threads + 1));
    let wrong = AtomicBool::new(false);
    let mut kept = 0u64;
    let mut res = Ok(());
    std::thread::scope(|s| {
        for t in 0..threads {
            let (passes, bufs, go, done, wrong) = (&passes, &bufs, &go, &done, &wrong);
            s.spawn(move || {
                let most = passes.iter().skip(t).step_by(threads).map(|p| p.2).max();
                for r in 0..rounds {
                    go.wait();
                    if let Some(&(lo, hi, count)) = passes.get(r * threads + t) {
                        let mut v = bufs[t].lock().unwrap_or_else(|e| e.into_inner());
                        if v.capacity() == 0 {
                            v.reserve_exact(usize::try_from(most.unwrap_or(0)).unwrap_or(0));
                        }
                        if !pass(g, lo, hi, count, &mut v) {
                            wrong.store(true, Ordering::Relaxed);
                        }
                    }
                    done.wait();
                }
            });
        }
        let mut bytes = Vec::with_capacity(1 << 20);
        for r in 0..rounds {
            go.wait();
            done.wait();
            for (t, buf) in bufs.iter().enumerate() {
                if res.is_err() || r * threads + t >= passes.len() {
                    break;
                }
                let v = buf.lock().unwrap_or_else(|e| e.into_inner());
                kept += v.len() as u64;
                for chunk in v.chunks(1 << 17) {
                    bytes.clear();
                    for x in chunk {
                        bytes.extend_from_slice(&x.to_le_bytes());
                    }
                    if let Err(e) = mk.write_all(&bytes) {
                        res = Err(e);
                        break;
                    }
                }
            }
        }
    });
    res?;
    if wrong.load(Ordering::Relaxed) {
        return Err(io::Error::other(
            "a pass found another number of k-mers than the count before it (a bug in og)",
        ));
    }
    Ok(Summary {
        added,
        kept,
        passes: passes.len(),
    })
}

/// The memory a build may use for its pass buffers, out of `mem` bytes in
/// all, and the most values one pass holds; an error when `mem` is too
/// small for this genome.
pub fn plan_for(mem: u64, genome: u64, threads: usize) -> Result<Plan, String> {
    let threads = threads.max(1) as u64;
    // Each histogram thread holds 8 MiB of counts; 128 MiB for the rest.
    let fixed = genome + threads * (8 << 20) + (128 << 20);
    let floor = threads * (8 << 20);
    match mem.checked_sub(fixed) {
        Some(b) if b >= floor => Ok(Plan {
            threads: threads as usize,
            pass_values: b / 8 / threads,
        }),
        // The value as given (the command line rounds GB to whole bytes), and
        // the minimum rounded up, so --mem of the number shown fits.
        _ => Err(format!(
            "--mem {} is too little for this reference on {threads} threads: it needs at least {:.1} GB",
            mem as f64 / 1e9,
            ((fixed + floor) as f64 / 1e8).ceil() / 10.0
        )),
    }
}

/// Builds `out/prefix.mk` and `out/prefix.mi` from `fasta` and `fasta.fai`.
pub fn build(
    fasta: &Path,
    out: &Path,
    prefix: &str,
    mem: u64,
    threads: usize,
) -> Result<Summary, Error> {
    let mut fai_path = fasta.as_os_str().to_owned();
    fai_path.push(".fai");
    let fai_path = PathBuf::from(fai_path);
    let fai_name = fai_path.display().to_string();
    let text = std::fs::read(&fai_path).map_err(|e| Error::Io(format!("{fai_name}: {e}")))?;
    // MEGAnE reads a contig name as bytes, and og shows names only in its
    // messages: a name that is not UTF-8 is read, with U+FFFD for its bad
    // bytes. The number columns are checked as digits either way.
    let text = String::from_utf8_lossy(&text);
    let fai = parse_fai(&text, &fai_name)?;
    let plan = plan_for(mem, Genome::bytes_for(&fai), threads).map_err(Error::Input)?;
    let g = Genome::load(fasta, &fai)?;
    build_genome(&g, out, prefix, plan)
}

/// Writes `out/prefix.mk` and `out/prefix.mi` for a loaded genome.
pub fn build_genome(g: &Genome, out: &Path, prefix: &str, plan: Plan) -> Result<Summary, Error> {
    std::fs::create_dir_all(out).map_err(|e| Error::Io(format!("{}: {e}", out.display())))?;
    let mk_path = out.join(format!("{prefix}.mk"));
    let mi_path = out.join(format!("{prefix}.mi"));
    let fail = |p: &Path| {
        let p = p.display().to_string();
        move |e: io::Error| Error::Io(format!("{p}: {e}"))
    };
    let mut mk = AtomicFile::create(&mk_path).map_err(fail(&mk_path))?;
    let s = write_kmers(g, plan, &mut mk).map_err(fail(&mk_path))?;
    let mut mi = AtomicFile::create(&mi_path).map_err(fail(&mi_path))?;
    writeln!(mi, "{}", s.kept).map_err(fail(&mi_path))?;
    mk.commit().map_err(fail(&mk_path))?;
    mi.commit().map_err(fail(&mi_path))?;
    Ok(s)
}

const USAGE: &str =
    "usage: og setup megane-kmers --fa FASTA --out DIR [--prefix P] [--mem GB] [--threads N]";

/// The default memory budget, in GB: GRCh38 fits it on 8 threads.
pub const DEFAULT_MEM_GB: f64 = 5.0;

/// `--mem` in bytes, rounded to the nearest byte: 1.15 GB is
/// 1,150,000,000 bytes, where truncating would give 1,149,999,999, and the
/// too-little message shows the value as given.
fn mem_bytes(gb: f64) -> u64 {
    (gb * 1e9).round() as u64
}

/// The command line; returns the exit code: 0, 1 on an I/O error, 2 on
/// input it refuses or a wrong command line.
pub fn cli(args: &[String]) -> u8 {
    let usage = |m: &str| {
        eprintln!("{USAGE}\nog setup megane-kmers: {m}");
        2
    };
    let (mut fa, mut out, mut prefix) = (None, None, None);
    let mut mem = DEFAULT_MEM_GB;
    let mut threads = std::thread::available_parallelism()
        .map_or(1, usize::from)
        .min(8);
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let mut val = || it.next().cloned();
        match a.as_str() {
            "--fa" => fa = val(),
            "--out" => out = val(),
            "--prefix" => prefix = val(),
            "--mem" => match val().and_then(|v| v.parse::<f64>().ok()) {
                Some(m) if m.is_finite() && m > 0.0 => mem = m,
                _ => return usage("--mem takes a number of GB above 0"),
            },
            "--threads" => match val().and_then(|v| v.parse::<usize>().ok()) {
                Some(n) if n > 0 => threads = n,
                _ => return usage("--threads takes a whole number above 0"),
            },
            s => return usage(&format!("unknown argument {s}")),
        }
    }
    let (Some(fa), Some(out)) = (fa, out) else {
        return usage("--fa and --out are required");
    };
    // MEGAnE's own default prefix: the FASTA's file name.
    let prefix = prefix.unwrap_or_else(|| {
        Path::new(&fa)
            .file_name()
            .map_or_else(|| fa.clone(), |n| n.to_string_lossy().into_owned())
    });
    if prefix.is_empty() || prefix.contains('/') {
        return usage("--prefix is a file name, without a /");
    }
    let start = std::time::Instant::now();
    match build(
        Path::new(&fa),
        Path::new(&out),
        &prefix,
        mem_bytes(mem),
        threads,
    ) {
        Ok(s) => {
            println!(
                "og setup megane-kmers: {} repeated 32-mers of {} in {prefix}.mk ({} passes on {threads} threads, {:.0} s)",
                s.kept,
                s.added,
                s.passes,
                start.elapsed().as_secs_f64()
            );
            0
        }
        Err(e) => {
            eprintln!("og setup megane-kmers: {}", e.message());
            e.code()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn contig(seq: &[u8]) -> Contig {
        let mut packed = vec![0u8; seq.len().div_ceil(2)];
        for (i, &b) in seq.iter().enumerate() {
            packed[i >> 1] |= CLASS[b as usize] << ((i & 1) << 2);
        }
        Contig {
            name: "t".into(),
            len: seq.len(),
            packed,
        }
    }

    fn values(seq: &[u8]) -> Vec<u64> {
        let mut v = Vec::new();
        each_value(&contig(seq), |x| v.push(x));
        v
    }

    #[test]
    fn the_first_base_is_in_the_top_bits_and_a_is_0_t_1_g_2_c_3() {
        let s = [b"C".as_slice(), &[b'A'; 31]].concat();
        assert_eq!(values(&s), [3 << 62, 0x5555_5555_5555_5556]);
        let s = [&[b'A'; 31][..], b"G"].concat();
        // Forward AAA..AG = 2; its reverse complement C TTT..T.
        assert_eq!(values(&s), [2, (3 << 62) | 0x1555_5555_5555_5555]);
    }

    #[test]
    fn a_palindrome_adds_one_value() {
        // ACGT repeated is its own reverse complement.
        let s = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        assert_eq!(values(s).len(), 1);
    }

    #[test]
    fn n_drops_windows_but_an_n_at_base_0_alone_keeps_the_first() {
        let tail = [b'C'; 32];
        let n0 = [b"N".as_slice(), &tail].concat();
        assert_eq!(values(&n0).len(), 4, "both windows kept, two values each");
        let n1 = [b"CN".as_slice(), &tail].concat();
        assert_eq!(values(&n1).len(), 2, "only the window after the N");
        let n00 = [b"NN".as_slice(), &tail].concat();
        assert_eq!(
            values(&n00).len(),
            2,
            "the last N of the first window is base 1"
        );
    }

    #[test]
    fn fai_lines_are_checked() {
        let ok = parse_fai("a\t40\t3\t20\t21\n", "x.fai").unwrap();
        assert_eq!(ok[0].line_bytes, 21);
        assert!(parse_fai("a\t40\t3\t20\n", "x.fai").is_err());
        assert!(parse_fai("a\t40\t3\t20\t1\n", "x.fai").is_err());
        assert!(parse_fai("a\t4x\t3\t20\t21\n", "x.fai").is_err());
        assert!(parse_fai("a\t40\t3\t20\t21\r\n", "x.fai").is_err());
    }

    #[test]
    fn the_mem_message_names_the_value_as_given_and_a_minimum_that_fits() {
        let e = plan_for(1_150_000_000, 2_000_000_000, 8).unwrap_err();
        assert!(e.starts_with("--mem 1.15 is too little"), "{e}");
        // 2.27 GB is needed: 2.3 fits, 2.2 would not.
        assert!(e.ends_with("it needs at least 2.3 GB"), "{e}");
        // 0.27 GB of genome + 8 + 128 + 8 MiB = 0.421 GB: rounded up, since
        // --mem 0.4 would not fit.
        let e = plan_for(250_000_000, 270_000_000, 1).unwrap_err();
        assert!(e.ends_with("at least 0.5 GB"), "{e}");
    }

    #[test]
    fn mem_is_rounded_to_whole_bytes_not_truncated() {
        assert_eq!(mem_bytes(1.15), 1_150_000_000);
        assert_eq!(mem_bytes(0.25), 250_000_000);
        assert_eq!(mem_bytes(5.0), 5_000_000_000);
        let e = plan_for(mem_bytes(1.15), 2_000_000_000, 8).unwrap_err();
        assert!(e.starts_with("--mem 1.15 is"), "{e}");
    }

    #[test]
    fn passes_ascend_and_cover_every_counted_prefix() {
        let h = [0, 5, 0, 3, 9, 1];
        assert_eq!(cut(&h, 8), [(1, 4, 8), (4, 5, 9), (5, 6, 1)]);
        assert_eq!(cut(&h, 100), [(1, 6, 18)]);
        assert!(cut(&[0, 0], 1).is_empty());
    }

    #[test]
    fn the_memory_plan_splits_what_is_left_over_threads() {
        let p = plan_for(5_000_000_000, 1_550_000_000, 8).unwrap();
        assert_eq!(p.threads, 8);
        assert!(p.pass_values > 40_000_000, "{p:?}");
        assert!(plan_for(1_000_000_000, 1_550_000_000, 8).is_err());
    }
}
