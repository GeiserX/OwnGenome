//! BGZF on many cores.
//!
//! BGZF (the format of `bgzip`, `.vcf.gz` and `.bam`) is gzip cut into
//! members of at most 64 KiB, each one saying how long it is. So the members
//! can be found without inflating them and inflated on several cores at
//! once, then handed back in order.
//!
//! - [`reader`] and [`open`] read what [`crate::open`] reads (plain, gzip,
//!   multi-member gzip, BGZF, NUL padding between members) and return the
//!   same bytes. A BGZF member is inflated by one of `workers` threads with
//!   zlib-rs; any other gzip member is inflated in order by the thread
//!   that reads the input, so a plain `gzip` file still reads correctly and
//!   still overlaps with the caller's work.
//! - [`Writer`] writes BGZF as `bgzip` lays it out (blocks of 65280 input
//!   bytes, the 28-byte end-of-file block), compressing blocks on `workers`
//!   threads. The bytes depend on the input only, never on `workers`.

use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::mem;
use std::num::NonZero;
use std::path::Path;
use std::sync::mpsc::{self, Receiver, SyncSender};
use std::sync::{Arc, Mutex, PoisonError};
use std::thread;

use flate2::bufread::GzDecoder;
use zlib_rs::{Deflate, DeflateFlush, Inflate, InflateFlush, Status};

/// The most a BGZF member holds, compressed or not.
const MAX_BLOCK: usize = 1 << 16;
/// Input bytes per written block, as `bgzip` and htslib cut them.
const BLOCK_INPUT: usize = 0xff00;
/// Header of a written block up to BSIZE: gzip, FEXTRA, MTIME 0, OS 255,
/// one extra subfield `BC` of two bytes.
const HEADER: [u8; 16] = [
    0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0,
];
/// The empty block that ends a BGZF file (SAM/BAM specification, 4.1.2).
pub const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0, 0x1b, 0, 3, 0, 0, 0, 0, 0, 0, 0,
    0, 0,
];
/// Raw deflate, as gzip carries it: no zlib header, a 32 KiB window.
const WINDOW_BITS: u8 = 15;
/// zlib's default level, the one `bgzip` uses.
const LEVEL: i32 = 6;
/// Decompressed bytes handed over at a time from a member that is not BGZF.
const CHUNK: usize = 1 << 16;

/// One result per member (or chunk), in input order.
type Slot = Receiver<io::Result<Vec<u8>>>;

/// The number of cores this process may use: a sensible `workers`.
pub fn workers() -> NonZero<usize> {
    thread::available_parallelism().unwrap_or(NonZero::<usize>::MIN)
}

/// Opens `path` (`-` is stdin) like [`crate::open`], inflating BGZF members
/// on `workers` threads.
pub fn open(
    path: impl AsRef<Path>,
    workers: NonZero<usize>,
) -> io::Result<Box<dyn BufRead + Send>> {
    let path = path.as_ref();
    if path == Path::new("-") {
        reader(io::stdin(), workers)
    } else {
        reader(File::open(path)?, workers)
    }
}

/// Wraps any reader like [`crate::reader`]: plain input is passed through,
/// gzip input is inflated, BGZF members on `workers` threads.
pub fn reader(
    mut inner: impl Read + Send + 'static,
    workers: NonZero<usize>,
) -> io::Result<Box<dyn BufRead + Send>> {
    // Read the magic byte by byte, as crate::reader does: a pipe may hand
    // over one byte at a time.
    let mut head = Vec::with_capacity(2);
    let mut byte = [0u8; 1];
    while head.len() < 2 {
        match inner.read(&mut byte) {
            Ok(0) => break,
            Ok(_) => head.push(byte[0]),
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    let gz = head == crate::GZIP_MAGIC;
    let input = io::Cursor::new(head).chain(inner);
    if !gz {
        return Ok(Box::new(BufReader::with_capacity(1 << 16, input)));
    }
    let n = workers.get();
    let jobs = pool(n, "og-bgzf-inflate", inflater, inflate)?;
    // Room for every worker's block plus as many again ready to be read.
    let (order_tx, order_rx) = mpsc::sync_channel(2 * n + 2);
    thread::Builder::new()
        .name("og-bgzf-read".into())
        .spawn(move || split(BufReader::with_capacity(1 << 17, input), &jobs, &order_tx))?;
    Ok(Box::new(Ordered {
        order: order_rx,
        buf: Vec::new(),
        pos: 0,
        done: false,
    }))
}

/// A BGZF member to inflate: its deflate data and gzip trailer.
struct InflateJob {
    body: Vec<u8>,
    reply: SyncSender<io::Result<Vec<u8>>>,
}

fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

fn inflater() -> Inflate {
    Inflate::new(false, WINDOW_BITS)
}

fn crc32(data: &[u8]) -> u32 {
    zlib_rs::crc32::crc32(0, data)
}

fn inflate(d: &mut Inflate, job: InflateJob) {
    let _ = job.reply.send(inflate_member(d, &job.body));
}

/// Inflates deflate data followed by the gzip trailer (CRC-32, ISIZE).
fn inflate_member(d: &mut Inflate, body: &[u8]) -> io::Result<Vec<u8>> {
    let (data, trailer) = body.split_at(body.len() - 8);
    let crc = u32::from_le_bytes(trailer[..4].try_into().expect("4 bytes"));
    let size = u32::from_le_bytes(trailer[4..].try_into().expect("4 bytes")) as usize;
    if size > MAX_BLOCK {
        return Err(invalid(format!(
            "BGZF block of {size} bytes, more than 65536"
        )));
    }
    let mut out = vec![0; size];
    d.reset(false);
    let status = d
        .decompress(data, &mut out, InflateFlush::Finish)
        .map_err(|e| invalid(format!("BGZF block: {}", e.as_str())))?;
    let n = d.total_out();
    if status != Status::StreamEnd || n != size as u64 {
        return Err(invalid(format!(
            "BGZF block does not inflate to the {size} bytes its trailer says"
        )));
    }
    if crc32(&out) != crc {
        return Err(invalid("BGZF block: CRC-32 does not match"));
    }
    Ok(out)
}

/// Starts `n` threads, each with its own `state`, running `run` on every job
/// sent to the returned sender. They stop when the sender is dropped.
fn pool<J, S>(
    n: usize,
    name: &str,
    state: fn() -> S,
    run: fn(&mut S, J),
) -> io::Result<SyncSender<J>>
where
    J: Send + 'static,
    S: 'static,
{
    let (tx, rx) = mpsc::sync_channel::<J>(n);
    let rx = Arc::new(Mutex::new(rx));
    for _ in 0..n {
        let rx = Arc::clone(&rx);
        thread::Builder::new().name(name.into()).spawn(move || {
            let mut s = state();
            loop {
                let job = rx.lock().unwrap_or_else(PoisonError::into_inner).recv();
                match job {
                    Ok(job) => run(&mut s, job),
                    Err(_) => return,
                }
            }
        })?;
    }
    Ok(tx)
}

/// Hands an already known result to the reader, in order. False once the
/// reader is gone.
fn ready(order: &SyncSender<Slot>, item: io::Result<Vec<u8>>) -> bool {
    let (tx, rx) = mpsc::sync_channel(1);
    let _ = tx.send(item);
    order.send(rx).is_ok()
}

/// The thread that reads the input: finds each member, sends BGZF members
/// to the workers and inflates other members itself, queueing one slot per
/// result in input order. An error ends the stream after the slots before it.
fn split<R: BufRead>(mut input: R, jobs: &SyncSender<InflateJob>, order: &SyncSender<Slot>) {
    loop {
        match member(&mut input, jobs, order) {
            Ok(true) => {}
            Ok(false) => return,
            Err(e) => {
                ready(order, Err(e));
                return;
            }
        }
    }
}

/// Reads one member. False at the end of the input or once the reader is gone.
fn member<R: BufRead>(
    input: &mut R,
    jobs: &SyncSender<InflateJob>,
    order: &SyncSender<Slot>,
) -> io::Result<bool> {
    // NUL bytes before a member are padding, as Python's gzip reads them.
    loop {
        let rest = input.fill_buf()?;
        if rest.is_empty() {
            return Ok(false);
        }
        let zeros = rest.iter().take_while(|&&c| c == 0).count();
        let more = zeros < rest.len();
        input.consume(zeros);
        if more {
            break;
        }
    }
    // A gzip member is at least 18 bytes, so 12 can always be read.
    let mut head = vec![0u8; 12];
    read_full(input, &mut head)?;
    if head[..2] != crate::GZIP_MAGIC {
        return Err(invalid(
            "not a gzip member after the end of the previous one",
        ));
    }
    let mut bsize = None;
    if head[2] == 8 && head[3] == 4 {
        let xlen = usize::from(u16::from_le_bytes([head[10], head[11]]));
        head.resize(12 + xlen, 0);
        read_full(input, &mut head[12..])?;
        bsize = bgzf_size(&head[12..]);
    }
    let Some(total) = bsize else {
        // Not BGZF: inflate the member here, replaying the bytes read so far.
        let mut dec = GzDecoder::new(io::Cursor::new(head).chain(&mut *input));
        loop {
            let mut chunk = vec![0; CHUNK];
            let n = read_full_or_end(&mut dec, &mut chunk)?;
            if n == 0 {
                return Ok(true);
            }
            chunk.truncate(n);
            if !ready(order, Ok(chunk)) {
                return Ok(false);
            }
        }
    };
    if total < head.len() + 8 + 2 {
        return Err(invalid(format!("BGZF block of {total} bytes is too short")));
    }
    let mut body = vec![0; total - head.len()];
    read_full(input, &mut body)?;
    let (reply, slot) = mpsc::sync_channel(1);
    if jobs.send(InflateJob { body, reply }).is_err() || order.send(slot).is_err() {
        return Ok(false);
    }
    Ok(true)
}

/// The total size of a BGZF member from its extra field: the `BC` subfield
/// of two bytes holds that size minus one. `None` when there is none.
fn bgzf_size(mut extra: &[u8]) -> Option<usize> {
    while extra.len() >= 4 {
        let len = usize::from(u16::from_le_bytes([extra[2], extra[3]]));
        let field = extra.get(4..4 + len)?;
        if extra[..2] == *b"BC" && len == 2 {
            return Some(usize::from(u16::from_le_bytes([field[0], field[1]])) + 1);
        }
        extra = &extra[4 + len..];
    }
    None
}

/// `read_exact`, with a truncated input reported as such.
fn read_full(r: &mut impl Read, buf: &mut [u8]) -> io::Result<()> {
    r.read_exact(buf).map_err(|e| {
        if e.kind() == io::ErrorKind::UnexpectedEof {
            io::Error::new(e.kind(), "the gzip input ends inside a member")
        } else {
            e
        }
    })
}

/// Fills `buf` unless the input ends first; returns the bytes read.
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

/// The caller's side: the slots in order, one decompressed buffer at a time.
struct Ordered {
    order: Receiver<Slot>,
    buf: Vec<u8>,
    pos: usize,
    done: bool,
}

impl Read for Ordered {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        let n = {
            let mut src = self.fill_buf()?;
            src.read(out)?
        };
        self.consume(n);
        Ok(n)
    }
}

impl BufRead for Ordered {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        while self.pos == self.buf.len() && !self.done {
            let Ok(slot) = self.order.recv() else {
                self.done = true;
                break;
            };
            match slot.recv() {
                Ok(Ok(buf)) => {
                    self.buf = buf;
                    self.pos = 0;
                }
                Ok(Err(e)) => {
                    self.done = true;
                    return Err(e);
                }
                Err(_) => {
                    self.done = true;
                    return Err(io::Error::other("a BGZF worker stopped"));
                }
            }
        }
        Ok(&self.buf[self.pos..])
    }

    fn consume(&mut self, amt: usize) {
        self.pos = (self.pos + amt).min(self.buf.len());
    }
}

/// A block of input to compress into one BGZF member.
struct DeflateJob {
    data: Vec<u8>,
    reply: SyncSender<io::Result<Vec<u8>>>,
}

/// Compressors at level 6 and at level 0, for a block that would not
/// shrink enough to fit.
fn compressors() -> (Deflate, Deflate) {
    (
        Deflate::new(LEVEL, false, WINDOW_BITS),
        Deflate::new(0, false, WINDOW_BITS),
    )
}

fn deflate(c: &mut (Deflate, Deflate), job: DeflateJob) {
    let _ = job.reply.send(deflate_block(c, &job.data));
}

/// Deflates all of `data` into `out`; the length written, or `None` when it
/// does not fit.
fn deflate_into(c: &mut Deflate, data: &[u8], out: &mut [u8]) -> io::Result<Option<usize>> {
    c.reset();
    let status = c
        .compress(data, out, DeflateFlush::Finish)
        .map_err(|e| io::Error::other(format!("BGZF deflate: {}", e.as_str())))?;
    Ok((status == Status::StreamEnd)
        .then(|| usize::try_from(c.total_out()).expect("at most 64 KiB")))
}

/// One BGZF member holding `data` (at most [`BLOCK_INPUT`] bytes).
fn deflate_block(c: &mut (Deflate, Deflate), data: &[u8]) -> io::Result<Vec<u8>> {
    let mut out = vec![0; MAX_BLOCK];
    let room = MAX_BLOCK - HEADER.len() - 2 - 8;
    let body = &mut out[HEADER.len() + 2..HEADER.len() + 2 + room];
    let n = match deflate_into(&mut c.0, data, body)? {
        Some(n) => n,
        // Stored blocks add 5 bytes per 65535: 65280 input bytes always fit.
        None => deflate_into(&mut c.1, data, body)?
            .ok_or_else(|| io::Error::other("a stored BGZF block does not fit"))?,
    };
    let total = HEADER.len() + 2 + n + 8;
    out.truncate(total);
    out[..HEADER.len()].copy_from_slice(&HEADER);
    let bsize = u16::try_from(total - 1).expect("a BGZF block fits in 64 KiB");
    out[HEADER.len()..HEADER.len() + 2].copy_from_slice(&bsize.to_le_bytes());
    out[total - 8..total - 4].copy_from_slice(&crc32(data).to_le_bytes());
    let size = u32::try_from(data.len()).expect("a BGZF block fits in 64 KiB");
    out[total - 4..].copy_from_slice(&size.to_le_bytes());
    Ok(out)
}

/// A BGZF writer that compresses blocks on several threads. Call
/// [`finish`](Writer::finish) to write the end-of-file block and get errors;
/// a writer dropped without it finishes the file and ignores errors.
pub struct Writer<W: Write> {
    inner: Option<W>,
    buf: Vec<u8>,
    pending: VecDeque<Slot>,
    jobs: SyncSender<DeflateJob>,
    max_pending: usize,
}

impl<W: Write> Writer<W> {
    pub fn new(inner: W, workers: NonZero<usize>) -> io::Result<Self> {
        let n = workers.get();
        Ok(Writer {
            inner: Some(inner),
            buf: Vec::with_capacity(BLOCK_INPUT),
            pending: VecDeque::new(),
            jobs: pool(n, "og-bgzf-deflate", compressors, deflate)?,
            max_pending: 2 * n + 2,
        })
    }

    /// Sends the buffered input off as one block.
    fn dispatch(&mut self) -> io::Result<()> {
        let data = mem::replace(&mut self.buf, Vec::with_capacity(BLOCK_INPUT));
        let (reply, slot) = mpsc::sync_channel(1);
        self.jobs
            .send(DeflateJob { data, reply })
            .map_err(|_| io::Error::other("a BGZF worker stopped"))?;
        self.pending.push_back(slot);
        while self.pending.len() >= self.max_pending {
            self.write_one()?;
        }
        Ok(())
    }

    /// Writes the oldest compressed block out.
    fn write_one(&mut self) -> io::Result<()> {
        let Some(slot) = self.pending.pop_front() else {
            return Ok(());
        };
        let block = slot
            .recv()
            .map_err(|_| io::Error::other("a BGZF worker stopped"))??;
        self.inner.as_mut().expect("not finished").write_all(&block)
    }

    /// Writes every block, then the end-of-file block, and returns the inner
    /// writer, flushed.
    pub fn finish(mut self) -> io::Result<W> {
        self.flush()?;
        let mut inner = self.inner.take().expect("finished once");
        inner.write_all(&EOF_BLOCK)?;
        inner.flush()?;
        Ok(inner)
    }
}

impl<W: Write> Write for Writer<W> {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        let n = data.len().min(BLOCK_INPUT - self.buf.len());
        self.buf.extend_from_slice(&data[..n]);
        if self.buf.len() == BLOCK_INPUT {
            self.dispatch()?;
        }
        Ok(n)
    }

    /// Ends the current block early (as htslib's `bgzf_flush` does) and
    /// writes every block out.
    fn flush(&mut self) -> io::Result<()> {
        if !self.buf.is_empty() {
            self.dispatch()?;
        }
        while !self.pending.is_empty() {
            self.write_one()?;
        }
        self.inner.as_mut().expect("not finished").flush()
    }
}

impl<W: Write> Drop for Writer<W> {
    fn drop(&mut self) {
        if self.inner.is_some()
            && self.flush().is_ok()
            && let Some(inner) = self.inner.as_mut()
        {
            let _ = inner.write_all(&EOF_BLOCK).and_then(|()| inner.flush());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::process::{Command, Stdio};

    fn n(k: usize) -> NonZero<usize> {
        NonZero::new(k).unwrap()
    }

    fn gz(data: &[u8]) -> Vec<u8> {
        let mut e = GzEncoder::new(Vec::new(), Compression::default());
        e.write_all(data).unwrap();
        e.finish().unwrap()
    }

    fn bgzf(data: &[u8], workers: usize) -> Vec<u8> {
        let mut w = Writer::new(Vec::new(), n(workers)).unwrap();
        w.write_all(data).unwrap();
        w.finish().unwrap()
    }

    /// Text that compresses, long enough for many blocks.
    fn text(len: usize) -> Vec<u8> {
        let mut out = Vec::with_capacity(len + 64);
        let mut i = 0u64;
        while out.len() < len {
            let x = i.wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 40;
            out.extend(
                format!("chr1\t{}\t.\tA\tG\t{}\tPASS\tDP={}\n", i * 7, x % 100, x).as_bytes(),
            );
            i += 1;
        }
        out.truncate(len);
        out
    }

    /// Bytes that do not compress, so a block falls back to stored deflate.
    fn noise(len: usize) -> Vec<u8> {
        let mut x = 0x2545_f491_4f6c_dd1du64;
        (0..len)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                x as u8
            })
            .collect()
    }

    fn read_par(input: Vec<u8>, workers: usize) -> io::Result<Vec<u8>> {
        let mut out = Vec::new();
        reader(io::Cursor::new(input), n(workers))?.read_to_end(&mut out)?;
        Ok(out)
    }

    fn read_seq(input: Vec<u8>) -> io::Result<Vec<u8>> {
        let mut out = Vec::new();
        crate::reader(io::Cursor::new(input))?.read_to_end(&mut out)?;
        Ok(out)
    }

    /// `gzip -dc` on `input`; `None` when gzip is not installed.
    fn gzip_dc(input: &[u8]) -> Option<(bool, Vec<u8>)> {
        let mut child = Command::new("gzip")
            .arg("-dc")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .ok()?;
        let mut stdin = child.stdin.take().unwrap();
        let data = input.to_vec();
        let feed = thread::spawn(move || {
            let _ = stdin.write_all(&data);
        });
        let out = child.wait_with_output().unwrap();
        feed.join().unwrap();
        Some((out.status.success(), out.stdout))
    }

    #[test]
    fn equals_gzip_dc_and_the_sequential_reader() {
        let big = text(1_000_000);
        let mut multi = gz(b"one\n");
        multi.extend(gz(b""));
        multi.extend(gz(b"two\n"));
        let mut mixed = bgzf(&big[..300_000], 3);
        mixed.extend(gz(&big[300_000..400_000]));
        mixed.extend(bgzf(&big[400_000..], 2));
        let mut two_files = bgzf(b"#header\n", 1);
        two_files.extend(bgzf(&big, 4));
        let cases: Vec<(&str, Vec<u8>, Vec<u8>)> = vec![
            ("plain gzip", gz(&big), big.clone()),
            ("multi-member gzip", multi, b"one\ntwo\n".to_vec()),
            ("BGZF", bgzf(&big, 4), big.clone()),
            ("BGZF of noise", bgzf(&noise(200_000), 2), noise(200_000)),
            ("BGZF, gzip, BGZF", mixed, big.clone()),
            (
                "two BGZF files one after the other",
                two_files,
                [&b"#header\n"[..], &big].concat(),
            ),
            ("empty BGZF", bgzf(b"", 1), Vec::new()),
            ("the EOF block alone", EOF_BLOCK.to_vec(), Vec::new()),
        ];
        for (name, input, want) in cases {
            assert_eq!(read_seq(input.clone()).unwrap(), want, "{name}: sequential");
            for w in [1, 2, 8] {
                assert_eq!(
                    read_par(input.clone(), w).unwrap(),
                    want,
                    "{name}: {w} workers"
                );
            }
            if let Some((ok, out)) = gzip_dc(&input) {
                assert!(ok, "{name}: gzip -dc failed");
                assert_eq!(out, want, "{name}: gzip -dc");
            }
        }
    }

    #[test]
    fn plain_input_and_padding() {
        for input in [b"plain\n".to_vec(), Vec::new(), vec![0x1f]] {
            assert_eq!(read_par(input.clone(), 4).unwrap(), input);
        }
        // NUL padding between and after members, as Python's gzip skips it.
        let mut padded = bgzf(b"one\n", 1);
        padded.extend([0u8; 16]);
        padded.extend(gz(b"two\n"));
        padded.extend([0u8; 40]);
        assert_eq!(read_par(padded.clone(), 2).unwrap(), b"one\ntwo\n");
        assert_eq!(read_seq(padded).unwrap(), b"one\ntwo\n");
    }

    #[test]
    fn damaged_input_is_an_error() {
        let good = bgzf(&text(200_000), 2);
        // Cut inside a block, inside a header, and garbage after a block.
        let mut cases = vec![
            ("cut inside a block", good[..good.len() / 2].to_vec()),
            (
                "cut inside a header",
                good[..good.len() - EOF_BLOCK.len() + 5].to_vec(),
            ),
            ("garbage after a block", [&good[..], b"junk"].concat()),
        ];
        let mut bad_crc = good.clone();
        let first_trailer = usize::from(u16::from_le_bytes([good[16], good[17]])) + 1 - 8;
        bad_crc[first_trailer] ^= 1;
        cases.push(("wrong CRC", bad_crc));
        let mut bad_data = good.clone();
        bad_data[30] ^= 0xff;
        cases.push(("damaged deflate data", bad_data));
        for (name, input) in cases {
            for w in [1, 4] {
                assert!(
                    read_par(input.clone(), w).is_err(),
                    "{name}: {w} workers read it"
                );
            }
            assert!(
                read_seq(input.clone()).is_err(),
                "{name}: the sequential reader read it"
            );
            if let Some((ok, _)) = gzip_dc(&input) {
                assert!(!ok, "{name}: gzip -dc read it");
            }
        }
    }

    #[test]
    fn writer_lays_out_blocks_like_bgzip() {
        let data = text(BLOCK_INPUT * 3 + 17);
        let one = bgzf(&data, 1);
        // The output depends on the input only.
        assert_eq!(bgzf(&data, 8), one);
        assert!(one.ends_with(&EOF_BLOCK));
        // Four data blocks: three full, one of 17 bytes, then EOF.
        let mut rest = &one[..];
        let mut sizes = Vec::new();
        while !rest.is_empty() {
            let total = bgzf_size(&rest[12..18]).unwrap();
            assert_eq!(rest[..16], HEADER);
            let isize = u32::from_le_bytes(rest[total - 4..total].try_into().unwrap());
            sizes.push(isize as usize);
            rest = &rest[total..];
        }
        assert_eq!(sizes, [BLOCK_INPUT, BLOCK_INPUT, BLOCK_INPUT, 17, 0]);
        // Small writes and a flush in the middle give the same content.
        let mut w = Writer::new(Vec::new(), n(3)).unwrap();
        for piece in data.chunks(1000) {
            w.write_all(piece).unwrap();
        }
        w.flush().unwrap();
        w.write_all(b"after a flush\n").unwrap();
        let out = w.finish().unwrap();
        assert_eq!(
            read_seq(out).unwrap(),
            [&data[..], b"after a flush\n"].concat()
        );
        // Dropped without finish: still a whole file.
        let mut buf = Vec::new();
        {
            let mut w = Writer::new(&mut buf, n(2)).unwrap();
            w.write_all(b"dropped\n").unwrap();
        }
        assert!(buf.ends_with(&EOF_BLOCK));
        assert_eq!(read_seq(buf).unwrap(), b"dropped\n");
    }

    #[test]
    fn stops_cleanly_when_dropped_early() {
        let input = bgzf(&text(2_000_000), 2);
        let mut r = reader(io::Cursor::new(input), n(4)).unwrap();
        let mut first = [0u8; 10];
        r.read_exact(&mut first).unwrap();
        drop(r);
    }
}
