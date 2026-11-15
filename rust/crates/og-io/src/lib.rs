//! Reading and writing files the way OwnGenome's Python tools do.
//!
//! - [`open`]: a path or `-` (stdin), plain or gzip, sniffed by its magic
//!   bytes. Multi-member gzip, BGZF included, is read to the end, and NUL
//!   padding after a member is skipped, as Python's `gzip` does.
//! - [`Lines`]: Python's text mode with `newline=None`: `\n`, `\r\n` and a
//!   lone `\r` each end a line and read as `\n`; UTF-8 strict or `replace`.
//! - [`AtomicFile`]: written beside the target and renamed over it on
//!   [`AtomicFile::commit`], so a reader never sees half a file.
//! - [`bgzf`]: the same reading with BGZF inflated on many cores, and a
//!   BGZF writer that compresses on many cores.
//! - [`tabix`]: the lines of a region of a bgzipped file, through its
//!   `.tbi` or `.csi` index, as `tabix FILE seq:beg-end` prints them.

pub mod bgzf;
pub mod tabix;

use flate2::bufread::GzDecoder;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Opens `path` for reading (`-` is stdin), decompressing gzip and BGZF.
pub fn open(path: impl AsRef<Path>) -> io::Result<Box<dyn BufRead>> {
    let path = path.as_ref();
    let inner: Box<dyn Read> = if path == Path::new("-") {
        Box::new(io::stdin().lock())
    } else {
        Box::new(File::open(path)?)
    };
    reader(inner)
}

/// Wraps any reader: gzip when it starts with the gzip magic, plain otherwise.
pub fn reader(mut inner: impl Read + 'static) -> io::Result<Box<dyn BufRead>> {
    // Read the magic byte by byte: a pipe may hand over one byte at a time.
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
    let gz = head == GZIP_MAGIC;
    let buf = BufReader::with_capacity(1 << 16, io::Cursor::new(head).chain(inner));
    Ok(if gz {
        Box::new(BufReader::with_capacity(
            1 << 16,
            Members {
                dec: Some(GzDecoder::new(buf)),
            },
        ))
    } else {
        Box::new(buf)
    })
}

/// Gzip members one after another, as Python's `gzip` reads them: NUL bytes
/// between or after members are padding, and the stream may end in them.
struct Members<R: BufRead> {
    /// The member being read; `None` once the input has ended.
    dec: Option<GzDecoder<R>>,
}

impl<R: BufRead> Read for Members<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        while let Some(dec) = self.dec.as_mut() {
            let n = dec.read(buf)?;
            if n > 0 {
                return Ok(n);
            }
            // This member has ended: skip padding, then the next one starts.
            let mut inner = self.dec.take().expect("a member").into_inner();
            loop {
                let rest = inner.fill_buf()?;
                if rest.is_empty() {
                    return Ok(0);
                }
                let zeros = rest.iter().take_while(|&&c| c == 0).count();
                let more = zeros < rest.len();
                inner.consume(zeros);
                if more {
                    break;
                }
            }
            self.dec = Some(GzDecoder::new(inner));
        }
        Ok(0)
    }
}

/// How invalid UTF-8 is read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Utf8 {
    /// `errors="strict"`: an error.
    Strict,
    /// `errors="replace"`: U+FFFD for each invalid sequence.
    Replace,
}

/// Python's universal-newline text reader over any [`BufRead`].
pub struct Lines<R> {
    inner: R,
    raw: Vec<u8>,
    skip_lf: bool,
    utf8: Utf8,
}

impl<R: BufRead> Lines<R> {
    pub fn new(inner: R, utf8: Utf8) -> Self {
        Lines {
            inner,
            raw: Vec::new(),
            skip_lf: false,
            utf8,
        }
    }

    /// Python's `readline()`: the next line with its ending turned into
    /// `\n` (none on a last line without one), or `false` at the end.
    pub fn read_line(&mut self, out: &mut String) -> io::Result<bool> {
        out.clear();
        self.raw.clear();
        // A CR ends this record immediately. Consume a following LF on the
        // next call, including when CR is the last byte of a buffer.
        if self.skip_lf {
            if self.inner.fill_buf()?.first() == Some(&b'\n') {
                self.inner.consume(1);
            }
            self.skip_lf = false;
        }
        let ended = loop {
            let bytes = self.inner.fill_buf()?;
            if bytes.is_empty() {
                if self.raw.is_empty() {
                    return Ok(false);
                }
                break false;
            }
            if let Some(k) = bytes.iter().position(|&b| b == b'\r' || b == b'\n') {
                self.raw.extend_from_slice(&bytes[..k]);
                self.skip_lf = bytes[k] == b'\r';
                self.inner.consume(k + 1);
                break true;
            }
            let n = bytes.len();
            self.raw.extend_from_slice(bytes);
            self.inner.consume(n);
        };
        // Consume the record before decoding so a strict UTF-8 error is
        // raised only once when a caller continues with the next record.
        let body = &self.raw;
        match self.utf8 {
            Utf8::Strict => out.push_str(
                std::str::from_utf8(body)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
            ),
            Utf8::Replace => out.push_str(&String::from_utf8_lossy(body)),
        }
        if ended {
            out.push('\n');
        }
        Ok(true)
    }
}

impl<R: BufRead> Iterator for Lines<R> {
    type Item = io::Result<String>;

    /// Each line without its ending, as `line.rstrip("\n")` after Python's
    /// universal-newline reading.
    fn next(&mut self) -> Option<Self::Item> {
        let mut s = String::new();
        match self.read_line(&mut s) {
            Ok(true) => {
                if s.ends_with('\n') {
                    s.pop();
                }
                Some(Ok(s))
            }
            Ok(false) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// A file written under a temporary name and renamed into place by
/// [`commit`](AtomicFile::commit). Dropped without a commit, the temporary
/// file is removed and the target is untouched.
pub struct AtomicFile {
    tmp: PathBuf,
    target: PathBuf,
    out: Option<BufWriter<File>>,
}

impl AtomicFile {
    pub fn create(target: impl AsRef<Path>) -> io::Result<Self> {
        let target = target.as_ref().to_path_buf();
        let name = target
            .file_name()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "no file name"))?;
        let mut tmp_name = std::ffi::OsString::from(".");
        tmp_name.push(name);
        tmp_name.push(format!(".tmp-{}", std::process::id()));
        let tmp = target.with_file_name(tmp_name);
        let out = BufWriter::new(File::create(&tmp)?);
        Ok(AtomicFile {
            tmp,
            target,
            out: Some(out),
        })
    }

    /// Flushes and renames over the target. On an error the temporary file
    /// is removed and the target is untouched.
    pub fn commit(mut self) -> io::Result<()> {
        let out = self.out.take().expect("committed once");
        let done = out
            .into_inner()
            .map_err(io::IntoInnerError::into_error)
            .and_then(|f| {
                drop(f);
                fs::rename(&self.tmp, &self.target)
            });
        if done.is_err() {
            let _ = fs::remove_file(&self.tmp);
        }
        done
    }
}

impl Write for AtomicFile {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.out.as_mut().expect("not committed").write(buf)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.out.as_mut().expect("not committed").flush()
    }
}

impl Drop for AtomicFile {
    fn drop(&mut self) {
        if self.out.take().is_some() {
            let _ = fs::remove_file(&self.tmp);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;

    fn lines(bytes: &[u8], utf8: Utf8) -> Vec<String> {
        let mut l = Lines::new(io::Cursor::new(bytes.to_vec()), utf8);
        let mut out = Vec::new();
        let mut s = String::new();
        while l.read_line(&mut s).unwrap() {
            out.push(s.clone());
        }
        out
    }

    #[test]
    fn universal_newlines() {
        assert_eq!(
            lines(b"a\r\nb\rc\n\rd", Utf8::Strict),
            ["a\n", "b\n", "c\n", "\n", "d"]
        );
        assert_eq!(lines(b"x\r", Utf8::Strict), ["x\n"]);
        assert_eq!(lines(b"", Utf8::Strict), Vec::<String>::new());
        assert_eq!(
            lines(b"\xff\xe2\x82\n", Utf8::Replace),
            ["\u{fffd}\u{fffd}\n"]
        );
        // A strict error on one line is raised once; reading goes on after it.
        let mut l = Lines::new(io::Cursor::new(b"\xff\nok\n".to_vec()), Utf8::Strict);
        let mut s = String::new();
        assert!(l.read_line(&mut s).is_err());
        assert!(l.read_line(&mut s).unwrap());
        assert_eq!(s, "ok\n");
        assert!(!l.read_line(&mut s).unwrap());
        let it: Vec<String> = Lines::new(io::Cursor::new(b"a\r\nb".to_vec()), Utf8::Strict)
            .map(Result::unwrap)
            .collect();
        assert_eq!(it, ["a", "b"]);
    }

    fn gz(data: &[u8]) -> Vec<u8> {
        let mut e = GzEncoder::new(Vec::new(), Compression::default());
        e.write_all(data).unwrap();
        e.finish().unwrap()
    }

    #[test]
    fn reads_plain_gzip_and_multi_member() {
        let mut two = gz(b"one\n");
        two.extend(gz(b"two\n"));
        // NUL padding between and after members, as Python's gzip skips it.
        let mut padded = gz(b"one\n");
        padded.extend([0u8; 16]);
        padded.extend(gz(b""));
        padded.extend(gz(b"two\n"));
        padded.extend([0u8; 12]);
        for (input, want) in [
            (b"plain\n".to_vec(), "plain\n"),
            (gz(b"zipped\n"), "zipped\n"),
            (two, "one\ntwo\n"),
            (padded, "one\ntwo\n"),
            (vec![0x1f], "\u{1f}"),
        ] {
            let mut s = String::new();
            reader(io::Cursor::new(input))
                .unwrap()
                .read_to_string(&mut s)
                .unwrap();
            assert_eq!(s, want);
        }
        // Anything else after a member is an error, as in Python.
        let mut junk = gz(b"one\n");
        junk.extend(b"\0\0junk");
        let mut s = String::new();
        let mut r = reader(io::Cursor::new(junk)).unwrap();
        assert!(r.read_to_string(&mut s).is_err());
    }

    #[test]
    fn atomic_write() {
        let dir = std::env::temp_dir().join(format!("og-io-test-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let target = dir.join("out.tsv");
        fs::write(&target, "old").unwrap();
        {
            let mut f = AtomicFile::create(&target).unwrap();
            f.write_all(b"half").unwrap();
        }
        assert_eq!(fs::read_to_string(&target).unwrap(), "old");
        let mut f = AtomicFile::create(&target).unwrap();
        f.write_all(b"new").unwrap();
        f.commit().unwrap();
        assert_eq!(fs::read_to_string(&target).unwrap(), "new");
        assert_eq!(fs::read_dir(&dir).unwrap().count(), 1);
        // A rename that fails (the target is a directory that is not empty)
        // leaves no temporary file behind.
        let busy = dir.join("busy");
        fs::create_dir_all(busy.join("x")).unwrap();
        let mut f = AtomicFile::create(&busy).unwrap();
        f.write_all(b"data").unwrap();
        assert!(f.commit().is_err());
        assert_eq!(fs::read_dir(&dir).unwrap().count(), 2);
        fs::remove_dir_all(&dir).unwrap();
    }
}

#[cfg(test)]
mod bounded_line_tests {
    use super::*;

    #[test]
    fn all_endings_and_buffer_edges() {
        for bytes in [
            b"a\nb\n".as_slice(),
            b"a\r\nb\r\n",
            b"a\rb\r",
            b"a\r\nb\rc\nd",
            b"tail",
            b"",
        ] {
            for capacity in 1..=8 {
                let mut input = Lines::new(
                    BufReader::with_capacity(capacity, io::Cursor::new(bytes)),
                    Utf8::Strict,
                );
                let mut got = Vec::new();
                let mut line = String::new();
                while input.read_line(&mut line).unwrap() {
                    got.push(line.clone());
                }
                let normalized = String::from_utf8(bytes.to_vec())
                    .unwrap()
                    .replace("\r\n", "\n")
                    .replace('\r', "\n");
                let want: Vec<_> = normalized.split_inclusive('\n').collect();
                assert_eq!(got, want, "capacity {capacity}: {bytes:?}");
            }
        }
    }

    #[test]
    fn cr_only_storage_is_bounded_by_one_line() {
        let input = b"a,b\r".repeat(1 << 18);
        let mut input = Lines::new(
            BufReader::with_capacity(64, io::Cursor::new(input)),
            Utf8::Strict,
        );
        let mut line = String::new();
        for _ in 0..1 << 18 {
            assert!(input.read_line(&mut line).unwrap());
            assert_eq!(line, "a,b\n");
            assert!(
                input.raw.capacity() <= 64,
                "buffer grew to {} bytes",
                input.raw.capacity()
            );
        }
        assert!(!input.read_line(&mut line).unwrap());
    }
}
