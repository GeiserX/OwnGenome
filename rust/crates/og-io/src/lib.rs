//! Reading and writing files the way OwnGenome's Python tools do.
//!
//! - [`open`]: a path or `-` (stdin), plain or gzip, sniffed by its magic
//!   bytes. Multi-member gzip, BGZF included, is read to the end.
//! - [`Lines`]: Python's text mode with `newline=None`: `\n`, `\r\n` and a
//!   lone `\r` each end a line and read as `\n`; UTF-8 strict or `replace`.
//! - [`AtomicFile`]: written beside the target and renamed over it on
//!   [`AtomicFile::commit`], so a reader never sees half a file.

use flate2::bufread::MultiGzDecoder;
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
        Box::new(BufReader::with_capacity(1 << 16, MultiGzDecoder::new(buf)))
    } else {
        Box::new(buf)
    })
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
    pos: usize,
    utf8: Utf8,
}

impl<R: BufRead> Lines<R> {
    pub fn new(inner: R, utf8: Utf8) -> Self {
        Lines {
            inner,
            raw: Vec::new(),
            pos: 0,
            utf8,
        }
    }

    /// Python's `readline()`: the next line with its ending turned into
    /// `\n` (none on a last line without one), or `false` at the end.
    pub fn read_line(&mut self, out: &mut String) -> io::Result<bool> {
        out.clear();
        if self.pos == self.raw.len() {
            self.raw.clear();
            self.pos = 0;
            if self.inner.read_until(b'\n', &mut self.raw)? == 0 {
                return Ok(false);
            }
        }
        let rest = &self.raw[self.pos..];
        // read_until stops after a \n or at the end, so a \r\n is never split.
        let (body, used) = match rest.iter().position(|&c| c == b'\n' || c == b'\r') {
            None => (rest, rest.len()),
            Some(k) if rest[k] == b'\r' && rest.get(k + 1) == Some(&b'\n') => (&rest[..k], k + 2),
            Some(k) => (&rest[..k], k + 1),
        };
        let ended = used > body.len();
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
        self.pos += used;
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

    /// Flushes and renames over the target.
    pub fn commit(mut self) -> io::Result<()> {
        let out = self.out.take().expect("committed once");
        out.into_inner().map_err(io::IntoInnerError::into_error)?;
        fs::rename(&self.tmp, &self.target)
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
        let mut l = Lines::new(io::Cursor::new(b"\xff\n".to_vec()), Utf8::Strict);
        assert!(l.read_line(&mut String::new()).is_err());
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
        for (input, want) in [
            (b"plain\n".to_vec(), "plain\n"),
            (gz(b"zipped\n"), "zipped\n"),
            (two, "one\ntwo\n"),
            (vec![0x1f], "\u{1f}"),
        ] {
            let mut s = String::new();
            reader(io::Cursor::new(input))
                .unwrap()
                .read_to_string(&mut s)
                .unwrap();
            assert_eq!(s, want);
        }
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
        fs::remove_dir_all(&dir).unwrap();
    }
}
