//! A zip file read as Python's `zipfile` reads it: the central directory
//! (ZIP64 included) gives the members in order, names in UTF-8 when the
//! member says so and CP437 otherwise, and a member is streamed out through
//! its stored or deflated bytes with its CRC-32 and size checked at the end.
//! Nothing is held in memory but the directory and one 64 KiB buffer.

use flate2::Crc;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::Path;
use zlib_rs::{Inflate, InflateFlush, Status};

const EOCD: u32 = 0x0605_4b50;
const EOCD64: u32 = 0x0606_4b50;
const LOCATOR64: u32 = 0x0706_4b50;
const CENTRAL: u32 = 0x0201_4b50;
const LOCAL: u32 = 0x0403_4b50;

/// The upper half of code page 437, as Python's `cp437` codec maps it.
const CP437_HIGH: &str = "ÇüéâäàåçêëèïîìÄÅÉæÆôöòûùÿÖÜ¢£¥₧ƒáíóúñÑªº¿⌐¬½¼¡«»░▒▓│┤╡╢╖╕╣║╗╝╜╛┐└┴┬├─┼╞╟╚╔╩╦╠═╬╧╨╤╥╙╘╒╓╫╪┘┌█▄▌▐▀αßΓπΣσµτΦΘΩδ∞φε∩≡±≥≤⌠⌡÷≈°∙·√ⁿ²■\u{a0}";

/// Why a zip file could not be read.
#[derive(Debug)]
pub enum ZipError {
    Io(io::Error),
    /// The file is not a zip, or a member is damaged or of a kind not read.
    Bad(String),
}

impl From<io::Error> for ZipError {
    fn from(e: io::Error) -> Self {
        ZipError::Io(e)
    }
}

impl fmt::Display for ZipError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ZipError::Io(e) => write!(f, "{e}"),
            ZipError::Bad(m) => f.write_str(m),
        }
    }
}

fn bad(msg: impl Into<String>) -> ZipError {
    ZipError::Bad(msg.into())
}

/// One member as the central directory describes it.
#[derive(Debug, Clone)]
pub struct Entry {
    /// The name, decoded, cut at a NUL as `zipfile` cuts it.
    pub name: String,
    flags: u16,
    method: u16,
    crc: u32,
    csize: u64,
    /// The size once unpacked.
    pub size: u64,
    local: u64,
}

impl Entry {
    /// A directory entry: its name ends in `/`.
    pub fn is_dir(&self) -> bool {
        self.name.ends_with('/')
    }
}

/// A raw deflate stream, inflated with zlib-rs (the same inflate og-io's
/// BGZF reader uses; it takes about two thirds of the CPU time flate2's
/// default backend takes on Exomiser's multi-gigabyte members).
struct Inflater<R> {
    src: R,
    z: Inflate,
    done: bool,
}

impl<R: BufRead> Inflater<R> {
    fn new(src: R) -> Self {
        Self {
            src,
            z: Inflate::new(false, 15),
            done: false,
        }
    }
}

impl<R: BufRead> Read for Inflater<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        while !self.done && !out.is_empty() {
            let input = self.src.fill_buf()?;
            let eof = input.is_empty();
            let (in0, out0) = (self.z.total_in(), self.z.total_out());
            let status = self
                .z
                .decompress(input, out, InflateFlush::NoFlush)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.as_str()))?;
            let used = usize::try_from(self.z.total_in() - in0).expect("fits the buffer");
            let made = usize::try_from(self.z.total_out() - out0).expect("fits the buffer");
            self.src.consume(used);
            if status == Status::StreamEnd {
                self.done = true;
            }
            if made > 0 {
                return Ok(made);
            }
            if eof {
                // A cut stream: what came out is short, and the size and
                // CRC check after it says so.
                return Ok(0);
            }
        }
        Ok(0)
    }
}

/// An open zip file and its members, in the order of the central directory.
pub struct Archive {
    file: File,
    entries: Vec<Entry>,
}

fn u16_at(b: &[u8], at: usize) -> Result<u16, ZipError> {
    b.get(at..at + 2)
        .map(|s| u16::from_le_bytes([s[0], s[1]]))
        .ok_or_else(|| bad("Truncated central directory"))
}

fn u32_at(b: &[u8], at: usize) -> Result<u32, ZipError> {
    b.get(at..at + 4)
        .map(|s| u32::from_le_bytes([s[0], s[1], s[2], s[3]]))
        .ok_or_else(|| bad("Truncated central directory"))
}

fn u64_at(b: &[u8], at: usize) -> Result<u64, ZipError> {
    b.get(at..at + 8)
        .and_then(|s| s.try_into().ok())
        .map(u64::from_le_bytes)
        .ok_or_else(|| bad("Truncated central directory"))
}

fn read_at(f: &mut File, at: u64, len: u64) -> Result<Vec<u8>, ZipError> {
    let len = usize::try_from(len).map_err(|_| bad("Truncated central directory"))?;
    f.seek(SeekFrom::Start(at))?;
    let mut b = vec![0; len];
    f.read_exact(&mut b).map_err(|e| match e.kind() {
        io::ErrorKind::UnexpectedEof => bad("Truncated file header"),
        _ => ZipError::Io(e),
    })?;
    Ok(b)
}

/// A member name as `zipfile` decodes it: UTF-8 when flag bit 11 is set,
/// CP437 otherwise; cut at the first NUL.
fn decode_name(raw: &[u8], flags: u16) -> Result<String, ZipError> {
    let mut s = if flags & 0x800 != 0 {
        String::from_utf8(raw.to_vec())
            .map_err(|_| bad(format!("a member name is not UTF-8: {raw:?}")))?
    } else {
        raw.iter()
            .map(|&b| {
                if b < 0x80 {
                    char::from(b)
                } else {
                    CP437_HIGH
                        .chars()
                        .nth(usize::from(b - 0x80))
                        .unwrap_or('\u{fffd}')
                }
            })
            .collect()
    };
    if let Some(nul) = s.find('\0') {
        s.truncate(nul);
    }
    Ok(s)
}

/// The ZIP64 extra field's values for the fields stored as all ones, in
/// the order size, compressed size, local header offset.
fn zip64_extra(extra: &[u8], want: [bool; 3]) -> Result<[Option<u64>; 3], ZipError> {
    let mut out = [None; 3];
    let mut i = 0;
    while i + 4 <= extra.len() {
        let id = u16_at(extra, i)?;
        let len = usize::from(u16_at(extra, i + 2)?);
        if id == 1 {
            let mut at = i + 4;
            for (k, w) in want.iter().enumerate() {
                if *w {
                    out[k] = Some(u64_at(extra, at).map_err(|_| bad("Corrupt extra field"))?);
                    at += 8;
                }
            }
            return Ok(out);
        }
        i += 4 + len;
    }
    Ok(out)
}

impl Archive {
    /// Opens `path` and reads its central directory.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, ZipError> {
        let mut f = File::open(path)?;
        let size = f.seek(SeekFrom::End(0))?;
        let tail_len = size.min(22 + 0xffff);
        let tail = read_at(&mut f, size - tail_len, tail_len)?;
        let at = (0..=tail.len().saturating_sub(22))
            .rev()
            .find(|&i| u32_at(&tail, i).ok() == Some(EOCD))
            .ok_or_else(|| bad("File is not a zip file"))?;
        let mut cd_size = u64::from(u32_at(&tail, at + 12)?);
        let mut cd_offset = u64::from(u32_at(&tail, at + 16)?);
        if at >= 20 && u32_at(&tail, at - 20)? == LOCATOR64 {
            let rec = read_at(&mut f, u64_at(&tail, at - 12)?, 56)?;
            if u32_at(&rec, 0)? != EOCD64 {
                return Err(bad("Corrupt zip64 end of central directory locator"));
            }
            cd_size = u64_at(&rec, 40)?;
            cd_offset = u64_at(&rec, 48)?;
        }
        let cd = read_at(&mut f, cd_offset, cd_size)?;
        let mut entries = Vec::new();
        let mut i = 0;
        // As zipfile does: read entries until the directory's size is used up.
        while i < cd.len() {
            if u32_at(&cd, i)? != CENTRAL {
                return Err(bad("Bad magic number for central directory"));
            }
            let flags = u16_at(&cd, i + 8)?;
            let name_len = usize::from(u16_at(&cd, i + 28)?);
            let extra_len = usize::from(u16_at(&cd, i + 30)?);
            let comment_len = usize::from(u16_at(&cd, i + 32)?);
            let name = cd
                .get(i + 46..i + 46 + name_len)
                .ok_or_else(|| bad("Truncated central directory"))?;
            let extra = cd
                .get(i + 46 + name_len..i + 46 + name_len + extra_len)
                .ok_or_else(|| bad("Truncated central directory"))?;
            let size = u64::from(u32_at(&cd, i + 24)?);
            let csize = u64::from(u32_at(&cd, i + 20)?);
            let local = u64::from(u32_at(&cd, i + 42)?);
            let [w_size, w_csize, w_local] = zip64_extra(
                extra,
                [
                    size == 0xffff_ffff,
                    csize == 0xffff_ffff,
                    local == 0xffff_ffff,
                ],
            )?;
            entries.push(Entry {
                name: decode_name(name, flags)?,
                flags,
                method: u16_at(&cd, i + 10)?,
                crc: u32_at(&cd, i + 16)?,
                csize: w_csize.unwrap_or(csize),
                size: w_size.unwrap_or(size),
                local: w_local.unwrap_or(local),
            });
            i += 46 + name_len + extra_len + comment_len;
        }
        Ok(Archive { file: f, entries })
    }

    /// Every member, in the order of the central directory.
    pub fn entries(&self) -> &[Entry] {
        &self.entries
    }

    /// The member called `name`; the last one when the name repeats, as in
    /// `zipfile`'s name table.
    pub fn find(&self, name: &str) -> Option<&Entry> {
        self.entries.iter().rev().find(|e| e.name == name)
    }

    /// Streams the content of `e` into `out` and checks its CRC-32 and size.
    pub fn copy_to(&mut self, e: &Entry, out: &mut impl Write) -> Result<u64, ZipError> {
        if e.flags & 1 != 0 {
            return Err(bad(format!(
                "File {} is encrypted, password required for extraction",
                og_filters::py_repr(&e.name)
            )));
        }
        let head = read_at(&mut self.file, e.local, 30)?;
        if u32_at(&head, 0)? != LOCAL {
            return Err(bad("Bad magic number for file header"));
        }
        let name_len = u64::from(u16_at(&head, 26)?);
        let extra_len = u64::from(u16_at(&head, 28)?);
        let local_name = read_at(&mut self.file, e.local + 30, name_len)?;
        if decode_name(&local_name, e.flags)? != e.name {
            return Err(bad(format!(
                "File name in directory {} and header {} differ.",
                og_filters::py_repr(&e.name),
                og_filters::py_repr(&String::from_utf8_lossy(&local_name))
            )));
        }
        self.file
            .seek(SeekFrom::Start(e.local + 30 + name_len + extra_len))?;
        let raw = BufReader::with_capacity(1 << 16, (&self.file).take(e.csize));
        let mut src: Box<dyn Read + '_> = match e.method {
            0 => Box::new(raw),
            8 => Box::new(Inflater::new(raw)),
            m => {
                return Err(bad(format!(
                    "That compression method is not supported ({m}): {}",
                    og_filters::py_repr(&e.name)
                )));
            }
        };
        let mut crc = Crc::new();
        let mut n = 0u64;
        let mut buf = vec![0u8; 1 << 16];
        loop {
            let k = match src.read(&mut buf) {
                Ok(0) => break,
                Ok(k) => k,
                Err(err) if err.kind() == io::ErrorKind::Interrupted => continue,
                Err(err)
                    if matches!(
                        err.kind(),
                        io::ErrorKind::InvalidInput | io::ErrorKind::InvalidData
                    ) =>
                {
                    return Err(bad(format!(
                        "Error -3 while decompressing data: {err}: {}",
                        og_filters::py_repr(&e.name)
                    )));
                }
                Err(err) => return Err(err.into()),
            };
            crc.update(&buf[..k]);
            n += k as u64;
            out.write_all(&buf[..k])?;
        }
        if n != e.size || crc.sum() != e.crc {
            return Err(bad(format!(
                "Bad CRC-32 for file {}",
                og_filters::py_repr(&e.name)
            )));
        }
        Ok(n)
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::DeflateEncoder;

    /// A zip of (name, content, deflate) members, written by hand; `zip64`
    /// adds a ZIP64 end record and stores every offset there.
    pub(crate) fn build(members: &[(&str, &[u8], bool)], zip64: bool) -> Vec<u8> {
        let mut out = Vec::new();
        let mut cd = Vec::new();
        for (name, data, deflate) in members {
            let mut c = Crc::new();
            c.update(data);
            let body = if *deflate {
                let mut enc = DeflateEncoder::new(Vec::new(), Compression::default());
                enc.write_all(data).unwrap();
                enc.finish().unwrap()
            } else {
                data.to_vec()
            };
            let method: u16 = if *deflate { 8 } else { 0 };
            let offset = out.len() as u32;
            out.extend_from_slice(&LOCAL.to_le_bytes());
            out.extend_from_slice(&[20, 0, 0, 0]);
            out.extend_from_slice(&method.to_le_bytes());
            out.extend_from_slice(&[0; 4]);
            out.extend_from_slice(&c.sum().to_le_bytes());
            out.extend_from_slice(&(body.len() as u32).to_le_bytes());
            out.extend_from_slice(&(data.len() as u32).to_le_bytes());
            out.extend_from_slice(&(name.len() as u16).to_le_bytes());
            out.extend_from_slice(&[0, 0]);
            out.extend_from_slice(name.as_bytes());
            out.extend_from_slice(&body);
            cd.extend_from_slice(&CENTRAL.to_le_bytes());
            cd.extend_from_slice(&[20, 0, 20, 0, 0, 0]);
            cd.extend_from_slice(&method.to_le_bytes());
            cd.extend_from_slice(&[0; 4]);
            cd.extend_from_slice(&c.sum().to_le_bytes());
            cd.extend_from_slice(&(body.len() as u32).to_le_bytes());
            cd.extend_from_slice(&(data.len() as u32).to_le_bytes());
            cd.extend_from_slice(&(name.len() as u16).to_le_bytes());
            cd.extend_from_slice(&[0; 12]);
            cd.extend_from_slice(&offset.to_le_bytes());
            cd.extend_from_slice(name.as_bytes());
        }
        let cd_offset = out.len() as u64;
        out.extend_from_slice(&cd);
        if zip64 {
            let eocd64 = out.len() as u64;
            out.extend_from_slice(&EOCD64.to_le_bytes());
            out.extend_from_slice(&44u64.to_le_bytes());
            out.extend_from_slice(&[45, 0, 45, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
            out.extend_from_slice(&(members.len() as u64).to_le_bytes());
            out.extend_from_slice(&(members.len() as u64).to_le_bytes());
            out.extend_from_slice(&(cd.len() as u64).to_le_bytes());
            out.extend_from_slice(&cd_offset.to_le_bytes());
            out.extend_from_slice(&LOCATOR64.to_le_bytes());
            out.extend_from_slice(&[0; 4]);
            out.extend_from_slice(&eocd64.to_le_bytes());
            out.extend_from_slice(&1u32.to_le_bytes());
        }
        out.extend_from_slice(&EOCD.to_le_bytes());
        out.extend_from_slice(&[0; 4]);
        let n = if zip64 { 0xffff } else { members.len() as u16 };
        out.extend_from_slice(&n.to_le_bytes());
        out.extend_from_slice(&n.to_le_bytes());
        out.extend_from_slice(&(cd.len() as u32).to_le_bytes());
        let off = if zip64 { 0xffff_ffff } else { cd_offset as u32 };
        out.extend_from_slice(&off.to_le_bytes());
        out.extend_from_slice(&[0, 0]);
        out
    }

    fn write_tmp(bytes: &[u8], name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("og-setup-zip-{}-{name}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("a.zip");
        std::fs::write(&p, bytes).unwrap();
        p
    }

    #[test]
    fn reads_stored_and_deflated_members_and_the_last_of_a_name() {
        for zip64 in [false, true] {
            let z = build(
                &[
                    ("a.txt", b"first", false),
                    ("d/", b"", false),
                    ("a.txt", b"second, deflated", true),
                ],
                zip64,
            );
            let p = write_tmp(&z, &format!("basic{zip64}"));
            let mut a = Archive::open(&p).unwrap();
            let names: Vec<_> = a.entries().iter().map(|e| e.name.clone()).collect();
            assert_eq!(names, ["a.txt", "d/", "a.txt"]);
            assert!(a.entries()[1].is_dir());
            let e = a.find("a.txt").unwrap().clone();
            let mut out = Vec::new();
            a.copy_to(&e, &mut out).unwrap();
            assert_eq!(out, b"second, deflated");
        }
    }

    #[test]
    fn a_flipped_byte_fails_the_crc_check() {
        let mut z = build(&[("a.txt", b"some content", false)], false);
        z[30 + 5 + 2] ^= 1;
        let p = write_tmp(&z, "crc");
        let mut a = Archive::open(&p).unwrap();
        let e = a.entries()[0].clone();
        let err = a.copy_to(&e, &mut Vec::new()).unwrap_err();
        assert_eq!(err.to_string(), "Bad CRC-32 for file 'a.txt'");
    }

    #[test]
    fn not_a_zip_is_refused() {
        let p = write_tmp(b"plain text, no zip here", "notzip");
        assert_eq!(
            Archive::open(&p).err().unwrap().to_string(),
            "File is not a zip file"
        );
    }

    #[test]
    fn names_are_cp437_unless_flagged_utf8() {
        assert_eq!(decode_name(&[b'a', 0x81, 0xff], 0).unwrap(), "aü\u{a0}");
        assert_eq!(decode_name("é".as_bytes(), 0x800).unwrap(), "é");
        assert_eq!(decode_name(b"ab\0cd", 0).unwrap(), "ab");
        assert_eq!(CP437_HIGH.chars().count(), 128);
    }
}
