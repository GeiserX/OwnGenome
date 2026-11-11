//! One member of a zip file (a jar), as Python's `zipfile` reads it: found
//! by name in the central directory (the last entry of a repeated name),
//! stored or deflated, its CRC checked. Only the end of the file, the
//! central directory and the member are read, not the whole jar.

use flate2::Crc;
use flate2::read::DeflateDecoder;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};

const EOCD: u32 = 0x0605_4b50;
const EOCD64: u32 = 0x0606_4b50;
const LOCATOR64: u32 = 0x0706_4b50;
const CENTRAL: u32 = 0x0201_4b50;
const LOCAL: u32 = 0x0403_4b50;

/// Why the member could not be read, worded as `zipfile` words it.
#[derive(Debug)]
pub enum ZipError {
    Io(io::Error),
    /// `BadZipFile` and its message.
    Bad(String),
    /// `KeyError`: no member of that name.
    Missing,
}

impl From<io::Error> for ZipError {
    fn from(e: io::Error) -> Self {
        ZipError::Io(e)
    }
}

fn bad(msg: &str) -> ZipError {
    ZipError::Bad(msg.to_owned())
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
        .map(|s| u64::from_le_bytes(s.try_into().unwrap_or_default()))
        .ok_or_else(|| bad("Truncated central directory"))
}

fn read_at(f: &mut File, at: u64, len: usize) -> Result<Vec<u8>, ZipError> {
    f.seek(SeekFrom::Start(at))?;
    let mut b = vec![0; len];
    f.read_exact(&mut b).map_err(|e| match e.kind() {
        io::ErrorKind::UnexpectedEof => bad("Truncated file header"),
        _ => ZipError::Io(e),
    })?;
    Ok(b)
}

/// A ZIP64 extra field value: the 64-bit one when the 32-bit is all ones.
fn zip64_extra(extra: &[u8], want: &[bool]) -> Result<Vec<Option<u64>>, ZipError> {
    let mut i = 0;
    while i + 4 <= extra.len() {
        let id = u16_at(extra, i)?;
        let len = usize::from(u16_at(extra, i + 2)?);
        if id == 1 {
            let mut at = i + 4;
            let mut out = Vec::new();
            for &w in want {
                if w {
                    out.push(Some(u64_at(extra, at)?));
                    at += 8;
                } else {
                    out.push(None);
                }
            }
            return Ok(out);
        }
        i += 4 + len;
    }
    Ok(want.iter().map(|_| None).collect())
}

/// The content of member `name` of the zip file at `path`.
pub fn member(path: &str, name: &str) -> Result<Vec<u8>, ZipError> {
    let mut f = File::open(path)?;
    let size = f.seek(SeekFrom::End(0))?;
    let tail_len = size.min(22 + 0xffff);
    let tail = read_at(
        &mut f,
        size - tail_len,
        usize::try_from(tail_len).unwrap_or(0),
    )?;
    let at = (0..=tail.len().saturating_sub(22))
        .rev()
        .find(|&i| u32_at(&tail, i).ok() == Some(EOCD))
        .ok_or_else(|| bad("File is not a zip file"))?;
    let mut count = u64::from(u16_at(&tail, at + 10)?);
    let mut cd_size = u64::from(u32_at(&tail, at + 12)?);
    let mut cd_offset = u64::from(u32_at(&tail, at + 16)?);
    if at >= 20 && u32_at(&tail, at - 20)? == LOCATOR64 {
        let eocd64 = u64_at(&tail, at - 12)?;
        let rec = read_at(&mut f, eocd64, 56)?;
        if u32_at(&rec, 0)? != EOCD64 {
            return Err(bad("Corrupt zip64 end of central directory locator"));
        }
        count = u64_at(&rec, 32)?;
        cd_size = u64_at(&rec, 40)?;
        cd_offset = u64_at(&rec, 48)?;
    }
    let cd = read_at(
        &mut f,
        cd_offset,
        usize::try_from(cd_size).map_err(|_| bad("Truncated central directory"))?,
    )?;
    let mut i = 0;
    let mut found = None;
    for _ in 0..count {
        if u32_at(&cd, i)? != CENTRAL {
            return Err(bad("Bad magic number for central directory"));
        }
        let name_len = usize::from(u16_at(&cd, i + 28)?);
        let extra_len = usize::from(u16_at(&cd, i + 30)?);
        let comment_len = usize::from(u16_at(&cd, i + 32)?);
        let this = cd
            .get(i + 46..i + 46 + name_len)
            .ok_or_else(|| bad("Truncated central directory"))?;
        if this == name.as_bytes() {
            // The last entry of a name wins, as in zipfile's NameToInfo.
            found = Some(i);
        }
        i += 46 + name_len + extra_len + comment_len;
    }
    let i = found.ok_or(ZipError::Missing)?;
    let flags = u16_at(&cd, i + 8)?;
    let method = u16_at(&cd, i + 10)?;
    let crc = u32_at(&cd, i + 16)?;
    let mut csize = u64::from(u32_at(&cd, i + 20)?);
    let mut usize_ = u64::from(u32_at(&cd, i + 24)?);
    let mut local = u64::from(u32_at(&cd, i + 42)?);
    let name_len = usize::from(u16_at(&cd, i + 28)?);
    let extra_len = usize::from(u16_at(&cd, i + 30)?);
    let extra = cd
        .get(i + 46 + name_len..i + 46 + name_len + extra_len)
        .unwrap_or_default();
    let big = [
        usize_ == 0xffff_ffff,
        csize == 0xffff_ffff,
        local == 0xffff_ffff,
    ];
    let wide = zip64_extra(extra, &big)?;
    usize_ = wide[0].unwrap_or(usize_);
    csize = wide[1].unwrap_or(csize);
    local = wide[2].unwrap_or(local);
    if flags & 1 != 0 {
        return Err(bad("the member is encrypted"));
    }
    let head = read_at(&mut f, local, 30)?;
    if u32_at(&head, 0)? != LOCAL {
        return Err(bad("Bad magic number for file header"));
    }
    let data = local + 30 + u64::from(u16_at(&head, 26)?) + u64::from(u16_at(&head, 28)?);
    let raw = read_at(
        &mut f,
        data,
        usize::try_from(csize).map_err(|_| bad("Truncated file header"))?,
    )?;
    let content = match method {
        0 => raw,
        8 => {
            let mut v = Vec::with_capacity(usize::try_from(usize_).unwrap_or(0));
            DeflateDecoder::new(raw.as_slice())
                .read_to_end(&mut v)
                .map_err(|e| ZipError::Bad(format!("Error -3 while decompressing data: {e}")))?;
            v
        }
        m => {
            return Err(ZipError::Bad(format!(
                "That compression method is not supported ({m})"
            )));
        }
    };
    let mut c = Crc::new();
    c.update(&content);
    if c.sum() != crc {
        return Err(ZipError::Bad(format!(
            "Bad CRC-32 for file {}",
            crate::py_repr(name)
        )));
    }
    Ok(content)
}
