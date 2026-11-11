//! The zip files of step 40: one deflated member each, stamped with a
//! date, so the same text on the same date gives the same bytes. The reader
//! takes what the strict reader needs: stored and deflated members, by name.

use std::fmt;
use std::io::{self, Read, Write};

use flate2::Crc;
use flate2::read::DeflateDecoder;
use zlib_rs::{Deflate, DeflateFlush, Status};

const LOCAL: u32 = 0x0403_4b50;
const CENTRAL: u32 = 0x0201_4b50;
const EOCD: u32 = 0x0605_4b50;
/// Version 2.0 (deflate), made on Unix, as Python's zipfile writes it.
const VERSION: u16 = 20;
const MADE_BY: u16 = (3 << 8) | VERSION;
/// A regular file, mode 0644, in the high half of the external attributes.
const EXTERNAL: u32 = 0o100_644 << 16;

/// A calendar date from 1980-01-01 to 2107-12-31, the range of a zip date.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DosDate {
    year: u16,
    month: u8,
    day: u8,
}

fn days_in_month(year: u16, month: u8) -> u8 {
    match month {
        2 if (year.is_multiple_of(4) && !year.is_multiple_of(100)) || year.is_multiple_of(400) => {
            29
        }
        2 => 28,
        4 | 6 | 9 | 11 => 30,
        _ => 31,
    }
}

impl DosDate {
    fn new(year: u16, month: u8, day: u8) -> Option<DosDate> {
        ((1980..=2107).contains(&year)
            && (1..=12).contains(&month)
            && (1..=days_in_month(year, month)).contains(&day))
        .then_some(DosDate { year, month, day })
    }

    /// `YYYY-MM-DD`.
    pub fn parse(s: &str) -> Option<DosDate> {
        let b = s.as_bytes();
        if b.len() != 10 || b[4] != b'-' || b[7] != b'-' {
            return None;
        }
        let num = |r: std::ops::Range<usize>| -> Option<u16> {
            let d = &s[r];
            d.bytes()
                .all(|c| c.is_ascii_digit())
                .then(|| d.parse().ok())?
        };
        let (y, m, d) = (num(0..4)?, num(5..7)?, num(8..10)?);
        DosDate::new(y, u8::try_from(m).ok()?, u8::try_from(d).ok()?)
    }

    /// The date `days` days after 1970-01-01 (the civil calendar, UTC).
    pub fn from_unix_days(days: u64) -> Option<DosDate> {
        // Howard Hinnant's civil_from_days, for days >= 0.
        let z = i64::try_from(days).ok()? + 719_468;
        let era = z / 146_097;
        let doe = z - era * 146_097;
        let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
        let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
        let mp = (5 * doy + 2) / 153;
        let day = doy - (153 * mp + 2) / 5 + 1;
        let month = if mp < 10 { mp + 3 } else { mp - 9 };
        let year = yoe + era * 400 + i64::from(month <= 2);
        DosDate::new(
            u16::try_from(year).ok()?,
            u8::try_from(month).ok()?,
            u8::try_from(day).ok()?,
        )
    }

    /// The MS-DOS date field; the time field is 0 (midnight).
    fn field(self) -> u16 {
        ((self.year - 1980) << 9) | (u16::from(self.month) << 5) | u16::from(self.day)
    }
}

impl fmt::Display for DosDate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:04}-{:02}-{:02}", self.year, self.month, self.day)
    }
}

fn too_big(what: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidInput,
        format!("{what} too large for a zip without ZIP64"),
    )
}

/// The text a [`Member`] collects before it hands it to the deflater.
const TEXT_CHUNK: usize = 1 << 16;

/// One zip member deflated as its text arrives, at zlib's default level:
/// only the compressed bytes and 64 KiB of text are held.
pub struct Member {
    d: Deflate,
    crc: Crc,
    /// The length of the text, which `Crc::amount` keeps only modulo 2^32.
    size: u64,
    text: Vec<u8>,
    packed: Vec<u8>,
    chunk: Vec<u8>,
}

impl Default for Member {
    fn default() -> Member {
        Member {
            d: Deflate::new(6, false, 15),
            crc: Crc::new(),
            size: 0,
            text: Vec::with_capacity(TEXT_CHUNK),
            packed: Vec::new(),
            chunk: vec![0; TEXT_CHUNK],
        }
    }
}

impl Member {
    /// Appends `data` to the member's text.
    pub fn write(&mut self, data: &[u8]) -> io::Result<()> {
        self.crc.update(data);
        self.size += data.len() as u64;
        if !self.text.is_empty() && self.text.len() + data.len() > TEXT_CHUNK {
            self.pump(DeflateFlush::NoFlush)?;
        }
        self.text.extend_from_slice(data);
        Ok(())
    }

    /// Deflates the text held; `Finish` also ends the stream.
    fn pump(&mut self, flush: DeflateFlush) -> io::Result<()> {
        let mut at = 0;
        loop {
            let (inp, done) = (self.d.total_in(), self.d.total_out());
            let status = self
                .d
                .compress(&self.text[at..], &mut self.chunk, flush)
                .map_err(|e| io::Error::other(format!("deflate: {}", e.as_str())))?;
            at += usize::try_from(self.d.total_in() - inp).map_err(|_| too_big("member"))?;
            let n = usize::try_from(self.d.total_out() - done).map_err(|_| too_big("member"))?;
            self.packed.extend_from_slice(&self.chunk[..n]);
            let finished = match flush {
                DeflateFlush::Finish => status == Status::StreamEnd,
                _ => at == self.text.len() && n < self.chunk.len(),
            };
            if finished {
                self.text.clear();
                return Ok(());
            }
        }
    }

    /// Writes the zip of this one member, named `name` and dated `date`.
    pub fn finish(mut self, out: &mut impl Write, name: &str, date: DosDate) -> io::Result<()> {
        self.pump(DeflateFlush::Finish)?;
        let packed = &self.packed;
        let csize = u32::try_from(packed.len()).map_err(|_| too_big("member"))?;
        let usize_ = u32::try_from(self.size).map_err(|_| too_big("member"))?;
        let crc = self.crc.sum();
        let name_len = u16::try_from(name.len()).map_err(|_| too_big("name"))?;
        let (time, day) = (0u16, date.field());

        let mut buf = Vec::with_capacity(2 * name.len() + 128);
        // Local header.
        buf.extend_from_slice(&LOCAL.to_le_bytes());
        for v in [VERSION, 0, 8, time, day] {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        for v in [crc, csize, usize_] {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf.extend_from_slice(&name_len.to_le_bytes());
        buf.extend_from_slice(&0u16.to_le_bytes());
        buf.extend_from_slice(name.as_bytes());
        // The member's data goes out as it is, never copied into the header buffer.
        out.write_all(&buf)?;
        out.write_all(packed)?;
        let central_at = u32::try_from(buf.len() + packed.len()).map_err(|_| too_big("zip"))?;
        buf.clear();
        // Central directory.
        buf.extend_from_slice(&CENTRAL.to_le_bytes());
        for v in [MADE_BY, VERSION, 0, 8, time, day] {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        for v in [crc, csize, usize_] {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        for v in [name_len, 0, 0, 0, 0] {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf.extend_from_slice(&EXTERNAL.to_le_bytes());
        buf.extend_from_slice(&0u32.to_le_bytes());
        buf.extend_from_slice(name.as_bytes());
        let central_len = u32::try_from(buf.len()).map_err(|_| too_big("zip"))?;
        // End of central directory.
        buf.extend_from_slice(&EOCD.to_le_bytes());
        for v in [0u16, 0, 1, 1] {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf.extend_from_slice(&central_len.to_le_bytes());
        buf.extend_from_slice(&central_at.to_le_bytes());
        buf.extend_from_slice(&0u16.to_le_bytes());
        out.write_all(&buf)
    }
}

/// Writes a zip of one member `name` holding `data`, deflated, dated `date`.
pub fn write_one(out: &mut impl Write, name: &str, data: &[u8], date: DosDate) -> io::Result<()> {
    let mut m = Member::default();
    m.write(data)?;
    m.finish(out, name, date)
}

fn u16_at(b: &[u8], at: usize) -> Result<u16, String> {
    b.get(at..at + 2)
        .map(|s| u16::from_le_bytes([s[0], s[1]]))
        .ok_or_else(|| "a zip structure runs past the end of the file".to_owned())
}

fn u32_at(b: &[u8], at: usize) -> Result<u32, String> {
    b.get(at..at + 4)
        .map(|s| u32::from_le_bytes([s[0], s[1], s[2], s[3]]))
        .ok_or_else(|| "a zip structure runs past the end of the file".to_owned())
}

fn at(v: u32) -> usize {
    usize::try_from(v).unwrap_or(usize::MAX)
}

/// Every member's name and content, in the order of the central directory.
/// Stored and deflated members are read and their CRC-32 checked; ZIP64,
/// encryption and other methods are refused.
pub fn read_all(b: &[u8]) -> Result<Vec<(String, Vec<u8>)>, String> {
    let lo = b.len().saturating_sub(22 + 0xffff);
    let eocd = (lo..=b.len().saturating_sub(22))
        .rev()
        .find(|&i| u32_at(b, i).ok() == Some(EOCD))
        .ok_or("no end-of-central-directory record")?;
    let count = u16_at(b, eocd + 10)?;
    let offset = u32_at(b, eocd + 16)?;
    if count == 0xffff || offset == 0xffff_ffff {
        return Err("a ZIP64 file".to_owned());
    }
    let mut p = at(offset);
    let mut out = Vec::with_capacity(usize::from(count));
    for _ in 0..count {
        if u32_at(b, p)? != CENTRAL {
            return Err(format!("no central directory entry at byte {p}"));
        }
        let flags = u16_at(b, p + 8)?;
        let method = u16_at(b, p + 10)?;
        let crc = u32_at(b, p + 16)?;
        let size = at(u32_at(b, p + 20)?);
        let name_len = usize::from(u16_at(b, p + 28)?);
        let extra_len = usize::from(u16_at(b, p + 30)?);
        let comment_len = usize::from(u16_at(b, p + 32)?);
        let local = at(u32_at(b, p + 42)?);
        let name = b
            .get(p + 46..p + 46 + name_len)
            .ok_or("a member name runs past the end of the file")?;
        let name = String::from_utf8_lossy(name).into_owned();
        if flags & 1 != 0 {
            return Err(format!("member {name} is encrypted"));
        }
        if u32_at(b, local)? != LOCAL {
            return Err(format!("member {name}: no local header at byte {local}"));
        }
        let data_at =
            local + 30 + usize::from(u16_at(b, local + 26)?) + usize::from(u16_at(b, local + 28)?);
        let raw = b
            .get(data_at..data_at.saturating_add(size))
            .ok_or_else(|| format!("member {name} runs past the end of the file"))?;
        let content = match method {
            0 => raw.to_vec(),
            8 => {
                let mut v = Vec::new();
                DeflateDecoder::new(raw)
                    .read_to_end(&mut v)
                    .map_err(|e| format!("member {name}: {e}"))?;
                v
            }
            m => return Err(format!("member {name}: compression method {m}")),
        };
        let mut c = Crc::new();
        c.update(&content);
        if c.sum() != crc {
            return Err(format!("member {name}: CRC-32 does not match"));
        }
        out.push((name, content));
        p += 46 + name_len + extra_len + comment_len;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dates() {
        assert_eq!(DosDate::from_unix_days(0), None);
        assert_eq!(
            DosDate::from_unix_days(3652).unwrap().to_string(),
            "1980-01-01"
        );
        assert_eq!(
            DosDate::from_unix_days(20402).unwrap().to_string(),
            "2025-11-10"
        );
        assert_eq!(
            DosDate::from_unix_days(11016).unwrap().to_string(),
            "2000-02-29"
        );
        assert_eq!(
            DosDate::parse("2024-02-29").unwrap().to_string(),
            "2024-02-29"
        );
        for bad in [
            "2023-02-29",
            "1979-12-31",
            "2108-01-01",
            "2026-1-01",
            "2026-01-0x",
            "+026-01-01",
        ] {
            assert_eq!(DosDate::parse(bad), None, "{bad}");
        }
        assert_eq!(
            DosDate::parse("2026-11-10").unwrap().field(),
            (46 << 9) | (11 << 5) | 10
        );
    }

    #[test]
    fn round_trip_and_reproducible() {
        let d = DosDate::parse("2026-11-10").unwrap();
        let text = b"# a\nrs1\t1\t100\tAG\n".repeat(50);
        let (mut a, mut b) = (Vec::new(), Vec::new());
        write_one(&mut a, "S_23andme_v5.txt", &text, d).unwrap();
        write_one(&mut b, "S_23andme_v5.txt", &text, d).unwrap();
        assert_eq!(a, b);
        assert_eq!(
            read_all(&a).unwrap(),
            [("S_23andme_v5.txt".to_owned(), text.clone())]
        );
        let mut c = Vec::new();
        write_one(
            &mut c,
            "S_23andme_v5.txt",
            &text,
            DosDate::parse("2026-11-11").unwrap(),
        )
        .unwrap();
        assert_ne!(a, c);
    }

    #[test]
    fn a_member_fed_line_by_line_holds_the_text() {
        let d = DosDate::parse("2026-11-10").unwrap();
        // Past TEXT_CHUNK several times, so the deflater runs mid-stream.
        let text: Vec<u8> = (0..40_000)
            .flat_map(|i| format!("rs{i}\t1\t{}\tAG\n", 1000 + 7 * i).into_bytes())
            .collect();
        let feed = || {
            let mut m = Member::default();
            for l in text.split_inclusive(|&b| b == b'\n') {
                m.write(l).unwrap();
            }
            let mut z = Vec::new();
            m.finish(&mut z, "x.txt", d).unwrap();
            z
        };
        let z = feed();
        assert_eq!(z, feed());
        // The local header's sizes: compressed, then the text's length.
        let csize = u32_at(&z, 18).unwrap() as usize;
        assert_eq!(z.len(), 30 + 5 + csize + 46 + 5 + 22);
        assert_eq!(u32_at(&z, 22).unwrap() as usize, text.len());
        assert_eq!(read_all(&z).unwrap(), [("x.txt".to_owned(), text)]);
    }

    #[test]
    fn a_damaged_member_is_refused() {
        let mut z = Vec::new();
        write_one(
            &mut z,
            "x.txt",
            b"hello hello hello\n",
            DosDate::parse("2026-01-01").unwrap(),
        )
        .unwrap();
        // Flip the stored CRC in the central directory.
        let central = z
            .windows(4)
            .position(|w| w == CENTRAL.to_le_bytes())
            .unwrap();
        z[central + 16] ^= 1;
        assert!(read_all(&z).unwrap_err().contains("CRC-32"));
        assert!(read_all(b"not a zip").is_err());
    }
}
