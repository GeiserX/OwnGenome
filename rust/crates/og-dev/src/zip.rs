//! The members of a zip file, by name: only what `zip-content` compares.
//! Stored and deflated members; ZIP64 and encryption are refused.

use flate2::read::DeflateDecoder;
use std::collections::BTreeMap;
use std::io::Read;

const EOCD: u32 = 0x0605_4b50;
const CENTRAL: u32 = 0x0201_4b50;
const LOCAL: u32 = 0x0403_4b50;

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

fn usize_at(b: &[u8], at: usize) -> Result<usize, String> {
    usize::try_from(u32_at(b, at)?).map_err(|e| e.to_string())
}

/// Every member's name and decompressed content. A name that appears twice
/// is an error: which one a reader takes depends on the reader.
pub fn members(b: &[u8]) -> Result<BTreeMap<String, Vec<u8>>, String> {
    // The end-of-central-directory record: 22 bytes plus a comment of up to
    // 65535 bytes, at the end of the file.
    let lo = b.len().saturating_sub(22 + 0xffff);
    let eocd = (lo..=b.len().saturating_sub(22))
        .rev()
        .find(|&i| u32_at(b, i).ok() == Some(EOCD))
        .ok_or("not a zip file (no end-of-central-directory record)")?;
    let count = u16_at(b, eocd + 10)?;
    let offset = u32_at(b, eocd + 16)?;
    if count == 0xffff || offset == 0xffff_ffff {
        return Err("a ZIP64 file, which zip-content does not read".to_owned());
    }
    let mut at = usize::try_from(offset).map_err(|e| e.to_string())?;
    let mut out = BTreeMap::new();
    for _ in 0..count {
        if u32_at(b, at)? != CENTRAL {
            return Err(format!("no central directory entry at byte {at}"));
        }
        let flags = u16_at(b, at + 8)?;
        let method = u16_at(b, at + 10)?;
        let size = usize_at(b, at + 20)?;
        let name_len = usize::from(u16_at(b, at + 28)?);
        let extra_len = usize::from(u16_at(b, at + 30)?);
        let comment_len = usize::from(u16_at(b, at + 32)?);
        let local = usize_at(b, at + 42)?;
        let name = b
            .get(at + 46..at + 46 + name_len)
            .ok_or("a member name runs past the end of the file")?;
        let name = String::from_utf8_lossy(name).into_owned();
        if flags & 1 != 0 {
            return Err(format!("member {name} is encrypted"));
        }
        if u32_at(b, local)? != LOCAL {
            return Err(format!("member {name}: no local header at byte {local}"));
        }
        let data =
            local + 30 + usize::from(u16_at(b, local + 26)?) + usize::from(u16_at(b, local + 28)?);
        let raw = b
            .get(data..data + size)
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
            m => {
                return Err(format!(
                    "member {name}: compression method {m}, which zip-content does not read"
                ));
            }
        };
        if out.insert(name.clone(), content).is_some() {
            return Err(format!("member {name} appears twice"));
        }
        at += 46 + name_len + extra_len + comment_len;
    }
    Ok(out)
}

/// A zip file of `files`, deflated when `deflate`, with `mtime` (DOS time)
/// on every member: for the self-test and the unit tests.
pub fn write(files: &[(&str, &[u8])], deflate: bool, mtime: u16) -> Vec<u8> {
    use flate2::{Compression, write::DeflateEncoder};
    use std::io::Write;
    let mut out = Vec::new();
    let mut central = Vec::new();
    for (name, content) in files {
        let data = if deflate {
            let mut e = DeflateEncoder::new(Vec::new(), Compression::best());
            e.write_all(content).expect("in memory");
            e.finish().expect("in memory")
        } else {
            content.to_vec()
        };
        let crc = {
            let mut h = flate2::Crc::new();
            h.update(content);
            h.sum()
        };
        let method: u16 = if deflate { 8 } else { 0 };
        let size = |n: usize| u32::try_from(n).expect("a small test file");
        let offset = size(out.len());
        let name_len = u16::try_from(name.len()).expect("a short name");
        let mut common = Vec::new();
        common.extend_from_slice(&20u16.to_le_bytes()); // version needed
        common.extend_from_slice(&0u16.to_le_bytes()); // flags
        common.extend_from_slice(&method.to_le_bytes());
        common.extend_from_slice(&mtime.to_le_bytes());
        common.extend_from_slice(&0x5a21u16.to_le_bytes()); // date
        common.extend_from_slice(&crc.to_le_bytes());
        common.extend_from_slice(&size(data.len()).to_le_bytes());
        common.extend_from_slice(&size(content.len()).to_le_bytes());
        common.extend_from_slice(&name_len.to_le_bytes());
        common.extend_from_slice(&0u16.to_le_bytes()); // extra length
        out.extend_from_slice(&LOCAL.to_le_bytes());
        out.extend_from_slice(&common);
        out.extend_from_slice(name.as_bytes());
        out.extend_from_slice(&data);
        central.extend_from_slice(&CENTRAL.to_le_bytes());
        central.extend_from_slice(&20u16.to_le_bytes()); // version made by
        central.extend_from_slice(&common);
        central.extend_from_slice(&[0; 6]); // comment length, disk, internal attributes
        central.extend_from_slice(&0u32.to_le_bytes()); // external attributes
        central.extend_from_slice(&offset.to_le_bytes());
        central.extend_from_slice(name.as_bytes());
    }
    let cd_offset = u32::try_from(out.len()).expect("small");
    let cd_size = u32::try_from(central.len()).expect("small");
    let n = u16::try_from(files.len()).expect("few members");
    out.extend_from_slice(&central);
    out.extend_from_slice(&EOCD.to_le_bytes());
    out.extend_from_slice(&[0; 4]);
    out.extend_from_slice(&n.to_le_bytes());
    out.extend_from_slice(&n.to_le_bytes());
    out.extend_from_slice(&cd_size.to_le_bytes());
    out.extend_from_slice(&cd_offset.to_le_bytes());
    out.extend_from_slice(&0u16.to_le_bytes());
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_back_stored_and_deflated_members() {
        let files: &[(&str, &[u8])] = &[("a.txt", b"hello hello hello"), ("d/b.tsv", b"x\ty\n")];
        for deflate in [false, true] {
            let m = members(&write(files, deflate, 7)).unwrap();
            assert_eq!(m.len(), 2);
            assert_eq!(m["a.txt"], b"hello hello hello");
            assert_eq!(m["d/b.tsv"], b"x\ty\n");
        }
    }

    #[test]
    fn refuses_what_is_not_a_zip() {
        assert!(
            members(b"plain text")
                .unwrap_err()
                .contains("not a zip file")
        );
        let twice: &[(&str, &[u8])] = &[("a", b"1"), ("a", b"2")];
        assert!(
            members(&write(twice, false, 0))
                .unwrap_err()
                .contains("appears twice")
        );
    }
}
