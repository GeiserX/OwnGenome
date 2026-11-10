//! The input stamp and the cache it keys.
//!
//! The stamp is a 64-bit FNV-1a over the format version, the og tree hash,
//! the stdin bytes and each input's path, size and modification time (or a
//! mark that it is absent, so a file that appears changes the stamp). It
//! detects change; it does not resist someone who wants a collision.
//!
//! The cache holds one document per run: `<cache-dir>/<run>/<stamp>.json.gz`.
//! Writing a new one removes the run's older stamps.

use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::time::UNIX_EPOCH;

const OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const PRIME: u64 = 0x0000_0100_0000_01b3;

/// FNV-1a, 64 bits.
#[derive(Debug, Clone, Copy)]
pub struct Fnv(u64);

impl Default for Fnv {
    fn default() -> Self {
        Fnv(OFFSET)
    }
}

impl Fnv {
    pub fn bytes(&mut self, data: &[u8]) {
        for &b in data {
            self.0 ^= u64::from(b);
            self.0 = self.0.wrapping_mul(PRIME);
        }
    }

    /// A field with its length in front, so two fields never read as one.
    pub fn field(&mut self, data: &[u8]) {
        self.bytes(&(data.len() as u64).to_le_bytes());
        self.bytes(data);
    }

    pub fn finish(self) -> u64 {
        self.0
    }
}

/// The stamp of one build, as 16 hex digits.
pub fn stamp(version: u32, tree_hash: &str, stdin: &[u8], inputs: &[PathBuf]) -> String {
    let mut h = Fnv::default();
    h.field(&version.to_le_bytes());
    h.field(tree_hash.as_bytes());
    h.field(stdin);
    for p in inputs {
        h.field(p.as_os_str().as_encoded_bytes());
        match fs::metadata(p) {
            Ok(m) => {
                let ns = m
                    .modified()
                    .ok()
                    .and_then(|t| t.duration_since(UNIX_EPOCH).ok())
                    .map_or(0, |d| d.as_nanos());
                h.field(b"present");
                h.field(&m.len().to_le_bytes());
                h.field(&ns.to_le_bytes());
            }
            Err(_) => h.field(b"absent"),
        }
    }
    format!("{:016x}", h.finish())
}

/// Where the document of `run` with `stamp` is cached.
pub fn cache_path(cache_dir: &Path, run: &str, stamp: &str) -> PathBuf {
    cache_dir.join(run).join(format!("{stamp}.json.gz"))
}

/// The cached document, if there is one.
pub fn lookup(path: &Path) -> Option<Vec<u8>> {
    fs::read(path).ok()
}

/// Writes `doc` at `path` atomically, then removes the run's other stamps.
pub fn store(path: &Path, doc: &[u8]) -> io::Result<()> {
    let dir = path.parent().expect("a cache path has a run folder");
    fs::create_dir_all(dir)?;
    let mut f = og_io::AtomicFile::create(path)?;
    f.write_all(doc)?;
    f.commit()?;
    let keep = path.file_name();
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let name = entry.file_name();
        let s = name.to_string_lossy();
        // Another writer's temporary file starts with a dot: never touch it.
        if Some(name.as_os_str()) != keep && !s.starts_with('.') && s.ends_with(".json.gz") {
            fs::remove_file(entry.path())?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fnv_known_values() {
        // The published FNV-1a 64 test vectors.
        let h = |s: &[u8]| {
            let mut f = Fnv::default();
            f.bytes(s);
            f.finish()
        };
        assert_eq!(h(b""), 0xcbf29ce484222325);
        assert_eq!(h(b"a"), 0xaf63dc4c8601ec8c);
        assert_eq!(h(b"foobar"), 0x85944171f73967e8);
    }

    #[test]
    fn stamp_follows_every_input() {
        let dir = std::env::temp_dir().join(format!("og-map-stamp-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let f = dir.join("in.txt");
        let inputs = [f.clone()];
        let absent = stamp(1, "t", b"{}", &inputs);
        fs::write(&f, "a").unwrap();
        let one = stamp(1, "t", b"{}", &inputs);
        assert_ne!(absent, one);
        assert_eq!(one, stamp(1, "t", b"{}", &inputs));
        assert_ne!(one, stamp(2, "t", b"{}", &inputs));
        assert_ne!(one, stamp(1, "u", b"{}", &inputs));
        assert_ne!(one, stamp(1, "t", b"{ }", &inputs));
        fs::write(&f, "ab").unwrap();
        assert_ne!(one, stamp(1, "t", b"{}", &inputs));
        assert_eq!(one.len(), 16);

        let cache = dir.join("cache");
        let old = cache_path(&cache, "DEMO-001", "0000000000000001");
        let new = cache_path(&cache, "DEMO-001", "0000000000000002");
        let other = cache_path(&cache, "DEMO-002", "0000000000000001");
        store(&old, b"old").unwrap();
        store(&other, b"other").unwrap();
        fs::write(old.with_file_name(".x.json.gz.tmp-1"), b"busy").unwrap();
        assert_eq!(lookup(&old).as_deref(), Some(&b"old"[..]));
        store(&new, b"new").unwrap();
        assert_eq!(lookup(&old), None);
        assert_eq!(lookup(&new).as_deref(), Some(&b"new"[..]));
        assert_eq!(lookup(&other).as_deref(), Some(&b"other"[..]));
        assert_eq!(fs::read_dir(new.parent().unwrap()).unwrap().count(), 2);
        fs::remove_dir_all(&dir).unwrap();
    }
}
