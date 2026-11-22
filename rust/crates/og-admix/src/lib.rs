//! What the fine-ancestry engines of step 26b share, so each engine is
//! built against one fixed interface:
//!
//! - [`panel`]: a phased panel VCF read into bit-packed haplotypes, for
//!   chosen people at chosen sites.
//! - [`map`]: a PLINK genetic map, cM by position and position by cM.
//! - [`geno`]: `model/ref.geno`, the reference people's genotypes, 2 bits
//!   each, with its writer and reader.
//! - [`rng`]: a seeded xoshiro256** generator, so a seed gives the same
//!   bytes on every machine.
//! - [`simulate`]: an admixed diploid genome with known truth, built from
//!   held-out people's haplotypes; every engine's truth tests start from it.
//! - [`paint`]: our own local-ancestry painting.
//! - [`project`]: the ancestry fractions projected on fixed group
//!   frequencies, in place of `admixture -P`, and its block bootstrap.

use std::fmt;
use std::io;
use std::path::Path;

pub mod geno;
pub mod map;
pub mod paint;
pub mod panel;
pub mod project;
pub mod rng;
pub mod simulate;

/// Why a reader or a writer stopped; the text is a whole sentence for an
/// `ERROR:` line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Error(pub String);

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for Error {}

pub type Result<T> = std::result::Result<T, Error>;

/// An I/O error's text without Rust's ` (os error N)`.
fn reason(e: &io::Error) -> String {
    let s = e.to_string();
    match s.rfind(" (os error ") {
        Some(i) if s.ends_with(')') => s[..i].to_owned(),
        _ => s,
    }
}

fn cannot_read(path: &Path, e: &io::Error) -> Error {
    Error(format!("cannot read {}: {}", path.display(), reason(e)))
}

fn cannot_write(path: &Path, e: &io::Error) -> Error {
    Error(format!("cannot write {}: {}", path.display(), reason(e)))
}

/// Python's `repr()` of a str, as the messages of og ancestry-fine quote a
/// value.
fn repr(s: &str) -> String {
    let quote = if s.contains('\'') && !s.contains('"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::with_capacity(s.len() + 2);
    out.push(quote);
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            c if (c as u32) < 0x20 || c as u32 == 0x7f => {
                out.push_str(&format!("\\x{:02x}", c as u32));
            }
            c => out.push(c),
        }
    }
    out.push(quote);
    out
}

/// True when the red controls are planted (feature `planted-red`).
#[cfg(test)]
pub(crate) const PLANTED: bool = cfg!(feature = "planted-red");

#[cfg(test)]
pub(crate) mod testutil {
    use std::path::PathBuf;

    /// A fresh directory under the temp dir, removed on drop.
    pub struct Dir(pub PathBuf);

    impl Dir {
        pub fn new(name: &str) -> Dir {
            let p = std::env::temp_dir().join(format!(
                "og-admix-{name}-{}-{:?}",
                std::process::id(),
                std::thread::current().id()
            ));
            let _ = std::fs::remove_dir_all(&p);
            std::fs::create_dir_all(&p).unwrap();
            Dir(p)
        }

        /// Writes `text` at `name` (folders made as needed); its path.
        pub fn put(&self, name: &str, text: &str) -> PathBuf {
            let p = self.0.join(name);
            if let Some(d) = p.parent() {
                std::fs::create_dir_all(d).unwrap();
            }
            std::fs::write(&p, text).unwrap();
            p
        }
    }

    impl Drop for Dir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
}
