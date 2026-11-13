//! UCSC chain files, for lifting single positions from one assembly to
//! another (`og ancient lift`).
//!
//! A chain file is a list of chains. Each starts with a header line,
//! `chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart
//! qEnd id`, followed by block lines `size dt dq` (the last one `size`
//! alone) and a blank line. The `t` side is the assembly a position is
//! lifted from (the source), the `q` side the one it is lifted to. Each
//! block is an ungapped stretch of `size` bases; `dt` and `dq` are the gaps
//! before the next block on each side.
//!
//! [`Chain::lift`] maps a position only when the blocks that cover it all
//! belong to one chain; a position in a gap, or covered by two chains, has
//! no single place and gives `None`.

use og_io::{Lines, Utf8};
use pycompat::num::{parse_int, py_isspace};
use std::collections::HashMap;
use std::fmt;
use std::io::{self, BufRead};

/// Why a chain file could not be read.
#[derive(Debug)]
pub enum ChainError {
    Io(io::Error),
    /// A chain whose source strand is not `+`: the chain's id (or its
    /// 0-based number when the header has no id) and the strand.
    SourceStrand {
        chain: String,
        strand: String,
    },
    /// A block line before any chain header, or after a chain's last block.
    Orphan(String),
    /// A line that cannot be used: its number (1-based) and why.
    Line {
        line: usize,
        why: String,
    },
}

impl fmt::Display for ChainError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ChainError::Io(e) => write!(f, "{e}"),
            ChainError::SourceStrand { chain, strand } => {
                write!(f, "chain {chain} has source strand {strand}")
            }
            ChainError::Orphan(line) => write!(f, "a block line outside a chain: {line}"),
            ChainError::Line { line, why } => write!(f, "line {line}: {why}"),
        }
    }
}

impl std::error::Error for ChainError {}

impl From<io::Error> for ChainError {
    fn from(e: io::Error) -> Self {
        ChainError::Io(e)
    }
}

/// One ungapped block: `[start, end)` on the source, 0-based.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Block {
    start: u32,
    end: u32,
    /// The chain's number in the file (0-based).
    chain: u32,
    /// The target sequence's name, an index into [`Chain::names`].
    name: u32,
    /// The block's first base on the target, counted on its strand.
    to_start: u32,
    /// The target strand, an index into [`Chain::names`].
    strand: u32,
    /// The target sequence's length.
    to_size: u32,
}

/// The blocks of one source sequence, sorted, with the furthest end of
/// each prefix, so a lookup stops where no earlier block can reach.
#[derive(Debug, Default)]
struct Source {
    blocks: Vec<Block>,
    reach: Vec<u32>,
}

/// A chain file, indexed by source sequence.
#[derive(Debug, Default)]
pub struct Chain {
    /// Target names and strands, interned.
    names: Vec<String>,
    sources: HashMap<String, Source>,
}

/// Where [`Chain::lift`] puts a position.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Lifted<'a> {
    /// The target sequence.
    pub chrom: &'a str,
    /// The 1-based position on it.
    pub pos: i64,
    /// The target strand as the file writes it (`+` or `-`).
    pub strand: &'a str,
}

/// The header fields a block line needs, kept while the chain's blocks run.
struct Open {
    source: String,
    at: i64,
    name: u32,
    to_at: i64,
    strand: u32,
    to_size: i64,
    chain: u32,
}

impl Chain {
    /// Reads a chain file (universal newlines; invalid UTF-8 is replaced).
    /// Fields are split on whitespace, so a header with tabs reads too.
    pub fn read(input: impl BufRead) -> Result<Chain, ChainError> {
        let mut lines = Lines::new(input, Utf8::Replace);
        let mut text = String::new();
        let mut out = Chain::default();
        let mut interned: HashMap<String, u32> = HashMap::new();
        let mut intern = |names: &mut Vec<String>, s: &str| -> u32 {
            if let Some(&i) = interned.get(s) {
                return i;
            }
            let i = names.len() as u32;
            names.push(s.to_owned());
            interned.insert(s.to_owned(), i);
            i
        };
        let mut cur: Option<Open> = None;
        let mut chains: u32 = 0;
        let mut n = 0;
        while lines.read_line(&mut text)? {
            n += 1;
            let c: Vec<&str> = text.split(py_isspace).filter(|f| !f.is_empty()).collect();
            let bad = |why: String| ChainError::Line { line: n, why };
            if c.is_empty() {
                cur = None;
                continue;
            }
            if c[0] == "chain" {
                if c.len() < 11 {
                    return Err(bad(format!(
                        "a chain header of {} fields, want 13",
                        c.len()
                    )));
                }
                let number = chains;
                chains += 1;
                if c[4] != "+" {
                    return Err(ChainError::SourceStrand {
                        chain: c
                            .get(12)
                            .map_or_else(|| number.to_string(), |s| (*s).to_owned()),
                        strand: c[4].to_owned(),
                    });
                }
                let num = |k: usize| coord(c[k], n, k + 1);
                cur = Some(Open {
                    source: c[2].to_owned(),
                    at: num(5)?,
                    name: intern(&mut out.names, c[7]),
                    to_size: num(8)?,
                    strand: intern(&mut out.names, c[9]),
                    to_at: num(10)?,
                    chain: number,
                });
                continue;
            }
            let Some(open) = cur.as_mut() else {
                return Err(ChainError::Orphan(
                    pycompat::num::py_strip(&text).to_owned(),
                ));
            };
            let size = coord(c[0], n, 1)?;
            let fit = |v: i64| {
                u32::try_from(v).map_err(|_| bad(format!("position {v} does not fit 32 bits")))
            };
            let block = Block {
                start: fit(open.at)?,
                end: fit(open.at + size)?,
                chain: open.chain,
                name: open.name,
                to_start: fit(open.to_at)?,
                strand: open.strand,
                to_size: fit(open.to_size)?,
            };
            out.sources
                .entry(open.source.clone())
                .or_default()
                .blocks
                .push(block);
            if c.len() >= 3 {
                open.at += size + coord(c[1], n, 2)?;
                open.to_at += size + coord(c[2], n, 3)?;
            } else {
                cur = None;
            }
        }
        for s in out.sources.values_mut() {
            let names = &out.names;
            s.blocks.sort_by(|a, b| {
                (a.start, a.end, a.chain)
                    .cmp(&(b.start, b.end, b.chain))
                    .then_with(|| names[a.name as usize].cmp(&names[b.name as usize]))
                    .then_with(|| a.to_start.cmp(&b.to_start))
                    .then_with(|| names[a.strand as usize].cmp(&names[b.strand as usize]))
                    .then_with(|| a.to_size.cmp(&b.to_size))
            });
            let mut far = 0;
            s.reach = s
                .blocks
                .iter()
                .map(|b| {
                    far = far.max(b.end);
                    far
                })
                .collect();
        }
        Ok(out)
    }

    /// How many blocks were read.
    pub fn blocks(&self) -> usize {
        self.sources.values().map(|s| s.blocks.len()).sum()
    }

    /// The 1-based position `pos` of source sequence `chrom` on the target,
    /// or `None` when no block covers it, or blocks of more than one chain
    /// do. On a `-` target strand the position is mirrored:
    /// `size - 1 - offset`, 0-based.
    pub fn lift(&self, chrom: &str, pos: i64) -> Option<Lifted<'_>> {
        let src = self.sources.get(chrom)?;
        let p = pos - 1;
        // Blocks that start at or before p, walked back while one of them
        // may still reach past p.
        let mut j = src.blocks.partition_point(|b| i64::from(b.start) <= p);
        let mut first: Option<&Block> = None;
        while j > 0 && i64::from(src.reach[j - 1]) >= p {
            j -= 1;
            let b = &src.blocks[j];
            if i64::from(b.end) >= p {
                match first {
                    None => first = Some(b),
                    Some(f) if f.chain != b.chain => return None,
                    Some(_) => {}
                }
            }
        }
        let b = first?;
        let mut q = i64::from(b.to_start) + (p - i64::from(b.start));
        let strand = self.names[b.strand as usize].as_str();
        if strand == "-" {
            q = i64::from(b.to_size) - 1 - q;
        }
        Some(Lifted {
            chrom: &self.names[b.name as usize],
            pos: q + 1,
            strand,
        })
    }
}

/// A whole-number field of line `line`, column `col`: Python's `int()`.
fn coord(field: &str, line: usize, col: usize) -> Result<i64, ChainError> {
    parse_int(field)
        .ok()
        .and_then(|v| i64::try_from(v).ok())
        .ok_or_else(|| ChainError::Line {
            line,
            why: format!("field {col} '{field}' is not a whole number"),
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The chain of tests/test_ancient_report.py: chr1 with a gap, chr2 on
    /// the minus strand, chr3 covered twice.
    const TEST: &str = "chain 100 chr1 100000 + 1000 1200 chr1 100000 + 2000 2200 1
100 50 50
50

chain 90 chr2 50000 + 500 600 chr2 10000 - 300 400 2
100

chain 80 chr3 50000 + 100 200 chr3 50000 + 100 200 3
100

chain 70 chr3 50000 + 150 250 chr3 50000 + 5150 5250 4
100
";

    fn chain() -> Chain {
        Chain::read(TEST.as_bytes()).unwrap()
    }

    fn lift(c: &Chain, chrom: &str, pos: i64) -> Option<(String, i64, String)> {
        c.lift(chrom, pos)
            .map(|l| (l.chrom.to_owned(), l.pos, l.strand.to_owned()))
    }

    fn at(chrom: &str, pos: i64, strand: &str) -> Option<(String, i64, String)> {
        Some((chrom.to_owned(), pos, strand.to_owned()))
    }

    #[test]
    fn blocks_move_by_their_offset() {
        let c = chain();
        assert_eq!(c.blocks(), 5);
        assert_eq!(lift(&c, "chr1", 1051), at("chr1", 2051, "+"));
        assert_eq!(lift(&c, "chr1", 1160), at("chr1", 2160, "+"));
        assert_eq!(lift(&c, "chr1", 1120), None, "a gap has no place");
        assert_eq!(lift(&c, "chr9", 1051), None, "an unknown sequence");
    }

    #[test]
    fn block_ends_are_half_open() {
        let c = chain();
        // The first block holds 0-based 1000..1100: 1-based 1001 to 1100.
        assert_eq!(lift(&c, "chr1", 1001), at("chr1", 2001, "+"));
        assert_eq!(lift(&c, "chr1", 1100), at("chr1", 2100, "+"));
        assert_eq!(lift(&c, "chr1", 1000), None);
        assert_eq!(lift(&c, "chr1", 1101), None);
        assert_eq!(lift(&c, "chr1", 1151), at("chr1", 2151, "+"));
        assert_eq!(lift(&c, "chr1", 1200), at("chr1", 2200, "+"));
        assert_eq!(lift(&c, "chr1", 1201), None);
        assert_eq!(lift(&c, "chr1", 0), None);
        assert_eq!(lift(&c, "chr1", -5), None);
    }

    #[test]
    fn minus_strand_mirrors() {
        // 0-based 510 is 310 on the reverse strand: 10000 - 1 - 310 = 9689.
        assert_eq!(lift(&chain(), "chr2", 511), at("chr2", 9690, "-"));
    }

    #[test]
    fn two_chains_over_a_position_give_none() {
        let c = chain();
        assert_eq!(lift(&c, "chr3", 171), None);
        assert_eq!(lift(&c, "chr3", 121), at("chr3", 121, "+"));
        assert_eq!(lift(&c, "chr3", 221), at("chr3", 5221, "+"));
    }

    #[test]
    fn a_long_block_far_back_still_counts() {
        // A long block of chain 1 starts long before a short block of
        // chain 2; a position inside both is covered twice.
        let text = "chain 1 c 9000 + 0 5000 c 9000 + 0 5000 1\n5000\n\n\
                    chain 1 c 9000 + 4000 4010 d 9000 + 0 10 2\n10\n";
        let c = Chain::read(text.as_bytes()).unwrap();
        assert_eq!(lift(&c, "c", 4005), None);
        assert_eq!(lift(&c, "c", 3000), at("c", 3000, "+"));
    }

    #[test]
    fn errors() {
        let e = Chain::read("chain 1 c 10 - 0 5 c 10 + 0 5 7\n5\n".as_bytes()).unwrap_err();
        assert_eq!(e.to_string(), "chain 7 has source strand -");
        let e = Chain::read("chain 1 c 10 - 0 5 c 10 + 0 5\n5\n".as_bytes()).unwrap_err();
        assert_eq!(e.to_string(), "chain 0 has source strand -");
        let e = Chain::read("  5 1 1 \n".as_bytes()).unwrap_err();
        assert_eq!(e.to_string(), "a block line outside a chain: 5 1 1");
        // After a block line with one field the chain is over.
        let e = Chain::read("chain 1 c 10 + 0 5 c 10 + 0 5 1\n5\n6\n".as_bytes()).unwrap_err();
        assert_eq!(e.to_string(), "a block line outside a chain: 6");
        let e = Chain::read("chain 1 c 10 + x 5 c 10 + 0 5 1\n".as_bytes()).unwrap_err();
        assert_eq!(e.to_string(), "line 1: field 6 'x' is not a whole number");
    }
}
