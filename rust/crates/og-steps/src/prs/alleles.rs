//! `score_alleles.tsv`: every effect and other allele of every written row,
//! `chrN<TAB>pos<TAB>allele`, unique, sorted by (chrom as a string, pos as
//! a number, allele as a string), the order Python's
//! `sorted(set, key=lambda t: (t[0], int(t[1]) if t[1].isdigit() else 0, t[2]))`
//! gives.
//!
//! Tens of millions of alleles pass through here, most of them shared by
//! several scores, so each is packed into one `u64`: the chromosome's rank
//! in string order (5 bits), the position (32 bits) and an id of the allele
//! string (27 bits). Each score's keys are sorted and deduplicated, then
//! merged into the running set in place; only the ids' order differs from
//! the strings' order, and that is fixed per position when the file is
//! written. A position that is not a plain decimal of at most 32 bits
//! without a leading zero (never seen in a PGS Catalog file) goes to a
//! short list sorted with the Python key instead.

use pycompat::text::isdigit;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::io::{self, Write};

/// The chromosomes a written row can be on, in the order of their names
/// with `chr` in front: `chr1 < chr10 < ... < chr19 < chr2 < ... < chrX`.
const CHROMS: [&str; 23] = [
    "1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "2", "20", "21", "22", "3",
    "4", "5", "6", "7", "8", "9", "X",
];
const ID_BITS: u32 = 27;
const POS_BITS: u32 = 32;
const MAX_ID: u32 = (1 << ID_BITS) - 1;

/// The rank of a chromosome name (without `chr`) in [`CHROMS`].
pub fn chrom_rank(chrom: &str) -> Option<u8> {
    CHROMS
        .iter()
        .position(|c| *c == chrom)
        .and_then(|i| u8::try_from(i).ok())
}

/// A position as packed, or `None` for one the short list keeps.
fn packed_pos(pos: &str) -> Option<u32> {
    let plain = !pos.is_empty()
        && pos.bytes().all(|b| b.is_ascii_digit())
        && (pos.len() == 1 || !pos.starts_with('0'));
    if plain { pos.parse().ok() } else { None }
}

/// Python's `int(pos) if pos.isdigit() else 0`. A digit string Python's
/// `int` cannot read (a superscript) sorts as 0 here; Python would stop.
fn py_pos_key(pos: &str) -> i128 {
    if isdigit(pos) {
        pycompat::num::parse_int(pos).unwrap_or(0)
    } else {
        0
    }
}

/// An allele the packed form cannot hold.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Odd {
    rank: u8,
    key: i128,
    pos: String,
    allele: String,
}

impl Odd {
    /// The Python key, then the position string (Python leaves such ties
    /// in the order of its set, which is not stable from run to run).
    fn cmp_key(&self, other: &Odd) -> Ordering {
        (self.rank, self.key, &self.allele, &self.pos).cmp(&(
            other.rank,
            other.key,
            &other.allele,
            &other.pos,
        ))
    }
}

/// The alleles of one score, added row by row.
#[derive(Debug, Default)]
pub struct ScoreAlleles {
    keys: Vec<u64>,
    odd: Vec<Odd>,
}

/// Every score's alleles so far, and the allele strings they use.
#[derive(Debug, Default)]
pub struct Alleles {
    ids: HashMap<String, u32>,
    names: Vec<String>,
    keys: Vec<u64>,
    odd: Vec<Odd>,
}

impl Alleles {
    /// Adds `allele` at `pos` of the chromosome of rank `rank` to `score`.
    pub fn push(&mut self, score: &mut ScoreAlleles, rank: u8, pos: &str, allele: &str) {
        if let Some(p) = packed_pos(pos)
            && let Some(id) = self.id(allele)
        {
            score.keys.push(
                (u64::from(rank) << (POS_BITS + ID_BITS))
                    | (u64::from(p) << ID_BITS)
                    | u64::from(id),
            );
            return;
        }
        score.odd.push(Odd {
            rank,
            key: py_pos_key(pos),
            pos: pos.to_owned(),
            allele: allele.to_owned(),
        });
    }

    fn id(&mut self, allele: &str) -> Option<u32> {
        if let Some(&id) = self.ids.get(allele) {
            return Some(id);
        }
        let id = u32::try_from(self.names.len())
            .ok()
            .filter(|&i| i <= MAX_ID)?;
        self.ids.insert(allele.to_owned(), id);
        self.names.push(allele.to_owned());
        Some(id)
    }

    /// Merges a finished score into the set: sorted, deduplicated, then
    /// merged from the back into the room made at the end of the set, so
    /// the peak is the set plus this score.
    pub fn add(&mut self, mut score: ScoreAlleles) {
        score.keys.sort_unstable();
        score.keys.dedup();
        let (old, new) = (self.keys.len(), score.keys.len());
        self.keys.reserve_exact(new);
        self.keys.resize(old + new, 0);
        let (mut i, mut j, mut k) = (old, new, old + new);
        while j > 0 {
            k -= 1;
            if i > 0 && self.keys[i - 1] > score.keys[j - 1] {
                i -= 1;
                self.keys[k] = self.keys[i];
            } else {
                j -= 1;
                self.keys[k] = score.keys[j];
            }
        }
        self.keys.dedup();
        self.odd.append(&mut score.odd);
    }

    /// Writes the set, one `chrN<TAB>pos<TAB>allele` line each, in Python's order.
    pub fn write(mut self, out: &mut impl Write) -> io::Result<()> {
        self.odd.sort_unstable_by(Odd::cmp_key);
        self.odd.dedup();
        let mut odd = self.odd.iter().peekable();
        let mut group: Vec<&str> = Vec::new();
        let mut rest = self.keys.as_slice();
        while let Some(&first) = rest.first() {
            let site = first >> ID_BITS;
            let n = rest.iter().take_while(|&&k| k >> ID_BITS == site).count();
            group.clear();
            group.extend(
                rest[..n]
                    .iter()
                    .map(|&k| self.names[(k & u64::from(MAX_ID)) as usize].as_str()),
            );
            group.sort_unstable_by(|a, b| b.cmp(a));
            rest = &rest[n..];
            let rank = (site >> POS_BITS) as u8;
            let pos = (site & ((1 << POS_BITS) - 1)) as u32;
            let pos_text = pos.to_string();
            for allele in &group {
                while let Some(o) = odd.next_if(|o| {
                    (o.rank, o.key, o.allele.as_str(), o.pos.as_str())
                        < (rank, i128::from(pos), *allele, pos_text.as_str())
                }) {
                    write_row(out, o.rank, &o.pos, &o.allele)?;
                }
                write_row(out, rank, &pos_text, allele)?;
            }
        }
        for o in odd {
            write_row(out, o.rank, &o.pos, &o.allele)?;
        }
        Ok(())
    }
}

fn write_row(out: &mut impl Write, rank: u8, pos: &str, allele: &str) -> io::Result<()> {
    writeln!(out, "chr{}\t{pos}\t{allele}", CHROMS[usize::from(rank)])
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The Python sort of the unique tuples, written out.
    fn python_order(rows: &[(&str, &str, &str)]) -> String {
        let mut set: Vec<(String, String, String)> = rows
            .iter()
            .map(|(c, p, a)| (format!("chr{c}"), (*p).to_owned(), (*a).to_owned()))
            .collect();
        set.sort();
        set.dedup();
        set.sort_by(|x, y| {
            (&x.0, py_pos_key(&x.1), &x.2, &x.1).cmp(&(&y.0, py_pos_key(&y.1), &y.2, &y.1))
        });
        set.iter()
            .map(|(c, p, a)| format!("{c}\t{p}\t{a}\n"))
            .collect()
    }

    fn ours(scores: &[&[(&str, &str, &str)]]) -> String {
        let mut all = Alleles::default();
        for rows in scores {
            let mut s = ScoreAlleles::default();
            for (c, p, a) in rows.iter() {
                all.push(&mut s, chrom_rank(c).expect("chrom"), p, a);
            }
            all.add(s);
        }
        let mut out = Vec::new();
        all.write(&mut out).expect("write");
        String::from_utf8(out).expect("utf8")
    }

    #[test]
    fn chromosomes_sort_as_strings_and_positions_as_numbers() {
        let rows = [
            ("2", "5", "A"),
            ("10", "100", "G"),
            ("1", "20", "T"),
            ("1", "3", "C"),
            ("X", "1", "A"),
            ("1", "20", "AT"),
            ("1", "20", "A"),
        ];
        assert_eq!(ours(&[&rows]), python_order(&rows));
        assert!(
            ours(&[&rows])
                .starts_with("chr1\t3\tC\nchr1\t20\tA\nchr1\t20\tAT\nchr1\t20\tT\nchr10\t")
        );
    }

    #[test]
    fn scores_merge_into_one_unique_set() {
        let a = [("1", "10", "A"), ("1", "10", "G"), ("3", "7", "C")];
        let b = [
            ("1", "10", "G"),
            ("1", "9", "T"),
            ("22", "1", "A"),
            ("3", "7", "C"),
        ];
        let mut both = a.to_vec();
        both.extend_from_slice(&b);
        assert_eq!(ours(&[&a, &b]), python_order(&both));
    }

    #[test]
    fn odd_positions_sort_with_the_python_key() {
        let rows = [
            ("1", "010", "A"),
            ("1", "9", "A"),
            ("1", "11", "A"),
            ("1", "abc", "A"),
            ("1", "4294967296", "C"),
            ("1", "4294967295", "C"),
            ("1", "0", "G"),
            ("1", "10", "A"),
            ("1", "\u{0661}2", "T"),
        ];
        assert_eq!(ours(&[&rows]), python_order(&rows));
    }

    /// The packed order against a brute-force sort of the strings, on many
    /// random rows over many scores (a small xorshift, so no crate is needed).
    #[test]
    fn random_rows_match_the_string_sort() {
        let alleles = [
            "A",
            "C",
            "G",
            "T",
            "AT",
            "GCC",
            "TA",
            "a",
            "N",
            "AAAAAAAAAAAAAAAA",
        ];
        let positions = [
            "1",
            "2",
            "9",
            "10",
            "11",
            "99",
            "100",
            "248956422",
            "4294967295",
            "007",
            "x",
        ];
        let mut seed = 0x9e37_79b9_7f4a_7c15_u64;
        let mut next = |n: usize| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed % n as u64) as usize
        };
        let mut scores: Vec<Vec<(&str, &str, &str)>> = Vec::new();
        for _ in 0..7 {
            let n = 1 + next(20_000);
            scores.push(
                (0..n)
                    .map(|_| {
                        (
                            CHROMS[next(CHROMS.len())],
                            positions[next(positions.len())],
                            alleles[next(alleles.len())],
                        )
                    })
                    .collect(),
            );
        }
        let refs: Vec<&[(&str, &str, &str)]> = scores.iter().map(Vec::as_slice).collect();
        let all: Vec<(&str, &str, &str)> = scores.concat();
        assert_eq!(ours(&refs), python_order(&all));
    }
}
