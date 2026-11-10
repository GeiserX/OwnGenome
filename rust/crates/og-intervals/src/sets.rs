//! Sorted interval sets, 0-based and half-open: `(start, end)` covers the
//! bases `start..end`, so `(10, 20)` holds ten bases and an interval with
//! `end == start` holds none.
//!
//! - [`Set`]: disjoint, merged intervals of one sequence; union,
//!   intersection, and point or range lookups by binary search.
//! - [`Sets`]: a [`Set`] per sequence, from a BED file such as the
//!   pseudoautosomal regions; [`Sets::contains_pos`] takes a VCF `POS`.
//! - [`Index`]: intervals kept as given (not merged) with a value each:
//!   which of them overlap a range, and whether one overlaps a range
//!   reciprocally (each covers at least a fraction of the other).

use crate::bed::{BedError, Rows, merge};
use std::collections::HashMap;
use std::io::BufRead;

/// Disjoint, ascending, non-empty half-open intervals.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Set(Vec<(u64, u64)>);

impl Set {
    /// The bases any of `rows` covers: sorted, overlapping and touching
    /// rows joined, empty rows dropped.
    pub fn new(rows: impl IntoIterator<Item = (u64, u64)>) -> Set {
        Set(merge(&rows.into_iter().collect::<Vec<_>>()))
    }

    /// The intervals.
    pub fn intervals(&self) -> &[(u64, u64)] {
        &self.0
    }

    /// Bases covered.
    pub fn bases(&self) -> u64 {
        self.0.iter().map(|(s, e)| e - s).sum()
    }

    /// The bases of either set.
    pub fn union(&self, other: &Set) -> Set {
        Set::new(self.0.iter().chain(&other.0).copied())
    }

    /// The bases of both sets.
    pub fn intersect(&self, other: &Set) -> Set {
        let (a, b) = (&self.0, &other.0);
        let (mut i, mut j) = (0, 0);
        let mut out = Vec::new();
        while i < a.len() && j < b.len() {
            let (s, e) = (a[i].0.max(b[j].0), a[i].1.min(b[j].1));
            if s < e {
                out.push((s, e));
            }
            if a[i].1 < b[j].1 { i += 1 } else { j += 1 }
        }
        Set(out)
    }

    /// Whether the 0-based base `pos` is covered.
    pub fn contains(&self, pos: u64) -> bool {
        let i = self.0.partition_point(|&(_, e)| e <= pos);
        self.0.get(i).is_some_and(|&(s, _)| s <= pos)
    }

    /// Whether any base of `[start, end)` is covered.
    pub fn overlaps(&self, start: u64, end: u64) -> bool {
        let i = self.0.partition_point(|&(_, e)| e <= start);
        self.0.get(i).is_some_and(|&(s, _)| s < end) && start < end
    }

    /// How many bases of `[start, end)` are covered.
    pub fn covered(&self, start: u64, end: u64) -> u64 {
        let i = self.0.partition_point(|&(_, e)| e <= start);
        self.0[i..]
            .iter()
            .take_while(|&&(s, _)| s < end)
            .map(|&(s, e)| e.min(end).saturating_sub(s.max(start)))
            .sum()
    }
}

/// A [`Set`] per sequence name.
#[derive(Debug, Clone, Default)]
pub struct Sets(HashMap<String, Set>);

impl Sets {
    /// Reads a BED file with [`Rows::read`]'s rules: whitespace-split,
    /// `#`, `track` and `browser` lines skipped, a bad coordinate an error.
    pub fn read(input: impl BufRead) -> Result<Sets, BedError> {
        let rows = Rows::read(input)?;
        Ok(Sets(
            rows.names()
                .map(|n| (n.to_owned(), Set(rows.merged(n).to_vec())))
                .collect(),
        ))
    }

    /// The set of `seq`, empty when the BED has no row on it.
    pub fn get(&self, seq: &str) -> Option<&Set> {
        self.0.get(seq)
    }

    /// Whether the 0-based base `pos` of `seq` is covered.
    pub fn contains(&self, seq: &str, pos: u64) -> bool {
        self.0.get(seq).is_some_and(|s| s.contains(pos))
    }

    /// Whether the 1-based position `pos` (a VCF `POS`) of `seq` is
    /// covered: a BED row `start end` covers `start < pos <= end`.
    pub fn contains_pos(&self, seq: &str, pos: u64) -> bool {
        pos > 0 && self.contains(seq, pos - 1)
    }
}

/// Intervals as given, each with a value, sorted by start.
#[derive(Debug, Clone)]
pub struct Index<T> {
    rows: Vec<(u64, u64, T)>,
    /// The longest interval, so a lookup knows how far back to look.
    longest: u64,
}

impl<T> Index<T> {
    /// Indexes `rows` (`start`, `end`, value); an empty row overlaps nothing.
    pub fn new(mut rows: Vec<(u64, u64, T)>) -> Index<T> {
        rows.sort_by_key(|r| (r.0, r.1));
        let longest = rows
            .iter()
            .map(|r| r.1.saturating_sub(r.0))
            .max()
            .unwrap_or(0);
        Index { rows, longest }
    }

    /// The rows sharing at least one base with `[start, end)`, by start.
    pub fn overlapping(&self, start: u64, end: u64) -> impl Iterator<Item = &(u64, u64, T)> {
        // A row ending after start began less than `longest` before it.
        let from = start.saturating_sub(self.longest);
        let lo = self.rows.partition_point(|r| r.0 < from);
        let hi = self.rows.partition_point(|r| r.0 < end);
        self.rows[lo..hi.max(lo)]
            .iter()
            .filter(move |r| r.1 > start && r.0 < r.1 && start < end)
    }

    /// The first row (by start) that overlaps `[start, end)` reciprocally:
    /// the shared bases are at least `frac` of each of the two lengths.
    pub fn reciprocal(&self, start: u64, end: u64, frac: f64) -> Option<&(u64, u64, T)> {
        let len = end.saturating_sub(start);
        self.overlapping(start, end).find(|r| {
            let shared = r.1.min(end) - r.0.max(start);
            reciprocal(shared, len, r.1 - r.0, frac)
        })
    }
}

/// Whether `shared` bases are at least `frac` of both `a` and `b`, compared
/// exactly as Python compares an int with a float.
pub fn reciprocal(shared: u64, a: u64, b: u64, frac: f64) -> bool {
    shared > 0 && shared as f64 >= frac * a as f64 && shared as f64 >= frac * b as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    /// xorshift64: the same rows on every run.
    struct Rng(u64);
    impl Rng {
        fn below(&mut self, n: u64) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0 % n
        }
        fn rows(&mut self, n: usize, span: u64, max_len: u64) -> Vec<(u64, u64)> {
            (0..n)
                .map(|_| {
                    let s = self.below(span);
                    (s, s + self.below(max_len + 1))
                })
                .collect()
        }
    }

    /// The covered bases of `rows` below `span`, one flag per base.
    fn bitmap(rows: &[(u64, u64)], span: u64) -> Vec<bool> {
        let mut m = vec![false; span as usize];
        for &(s, e) in rows {
            for p in s..e.min(span) {
                m[p as usize] = true;
            }
        }
        m
    }

    #[test]
    fn set_operations_equal_a_base_by_base_model() {
        let mut rng = Rng(42);
        for round in 0..300 {
            let span = 400;
            let n = rng.below(12) as usize;
            let ra = rng.rows(n, span - 40, 40);
            let n = rng.below(12) as usize;
            let rb = rng.rows(n, span - 40, 40);
            let (a, b) = (Set::new(ra.clone()), Set::new(rb.clone()));
            let (ma, mb) = (bitmap(&ra, span), bitmap(&rb, span));
            // Disjoint, ascending, not touching, none empty.
            for s in [&a, &b] {
                assert!(s.0.iter().all(|(x, y)| x < y));
                assert!(s.0.windows(2).all(|w| w[0].1 < w[1].0), "{round}: {s:?}");
            }
            let mu: Vec<bool> = ma.iter().zip(&mb).map(|(x, y)| *x || *y).collect();
            let mi: Vec<bool> = ma.iter().zip(&mb).map(|(x, y)| *x && *y).collect();
            assert_eq!(bitmap(a.union(&b).intervals(), span), mu, "{round} union");
            assert_eq!(
                bitmap(a.intersect(&b).intervals(), span),
                mi,
                "{round} intersect"
            );
            assert_eq!(a.intersect(&b), b.intersect(&a));
            assert_eq!(a.union(&b), b.union(&a));
            assert_eq!(a.bases(), ma.iter().filter(|x| **x).count() as u64);
            for p in 0..span {
                assert_eq!(a.contains(p), ma[p as usize], "{round} contains {p}");
            }
            for _ in 0..50 {
                let s = rng.below(span);
                let e = s + rng.below(span - s);
                let want = ma[s as usize..e as usize].iter().filter(|x| **x).count() as u64;
                assert_eq!(a.covered(s, e), want, "{round} covered [{s},{e})");
                assert_eq!(a.overlaps(s, e), want > 0, "{round} overlaps [{s},{e})");
            }
        }
    }

    #[test]
    fn index_queries_equal_a_full_scan() {
        let mut rng = Rng(7);
        for round in 0..200 {
            let n = rng.below(40) as usize;
            let rows: Vec<(u64, u64, usize)> = rng
                .rows(n, 1_000, if round % 3 == 0 { 600 } else { 30 })
                .into_iter()
                .enumerate()
                .map(|(i, (s, e))| (s, e, i))
                .collect();
            let idx = Index::new(rows.clone());
            for _ in 0..60 {
                let s = rng.below(1_100);
                let e = s + rng.below(300);
                let mut want: Vec<_> = rows
                    .iter()
                    .filter(|r| r.0 < r.1 && r.0 < e && r.1 > s && s < e)
                    .collect();
                want.sort_by_key(|r| (r.0, r.1));
                let got: Vec<_> = idx.overlapping(s, e).collect();
                let key = |v: &[&(u64, u64, usize)]| v.iter().map(|r| r.2).collect::<Vec<_>>();
                let mut g = key(&got);
                let mut w = key(&want);
                g.sort_unstable();
                w.sort_unstable();
                assert_eq!(g, w, "{round} [{s},{e})");
                for frac in [0.5, 0.8, 1.0] {
                    let any = rows.iter().any(|r| {
                        let shared = r.1.min(e).saturating_sub(r.0.max(s));
                        shared > 0 && shared as f64 >= frac * (e - s).max(r.1 - r.0) as f64
                    });
                    assert_eq!(
                        idx.reciprocal(s, e, frac).is_some(),
                        any,
                        "{round} [{s},{e}) {frac}"
                    );
                }
            }
        }
    }

    #[test]
    fn reciprocal_overlap_is_measured_on_the_longer_call() {
        let idx = Index::new(vec![(100, 200, "a")]);
        // 50 shared bases: half of the 100-base row and of a 100-base call.
        assert!(idx.reciprocal(150, 250, 0.5).is_some());
        // The same 50 bases are only a third of a 150-base call.
        assert!(idx.reciprocal(150, 300, 0.5).is_none());
        // Inside the row: all of a 40-base call, but 40% of the row.
        assert!(idx.reciprocal(120, 160, 0.5).is_none());
        // Touching is not overlapping: half-open ends share no base.
        assert!(idx.reciprocal(200, 300, 0.0).is_none());
        assert!(idx.overlapping(200, 300).next().is_none());
        assert!(idx.overlapping(0, 100).next().is_none());
        assert_eq!(idx.overlapping(199, 200).count(), 1);
        assert_eq!(idx.overlapping(99, 101).count(), 1);
        assert!(reciprocal(1, 2, 2, 0.5));
        assert!(!reciprocal(0, 0, 0, 0.0));
    }

    #[test]
    fn par_rows_are_half_open_in_vcf_positions() {
        // GRCh38's PAR1 on chrX, as assets/par_grch38.bed writes it.
        let par = Sets::read("chrX\t10000\t2781479\nchrY\t10000\t2781479\n".as_bytes()).unwrap();
        // BED start 10000 is the 1-based POS 10001: POS 10000 is outside.
        assert!(!par.contains_pos("chrX", 10_000));
        assert!(par.contains_pos("chrX", 10_001));
        // BED end 2781479 is the last 1-based POS inside.
        assert!(par.contains_pos("chrX", 2_781_479));
        assert!(!par.contains_pos("chrX", 2_781_480));
        assert!(par.contains("chrX", 10_000));
        assert!(!par.contains("chrX", 2_781_479));
        assert!(!par.contains_pos("chrX", 0));
        assert!(!par.contains_pos("chr1", 20_000));
        assert_eq!(par.get("chrY").unwrap().bases(), 2_771_479);
    }

    #[test]
    fn empty_sets_and_ranges() {
        let e = Set::default();
        assert!(!e.contains(0));
        assert!(!e.overlaps(0, u64::MAX));
        let a = Set::new([(5, 10), (10, 12), (3, 3)]);
        assert_eq!(a.intervals(), &[(5, 12)]);
        assert!(!a.overlaps(7, 7));
        assert_eq!(a.covered(7, 7), 0);
        assert_eq!(a.intersect(&e), e);
        assert_eq!(a.union(&e), a);
        assert!(Index::<()>::new(vec![]).reciprocal(0, 10, 0.5).is_none());
    }
}
