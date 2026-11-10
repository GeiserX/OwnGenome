//! What a closed result could reveal, as tags in the document.
//!
//! The document is the same for every reader: nothing is left out of it.
//! Instead, every row that could show a result the page keeps closed
//! carries that result's id, and the page hides a tagged row with the same
//! rule it hides the result's card with. Whether a page is closed is the
//! page's rule alone; this module decides only what could be revealed.
//!
//! 1. Closable findings: level `act` (closed on a child's page and while
//!    the family's answers are not known), and every gated finding (`gate`
//!    set, APOE), which the page keeps behind its warning on every page,
//!    an adult's too.
//! 2. Veiled spans: each closable finding's place on the map: its position,
//!    its span, or its gene's span when it has no position.
//! 3. Row tags: a row of a layer that can be a finding's own source (a
//!    structural variant, a mosaic call, a ClinVar hit) that overlaps a
//!    veiled span gets that finding's index in `dict.veil`, else -1. A row
//!    holds one tag: when several spans overlap it, the first in
//!    `dict.veil`, where gated findings come first (closed on every page),
//!    then the act ones in the page's order. A finding's own row is tagged
//!    the same way, so an open finding over a closable one's span is
//!    tagged. Runs of homozygosity and bands are never a finding's source
//!    and carry no tag.
//! 4. Counts, totals and density bins cover untagged rows only, so no number on
//!    the map counts a closed row.

use crate::facts::Finding;
use crate::reference::{CHROMS, Reference};

/// The tag of a row no closed result could reveal.
pub const NONE: i64 = -1;

/// True for a finding a page can keep closed.
pub fn closable(f: &Finding) -> bool {
    f.level == "act" || f.gate.is_some()
}

/// A row on the map: its chromosome and its 0-based, end-exclusive span.
pub trait Span {
    /// Index into [`CHROMS`].
    fn chrom(&self) -> usize;
    fn span(&self) -> (u64, u64);
}

/// A point (`s == e`, an insertion or a breakend) still covers its base.
fn nonempty((s, e): (u64, u64)) -> (u64, u64) {
    (s, e.max(s + 1))
}

/// The veiled spans of one person's closable findings, by chromosome.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Veil {
    /// `dict.veil`: the closable findings' ids, gated ones first, then in
    /// the page's order.
    pub ids: Vec<String>,
    /// Per chromosome: (s, e, index into `ids`), in `ids` order.
    spans: Vec<Vec<(u64, u64, u32)>>,
}

impl Veil {
    /// The veil of the closable findings in `placed`, each given as (its
    /// key, id, chromosome, span). `ids` follow the keys, smallest first.
    pub fn new<'a, K: Ord>(
        placed: impl IntoIterator<Item = (K, &'a str, usize, (u64, u64))>,
    ) -> Veil {
        let mut all: Vec<_> = placed.into_iter().collect();
        all.sort_by(|a, b| a.0.cmp(&b.0));
        let mut v = Veil {
            ids: Vec::with_capacity(all.len()),
            spans: vec![Vec::new(); CHROMS.len()],
        };
        // In key order, so the first span a row meets is the first finding
        // by key.
        for (k, (_, id, chrom, span)) in all.into_iter().enumerate() {
            v.ids.push(id.to_owned());
            let (s, e) = nonempty(span);
            v.spans[chrom].push((s, e, k as u32));
        }
        v
    }

    /// The index in `ids` of the first closable finding whose span `row`
    /// overlaps, if any.
    pub fn tag(&self, chrom: usize, span: (u64, u64)) -> Option<u32> {
        let (s, e) = nonempty(span);
        self.spans
            .get(chrom)?
            .iter()
            .find(|&&(vs, ve, _)| s < ve && vs < e)
            .map(|x| x.2)
    }

    /// The tag of every row, in order.
    pub fn tags<R: Span>(&self, rows: &[R]) -> Vec<Option<u32>> {
        rows.iter().map(|r| self.tag(r.chrom(), r.span())).collect()
    }
}

/// A tag column as the document writes it: the index, or [`NONE`].
pub fn column(tags: &[Option<u32>]) -> Vec<i64> {
    tags.iter().map(|t| t.map_or(NONE, i64::from)).collect()
}

/// The rows no closed result could reveal.
pub fn untagged(tags: &[Option<u32>]) -> usize {
    tags.iter().filter(|t| t.is_none()).count()
}

/// The count of untagged rows per `size` bin, by the bin of each row's
/// start, for every chromosome with at least one untagged row. A tagged
/// row is never counted, so no bin points at a closed result.
pub fn bins<R: Span>(
    rows: &[R],
    tags: &[Option<u32>],
    reference: &Reference,
    size: u64,
) -> Vec<(usize, Vec<u64>)> {
    let mut per: Vec<Option<Vec<u64>>> = vec![None; CHROMS.len()];
    for (r, t) in rows.iter().zip(tags) {
        if t.is_some() {
            continue;
        }
        let c = r.chrom();
        let n = reference.chroms[c].len.div_ceil(size) as usize;
        let b = per[c].get_or_insert_with(|| vec![0; n]);
        let i = ((r.span().0 / size) as usize).min(n.saturating_sub(1));
        b[i] += 1;
    }
    per.into_iter()
        .enumerate()
        .filter_map(|(c, b)| b.map(|b| (c, b)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    struct R(usize, u64, u64);
    impl Span for R {
        fn chrom(&self) -> usize {
            self.0
        }
        fn span(&self) -> (u64, u64) {
            (self.1, self.2)
        }
    }

    #[test]
    fn half_open_overlap_and_the_first_id() {
        // Finding b comes first on the page though it was placed second.
        let v = Veil::new([(1, "a", 0, (100, 101)), (0, "b", 0, (50, 200))]);
        assert_eq!(v.ids, ["b", "a"]);
        assert_eq!(v.tag(0, (100, 101)), Some(0));
        assert_eq!(v.tag(0, (0, 50)), None);
        assert_eq!(v.tag(0, (200, 300)), None);
        assert_eq!(v.tag(0, (199, 199)), Some(0));
        assert_eq!(v.tag(1, (100, 101)), None);
        let rows = [R(0, 10, 20), R(0, 60, 61), R(1, 0, 5)];
        let tags = v.tags(&rows);
        assert_eq!(column(&tags), [NONE, 0, NONE]);
        assert_eq!(untagged(&tags), 2);
    }
}
