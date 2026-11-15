//! A PLINK genetic map (`chromosome id cM position`, the layout of Beagle's
//! maps in `maps/chrN.map`): the genetic position of a base position, and
//! back.
//!
//! Between two markers the map is linear; before the first marker and after
//! the last it is flat, at the end marker's cM.

use crate::{Error, Result, cannot_read, repr};
use std::path::Path;

/// Markers sorted by position, with cM never going down.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct GeneticMap {
    pos: Vec<i64>,
    cm: Vec<f64>,
}

impl GeneticMap {
    /// Reads a PLINK `.map`: the third column is cM, the fourth the
    /// position, split on whitespace; a line with fewer than 4 columns is
    /// skipped. gzip is read through.
    ///
    /// A position that is not a whole number, a cM that is not a finite
    /// number, or a cM that goes down from one position to the next is an
    /// error: a map like that cannot place a base or a crossover.
    pub fn read(path: &Path) -> Result<GeneticMap> {
        let input = og_io::open(path).map_err(|e| cannot_read(path, &e))?;
        let mut lines = og_io::Lines::new(input, og_io::Utf8::Strict);
        let mut line = String::new();
        let mut markers = Vec::new();
        let mut n = 0u64;
        while lines
            .read_line(&mut line)
            .map_err(|e| cannot_read(path, &e))?
        {
            n += 1;
            let text = line.strip_suffix('\n').unwrap_or(&line);
            let c: Vec<&str> = text
                .split(pycompat::num::py_isspace)
                .filter(|w| !w.is_empty())
                .collect();
            if c.len() < 4 {
                continue;
            }
            let bad = |what: &str, v: &str| {
                Error(format!("{}: line {n}: {what}: {}", path.display(), repr(v)))
            };
            let pos = pycompat::num::parse_int(c[3])
                .ok()
                .and_then(|p| i64::try_from(p).ok())
                .ok_or_else(|| bad("position is not a whole number", c[3]))?;
            let cm =
                pycompat::num::parse_float(c[2]).map_err(|_| bad("cM is not a number", c[2]))?;
            if !cm.is_finite() {
                return Err(bad("cM is not a finite number", c[2]));
            }
            markers.push((pos, cm));
        }
        GeneticMap::from_markers(markers).map_err(|e| Error(format!("{}: {e}", path.display())))
    }

    /// A map from (position, cM) markers in any order. The cM must be
    /// finite and must not go down as the position goes up.
    pub fn from_markers(mut markers: Vec<(i64, f64)>) -> Result<GeneticMap> {
        if let Some(&(p, c)) = markers.iter().find(|(_, c)| !c.is_finite()) {
            return Err(Error(format!("cM {c} at position {p} is not finite")));
        }
        markers.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
        if let Some(w) = markers.windows(2).find(|w| w[1].1 < w[0].1) {
            return Err(Error(format!(
                "cM goes down from {} at position {} to {} at position {}",
                w[0].1, w[0].0, w[1].1, w[1].0
            )));
        }
        let (pos, cm) = markers.into_iter().unzip();
        Ok(GeneticMap { pos, cm })
    }

    /// The number of markers.
    pub fn len(&self) -> usize {
        self.pos.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pos.is_empty()
    }

    /// The markers, (position, cM), in position order.
    pub fn markers(&self) -> impl Iterator<Item = (i64, f64)> + '_ {
        self.pos.iter().copied().zip(self.cm.iter().copied())
    }

    /// The genetic position of `pos` in cM: linear between the markers
    /// around it, flat beyond the ends, 0 on an empty map. At a position
    /// with two markers, the last one counts.
    pub fn cm_at(&self, pos: i64) -> f64 {
        let (Some(&p_first), Some(&p_last)) = (self.pos.first(), self.pos.last()) else {
            return 0.0;
        };
        if pos < p_first {
            return self.cm[0];
        }
        if pos >= p_last {
            return self.cm[self.cm.len() - 1];
        }
        // The last marker at or before pos, and the one after it.
        let hi = self.pos.partition_point(|&p| p <= pos);
        let (p0, c0, p1, c1) = (self.pos[hi - 1], self.cm[hi - 1], self.pos[hi], self.cm[hi]);
        if p1 > p0 {
            c0 + (c1 - c0) * (pos - p0) as f64 / (p1 - p0) as f64
        } else {
            c0
        }
    }

    /// The base position (not rounded) at genetic position `cm`: the
    /// inverse of [`cm_at`](GeneticMap::cm_at) where the map rises, the
    /// first position of a flat stretch, the end positions beyond the ends,
    /// 0 on an empty map.
    pub fn pos_at(&self, cm: f64) -> f64 {
        let (Some(&c_first), Some(&c_last)) = (self.cm.first(), self.cm.last()) else {
            return 0.0;
        };
        if cm <= c_first {
            return self.pos[0] as f64;
        }
        if cm > c_last {
            return self.pos[self.pos.len() - 1] as f64;
        }
        // The first marker at or above cm; it exists, and is not the first.
        let hi = self.cm.partition_point(|&c| c < cm);
        if self.cm[hi] == cm {
            return self.pos[hi] as f64;
        }
        // Here cm[hi - 1] < cm < cm[hi], so the stretch rises.
        let (p0, c0, p1, c1) = (self.pos[hi - 1], self.cm[hi - 1], self.pos[hi], self.cm[hi]);
        p0 as f64 + (p1 - p0) as f64 * (cm - c0) / (c1 - c0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::Dir;

    /// Three markers: 1 cM per Mb, then 3 cM per Mb, written out of order.
    const MAP: &str = "chr1 . 3.0 2000000\nchr1 . 0.0 1000000\n\nchr1 . 6.0 3000000\n";

    /// The same markers with the cM column shuffled across them: what a
    /// map read from the wrong column or the wrong file looks like.
    const SHUFFLED: &str = "chr1 . 6.0 2000000\nchr1 . 3.0 1000000\n\nchr1 . 0.0 3000000\n";

    #[test]
    fn interpolates_hand_values_at_both_ends_and_between() {
        let d = Dir::new("map");
        let text = if crate::PLANTED { SHUFFLED } else { MAP };
        let map = GeneticMap::read(&d.put("chr1.map", text)).unwrap();
        assert_eq!(map.len(), 3);
        // Before the first marker, at it, between, at the last, after it.
        assert_eq!(map.cm_at(1), 0.0);
        assert_eq!(map.cm_at(1_000_000), 0.0);
        assert_eq!(map.cm_at(1_500_000), 1.5);
        assert_eq!(map.cm_at(2_000_000), 3.0);
        assert_eq!(map.cm_at(2_250_000), 3.75);
        assert_eq!(map.cm_at(3_000_000), 6.0);
        assert_eq!(map.cm_at(250_000_000), 6.0);
        // And back.
        assert_eq!(map.pos_at(-1.0), 1_000_000.0);
        assert_eq!(map.pos_at(0.0), 1_000_000.0);
        assert_eq!(map.pos_at(1.5), 1_500_000.0);
        assert_eq!(map.pos_at(3.75), 2_250_000.0);
        assert_eq!(map.pos_at(6.0), 3_000_000.0);
        assert_eq!(map.pos_at(60.0), 3_000_000.0);
        for pos in [1_000_000, 1_234_567, 2_000_001, 2_999_999] {
            assert!((map.pos_at(map.cm_at(pos)) - pos as f64).abs() < 1e-6);
        }
    }

    #[test]
    fn a_shuffled_cm_column_is_refused() {
        let d = Dir::new("map-shuffled");
        let err = GeneticMap::read(&d.put("chr1.map", SHUFFLED)).unwrap_err();
        assert!(
            err.0.ends_with(
                "chr1.map: cM goes down from 6 at position 2000000 to 0 at position 3000000"
            ),
            "{err}"
        );
        // Lines out of order are only sorted: the map is the same.
        let a = GeneticMap::read(&d.put("a.map", MAP)).unwrap();
        let sorted = "chr1 . 0.0 1000000\nchr1 . 3.0 2000000\nchr1 . 6.0 3000000\n";
        let b = GeneticMap::read(&d.put("b.map", sorted)).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn flat_stretches_and_empty_maps() {
        let map =
            GeneticMap::from_markers(vec![(100, 1.0), (200, 1.0), (300, 2.0), (300, 2.5)]).unwrap();
        assert_eq!(map.cm_at(150), 1.0);
        // Two markers at one position: the higher cM counts from there on.
        assert_eq!(map.cm_at(300), 2.5);
        // A flat stretch maps back to its first position.
        assert_eq!(map.pos_at(1.0), 100.0);
        assert_eq!(map.pos_at(1.5), 250.0);
        // A first position with two markers: the last one counts there too.
        let map = GeneticMap::from_markers(vec![(100, 1.0), (100, 2.0), (200, 3.0)]).unwrap();
        assert_eq!(map.cm_at(99), 1.0);
        assert_eq!(map.cm_at(100), 2.0);
        assert_eq!(map.cm_at(150), 2.5);
        // An inner and a trailing flat stretch map back to their first
        // position; a cM between markers still interpolates.
        let map =
            GeneticMap::from_markers(vec![(0, 0.0), (10, 1.0), (20, 1.0), (30, 2.0), (40, 2.0)])
                .unwrap();
        assert_eq!(map.pos_at(1.0), 10.0);
        assert_eq!(map.pos_at(2.0), 30.0);
        assert_eq!(map.pos_at(0.5), 5.0);
        assert_eq!(map.pos_at(1.5), 25.0);
        assert_eq!(map.pos_at(3.0), 40.0);
        let empty = GeneticMap::default();
        assert!(empty.is_empty());
        assert_eq!(empty.cm_at(5), 0.0);
        assert_eq!(empty.pos_at(5.0), 0.0);
    }

    #[test]
    fn bad_values_name_the_line() {
        let d = Dir::new("map-bad");
        let p = d.put("x.map", "chr1 . 0.0 10\nchr1 . 1.0 2x\n");
        let err = GeneticMap::read(&p).unwrap_err().0;
        assert!(
            err.ends_with("x.map: line 2: position is not a whole number: '2x'"),
            "{err}"
        );
        let p = d.put("y.map", "chr1 . one 10\n");
        let err = GeneticMap::read(&p).unwrap_err().0;
        assert!(
            err.ends_with("y.map: line 1: cM is not a number: 'one'"),
            "{err}"
        );
        let p = d.put("z.map", "chr1 . nan 10\n");
        let err = GeneticMap::read(&p).unwrap_err().0;
        assert!(
            err.ends_with("z.map: line 1: cM is not a finite number: 'nan'"),
            "{err}"
        );
        let err = GeneticMap::read(&d.0.join("missing.map")).unwrap_err().0;
        assert!(err.starts_with("cannot read "), "{err}");
    }
}
